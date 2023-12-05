import json
import os
from random import randint
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import csv
import cnapy.core
from cnapy.core import efm_computation, Scenario
import itertools
from collections import defaultdict
from typing import Dict, Tuple, List
from collections import Counter
import gurobipy
import numpy as np
import cobra
from cobra.util.array import create_stoichiometric_matrix
from cobra.core.dictlist import DictList
from optlang.symbolics import Zero, Add
from qtpy.QtWidgets import QMessageBox
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis




def save_efm_outputs(efm_results, output_file, cobra_model):
    efm_data = []
    efms = efm_results[0]
    reac_id = np.array(efms.reac_id)
    fv_mat = efms.fv_mat
    for i in range(reac_id.shape[0]):
        flux = fv_mat[:, i]
        flux = flux[0]
        reaction = cobra_model.reactions.get_by_id(reac_id[i])
        lb = reaction.lower_bound
        ub = reaction.upper_bound
        efm_data.append([reac_id[i], reaction.id, "", flux, lb, ub])
    df = pd.DataFrame(efm_data, columns=['Id', 'Name', 'Scenario', 'Flux', 'LB', 'UB'])
    df.to_csv(output_file, index=False)


def save_outputs(cobra_py_model, output_file):
    cobra_py_model.optimize()
    fva_data = []
    for reaction in cobra_py_model.reactions:
        flux_range = cobra.flux_analysis.flux_variability_analysis(cobra_py_model, reaction, fraction_of_optimum=0.9)
        flux_min = flux_range.minimum
        flux_max = flux_range.maximum
        fva_data.append([reaction.id, reaction.name, '', flux_min, flux_max])

    df = pd.DataFrame(fva_data, columns=['Id', 'Name', 'Scenario', 'Flux', 'LB', 'UB'])
    df.to_csv(output_file, index=False)


def save_fva_outputs(fva_results, output_file):
    fva_data = []
    for reaction_id, reaction_data in fva_results.iterrows():
        flux = reaction_data['maximum'] if reaction_data['maximum'] == reaction_data['minimum'] else ''
        lb = reaction_data['minimum']
        ub = reaction_data['maximum']
        fva_data.append([reaction_id, reaction_id, '', flux, lb, ub])

    df = pd.DataFrame(fva_data, columns=['Id', 'Name', 'Scenario', 'Flux', 'LB', 'UB'])
    df.to_csv(output_file, index=False)


def load_scenario_from_csv(file_path):
    scenario = Scenario()
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            reaction_id = row['Reaction']
            lower_bound = float(row['LowerBound'])
            upper_bound = float(row['UpperBound'])
            scenario[reaction_id] = (lower_bound, upper_bound)
    print(scenario)  # Add this line for debugging
    return scenario


def save_sum_results(result, output_file):
    with open(output_file, 'w') as file:
        file.write(str(result))


def run():
    # Check if 'model.xml' exists, and if so, read it
    # Check if 'model.xml' exists, and if so, read it
    if os.path.isfile('out/model.xml'):
        model_file = 'out/model.xml'
        model = cobra.io.read_sbml_model(model_file)
    # If 'model.xml' doesn't exist, check if 'model.sbml' exists, and if so, read it
    elif os.path.isfile('out/model.sbml'):
        model_file = 'out/model.sbml'
        model = cobra.io.read_sbml_model(model_file)
    else:
        raise FileNotFoundError("No model file (model.xml or model.sbml) found")



    # Load scenario data from CSV
    scenario_file = 'out/scenario.csv'
    scen_values = load_scenario_from_csv(scenario_file)

    # Apply the scenario data to the model
    for reaction_id, (lb, ub) in scen_values.items():
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.lower_bound = lb
        reaction.upper_bound = ub

    # Perform EFM analysis using the loaded scenario values
    efm_results = cnapy.core.efm_computation(model, {}, True)

    # Part 1: Update nodes and edges using the latest efm_results
    nodes = []
    edges = []

    # Add nodes for metabolites
    for metabolite in model.metabolites:
        node = {
            'id': metabolite.id,
            'label': metabolite.id,
            'title': metabolite.name,
            'color': 'lightblue',
            'shape': 'dot'
        }
        nodes.append(node)

    # Add nodes for reactions with flux and gene information
    for i, reaction in enumerate(model.reactions):
        flux_value = efm_results[0].fv_mat[:, i][0]
        gene_info = model.reactions.get_by_id(reaction.id).gene_reaction_rule

        # Get reaction bounds
        lower_bound = reaction.lower_bound
        upper_bound = reaction.upper_bound

        node_label = f"Reaction: {reaction.id}<br>Flux: {flux_value:.2f}<br>"
        node_label += f"LowerBound: {lower_bound}<br>"
        node_label += f"UpperBound: {upper_bound}"

        node = {
            'id': reaction.id,
            'label': reaction.id,
            'title': node_label,
            'color': 'lightgreen',
            'shape': 'box'
        }
        nodes.append(node)

        # Add edges for metabolite-reaction connections
        for metabolite, coefficient in reaction.metabolites.items():
            edge = {
                'from': metabolite.id,
                'to': reaction.id,
                'title': f"Coefficient: {coefficient:.2f}"
            }
            edges.append(edge)

    # Write nodes to nodes.csv file
    with open('CNApy_outputs/nodes.csv', 'w', newline='') as nodes_file:
        fieldnames = ['id', 'label', 'title', 'color', 'shape']
        writer = csv.DictWriter(nodes_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(nodes)

    # Write edges to edges.csv file
    with open('CNApy_outputs/edges.csv', 'w', newline='') as edges_file:
        fieldnames = ['from', 'to', 'title']
        writer = csv.DictWriter(edges_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(edges)

    print("CNApy_outputs/nodes.csv and edges.csv files have been created successfully.")

    # Save EFM analysis results
    save_efm_outputs(efm_results, 'CNApy_outputs/efm_flux_results.csv', model)


    optimize_efm = model.optimize()
    summary_efm = model.summary()

    # Save summary_fva as text files
    save_sum_results(summary_efm, 'CNApy_outputs/summary_efm.txt')

    # Save optimize_fva as a CSV file
    optimize_efm_df = pd.DataFrame({'Reaction': optimize_efm.fluxes.index,
                                    'Flux': optimize_efm.fluxes.values,
                                    'Reduced Costs': optimize_efm.reduced_costs.values})
    optimize_efm_df.to_csv('CNApy_outputs/optimize_efm.csv', index=False)

    

    # Part 2: FVA analysis
    fva_scenario = load_scenario_from_csv('out/scenario.csv')
    appdata = initialize_cna()
    appdata.scen_values_set_multiple(list(fva_scenario.keys()), list(fva_scenario.values()))
    appdata.project.load_scenario_into_model(model, appdata)

    fva_results = cobra.flux_analysis.flux_variability_analysis(model, model.reactions)
    save_fva_outputs(fva_results, 'CNApy_outputs/fva_results.csv')

    optimize_fva = model.optimize()
    summary_fva = model.summary()

    # Save summary_fva as text files
    save_sum_results(summary_fva, 'CNApy_outputs/summary_fva.txt')

    # Save optimize_fva as a CSV file
    optimize_fva_df = pd.DataFrame({'Reaction': optimize_fva.fluxes.index,
                                    'Flux': optimize_fva.fluxes.values,
                                    'Reduced Costs': optimize_fva.reduced_costs.values})
    optimize_fva_df.to_csv('CNApy_outputs/optimize_fva.csv', index=False)

    # Part 3: pFBA analysis
    pfba_scenario = load_scenario_from_csv('out/scenario.csv')
    appdata.scen_values_set_multiple(list(pfba_scenario.keys()), list(pfba_scenario.values()))
    appdata.project.load_scenario_into_model(model, appdata)

    pfba_solution = cobra.flux_analysis.pfba(model)
    pfba_df = pd.DataFrame(pfba_solution.fluxes.items(), columns=['Reaction', 'Flux'])
    pfba_df.to_csv('CNApy_outputs/pfba_results.csv', index=False)

    optimize_pfba = model.optimize()
    summary_pfba = model.summary()

    # Save summary_fva as text files
    save_sum_results(summary_pfba, 'CNApy_outputs/summary_pfba.txt')

    # Save optimize_fva as a CSV file
    optimize_pfba_df = pd.DataFrame({'Reaction': optimize_pfba.fluxes.index,
                                    'Flux': optimize_pfba.fluxes.values,
                                    'Reduced Costs': optimize_pfba.reduced_costs.values})
    optimize_pfba_df.to_csv('CNApy_outputs/optimize_pfba.csv', index=False)

    # Part 4: FBA analysis
    fba_scenario = load_scenario_from_csv('out/scenario.csv')
    appdata.scen_values_set_multiple(list(fba_scenario.keys()), list(fba_scenario.values()))
    appdata.project.load_scenario_into_model(model, appdata)

    FBA_result = model.optimize()
    fba_df = pd.DataFrame({'Reaction': FBA_result.fluxes.index, 'Flux': FBA_result.fluxes.values})
    fba_df.to_csv('CNApy_outputs/fba_results.csv', index=False)

    optimize_fba = model.optimize()
    summary_fba = model.summary()

    # Save summary_fva as text files
    save_sum_results(summary_fba, 'CNApy_outputs/summary_fba.txt')

    # Save optimize_fva as a CSV file
    optimize_fba_df = pd.DataFrame({'Reaction': optimize_fba.fluxes.index,
                                    'Flux': optimize_fba.fluxes.values,
                                    'Reduced Costs': optimize_fba.reduced_costs.values})
    optimize_fba_df.to_csv('CNApy_outputs/optimize_fba.csv', index=False)

    # network data

    

def initialize_cna():
    class CNAData:
        def __init__(self):
            self.project = initialize_project()

        def scen_values_set_multiple(self, keys, values):
            self.scen_values = dict(zip(keys, values))

    return CNAData()


def initialize_project():
    class Project:
        def load_scenario_into_model(self, model, appdata):
            for reaction_id, (lb, ub) in appdata.scen_values.items():
                reaction = model.reactions.get_by_id(reaction_id)
                reaction.lower_bound = lb
                reaction.upper_bound = ub

    return Project()


# Run the code
run()

