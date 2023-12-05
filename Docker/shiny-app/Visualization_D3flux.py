import csv
from cobra.core import Metabolite, Reaction, Model
from cobra.io import write_sbml_model

# Create a new model
model = Model('simple_model')

# Initialize empty dictionaries for metabolites and reactions
metabolites = {}
reactions = {}

# Read the CSV file
with open('out/metabolites.csv', mode='r', newline='', encoding='utf-8') as file:
    reader = csv.DictReader(file)
    for row in reader:
        if row['type'] == 'metabolite':
            # Create metabolite with a default compartment and add to the metabolites dictionary
            metabolites[row['name']] = Metabolite(row['name'], compartment='c')
        elif row['type'] == 'reaction':
            # Create reaction and add to the reactions dictionary
            reactions[row['name']] = Reaction(row['name'])

# Add the metabolites and reactions to the model
model.add_metabolites(metabolites.values())
model.add_reactions(reactions.values())

# Re-read the CSV file to build the reaction from string as the metabolites should be added first
with open('out/metabolites.csv', mode='r', newline='', encoding='utf-8') as file:
    reader = csv.DictReader(file)
    for row in reader:
        if row['type'] == 'reaction':
            # Build reaction from string
            model.reactions.get_by_id(row['name']).build_reaction_from_string(row['equation'])

# Save the model as an SBML file
write_sbml_model(model, "out/model.xml")

from d3flux import flux_map
html = flux_map(model, display_name_format=lambda x: str(x.id),
         flux_dict={rxn.id: None for rxn in model.reactions}, figsize=(1080, 1080))

from mfo import render_html_to_image
raw_html = '<script src="https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.6/require.min.js"></script>' + html.__html__()
with open("output.html", "w") as f:
    f.write(raw_html)
render_html_to_image()
