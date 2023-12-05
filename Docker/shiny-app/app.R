#########################################
##### Cat-E Bioinformatics WepApp #######
#### Written by Rana Salihoglu       ####
#### Update:14_Nov_2023              ####
#### Ⓒ                              ####
#########################################

library("gatom")
library("mwcsr")
library("shiny")
library("bs4Dash")
library("thematic")
library("waiter")
library("dplyr")
library("billboarder")
library("readxl")
library("shinyWidgets")
library("data.table")
library("visNetwork")
library("png")
library("plotly")
library("biomaRt")
#library("TCGAbiolinks")
library("shinycssloaders")
library("ggplot2")
library("ggplotify")
library("DT")
library("RColorBrewer")
library("survival")
library("survminer")
library("limma")
library("edgeR")
library("reshape2")
library("maftools")
library("g3viz")
library("genefilter")
library("SummarizedExperiment")
library("imsig")
library("DiffCorr")
library("igraph")
library("reticulate")
library("shinyjs")
library("htmlwidgets")
library("shinyCyJS")
library("NGLVieweR")
library("httr")
library("colourpicker")
require("topGO")
require("org.Hs.eg.db")
library("shinyBS")
library("rintrojs")
library("stringr")
library("ggvenn")
library("clusterProfiler")
library("pathview")
library("GEOquery")
library("heatmaply")


source('setShadow.R')
source('build-jimena-app.R')
options(shiny.maxRequestSize = 90*1024^2)

# Import Python modules
rdkit <- import("rdkit")
urllib <- import("urllib.request")
json <- import("json")



thematic_shiny()

# toast options
toastOpts <- list(
  autohide = TRUE,
  icon = "fas fa-home",
  close = FALSE,
  position = "bottomRight"
)





ui = dashboardPage(
  
  
  preloader = list(html = tagList(spin_2(), "Loading ..."), color = "#34495E"),
  #dark = TRUE,
  #help = TRUE,
  fullscreen = TRUE,
  scrollToTop = TRUE,
  dashboardHeader(
    tags$style(HTML("
    .header-logo-container {
      display: flex;
      align-items: center;
      justify-content: center;
    }
  ")),
    
    title = tags$div(
      class = "header-logo-container",
      tags$img(src = "logo.png",  # Adjust the path as needed
               class = "brand-image img-square elevation-3",
               style = "height: 75px; width: 140px;"
      )
    ),
    
    titleWidth = 300,  # Adjust the width as needed
    
    leftUi = tagList(
      
    )
  )
  
  
  ,
  
  tags$head(
    tags$script(
      HTML(
        '
        function updateBounds(input) {
            var reactionId = input.id.replace(\'_lower_bound\', \'\').replace(\'_upper_bound\', \'\');
            var inputType = input.id.endsWith(\'_lower_bound\') ? \'lower_bound\' : \'upper_bound\';
            var boundValue = parseFloat(input.value);
            var reactionNode = network.findNode(reactionId);
            if (reactionNode !== null) {
                reactionNode.title = reactionNode.title.replace(new RegExp(inputType.charAt(0).toUpperCase() + inputType.slice(1) + \': (-?\\d*\\.?\\d+)\'), inputType.charAt(0).toUpperCase() + inputType.slice(1) + \': \' + boundValue);
                network.stabilize();
            }
            form_values[reactionId][inputType] = boundValue;
            console.log(form_values);
        }
        '
      )
    )
  ),
  
  
  sidebar = dashboardSidebar(
    introjsUI(),
    fixed = TRUE,
    skin = "light",
    status = "primary",
    id = "sidebar",
    setShadow(id = "home_human"),
    setShadow(id = "home_donut"),
    setShadow(id = "cancers"),
    setShadow(id = "analyTab"),
    setShadow(id = "plotS"),
    setShadow(id = "models"),
    setShadow(id = "pertbox"),
    setShadow(id = "simjim"),
    setShadow(id = "jimedit"),
    setShadow(id = "filtered_onc"),
    setShadow(id = "PPI_onc"),
    setShadow(id = "Enrich"),
    setShadow(id = "Enrich2"),
    setShadow(id = "metabolic1"),
    setShadow(id = "metabolic2"),
    setShadow(id = "metabolic3"),
    setShadow(id = "flux1"),
    setShadow(id = "flux3"),
    setShadow(id = "flux4"),
    setShadow(id = "SNV"),
    setShadow(id = "SNVbox1"),
    setShadow(id = "SNVbox2"),
    setShadow(id = "infcell"),
    setShadow(id = "celltypeS"),
    setShadow(id = "cellMarkers"),
    setShadow(id = "Jimena_table"),
    setShadow(id = "drugTissue"),
    setShadow(id = "drugDesc"),
    setShadow(id = "drugStruc"),
    setShadow(id = "drugDesTab"),
    setShadow(id = "drugDesTab2"),
    setShadow(id = "drugDesTab3"),
    setShadow(id = "drugGenes"),
    setShadow(id = "members"),
    setShadow(id = "linksForHotlist"),
    setShadow(id = "Imsig"),
    setShadow(id = "ImsigBar"),
    setShadow(id = "ImsigRes"),
    setShadow(id = "Survival"),
    setShadow(id = "Survival2"),
    setShadow(id = "CarTarget"),
    setShadow(id = "BispecificA"),
    setShadow(id = "cytostaticBox"),
    setShadow(id = "ovirusTherapyBox"),
    setShadow(id = "check1"),
    setShadow(id = "check2"),
    
    sidebarMenu(
      id = "tabs",
      flat = FALSE,
      compact = FALSE,
      childIndent = TRUE,
      menuItem(
        "Home",
        tabName = "tab_home",
        icon = icon("home")
      ),
      sidebarHeader("Cards"),
      menuItem(
        "Oncolytic virus",
        tabName = "v_species",
        icon = icon("viruses")
      ),
      menuItem(
        "Analyses",
        tabName = "cancer_type",
        icon = icon("disease"),
        startExpanded = FALSE,
        menuSubItem(
          text = HTML(
            paste(
              "Expression",
              dashboardBadge(
                # "new",
                position = "right",
                color = "danger"
              )
            )
          ),
          tabName = "DeG",
          icon = icon("circle")
        ),
        menuSubItem(
          text = HTML(
            paste(
              "Metabolic",
              dashboardBadge(
                # "new",
                position = "right",
                color = "danger"
              )
            )
          ),
          tabName = "metabolic_tab",
          icon = icon("circle")
        ),
        menuSubItem(
          text = HTML(
            paste(
              "Flux",
              dashboardBadge(
                "new",
                position = "right",
                color = "danger"
              )
            )
          ),
          tabName = "flux_tab",
          icon = icon("circle")
        ),
        menuSubItem(
          text = HTML(
            paste(
              "Enrichment",
              dashboardBadge(
                #  "!",
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "analy2",
          icon = icon("circle")
        ),
        
        
        menuSubItem(
          text = HTML(
            paste(
              "Survival",
              dashboardBadge(
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "analy3",
          icon = icon("circle")
        ),
        
        menuSubItem(
          text = HTML(
            paste(
              "Immune Signature",
              dashboardBadge(
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "analy6",
          icon = icon("circle")
        ),
        menuSubItem(
          text = HTML(
            paste(
              "Single Nucleotide Variation",
              dashboardBadge(
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "analy4",
          icon = icon("circle")
        )
        # menuSubItem(
        #    text = HTML(
        #     paste(
        #      "Copy Number Variation",
        #      dashboardBadge(
        #        position = "right",
        #        color = "success"
        #      )
        #    )
        #  ),
        #  tabName = "analy5"
        #)
      ),
      
      menuItem(
        "Cell line",
        tabName = "cell_line",
        icon = icon("syringe")
      ),
      menuItem(
        "Cell type",
        tabName = "cell_type",
        icon = icon("syringe")
      ),
      menuItem(
        "Jimena",
        tabName = "cancerModel",
        icon = icon("frog")
      ),
      menuItem(
        "Protein Structure",
        tabName = "strucPro",
        icon = icon("dna")
      ),
      
      sidebarHeader("Treatment"),
      menuItem(
        text = "Drug",
        icon = icon("pills"),
        startExpanded = FALSE,
        menuSubItem(
          text = HTML(
            paste(
              "Description",
              dashboardBadge(
                #"new",
                position = "right",
                color = "danger"
              )
            )
          ),
          tabName = "drug2"
          #icon = icon("circle")
        ),
        menuSubItem(
          text = HTML(
            paste(
              "Tissues",
              dashboardBadge(
                # "new",
                position = "right",
                color = "danger"
              )
            )
          ),
          tabName = "drug1"
          # icon = icon("circle")
        ),
        menuSubItem(
          text = HTML(
            paste(
              "Genes",
              dashboardBadge(
                #"!",
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "geneDrug",
          
        ),
        
        menuSubItem(
          text = HTML(
            paste(
              "STITCH",
              dashboardBadge(
                #"!",
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "stitch",
          
        )
      ),
      
      menuItem(
        text = "Immune Modulatory",
        icon = icon("pills"),
        startExpanded = FALSE,
        
        menuSubItem(
          text = HTML(
            paste(
              "Bispesific Antibody",
              dashboardBadge(
                # "!",
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "bispecificTab"
        ),
        
        menuSubItem(
          text = HTML(
            paste(
              "CAR-T cell Therapy",
              dashboardBadge(
                "new",
                position = "right",
                color = "danger"
              )
            )
          ),
          tabName = "carTab"
          #icon = icon("circle")
        ),
        menuSubItem(
          text = HTML(
            paste(
              "Checkpoint",
              dashboardBadge(
                # "!",
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "checkPoinTab"
        ),
        
        
        menuSubItem(
          text = HTML(
            paste(
              "Cytostatic Therapy",
              dashboardBadge(
                # "!",
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "cytostaticTab"
        ),
        
        menuSubItem(
          text = HTML(
            paste(
              "Oncolytic virus Therapy",
              dashboardBadge(
                # "!",
                position = "right",
                color = "success"
              )
            )
          ),
          tabName = "ovirusTherapyTab"
        )
        
      ),
      
      sidebarHeader("Data sources"),
      menuItem("Database",
               tabName = "DataSour",
               icon = icon("database")
      ),
      menuItem("Hotlist",
               tabName = "hotList",
               icon = icon("database")
      ),
      menuItem("About us",
               tabName = "usAbout",
               icon = icon("people-group")
      ),
      actionBttn(
        inputId = "help3",
        label = "Help",
        style = "unite", 
        color = "success",
        icon = icon("info-circle", verify_fa = FALSE)
      )
      
    )
  ),
  body = dashboardBody(
    tabItems(
      ###### home tab UI #################
      tabItem(
        introjsUI(),
        useShinyjs(),
        
        tabName = "tab_home",
        fluidRow(
          column(
            width = 10,  # Adjust the width as needed
            offset = 1, # Adjust the offset to center the box
            box(
              id = "home_donut",
              title = "", 
              width = 10,
              solidHeader = TRUE,
              status = "gray-dark", 
              closable = FALSE,
              maximizable = TRUE, 
              collapsible = TRUE,
              billboarderOutput("pie", height = 700)
            )),
          column(width = 1,
                 uiOutput("home1"),
                 uiOutput("home2"),
                 uiOutput("home3")
          ),
          
          column(
            width = 10,  # Adjust the width as needed
            offset = 1, # Adjust the offset to center the box
            box(
              id = "home_human",
              title = "Cases", 
              width = 10,
              solidHeader = TRUE,
              status = "gray-dark", 
              closable = FALSE,
              maximizable = TRUE, 
              collapsible = TRUE,
              #imageOutput("humanBody", width = 700, height = 700),
              plotlyOutput("plotCase", height = 700)
            )
          ),
          column(width = 1,
                 uiOutput("home4"),
                 uiOutput("home5")
          ))
        
        
      ),
      ################## oncolytic virus tab UI ##############################    
      
      tabItem(
        tabName = "v_species",
        
        box(
          id = "filtered_onc",
          title = "Filtered by selected oncolytic virus", 
          closable = TRUE, 
          maximizable = TRUE,
          width = 12,
          status = "warning", 
          solidHeader = FALSE, 
          collapsible = TRUE,
          selectInput(
            inputId = "onc_spec", 
            label = "Oncolytic viruses", 
            choices = c("",
                        "Adenovirus",
                        "Herpes simplex virus",
                        "Measles virus",
                        "Canine parvovirus",
                        "Vaccinia virus",
                        "Alphavirus",
                        "Reovirus",
                        "Avian influenza A virus",
                        "Sindbis virus",
                        "Vesicular stomatitis virus",
                        "Parovirus",
                        "Newcastle disease virus",
                        "Enterovirus",
                        "Respiratory synctial virus",
                        "Mumps virus",
                        "Measles and mumps virus",
                        "Tanapoxvirus",
                        "Poxvirus",
                        "Sendai virus",
                        "Semliki forest virus",
                        "Maraba virus",
                        "Myxoma virus",
                        "Bovine herpesvirus",
                        "Avian orthoreovirus"
                        
            ), selected = 1),
          
          DT::dataTableOutput("virus")
        ),
        
        box(
          id = "PPI_onc",
          title = "Human and Oncolytic virus Protein-Protein Interaction", 
          closable = TRUE, 
          maximizable = TRUE,
          width = 12,
          status = "warning", 
          solidHeader = FALSE, 
          collapsible = TRUE,
          selectInput(
            inputId = "virusID", 
            label = "Oncolytic virus species", 
            choices = c("",
                        "Human adenovirus A serotype 12 (HAdV-12)"= "28282",
                        "Human adenovirus A serotype 31"="10529",
                        "Human adenovirus B serotype 3"="45659",
                        "Human adenovirus C serotype 2 (HAdV-2)"="10515",
                        "Human adenovirus C serotype 5 (HAdV-5)"="28285",
                        "Human adenovirus D serotype 9"="10527",
                        "Human adenovirus E"="130308",
                        "Human adenovirus 36" = "46936",
                        "Human herpesvirus 1" ="10298",
                        "Human herpesvirus 1 (strain 17)" ="10299",
                        "Human herpesvirus 8" ="868565",
                        "Human herpesvirus 8 (HHV-8)" = "37296",
                        "Human herpesvirus 8 type M" = "435895",
                        "Human herpesvirus" ="10335,10304,10315,10310,154633,10306,10370,10301,10308,10369,36351",
                        "Murine adenovirus A serotype 1 (MAdV-1)"="10530",
                        #"Parovirus"=" 10788",
                        "Newcastle disease virus" = "652953",
                        #"Enterovirus"="39054",
                        "Respiratory synctial virus"="11259,11250,11260,79692",
                        "Mumps virus"="11167,11168,11171,11173",
                        "Measles virus"="645098,11234,36408,132487,11235,11236,11239",
                        "Tanapoxvirus"="132475",
                        "Poxvirus" = "10243,32606,10269",
                        "Reovirus type 1 (strain Lang) (T1L)" ="10884",
                        "Reovirus type 3 (strain Dearing)" ="10886",
                        "Sendai virus"="11195,11198",
                        "Semliki forest virus"="11033",
                        "Sindbis virus" = "11034,31699",
                        "Vesicular stomatitis virus" = "11276,434490,434489,434488,11278,11284,11285",
                        #"Maraba virus",
                        "Myxoma virus"="31530",
                        "Bovine herpesvirus"="10320,10323,31518,10324,10385,79889",
                        "Avian orthoreovirus" = "11862,160750,269446,11866,11879",
                        #"Avian " ="642261",
                        "Vaccinia virus" = '10245',
                        "Vaccinia virus Copenhagen" = '10249',
                        "Vaccinia virus WR" = '10254',
                        "Vaccinia virus L-IPV" = '31531',
                        "Vaccinia virus Ankara" = '126794',
                        "All Vaccinia virus strain"='10245,10249,10254,31531,126794',
                        "Vaccinia virus GLV-1h68" = "1",
                        "Vaccinia virus VVΔTKΔN1L" ="10"
                        
            )),
          
          actionBttn(
            inputId = "runVisPPI",
            label = "Search", 
            style = "gradient",
            color = "danger",
            icon = icon("circle-play")
          ),
          tags$style(type="text/css",
                     ".shiny-output-error { visibility: hidden; }",
                     ".shiny-output-error:before { visibility: hidden; }"
          ),
          
          visNetworkOutput("viRusNetwork", width = "100%", height = "700px"),
          DT::dataTableOutput(outputId = "PPIdata")
        )),
      ############### cancer tab- Differential Expression Analysis UI ##########################
      
      tabItem(
        tabName = "DeG",
        fluidRow(
          box(
            id = "cancers",
            title = "", 
            width = 3,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            
            conditionalPanel(
              condition = "input.dataChoice == 'TCGA-GTEX data' || input.compVen",
              
            selectInput(
              inputId = "CancerType", 
              label = "Select Cancer Type Dataset", 
              choices = c("",
                          #"Adrenocortical carcinoma (ACC)"= "TCGA-ACC",
                          "Bladder Urothelial Carcinoma (BLCA)"="TCGA-BLCA",
                          "Breast invasive carcinoma (BRCA)"="TCGA-BRCA",
                          "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)"="TCGA-CESC",
                          "Cholangio carcinoma (CHOL)"="TCGA-CHOL",
                          "Colon adenocarcinoma (COAD)"="TCGA-COAD",
                          "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)"="TCGA-DLBC",
                          "Esophageal carcinoma (ESCA)"="TCGA-ESCA",
                          "Glioblastoma multiforme (GBM)"="TCGA-GBM",
                          "Head and Neck squamous cell carcinoma (HNSC)"="TCGA-HNSC",
                          "Kidney Chromophobe (KICH)"= "TCGA-KICH",
                          "Kidney renal clear cell carcinoma (KIRC)"="TCGA-KIRC",
                          "Kidney renal papillary cell carcinoma (KIRP)"="TCGA-KIRP",
                          "Acute Myeloid Leukemia (LAML)"="TCGA-LAML",
                          "Brain Lower Grade Glioma (LGG)"="TCGA-LGG",
                          "Liver hepatocellular carcinoma (LIHC)"="TCGA-LIHC",
                          "Lung adenocarcinoma (LUAD)"="TCGA-LUAD",
                          "Lung squamous cell carcinoma (LUSC)"="TCGA-LUSC",
                          "Mesothelioma (MESO)"="TCGA-MESO",
                          "Ovarian serous cystadenocarcinoma (OV)"="TCGA-OV",
                          "Pancreatic adenocarcinoma (PAAD)"="TCGA-PAAD",
                          "Pheochromocytoma and Paraganglioma (PCPG)"="TCGA-PCPG",
                          "Prostate adenocarcinoma (PRAD)"="TCGA-PRAD",
                          "Rectum adenocarcinoma (READ)"="TCGA-READ",
                          "Sarcoma (SARC)"="TCGA-SARC",
                          "Skin Cutaneous Melanoma (SKCM)"="TCGA-SKCM",
                          "Stomach adenocarcinoma (STAD)"="TCGA-STAD",
                          "Testicular Germ Cell Tumors (TGCT)"="TCGA-TGCT",
                          "Thyroid carcinoma (THCA)"= "TCGA-THCA",
                          "Thymoma (THYM)"= "TCGA-THYM",
                          "Uterine Corpus Endometrial Carcinoma (UCEC)"= "TCGA-UCEC",
                          "Uterine Carcinosarcoma (UCS)"= "TCGA-UCS"
                          #"Uveal Melanoma (UVM)"= "TCGA-UVM"
              ), selected = 1)),
            
            conditionalPanel(
              condition = "input.compVen==0",
            prettyRadioButtons(
              inputId = "dataChoice",
              label = "Choose Data Source:",
              choices = c("TCGA-GTEX data","Use own data", "Use GEO data"),
              selected = "TCGA-GTEX data",
              icon = icon("check"),
              bigger = TRUE,
              status = "info",
              animation = "jelly"
            )),
            
            
            
            conditionalPanel(
              condition = "input.dataChoice == 'Use GEO data'&& input.compVen==0",
              searchInput(
                inputId = "GEOid",
                placeholder = "Search GEO id (Example: GSE33126)",
                btnSearch = icon("search", verify_fa = FALSE), 
                btnReset = icon("trash")
              ),
              textInput("assignGroup", "Assign to group: ")),
            
            awesomeCheckbox(
              inputId = "compVen",
              label = "Create Venn Diagram", 
              value = FALSE,
              status = "info",
            ),
            
            conditionalPanel(
              condition = "input.dataChoice == 'Use own data' || input.compVen==1",
              
              fileInput(inputId = "userDif", "Upload DEG data", 
                        accept = c("text/csv","text/comma-separated-values,text/plain", ".csv")),
              
              
              pickerInput('geneColumn','Gene column', ""),
              pickerInput('logFColumn','log2FC column', ""),
              pickerInput('pAdColumn','pAdjust column', "") 
            ),
            numericInput("logFCcut", label = "|Log2FC| Cut off", value = 0),
            numericInput("qvalcut", label = "q-value Cut off",0.01,
                         min = 0, max = 1, step =0.001),
            
            
            conditionalPanel(
            condition = "input.dataChoice == 'TCGA-GTEX data' && input.compVen==0",
              
            prettyRadioButtons("AnalTy", label = "Select Differential Methods", 
                               choices = c("ANOVA"="edgeR", "Limma"="limma"),icon = icon("check"), 
                               bigger = TRUE,
                               status = "info",
                               animation = "jelly")),
            actionBttn(
              inputId = "runAnalys",
              label = "RUN", 
              style = "gradient",
              color = "danger",
              icon = icon("circle-play")
            )
            
          ),
          
          box(
            id = "analyTab",
            title = "", 
            width = 9,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            conditionalPanel(
              condition = "input.dataChoice == 'Use GEO data' && input.compVen==0",
              dataTableOutput('groupsample')),
            
            DT::dataTableOutput(outputId = "ExTab"),
            
            conditionalPanel(
              condition = "input.dataChoice == 'Use GEO data' && input.compVen==0 && input.runAnalys >0",
              
              DT::dataTableOutput(outputId = "deg")),
            
            conditionalPanel(
              condition = "input.compVen==1 && input.runAnalys >0",
              style = "display: none;",
              withSpinner(plotOutput("vennPlot", click = "vennPlot_click")),
              DTOutput("elementTable"))
          ),
          
          
          box(
            id = "plotS",
            title = "", 
            width = 12,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            conditionalPanel(
              condition = "input.dataChoice == 'TCGA-GTEX data' && input.compVen==0 && input.runAnalys >0",
            selectInput(
              inputId = "PloType", 
              label = "Select Plot Type", 
              choices = c("","Volcano plot"="plot1",
                          "Violin"="plot3",
                          "PCA Analysis Plot"="plot4"))),
            
            conditionalPanel(
              condition = "input.dataChoice == 'Use own data' && input.compVen==0 && input.runAnalys >0",
              selectInput(
                inputId = "PloType2", 
                label = "Select Plot Type", 
                choices = c("","Volcano plot"="plot1"))),
            
            conditionalPanel(
              condition = "input.dataChoice == 'Use GEO data' && input.compVen==0 && input.runAnalys >0",
              selectInput(
                inputId = "PloType3", 
                label = "Select Plot Type", 
                choices = c("","Volcano plot"="plot1",
                            "PCA plot" = "GEOpca",
                            "Heatmap" = "GEOheatmap"))
              ),
            
            conditionalPanel(
              condition = "input.runAnalys > 0 && input.compVen == 0 && (input.PloType == 'plot1' || input.PloType2 == 'plot1' || input.PloType3 == 'plot1')",
              plotlyOutput("erciyes", height = 700)
            ),
            
            conditionalPanel(
              condition = "input.PloType == 'plot3' && input.compVen==0",
              uiOutput("geneViolin"),
              plotlyOutput("violin",height = 700)
            ),
            conditionalPanel(
              condition = "input.PloType3 == 'GEOheatmap' && input.compVen==0",
              plotlyOutput("GEO_heatmap",height = 700)
            ),
            conditionalPanel(
              condition = "input.PloType3 == 'GEOpca' && input.compVen==0",
              plotlyOutput("GEO_pca",height = 700)
            ),
            conditionalPanel(
              condition = "input.PloType == 'plot4' && input.compVen==0",
              plotlyOutput("PCAnalysis",height = 700)
            )
          )
        )
      ),
      
      
      ################## Enrichment Analysis UI ########################
      tabItem(
        tabName = "analy2",
        tabsetPanel(
          tabPanel("GO",
                   br(),
        fluidRow(
          box(
            id = "Enrich",
            title = "", 
            width = 3,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            
            conditionalPanel(
              condition = "input.userinputEnrich == 0",
              
            selectInput(
              inputId = "typeForEnrc", 
              label = "Select Cancer Type Dataset", 
              choices = c("",
                          "Adrenocortical carcinoma (ACC)"= "TCGA-ACC",
                          "Bladder Urothelial Carcinoma (BLCA)"="TCGA-BLCA",
                          "Breast invasive carcinoma (BRCA)"="TCGA-BRCA",
                          "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)"="TCGA-CESC",
                          "Cholangio carcinoma (CHOL)"="TCGA-CHOL",
                          "Colon adenocarcinoma (COAD)"="TCGA-COAD",
                          "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)"="TCGA-DLBC",
                          "Esophageal carcinoma (ESCA)"="TCGA-ESCA",
                          "Glioblastoma multiforme (GBM)"="TCGA-GBM",
                          "Head and Neck squamous cell carcinoma (HNSC)"="TCGA-HNSC",
                          "Kidney Chromophobe (KICH)"= "TCGA-KICH",
                          "Kidney renal clear cell carcinoma (KIRC)"="TCGA-KIRC",
                          "Kidney renal papillary cell carcinoma (KIRP)"="TCGA-KIRP",
                          "Acute Myeloid Leukemia (LAML)"="TCGA-LAML",
                          "Brain Lower Grade Glioma (LGG)"="TCGA-LGG",
                          "Liver hepatocellular carcinoma (LIHC)"="TCGA-LIHC",
                          "Lung adenocarcinoma (LUAD)"="TCGA-LUAD",
                          "Lung squamous cell carcinoma (LUSC)"="TCGA-LUSC",
                          "Mesothelioma (MESO)"="TCGA-MESO",
                          "Ovarian serous cystadenocarcinoma (OV)"="TCGA-OV",
                          "Pancreatic adenocarcinoma (PAAD)"="TCGA-PAAD",
                          "Pheochromocytoma and Paraganglioma (PCPG)"="TCGA-PCPG",
                          "Prostate adenocarcinoma (PRAD)"="TCGA-PRAD",
                          "Rectum adenocarcinoma (READ)"="TCGA-READ",
                          "Sarcoma (SARC)"="TCGA-SARC",
                          "Skin Cutaneous Melanoma (SKCM)"="TCGA-SKCM",
                          "Stomach adenocarcinoma (STAD)"="TCGA-STAD",
                          "Testicular Germ Cell Tumors (TGCT)"="TCGA-TGCT",
                          "Thyroid carcinoma (THCA)"= "TCGA-THCA",
                          "Thymoma (THYM)"= "TCGA-THYM",
                          "Uterine Corpus Endometrial Carcinoma (UCEC)"= "TCGA-UCEC",
                          "Uterine Carcinosarcoma (UCS)"= "TCGA-UCS"
                          #"Uveal Melanoma (UVM)"= "TCGA-UVM"
              ), selected = 1)),
            
            awesomeCheckbox(
              inputId = "userinputEnrich",
              label = "Use own data", 
              value = FALSE,
              status = "info",
            ),
            
            conditionalPanel(
              condition = "input.userinputEnrich",
              
              fileInput(inputId = "userEnrich", "Upload DEG data", 
                        accept = c("text/csv","text/comma-separated-values,text/plain", ".csv")),
              pickerInput('geneColumnE','Gene column', ""),
              pickerInput('logFColumnE','log2FC column', ""),
              pickerInput('pAdColumnE','pAdjust column', "")
            ),
            
            numericInput("logFCcut", label = "|Log2FC| Cut off", value = 0),
            numericInput("qvalcut", label = "q-value Cut off",0.01,
                         min = 0, max = 1, step =0.0000001),
            actionBttn(
              inputId = "runEnrich",
              label = "RUN", 
              style = "gradient",
              color = "danger",
              icon = icon("circle-play")
            )
          ),
          box(
            id = "Enrich2",
            title = "",
            width = 9,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            selectInput("EnrichType",label = "Enrichment category",
                        choices = c("Biological Process" = "BPplot1",
                                    "Cellular Component"= "CCplot2",
                                    "Molecular Function" = "MFplot3"), selected = 1),
            conditionalPanel(
              condition = "input.runEnrich > 0",
              style = "display: none;",
              withSpinner(plotlyOutput(outputId = "whichenrich")),
            withSpinner(DT::dataTableOutput("TablEnrich")))
          )
        )),
        tabPanel("KEGG",
                 br(),
                 fluidRow(
                   # Left side box
                   box(
                     title = "Parameters",
                     width = 3,
                     status = "primary", 
                     closable = FALSE,
                     maximizable = TRUE, 
                     collapsible = TRUE,
                     prettyRadioButtons(
                       inputId = "KEGGinput",
                       label = "Choose:", 
                       choices = c("Cancer Data", "PPI Data"),
                       icon = icon("check"), 
                       bigger = TRUE,
                       status = "info",
                       animation = "jelly"
                     ),
                     selectInput(
                       inputId = "typeforKEGG", 
                       label = "Select Cancer Type Dataset",
                       choices = c("",
                                   "Adrenocortical carcinoma (ACC)"= "TCGA-ACC",
                                   "Bladder Urothelial Carcinoma (BLCA)"="TCGA-BLCA",
                                   "Breast invasive carcinoma (BRCA)"="TCGA-BRCA",
                                   "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)"="TCGA-CESC",
                                   "Cholangio carcinoma (CHOL)"="TCGA-CHOL",
                                   "Colon adenocarcinoma (COAD)"="TCGA-COAD",
                                   "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)"="TCGA-DLBC",
                                   "Esophageal carcinoma (ESCA)"="TCGA-ESCA",
                                   "Glioblastoma multiforme (GBM)"="TCGA-GBM",
                                   "Head and Neck squamous cell carcinoma (HNSC)"="TCGA-HNSC",
                                   "Kidney Chromophobe (KICH)"= "TCGA-KICH",
                                   "Kidney renal clear cell carcinoma (KIRC)"="TCGA-KIRC",
                                   "Kidney renal papillary cell carcinoma (KIRP)"="TCGA-KIRP",
                                   "Acute Myeloid Leukemia (LAML)"="TCGA-LAML",
                                   "Brain Lower Grade Glioma (LGG)"="TCGA-LGG",
                                   "Liver hepatocellular carcinoma (LIHC)"="TCGA-LIHC",
                                   "Lung adenocarcinoma (LUAD)"="TCGA-LUAD",
                                   "Lung squamous cell carcinoma (LUSC)"="TCGA-LUSC",
                                   "Mesothelioma (MESO)"="TCGA-MESO",
                                   "Ovarian serous cystadenocarcinoma (OV)"="TCGA-OV",
                                   "Pancreatic adenocarcinoma (PAAD)"="TCGA-PAAD",
                                   "Pheochromocytoma and Paraganglioma (PCPG)"="TCGA-PCPG",
                                   "Prostate adenocarcinoma (PRAD)"="TCGA-PRAD",
                                   "Rectum adenocarcinoma (READ)"="TCGA-READ",
                                   "Sarcoma (SARC)"="TCGA-SARC",
                                   "Skin Cutaneous Melanoma (SKCM)"="TCGA-SKCM",
                                   "Stomach adenocarcinoma (STAD)"="TCGA-STAD",
                                   "Testicular Germ Cell Tumors (TGCT)"="TCGA-TGCT",
                                   "Thyroid carcinoma (THCA)"= "TCGA-THCA",
                                   "Thymoma (THYM)"= "TCGA-THYM",
                                   "Uterine Corpus Endometrial Carcinoma (UCEC)"= "TCGA-UCEC",
                                   "Uterine Carcinosarcoma (UCS)"= "TCGA-UCS"
                                   #"Uveal Melanoma (UVM)"= "TCGA-UVM"
                       ), selected = 1),
                     
                     conditionalPanel(
                       condition = "input.KEGGinput == 'PPI Data'",
                       selectInput(
                         inputId = "virusID2", 
                         label = "Oncolytic virus species",
                         choices = c("",
                                     "Human adenovirus A serotype 12 (HAdV-12)"= "28282",
                                     "Human adenovirus A serotype 31"="10529",
                                     "Human adenovirus B serotype 3"="45659",
                                     "Human adenovirus C serotype 2 (HAdV-2)"="10515",
                                     "Human adenovirus C serotype 5 (HAdV-5)"="28285",
                                     "Human adenovirus D serotype 9"="10527",
                                     "Human adenovirus E"="130308",
                                     "Human adenovirus 36" = "46936",
                                     "Human herpesvirus 1" ="10298",
                                     "Human herpesvirus 1 (strain 17)" ="10299",
                                     "Human herpesvirus 8" ="868565",
                                     "Human herpesvirus 8 (HHV-8)" = "37296",
                                     "Human herpesvirus 8 type M" = "435895",
                                     "Human herpesvirus" ="10335,10304,10315,10310,154633,10306,10370,10301,10308,10369,36351",
                                     "Murine adenovirus A serotype 1 (MAdV-1)"="10530",
                                     #"Parovirus"=" 10788",
                                     "Newcastle disease virus" = "652953",
                                     "Enterovirus"="39054",
                                     "Respiratory synctial virus"="11259,11250,11260,79692",
                                     "Mumps virus"="11167,11168,11171,11173",
                                     "Measles virus"="645098,11234,36408,132487,11235,11236,11239",
                                     "Tanapoxvirus"="132475",
                                     "Poxvirus" = "10243,32606,10269",
                                     "Reovirus type 1 (strain Lang) (T1L)" ="10884",
                                     "Reovirus type 3 (strain Dearing)" ="10886",
                                     "Sendai virus"="11195,11198",
                                     "Semliki forest virus"="11033",
                                     "Sindbis virus" = "11034,31699",
                                     "Vesicular stomatitis virus" = "11276,434490,434489,434488,11278,11284,11285",
                                     #"Maraba virus",
                                     "Myxoma virus"="31530",
                                     "Bovine herpesvirus"="10320,10323,31518,10324,10385,79889",
                                     "Avian orthoreovirus" = "11862,160750,269446,11866,11879",
                                     #"Avian " ="642261",
                                     "Vaccinia virus" = '10245',
                                     "Vaccinia virus Copenhagen" = '10249',
                                     "Vaccinia virus WR" = '10254',
                                     "Vaccinia virus L-IPV" = '31531',
                                     "Vaccinia virus Ankara" = '126794',
                                     "All Vaccinia virus strain"='10245,10249,10254,31531,126794',
                                     "Vaccinia virus GLV-1h68" = "1",
                                     "Vaccinia virus VVΔTKΔN1L" ="10"
                                     
                         )
                       ),
                       actionBttn(
                         inputId = "runVisPPI2",
                         label = "Network", 
                         style = "gradient",
                         color = "danger",
                         icon = icon("circle-play")
                       ) 
                     ),
                     br(),
                     numericInput("qvalCut", "qValue", 0.05, min = 0, max = 1, step = 0.01),
                     numericInput("pvalCut", "pValue", 0.05, min = 0, max = 1, step = 0.01),
                     br(),
                     actionBttn(
                       inputId = "runKEGG",
                       label = "Run KEGG", 
                       style = "gradient",
                       color = "danger",
                       icon = icon("circle-play")
                     )
                   ),
                   # Right side box
                   box(
                     title = "Editable Network",
                     width = 9,
                     status = "primary", 
                     closable = FALSE,
                     maximizable = TRUE, 
                     collapsible = TRUE,
                     conditionalPanel(
                       condition = "input.runVisPPI2 > 0",
                       style = "display: none;",
                     withSpinner(visNetworkOutput("editable_network", height = 600)))
                   ),
                   box(
                     title = "KEGG Enrichment Output",
                     width = 12,
                     status = "primary", 
                     closable = FALSE,
                     maximizable = TRUE, 
                     collapsible = TRUE,
                     conditionalPanel(
                       condition = "input.runKEGG > 0",
                       style = "display: none;",
                     withSpinner(DT::dataTableOutput("kk_df")))
                   ),
                   
                   box(title = "Search Pathway",
                     status = "primary", collapsible = TRUE, width = 12,id = "pathview1",
                     searchInput(
                       inputId = "PathviEw",
                       placeholder = "Input pathway",
                       btnSearch = icon("search", verify_fa = FALSE), 
                       btnReset = icon("trash"),
                       width = "100%",
                       
                     ),
                     tags$div(class = "clearBoth"),
                     
                     imageOutput(outputId = "pathway1", inline = T),
                     br(),
                     br(),
                     # Place the download buttons in a separate column
                     downloadBttn('downloadPathviewPng', 'Download .png', style = "bordered", color = "primary"),
                     downloadBttn('downloadPathviewPdf', 'Download .pdf', style = "bordered", color = "primary")
                      ))
        )
        
        )),
      
      ################## Metabolic Analysis UI ########################
      tabItem(
        tabName = "metabolic_tab",
        fluidRow(
          box(
            id = "metabolic1",
            title = "", 
            width = 3,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            
            selectInput(
              inputId = "typeForMetabolic", 
              label = "Select Cancer Type Dataset", 
              choices = c("",
                          "Adrenocortical carcinoma (ACC)"= "TCGA-ACC",
                          "Bladder Urothelial Carcinoma (BLCA)"="TCGA-BLCA",
                          "Breast invasive carcinoma (BRCA)"="TCGA-BRCA",
                          "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)"="TCGA-CESC",
                          "Cholangio carcinoma (CHOL)"="TCGA-CHOL",
                          "Colon adenocarcinoma (COAD)"="TCGA-COAD",
                          "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)"="TCGA-DLBC",
                          "Esophageal carcinoma (ESCA)"="TCGA-ESCA",
                          "Glioblastoma multiforme (GBM)"="TCGA-GBM",
                          "Head and Neck squamous cell carcinoma (HNSC)"="TCGA-HNSC",
                          "Kidney Chromophobe (KICH)"= "TCGA-KICH",
                          "Kidney renal clear cell carcinoma (KIRC)"="TCGA-KIRC",
                          "Kidney renal papillary cell carcinoma (KIRP)"="TCGA-KIRP",
                          "Acute Myeloid Leukemia (LAML)"="TCGA-LAML",
                          "Brain Lower Grade Glioma (LGG)"="TCGA-LGG",
                          "Liver hepatocellular carcinoma (LIHC)"="TCGA-LIHC",
                          "Lung adenocarcinoma (LUAD)"="TCGA-LUAD",
                          "Lung squamous cell carcinoma (LUSC)"="TCGA-LUSC",
                          "Mesothelioma (MESO)"="TCGA-MESO",
                          "Ovarian serous cystadenocarcinoma (OV)"="TCGA-OV",
                          "Pancreatic adenocarcinoma (PAAD)"="TCGA-PAAD",
                          "Pheochromocytoma and Paraganglioma (PCPG)"="TCGA-PCPG",
                          "Prostate adenocarcinoma (PRAD)"="TCGA-PRAD",
                          "Rectum adenocarcinoma (READ)"="TCGA-READ",
                          "Sarcoma (SARC)"="TCGA-SARC",
                          "Skin Cutaneous Melanoma (SKCM)"="TCGA-SKCM",
                          "Stomach adenocarcinoma (STAD)"="TCGA-STAD",
                          "Testicular Germ Cell Tumors (TGCT)"="TCGA-TGCT",
                          "Thyroid carcinoma (THCA)"= "TCGA-THCA",
                          "Thymoma (THYM)"= "TCGA-THYM",
                          "Uterine Corpus Endometrial Carcinoma (UCEC)"= "TCGA-UCEC",
                          "Uterine Carcinosarcoma (UCS)"= "TCGA-UCS"
                          #"Uveal Melanoma (UVM)"= "TCGA-UVM"
              ), selected = 1),
            
            awesomeCheckbox(
              inputId = "userinput",
              label = "Use own data", 
              value = FALSE,
              status = "info",
            ),
            
            conditionalPanel(
              condition = "input.userinput",
              
              selectInput("organism",label = "Select Organism",
                          choices = c("Human"="hsa",
                                      "Mouse"="mmu",
                                      "Arabidopsis"="ath",
                                      "Yeast"="sce"
                          ), selected = 1),
              
              fileInput(inputId = "geneDE", "DE Genomic Data", 
                        accept = c("text/csv","text/comma-separated-values,text/plain", ".csv")),
              fileInput(inputId = "metDE", "DE Metabolic Data", 
                        accept = c("text/csv","text/comma-separated-values,text/plain", ".csv")),
              
              
            ),
            selectInput("network",label = "Select Network",
                        choices = c("KEGG network"="kegg",
                                    "Rhea network"="rhea",
                                    "Rhea lipid subnetwork"="lipidomic"),
                        selected = 1),
            selectInput("nodesAs",label = "Select Topology",
                        choices = c("metabolites"="metabolites",
                                    "atoms"="atoms"),
                        selected = 1),
            
            numericInput("kgene",
                         label=HTML("Number of positive genes"),
                         max=100, min=0, value=50, step=1),
            
            actionBttn(
              inputId = "runmetabolic",
              label = "RUN", 
              style = "gradient",
              color = "primary",
              icon = icon("circle-play")
            ) 
          ),
          tabBox(
            id = "metabolic2",
            title = "",
            width = 9,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            tabPanel("Data",
                     DT::dataTableOutput("metabolicgenedf")),
            tabPanel("Pathway",
                     DT::dataTableOutput("pathwaysTable")),
            tabPanel("Reaction",
                     DT::dataTableOutput("reactionTable"))
          )),
        box(
          id = "metabolic3",
          title = "",
          width = 12,
          status = "primary", 
          closable = FALSE,
          maximizable = TRUE, 
          collapsible = TRUE,
          visNetworkOutput("metabolicnetwork", height = 600),
          downloadLink('downloadNetwork', 'Download network')
          #downloadButton("downloadModule", "XGMML")
        )
      ),
      
      ################## Metabolic Flux Analysis with CNApy UI ########################
      tabItem(
        tabName = "flux_tab",
        fluidRow(
          box(
            id = "flux1",
            title = "", 
            width = 6,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            selectInput(
              inputId = "typeForFlux", 
              label = "Select Metabolic Flux Analysis", 
              choices = c("",
                          "Flux Balance Analysis"= "FBA",
                          "Parsimonious Flux Balance Analysis"="pFBA",
                          "Flux Variability Analysis"="FVA",
                          "Elementary Flux Modes"="EFM"
              ), selected = 1),
            
            prettyRadioButtons(
              inputId = "selectSBML",
              label = "Choose:", 
              choices = c("Upload own SBML/XML", "Create SBML"),
              icon = icon("check"), 
              bigger = TRUE,
              status = "info",
              animation = "jelly"
            ),
            conditionalPanel(
              condition = 'input.selectSBML == "Upload own SBML/XML"',
            fileInput(inputId = "fluxsbml", label = "XML file",
                      accept = c(".xml",".sbml"))),
            
            conditionalPanel(
              condition = 'input.selectSBML == "Create SBML"',
              br(),
              fileInput(inputId = "sbmlCsv", label = "Upload CSV file", accept = c(".csv")),
              DT::dataTableOutput("ediTableSbml"),
              actionBttn(
                inputId = "createSBML",
                label = "Create SBML", 
                style = "gradient",
                color = "primary",
                icon = icon("circle-play")
              )
            ),
            
            
            br(),
            DT::dataTableOutput("editableTableScenario"),
            actionBttn("addRowBtn", "Add Row"),  # Add Row button
            actionBttn("removeRowBtn", "Remove Row"),  # Remove Row button
            br(),
            br(),
            actionBttn(
              inputId = "runflux",
              label = "RUN", 
              style = "gradient",
              color = "primary",
              icon = icon("circle-play")
            )),
         
          tabBox(
            id = "flux3",
            title = "",
            #width = 6,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            tabPanel("Created Model",
                     tags$head(
                       tags$script(src = "https://unpkg.com/panzoom@9.4.0/dist/panzoom.min.js")
                     ),
                     
                     
                     div(
                       imageOutput("model_sbml", width = "140%", height = "80vh"),
                       style = "border: 2px solid #FFFFFF; padding: 10px; width: 90%; margin: 0 auto; max-width: 100%; max-height: 100%; overflow: hidden;"
                     ),
                     
                     tags$script(
                       HTML('panzoom($("#model_sbml")[0])')
                     ),
                     downloadButton("downSBMLfig", label = "Download SVG", style = "color: white; background-color: #4AA0F6;")
                     
                     #downloadButton("downSBMLfig", label = "Download SVG")
            ),
            tabPanel("Flux Result",
                     DT::dataTableOutput("fluxOutput")),
            tabPanel("Optimize",
                     DT::dataTableOutput("fluxOptimize")),
            tabPanel("Summary",
                     DT::dataTableOutput("fluxSummary")))),
          
          
          box(
            id = "flux4",
            title = "",
            width = 12,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            visNetworkOutput("fluxnetwork", height = 600),
            
            downloadButton("downloadModule2", "XGMML"))
      ),
      ########################## Immune Signature Analysis UI ###############################################
      
      tabItem(
        tabName = "analy6",
        fluidRow(
          box(width = 3, id = "Imsig",status = "primary", 
              
                              selectInput(
                                inputId = "typeForImsig", 
                                label = "Select Cancer Type Dataset", 
                                choices = c("",
                                            "Adrenocortical carcinoma (ACC)"= "TCGA-ACC",
                                            "Bladder Urothelial Carcinoma (BLCA)"="TCGA-BLCA",
                                            "Breast invasive carcinoma (BRCA)"="TCGA-BRCA",
                                            "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)"="TCGA-CESC",
                                            "Cholangio carcinoma (CHOL)"="TCGA-CHOL",
                                            "Colon adenocarcinoma (COAD)"="TCGA-COAD",
                                            "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)"="TCGA-DLBC",
                                            "Esophageal carcinoma (ESCA)"="TCGA-ESCA",
                                            "Glioblastoma multiforme (GBM)"="TCGA-GBM",
                                            "Head and Neck squamous cell carcinoma (HNSC)"="TCGA-HNSC",
                                            "Kidney Chromophobe (KICH)"= "TCGA-KICH",
                                            "Kidney renal clear cell carcinoma (KIRC)"="TCGA-KIRC",
                                            "Kidney renal papillary cell carcinoma (KIRP)"="TCGA-KIRP",
                                            "Acute Myeloid Leukemia (LAML)"="TCGA-LAML",
                                            "Brain Lower Grade Glioma (LGG)"="TCGA-LGG",
                                            "Liver hepatocellular carcinoma (LIHC)"="TCGA-LIHC",
                                            "Lung adenocarcinoma (LUAD)"="TCGA-LUAD",
                                            "Lung squamous cell carcinoma (LUSC)"="TCGA-LUSC",
                                            "Mesothelioma (MESO)"="TCGA-MESO",
                                            "Ovarian serous cystadenocarcinoma (OV)"="TCGA-OV",
                                            "Pancreatic adenocarcinoma (PAAD)"="TCGA-PAAD",
                                            "Pheochromocytoma and Paraganglioma (PCPG)"="TCGA-PCPG",
                                            "Prostate adenocarcinoma (PRAD)"="TCGA-PRAD",
                                            "Rectum adenocarcinoma (READ)"="TCGA-READ",
                                            "Sarcoma (SARC)"="TCGA-SARC",
                                            "Skin Cutaneous Melanoma (SKCM)"="TCGA-SKCM",
                                            "Stomach adenocarcinoma (STAD)"="TCGA-STAD",
                                            "Testicular Germ Cell Tumors (TGCT)"="TCGA-TGCT",
                                            "Thyroid carcinoma (THCA)"= "TCGA-THCA",
                                            "Thymoma (THYM)"= "TCGA-THYM",
                                            "Uterine Corpus Endometrial Carcinoma (UCEC)"= "TCGA-UCEC",
                                            "Uterine Carcinosarcoma (UCS)"= "TCGA-UCS"
                                            #"Uveal Melanoma (UVM)"= "TCGA-UVM"
                                ), selected = 1),
                             
                      
                              actionBttn(
                                inputId = "runImSig",
                                label = "RUN", 
                                style = "gradient",
                                color = "danger",
                                icon = icon("circle-play"))
          ),
          box(width = 9,id = "imSigBar",status = "primary", 
              tabsetPanel(
                tabPanel("Bar Plot",
                         fluidRow(
                           column(4,
                                  br(),
                                  selectInput(
                                    inputId = "ImsigType", 
                                    label = "Immune Cells", 
                                    choices = c("B cells" = "B_cells",
                                                "Interferon" = "Interferon",
                                                "Macrophages" = "Macrophages",
                                                "Monocytes" = "Monocytes",
                                                "Neutrophils" = "Neutrophils",
                                                "NK cells" = "NK_cells",
                                                "Plasma cells" = "Plasma_cells",
                                                "Proliferation" = "Proliferation",
                                                "T cells" = "T_cells",
                                                "Translation" = "Translation"), selected = 1)),
                           plotlyOutput("imsigplot1"))),
                tabPanel("Network",
                         fluidRow(visNetworkOutput("imsigplot2"),
                                  downloadLink('downloadNetwork2', 'Download network')))
              )   
          ),
          box(width = 12, id = "imSigRes",status = "primary", 
              tabsetPanel(
                tabPanel("Analysis Result",
                         br(),
                         fluidRow(DT::dataTableOutput("info_imsig"))),
                tabPanel("Signature Genes",
                         br(),
                         fluidRow(DT::dataTableOutput("info_imsig2"))),
                tabPanel("Result of Signature Genes",
                         br(),
                         fluidRow(DT::dataTableOutput("info_imsig_out")))
                
                
                
              )    
          )
        )
      ),
      
      #########################################################################
      
      
      ########################## Survival Analysis UI ###############################################
      
      tabItem(
        tabName = "analy3",
        fluidRow(
          box(width = 12,id = "Survival",status = "primary",
              fluidRow(width = 12,
                       closable = FALSE,
                       column(width=6,
                              selectInput(
                                inputId = "typeForSurv", 
                                label = "Select Cancer Type Dataset", 
                                choices = c("",
                                            "Adrenocortical carcinoma (ACC)"= "TCGA-ACC",
                                            "Bladder Urothelial Carcinoma (BLCA)"="TCGA-BLCA",
                                            "Breast invasive carcinoma (BRCA)"="TCGA-BRCA",
                                            "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)"="TCGA-CESC",
                                            "Cholangio carcinoma (CHOL)"="TCGA-CHOL",
                                            "Colon adenocarcinoma (COAD)"="TCGA-COAD",
                                            "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)"="TCGA-DLBC",
                                            "Esophageal carcinoma (ESCA)"="TCGA-ESCA",
                                            "Glioblastoma multiforme (GBM)"="TCGA-GBM",
                                            "Head and Neck squamous cell carcinoma (HNSC)"="TCGA-HNSC",
                                            "Kidney Chromophobe (KICH)"= "TCGA-KICH",
                                            "Kidney renal clear cell carcinoma (KIRC)"="TCGA-KIRC",
                                            "Kidney renal papillary cell carcinoma (KIRP)"="TCGA-KIRP",
                                            "Acute Myeloid Leukemia (LAML)"="TCGA-LAML",
                                            "Brain Lower Grade Glioma (LGG)"="TCGA-LGG",
                                            "Liver hepatocellular carcinoma (LIHC)"="TCGA-LIHC",
                                            "Lung adenocarcinoma (LUAD)"="TCGA-LUAD",
                                            "Lung squamous cell carcinoma (LUSC)"="TCGA-LUSC",
                                            "Mesothelioma (MESO)"="TCGA-MESO",
                                            "Ovarian serous cystadenocarcinoma (OV)"="TCGA-OV",
                                            "Pancreatic adenocarcinoma (PAAD)"="TCGA-PAAD",
                                            "Pheochromocytoma and Paraganglioma (PCPG)"="TCGA-PCPG",
                                            "Prostate adenocarcinoma (PRAD)"="TCGA-PRAD",
                                            "Rectum adenocarcinoma (READ)"="TCGA-READ",
                                            "Sarcoma (SARC)"="TCGA-SARC",
                                            "Skin Cutaneous Melanoma (SKCM)"="TCGA-SKCM",
                                            "Stomach adenocarcinoma (STAD)"="TCGA-STAD",
                                            "Testicular Germ Cell Tumors (TGCT)"="TCGA-TGCT",
                                            "Thyroid carcinoma (THCA)"= "TCGA-THCA",
                                            "Thymoma (THYM)"= "TCGA-THYM",
                                            "Uterine Corpus Endometrial Carcinoma (UCEC)"= "TCGA-UCEC",
                                            "Uterine Carcinosarcoma (UCS)"= "TCGA-UCS"
                                            #"Uveal Melanoma (UVM)"= "TCGA-UVM"
                                ), selected = 1)),
                       column(width=6,
                              selectInput(inputId = "inVar3", label = "Clinical Feature", choices = c(
                                
                                "synchronous_malignancy",
                                "ajcc_pathologic_stage",
                                "days_to_diagnosis",
                                "days_to_last_follow_up",
                                "tissue_or_organ_of_origin",
                                "primary_diagnosis",
                                "prior_malignancy",
                                "ajcc_pathologic_t",
                                "morphology",
                                "ajcc_pathologic_n",
                                "icd_10_code",
                                "site_of_resection_or_biopsy",
                                "race",
                                "gender",
                                "ethnicity",
                                "vital_status",
                                "age_at_index",
                                "treatments_pharmaceutical_treatment_type",
                                "treatments_pharmaceutical_treatment_or_therapy",
                                "treatments_radiation_treatment_type",
                                "treatments_radiation_treatment_or_therapy"
                              ))
                              
                       ),
                       actionBttn(
                         inputId = "runSurv",
                         label = "Plot", 
                         style = "gradient",
                         color = "danger",
                         icon = icon("circle-play")
                       )),
              conditionalPanel(
                condition = "input.runSurv > 0",
                style = "display: none;",
                withSpinner(plotlyOutput("SurvivalAn", height = 600)))
          ),
          box(width = 12,id = "SurvivalTable",
              status = "primary",
              closable = FALSE,
              conditionalPanel(
                condition = "input.runSurv > 0",
                style = "display: none;",
                withSpinner(dataTableOutput("survTable")))
          ),
          
          
          
          box(width = 12,id = "Survival2",status = "primary", 
              fluidRow(width = 12,
                       closable = FALSE,
                       maximizable = TRUE, 
                       collapsible = TRUE,
                       column(width=4,
                              prettyRadioButtons(
                                inputId = "IdGene",
                                label = "Search:", 
                                choices = c("Gene Name" = "gene_name", "Ensembl ID"="gene_id"),
                                icon = icon("check"), 
                                bigger = TRUE,
                                status = "info",
                                animation = "jelly"
                              )),
                       column(width=8,
                              uiOutput('inVarGene')     
                       ),
                       actionBttn(
                         inputId = "runSurv2",
                         label = "Plot", 
                         style = "gradient",
                         color = "danger",
                         icon = icon("circle-play")
                       ),
                       column(width=12,
                              conditionalPanel(
                                condition = "input.runSurv2 > 0",
                                style = "display: none;",
                                withSpinner(plotlyOutput("geneSurv",height = 600))))        
              )))
      ),
      
      #########################################################################
      
      ########################## SNV Analysis UI ###############################################
      tabItem(
        tabName = "analy4",
        fluidRow(
          box(
            id = "SNV",
            title = "", 
            width = 12,
            status = "primary", 
            closable = FALSE,
            maximizable = FALSE, 
            collapsible = FALSE,
            fluidRow(
              column(3,
                     selectInput(
                       inputId = "typeForSNV", 
                       label = "Cancer Type Dataset", 
                       choices = c("",
                                   "Adrenocortical carcinoma (ACC)"= "TCGA-ACC",
                                   "Bladder Urothelial Carcinoma (BLCA)"="TCGA-BLCA",
                                   "Breast invasive carcinoma (BRCA)"="TCGA-BRCA",
                                   "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)"="TCGA-CESC",
                                   "Cholangio carcinoma (CHOL)"="TCGA-CHOL",
                                   "Colon adenocarcinoma (COAD)"="TCGA-COAD",
                                   "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)"="TCGA-DLBC",
                                   "Esophageal carcinoma (ESCA)"="TCGA-ESCA",
                                   "Glioblastoma multiforme (GBM)"="TCGA-GBM",
                                   "Head and Neck squamous cell carcinoma (HNSC)"="TCGA-HNSC",
                                   "Kidney Chromophobe (KICH)"= "TCGA-KICH",
                                   "Kidney renal clear cell carcinoma (KIRC)"="TCGA-KIRC",
                                   "Kidney renal papillary cell carcinoma (KIRP)"="TCGA-KIRP",
                                   "Acute Myeloid Leukemia (LAML)"="TCGA-LAML",
                                   "Brain Lower Grade Glioma (LGG)"="TCGA-LGG",
                                   "Liver hepatocellular carcinoma (LIHC)"="TCGA-LIHC",
                                   "Lung adenocarcinoma (LUAD)"="TCGA-LUAD",
                                   "Lung squamous cell carcinoma (LUSC)"="TCGA-LUSC",
                                   "Mesothelioma (MESO)"="TCGA-MESO",
                                   "Ovarian serous cystadenocarcinoma (OV)"="TCGA-OV",
                                   "Pancreatic adenocarcinoma (PAAD)"="TCGA-PAAD",
                                   "Pheochromocytoma and Paraganglioma (PCPG)"="TCGA-PCPG",
                                   "Prostate adenocarcinoma (PRAD)"="TCGA-PRAD",
                                   "Rectum adenocarcinoma (READ)"="TCGA-READ",
                                   "Sarcoma (SARC)"="TCGA-SARC",
                                   "Skin Cutaneous Melanoma (SKCM)"="TCGA-SKCM",
                                   "Stomach adenocarcinoma (STAD)"="TCGA-STAD",
                                   "Testicular Germ Cell Tumors (TGCT)"="TCGA-TGCT",
                                   "Thyroid carcinoma (THCA)"= "TCGA-THCA",
                                   "Thymoma (THYM)"= "TCGA-THYM",
                                   "Uterine Corpus Endometrial Carcinoma (UCEC)"= "TCGA-UCEC",
                                   "Uterine Carcinosarcoma (UCS)"= "TCGA-UCS"
                                   #"Uveal Melanoma (UVM)"= "TCGA-UVM"
                       ), selected = 1),
                    
              ),
              column(3,
                     selectInput(
                       inputId = "PloTsnv", 
                       label = "Plot Type", 
                       choices = c("","Oncoplot"="snvplot1",
                                   "Heatmap"="snvplot2",
                                   "Lollipop Plot"="snvplot3"))),
              column(3,
                     numericInput(inputId =  "GeneSize",label =  "Top Gene Size", 5,min = 1, max = 60),
              ),
              #column(2, selectInput(inputId = "mutation",
              #                      label = "Mutation type",
              #                      choices = c("Somatic driver mutation" = "TRUE",
              #                                  "All" = "FALSE"))),
              
              column(3,
                     br(),
                     actionBttn(
                       inputId = "runSNV",
                       label = "Plot", 
                       style = "unite",
                       color = "danger",
                       icon = icon("circle-play"))
              )
            )),
          box(
            id = "SNVbox1",
            title = "",
            width = 12,
            status = "primary", 
            closable = FALSE,
            collapsible = TRUE,
            conditionalPanel(
              condition = "input.PloTsnv == 'snvplot3'",
              uiOutput('inVar5'),
              g3LollipopOutput(outputId = "lollipop", width = 1550, height = 600)
            ),
            conditionalPanel(
              condition = "input.PloTsnv == 'snvplot1'",
              div(
                style = "display: flex; justify-content: center; align-items: center; height: 100%;",
                div(
                  style = "text-align: center;",
                  plotOutput(outputId = "oncoplot1", width = 1610, height = 600)
                )
              )
            ),
            conditionalPanel(
              condition = "input.PloTsnv == 'snvplot2'",
              plotlyOutput(outputId = "SNVheatmap", height = 600)
            ),
            downloadButton(outputId = "downloadPlot", label = "Download")
            #)
            
          ),
          box(
            id = "SNVbox2",
            title = "",
            width = 12,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            conditionalPanel(
              condition = "input.PloTsnv == 'snvplot3'",
              withSpinner(dataTableOutput("SNVtable"))
            ),
            conditionalPanel(
              condition = "input.PloTsnv == 'snvplot2'",
              style = "display: none;",
              withSpinner(dataTableOutput("SNVheatTable"))
              #withSpinner(dataTableOutput("SNVpathway"))
            )
            #plotlyOutput("geneSurv")
          )
        )
      ),
      
      
      #########################################################################
      ########################## CNV Analysis UI ###############################################
      
      tabItem(
        tabName = "analy5",
        fluidRow(
          box(
            id = "CNV",
            title = "", 
            width = 12,
            status = "primary", 
            closable = FALSE,
            maximizable = FALSE, 
            collapsible = FALSE,
            fluidRow(
              column(3,
                     selectInput(
                       inputId = "typeForCNV", 
                       label = "Cancer Type Dataset", 
                       choices = c("",
                                   "Adrenocortical carcinoma (ACC)"= "TCGA-ACC",
                                   "Bladder Urothelial Carcinoma (BLCA)"="TCGA-BLCA",
                                   "Breast invasive carcinoma (BRCA)"="TCGA-BRCA",
                                   "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)"="TCGA-CESC",
                                   "Cholangio carcinoma (CHOL)"="TCGA-CHOL",
                                   "Colon adenocarcinoma (COAD)"="TCGA-COAD",
                                   "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)"="TCGA-DLBC",
                                   "Esophageal carcinoma (ESCA)"="TCGA-ESCA",
                                   "Glioblastoma multiforme (GBM)"="TCGA-GBM",
                                   "Head and Neck squamous cell carcinoma (HNSC)"="TCGA-HNSC",
                                   "Kidney Chromophobe (KICH)"= "TCGA-KICH",
                                   "Kidney renal clear cell carcinoma (KIRC)"="TCGA-KIRC",
                                   "Kidney renal papillary cell carcinoma (KIRP)"="TCGA-KIRP",
                                   "Acute Myeloid Leukemia (LAML)"="TCGA-LAML",
                                   "Brain Lower Grade Glioma (LGG)"="TCGA-LGG",
                                   "Liver hepatocellular carcinoma (LIHC)"="TCGA-LIHC",
                                   "Lung adenocarcinoma (LUAD)"="TCGA-LUAD",
                                   "Lung squamous cell carcinoma (LUSC)"="TCGA-LUSC",
                                   "Mesothelioma (MESO)"="TCGA-MESO",
                                   "Ovarian serous cystadenocarcinoma (OV)"="TCGA-OV",
                                   "Pancreatic adenocarcinoma (PAAD)"="TCGA-PAAD",
                                   "Pheochromocytoma and Paraganglioma (PCPG)"="TCGA-PCPG",
                                   "Prostate adenocarcinoma (PRAD)"="TCGA-PRAD",
                                   "Rectum adenocarcinoma (READ)"="TCGA-READ",
                                   "Sarcoma (SARC)"="TCGA-SARC",
                                   "Skin Cutaneous Melanoma (SKCM)"="TCGA-SKCM",
                                   "Stomach adenocarcinoma (STAD)"="TCGA-STAD",
                                   "Testicular Germ Cell Tumors (TGCT)"="TCGA-TGCT",
                                   "Thyroid carcinoma (THCA)"= "TCGA-THCA",
                                   "Thymoma (THYM)"= "TCGA-THYM",
                                   "Uterine Corpus Endometrial Carcinoma (UCEC)"= "TCGA-UCEC",
                                   "Uterine Carcinosarcoma (UCS)"= "TCGA-UCS"
                                   #"Uveal Melanoma (UVM)"= "TCGA-UVM"
                       ), selected = 1)),
              column(3,
                     selectInput(
                       inputId = "PloTcnv", 
                       label = "Plot Type", 
                       choices = c("","Oncoplot"="cnvplot1",
                                   "Interactive Oncoplot"="cnvplot2",
                                   "Lollipop Plot"="cnvplot3",
                                   "Somatic Interaction" = "cnvplot4"))),
              column(2,
                     numericInput(inputId =  "GeneSize2",label =  "Top Gene Size", 10,min = 5, max = 60),
              ),
              
              column(2,
                     actionBttn(
                       inputId = "runCNV",
                       label = "Plot", 
                       style = "gradient",
                       color = "danger",
                       icon = icon("circle-play"))
              )
            )),
          box(
            id = "CNVbox1",
            title = "",
            width = 12,
            status = "primary", 
            closable = FALSE,
            collapsible = TRUE
            
          ),
          box(
            id = "SNVbox2",
            title = "",
            width = 12,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE
            #conditionalPanel(
            # condition = "input.runSNV > 0",
            #style = "display: none;",
            #withSpinner(dataTableOutput("SNVtable"))
            # )
            #plotlyOutput("geneSurv")
          )
        )
      ),
      
      #########################################################################
      
      ################# cell line tab UI ###############################
      tabItem(
        tabName = "cell_line",
        fluidRow(
          
          box(
            id = "infcell",
            title = "Cell line information", 
            width = 12,
            status = "teal", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            searchInput(
              inputId = "searchCell",
              label = "Click search icon to update or hit 'Enter'", 
              placeholder = "Search cell line",
              btnSearch = icon("search"), 
              btnReset = icon("trash"),
              width = "100%"
              
              
            ),
            #tags$style(HTML('table.dataTable tr:nth-child(even) {background-color: #65646C !important;}')),
            #tags$style(HTML('table.dataTable tr:nth-child(odd) {background-color: #565661 !important;}')),
            #tags$style(HTML('table.dataTable th {background-color: white !important;}')),
            
            DT::dataTableOutput("info_celline")
            
          )
        )
        
      ),
      
      ################# cell type tab UI ###############################
      tabItem(
        tabName = "cell_type",
        fluidRow(
          
          box(
            id = "celltypeS",
            title = "Filtered by selected cell type", 
            closable = TRUE, 
            maximizable = TRUE,
            width = 12,
            status = "success", 
            solidHeader = FALSE, 
            collapsible = TRUE,
            selectInput(
              inputId = "tissueType", 
              label = "Select Tissue For Cell Type Samples", 
              choices = c("",
                          "Cultured embryonic stem cells", "Lung mesenchyme", "Kidney cortex", "Cervical and lumbar spinal cord", "Ventral midbrain", "Colon", "Caudal ganglionic eminence", "Medial ganglionic eminence", "Cortex", "Subcortex", "Lung progenitors", "Bone marrow stem cells", "Tonsil", "Mesenteric", "Cd14+ monocytes", "Embryonic kidney cortex", "Small intestine", "Aortic leukocytes", "Subventricular zone", "Bone marrow", "Alveolar rhabdomyosarcoma", "Trigeminal neurons", "Peripheral blood mononuclear cells", "Bladder", "Mammary gland", "Brain", "Embryo", "Fetal brain", "Fetal intestine", "Fetal kidney", "Fetal gonad", "Fetal liver", "Fetal lung", "Fetal stomach", "Kidney", "Lung", "Liver", "Mesenchymal stem cell", "Muscle", "Calvaria", "Rib", "Skin", "Ovary", "Pancreas", "Blood", "Placenta", "Stomach", "Spleen", "Testis", "Trophoblast", "Thymus", "Uterus", "Spinal cord", "Dorsal root ganglion", "Ventral striatum", "Hippocampus (ca1)", "Dentate gyrus", "Neo-cortex", "Medulla oblongata", "Thalamus", "Hippocampus", "Dorsal striatum", "Cerebellum", "Hypothalamus", "Amygdala", "Pons", "Olfactory bulb?", "Dorsal midbrain", "Mesencephalon?", "Cortex 1", "Olfactory bulb", "Cortex 3", "Cortex 2", "Somatosensory cortex", "Enteric nervous system", "Sympathetic nervous system", "Neurons", "Cortex, hippocampus and subventricular zone", "Mammary epithelial cells", "Left Ventricle", "Intestinal crypt epithelium", "Microglia", "Human embryo forebrain", "Airway submucosal glandular epithelial cells", "Epithelial cells in airway luminal surface", "Neonatal mouse stomach organoid", "Stomach organoid", "Neonatal mouse stomach explants", "Lks+ hematopoietic stem cells", "Arcuate-median eminence", "Lsk cells from bone marrow", "Mammary gland epithelium cell", "Peripheral blood mononuclear cell", "Pulmonary alveolar type i cells", "Pulmonary alveolar type ii cells", "Hepatocellular carcinoma", "B cells", "Lymph node", "Melanoma", "Kidney organoids", "Ovarian tumor", "Renal cortex", "Embryonic kidney", "Whole kidney", "Tumor", "Embryonic stem cell", "Heart", "Embryoid body", "Lacrimal gland", "Umbilical vein endothelial cells", "Breast epithelium", "Lateral geniculate nucleus", "Retina", "Medial amygdala", "Adipose", "Lateral ganglionic eminence", "Neuronal epithelium", "Aorta", "Heart valve", "Endothelial cells", "Retinal ganglion cells", "Substantia nigra-ventral tegmental area", "Breast cancer cell line", "Pancreatic islets", "Embryonic fibroblasts", "Circulating tumor cells in hepatocellular carcinoma", "Lung distal alveolar stromal cells", "Lung proximal airway stromal cells", "Oligodendrocyte cultures", "ES-derived embryoid body", "CD34+ cord blood cells", "NK cells", "Neural stem cells from subventricular zone", "Embryonic mouse mammary gland", "Cortex CA1", "Mammary epithelium", "In vitro differentiation of mESCs into motor neurons", "Striatum", "Nasal airway epithelium", "White adipose tissue", "Dermis", "Cortical organoids", "Cardiac progenitor cells", "Leukocytes from aortas", "Foam cells from aortas", "Lung airway epithelial cells", "Cardiac tissue", "T cells", "Epidermal innate lymphoid cells", "Dermal innate lymphoid cells", "Subcutaneous innate lymphoid cells", "Colon (Ulcerative Colitis)", "Hepatocyte-derived liver progenitor-like cells", "Hippocampus (24h Sham surgery)", "Hippocampus (injury)", "Skeletal muscle macrophages from hindlimb", "Merkel Cell Carcinoma", "PBMC", "Embryonic skin", "Olfactory epithelium", "NK cells from spleen", "NK cells from blood", "Bone marrow from melanoma tumor-bearing mice", "Umbilical cord blood", "Immortalized gonadotrope cells", "Splenocytes enriched for B cells", "Pituitary gland", "SVZ-derived neural stem cells", "Aortic cells", "CD45+ cells", "Testicular cells", "Testicle", "Decidua", "Intestinal organoid epithelial cells", "CD24+ CD45- neuronal cells", "Submandibular Gland", "Mouse embryonic fibroblasts", "Cochlea", "Whole heart", "Intraembryonic tissues", "Yolk sack", "Blood CD8+ T lymphocytes", "Cortical cells", "Preoptic region of the hypothalamus", "Periaortic lymph nodes", "Mesenteric lymph nodes", "Heart macrophages and dendritic cells", "Diencephalon", "Ear skin", "Peripheral lymph node", "Alpha cells", "Beta cells", "Pancreatic progenitor cells", "Trachea", "Tongue", "Marrow", "Mammary tissue", "Brain (neurons)", "Fat", "ES-derived kidney organoid", "Prostate", "Mandible", "Lymphoblastoid cell line", "Embryonic thymus", "Lung endoderm", "Monocyte-derived macrophages", "Leukemic cells", "Embryonic stem cell line", "Macrophage from E14.5", "Kaposi's sarcoma", "Thymocytes", "Prostate basal cells", "Nerve", "Sub-ventricular zone", "Visual cortex", "Intestinal epithelium", "Peritoneal cavity cells", "Muscle satellite cells", "Bone marrow mesenchymal and endothelial cells", "Inflamed area of colon", "Unknown", "Stria terminalis (brain)", "Substantia nigra", "Bed nucleus of the stria terminalis", "Cord blood", "Tonsil CD14+ cells", "Sorted myeloid cells", "Epidermal cells", "Hind limb muscles", "Tonsil CD4 T cells", "Commensal-specific CD8+ T cells", "Gonadal adipose tissue", "Lungs", "CD45+ dermal immune cells", "Distal small intestine", "Embryonic hindbrain", "Hindbrain", "Dura mater", "Choroid plexus", "Subdural meninges"
                          
              ), selected = 1),
            
            DT::dataTableOutput("out_celltype")
            
          ),
          
          
          box(
            id = "cellMarkers",
            title = "", 
            closable = TRUE, 
            maximizable = TRUE,
            width = 12,
            status = "success", 
            solidHeader = FALSE, 
            collapsible = TRUE,
            selectInput(
              inputId = "cellType", 
              label = "Select Cell Type for Gene Markers", 
              choices = c("","Normal cell","Cancer cell"
                          
              ), selected = 1),
            
            DT::dataTableOutput("out_cellMarker")
            
          ) 
          
        )
        
      ),
      
      ################# Jimena tab UI ###############################
      tabItem(
        tabName = "cancerModel",
        fluidRow(
          box(
            id = "models",
            title = "", 
            #width = 7,
            status = "navy", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            prettyRadioButtons(
              inputId = "choseJimena",
              label = "Choose:", 
              choices = c("Create data" = "createJim","Use example data" = "exampleJim", "Upload your data"="userJim"),
              icon = icon("check"), 
              bigger = TRUE,
              status = "info",
              animation = "jelly"
            ),
            
            conditionalPanel(condition = "input.choseJimena=='createJim'",
                             searchInput("searchProtein", h4("Search Protein"), placeholder = "Enter gene symbol",btnSearch = icon("search", verify_fa = FALSE), 
                                         btnReset = icon("trash")),
                             DTOutput("source_table"),
                             br(),
                             br()
                             
            ),
            
            conditionalPanel(condition = "input.choseJimena=='userJim'",
                             
                             fileInput(inputId = "file1", h4("Choose Txt File"), 
                                       accept = c("text/csv","text/comma-separated-values,text/plain", ".csv"))),
            withSpinner(DT::dataTableOutput(outputId = "whichjim")),
            actionBttn(
              inputId = "runGraphml",
              label = "Convert Graphml",
              style = "fill", 
              color = "warning"
            ),
            actionBttn(
              inputId = "runJimena",
              label = "Run Jimena",
              style = "gradient", 
              color = "primary",
              icon = icon("glyphicon glyphicon-play",lib = "glyphicon")
            )),
          
          box(
            id = "targeTable",
            title = "", 
            #width = 7,
            status = "navy", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            DTOutput("target_table"),
            actionBttn("remove_button", "Remove Selected Rows"),
            br(),
            br()
          ),
          box(
            width = 6,
            status = "navy", solidHeader = TRUE, collapsible = TRUE, closable = TRUE, id = "pertbox",
            title = HTML("Perturbation setting", 
                         # as.character(actionLink(inputId = "info9", 
                         #                         label = "Step 2", 
                         #                         icon = icon("info-circle", verify_fa = FALSE,"<font size='5'>"),))
            ), 
            
            uiOutput('inVar2'),
            
            numericInput("strt", "Start", 0,
                         min = 0, max = 1000),
            numericInput("end", "End", 1000,
                         min = 1, max = 1000),
            numericInput("val", "Value", 1,
                         min = 0, max = 1, step =0.01),
            
            
            actionBttn(inputId = "add_pert",
                       label = "Add New Perturbation",
                       color = "success",
                       style = "unite"),
            
            actionBttn(inputId = "remove_per",
                       label = "Remove Perturbation",
                       color = "danger",
                       style = "unite")
            
          ),
          
          box(width = 6,id = "jimedit",status = "navy",
              DT::dataTableOutput("shiny_table")
          ),
          box(status = "navy", collapsible = TRUE, width = 12,id = "simjim",
              
              #conditionalPanel(
              # condition = "input.runJimena > 0",
              #style = "display: none;",
              withSpinner(plotlyOutput("jimenaResult")),
              #),
              br()),
          box( status = "navy",collapsible = TRUE, id = "jimena_table",width = 12,
               hr(),withSpinner(DT::dataTableOutput("jimenaoutput_tb")),data.step = 1)
        )
        
      ),
      
      
      ################# drugs-tissue tab ###############################
      tabItem(
        tabName = "drug1",
        fluidRow(
          box(
            id = "drugTissue",
            title = "", 
            width = 12,
            status = "purple", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            selectInput(
              inputId = "tissue", 
              label = "Select Tissue", 
              choices = c("",
                          "Acute myeloid leukemia (LAML)" = "LAML",
                          "Acute lymphoblastic leukemia (ALL)" = "ALL",
                          "Adrenocortical carcinoma (ACC)" = "ACC",
                          "Bladder urothelial carcinoma (BLCA)" = "BLCA",
                          "Brain lower grade glioma (LGG)" ="LGG",
                          "Breast invasive carcinoma (BRCA)" = "BRCA",
                          "Colon and rectum adenocarcinoma (COAD/READ)" ="COREAD",
                          "Chronic myelogenous leukemia (ALL)" = "ALL",
                          "Chronic lymphocytic leukemia (CLL)" = "CLL",
                          "Cervical squamous cell carcinoma (CESC)" ="CESC",
                          "Esophageal carcinoma (ESCA)" = "ESCA",
                          "Glioblastoma multiforme (GBM)"="GBM",
                          "Head and neck squamous cell carcinoma (HNSC)"="HNSC",
                          "Kidney renal clear cell carcinoma (KIRC)"="KIRC",
                          "Liver hepatocellular carcinoma (LIHC)"="LIHC",
                          "Lung adenocarcinoma (LUAD)"="LUAD",
                          "Lung squamous cell carcinoma (LUSC)"="LUSC",
                          "Mesothelioma (MESO)" = "MESO",
                          "Medulloblastoma (MB)" = "MB",
                          "Multiple myeloma (MM)" = "MM",
                          "Neuroblastoma (NB)" ="NB",
                          "Lymphoid neoplasm diffuse large B-cell lymphoma (DLBC)"="DLBC",
                          "Ovarian secrous cystadenocarcinoma (OV)"="OV",
                          "Pancreatic adenocarcinoma (PAAD)"="PAAD",
                          "Prostate adenocarcinoma (PRAD)"="PRAD",
                          "Small cell lung cancer (SCLC)" = "SCLC",
                          "Skin cutaneous melanoma (SKCM)"="SKCM",
                          "Stomach adenocarcinoma (STAD)"="STAD",
                          "Thyroid carcinoma (THCA)"="THCA",
                          "Uterine corpus endometrial carcinoma (UCEC)"="UCEC",
                          "NO TCGA classification"="UNCLASSIFIED"
              ), selected = 1),
            
            
            DT::dataTableOutput("drug")
            
          )
        )
        
      ),
      #############################################################################
      
      ################# drugs-description tab ###############################
      tabItem(
        tabName = "drug2",
        tags$head(
          tags$style(HTML("
      .small-img {
        max-width: 300px;
        max-height: 300px;
        #border: 1px solid gray;
        padding: 2px;
        margin: 2px;
      }
      
      .popover {
        max-width: 600px;
        background: none;
        width: 600px !important;
        height: 500px !important;
      }
      

      .img-popover {
        max-width: 600px;
        max-height: 600px;
      }
    "))
        ),
        fluidRow(
          column(5, 
                 box(
                   id = "drugDesc",
                   title = "", 
                   width = 12,
                   status = "purple", 
                   closable = FALSE,
                   maximizable = TRUE, 
                   collapsible = TRUE,
                   prettyRadioButtons(
                     inputId = "choseDrug",
                     label = "Choose:", 
                     choices = c("Drug" = "searchDrug","Indication" = "searchIndic", "Smiles"="searchSmi"),
                     icon = icon("check"), 
                     bigger = TRUE,
                     status = "info",
                     animation = "jelly"
                   ),
                   
                   conditionalPanel(condition = "input.choseDrug=='searchDrug'",
                                    searchInput(
                                      inputId = "searchDrugName",
                                      label = "Click search icon to update or hit 'Enter'", 
                                      placeholder = "Search drug name",
                                      btnSearch = icon("search"), 
                                      btnReset = icon("trash")
                                      #width = "100%",
                                      
                                    )),
                   
                   conditionalPanel(condition = "input.choseDrug=='searchIndic'",
                                    searchInput(
                                      inputId = "searchIndication",
                                      label = "Click search icon to update or hit 'Enter'", 
                                      placeholder = "Search indication",
                                      btnSearch = icon("search"), 
                                      btnReset = icon("trash")
                                      #width = "100%",
                                    )),
                   conditionalPanel(condition = "input.choseDrug=='searchSmi'",
                                    searchInput(
                                      inputId = "searchSmile",
                                      label = "Click search icon to update or hit 'Enter'", 
                                      placeholder = "Search indication",
                                      btnSearch = icon("search"), 
                                      btnReset = icon("trash")
                                      #width = "100%",
                                      
                                    ))
                   
                   
                   
                 ),
                 
                 conditionalPanel(condition = "input.choseDrug=='searchDrug'",
                                  box(
                                    id = "drugStruc",
                                    title = "Structure of drug", 
                                    width = 12,
                                    status = "purple", 
                                    closable = FALSE,
                                    maximizable = TRUE, 
                                    collapsible = TRUE,
                                    uiOutput("my_output"),
                                    verbatimTextOutput("drug_info")))
          ),
          
          
          column(7, 
                 conditionalPanel(condition = "input.choseDrug=='searchDrug'",
                                  box(
                                    id = "drugDesTab",
                                    title = "", 
                                    width = 12,
                                    status = "purple", 
                                    closable = FALSE,
                                    maximizable = TRUE, 
                                    collapsible = TRUE,
                                    DT::dataTableOutput("drug_desc"))
                 ),
                 conditionalPanel(condition = "input.choseDrug=='searchIndic'",
                                  box(
                                    id = "drugDesTab2",
                                    title = "", 
                                    width = 12,
                                    status = "purple", 
                                    closable = FALSE,
                                    maximizable = TRUE, 
                                    collapsible = TRUE,
                                    DT::dataTableOutput("drug_indic"))
                 ),
                 conditionalPanel(condition = "input.choseDrug=='searchSmi'",
                                  box(
                                    id = "drugDesTab3",
                                    title = "", 
                                    width = 12,
                                    status = "purple", 
                                    closable = FALSE,
                                    maximizable = TRUE, 
                                    collapsible = TRUE,
                                    DT::dataTableOutput("drug_sim"))
                 )
                 
          ))
        
        
      ),
      #############################################################################
      
      ################# drugs-genes tab ###############################
      tabItem(
        tabName = "geneDrug",
        fluidRow(
          box(
            id = "drugGenes",
            title = "", 
            width = 12,
            status = "info", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            searchInput(
              inputId = "searchBttn",
              label = "Click search icon to update or hit 'Enter'", 
              placeholder = "Search gene",
              btnSearch = icon("search"), 
              btnReset = icon("trash"),
              width = "100%"
            ),
            DT::dataTableOutput("searchgene")
          )
        )
        
      ),
      
      ################# STITCH tab ###############################
      tabItem(
        tabName = "stitch",
        fluidRow(
          box(
            id = "stitchWeb",
            title = "", 
            width = 12,
            status = "info", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            uiOutput("stitch_web") 
          )
        ) 
      ),
      
      ################# About us tab ###############################
      tabItem(
        tabName = "usAbout",
        fluidRow(
          box(
            id = "members",
            title = "", 
            width = 12,
            status = "info", 
            closable = FALSE,
            maximizable = FALSE, 
            collapsible = FALSE,
            #style = "background-color: #80CDF7;", 
            fluidRow(
              column(width = 4, align = "center",
                     imageOutput("person1"),
                     HTML('<div class="details" style="font-size: 24px; font-family: monospace; margin-top: -100px;">
                  <p><strong>Thomas Dandekar</strong></p>
                  <p>Leitung</p>
                </div>')
              ),
              column(width = 4, align = "center",
                     imageOutput("person2"),
                     HTML('<div class="details" style="font-size: 24px; font-family: monospace;  margin-top: -100px;">
                  <p><strong>Elena Bencurova</strong></p>
                  <p>Molekularbiologie</p>
                </div>')
              ),
              column(width = 4, align = "center",
                     imageOutput("person3"),
                     HTML('<div class="details" style="font-size: 24px; font-family: monospace; margin-top: -100px;">
                  <p><strong>Rana Salihoglu</strong></p>
                  <p>Programmierung</p>
                </div>')
              )
            ))
        ) 
      ),
      
      ################# Uses Database tab ###############################
      tabItem(
        tabName = "DataSour",
        fluidRow(
          box(
            id = "databaseinf",
            title = "", 
            width = 12,
            status = "info", 
            closable = FALSE,
            maximizable = FALSE, 
            collapsible = FALSE,
            style = "background-color: #FFFFFF;",
            
            #HTML('<a href="http://www.ckttdb.org/" target="_blank">CKTTDB</a>'),br(),br(),
            
            HTML('<a href="https://webs.iiitd.edu.in/raghava/ovirustdb/" target="_blank">OvirusTdb</a>'),br(),br(),
            HTML('<a href="https://www.ebi.ac.uk/intact/home" target="_blank">IntAct</a>'),br(),br(),
            HTML('<a href="https://depmap.org/portal/" target="_blank">depmap</a>'),br(),br(),
            HTML('<a href="http://117.50.127.228/CellMarker/index.html/" target="_blank">CellMarker</a>'),br(),br(),
            HTML('<a href="https://panglaodb.se/" target="_blank">PanglaoDB</a>'),br(),br(),
            HTML('<a href="https://go.drugbank.com/" target="_blank">DRUGBANK</a>'),br(),br(),
            HTML('<a href="https://www.cancerrxgene.org/" target="_blank">cancerrxgene</a>'),br(),br(),
            HTML('<a href="https://pubchem.ncbi.nlm.nih.gov/" target="_blank">pubchem</a>'),br(),br(),
            HTML('<a href="https://www.dgidb.org/" target="_blank">dgidb</a>'),br(),br(),
            HTML('<a href="http://stitch.embl.de/" target="_blank">STITCH</a>'),br(),br(),
            HTML('<a href="https://classic.clinicaltrials.gov/ct2/home" target="_blank">clinicalTrials.gov</a>'),br(),br(),
            HTML('<a href="http://www.ckttdb.org/" target="_blank">CKTTDB</a>'),br(),br() 
          )
        )),
      
      ################# Hotlist tab ###############################
      tabItem(
        tabName = "hotList",
        fluidRow(
          tags$head(
            tags$style(
              HTML('
        .link-description {
          font-family: Arial, sans-serif; /* Change font type */
          font-size: 16px; /* Change font size */
        }
        
        .link {
          font-family: Arial, sans-serif; /* Change font type */
          font-size: 16px; /* Change font size */
          color: #007BFF; /* Change link color */
        }
      '))
          ),
          box(
            id = "linksForHotlist",
            title = "", 
            width = 12,
            status = "info", 
            closable = FALSE,
            maximizable = FALSE, 
            collapsible = FALSE,
            style = "background-color: #FFFFFF;",
            HTML('<p class="link-description"><a class="link" href="https://www.cochranelibrary.com" target="_blank">Cochrane Library</a> - The Cochrane Library is a reputable online resource for evidence-based healthcare research and systematic reviews.</p>'),br(),br(),
            HTML('<p class="link-description"><a class="link" href="https://drumpid.bioapps.biozentrum.uni-wuerzburg.de/compounds/index.php" target="_blank">DrumpID</a> - The DrumPID database is a comprehensive resource offering customized information on drugs and their protein networks, including indications, targets, and side-targets, making it invaluable for drug development, predicting side-effects, and studying structure-activity relationships.</p>'),br(),br(),
            HTML('<p class="link-description"><a class="link" https://europepmc.org" target="_blank">Europe PMC</a> - Europe PMC is a valuable online platform for accessing a wide range of biomedical literature and research articles, facilitating easy access to scientific information for researchers and healthcare professionals.</p>'),br(),br(),
            HTML('<p class="link-description"><a class="link" href="https://www.novoprolabs.com/tools/convert-peptide-to-smiles-string" target="_blank">Convert Peptide to Smiles</a> - The NovoPro Labs peptide-to-SMILES string conversion tool provides a convenient online resource for converting peptide sequences into simplified molecular input line entry system (SMILES) notation for chemical structure representation and analysis.</p>'),br(),br()
          )
        )
      ),
      
      
      
      ################# CAR-T cell therapy tab ###############################
      tabItem(
        tabName = "carTab",
        fluidRow(
          box(
            id = "CarTarget",
            title = "", 
            width = 12,
            status = "purple", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            selectInput(
              inputId = "carDisease", 
              label = "Targets with CAR T-Cell Therapy by Disease", 
              choices = c("",
                          "Solid Tumor",
                          "Brain Cancer",
                          "Spinal cord",
                          "Myelodysplastic Syndrome",
                          "Leukaemia",
                          "Sarcoma",
                          "Acute Myeloid Leukaemia",
                          "Precursor B-lymphoblastic Neoplasm",
                          "Lymphoma",
                          "Follicular Lymphoma",
                          "Diffuse large B-cell Lymphoma",
                          "Mature B-cell Leukaemia",
                          "Multiple Myeloma",
                          "Mature B-cell Lymphoma",
                          "B-cell Lymphoma",
                          "Mature T-cell Lymphoma",
                          "Mycosis Fungoides",
                          "Hodgkin Lymphoma",
                          "Immunoproliferative Disorder",
                          "Malignant Haematopoietic Neoplasm",
                          "Osteosarcoma",
                          "Synovial Sarcoma",
                          "Nasopharyngeal Cancer",
                          "Esophageal Cancer",
                          "Stomach Cancer",
                          "Colon Cancer",
                          "Colorectal Cancer",
                          "Pancreatic Cancer",
                          "Liver Cancer",
                          "Lung Cancer",
                          "Pleural Mesothelioma",
                          "Melanoma",
                          "Mesothelin Positive Tumor",
                          "Peritoneal Cancer",
                          "Breast Cancer",
                          "Ovarian Cancer",
                          "Fallopian Tube Cancer",
                          "Endometrial Cancer",
                          "Cervical Cancer",
                          "Testicular Cancer",
                          "Prostate Cancer",
                          "Renal Cell Carcinoma",
                          "Bladder Cancer",
                          "Thyroid Cancer",
                          "Adrenal Cancer",
                          "Adenocarcinoma",
                          "Metastatic Malignant Neoplasm",
                          "Metastatic Tumor",
                          "Metastatic Pleura Neoplasm",
                          "Metastatic Peritoneum Neoplasm",
                          "Immune System Disease",
                          "Lupus Erythematosus",
                          "Lymphatic Vessel/Lymph Nodes Disorder",
                          "Ascites"
              ), selected = 1),
            DT::dataTableOutput("car_Therapy") 
          )
        ) 
      ),
      
      
      ##########################################################################
      
      ################# Bispecific therapy tab ###############################
      tabItem(
        tabName = "bispecificTab",
        fluidRow(
          box(
            id = "BispecificA",
            title = "", 
            width = 12,
            status = "purple", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            DT::dataTableOutput("bispecific_Therapy") 
          )
        ) 
      ),
      
      #############################################################################
      ################# Cytostatic therapy tab ###############################
      tabItem(
        tabName = "cytostaticTab",
        fluidRow(
          box(
            id = "cytostaticBox",
            title = "", 
            width = 12,
            status = "purple", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            DT::dataTableOutput("cytostatic_Therapy")
          )
        ) 
      ),
      
      #############################################################################
      ################# Oncolytic Virus therapy tab ###############################
      tabItem(
        tabName = "ovirusTherapyTab",
        fluidRow(
          box(
            id = "ovirusTherapyBox",
            title = "", 
            width = 12,
            status = "purple", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            DT::dataTableOutput("ovirus_Therapy")
          )
        )  
      ),
      
      #############################################################################
      
      ################# CheckPoint Inhbitors tab ###############################
      tabItem(
        tabName = "checkPoinTab",
        fluidRow(
          box(
            id = "check1",
            title = "", 
            width = 12,
            status = "purple", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            DT::dataTableOutput("checkpoint_target")
          ),
          
          box(
            id = "check2",
            title = "", 
            width = 12,
            status = "purple", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            DT::dataTableOutput("checkpoint_modulator")
          ) 
        )
      ),
      #############################################################################
      
      ###### Protein structure tab #################
      tabItem(
        tags$head(
          tags$script(src = "https://3Dmol.org/build/3Dmol-min.js"),
          tags$script(src = "https://3Dmol.org/build/3Dmol.ui-min.js")
        ),
        tabName = "strucPro",
        fluidRow(
          box(
            id = "alphafold2",
            title = "Options", 
            width = 2,
            solidHeader = TRUE,
            status = "primary", 
            closable = FALSE,
            textInput("pdb_name", "Enter PDB name (UniProt ID):", value = "Q5VSL9"),
            actionBttn(
              inputId = "download_btn",
              label = "Search", 
              style = "stretch",
              color = "warning"
            ),
            br(),
            textInput("selection", "Selection", "1-20"),
            selectInput("type_structure", "Structure", c( "cartoon","hide", "ball+stick", "ribbon", "backbone", "licorice", "spacefill", "line", "contact",
                                                          "helixorient", "hyperball", "rocket")),
            selectInput("type_surface", "Surface", c("hide", "electrostatic", "uniform", "hydrobhobicity", "element", "bfactor", "residueindex", "random", "resname",
                                                     "sstruc")),
            selectInput("type_ligand", "Ligand", c("hide","ball+stick", "spacefill", "line", "contact",
                                                   "helixorient", "hyperball", "rocket")),
            colourInput("col", "Select colour", "yellow"),
            colourInput("backgroundColor", "Background color", "black"),
            #selectInput("backgroundColor", "Background Color", c("black", "white", "grey")),
            actionButton("add", "Add"),
            actionButton("remove", "Remove"),
            actionBttn(
              inputId = "screenshot",
              label = "Save", 
              style = "minimal",
              color = "primary"
            ),
            
            br(),
            br(),
            materialSwitch(
              inputId = "fullscreen",
              label = "Fullscreen", 
              value = TRUE,
              status = "success"
            ),
            
            prettyRadioButtons(
              inputId = "animat",
              label = "Animation", 
              choices = c("None"="none", "Spin"="spin", "Rock"="rock"),
              icon = icon("check"), 
              bigger = TRUE,
              status = "info",
              animation = "jelly"
            )
          ),
          box(
            id = "alphafold",
            title = "AlphaFold Protein Structure Prediction", 
            width = 10,
            solidHeader = TRUE,
            status = "primary", 
            closable = FALSE,
            maximizable = TRUE, 
            collapsible = TRUE,
            sidebarPanel(
              h4("Legend"),
              uiOutput("legend")
            ),
            #textOutput("info2"),
            NGLVieweROutput("structure", width = "100%", height = "800px")
          ),
          
          column(
            width = 11,  # Adjust the width as needed
            offset = 1, # Adjust the offset to center the box
            box(
              id = "2D3D",
              title = "3D Protein Structure", 
              width = 11,
              solidHeader = TRUE,
              status = "success", 
              closable = FALSE,
              maximizable = TRUE, 
              collapsible = TRUE,
              br(),
              HTML('<p class="link-description"><a class="link" href="https://www.rcsb.org/" target="_blank">Find your pdb code</a> </p>'),br(),br(),
              
              tags$div(
                style = "height: 700px; 
              width: 1250px; 
              position: relative ",
                class = "viewer_3Dmoljs", 
                "data-pdb" = "1WQ9", 
                "data-backgroundcolor" = "0xffffff", 
                "data-style" = "stick", 
                "data-ui" = "true"
              )
            ))
        ) 
      )
      
      #############################################################################
    )
  ),
  
  
  
  ##########____________________________________________________________________
  
  controlbar = dashboardControlbar(
    id = "controlbar",
    skin = "light",
    pinned = TRUE,
    overlay = FALSE,
    controlbarMenu(
      id = "controlbarMenu",
      type = "pills",
      controlbarItem(
        "Skin",
        skinSelector()
      )
    )
  ),
  footer = dashboardFooter(
    fixed = FALSE,
    left = a(
      href = "",
      target = ""
    ),
    right = "2023"
  ),
  title = "Cat-E-Cancer Theraphy Explorer"
)



#### S E R V E R ########################

server = function(input, output, session) {
  options(shiny.maxRequestSize=90*1024^2)
  useAutoColor()
  source("connection-mysql.R",local = TRUE)
  source("mysql_queries.R",local = TRUE)
  
  ovirus <- read_excel("data/ovirusTDB.xlsx")
  ovirus2 <- reactiveValues(data= ovirus)
  
  #### home piechart ######
  output$pie <-  renderBillboarder({
    for_donut <- ovirus2$data
    
    y <-for_donut %>% 
      group_by(Virus_Name) %>% 
      summarise( percent = 100 * n() / nrow(for_donut))
    # Main title
    billboarder() %>% 
      bb_donutchart(y) %>% bb_legend(position = 'right')%>%
      bb_donut(title = "Oncolytic Virus") %>% 
      bb_title(
        #text = "Main title",
        padding = list(top = 20, right = 10, bottom = 10, left = 10), 
        position = "left"
      ) %>% 
      bb_add_style(
        ".bb-chart-arcs-title" = "font-size: 34px; fill: black;", 
        ".bb-title" = "font-size: 42px; fill: green;"
      )
  })
  
  output$home1 <- renderUI({
    div(
      img(
        src = "cancer1.png", # Change this to the path of your image
        class = "small-img",
        id = "img1"
      ),
      tags$script(HTML('$("#img1").popover({trigger: "hover", html: true, content: "<img src=cancer1.png class=img-popover>"});'))
    )
  })
  
  output$home2 <- renderUI({
    div(
      img(
        src = "cancer2.png", # Change this to the path of your image
        class = "small-img",
        id = "img1"
      ),
      tags$script(HTML('$("#img1").popover({trigger: "hover", html: true, content: "<img src=cancer2.png class=img-popover>"});'))
    )
  })
  
  output$home3 <- renderUI({
    div(
      img(
        src = "cancer3.png", # Change this to the path of your image
        class = "small-img",
        id = "img1"
      ),
      tags$script(HTML('$("#img1").popover({trigger: "hover", html: true, content: "<img src=cancer3.png class=img-popover>"});'))
    )
  })
  
  output$home4 <- renderUI({
    div(
      img(
        src = "viruscancer.png", # Change this to the path of your image
        class = "small-img",
        id = "img1"
      ),
      tags$script(HTML('$("#img1").popover({trigger: "hover", html: true, content: "<img src=viruscancer.png class=img-popover>"});'))
    )
  })
  
  output$home5 <- renderUI({
    div(
      img(
        src = "cancer4.png", # Change this to the path of your image
        class = "small-img",
        id = "img1"
      ),
      tags$script(HTML('$("#img1").popover({trigger: "hover", html: true, content: "<img src=cancer4.png class=img-popover>"});'))
    )
  })
  
  ######## Help button #####
  ##########################
  source("HelpInfo.R")
  
  observeEvent(
    eventExpr = input$help3,
    if(input$tabs > 0){
      
      handlerExpr = {
        introjs(session, 
                options = list(
                  "showBullets"="false", "showProgress"="true", 
                  "showStepNumbers"="false","nextLabel"="Next","prevLabel"="Prev","skipLabel"="Skip",
                  steps=helptext()[tab == input$tabs]
                ),events = list(onbeforechange = readCallback("switchTabs"))
        )
      }
    })
  
  
  ##########################################
  #########################################
  
  
  output$plotCase <- renderPlotly({
    dframe <- read.csv("data/cases.csv")
    colors34 <- brewer.pal(12, "Paired")
    colors34 <- rep(colors34, length.out = 34)
    
    # View the custom palette
    plot_ly(
      data = dframe,
      x = ~ Cases,
      y = ~ Project,
      type = 'bar',
      color = ~ Project,
      colors = colors34
    )
  })
  
  
  #_____________________________________________________________________________
  
  
  ######### oncolytic virus table
  output$virus <- DT::renderDataTable(server=FALSE,{
    if (input$onc_spec > 0) {
      filtered_oncv <- ovirus %>% filter(Virus_Name == input$onc_spec)
      
      DT::datatable(filtered_oncv,
                    style="bootstrap4",
                    extensions = c('Buttons','Responsive','Scroller'),
                    options = list(
                      dom = 'Bfrtip',
                      scrollY = 400,
                      ordering=TRUE,
                      #extensions = c('Buttons','Responsive','Scroller'), 
                      #options = list(
                      # dom = 'Bfrtip',
                      #lengthChange = FALSE,
                      #deferRender = TRUE,
                      #ordering = TRUE,
                      #scrollY = 400,
                      #scroller = TRUE,
                      buttons = list(
                        list(extend = "collection", buttons = c('excel','csv'),
                             text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                        list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'tt'))
                    ))
    }
  })
  
  
  ############# Virus VisNetwork ###########################################
  
  observeEvent(input$runVisPPI,{
    output$viRusNetwork <- renderVisNetwork({
      if (input$virusID > 0){
        print(input$virusID)
        edges <-e$data
        
        colnames(edges) <- c("from","Gene_1","Desc_1","tax_from","organism_1", "link","to","Gene_2","Desc_2", "tax_to","organism_2")
        
        # Find the indices of empty elements in Column1
        empty_indices <- which(is.na(edges$Gene_1) | edges$Gene_1 == "")
        
        # Fill in the empty elements in Column1 with data from Column2
        edges$Gene_1[empty_indices] <- edges$from[empty_indices]
        
        # Find the indices of empty elements in Column1
        empty_indices <- which(is.na(edges$Gene_2) | edges$Gene_2 == "")
        
        # Fill in the empty elements in Column1 with data from Column2
        edges$Gene_2[empty_indices] <- edges$to[empty_indices]
        
        nodevistb <- edges
        
        colnames(nodevistb) <- c("From_A","Gene_A","Description_A","Taxonomy_A","Organism_A", "Interaction_Type","To_B","Gene_B","Description_B", "Taxonomy_B","Organism_B")
        nodevistb <- nodevistb[,c(1,7,4,2,8,5,10,6,11,3,9)]
        
        Edges <- as.data.frame(nodevistb$Gene_A)
        Edges$to <- nodevistb$Gene_B
        Edges$link <- nodevistb$Interaction_Type
        colnames(Edges) <- c("from","to","link")
        
        nodes1 <- as.data.frame(nodevistb$From_A)
        nodes1$tax <- nodevistb$Taxonomy_A
        #print(nodevistb$Taxonomy_A)
        nodes1$group <- nodevistb$Organism_A
        nodes1$genes <- nodevistb$Gene_A
        colnames(nodes1) <- c("id","tax","group","genes")
        
        nodes2 <- as.data.frame(nodevistb$To_B)
        nodes2$tax <- nodevistb$Taxonomy_B
        nodes2$group <- nodevistb$Organism_B
        nodes2$genes <- nodevistb$Gene_B
        colnames(nodes2) <- c("id","tax","group","genes")
        
        nodev1 <- rbind(nodes1,nodes2)
        nodev1 <- distinct(nodev1)
        nodev2 <- as.data.frame(nodev1)
        nodev2$id <- nodev1$genes
        nodev2$label <- nodev1$genes
        nodev2 <- nodev2%>% mutate(font.size = ifelse(id=="GL068-",34,
                                                      ifelse(id=="GL121-",34,
                                                             ifelse(id=="GL245-",34,
                                                                    ifelse(id=="gusA",27,
                                                                           ifelse(id=="lacZ",27,
                                                                                  ifelse(id=="RUC-GFP",27,
                                                                                         16)))))))
        
        for_symbol <- NULL
        for_symbol$id <- as.data.frame(nodev1$id)
        for_symbol$gene <- as.data.frame(nodev1$genes)
        for_symbol <- as.data.frame(for_symbol)
        colnames(for_symbol) <- c("id","symbol")
        
        nodev2$title <- paste0("Gene :",nodev1$genes,
                               "<br>UniprotKB :",nodev1$id,
                               "<br>Species :",nodev1$group,
                               "<br>Taxonomy :",nodev1$tax)
        
        
        
        nodes <- as.data.frame(nodev2)
        nodes <- subset(nodes, select = -c(genes))
        nodes <- nodes %>% group_by(id) %>% filter (! duplicated(id)) 
        
        edges <- as.data.frame(Edges)
        edges<-edges %>% mutate(color = ifelse(link=="activation","#5ebf87",
                                               ifelse(link=="inhibition","#e0392b",
                                                      ifelse(link=="MI:0194","#d0d676",
                                                             ifelse(link=="MI:0203","#d6b076",
                                                                    ifelse(link=="MI:0217","#8ec4b7",
                                                                           ifelse(link=="MI:0220","#8eb5c4",
                                                                                  ifelse(link=="MI:0403","#998ec4",
                                                                                         ifelse(link=="MI:0407","#bd8ec4",
                                                                                                ifelse(link=="MI:0414","#c48ea3",
                                                                                                       ifelse(link=="MI:0570","#d3815f",
                                                                                                              ifelse(link=="MI:0844","#bfd35f",
                                                                                                                     ifelse(link=="MI:0882","#76d35f",
                                                                                                                            ifelse(link=="MI:0914","#a994bf",
                                                                                                                                   ifelse(link=="MI:0915","#797e85",
                                                                                                                                          ifelse(link=="MI:0945","#5fcfd3",
                                                                                                                                                 ifelse(link=="MI:1126","#5f94d3",
                                                                                                                                                        ifelse(link=="MI:2364","#d35fcc",
                                                                                                                                                               "#9c95a2"))))))))))))))))))
        
        edges<-edges %>% mutate(arrows.to.type = ifelse(link=="activation","arrow",
                                                        ifelse(link=="inhibition","bar", "image")))
        edges<-edges %>% mutate(title = edges$link)
        edges <- distinct(edges)
        edges$id <- row.names(edges)
        
        ledges <-edges[, c("color", "title","arrows.to.type")]
        ledges <- distinct(ledges)
        
        colnames(ledges) <- c("color", "label","arrows")
        
        PPI_network <- reactive({
          visNetwork(nodes,edges,height = "800px", width = "120%") %>%visIgraphLayout(
            layout = "layout_nicely")%>% 
            visGroups(groupname ="Human adenovirus A serotype 12", color =list(highlight ="#C0392B",background = "#80C1F7",border ="#80C1F7")) %>%
            visGroups(groupname = "Human adenovirus A serotype 31", color =list(highlight ="#C0392B",background = "#E57373",border ="#E57373")) %>%
            visGroups(groupname = "Human adenovirus B serotype 3",color =list(highlight ="#C0392B",background = "#F06292",border ="#F06292")) %>%
            visGroups(groupname = "Human adenovirus C serotype 2",color =list(highlight ="#C0392B",background = "#D77DC1",border ="#D77DC1")) %>%
            visGroups(groupname = "Human adenovirus C serotype 5",color =list(highlight ="#C0392B",background = "#BF7785",border ="#BF7785")) %>%
            visGroups(groupname = "Human adenovirus D serotype 9",color =list(highlight ="#C0392B",background = "#C9D36F",border ="#C9D36F")) %>%
            visGroups(groupname = "Murine adenovirus A serotype 1", color =list(highlight ="#C0392B",background = "#BA68C8",border ="#BA68C8")) %>%
            visGroups(groupname = "Vesicular stomatitis Indiana virus (strain Orsay)", color =list(highlight ="#C0392B",background = "#9575CD",border ="#9575CD")) %>%
            visGroups(groupname = "Newcastle disease virus (strain Chicken/United States/B1/48)", color =list(highlight ="#C0392B", background ="#64B5F6",border ="#64B5F6")) %>%
            visGroups(groupname = "Enterovirus", color =list(highlight ="#C0392B",background = "#4FC3F7",border ="#4FC3F7")) %>%
            visGroups(groupname = "Respiratory synctial virus", color =list(highlight ="#C0392B",background = "#4DD0E1",border ="#4DD0E1")) %>%
            visGroups(groupname = "Mumps virus (strain Enders)", color =list(highlight ="#C0392B",background = "#5FE3AD",border ="#5FE3AD")) %>%
            visGroups(groupname = "Measles virus (strain Ichinose-B95a)", color =list(highlight ="#C0392B",background = "#F1948A",border ="#F1948A")) %>%
            visGroups(groupname = "Yaba-like disease virus", color =list(highlight ="#C0392B",background = "#1ABC9C",border ="#1ABC9C")) %>%
            visGroups(groupname = "Poxvirus", color =list(highlight ="#C0392B",background = "#47CFC3",border ="#47CFC3")) %>%
            visGroups(groupname = "Sendai virus (strain Fushimi)", color =list(highlight ="#C0392B",background = "#4A7DF4",border ="#4A7DF4")) %>%
            visGroups(groupname = "Semliki forest virus", color =list(highlight ="#C0392B",background = "#4DB6AC",border ="#4DB6AC")) %>%
            visGroups(groupname = "Homo sapiens", color =list(highlight ="#C0392B",background = "#7986CB",border ="#7986CB")) %>%
            visGroups(groupname = "Myxoma virus", color =list(highlight ="#C0392B",background = "#81C784",border ="#81C784")) %>%
            visGroups(groupname = "Bovine herpesvirus", color =list(highlight ="#C0392B",background = "#8484b5",border ="#8484b5")) %>%
            visGroups(groupname = "Avian orthoreovirus", color =list(highlight ="#C0392B",background = "#1976D2",border ="#1976D2")) %>%
            visGroups(groupname = "Vaccinia virus GLV-1h68", color =list(highlight ="#C0392B",background = "#1976D2",border ="#1976D2")) %>%
            
            
            visGroups(groupname = "insertion", color = "#db3058", shape = "triangle") %>%
            visGroups(groupname = "Deletion", color = "#db3058", shape = "triangleDown") %>%
            visNodes(shape = input$nodeShape, size=45, shadow = TRUE)%>%
            visEdges(arrows = "to",label = edges$link, width =4, shadow = TRUE, smooth = list(enabled = TRUE, roundness= 0.1)) %>%
            visInteraction(multiselect = TRUE, navigationButtons = TRUE, hover = TRUE)%>%
            visExport(type = "pdf")%>%
            visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_selection', nodes.nodes);
                ;}")%>%
            
            visLegend(addEdges = ledges, main = list(text = "Legend",
                                                     style = "font-family:Comic Sans MS;color: #212f3c;font-size:20px;text-align:center;"),
                      
                      width = 0.20, stepY= 53, position = "right")%>%
            visOptions(manipulation = T)
        })
        
        output$downloadNetwork <- downloadHandler(
          filename = function() {
            paste('network-', Sys.Date(), '.html', sep='')
          },
          content = function(con) {
            PPI_network() %>% visSave(con, selfcontained = TRUE)
          })
        
        PPI_network()
      }
    })
  
   
   output$PPIdata <- renderDataTable({
     edges <-e$data
     colnames(edges) <- c("From","Gene1","Desc1","Tax-From","Organism1", "Link","To","Gene2","Desc2", "Tax-To","Organism2")
     
     DT::datatable(edges,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  style="bootstrap4",
                  options = list(
                    dom = 'Bfrtip',
                    scrollY = 400,
                    ordering=TRUE,
                    buttons = list(
                      list(extend = "collection", buttons = c('excel','csv'),
                           text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                      list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_PPI'))
                  ))
   })
  })
  
  
  
  ############# Differential Expression Analysis ##################
  data <- reactiveValues(cleanDeg = NULL,
                         DEGsData = NULL,
                         GEO_gse = NULL,
                         sampleInfo = NULL)
  values <- reactiveValues(
    df_data = data.frame(),
    results = data.frame(),
    gset = data.frame(),
    selected_rows = NULL
  )
  
  
  reactives <- reactiveValues(userData = NULL)
  
  observeEvent(input$dataChoice, {
    DEGsData <- NULL
    cleanDeg <- NULL
    shinyjs::reset("userDif")
    output$ExTab <- renderDataTable({NULL}) 
    output$erciyes <- renderPlotly({NULL})
    output$groupsample <- renderDataTable({NULL})
    
  })
  
  
  observeEvent(input$userDif, {
    req(input$userDif)  # Ensure a file is uploaded
    df <- read.csv(input$userDif$datapath)
    reactives$userData <- df
    
    # Update the pickerInput choices based on the uploaded data columns
    updatePickerInput(session, inputId = 'geneColumn', label = 'Select Gene Column', choices = colnames(df))
    updatePickerInput(session, inputId = 'logFColumn', label = 'Select log2FC Column', choices = colnames(df))
    updatePickerInput(session, inputId = 'pAdColumn', label = 'Select pAdj Column', choices = colnames(df))
  })
  
  user_data <- reactive({
    req(reactives$userData)
    df <- reactives$userData
    df
  })
  
  output$UserData <- renderDataTable({
    user_data()
  })
  
  
  observe({
    if (input$dataChoice != 'Use own data' && input$dataChoice != 'Use GEO data' && input$CancerType > 0 && input$AnalTy == "edgeR" ) {
      expr_df <- read.csv(file = paste0("TCGA_ANOVA_OUTPUT/",input$CancerType, ".csv"))
      
      outVar6 <- reactive({
        vars4 <- as.data.frame(expr_df$gene_name)
        vars4 <- vars4[-1,]
        return(vars4)
      })
      
      output$geneViolin <- renderUI({
        selectInput(inputId = "inVar6", label = "Search Gene", choices = outVar6())
      })
      
    } else if (input$dataChoice != 'Use own data' && input$dataChoice != 'Use GEO data' && input$CancerType > 0 && input$AnalTy == "limma") {
      expr_df <- read.csv(file = paste0("TCGA_LIMMA_OUTPUT/",input$CancerType, ".csv"))
      
      outVar7 <- reactive({
        vars4 <- as.data.frame(expr_df$hgnc_symbol)
        vars4 <- vars4[-1,]
        return(vars4)
      })
      
      output$geneViolin <- renderUI({
        selectInput(inputId = "inVar7", label = "Search Gene", choices = outVar7())
      })
    }
  })
  
  
  observe({
    if (input$dataChoice == 'Use GEO data' && input$compVen == 0) {
      req(input$GEOid)
      print(input$GEOid)
      
      
      withProgress(message = "Running please wait ...",{
      tryCatch({
        gset <- getGEO(input$GEOid, GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL1977", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        data$GEO_gse <- gset
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        exprs(gset) <- log2(exprs(gset))
        values$gset <- gset
        sampleInfo <- pData(gset)
        # Assuming sampleInfo is your data frame
        rownames(sampleInfo) <- NULL
        
        #sampleInfo <- select(sampleInfo, geo_accession, title, source_name_ch1, source_name_ch2)
        # Assuming sampleInfo is your data frame
        sampleInfo <- sampleInfo[, 1:12]
        
        sampleInfo$Group <- NA
        data$sampleInfo <- sampleInfo
        
        
        output$groupsample = renderDataTable({
          # Move the "Group" column to the first position
          data$sampleInfo <- data$sampleInfo[, c("Group", setdiff(names(data$sampleInfo), "Group"))]
          
          datatable(
            data$sampleInfo,
            selection = list(target = 'row', mode = 'multiple'),
            options = list(
              rowCallback = JS(
                'function(row, data, index, api) {
      var group = data[0];  // Assuming the grouping column is at index 4
      if (group !== null) {
        $(row).css("background-color", group);
      }
    }'
              ),
              rownames = FALSE,  # Remove row names
              paging = TRUE,
              scrollY = 500,     # Set the height of the scrollable area
              scroller = TRUE,   # Enable the Scroller extension
              deferRender = TRUE,
              autoWidth = FALSE, # Disable automatic column width calculation
              scrollX = TRUE     # Enable horizontal scrolling
            )
          )
        })
        
        
        
        
        observeEvent(input$groupsample_rows_selected, {
          selected_rows <- isolate(input$groupsample_rows_selected)  # Get the current selection
          current_selected <- values$selected_rows
          
          if (!is.null(current_selected)) {
            # Append the current selection to the existing selected rows
            values$selected_rows <- c(current_selected, selected_rows)
          } else {
            # If it's the first selection, initialize the list
            values$selected_rows <- selected_rows
          }
          
          s <- input$groupsample_rows_selected
          data$sampleInfo[s, "Group"] <- input$assignGroup
          replace_color <- JS(
            'function(row, data, index, api) {
         var group = data[0];  // Assuming the grouping column is at index 4
         if (group !== null) {
           $(row).css("background-color", group);
         }
       }'
          )
          proxy = dataTableProxy('groupsample')
          replaceData(proxy, data$sampleInfo, resetPaging = FALSE)
          addRowSelect = c(
            JS("table.on('select', function(e, dt, type, indexes) {Shiny.setInputValue('groupsample_rows_selected', table.rows({selected: true}).indexes().toArray());});"),
            replace_color
          )
          appendContent(proxy, selected = TRUE, data = addRowSelect)
        })
      
      
      # Display a success notification
      showNotification("Data fetched successfully!", type = "message")
      }, error = function(e) {
        # Handle the error and display an error notification
        showNotification("An error occurred. Please search for GEO ID again.", type = "error")
      })
      
      })
    }
    
  })
  
 
  
  observeEvent(input$runAnalys, {
    if (input$dataChoice == 'Use GEO data' && ((length(unique(data$sampleInfo$Group)) > 3) || (length(unique(data$sampleInfo$Group)) < 3))) {
      
      ask_confirmation(
        "form3",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please determine 2 groups! "
          ))
      ) 
      return()
    }
   
    
    if(input$dataChoice != 'Use own data' && input$dataChoice != 'Use GEO data' && input$CancerType < 0) {
      ask_confirmation(
        "form3",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select the cancer type dataset or upload own data! "
          ))
      ) 
      
    } else {
      
      if(input$dataChoice == 'Use own data'&& input$compVen ==0){
        
        user_df <- as.data.frame(reactives$userData)
        logFC <- (user_df[input$logFColumn])
        GeneCol <- (user_df[input$geneColumn])
        pAdjCol <- (user_df[input$pAdColumn])
        DEGsData <- data.frame(GeneCol, logFC, pAdjCol)
        colnames(DEGsData) <- c("geneSymbol", "LogFC", "P.Value")
        filtered_df <- filter(DEGsData, logFC > input$logFCcut | logFC < -input$logFCcut)
        
        cleanDeg <- filtered_df %>% filter(P.Value < input$qvalcut)
        data$DEGsData <- DEGsData
        
      }
      
      if(input$compVen ==1 && !is.null(reactives$userData) && input$CancerType > 0){
        
        data_file <- reactive({
          file_name <- paste0(input$CancerType)
          file_path <- paste0("TCGA_LIMMA_OUTPUT/", file_name)
          return(file_path)
        })
        
        DEGsData <- read.csv(file= paste0(data_file(),".csv"))
        
        head(DEGsData)
        DEGsData <- DEGsData[, -1]
        colnames(DEGsData) <- c("Ensembl", "logFC", "AverExpr", "t", "Pvalue", "adjPval","B", "geneSymbol")
        head(DEGsData)
        data$DEGsData <- DEGsData
        filtered_df <- filter(DEGsData, logFC > input$logFCcut | logFC < -input$logFCcut)
        
        selectedDeg <- filtered_df %>% filter(adjPval < input$qvalcut)
        
        user_df <- as.data.frame(reactives$userData)
       
        
        logFC <- (user_df[input$logFColumn])
        GeneCol <- (user_df[input$geneColumn])
        pAdjCol <- (user_df[input$pAdColumn])
        DEGsData <- data.frame(GeneCol, logFC, pAdjCol)
        colnames(DEGsData) <- c("geneSymbol", "LogFC", "P.Value")
        data$DEGsData <- DEGsData
        filtered_df <- filter(DEGsData, logFC > input$logFCcut | logFC < -input$logFCcut)
        
        userDeg <- filtered_df %>% filter(P.Value < input$qvalcut)
        
        x <- list(SELECTED_DATA = selectedDeg$geneSymbol, USER_DATA = userDeg$geneSymbol)
        
        venn_plot <- ggvenn(x, stroke_size = 0.5, set_name_size = 4)
        
        output$vennPlot <- renderPlot({
          print(venn_plot)
        })
        
        observeEvent(input$vennPlot_click, {
          click_data <- input$vennPlot_click
          if (is.null(click_data)) return()
          
          x_coord <- click_data$x
          y_coord <- click_data$y
          
          selected_elements <- determine_selected_elements(x_coord, y_coord)
          
          output$elementTable <- renderDT({
            selected_data <- data.frame(Element = selected_elements)
            datatable(selected_data)
          })
        })
        
        determine_selected_elements <- function(x_coord, y_coord) {
          # Determine which region was clicked
          clicked_region <- determine_clicked_region(x_coord, y_coord)
          
          if (clicked_region == "SELECTED_DATA") {
            return(setdiff(selectedDeg$geneSymbol, userDeg$geneSymbol))
          } else if (clicked_region == "USER_DATA") {
            return(setdiff(userDeg$geneSymbol, selectedDeg$geneSymbol))  # Use setdiff the other way around
          } else if (clicked_region == "Intersection") {
            return(intersect(selectedDeg$geneSymbol, userDeg$geneSymbol))
          } else {
            return(NULL)
          }
        }
        
        
        determine_clicked_region <- function(x_coord, y_coord) {
          # Implement logic to determine which region was clicked based on coordinates
          if (x_coord < 0.5 && y_coord < 0.5) {
            return("SELECTED_DATA")
          } else if (x_coord > 0.5 && y_coord < 0.5) {
            return("USER_DATA")
          } else if (y_coord > 0.5) {
            return("Intersection")
          } else {
            return(NULL)
          }
          
        }
      } 
      
      if (input$compVen ==0 && input$dataChoice != 'Use own data' && input$dataChoice != 'Use GEO data' && input$CancerType > 0 && input$AnalTy == "edgeR") {
        data_file <- reactive({
          file_name <- paste0(input$CancerType)
          file_path <- paste0("TCGA_ANOVA_OUTPUT/", file_name)
          print(file_path)
          return(file_path)
          
        })
        
        DEGsData <- read.csv(file= paste0(data_file(),".csv"))
        # Check if column names match the first set of names
        if (all(colnames(DEGsData) %in% c("X", "logFC", "logCPM", "LR", "PValue", "FDR", "gene_name", "gene_type"))) {
          DEGsData <- DEGsData[, -ncol(DEGsData)]
          colnames(DEGsData) <- c("Ensembl", "logFC", "logCPM", "LR", "Pvalue", "FDR", "geneSymbol")
          
        }
        
        # Check if column names match the second set of names
        if (all(colnames(DEGsData) %in% c("X", "genes", "logFC", "logCPM", "LR", "PValue", "FDR", "hgnc_symbol"))) {
          DEGsData <- DEGsData[, -1]
          colnames(DEGsData) <- c("Ensembl", "logFC", "logCPM", "LR", "Pvalue", "FDR", "geneSymbol")
          
        }
        
        data$DEGsData <- DEGsData
        filtered_df <- filter(DEGsData, logFC > input$logFCcut | logFC < -input$logFCcut)
        cleanDeg <- filtered_df %>% filter(FDR < input$qvalcut) %>%
          filter(FDR != 0)
        
      }
      
      if (is.null(reactives$userData) && input$dataChoice != 'Use own data' && input$dataChoice != 'Use GEO data' &&  input$CancerType > 0 && input$AnalTy == "limma" && input$compVen ==0) {
        data_file <- reactive({
          file_name <- paste0(input$CancerType)
          file_path <- paste0("TCGA_LIMMA_OUTPUT/", file_name)
          return(file_path)
        })
        
        DEGsData <- read.csv(file= paste0(data_file(),".csv"))
        
        head(DEGsData)
        DEGsData <- DEGsData[, -1]
        colnames(DEGsData) <- c("Ensembl", "logFC", "AverExpr", "t", "Pvalue", "adjPval","B", "geneSymbol")
        head(DEGsData)
        data$DEGsData <- DEGsData
        filtered_df <- filter(DEGsData, logFC > input$logFCcut | logFC < -input$logFCcut)
        
        cleanDeg <- filtered_df %>% filter(adjPval < input$qvalcut)
        
      }
      
      if (input$dataChoice == 'Use GEO data' && input$compVen == FALSE) {
       
          gset <-values$gset
          # make proper column names to match toptable 
          fvarLabels(gset) <- make.names(fvarLabels(gset))
          
          num_samples <- ncol(gset)
          
          # Initialize gsms with "X" for all rows
          gsms <- rep("X", num_samples)
          
          selected_rows <- values$selected_rows
          
          # Extract the "Group" column from data$sampleInfo
          group_column <- as.character(data$sampleInfo$Group)
          
          # Update the selected rows in gsms based on the selected groups
          gsms[selected_rows] <- group_column[selected_rows]
          
          # Convert it into a factor
          gsms <- as.factor(gsms)
          
          # Extract the levels
          
          groups <- levels(gsms)[levels(gsms) != "X"]
          gsms <- factor(ifelse(gsms == "1", "0", gsms))
          
          # Set the levels of the factor
          levels(gsms) <- c("0", "1", "X")
          
          sel <- which(gsms != "X")
          
          gsms <- gsms[sel]
          gset <- gset[ ,sel]
          
          # log2 transformation
          ex <- exprs(gset)
          qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
          LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0)
          if (LogC) { ex[which(ex <= 0)] <- NaN
          exprs(gset) <- log2(ex) }
          
          # assign samples to groups and set up design matrix
          gs <- factor(gsms)
          
          design <- model.matrix(~gs + 0, gset)
          colnames(design) <- groups
          
          gset <- gset[complete.cases(exprs(gset)), ] # skip missing values
          
          fit <- lmFit(gset, design)  # fit linear model
          
          # set up contrasts of interest and recalculate model coefficients
          cts <- paste(groups[1], groups[2], sep="-")
          cont.matrix <- makeContrasts(contrasts=cts, levels=design)
          fit2 <- contrasts.fit(fit, cont.matrix)
          
          # compute statistics and table of top significant genes
          fit2 <- eBayes(fit2)
          #tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
          full_results <- topTable(fit2, number=Inf)
          #full_results <- tibble::rownames_to_column(full_results,"ID")
          
          
          tT <- subset(full_results, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol"))
          
          # Assuming tT is your data frame
          colnames(tT)[colnames(tT) == "Gene.symbol"] <- "Symbol"
          rownames(tT) <- NULL
          
          # Set the reactive value
          data$DEGsData <- as.data.frame(tT)
          
          filtered_df <- filter(tT, logFC > input$logFCcut | logFC < -input$logFCcut)
          cleanDeg <- filter(filtered_df, adj.P.Val < input$qvalcut)
          
      }
      
      withProgress(message = "Running ...",{
        
        output$ExTab <-  DT::renderDataTable(server=FALSE,{
          
          #cleanDeg<- cleanDeg[,-1]
          data$cleanDeg <- cleanDeg
          cleanDeg <- data.table(cleanDeg)
          DT::datatable(cleanDeg,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  style="bootstrap4",
                  options = list(
                    dom = 'Bfrtip',
                    scrollY = 400,
                    ordering=TRUE,
                    buttons = list(
                      list(extend = "collection", buttons = c('excel','csv'),
                           text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                      list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_DEG_limma'))
                  ))
         
        }) 
      }) 
    }
 
  #___ Volcano Plot ____#
  output$erciyes <- renderPlotly({
    req(input$PloType == "plot1" || input$PloType2 == "plot1" || input$PloType3 == "plot1" )
    req(data$DEGsData)
    
    cleanDeg <- data$DEGsData
    cleanDeg$diffexpressed <- "Insignificant"
    
    if (input$dataChoice == 'Use own data' && input$PloType2 == "plot1") {
     
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      cleanDeg$diffexpressed[cleanDeg$LogFC > input$logFCcut & cleanDeg$P.Value < input$qvalcut] <- "Up-regulated"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      cleanDeg$diffexpressed[cleanDeg$LogFC < -input$logFCcut & cleanDeg$P.Value < input$qvalcut] <- "Down-regulated"
      
      cleanDeg$negLogpAdj <- -log10(cleanDeg$P.Value)
      
      volc <- ggplot(cleanDeg,aes(label = geneSymbol,  x=LogFC, y=negLogpAdj, color=diffexpressed)) +
        geom_point() +
        coord_cartesian() +
        scale_color_manual(values=c("blue", "black", "red")) +
        ylab("-log10 FDR") +
        xlab("logFC")+
        geom_vline(xintercept=c(-input$logFCcut, input$logFCcut), linetype='dashed',col="gray") +
        geom_hline(yintercept= -log10(input$qvalcut), linetype='dashed',col="gray")+
        theme_bw()
      ggplotly(volc)
      
    }
    
    if (input$dataChoice == 'TCGA-GTEX data' && input$PloType == "plot1") {
      if (input$AnalTy == "edgeR") {
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      cleanDeg$diffexpressed[cleanDeg$logFC > input$logFCcut & cleanDeg$FDR < input$qvalcut] <- "Up-regulated"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      cleanDeg$diffexpressed[cleanDeg$logFC < -input$logFCcut & cleanDeg$FDR < input$qvalcut] <- "Down-regulated"
      
      #convert FDR to -log(FDR)
      cleanDeg$negLogpAdj <- -log10(cleanDeg$FDR)
      
      volc <- ggplot(cleanDeg,aes(label = geneSymbol,  x=logFC, y=negLogpAdj, color=diffexpressed)) +
        geom_point() +
        coord_cartesian() +
        scale_color_manual(values=c("blue", "black", "red")) +
        ylab("-log10 FDR") +
        xlab("logFC")+
        geom_vline(xintercept=c(-input$logFCcut, input$logFCcut), linetype='dashed',col="gray") +
        geom_hline(yintercept= -log10(input$qvalcut), linetype='dashed',col="gray")+
        theme_bw()
      ggplotly(volc)
      } else if (input$AnalTy == "limma") {
        
        # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
        cleanDeg$diffexpressed[cleanDeg$logFC > input$logFCcut & cleanDeg$adjPval < input$qvalcut] <- "Up-regulated"
        # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
        cleanDeg$diffexpressed[cleanDeg$logFC < -input$logFCcut & cleanDeg$adjPval < input$qvalcut] <- "Down-regulated"
        
        #convert FDR to -log(FDR)
        cleanDeg$negLogpAdj <- -log10(cleanDeg$adjPval)
        
        volc <- ggplot(cleanDeg,aes(label = geneSymbol,  x=logFC, y=negLogpAdj, color=diffexpressed)) +
          geom_point() +
          coord_cartesian() +
          scale_color_manual(values=c("blue", "black", "red")) +
          ylab("-log10 adj.P.Val") +
          xlab("logFC")+
          geom_vline(xintercept=c(-input$logFCcut, input$logFCcut), linetype='dashed',col="gray") +
          geom_hline(yintercept= -log10(input$qvalcut), linetype='dashed',col="gray")+
          theme_bw()
        ggplotly(volc)
      }
      
    }
    if (input$dataChoice == 'Use GEO data' && input$PloType3 == "plot1") {
      
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      cleanDeg$diffexpressed[cleanDeg$logFC > input$logFCcut & cleanDeg$adj.P.Val < input$qvalcut] <- "Up-regulated"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      cleanDeg$diffexpressed[cleanDeg$logFC < -input$logFCcut & cleanDeg$adj.P.Val < input$qvalcut] <- "Down-regulated"
      
      #convert FDR to -log(FDR)
      cleanDeg$negLogpAdj <- -log10(cleanDeg$adj.P.Val)
      
      volc <- ggplot(cleanDeg,aes(label = Symbol,  x=logFC, y=negLogpAdj, color=diffexpressed)) +
        geom_point() +
        coord_cartesian() +
        scale_color_manual(values=c("blue", "black", "red")) +
        ylab("-log10 adj.P.Val") +
        xlab("logFC")+
        geom_vline(xintercept=c(-input$logFCcut, input$logFCcut), linetype='dashed',col="gray") +
        geom_hline(yintercept= -log10(input$qvalcut), linetype='dashed',col="gray")+
        theme_bw()
      ggplotly(volc)
    }
    if (exists("volc")) {
      return(volc)
    }
  })
  
  
  output$GEO_pca <- renderPlotly({
    req(input$PloType3 == "GEOpca" )
    if (input$dataChoice == 'Use GEO data' && input$compVen == 0) {
      gse <- data$GEO_gse
      sampleInfo <- data$sampleInfo
      
      pca <- prcomp(t(exprs(gse)))
      
        pca_plot <- plot_ly(
          data = cbind(sampleInfo, pca$x),
          x = ~PC1,
          y = ~PC2,
          color = ~Group,
          text = ~paste("Patient", geo_accession),
          type = 'scatter',
          mode = 'markers',
          marker = list(size = 8),
          textposition = "bottom center"
        )
        
        pca_plot <- pca_plot %>% layout(
          xaxis = list(title = "PC1"),
          yaxis = list(title = "PC2"),
          title = "Interactive PCA Plot"
        )
        
        pca_plot
    }
      })
  
  
  
    output$GEO_heatmap <- renderPlotly({
      req(input$PloType3 == "GEOheatmap" )
      if(input$dataChoice == 'Use GEO data' && input$compVen == 0){
      
        #gse <- as.data.frame(data$GEO_gse)
        corMatrix <- cor(exprs(data$GEO_gse), use = "c")
        heatmap <- heatmaply(corMatrix)
        heatmap
      }
      
    })
  
  
  ## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX
  
  output$violin <- renderPlotly({
    if(input$PloType =="plot3" && input$AnalTy == "edgeR"){
      req(input$inVar6)
      
      load(file = paste0("TCGA_DATA/",input$CancerType, ".rda"))
      tcga_data <- data
      expr_df <- read.csv(file = paste0("TCGA_ANOVA_OUTPUT/", input$CancerType, ".csv"))
      expr_df <- expr_df[, -ncol(expr_df)]
      colnames(expr_df) <- c("gene_id", "logFC", "logCPM", "LR", "Pvalue", "FDR", "gene_name")
      # Get the expression values for the selected gene
      col_names <- colnames(tcga_data)
      row_names <- rownames(tcga_data)
      d_mat <- data@assays@data@listData[["fpkm_uq_unstrand"]]
      rownames(d_mat) <- row_names
      colnames(d_mat) <- col_names
      
      # Extract the part before the "." in the row names of data
      data_row_names <- sub("\\..*", "", rownames(d_mat))
      rownames(d_mat) <- data_row_names
      
      # Extract the gene IDs from "expr_df" row names
      deg_gene_ids <- expr_df$gene_id
      
      # Subset "d_mat" based on DEG gene IDs
      sub_d_mat <- subset(d_mat, rownames(d_mat) %in% deg_gene_ids)
      sub_d_mat<- as.data.frame(sub_d_mat)
      sub_d_mat$gene_id <- rownames(sub_d_mat)
      
      
      merged_data <- merge(sub_d_mat, expr_df, by = "gene_id", all.x = TRUE)
      
      cols_to_drop <- c("logFC", "LogCPM","LR","Pvalue","FDR","gene_id")  # Replace column_name1, column_name2, ... with the names of columns you want to drop
      filtered_data <- merged_data[, !(names(merged_data) %in% cols_to_drop)]
      
      # Get the "TP" and "NT" sample labels from tcga_data
      sample_labels <- tcga_data@colData@listData[["shortLetterCode"]]
      
      # Split d_mat into cancer and normal samples
      cancer_samples <- d_mat[, sample_labels == "TP"]
      cancer_samples <- as.data.frame(cancer_samples)
      cancer_samples$gene_id <- rownames(cancer_samples)
      
      selected_expr_df <- expr_df[, c("gene_id", "gene_name"), drop = FALSE]
      cancer_samples_data <- merge(cancer_samples, selected_expr_df, by = "gene_id", all.x = TRUE)
      
      normal_samples <- d_mat[, sample_labels == "NT"]
      normal_samples <- as.data.frame(normal_samples)
      normal_samples$gene_id <- rownames(normal_samples)
      normal_samples_data <- merge(normal_samples, selected_expr_df, by = "gene_id", all.x = TRUE)
      
      # Assuming "gene_symbol_of_interest" is the HGNC symbol (gene symbol) for the gene of interest
      gene_symbol_of_interest <- input$inVar6
      
      # Filter cancer and normal samples for the gene of interest
      gene_cancer_data <- cancer_samples_data %>%
        filter(gene_name == gene_symbol_of_interest)
      
      
      
      gene_normal_data <- normal_samples_data  %>%
        filter(gene_name == gene_symbol_of_interest)
      
      gene_cancer_data$gene_name <- NULL
      gene_normal_data$gene_name <- NULL
      gene_cancer_data$gene_id <- NULL
      gene_normal_data$gene_id <- NULL
      
      gene_cancer_data <- t(gene_cancer_data)
      gene_cancer_data <-as.data.frame(gene_cancer_data)
      
      gene_normal_data <- t(gene_normal_data)
      gene_normal_data <-as.data.frame(gene_normal_data)
      
      fig_normal <- gene_normal_data |>
        plot_ly(
          x = ~rep(1, nrow(gene_normal_data)),
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = ~paste("Sample: ", rownames(gene_normal_data), "<br>Expression: ", V1),
          hoverinfo = "text",  # Show only the custom text on hover
          type = 'violin',
          name = "Normal",
          showlegend = FALSE
        ) |>
        add_boxplot(
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = rownames(gene_normal_data),  # Show row names as labels on hover
          hoverinfo = "text+y",  # Show row names and y-values on hover
          jitter = 0.4,  # Adjust the jitter to control the spread of dots
          marker = list(size = 6, symbol = "dot"),
          name = "Normal",
          showlegend = FALSE,
          hovertemplate = "Sample: %{text}<br>Expression: %{y}<extra></extra>"  # Hover template for the boxplot trace
        ) |>
        add_markers(
          x = ~rep(1, nrow(gene_normal_data)),  # Position the dots at x=1 (centered within the violin)
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = rownames(gene_normal_data),  # Show row names as labels on hover
          hoverinfo = "text",  # Show row names on hover
          marker = list(size = 8, symbol = "dot", opacity = 0.5, line = list(width = 0)),  # Customize the dot appearance
          name = "Normal",
          showlegend = FALSE,
          hovertemplate = "Sample: %{text}<br>Expression: %{y}<extra></extra>"  # Hover template for the dot plot trace
        ) |>
        layout(yaxis = list(zeroline = FALSE, title = "Gene Expression (fpkm_uq_unstrand)"), 
               xaxis = list(
                 title = "Sample",
                 showticklabels = FALSE,  # This removes the tick labels on the x-axis
                 tickvals = "",  # This sets the tick values on the x-axis to empty
                 ticktext = ""   # This sets the tick text on the x-axis to empty,
               ),
               title = input$inVar6
        )
      fig_cancer <- gene_cancer_data |>
        plot_ly(
          x = ~rep(2, nrow(gene_cancer_data)),  # Position the cancer violin plot to the right of the normal violin plot
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = ~paste("Sample: ", rownames(gene_cancer_data), "<br>Expression: ", V1),
          hoverinfo = "text",  # Show only the custom text on hover
          type = 'violin',
          name = "Cancer",
          showlegend = FALSE
        ) |>
        add_boxplot(
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = rownames(gene_cancer_data),  # Show row names as labels on hover
          hoverinfo = "text+y",  # Show row names and y-values on hover
          jitter = 0.4,  # Adjust the jitter to control the spread of dots
          marker = list(size = 6, symbol = "dot"),
          name = "Cancer",
          showlegend = FALSE,
          hovertemplate = "Sample: %{text}<br>Expression: %{y}<extra></extra>"  # Hover template for the boxplot trace
        ) |>
        add_markers(
          x = ~rep(2, nrow(gene_cancer_data)),  # Position the dots at x=2 (next to the normal violin)
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = rownames(gene_cancer_data),  # Show row names as labels on hover
          hoverinfo = "text",  # Show row names on hover
          marker = list(size = 8, symbol = "dot", opacity = 0.5, line = list(width = 0)),  # Customize the dot appearance
          name = "Cancer",
          showlegend = FALSE,
          hovertemplate = "Sample: %{text}<br>Expression: %{y}<extra></extra>"  # Hover template for the dot plot trace
        )
      
      # Combine the two violin plots side by side
      combined_violin_plot <- subplot(
        fig_normal,
        fig_cancer,
        nrows = 1,
        shareY = TRUE,
        margin = 0.05  # Adjust the margin as needed
      )
      
      combined_violin_plot
    }
    
    else if(input$PloType =="plot3" && input$AnalTy == "limma"){
      req(input$inVar7)
      load(file = paste0("TCGA_DATA/",input$CancerType, ".rda"))
      tcga_data <- data
      expr_df <- read.csv(file = paste0("TCGA_LIMMA_OUTPUT/", input$CancerType, ".csv"))
      
      # Get the expression values for the selected gene
      col_names <- colnames(tcga_data)
      row_names <- rownames(tcga_data)
      d_mat <- data@assays@data@listData[["fpkm_uq_unstrand"]]
      rownames(d_mat) <- row_names
      colnames(d_mat) <- col_names
      
      # Extract the part before the "." in the row names of data
      data_row_names <- sub("\\..*", "", rownames(d_mat))
      rownames(d_mat) <- data_row_names
      
      # Extract the gene IDs from "expr_df" row names
      deg_gene_ids <- expr_df$hgnc_symbol
      
      # Subset "d_mat" based on DEG gene IDs
      sub_d_mat <- subset(d_mat, rownames(d_mat) %in% deg_gene_ids)
      sub_d_mat<- as.data.frame(sub_d_mat)
      sub_d_mat$gene_id <- rownames(sub_d_mat)
      merged_data <- merge(sub_d_mat, expr_df, by = "gene_id", all.x = TRUE)
      cols_to_drop <- c("X", "logFC", "AveExpr","t","P.Value","adj.P.Val","B","gene_id")  # Replace column_name1, column_name2, ... with the names of columns you want to drop
      filtered_data <- merged_data[, !(names(merged_data) %in% cols_to_drop)]
      
      # Get the "TP" and "NT" sample labels from tcga_data
      sample_labels <- tcga_data@colData@listData[["shortLetterCode"]]
      
      # Split d_mat into cancer and normal samples
      cancer_samples <- d_mat[, sample_labels == "TP"]
      cancer_samples <- as.data.frame(cancer_samples)
      cancer_samples$gene_id <- rownames(cancer_samples)
      
      selected_expr_df <- expr_df[, c("gene_id", "hgnc_symbol"), drop = FALSE]
      cancer_samples_data <- merge(cancer_samples, selected_expr_df, by = "gene_id", all.x = TRUE)
      
      normal_samples <- d_mat[, sample_labels == "NT"]
      normal_samples <- as.data.frame(normal_samples)
      normal_samples$gene_id <- rownames(normal_samples)
      normal_samples_data <- merge(normal_samples, selected_expr_df, by = "gene_id", all.x = TRUE)
      
      # Assuming "gene_symbol_of_interest" is the HGNC symbol (gene symbol) for the gene of interest
      gene_symbol_of_interest <- input$inVar7
      
      # Filter cancer and normal samples for the gene of interest
      gene_cancer_data <- cancer_samples_data %>%
        filter(hgnc_symbol == gene_symbol_of_interest)
      
      gene_normal_data <- normal_samples_data  %>%
        filter(hgnc_symbol == gene_symbol_of_interest)
      
      gene_cancer_data$hgnc_symbol <- NULL
      gene_normal_data$hgnc_symbol <- NULL
      gene_cancer_data$gene_id <- NULL
      gene_normal_data$gene_id <- NULL
      
      gene_cancer_data <- t(gene_cancer_data)
      gene_cancer_data <-as.data.frame(gene_cancer_data)
      
      gene_normal_data <- t(gene_normal_data)
      gene_normal_data <-as.data.frame(gene_normal_data)
      
      # Assuming you have already loaded the necessary libraries (e.g., dplyr, plotly)
      fig_normal <- gene_normal_data |>
        plot_ly(
          x = ~rep(1, nrow(gene_normal_data)),
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = ~paste("Sample: ", rownames(gene_normal_data), "<br>Expression: ", V1),
          hoverinfo = "text",  # Show only the custom text on hover
          type = 'violin',
          name = "Normal",
          showlegend = FALSE
        ) |>
        add_boxplot(
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = rownames(gene_normal_data),  # Show row names as labels on hover
          hoverinfo = "text+y",  # Show row names and y-values on hover
          jitter = 0.4,  # Adjust the jitter to control the spread of dots
          marker = list(size = 6, symbol = "dot"),
          name = "Normal",
          showlegend = FALSE,
          hovertemplate = "Sample: %{text}<br>Expression: %{y}<extra></extra>"  # Hover template for the boxplot trace
        ) |>
        add_markers(
          x = ~rep(1, nrow(gene_normal_data)),  # Position the dots at x=1 (centered within the violin)
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = rownames(gene_normal_data),  # Show row names as labels on hover
          hoverinfo = "text",  # Show row names on hover
          marker = list(size = 8, symbol = "dot", opacity = 0.5, line = list(width = 0)),  # Customize the dot appearance
          name = "Normal",
          showlegend = FALSE,
          hovertemplate = "Sample: %{text}<br>Expression: %{y}<extra></extra>"  # Hover template for the dot plot trace
        ) |>
        layout(yaxis = list(zeroline = FALSE, title = "Gene Expression (fpkm_uq_unstrand)"), 
               xaxis = list(
                 title = "Sample",
                 showticklabels = FALSE,  # This removes the tick labels on the x-axis
                 tickvals = "",  # This sets the tick values on the x-axis to empty
                 ticktext = ""   # This sets the tick text on the x-axis to empty,
               ),
               title = input$inVar4
        )
      fig_cancer <- gene_cancer_data |>
        plot_ly(
          x = ~rep(2, nrow(gene_cancer_data)),  # Position the cancer violin plot to the right of the normal violin plot
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = ~paste("Sample: ", rownames(gene_cancer_data), "<br>Expression: ", V1),
          hoverinfo = "text",  # Show only the custom text on hover
          type = 'violin',
          name = "Cancer",
          showlegend = FALSE
        ) |>
        add_boxplot(
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = rownames(gene_cancer_data),  # Show row names as labels on hover
          hoverinfo = "text+y",  # Show row names and y-values on hover
          jitter = 0.4,  # Adjust the jitter to control the spread of dots
          marker = list(size = 6, symbol = "dot"),
          name = "Cancer",
          showlegend = FALSE,
          hovertemplate = "Sample: %{text}<br>Expression: %{y}<extra></extra>"  # Hover template for the boxplot trace
        ) |>
        add_markers(
          x = ~rep(2, nrow(gene_cancer_data)),  # Position the dots at x=2 (next to the normal violin)
          y = ~V1,  # Assuming the gene expression data is in the first column
          text = rownames(gene_cancer_data),  # Show row names as labels on hover
          hoverinfo = "text",  # Show row names on hover
          marker = list(size = 8, symbol = "dot", opacity = 0.5, line = list(width = 0)),  # Customize the dot appearance
          name = "Cancer",
          showlegend = FALSE,
          hovertemplate = "Sample: %{text}<br>Expression: %{y}<extra></extra>"  # Hover template for the dot plot trace
        )
      
      # Combine the two violin plots side by side
      combined_violin_plot <- subplot(
        fig_normal,
        fig_cancer,
        nrows = 1,
        shareY = TRUE,
        margin = 0.05  # Adjust the margin as needed
      )
      
      combined_violin_plot
    }
  })
  
  
  output$PCAnalysis <- renderPlotly({
    if(input$PloType =="plot4" && input$AnalTy == "edgeR"){
    load(file = paste0("TCGA_DATA/", input$CancerType, ".rda"))
    tcga_data <- data
    expr_df <- read.csv(file = paste0("TCGA_ANOVA_OUTPUT/",  input$CancerType, ".csv"))
    colnames(expr_df)[colnames(expr_df) == "X"] <- "gene_id"
    
    # Get the expression values for the selected gene
    col_names <- colnames(tcga_data)
    row_names <- rownames(tcga_data)
    d_mat <- data@assays@data@listData[["fpkm_uq_unstrand"]]
    rownames(d_mat) <- row_names
    colnames(d_mat) <- col_names
    
    # Extract the part before the "." in the row names of data
    data_row_names <- sub("\\..*", "", rownames(d_mat))
    rownames(d_mat) <- data_row_names
    
    cols_to_drop <- c("logFC", "LogCPM","LR","Pvalue","FDR","gene_id")
    filtered_data <- d_mat[, !(names(d_mat) %in% cols_to_drop)]
    
    
    # Get the "TP" and "NT" sample labels from tcga_data
    sample_labels <- colData(tcga_data)$shortLetterCode
    # Split d_mat into cancer and normal samples
    cancer_samples <- d_mat[, sample_labels == "TP"]
    cancer_samples <- as.data.frame(cancer_samples)
    cancer_samples$gene_id <- rownames(cancer_samples)
    
    cancer_samples_data <- merge(cancer_samples, expr_df, by = "gene_id", all.y = TRUE)
    
    normal_samples <- d_mat[, sample_labels == "NT"]
    
    normal_samples <- as.data.frame(normal_samples)
    normal_samples$gene_id <- rownames(normal_samples)
    normal_samples_data <- merge(normal_samples, expr_df, by = "gene_id", all.y = TRUE)
    all_samples_data <- merge(normal_samples_data, cancer_samples_data, by = "gene_id", all = TRUE)
    # Sort the genes based on a relevant metric, such as log fold change or p-value
    sorted_genes <- all_samples_data[order(all_samples_data$logFC.x, decreasing = TRUE), ]
    
    # Select the top 200 genes
    
    top_200_genes <- head(sorted_genes, 200)
    top_200_genes$Regulation <- "Up"
    down_200_genes <- tail(sorted_genes, 200)
    down_200_genes$Regulation <- "Down"
    selected_genes <- rbind(top_200_genes, down_200_genes)
   
    cols_to_drop2 <- c("gene_id","gene_name.y","logFC.x","logFC.y","logCPM.y",
                       "logCPM.x","LR.y","PValue.x","PValue.y","FDR.y","FDR.x",
                       "LR.x","gene_type.y","gene_id.y","gene_name.y","gene_name.x","gene_type.x")  # Replace column_name1, column_name2, ... with the names of columns you want to drop
    filtered_400_genes <- selected_genes[, !(names(selected_genes) %in% cols_to_drop2)]
    symbol <- as.data.frame(selected_genes$hgnc_symbol.y)
    rownames(filtered_400_genes) <- filtered_400_genes$gene_id
    # Drop the gene_id column
    filtered_400_genes$gene_id <- NULL
    filtered_400_genes$hgnc_symbol.y <- NULL
    #filtered_400_genes_transposed <- t(filtered_400_genes)
    filtered_400_genes_transposed <- filtered_400_genes
    # Create a "Sample" column from row names
    Sample <- rownames(filtered_400_genes_transposed)
    
    # Add the "Sample" column to the transposed dataframe
    filtered_400_genes_transposed <- as.data.frame(filtered_400_genes_transposed) 
    filtered_400_genes_transposed$Sample<- Sample
    
    sample_column_names <- c(colnames(normal_samples), colnames(cancer_samples))
    
    # Create a data frame to hold the sample information
    sample_info <- data.frame(
      Sample = sample_column_names,
      Type = ifelse(sample_column_names %in% colnames(normal_samples), "Normal", "Cancer")
    )
    
    filtered_400_genes_transposed <- merge(filtered_400_genes_transposed, sample_info, by = "Sample", all.x = TRUE)
    rownames(filtered_400_genes_transposed) <- filtered_400_genes_transposed$Sample
    filtered_400_genes_transposed$Sample<- NULL
    
    
    X <- subset(filtered_400_genes_transposed, select = -c(Type))
    
    X <- subset(X, select = -c(Regulation))
    
    prin_comp <- prcomp(X, rank. = 3)
    
    components <- prin_comp[["x"]]
    components <- data.frame(components)
    components$PC2 <- -components$PC2
    components$PC3 <- -components$PC3
    components = cbind(components, filtered_400_genes_transposed$Type)
    
    tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
    tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)
    
    fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3,
                   color = filtered_400_genes_transposed$Regulation, colors = c('#EF553B', '#636EFA')) %>%
      add_markers(size = 10,text= rownames(X),marker = list(size = 10)) %>%
      layout(
        scene = list(bgcolor = "#e5ecf6",
                     width = 400,  # Adjust the width as needed
                     height = 400  # Adjust the height as needed)
        ))
    fig
    }
    else if(input$PloType =="plot4" && input$AnalTy == "limma"){
      load(file = paste0("TCGA_DATA/", input$CancerType, ".rda"))
      tcga_data <- data
      expr_df <- read.csv(file = paste0("TCGA_LIMMA_OUTPUT/",  input$CancerType, ".csv"))
      # Get the expression values for the selected gene
      col_names <- colnames(tcga_data)
      row_names <- rownames(tcga_data)
      d_mat <- data@assays@data@listData[["fpkm_uq_unstrand"]]
      rownames(d_mat) <- row_names
      colnames(d_mat) <- col_names
      
      # Extract the part before the "." in the row names of data
      data_row_names <- sub("\\..*", "", rownames(d_mat))
      rownames(d_mat) <- data_row_names
      
      cols_to_drop <- c("X", "logFC", "AveExpr","t","P.Value","adj.P.Val","B","gene_id")  # Replace column_name1, column_name2, ... with the names of columns you want to drop
      filtered_data <- d_mat[, !(names(d_mat) %in% cols_to_drop)]
      
      # Get the "TP" and "NT" sample labels from tcga_data
      sample_labels <- colData(tcga_data)$shortLetterCode
      # Split d_mat into cancer and normal samples
      cancer_samples <- d_mat[, sample_labels == "TP"]
      cancer_samples <- as.data.frame(cancer_samples)
      cancer_samples$gene_id <- rownames(cancer_samples)
      
      cancer_samples_data <- merge(cancer_samples, expr_df, by = "gene_id", all.y = TRUE)
      
      normal_samples <- d_mat[, sample_labels == "NT"]
      normal_samples <- as.data.frame(normal_samples)
      normal_samples$gene_id <- rownames(normal_samples)
      normal_samples_data <- merge(normal_samples, expr_df, by = "gene_id", all.y = TRUE)
      all_samples_data <- merge(normal_samples_data, cancer_samples_data, by = "gene_id", all = TRUE)
      # Sort the genes based on a relevant metric, such as log fold change or p-value
      sorted_genes <- all_samples_data[order(all_samples_data$logFC.x, decreasing = TRUE), ]
      
      # Select the top 200 genes
      
      top_200_genes <- head(sorted_genes, 200)
      top_200_genes$Regulation <- "Up"
      down_200_genes <- tail(sorted_genes, 200)
      down_200_genes$Regulation <- "Down"
      selected_genes <- rbind(top_200_genes, down_200_genes)
      
      cols_to_drop2 <- c("X.x","X.y","logFC.x","logFC.y", "AveExpr.x","AveExpr.y",
                         "t.x","t.y","P.Value.x","P.Value.y","adj.P.Val.x","adj.P.Val.y",
                         "B.x","B.y","gene_id.y","hgnc_symbol.x")  # Replace column_name1, column_name2, ... with the names of columns you want to drop
      filtered_400_genes <- selected_genes[, !(names(selected_genes) %in% cols_to_drop2)]
      symbol <- as.data.frame(selected_genes$hgnc_symbol.y)
      rownames(filtered_400_genes) <- filtered_400_genes$gene_id
      # Drop the gene_id column
      filtered_400_genes$gene_id <- NULL
      filtered_400_genes$hgnc_symbol.y <- NULL
      #filtered_400_genes_transposed <- t(filtered_400_genes)
      filtered_400_genes_transposed <- filtered_400_genes
      # Create a "Sample" column from row names
      Sample <- rownames(filtered_400_genes_transposed)
      
      # Add the "Sample" column to the transposed dataframe
      filtered_400_genes_transposed <- as.data.frame(filtered_400_genes_transposed) 
      filtered_400_genes_transposed$Sample<- Sample
      
      sample_column_names <- c(colnames(normal_samples), colnames(cancer_samples))
      
      # Create a data frame to hold the sample information
      sample_info <- data.frame(
        Sample = sample_column_names,
        Type = ifelse(sample_column_names %in% colnames(normal_samples), "Normal", "Cancer")
      )
      
      filtered_400_genes_transposed <- merge(filtered_400_genes_transposed, sample_info, by = "Sample", all.x = TRUE)
      rownames(filtered_400_genes_transposed) <- filtered_400_genes_transposed$Sample
      filtered_400_genes_transposed$Sample<- NULL
      
      X <- subset(filtered_400_genes_transposed, select = -c(Type))
      X <- subset(X, select = -c(Regulation))
      
      prin_comp <- prcomp(X, rank. = 3)
      
      components <- prin_comp[["x"]]
      components <- data.frame(components)
      components$PC2 <- -components$PC2
      components$PC3 <- -components$PC3
      components = cbind(components, filtered_400_genes_transposed$Type)
      
      tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
      tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)
      
      fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3,
                     color = filtered_400_genes_transposed$Regulation, colors = c('#EF553B', '#636EFA')) %>%
        add_markers(size = 10,text= rownames(X),marker = list(size = 10)) %>%
        layout(
          scene = list(bgcolor = "#e5ecf6",
                       width = 400,  # Adjust the width as needed
                       height = 400  # Adjust the height as needed)
          ))
      fig
      
    }
  })
  })
  
  
  ##############################################################################
  ############### Metabolic Analysis ###########################################
  observeEvent(input$runmetabolic,{
    if((input$userinput ==FALSE && input$typeForMetabolic< 0)){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please make sure you select the parameters completely and correctly."
          ))
      )}else {
        
        withProgress(message = "Running please wait ...",{
        conf <- config::get(file="config.yml", use_parent = FALSE)
        
        ######## function
        lazyReadRDS <- function(name, path, envir=.GlobalEnv) {
          if (!name %in% ls(envir=envir)) {
            message(paste0("No ", name, ", loading"))
            if (startsWith(path, "http://") || startsWith(path, "https://")) {
              res <- readRDS(url(path))
            } else {
              res <- readRDS(path)
            }
            assign(name, res, envir=envir)
            message("Done")
            return(res)
          } else {
            return(get(name, envir=envir))
          }
        }
        
        normalizeName <- function(x) {
          gsub("[^a-z0-9]", "", tolower(x))
        }
        
        necessary.gene.de.fields <- c("ID", "pval", "log2FC", "baseMean")
        
        necessary.met.de.fields <- c("ID", "pval", "log2FC")
        
        vector2html <- function(v) {
          paste0("<ul>\n",
                 paste("<li>", names(v), ": ", v, "</li>\n", collapse=""),
                 "</ul>\n")
        }
        #######################
        networks <- list(
          "kegg" = "network.kegg",
          "rhea" = "network.rhea",
          "lipidomic" = "network.rhea.lipids"
        )
        
        networkPaths <- list(
          "kegg" = conf$path.to.kegg.network,
          "rhea" = conf$path.to.rhea.network,
          "lipidomic" = conf$path.to.lipid.network
        )
        
        annotations <- list(
          "mmu" = "org.Mm.eg.gatom.anno",
          "hsa" = "org.Hs.eg.gatom.anno",
          "ath" = "org.At.tair.gatom.anno",
          "sce" = "org.Sc.sgd.gatom.anno"
        )
        
        annotationPaths <- list(
          "mmu" = conf$path.to.org.Mm.eg.gatom.anno,
          "hsa" = conf$path.to.org.Hs.eg.gatom.anno,
          "ath" = conf$path.to.org.At.tair.gatom.anno,
          "sce" = conf$path.to.org.Sc.sgd.gatom.anno
        )
        
        annotationsRhea <- list(
          "mmu" = "gene2reaction.rhea.mmu.eg",
          "hsa" = "gene2reaction.rhea.hsa.eg",
          "ath" = "gene2reaction.rhea.ath.eg",
          "sce" = "gene2reaction.rhea.sce.eg"
        )
        
        annotationRheaPaths <- list(
          "mmu" = conf$path.to.gene2reaction.rhea.mmu.eg,
          "hsa" = conf$path.to.gene2reaction.rhea.hsa.eg,
          "ath" = conf$path.to.gene2reaction.rhea.ath.eg,
          "sce" = conf$path.to.gene2reaction.rhea.sce.eg
        )
        
        getNetwork <- reactive({
          res <- lazyReadRDS(networks[[input$network]],
                             path = networkPaths[[input$network]])
          res
        })
        
        getAnnotation <- reactive({
          res <- lazyReadRDS(annotations[[input$organism]],
                             path=annotationPaths[[input$organism]])
          res
        })
        
        if (input$userinput ==TRUE && input$typeForMetabolic< 0) {
          geneDEInputRaw <- reactive({
            if (!is.null(input$geneDE)) {
              
              inFile2 <- input$geneDE
              print(inFile2)
              if (is.null(inFile2)) return(NULL)   
              dff <- read.csv(inFile2$datapath,header = TRUE)
              print(head(dff))
              dff
            }
            
            
            if (!is.null(input$geneDE)) {
              
              path <- input$geneDE$datapath
              deName <- input$geneDE$name
            }
            
            
            if (grepl("xlsx?$", input$geneDE$name)){
              res <- read.xlsx(path)
            } else {
              res <- fread(path, colClasses="character")
            }
            attr(res, "name") <- deName
            
            res
          })
          
          metDEInputRaw <- reactive({
            
            if (is.null(input$metDE)) {
              # User has not uploaded a file yet
              return(NULL)
            }
            
            if (!is.null(input$gmetDE)) {
              
              inFile3 <- input$metDE
              print(inFile3)
              if (is.null(inFile3)) return(NULL)   
              dff2 <- read.csv(inFile3$datapath,header = TRUE)
              print(head(dff2))
              dff2
            }
            
            
            if (grepl("xlsx?$", input$metDE$name)){
              res <- read.xlsx(input$metDE$datapath)
            } else {
              res <- fread(input$metDE$datapath, colClasses="character")
            }
            attr(res, "name") <- input$metDE$name
            
            res
          })
        }
        
        if (input$userinput ==FALSE && input$typeForMetabolic> 0) {
          geneDEInputRaw <- reactive({
            genedata <- read.csv(file= paste0("TCGA_LIMMA_OUTPUT/",input$typeForMetabolic, ".csv"))
            
            genedata <- as.data.frame(genedata)
            genedata <- genedata[c('hgnc_symbol', 'logFC','P.Value')]
            print(head(genedata))
            colnames(genedata) <- c("Symbol","logFC","PValue")
            genedata
            
          })
          
          metDEInputRaw <- reactive({
            
            if (is.null(input$metDE)) {
              # User has not uploaded a file yet
              return(NULL)
            }
          })
          
        }
        
        geneDEMeta <- reactive({
          gene.de.raw <- geneDEInputRaw()
          
          if (is.null(gene.de.raw)) {
            return(NULL)
          }
          
          org.gatom.anno <- getAnnotation()
          gene.de.meta <- getGeneDEMeta(gene.de.raw, org.gatom.anno)
          gene.de.meta
        })
        
        geneDEInput <- reactive({
          gene.de.raw <- geneDEInputRaw()
          
          if (is.null(gene.de.raw)) {
            return(NULL)
          }
          
          org.gatom.anno <- getAnnotation()
          gene.de.meta <- geneDEMeta()
          
          res <- prepareDE(gene.de.raw, gene.de.meta)
          # res[, signalRank := NULL] # todo
          
          if (!all(necessary.gene.de.fields %in% names(res))) {
            
            stop(paste0("Genomic differential expression data should contain at least these fields: ",
                        paste(necessary.gene.de.fields, collapse=", ")))
          }
          
          attr(res, "name") <- attr(gene.de.raw, "name")
          
          res
        })
        
        geneIdsType <- reactive({
          gene.de.meta <- geneDEMeta()
          
          if (is.null(gene.de.meta)) {
            return(NULL)
          }
          
          res <- gene.de.meta$idType
          
          if (length(res) != 1) {
            stop("Can't determine type of IDs for genes")
          }
          res
        })
        
        getMetDb <- reactive({
          met.de.raw <- metDEInputRaw()
          if (is.null(met.de.raw)) {
            return(NULL)
          }
          
          if ((input$network == "lipidomic")) {
            met.lipid.db <- lazyReadRDS(name = "met.lipid.db",
                                        path = conf$path.to.met.lipid.db)
            met.db <- met.lipid.db
          } else if (input$network == "kegg"){
            met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                       path = conf$path.to.met.kegg.db)
            met.db <- met.kegg.db
          } else {
            met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                       path = conf$path.to.met.rhea.db)
            met.db <- met.rhea.db
          }
          
          met.db
        })
        
        metDEMeta <- reactive({
          met.de.raw <- metDEInputRaw()
          if (is.null(met.de.raw)) {
            return(NULL)
          }
          
          met.db <- getMetDb()
          
          met.de.meta <- getMetDEMeta(met.de.raw, met.db)
          met.de.meta
        })
        
        metDEInput <- reactive({
          met.de.raw <- metDEInputRaw()
          if (is.null(met.de.raw)) {
            return(NULL)
          }
          
          met.db <- getMetDb()
          met.de.meta <- metDEMeta()
          res <- prepareDE(met.de.raw, met.de.meta)
          
          if (!all(necessary.met.de.fields %in% names(res))) {
            
            stop(paste0("Metabolic differential expression data should contain at least these fields: ",
                        paste(necessary.met.de.fields, collapse=", ")))
          }
          
          attr(res, "name") <- attr(met.de.raw, "name")
          res
          
        })
        
        metIdsType <- reactive({
          met.de.meta <- metDEMeta()
          if (is.null(met.de.meta)) {
            return(NULL)
          }
          
          res <- met.de.meta$idType
          
          if (length(res) != 1) {
            stop("Can't determine type of IDs for metabolites")
          }
          res
        })
        
        metDETag <- reactive({
          metData <- metDEInput()
          tag <- attr(metData, "name")
          tag <- gsub("\\.gz$", "", tag)
          tag <- gsub("\\.([ct]sv|txt|xlsx)$", "", tag)
          tag
        })
        
        experimentTag <- reactive({
          geneData <- geneDEInput()
          metData <- metDEInput()
          name <- if (!is.null(geneData)) attr(geneData, "name") else attr(metData, "name")
          tag <- name
          tag <- gsub("\\.gz$", "", tag)
          tag <- gsub("\\.([ct]sv|txt|xlsx)$", "", tag)
          tag
        })
        
        gInput <- reactive({
          org.gatom.anno <- getAnnotation()
          
          if(!is.null(input$metDE)){
            met.de <- isolate(metDEInput())
          }else{
            gene.de <- isolate(geneDEInput())
            
          }
          
          network <- isolate(getNetwork())
          gene.ids <- isolate(geneIdsType())
          met.ids <- isolate(metIdsType())
          tag <- isolate(experimentTag())
          
          tryCatch({
            
            topology <- isolate(input$nodesAs)
            org.gatom.anno <- isolate(getAnnotation())
            keepReactionsWithoutEnzymes <- FALSE
            gene2reaction.extra <- NULL
            
            
            if ((isolate(input$network) == "lipidomic")) {
              met.lipid.db <- lazyReadRDS(name = "met.lipid.db",
                                          path = conf$path.to.met.lipid.db)
              met.db <- met.lipid.db
              
              gene2reaction.extra <- (fread(annotationRheaPaths[[input$organism]], colClasses="character"))[gene != "-"]
              
              if (met.ids == "Species") {
                met.de$SpecialSpeciesLabelColumn <- met.de$ID
              }
              
              topology <- "metabolites"
              keepReactionsWithoutEnzymes <- TRUE
            } else if (isolate(input$network) == "kegg"){
              met.kegg.db <- lazyReadRDS(name = "met.kegg.db",
                                         path = conf$path.to.met.kegg.db)
              met.db <- met.kegg.db
            } else {
              met.rhea.db <- lazyReadRDS(name = "met.rhea.db",
                                         path = conf$path.to.met.rhea.db)
              met.db <- met.rhea.db
            }
            
            if((is.null(input$metDE) && input$userinput ==FALSE && input$typeForMetabolic> 0) || (is.null(input$metDE) && input$userinput ==TRUE && input$typeForMetabolic< 0) ){
              
              g <- makeMetabolicGraph(network=network,
                                      topology=topology,
                                      org.gatom.anno=org.gatom.anno,
                                      gene.de=gene.de,
                                      met.db=met.db,
                                      met.de=NULL,
                                      #met.to.filter=fread(system.file("mets2mask.lst",
                                       #                               package="gatom"))$ID,
                                      keepReactionsWithoutEnzymes = keepReactionsWithoutEnzymes,
                                      gene2reaction.extra = gene2reaction.extra)
              
            }
            else{
              g <- makeMetabolicGraph(network=network,
                                      topology=topology,
                                      org.gatom.anno=org.gatom.anno,
                                      gene.de=gene.de,
                                      met.db=met.db,
                                      met.de=met.de,
                                      #met.to.filter=fread(system.file("mets2mask.lst",
                                       #                               package="gatom"))$ID,
                                      keepReactionsWithoutEnzymes = keepReactionsWithoutEnzymes,
                                      gene2reaction.extra = gene2reaction.extra)
              
            }
            
            attr(g, "tag") <- tag
            g$organism <- isolate(input$organism)
            g
            
          })
        })
        
        g <- isolate(gInput())
        
        if((is.null(input$metDE) && input$userinput ==FALSE && input$typeForMetabolic> 0) || (is.null(input$metDE) && input$userinput ==TRUE && input$typeForMetabolic< 0) ){
          gs <- scoreGraph(g, k.gene = input$kgene, k.met = NULL)
        }
        
        else {
          gs <- scoreGraph(g, k.gene = input$kgene, k.met = input$kgene)
          
        }
        
        solver <- rnc_solver()
        res <- solve_mwcsp(solver, gs)
        m <- res$graph
        
        #print(head(E(m)$label))
        #print(head(V(m)$label))
        data <- as_long_data_frame(m)
        data <- data[!duplicated(data[, c("from", "to", "label")]), ]
        data <- distinct(data)
        
        # Get unique labels from "from_label" and "to_label" columns
        node_labels <- unique(c(as.character(data$from_label),as.character(data$to_label)))
        node_labels <- as.data.frame(node_labels)
        nodes_df <- distinct(node_labels)
        colnames(nodes_df) <- "id"
        nodes_df$label <- nodes_df$id
        
        # Create edges data frame
        edges_df <- data.frame(from = data$from_label,
                               to = data$to_label,
                               label = data$label,
                               log2FC = data$log2FC)  # Include log2FC column
        
        edges_df <- distinct(edges_df)
        
        # Set color and width of edges based on log2FC
        edges_df$color <- ifelse(edges_df$log2FC > 0, "green", "red")  # Set color based on log2FC
        edges_df$width <- abs(edges_df$log2FC)*3  # Set width based on absolute value of log2FC
        
      
        m.ext <- addHighlyExpressedEdges(m, gs)
        m1 <- connectAtomsInsideMetabolite(m.ext)
        m2 <- collapseAtomsIntoMetabolites(m.ext)
       
        
        
        # Render the network
        output$metabolicnetwork <- renderVisNetwork({
          Metabolic_network <- reactive({
            visNetwork(nodes_df, edges_df, main = "Network Visualization") %>%
              visNodes(shape = "dot" ,color = "#8CBEF3", size = 15 ) %>%
              visEdges(color = "color", width = "width") %>%
              visOptions(manipulation = TRUE)
          })
          
          output$downloadNetwork <- downloadHandler(
            filename = function() {
              paste('network-', Sys.Date(), '.html', sep='')
            },
            content = function(con) {
              Metabolic_network() %>% visSave(con, selfcontained = TRUE)
            })
          
          Metabolic_network() 
        })
        
        
        
        output$pathwaysTable <- DT::renderDataTable(server=FALSE,{
          org.gatom.anno <- getAnnotation()
          
          foraRes <- fgsea::fora(pathways=org.gatom.anno$pathways,
                                 genes=E(m)$gene,
                                 universe=unique(E(g)$gene),
                                 minSize=5)
          res <- foraRes[padj < 0.05]
          DT::datatable(res,
                        extensions = c('Buttons','Responsive','Scroller'), 
                        style="bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering=TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel','csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E-gatom'))
                        ),escape = FALSE)
        })
        
        
        output$metabolicgenedf <- DT::renderDataTable(server=FALSE,{
          df_gene.de <- isolate(geneDEInput())
          DT::datatable(df_gene.de,
                        extensions = c('Buttons','Responsive','Scroller'), 
                        style="bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering=TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel','csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'tt'))
                        ),rownames = FALSE,escape = FALSE)
        })
        
        output$reactionTable <- DT::renderDataTable(server=FALSE,{
          
          DT::datatable(data,
                        extensions = c('Buttons','Responsive','Scroller'), 
                        style="bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering=TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel','csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'tt'))
                        ),rownames = FALSE,escape = FALSE)
        })
        })
      }
    
  })
  
  
  
  
  ##############################################################################
  ############### Metabolic Flux Analysis with CNApy ###########################
  
  # Render the DataTable when "Create SBML" is selected
  observeEvent(input$selectSBML, {
    if (input$selectSBML == "Create SBML") {
      forsbml_df <- reactiveValues(
        forsbml = data.frame(
          type = character(),
          name = character(),
          equation = character(),
          stringsAsFactors = FALSE
        )
      )
      
      
      # Load the existing data from the CSV file if it exists
      
        observeEvent(input$sbmlCsv, {
          req(input$sbmlCsv)  # Ensure a file is uploaded
          df <- read.csv(input$sbmlCsv$datapath)
          forsbml_df$forsbml <- df})
        
        observe({
          req(input$sbmlCsv)  # Ensure a file is uploaded
          
          # Define the path where you want to save the file
          output_path <- normalizePath("out/metabolites.csv")
          
          # Save the uploaded file to the specified path
          if (!is.null(input$sbmlCsv$datapath)) {
            file.copy(input$sbmlCsv$datapath, output_path, overwrite = TRUE)
          }
        })
        
     
      
      output$ediTableSbml <- DT::renderDataTable(server = FALSE, {
        datatable(
          forsbml_df$forsbml,
          editable = TRUE,
          style = "bootstrap4",
          rownames = FALSE,
          extensions = c('Buttons', 'Responsive', 'Scroller'),
          options = list(
            dom = 'Bfrtip',
            scrollY = "400px",  # Change the scrollY value to your desired height in pixels
            paging = FALSE,  # Disable paging
            ordering = TRUE,
            
            buttons = list(
              list(extend = "collection", buttons = c('excel', 'csv'),
                   text = '<span class="glyphicon glyphicon-download-alt"></span>'),
              list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_for_sbml')
            )
          )
        )
      })
     
      
     
      
      observeEvent(input$createSBML,{
        
        if(is.null(input$sbmlCsv) && !(file.exists("model.xml") || file.exists("model.sbml")) && input$selectSBML =='Create SBML'){
          ask_confirmation(
            "form1",
            title = " ",
            type = "info",
            btn_labels = c("Cancel", "OK"),
            allowEscapeKey = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            text =  
              div(HTML("Oops!", "When creating the sbml file, please make sure that the required csv file contains the correct reactions and interactions."
              ))
          )}
        
        
        
        withProgress(message = "Running please wait ...",{
            
      # Run Python script in a separate process
      command <- "/opt/conda/bin/python"
      script <- "Visualization_D3flux.py"
      
      system2(command, args = script)
      
      
      reactives <- reactiveValues()
      
      observe({
        model_sbml_svg <- list(
          list(src = "untitled.svg")
        )
        output$model_sbml <- renderImage({
          list(src = model_sbml_svg[[1]]$src)
        }, deleteFile = FALSE)
      })
      
      output$model_sbml <- renderImage(reactives$graph, deleteFile = FALSE)
      
      # Add a download button for the SVG file
      output$downSBMLfig <- downloadHandler(
        filename = function() {
          "untitled.svg"
        },
        content = function(file) {
          file.copy("untitled.svg", file)
        }
      )
      
      })
      })
        
    }
  })
  
  
  
  observeEvent(input$fluxsbml, {
    file <- input$fluxsbml
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    
    # Validate that the file extension is either "xml" or "sbml"
    validate(need(ext %in% c("xml", "sbml"), "Please upload an XML or SBML file"))
    
    # Determine the desired file extension
    desired_ext <- if (ext == "xml") {
      "xml"
    } else if (ext == "sbml") {
      "sbml"
    }
    
    # Set the file name to "model.xml" or "model.sbml" in the "out/" directory
    save_path <- file.path("out", paste0("model.", desired_ext))
    
    # Copy the uploaded file to the desired path
    file.copy(from = file$datapath, to = save_path, overwrite = TRUE)
  })
  
  
  # Define the reactive values to track the rows
  rv <- reactiveValues(
    scenario = data.frame(
      Reaction = character(),         # Empty character vector for Reaction column
      LowerBound = numeric(),     # Empty numeric vector for Lower Bound column
      UpperBound = numeric(),     # Empty numeric vector for Upper Bound column
      stringsAsFactors = FALSE       # Avoid converting strings to factors
    )
  )
  
  # Load the existing data from the CSV file if it exists
  if (!file.exists("out/scenario.csv")) {
    initial_data <- data.frame(
      Reaction = character(),
      LowerBound = numeric(),
      UpperBound = numeric(),
      stringsAsFactors = FALSE
    )
    write.csv(initial_data, "out/scenario.csv", row.names = FALSE)
  }
  
 
  # Load the existing data from the CSV file if it exists
  if (file.exists("out/scenario.csv")) {
    rv$scenario <- read.csv("out/scenario.csv", stringsAsFactors = FALSE)
  }
  
  output$editableTableScenario <- DT::renderDataTable(server = FALSE, {
    datatable(
      rv$scenario,
      editable = TRUE,
      style = "bootstrap4",
      extensions = c('Buttons', 'Responsive', 'Scroller'),
      options = list(
        dom = 'Bfrtip',
        scrollY = "400px",  # Change the scrollY value to your desired height in pixels
        paging = FALSE,  # Disable paging
        ordering = TRUE,
        buttons = list(
          list(extend = "collection", buttons = c('excel', 'csv'),
               text = '<span class="glyphicon glyphicon-download-alt"></span>'),
          list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_scenario')
        )
      )
    )
  })
  
  
  # Add event for the "Add Row" button
  observeEvent(input$addRowBtn, {
    new_row <- data.frame(Reaction = "", LowerBound = 0, UpperBound = 0)
    rv$scenario <- rbind(rv$scenario, new_row)
  })
  
  # Add event for the "Remove Row" button
  observeEvent(input$removeRowBtn, {
    if (nrow(rv$scenario) > 0) {
      rv$scenario <- rv$scenario[-nrow(rv$scenario), ]
    }
  })
  
  # Automatically save the edited table
  observeEvent(input$editableTableScenario_cell_edit, {
    edited_data <- input$editableTableScenario_cell_edit
    row_index <- edited_data$row
    column_name <- edited_data$col
    new_value <- edited_data$value
    
    # Update the dataframe with the edited value
    rv$scenario[row_index, column_name] <- new_value
    
    # Save the edited dataframe
    write.csv(rv$scenario, "out/scenario.csv", row.names = FALSE)
  })
  
  
  observeEvent(input$runflux,{
    
    if((is.null(input$fluxsbml) && input$typeForFlux< 0 && input$selectSBML =='Upload own SBML/XML')){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please check the xml file you created in CellDesigner and select the flux analysis you are interested in."
          ))
      )}else {
        
        withProgress(message = "Running please wait ...",{
        
        
        # Run Python script in a separate process
        
        command <- "/opt/conda/envs/cnapy-1.1.8/bin/python"
        script <- "Metabolic_Flux_RS.py"
        
        system2(command, args = script)
        
        if(is.null(input$sbmlCsv) && input$selectSBML =='Upload own SBML/XML' && !(file.exists("CNApy_outputs/nodes.csv"))){
          ask_confirmation(
            "form1",
            title = " ",
            type = "info",
            btn_labels = c("Cancel", "OK"),
            allowEscapeKey = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            text =  
              div(HTML("Oops!", "Something went wrong reading the SBML model. Most likely the SBML model is not valid.Please creat it in CellDesigner"
              ))
          )}else{
        
        # Read nodes.csv and edges.csv data
        nodes <- fread("CNApy_outputs/nodes.csv")
        edges <- fread("CNApy_outputs/edges.csv")
       
        
        
        output$fluxnetwork <- renderVisNetwork({
          visNetwork(nodes, edges) %>%
            visEdges(arrows = "from") %>% 
            visOptions(manipulation = TRUE) # Enable manipulation of nodes
        })
        
        output$fluxOutput <- DT::renderDataTable(server=FALSE,{
          req(input$typeForFlux)
          
          file_flux <- switch(input$typeForFlux,
                              "FBA" = "CNApy_outputs/fba_results.csv",
                              "pFBA" = "CNApy_outputs/pfba_results.csv",
                              "FVA" = "CNApy_outputs/fva_results.csv",
                              "EFM" = "CNApy_outputs/efm_flux_results.csv")
          
          
          # Read the selected result CSV file
          flux_results <- read.csv(file_flux)
          
          DT::datatable(flux_results,
                        extensions = c('Buttons','Responsive','Scroller'), 
                        style = "bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering = TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel','csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_flux_result'))
                        ),rownames = FALSE,escape = FALSE)
        })
        
        
        output$fluxOptimize <- DT::renderDataTable(server=FALSE,{
          req(input$typeForFlux)
          
          flux_optimize <- switch(input$typeForFlux,
                                  "FBA" = "CNApy_outputs/optimize_fba.csv",
                                  "pFBA" = "CNApy_outputs/optimize_pfba.csv",
                                  "FVA" = "CNApy_outputs/optimize_fva.csv",
                                  "EFM" = "CNApy_outputs/optimize_efm.csv")
          
          
          # Read the selected result CSV file
          optimize_results <- read.csv(flux_optimize)
          DT::datatable(optimize_results,
                        extensions = c('Buttons','Responsive','Scroller'), 
                        style="bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering=TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel','csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_optimize_result'))
                        ),rownames = FALSE, escape = FALSE)
        })
        
        output$fluxSummary <- DT::renderDataTable(server = FALSE, {
          req(input$typeForFlux)
          
          flux_summary <- switch(input$typeForFlux,
                                 "FBA" = "CNApy_outputs/summary_fba.txt",
                                 "pFBA" = "CNApy_outputs/summary_pfba.txt",
                                 "FVA" = "CNApy_outputs/summary_fva.txt",
                                 "EFM" = "CNApy_outputs/summary_efm.txt")
          
          # Read the selected result text file
          summary_results <- readLines(flux_summary)
          
          # Convert text data to data frame
          summary_results <- data.frame(Content = summary_results)
          
          DT::datatable(summary_results,
                        extensions = c('Buttons','Responsive','Scroller'), 
                        style = "bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering = TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel','csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_summary_result'))
                        ),rownames = FALSE,escape = FALSE)
        })
          }
        
        })
        
      } 
    
    
  })
  ##############################################################################
  ############### Survival Analysis #####################################
  observeEvent(input$runSurv,{
    if(input$typeForSurv<0){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select cancer type."
          ))
      )}else {
        
        
        clin2 <- readRDS(file= paste0("TCGA_CLIN_OUTPUT/",input$typeForSurv,".rds"))
        clin2$vital_status <- ifelse(clin2$vital_status == "Alive", 0, 1)
        
        # Use the reactive expression to create the plot
        output$SurvivalAn <- renderPlotly({
          custom_palette <- brewer.pal(9, "Set1") # You can adjust the number of colors as needed
          
          #surv <- survfit(as.formula(paste("Surv(days_to_death, status) ~", ",input$inVar3,")), data = clin2)
          fit <- eval(parse(text = paste0("survfit(Surv(days_to_death, vital_status) ~ ", input$inVar3, ", data = clin2)")))
          # Create the survival plot using ggplot2 and survminer
          p1 <-ggsurvplot(fit, data = clin2, pval = TRUE, palette = custom_palette)
          plotly::ggplotly(p1[[1]])
        })
        output$survTable <- DT::renderDataTable(server=FALSE,{
          
          DT::datatable(clin2,
                        extensions = c('Buttons', 'Responsive', 'Scroller'),
                        style = "bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering = TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel', 'csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_clinical'))
                        ))
        })
      }
    
  })
  
  ######________________________________________________________________________
  
  expression <- reactiveValues(data = NULL)
  
  observe({
    if (input$runSurv > 0 && input$typeForSurv > 0 && input$IdGene == "gene_name") {
      expr_df <- read.csv(file = paste0("TCGA_LIMMA_OUTPUT/", input$typeForSurv, ".csv"))
      
      outVar4 <- reactive({
        vars4 <- as.data.frame(expr_df$hgnc_symbol)
        vars4 <- vars4[-1,]
        return(vars4)
      })
      
      output$inVarGene <- renderUI({
        selectInput(inputId = "inVar4", label = "Search Gene", choices = outVar4())
      })
      
    } else if (input$typeForSurv > 0 && input$IdGene == "gene_id") {
      expr_df <- read.csv(file = paste0("TCGA_LIMMA_OUTPUT/", input$typeForSurv, ".csv"))
      
      outVar4 <- reactive({
        vars4 <- as.data.frame(expr_df$gene_id)
        vars4 <- vars4[-1,]
        return(vars4)
      })
      
      output$inVarGene <- renderUI({
        selectInput(inputId = "inVar4", label = "Search Gene", choices = outVar4())
      })
    }
  })
  
  observeEvent(input$runSurv2, {
    print(input$inVar4)
    
    output$geneSurv <- renderPlotly({
      print(input$typeForSurv)
      
      load(file = paste0("TCGA_DATA/", input$typeForSurv, ".rda"))
      tcga_data <- data
      expr_df <- read.csv(file = paste0("TCGA_LIMMA_OUTPUT/", input$typeForSurv, ".csv"))
      
      # Extract clinical data
      clinical <- tcga_data@colData
      
      # We are only interested in the "Primary solid Tumor" cases for survival
      primary_tumor_mask <- clinical@listData[["definition"]] == "Primary solid Tumor"
      clin_df <- clinical[primary_tumor_mask, c("patient", "vital_status", "days_to_death", "days_to_last_follow_up", "gender", "ajcc_pathologic_stage")]
      clin_df$deceased <- clin_df$vital_status == "Dead"
      
      # Create an "overall survival" variable that is equal to days_to_death
      # for dead patients, and to days_to_last_follow_up for patients who
      # are still alive
      clin_df$overall_survival <- ifelse(clin_df$deceased,
                                         clin_df$days_to_death,
                                         clin_df$days_to_last_follow_up)
      
      # 5.1 Kaplan-Meier plots
      Surv(clin_df$overall_survival, clin_df$deceased)
      Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender
      fit <- survfit(Surv(overall_survival, deceased) ~ gender, data = clin_df)
      ggsurvplot(fit, data = clin_df, pval = T, risk.table = T, risk.table.col = "strata")
      
      # Use grep() to search for the selected gene (either gene_name or gene_id) in the respective column
      if (input$IdGene == "gene_name") {
        matches <- grep(input$inVar4, expr_df$hgnc_symbol)
      } else if (input$IdGene == "gene_id") {
        matches <- grep(input$inVar4, expr_df$gene_id)
      } else {
        # Handle the case when no option is selected or an invalid option is provided
        # You may choose to display a message or handle it differently based on your requirement.
        # For example, you can set matches to an empty vector if nothing is selected:
        matches <- integer(0)
      }
      
      # Check if there are matches before proceeding
      if (length(matches) > 0) {
        # Use subset() to get the rows where the selected gene was found
        if (input$IdGene == "gene_name") {
          subset_df <- subset(expr_df, grepl(input$inVar4, hgnc_symbol))
        } else if (input$IdGene == "gene_id") {
          subset_df <- subset(expr_df, grepl(input$inVar4, gene_id))
        }
        
        # Continue with the rest of the code
        # Get the gene_id and gene_name from the subset_df
        gene_id <- subset_df$gene_id
        gene_name <- subset_df$hgnc_symbol
        
        # Get the expression values for the selected gene
        col_names <- colnames(tcga_data)
        row_names <- rownames(tcga_data)
        d_mat <- tcga_data@assays@data@listData[["unstranded"]]
        rownames(d_mat) <- row_names
        colnames(d_mat) <- col_names
        
        # Extract the part before the "." in the row names of data
        data_row_names <- sub("\\..*", "", rownames(d_mat))
        rownames(d_mat) <- data_row_names
        d_mat <- as.matrix(t(d_mat))
        
        # Get the expression values for the selected gene
        clin_df$gene_value <- d_mat[rownames(clin_df), gene_id]
        
        # Find the median value of the gene and print it
        median_value <- median(clin_df$gene_value)
        
        # Divide patients into two groups, up and down-regulated.
        # If the patient's expression is greater or equal to the median, we put it
        # among the "up-regulated", otherwise among the "down-regulated"
        clin_df$gene <- ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
        
        # We can fit a survival model, like we did in the previous section
        fit <- survfit(Surv(overall_survival, deceased) ~ gene, data = clin_df)
        
        # We can extract the survival p-value and print it
        pval <- surv_pvalue(fit, data = clin_df)$pval
        custom_palette <- brewer.pal(12, "Dark2")
        p1 <- ggsurvplot(fit, data = clin_df, pval = T, palette = custom_palette, title = paste(input$inVar4))
        plotly::ggplotly(p1[[1]])
      }
    })
  })
  ###############################################################################
  ############### Immune Signature Analysis #####################################
  
  
  observeEvent(input$runImSig, {
    if(input$typeForImsig < 0){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select cancer type!"
          ))
      )}
    
    if (input$typeForImsig > 0) {
      observe({
        DEGs <- read.csv(file = paste0("TCGA_ANOVA_OUTPUT/", input$typeForImsig, ".csv"))
        colnames(DEGs)[colnames(DEGs) == "X"] <- "mRNA"
        colnames(DEGs)[colnames(DEGs) == "gene_name"] <- "geneSymbol"
        unique_gene_symbols <- unique(DEGs$geneSymbol)
        DEGs <- subset(DEGs, !duplicated(geneSymbol))
        
        load(file = paste0("TCGA_DATA/", input$typeForImsig, ".rda"))
        
        gene_expr <- assay(data)
        
        sample_info <- data@colData@listData[["definition"]]
        sample_info <- as.data.frame(sample_info)
        sample_name <- data@colData@rownames
        sample_name <- as.data.frame(sample_name)
        merge_sampleinfo <- cbind(sample_info, sample_name)
        df_tumor <- merge_sampleinfo[grep("Tumor", merge_sampleinfo$sample_info),]
        
        mRNA <- as.data.frame(rownames(gene_expr))
        mRNA <- gsub("\\..*", "", mRNA$`rownames(gene_expr)`)
        mRNA <- as.data.frame(mRNA)
        colnames(mRNA) <- "mRNA"
        gene_expr <- cbind(gene_expr, mRNA)
        
        merged_data <- merge(DEGs, gene_expr, by = "mRNA", all.x = TRUE)
        mydata_unique <- unique(merged_data$geneSymbol)
        mydata_unique <- as.data.frame(mydata_unique)
        mydata_unique <- na.omit(mydata_unique)
        colnames(mydata_unique) <- "geneSymbol"
        
        merged_data <- merged_data[, -2]
        merged_data <- merge(mydata_unique, merged_data, by = "geneSymbol")
        mydata_filtered <- merged_data %>%
          filter(geneSymbol != "")
        
        # Find duplicated geneSymbol
        duplicated_genes <- duplicated(mydata_filtered$geneSymbol)
        
        # Add suffix to duplicated geneSymbol to make them unique
        duplicated_genes_suffix <- ifelse(duplicated_genes, paste0(mydata_filtered$geneSymbol, ".1"), mydata_filtered$geneSymbol)
        
        # Set the row names of mydata_filtered to geneSymbol column with suffix added
        rownames(mydata_filtered) <- duplicated_genes_suffix
        
        mydata_filtered <- mydata_filtered[, -c(1, 2, 3, 4, 5, 6, 7, 8, 9)]
        exp <- mydata_filtered
        
        # Assuming df_tumor$sample_name contains the column names you want to select from mydata_filtered
        selected_columns <- df_tumor$sample_name
        
        # Check if the selected column names are present in mydata_filtered
        matching_columns <- intersect(selected_columns, colnames(mydata_filtered))
        
        # Check for missing columns
        missing_columns <- setdiff(selected_columns, matching_columns)
        
        if (length(missing_columns) > 0) {
          print(paste("The following columns are not found in mydata_filtered:", paste(missing_columns, collapse = ", ")))
        } 
        # If all selected columns are present, proceed with the selection
        mydata_filtered_sub <- mydata_filtered[, matching_columns]
        
        ### update imsig function
        pp_imsig <- function(exp, sig, r = 0.7, sort = TRUE, sort_by = 'T cells') {
          sig_fg <- sig[which(as.character(sig$gene) %in% rownames(exp)),]
          avg_data <- exp[row.names(exp) %in% as.character(sig_fg$gene),]
          cc <- data.frame(matrix(nrow = ncol(exp)))
          cc <- cc[, -1]
          for (i in levels(sig_fg$cell)) {
            s <- sig_fg[sig_fg$cell %in% i,]
            e <- avg_data[as.character(s$gene),]
            e_avg <- data.frame(colMeans(e, na.rm = TRUE))
            colnames(e_avg) <- i
            cc <- cbind(cc, e_avg)
          }
          if (sort == TRUE) {
            cc <- cc[sort.list(unlist(data.frame(cc[, sort_by]))),]
          }
          return(cc)
        }
        
        imsig_output <- pp_imsig(exp = mydata_filtered_sub, sig = sig, r = 0.7)
        imsig_output_bar <- pp_imsig(exp = mydata_filtered, sig = sig, r = 0.7)
        
        # merge with sample info data frame by row names
        data_for_barplot <- merge(imsig_output, df_tumor, by.x = "row.names", by.y = "sample_name", all.x = TRUE)
        merge_data <- merge(imsig_output_bar, merge_sampleinfo, by.x = "row.names", by.y = "sample_name", all.x = TRUE)
        rownames(merge_data) <- merge_data$Row.names
        # remove the row.names column
        merge_data <- merge_data[, -1]
        colnames(merge_data) <- c("B_cells", "Interferon", "Macrophages", "Monocytes", "Neutrophils",
                                  "NK_cells", "Plasma_cells", "Proliferation", "T_cells", "Translation", "Sample_Type")
        
        # Set sample_info as factor variable for ordering the bars in plot
        merge_data$Sample_Type <- factor(merge_data$Sample_Type, levels = unique(merge_data$Sample_Type))
        merge_data <- merge_data[order(merge_data$Sample_Type),]
        
        # selectGene <- feature_select(mydata_filtered_sub)
        selectGene <- sig[which(sig$gene %in% rownames(mydata_filtered)),]
        
        # sample_colors <- ifelse(merge_sampleinfo$Sample_Type == "Primary solid Tumor", "#E74C3C", "#1ABC9C")
        
        output$imsigplot1 <- renderPlotly({
          plot_ly(
            data = merge_data,
            y = ~(get(input$ImsigType)),
            type = 'bar',
            hovertemplate = ~paste(
              "<b>Sample type:</b>", Sample_Type,
              "<br><b>",input$ImsigType,":</b>", input$ImsigType,
              "<br><b>Sample name:</b>", rownames(merge_data)
            ),
            color = ~Sample_Type
          ) %>%
            layout(
              barmode = 'stack',
              yaxis = list(title = input$ImsigType)
            )
        })
        imsig <-reactiveValues(result=NULL)
        download_plots <- reactiveValues()
        output$imsigplot2 <- renderVisNetwork({
          
          use_it <- merge(selectGene, mydata_filtered_sub, by.x = "gene", by.y = "row.names", all.x = TRUE)
          rownames(use_it) <- use_it$gene
          use_it <- use_it[, -c(1, 2)]
          
          # g <- Reduce(intersect, list(as.character(row.names(exp), sig$gene)))
          g <- intersect(row.names(use_it), as.character(sig$gene))
          exp <- use_it[as.character(g),]
          print(str(exp))
          print(summary(exp))
          
          # g <- intersect(row.names(use_it), as.character(sig$gene))
          sig <- sig[sig$gene %in% as.character(g),]
          
          cor_data <- cor(t(exp), method = "pearson")
          cor_data[cor_data <= 0.7] = 0
          diag(cor_data) = 0
          # fg <- feature_select(use_it, 0.7)
          # cor_data <- cor_data[fg,fg]
          network <- graph_from_adjacency_matrix(cor_data, weighted = T, mode = "undirected", diag = F)
          
          edges_d3 <- as.data.frame(get.edgelist(network))
          nodes <- unique(c(edges_d3$V1, edges_d3$V2))
          
          edges_d3 <- data.frame(from = edges_d3$V1,
                                 to = edges_d3$V2)
          
          # edges_d3 <- data.frame(from = edges_d3$V1,
          #                      to = edges_d3$V2,
          #                     width = edge.betweenness(network,e=E(network)))
          nodes_info <- sig[sig$gene %in% as.character(nodes),]
          nodes_d3 <- data.frame(id = nodes_info$gene,
                                 label = nodes_info$gene,
                                 group = nodes_info$cell)
          nodes_d4 <- data.frame(nodes_d3)
          imsig$result <- nodes_d4
         
          imsig_network <- reactive({
            visNetwork(nodes = nodes_d3, edges = edges_d3) %>%
              # Create legend with title 'Legend'
              visLegend(main = "Legend") %>%
              visPhysics(enabled = FALSE) %>%
              visIgraphLayout(layout = "layout_nicely")
            
          })
          output$downloadNetwork2 <- downloadHandler(
            filename = function() {
              paste('Imsig_network-', Sys.Date(), '.html', sep='')
            },
            content = function(con) {
              imsig_network() %>% visSave(con, selfcontained = TRUE)
            })
          
          imsig_network() 
          
        })
        
        
        
        
        output$info_imsig <- DT::renderDataTable(server = FALSE, {
          
          DT::datatable(imsig_output_bar,
                        extensions = c('Buttons', 'Responsive', 'Scroller'),
                        style = "bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering = TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel', 'csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_imSig'))
                        ))
          
        })
        
        output$info_imsig2 <- DT::renderDataTable(server = FALSE, {
          
          DT::datatable(sig,
                        extensions = c('Buttons', 'Responsive', 'Scroller'),
                        style = "bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering = TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel', 'csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_sig'))
                        ))
          
        })
        
        output$info_imsig_out <- DT::renderDataTable(server = FALSE, {
          imsigRes <- imsig$result
          
          DT::datatable(imsigRes,
                        extensions = c('Buttons', 'Responsive', 'Scroller'),
                        style = "bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering = TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel', 'csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_sig'))
                        ))
          
        })
        
        
        
      })
    }
  })
  
  
  ##############################################################################
  ################## Single Nucleotide Variations Analysis ##################
  mutation <- reactiveValues(data = NULL)
  gene_data <- reactiveVal(NULL)
  
  observe({
    
    if (input$typeForSNV > 0) {
      maf <- read.csv(file = paste0("TCGA_MUT_OUTPUT/", input$typeForSNV, ".csv"))
      mutation$data <- as.data.frame(maf)
      my_table <- table(maf$Hugo_Symbol)
      my_table_sorted <- sort(my_table, decreasing = TRUE)
      my_table_sorted <- as.data.frame(my_table_sorted)
      #my_table_sorted <- head(my_table_sorted,input$GeneSize)  # Assuming you want the top 10 genes
      
      
      outVar <- reactive({
        vars <- my_table_sorted$Var1  # Assuming the gene names are in the first column
        return(vars)
      })
      
      output$inVar5 <- renderUI({
        selectInput(inputId = "inVar5", label = "Gene", choices = outVar())
      })
    }
  })
  
  observeEvent(input$inVar5, {
    # Update the selected_gene reactiveVal based on the selected gene
    maf <- mutation$data
    maf <- maf %>% read.maf
    selected_gene <- subset(maf@data, Hugo_Symbol == input$inVar5)
    selected_gene <- subset(selected_gene, select = -c(1:2))
    df <- selected_gene[complete.cases(selected_gene$Protein_position), ]
    df$Protein_position <- sub("/.*", "", df$Protein_position)
    gene_data(df)  # Using parentheses to access the value inside the reactiveVal
  })
  
  
  observeEvent(input$runSNV,{
    if(input$typeForSNV < 0){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select cancer and plot type!"
          ))
      )}else{
    #maf <- read.csv(file= paste0("TCGA_MUT_OUTPUT/",input$typeForSNV, ".csv"))
    maf <-  mutation$data
    maf <- maf %>% read.maf
    selected_gene <- subset(maf@data, Hugo_Symbol == input$inVar5)
    selected_gene <- subset(selected_gene, select = -c(1:2))
    df <- selected_gene[complete.cases(selected_gene$Protein_position), ]
    df$Protein_position <- sub("/.*", "", df$Protein_position)
    
    output$lollipop <-renderG3Lollipop({
      req(gene_data())  # Ensure that the selected_gene data is available
      df <- gene_data()
      
      rs.options <-g3Lollipop.options(
        chart.width = 1500,
        chart.type = "circle",
        # lollipop track options
        lollipop.track.height = 300,
        lollipop.track.background = "transparent",
        lollipop.pop.max.size = 4,
        lollipop.pop.min.size = 4,
        # set larger than lollipop.pop.max.size to turn off pop info
        lollipop.pop.info.limit = 4.1,
        # y-axis label
        lollipop.line.color = "grey",
        lollipop.line.width = 0.5,
        lollipop.circle.color = "grey",
        lollipop.circle.width = 0.5,
        lollipop.color.scheme = "bottlerocket2",
        y.axis.line.color = "transparent",
        # domain annotation bar
        anno.bar.fill = "#969696",
        anno.bar.margin = list(top = 4, bottom = 8),
        # domain track options
        domain.margin = list(top = 2, bottom = 6),
        domain.color.scheme = "bottlerocket1",
        domain.text.font = "normal 12px Arial",
        domain.text.color = "white",
        # highlight text
        highlight.text.angle = 45,
        # disable brush
        brush = FALSE,
        # disable legend
        legend = TRUE
      )
      

      g3Lollipop(df,
                 gene.symbol = input$inVar5,
                 btn.style = "gray", # gray-style chart download buttons
                 plot.options = rs.options,
                 factor.col = "Variant_Classification",
                 aa.pos.col = "Protein_position",
                 output.filename = "Cat-E_SNV")
      
    })
    
    output$Vaf <-renderPlotly({
      plotVaf(maf = maf)
    })
    output$Somatic <-renderPlotly({
      somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))
    })
    
    
    download_plots <- reactiveValues()
    # Define a reactive value to store the output_df
    output_df_reactive <- reactiveVal()
    
   
    output$oncoplot1 <- renderPlot({
      # Store the plot in a variable
     oncoplot(maf = maf, top = input$GeneSize, removeNonMutated = TRUE, draw_titv = TRUE, fontSize = 1.0)
      
    })
    
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste('myplots.png', sep='')
      },
      content = function(file) {
        
          png(file, width = 1200, height = 500, res = NA)
          oncoplot(maf = maf, top = input$GeneSize, removeNonMutated = TRUE, draw_titv = TRUE, fontSize = 1.0)
          dev.off()
       
      },
      contentType = 'image/png'
    )
    
    
   
   
    
    # Update the SNVtable when the selected gene changes
    output$SNVtable <- DT::renderDataTable(server = FALSE, {
      req(gene_data())  # Ensure that the selected_gene data is available
      df <- gene_data()
      DT::datatable(df,
                    extensions = c('Buttons','Responsive','Scroller'), 
                    style="bootstrap4",
                    options = list(
                      dom = 'Bfrtip',
                      scrollY = 400,
                      ordering=TRUE,
                      buttons = list(
                        list(extend = "collection", buttons = c('excel','csv'),
                             text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                        list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_SNV'))
                    ))
    })
    
    
    output$SNVheatmap <- renderPlotly({
      maf_df <- mutation$data
      select <- head(maf@gene.summary,input$GeneSize)
      top_genes <- as.data.frame(select$Hugo_Symbol)
      
      
      # Filter the maf data for the top genes
      maf_filtered <- maf_df[maf_df$Hugo_Symbol %in% top_genes$`select$Hugo_Symbol`, ]
      
      
      heatmap_data <- maf_filtered %>%
        group_by(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
        summarise(Count = n()) %>%
        arrange(match(Hugo_Symbol, select$Hugo_Symbol))
      
      
      heatmap_plot <-ggplot(heatmap_data, aes(x = Tumor_Sample_Barcode, y = Hugo_Symbol, fill = Variant_Classification)) +
        geom_tile() +
        # scale_fill_viridis(discrete = TRUE) +  # Using the viridis color palette
        labs(
          x = "Tumor Sample Barcode",
          y = "Hugo Symbol")+
        theme(legend.position="bottom",
              axis.title.x = element_blank(),
              axis.text.x = element_blank()   # Remove x-axis tick labels
        ) 
      
      
      heatmap_plotly <- ggplotly(heatmap_plot)
      
      # Create a barplot (histogram)
      
      barplot_plot <- ggplot(heatmap_data, aes(x = Tumor_Sample_Barcode, y = Count)) +  ylab(NULL) +
        geom_bar(stat = "identity") +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),  # Remove panel background
          plot.background = element_rect(fill = "transparent", color = NA)   # Remove plot background
        )
      
      
      # Combine both plots into a single interactive plot
      combined_plot <- subplot(
        barplot_plot,
        heatmap_plot,
        nrows = 2)
      
      # Make the plot interactive using plotly
      combined_plot <- ggplotly(combined_plot)
      
      # Display the interactive plot
      combined_plot
      
    })
    
    
    
    output$SNVheatTable <- DT::renderDataTable(server = FALSE, {
      select <- head(maf@gene.summary,input$GeneSize)
      
      DT::datatable(select,
                    extensions = c('Buttons','Responsive','Scroller'), 
                    style="bootstrap4",
                    options = list(
                      dom = 'Bfrtip',
                      scrollY = 400,
                      ordering=TRUE,
                      buttons = list(
                        list(extend = "collection", buttons = c('excel','csv'),
                             text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                        list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_snv'))
                    ))
    })
      }
    
  })
  ###########################################################################
  ###########################################################################
  
  ########### Enrichment Analysis ################################
  reactives <- reactiveValues(userData = NULL)
  
  observeEvent(input$userEnrich, {
    req(input$userEnrich)  # Ensure a file is uploaded
    df <- read.csv(input$userEnrich$datapath)
    reactives$userData <- df
    
    # Update the pickerInput choices based on the uploaded data columns
    updatePickerInput(session, inputId = 'geneColumnE', label = 'Select Gene Column', choices = colnames(df))
    updatePickerInput(session, inputId = 'logFColumnE', label = 'Select log2FC Column', choices = colnames(df))
    updatePickerInput(session, inputId = 'pAdColumnE', label = 'Select pAdj Column', choices = colnames(df))
  })
  
  user_data <- reactive({
    req(reactives$userData)
    df <- reactives$userData
    df
  })
  
  output$UserData <- renderDataTable({
    user_data()
  })
  observeEvent(input$runEnrich, {
    if(input$typeForEnrc<0 && input$userinputEnrich ==FALSE){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select cancer type data GO analysis."
          ))
      )}else {
        
        withProgress(message = 'I need about  1 minute to finish complete  Enrichment analysis...', {
          
          if(input$typeForEnrc > 0 && input$userinputEnrich ==FALSE){
            data_file <- reactive({
              file_name <- paste0(input$typeForEnrc)
              file_path <- paste0("TCGA_ANOVA_OUTPUT/", file_name)
              return(file_path)
              
            })
            DEGsData <- read.csv(file= paste0(data_file(),".csv"))
            
            filtered_df <- filter(DEGsData, logFC > input$logFCcut | logFC < -input$logFCcut)
            cleanDeg <- filtered_df %>% filter(FDR < input$qvalcut)
            
            
          }
          
          if(input$userinputEnrich > 0 && input$typeForEnrc < 0 ){
            
            user_df <- as.data.frame(reactives$userData)
            logFC <- (user_df[input$logFColumnE])
            GeneCol <- (user_df[input$geneColumnE])
            pAdjCol <- (user_df[input$pAdColumnE])
            DEGsData <- data.frame(GeneCol, logFC, pAdjCol)
            colnames(DEGsData) <- c("gene_name", "LogFC", "PValue")
            data$DEGsData <- DEGsData
            filtered_df <- filter(DEGsData, logFC > input$logFCcut | logFC < -input$logFCcut)
            cleanDeg <- filtered_df %>% filter(PValue < input$qvalcut)
           
            
          }
          
          print(input$typeForEnrc)
          
          
          #system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",GeneList))
          # Assuming genes is the dataframe with columns "gene_name" and "PValue"
          unique_genes <- cleanDeg %>%
            distinct(gene_name, .keep_all = TRUE)
          
          # Extract the gene names as a character vector
          gene_names <- unique_genes$gene_name
          
          # Extract the P-values as a numeric vector
          p_values <- unique_genes$PValue
          
          # Create a named list with gene names as names and P-values as values
          gene_p_values <- setNames(p_values, gene_names)
          
          # Define a vector with the ontology types (BP, MF, CC)
          ontology_types <- c("BP", "MF", "CC")
          
          # Create an empty list to store the results for each ontology type
          all_goEnrichment <- list()
          
          selection <- function(allScore) {
            return(allScore < input$qvalcut)
          }
          # Loop through each ontology type
          selected_enrichment <- reactiveVal(NULL)
          for (ontology_type in ontology_types) {
            # Create the GOdata object for the current ontology type
            allGO2genes <- annFUN.org(whichOnto = ontology_type, feasibleGenes = NULL, mapping = "org.Hs.eg.db", ID = "symbol")
            GOdata <- new("topGOdata",
                          ontology = ontology_type,
                          allGenes = gene_p_values,
                          annot = annFUN.GO2genes,
                          GO2genes = allGO2genes,
                          geneSel = selection,
                          nodeSize = 10)
            
            # Run the test and obtain the enrichment results
            results.ks <- runTest(GOdata, algorithm = "classic", statistic = "ks")
            goEnrichment <- GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 20)
            goEnrichment$KS <- as.numeric(goEnrichment$KS)
            
            # Filter for significant results (KS < 0.05)
            goEnrichment <- goEnrichment[goEnrichment$KS < input$qvalcut, ]
            
            # Store the results for the current ontology type in the list
            all_goEnrichment[[ontology_type]] <- goEnrichment[, c("GO.ID", "Term", "KS")]
            
          }
          
          # Assuming you have already defined all_goEnrichment with separate data frames for each ontology type
          observe({
            if (input$EnrichType == "BPplot1") {
              #selected_enrichment(NULL)
              selected_enrichment(all_goEnrichment$BP)
            } 
            if (input$EnrichType == "MFplot3") {
              #selected_enrichment(NULL)
              selected_enrichment(all_goEnrichment$MF)
            }
            if (input$EnrichType == "CCplot2") {
              #selected_enrichment(NULL)
              selected_enrichment(all_goEnrichment$CC)
              
            }
          })
          
          
          # Create the reactive plot
          output$whichenrich <- renderPlotly({
            # Get the selected enrichment data for the plot
            goEnrichment <- selected_enrichment()
            
            if (!is.null(goEnrichment)) {
              # Calculate the negative log10 of KS values
              goEnrichment$KS_log10 <- -log10(goEnrichment$KS)
              
              # Define colors based on significance level (KS value)
              colors <- ifelse(goEnrichment$KS_log10 < 5, '#FF5733', ifelse(goEnrichment$KS_log10 < 10, '#FFC300', '#008080'))
              
              # Concatenate GO.ID and Term columns for bar labels
              goEnrichment$Label <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ":")
              
              # Create the plot using plot_ly
              pl <- plot_ly(goEnrichment, x = ~KS_log10, y = ~Label, type = 'bar', marker = list(color = colors)) %>%
                layout(
                  xaxis = list(
                    title = "Enrichment (-log10(p-value))",
                    tickvals = round(seq(0, max(goEnrichment$KS_log10), by = 2), 1)
                  ),
                  yaxis = list(
                    title = "GO Term",
                    tickfont = list(size = 12),
                    tickangle = 0
                  ),
                  title = paste("GO Enrichment Results"),
                  showlegend = FALSE
                )
              
              pl
            }
          })
          
          output$TablEnrich <- DT::renderDataTable(server=FALSE,{
            
            if (input$EnrichType > 0) {
              goEnrichment <- selected_enrichment()
              DT::datatable(goEnrichment,
                            extensions = c('Buttons','Responsive','Scroller'), 
                            style="bootstrap4",
                            options = list(
                              dom = 'Bfrtip',
                              scrollY = 400,
                              ordering=TRUE,
                              buttons = list(
                                list(extend = "collection", buttons = c('excel','csv'),
                                     text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                                list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_enrich'))
                            ))
            }
          })
          
        })
      }
  })
  
  
  #### enrich KEGG
  v <- reactiveValues(
    data = NULL,
    edgevis = data.frame(),
    for_symbol = data.frame(),
    current_nodev1 = NULL,
    node_list = data.frame() # Add a new reactiveValues object for the list of nodes
  )
  
  userSigDF <- reactiveValues(SiGeneListFC= data.frame(),
                              SiGeneList = data.frame())
  
  observeEvent(input$runVisPPI2, {
    if(input$typeforKEGG<0 && input$KEGGinput != "Cancer Data" && input$virusID2 <0){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Before the click Network button, please select cancer type data and select virus."
          ))
      )}else{
    
    
    observe({
      edges <-e$data
      
      colnames(edges) <- c("from","Gene_1","Desc_1","tax_from", "organism_1","link","to","Gene_2","Desc_2", "tax_to","organism_2")
      
      # Find the indices of empty elements in Column1
      empty_indices <- which(is.na(edges$Gene_1) | edges$Gene_1 == "")
      #print(empty_indices)
      # Fill in the empty elements in Column1 with data from Column2
      edges$Gene_1[empty_indices] <- edges$from[empty_indices]
      
      # Find the indices of empty elements in Column1
      empty_indices <- which(is.na(edges$Gene_2) | edges$Gene_2 == "")
      #print(empty_indices)
      # Fill in the empty elements in Column1 with data from Column2
      edges$Gene_2[empty_indices] <- edges$to[empty_indices]
      
      nodevistb <- edges
      
      colnames(nodevistb) <- c("From_A","Gene_A","Description_A","Taxonomy_A","Organism_A", "Interaction_Type","To_B","Gene_B","Description_B", "Taxonomy_B","Organism_B")
      nodevistb <- nodevistb[,c(1,7,4,2,8,5,10,6,11,3,9)]
      v$edgevis <- as.data.frame(nodevistb)
      
      Edges <- as.data.frame(nodevistb$Gene_A)
      Edges$to <- nodevistb$Gene_B
      Edges$link <- nodevistb$Interaction_Type
      colnames(Edges) <- c("from","to","link")
      
      nodes1 <- as.data.frame(nodevistb$From_A)
      nodes1$tax <- nodevistb$Taxonomy_A
      nodes1$group <- nodevistb$Organism_A
      nodes1$genes <- nodevistb$Gene_A
      colnames(nodes1) <- c("id","tax","group","genes")
      
      
      nodes2 <- as.data.frame(nodevistb$To_B)
      nodes2$tax <- nodevistb$Taxonomy_B
      nodes2$group <- nodevistb$Organism_B
      nodes2$genes <- nodevistb$Gene_B
      colnames(nodes2) <- c("id","tax","group","genes")
     
      
      nodev1 <- rbind(nodes1,nodes2)
      nodev1 <- distinct(nodev1)
      
      v$current_nodev1 <- nodev1
      
      nodev2 <- as.data.frame(nodev1)
      
      nodev2$id <- nodev1$genes
      nodev2$label <- nodev1$genes
      nodev2 <- nodev2%>% mutate(font.size = ifelse(id=="GL068-",34,
                                                    ifelse(id=="GL121-",34,
                                                           ifelse(id=="GL245-",34,
                                                                  ifelse(id=="gusA",27,
                                                                         ifelse(id=="lacZ",27,
                                                                                ifelse(id=="RUC-GFP",27,
                                                                                       16)))))))
      
      for_symbol <- NULL
      for_symbol$id <- as.data.frame(nodev1$id)
      for_symbol$gene <- as.data.frame(nodev1$genes)
      for_symbol <- as.data.frame(for_symbol)
      colnames(for_symbol) <- c("id","symbol")
      v$for_symbol <- for_symbol
      
      nodev2$title <- paste0("Gene :",nodev1$genes,
                             "<br>UniprotKB :",nodev1$id,
                             "<br>Species :",nodev1$group,
                             "<br>Taxonomy :",nodev1$tax)
      
      
      
      nodes <- as.data.frame(nodev2)
      
      #v$data<- as.data.frame(nodes)
      
      nodes <- subset(nodes, select = -c(genes))
      nodes <- nodes %>% group_by(id) %>% filter (! duplicated(id)) 
      
      edges <- as.data.frame(Edges)
      
      edges<-edges %>% mutate(color = ifelse(link=="activation","#5ebf87",
                                             ifelse(link=="inhibition","#e0392b",
                                                    ifelse(link=="MI:0194","#d0d676",
                                                           ifelse(link=="MI:0203","#d6b076",
                                                                  ifelse(link=="MI:0217","#8ec4b7",
                                                                         ifelse(link=="MI:0220","#8eb5c4",
                                                                                ifelse(link=="MI:0403","#998ec4",
                                                                                       ifelse(link=="MI:0407","#bd8ec4",
                                                                                              ifelse(link=="MI:0414","#c48ea3",
                                                                                                     ifelse(link=="MI:0570","#d3815f",
                                                                                                            ifelse(link=="MI:0844","#bfd35f",
                                                                                                                   ifelse(link=="MI:0882","#76d35f",
                                                                                                                          ifelse(link=="MI:0914","#a994bf",
                                                                                                                                 ifelse(link=="MI:0915","#797e85",
                                                                                                                                        ifelse(link=="MI:0945","#5fcfd3",
                                                                                                                                               ifelse(link=="MI:1126","#5f94d3",
                                                                                                                                                      ifelse(link=="MI:2364","#d35fcc",
                                                                                                                                                             "#9c95a2"))))))))))))))))))
      
      edges<-edges %>% mutate(arrows.to.type = ifelse(link=="activation","arrow",
                                                      ifelse(link=="inhibition","bar", "image")))
      edges<-edges %>% mutate(title = edges$link)
      edges <- distinct(edges)
      edges$id <- row.names(edges)
      
      
      ledges <-edges[, c("color", "title","arrows.to.type")]
      ledges <- distinct(ledges)
      
      colnames(ledges) <- c("color", "label","arrows")
      
      gene_lst <- reactiveVal(nodes)
      graph_data$nodes <- nodes
      
      graph_data$edges <- edges
      graph_data$ledges <-ledges
      v$data <- nodes
      
    })
    
    
    graph_data = reactiveValues(
      nodes = NULL,
      edges =NULL,
      ledges=NULL
    )
    
    
    output$editable_network <- renderVisNetwork({
      
      PPI_network <- reactive({
        visNetwork(graph_data$nodes, graph_data$edges,height = "800px", width = "120%") %>%
          visIgraphLayout(layout = "layout_nicely")%>%
          visPhysics(solver = "repulsion")%>%
          #visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -10))%>%
          visGroups(groupname = "Homo sapiens", color =list(highlight ="#C0392B",background = "#7986CB",border ="#7986CB"))%>% 
          visEdges(arrows = "to",label = graph_data$edges$link, width =4, shadow = TRUE, smooth = list(enabled = TRUE, roundness= 0.1)) %>%
          visInteraction(multiselect = TRUE, navigationButtons = TRUE, hover = TRUE)%>%
          visExport(type = "pdf")%>%
          visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_selection', nodes.nodes);
                ;}")%>%
          
          visLegend(addEdges = graph_data$ledges, main = list(text = "Legend",
                                                              style = "font-family:Comic Sans MS;color: #212f3c;font-size:20px;text-align:center;"),
                    
                    width = 0.20, stepY= 53, position = "right")%>%
          visOptions(manipulation = T)
      })
      
      output$downloadNetwork <- downloadHandler(
        filename = function() {
          paste('network-', Sys.Date(), '.html', sep='')
        },
        content = function(con) {
          PPI_network() %>% visSave(con, selfcontained = TRUE)
        })
      
      PPI_network() 
    })
    
    
    observeEvent(input$editable_network_graphChange, {
      
      # If the user added a node, add it to the data frame of nodes.
      if(input$editable_network_graphChange$cmd == "addNode") {
        temp = bind_rows(
          graph_data$nodes,
          data.frame(id = input$editable_network_graphChange$id,
                     label = input$editable_network_graphChange$label,
                     stringsAsFactors = F)
        )
        graph_data$nodes = temp
        
      }
      # If the user added an edge, add it to the data frame of edges.
      else if(input$editable_network_graphChange$cmd == "addEdge") {
        temp = bind_rows(
          graph_data$edges,
          data.frame(id = input$editable_network_graphChange$id,
                     from = input$editable_network_graphChange$from,
                     to = input$editable_network_graphChange$to,
                     stringsAsFactors = F)
        )
        graph_data$edges = temp
      }
      # If the user edited a node, update that record.
      else if(input$editable_network_graphChange$cmd == "editNode") {
        temp = graph_data$nodes
        temp$label[temp$id == input$editable_network_graphChange$id] = input$editable_network_graphChange$label
        graph_data$nodes = temp
      }
      # If the user edited an edge, update that record.
      else if(input$editable_network_graphChange$cmd == "editEdge") {
        temp = graph_data$edges
        temp$from[temp$id == input$editable_network_graphChange$id] = input$editable_network_graphChange$from
        temp$to[temp$id == input$editable_network_graphChange$id] = input$editable_network_graphChange$to
        graph_data$edges = temp
      }
      # If the user deleted something, remove those records.
      else if(input$editable_network_graphChange$cmd == "deleteElements") {
        for(node.id in input$editable_network_graphChange$nodes) {
          temp = graph_data$nodes
          temp = temp[temp$id != node.id,]
          graph_data$nodes = temp
        }
        for(edge.id in input$editable_network_graphChange$edges) {
          temp = graph_data$edges
          temp = temp[temp$id != edge.id,]
          graph_data$edges = temp
        }
      }
    })
    
    #############
    
    vis_cho <- reactive({
      vis_node_list <- as.data.frame(graph_data$nodes)
      vis_node <- as.data.frame(vis_node_list$id)
      #colnames(vis_node) <- "Node List"
      print(vis_node)
      v$updateDat <- vis_node
      return(vis_node)
    })
    
    
    
    observe({
      updateSelectizeInput(session, 'SearchVisNode', choices = vis_cho())
    })
    
    observeEvent(input$btn1,{
      if(input$SearchVisNode > 0){
        isolate({
          nodes_visSearch <- as.data.frame(v$data)
          current_vis_node <- nodes_visSearch[grep(input$SearchVisNode, nodes_visSearch$label), "id"]
          visNetworkProxy("editable_network") %>% visSelectNodes(id  = current_vis_node)
        })
      }
    })
    
    #, smooth = list(enabled = TRUE, roundness= 0.1)
    # display go info by selected nodes
    observeEvent(input$current_node_selection, {
      edgeslist_for_current <-  v$current_nodev1
      result <- edgeslist_for_current[edgeslist_for_current$genes == input$current_node_selection, "id"]
      print(result)
      
      
      
      # get position info
      observeEvent(input$store_position, {
        visNetworkProxy("editable_network") %>% visGetPositions()
      })
      
      # format positions
      nodes_positions <- reactive({
        positions <- input$network_positions
        if(!is.null(positions)){
          nodes_positions <- do.call("rbind", lapply(positions, function(x){ data.frame(x = x$x, y = x$y)}))
          nodes_positions$id <- names(positions)
          nodes_positions
        } else {
          NULL
        }
      })
      
    })
      }
  })
  
  
  
  
  observeEvent(input$runKEGG, {
    if(input$typeforKEGG<0 && input$KEGGinput != "PPI Data"){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select cancer type data for KEGG analysis."
          ))
      )
      return()
      }
    if(input$typeforKEGG<0 && input$KEGGinput == "Cancer Data"){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select cancer type data for KEGG analysis."
          ))
      )
      return()
      }
    if(input$typeforKEGG<0 && input$KEGGinput != "Cancer Data"){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select cancer type data and select virus."
          ))
        )
      return()
      }else {
     
      if (input$KEGGinput == "PPI Data") {
        
        clean_data <- v$for_symbol
        genes_symbolup <- v$updateDat
        colnames(genes_symbolup) <- "symbol"
        
        ##### warning part
        spsComps::shinyCatch({
          if (suppressWarnings(
            nrow(ids <- bitr(genes_symbolup$symbol, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db"))<1
          )
          )
            
            stop("Please select correct organisms")
          
          message("Continues...")
          
        }, blocking_level = "error")
        
       
        DEG_data <- read.csv(file = paste0("TCGA_LIMMA_OUTPUT/", input$typeforKEGG, ".csv"))
        
        colnames(DEG_data)[colnames(DEG_data) == "hgnc_symbol"] <- "SYMBOL"
        
        selected_DEG_data <- merge(DEG_data, ids, by = "SYMBOL")
        
        
      }
      if (input$KEGGinput == "Cancer Data") {
        selected_DEG_data <- read.csv(file = paste0("TCGA_LIMMA_OUTPUT/", input$typeforKEGG, ".csv"))
        
        colnames(selected_DEG_data)[colnames(selected_DEG_data) == "hgnc_symbol"] <- "SYMBOL"
        
        spsComps::shinyCatch({
          if (suppressWarnings(
            nrow(ids <- bitr(selected_DEG_data$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID","UNIPROT"), OrgDb="org.Hs.eg.db"))<1
          )
          )
            
            stop("Please select correct organisms")
          
          message("Continues...")
          
        }, blocking_level = "error")
        selected_DEG_data <- merge(selected_DEG_data, ids, by = "SYMBOL")
      }
      
      
      geneSig <- selected_DEG_data$logFC
      names(geneSig) <- selected_DEG_data$ENTREZID
      geneSig <- na.omit(geneSig)
      userSigDF$SiGeneListFC <- geneSig
      userSigDF$SiGeneList <- names(geneSig)
      
      observe({
        
        id_cnvrt <- bitr_kegg(userSigDF$SiGeneList, fromType = 'ncbi-geneid', toType = 'kegg', organism = "hsa")
        
        userSigDF$SiGeneList <- id_cnvrt$kegg
        
        ekk <- enrichKEGG(
          gene = userSigDF$SiGeneList,
          organism = "hsa",
          minGSSize = 1,
          maxGSSize = 500,
          pAdjustMethod = "fdr",  # Adjusted for the example
          pvalueCutoff = input$pvalCut,
          qvalueCutoff = input$qvalCut
        )
        
        kk_df <- data.frame(ekk)
        row.names(kk_df) <- NULL
        
        output$kk_df <- DT::renderDataTable({
          DT::datatable(kk_df,
                        extensions = c('Buttons','Responsive','Scroller'), 
                        style="bootstrap4",
                        options = list(
                          dom = 'Bfrtip',
                          scrollY = 400,
                          ordering=TRUE,
                          buttons = list(
                            list(extend = "collection", buttons = c('excel','csv'),
                                 text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                            list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_kegg'))
                        ), escape = FALSE)
        })
         
      })
      
      keggFig <- reactiveValues()
      pathviewReactive <- eventReactive(input$PathviEw,{
        withProgress(message = 'Plotting Pathview Pathway...', {
          
          isolate({
            setProgress(value = 0.7, detail = paste0("Pathview ID ",input$PathviEw," ..."))
            dme <- pathview(gene.data  =userSigDF$SiGeneListFC,
                            pathway.id = input$PathviEw,
                            species    = "hsa",
                            limit      = list(gene=max(abs(userSigDF$SiGeneListFC)), cpd=1))
            file.copy(paste0(input$PathviEw,".pathview.png"),paste0(input$PathviEw))
            
            #setProgress(value = 0.7, detail = paste0("Pathview ID ",input$PathviEw," generating pdf ..."))
            dmePdf <- pathview(gene.data  =userSigDF$SiGeneListFC,
                               pathway.id = input$PathviEw,
                               species    = "hsa", kegg.native = F,
                               limit      = list(gene=max(abs(userSigDF$SiGeneListFC)), cpd=1))
            
            keggFig$imagePath = paste0(input$PathviEw,".pathview.")
            return(list(
              src = paste0(input$PathviEw),
              filetype = tempfile("image/png"),
              alt = "pathview image"
            ))
          })
          
        })
      })
      
      output$pathway1  = renderImage({
        
        return(pathviewReactive())
        
      })
      
      
      
      output$downloadPathviewPng <- downloadHandler(
        filename = function()  {paste0(keggFig$imagePath,"png")},
        content = function(file) {
          file.copy(paste0(getwd(),'/',keggFig$imagePath,"png"), file)
        }
      )
      
      output$downloadPathviewPdf <- downloadHandler(
        filename = function()  {paste0(keggFig$imagePath,"pdf")},
        content = function(file) {
          file.copy(paste0(getwd(),'/',keggFig$imagePath,"pdf"), file)
        }
      )
        }
    
    
  })
  
  
  
  ######### celline info tab ################################################
  
  output$info_celline <- DT::renderDataTable(server=FALSE,{
    
    if (input$searchCell > 0) {
      info_cell <- read.csv("data/depmap_sample_info.csv")
      info_cell <- data.table(info_cell)
      
      searched_celline <- info_cell[cell_line_name %ilike% input$searchCell]
      
      DT::datatable(searched_celline,
                    extensions = c('Buttons','Responsive','Scroller'), 
                    style="bootstrap4",
                    options = list(
                      dom = 'Bfrtip',
                      scrollY = 400,
                      ordering=TRUE,
                      buttons = list(
                        list(extend = "collection", buttons = c('excel','csv'),
                             text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                        list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_celline'))
                    ))
    }
  })
  
  ######### drug-tissue and celline ################################################
  
  output$drug <- DT::renderDataTable(server=FALSE,{
    
    if (input$tissue > 0) {
      drugs <- read.csv("data/drug-celline_.csv")
      
      filtered_tissue <- drugs %>% filter(TCGA.classification == input$tissue)
      
      DT::datatable(filtered_tissue,
                    extensions = c('Buttons','Responsive','Scroller'), 
                    style="bootstrap4",
                    options = list(
                      dom = 'Bfrtip',
                      scrollY = 400,
                      ordering=TRUE,
                      buttons = list(
                        list(extend = "collection", buttons = c('excel','csv'),
                             text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                        list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_drug_celline'))
                    ))
    }
  })
  
  
  ############################### cell type form panglaoDB ################
  
  output$out_celltype <- DT::renderDataTable(server=FALSE,{
    
    if (input$tissueType > 0) {
      
      getSampleComposition(
        tissue = input$tissueType)
    }
  })
  
  output$out_cellMarker <- DT::renderDataTable(server=FALSE,{
    
    if (input$cellType > 0) {
      cell_type <- read_excel("data/Cell_marker_All.xlsx")
      filtered_cell_type <- cell_type %>% filter(cell_type == input$cellType)
      DT::datatable(filtered_cell_type,
                    style="bootstrap4",
                    extensions = c('Buttons','Responsive','Scroller'),
                    options = list(
                      dom = 'Bfrtip',
                      scrollY = 400,
                      ordering=TRUE,
                      buttons = list(
                        list(extend = "collection", buttons = c('excel','csv'),
                             text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                        list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_celltype'))
                    ))
      
      
    }
  })
  
  ############# JIMENA TAB ##############################################
  #____ graphml converter __#
  
  #### uploaded txt file in box
  
  output$whichjim <-  DT::renderDataTable(server=FALSE,{
    if (input$choseJimena == "userJim") {
      shinyjs::hide("targeTable")
      
      withProgress(message = "Importing your file...",{
        observe({
          
          inFile <- input$file1
          print(inFile)
          if (is.null(inFile)) return(NULL)   
          dff <- read.csv(inFile$datapath,header = TRUE)
          print(dff)
          write.csv(dff, "input.txt",row.names = FALSE,quote=FALSE)
          
        }) 
      })
      
      user_txt <- read.csv("input.txt",header = TRUE,quote=FALSE)
      
      DT::datatable(user_txt,
                    extensions = c('Responsive','Scroller'), 
                    options = list(
                      lengthChange = FALSE,
                      deferRender = TRUE,
                      ordering = TRUE,
                      scrollY = 200,
                      scroller = TRUE))
    }
    else if(input$choseJimena == "exampleJim") {
      
      shinyjs::hide("targeTable")
      test_file <- read.table("testJimena.txt",header = TRUE)
      
      write.table(test_file, "input.txt",row.names = FALSE,quote=FALSE)
      DT::datatable(test_file,
                    extensions = c('Responsive','Scroller'), 
                    options = list(
                      lengthChange = FALSE,
                      deferRender = TRUE,
                      ordering = TRUE,
                      scrollY = 200,
                      scroller = TRUE
                    ))}
    source_data <- reactiveVal(data.frame(node1 = character(0), label = character(0), node2 = character(0)))
    target_data <- reactiveVal(data.frame(node1 = character(0), label = character(0), node2 = character(0)))
    
    # Function to update source table based on search input
    updateSourceTable <- function(gene) {
      # Replace this with your code to fetch data based on the gene symbol
      # For this example, we generate some sample data
      req(e$data)
      edges <- e$data
      
      colnames(edges) <- c("node1", "label", "node2")
      source_data(edges)
    }
    
    # Initialize source table
    observe({
      if (input$choseJimena == "createJim" && input$searchProtein > 0) {
        updateSourceTable(input$searchProtein)
        shinyjs::show("targeTable")
      }
    })
    
    # Render source table
    output$source_table <- renderDT({
      datatable(source_data(), selection = 'single')
    })
    
    # Render target table
    output$target_table <- renderDT({
      datatable(target_data(), selection = 'multiple')
    })
    
    # JavaScript code to capture selected rows and add them to the target table
    observeEvent(input$source_table_rows_selected, {
      selected_rows <- input$source_table_rows_selected
      if (length(selected_rows) > 0) {
        row_data <- source_data()[selected_rows, , drop = FALSE]
        current_target_data <- target_data()
        new_target_data <- rbind(current_target_data, row_data)
        target_data(new_target_data)
      }
    })
    
    # JavaScript code to remove selected rows from the target table
    observeEvent(input$remove_button, {
      selected_rows <- input$target_table_rows_selected
      if (length(selected_rows) > 0) {
        current_target_data <- target_data()
        new_target_data <- current_target_data[-selected_rows, , drop = FALSE]
        target_data(new_target_data)
      }
    })
    
    # Write target_data to a file (You can modify this part as needed)
    observe({
      if (input$choseJimena == "createJim" && input$searchProtein > 0) {
        write.table(target_data(), "input.txt", row.names = FALSE, quote = FALSE)
      }
    })
    
  })
  
  
  #### run button for jimena 
  observeEvent(input$runGraphml,{
    withProgress(message = "Converting graphml...",{
      system(paste("perl txt_to_graphml.pl"))
    })
  })
  ##################################################
  daf <- reactiveVal()
  add_per <- reactiveValues(data=NULL)
  result_val <- reactiveVal()
  add_per <- reactiveValues(NodeInfo=NULL)
  add_per <- reactiveValues(NodeIndex=NULL)
  add_per <- reactiveValues(remained_nodes=NULL)
  remain_per<- reactiveValues(NodeIndex=NULL)
  remain_per<- reactiveValues(data_name=NULL)
  add_parameter <- reactiveValues()
  remain_parameter <- reactiveValues()
  squad_nodes <- reactiveValues(data=NULL)
  
  jimenaout <- reactiveValues(data = NULL)
  observeEvent(input$runJimena,{
    if(input$runGraphml == FALSE) {
      ask_confirmation(
        "form3",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please convert your txt file to graphml format before running Jimena."
          ))
      )
      
    }
    
    withProgress(message = 'Please wait for Jimena calculation...', {
      req(input$runGraphml)
      GRAPHML_PATH <- "input.graphml"
      
      system("rm jimena_time_series_data.csv")
      system(paste("java -classpath jimena-app/jimena.jar:jimena-app/  App ", GRAPHML_PATH, sep=" "))
      
      spsComps::shinyCatch({
        if (nrow(suppressWarnings(jimena_output <- read.csv("jimena_time_series_data.csv")))<1){
          stop("Upload correct data!")
        }else{
          message("Running...")
          
        }
        
      }, blocking_level = "error")
      
      jimena_output_long <- melt(jimena_output, id.vars = "time")   
      
      jimenaout$data <-jimena_output
      
      accumulate_by <- function(dat, var) {
        var <- lazyeval::f_eval(var, dat)
        lvls <- plotly:::getLevels(var)
        dats <- lapply(seq_along(lvls), function(x) {
          cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
        })
        dplyr::bind_rows(dats)
      }
      
      fig <- jimena_output_long %>% accumulate_by(~time)
      jimena_time_series <- 
        
        fig%>%plot_ly()%>%
        add_lines(x =~time,  y = ~value, split = ~variable)
      jimena_time_series <- jimena_time_series %>% layout(
        xaxis = list(
          title = "Time",
          zeroline = F
        ),
        yaxis = list(
          title = "Activation",
          zeroline = F
        )
      )
      
      result_val(jimena_time_series)
      
      add_per$data <- data.frame(jimena_output)
      Snode <- as.data.frame(colnames(add_per$data))
      squad_nodes$data <- Snode[-1,]
      
      print(remain_per$data_name)
      
      outVar <- reactive({
        vars <- as.data.frame(colnames(add_per$data))
        vars <- vars[-1,]
        print(vars)
        
        return(vars)
        
      })
      
      output$inVar2 <- renderUI({
        selectInput(inputId = "inVar2", label = h4("Add Perturbation"), choices =  outVar())
        
        
      })
      
    })
    
    
  }) 
  
  output$jimenaResult <-renderPlotly(
    result_val())
  
  
  #___ JIMENA ADD PERTURBATION PART ___#
  
  
  observeEvent(input$add_pert, {
    req(input$runJimena)
    t = rbind(data.frame(Nodes = input$inVar2,
                         Start = input$strt,
                         End = input$end,
                         Value = input$val), daf())
    t <- distinct(t)
    daf(t)
    all_node_ <- as.data.frame(colnames(add_per$data))
    all_node_ <- all_node_[-1,]
    all_node_ <- as.data.frame(all_node_)
    node_length <- nrow(all_node_)
    all_node_$index <- rep(0:(node_length-1))
    colnames(all_node_) <- c("Nodes", "NodeIndex")
    
    
    
    add_per$NodeInfo <-all_node_
    r <- all_node_[which(all_node_$Nodes==input$inVar2),]
    remain_parameter<- r
    remain_per$data_name <- remain_parameter$Nodes
    
    print(r$NodeIndex)
    add_per$NodeIndex <- r$NodeIndex
    add_parameter <- rbind(data.frame(    Index = 0,
                                          Nodes = add_per$NodeIndex,
                                          Start = input$strt,
                                          End = input$end,
                                          Value = input$val))
    
    add_parameter <- distinct(add_parameter)
    csv_fname = "jimena-app/ParameterInputs.csv"
    write.table(add_parameter, file = csv_fname, sep = ",",
                append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    
  })
  
  
  
  observeEvent(input$remove_per, {
    req(input$runJimena)
    t = daf()
    print(t)
    print(nrow(t))
    
    if (!is.null(input$shiny_table_rows_selected)) {
      t <- t[-as.numeric(input$shiny_table_rows_selected),]
      add_per$remained_nodes <- t
      
    }
    daf(t)
    
    remain <- as.data.frame(add_per$remained_nodes)
    
    node_info <- add_per$NodeInfo
    
    
    last_remain_info <- subset(node_info, (Nodes %in% remain$Nodes))
    
    if (nrow(last_remain_info)>0){
      
      print(last_remain_info)
      merged_part <- merge(last_remain_info, remain, by= "Nodes")
      
      
      remain_parameter <- rbind(data.frame(    Index = 0,
                                               Nodes = merged_part$NodeIndex,
                                               Start = merged_part$Start,
                                               End = merged_part$End,
                                               Value = merged_part$Value,
                                               Fix = NA))
      
      remain_parameter <- distinct(remain_parameter)
      colnames(remain_parameter) <- c(0,0.05,10,1,0.1,0.5)
      
      print(remain_parameter)
      csv_fname = "jimena-app/ParameterInputs.csv"
      write.table(remain_parameter, file = csv_fname, sep = ",",
                  append = FALSE, quote = FALSE,
                  col.names = TRUE, row.names = FALSE)
    }
    
    if (nrow(last_remain_info)==0){
      
      remain_parameter <- rbind(data.frame(    A = 0,
                                               B = 0.05,
                                               C = 10,
                                               D = 1,
                                               E = 0.1,
                                               G = 0.5))
      
      print(remain_parameter)
      csv_fname = "jimena-app/ParameterInputs.csv"
      write.table(remain_parameter, file = csv_fname, sep = ",",
                  append = FALSE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }
    
  })
  
  output$shiny_table <- renderDataTable(
    daf(), selection = 'single', extensions = c('Buttons','Responsive','Scroller'), 
    options = list(
      dom = 'Bfrtip',
      lengthChange = FALSE,
      deferRender = TRUE,
      ordering = TRUE,
      scrollY = 100,
      scroller = TRUE,
      buttons = 
        list( list(
          extend = 'collection',
          
          buttons = list(list(extend='csv',
                              filename = 'CatE_perturbation_table'),
                         list(extend='excel',
                              filename = 'CatE_perturbation_table'),
                         list(extend='pdf',
                              filename= 'CatE_perturbation_table')),
          text = '<span class="glyphicon glyphicon-download-alt"></span>'
        )),
      
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#0080F1', 'color': '#fff'});",
        "}")))
  
  
  
  #____ Jimena output Data Table ____#
  
  output$jimenaoutput_tb <- DT::renderDataTable(server = FALSE,{
    DT::datatable(jimenaout$data,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  options = list(
                    dom = 'Bfrtip',
                    lengthChange = FALSE,
                    deferRender = TRUE,
                    ordering = TRUE,
                    scrollY = 400,
                    scroller = TRUE,
                    buttons = 
                      list( list(
                        extend = 'collection',
                        
                        buttons = list(list(extend='csv',
                                            filename = 'CatE_jimenaout_table'),
                                       list(extend='excel',
                                            filename = 'CatE_jimenaout_table'),
                                       list(extend='pdf',
                                            filename= 'CatE_jimenaout_table')),
                        text = '<span class="glyphicon glyphicon-download-alt"></span>'
                      )),
                    
                    initComplete = JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#0080F1', 'color': '#fff'});",
                      "}")))})
  
  
  
  #########################################################################
  ####### AlphaFold Protein Structure Prediction ##########################
  
  
  
  # Define color bands for pLDDT visualization
  PLDDT_BANDS <- list(
    list(start = 0, end = 50, color = '#FF7D45', label = "Very low (pLDDT < 50)"),
    list(start = 50, end = 70, color = '#FFDB13', label = "Low (70 > pLDDT > 50)"),
    list(start = 70, end = 90, color = '#65CBF3', label = "Confident (90 > pLDDT > 70)"),
    list(start = 90, end = 100, color = '#0053D6', label = "Very high (pLDDT > 90)")
  )
  
  observeEvent(input$download_btn, {
    pdb_name <- input$pdb_name
    
    # Define the URL of the PDB file
    pdb_url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", pdb_name, "-F1-model_v4.pdb")
    
    # Send a GET request to download the PDB file
    response <- GET(pdb_url)
    
    # Check if the request was successful (status code 200)
    if (status_code(response) == 200) {
      # Save the content to a file
      file_path <- paste0("www/",pdb_name, ".pdb")
      file_content <- content(response, "text")
      writeLines(file_content, file_path)
      print("PDB file downloaded successfully.")
      
      # Parse the PDB file to extract pLDDT scores and atom information
      pdb_data <- readLines(file_path)
      plddt_scores <- numeric(length(pdb_data))
      atom_info <- character(length(pdb_data))
      
      for (i in seq_along(pdb_data)) {
        line <- pdb_data[i]
        if (substr(line, 1, 4) == "ATOM" || substr(line, 1, 6) == "HETATM") {
          plddt_scores[i] <- as.numeric(substr(line, 61, 66))
          atom_info[i] <- paste0(substr(line, 13, 16), substr(line, 18, 20), ":", substr(line, 22, 26), ":", plddt_scores[i])
        }
      }
      
      # Apply a cartoon style with color based on pLDDT score
      output$structure <- renderNGLVieweR({
        NGLVieweR(file_path) %>%
          addRepresentation(input$type_structure,
                            param = list(
                              name = "cartoon",
                              colorScheme = "bfactor",
                              colorValue = plddt_scores,
                              colorOptions = list(palette = sapply(PLDDT_BANDS, function(b) b$color))
                            )
          ) %>%
          stageParameters(backgroundColor = input$backgroundColor) %>%
          setQuality("high") %>%
          setFocus(0) %>%
          setSpin(FALSE) %>%
          addRepresentation("ball+stick", param = list(
            name = "ball+stick", sele =  isolate(input$selection),
            colorScheme = "element", colorValue = input$col
          ))
        
        #%>%
        # zoomMove(
        #  center = "27:B",
        # zoom = "27:B",
        #z_offSet = -20
        #)
      })
    }
  })
  
  output$info <- renderText({
    input$structure_selection
  })
  
  output$info2 <- renderText({
    input$structure_sequence
  })
  
  observe({
    
    
    if(input$type_surface  != "hide") {
      NGLVieweR_proxy("structure") %>%
        addSelection(
          "surface",
          param = list(
            name = "surface",
            colorValue = isolate(input$type_surface),
            opacity = 0.2
            
          )
        )
    }
  })
  
  observeEvent(input$screenshot, {
    NGLVieweR_proxy("structure") %>%
      snapShot("Snapshot", param = list(
        antialias = TRUE,
        trim = TRUE,
        transparent = TRUE,
        scale = 1
      ))
  })
  
  observeEvent(input$fullscreen, {
    NGLVieweR_proxy("structure") %>%
      updateFullscreen()
  })
  
  
  observeEvent(input$remove, {
    NGLVieweR_proxy("structure") %>%
      removeSelection("sel1")
  })
  
  output$legend <- renderUI({
    tags$div(
      style = "text-align: center;",
      lapply(PLDDT_BANDS, function(b) {
        tags$div(
          style = paste0(
            "display: flex; align-items: center; margin-bottom: 5px;"
          ),
          tags$div(
            style = paste0(
              "width: 20px; height: 20px; background-color: ", b$color, "; margin-right: 5px;"
            )
          ),
          b$label
        )
      })
    )
  })
  
  # Highlight relevant part of protein structure on hover
  observeEvent(input$structure_selection, {
    NGLVieweR_proxy("structure") %>%
      removeSelection("hover") %>%
      addSelection(
        "line",
        param = list(
          name = "hover",
          sele = input$structure_selection,
          color = "red",
          opacity = 0.8
        )
      )
  })
  
  observe({
    
    if(input$animat == "spin")
    {
      NGLVieweR_proxy("structure") %>%
        updateSpin(TRUE)
    }
    if(input$animat == "rock"){
      NGLVieweR_proxy("structure") %>%
        updateRock(TRUE)
    }else if(input$animat == "none")
    {
      NGLVieweR_proxy("structure") %>%
        
        updateRock(FALSE) %>%
        updateSpin(FALSE)
      
    }
  })
  
  
  #########################################################################
  #########################################################################
  ############### DRUG DESCRIPTION PART ###################################
  # Draw an SVG image of the drug molecule structure
  
  observeEvent(input$searchDrugName, {
    if(input$searchDrugName > 0){
      
      
      drug_name <- input$searchDrugName
      print(drug_name)
      # Search for the InChI string and other information of the drug using the PubChem database
      url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/InChI,MolecularFormula/json", drug_name)
      response <- urllib$urlopen(url)
      content <- response$read()
      data <- json$loads(content$decode("utf-8"))
      inchi_string <- data$PropertyTable$Properties[[1]]$InChI
      cid <- data$PropertyTable$Properties[[1]]$CID
      print(cid)
      molecular_formula <- data$PropertyTable$Properties[[1]]$MolecularFormula
      if (is.null(molecular_formula) || molecular_formula == "") {
        molecular_formula <- "Not available"
      }
      
      # Get the synonyms and molecular weight information of the drug using the PubChem database
      synonyms_url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/synonyms/json", drug_name)
      synonyms_response <- urllib$urlopen(synonyms_url)
      synonyms_content <- synonyms_response$read()
      synonyms_data <- json$loads(synonyms_content$decode("utf-8"))
      synonyms <- synonyms_data$InformationList$Information[[1]]$Synonym
      
      molecular_weight_url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/MolecularWeight/json", drug_name)
      molecular_weight_response <- urllib$urlopen(molecular_weight_url)
      molecular_weight_content <- molecular_weight_response$read()
      molecular_weight_data <- json$loads(molecular_weight_content$decode("utf-8"))
      molecular_weight <- molecular_weight_data$PropertyTable$Properties[[1]]$MolecularWeight
      
      # Get the chemical structure of the drug using RDKit
      mol <- rdkit$Chem$MolFromInchi(inchi_string)
      
      # Generate the Canonical SMILES, Standard InChI, and Standard InChI Key information
      canonical_smiles <- rdkit$Chem$MolToSmiles(mol)
      standard_inchi <- rdkit$Chem$MolToInchi(mol)
      standard_inchikey <- rdkit$Chem$InchiToInchiKey(standard_inchi)
      
      observe({
        # Draw an SVG image of the drug molecule structure
        #output$drug_image <- renderImage({
        
        drawer <- rdkit$Chem$Draw$rdMolDraw2D$MolDraw2DSVG(as.integer(500), as.integer(500))
        drawer$DrawMolecule(mol)
        drawer$FinishDrawing()
        svg <- drawer$GetDrawingText()
        
        # Modify the SVG to set a transparent background
        svg <- gsub("<rect style='opacity:1.0;fill:#FFFFFF;stroke:none'", "<rect style='opacity:0.0;fill-opacity:0;stroke-opacity:0'", svg)
        
        # Write the SVG to a file
        svg_file_name <- paste0("www/",drug_name, ".svg")
        
        writeLines(svg, con = svg_file_name)
        
        
        # Return a list containing the filename and the raw SVG data
        list(src = svg_file_name, width = 460, height = 460)
      })
      
      
      
      output$drug_desc <- DT::renderDataTable(server=FALSE,{
        
        drug_desc <- read.csv("data/drugDescriptionHyperlink.csv")
        drug_desc <- data.table(drug_desc)
        #drug_desc$DrugBank_Accession_Number <- sprintf(
        # '<a href="https://go.drugbank.com/drugs/%s" target="_blank">%s</a>',
        #drug_desc$DrugBank_Accession_Number,
        #drug_desc$DrugBank_Accession_Number
        #)
        
        #searched_genes <- drugsgenes %>% filter(gene_name %like% input$searchBttn)
        
        searched_drug <- drug_desc[Generic_Name %ilike% input$searchDrugName]
        DT::datatable(searched_drug,
                      extensions = c('Buttons','Responsive','Scroller'), 
                      style="bootstrap4",
                      options = list(
                        dom = 'Bfrtip',
                        scrollY = 400,
                        ordering=TRUE,
                        buttons = list(
                          list(extend = "collection", buttons = c('excel','csv'),
                               text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                          list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_drug'))
                      ),escape = FALSE)
      })
      
      # Update the drug information in the UI
      output$drug_info <- renderPrint({
        cat("InChI string:", inchi_string, "\n",
            "PubChem CID:", cid, "\n",
            "Canonical SMILES:", canonical_smiles, "\n",
            "Standard InChI:", standard_inchi, "\n",
            "Standard InChI Key:", standard_inchikey, "\n",
            "Molecular Formula:", molecular_formula, "\n",
            "Synonyms:", synonyms, "\n",
            "Molecular Weight:", molecular_weight)
      })
      
      output$my_output <- renderUI({
        lapply(1, function(i) {
          div(
            img(
              src = paste0(input$searchDrugName, ".svg"),
              class = "small-img",
              id = paste0("img", i)
            ),
            tags$script(HTML(paste0('$("#img', i, '").popover({trigger: "hover", html: true, content: "<img src= ', if (!is.null(input$searchDrugName)) input$searchDrugName else '', '.svg class=img-popover>"});')))
          )
        })
      })
      
      
    }
  })
  
  observeEvent(input$searchIndication, {
    if (input$searchIndication > 0) {
      drug_desc <- read.csv("data/drugDescriptionHyperlink.csv")
      drug_desc <- data.table(drug_desc)
      
      searched_indic <- drug_desc[Indication %ilike% input$searchIndication]
      
      output$drug_indic <- DT::renderDataTable({
        DT::datatable(searched_indic,
                      extensions = c('Buttons','Responsive','Scroller'), 
                      style="bootstrap4",
                      options = list(
                        dom = 'Bfrtip',
                        scrollY = 400,
                        ordering=TRUE,
                        buttons = list(
                          list(extend = "collection", buttons = c('excel','csv'),
                               text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                          list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_drug'))
                      ), escape = FALSE)
      })
    }
  })
  observeEvent (input$choseDrug == "searchSmi", {
    
    if(input$searchSmile > 0){
      
      output$drug_sim <- DT::renderDataTable(server=FALSE,{
        
        drug_smi <- read.csv("data/drugDescriptionHyperlink.csv")
        drug_smi <- data.table(drug_smi)
        smiles_input <- input$searchSmile
        # Use str_detect to identify rows that contain target_smiles
        rows_with_target <- which(str_detect(drug_smi$Smiles, fixed(smiles_input)))
        
        # Extract the rows from drug_smi
        drug_smi_with_target <- drug_smi[rows_with_target,]
        
        DT::datatable(drug_smi_with_target,
                      extensions = c('Buttons','Responsive','Scroller'), 
                      style="bootstrap4",
                      options = list(
                        dom = 'Bfrtip',
                        scrollY = 400,
                        ordering=TRUE,
                        buttons = list(
                          list(extend = "collection", buttons = c('excel','csv'),
                               text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                          list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_drug'))
                      ),escape = FALSE)
      })
    }
  })
  
  
  #########################################################################
  ######### drug and genes ################################################
  
  output$searchgene <- DT::renderDataTable(server=FALSE,{
    
    if (input$searchBttn > 0) {
      drugsgenes <- read.csv("data/dgidb_interaction.csv")
      drugsgenes <- data.table(drugsgenes)
      
      #searched_genes <- drugsgenes %>% filter(gene_name %like% input$searchBttn)
      
      searched_genes <- drugsgenes[gene_name %ilike% input$searchBttn]
      DT::datatable(searched_genes,
                    extensions = c('Buttons','Responsive','Scroller'), 
                    style="bootstrap4",
                    options = list(
                      dom = 'Bfrtip',
                      scrollY = 400,
                      ordering=TRUE,
                      buttons = list(
                        list(extend = "collection", buttons = c('excel','csv'),
                             text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                        list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_gene_drug'))
                    ))
    }
  })
  
  
  #################### STITCH WEBSITE ########################################
  
  output$stitch_web <- renderUI({
    tags$iframe(src = "http://stitch.embl.de/cgi/input.pl?UserId=UGnzcX9PnC89&sessionId=9M3znCRZvuNh",
                width = "100%", height = "800px")
  })
  
  ############################################################################
  #############################################################################
  #################### CAR-T cell Therapy Tab #################################
  output$car_Therapy <- DT::renderDataTable(server=FALSE,{
    
    if (input$carDisease > 0) {
      cartdata <- read.csv("data/CAR-TcellTherapy_hyper.csv")
      cartdata <- data.table(cartdata)
      
      #searched_genes <- drugsgenes %>% filter(gene_name %like% input$searchBttn)
      
      searched_therapy <- cartdata[Conditions %ilike% input$carDisease]
      DT::datatable(searched_therapy,
                    extensions = c('Buttons','Responsive','Scroller'), 
                    style="bootstrap4",
                    options = list(
                      dom = 'Bfrtip',
                      scrollY = 400,
                      ordering=TRUE,
                      buttons = list(
                        list(extend = "collection", buttons = c('excel','csv'),
                             text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                        list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_CarT'))
                    ),escape = FALSE)
    }
  })
  
  #############################################################################
  #################### Bispecific Therapy Tab #################################
  output$bispecific_Therapy <- DT::renderDataTable(server=FALSE,{
    
    
    bispecificdata <- read.csv("data/CytostaticTherapy_hyper.csv")
    bispecificdata <- data.table(bispecificdata)
    
    DT::datatable(bispecificdata,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  style="bootstrap4",
                  options = list(
                    dom = 'Bfrtip',
                    scrollY = 400,
                    ordering=TRUE,
                    buttons = list(
                      list(extend = "collection", buttons = c('excel','csv'),
                           text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                      list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_Cytostatic'))
                  ),escape = FALSE)
  })
  
  
  #############################################################################
  #################### Cytostatic Therapy Tab #################################
  output$cytostatic_Therapy <- DT::renderDataTable(server=FALSE,{
    
    
    cytostaticdata <- read.csv("data/CytostaticTherapy_hyper.csv")
    cytostaticdata <- data.table(cytostaticdata)
    
    #searched_genes <- drugsgenes %>% filter(gene_name %like% input$searchBttn)
    
    
    DT::datatable(cytostaticdata,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  style="bootstrap4",
                  options = list(
                    dom = 'Bfrtip',
                    scrollY = 400,
                    ordering=TRUE,
                    buttons = list(
                      list(extend = "collection", buttons = c('excel','csv'),
                           text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                      list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_cyto'))
                  ),escape = FALSE)
  })
  
  #############################################################################
  #################### Oncolytic Virus Therapy Tab #################################
  output$ovirus_Therapy <- DT::renderDataTable(server=FALSE,{
    
    
    ovirusTherapdata <- read.csv("data/oncolyticvirusTherapy_hyper.csv")
    ovirusTherapdata <- data.table(ovirusTherapdata)
    
    DT::datatable(ovirusTherapdata,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  style="bootstrap4",
                  options = list(
                    dom = 'Bfrtip',
                    scrollY = 400,
                    ordering=TRUE,
                    buttons = list(
                      list(extend = "collection", buttons = c('excel','csv'),
                           text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                      list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_onco'))
                  ),escape = FALSE)
  })
  ############################################################################
  #################### CheckPoint Inhibitor Therapy Tab #################################
  output$checkpoint_target <- DT::renderDataTable(server=FALSE,{
    data <- read_excel("data/CKTTDB.xlsx", sheet = "Targets")
    
    data <- data.table(data)
    DT::datatable(data,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  style="bootstrap4",
                  options = list(
                    dom = 'Bfrtip',
                    scrollY = 400,
                    ordering=TRUE,
                    buttons = list(
                      list(extend = "collection", buttons = c('excel','csv'),
                           text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                      list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_chekpoint'))
                  ))
  })
  
  output$checkpoint_modulator <- DT::renderDataTable(server=FALSE,{
    data <- read_excel("data/CKTTDB.xlsx", sheet = "Modulators")
    
    data <- data.table(data)
    DT::datatable(data,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  style="bootstrap4",
                  options = list(
                    dom = 'Bfrtip',
                    scrollY = 400,
                    ordering=TRUE,
                    buttons = list(
                      list(extend = "collection", buttons = c('excel','csv'),
                           text = '<span class="glyphicon glyphicon-download-alt"></span>'),
                      list(extend = 'pdf', pageSize = 'A3', orientation = 'landscape', filename = 'Cat-E_checkpoint'))
                  ))
  })
  
  ################### L O G O #########################################################
  
  #### About us server #########################################################
  
  person_info <- list(
    list(src = "images/person1.png"),
    list(src = "images/person2.png"),
    list(src = "images/person3.png")
    # Add more person information as needed
  )
  
  # Render images and information
  output$person1 <- renderImage({
    list(src = person_info[[1]]$src, alt = "Person 1", width = "250px", style = "margin-bottom: 50px;")
  }, deleteFile = FALSE)
  
  output$person2 <- renderImage({
    list(src = person_info[[2]]$src, alt = "Person 2", width = "250px", style = "margin-bottom: 50px;")
  }, deleteFile = FALSE)
  
  output$person3 <- renderImage({
    list(src = person_info[[3]]$src, alt = "Person 3", width = "250px", style = "margin-bottom: 50px;")
  }, deleteFile = FALSE)
  ##############################################################################
  ###################################
  # app button --------------------------------------------------------------
  output$btnVal <- renderText(input$myAppButton)
  
  remove_files_out <- function() {
    files_out <- list.files("out")
    for (file in files_out) {
      file_path <- file.path("out", file)
      if (!file.info(file_path)$isdir) { # Check if file is not a directory
        file.remove(file_path)
      }
    }
  }
  
  remove_files_cnapy <- function() {
    files_cnapy <- list.files("CNApy_outputs")
    for (file in files_cnapy) {
      file_path_cnapy <- file.path("CNApy_outputs", file)
      if (!file.info(file_path_cnapy)$isdir) { # Check if file is not a directory
        file.remove(file_path_cnapy)
      }
    }
  }
  
  remove_files <- function() {
    files <- list.files("www")
    for (file in files) {
      file_path <- file.path("www", file)
      if (!file.info(file_path)$isdir && file != "logo.png" && file != "cancer1.png"&& file != "cancer2.png"&& file != "cancer3.png"&& file != "viruscancer.png"&& file != "cancer4.png") { # Check if file is not a directory and not "logo.png"
        file.remove(file_path)
      }
    }
  }
  
  
  session$onSessionEnded(function() {
    stopApp()
    system("rm jimena_time_series_data.csv")
    system("rm hsa*")
    system("rm untitled.svg")
    system("rm output.html")
    source("killConnection.R", local = TRUE)
    killDbConnections()
    remove_files()
    remove_files_cnapy()
    remove_files_out()
    clean_parameter <- rbind(data.frame(    A = 0,
                                            B = 0.05,
                                            C = 10,
                                            D = 1,
                                            E = 0.1,
                                            G = 0.5))
    
    print(clean_parameter)
    csv_fname = "jimena-app/ParameterInputs.csv"
    write.table(clean_parameter, file = csv_fname, sep = ",",
                append = FALSE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
  })
  
  
  ############################################################################
  
  
}

shinyApp(ui, server)
