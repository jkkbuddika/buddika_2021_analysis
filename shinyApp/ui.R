# The ui.R file creates the user interface of the Shiny App.
# This Shiny App was written by Kasun Buddika to distribute data from Buddika et al., 2021.

## Load required libraries
library(DT)
library(ggiraph)
library(tidyverse)
library(shinycustomloader)
library(shinythemes)

navbarPage(
  theme = shinytheme("sandstone"),
  "Buddika et al. 2021 Data",
  
  ## Home page
  tabPanel(
    "Home",
    wellPanel(
      "This Shiny App has been written and deployed by Kasun Buddika to distribute 
      data presented in",
      a(href = "https://www.biorxiv.org/content/10.1101/2020.06.27.175398v1.abstract",
        "Buddika et al. 2021", style = "color: darkgreen"), ", a study that is currently in review.",
      "This app provides access to, results of the differential gene expression (DGE) analyses 
        performed in the study, volcano plot visualization of DGE data and data related to the 
        genetic screen described in the manuscript. Read the abstract given below."
    ),
    hr(),
    wellPanel(
      h4(strong("Coordinate transcriptional and post-transcriptional repression of 
                pro-differentiation genes maintains intestinal stem cell identity")),
      "The role of Processing bodies (P-bodies), key sites of post-transcriptional control, 
      in adult   stem cells remains poorly understood. Here, we report that adult Drosophila 
      intestinal stem cells, but not surrounding differentiated cells such as absorptive 
      Enterocytes (ECs), harbor P-bodies that contain Drosophila orthologs of mammalian P-body
      components DDX6, EDC3, EDC4 and LSM14A/B. A targeted RNAi screen in intestinal progenitor 
      cells identified 39 previously known and 64 novel P-body regulators, including Patr-1, 
      a gene necessary for P-body assembly. Loss of Patr-1-dependent P-bodies leads to a loss 
      of stem cells that is associated with inappropriate translation and expression of EC-fate 
      gene nubbin. Transcriptomic analysis of progenitor cells identifies a cadre of such weakly
      transcribed pro-differentiation transcripts that are elevated after P-body loss. Altogether, 
      this study identifies a coordinated P-body dependent, translational and transcriptional 
      repression program that maintains a defined set of in vivo stem cells in a state primed 
      for differentiation.",
      style = "background: #F0F8FF"
    )
  ),
  
  ## RNA-seq data page
  tabPanel(
    "RNA-seq Data",
    fluidPage(
      wellPanel(
        "This page provides differential gene expression datatables associated with",
        a(href = "https://www.biorxiv.org/content/10.1101/2020.06.27.175398v1.abstract",
          "Buddika et al. 2021", style = "color: #DC143C"), "(in review).",
        p("Use the scroll down menu to change the dataset you wish to access. The study present 
          data generated in the current study and reanalyzed date from",
          a(href = "https://www.embopress.org/doi/full/10.15252/embj.201489072",
            "Korzelius et al. 2014", style = "color: #DC143C"), "and",
          a(href = "https://elifesciences.org/articles/48160#content",
            "Khaminets et al. 2020", style = "color: #DC143C"), "available from public databases."),
        p("Select the dataset and click the download button to download the data to your computer"),
        style = "background: #DCDCDC"
      ),
      selectInput(
        "select_dataset",
        label = "Select a Dataset",
        choices = NULL),
      downloadButton(
        "download_datatables",
        "Download Data"
      ),
      hr(),
      DTOutput("dgeData"))),
  
  ## Volcano plot page
  tabPanel(
    "Plots",
    fluidPage(
      wellPanel(
        "This page visualises differential gene expression data associated with",
        a(href = "https://www.biorxiv.org/content/10.1101/2020.06.27.175398v1.abstract",
          "Buddika et al. 2021", style = "color: goldenrod2"), "(in review) in an
        interactive volcano plot.",
        p("Use the scroll down menu to change the dataset you wish to visualize. The study present 
          data generated in the current study and reanalyzed date from",
          a(href = "https://www.embopress.org/doi/full/10.15252/embj.201489072",
            "Korzelius et al. 2014", style = "color: goldenrod2"), "and",
          a(href = "https://elifesciences.org/articles/48160",
            "Khaminets et al. 2020", style = "color: goldenrod2"), "available from public databases."),
        p("Select the dataset and hover over the plot to identify represented genes. Click the 
          download icon to download the plot to your computer."),
        style = "background: #FFEFFF"
      ),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "select_dataset_plot",
            label = "Select a Dataset",
            choices = NULL),
          style = "background: #B2D3C2"),
        
        mainPanel(
          withLoader(girafeOutput("plot")))
        ))),
  
  ## The RNAi screen page
  tabPanel(
    "Screen",
    fluidPage(
      wellPanel(
        "This page provides results of the targeted P-body genetic screen described in",
        a(href = "https://www.biorxiv.org/content/10.1101/2020.06.27.175398v1.abstract",
          "Buddika et al. 2021", style = "color: #1E90FF"), 
        "manuscript. The P-body morphology was assessed using two key markers", 
        strong("TRAL"), "and", strong("PATR-1."), "The effect on IPSGs described in",
        a(href = "https://jcs.biologists.org/content/133/10/jcs243451.abstract",
          "Buddika et al. 2020", style = "color: #1E90FF"),
        "was assessed using canonical stress granule markers",
        strong("FMRP"), "and", strong("ROX8."),
        "Note that not all genes were secreened with IPSG markers.",
        style = "background: #F9D7D7"
      ),
      p("Click the download button below to download the data to your machine. 
        Abbreviations used in the table include:", 
        strong("B = Big"), ",", strong("D = Diffuse"), ",",
        strong("S = Small"), ",", strong("DB = Diffuse + Big"), "and",
        strong("N = Normal.")), br(),
        downloadButton(
          "download_screenData",
          "Download Data"
        ),
        hr(),
        DTOutput("screenData"))),
  collapsible = TRUE
  )



