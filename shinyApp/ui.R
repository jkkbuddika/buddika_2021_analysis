library(DT)
library(ggiraph)
library(tidyverse)
library(shinycustomloader)
library(shinythemes)

navbarPage(
  theme = shinytheme("sandstone"),
  "Buddika et al. 2021 Data",
  
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
    wellPanel(
      h4(strong("Coordinate transcriptional and post-transcriptional repression of pro-differentiation genes maintains intestinal stem cell identity")),
      "The role of Processing bodies (P-bodies), key sites of post-transcriptional control, in adult stem cells remains poorly understood. Here, we report that adult Drosophila intestinal stem cells, but not surrounding differentiated cells such as absorptive Enterocytes (ECs), harbor P-bodies that contain Drosophila orthologs of mammalian P-body components DDX6, EDC3, EDC4 and LSM14A/B. A targeted RNAi screen in intestinal progenitor cells identified 39 previously known and 64 novel P-body regulators, including Patr-1, a gene necessary for P-body assembly. Loss of Patr-1-dependent P-bodies leads to a loss of stem cells that is associated with inappropriate translation and expression of EC-fate gene nubbin. Transcriptomic analysis of progenitor cells identifies a cadre of such weakly transcribed pro-differentiation transcripts that are elevated after P-body loss. Altogether, this study identifies a coordinated P-body dependent, translational and transcriptional repression program that maintains a defined set of in vivo stem cells in a state primed for differentiation.",
      style = "background: #F0F8FF"
    )
  ),
  
  tabPanel(
    "RNA-seq Data",
    fluidPage(
      wellPanel(
        "This page provides differential gene expression datatables associated with",
        a(href = "https://www.biorxiv.org/content/10.1101/2020.06.27.175398v1.abstract",
          "Buddika et al. 2021", style = "color: #DC143C"), "currently in review.",
        p("Use the scroll down menu to change the dataset you wish to access. The study present 
          data generated in the current study and reanalyzed date from",
          a(href = "https://www.embopress.org/doi/full/10.15252/embj.201489072",
            "Korzelius et al. 2014", style = "color: #DC143C"), "and",
          a(href = "https://elifesciences.org/articles/48160#content",
            "Khaminets et al. 2020", style = "color: #DC143C"), "available from public databases."),
        p("Select the dataset and click the download button to download the data to your machine."),
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
  
  tabPanel(
    "Plots",
    fluidPage(
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
  
  tabPanel(
    "Screen",
    fluidPage(
      wellPanel(
        "This page provides results of the targeted P-body genetic screen described in",
        a(href = "https://www.biorxiv.org/content/10.1101/2020.06.27.175398v1.abstract",
          "Buddika et al. 2021", style = "color: #1E90FF"), "manuscript.",
        p("Click the download button below to download the data to your machine."),
        style = "background: #F9D7D7"
      ),
      downloadButton(
        "download_screenData",
        "Download Data"
      ),
      hr(),
      DTOutput("screenData"))),
  collapsible = TRUE
  )



