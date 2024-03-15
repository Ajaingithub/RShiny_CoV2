## Please use the RunApp in the RStudio to launch the App.
## Using this app you can identify your gene of interest in SARS-CoV-2 specific CD4 and CD8 T cells.
## It has 4 Tabs:
## 1. UMAP: UMAP and Violin Plots for CD4 and CD8. All the plots are downloadable in pdf format
## 2. DEG: To perform differential genes between multiple clusters, samples, vaccine etc. All the tables are downloadable in csv format
## 3. DEP: To perform differential protein between multiple clusters, samples, vaccine etc. All the tables are downloadable in csv format
## 4. Gene Score: This tab perform the combine gene score UMAP and Violin plot. All the plots are downloadable in pdf format.

library(shiny)
library(Seurat)
library(session)
library(shinythemes)
library(shinydashboard)
library(shinyWidgets)
library(ggplot2)
library(UCell)

CD4 <- readRDS("./Data/CD4_modality_integrate_cluster_0.8_Diet.RDS")
CD8 <- readRDS("./Data/CD8_modality_integrate_cluster_0.8_Diet.RDS")

source("./Data/FeaturePlotSingle.R")
source("./Data/FeaturePlotSingle_ADT.R")
source("./Data/express_cell_front.R")
source("./Data/cell_express_front_splitted.R")
source("./Data/cell_express_front_splitted_Samples.R")

CD4$Sample_cluster <- paste(CD4$orig.ident,CD4$seurat_clusters,sep = "_")
CD4$vaccine_cluster <- paste(CD4$vaccine,CD4$seurat_clusters,sep = "_")

CD8$Sample_cluster <- paste(CD8$orig.ident,CD8$seurat_clusters,sep = "_")
CD8$vaccine_cluster <- paste(CD8$vaccine,CD8$seurat_clusters,sep = "_")

shinyApp(
  ui = navbarPage("scCITESeq SARS CoV2", theme = shinytheme("flatly"),
                  tabPanel("UMAP",
                           fluidRow(
                             column(2,img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Mayo_Clinic_logo.svg/220px-Mayo_Clinic_logo.svg.png", 
                                          height = 70, width = 70)),
                             column(8, align = "center",h1("Single Cell CiteSeq SARS CoV 2"))
                           ),
                           fluidRow(
                             column(6,wellPanel(selectInput(inputId = "Protein",label = "Choose a Protein",choices = row.names(CD4@assays$ADT@data)),
                                                actionButton(inputId = "go",label = "Update"))),
                             column(6,wellPanel(selectInput(inputId = "Gene",label = "Choose a Gene",choices = row.names(CD4@assays$RNA@data)),
                                                actionButton(inputId = "go2",label = "Update")))),
                           br(),
                           column(12, align = "center",h3("CD4 UMAP RNA and ADT")),
                           column(3,plotOutput("UMAP_wt_CD4_ADT")),
                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                           column(3,tags$img(height = 400, width = 400, src = "CD4_bubble_plot.png")),
                           column(3,plotOutput("UMAP_wt_CD4_RNA")),
                           br(),
                           column(6,plotOutput("UMAP_wt_CD4_vaccine_ADT")),
                           column(6,plotOutput("UMAP_wt_CD4_vaccine_RNA")),
                           br(),
                           column(6,plotOutput("UMAP_wt_CD4_sample_ADT")),
                           column(6,plotOutput("UMAP_wt_CD4_sample_RNA")),
                           br(),
                           downloadButton("downloadPlot_CD4_ADTwt", "Download UMAP ADT CD4 weighted"),
                           downloadButton("downloadPlot_CD4_ADTvln","Download ADT CD4 Violin"),
                           downloadButton("downloadPlot_CD4_RNAwt", "Download UMAP RNA CD4 weighted"),
                           downloadButton("downloadPlot_CD4_RNAvln","Download RNA CD4 Violin"),
                           br(),
                           column(12, align = "center",h3("CD8 UMAP RNA and ADT")),
                           column(3,plotOutput("UMAP_wt_CD8_ADT")),
                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png")),
                           column(3,tags$img(height = 400, width = 400, src = "CD8_bubble_plot.png")),
                           column(3,plotOutput("UMAP_wt_CD8_RNA")),
                           br(),
                           column(6,plotOutput("UMAP_wt_CD8_vaccine_ADT")),
                           column(6,plotOutput("UMAP_wt_CD8_vaccine_RNA")),
                           br(),
                           column(6,plotOutput("UMAP_wt_CD8_sample_ADT")),
                           column(6,plotOutput("UMAP_wt_CD8_sample_RNA")),
                           br(),
                           downloadButton("downloadPlot_CD8_ADTwt", "Download UMAP ADT CD8 weighted"),
                           downloadButton("downloadPlot_CD8_ADTvln","Download ADT CD8 Violin"),
                           downloadButton("downloadPlot_CD8_RNAwt", "Download UMAP RNA CD8 weighted"),
                           downloadButton("downloadPlot_CD8_RNAvln","Download RNA CD8 Violin"),
                           br(),
                           fluidRow(
                             column(12, align = "center",h3("CD4 Vln Plot")),
                             column(6,plotOutput("VlnPlot_CD4_ADT")),
                             column(6,plotOutput("VlnPlot_CD4_RNA"))
                           ),
                           fluidRow(
                             column(12, align = "center",h3("CD8 Vln Plot")),
                             column(6,plotOutput("VlnPlot_CD8_ADT")),
                             column(6,plotOutput("VlnPlot_CD8_RNA"))
                           )
                  ),
                  tabPanel("DEG",
                           dashboardPage(
                             dashboardHeader(title = "RNA Differential"),
                             dashboardSidebar(
                               sidebarMenu(
                                 menuItem("Sample", tabName = "sample", icon = icon("database")),
                                 menuItem("Cluster", tabName = "cluster", icon = icon("database")),
                                 menuItem("Sample_Cluster", tabName = "sample_Cluster", icon = icon("database")),
                                 menuItem("Vaccine", tabName = "vaccine", icon = icon("database")),
                                 menuItem("Vaccine_cluster", tabName = "vaccine_cluster", icon = icon("database"))
                               )
                             ),
                             dashboardBody(
                               tabItems(
                                 tabItem(tabName = "sample",
                                         column(6,wellPanel(selectInput("Sample1", "Sample1",choices = sort(unique(CD4@meta.data$orig.ident)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("Sample2","Sample2",choices = unique(CD4@meta.data$orig.ident),
                                                                        multiple = TRUE),actionButton(inputId = "go3",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_Sample")),
                                           downloadButton("CD4_Sample_download", "Download CD4 marker Data Table"),
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_Sample")),
                                           downloadButton("CD8_Sample_download", "Download CD8 marker Data Table")
                                         ),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 ),
                                 tabItem(tabName = "cluster",
                                         column(6,wellPanel(selectInput("CD4_Cluster1","CD4 Cluster 1",choices = levels(CD4),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("CD4_Cluster2","CD4 Cluster 2",choices = levels(CD4),
                                                                        multiple = TRUE),actionButton(inputId = "CD4_go4",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_Cluster")),
                                           downloadButton("CD4_Cluster_download", "Download CD4 marker Data Table")
                                         ),
                                         column(6,wellPanel(selectInput("CD8_Cluster1","CD8 Cluster 1",choices = levels(CD8),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("CD8_Cluster2","CD8 Cluster 2",choices = levels(CD8),
                                                                        multiple = TRUE),actionButton(inputId = "CD8_go4",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_Cluster")),
                                           downloadButton("CD8_Cluster_download", "Download CD8 marker Data Table")
                                         ),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 ),
                                 tabItem(tabName = "sample_Cluster",
                                         column(6,wellPanel(selectInput("Sample_CD4_Cluster1","CD4 Sample Cluster1",choices = sort(unique(CD4$Sample_cluster)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("Sample_CD4_Cluster2","CD4 Sample Cluster2",choices = sort(unique(CD4$Sample_cluster)),
                                                                        multiple = TRUE),actionButton(inputId = "CD4_go6",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_Sample_Cluster")),
                                           downloadButton("CD4_sample_cluster_download", "Download CD4 Sample Cluster Data Table")
                                         ),
                                         column(6,wellPanel(selectInput("Sample_CD8_Cluster1","CD8 Sample Cluster1",choices = sort(unique(CD8$Sample_cluster)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("Sample_CD8_Cluster2","CD8 Sample Cluster2",choices = sort(unique(CD8$Sample_cluster)),
                                                                        multiple = TRUE),actionButton(inputId = "CD8_go6",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_Sample_Cluster")),
                                           downloadButton("CD8_sample_cluster_download", "Download CD8 Sample Cluster Data Table")
                                         ),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 ),
                                 tabItem(tabName = "vaccine",
                                         column(6,wellPanel(selectInput("Vaccine1","vaccine1",choices = unique(CD4@meta.data$vaccine),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("Vaccine2","vaccine2",choices = unique(CD4@meta.data$vaccine),
                                                                        multiple = TRUE),actionButton(inputId = "go5",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_vaccine")),
                                           downloadButton("CD4_vaccine_download", "Download CD4 Vaccine Data Table"),
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_vaccine")),
                                           downloadButton("CD8_vaccine_download", "Download CD8 Vaccin Data Table")
                                         ),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 ),
                                 tabItem(tabName = "vaccine_cluster",
                                         column(6,wellPanel(selectInput("CD4_Vaccine_Cluster1","CD4 Vaccine Cluster1",choices = sort(unique(CD4@meta.data$vaccine_cluster)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("CD4_Vaccine_Cluster2","CD4 Vaccine Cluster2",choices = sort(unique(CD4@meta.data$vaccine_cluster)),
                                                                        multiple = TRUE),actionButton(inputId = "CD4_go7",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_vaccine_cluster")),
                                           downloadButton("CD4_vaccine_cluster_download", "Download CD4 Vaccine Cluster Data Table")
                                         ),
                                         column(6,wellPanel(selectInput("CD8_Vaccine_Cluster1","CD8 Vaccine Cluster1",choices = sort(unique(CD8@meta.data$vaccine_cluster)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("CD8_Vaccine_Cluster2","CD8 Vaccine Cluster2",choices = sort(unique(CD8@meta.data$vaccine_cluster)),
                                                                        multiple = TRUE),actionButton(inputId = "CD8_go7",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_vaccine_cluster")),
                                           downloadButton("CD8_vaccine_cluster_download", "Download CD8 Vaccine Cluster Data Table")
                                         ),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 )
                               )
                             )
                           )
                  ),
                  tabPanel("DEP",
                           dashboardPage(
                             dashboardHeader(title = "Protein Differential"),
                             dashboardSidebar(
                               sidebarMenu(
                                 menuItem("Sample", tabName = "sample_dep", icon = icon("database")),
                                 menuItem("Cluster", tabName = "cluster_dep", icon = icon("database")),
                                 menuItem("Sample_Cluster", tabName = "sample_Cluster_dep", icon = icon("database")),
                                 menuItem("Vaccine", tabName = "vaccine_dep", icon = icon("database")),
                                 menuItem("Vaccine_cluster", tabName = "vaccine_cluster_dep", icon = icon("database"))
                               )
                             ),
                             dashboardBody(
                               tabItems(
                                 tabItem(tabName = "sample_dep",
                                         column(6,wellPanel(selectInput("Sample1_dep", "Sample1",choices = sort(unique(CD4@meta.data$orig.ident)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("Sample2_dep","Sample2",choices = unique(CD4@meta.data$orig.ident),
                                                                        multiple = TRUE),actionButton(inputId = "go3_dep",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_Sample_dep")),
                                           downloadButton("CD4_Sample_download_dep", "Download CD4 marker DEP Data Table"),
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_Sample_dep")),
                                           downloadButton("CD8_Sample_download_dep", "Download CD8 marker DEP Data Table")
                                         ),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 ),
                                 tabItem(tabName = "cluster_dep",
                                         column(6,wellPanel(selectInput("CD4_Cluster1_dep","CD4 Cluster 1",choices = levels(CD4),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("CD4_Cluster2_dep","CD4 Cluster 2",choices = levels(CD4),
                                                                        multiple = TRUE),actionButton(inputId = "CD4_go4_dep",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_Cluster_dep")),
                                           downloadButton("CD4_Cluster_download_dep", "Download CD4 marker Data Table")
                                         ),
                                         column(6,wellPanel(selectInput("CD8_Cluster1_dep","CD8 Cluster 1",choices = levels(CD8),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("CD8_Cluster2_dep","CD8 Cluster 2",choices = levels(CD8),
                                                                        multiple = TRUE),actionButton(inputId = "CD8_go4_dep",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_Cluster_dep")),
                                           downloadButton("CD8_Cluster_download_dep", "Download CD8 marker Data Table")
                                         ),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 ),
                                 tabItem(tabName = "sample_Cluster_dep",
                                         column(6,wellPanel(selectInput("Sample_CD4_Cluster1_dep","CD4 Sample Cluster1",choices = sort(unique(CD4$Sample_cluster)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("Sample_CD4_Cluster2_dep","CD4 Sample Cluster2",choices = sort(unique(CD4$Sample_cluster)),
                                                                        multiple = TRUE),actionButton(inputId = "CD4_go6_dep",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_Sample_Cluster_dep")),
                                           downloadButton("CD4_sample_cluster_download_dep", "Download CD4 Sample Cluster Data Table")),
                                         column(6,wellPanel(selectInput("Sample_CD8_Cluster1_dep","CD8 Sample Cluster1",choices = sort(unique(CD8$Sample_cluster)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("Sample_CD8_Cluster2_dep","CD8 Sample Cluster2",choices = sort(unique(CD8$Sample_cluster)),
                                                                        multiple = TRUE),actionButton(inputId = "CD8_go6_dep",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_Sample_Cluster_dep")),
                                           downloadButton("CD8_sample_cluster_download_dep", "Download CD8 Sample Cluster Data Table")),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 ),
                                 tabItem(tabName = "vaccine_dep",
                                         column(6,wellPanel(selectInput("Vaccine1_dep","vaccine1",choices = unique(CD4@meta.data$vaccine),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("Vaccine2_dep","vaccine2",choices = unique(CD4@meta.data$vaccine),
                                                                        multiple = TRUE),actionButton(inputId = "go5_dep",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_vaccine_dep")),
                                           downloadButton("CD4_vaccine_download_dep", "Download CD4 Vaccine Data Table"),
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_vaccine_dep")),
                                           downloadButton("CD8_vaccine_download_dep", "Download CD8 Vaccin Data Table")),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 ),
                                 tabItem(tabName = "vaccine_cluster_dep",
                                         column(6,wellPanel(selectInput("CD4_Vaccine_Cluster1_dep","CD4 Vaccine Cluster1",choices = sort(unique(CD4@meta.data$vaccine_cluster)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("CD4_Vaccine_Cluster2_dep","CD4 Vaccine Cluster2",choices = sort(unique(CD4@meta.data$vaccine_cluster)),
                                                                        multiple = TRUE),actionButton(inputId = "CD4_go7_dep",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD4 Data Table")),
                                           column(12,dataTableOutput("CD4_vaccine_cluster_dep")),
                                           downloadButton("CD4_vaccine_cluster_download_dep", "Download CD4 Vaccine Cluster Data Table")
                                         ),
                                         column(6,wellPanel(selectInput("CD8_Vaccine_Cluster1_dep","CD8 Vaccine Cluster1",choices = sort(unique(CD8@meta.data$vaccine_cluster)),
                                                                        multiple = TRUE))),
                                         column(6,wellPanel(selectInput("CD8_Vaccine_Cluster2_dep","CD8 Vaccine Cluster2",choices = sort(unique(CD8@meta.data$vaccine_cluster)),
                                                                        multiple = TRUE),actionButton(inputId = "CD8_go7_dep",label = "Update"))),
                                         fluidRow(
                                           column(12, align = "center",h3("CD8 Data Table")),
                                           column(12,dataTableOutput("CD8_vaccine_cluster_dep")),
                                           downloadButton("CD8_vaccine_cluster_download_dep", "Download CD8 Vaccine Cluster Data Table")),
                                         fluidRow(
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                                           column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                                         )
                                 )
                               )
                             )
                           )
                  ),
                  tabPanel("GeneScore",
                           fluidRow(
                             column(2,img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Mayo_Clinic_logo.svg/220px-Mayo_Clinic_logo.svg.png", height = 70, width = 70)),
                             column(8, align = "center",h1("GeneScore Analysis"))
                           ),
                           fluidRow(
                             column(6,fileInput("upload_GeneScoregene1", "Upload Upregulated Genes")),
                             column(6,fileInput("upload_GeneScoregene2", "Upload Downregulated Genes"), actionButton(inputId = "upload_GeneScorego",label = "Update"))
                           ),
                           fluidRow(
                             column(6, plotOutput("upload_GeneScore_CD4_RNA")),
                             column(6, plotOutput("upload_GeneScore_CD8_RNA"))
                           ),
                           fluidRow(
                             column(12, h3("These Genes were absent from CD4 or CD8")),
                             column(12,dataTableOutput("absent_gene_table"))
                           ),
                           fluidRow(
                             column(3,tags$img(height = 400, width = 400, src = "CD4_Sample_split_vaccine.png")),
                             column(3,tags$img(height = 400, width = 400, src = "CD4_new_seurat_cluster.png")),
                             column(3,tags$img(height = 400, width = 400, src = "CD8_Sample_split_vaccine.png")),
                             column(3,tags$img(height = 400, width = 400, src = "CD8_new_seurat_cluster.png"))
                           ),
                           fluidRow(
                             column(6,wellPanel(selectizeInput(inputId = "GeneScoregene1",label = "Upregulated Genes", choices = row.names(CD4@assays$RNA@data), multiple = TRUE))),
                             column(6,wellPanel(selectizeInput(inputId = "GeneScoregene2",label = "Downregulated Genes", choices = row.names(CD4@assays$RNA@data), multiple = TRUE)),
                                    actionButton(inputId = "GeneScorego",label = "Update"))),
                           fluidRow(
                             column(6, plotOutput("GeneScore_CD4_RNA")),
                             column(6, plotOutput("GeneScore_CD8_RNA"))
                           )
                  )
  ),
  server = function(input, output){
    CD4_data4 <- eventReactive(input$go, {
      DefaultAssay(CD4) <- "MAGIC_ADT"
      featureplot_front(CD4, feature = input$Protein, reduction = "wnn.umap", x = "wnnUMAP_1", y = "wnnUMAP_2",
                        color = "darkgreen") + ggtitle(paste(input$Protein))
    })
    output$UMAP_wt_CD4_ADT <- renderPlot({
      CD4_data4()
    })

    CD8_data4 <- eventReactive(input$go, {
      DefaultAssay(CD8) <- "MAGIC_ADT"
      featureplot_front(CD8, feature = input$Protein, reduction = "wnn.umap", x = "wnnUMAP_1", y = "wnnUMAP_2",
                        color = "darkgreen") + ggtitle(paste(input$Protein))
    })
    output$UMAP_wt_CD8_ADT <- renderPlot({
      CD8_data4()
    })

    output$downloadPlot_CD4_ADTwt <- downloadHandler(
      filename = function() { paste(input$Protein,"wUMAP_CD4_protein_impute.pdf", sep='_') },
      content = function(file) {
        #device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = CD4_data4(), device = "pdf", width = 5, height = 5)
      })

    output$downloadPlot_CD8_ADTwt <- downloadHandler(
      filename = function() { paste(input$Protein,"wUMAP_CD8_protein_impute.pdf", sep='_') },
      content = function(file) {
        #device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = CD8_data4(), device = "pdf", width = 5, height = 5)
      })

    CD4_data5 <- eventReactive(input$go, {
      DefaultAssay(CD4) <- "ADT"
      VlnPlot(CD4, input$Protein, pt.size = 0.00001, split.by = "vaccine") + ggplot2::scale_fill_manual(values = c('gray60',"maroon1", 'lightseagreen'))
    })

    output$VlnPlot_CD4_ADT <- renderPlot({
      CD4_data5()
    })

    CD8_data5 <- eventReactive(input$go, {
      DefaultAssay(CD8) <- "ADT"
      VlnPlot(CD8, input$Protein, pt.size = 0.00001, split.by = "vaccine") + ggplot2::scale_fill_manual(values = c('gray60',"maroon1", 'lightseagreen'))
    })

    output$VlnPlot_CD8_ADT <- renderPlot({
      CD8_data5()
    })

    output$downloadPlot_CD4_ADTvln <- downloadHandler(
      filename = function() { paste(input$Protein,"vln_protein.pdf", sep='_') },
      content = function(file) {
        #device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = CD4_data5(), device = "pdf", width = 9, height = 7)
      })

    output$downloadPlot_CD8_ADTvln <- downloadHandler(
      filename = function() { paste(input$Protein,"vln_protein.pdf", sep='_') },
      content = function(file) {
        #device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = CD8_data5(), device = "pdf", width = 9, height = 7)
      })

    CD4_data3 <- eventReactive(input$go2, {
      DefaultAssay(CD4) <- "MAGIC_RNA"
      featureplot_front(CD4, feature = input$Gene, reduction = "wnn.umap", x = "wnnUMAP_1", y = "wnnUMAP_2",
                        color = "blue") + ggtitle(paste(input$Gene))
    })

    output$UMAP_wt_CD4_RNA <- renderPlot({
      CD4_data3()
    })

    CD8_data3 <- eventReactive(input$go2, {
      DefaultAssay(CD8) <- "MAGIC_RNA"
      featureplot_front(CD8, feature = input$Gene, reduction = "wnn.umap", x = "wnnUMAP_1", y = "wnnUMAP_2",
                        color = "blue") + ggtitle(paste(input$Gene))
    })

    output$UMAP_wt_CD8_RNA <- renderPlot({
      CD8_data3()
    })

    output$downloadPlot_CD4_RNAwt <- downloadHandler(
      filename = function() { paste(input$Gene,"wtUMAP_CD4__RNA.pdf", sep='_') },
      content = function(file) {
        #device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = CD4_data3(), device = "pdf", width = 5, height = 5)
      })

    output$downloadPlot_CD8_RNAwt <- downloadHandler(
      filename = function() { paste(input$Gene,"wtUMAP_CD8_RNA.pdf", sep='_') },
      content = function(file) {
        #device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = CD8_data3(), device = "pdf", width = 5, height = 5)
      })

    CD4_data6 <- eventReactive(input$go2, {
      DefaultAssay(CD4) <- "RNA"
      VlnPlot(CD4, input$Gene, pt.size = 0.00001, split.by = "vaccine") + ggplot2::scale_fill_manual(values = c('gray60',"maroon1", 'lightseagreen'))
    })
    output$VlnPlot_CD4_RNA <- renderPlot({
      CD4_data6()
    })

    CD8_data6 <- eventReactive(input$go2, {
      DefaultAssay(CD8) <- "RNA"
      VlnPlot(CD8, input$Gene, pt.size = 0.00001, split.by = "vaccine") + ggplot2::scale_fill_manual(values = c('gray60',"maroon1", 'lightseagreen'))
    })
    output$VlnPlot_CD8_RNA <- renderPlot({
      CD8_data6()
    })

    output$downloadPlot_CD4_RNAvln <- downloadHandler(
      filename = function() { paste(input$Protein,"CD4_vln_RNA.pdf", sep='_') },
      content = function(file) {
        #device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = CD4_data6(), device = "pdf", width = 9, height = 7)
      })

    output$downloadPlot_CD8_RNAvln <- downloadHandler(
      filename = function() { paste(input$Protein,"CD8_vln_RNA.pdf", sep='_') },
      content = function(file) {
        #device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = CD8_data6(), device = "pdf", width = 9, height = 7)
      })

    ## Vaccine UMAP
    CD4_vaccine_data1 <- eventReactive(input$go,{
      DefaultAssay(CD4) <- "MAGIC_ADT"
      featureplot_front_split(CD4,input$Protein, reduction = "wnn.umap",x = "wnnUMAP_1", y = "wnnUMAP_2",
                              color = "darkgreen") + ggtitle(paste(input$Protein))
    })

    output$UMAP_wt_CD4_vaccine_ADT <- renderPlot({
      CD4_vaccine_data1()
    })

    CD8_vaccine_data1 <- eventReactive(input$go,{
      DefaultAssay(CD8) <- "MAGIC_ADT"
      featureplot_front_split(CD8,input$Protein, reduction = "wnn.umap",x = "wnnUMAP_1", y = "wnnUMAP_2",
                              color = "darkgreen") + ggtitle(paste(input$Protein))
    })

    output$UMAP_wt_CD8_vaccine_ADT <- renderPlot({
      CD8_vaccine_data1()
    })

    ### RNA
    CD4_vaccine_data2 <- eventReactive(input$go2,{
      DefaultAssay(CD4) <- "MAGIC_RNA"
      featureplot_front_split(CD4, feature = input$Gene, reduction = "wnn.umap", x = "wnnUMAP_1", y = "wnnUMAP_2",
                        color = "blue") + ggtitle(paste(input$Gene))
    })

    output$UMAP_wt_CD4_vaccine_RNA <- renderPlot({
      CD4_vaccine_data2()
    })

    CD8_vaccine_data2 <- eventReactive(input$go2,{
      DefaultAssay(CD8) <- "MAGIC_RNA"
      featureplot_front_split(CD8, feature = input$Gene, reduction = "wnn.umap", x = "wnnUMAP_1", y = "wnnUMAP_2",
                        color = "blue") + ggtitle(paste(input$Gene))
    })

    output$UMAP_wt_CD8_vaccine_RNA <- renderPlot({
      CD8_vaccine_data2()
    })

    ## Sample UMAP
    CD4_sample_data1 <- eventReactive(input$go,{
      DefaultAssay(CD4) <- "MAGIC_ADT"
      featureplot_front_split_sample(CD4,input$Protein, reduction = "wnn.umap",x = "wnnUMAP_1", y = "wnnUMAP_2",
                              color = "darkgreen",nrow=3, ncol=5) + ggtitle(paste(input$Protein))
      # FeaturePlot(CD4,features = input$Protein, reduction = "wnn.umap",split.by = "orig.ident", cols = c("lightgrey","darkgreen"), pt.size = 0.000001) +
      #   patchwork::plot_layout(ncol = 5, nrow = 3)
      # # & theme(legend.position = c(0.9,0.2))
    })

    output$UMAP_wt_CD4_sample_ADT <- renderPlot({
      CD4_sample_data1()
    })

    CD8_sample_data1 <- eventReactive(input$go,{
      DefaultAssay(CD8) <- "MAGIC_ADT"
      featureplot_front_split_sample(CD8,input$Protein, reduction = "wnn.umap",x = "wnnUMAP_1", y = "wnnUMAP_2",
                                     color = "darkgreen",nrow=3, ncol=5) + ggtitle(paste(input$Protein))
      # FeaturePlot(CD8,features = input$Protein, reduction = "wnn.umap", split.by = "orig.ident", cols = c("lightgrey","darkgreen"), pt.size = 0.000001) +
      #   patchwork::plot_layout(ncol = 5, nrow = 3)
      # & theme(legend.position = c(0.9,0.2))
    })

    output$UMAP_wt_CD8_sample_ADT <- renderPlot({
      CD8_sample_data1()
    })

    ### RNA
    CD4_sample_data2 <- eventReactive(input$go2,{
      DefaultAssay(CD4) <- "MAGIC_RNA"
      featureplot_front_split_sample(CD4,input$Gene, reduction = "wnn.umap",x = "wnnUMAP_1", y = "wnnUMAP_2",
                                     color = "blue",nrow=3, ncol=5) + ggtitle(paste(input$Gene))
      # FeaturePlot(CD4,features = input$Gene, reduction = "wnn.umap", split.by = "orig.ident", cols = c("lightgrey","blue"), pt.size = 0.000001) +
      #   patchwork::plot_layout(ncol = 5, nrow = 3)
      # & theme(legend.position = c(0.9,0.2))
    })

    output$UMAP_wt_CD4_sample_RNA <- renderPlot({
      CD4_sample_data2()
    })

    CD8_sample_data2 <- eventReactive(input$go2,{
      DefaultAssay(CD8) <- "MAGIC_RNA"
      featureplot_front_split_sample(CD8,input$Gene, reduction = "wnn.umap",x = "wnnUMAP_1", y = "wnnUMAP_2",
                                     color = "blue",nrow=3, ncol=5) + ggtitle(paste(input$Gene))
      # FeaturePlot(CD8,features = input$Gene, reduction = "wnn.umap", split.by = "orig.ident", cols = c("lightgrey","blue"), pt.size = 0.000001) +
      #   patchwork::plot_layout(ncol = 5, nrow = 3)
      # & theme(legend.position = c(0.9,0.2))
    })

    output$UMAP_wt_CD8_sample_RNA <- renderPlot({
      CD8_sample_data2()
    })


    ###### DEG #####
    # Samples
    CD4_df_data1 <- eventReactive(input$go3, {
      DefaultAssay(CD4) <- "RNA"
      fm <- FindMarkers(CD4, ident.1 = input$Sample1, ident.2 = input$Sample2, group.by = "orig.ident", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })
    output$CD4_Sample <- renderDataTable({
      CD4_df_data1()
    })

    CD8_df_data1 <- eventReactive(input$go3, {
      DefaultAssay(CD8) <- "RNA"
      fm <- FindMarkers(CD8, ident.1 = input$Sample1, ident.2 = input$Sample2, group.by = "orig.ident", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })
    output$CD8_Sample <- renderDataTable({
      CD8_df_data1()
    })

    output$CD4_Sample_download <- downloadHandler(
      filename = function(){paste(input$Sample1,"vs",input$Sample2,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data1(), fname)
      }
    )

    output$CD8_Sample_download <- downloadHandler(
      filename = function(){paste(input$Sample1,"vs",input$Sample2,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data1(), fname)
      }
    )

    # Cluster
    CD4_df_data2 <- eventReactive(input$CD4_go4, {
      DefaultAssay(CD4) <- "RNA"
      fm <- FindMarkers(CD4, ident.1 = input$CD4_Cluster1, ident.2 = input$CD4_Cluster2, group.by = "seurat_clusters", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD4_Cluster <- renderDataTable({
      CD4_df_data2()
    })

    CD8_df_data2 <- eventReactive(input$CD8_go4, {
      DefaultAssay(CD8) <- "RNA"
      fm <- FindMarkers(CD8, ident.1 = input$CD8_Cluster1, ident.2 = input$CD8_Cluster2, group.by = "seurat_clusters", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD8_Cluster <- renderDataTable({
      CD8_df_data2()
    })

    output$CD4_Cluster_download <- downloadHandler(
      filename = function(){paste(input$CD4_Cluster1,"vs",input$CD4_Cluster2,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data2(), fname)
      }
    )

    output$CD8_Cluster_download <- downloadHandler(
      filename = function(){paste(input$CD8_Cluster1,"vs",input$CD8_Cluster2,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data2(), fname)
      }
    )

    # Sample Cluster
    CD4_df_data3 <- eventReactive(input$CD4_go6, {
      DefaultAssay(CD4) <- "RNA"
      fm <- FindMarkers(CD4, ident.1 = input$Sample_CD4_Cluster1, ident.2 = input$Sample_CD4_Cluster2, group.by = "Sample_cluster", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD4_Sample_Cluster <- renderDataTable({
      CD4_df_data3()
    })

    CD8_df_data3 <- eventReactive(input$CD8_go6, {
      DefaultAssay(CD8) <- "RNA"
      fm <- FindMarkers(CD8, ident.1 = input$Sample_CD8_Cluster1, ident.2 = input$Sample_CD8_Cluster2, group.by = "Sample_cluster", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD8_Sample_Cluster <- renderDataTable({
      CD8_df_data3()
    })

    output$CD4_sample_cluster_download <- downloadHandler(
      filename = function(){paste(input$Sample_CD4_Cluster1,"vs",input$Sample_CD4_Cluster2,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data3(), fname)
      }
    )

    output$CD8_sample_cluster_download <- downloadHandler(
      filename = function(){paste(input$Sample_CD8_Cluster1,"vs",input$Sample_CD8_Cluster2,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data3(), fname)
      }
    )

    # vaccine
    CD4_df_data4 <- eventReactive(input$go5, {
      DefaultAssay(CD4) <- "RNA"
      fm <- FindMarkers(CD4, ident.1 = input$Vaccine1, ident.2 = input$Vaccine2, group.by = "vaccine", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD4_vaccine <- renderDataTable({
      CD4_df_data4()
    })

    CD8_df_data4 <- eventReactive(input$go5, {
      DefaultAssay(CD8) <- "RNA"
      fm <- FindMarkers(CD8, ident.1 = input$Vaccine1, ident.2 = input$Vaccine2, group.by = "vaccine", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD8_vaccine <- renderDataTable({
      CD8_df_data4()
    })

    output$CD4_vaccine_download <- downloadHandler(
      filename = function(){paste(input$Vaccine1,"vs",input$Vaccine2,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data4(), fname)
      }
    )

    output$CD8_vaccine_download <- downloadHandler(
      filename = function(){paste(input$Vaccine1,"vs",input$Vaccine2,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data4(), fname)
      }
    )

    # Vaccine Cluster
    CD4_df_data5 <- eventReactive(input$CD4_go7, {
      DefaultAssay(CD4) <- "RNA"
      fm <- FindMarkers(CD4, ident.1 = input$CD4_Vaccine_Cluster1, ident.2 = input$CD4_Vaccine_Cluster2, group.by = "vaccine_cluster", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD4_vaccine_cluster <- renderDataTable({
      CD4_df_data5()
    })

    CD8_df_data5 <- eventReactive(input$CD8_go7, {
      DefaultAssay(CD8) <- "RNA"
      fm <- FindMarkers(CD8, ident.1 = input$CD8_Vaccine_Cluster1, ident.2 = input$CD8_Vaccine_Cluster2, group.by = "vaccine_cluster", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD8_vaccine_cluster <- renderDataTable({
      CD8_df_data5()
    })

    output$CD4_vaccine_cluster_download <- downloadHandler(
      filename = function(){paste(input$CD4_Vaccine_Cluster1,"vs",input$CD4_Vaccine_Cluster2,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data5(), fname)
      }
    )

    output$CD8_vaccine_cluster_download <- downloadHandler(
      filename = function(){paste(input$CD8_Vaccine_Cluster1,"vs",input$CD8_Vaccine_Cluster2,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data5(), fname)
      }
    )

    ###### DEP #####
    # Samples
    CD4_df_data1_dep <- eventReactive(input$go3_dep, {
      DefaultAssay(CD4) <- "ADT"
      fm <- FindMarkers(CD4, ident.1 = input$Sample1_dep, ident.2 = input$Sample2_dep, group.by = "orig.ident", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })
    output$CD4_Sample_dep <- renderDataTable({
      CD4_df_data1_dep()
    })

    CD8_df_data1_dep <- eventReactive(input$go3_dep, {
      DefaultAssay(CD8) <- "ADT"
      fm <- FindMarkers(CD8, ident.1 = input$Sample1_dep, ident.2 = input$Sample2_dep, group.by = "orig.ident", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })
    output$CD8_Sample_dep <- renderDataTable({
      CD8_df_data1_dep()
    })

    output$CD4_Sample_download_dep <- downloadHandler(
      filename = function(){paste(input$Sample1_dep,"vs",input$Sample2_dep,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data1_dep(), fname)
      })

    output$CD8_Sample_download_dep <- downloadHandler(
      filename = function(){paste(input$Sample1_dep,"vs",input$Sample2_dep,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data1_dep(), fname)
      }
    )

    # Cluster
    CD4_df_data2_dep <- eventReactive(input$CD4_go4_dep, {
      DefaultAssay(CD4) <- "ADT"
      fm <- FindMarkers(CD4, ident.1 = input$CD4_Cluster1_dep, ident.2 = input$CD4_Cluster2_dep, group.by = "seurat_clusters", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD4_Cluster_dep <- renderDataTable({
      CD4_df_data2_dep()
    })

    CD8_df_data2_dep <- eventReactive(input$CD8_go4_dep, {
      DefaultAssay(CD8) <- "ADT"
      fm <- FindMarkers(CD8, ident.1 = input$CD8_Cluster1_dep, ident.2 = input$CD8_Cluster2_dep, group.by = "seurat_clusters", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD8_Cluster_dep <- renderDataTable({
      CD8_df_data2_dep()
    })

    output$CD4_Cluster_download_dep <- downloadHandler(
      filename = function(){paste(input$CD4_Cluster1_dep,"vs",input$CD4_Cluster2_dep,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data2_dep(), fname)
      }
    )

    output$CD8_Cluster_download_dep <- downloadHandler(
      filename = function(){paste(input$CD8_Cluster1_dep,"vs",input$CD8_Cluster2_dep,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data2_dep(), fname)
      }
    )

    # Sample Cluster
    CD4_df_data3_dep <- eventReactive(input$CD4_go6_dep, {
      DefaultAssay(CD4) <- "ADT"
      fm <- FindMarkers(CD4, ident.1 = input$Sample_CD4_Cluster1_dep, ident.2 = input$Sample_CD4_Cluster2_dep, group.by = "Sample_cluster", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD4_Sample_Cluster_dep <- renderDataTable({
      CD4_df_data3_dep()
    })

    CD8_df_data3_dep <- eventReactive(input$CD8_go6_dep, {
      DefaultAssay(CD8) <- "ADT"
      fm <- FindMarkers(CD8, ident.1 = input$Sample_CD8_Cluster1_dep, ident.2 = input$Sample_CD8_Cluster2_dep, group.by = "Sample_cluster", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD8_Sample_Cluster_dep <- renderDataTable({
      CD8_df_data3_dep()
    })

    output$CD4_sample_cluster_download_dep <- downloadHandler(
      filename = function(){paste(input$Sample_CD4_Cluster1_dep,"vs",input$Sample_CD4_Cluster2_dep,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data3_dep(), fname)
      }
    )

    output$CD8_sample_cluster_download_dep <- downloadHandler(
      filename = function(){paste(input$Sample_CD8_Cluster1_dep,"vs",input$Sample_CD8_Cluster2_dep,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data3_dep(), fname)
      }
    )

    # vaccine
    CD4_df_data4_dep <- eventReactive(input$go5_dep, {
      DefaultAssay(CD4) <- "ADT"
      fm <- FindMarkers(CD4, ident.1 = input$Vaccine1_dep, ident.2 = input$Vaccine2_dep, group.by = "vaccine", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD4_vaccine_dep <- renderDataTable({
      CD4_df_data4_dep()
    })

    CD8_df_data4_dep <- eventReactive(input$go5_dep, {
      DefaultAssay(CD8) <- "ADT"
      fm <- FindMarkers(CD8, ident.1 = input$Vaccine1_dep, ident.2 = input$Vaccine2_dep, group.by = "vaccine", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD8_vaccine_dep <- renderDataTable({
      CD8_df_data4_dep()
    })

    output$CD4_vaccine_download_dep <- downloadHandler(
      filename = function(){paste(input$Vaccine1_dep,"vs",input$Vaccine2_dep,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data4_dep(), fname)
      }
    )

    output$CD8_vaccine_download_dep <- downloadHandler(
      filename = function(){paste(input$Vaccine1_dep,"vs",input$Vaccine2_dep,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data4_dep(), fname)
      }
    )

    # Vaccine Cluster
    CD4_df_data5_dep <- eventReactive(input$CD4_go7_dep, {
      DefaultAssay(CD4) <- "ADT"
      fm <- FindMarkers(CD4, ident.1 = input$CD4_Vaccine_Cluster1_dep, ident.2 = input$CD4_Vaccine_Cluster2_dep, group.by = "vaccine_cluster", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD4_vaccine_cluster_dep <- renderDataTable({
      CD4_df_data5_dep()
    })

    CD8_df_data5_dep <- eventReactive(input$CD8_go7_dep, {
      DefaultAssay(CD8) <- "ADT"
      fm <- FindMarkers(CD8, ident.1 = input$CD8_Vaccine_Cluster1_dep, ident.2 = input$CD8_Vaccine_Cluster2_dep, group.by = "vaccine_cluster", min.pct = 0, logfc.threshold = 0.1)
      fm$Gene <- rownames(fm)
      fm
    })

    output$CD8_vaccine_cluster_dep <- renderDataTable({
      CD8_df_data5_dep()
    })

    output$CD4_vaccine_cluster_download_dep <- downloadHandler(
      filename = function(){paste(input$CD4_Vaccine_Cluster1_dep,"vs",input$CD4_Vaccine_Cluster2_dep,"CD4_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD4_df_data5_dep(), fname)
      })

    output$CD8_vaccine_cluster_download_dep <- downloadHandler(
      filename = function(){paste(input$CD8_Vaccine_Cluster1_dep,"vs",input$CD8_Vaccine_Cluster2_dep,"CD8_markers.csv", sep='_')},
      content = function(fname){
        write.csv(CD8_df_data5_dep(), fname)
      })

    ### GeneScore Analysis ####
    GeneScore_analysis_CD4 <- eventReactive(input$GeneScorego, {
      DefaultAssay(CD4) <- "RNA"
      markers <- list()
      markers$GeneScore <- c(input$GeneScoregene1,paste(input$GeneScoregene2,"-",sep = ""))
      CD4 <- AddModuleScore_UCell(CD4, features = markers)
      signature.names <- paste0(names(markers), "_UCell")
      p1 <- VlnPlot(CD4, features = signature.names, group.by = "seurat_clusters")
      p2 <- FeaturePlot(CD4, features = signature.names, reduction = "wnn.umap")
      p1|p2
    })

    output$GeneScore_CD4_RNA <- renderPlot({
      GeneScore_analysis_CD4()
    })

    GeneScore_analysis_CD8 <- eventReactive(input$GeneScorego, {
      DefaultAssay(CD8) <- "RNA"
      markers <- list()
      markers$GeneScore <- c(input$GeneScoregene1,paste(input$GeneScoregene2,"-",sep = ""))
      CD8 <- AddModuleScore_UCell(CD8, features = markers)
      signature.names <- paste0(names(markers), "_UCell")
      p1 <- VlnPlot(CD8, features = signature.names, group.by = "seurat_clusters")
      p2 <- FeaturePlot(CD8, features = signature.names, reduction = "wnn.umap")
      p1|p2
    })

    output$GeneScore_CD8_RNA <- renderPlot({
      GeneScore_analysis_CD8()
    })

    upload_GeneScore_analysis_CD4 <- eventReactive(input$upload_GeneScorego, {
      DefaultAssay(CD4) <- "RNA"
      markers <- list()
      path_gene1 = read.table(input$upload_GeneScoregene1$datapath, header = FALSE, sep = "\t")
      path_gene2 = read.table(input$upload_GeneScoregene2$datapath, header = FALSE, sep = "\t")
      markers$GeneScore <- c(path_gene1[,1],paste(path_gene2[,1],"-",sep = ""))
      CD4 <- AddModuleScore_UCell(CD4, features = markers)
      signature.names <- paste0(names(markers), "_UCell")
      p1 <- VlnPlot(CD4, features = signature.names, group.by = "seurat_clusters")
      p2 <- FeaturePlot(CD4, features = signature.names, reduction = "wnn.umap") + ggtitle("CD4")
      p1|p2
    })

    output$upload_GeneScore_CD4_RNA <- renderPlot({
      upload_GeneScore_analysis_CD4()
    })

    upload_GeneScore_analysis_CD8 <- eventReactive(input$upload_GeneScorego, {
      DefaultAssay(CD8) <- "RNA"
      markers <- list()
      path_gene1 = read.table(input$upload_GeneScoregene1$datapath, header = FALSE, sep = "\t")
      path_gene2 = read.table(input$upload_GeneScoregene2$datapath, header = FALSE, sep = "\t")
      markers$GeneScore <- c(path_gene1[,1],paste(path_gene2[,1],"-",sep = ""))
      CD8 <- AddModuleScore_UCell(CD8, features = markers)
      signature.names <- paste0(names(markers), "_UCell")
      p1 <- VlnPlot(CD8, features = signature.names, group.by = "seurat_clusters")
      p2 <- FeaturePlot(CD8, features = signature.names, reduction = "wnn.umap") + ggtitle("CD8")
      p1|p2
    })

    output$upload_GeneScore_CD8_RNA <- renderPlot({
      upload_GeneScore_analysis_CD8()
    })

    upload_GeneScore_gene <- eventReactive(input$upload_GeneScorego, {
      path_gene1 = read.table(input$upload_GeneScoregene1$datapath, header = FALSE, sep = "\t")
      path_gene2 = read.table(input$upload_GeneScoregene2$datapath, header = FALSE, sep = "\t")
      GeneScore_genes = c(path_gene1[,1],path_gene2[,1])
      CD4_genes = rownames(CD4@assays$RNA@counts)
      CD8_genes = rownames(CD8@assays$RNA@counts)
      `%notin%` <- Negate(`%in%`)
      CD4_GeneScore_genes_absent = GeneScore_genes[which(GeneScore_genes %notin% CD4_genes)]
      CD8_GeneScore_genes_absent = GeneScore_genes[which(GeneScore_genes %notin% CD8_genes)]
      genes = as.data.frame(c(CD4_GeneScore_genes_absent,CD8_GeneScore_genes_absent))
      colnames(genes) = "Genes"
      return(genes)
    })

    output$absent_gene_table <- renderDataTable({
      upload_GeneScore_gene()
    })
  }
)

