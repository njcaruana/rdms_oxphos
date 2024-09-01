
#library(BiocManager)
#getOption("repos")

#if (!require("pacman")) install.packages("pacman")
#pacman::p_load("shinymanager","tidyverse","vroom","DT", "ggplot2", "ggpubr", "ggiraph", "plotly","EnhancedVolcano","RColorBrewer", "shinyjs")

library(vroom)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DT)
library(ggiraph)
library(plotly)
library(EnhancedVolcano)
library(RColorBrewer)
library(shinyjs)
library(shinymanager)
library(gtools)


datalist <- list.files(path = "data/rca-data/",pattern=NULL, full.names = FALSE)
datalist <- gsub('.csv', '', datalist)
#datalist <- gsub('\\.', '-', datalist)
#datalist <- sort(datalist, decreasing = TRUE )
datalist <- mixedsort(datalist)

#inactivity script to timeout login page after 2 mins
inactivity <- "function idleTimer() {
var t = setTimeout(logout, 6000000);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
window.close();  //close the window
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, 6000000);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();"

# data.frame with credentials info
# credentials <- data.frame(
#   user = c("RDMS"),
#   password = c("Mito1!23"),
#   stringsAsFactors = FALSE
# )

#Function for RCA stats
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  ) 
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Function for calculating the outliers within the uploaded data.
outlier_fun <- function(m){
  m2 <- as.data.frame(sapply(m, function(m) (abs(m-mean(m))/sd(m))))
  colnames(m2) <- colnames(m)
  rownames(m2) <- rownames(m)
  zscore <- m2 %>% mutate(outliercheck = !rowSums(m2>3), )
  zscore$uniprot <- rownames(zscore)
  zscore = zscore %>% filter(outliercheck == TRUE)
  zscore = zscore[c("uniprot")]
  return(zscore)
}

#Function for loading page
# load_data <- function() {
#   Sys.sleep(2)
#   hide("loading_page")
#   show("main_content")
# }

#linebreaks
linebreaks <- function(n){HTML(strrep(br(), n))}


ui <- shinyUI(#secure_app(head_auth = tags$script(inactivity),
                         fluidPage(
                           tags$head(includeHTML("google-analytics.html")),
                           useShinyjs(),
                           theme = bslib::bs_theme(bootswatch = "lux"),
                           br(),
                           img(src='RDMSExplorer_logo.png', align = "center",height='150px'),
                           br(),
                           br(),
                           #em(strong("Exploration of Rare Disease Proteomics Data"), style = "font-size:20px;"),
                           #br(),
                           #br(),
                           tags$footer(strong("RDMSExplorer is a dynamic website which is constantly being updated. Citation: Untargeted proteomics enables ultra-rapid variant prioritization in mitochondrial and other rare diseases. Hock et al. 2024. Current Medrxiv at: https://shorturl.at/s85s7"),
                                       # strong() = bold
                                       align = "left",
                                       style = "
                           position:fixed;
                           bottom:11.5px;
                           width:100%;
                           height:40px;
                           color: white;
                           padding: 10px;
                           opacity: 0.7;
                           font-size: 12px;
                           background-color: #2A6C8C;
                           z-index: 100;
              "),
                           tabsetPanel(
                             
                             tabPanel("Project Overview",
                                      fluidPage(fluidRow(class = "text-center",
                                                         column(width = 12,
                                                                column(5),
                                                                
                                                                br(),
                                                                br(),
                                                                HTML('<footer> <img src="RDMSopener.png", align = "center", height = "950cm"</img> </footer>'),
                                                                linebreaks(7),
                                                                HTML('<footer> <img src="RDMS_logo-opacity.png", align = "center", height = "80px"</img> </footer>'),
                                                                br())
                                      )
                                      )
                             ),
                             #Tab Panel for Patient Plots with uploaded data. 
                             tabPanel("RCA Plots",
                                      br(),
                                     # h2(strong("Exploration of Data"), style = "font-size:20px;"),
                                      fluidRow(
                                        column(12, align = "center",
                                               br(),
                                               selectInput('inputFile', label='Select or search for data:',
                                                           choice= datalist, width = "300px"),
                                               br()
                                        )),
                                      fluidRow(column(12, align = "center",
                                                      #textOutput("mitofactor"),
                                                      br(),
                                                      girafeOutput('p_RCAPlot'),
                                                      br()
                                      )), 
                                      fluidRow(
                                        column(12, align = "center",
                                               br(),
                                               p("ns: p > 0.05, *: p <= 0.05, **: p <= 0.01, ***: p <= 0.001, ****: p <= 0.0001"),
                                        )),
                                      fluidRow(
                                        column(12 , align = "center",
                                               br(),
                                               tableOutput("p_rcadata"),
                                               br(),
                                               br(),
                                               p("Abbreviations for samples: VC - Validation Cohort. KC - Knock Out Cohort. UDP - Undiagnosed Patient. SC - Supporting Cohort."),
                                               br(),
                                               br(),
                                               p("For Relative Complex Abundance (RCA) plot, 
                                   MS2 intensity from proteins belonging to OXPHOS 
                                   complexes were plotted in R calculating 
                                   the difference between the sample of interest and control MS2 intensities for each subunit. 
                                   The mean and standard deviation were then calculated, along with the confidence 
                                   interval based on the t-statistic for each complex 
                                   (calculated from the difference between the control and sample of interest from each subunit). 
                                   A paired t-test calculated significance between the control and sample of interest relative abundances for each complex." ),
                                               br(),
                                               br(),
                                               br()
                                        ))
                             ),
                             
                             tabPanel("RCA Data",
                                      fluidPage(
                                        br(),
                                        #h2(strong("Data for Selected Sample:"), style = "font-size:20px;"),
                                        mainPanel(
                                          dataTableOutput('patientdatatable'),
                                          linebreaks(7)
                                        )
                                      )
                             ),
                             tabPanel("Volcano Plot",
                                      br(),
                                      #h2(strong("Volcano Plot"), style = "font-size:20px;"),
                                      fluidRow(
                                        column(3,align = "center",offset = 4,
                                               br(),
                                               selectInput('inputFile_vol', label='select a file:',
                                                           choice=datalist),
                                               br()
                                               )),
                                      fluidRow(column(3,align = "center", offset = 4,
                                                      plotOutput('volcanoPlot',click='plot_click', brush='plot_brush', width = "100%"),
                                                      br()
                                                      )),
                                      fluidRow(column(11, align = "center",
                                                      tableOutput('clickedPoints'),
                                                      tableOutput('brushedPoints'),
                                                      linebreaks(7)
                                                      ))
                             ),
                             ### FILE UPLOAD COMPARISON ###    
                             tabPanel("Upload User Data",
                                      br(),
                                      h2(strong("Upload your own file"), style = "font-size:20px;"),
                                      sidebarLayout(
                                        sidebarPanel(
                                          fileInput('file1', 'Choose CSV File',
                                                    accept=c('text/csv', 
                                                             'text/comma-separated-values,text/plain', 
                                                             '.csv')),
                                          selectInput('inputFile_patientcomp', label='Select a sample for comparison correlation:',choice=datalist),
                                          tags$br(),
                                          downloadButton("examplefile", label = "Download an example file")
                                          
                                        ),
                                        mainPanel(
                                          tabsetPanel(
                                            tabPanel(
                                              title = "Correlation Plot",
                                              plotlyOutput(outputId = "corrplot"),
                                              selectizeInput(inputId = "Names", label = "Search for a Protein", choices = NULL, multiple = TRUE)
                                            ),
                                            tabPanel(
                                              title = "Data",
                                              dataTableOutput("contents"),
                                              linebreaks(7)
                                            ),
                                            
                                            tabPanel(
                                              title = "RCA",
                                              br(),
                                              downloadButton("SaveRCA", label = "Download Normalised Plot as PDF"),
                                              fluidRow(splitLayout(cellArgs = list(style = "padding: .7em 5.2%"),
                                                                   textOutput("mitofactor2"))),
                                              fluidRow(splitLayout(cellArgs = list(style = "padding: 30px"),cellWidths = 600,
                                                                   girafeOutput("rcaplot_unnorm", width = "100%", height = "100%" ),girafeOutput("rcaplot",width = "100%", height = "100%") )),
                                              fluidRow(splitLayout(cellArgs = list(style = "padding: 25px"),cellWidths = 700, tableOutput("up_rcadata_unnorm"), tableOutput("up_rcadata"))),
                                              br(),
                                              textOutput("signif_up"),
                                              linebreaks(7)
                                              #need to add a new girafeOutput("extra_complexes")
                                            ),
                                            tabPanel(
                                              title = "Relative Mitochondrial Factor",
                                              br(),
                                              br(),
                                              # fluidRow(splitLayout(cellArgs = list(style = "padding: .7em 15%"),
                                              # textOutput("mitofactor2"))),
                                              br(),
                                              fluidRow(splitLayout(cellWidths = 600,
                                                                   girafeOutput("mitofactorplot"))),
                                              linebreaks(7)
                                              
                                            )
                                          )
                                        )
                                      )
                             ),
                             tabPanel("Acknowledgments",
                                      fluidPage(
                                        br(),
                                        h2(strong("Acknowledgments"), style = "font-size:20px;"),
                                        br(),
                                        mainPanel(
                                          h3("Citation",style = "font-size:12px;"),
                                          br(),
                                          HTML('<div class="container-fluid">
                                 <label for="text"> Untargeted proteomics enables ultra-rapid variant prioritization in mitochondrial and other rare diseases. Hock et al. 2024. Current Medrxiv Citation at: https://shorturl.at/s85s7 </label>
                                 </div>'),
                                          #HTML('<span class="__dimensions_badge_embed__" data-doi="10.1001/jama.2016.9797"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>'),
                                          HTML('<hr>'),
                                          h3("Funding", style = "font-size:12px;"),
                                          HTML('<div class="container-fluid">
                                 <label for="text"> This research was supported by Australian National Health and 
                                 Medical Research Council (NHMRC) Project and Ideas grants (1140906 to DAS; 1164479 to DRT; 2010939 to MTR), Investigator Fellowships 
                                 (2009732 to DAS, 2010149 to LuEF and 1155244 to DRT) and a 
                                 Principal Research Fellowship (1155244 to DRT) along with funding by Australian 
                                 Genomics Health Alliance (Australian Genomics) NHMRC Targeted Call for Research grant GNT1113531. 
                                 Additional support came from the Australian Medical Research Future Fund Genomics Health Futures 
                                 Mission (2007959 to DRT) and the US Department of Defense Congressionally Directed Medical Research 
                                 Programs (PR170396 to DRT). We thank the Mito Foundation for the provision of instrumentation through 
                                 research equipment grants to DAS and DHH. Additionally, LuEF acknowledges support from the Mito Foundation. 
                                 Work at the MCRI is supported through the Victorian Government’s Operational Infrastructure Support Program. 
                                 RWT is funded by the Wellcome Centre for Mitochondrial Research (203105/Z/16/Z), the Mitochondrial Disease Patient 
                                 Cohort (UK) (G0800674), the Medical Research Council (MR/W019027/1), the Lily Foundation, Mito Foundation, the 
                                 Pathological Society, the UK NIHR Biomedical Research Centre for Ageing and Age-related disease award to 
                                 the Newcastle upon Tyne Foundation Hospitals NHS Trust, LifeArc and the UK NHS Highly Specialised Service 
                                 for Rare Mitochondrial Disorders of Adults and Children.
  </label>
                                 </div>'),
                                          HTML('<hr>'),
                                          h3("Ethics Statment", style = "font-size:12px;"),
                                          HTML('<div class="container-fluid">
                                 <label for="text"> This study was conducted in accordance with the revised Declaration of Helsinki 
                                 and following the Australian National Health and Medical Research Council statement of ethical conduct 
                                 in research involving humans. Samples were obtained after receiving written, informed consent for 
                                 diagnostic or research investigations from the respective responsible human ethics institutional
                                 review boards under HREC/89419/RCHM-2022, HREC/RCH/34228, HREC/16/MH/251, HREC/82160/RCHM-2022, HREC/RCH/34183 
                                 and REC ref 2002/205 by the Newcastle and North Tyneside Local Research Ethics Committee. </label>
                                 </div>'),
                                          HTML('<hr>'),
                                          h3("Credits", style = "font-size:12px;"),
                                          HTML('<div class="container-fluid">
                                 <label for="text"> Website Conceptualisation: David Stroud, Daniella Hock, Nikeisha Caruana <br>
                                                    Website Design and Programming: Nikeisha Caruana<br>
                                                    Beta Testing: Liana Semcesen, Tanavi Sharma, Daniella Hock<br>
                                                    Datasets: Obtained from all sources acknolwedged within the publication and dataset tables<br> 
                                                    RDMSExplorer is a live database regularly updated with new data. </label>
                                 </div>'),
                                          HTML('<hr>'),
                                          h3("Copyright", style = "font-size:12px;"),
                                          HTML('<div class="container-fluid">
                                 <label for="text"> RDMSExplorer content and code are published under the Creative Commons Attribution-NonCommercial 4.0 International </label>
                                 </div>'),
                                        )
                                      )
                             )
                           )
                         )
)#)

server <- shinyServer(function(input, output, session) {
  
  ### PASSWORD ###
  result_auth <- secure_server(check_credentials = check_credentials(credentials))
  
  output$res_auth <- renderPrint({
    reactiveValuesToList(result_auth)
  })
  
  ################# PATIENT DATA #################  
  ## Main Page Summary ##
  
  ## Data ##
  p_dataraw <- reactive({ 
    filename <- paste0("data/rca-data/",input$inputFile,".csv")
    read.table(file=filename, header=T, sep=',')
  })
  
  ## Normalisation with Mean ##
  
  mitodata <- reactive({
    mitodata <- p_dataraw()[c(1:6)]
    mitodata <- mitodata[which(mitodata$mito == "mito"), ]
  })
  
  normfactor_patient <- reactive({
    rca_patient <- mitodata()[c('uniprot','mean_patient')] %>% 
      gather(samples, values, -"uniprot") %>% 
      group_by(samples) %>% 
      dplyr::summarize(mitomean_pat = mean(values, na.rm = TRUE))
  })
  
  patient <- reactive({
    normfactor_patient()$mitomean_pat
  })
  
  meanpatient_data <- reactive({
    meanpatient_data <- p_dataraw()[c('uniprot','genename','mean_patient','complex')] %>%
      gather(samples, values, -"uniprot", -"complex", -"genename") %>% 
      group_by(uniprot, samples, complex, genename) %>% 
      mutate(norm = values/patient()) %>% 
      separate(samples, into = c("mean", "type")) 
  })
  
  normfactor_ctrl <- reactive({
    rca_ctrl <- mitodata()[c('uniprot','mean_ctrl')] %>% 
      gather(samples, values, -"uniprot") %>% 
      group_by(samples) %>% 
      dplyr::summarize(mitomean_ctrl = mean(values, na.rm = TRUE))
  })
  
  ctrl <- reactive({
    normfactor_ctrl()$mitomean_ctrl
  })
  
  meanctrl_data <- reactive({
    meanctrl_data <- p_dataraw()[c('uniprot','mean_ctrl','complex', 'genename')] %>%
      gather(samples, values, -"uniprot", -"complex", -"genename") %>% 
      group_by(uniprot, samples, complex, genename) %>% 
      mutate(norm = values/ctrl()) %>% 
      separate(samples, into = c("mean", "type")) 
  })
  
  mergedmean <- reactive({
    rbind(meanpatient_data(),meanctrl_data())
  })
  
  ## Normalised Stats Calculation ##
  
  normed_data<- reactive({ 
    merged_mean <- mergedmean()[-which(mergedmean()$complex == ""), ]
    merged_mean <- merged_mean %>% 
      dplyr::group_by(uniprot, type, complex, genename) %>% 
      dplyr::summarize(log10 = log10(norm)) %>% 
      group_by(uniprot) %>% 
      dplyr::mutate(Diff = log10 - log10[type == 'ctrl']) %>% 
      dplyr::group_by(uniprot) %>% 
      dplyr::mutate(Ratio = log10 / log10[type == 'ctrl']) %>% 
      dplyr::group_by(uniprot) %>% 
      dplyr::mutate(RCA = 10^Diff)
  })
  
  normed_data2 <- reactive({ 
    normed_data2 <- normed_data() %>% 
      filter(type != "ctrl") %>%
      filter(complex != "Ribo")
  })
  
  summarystats_norm <- reactive({ 
    stats_norm <- summarySE(normed_data2(), measurevar = "RCA", groupvars=c("type", "complex"))
  })
  
  output$p_RCAPlot <- renderGirafe({
    girafe(ggobj = ggplot(normed_data2(), aes(complex, RCA)) + 
             geom_jitter_interactive(aes(tooltip = genename),fill="lightgrey", colour = "grey35", pch=21, width = 0.2, size = 6) +
             geom_errorbar(data = summarystats_norm(), aes(ymin=RCA-sd, ymax=RCA+sd), size = 1, width=.05, colour = "black") +
             geom_hline(yintercept = 1, linetype = "dotted") + 
             stat_summary(fun=mean, geom="crossbar", color="black", size = .5, width = .3) +
             xlab("Complex") +
             ylab("Relative Complex Abundance") +
             scale_y_continuous(limits = c(0.2,4), labels = scales::percent) +
             theme_classic()
    )
  })
  
  
  
  # RCA Stats Table #
  output$p_rcadata <- renderTable({
    p_rcadata <- compare_means(log10 ~ type, data = normed_data(), group.by = "complex", method = "t.test", paired = TRUE)
    p_rcadata <- merge(p_rcadata, summarystats_norm(), by = "complex")
    p_rcadata <- p_rcadata %>% mutate(`RCA ` = sprintf("%1.0f%%", 100*RCA))
    p_rcadata <- p_rcadata[c(1,7,8,11,13,15,16)]
  })
  
  # Significance Text #
  output$signif <- renderText({
  })
  
  ## RCA Patient Data Table ##
  patientdata <- reactive({ 
    patientdata <- normed_data2()[c('genename','RCA')]
    patientdata <- merge(patientdata,p_dataraw(), by = c("genename"))
    patientdata <- patientdata[c("genename","uniprot","complex","mean_ctrl","mean_patient","RCA")]
  })
  
  output$patientdatatable <- renderDataTable({
    datatable(patientdata()) %>% formatRound(c("mean_ctrl", "mean_patient"), 0) %>% formatPercentage(c("RCA"), 0)
  })
  
  
  ## Stats Calculation Unnormalised Data ##
  p_data <- reactive({ 
    p_data <- p_dataraw()[c(1:4)]
    p_data <- p_data[-which(p_data$complex == ""), ]
    p_data <- p_data %>%
      gather(samples, values, -"genename", -"complex") %>% 
      dplyr::group_by(genename, samples, complex) %>% 
      dplyr::summarize(log10 = log10(values)) %>% 
      dplyr::group_by(genename) %>% 
      dplyr::mutate(Diff = log10 - log10[samples == 'mean_ctrl']) %>%
      dplyr::group_by(genename) %>% 
      dplyr::mutate(Ratio = log10 / log10[samples == 'mean_ctrl']) %>% 
      dplyr::group_by(genename) %>% 
      dplyr::mutate(RCA = 10^Diff)
  })
  
  p_data2 <- reactive({ 
    p_data2 <- p_data() %>% filter(samples != "mean_ctrl")
  })
  
  summarystats <- reactive({ 
    stats <- summarySE(p_data2(), measurevar = "RCA", groupvars=c("samples", "complex"))
  })
  
  
  # output$mitofactor <- renderText({
  #     "Relative Mitochondrial Factor:"
  # })
  
  
  ## Volcano Plot ##
  volcanodata <-  reactive({
    filename <- paste0("data/volcano-data/",input$inputFile_vol,".csv")
    read.table(file=filename, header=T, sep=',')
  })
  
  colScale <- reactive({
    myColors <- brewer.pal(5,"Set1")
    names(myColors) <- levels(volcanodata()$MitoList)
    colScale <- scale_colour_manual(name = "MitoList", values = myColors)
  })
  
  #plot it normally with ggplot:
  output$volcanoPlot <- renderPlot({ 
    ggplot(volcanodata(),aes(x=log2FC_patient, y=negLogPval, color=MitoList)) +
      geom_point() +
      coord_cartesian() +
      ylab("-log10 P value") +
      xlab("log2 fold change") +
      colScale() + 
      theme(legend.position = "bottom",
            text = element_text(size = 16),
            legend.text = element_text(size = 16)) 
  })
  
  brushed <- reactive({
    # We need to tell it what the x and y variables are:
    brushedPoints(volcanodata(), input$plot_brush, xvar = "log2FC_patient", yvar = "negLogPval")
  })
  
  output$brushedPoints <- renderTable({
    brushed()
  }, rownames = T)
  
  
  
  
  ######### DATA UPLOAD #########  
  
  #Patient data (used from Volcano plots - has all data not just mitos)
  data_patient <- reactive({
    filename_patientcomp <- paste0("data/volcano-data/",input$inputFile_patientcomp,".csv")
    read.table(file=filename_patientcomp, header=T, sep=',')
  })
  
  #read in custom mitocarta database. 
  mitocarta <- reactive({ 
    read.csv(file="data/MitoCarta3-Dani.csv")
  })
  
  #patient data upload
  data_upload <- reactive({ 
    req(input$file1)
    inFile <- input$file1 
    df <- read.csv(inFile$datapath)
    return(df)
  })
  #log data for further analysis 
  datalog <- reactive({
    df <- data_upload() %>% 
      mutate(log_P = log2(mean_patient)) %>% 
      mutate(log_C = log2(mean_ctrl)) %>% 
      mutate(log2FC = log_P-log_C) %>%
      select(-c(uniprot))
  })
  
  #Need to annotate complexes and Other, then duplicate complexes for analysis and graphing.
  #Also need to indicate which are mitos and which are not.
  #The issue is that there is no NonMitos in the complexsub_merge table. 
  
  #subset mitocarta
  complexsubs <- reactive({
    mitocarta()[c('uniprot', 'genename','MitoList','complex')]
  })
  
  #add in complexes [need to ID which are mitos in mitocarta, change + to Mito, and then - to NonMito ]
  complexsub_merge <- reactive({
    complexmerge <- merge(datalog(), complexsubs(), by = c("genename"))
  })
  
  
  
  #Data for outlier removal (just the single set ID'd by 'Other')
  # data_outliers <- reactive({
  #   data2 <- complexsub_merge() %>% filter(complex == "Other")
  #   rownames(data2) <- data2$uniprot
  #   return(data2)
  # })
  
  #Example File Download
  # The requested dataset
  exampledataoutput <- reactive({
    filename <- paste0("data/examplefile_RDMSExplorer.csv")
    read.table(file=filename, header=T, sep=',')
  })
  
  output$examplefile <- downloadHandler(
    filename = function(file) {
      "examplefile.csv"
    },
    
    content = function(file) {
      write.csv(exampledataoutput(), file,row.names = FALSE)
    }
  )
  
  ### OUTLIER TEST and warning for upload ###
  
  #Removal of Outliers for Mean calculation. Need to calculate before any display in table or RCA.
  #Used for both Norm and Unnorm data. But the other data needs two sets of complexes doesnt it? No it doesnt,
  #You only need the mean overall value for the control and the patient and then you are going to apply it to the
  #original data uploaded because youre only looking at complexes for the RCA (unless the issues are with the complexes?)
  
  # rmoutlier_data <- reactive({
  #   data2 <- data_outliers()[c('mean_ctrl','mean_patient')]
  #   outlierrm <- outlier_fun(data2)
  #   data_nooutlier <- merge(data_outliers(), outlierrm, by = "uniprot")
  #   return(data_nooutlier)
  # })
  # rmprots <- reactive({
  #   ogdata <- data_outliers()$genename
  #   newdata <- rmoutlier_data()$genename
  #   setdiff(ogdata,newdata)
  # })
  # rmprots2 <- reactive({
  #   toString(rmprots())
  # })
  # observeEvent(input$file1, {
  #   showModal(modalDialog(
  #     title = "Outliers Detected",
  #     print(paste0("Warning: Mitochondrial proteins removed for  mean calculation for mitochondrial normalisation
  #                  due to high variation. These proteins are: ", rmprots2(), ". These proteins are only removed for the mean
  #                  calculation. We recommend checking these outliers in your original data to ensure they make sense.")),
  #     easyClose = TRUE
  #   ))
  #})
  
  ### CALCULATION OF MEANS FOR PATIENT AND CONTROL FOR NORMALISATION ###
  
  ##change this to rmoutlier_data for outlier calculation in 2.0
  
  cleandata <- reactive({
    complexsub_merge()
  })
  
  #data without nonmitos
  data_upload2 <- reactive({
    data_upload2 <- complexsub_merge() %>% filter(MitoList != "NonMito")
  })
  
  
  #calculate mito mean for patients
  patient_mean <- reactive({
    mean_pat <- data_upload2()[c('uniprot', 'mean_patient')] %>%
      gather(samples, values, -"uniprot") %>%
      group_by(samples) %>%
      dplyr::summarize(mean_size = mean(values, na.rm = TRUE))
    mitos_pat_mean <- mean_pat %>%
      group_by(samples) %>%
      dplyr::summarize(mitomean_pat = mean(mean_size, na.rm = TRUE))
    patientmean <- mitos_pat_mean$mitomean_pat
    return(patientmean)
  })

  #apply mean to all patient data
  patientmean_norm <- reactive({
    p_mitomean <- data_upload2()[c('uniprot','genename','mean_patient','complex','MitoList')]
    p_mitomean<- p_mitomean %>%
      gather(samples, values, -"genename", -"uniprot", -"complex", -"MitoList") %>%
      group_by(uniprot, genename, samples, complex, MitoList) %>%
      mutate(norm = values/patient_mean()) %>%
      separate(samples, into = c("mean", "type"))
  })
  
  #calculate mito mean for controls
  control_mean <- reactive({
    mean_ctrl <- data_upload2()[c('uniprot', 'mean_ctrl')] %>%
      gather(samples, values, -"uniprot") %>%
      group_by(samples) %>%
      dplyr::summarize(mean_size = mean(values, na.rm = TRUE))
    mitos_ctrl_mean <- mean_ctrl %>%
      group_by(samples) %>%
      dplyr::summarize(mitomean_ctrl = mean(mean_size, na.rm = TRUE))
    ctrlmean <- mitos_ctrl_mean$mitomean_ctrl
    return(ctrlmean)
  })
  
  #apply mean to all control data
  ctrlmean_norm <- reactive({
    c_mitomean <- data_upload2()[c('genename','uniprot','mean_ctrl','complex','MitoList')]
    c_mitomean<- c_mitomean %>%
      gather(samples, values,-"genename", -"uniprot", -"complex", -"MitoList") %>%
      group_by(uniprot, genename, samples, complex, MitoList) %>%
      mutate(norm = values/control_mean()) %>%
      separate(samples, into = c("mean", "type"))
  })
  
  #merge datasets
  mergedmean_normed <- reactive({
    rbind(ctrlmean_norm(),patientmean_norm())
  })
  
  
  # genename_symbols <- reactive({
  #   symbols <- data_upload2() %>%
  #     filter(complex == "Other") %>% select(genename, uniprot)
  # })
  
  genename_symbols <- reactive({
    symbols <- data_upload2() %>%
      select(genename, uniprot)
  })
  
  ### NORMALISED DATA STATS AND  RCA GRAPH ###
  
  
  # #Run the below for the stats and then output the graph. Need to log first
  up_data_normed <- reactive({
    up_data_normed <- mergedmean_normed() %>%
      dplyr::group_by(uniprot, genename, type, complex, MitoList) %>%
      dplyr::summarize(log10 = log10(norm)) %>%
      dplyr::group_by(uniprot, genename) %>%
      dplyr::mutate(Diff = log10 - log10[type == 'ctrl']) %>%
      dplyr::group_by(uniprot, genename) %>%
      dplyr::mutate(Ratio = log10 / log10[type == 'ctrl']) %>%
      dplyr::group_by(uniprot, genename) %>%
      dplyr::mutate(relativeabund = 10^Diff)
  })
  
  
  # symbolmerge <- reactive({
  #   merge(up_data_normed(), genename_symbols(), by = "uniprot")
  # })
  
  
  #maybe this is the section that can be changed to be dynamic? Click boxes turn on/off subunits. Keep CI-CV as standard output, 
  #but then have the option to turn all check boxes on and off?
  #Do we even want the unnormalised graph 
  
  up_data_normed2 <- reactive({
    up_data_normed2 <- up_data_normed() %>%
      filter(type != "ctrl") %>%
      filter(complex != "Other") 
  })
  
  summarystats_up_normed <- reactive({
    stats_up_normed <- summarySE(up_data_normed2(), measurevar = "relativeabund", groupvars=c("type", "complex"))
  })
  
  output$rcaplot <- renderGirafe({
    girafe(ggobj = ggplot(up_data_normed2(), aes(complex, relativeabund)) +
             geom_jitter_interactive(aes(tooltip = genename),fill="lightgrey", colour = "gray35", pch=21, width = 0.2, size = 6) +
             geom_errorbar(data = summarystats_up_normed(), aes(ymin=relativeabund-sd, ymax=relativeabund+sd), size = 1, width=.1, colour = "black") +
             geom_hline(yintercept = 1, linetype = "dotted") +
             stat_summary(fun=mean, geom="crossbar", color="black", size = .8, width = .3) +
             ggtitle("Normalised") +
             xlab("Complex") +
             ylab("Relative Complex Abundance") +
             scale_y_continuous(labels = scales::percent) +
             theme_classic()
    )
  })
  
  output$up_rcadata <- renderTable({
    up_rcadata <- compare_means(log10 ~ type, data = up_data_normed(), group.by = "complex", method = "t.test", paired = TRUE)
    up_rcadata <- merge(up_rcadata, summarystats_up_normed(), by = "complex")
    up_rcadata <- up_rcadata %>% mutate(relativeabund = sprintf("%1.0f%%", 100*relativeabund))
    up_rcadata <- up_rcadata[c(1,7,8,11,13,15,12)]
  })
  
  pdfplot <- reactive({
    ggplot(up_data_normed2(), aes(complex, relativeabund)) +
      geom_jitter_interactive(aes(tooltip = genename),fill="lightgrey", colour = "gray35", pch=21, width = 0.2, size = 4) +
      geom_errorbar(data = summarystats_up_normed(), aes(ymin=relativeabund-sd, ymax=relativeabund+sd), size = 0.7, width=.2, colour = "black") +
      geom_hline(yintercept = 1, linetype = "dotted") +
      stat_summary(fun=mean, geom="crossbar", color="black", size = .4, width = .5) +
      theme_classic()+
      ggtitle("Normalised") +
      xlab("Complex") +
      ylab("Relative Complex Abundance") +
      scale_y_continuous(labels = scales::percent)
  })
  
  #output plot
  output$SaveRCA <- downloadHandler(
    filename = function(file) {
      "RCAplot.pdf"
    },
    content = function(file) {
      ggsave(file, plot = pdfplot(), width = 150, height = 150, units = "mm", device = "pdf")
    }
  )
  
  
  # 
  # 
  ### STATS and DATA UNNORMED DATA ###
  
  up_data <- reactive({
    up_data <- complexsub_merge()[c('mean_ctrl','mean_patient','complex','genename', 'MitoList')]
    up_data <- up_data[-which(up_data$complex == "Other" | up_data$MitoList == "NonMito"), ]
    up_data <- up_data %>%
      gather(samples, values, -"genename", -"complex", -"MitoList") %>%
      dplyr::group_by(genename, samples, complex, MitoList) %>%
      dplyr::summarize(log10 = log10(values)) %>%
      dplyr::group_by(genename) %>%
      dplyr::mutate(Diff = log10 - log10[samples == 'mean_ctrl']) %>%
      dplyr::group_by(genename) %>%
      dplyr::mutate(Ratio = log10 / log10[samples == 'mean_ctrl']) %>%
      dplyr::group_by(genename) %>%
      dplyr::mutate(relativeabund = 10^Diff)
  })
  up_data2 <- reactive({
    up_data2 <- up_data() %>% 
      filter(samples != "mean_ctrl") %>%
      filter(complex != "Ribo")
  })
  summarystats_up <- reactive({
    stats_up <- summarySE(up_data2(), measurevar = "relativeabund", groupvars=c("samples", "complex"))
  })
  
  
  
  ### RCA PLOT ###
  output$rcaplot_unnorm <- renderGirafe({
    girafe(ggobj = ggplot(up_data2(), aes(complex, relativeabund)) +
             geom_jitter_interactive(aes(tooltip = genename),fill="lightgrey", colour = "black", pch=21, width = 0.2, size = 6) +
             geom_errorbar(data = summarystats_up(), aes(ymin=relativeabund-ci, ymax=relativeabund+ci), size = 1, width=.1, colour = "black") +
             geom_hline(yintercept = 1, linetype = "dotted") +
             stat_summary(fun=mean, geom="crossbar", color="black", size = .8, width = .3) +
             ggtitle("Unnormalised") +
             xlab("Complex") +
             ylab("Relative Complex Abundance") +
             scale_y_continuous(labels = scales::percent) +
             theme_classic()
    )
  })
  # 
  # DATA UNDER RCA PLOTS #
  
  rcadata_unnorm <- reactive({
    up_rcadata <- compare_means(log10 ~ samples, data = up_data(), group.by = "complex", method = "t.test", paired = TRUE)
    up_rcadata <- merge(up_rcadata, summarystats_up(), by = "complex")
    up_rcadata <- up_rcadata %>% mutate(relativeabund = sprintf("%1.0f%%", 100*relativeabund))
    up_rcadata <- up_rcadata[c(1,7,8,11,13,15,12)]
  })
  
  output$up_rcadata_unnorm <- renderTable({
    rcadata_unnorm()
  })
  
  output$signif_up <- renderText({
    "ns: p > 0.05, *: p <= 0.05, **: p <= 0.01, ***: p <= 0.001, ****: p <= 0.0001"
  })
  
  
  
  
  
  ### Extra graph of MITO FACTOR ###
  
  mitofactor_unnorm <- reactive({
    unnormfactor_data <- up_data() %>% 
      filter(samples != "mean_ctrl") %>%
      filter(MitoList == "Mito") %>%
      mutate(MitoList = str_replace(MitoList, "Mito", "Unnorm"))
  })
  
  mitofactor_norm <- reactive({
    normfactor_data <- up_data_normed() %>%
      filter(type != "ctrl") %>%
      filter(complex != "Other") %>%
      filter(MitoList == "Mito") %>%
      mutate(MitoList = str_replace(MitoList, "Mito", "Norm"))
  })
  
  mitofactortotal <- reactive({
    mitofactortotal <- rbind(mitofactor_unnorm(), mitofactor_norm())
    mitofactortotal <- mitofactortotal %>% mutate(MitoList = factor(MitoList, levels = c("Unnorm", "Norm"), ordered = TRUE))
  })
  
  summarystats_mitototal <- reactive({
    stats_up <- summarySE(mitofactortotal(), measurevar = "relativeabund", groupvars=c("samples", "MitoList"))
  })
  
  output$mitofactorplot <- renderGirafe({
    girafe(ggobj = ggplot(mitofactortotal(), aes(MitoList, relativeabund)) +
             geom_jitter_interactive(aes(tooltip = genename),fill="lightgrey", colour = "gray35", pch=21, width = 0.1, size = 6) +
             geom_errorbar(data = summarystats_mitototal(), aes(ymin=relativeabund-ci, ymax=relativeabund+ci), size = 1, width=.1, colour = "black") +
             geom_hline(yintercept = 1, linetype = "dotted") +
             stat_summary(fun=mean, geom="crossbar", color="black", size = .8, width = .3) +
             ylab("Relative Complex Abundance") +
             scale_y_continuous(labels = scales::percent) +
             coord_fixed(1.75) +
             theme_classic()
    )
  })
  
  # Mito Factor Calculation #
  
  mitofactor <- reactive({
    mitofactorcalc <- summarystats_mitototal() %>% mutate(mitofactor = sum(relativeabund[MitoList == 'Norm']) - sum(relativeabund[MitoList == 'Unnorm']))
    mitofactorcalc <- mitofactorcalc %>% mutate(mitofactor = sprintf("%1.0f%%", 100*mitofactor))
    mitofactorcalc <- mitofactorcalc$mitofactor[[1]]
  })
  
  output$mitofactor2 <- renderText({
    print(paste0("Relative Mitochondrial Factor: ", mitofactor()))
  })

  
  
  # 
  #### DATA TABLE TAB ####
  
  mergeddata <- reactive({
    mdata <- merge(complexsub_merge(), data_patient(), by = c("genename", "uniprot", "complex"))
    mdata <- mdata %>% arrange(match(complex, c( "CI", "CII", "CIII","CIV","CV","mtSSU","mtLSU","Other")))
  })
  
  output$contents <- renderDataTable({
    dtable <- mergeddata()[c("genename","uniprot","complex","log2FC","log2FC_patient")]
  })

  
  #### CORRELATION PLOT ####
  
  
  observe({
    updateSelectizeInput(session, "Names", choices = mergeddata()$genename)
  })
  
  
  ##function that creates plotly dictionary object for single annotation
  listify_func <- function(x, y, text){
    return(list(
      x = x,
      y = y,
      text = as.character(text),
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 7,
      arrowcolor = "#ff9900",
      font = list(color = "#000000", size = 10),
      ax = runif(1, 1, 90),
      ay = -runif(1, 1, 90),
      bgcolor = "#f5f5f5",
      bordercolor = "#b3b3b3"
    ))
  }
  
  ##creating reactive object containing the subsetted data:
  highlighted <- eventReactive(input$Names,{
    data_sub <- mergeddata() %>% dplyr::filter(mergeddata()$genename %in% input$Names)
  })
  
  ##creating a recursive list of annotation lists for selected genes
  annotation_list <- reactiveValues(data = NULL)
  observeEvent(input$Names,{
    x <- highlighted()$log2FC
    y <- highlighted()$log2FC_patient
    text <- highlighted()$genename
    ##create dataframe with relative x, y and text values to create
    ##annotations:
    df <- data.frame(x = x, y = y, text = text)
    ##create matrix of lists defining each annotation
    ma <- mapply(listify_func, df$x, df$y, df$text)
    if(length(ma) > 0){
      ##convert matrix to list of lists:
      annotation_list$data <- lapply(seq_len(ncol(ma)), function(i) ma[,i])
    }
  })
  ##if nothing is selected, clear recursive list (i.e. remove annotations):
  observe({
    if(is.null(input$Names)){
      annotation_list$data <- list()
    }
  })
  
  
  output$corrplot <- renderPlotly({
    p <- plot_ly(mergeddata(),
                 x = ~log2FC,
                 y = ~log2FC_patient,
                 type = "scatter",
                 mode = "markers",
                 color = ~factor(complex, levels = c("NonMito", "Other","CI", "CII", "CIII","CIV","CV")),
                 hoverinfo = "text",
                 hovertext = ~paste("Name: ",genename)) %>%
      layout(yaxis = list( title= ~paste(input$inputFile_patientcomp,"log2FC")),
             xaxis = list( title = "uploaded_log2FC"))
  })
  
  
  ##proxy updating recursive list for annotations
  observeEvent(annotation_list$data,{
    plotlyProxy("corrplot", session) %>%
      plotlyProxyInvoke(
        "relayout",
        list(annotations = annotation_list$data)
      )
  })
  
  
})

shinyApp(ui, server)