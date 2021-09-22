#---------------------------------------------------------------------------------------
# Archivo: geneticae_shiny_app
# Descripcion: este archivo genera una aplicación Shiny que permite analizar datos
#              provenientes de ensayos multiambientales. Análisis descriptivos, ANOVA,
#              y biplots se pueden realizar. Estos últimos se obtienen del paquete
#              geneticae.
# Autora: Julia Angelini, Marcos Prunello, Gerardo Cervigni
#---------------------------------------------------------------------------------------

#-----------
# Paquetes
#-----------

rm(list=ls(all=TRUE))
library(shiny)
library(ggplot2)
library(car)
library(shinythemes)
library(DT)
library(geneticae)
library(dplyr)
library(kableExtra)
library(ggcorrplot)
library(shinyWidgets)
library(ggthemes)
library(plotly)
library(shinyhelper)
library(tidyr)

#---------------------
# Funciones auxiliares
#---------------------

# Funcion que permite importar un conjunto de datos
source("csvFile.R")

# Funcion para realizar un gráfico de probabilidad normal
gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))

  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }

  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE


  p <- ggplot(df, aes(x=z, y=ord.x)) +
    # geom_point() +
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) +
    theme_few() + xlab("theorical quantile") + ylab("sample quantile")


  print(p)
}


#-------------
# Define UI
#-------------

ui <- navbarPage(
  title=HTML('<span style="font-family: Lobste, cursive;font-size:175%;color:white;font-weight:bold;">Geneticae APP</span></a>'),
  theme = shinytheme("united"),
  # Data
  tabPanel("Data",
           tabsetPanel(
                # Import a dataset
                 tabPanel("Upload data", icon = icon("table"),
                          sidebarPanel(
                                csvFileInput("datafile", "Upload data (.csv format)"), br(),
                                tags$p("Select the column name that contains:"),
                                uiOutput("select_gen"),
                                uiOutput("select_env"),
                                uiOutput("select_rep"),
                                uiOutput("select_pheno"),
                                width = 3
                                ),
                          mainPanel(dataTableOutput("table"), width=8)
                        ),

                 # Examples dataset
                 tabPanel(strong("Example datasets"), icon = icon("table"),
                          sidebarPanel(br(),
                                prettyRadioButtons(
                                   inputId = "exampledataset",
                                   label = "Example data",
                                   choices = c("without repetitions", "with repetitions"),
                                   icon = icon("check"),
                                   status = "info",
                                   animation = "rotate"
                                ),
                                actionButton("show_exampledataset", "Show"), br(), br(),
                                textInput("Filename1", "Choose file name", value = "Example dataset"),
                                downloadButton("downloadexample", "Download"),
                                width = 3
                                ),
                          mainPanel(htmlOutput("show_example"), width=8)
                          )
             )
  ),
  # Descriptive Analysis
  tabPanel("Descriptive analysis",
           tabsetPanel(
                # Boxplot
                tabPanel("Boxplot", icon = icon("bar-chart-o"),
                         sidebarLayout(
                              sidebarPanel(
                                  prettyRadioButtons(
                                      inputId = "Var",
                                      label = "Graph for...",
                                      choices = c("Environment", "Genotype"),
                                      icon = icon("check"),
                                      status = "info",
                                      animation = "rotate"
                                      )%>%
                                    helper(type = "inline",
                                           title = "Inline Help",
                                           content = c("Choose if you want to compare the phenotypic trait between genotypes or between environments."),
                                           size = "m"),
                                  pickerInput(
                                      inputId = "fillcol",
                                      label = "Fill color",
                                      choices = c("orange","dimgrey", "red","lightpink","white")
                                  ),
                                  textInput("axisx", "X-axis label:", value = "Environment"),
                                  textInput("axisy", "Y-axis label:", value = "Yield"),

                                  actionButton("do", "Run"),
                                  br(), br(),

                                  textInput("Filename_box", "File name", value = "Boxplot")%>%
                                    helper(type = "inline",
                                           title = "Inline Help",
                                           content = c("An <b>.html</b> file will be downloaded."),
                                           size = "m"),
                                  downloadButton("download_box", "Download"),
                                  width=3
                                  ),
                        mainPanel(br(),plotlyOutput("boxplot"), width=6)
                        )
                  ),
                # Correlation Plot
                tabPanel("Correlation plot", icon = icon("bar-chart-o"),
                          sidebarLayout(
                              sidebarPanel(
                                  prettyRadioButtons(
                                      inputId = "Var2",
                                      label = "Graph for...",
                                      choices = c("Environment", "Genotype"),
                                      icon = icon("check"),
                                      status = "info",
                                      animation = "rotate"
                                      )%>%
                                    helper(type = "inline",
                                           title = "Inline Help",
                                           content = c("Choose if you want to plot the phenotypic trait correlations between genotypes or environments."),
                                           size = "m"),
                                  prettyRadioButtons(
                                      inputId = "corrType",
                                      label = "Correlation index",
                                      choices = c("Spearman", "Pearson"),
                                      icon = icon("check"),
                                      status = "info",
                                      animation = "rotate"
                                      )%>%
                                    helper(type = "inline",
                                           title = "Inline Help",
                                           content = c("Choose a correlation coefficient to be computed."),
                                           size = "m"),
                                  actionButton("do_corrplot", "Run"), br(), br(),
                              br(),
                              textInput("Filename3", "File name", value = "Correlation plot")%>%
                                helper(#icon = "exclamation",
                                       #colour = "red",
                                       type = "inline",
                                       title = "Inline Help",
                                       content = c("A <b>.png</b> file will be downloaded."),
                                       size = "s"),
                              downloadButton("downloadcorplot", "Download"),
                              width=3
                              ),
                            mainPanel(br(), plotOutput("corplot"), width=6))
             ),
             # Correlation Matrix
             tabPanel("Correlation matrix", icon = icon("table"),
                        sidebarLayout(
                            sidebarPanel(
                                prettyRadioButtons(
                                    inputId = "Var3",
                                    label = "Calcullate for...",
                                    choices = c("Environment", "Genotype"),
                                    icon = icon("check"),
                                    status = "info",
                                    animation = "rotate"
                                    )%>%
                                  helper(type = "inline",
                                         title = "Inline Help",
                                         content = c("Choose if you want to estimate the phenotypic trait correlations between genotypes or environments."),
                                         size = "m"),
                                prettyRadioButtons(
                                    inputId = "correlType",
                                    label = "Correlation index",
                                    choices = c("Spearman","Pearson"),
                                    icon = icon("check"),
                                    status = "info",
                                    animation = "rotate"
                                    )%>%
                                  helper(type = "inline",
                                         title = "Inline Help",
                                         content = c("Choose a correlation coefficient to be computed."),
                                         size = "m"),
                                actionButton("do_corrmat", "Run"), br(), br(),
                                width=3),
                            mainPanel(br(), htmlOutput("cormat"))
                            )
             ),
             # Interaction Plot
             tabPanel("Interaction Plot", icon = icon("bar-chart-o"),
                          sidebarLayout(
                              sidebarPanel(
                                  prettyRadioButtons(
                                    inputId = "Var4",
                                    label = "Graph for...",
                                    choices = c("Environment", "Genotype"),
                                    icon = icon("check"),
                                    status = "info",
                                    animation = "rotate"
                                    )%>%
                                    helper(type = "inline",
                                           title = "Inline Help",
                                           content = c("Choose if you want to plot the phenotypic trait between genotypes or environments."),
                                           size = "m"),
                              textInput("axisx_int", "X-axis label:", value = "Environment"),
                              textInput("axisy_int", "Y-axis label:", value = "Yield"),
                              actionButton("do_int", "Run"), br(), br(),
                              textInput("Filename_int", "File name", value = "Interaction plot")%>%
                                helper(type = "inline",
                                       title = "Inline Help",
                                       content = c("An <b>.html</b> file will be downloaded."),
                                       size = "m"),
                              downloadButton("download_int", "Download"),
                              width=3),
                          mainPanel(br(), plotlyOutput("int"), width=6)
                          )
             )
       )
 ),
  # GGE biplot
  tabPanel("GGE Biplot", icon = icon("bar-chart-o"),
           sidebarLayout(
                sidebarPanel(
                  prettyRadioButtons(
                    inputId = "SVP",
                    label = "SVP type",
                    choices = c("Symmetrical","Genotype-focused", "Environment-focused"),
                    icon = icon("check"),
                    status = "info",
                    animation = "rotate"
                  )%>%
                    helper(type = "inline",
                           title = "Inline Help",
                           content = c("The SVD method must be selected, however, this choice does not alter the relationships or relative interactions between genotypes and environments, although the appearance of the biplot will be different.
                                       Symmetrical option allows comparison for both genotypes and environments; Genotype-Focused displays the interrelationship among genotypes more accurately than any other method does, and
                                       Environment-focused is most informative of interrelationships among environments."),
                           size = "m"),
                  prettyRadioButtons(
                         inputId = "plotType",
                         label = "Plot type",
                         choices = c("Biplot",
                                     "Selected Environment",
                                     "Selected Genotype",
                                     "Comparison of Genotype",
                                     "Which Won Where/What",
                                     "Mean vs. Stability",
                                     "Ranking Genotypes",
                                     "Relationship Among Environments",
                                     "Ranking Environments"),
                         icon = icon("check"),
                         status = "info",
                         animation = "rotate"
                     )%>%
                    helper(type = "inline",
                           title = "Inline Help",
                           content = c("Several GGE biplots views can be obtained: basic biplot, Selected Environment to identify the most suitable cultivars for a particular environment of interest,
                                       Selected Genotype to determine which is the most suitable environment for a genotype; Comparison of Genotype to compare two cultivars;
                                       Which Won Where/What allow the identification of the best cultivar in each environment or mega-environment;
                                       Mean vs. Stability to visualize mean yield and stability of genotypes in yield units per se;
                                       Ranking Genotypes compares the cultivars to the “ideal” one with the highest yield and absolute stability;
                                       Relationship Among Environments to understand the interrelationships between the environments and
                                       Ranking Environments to classify the environments with respect to the ideal one."),
                           size = "m"),
                  # Esto solo se muestra si distr_media se elige en Normal:
                  conditionalPanel(condition = 'input.plotType == "Selected Environment"',
                                   pickerInput('SelectedE','Environment:',
                                               options = list(`actions-box` = TRUE, size = 10), multiple = FALSE,
                                               choices = NULL)%>%
                                     helper(type = "inline",
                                            title = "Inline Help",
                                            content = c("Select the environment of interest"),
                                            size = "m")
                                   ),
                  conditionalPanel(condition = 'input.plotType == "Mean vs. Stability" ||
                                                input.plotType == "Ranking Genotypes"  ||
                                                input.plotType == "Relationship Among Environments" ||
                                                input.plotType == "Ranking Environments"',
                                   pickerInput('ME','Environments inside the mega-environment:',
                                               options = list(`actions-box` = TRUE, size = 10), multiple = TRUE,
                                               choices = NULL)%>%
                                     helper(type = "inline",
                                            title = "Inline Help",
                                            content = c("Select the environments inside the mega-environment to analyze"),
                                            size = "m")
                                  ),
                  conditionalPanel(condition = 'input.plotType == "Selected Genotype"',
                                   pickerInput('SelectedG','Genotype:',
                                               options = list(`actions-box` = TRUE, size = 10), multiple = FALSE,
                                               choices = NULL)%>%
                                     helper(type = "inline",
                                            title = "Inline Help",
                                            content = c("Select the genotype of interest"),
                                            size = "m")
                                  ),

                  conditionalPanel(condition = 'input.plotType == "Comparison of Genotype"',
                                   pickerInput('SelectedG1','Genotype 1:',
                                               options = list(`actions-box` = TRUE, size = 10), multiple = FALSE,
                                               choices = NULL)%>%
                                     helper(type = "inline",
                                            title = "Inline Help",
                                            content = c("Select the genotypes to be compared"),
                                            size = "m"),
                                   pickerInput('SelectedG2','Genotype 2',
                                               options = list(`actions-box` = TRUE, size = 10), multiple = FALSE,
                                               choices = NULL)
                                   ),
                     hr(),
                     materialSwitch(
                       inputId = "footnote",
                       label = "Footnote",
                       value = TRUE,
                       status = "primary"
                     ),
                     materialSwitch(
                       inputId = "title",
                       label = "Title",
                       value = TRUE,
                       status = "primary"
                     ),
                    materialSwitch(
                      inputId = "axislabels",
                      label = "Axislabels",
                      value = TRUE,
                      status = "primary"
                    ),
                    materialSwitch(
                    inputId = "axes",
                    label = "Axes",
                    value = TRUE,
                    status = "primary"
                    ),
                     hr(),
                     pickerInput(
                       inputId = "colgen",
                       label = "Genotype color",
                       choices = c("dimgrey", "red", "black")
                     ),
                     pickerInput(
                       inputId = "colenv",
                       label = "Environment color",
                       choices = c("dimgrey", "red", "black")
                     ),
                     pickerInput(
                       inputId = "colsegment",
                       label = "Segment color",
                       choices = c("dimgrey", "red", "black")
                      ),
                    sliderInput("sizeGen",
                              "Genotype marker size:",
                              min = 0,  max = 10,  value = 4),
                    sliderInput("sizeEnv",
                              "Environment marker size:",
                              min = 0,  max = 10,  value = 4),
                     actionButton("do_GGE", "Run"), br(), br(),
                    br(), br(),
                    textInput("Filename_GGE", "File name", value = "GGE Biplot"),
                    downloadButton("download_gge", "Download"),
                    width = 3
              ),
              mainPanel(plotOutput("plotGGE"))
      )
  ),

  #AMMI biplot
  tabPanel("AMMI Biplot", icon = icon("bar-chart-o"),
           sidebarLayout(
                sidebarPanel(
                   prettyRadioButtons(
                       inputId = "robustType",
                       label = "Plot type",
                       choices = c("AMMI","rAMMI", "hAMMI", "gAMMI", "lAMMI", "ppAMMI"),
                       icon = icon("check"),
                       status = "info",
                       animation = "rotate"),
                 hr(),
                   materialSwitch(
                       inputId = "foot",
                       label = "Footnote",
                       value = TRUE,
                       status = "primary"
                   ),
                   materialSwitch(
                       inputId = "tit",
                       label = "Title",
                       value = TRUE,
                       status = "primary"
                   ),
                 materialSwitch(
                   inputId = "axislabels_AMMI",
                   label = "Axislabels",
                   value = TRUE,
                   status = "primary"
                 ),
                 materialSwitch(
                   inputId = "axes_AMMI",
                   label = "Axes",
                   value = TRUE,
                   status = "primary"
                 ),
               hr(),
                   pickerInput(
                       inputId = "colorgen",
                       label = "Genotype color",
                       choices = c("dimgrey", "red", "black")
                   ),
                   pickerInput(
                       inputId = "colorenv",
                       label = "Environment color",
                       choices = c("dimgrey", "red", "black")
                   ),
               sliderInput("sizeGen_AMMI",
                           "Genotype marker size:",
                           min = 0,  max = 10,  value = 4),
               sliderInput("sizeEnv_AMMI",
                           "Environment marker size:",
                           min = 0,  max = 10,  value = 4),
               actionButton("do_AMMI", "Run"), br(), br(),
               textInput("Filename_AMMI", "File name", value = "AMMI Biplot"),
               downloadButton("download_ammi", "Download"),
               width = 3
              ),
            mainPanel( plotOutput("plotAMMI"))
        )
  ),
  # Help
  navbarMenu("Help",icon = icon("question-circle"),
             tabPanel("Getting Started",
                      mainPanel(
                        tags$h2(strong("Getting Started")),
                        br(),
                        tags$h4(strong("Motivation")),
                        tags$div(
                            tags$p("Understanding the relationship between crops performance and environment is a
                                    key problem for plant breeders and geneticists. In advanced stages of breeding
                                    programs, in which few genotypes are evaluated, multi-environment trials (MET)
                                    are one of the most used experiments. Such studies test a number of genotypes in
                                    multiple environments in order to identify the superior genotypes according to
                                    their performance. In these experiments, crop performance is modeled as a
                                    function of genotype (G), environment (E) and genotype-environment interaction
                                    (GEI). The presence of GEI generates differential genotypic responses in the
                                    different environments (Angelini et al., 2019; Crossa, 1990; Kang and Magari,
                                    1996). Therefore appropriate statistical methods should be used to obtain an
                                    adequate GEI analysis, which is essential for plant breeders (Giauffret et al.,
                                    2000).",
                                   align = "justify"),
                            tags$p("The average performance of genotypes through different environments can only be
                                    considered in the absence of GEI (Yan and Kang, 2003). However, GEI is almost
                                    always present and the comparison of the mean performance between genotypes is
                                    not enough. The most widely used methods to analyze MET data are based on
                                    regression models, analysis of variance (ANOVA) and multivariate techniques. In
                                    particular, two statistical models are widely used among plant breeders as they
                                    provide useful graphical tools for the study of GEI: the Additive Main effects
                                    and Multiplicative Interaction model (AMMI) (Kempton, 1984; Gauch, 1988) and the
                                    Site Regression Model (SREG) (Cornelius et al., 1996; Gauch and Zobel, 1997).
                                    However, these models are not always efficient enough to analyze MET data
                                    structure of plant breeding programs. They present serious limitations in the
                                    presence of atypical observations and missing values, which occur very
                                    frequently. To overcome this, several imputation alternatives and a robust AMMI
                                    were recently proposed in literature.",
                                   align = "justify"),
                            tags$p("Although there are R packages which tackle different aspects of MET data
                                    analysis, there aren't any packages capable of performing all the steps that
                                    need to be considered. The geneticae package was created to gather in one
                                    place the most useful functions for this type of analysis and it also implements
                                    new methodology which can be found in recent literature. More importantly,
                                    geneticae is the first package to implement the robust AMMI model and new
                                    imputation methods not available before. In addition, there is no need to
                                    preprocess the data to use the `geneticae` package, as it the case of some
                                    previous packages which require a data frame or matrix containing genotype by
                                    environment  means  with  the genotypes in rows and the environments in columns.
                                    In this package, data in long format is required. Genotypes, environments,
                                    repetitions (if any) and phenotypic traits of interest can be presented in any
                                    order and there is no restriction on columns names. Also, extra information that
                                    will not be used in the analysis may be present in the dataset. Finally,
                                    `geneticae` offers a wide variety of options to customize the biplots, which are
                                    part of the graphical output of these methods.",
                                    align = "justify"),
                            tags$p("The goal of the Geneticae Shiny Web APP is to provide a graphical user interface
                                   for the R package, so that it can be used by breeders and analysts with no previous
                                   experience in R programming. It is an interactive, noncommercial and open source
                                   software, offering a free alternative to available commercial software to analize METs.",
                                   align = "justify")
                        ),
                        tags$br(),
                        tags$h4(strong("Small example"),
                                align = "justify"),
                        tags$p("If you are just getting started with Geneticae APP we recommend visiting and exploring
                                  the examples throughout the tutorial. Here we present a small example.",
                                align = "justify"),
                        tags$p("The dataset yan.winterwheat available as example dataset has information about the yield
                                of 18 winter wheat varieties grown in nine environments in Ontario at 1993. The GGE biplot
                                visually addresses many issues relative to genotype and test environment evaluation. The
                                ", em (" GGE Biplot "), "tab allows to builds several GGE biplots views. The basic one is produced
                                by default. If Which Won Where/What is indicate in plot type, and the other options are
                                left by default the polygonal view of the GGE biplots is provides (Figure 1). This is
                                an effective way to visualize the which-won-where pattern of MET data. Cultivars in the
                                vertices of the polygon (Fun,Zav, Ena, Kat and Luc) are those with the longest vectors,
                                in their respective directions, which is a measure of the ability to respond to
                                environments. The vertex cultivars are, therefore, among the most responsive cultivars;
                                all other cultivars are less responsive in their respective directions.",
                               align = "justify"),
                       tags$p("The dotted lines are perpendicular to the polygon sides and divide the biplot into
                               mega-environments, each of which has a vertex cultivar, which is the one with the
                               highest yield (phenotype) in all environments found in it. OA93 and KE93 are in the
                               same sector, separated from the rest of the biplot by two perpendicular lines, and
                               Zav is the highest-yielding cultivar in this sector. Fun is the highest-yielding
                               cultivar in its sector, which contains seven environments, namely, EA93, BH93, HW93,
                               ID93, WP93, NN93, and RN93. No environments fell in the sectors with Ena, Kat, and
                               Luc as vertex cultivars. This indicates that these vertex cultivars were not the
                               best in any of the test environments. Moreover, these cultivars were the poorest
                               in some or all of the environments.",
                               align = "justify"),
                       tags$div(
                       tags$img(src="GGE_WhichWonWhereWhat.png", height="50%", width="50%", align = "middle"),
                       tags$h5("Figure 1: polygon view of the GGE biplot, showing which cultivars presented
                                      highest yield in each environment. The scaling method used is symmetrical
                                      singular value partitioning (by default). The 78% of G + GE variability is
                                      explained by the first two multiplicative terms. Cultivars are shown in
                                      lowercase and environments in uppercase."),
                       align = "center"
                       ),
                       br(),
                       tags$br(),
                       tags$h4(strong("Authors"),
                               align = "justify"
                               ),
                       tags$ul(
                         tags$li("Julia Angelini"),
                         tags$li("Marcos Prunello"),
                         tags$li("Gerardo Cervigni"),
                         align = "justify"
                         ),
                       tags$h5(strong("Maintainer:"), "Julia Angelini"),
                       a(actionButton(inputId = "email1",
                                      label = "Contact Maintainer",
                                      icon = icon("envelope",
                                                   lib = "font-awesome")
                                       ),
                          href="mailto: jangelini_93@hotmail.com.com"),
                       width = 11),
                      br(),br()
             ),
             tabPanel("How to use the APP",
                    navlistPanel(widths = c(3, 9),
                          tabPanel("Data",
                                 mainPanel(
                                   tags$h4(strong("Preparing a data file for the Geneticae APP")),
                                   tags$div(
                                       tags$p("Geneticae APP uses data in .csv format, delimited by commas or semicolons, with
                                               the columns names in the first row of the file (heading). Data in long format
                                               is required, i.e. each row corresponds to one observation and each column to one
                                               variable (genotype, environment, repetition (if any) and the observed phenotype).
                                               If each genotype has been evaluated more than once at each environment, the phenotypic
                                               mean required by SREG and AMMI model for each combination of genotype and environment
                                               is internally calculated and then the model is estimated. Extra variables that will
                                               not be used in the analysis may be present in the dataset. Missing values are not allowed.",
                                              align = "justify"),
                                       tags$p("Two datasets are available in the APP:", align = "justify"),
                                       tags$ol(
                                         tags$li(em("plrv"), " dataset (de Mendiburu, 2020): study about resistance to PLRV (Patato Leaf Roll
                                                 Virus) causing leaf curl. 28 genotypes were experimented at 6 locations in Peru.
                                                 Each clone was evaluated three times in each environment, and yield, plant weight
                                                 and plot were registered.", align = "justify"),
                                         br(),
                                         tags$li(em("yanwinterwheat"), " dataset (Wright, 2020): yield of 18 winter wheat varieties
                                                 grown in nine environments in Ontario at 1993. Although four blocks or replicas in each
                                                 environment were performed in the experiment, only yield mean for each variety and
                                                 environment combination was available.",
                                                 align = "justify")
                                       ),
                                       tags$p("They are available in the tab", em("The data -> Example datasets"), ", and can be downloaded in
                                              .cvs format (Figure 2). The ", em("yanwinterwheat"), "dataset will be used to illustrate the methodology
                                              included in Geneticae APP to analyse MET data.",
                                              align = "justify")
                                   ),
                                   br(),
                                   tags$div(
                                   tags$img(src="dataexamples.png", height="100%", width="100%"),
                                   tags$h5(strong("Figure 2: (A)"),em("Plrv"), "dataset", strong("(B) "),em("yanwinterwheat"), "dataset"),
                                   align="center"),
                                   br(),br(),
                                   tags$div(
                                     tags$h4(strong("Loading a dataset into the APP")),
                                     tags$p("The dataset to be analyzed must be loaded in the ", em (" Data -> User data "), "tab.
                                            For example, to import the ",em("yanwinterwheat"), "dataset, the .csv file must be uploaded,
                                            indicating that it is delimited by comma, that the first row contains the names of
                                            each variable (header) and also the column names for genotype, environment and
                                            phenotypic trait information (gen, env and yield in this case, see Figure 3).
                                            If repetitions are available, specify the name of that column, otherwise, do not
                                            indicate anything.",
                                            align = "justify"),
                                     tags$div(
                                      tags$img(src="data.png", height="60%", width="60%"), br(),
                                      tags$h5(strong("Figure 3"),"Loading", em("yanwinterwheat"), "dataset to Geneticae APP"),br(),br(),
                                      align="center"
                                      ),
                                   ),
                                   width = 11
                                   )
                                 ),
                          tabPanel("Descriptive analysis",
                                 mainPanel(
                                   tags$h4(strong("Descriptive analysis of the dataset")),
                                   tags$p("Any study should start with a descriptive analysis of the dataset. The", em ("Descriptive Analysis"),
                                          "tab provides some tools for this. Boxplots, correlation plots and matrix
                                          and interaction plots may be obtained.",
                                          align = "justify"),
                                   tags$p("A boxplot comparing the quantitative trait across environments (Figure 4) or
                                          across genotypes (Figure 5) may be one of the plots of interest. The summary measures used
                                          for its construction are shown interactively by moving the mouse within the figure panel.
                                          In addition, it can be downloaded as an interactive file (.html) as well as a .png file,
                                          by clicking on the Download button or on the camera that appears in the figure,
                                          respectively (Figure 4). Some aspects of the graph can be customized by the user,
                                          such as box color and axes names.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                    tags$img(src="boxplots.png", height="100%", width="100%"), br(),
                                    tags$h5(strong("Figure 4: "),"boxplot of" , strong("(A) "), "genotypes and ", strong("(B) "),
                                            "environments for", em("yanwinterwheat"), "dataset"),
                                    align="center"),
                                   br(),
                                   br(),
                                   tags$p("Pearson or Spearman correlation coefficients between genotypes can be shown as a plot
                                          or a matrix (Figures 5). Positive correlations are shown in blue and negative in red, while
                                          the intensity of the color and the size of the circle are proportional to the correlation
                                          coefficients. The correlation plot can be downloaded in .png format. High correlations are
                                          observed between the yield of the genotypes studied.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                    tags$h5(strong("Poner captura de matriz y graf de correlacion, no de los graficos los dos graficos
                                                   de correlacion... o sea por ej para los genotipos ambas formas de representar")),
                                    tags$img(src="corrplot.png", height="100%", width="100%"), br(),
                                    tags$h5(strong("Figure 5:"),"Correlation (A) plot and (B) matrix between genotypes",
                                            em("yanwinterwheat"), "dataset"),
                                    align="center"),
                                   br(),
                                   br(),
                                   tags$p("Since GEI generates differential genotypic responses in different environments, which
                                          complicates the task of selecting the best cultivars, an interaction plot may be of interest.
                                          The change in genotypic effect across environments are shown in Figure 6-A, while the change
                                          in the environmental effect through genotypes in Figure 6-B. It is also an interactive
                                          plot, and therefore, it is possible to download it in .HTML or .png
                                          formats with the Download button or clicking on the camera, respectively. Additionally,
                                          axes names can be customized by the user. In this example inconsistencies in the performance
                                          of genotypes in different environments can be seen.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                    tags$img(src="int_plot.png", height="100%", width="100%"), br(),
                                    tags$h5(strong("Figure 6:"),"Interaction plot for (A) environments through genotypes and (B)
                                            genotypes through environments of", em("yanwinterwheat"), "dataset"),
                                    align="center"),
                                   br(),
                                   width = 10
                                 )
                         ),
                         tabPanel("GGE Biplot",
                                 mainPanel(
                                   h4(strong("Site regression model")),
                                   tags$p("The Site Regression model (SREG, also called genotype plus genotype-by-environment
                                           model or GGE model) is another powerful tool for the analysis and interpretation of MET data
                                           in breeding programs. In this case, an ANOVA is performed to obtain estimates for the additive
                                           main effects of environments and a SVD is performed on the residuals matrix in order to explore
                                           patterns related to genotype (G) and GEI. The result of the first two multiplicative terms of the
                                           SVD is often presented in a biplot called genotype plus genotype x environment interaction (GGE)
                                           (Yan et al., 2000).  The GGE biplot addresses many issues relative to genotype and test environment
                                           evaluation. Considering the average performance of each genotype, this plot can be used to evaluate
                                           specific and general adaptation. In addition, environments can be visually grouped according to their
                                           ability to discriminate among genotypes and their representativeness of other test environments.
                                           GGE biplot reveals the which-won-where pattern and allows to recommend specific genotypes for
                                           each environment (Yan and Tinker, 2005).",
                                          align = "justify"),
                                   tags$p("The ", em (" GGE Biplot "), "tab builds several GGE biplots views, in which cultivars are
                                           shown in lowercase and environments in uppercase. Since the model requires a single observation
                                           for each combination of genotype and environment, if there are repetitions, the phenotypic average
                                           value is automatically calculated before fitting the model. Missing values are not allowed.",
                                          align = "justify"),
                                   tags$p("The SVD method must be selected, however, this choice does not alter the relationships or
                                           relative interactions
                                           between genotypes and environments, although the appearance of the biplot will be different (Yan, 2002).
                                           Symmetrical option allows comparison for both genotypes and environments; Genotype-Focused displays
                                           the interrelationship among genotypes more accurately than any other method does, and Environment-focused
                                           is most informative of interrelationships among environments.
                                          A footnote indicating that the centering method is tester-center to obtain GGE bioplot, no scaling
                                          is applied to the data, SVD method selected by the user and the the percentage of G + GEI variation
                                          explained by the two axes can be added. The title graph, axes and axis names can be configure to
                                          appear or not. Also, the genotype andenvironments marker color and size an be customized by the user.",
                                          align = "justify"),
                                   tags$p("A basic GGE biplot is produced by default (Figure 7). In this example the 78% of G + GE variability
                                          is explained by the fist two multiplicative terms. The angles between genotypes markers and environments
                                          vectors are considered to understand this plot. Thus, for example, Kat performs below the average in all
                                          environments, as it has an angle greater than 90° with all environments. On the other hand, Fun presents
                                          an above-average performance in all locations except OA93 and KE93, as indicated by the acute angles.
                                          The length of the environment vectors is a measure of the environment’s ability to discriminate between crops.",
                                         align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$img(src="biplot_GGE.png", height="100%", width="100%"), br(),
                                     tags$h5(strong("Figure 7:"),"GGE biplot based on yield data of 1993 Ontario winter wheat performance trials.
                                             The scaling method used is symmetrical singular value partitioning (by default). The 78% of G + GE
                                             variability is explained by the first two multiplicative terms. Cultivars are shown in lowercase and
                                             environments in uppercase."),
                                     align="center"),
                                   br(),
                                   tags$p("Breeders usually want to identify the most suitable cultivars for a particular environment
                                   of interest, i.e., OA93. To do this with GGE biplots, Yan and Hunt (2002) suggest drawing a line that
                                   passes through the environment marker and the biplot origin, which may be referred to as the OA93 axis.
                                   The performance of the cultivars in this particular environment can be ranked projecting them onto this axis.
                                   This can be done by setting Selected Environment in plot type and providing the name of the environment (OA93)
                                   in Environment selected (Figure 8). Thus, at OA93, the highest-yielding cultivar was Zav, and the lowest-yielding cultivar
                                   was Luc. The line that passes through the biplot origin and is perpendicular to the OA93 axis separates genotypes
                                   that yielded above and below the mean in this environment.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$h5(strong("Figure 8:"),"comparison of cultivar performance in a selected environment (OA93). The scaling
                                             method used is symmetrical singular value partitioning (by default). The 78% of G + GE variability is
                                             explained by the first two multiplicative terms. Cultivars are shown in lowercase and environments
                                             in uppercase."),
                                     align="center"),
                                   br(),
                                   tags$p("Another goal of plant breeders is to determine which is the most suitable environment for
                                   a genotype. Yan and Hunt (2002) suggest plotting a line that passes through the origin and a cultivar
                                   marker, i.e., Kat. To obtain this GGE biplots view the option Selected Genotype in plot type and
                                   the name of the genotype of interest in Genotype selected must be indicated (Figure 9). Environments
                                   are classified along the genotype axis in the direction indicated by the arrow. The perpendicular axis
                                   separates the environments in which the cultivar presented a performance below or above the average.
                                   In this example, Kat presented a performance below the average in all the environments studied.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$h5(strong("Figure 9:"),"comparison of the performance of cultivar Luc in different environments.
                                     The scaling method used is symmetrical singular value partitioning (by default). The 78% of G + GE
                                     variability is explained by the first two multiplicative terms. Cultivars are shown in lowercase and
                                     environments in uppercase." ),
                                     align="center"),
                                   br(),
                                   tags$p("It is also possible to compare two cultivars, i.e. Kat and Cas, linking them with a line and a
                                   segment perpendicular to it. To obtain this GGE biplots view the option Comparison of Genotype in plot type
                                   and the genotypes to be compared Kat and Cas must be indicated in Genotype 1 selected and To be compared with
                                   Genotype 2 selected (Figure 10). Cas was more yielding than Kat in all environments as they all are in the
                                   same side of the perpendicular line as Cas.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$h5(strong("Figure 10:"),"comparison of the cultivars Kat and Cas. The scaling method used is symmetrical
                                             singular value partitioning (by default). The 78% of G + GE variability is explained by the first
                                             two multiplicative terms. Cultivars are shown in lowercase and environments in uppercase. " ),
                                     align="center"),
                                   br(),
                                   tags$p("If Which Won Where/What is indicate in plot type, and the other options are
                                            left by default the polygonal view of the GGE biplots is provides (Figure 11).This view of the GGE
                                            biplots provides an effective way to visualize the which-won-where pattern of MET data, allowing
                                            the identification of the best cultivar in each environment or mega-environment. Cultivars in the vertices
                                            of the polygon (Fun,Zav, Ena, Kat and Luc) are  those with the longest vectors, in their respective directions,
                                            which is a measure of the ability to respond to environments. The vertex cultivars are, therefore, among the
                                            most responsive cultivars; all other cultivars are less responsive in their respective directions.",
                                          align = "justify"),
                                   br(),
                                   tags$p("The dotted lines are perpendicular to the polygon sides and divide the biplot into mega-environments,
                                           each of which has a vertex cultivar, which is the one with the highest yield (phenotype) in all environments
                                           found in it. OA93 and KE93 are in the same sector, separated from the rest of the biplot by two perpendicular
                                           lines, and Zav is the highest-yielding cultivar in this sector. Fun is the highest-yielding cultivar in its sector,
                                           which contains seven environments, namely, EA93, BH93, HW93, ID93, WP93, NN93, and RN93. No environments fell in
                                           the sectors with Ena, Kat, and Luc as vertex cultivars. This indicates that these vertex cultivars were not the best
                                           in any of the test environments. Moreover, these cultivars were the poorest in some or all of the environments.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$img(src="GGE_WhichWonWhereWhat.png", height="100%", width="100%"), br(),
                                     tags$h5(strong("Figure 11:"),"polygon view of the GGE biplot, showing which cultivars presented highest yield in
                                             each environment. The scaling method used is symmetrical singular value partitioning (by default).
                                             The 78% of G + GE variability is explained by the first two multiplicative terms. Cultivars are shown
                                             in lowercase and environments in uppercase. "),
                                     align="center"),
                                   br(),
                                   strong("DE ACA PARA ABAJO NO CORREN PORQUE NO ME SALE QUE INDIQUEN UN MEGA AMBIENTE, ES DECIR QUE COMPLETEN
                                          LOS NOMBRES DE LOS AMBIENTES EN UN MEGA AMBIENTE Y SE ANALICE SOLO ESE SUBCONJUNTO DE DATOS"),
                                   br(),
                                   tags$p("Selecting cultivars within each mega-environments is an issue among plant breeders. Figure 12 clearly suggests that
                                          Zav is the best cultivar for OA93 and KE93, and Fun is the best cultivar for the other locations. However, breeders
                                          do not select a single cultivar in each megaenvironment. Instead, they evaluate all cultivars in order to get an
                                          idea of their performance (yield and stability).",
                                          align = "justify"),
                                   br(),
                                   tags$p("In the GGE biplot it is also possible to visualize mean yield and stability of genotypes in yield units per se
                                          (Figures 12 and 13). The GGE biplot based on genotype-focused scaling, obtained indicating the option row in SVP,
                                          provides an useful way to visualize both mean performance and stability of the tested genotypes. This is because
                                          the unit of both axes for the genotypes is the original unit of the data.",
                                          align = "justify"),
                                   br(),
                                   tags$p("Visualization of the mean and stability of genotypes is achieved by drawing an average environment coordinate (AEC).
                                          For example, Figure 12 shows the AEC for the mega-environment composed of he environments BH93, EA93, HW93, ID93, NN93,
                                          RN93, WP93. The abscissa represents the G effect, thus, the cultivars are ranked along the AEC abscissa. Cultivar Fun
                                          was clearly the highest-yielding cultivar, on average, in this mega-environment, followed by Cas and Har,and Kat was
                                          the poorest. The AEC ordinate approximate the GEI associated with each genotype, which is a measure of the variability
                                          or instability of the genotype. Rub and Dia are more variable and less stable than other cultivars, by the contrary,
                                          Cas, Zav, Reb, Del, Ari, and Kar, were more stable.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$h5(strong("Figure 12:"),"average environment view of the GGE biplot based on genotype-focused scaling, showing mean
                                             yield and stability of genotypes. Cultivars are shown in lowercase and environments in uppercase."),
                                     align="center"),
                                   br(),
                                   tags$p("Figure 13 compares the cultivars to the “ideal” one with the highest yield and absolute stability.
                                          This ideal cultivar is represented by a small circle and is used as a reference, as it rarely exists.
                                          The distance between cultivars and the ideal one can be used as a measure of convenience. Concentric
                                          circles help to visualize these distances. In the example, Fun is the closest one to the ideal crop,
                                          and therefore the most desirable one, followed by Cas and Hay, which in turn are followed by Rum, Ham,
                                          Rub, Zav, Del and Reb, etc.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$h5(strong("Figure 13:")," Classification of genotypes with respect to the ideal genotype. Genotype-focused
                                             scaling is used. Cultivars are shown in lowercase and environments in uppercase."),
                                     align="center"),
                                   br(),
                                   tags$p(" Although METs are performed to study cultivars, they are equally useful for the analysis of the
                                          environments. This includes several aspects: (i) evaluating whether the target region belongs
                                          to one or more megaenvironments; (ii) identifying better test environments; (iii) detecting
                                          redundant environments that do not provide additional information on cultivars; and (iv) determining
                                          environments that can be used for indirect selection. To obtain GGE biplots for comparing environments
                                          the environment-focused scaling should be used as is most informative of interrelationships among them
                                          (Figure 14 and 15). This is obtained indicatig column in SVP.",
                                          align = "justify"),
                                   br(),
                                   tags$p("In Figure 14 environments are connected to the origin through vectors, allowing us to understand the
                                          interrelationships between. The coefficient of correlation between two environments it is approximated
                                          by the cosine of the angle formed by the respective vectors. In this example the relation between the
                                          environments for the mega-environment with BH93, EA93, HW93, ID93, NN93, RN93 and WP93 is considered.
                                          The angle between the vectors for the environments NN93 and WP93 is approximately 10º; therefore, they
                                          are closely related; while RN93 and OA93 present a weak negative correlation since the angle is slightly
                                          greater than 90º. The cosine of the angles does not translate precisely into coefficients of correlation,
                                          since the biplot does not explain all the variability in the dataset. However, they are informative enough
                                          to understand the interrelationship between test environments.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$h5(strong("Figure 14:")," Relationship between environments. Environment-focused scaling is used. "),
                                     align="center"),
                                   br(),
                                   tags$p("Discrimination ability as well as representativeness with respect to the target environment are
                                          fundamental measures for an environment. An ideal test environment should be both discriminating and
                                          representative. If it does not have the ability to discriminate, it does not provide information on
                                          cultivars and is therefore of no use. At the same time, if it is not representative, not only does it
                                          lack usefulness but it can also provide biased information on the evaluated cultivars.",
                                          align = "justify"),

                                   tags$p("To visualize these measurements, an average environment coordinate is defined and the center of
                                          a set of concentric circles represents the ideal environment. Figure 15 shows the GGE biplots view
                                          for the mega-environment with BH93, EA93, HW93, ID93, NN93, RN93 and WP93. The angle between the
                                          vector of an environment and the AEC provides a measure of representativeness. Therefore, EA93 and
                                          ID93 are the most representative, while RN93 and BH93 are the least representative of the average
                                          environment, when the mega-environment is analyzed. On the other hand, an environment to be
                                          discriminative must be close to the ideal environment. HW93 is the closest to ideal environment
                                          and therefore the most desirable of the mega-environment, followed by EA93 and ID93. By the contrary,
                                          RN93 and BH93 were the least desirable test environments of this mega-environment.",
                                          align = "justify"),
                                   br(),
                                   tags$div(
                                     tags$h5(strong("Figure 15:"),"classification of environments with respect to the ideal environment.
                                             Environment-focused scaling is used."),
                                     align="center"),
                                   br(),
                                   width = 10
                                 )
                         ),
                         tabPanel("AMMI Biplot",
                                 mainPanel(
                                   h4(strong("Additive main effects and multiplicative interaction model")),
                                   tags$p("The AMMI model (Gauch, 1988) is widely used to analyse the effect of GEI. This model includes
                                          two stages. First, an ANOVA is performed to obtain estimates for the additive main effects of
                                          environments and genotypes. Secondly, the residuals from the ANOVA are arranged in a matrix with
                                          genotypes in the rows and environments in the columns and a singular value decomposition (SVD)
                                          is applied in order to explore patterns related to GEI, still present in the residuals. The result
                                          of the first two multiplicative terms of the SVD is often presented in a biplot called GE and
                                          represents a two-rank approximation of GEI effects.",
                                          align = "justify"),


                                   tags$p("The ", em (" AMMI Biplot "), "tab builds GE biplot, in which cultivars are shown in lowercase and
                                          environments in uppercase. Since either the clasic and the robust alternatives requires a single
                                          observation for each combination of genotype and environment, if there are repetitions, the phenotypic
                                          average value is automatically calculated before fitting the model. Missing values are not allowed.
                                          As in GGE biplot, a footnote indicating that the percentage of GEI variation explained by the two axes, title graph,
                                          axes and axis names can be configure to appear or not. Also, the genotype and environments marker
                                          color and size an be customized by the user.",
                                          align = "justify"),


                                  tags$p("The GE biplot obtained using AMMI clasic for yan.winterwheat dataset is obtain by default (Figure 16).
                                          In this example, BH93, KE93 and OA93 are the environments that contributethe most to the
                                          interaction as their vectors are the longest ones. The cultivars m12 and Kat present similar
                                          interaction patterns (their markers are close to each other in the biplot) and they are very
                                          different from Ann and Aug, for example. The closeness between the cultivar Dia and the
                                          environment BH93 indicates a strong positive association between them, which means that BH93
                                          is a extremely favorable environment for that genotype. As OA93 and Luc markers are opposite,
                                          this environment is considerably unfavorable for that genotype. Finally, Cas and Reb are close
                                          to the origin, which means that they adapt equally to all environments.",
                                          align = "justify"),
                                  br(),
                                  tags$div(
                                   tags$img(src="AMMI_S.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 16:"),"GE biplot based on yield data of 1993 Ontario winter wheat performance
                                           trials. The 71.66% of GE variability is explained by the first two multiplicative terms.
                                           Cultivars are shown in lowercase and environments in uppercase. "),
                                   align="center"),
                                   br(),
                                   br(), br(),
                                   h4(strong("Robust AMMI models")),
                                   tags$p("The AMMI model, in its standard form, assumes that no outliers are present in the data. To overcome the problem
                                         of data contamination with outlying observations, Rodrigues et al. (2016) proposed five robust AMMI models, which
                                         can be obtained in two stages: (i) fitting a robust regression model with an M-Huber estimator (Huber, 1981) to
                                         replace the ANOVA model; and (ii) using a robust SVD or principal components analysis (PCA) procedure to replace
                                         the standard SVD. Until now, robust AMMI models were not available in any R package. All robust biplots proposed
                                         by Rodrigues et al. (2016) can be obtained using rAMMI(). The argument type can be used to specify the type of
                                         model to be fitted:",
                                   align = "justify"),
                                   br(),
                                   tags$ol(
                                     tags$li("rAMMI"),
                                     tags$li("hAMMI"),
                                     tags$li("gAMMI"),
                                     tags$li("lAMMI"),
                                     tags$li("ppAMMI")
                                   ),
                                   tags$p("Since the sample yan.winterwheat dataset does not present outliers, the conclusions obtained
                                           with robust biplots will not differ from those made with the classic biplot (Rodrigues et al., 2016).
                                           Thus, no interpretation of the robust biplots is presented in this tutorial.",
                                          align = "justify"),
                                   br(),
                                   width = 10
                                 )
                        )
                  )
              ),
             tabPanel("Video-Tutorial",
                      mainPanel(
                        tags$div(
                        tags$h4(strong("Tutorial video: how to use the Geneticae Shiny APP")),
                        HTML('<iframe width="800" height="400" src="https://www.youtube.com/embed/WAODnJDIZ-s" frameborder="0" allowfullscreen></iframe>'),
                        style="text-align: center;")
                      )
             )
  )
)



#------------------------
# Define funcion server
#------------------------

server <- function(input, output, session) {
  observe_helpers(withMathJax = TRUE)


####################
# Load user dataset
####################

datafile <- callModule(csvFile, "datafile",
                         stringsAsFactors = FALSE)

# Show dataset
output$table <- renderDataTable({
dat <- as.data.frame(datafile())
datatable(dat, editable = "cell", class = 'cell-border stripe', rownames = F) %>%
      formatStyle(1:ncol(dat),
                  # background = styleColorBar(c(0, 1), 'lightgray'),
                  backgroundSize = '98% 88%',
                  backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'center', color = "black",
                  backgroundColor = "white",
                  fontWeight = "bold", textAlign = 'center',
                  striped = FALSE)
})

output$select_gen <- renderUI({
checkboxGroupInput(inputId = "select_gen",
                   label = "Genotypes",
                   choices = names(datafile()))
})

output$select_env <- renderUI({
  checkboxGroupInput(inputId = "select_env",
                     label = "Environments",
                     choices = names(datafile()))
})

output$select_rep <- renderUI({
  checkboxGroupInput(inputId = "select_rep",
                     label = "Repetitions",
                     choices = names(datafile()))


})

output$select_pheno <- renderUI({
  checkboxGroupInput(inputId = "select_pheno",
                     label = "Phenotype",
                     choices = names(datafile()))
})


####################
# Examples dataset
####################

# Data without replication
example_withoutrep <- reactive({
  data(yan.winterwheat)
  data_withoutrep <- yan.winterwheat
})

# Data with replication
example_withrep <- reactive({
  data(plrv)
  data_withrep <- plrv
})


# Show example dataset
exampleInput <- eventReactive(input$show_exampledataset, {

    if(input$exampledataset == "without repetitions"){
      data_example <- example_withoutrep()
    }else{
      data_example <- example_withrep()
    }

  knitr::kable(data_example, font_size = 12) %>%
    kable_styling(position = "center")
})

output$show_example <- renderText({
  print(exampleInput())
})


#
# Downloadable csv of sample dataset
output$downloadexample <- downloadHandler(
  filename = function() {
    paste(input$Filename1, '.csv', sep = '')
  },
  content = function(file) {
    write.csv(example_withoutrep(), file, row.names = FALSE, quote = FALSE)
  }
)



####################
# Preprocessing user dataset
####################

sinrep <-reactive({
  data <- datafile() %>%
          group_by(!!sym(input$select_gen), !!sym(input$select_env)) %>%
          summarise(mean_resp = mean(!!sym(input$select_pheno)))%>%
          spread(!!sym(input$select_env), mean_resp) %>%
          as.data.frame()

  rownames(data) <- data[,1]
  data <- data[,-1]
})

dataset <- reactive({
  datos <- datafile() %>%
    mutate(gen = factor(!!sym(input$select_gen)),
           env = factor(!!sym(input$select_env)),
           pheno = !!sym(input$select_pheno))%>%
    as.data.frame()
})

####################
# Descriptive analysis
####################

# Boxplot
boxplotInput <- eventReactive( input$do, {
    if(input$Var == "Environment"){
      boxp_ <-  ggplot(dataset(), aes(x = env, y = pheno))
    } else{
      boxp_ <-  ggplot(dataset(), aes(x = gen, y = pheno))
    }

  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:30) {
                   incProgress(1/20)
                   Sys.sleep(0.30)
                 }
               })

    boxp <- boxp_ +
      geom_boxplot(fill = input$fillcol)+
      labs(y = input$axisy, x = input$axisx) +
      theme_few()+
      theme(text=element_text(family = "Times", size = 9),
            axis.text = element_text(family = "Times", size = 9, angle = 90),
            axis.title = element_text(family = "Times", size = 12))


    ggplotly(boxp) %>%
      config(displaylogo = FALSE,
             modeBarButtonsToRemove = list(
               'sendDataToCloud',
               'autoScale2d',
               'resetScale2d',
               'hoverClosestCartesian',
               'hoverCompareCartesian',
               'zoom2d',
               'pan2d',
               'select2d',
               'lasso2d',
               'zoomIn2d',
               'zoomOut2d',
               'toggleSpikelines'
             )
      )

})

output$boxplot <- renderPlotly({
      print(boxplotInput())
})

output$download_box <- downloadHandler(
  filename = function() {
    paste(input$Filename_box,'.html',sep='')
    },
  content = function(file) {
    htmlwidgets::saveWidget(as_widget(boxplotInput()), file)
  }
)

# Correlation plot
corrInput <- reactive({

  if(input$Var2 == "Environment"){
    M_graf <- cor(sinrep(), method = tolower(input$corrType))
  }else{
    M_graf <- cor(t(sinrep()), method = tolower(input$corrType))
  }

})

corplotInput <- eventReactive( input$do_corrplot, {
  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:20) {
                   incProgress(1/20)
                   Sys.sleep(0.30)
                 }
               })
    ggcorrplot(corrInput(),
               digits = 2, type="upper",
               outline.color = "white",
               ggtheme = ggthemes::theme_few(),
               method = "circle",
               colors = c("darkred", "white", "darkblue")) +
      theme(text = element_text(family = "Times", size = 9),
            axis.text.x = element_text(family = "Times", size = 9,angle = 90),
            axis.text.y = element_text(family = "Times", size = 9),
            axis.title = element_text(family = "Times", size = 9))
})

output$corplot <- renderPlot({
  print(corplotInput())
})

output$downloadcorplot <- downloadHandler(
  filename = function() {
    paste(input$Filename3, '.png', sep='')
    },
  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = corplotInput(), device = device)
  }
  )

# Correlation matrix
cormatInput <- eventReactive( input$do_corrmat, {
  withProgress(message = 'Estimation in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:20) {
                   incProgress(1/20)
                   Sys.sleep(0.30)
                 }
               })
  if(input$Var3 == "Environment"){

    M <- cor(sinrep(), method = tolower(input$correlType))
  }else{
    M <- cor(t(sinrep()), method = tolower(input$correlType))
  }

    knitr::kable(M,"html", digits = 2,full_width = F, font_size = 12) %>%
    kable_styling(position = "center")  %>%
    column_spec(1:ncol(M), bold = T, color = "black") %>%
    row_spec(1:nrow(M), bold = T, color = "black")
})



output$cormat<-renderText({
print(cormatInput())
})

# Interaction plot
datos_summ <- reactive({
  datos <-
    datafile() %>%
    group_by(!!sym(input$select_gen), !!sym(input$select_env)) %>%
    summarise(y = mean(!!sym(input$select_pheno)))

  datos <-data.frame(datos)

})

interacInput <-eventReactive( input$do_int, {
  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:20) {
                   incProgress(1/20)
                   Sys.sleep(0.30)
                 }
               })

  if(input$Var4 == "Environment"){
    intp <- ggplot2::ggplot(data = datos_summ(), aes(x = !!sym(input$select_env), y = y, colour = !!sym(input$select_gen), group = !!sym(input$select_gen)))
  }else{
    intp <- ggplot2::ggplot(data = datos_summ(), aes(x = !!sym(input$select_gen), y = y, colour = !!sym(input$select_env), group = !!sym(input$select_env)))
  }

  intp <- intp +
    stat_summary(fun = mean, geom = "point") +
    stat_summary(fun = mean,geom = "line")+
    labs(y = input$axisy_int, x=input$axisx_int)+
    theme_few() +
    theme(text=element_text(family = "Times", size = 9),
          axis.text=element_text(family = "Times", size = 9),
          axis.title=element_text(family = "Times", size = 9))

  ggplotly(intp) %>%
    config(displaylogo = FALSE,
           modeBarButtonsToRemove = list(
             'sendDataToCloud',
             'autoScale2d',
             'resetScale2d',
             'hoverClosestCartesian',
             'hoverCompareCartesian',
             'zoom2d',
             'pan2d',
             'select2d',
             'lasso2d',
             'zoomIn2d',
             'zoomOut2d',
             'toggleSpikelines'
           )
    )
})

output$int <- renderPlotly({
  print(interacInput())
})

output$download_int <- downloadHandler(
  filename = function() {
    paste(input$Filename_int,'.html',sep='')
    },
  content = function(file) {
    htmlwidgets::saveWidget(as_widget(interacInput()), file)
  }
)


####################
#   GGE Biplot
####################

# Environments inside mega-environment
choices_ME <- reactive({
  datafile() %>%
    select(input$select_env) %>%
    unique()
})

observe({updatePickerInput(session, 'ME', choices = choices_ME())})

datafile_ME <- reactive({
    datafile() %>%
    filter(!!sym(input$select_env) %in% input$ME) %>%
    as.data.frame()
})

# Selection of environment
choices_E <- reactive({
  datafile() %>%
    select(input$select_env) %>%
    unique()
})

observe({updatePickerInput(session, 'SelectedE', choices = choices_E())})

# Selection of genotype
choices_G <- reactive({
  datafile() %>%
    select(input$select_gen) %>%
    unique()
})

observe({updatePickerInput(session, 'SelectedG', choices = choices_G())})

# Selection of genotype 1
choices_G1 <- reactive({
  datafile() %>%
    select(input$select_gen) %>%
    unique()
})

observe({updatePickerInput(session, 'SelectedG1', choices = choices_G1())})

# Selection of genotype 2
choices_G2 <- reactive({
  datafile() %>%
    select(input$select_gen) %>%
    unique()
})

observe({updatePickerInput(session, 'SelectedG2', choices = choices_G2())})


modelInput <- reactive({

  if(input$plotType %in% c("Mean vs. Stability","Ranking Genotypes","Relationship Among Environments","Ranking Environments")){
    if(input$SVP=="Genotype-focused"){
      geneticae::GGEmodel(datafile_ME(), genotype = input$select_gen, environment = input$select_env, rep = input$select_rep, response = input$select_pheno,
                          SVP = "row", centering = "tester", scaling = "none")
    } else if(input$SVP=="Environment-focused"){
      geneticae::GGEmodel(datafile_ME(), genotype = input$select_gen, environment = input$select_env, rep = input$select_rep, response = input$select_pheno,
                          SVP = "column", centering = "tester", scaling = "none")
    }else{
      geneticae::GGEmodel(datafile_ME(), genotype = input$select_gen, environment = input$select_env, rep = input$select_rep, response = input$select_pheno,
                          SVP = "symmetrical", centering = "tester", scaling = "none")
    }

  }else{
    if(input$SVP=="Genotype-focused"){
      geneticae::GGEmodel(datafile(), genotype = input$select_gen, environment = input$select_env, rep = input$select_rep, response = input$select_pheno,
                          SVP = "row", centering = "tester", scaling = "none")
    } else if(input$SVP=="Environment-focused"){
      geneticae::GGEmodel(datafile(), genotype = input$select_gen, environment = input$select_env, rep = input$select_rep, response = input$select_pheno,
                          SVP = "column", centering = "tester", scaling = "none")
    }else{
      geneticae::GGEmodel(datafile(), genotype = input$select_gen, environment = input$select_env, rep = input$select_rep, response = input$select_pheno,
                          SVP = "symmetrical", centering = "tester", scaling = "none")
    }
  }
})

plotGGEInput <- eventReactive( input$do_GGE, {
  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:25) {
                   incProgress(1/25)
                   Sys.sleep(0.25)
                 }
               })

  if(input$plotType=="Selected Environment"){

    p <-  geneticae::GGEPlot(modelInput(), type = input$plotType, selectedE = input$SelectedE, footnote = input$footnote,
                           colGen = input$colgen, colEnv = input$colenv, colSegment = input$colsegment,
                           titles = input$title, sizeGen = as.numeric(input$sizeGen), sizeEnv = as.numeric(input$sizeEnv),
                           axislabels = input$axislabels, axes = input$axes)
  } else if(input$plotType=="Selected Genotype"){

    p <-  geneticae::GGEPlot(modelInput(), type = input$plotType, selectedG = input$SelectedG, footnote = input$footnote,
                             colGen = input$colgen, colEnv = input$colenv, colSegment = input$colsegment,
                             titles = input$title, sizeGen = as.numeric(input$sizeGen), sizeEnv = as.numeric(input$sizeEnv),
                             axislabels = input$axislabels, axes = input$axes)

  } else if(input$plotType=="Comparison of Genotype"){

    p <-  geneticae::GGEPlot(modelInput(), type = input$plotType, selectedG1 = input$SelectedG1, selectedG2 = input$SelectedG2,
                             footnote = input$footnote, colGen = input$colgen, colEnv = input$colenv, colSegment = input$colsegment,
                             titles = input$title, sizeGen = as.numeric(input$sizeGen), sizeEnv = as.numeric(input$sizeEnv),
                             axislabels = input$axislabels, axes = input$axes)
  }else{
    p <- geneticae::GGEPlot(modelInput(), type = input$plotType, footnote = input$footnote,
                             colGen = input$colgen, colEnv = input$colenv, colSegment = input$colsegment,
                             titles = input$title, sizeGen = as.numeric(input$sizeGen), sizeEnv = as.numeric(input$sizeEnv),
                             axislabels = input$axislabels, axes = input$axes)
  }

})

output$plotGGE <- renderPlot({
  print(plotGGEInput())
})

output$download_gge <- downloadHandler(
  filename = function() {
    paste(input$Filename_GGE, '.png', sep='')
    },
  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = plotGGEInput(), device = device)
  }
  )

####################
#   GE Biplot
####################
plotAMMIInput <- eventReactive( input$do_AMMI, {
  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:20) {
                   incProgress(1/20)
                   Sys.sleep(0.30)
                 }
               })
  if(!is.null(input$select_rep)){
  p <- geneticae::rAMMI(dataset(), genotype = input$select_gen, environment = input$select_env, rep = input$select_rep, response = input$select_pheno,
                         type = input$robustType, footnote = input$foot, colGen = input$colorgen, colEnv = input$colorenv,
                         titles = input$tit, sizeGen = as.numeric(input$sizeGen_AMMI), sizeEnv = as.numeric(input$sizeEnv_AMMI),
                         axislabels = input$axislabels_AMMI,axes = input$axes_AMMI)
  }else{
  p <- geneticae::rAMMI(dataset(), genotype = input$select_gen, environment= input$select_env, rep = NULL, response = input$select_pheno,
                         type = input$robustType, footnote = input$foot, colGen = input$colorgen, colEnv = input$colorenv,
                         titles = input$tit, sizeGen=as.numeric(input$sizeGen_AMMI), sizeEnv = as.numeric(input$sizeEnv_AMMI),
                         axislabels = input$axislabels_AMMI, axes = input$axes_AMMI)

  }

  })


output$plotAMMI <- renderPlot({

  print(plotAMMIInput())
})

output$download_ammi <- downloadHandler(
  filename = function() {
    paste(input$Filename_AMMI, '.png', sep='')
    },
  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = plotAMMIInput(), device = device)
  }
  )

}

#-----------------------
# Create Shiny app
#-----------------------

shinyApp(ui, server)




