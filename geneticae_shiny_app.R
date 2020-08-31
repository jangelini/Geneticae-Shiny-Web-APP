#---------------------------------------------------------------------------------------
# Archivo: geneticae_shiny_app
# Descripcion: este archivo genera una aplicación Shiny que permite analizar datos
#              provenientes de ensayos multiambientales. Análisis descriptivos, ANOVA,
#              y biplots se pueden realizar. Estos últimos se obtienen del paquete
#              geneticae.
# Autora: Julia Angelini
#---------------------------------------------------------------------------------------

#-----------
# Paquetes
#-----------

rm(list=ls(all=TRUE))
library(shiny)
library(ggplot2)
library(tidyr)
# library(corrplot)
library(shinythemes)
library(DT)
# library(formattable)
library(geneticae)
library(dplyr)
library(kableExtra)
# library(RColorBrewer)
library(ggfortify)
library(ggcorrplot)
library(shinyWidgets)
library(ggthemes)
library(plotly)
library(shinyhelper)

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
  tabPanel("The data",
           tabsetPanel(
                 tabPanel("User data", icon = icon("table"),
                          sidebarPanel(
                                csvFileInput("datafile", "User data (.csv format)"), br(),
                                tags$p("Select the column name that containts:"),
                                uiOutput("select_gen"),
                                uiOutput("select_env"),
                                uiOutput("select_rep"),
                                uiOutput("select_pheno"),
                                width = 3
                                ),
                          mainPanel(dataTableOutput("table"), width=8)
                        ),
                 tabPanel(strong("Examples dataset"), icon = icon("table"),
                          sidebarPanel(br(), strong("Example without repetitions"), br(),
                                actionButton("example_withoutrep", "Show"), br(), br(),
                                strong("Download sample dataset"),
                                textInput("Filename1", "File name", value = "Example without repetitions"),
                                downloadButton("downloadexample_withoutrep", "without repetitions"),
                                br(), br(), br(),
                                strong("Example with repetitions"),
                                br(),
                                actionButton("example_withrep", "Show"), br(), br(),
                                strong("Download sample dataset"),
                                textInput("Filename2", "File name", value = "Example with repetitions"),
                                downloadButton("downloadexample_withrep", "with repetitions"),
                                width = 3
                                ),
                          mainPanel(htmlOutput("show_example"), width=8)
                          )
             )
  ),
  tabPanel("Descriptive analysis",
           tabsetPanel(
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
                                           content = c("Select if you want to plot the phenotypic trait across genotypes or environments"),
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
                                           content = c("a <b>.HTML</b> file will be download"),
                                           size = "m"),
                                  downloadButton("download_box", "Download"),
                                  width=3
                                  ),
                        mainPanel(br(),plotlyOutput("boxplot"), width=6)
                        )
                  ),
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
                                           content = c("Select if you want to plot the phenotypic trait correlations between genotypes or environments"),
                                           size = "m"),
                                  prettyRadioButtons(
                                      inputId = "corrType",
                                      label = "Correlation index",
                                      choices = c("spearman", "pearson"),
                                      icon = icon("check"),
                                      status = "info",
                                      animation = "rotate"
                                      )%>%
                                    helper(type = "inline",
                                           title = "Inline Help",
                                           content = c("Indicate which correlation coefficient is to be computed"),
                                           size = "m"),
                                  actionButton("do_corrplot", "Run"), br(), br(),
                              br(),
                              textInput("Filename3", "File name", value = "Correlation plot")%>%
                                helper(icon = "exclamation",
                                       colour = "red",
                                       type = "inline",
                                       title = "Inline Help",
                                       content = c("a <b>.png</b> file will be download"),
                                       size = "s"),
                              downloadButton("downloadcorplot", "Download"),
                              width=3
                              ),
                            mainPanel(br(), plotOutput("corplot"), width=6))
             ),
             tabPanel("Correlation matrix", icon = icon("table"),
                        sidebarLayout(
                            sidebarPanel(
                                prettyRadioButtons(
                                    inputId = "Var3",
                                    label = "Obtain for...",
                                    choices = c("Environment", "Genotype"),
                                    icon = icon("check"),
                                    status = "info",
                                    animation = "rotate"
                                    )%>%
                                  helper(type = "inline",
                                         title = "Inline Help",
                                         content = c("Select if you want to estimate the phenotypic trait correlations between genotypes or environments"),
                                         size = "m"),
                                prettyRadioButtons(
                                    inputId = "correlType",
                                    label = "Correlation index",
                                    choices = c("spearman","pearson"),
                                    icon = icon("check"),
                                    status = "info",
                                    animation = "rotate"
                                    )%>%
                                  helper(type = "inline",
                                         title = "Inline Help",
                                         content = c("Indicate which correlation coefficient is to be computed"),
                                         size = "m"),
                                actionButton("do_corrmat", "Run"), br(), br(),
                                width=3),
                            mainPanel(br(), htmlOutput("cormat"))
                            )
             ),
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
                                           content = c("Select if you want to plot the phenotypic trait between genotypes or environments"),
                                           size = "m"),
                              textInput("axisx_int", "X-axis label:", value = "Environment"),
                              textInput("axisy_int", "Y-axis label:", value = "Yield"),
                              actionButton("do_int", "Run"), br(), br(),
                              textInput("Filename_int", "File name", value = "Interaction plot")%>%
                                helper(type = "inline",
                                       title = "Inline Help",
                                       content = c("a <b>.HTML</b> file will be download"),
                                       size = "m"),
                              downloadButton("download_int", "Download"),
                              width=3),
                          mainPanel(br(), plotlyOutput("int"), width=6)
                          )
             )
       )
 ),
tabPanel("Analysis of variance",
         tabsetPanel(
              tabPanel("ANOVA" , icon = icon("table"), br(),br(),
                      mainPanel( htmlOutput("Anova"), verbatimTextOutput("Anova2"))
                      ),
             tabPanel("Check Normality", icon = icon("bar-chart-o"),
                      sidebarLayout(
                           sidebarPanel( br(),
                                strong("Histogram residuals"),
                                pickerInput(
                                    inputId = "fillhist",
                                    label = "Fill color",
                                    choices = c("orange","dimgrey", "red","lightpink","white", "blue")
                                ),
                                textInput("Filename_hist", "File name", value = "Histogram"),
                                downloadButton("downloadHist", "Download"),
                                br(), br(),
                                strong("QQ-Plot residuals"),
                                pickerInput(
                                    inputId = "fillqq",
                                    label = "Fill color",
                                    choices = c("orange","dimgrey", "red","lightpink","white", "blue")
                                ),
                                textInput("Filename_qq", "File name", value = "QQPlot"),
                                downloadButton("downloadqqnorm", "Download"), br(),br(),
                                strong("Shapiro-Wilks normality test"),
                                actionButton("do_test", "Run"),
                                width=3
                          ),
                          mainPanel(plotlyOutput("residual_hist"), br(),
                                    plotlyOutput("residual_qqnorm"), br(),
                                    verbatimTextOutput("residual_test")
                          )
                    )
             ),
             tabPanel("Check Homoscedasticity", icon = icon("bar-chart-o"),
                      sidebarLayout(
                          sidebarPanel( br(),
                              pickerInput(
                                  inputId = "fillhomo",
                                  label = "Fill color",
                                  choices = c("orange","dimgrey", "red","lightpink","white", "blue")
                              ),
                              textInput("Filename_homo", "File name", value = "Check Homocedasticity"),
                              downloadButton("downloadlinearity", "Download"), br(), br(),
                              strong("Levene test for environments"),
                              actionButton("do_leve_env", "Run"), br(), br(),
                              strong("Levene test for genotypes"),
                              actionButton("do_leve_gen", "Run"),
                              width=3
                         ),
                      mainPanel(plotOutput("residual_linearity"))
                      )
             ),
             tabPanel("Outliers", icon = icon("bar-chart-o"),
                      sidebarLayout(
                          sidebarPanel( br(),
                              pickerInput(
                                  inputId = "fillout",
                                  label = "Fill color",
                                  choices = c("orange","dimgrey", "red","lightpink","white", "blue")
                              ),
                              textInput("Filename_out", "File name", value = "Outliers"),
                              downloadButton("downloadoutliers", "Download"),
                              width=3
                          ),
                          mainPanel(plotOutput("residual_outliers"))
                      )
              )
        )
  ),
  tabPanel("GGE Biplot", icon = icon("bar-chart-o"),
           sidebarLayout(
                sidebarPanel(
                  prettyRadioButtons(
                    inputId = "SVP",
                    label = "SVP type",
                    choices = c("symmetrical","row", "column","dual"),
                    icon = icon("check"),
                    status = "info",
                    animation = "rotate"
                  ),
                  prettyRadioButtons(
                    inputId =  "center",
                    label = "center type",
                    choices = c("tester", "global","double"),
                    icon = icon("check"),
                    status = "info",
                    animation = "rotate"
                  ),
                  prettyRadioButtons(
                    inputId = "scale",
                    label = "scale type",
                    choices = c("none","sd"),
                    icon = icon("check"),
                    status = "info",
                    animation = "rotate"
                  ),
                     prettyRadioButtons(
                         inputId = "plotType",
                         label = "Plot type",
                         choices = c("Biplot","Relationship Among Environments",
                                     "Which Won Where/What","Discrimination vs. representativeness",
                                     "Ranking Environments","Mean vs. Stability","Ranking Genotypes"),
                         icon = icon("check"),
                         status = "info",
                         animation = "rotate"
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
                      selectInput("sizeGen", "Size Gen marker:", choices = c(4,0,1,2,3,5,6,7,8,9,10)),
                      selectInput("sizeEnv", "Size Env marker:", choices = c(4,0,1,2,3,5,6,7,8,9,10)),
                     actionButton("do_GGE", "Run"), br(), br(),

                    br(), br(),
                    textInput("Filename_GGE", "File name", value = "GGE Biplot"),
                    downloadButton("download_gge", "Download"),
                    width = 3
              ),
              mainPanel(plotOutput("plotGGE"))
      )
  ),
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
               selectInput("sizeGen_AMMI", "Size Gen marker:", choices = c(4,0,1,2,3,5,6,7,8,9,10)),
               selectInput("sizeEnv_AMMI", "Size Env marker:", choices = c(4,0,1,2,3,5,6,7,8,9,10)),
               actionButton("do_AMMI", "Run"), br(), br(),
               textInput("Filename_AMMI", "File name", value = "AMMI Biplot"),
               downloadButton("download_ammi", "Download"),
               width = 3
              ),
            mainPanel( plotOutput("plotAMMI"))
        )
  ),

  navbarMenu("Help",icon = icon("question-circle"),
             tabPanel("General",
                      mainPanel(
                        tags$h4(strong("Motivation")),
                        tags$div(
                            tags$p("Understanding the relationship between crops performance and environment is a key problem for plant breeders
                                and geneticists. In advanced stages of breeding programs, multi-environmental
                                trials (MET) is one of the most common experiments. They are conducted by testing a number of genotypes across
                                multiple environments, allowing the identification of superior genotypes. Crop performance, is a function of genotype
                                (G), environment (E) and genotype x environment interaction (GEI). METs are essential due to the presence of GEI
                                which generates differential genotypic responses in the different environments evaluated (Crossa et al., 1990;
                                Cruz Medina, 1992; Kang and Magari, 1996). This is a particular problem for plant breeders (Giauffret et al., 2000),
                                therefore appropriate statistical methods should be used to obtain an adequate GEI analysis.", align = "justify"),
                           tags$p("Nowadays, computer programs have become a basic tool for data analysis. Currently, R is one of the most used programs due to
                                its power and free distribution. This software consists of a large number of packages which can perform different types of analysis.
                                In particular, geneticae package was developed by our group and offers functions for METs analysis. However, R has a complex syntax, and therefore is not friendly
                                for those who do not know the R programming language. Frequently, breeders use programs with graphical user interface
                                to perform statistical analysis but not all of them allow all the analyzes of interest to be carried out and therefore, several must be used to fulfill an objective.
                                This generates the need to have different programs to perform different analyzes, to know the data formats required by each one and to understand the different types
                                  of results that they generate.", align = "justify"),
                            tags$p("The objective of Geneticae APP is to create a graphical interface of geneticae package, removing the obstacle of R programming language complexity.", align = "justify")
                        ), tags$br(),
                        tags$h4(strong("About Geneticae APP")),
                        tags$div(
                            tags$p("Geneticae APP is a statistical Shiny Web APP for phenotypic analyses in the plant breeding context, developed by Julia Angelini, Marcos Prunello
                               and Gerardo Cervigni, as part of the final work to obtain Specialist in Bioinformatics degree of the first author.", align = "justify"),
                            tags$p("Is an interactive, noncommercial and open source software, offering a free alternative to available commercial
                              software to analize METs.", align = "justify")
                        ), tags$br(),
                        tags$h4(strong("Citation and contact"), align = "justify"),
                        tags$h5("..............", align = "justify"), br(),
                        tags$h5("Maintainer: Julia Angelini"),
                        a(actionButton(inputId = "email1", label = "Contact Maintainer",
                                       icon = icon("envelope", lib = "font-awesome")),
                          href="mailto: jangelini_93@hotmail.com.com"),
                        width = 11)
             ),
             tabPanel("How to use the APP",
                    navlistPanel(widths = c(3, 9),
                          tabPanel("The data",
                                 mainPanel(
                                   tags$h4(strong("Preparing data file for Geneticae APP")),
                                   tags$div(
                                       tags$p("Geneticae APP allows data in .csv format, delimited by commas or semiclons; as well as a header in first row of the file.
                                        Since this application connects to R and uses geneticae package functions, the format required by it must be respected.
                                        Observations must be in rows and variables (genotypes, environments, repetitions (if any) and the observed phenotype) in the columns.
                                        Other variables that will not be used in the analysis can be included, since when the data set is loaded, the names of the columns corresponding
                                        to the genotype, environment, repetition (if any) and phenotype of interest must be indicated.", align = "justify"),
                                       tags$p("Two sample datasets, ", em("plrv"), "and", em("yanwinterwheat"), ", included in geneticae R package, are also available in the APP to be downloaded and run the
                                       the examples (Figure 1, 2).", align = "justify"),
                                       tags$ol(
                                         tags$li(em("plrv"), " dataset (de Mendiburu, 2020): yield, plant weight and plot of
                                                 28 genotypes in 6 locations in Peru, in order to study resistance to PLRV (Patato Leaf Roll Virus)
                                                 that causes leaf roll. Each clone was evaluated three times in each environment.", align = "justify"),
                                         br(),
                                         tags$li(em("yanwinterwheat"), " dataset (Wright, 2018): yield of 18 varieties of winterwheat grown
                                                 in nine environments in Ontario in 1993. Although the experiment had four blocks or replications, only the
                                                 yield average for every combination of variety and environment is available.", align = "justify")
                                       )
                                   ), br(),
                                   tags$img(src="Exampledatasets_withrep.png", height="100%", width="100%"),
                                   tags$h5(strong("Figure 1:"),em("Plrv"), "dataset"),br(),br(),
                                   tags$img(src="Exampledatasets_withoutrep.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 2"),em("yanwinterwheat"), "dataset"),br(),br(),

                                   tags$div(
                                     tags$h4(strong("Loading a dataset into the APP")),
                                     tags$p("To load a dataset, for example", em("yanwinterwheat"), ", it must indicated that the .csv file is delimited by
                                            semicolons, the header contains the names of each variable and also the name of the column
                                            containing genotype, environments and phenotypic information (Figure 3). If there are
                                            repeats available, as is the case for the plrv dataset, the column name must also be specified.
                                            As an example, this dataset will be used to show all the analyzes that can be performed using the Geneticae APP.", align = "justify"),
                                     tags$img(src="data.png", height="100%", width="100%"), br(),
                                     tags$h5(strong("Figure 3"),"Loading", em("yanwinterwheat"), "dataset"),br(),br(),
                                   ),
                                   width = 11
                                   )
                                 ),
                          tabPanel("Descriptive analysis",
                                 mainPanel(
                                   tags$h4(strong("Descriptive analysis of the dataset")),
                                   tags$p("Any study should start with a descriptive analysis of the dataset, the Descriptive Analysis tab provides some tools for that first step.
                                          One of the graphs of interest is a boxplot that compares the quantitative trait across environments (Figure 4) or across genotypes (Figure 5).
                                          This graphic is interactive, allowing to obtain the summary measures used for its construction by moving the mouse within it. In addition,
                                          they can be downloaded in interactive format (.HTML) as well as in .png format by clicking on the Download button and on the camera that
                                          appears in the graphic respectively (Figures 4 and 5).  The user can customize some aspects of the graph, such as the color of the
                                          box and the names of the axes.", align = "justify"),
                                   br(),
                                   tags$img(src="Boxplot_genotypes.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 4:"),"Boxplot of genotypes for", em("yanwinterwheat"), "dataset"),
                                   br(),
                                   br(),
                                   tags$img(src="Boxplot_environment.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 5:"),"Boxplot of environments for", em("yanwinterwheat"), "dataset"),
                                   br(),
                                   tags$p("The correlation between genotypes or between environments, for which both the correlation plot and
                                   matrix can be obteined using Pearson or Spearman method (Figures 6 and 7). The positive correlations are shown in blue in the correlogram
                                          and negative in red, the intensity of the color and the size of the circle are proportional to the correlation coefficients (Figure 6).", align = "justify"),
                                   br(),
                                   tags$img(src="corr_env.png", height="100%", width="100%"), br(),
                                   br(),
                                   tags$h5(strong("Figure 6:"),"Correlation plot between genotypes of", em("yanwinterwheat"), "dataset"),
                                   br(),
                                   tags$img(src="corr_env.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 7:"),"Correlation matrix between environments of", em("yanwinterwheat"), "dataset"),
                                   br(),
                                   tags$p("Since inconsistencies in the performance of genotypes in different environments complicate plant breeders job,
                                          an interaction plot may be of interest. The change in genotypic effect across environments are shown in Figure 8,
                                          while the change in the environmental effect through genotypes in Figure 9. It is a
                                          interactive graph, and therefore, it is possible to download in interactive (.HTML) as well as in .png format from the Download button and
                                          when clicking on the camera, respectively (Figure 8 y 9). Additionally, the names of the axes can be customized by the user.", align = "justify"),
                                   br(),
                                   tags$img(src="int_plotenv.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 8:"),"Interaction plot for environments through genotypes of", em("yanwinterwheat"), "dataset"),
                                   br(),
                                   br(),
                                   tags$img(src="int_plot.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 9:"),"Interaction plot for genotypes through environments of", em("yanwinterwheat"), "dataset"),
                                   br(),
                                   width = 10
                                 )
                         ),
                         tabPanel("Analysis of variance",
                                 mainPanel(
                                   h4(strong("ANOVA")),
                                   tags$p("The phenotypic (P) consist of environmental (E), genotypic (G) and gentype-by-environment
                                          interaction (GEI) variation: P=E+G+GEI. If only the effects of G and E are significant
                                          (that is, there is no GEI effect), the interaction must be ignored and the biplots are meaningless.
                                          If the genotypes were evaluated only once in each environment, that is, there are no repetitions, the interaction
                                          could not be tested.",
                                          "Therefore, in ", em("yanwinterwheat"), "dataset, only the genotypic and environmental effects can be tested (Figure 10).", align = "justify"),
                                   br(),
                                   tags$img(src="ANOVA.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 10:"),"Analysis of variance for", em("yanwinterwheat"), "dataset"),
                                   br(),br(),
                                   tags$p("The validity of the ANOVA conclusions depends on whether the errors have a normal distribution with a zero mean and constant variance.
                                        Three tabs: Check normality, Check homocedasticity and Outliers allow verifying the above assumptions.The graphics to verify the assumptions
                                        required by the ANOVA technique can be downloaded using the Download button in the corresponding tab. The normality  can be verified graphically with a histogram and a normal probability
                                        graph (Figure 12). In addition, it can be tested with the shapiro-wilks test on the ANOVA residues (Figure 13).", align = "justify"),
                                   br(),
                                   tags$img(src="Normalidad.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 12:"),"Normality errors asumption for", em("yanwinterwheat"), "dataset"),
                                   br(),br(),
                                   tags$img(src="Normalidad2.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 13:"),"Shapiro-wilks test for", em("yanwinterwheat"), "dataset"),
                                   br(),
                                   tags$p("Homoskedasticity can be tested with residuals vs. predicted values plot, as well as with the levene test (Figure 14).
                                        This test can be seen for both genotypes and environments (Figure 15).", align = "justify"),
                                   br(),
                                   tags$img(src="Homocedasticidad.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 14:"),"Homoskedasticity of errors for", em("yanwinterwheat"), "dataset"),
                                   br(),br(),
                                   tags$img(src="levene.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 15:"),"Levene test for", em("yanwinterwheat"), "dataset"),
                                   br(),br(),
                                   tags$p("Finally, ANOVA is not robust in the presence of atypical observations, therefore graphs are included to detect if outliers (Figure 16).", align = "justify"),
                                   br(),
                                   tags$img(src="Outliers.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 16:"),"Outliers in", em("yanwinterwheat"), "dataset"),
                                   width = 10
                                 )),
                         tabPanel("GGE Biplot",
                                 mainPanel(
                                   h4(strong("Site regression model")),
                                   tags$p("Site regression model (SREG) explore the joint effect of genotype plus genotype-by-environment performing a singular value decomposition
                                      (SVD) on the matrix of residuals from a one-way ANOVA with fixed effects for environments.
                                      The result of the first two multiplicative terms of the SVD is often presented in a biplot
                                      called genotype plus genotype x environment interaction (GGE) (Yan et al., 2000). Such biplot
                                      represents an approximation of the G+GE effects. Increasingly, plant breeders have
                                      found GGE biplots as useful tools in mega environment (ME) analysis (Yan et al., 2001; Yan
                                      and Rajcan, 2002), and genotype and environment evaluation (Bhan et al., 2005; Kang et al., 2006;
                                      Yan et al., 2007). ", align = "justify"),
                                      br(),
                                    h4(strong("GGE biplot")),
                                   tags$p("The GGE biplot visually addresses many issues relative to genotype and test
                                      environment evaluation:", align = "justify"),
                                      tags$ol(
                                        tags$li("Basic biplot."),
                                        tags$li("Relationship Among Environments."),
                                        tags$li("Which won where/what: Identifying the best cultivar in each environment."),
                                        tags$li("Discrimination vs. representativeness."),
                                        tags$li("Ranking environments: Ranking environments with respect to the ideal environment."),
                                        tags$li("Mean vs. stability: Evaluating cultivars based on both average yield and stability."),
                                        tags$li("Ranking gentoypes: Ranking genotypes with respect to the ideal genotype.")
                                      ),
                                   br(),
                                   h4(strong("Some details")),
                                   tags$p("Since the model requires a single observation for each combination of genotype and environment,
                                          if there are repetitions, the phenotypic average value is automatically calculated
                                          before fitting the model. Missing values are not allowed Missing values are NOT allowed.
                                          The centering, SVD, and scaling method must be selected. By default the data
                                          is centered by G and GE, giving rise to the SREG model, another option in this argument will give
                                          rise to a different model. The choice of the SVD method does not alter the relationships or relative
                                          interactions between genotypes and environments, although the appearance of the biplot will be different (Yan, 2002)
                                          . Finally, different biplots can be generated depending on the
                                          scaling method, more information is available in Yan and Kang (2003) and Yan and Tinker (2006).
                                          In these graphs the cultivars are shown in lower case, to differentiate them from the environments,
                                          which are in upper case. The centering, scaling and SVD method used, as well
                                          as the G + GE percentage explained by the two axes can be added in the graph as a footnote to the biplot
                                          s well as the title of the graph, and the name of the axes and if we want them not to appear", align = "justify"),
                                   br(),
                                   tags$img(src="biplot_GGE.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 17:"),"Basic GGE biplot for", em("yanwinterwheat"), "dataset"),
                                   br(),
                                   width = 10
                                 )
                         ),
                         tabPanel("AMMI Biplot",
                                 mainPanel(
                                   h4(strong("Additive main effects and multiplicative interaction model")),
                                   tags$p("The additive main effects and multiplicative interaction (AMMI)
                                    model (Gauch, 1988, 1992) is one of the most widely used tools to
                                    analyse and structure GEI. This model works under a fixed-model
                                    framework and is fit in two stages. First, the main effects of the
                                    model are estimated using the additive two-way analysis of variance (ANOVA) by least squares.
                                    Then, the singular value decomposition (SVD) is applied to the residuals from the ANOVA,
                                    i.e. to the interaction, to obtain the estimates for the multiplicative terms of the AMMI model.", align = "justify"),
                                   br(),
                                   tags$img(src="AMMI_S.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 18:"),"AMMI biplot for", em("yanwinterwheat"), "dataset"),
                                   br(), br(),
                                   h4(strong("Robust AMMI models")),
                                   tags$p(" The AMMI model, in its standard form, implicitly
                                      assumes equal weights for all entries of the two-way dataset and
                                      that no outliers (leverage points) are present in the data. Since, as in many other real-life
                                      studies the distribution of these data is usually not normal due to the presence of outlying observations,
                                      either resulting from measurement errors or sometimes from individual intrinsic characteristics,
                                      robust SVD methods were proposed by Rodrigues et al. (2016) porposed five robust AMMI models: ", align = "justify"),
                                   tags$ol(
                                     tags$li("rAMMI"),
                                     tags$li("hAMMI"),
                                     tags$li("gAMMI"),
                                     tags$li("lAMMI"),
                                     tags$li("ppAMMI")
                                   ),
                                   br(),
                                   h4(strong("Some details")),
                                   tags$p("Since the model requires a single observation for each combination of genotype and environment,
                                          if there are repetitions, the phenotypic average value is automatically calculated
                                          before fitting the model. Missing values are not allowed Missing values are NOT allowed.", align = "justify"),
                                   br(),
                                   tags$img(src="rAMMI_S.png", height="100%", width="100%"), br(),
                                   tags$h5(strong("Figure 19:"),"rAMMI biplot for", em("yanwinterwheat"), "dataset"),

                                   width = 10
                                 )
                        )
                  )
              ),
             tabPanel("Video-Tutorial",
                      mainPanel(
                        tags$div(
                        tags$h4(strong("Video Tutorial of Geneticae Shiny APP")),
                        HTML('<iframe width="900" height="500" src="https://www.youtube.com/embed/WAODnJDIZ-s" frameborder="0" allowfullscreen></iframe>'),
                        style="text-align: center;")
                      )
             )
  )
)



#------------------------
# Definir funci?n server
#------------------------

server <- function(input, output, session) {
  observe_helpers(withMathJax = TRUE)
# Seleccion de archivo .csv
datafile <- callModule(csvFile, "datafile",
                         stringsAsFactors = FALSE)

# Setear el tipo de letra para graficos
# windowsFonts(Times = windowsFont("Arial Unicode MS"))


# Mostrar el conjunto de datos

output$table <- renderDataTable({
dat<-as.data.frame(datafile())
datatable(dat, editable="cell", class = 'cell-border stripe', rownames = F) %>%
      formatStyle(1:ncol(dat),
                  # background = styleColorBar(c(0, 1), 'lightgray'),
                  backgroundSize = '98% 88%',
                  backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'center', color = "black",
                  backgroundColor="white",
                  fontWeight="bold", textAlign = 'center',
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

example_withoutrep <- reactive({
  # Data without replication
  data(yan.winterwheat)
  data_withoutrep <-yan.winterwheat
})

example_withrep <- reactive({
  # Data with replication
  data(plrv)
  data_withrep <- plrv
})

# Show example dataset
observeEvent( input$example_withoutrep,

output$show_example<-renderText({

  knitr::kable(example_withoutrep(), font_size = 12,
               caption = "Dataset example without repetitions") %>%
    kable_styling(position = "center")


})
)

observeEvent(input$example_withrep,
output$show_example<-renderText({
  knitr::kable(example_withrep(), font_size = 12,
               caption = "Dataset example with repetitions") %>%
    kable_styling(position = "center")
})
)

 # Downloadable csv of sample dataset ----

output$downloadexample_withoutrep<- downloadHandler(
  filename = function() {
    paste(input$Filename1, '.csv', sep = '')
  },
  content = function(file) {
    write.csv(example_withoutrep(), file, row.names = FALSE)
  }
)



output$downloadexample_withrep<- downloadHandler(
  filename = function() {
    paste(input$Filename2, '.csv', sep = '')
  },
  content = function(file) {
    write.csv(example_withrep(), file, row.names = FALSE)
  }
)


sinrep<-reactive({
  data<- datafile() %>%
          group_by(!!sym(input$select_gen), !!sym(input$select_env)) %>%
          summarise(mean_resp=mean(!!sym(input$select_pheno)))%>%
          spread(!!sym(input$select_env), mean_resp) %>%
          as.data.frame()

  rownames(data)<-data[,1]
  data<-data[,-1]
})
# -------------------------------------------------------
#    No anda por eso no hacen matrz de corr
# -------------------------------------------------------



dataset <- reactive({
  datos <-datafile() %>%
    mutate(gen=factor(!!sym(input$select_gen)),
           env=factor(!!sym(input$select_env)),
           pheno=!!sym(input$select_pheno))%>%
    as.data.frame()
})


# Descriptive analysis

# Boxplot
boxplotInput <- eventReactive( input$do, {
    if(input$Var=="Environment"){
      boxp_ <-  ggplot(dataset(), aes(x=env, y=pheno))
    } else{
      boxp_ <-  ggplot(dataset(), aes(x=gen, y=pheno))
    }

  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:30) {
                   incProgress(1/20)
                   Sys.sleep(0.30)
                 }
               })

    boxp<-boxp_ +
      geom_boxplot(fill=input$fillcol)+
      labs(y = input$axisy, x=input$axisx) +
      theme_few()+
      theme(text=element_text(family="Times", size=9),
            axis.text=element_text(family="Times", size=9,angle = 90),
            axis.title=element_text(family="Times", size=12))


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

corrInput <- reactive({
  if(input$Var3=="Environment"){
    M<-cor(sinrep(), method =input$correlType)
  }else{
    M<-cor(t(sinrep()), method =input$correlType)
  }
})

# Correlation plot
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
      theme(text=element_text(family="Times", size=9),
            axis.text.x=element_text(family="Times", size=9,angle = 90),
            axis.text.y=element_text(family="Times", size=9),
            axis.title=element_text(family="Times", size=9))
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
  if(input$Var3=="Environment"){
    M<-cor(sinrep(), method =input$correlType)
  }else{
    M<-cor(t(sinrep()), method =input$correlType)
  }

  knitr::kable(corrInput(),"html", digits = 2,full_width = F, font_size = 12) %>%
    kable_styling(position = "center")  %>%
    column_spec(1:ncol(corrInput()), bold = T, color = "black") %>%
    row_spec(1:nrow(corrInput()), bold = T, color = "black")
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

  if(input$Var4=="Environment"){
    intp<-ggplot2::ggplot(data = datos_summ(), aes(x = !!sym(input$select_env), y = y, colour = !!sym(input$select_gen), group = !!sym(input$select_gen)))
  }else{
    intp<-ggplot2::ggplot(data = datos_summ(), aes(x = !!sym(input$select_gen), y = y, colour = !!sym(input$select_env), group = !!sym(input$select_env)))
  }

  intp <- intp +
    stat_summary(fun.y = mean, geom = "point") +
    stat_summary(fun.y = mean,geom = "line")+
    labs(y = input$axisy_int, x=input$axisx_int)+
    theme_few() +
    theme(text=element_text(family="Times", size=9),
          axis.text=element_text(family="Times", size=9),
          axis.title=element_text(family="Times", size=9))

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



AnovaInput <- reactive({
  if(!is.null(input$select_rep)){
    mod<-lm(as.numeric(pheno) ~ as.factor(env) * as.factor(gen), data = dataset())
  }else{
    mod<-lm(as.numeric(pheno) ~ as.factor(env) + as.factor(gen), data = dataset())
  }
})

output$Anova <- renderText({
  withProgress(message = 'Estimation in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:20) {
                   incProgress(1/20)
                   Sys.sleep(0.25)
                 }
               })
  if(!is.null(input$select_rep)){
    anova<-as.data.frame(anova(AnovaInput()))
    row.names(anova)[1:3] = c("Environment", "Genotype", "Interaction")
    options(knitr.kable.NA = '-')
  }else{
    anova<-as.data.frame(anova(AnovaInput()))
    row.names(anova)[1:2] = c("Environment", "Genotype")
    options(knitr.kable.NA = '-')
  }
  knitr::kable(anova, full_width = F, font_size = 16, col.names = c("Df", "Sum Sq", "Mean Sq", "F-Value", "P-
          Value"), caption = 'Two-way ANOVA table') %>%
    kable_styling(position = "center") %>%
    column_spec(1, bold = T, color = "black")
})

output$Anova2 <- renderPrint({

  if(!is.null(input$select_rep)){
    showModal(modalDialog(
      "The interaction effect can be tested since there are repetitions in the data set",
      footer = list(
        actionButton("ok", "OK")
      )
    ))
  }else{
    showModal(modalDialog(
      "The interaction effect can not be tested since there aren´t repetitions in the data set",
      footer = list(
        actionButton("ok", "OK")
      )
    ))
  }

  observeEvent(input$ok,
               removeModal()
  )


})



# Histograma de residuos
histInput <- reactive({
  anova<-aov(AnovaInput())
  residuals<-as.data.frame(anova$residuals)
  histogram <- ggplot(residuals, aes(x=anova$residuals)) +
    geom_histogram(aes(y=..density..), colour="black", fill=input$fillhist) +
    xlab("Residuals") +
    theme_few()

  ggplotly(histogram) %>%
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


output$residual_hist <- renderPlotly({
  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:20) {
                   incProgress(1/20)
                   Sys.sleep(0.30)
                 }
               })
  print(histInput())
  })


output$downloadHist <- downloadHandler(
  filename = function() {
    paste(input$Filename_hist,'.html', sep='')
    },
  content = function(file) {
    htmlwidgets::saveWidget(as_widget(histInput()), file)
}
)




# qq-plot residuos
qqnormInput <- reactive({
  anova<-aov(AnovaInput())
  residuals<-as.data.frame(anova$residuals)
  qqplot <-gg_qq(residuals(anova))+
   geom_point(colour = input$fillqq)



  ggplotly(qqplot) %>%
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


output$residual_qqnorm <- renderPlotly({
  print(qqnormInput())
  })


output$downloadqqnorm <- downloadHandler(
  filename = function() {
    paste(input$Filename_qq,'.html', sep='')
    },
  content = function(file) {
    htmlwidgets::saveWidget(as_widget(histInput()), file)
}
)



# Prueba de normalidad de residuos
observeEvent(input$do_test,{
output$residual_test <- renderText({
  anova<-aov(AnovaInput())
  shapiro<-shapiro.test(residuals(anova))

  if(shapiro$p.value>0.05){
  sendSweetAlert(
    session = session,
    title = "Normality is verify!!",
    text = paste("Shapiro-Wilk statistic:", round(shapiro$statistic,2),"p-value:",round(shapiro$p.value,4),
                 "\n", "5% significance level"),
    type = "success"
  )
  }else{
    sendSweetAlert(
      session = session,
      title = "Normality is not verify",
      text = paste("Shapiro-Wilk statistic:", round(shapiro$statistic,2),"p-value:",round(shapiro$p.value,4),
                   "\n\n", "5% significance level"),
      type = "error"
    )
}
  })
})

# Grafico de linealidad
linearityInput <- reactive({
  autoplot(AnovaInput(),which =1, ncol = 1, label.size = 3,colour = input$fillhomo, smooth.colour = "black") + theme_few() +
    theme(text=element_text(family="Times", size=12),
          axis.text=element_text(family="Times", size=12),
          axis.title=element_text(family="Times", size=12))
})

output$residual_linearity <- renderPlot({
  print(linearityInput())
})



output$downloadlinearity<- downloadHandler(
  filename = function() {
    paste(input$Filename_homo, '.png', sep='')
    },
  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = linearityInput(), device = device)
  }
  )


observeEvent(input$do_leve_gen,{
    levene<-leveneTest(pheno ~ as.factor(gen), center=mean, data = dataset())

    if(levene$"Pr(>F)"[1]>0.05){
      sendSweetAlert(
        session = session,
        title = "Homogeneity of variance between genotype is verify!!",
        text = paste("Levene statistic:", round(levene$"F value"[1],2),"p-value:",round(levene$"Pr(>F)"[1],4),
                     "\n", "5% significance level"),
        type = "success"
      )
    }else{
      sendSweetAlert(
        session = session,
        title = "Homogeneity of variance between genotype is not verify!!",
        text = paste("Levene statistic:", round(levene$"F value"[1],2),"p-value:",round(levene$"Pr(>F)"[1],4),
                     "\n", "5% significance level"),
        type = "error"
      )
    }
})



observeEvent(input$do_leve_env,{
    levene<-leveneTest(pheno ~ as.factor(env), center=mean, data = dataset())

    if(levene$"Pr(>F)"[1]>0.05){
      sendSweetAlert(
        session = session,
        title = "Homogeneity of variance between environments is verify!!",
        text = paste("Levene statistic:", round(levene$"F value"[1],2),"p-value:",round(levene$"Pr(>F)"[1],4),
                     "\n", "5% significance level"),
        type = "success"
      )
    }else{
      sendSweetAlert(
        session = session,
        title = "Homogeneity of variance between environments is verify!!",
        text = paste("Levene statistic:", round(levene$"F value"[1],2),"p-value:",round(levene$"Pr(>F)"[1],4),
                     "\n", "5% significance level"),
        type = "error"
      )
    }
})




# Grafico de outliers
outliersInput <- reactive({
  p <-autoplot(AnovaInput(), which =4, ncol = 1, label.size = 3,
               colour = input$fillout, smooth.colour = "black") + theme_few() +
    theme(text=element_text(family="Times", size=12),
          axis.text=element_text(family="Times", size=12),
          axis.title=element_text(family="Times", size=12))
})

output$residual_outliers <- renderPlot({
  print(outliersInput())
})

output$downloadoutliers <- downloadHandler(
  filename = function() {
    paste(input$Filename_out, '.png', sep='')
    },
  content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = width, height = height,
                     res = 300, units = "in")
    }
    ggsave(file, plot = outliersInput(), device = device)
  }
  )



# GGE Biplot
# output$genot <- renderUI({
#   checkboxGroupInput(inputId = "select_gen_",
#                      label = "Genotypes",
#                      choices = names(datafile()))
# })


modelInput <- reactive({
  geneticae::GGEmodel(datafile(), genotype= input$select_gen, environment= input$select_env, rep = input$select_rep, response = input$select_pheno,
                      SVP = input$SVP, centering = input$center, scaling = input$scale)
})



plotGGEInput <- eventReactive( input$do_GGE, {
  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:25) {
                   incProgress(1/25)
                   Sys.sleep(0.25)
                 }
               })


  if(!is.null(input$select_rep)){

    p <-  geneticae::GGEPlot(modelInput(), type=input$plotType, footnote = input$footnote,
                           colGen =input$colgen, colEnv =input$colenv,colSegment=input$colsegment,
                           titles = input$title, sizeGen=as.numeric(input$sizeGen), sizeEnv = as.numeric(input$sizeEnv),
                           axislabels=input$axislabels,axes=input$axes)
  }else{
    p <-  geneticae::GGEPlot(modelInput(), type=input$plotType, footnote = input$footnote,
                             colGen =input$colgen, colEnv =input$colenv,colSegment=input$colsegment,
                             titles = input$title, sizeGen=as.numeric(input$sizeGen), sizeEnv = as.numeric(input$sizeEnv),
                             axislabels=input$axislabels,axes=input$axes)
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




# AMMI biplot
plotAMMIInput <- eventReactive( input$do_AMMI, {
  withProgress(message = 'Graphic in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:20) {
                   incProgress(1/20)
                   Sys.sleep(0.30)
                 }
               })
  if(!is.null(input$select_rep)){
  p <-  geneticae::rAMMI(dataset(), genotype= input$select_gen, environment= input$select_env, rep = input$select_rep, response = input$select_pheno,
                         type=input$robustType, footnote = input$foot, colGen =input$colorgen, colEnv =input$colorenv,
                         titles = input$tit, sizeGen=as.numeric(input$sizeGen_AMMI), sizeEnv = as.numeric(input$sizeEnv_AMMI),
                         axislabels=input$axislabels_AMMI,axes=input$axes_AMMI)
  }else{
  p <-  geneticae::rAMMI(dataset(), genotype= input$select_gen, environment= input$select_env, rep = NULL, response = input$select_pheno,
                         type=input$robustType, footnote = input$foot, colGen =input$colorgen, colEnv =input$colorenv,
                         titles = input$tit, sizeGen=as.numeric(input$sizeGen_AMMI), sizeEnv = as.numeric(input$sizeEnv_AMMI),
                         axislabels=input$axislabels_AMMI,axes=input$axes_AMMI)

  }

  })


output$plotAMMI <- renderPlot({

  print(plotAMMIInput())
})

output$download_ammi<- downloadHandler(
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




