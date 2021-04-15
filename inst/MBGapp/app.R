#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(geoR)
library(ggplot2)
library(magrittr)
library(dplyr)
library(readr)
library(tidyr)
library(tmap)
library(sf)
library(leaflet)
library(rgdal)
library(shinyjs)
require(PrevMap)
require(splancs)
require(grDevices)

options(shiny.maxRequestSize = 30*1024^2)
# jsCode <- "shinyjs.hideSidebar = function(params){$('body').addClass('sidebar-collapse');}"


########### useful functions to deal with variogram ###############
variog_envelope <- function (geodata, coords = geodata$coords, data = geodata$data, 
                             obj.variog, nsim = 999, save.sim = FALSE, messages) 
{
    call.fc <- match.call()
    if (missing(geodata)) 
        geodata <- list(coords = coords, data = data)
    if (missing(messages)) 
        messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                             TRUE, getOption("geoR.messages")))
    else messages.screen <- messages
    obj.variog$v <- NULL
    if ((is.matrix(data) | is.data.frame(data))) 
        if (ncol(data) > 1) 
            stop("envelops can be computed for only one data set at once")
    if (!is.null(obj.variog$estimator.type)) 
        estimator.type <- obj.variog$estimator.type
    else estimator.type <- "classical"
    if (abs(obj.variog$lambda - 1) > 1e-04) {
        if (abs(obj.variog$lambda) < 1e-04) 
            data <- log(data)
        else data <- ((data^obj.variog$lambda) - 1)/obj.variog$lambda
    }
    xmat <- unclass(trend.spatial(trend = obj.variog$trend, geodata = geodata))
    if (obj.variog$trend != "cte") {
        if (is.vector(data)) {
            data <- lm(data ~ xmat + 0)$residuals
            names(data) <- NULL
        }
        else {
            only.res <- function(y, x) {
                lm(y ~ xmat + 0)$residuals
            }
            data <- apply(data, 2, only.res, x = xmat)
        }
    }
    if (messages.screen) 
        cat(paste("variog.env: generating", nsim, "simulations by permutating data values\n"))
    simula <- list(coords = coords)
    n.data <- length(data)
    perm.f <- function(i, data, n.data) {
        return(data[sample(1:n.data)])
    }
    simula$data <- apply(as.matrix(1:nsim), 1, perm.f, data = data, 
                         n.data = n.data)
    if (messages.screen) 
        cat(paste("variog.env: computing the empirical variogram for the", 
                  nsim, "simulations\n"))
    nbins <- length(obj.variog$bins.lim) - 1
    if (obj.variog$direction == "omnidirectional") {
        bin.f <- function(sim) {
            cbin <- vbin <- sdbin <- rep(0, nbins)
            temp <- .C("binit", as.integer(obj.variog$n.data), 
                       as.double(as.vector(coords[, 1])), as.double(as.vector(coords[, 
                                                                                     2])), as.double(as.vector(sim)), as.integer(nbins), 
                       as.double(as.vector(obj.variog$bins.lim)), as.integer(estimator.type == 
                                                                                 "modulus"), as.double(max(obj.variog$u)), as.double(cbin), 
                       vbin = as.double(vbin), as.integer(FALSE), as.double(sdbin), 
                       PACKAGE = "geoR")$vbin
            return(temp)
        }
        simula.bins <- apply(simula$data, 2, bin.f)
    }
    else {
        variog.vbin <- function(x, ...) {
            variog(geodata = geodata, 
                   data = x, uvec = obj.variog$uvec, estimator.type = obj.variog$estimator.type, 
                   nugget.tolerance = obj.variog$nugget.tolerance, max.dist = obj.variog$max.dist, 
                   pairs.min = obj.variog$pairs.min, direction = obj.variog$direction, 
                   tolerance = obj.variog$tolerance, messages.screen = FALSE,...)$v
        }
        simula.bins <- apply(simula$data, 2, variog.vbin)
    }
    simula.bins <- simula.bins[obj.variog$ind.bin, ]
    if (save.sim == FALSE) 
        simula$data <- NULL
    if (messages.screen) 
        cat("variog.env: computing the envelops\n")
    limits <- apply(simula.bins, 1, quantile, prob = c(0.025, 0.975))
    res.env <- list(u = obj.variog$u, v.lower = limits[1, ], 
                    v.upper = limits[2, ])
    if (save.sim) 
        res.env$simulations <- simula$data
    res.env$call <- call.fc
    oldClass(res.env) <- "variogram.envelope"
    return(res.env)
}



thr.var <- function (x, max.dist, scaled = FALSE, ...) 
{
    my.l <- list()
    if (missing(max.dist)) {
        my.l$max.dist <- x$max.dist
        if (is.null(my.l$max.dist)) 
            stop("argument max.dist needed for this object")
    }
    else my.l$max.dist <- max.dist
    if (any(x$cov.model == c("matern", "powered.exponential", 
                             "cauchy", "gencauchy", "gneiting.matern"))) 
        my.l$kappa <- x$kappa
    else kappa <- NULL
    if (is.vector(x$cov.pars)) 
        my.l$sill.total <- x$nugget + x$cov.pars[1]
    else my.l$sill.total <- x$nugget + sum(x$cov.pars[, 1])
    my.l$nugget <- x$nugget
    my.l$cov.pars <- x$cov.pars
    my.l$cov.model <- x$cov.model
    if (scaled) {
        if (is.vector(x$cov.model)) 
            my.l$cov.pars[1] <- my.l$cov.pars[1]/my.l$sill.total
        else my.l$cov.pars[, 1] <- my.l$cov.cov.pars[, 1]/my.l$sill.total
        my.l$sill.total <- 1
    }
    gamma.f <- function(x, my.l) {
        if (any(my.l$cov.model == c("linear", "power"))) 
            return(my.l$nugget + my.l$cov.pars[1] * (x^my.l$cov.pars[2]))
        else return(my.l$sill.total - cov.spatial(x, cov.model = my.l$cov.model, 
                                                  kappa = my.l$kappa, cov.pars = my.l$cov.pars))
    }
    dd <- gamma.f(x= seq(0, my.l$max.dist, length.out = 101), my.l = my.l)
    # curve(gamma.f(x, my.l = my.l), from = 0, to = my.l$max.dist, 
    #       add = TRUE, ...)
    return(dd)
}


# Calculate and plot the variogram
ggvario <- function(coords, 
                    data, 
                    bins = 15, 
                    maxdist = max(dist(coords))/3, 
                    uvec = NULL, 
                    nsim = 999,
                    color = "royalblue1", 
                    xlab = "distance", 
                    show_nbins = F, envelop=1, cov.model="matern", fix.kappa=T) {
    require(geoR)
    res <- list()
    coords <- as.matrix(coords)
    min_dist <- min(dist(coords))
    if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
    empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
    if(envelop ==1){
        
        ### plot variogram alone 
        
        dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                              nbinns = empvario$n)
        p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbinns)) +
            geom_point(col = "black", fill = color, shape = 21, size = 3) +
            scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                               breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
            scale_y_continuous(name = "semivariance",
                               #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1),
                               limits = c(0,  max(dfvario$empirical))) +
            ggtitle("Empirical semivariogram") 
        # theme_classic()
        p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
        if(show_nbins){
            res[["pl"]] <- p2 
        } else {
            res[["pl"]] <- p1
        }
    } else if (envelop == 2){
        
        ## plot variogram  and the therectical variogram line 
        dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                              nbins = empvario$n)
        
        vari.fit <- variofit(vario = empvario, ini.cov.pars=c(mean(dfvario$empirical), variog_extent/3), cov.model = cov.model,
                             fix.nug=F, nugget = 0, fix.kappa = fix.kappa)
        
        p1 <- ggplot() +
            geom_point(data= dfvario, aes(y = empirical, x = distance, label = nbins), col = "black", fill = color, shape = 21, size = 3) +
            scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                               breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
            scale_y_continuous(name = "semivariance",
                               #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1),
                               limits = c(0,  max(dfvario$empirical))) +
            ggtitle("Empirical semivariogram") +
            geom_line(data = data.frame(xx= seq(0, maxdist, length.out = 101), yy = thr.var(vari.fit)), aes(x = xx, y = yy))
        # theme_classic()
        p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
        if(show_nbins){
            res[["pl"]] <- p2 
        } else {
            res[["pl"]] <- p1
        }
        res[["summ"]] <- vari.fit
    }else{
        ### plot the Monte Carlo envelope and variogram
        
        envmc <- variog_envelope(coords = coords, data = data, 
                                 obj.variog = empvario, nsim = nsim, messages = F)
        dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                              lowemp = envmc$v.lower, upemp = envmc$v.upper, 
                              nbins = empvario$n)
        p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
            geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = .3) +
            geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
            scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                               breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
            scale_y_continuous(name = "semivariance", 
                               #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1), 
                               limits = c(0, max(dfvario$upemp, dfvario$empirical))) +
            ggtitle("Empirical semivariogram") 
        # theme_classic()
        p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
        if(show_nbins){
            res[["pl"]] <- p2 
        } else {
            res[["pl"]] <- p1
        }
    }
    res
}




# emplogit<-function(p,N){
#     top=p*N+0.5
#     bottom=N*(1-p)+0.5
#     return(log(top/bottom))
# }


lonlat2UTM = function(lonlat) {
    utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
    if(lonlat[2] > 0) {
        utm + 32600
    } else{
        utm + 32700
    }
}


### convert to metres using 3857


create_labels <- function(x, greater = F, smaller = F) {
    n <- length(x)
    x <- gsub(" ", "", format(x))
    labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
    if (greater) {
        labs[length(labs)] <- paste("\u2265", x[n - 1])
    }
    if (smaller) {
        labs[1] <- paste("<", x[2])
    }
    
    return(labs)
}


# Convert epsg to epsg KM
epsgKM <- function(x) {
    crs <- st_crs(x)
    proj4KM <- gsub(pattern = "+.units=m", replacement = "+units=km", 
                    crs$proj4string)
    return(proj4KM)
}

############################## The begining of the APP ########################################################

# Define UI for application that draws a histogram
ui <- fluidPage(
    useShinyjs(),
    # extendShinyjs(text = jsCode, functions = c("hideSidebar")),
    # img(src='chicas_logo.png', align = "right"),
    # # Application title
    titlePanel(title=div(img(src="chicas_logo.png", align = "right", height = 30, width = 100), "Model-based geostatistics")),
    
    
    # Sidebar with a slider input the data and the shapefile
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "mbgdata", label = "Upload the data (csv file):"),
            radioButtons("maptype", "Do you know the projection of the location?:", 
                         c("Yes" = "view",
                           "No" = "plot"), selected = "plot"),
            conditionalPanel(condition = "input.maptype=='view'",
                             numericInput("crs", "Coordinate reference system (default = 4326):", 4326, min = 1, max = 100000)
                             ),
            fileInput(inputId = "mbgshp", label = "Upload the shapefile (optional):", 
                      accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj"), multiple=TRUE),   
            selectInput("datatype", 'Choose the data type', 
                        choices=c("Continuous data" ='continuous', "Prevalence data" = 'prevalence', "Count data" = 'count'), 
                        selected = NULL),  
            
            # radioButtons("maptype", "Map mode:", 
            #              c("Interactive viewing" = "view",
            #                "Static plotting" = "plot")),
            
            selectInput(
                inputId = "xaxis",
                label = "X-coordinate",
                choices = "",
            ),
            selectInput(
                inputId = "yaxis",
                label = "Y-coordinate",
                choices = "",
            ),
            conditionalPanel(condition = "input.datatype=='continuous'",
                             selectInput(
                                 inputId = "y",
                                 label = "Continuous outcome",
                                 choices = ""
                             )
            ),
            conditionalPanel(condition = "input.datatype=='prevalence'",
                             selectInput(
                                 inputId = "p",
                                 label = "Postives",
                                 choices = ""
                             ),
                             selectInput(
                                 inputId = "m",
                                 label = "Total Examined",
                                 choices = ""
                             )
            ),
            conditionalPanel(condition = "input.datatype=='count'",
                             selectInput(
                                 inputId = "c",
                                 label = "Count",
                                 choices = ""
                             ),
                             selectInput(
                                 inputId = "e",
                                 label = "Offset",
                                 choices = ""
                             )
                             
            ),
            selectInput(
                inputId = "D",
                label = "Covariate(s)",
                choices = "",
                multiple = T
            ),
            
            conditionalPanel(condition = "input.tabselected==1",
                             conditionalPanel(condition = "input.datatype=='continuous'",
                                              radioButtons("transformcont", "Choose outcome transformation", c("No-transform" = "identity", 
                                                                                                               "Log-transform" = "log"))
                                              
                                              
                             ),
                             conditionalPanel(condition = "input.datatype=='prevalence'",
                                              radioButtons("transformprev", "Choose outcome transformation", c("No-transform" = "identity", 
                                                                                                               "Log-transform" = "log", 
                                                                                                               "Logit-transform"= "logit"))
                                              
                                              
                             ),
                             conditionalPanel(condition = "input.datatype=='count'",
                                              radioButtons("transformcnt", "Choose outcome transformation", c("No-transform" = "identity", 
                                                                                                              "Log-transform" = "log"))
                                              
                                              
                             ),
                             radioButtons("transformcov", "Choose covariate transformation ", c("No-transform" = "identity",
                                                                                                "Log-transform" = "log", 
                                                                                                "square-root transform"= "sqrt"))
                             
                             
            ),
            conditionalPanel(condition = "input.tabselected==2",
                             
                             sliderInput(inputId = "nbins",
                                         label = "Number of bins:",
                                         min = 0,
                                         max = 50,
                                         value = 15, step=1), 
                             
                             sliderInput(inputId = "dist",
                                         label = "Distance:",
                                         min = 0,
                                         max = 100,
                                         value = 70, step=1), 
                             
                             # actionButton("change", "Change slider max value"),
                             selectInput("functions", 'Correlation functions', 
                                         choices=c("matern" = "matern", "exponential" = "exponential", "gaussian" = "gaussian", 
                                                   "spherical" = "spherical", "circular" = "circular", 
                                                   "cubic" = "cubic", "wave" ="wave", 
                                                   "powered.exponential" = "powered.exponential", "cauchy" = "cauchy", 
                                                   "gneiting" = "gneiting", 
                                                   "pure.nugget" = "pure.nugget"), 
                                         selected = NULL),  
                             radioButtons("envelop", "Choose plot", c("Variogram only" = "vario", 
                                                                      "Variogram with envelope"= "varioEnve",
                                                                      "Fit theorectical variogram" = "varifit")),
                             conditionalPanel(condition = "input.envelop =='varioEnve'", 
                                              numericInput("npermute", "Number of permutation", 999)),
                             
                             shiny::actionButton(inputId='ab1', label="Learn More",
                                                 icon = icon("th"),
                                                 onclick ="window.open('https://olatunjijohnson.shinyapps.io/variogshiny', '_blank')"),
                             
                             # tags$a(href="https://olatunjijohnson.shinyapps.io/variogshiny/", "Learn More!", target="_blank"),
                             
                             
                             #### This part helps to hide the error 
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             )
                             
                             
            ),
            conditionalPanel(condition = "input.tabselected==3",
                             conditionalPanel(condition = "input.datatype =='prevalence'",
                                              radioButtons("fitlinear", "Do you want to fit linear model?:", 
                                                           c("Yes" = "linearmodel",
                                                             "No" = "binomialmodel"), selected = "binomialmodel")),
                             
                             numericInput("phi", "Intial value of scale parameter", 50),
                             numericInput("nu", "Intial value of relative variance of the nugget effect", 0.1),
                             numericInput("kappa", "Value of kappa", 0.5),
                             conditionalPanel(condition = "input.datatype !='continuous' & input.fitlinear=='binomialmodel'",
                                              actionButton("AdvOption", "Advance options")),
                             conditionalPanel(condition = "input.AdvOption & input.fitlinear=='binomialmodel'", 
                                              numericInput("mcmcNsim", "Number of simulation", 1000),
                                              numericInput("mcmcNburn", "Number of burn-in", 200),
                                              numericInput("mcmcNthin", "Number of thinning", 8)
                                              ),
                             actionButton("ShowEst", "Show the result summary", icon = icon("fas fa-running")),
                             actionButton("gotab", "show table"),
                             
                             #### This part helps to hide the error 
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             )
                             
                             
            ),
            conditionalPanel(condition = "input.tabselected==4",
                             numericInput("resolution", label = "Spatial resolution", value = 10, min = 0, max = 100000),
                             fileInput(inputId = "gridpreddata", label = "Upload the predictive grid (optional):"),
                             fileInput(inputId = "predictorsdata", label = "Upload the predictors"),
                             
                             conditionalPanel(condition = "input.datatype=='continuous'",
                                              radioButtons(inputId = "predtomapcont", label = "Choose map", 
                                                           choices = c("Mean Outcome" = "meann", 
                                                                       "Standard error" = "sdd",
                                                                       "Exceedance probability" = "exprob",
                                                                       "Quantile" = "quant"))
                             ),
                             conditionalPanel(condition = "input.datatype=='prevalence'",
                                              radioButtons(inputId = "predtomapprev", label = "Choose map", 
                                                           choices = c("Mean prevalence" = "meann", 
                                                                       "Standard error" = "sdd",
                                                                       "Exceedance probability" = "exprob",
                                                                       "Quantile" = "quant"))
                                              
                             ),
                             conditionalPanel(condition = "input.datatype=='count'",
                                              radioButtons("predtomapcount", "Choose map", 
                                                           choices = c("Mean" = "meann", 
                                                                       "Standard error"= "sdd",
                                                                       "Exceedance probability" = "exprob",
                                                                       "Quantile" = "quant"))
                                              
                                              
                             ),
                             conditionalPanel(condition = "input.predtomapcount=='exprob' | input.predtomapprev=='exprob' | input.predtomapcont=='exprob'",
                                              sliderInput(inputId = "threshold",
                                                          label = "Exceedance probability threshold:",
                                                          min = 0,
                                                          max = 1,
                                                          value = 0.5, step=0.01,
                                                          animate=animationOptions(interval = 1000,loop=F))
                             ),
                             conditionalPanel(condition = "input.predtomapcount=='quant' | input.predtomapprev=='quant' | input.predtomapcont=='quant'",
                                              sliderInput(inputId = "quantprob",
                                                          label = "Probability:",
                                                          min = 0,
                                                          max = 1,
                                                          value = 0.5, step=0.05)
                             ),
                             
                             actionButton("ShowPred", "Map the prediction", icon = icon("fas fa-running")),
                             
                             #### This part helps to hide the error 
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             )
                             
                             
            ),
            conditionalPanel(condition = "input.tabselected==5"),
        ),
        # Show a map and plot  of the data
        mainPanel(
            
            tabsetPanel(type="pills", 
                        tabPanel("Explore", value = 1,       
                                 conditionalPanel(condition = "input.maptype == 'view'",
                                                  leafletOutput(outputId = "map")),
                                 conditionalPanel(condition = "input.maptype == 'plot'",
                                                  plotOutput(outputId = "map2")),
                                 # h3("Scatter plot of the outcome and the covariate"),
                                 plotOutput(outputId ="Plot")),
                        tabPanel("Variogram", value = 2,       
                                 plotOutput(outputId ="variogplot"),
                                 # h3("Summary of estimate covariance parameter"),
                                 verbatimTextOutput(outputId ="summary")),
                        tabPanel("Estimation", value = 3,       
                                 verbatimTextOutput(outputId ="estsummary"),
                                 htmlOutput("tab")), 
                        tabPanel("Prediction", value = 4,    
                                 conditionalPanel(condition = "input.maptype == 'view'",
                                                  leafletOutput(outputId = "predmap", height=800)),
                                 conditionalPanel(condition = "input.maptype == 'plot'",
                                                  plotOutput(outputId = "predmap2", height=800))),
                        
                        tabPanel("Report", value = 5,
                                 conditionalPanel(condition = "input.maptype == 'view'",
                                                  downloadButton("report2", "Download report")),
                                 conditionalPanel(condition = "input.maptype == 'plot'",
                                                  downloadButton("report", "Download report")),
                                 
                                 
                                 HTML("<br>"),
                                 helpText("Download the report."),
                                 checkboxGroupInput("whattoshow", "What to show in the report", inline=F,
                                                    c("Map of the outcome" = "fig1",
                                                      "Scatter plot of outcome and covariate" = "fig2",
                                                      "Variogram plot" = "fig3",
                                                      "Summary of the parameter" = "fig4",
                                                      "Prediction map" = "fig5"),
                                 selected = c("fig1")),
                                 
                                 # checkboxGroupInput("rtables", "Tables summary", inline=TRUE,
                                 #                    c("Population" = "Population", "Observed" = "Observed", "Expected" = "Expected", "SIR" = "SIR",
                                 #                      "Risk" = "Risk", "2.5 percentile" = "LowerLimitCI", "97.5 percentile" = "UpperLimitCI"),
                                 #                    selected = c("Population", "Observed", "Expected", "SIR", "Risk", "LowerLimitCI", "UpperLimitCI")),
                                 
                                 
                                 HTML("<br>")),
                        
                        id="tabselected"
            )
            
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    
    ##### hide some sidebars
    observeEvent(input$tabselected, {
        if(input$tabselected == 1){
            shinyjs::show(id = "mbgdata")
            shinyjs::show(id = "crs")
            shinyjs::show(id = "mbgshp")
            shinyjs::show(id = "datatype")
            shinyjs::show(id = "maptype")
            shinyjs::show(id = "xaxis")
            shinyjs::show(id = "yaxis")
            shinyjs::show(id = "y")
            shinyjs::show(id = "p")
            shinyjs::show(id = "m")
            shinyjs::show(id = "c")
            shinyjs::show(id = "e")
            shinyjs::show(id = "D")
        }else if (input$tabselected == 2){
            shinyjs::hide(id = "mbgdata")
            shinyjs::hide(id = "crs")
            shinyjs::hide(id = "mbgshp")
            shinyjs::hide(id = "datatype")
            shinyjs::hide(id = "maptype")
            shinyjs::show(id = "xaxis")
            shinyjs::show(id = "yaxis")
            shinyjs::show(id = "y")
            shinyjs::show(id = "p")
            shinyjs::show(id = "m")
            shinyjs::show(id = "c")
            shinyjs::show(id = "e")
            shinyjs::show(id = "D")
        }else if (input$tabselected == 3){
            shinyjs::hide(id = "mbgdata")
            shinyjs::hide(id = "crs")
            shinyjs::hide(id = "mbgshp")
            shinyjs::hide(id = "datatype")
            shinyjs::hide(id = "maptype")
            shinyjs::show(id = "xaxis")
            shinyjs::show(id = "yaxis")
            shinyjs::show(id = "y")
            shinyjs::show(id = "p")
            shinyjs::show(id = "m")
            shinyjs::show(id = "c")
            shinyjs::show(id = "e")
            shinyjs::show(id = "D")
        }else if (input$tabselected == 4){
            shinyjs::hide(id = "mbgdata")
            shinyjs::hide(id = "crs")
            shinyjs::hide(id = "mbgshp")
            shinyjs::hide(id = "datatype")
            shinyjs::hide(id = "maptype")
            shinyjs::hide(id = "xaxis")
            shinyjs::hide(id = "yaxis")
            shinyjs::hide(id = "y")
            shinyjs::hide(id = "p")
            shinyjs::hide(id = "m")
            shinyjs::hide(id = "c")
            shinyjs::hide(id = "e")
            shinyjs::hide(id = "D")
        }else if (input$tabselected == 5){
            shinyjs::hide(id = "mbgdata")
            shinyjs::hide(id = "crs")
            shinyjs::hide(id = "mbgshp")
            shinyjs::hide(id = "datatype")
            shinyjs::hide(id = "maptype")
            shinyjs::hide(id = "xaxis")
            shinyjs::hide(id = "yaxis")
            shinyjs::hide(id = "y")
            shinyjs::hide(id = "p")
            shinyjs::hide(id = "m")
            shinyjs::hide(id = "c")
            shinyjs::hide(id = "e")
            shinyjs::hide(id = "D")
        }
    })
    # Upload the data
    data_all <- reactive({
        req(input$mbgdata)
        dff <- input$mbgdata
        if (is.null(dff))
            return(NULL)
        if(grepl("\\.rds$", dff)){
            x <- readRDS(dff$datapath)
            x
        }else{
            x <- read_csv(dff$datapath)
            x$XXX <- 1
            x$YYY <- 1
            if(input$datatype== "prevalence") x$emplogit <- 1
            x
        }
    })
    
    
    # Upload the shapefile  
    map_all <- reactive({
        shpdf <- input$mbgshp
        if(is.null(shpdf)){
            return()
        }
        previouswd <- getwd()
        uploaddirectory <- dirname(shpdf$datapath[1])
        setwd(uploaddirectory)
        for(i in 1:nrow(shpdf)){
            file.rename(shpdf$datapath[i], shpdf$name[i])
        }
        setwd(previouswd)
        
        #map <- readShapePoly(paste(uploaddirectory, shpdf$name[grep(pattern="*.shp", shpdf$name)], sep="/"),  delete_null_obj=TRUE)
        #reads the file that finishes with .shp using $ at the end: grep(pattern="*.shp$", shpdf$name)
        map <- readOGR(paste(uploaddirectory, shpdf$name[grep(pattern="*.shp$", shpdf$name)], sep="/"))#,  delete_null_obj=TRUE)
        map <- st_transform(st_as_sf(map), crs=4326)
        # map <- st_as_sf(spTransform(map, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))
        
        map
        
    })
    
    
    # Upload the grid.locations 
    gridpred <- reactive({
        req(input$gridpreddata)
        dff <- input$gridpreddata
        if (is.null(dff))
            return(NULL)
        if(grepl("\\.rds$", dff)){
            x <- readRDS(dff$datapath)
            x
        }else{
            x <- read_csv(dff$datapath)
            x
        }
    })
    
    # Upload the predictors
    predictors <- reactive({
        req(input$predictorsdata)
        dff <- input$predictorsdata
        if (is.null(dff))
            return(NULL)
        if(grepl("\\.rds$", dff)){
            x <- readRDS(dff$datapath)
            x
        }else{
            x <- read_csv(dff$datapath)
            x
        }
    })
    
    
    
    # Update the choices when the data is uploaded
    # note that I can update the label in the below
    observe({
        df2 <- data_all()
        updateVarSelectInput(session, "xaxis",  data = df2)
        updateVarSelectInput(session, "yaxis",  data = df2)
        updateVarSelectInput(session, "y", data = df2)
        updateVarSelectInput(session, "D",  data = df2)
        updateVarSelectInput(session, "p", data = df2)
        updateVarSelectInput(session, "m", data = df2)
        updateVarSelectInput(session, "c", data = df2)
        updateVarSelectInput(session, "e", data = df2)
    })
    
    # Change the maximum distance of the variogram 
    
    # observeEvent(input$change,{
    #     df2 <- data_all()
    #     dummy_coords <- data.frame(df2[, c(input$xaxis,input$yaxis)])
    #     dummy_coords <- dummy_coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>%
    #         st_transform(., crs=epsgKM(as.numeric(lonlat2UTM(dummy_coords[1,]))))
    #     variog_extent <- max(dist(cbind(st_coordinates(dummy_coords))), na.rm=T)
    #     updateSliderInput(session, "dist", min = 0, max = variog_extent,
    #                       value = variog_extent/3,
    #                       step = round(variog_extent/100)+1)
    # })
    
    observeEvent(input$tabselected == 2, {
        df2 <- data_all()
        dummy_coords <- data.frame(df2[, c(input$xaxis,input$yaxis)])
        if(input$maptype == 'plot'){
            dummy_coords <- dummy_coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
            variog_extent <<- max(dist(cbind(st_coordinates(dummy_coords))), na.rm=T)
            updateSliderInput(session, "dist", min = 0, max = variog_extent,
                              value = variog_extent/3,
                              step = round(variog_extent/100)+1)
        }else{
            dummy_coords <- dummy_coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>%
                st_transform(., crs=epsgKM(as.numeric(lonlat2UTM(dummy_coords[1,]))))
            variog_extent <<- max(dist(cbind(st_coordinates(dummy_coords))), na.rm=T)
            updateSliderInput(session, "dist", label = "Distance (km)", min = 0, max = variog_extent,
                              value = variog_extent/3,
                              step = round(variog_extent/100)+1)
        }
        
    })
    
    
    
    # observeEvent(input$change,{
    #   updateSliderInput(session, "dist", max = 50000, step = round(50000/50))
    # })
    
    flname <- reactive({
        if(input$maptype == 'view'){
            filename <- "report.html"
            filename
        }else{
            filename <- "report.pdf"
            filename
        }
    })
    
    
    explore_map_lf <- reactive({
        df <- data_all()
        if(input$datatype=='continuous'){
            if(is.null(input$mbgshp)){
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=input$crs)
                mapdata <- st_transform(mapdata, crs=4326)
                # brks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
                # labs <- create_labels(brks, greater = F)
                # pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(mapdata) +
                        tm_symbols(col=input$y, style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1) +
                        # tm_shape(shp, is.master = T) +
                        # tm_borders(col="black") +
                        tm_layout()
                )
                l
                # l <- tmap::tm_shape(mapdata) +
                #         tm_symbols(col=input$y, style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1) +
                #         # tm_shape(shp, is.master = T) +
                #         # tm_borders(col="black") +
                #         tm_layout()
                # l
            }else{
                shp <- map_all()
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=crs(shp))
                mapdata <- st_transform(mapdata, crs=4326)
                
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(mapdata) +
                        tm_symbols(col=input$y, style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1) +
                        tm_shape(shp, is.master = T) +
                        tm_borders(col="black") +
                        tm_layout()
                )
                l
            }
        }else if (input$datatype=='prevalence'){
            if(is.null(input$mbgshp)){
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=input$crs)
                mapdata <- st_transform(mapdata, crs=4326)
                new_dat <- data.frame(df[, c(input$p, input$m, input$D), drop=FALSE])
                
                
                # mapdata["Emplogit"] <- log((df[,input$p] + 0.5)/(df[, input$m] - df[,input$p] + 0.5))
                mapdata[,"Prevalence"] <- df[,input$p]/df[, input$m] 
                # brks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
                # labs <- create_labels(brks, greater = F)
                # pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(mapdata) +
                        tm_symbols(col="Prevalence", style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1, 
                                   title.col = "Empirical prevalence") +
                        # tm_shape(shp, is.master = T) +
                        # tm_borders(col="black") +
                        tm_layout()
                )
                l
            }else{
                shp <- map_all()
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=crs(shp))
                mapdata <- st_transform(mapdata, crs=4326)
                mapdata[,"Prevalence"] <- df[,input$p]/df[, input$m] 
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(mapdata) +
                        tm_symbols(col="Prevalence", style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1, 
                                   title.col = "Empirical prevalence") +
                        tm_shape(shp, is.master = T) +
                        tm_borders(col="black") +
                        tm_layout()
                )
                l
            }
        }else{
            if(is.null(input$mbgshp)){
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=input$crs)
                mapdata <- st_transform(mapdata, crs=4326)
                # brks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
                # labs <- create_labels(brks, greater = F)
                # pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                mapdata[,"incidence"] <- df[,input$c]/df[, input$e]
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(mapdata) +
                        tm_symbols(col="incidence", style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1,
                                   title.col = "Incidence") +
                        # tm_shape(shp, is.master = T) +
                        # tm_borders(col="black") +
                        tm_layout()
                )
                l
            }else{
                shp <- map_all()
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=crs(shp))
                mapdata <- st_transform(mapdata, crs=4326)
                mapdata[,"incidence"] <- df[,input$c]/df[, input$e]
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(mapdata) +
                        tm_symbols(col="incidence", style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1,
                                   title.col = "Incidence") +
                        tm_shape(shp, is.master = T) +
                        tm_borders(col="black") +
                        tm_layout()
                )
                l
            }
        }
    })
    
    
    output$map <- renderLeaflet({
        if (is.null(explore_map_lf())) return(NULL)
        explore_map_lf()
    })
    
    ################### mapping without interactive map ###################################
    
    explore_map_st <- reactive({
        df <- data_all()
        if(input$datatype=='continuous'){
            if(is.null(input$mbgshp)){
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
                # mapdata <- st_transform(mapdata)
                # brks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
                # labs <- create_labels(brks, greater = F)
                # pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tm_shape(mapdata) +
                    tm_symbols(col=input$y, style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1) +
                    # tm_shape(shp, is.master = T) +
                    # tm_borders(col="black") +
                    tm_layout()
                l
            }else{
                shp <- map_all()
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
                # mapdata <- st_transform(mapdata)
                
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tm_shape(mapdata) +
                    tm_symbols(col=input$y, style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1) +
                    tm_shape(shp, is.master = T) +
                    tm_borders(col="black") +
                    tm_layout()
                l
            }
        }else if (input$datatype=='prevalence'){
            if(is.null(input$mbgshp)){
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
                # mapdata <- st_transform(mapdata)
                new_dat <- data.frame(df[, c(input$p, input$m, input$D), drop=FALSE])
                
                
                # mapdata["Emplogit"] <- log((df[,input$p] + 0.5)/(df[, input$m] - df[,input$p] + 0.5))
                mapdata[,"Prevalence"] <- df[,input$p]/df[, input$m] 
                # brks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
                # labs <- create_labels(brks, greater = F)
                # pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tm_shape(mapdata) +
                    tm_symbols(col="Prevalence", style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1, 
                               title.col = "Empirical prevalence") +
                    # tm_shape(shp, is.master = T) +
                    # tm_borders(col="black") +
                    tm_layout()
                l
            }else{
                shp <- map_all()
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
                # mapdata <- st_transform(mapdata)
                mapdata[,"Prevalence"] <- df[,input$p]/df[, input$m] 
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tm_shape(mapdata) +
                    tm_symbols(col="Prevalence", style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1, 
                               title.col = "Empirical prevalence") +
                    tm_shape(shp, is.master = T) +
                    tm_borders(col="black") +
                    tm_layout()
                l
            }
        }else{
            if(is.null(input$mbgshp)){
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
                # mapdata <- st_transform(mapdata)
                # brks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
                # labs <- create_labels(brks, greater = F)
                # pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                mapdata[,"incidence"] <- df[,input$c]/df[, input$e]
                l <- tmap::tm_shape(mapdata) +
                    tm_symbols(col="incidence", style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1,
                               title.col = "Incidence") +
                    # tm_shape(shp, is.master = T) +
                    # tm_borders(col="black") +
                    tm_layout()
                l
            }else{
                shp <- map_all()
                mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
                # mapdata <- st_transform(mapdata, crs=4326)
                mapdata[,"incidence"] <- df[,input$c]/df[, input$e]
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = 5, contrast = c(0, 1), plot = F)
                l <- tmap::tm_shape(mapdata) +
                    tm_symbols(col="incidence", style="equal", alpha=0.5, size=0.2, palette="-RdYlBu", contrast=1,
                               title.col = "Incidence") +
                    tm_shape(shp, is.master = T) +
                    tm_borders(col="black") +
                    tm_layout()
                l
            }
        }
    })
    
    output$map2 <- renderPlot({
        if (is.null(explore_map_st())) return(NULL)
        explore_map_st()
    })
    
    
    ######### Plotting the relationship between the outcomes and the covariates 
    
    scatter_ass_plot <- reactive({
        df <- data_all()
        
        func <- switch(input$transformcov,
                       log=log,
                       sqrt=sqrt,
                       identity)
        
        if(input$datatype=='continuous'){
            new_dat <- data.frame(df[, c(input$y, input$D), drop=FALSE])
            if (input$transformcont == "log"){
                toExclude <- names(new_dat)[1]
                new_dat[,toExclude] <- log(new_dat[,toExclude])
                new_dat2 <- gather(data = new_dat, key, value, -toExclude)
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = toExclude)) + 
                    facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="", y=paste0("Log-", toExclude))
                pp
            }else{
                toExclude <- names(new_dat)[1]
                # new_dat[,toExclude] <- func(new_dat[,toExclude])
                new_dat2 <- gather(data = new_dat, key, value, -toExclude)
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = toExclude)) + 
                    facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="")
                pp
                
            }
        }else if (input$datatype=='prevalence'){
            
            if(input$transformprev == "logit"){
                new_dat <- data.frame(df[, c(input$p, input$m, input$D), drop=FALSE])
                new_dat[,"Emplogit"] <- log((new_dat[,input$p] + 0.5)/(new_dat[, input$m] - new_dat[,input$p] + 0.5))
                new_dat2 <- gather(data = new_dat[, -c(1,2)], key, value, -Emplogit)
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "Emplogit")) + 
                    facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + 
                    labs(x="", y=paste0("Emp-logit prevalence"))
                pp
            }else if (input$transformprev == "log"){
                new_dat <- data.frame(df[, c(input$p, input$m, input$D), drop=FALSE])
                new_dat[,"logprev"] <- log((new_dat[,input$p])/(new_dat[, input$m]))
                new_dat2 <- gather(data = new_dat[, -c(1,2)], key, value, -logprev)
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "logprev")) + 
                    facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="", y=paste0("Log-prevalence"))
                pp
            }else{
                new_dat <- data.frame(df[, c(input$p, input$m, input$D), drop=FALSE])
                new_dat[,"pprev"] <- as.numeric((new_dat[,input$p])/(new_dat[, input$m]))
                new_dat2 <- gather(data = new_dat[, -c(1,2)], key, value, -pprev)
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "pprev")) + 
                    facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="", y=paste0("Prevalence"))
                pp
                
            }
        } else{
            
            if (input$transformcnt == "log"){
                new_dat <- data.frame(df[, c(input$c, input$e, input$D), drop=FALSE])
                new_dat[,"logincidence"] <- log((new_dat[,input$c])/(new_dat[, input$e]))
                new_dat2 <- gather(data = new_dat[, -c(1,2)], key, value, -logincidence)
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "logincidence")) + 
                    facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="", y=paste0("Log-incidence"))
                pp
            }else{
                new_dat <- data.frame(df[, c(input$c, input$e, input$D), drop=FALSE])
                new_dat[,"iincidence"] <- (new_dat[,input$c])/(new_dat[, input$e])
                new_dat2 <- gather(data = new_dat[, -c(1,2)], key, value, -iincidence)
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "iincidence")) + 
                    facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="", y=paste0("Incidence"))
                pp
                
            }
        }
        
    })
    
    
    output$Plot <- renderPlot({
        if (is.null(scatter_ass_plot())) return(NULL)
        scatter_ass_plot()
    })
    
    var_plot_sum <- reactive({
        df <- data_all()
        envestatus <- switch(input$envelop, vario=1,  varifit = 2, varioEnve = 3)
        if(input$datatype=='continuous'){
            if(is.null(input$D)){
                # xmat <- as.matrix(cbind(rep(intercept=1, nrow(df))))
                # fml <- as.formula(paste("Prevalence ~ ", paste(colnames(xmat), collapse= "+"), paste0("+ (1|ID)")))
                fml <- as.formula(paste(paste0(input$y, " ~ 1")))
                # temp.fit <- glmer(formula = fml, data = dat, family = gaussian, control=lmerControl(check.nobs.vs.nlev="ignore"))
                temp.fit <- lm(formula = fml, data = df)
                # temp.fit <- lm(as.matrix(df[, input$y]) ~ xmat + 0)
                # temp.fit <- glmer(formula =  as.matrix(df[, input$y]) ~ xmat + 0)
                # beta.ols <- temp.fit$coeff
                residd <- resid(temp.fit)
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                #                 data = residd, max.dist = input$dist)
                
                
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
                if(input$maptype == 'plot'){
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
                }else{
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                        st_transform(., crs=utmcode)
                }
                if(input$functions=="matern"){
                    fix.kappa = FALSE
                } else {
                    fix.kappa = TRUE
                }
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions, fix.kappa=fix.kappa, bins = input$nbins, nsim = input$npermute)
                # if(envestatus == 2) vari <<- plo$summ
                plo$utmcode <- utmcode
                plo
            }
            else{
                fml <- as.formula(paste(paste0(input$y, " ~ ", paste(input$D, collapse= "+"))))
                # xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                temp.fit <- lm(formula = fml, data = df)
                # beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                #                 data = residd, max.dist = input$dist)
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
                if(input$maptype == 'plot'){
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
                }else{
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                        st_transform(., crs=utmcode)
                }
                if(input$functions=="matern") fix.kappa = FALSE else fix.kappa = TRUE
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions, fix.kappa=fix.kappa, bins = input$nbins, nsim = input$npermute)
                # if(envestatus== 2) vari <<- plo$summ
                plo$utmcode <- utmcode
                plo
            }
            
        } else if(input$datatype=='prevalence'){
            if(is.null(input$D)){
                # xmat <- as.matrix(cbind(rep(1, nrow(df))))
                # logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                # temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
                
                # fml <- as.formula(paste(paste0("cbind(", input$m, "-", input$p, ",", input$m, ")~ 1")))
                # temp.fit <- glm(formula = fml, data = df, family = binomial)
                
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # residd <- residuals(temp.fit)
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
                if(input$maptype == 'plot'){
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
                }else{
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                        st_transform(., crs=utmcode)
                }
                if(input$functions=="matern") fix.kappa = FALSE else fix.kappa = TRUE
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions, fix.kappa=fix.kappa, bins = input$nbins, nsim = input$npermute)
                
                ### get the utmcode used 
                plo$utmcode <- utmcode
                plo
            } else{
                xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
                # 
                
                # fml <- as.formula(paste(paste0("cbind(", input$m, "-", input$p, ",", input$m, ") ~ ", paste(input$D, collapse= "+"))))
                # temp.fit <- glm(formula = fml, data = df, family = binomial)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # residd <- residuals(temp.fit)
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
                if(input$maptype == 'plot'){
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
                }else{
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                        st_transform(., crs=utmcode)
                }
                if(input$functions=="matern") fix.kappa = FALSE else fix.kappa = TRUE
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions, fix.kappa=fix.kappa, bins = input$nbins, nsim = input$npermute)
                
                #### get the utm code used
                plo$utmcode <- utmcode
                plo
            }
        }else{
            if(is.null(input$D)){
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                logc <- log((df[, input$c]+1)/(df[, input$e]))
                temp.fit <- lm(as.matrix(logc) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
                if(input$maptype == 'plot'){
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
                }else{
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                        st_transform(., crs=utmcode)
                }
                if(input$functions=="matern") fix.kappa = FALSE else fix.kappa = TRUE
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions, fix.kappa=fix.kappa, bins = input$nbins, nsim = input$npermute)
                
                #### get the utm code used
                plo$utmcode <- utmcode
                plo
            } else{
                xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                logc <- log((df[, input$c]+1)/(df[, input$e]))
                temp.fit <- lm(as.matrix(logc) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
                if(input$maptype == 'plot'){
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
                }else{
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                        st_transform(., crs=utmcode)
                }
                if(input$functions=="matern") fix.kappa = FALSE else fix.kappa = TRUE
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions, fix.kappa=fix.kappa, bins = input$nbins, nsim = input$npermute)
                
                #### get the utm code used
                plo$utmcode <- utmcode
                plo
            }
        }
    })
    
    
    output$variogplot <- renderPlot({
        if (is.null(var_plot_sum())) return(NULL)
        var_plot_sum()$pl
        # myvariogramplot(vario)
    })
    
    output$summary <- renderPrint({
        if (is.null(var_plot_sum())) return(NULL)
        var_plot_sum()$summ
    })
    
    
    model.fit <- eventReactive(input$ShowEst, {
        df <- data_all()
        if(input$datatype=='continuous'){
            if(is.null(input$D)){
                fml <- as.formula(paste(paste0(input$y, " ~ 1")))
            } else{
                fml <- as.formula(paste(paste0(input$y, " ~ ", paste(input$D, collapse= "+"))))
            }
            coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
            # utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
            utmcode <- var_plot_sum()$utmcode
            if(input$maptype == 'plot'){
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
            }else{
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                    st_transform(., crs=utmcode)
            }
            coords <- st_coordinates(coords)
            df[, "XXX"] <- coords[, "X"]
            df[, "YYY"] <- coords[, "Y"]
            if(input$nu > 0){
                fit.MLE <- linear.model.MLE(formula = fml,coords=as.formula(paste("~", paste(c("XXX", "YYY"), collapse= "+"))),
                                            data=df, start.cov.pars=c(input$phi, input$nu),
                                            kappa=input$kappa, messages = F, method = "nlminb")
            }else{
                fit.MLE <- linear.model.MLE(formula = fml,coords=as.formula(paste("~", paste(c("XXX", "YYY"), collapse= "+"))),
                                            data=df, start.cov.pars=c(input$phi),
                                            kappa=input$kappa, messages = F, method = "nlminb", fixed.rel.nugget = TRUE)
            }
            
            fit.MLE$fml <- fml
            fit.MLE$utmcode <- utmcode
            fit.MLE
        } else if(input$datatype=='prevalence'){
            
            if(input$fitlinear == "binomialmodel"){
                if(is.null(input$D)){
                    fml <- as.formula(paste(paste0(input$p, " ~ 1")))
                    xmat <- as.matrix(cbind(rep(1, nrow(df))))
                } else{
                    fml <- as.formula(paste(paste0(input$p, " ~ ", paste(input$D, collapse= "+"))))
                    xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                }
                control.mcmc <- control.mcmc.MCML(n.sim=input$mcmcNsim,burnin=input$mcmcNburn,thin=input$mcmcNthin)
                ####
                logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
                # 
                
                # fml <- as.formula(paste(paste0("cbind(", input$m, "-", input$p, ",", input$m, ") ~ ", paste(input$D, collapse= "+"))))
                # temp.fit <- glm(formula = fml, data = df, family = binomial)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # par0 <- c(beta.ols, var(residd), input$phi, input$nu*var(residd))
                ##### conversion to utm
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                # utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
                utmcode <- var_plot_sum()$utmcode
                if(input$maptype == 'plot'){
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
                }else{
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                        st_transform(., crs=utmcode)
                }
                coords <- st_coordinates(coords)
                df[, "XXX"] <- coords[, "X"]
                df[, "YYY"] <- coords[, "Y"]
                ########
                if(input$nu > 0){
                    par0 <- c(beta.ols, var(residd), input$phi, input$nu*var(residd))
                    fit.MCML <- binomial.logistic.MCML(formula = fml,
                                                       coords=as.formula(paste("~", paste(c("XXX", "YYY"), collapse= "+"))),
                                                       data=df, start.cov.pars=c(input$phi, input$nu), units.m= as.formula(paste("~", input$m)),
                                                       kappa=input$kappa, messages = F, method = "nlminb", control.mcmc = control.mcmc, par0= par0)
                    
                } else{
                    par0 <- c(beta.ols, var(residd), input$phi)
                    fit.MCML <- binomial.logistic.MCML(formula = fml,
                                                       coords=as.formula(paste("~", paste(c("XXX", "YYY"), collapse= "+"))),
                                                       data=df, start.cov.pars=c(input$phi), units.m= as.formula(paste("~", input$m)),
                                                       kappa=input$kappa, messages = F, method = "nlminb", 
                                                       control.mcmc = control.mcmc, par0= par0, fixed.rel.nugget = TRUE)
                    
                }
                fit.MCML$fml <- fml
                fit.MCML$utmcode <- utmcode
                fit.MCML
            }else if(input$fitlinear == "linearmodel"){
                emplogit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                df[, "emplogit"] <- emplogit
                if(is.null(input$D)){
                    fml <- as.formula(paste(paste0("emplogit", " ~ 1")))
                } else{
                    fml <- as.formula(paste(paste0("emplogit", " ~ ", paste(input$D, collapse= "+"))))
                }
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                # utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
                utmcode <- var_plot_sum()$utmcode
                if(input$maptype == 'plot'){
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
                }else{
                    coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                        st_transform(., crs=utmcode)
                }
                coords <- st_coordinates(coords)
                df[, "XXX"] <- coords[, "X"]
                df[, "YYY"] <- coords[, "Y"]
                if(input$nu > 0){
                    fit.MCML <- linear.model.MLE(formula = fml,coords=as.formula(paste("~", paste(c("XXX", "YYY"), collapse= "+"))),
                                                data=df, start.cov.pars=c(input$phi, input$nu),
                                                kappa=input$kappa, messages = F, method = "nlminb")
                }else{
                    fit.MCML <- linear.model.MLE(formula = fml,coords=as.formula(paste("~", paste(c("XXX", "YYY"), collapse= "+"))),
                                                data=df, start.cov.pars=c(input$phi),
                                                kappa=input$kappa, messages = F, method = "nlminb", fixed.rel.nugget = TRUE)
                }
                
                fit.MCML$fml <- fml
                fit.MCML$utmcode <- utmcode
                fit.MCML
            }
            
            
        }else{
            if(is.null(input$D)){
                fml <- as.formula(paste(paste0(input$c, " ~ 1")))
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
            } else{
                fml <- as.formula(paste(paste0(input$c, " ~ ", paste(input$D, collapse= "+"))))
                xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
            }
            control.mcmc <- control.mcmc.MCML(n.sim=input$mcmcNsim,burnin=input$mcmcNburn,thin=input$mcmcNthin)
            ####
            logit <- log((df[, input$c]+1)/ (df[, input$e]))
            temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
            # 
            
            # fml <- as.formula(paste(paste0("cbind(", input$m, "-", input$p, ",", input$m, ") ~ ", paste(input$D, collapse= "+"))))
            # temp.fit <- glm(formula = fml, data = df, family = binomial)
            beta.ols <- temp.fit$coeff
            residd <- temp.fit$residuals
            # par0 <- c(beta.ols, var(residd), input$phi, input$nu*var(residd))
            ##### conversion to utm
            coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
            # utmcode <- epsgKM(as.numeric(lonlat2UTM(coords[1,])))
            utmcode <- var_plot_sum()$utmcode
            if(input$maptype == 'plot'){
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis))
            }else{
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                    st_transform(., crs=utmcode)
            }
            coords <- st_coordinates(coords)
            df[, "XXX"] <- coords[, "X"]
            df[, "YYY"] <- coords[, "Y"]
            ########
            if(input$nu > 0){
                par0 <- c(beta.ols, var(residd), input$phi, input$nu*var(residd))
                fit.MCML <- poisson.log.MCML(formula = fml,
                                             coords=as.formula(paste("~", paste(c("XXX", "YYY"), collapse= "+"))),
                                             data=df, start.cov.pars=c(input$phi, input$nu), units.m= as.formula(paste("~", input$e)),
                                             kappa=input$kappa, messages = F, method = "nlminb", control.mcmc = control.mcmc, par0= par0)
                
            } else{
                par0 <- c(beta.ols, var(residd), input$phi)
                fit.MCML <- poisson.log.MCML(formula = fml,
                                             coords=as.formula(paste("~", paste(c("XXX", "YYY"), collapse= "+"))),
                                             data=df, start.cov.pars=c(input$phi), units.m= as.formula(paste("~", input$e)),
                                             kappa=input$kappa, messages = F, method = "nlminb", 
                                             control.mcmc = control.mcmc, par0= par0, fixed.rel.nugget = TRUE)
                
            }
            fit.MCML$fml <- fml
            fit.MCML$utmcode <- utmcode
            fit.MCML
        }
    })
    
    output$estsummary <- renderPrint({
        if (is.null(model.fit())) return(NULL)
        summary(model.fit(), log.cov.pars = F)
    })
    
    output$tab <- renderTable({
        if (is.null(model.fit())) return(NULL)
        if(input$gotab > 0 ){
            create_tab <- function(x) {
                pars <- as.numeric(x$estimate)
                estimates <- PrevMap:::coef.PrevMap(x) 
                n <- length(estimates)
                ncov <- n - (which(grepl(pattern = paste0("sigma", collapse = "|"),  names(estimates)))-1)
                p <-  n - ncov
                estimates <- round(estimates, 4)
                # print(estimates)
                
                
                se <- sqrt(diag(x$covariance))
                ci_up <- pars + qnorm(0.975) * se
                ci_low <- pars - qnorm(0.975) * se
                
                #### delta method
                # se <- msm::deltamethod(g = sapply(c(sapply(1:p, function(x) paste0("~x",x)), 
                #                                     sapply((p+1):(n-1), function(x) paste0("~exp(x",x, ")")), 
                #                                     paste0("~exp(x",n, "/", "x", p+1, ")")), 
                #                                   function(y) as.formula(y)),
                #                        mean = x$estimate, cov = x$covariance)
                # ci_up <- estimates + qnorm(0.975) * se
                # ci_low <- estimates - qnorm(0.975) * se
                
                
                ci <- cbind(ci_low, ci_up) 
                ci[p + ncov, ] <- ci[p + ncov, ] + ci[p + 1, ]
                ci[(p + 1):length(estimates), ] <- exp(ci[(p + 1):length(estimates), ])
                ci <- round(ci, 4)
                ci <- apply(ci, 1, function(x) paste0("(", x[1], ", ", x[2], ")"))
                tab <- data.frame(Parameter = names(estimates), Estimate = as.character(estimates), CI = ci)
                return(tab)
            }
            
            fit_tab <- create_tab(model.fit())
            fit_tab
            # aaa <- PrevMap:::summary.PrevMap(model.fit(), log.cov.pars = F)
            # tibble::rownames_to_column(data.frame(aaa$coefficients), "covariates")
        }else{
            return(NULL)
        }
        
    })
    
    
    pred.fit <- eventReactive(input$ShowPred, {
        fit <- model.fit()
        coords <- fit$coords
        if(is.null(input$gridpreddata)){
            if(is.null(input$mbgshp)){
                poly <- coords[chull(coords),]
                gridpred <- gridpts(poly, xs = input$resolution, ys = input$resolution)
            }else{
                ##### make grid inside the polygon
                shp <- map_all()
                shp <- shp %>% st_transform(., crs=fit$utmcode)
                poly <- st_coordinates(shp)[,1:2]
                gridpred <- gridpts(poly, xs = input$resolution, ys = input$resolution)
            }
        }else{
            gridpred <- gridpred()
        }
        if (is.null(input$predictorsdata)){
            predictors <- NULL
        } else{
            predictors <- data.frame(predictors())
        }
        fml <<- fit$fml
        if(input$datatype=='continuous'){
            pred.mle <- spatial.pred.linear.MLE(
                object=fit,
                grid.pred=gridpred,
                predictors = predictors,
                predictors.samples = NULL,
                type = "marginal",
                scale.predictions = c("logit", "odds"),
                quantiles = c(0.025, 0.975),
                n.sim.prev = 1000,
                standard.errors = FALSE,
                thresholds = NULL,
                scale.thresholds = NULL,
                messages = TRUE,
                include.nugget = FALSE
            )
            res_df <- data.frame(pred.mle$grid, pred.mle$samples)
            res_df
        } else if(input$datatype=='prevalence'){
            if(input$fitlinear == "binomialmodel"){
                control.mcmc <- control.mcmc.MCML(n.sim=input$mcmcNsim,burnin=input$mcmcNburn,thin=input$mcmcNthin)
                # print(head(predictors))
                pred.mle <- spatial.pred.binomial.MCML(
                    object=fit,
                    grid.pred=gridpred,
                    predictors = predictors,
                    control.mcmc = control.mcmc,
                    type = "marginal",
                    scale.predictions = c("logit", "prevalence", "odds"),
                    quantiles = c(0.025, 0.975),
                    standard.errors = FALSE,
                    thresholds = NULL,
                    scale.thresholds = NULL,
                    plot.correlogram = FALSE,
                    messages = TRUE
                )
                res_df <- data.frame(pred.mle$grid, plogis(pred.mle$samples))
                res_df
            } else if(input$fitlinear == "linearmodel"){
                pred.mle <- spatial.pred.linear.MLE(
                    object=fit,
                    grid.pred=gridpred,
                    predictors = predictors,
                    predictors.samples = NULL,
                    type = "marginal",
                    scale.predictions = c("logit", "odds"),
                    quantiles = c(0.025, 0.975),
                    n.sim.prev = 1000,
                    standard.errors = FALSE,
                    thresholds = NULL,
                    scale.thresholds = NULL,
                    messages = TRUE,
                    include.nugget = FALSE
                )
                res_df <- data.frame(pred.mle$grid, 1/(1+exp(-pred.mle$samples)))
                res_df
            }

        }else{
            control.mcmc <- control.mcmc.MCML(n.sim=input$mcmcNsim,burnin=input$mcmcNburn,thin=input$mcmcNthin)
            # print(head(predictors))
            pred.mle <- spatial.pred.poisson.MCML(
                object=fit,
                grid.pred=gridpred,
                predictors = predictors,
                control.mcmc = control.mcmc,
                type = "marginal",
                scale.predictions = c("log", "exponential"),
                quantiles = c(0.025, 0.975),
                standard.errors = FALSE,
                thresholds = NULL,
                scale.thresholds = NULL,
                plot.correlogram = FALSE,
                messages = TRUE
            )
            res_df <- data.frame(pred.mle$grid, exp(pred.mle$samples))
            res_df
        }
    })
    
    ##############################################################################
    ### Producing the map with leaft where the projection is supplied ###
    ###################################################################################
    pred_map_lf <- reactive({
        if(input$datatype=='continuous'){
            if (is.null(pred.fit())) return(NULL)
            all_df <- pred.fit()
            if(input$predtomapcont == "meann"){
                ras_dff <- data.frame(all_df[, 1:2], prevalence = apply(all_df[, - c(1:2)], 1, mean))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="prevalence", style="quantile", title = "Mean outcome",
                                  alpha=0.5, palette="-RdYlBu", contrast=1) +
                        tm_layout())
                l
                
            }else if(input$predtomapcont == "sdd"){
                ras_dff <- data.frame(all_df[, 1:2], stderror = apply(all_df[, - c(1:2)], 1, sd))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="stderror", style="quantile", title = "Standard error",
                                  alpha=0.5, palette="-RdYlBu", contrast=1) +
                        tm_layout())
                l
            }else if(input$predtomapcont == "exprob"){
                ras_dff <- data.frame(all_df[, 1:2], exprob = apply(all_df[, - c(1:2)], 1, function(x) mean(x>input$threshold)))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                brks <- c(0,  0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1)
                labs <- create_labels(brks, greater = F)
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="exprob", style="fixed", alpha=0.5, palette=pal, contrast=1, labels = labs, breaks = brks, 
                                  title = paste0("Ex-prob ", input$threshold*100, "%")) +
                        tm_layout())
                l
            }else if(input$predtomapcont == "quant"){
                ras_dff <- data.frame(all_df[, 1:2], quantile = apply(all_df[, - c(1:2)], 1, 
                                                                      function(x) quantile(x = x, probs= input$quantprob)))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <<- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="quantile", style="quantile", alpha=0.5, palette="-RdYlBu", contrast=1,
                                  title = paste0(input$quantprob*100, "%", " Quantile")) +
                        tm_layout())
                l
            }
        }else if(input$datatype=='prevalence'){
            if (is.null(pred.fit())) return(NULL)
            all_df <- pred.fit()
            if(input$predtomapprev == "meann"){
                ras_dff <- data.frame(all_df[, 1:2], prevalence = apply(all_df[, - c(1:2)], 1, mean))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="prevalence", style="quantile", title = "Prevalence",
                                  alpha=0.5, palette="-RdYlBu", contrast=1) +
                        tm_layout())
                l
                
            }else if(input$predtomapprev == "sdd"){
                ras_dff <- data.frame(all_df[, 1:2], stderror = apply(all_df[, - c(1:2)], 1, sd))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="stderror", style="quantile", title = "Standard error",
                                  alpha=0.5, palette="-RdYlBu", contrast=1) +
                        tm_layout())
                l
            }else if(input$predtomapprev == "exprob"){
                ras_dff <- data.frame(all_df[, 1:2], exprob = apply(all_df[, - c(1:2)], 1, function(x) mean(x>input$threshold)))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                brks <- c(0,  0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1)
                labs <- create_labels(brks, greater = F)
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="exprob", style="fixed", alpha=0.5, palette=pal, contrast=1, labels = labs, breaks = brks, 
                                  title = paste0("Ex-prob ", input$threshold*100, "%")) +
                        tm_layout())
                l
            }else if(input$predtomapprev == "quant"){
                ras_dff <- data.frame(all_df[, 1:2], quantile = apply(all_df[, - c(1:2)], 1, 
                                                                      function(x) quantile(x = x, probs= input$quantprob)))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="quantile", style="quantile", title = paste0(input$quantprob*100, "%", " Quantile"),
                                  alpha=0.5, palette="-RdYlBu", contrast=1) +
                        tm_layout())
                l
            }
            
            
        }else{
            if (is.null(pred.fit())) return(NULL)
            all_df <- pred.fit()
            if(input$predtomapcount == "meann"){
                ras_dff <- data.frame(all_df[, 1:2], incidence = apply(all_df[, - c(1:2)], 1, mean))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="incidence", style="quantile", title = "Incidence",
                                  alpha=0.5, palette="-RdYlBu", contrast=1) +
                        tm_layout())
                l
                
            }else if(input$predtomapcount == "sdd"){
                ras_dff <- data.frame(all_df[, 1:2], stderror = apply(all_df[, - c(1:2)], 1, sd))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="stderror", style="quantile", title = "Standard error",
                                  alpha=0.5, palette="-RdYlBu", contrast=1) +
                        tm_layout())
                l
            }else if(input$predtomapcount == "exprob"){
                ras_dff <- data.frame(all_df[, 1:2], exprob = apply(all_df[, - c(1:2)], 1, function(x) mean(x>input$threshold)))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                brks <- c(0,  0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1)
                labs <- create_labels(brks, greater = F)
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="exprob", style="fixed", alpha=0.5, palette=pal, contrast=1, labels = labs, breaks = brks, 
                                  title = paste0("Ex-prob ", input$threshold*100, "%")) +
                        tm_layout())
                l
            }else if(input$predtomapcount == "quant"){
                ras_dff <- data.frame(all_df[, 1:2], quantile = apply(all_df[, - c(1:2)], 1, 
                                                                      function(x) quantile(x = x, probs= input$quantprob)))
                pred.raster <- raster::rasterFromXYZ(ras_dff, 
                                                     crs = var_plot_sum()$utmcode)
                l <- tmap::tmap_leaflet(
                    tmap::tm_shape(pred.raster) +
                        tm_raster(col="quantile", style="quantile", title = paste0(input$quantprob*100, "%", " Quantile"),
                                  alpha=0.5, palette="-RdYlBu", contrast=1) +
                        tm_layout())
                l
            }
            
        }
    })
    
    output$predmap <- renderLeaflet({
        if (is.null(pred_map_lf())) return(NULL)
        pred_map_lf()
    })
    
    ##############################################################################
    ### Removing the basemap where the projection is not supplied ###
    ###################################################################################
    pred_map_st <- reactive({
        if(input$datatype=='continuous'){
            if (is.null(pred.fit())) return(NULL)
            all_df <- pred.fit()
            if(input$predtomapcont == "meann"){
                ras_dff <- data.frame(all_df[, 1:2], outcome = apply(all_df[, - c(1:2)], 1, mean))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="outcome", style="quantile", title = "Mean outcome",
                              alpha=0.5, palette="-RdYlBu", contrast=1) +
                    tm_layout()
                l
                
            }else if(input$predtomapcont == "sdd"){
                ras_dff <- data.frame(all_df[, 1:2], stderror = apply(all_df[, - c(1:2)], 1, sd))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="stderror", style="quantile", title = "Standard error",
                              alpha=0.5, palette="-RdYlBu", contrast=1) +
                    tm_layout()
                l
            }else if(input$predtomapcont == "exprob"){
                ras_dff <- data.frame(all_df[, 1:2], exprob = apply(all_df[, - c(1:2)], 1, function(x) mean(x>input$threshold)))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                brks <- c(0,  0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1)
                labs <- create_labels(brks, greater = F)
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="exprob", style="fixed", alpha=0.5, palette=pal, contrast=1, labels = labs, breaks = brks, 
                              title = paste0("Ex-prob ", input$threshold*100, "%")) +
                    tm_layout()
                l
            }else if(input$predtomapcont == "quant"){
                ras_dff <- data.frame(all_df[, 1:2], quantile = apply(all_df[, - c(1:2)], 1, 
                                                                      function(x) quantile(x = x, probs= input$quantprob)))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="quantile", style="quantile", alpha=0.5, palette="-RdYlBu", contrast=1,
                              title = paste0(input$quantprob*100, "%", " Quantile")) +
                    tm_layout()
                l
            }
        }else if(input$datatype=='prevalence'){
            if (is.null(pred.fit())) return(NULL)
            all_df <- pred.fit()
            if(input$predtomapprev == "meann"){
                ras_dff <- data.frame(all_df[, 1:2], prevalence = apply(all_df[, - c(1:2)], 1, mean))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="prevalence", style="quantile", title = "Prevalence",
                              alpha=0.5, palette="-RdYlBu", contrast=1) +
                    tm_layout()
                l
                
            }else if(input$predtomapprev == "sdd"){
                ras_dff <- data.frame(all_df[, 1:2], stderror = apply(all_df[, - c(1:2)], 1, sd))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="stderror", style="quantile", title = "Standard error",
                              alpha=0.5, palette="-RdYlBu", contrast=1) +
                    tm_layout()
                l
            }else if(input$predtomapprev == "exprob"){
                ras_dff <- data.frame(all_df[, 1:2], exprob = apply(all_df[, - c(1:2)], 1, function(x) mean(x>input$threshold)))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                brks <- c(0,  0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1)
                labs <- create_labels(brks, greater = F)
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="exprob", style="fixed", alpha=0.5, palette=pal, contrast=1, labels = labs, breaks = brks, 
                              title = paste0("Ex-prob ", input$threshold*100, "%")) +
                    tm_layout()
                l
            }else if(input$predtomapprev == "quant"){
                ras_dff <- data.frame(all_df[, 1:2], quantile = apply(all_df[, - c(1:2)], 1, 
                                                                      function(x) quantile(x = x, probs= input$quantprob)))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="quantile", style="quantile", title = paste0(input$quantprob*100, "%", " Quantile"),
                              alpha=0.5, palette="-RdYlBu", contrast=1) +
                    tm_layout()
                l
            }
            
            
        }else{
            if (is.null(pred.fit())) return(NULL)
            all_df <- pred.fit()
            if(input$predtomapcount == "meann"){
                ras_dff <- data.frame(all_df[, 1:2], incidence = apply(all_df[, - c(1:2)], 1, mean))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="incidence", style="quantile", title = "Incidence",
                              alpha=0.5, palette="-RdYlBu", contrast=1) +
                    tm_layout()
                l
                
            }else if(input$predtomapcount == "sdd"){
                ras_dff <- data.frame(all_df[, 1:2], stderror = apply(all_df[, - c(1:2)], 1, sd))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="stderror", style="quantile", title = "Standard error",
                              alpha=0.5, palette="-RdYlBu", contrast=1) +
                    tm_layout()
                l
            }else if(input$predtomapcount == "exprob"){
                ras_dff <- data.frame(all_df[, 1:2], exprob = apply(all_df[, - c(1:2)], 1, function(x) mean(x>input$threshold)))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                brks <- c(0,  0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1)
                labs <- create_labels(brks, greater = F)
                pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="exprob", style="fixed", alpha=0.5, palette=pal, contrast=1, labels = labs, breaks = brks, 
                              title = paste0("Ex-prob ", input$threshold*100, "%")) +
                    tm_layout()
                l
            }else if(input$predtomapcount == "quant"){
                ras_dff <- data.frame(all_df[, 1:2], quantile = apply(all_df[, - c(1:2)], 1, 
                                                                      function(x) quantile(x = x, probs= input$quantprob)))
                pred.raster <- raster::rasterFromXYZ(ras_dff)
                l <- tmap::tm_shape(pred.raster) +
                    tm_raster(col="quantile", style="quantile", title = paste0(input$quantprob*100, "%", " Quantile"),
                              alpha=0.5, palette="-RdYlBu", contrast=1) +
                    tm_layout()
                l
            }
            
        }
    })
    
    ##################### Remove the basemap ##########################
    output$predmap2 <- renderPlot({
        if (is.null(pred_map_st())) return(NULL)
        pred_map_st()
    })
    

    
    ########################################################
    # The tab for downloading the result
    ########################################################
    
    params_func <- reactive({
        #### set which map to plot 
        if(input$maptype == 'view'){
            print(input$whattoshow)
            if(any(input$whattoshow == 'fig1')){
                exploremap = explore_map_lf()
            }else{
                exploremap = NULL
            }
            
        }else{
            print(input$whattoshow)
            if(any(input$whattoshow == 'fig1')){
                exploremap = explore_map_st()
            }else{
                exploremap = NULL
            }
        }
        
        ### set for scatter plot of the association
        if(any(input$whattoshow == 'fig2')){
            scatterplot = scatter_ass_plot()
        }else{
            scatterplot = NULL
        }
        
        # set for variogram
        if(any(input$whattoshow == 'fig3')){
            varplot  = var_plot_sum()$pl
        }else{
            varplot  = NULL
        }
        
        ### set the parameter estimate
        if(any(input$whattoshow == 'fig4')){
            parasumm = summary(model.fit(), log.cov.pars = F)
        }else{
            parasumm = NULL
        }
        
        
        
        ##set which one to map
        if(input$maptype == 'view'){
            if(any(input$whattoshow == 'fig5')){
                pred_map = pred_map_lf()
            }else{
                pred_map = NULL
            }
            
        }else{
            if(any(input$whattoshow == 'fig5')){
                pred_map = pred_map_st()
            }else{
                pred_map = NULL
            }
        }
        
        
        params <- list(nameofanalysis = input$datatype,  exploremap = exploremap,  scatterplot = scatterplot,
                       varplot  = varplot, parasumm = parasumm, predmap= pred_map)
        params
    })
    
    flname <- reactive({
        if(input$maptype == 'view'){
            filename <- "report.html"
            filename
        }else{
            filename <- "report.pdf"
            filename
        }
    })
    
    ################### for pdf format ###############################
    output$report <- downloadHandler(
        # For PDF output, change this to "report.pdf"
        filename = flname(),
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path(tempdir(), "report.Rmd")
            file.copy("report.Rmd", tempReport, overwrite = TRUE)
            # Set up parameters to pass to Rmd document
            params <- params_func()

            # print(params)
            tempReport <- file.path("report.Rmd")
            
            # Knit the document, passing in the `params` list, and eval it in a
            # child of the global environment (this isolates the code in the document
            # from the code in this app).
            rmarkdown::render(tempReport,  output_file = file,
                              params = params,
                              envir = new.env(parent = globalenv())
                              
            )
        }
    )
    
    ########### for html format ####################
    output$report2 <- downloadHandler(
        # For PDF output, change this to "report.pdf"
        filename = flname(),
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path(tempdir(), "report2.Rmd")
            file.copy("report2.Rmd", tempReport, overwrite = TRUE)
            # Set up parameters to pass to Rmd document
            params <- params_func()
            
            # print(params)
            tempReport <- file.path("report2.Rmd")
            
            # Knit the document, passing in the `params` list, and eval it in a
            # child of the global environment (this isolates the code in the document
            # from the code in this app).
            rmarkdown::render(tempReport,  output_file = file,
                              params = params,
                              envir = new.env(parent = globalenv())
                              
            )
        }
    )
    ########################################################
    # END download report
    ########################################################
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
