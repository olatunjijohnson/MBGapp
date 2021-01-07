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
library(tidyverse)
library(tmap)
library(sf)
library(leaflet)
library(rgdal)


########### useful functions to deal with variogram ###############
variog_envelope <- function (geodata, coords = geodata$coords, data = geodata$data,
                             obj.variog, nsim = 99, save.sim = FALSE, messages)
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

# Calculate and plot the variogram
ggvario <- function(coords,
                    data,
                    bins = 15,
                    maxdist = max(dist(coords))/3,
                    uvec = NULL,
                    nsim = 999,
                    color = "royalblue1",
                    xlab = "distance",
                    show_nbins = F, envelop=F) {
    require(geoR)
    coords <- as.matrix(coords)
    min_dist <- min(dist(coords))
    if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
    empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
    if(envelop==FALSE){
        dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                              nbins = empvario$n)
        p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
            geom_point(col = "black", fill = color, shape = 21, size = 3) +
            scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                               breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
            scale_y_continuous(name = "semivariance",
                               #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1),
                               limits = c(0,  max(dfvario$empirical))) +
            ggtitle("Empirical semivariogram")
        # theme_classic()
        p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
        if(show_nbins) p2 else p1
    }else{
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
        if(show_nbins) p2 else p1
    }
}

############################## The begining of the APP ########################################################

# Define UI for application that draws a histogram
ui <- fluidPage(


    # Application title
    titlePanel("Model-based geostatistics"),

    # Sidebar with a slider input the data and the shapefile
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "mbgdata", label = "Upload a csv file:"),
            numericInput("crs", "coordinate reference system:", 4326, min = 1, max = 100000),
            fileInput(inputId = "mbgshp", label = "Upload a shape file (optional):"),
            selectInput(
                inputId = "xaxis",
                label = "X axis",
                choices = "",
            ),
            selectInput(
                inputId = "yaxis",
                label = "Y axis",
                choices = "",
            ),

            conditionalPanel(condition = "input.tabselected==1",
                             #### select the longitude and latitude


                             selectInput(
                                 inputId = "y",
                                 label = "variable to map",
                                 choices = ""
                             ),


                             selectInput(
                                 inputId = "D",
                                 label = "covariate",
                                 choices = "",
                                 multiple = T
                             )

            ),
            conditionalPanel(condition = "input.tabselected==2",

                             selectInput("datatype", 'choose the type of data',
                                         choices=c("Continuous data" ='continuous', "Prevalence data" = 'prevalence', "Count data" = 'count'),
                                         selected = NULL),

                             conditionalPanel(condition = "input.datatype=='continuous'",
                                              selectInput(
                                                  inputId = "yy",
                                                  label = "continuous variable",
                                                  choices = ""
                                              )
                             ),
                             conditionalPanel(condition = "input.datatype=='prevalence'",
                                              selectInput(
                                                  inputId = "p",
                                                  label = "Number of postives",
                                                  choices = ""
                                              ),
                                              selectInput(
                                                  inputId = "m",
                                                  label = "Number Examined",
                                                  choices = ""
                                              )
                             ),
                             conditionalPanel(condition = "input.datatype=='count'",
                                              selectInput(
                                                  inputId = "c",
                                                  label = "count variable",
                                                  choices = ""
                                              ),
                                              selectInput(
                                                  inputId = "e",
                                                  label = "offset",
                                                  choices = ""
                                              )

                             ),
                             selectInput(
                                 inputId = "DD",
                                 label = "covariates",
                                 choices = "",
                                 multiple = T
                             ),
                             sliderInput(inputId = "dist",
                                         label = "Distance:",
                                         min = 0,
                                         max = 5,
                                         value = 3, step=1),

                             actionButton("change", "Change slider max value"),
                             br(),
                             br(),
                             radioButtons("envelop", "Choose plot", c("Variogram only" = "vario",
                                                                      "Variogram with envelope"= "varioEnve")),
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             )


            ),
        ),

        # Show a map and plot  of the data
        mainPanel(
            tabsetPanel(type="tabs",
                        tabPanel("Explore", value = 1,
                                 leafletOutput(outputId = "map"),
                                 plotOutput(outputId ="Plot")),
                        tabPanel("Variogram", value = 2,
                                 plotOutput(outputId ="variogplot")), id="tabselected"
            )

        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    data_all <- reactive({
        req(input$mbgdata)
        dff <- input$mbgdata
        if (is.null(dff))
            return(NULL)
        x <- read_csv(dff$datapath)
        x
    })

    map_all <- reactive({
        req(input$mbgshp)
        dff <- input$mbgshp
        if (is.null(dff))
            return(NULL)
        if(grepl("\\.rds$", dff)){
            x <- readRDS(dff$datapath)
            x
        }else{
            x <- rgdal::readOGR(dff$datapath)
            x
        }

    })
    observe({
        df2 <- data_all()
        updateVarSelectInput(session, "xaxis", label = "select the x axis", data = df2)
        updateVarSelectInput(session, "yaxis", label = "select the y axis", data = df2)
        updateVarSelectInput(session, "y", label = "select the variable", data = df2)
        updateVarSelectInput(session, "D", label = "covariate", data = df2)
        updateVarSelectInput(session, "yy", label = "select the variable", data = df2)
        updateVarSelectInput(session, "p", data = df2)
        updateVarSelectInput(session, "m", data = df2)
        updateVarSelectInput(session, "c", data = df2)
        updateVarSelectInput(session, "e", data = df2)
        updateVarSelectInput(session, "DD", label = "covariates", data = df2)
    })

    observeEvent(input$change,{
        df2 <- data_all()
        updateSliderInput(session, "dist", min = 0, max = max(dist(cbind(df2[,c(input$xaxis, input$yaxis), drop=FALSE])), na.rm=T),
                          value = 4,
                          step = 1)
    })


    # observeEvent(input$change,{
    #   updateSliderInput(session, "dist", max = 50000, step = round(50000/50))
    # })

    output$map <- renderLeaflet({
        df <- data_all()
        if(is.null(input$mbgshp)){
            mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=input$crs)
            mapdata <- st_transform(mapdata, crs=4326)
            l <- tmap::tmap_leaflet(
                tmap::tm_shape(mapdata) +
                    tm_symbols(col=input$y, style="equal") +
                    # tm_shape(shp, is.master = T) +
                    # tm_borders(col="black") +
                    tm_layout()
            )
        }else{
            shp <- map_all()
            mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=crs(shp))
            shp <- st_transform(shp, crs=4326)
            mapdata <- st_transform(mapdata, crs=4326)


            l <- tmap::tmap_leaflet(
                tmap::tm_shape(mapdata) +
                    tm_symbols(col=input$y, style="equal") +
                    tm_shape(shp, is.master = T) +
                    tm_borders(col="black") +
                    tm_layout()
            )
        }





    })

    output$Plot <- renderPlot({
        df <- data_all()
        new_dat <- data.frame(df[, c(input$y, input$D), drop=FALSE])
        toExclude <- names(new_dat)[1]
        new_dat2 <- gather(data = new_dat, key, value, -toExclude)
        ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = toExclude)) +
            facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="")
    })


    output$variogplot <- renderPlot({
        df <- data_all()
        envestatus <- switch(input$envelop, vario=FALSE, varioEnve = TRUE)
        if(input$datatype=='continuous'){
            if(is.null(input$DD)){
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                temp.fit <- lm(as.matrix(df[, input$yy]) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]),
                #                 data = residd, max.dist = input$dist)
                ggvario(coords = cbind(df[, c(input$xaxis,input$yaxis)]), data=residd, maxdist = input$dist, envelop = envestatus)

            } else{
                xmat <- as.matrix(cbind(1, df[, input$DD, drop=FALSE]))
                temp.fit <- lm(as.matrix(df[, input$yy]) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]),
                #                 data = residd, max.dist = input$dist)
                ggvario(coords = cbind(df[, c(input$xaxis,input$yaxis)]), data=residd, maxdist = input$dist, envelop = envestatus)
            }

        } else if(input$datatype=='prevalence'){
            if(is.null(input$DD)){
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]),
                #                 data = residd, max.dist = input$dist)
                ggvario(coords = cbind(df[, c(input$xaxis,input$yaxis)]), data=residd, maxdist = input$dist, envelop = envestatus)
            } else{
                xmat <- as.matrix(cbind(1, df[, input$DD, drop=FALSE]))
                logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]),
                #                 data = residd, max.dist = input$dist)
                ggvario(coords = cbind(df[, c(input$xaxis,input$yaxis)]), data=residd, maxdist = input$dist, envelop = envestatus)
            }
        }else{
            if(is.null(input$DD)){
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                logc <- log((df[, input$c])/(df[, input$e]))
                temp.fit <- lm(as.matrix(logc) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]),
                #                 data = residd, max.dist = input$dist)
                ggvario(coords = cbind(df[, c(input$xaxis,input$yaxis)]), data=residd, maxdist = input$dist, envelop = envestatus)
            } else{
                xmat <- as.matrix(cbind(1, df[, input$DD, drop=FALSE]))
                logc <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(as.matrix(logc) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]),
                #                 data = residd, max.dist = input$dist)
                ggvario(coords = cbind(df[, c(input$xaxis,input$yaxis)]), data=residd, maxdist = input$dist, envelop = envestatus)
            }
        }
        # myvariogramplot(vario)
    })



}

# Run the application
shinyApp(ui = ui, server = server)
