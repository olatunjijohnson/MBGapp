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
                    show_nbins = F, envelop=1, cov.model="matern") {
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
        
        vari.fit <- variofit(vario = empvario, ini.cov.pars=c(mean(dfvario$empirical), 1000), cov.model = cov.model,
                             fix.nug=F, nugget = 0, fix.kappa = F, kappa = 0.5)
        
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




############################## The begining of the APP ########################################################

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # img(src='chicas_logo.png', align = "right"),
    # # Application title
    titlePanel(title=div(img(src="chicas_logo.png", align = "right", height = 30, width = 100), "Model-based geostatistics")),
    
    
    # Sidebar with a slider input the data and the shapefile
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "mbgdata", label = "Upload a csv file:"),
            numericInput("crs", "coordinate reference system (optional):", 4326, min = 1, max = 100000),
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
                                 label = "select outcome",
                                 choices = ""
                             ), 
                             
                             
                             selectInput(
                                 inputId = "D",
                                 label = "covariates",
                                 choices = "",
                                 multiple = T
                             ),
                             radioButtons("transform", "Choose outcome transformation", c("No-transform" = "identity", 
                                                                                          "Log-transform" = "log", 
                                                                      "Logit-transform"= "logit")),
                             conditionalPanel(condition = "input.transform == 'logit'",
                                              selectInput(
                                                  inputId = "pp",
                                                  label = "Number of postives",
                                                  choices = ""
                                              ),
                                              selectInput(
                                                  inputId = "mm",
                                                  label = "Number Examined",
                                                  choices = ""
                                              )
                             ),
                             radioButtons("transformcov", "Choose covariate transformation ", c("No-transform" = "identity",
                                                                                               "Log-transform" = "log", 
                                                                                          "square-root transform"= "sqrt"))
                             
                             
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
                                 h2("Scatter plot of the outcome and the covariate"),
                                 plotOutput(outputId ="Plot")),
                        tabPanel("Variogram", value = 2,       
                                 plotOutput(outputId ="variogplot"),
                                 h2("Summary of estimate covariance parameter"),
                                 verbatimTextOutput(outputId ="summary")), id="tabselected"
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
        if(grepl("\\.rds$", dff)){
            x <- readRDS(dff$datapath)
            x
        }else{
            x <- read_csv(dff$datapath)
            x
        }
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
        updateVarSelectInput(session, "y", label = "select outcome", data = df2)
        updateVarSelectInput(session, "D", label = "covariate", data = df2)
        updateVarSelectInput(session, "yy", label = "select outcome", data = df2)
        updateVarSelectInput(session, "p", data = df2)
        updateVarSelectInput(session, "m", data = df2)
        updateVarSelectInput(session, "c", data = df2)
        updateVarSelectInput(session, "e", data = df2)
        updateVarSelectInput(session, "DD", label = "covariates", data = df2)
        updateVarSelectInput(session, "pp", data = df2)
        updateVarSelectInput(session, "mm", data = df2)
        
        
    })
    
    observeEvent(input$change,{
        df2 <- data_all()
        dummy_coords <- data.frame(df2[, c(input$xaxis,input$yaxis)])
        dummy_coords <- dummy_coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>%
            st_transform(., crs=as.numeric(lonlat2UTM(dummy_coords[1,])))
        variog_extent <- max(dist(cbind(st_coordinates(dummy_coords))), na.rm=T)
        updateSliderInput(session, "dist", min = 0, max = variog_extent/3,
                          value = variog_extent/10,
                          step = round(variog_extent/200)+1)
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
        
        func <- switch(input$transformcov,
               log=log,
               sqrt=sqrt,
               identity)
        
        
        if(input$transform == "logit"){
            new_dat <- data.frame(df[, c(input$pp, input$mm, input$D), drop=FALSE])
            
            # emplogit<-function(new_dat){
            #     Y <- new_dat[, input$pp]
            #     N <- new_dat[, input$mm]
            #     top=Y+0.5
            #     bottom=N-Y+0.5
            #     return(log(top/bottom))
            # }
            
            new_dat["Emplogit"] <- log((new_dat[,input$pp] + 0.5)/(new_dat[, input$mm] - new_dat[,input$pp] + 0.5))
            new_dat2 <- gather(data = new_dat[, -c(1,2)], key, value, -Emplogit)
            new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
            ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "Emplogit")) + 
                facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="", y=paste0("Emp-logit prevalence"))
        }else if (input$transform == "log"){
            new_dat <- data.frame(df[, c(input$y, input$D), drop=FALSE])
            toExclude <- names(new_dat)[1]
            new_dat[,toExclude] <- log(new_dat[,toExclude])
            new_dat2 <- gather(data = new_dat, key, value, -toExclude)
            new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
            ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = toExclude)) + 
                facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="", y=paste0("Log-", toExclude))
        }else{
            new_dat <- data.frame(df[, c(input$y, input$D), drop=FALSE])
            toExclude <- names(new_dat)[1]
            # new_dat[,toExclude] <- func(new_dat[,toExclude])
            new_dat2 <- gather(data = new_dat, key, value, -toExclude)
            new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
            ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = toExclude)) + 
                facet_wrap(facets = ~key, scales = "free_x") + geom_point() + geom_smooth() + labs(x="")
            
        }
        
    })
    
    var_plot_sum <- reactive({
        df <- data_all()
        envestatus <- switch(input$envelop, vario=1,  varifit = 2, varioEnve = 3)
        if(input$datatype=='continuous'){
            if(is.null(input$DD)){
                # xmat <- as.matrix(cbind(rep(intercept=1, nrow(df))))
                # fml <- as.formula(paste("Prevalence ~ ", paste(colnames(xmat), collapse= "+"), paste0("+ (1|ID)")))
                fml <- as.formula(paste(paste0(input$yy, " ~ 1")))
                # temp.fit <- glmer(formula = fml, data = dat, family = gaussian, control=lmerControl(check.nobs.vs.nlev="ignore"))
                temp.fit <- lm(formula = fml, data = df)
                # temp.fit <- lm(as.matrix(df[, input$yy]) ~ xmat + 0)
                # temp.fit <- glmer(formula =  as.matrix(df[, input$yy]) ~ xmat + 0)
                # beta.ols <- temp.fit$coeff
                residd <- resid(temp.fit)
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                #                 data = residd, max.dist = input$dist)
                
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                    st_transform(., crs=as.numeric(lonlat2UTM(coords[1,])))
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions)
                # if(envestatus == 2) vari <<- plo$summ
                plo
            } else{
                fml <- as.formula(paste(paste0(input$yy, " ~ ", paste(input$DD, collapse= "+"))))
                # xmat <- as.matrix(cbind(1, df[, input$DD, drop=FALSE]))
                temp.fit <- lm(formula = fml, data = df)
                # beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                # vario <- variog(coords = cbind(df[, c(input$xaxis,input$yaxis)]), 
                #                 data = residd, max.dist = input$dist)
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                    st_transform(., crs=as.numeric(lonlat2UTM(coords[1,])))
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions)
                # if(envestatus== 2) vari <<- plo$summ
                plo
            }
            
        } else if(input$datatype=='prevalence'){
            if(is.null(input$DD)){
                # xmat <- as.matrix(cbind(rep(1, nrow(df))))
                # logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                # temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
                
                fml <- as.formula(paste(paste0("cbind(", input$m, "-", input$p, ",", input$m, ")~ 1")))
                temp.fit <- glm(formula = fml, data = df, family = binomial)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                    st_transform(., crs=as.numeric(lonlat2UTM(coords[1,])))
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions)
                plo
            } else{
                # xmat <- as.matrix(cbind(1, df[, input$DD, drop=FALSE]))
                # logit <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                # temp.fit <- lm(as.matrix(logit) ~ xmat + 0)
                # 
                
                fml <- as.formula(paste(paste0("cbind(", input$m, "-", input$p, ",", input$m, ") ~ ", paste(input$DD, collapse= "+"))))
                temp.fit <- glm(formula = fml, data = df, family = binomial)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                    st_transform(., crs=as.numeric(lonlat2UTM(coords[1,])))
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions)
                plo
            }
        }else{
            if(is.null(input$DD)){
                xmat <- as.matrix(cbind(rep(1, nrow(df))))
                logc <- log((df[, input$c])/(df[, input$e]))
                temp.fit <- lm(as.matrix(logc) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                    st_transform(., crs=as.numeric(lonlat2UTM(coords[1,])))
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions)
                plo
            } else{
                xmat <- as.matrix(cbind(1, df[, input$DD, drop=FALSE]))
                logc <- log((df[, input$p] + 0.5)/ (df[, input$m] - df[, input$p] + 0.5))
                temp.fit <- lm(as.matrix(logc) ~ xmat + 0)
                beta.ols <- temp.fit$coeff
                residd <- temp.fit$residuals
                coords <- data.frame(df[, c(input$xaxis,input$yaxis)])
                coords <- coords %>% st_as_sf(., coords=c(input$xaxis, input$yaxis), crs= input$crs) %>% 
                    st_transform(., crs=as.numeric(lonlat2UTM(coords[1,])))
                plo <- ggvario(coords = st_coordinates(coords), data=residd, maxdist = input$dist, envelop = envestatus,
                               cov.model=input$functions)
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
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
