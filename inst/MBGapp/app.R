#
# This is a Shiny web application for Model-based geostatistics. You can run the application by clicking
# the 'Run App' button above.
#
#
library(shiny)
library(geoR)
library(ggplot2)
library(magrittr)
library(dplyr)
library(readr)
library(tidyr)
library(sf)
library(leaflet)
# leafem: attached lazily via leafem:: (only addStarsImage is used) to cut startup time
# tidyterra: used only via tidyterra::geom_spatraster, so not attached at startup
library(shinyjs)
library(RiskMap)
library(terra)
require(grDevices)
library(splines)
# httr2: used only via httr2:: (Groq helper), so not attached at startup

options(shiny.maxRequestSize = 30*1024^2)
# jsCode <- "shinyjs.hideSidebar = function(params){$('body').addClass('sidebar-collapse');}"

########### Groq LLM helper ###############
call_groq <- function(prompt, api_key, model = "llama-3.3-70b-versatile", max_tokens = 700) {
    if (is.null(api_key) || nchar(trimws(api_key)) == 0) return(NULL)
    tryCatch({
        resp <- httr2::request("https://api.groq.com/openai/v1/chat/completions") |>
            httr2::req_headers(
                Authorization  = paste("Bearer", trimws(api_key)),
                `Content-Type` = "application/json"
            ) |>
            httr2::req_body_json(list(
                model       = model,
                messages    = list(
                    list(role = "system",
                         content = paste("You are a statistical analyst writing a scientific report.",
                                         "Explain geostatistical results clearly for a public health or",
                                         "environmental science audience. Avoid jargon. Use 2-3 concise paragraphs.")),
                    list(role = "user", content = prompt)
                ),
                max_tokens  = as.integer(max_tokens),
                temperature = 0.3
            )) |>
            httr2::req_timeout(30) |>
            httr2::req_perform()
        httr2::resp_body_json(resp)$choices[[1]]$message$content
    }, error = function(e) {
        paste0("[AI generation error: ", conditionMessage(e), "]")
    })
}

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


#### create formula
create_formula <- function(y, covars) {
    formula(paste(y, paste(covars, collapse = "+"), sep = "~"))
}

return_formula <- function(y, covars, nl_terms){
    covars <- covars
    fml <- create_formula(y = y, covars = covars)
    nl_terms <- nl_terms
    if(!is.null(nl_terms)) {
        id_rm <- sapply(covars,
                        function(x) any(grepl(paste0("\\b", x, "\\b"), nl_terms)))
        # covars[id_rm] <- strsplit(nl_terms, "\\+")[[1]]
        covars[id_rm] <- nl_terms
        fml <- create_formula(y = y, covars = covars)
    }
    return(fml)
}

return_formula2 <- function(y, covars, nl_terms){
    covars <- covars
    fml <- create_formula(y = y, covars = covars)
    nl_terms <- nl_terms
    if(!is.null(nl_terms)) {
        id_rm <- sapply(covars,
                        function(x) any(grepl(paste0("\\b", x, "\\b"), nl_terms)))
        covars[id_rm] <- strsplit(nl_terms, "\\+")[[1]]
        # covars[id_rm] <- nl_terms
        fml <- create_formula(y = y, covars = covars)
    }
    return(fml)
}

drop_formula_term <- function(the_formula, var_name) {
    f <- list()
    for(i in 1:length(var_name)){
        var_position <- grep(paste0("\\b", var_name[i], "\\b"), attr(terms(the_formula), "term.labels"))
        f0 <- update(the_formula,drop.terms(terms(the_formula), -var_position,
                                            keep.response=TRUE))
        ff <- gsub(pattern = paste0("\\b", all.vars(f0[[3]]), "\\b"), replacement = "x",  x=attr(terms(f0),"term.labels")[1])
        f[[var_name[i]]] <- create_formula("y", ff)

    }
    return(f)
}

############################## The begining of the APP ########################################################

# Define UI for application that draws a histogram
ui <- fluidPage(
    useShinyjs(),
    # extendShinyjs(text = jsCode, functions = c("hideSidebar")),
    # img(src='chicas_logo.png', align = "right"),
    # # Application title
    titlePanel(title=div(
        img(src="chicas_logo.png",          align="right", height=40, width=130, style="margin-left:8px;"),
        img(src="manchester_logo_big.gif",  align="right", height=40, style="margin-left:8px;"),
        "Model-based geostatistics"
    )),


    # Sidebar with a slider input the data and the shapefile
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "mbgdata", label = "Upload the data (csv file):"),
            radioButtons("maptype", "Do you know the projection of the location?:",
                         c("Yes" = "view",
                           "No" = "plot"), selected = "view"),
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
                label = "X-coordinate/Longitude",
                choices = "",
            ),
            selectInput(
                inputId = "yaxis",
                label = "Y-coordinate/Latitude",
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
                                 label = "Positives",
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
                                                                                                "square-root transform"= "sqrt",
                                                                                                "Personalised transformation" = "I"))



            ),
            conditionalPanel(condition = "input.transformcov=='I' | input.AdvOption",
                             textInput(inputId = "nl_terms",
                                       label = "Covariate Non-linear function:", value = ""),
                             conditionalPanel(condition = "input.transformcov=='I'",
                                              actionButton("showNL", "show fitted line"))),



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
                                                 onclick ="window.open('https://olatunjijohnson-variogramapp.hf.space', '_blank')"),

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

                             numericInput("phi", "Initial value of scale parameter", 50),
                             selectInput("includenugget", "Include the nugget effect", choices = c("Yes" = 1, "No" = 0)),
                             conditionalPanel(condition = "input.includenugget==1",
                                              numericInput("nu", "Initial value of relative variance of the nugget effect", 0.1)),
                             numericInput("kappa", "Value of kappa", 0.5),
                             actionButton("AdvOption", "Advanced options"),
                             # conditionalPanel(condition = "input.datatype !='continuous' & input.fitlinear=='binomialmodel'",
                             #                  actionButton("AdvOption", "Advanced options")),
                             conditionalPanel(condition = "(input.AdvOption & input.datatype =='prevalence' & input.fitlinear=='binomialmodel') | (input.AdvOption & input.datatype=='count')",
                                              numericInput("mcmcNsim", "Number of simulation", 5000),
                                              numericInput("mcmcNburn", "Number of burn-in", 1000),
                                              numericInput("mcmcNthin", "Number of thinning", 4)
                                              ),
                             uiOutput("model_selector_ui"),
                             actionButton("ShowEst", "Show the result summary", icon = icon("fas fa-running")),

                             #### This part helps to hide the error
                             tags$style(type="text/css",
                                        ".shiny-output-error { visibility: hidden; }",
                                        ".shiny-output-error:before { visibility: hidden; }"
                             )


            ),
            conditionalPanel(condition = "input.tabselected==4",
                             numericInput("resolution", label = "Spatial resolution (km)", value = 10, min = 0, max = 100000),
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
            conditionalPanel(condition = "input.tabselected==5",
                             p("Select the sections to include in your report, then click", strong("Download report"), "to generate a PDF.")
            ),
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
                                 h4("Model summary"),
                                 verbatimTextOutput(outputId ="estsummary"),
                                 hr(),
                                 h4("Parameter table"),
                                 tableOutput("tab")),
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
                                 helpText("Select sections to include, optionally add AI explanations, then download."),
                                 checkboxGroupInput("whattoshow", "What to show in the report", inline = FALSE,
                                                    c("Map of the outcome"               = "fig1",
                                                      "Scatter plot of outcome vs covariate" = "fig2",
                                                      "Variogram plot"                    = "fig3",
                                                      "Summary of the parameters"         = "fig4",
                                                      "Prediction map"                    = "fig5"),
                                                    selected = c("fig1")),

                                 hr(),
                                 h4("AI-generated explanations (Groq)"),
                                 uiOutput("llm_settings_ui"),
                                 uiOutput("llm_preview_ui"),
                                 HTML("<br>")),

                        id="tabselected"
            ),

            tags$hr(),
            tags$p(
                tags$small(
                    tags$strong("Authors: "),
                    "Olatunji Johnson (University of Manchester), ",
                    "Claudio Fronterre (University of Birmingham) and ",
                    "Emanuele Giorgi (University of Birmingham)."
                )
            )

        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    has_inla <- requireNamespace("INLA", quietly=TRUE)

    # ---- LLM reactive state ----
    llm_rv <- reactiveValues(
        data_text  = "",
        variog_text = "",
        est_text   = "",
        pred_text  = "",
        status     = ""
    )

    output$llm_status <- renderText({ llm_rv$status })

    output$llm_preview_ui <- renderUI({
        sections <- list()
        if (nchar(llm_rv$data_text) > 0) {
            sections <- c(sections, list(
                tags$h5("AI explanation — map of outcome:"),
                textAreaInput("edit_llm_data", NULL,
                              value = llm_rv$data_text, rows = 5, width = "100%")
            ))
        }
        if (nchar(llm_rv$variog_text) > 0) {
            sections <- c(sections, list(
                tags$h5("AI explanation — variogram:"),
                textAreaInput("edit_llm_variog", NULL,
                              value = llm_rv$variog_text, rows = 5, width = "100%")
            ))
        }
        if (nchar(llm_rv$est_text) > 0) {
            sections <- c(sections, list(
                tags$h5("AI explanation — model estimates:"),
                textAreaInput("edit_llm_est", NULL,
                              value = llm_rv$est_text, rows = 5, width = "100%")
            ))
        }
        if (nchar(llm_rv$pred_text) > 0) {
            sections <- c(sections, list(
                tags$h5("AI explanation — prediction map:"),
                textAreaInput("edit_llm_pred", NULL,
                              value = llm_rv$pred_text, rows = 5, width = "100%")
            ))
        }
        if (length(sections) > 0) {
            tagList(
                hr(),
                tags$p(tags$strong("Preview and edit AI explanations"),
                       " — the text below will be included in your downloaded report. Edit freely."),
                tagList(sections)
            )
        }
    })

    output$llm_settings_ui <- renderUI({
        key_set <- nchar(Sys.getenv("GROQ_API_KEY")) > 0
        model_choices <- c(
            "Llama 3.3 70B — best quality" = "llama-3.3-70b-versatile",
            "Llama 3.1 8B — fastest"       = "llama-3.1-8b-instant",
            "Mixtral 8x7B — balanced"      = "mixtral-8x7b-32768",
            "Gemma 2 9B"                   = "gemma2-9b-it"
        )
        if (key_set) {
            wellPanel(
                tags$p(icon("circle-check", style = "color:#28a745"),
                       tags$strong(" AI ready."),
                       " API key is pre-configured — no setup needed."),
                selectInput("groq_model", "AI model:", choices = model_choices,
                            selected = "llama-3.3-70b-versatile"),
                actionButton("gen_llm", "Generate AI explanations",
                             icon = icon("robot"), class = "btn-info"),
                tags$br(), tags$br(),
                tags$em(textOutput("llm_status", inline = TRUE))
            )
        } else {
            wellPanel(
                passwordInput("groq_api_key",
                              label = tags$span("Groq API key",
                                                tags$small(tags$a(" (get free key)",
                                                                  href   = "https://console.groq.com",
                                                                  target = "_blank"))),
                              placeholder = "gsk_..."),
                selectInput("groq_model", "AI model:", choices = model_choices,
                            selected = "llama-3.3-70b-versatile"),
                helpText("Your key is used only for this session and never stored."),
                actionButton("gen_llm", "Generate AI explanations",
                             icon = icon("robot"), class = "btn-info"),
                tags$br(), tags$br(),
                tags$em(textOutput("llm_status", inline = TRUE))
            )
        }
    })

    observeEvent(input$gen_llm, {
        api_key <- trimws(
            if (!is.null(input$groq_api_key) && nchar(trimws(input$groq_api_key)) > 0)
                input$groq_api_key
            else
                Sys.getenv("GROQ_API_KEY")
        )
        if (nchar(api_key) == 0) {
            llm_rv$status <- "Please enter a Groq API key (or set GROQ_API_KEY env var)."
            return()
        }
        model_id <- input$groq_model
        what     <- input$whattoshow
        df       <- tryCatch(data_all(), error = function(e) NULL)

        total_tasks <- sum(c(
            "fig1" %in% what && !is.null(df),
            "fig3" %in% what && !is.null(tryCatch(var_plot_sum(), error=function(e) NULL)),
            "fig4" %in% what && !is.null(tryCatch(model.fit(),   error=function(e) NULL)),
            "fig5" %in% what && !is.null(tryCatch(pred.fit(),    error=function(e) NULL))
        ))
        if (total_tasks == 0) {
            llm_rv$status <- "No sections with data are selected. Load data and run the analysis first."
            return()
        }

        done <- 0
        llm_rv$status <- paste0("Generating 0 / ", total_tasks, " explanations...")

        # --- Data / map ---
        if ("fig1" %in% what && !is.null(df)) {
            datatype <- input$datatype
            n_obs    <- nrow(df)
            prompt <- switch(
                datatype,
                continuous = {
                    yc  <- tryCatch(input$y, error = function(e) "outcome")
                    rng <- round(range(df[[yc]], na.rm = TRUE), 3)
                    paste0("A geostatistical survey recorded ", n_obs, " observations of the ",
                           "continuous variable '", yc, "' (range: ", rng[1], " to ", rng[2], "). ",
                           "Write 2-3 paragraphs for a scientific report explaining what the spatial ",
                           "map of these observations shows and why geostatistical methods are appropriate.")
                },
                prevalence = {
                    pc <- tryCatch(input$p, error = function(e) "positives")
                    mc <- tryCatch(input$m, error = function(e) "examined")
                    emp <- if (!is.null(df[[pc]]) && !is.null(df[[mc]]))
                        round(df[[pc]] / df[[mc]] * 100, 1) else NULL
                    rng_pct <- if (!is.null(emp)) round(range(emp, na.rm = TRUE), 1) else c(NA, NA)
                    paste0("A disease prevalence survey at ", n_obs, " locations recorded '",
                           pc, "' positives out of '", mc, "' individuals examined. ",
                           "Empirical prevalence ranges from ", rng_pct[1], "% to ", rng_pct[2], "%. ",
                           "Write 2-3 paragraphs for a scientific report explaining what the map of ",
                           "observed prevalence shows and why model-based geostatistics is used.")
                },
                count = {
                    cc <- tryCatch(input$c, error = function(e) "counts")
                    ec <- tryCatch(input$e, error = function(e) "offset")
                    paste0("A geostatistical survey recorded case counts ('", cc, "') with an ",
                           "exposure offset ('", ec, "') at ", n_obs, " locations. ",
                           "Write 2-3 paragraphs for a scientific report explaining what the spatial ",
                           "map of counts shows and why a geostatistical model is appropriate.")
                },
                paste0("A geostatistical dataset with ", n_obs, " spatial observations. ",
                       "Write 2-3 paragraphs explaining what the map shows.")
            )
            llm_rv$data_text  <- call_groq(prompt, api_key, model_id)
            done <- done + 1
            llm_rv$status <- paste0("Generated ", done, " / ", total_tasks, " explanations...")
        }

        # --- Variogram ---
        if ("fig3" %in% what) {
            vs <- tryCatch(var_plot_sum(), error = function(e) NULL)
            if (!is.null(vs)) {
                summ_text <- if (!is.null(vs$summ))
                    paste(capture.output(vs$summ), collapse = "\n")
                else "Variogram summary not available."
                fn <- if (!is.null(input$functions)) input$functions else "Matern"
                prompt <- paste0(
                    "An empirical variogram was fitted using the '", fn, "' correlation function. ",
                    "Variogram output:\n\n", summ_text,
                    "\n\nWrite 2-3 paragraphs for a scientific report: interpret the spatial range ",
                    "(phi parameter, in km) and nugget effect in plain language for a public health audience, ",
                    "and explain what the degree of spatial correlation implies for the study area."
                )
                llm_rv$variog_text <- call_groq(prompt, api_key, model_id)
                done <- done + 1
                llm_rv$status <- paste0("Generated ", done, " / ", total_tasks, " explanations...")
            }
        }

        # --- Parameter estimates ---
        if ("fig4" %in% what) {
            fit <- tryCatch(model.fit(), error = function(e) NULL)
            if (!is.null(fit)) {
                summ_text <- paste(capture.output(summary(fit)), collapse = "\n")
                datatype  <- input$datatype
                prompt    <- paste0(
                    "A model-based geostatistical model for ", datatype, " data was fitted. ",
                    "Model summary:\n\n", summ_text,
                    "\n\nWrite 2-3 paragraphs for a scientific report: interpret the regression ",
                    "coefficients, spatial parameters (range phi, partial sill sigma2), and nugget ",
                    "effect using plain language for a public health or environmental science audience."
                )
                llm_rv$est_text <- call_groq(prompt, api_key, model_id)
                done <- done + 1
                llm_rv$status <- paste0("Generated ", done, " / ", total_tasks, " explanations...")
            }
        }

        # --- Prediction ---
        if ("fig5" %in% what) {
            pf <- tryCatch(pred.fit(), error = function(e) NULL)
            if (!is.null(pf)) {
                vals  <- pf[, -c(1, 2), drop = FALSE]
                means <- rowMeans(vals, na.rm = TRUE)
                mn    <- round(mean(means, na.rm = TRUE), 4)
                rng   <- round(range(means, na.rm = TRUE), 4)
                datatype <- input$datatype
                map_what <- switch(datatype,
                    continuous = tryCatch(input$predtomapcont, error = function(e) "mean"),
                    prevalence = tryCatch(input$predtomapprev, error = function(e) "mean prevalence"),
                    count      = tryCatch(input$predtomapcount, error = function(e) "mean rate"),
                    "predicted surface"
                )
                outcome_label <- switch(datatype,
                    prevalence = "prevalence", count = "disease rate", "outcome")
                prompt <- paste0(
                    "A spatial prediction map (", map_what, ") was generated from a ",
                    datatype, " geostatistical model. The predicted ", outcome_label,
                    " values range from ", rng[1], " to ", rng[2],
                    " (spatial mean: ", mn, "). ",
                    "Write 2-3 paragraphs for a scientific report: describe the spatial pattern ",
                    "visible in the prediction map, highlight areas of high or low ", outcome_label,
                    ", and discuss practical or public health implications."
                )
                llm_rv$pred_text <- call_groq(prompt, api_key, model_id)
                done <- done + 1
                llm_rv$status <- paste0("Generated ", done, " / ", total_tasks, " explanations...")
            }
        }

        llm_rv$status <- paste0(
            "Done — ", done, " AI explanation(s) generated. ",
            "Review and edit the text below, then download the report."
        )
    })
    # ---- end LLM section ----

    output$model_selector_ui <- renderUI({
        choices <- c("RiskMap (MCMC)" = "riskmap")
        if(has_inla) {
            choices <- c(choices, "INLA (fast Bayes)" = "inla")
        } else {
            choices <- c(choices, "INLA — not installed" = "inla_na")
        }
        choices <- c(choices, "Stan — coming soon" = "stan")
        radioButtons("backend", "Fitting method:",
                     choices  = choices,
                     selected = "riskmap",
                     inline   = FALSE)
    })

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
            shinyjs::show(id="nl_terms")
            shinyjs::show(id="showNL")
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
            shinyjs::show(id="nl_terms")
            shinyjs::hide(id="showNL")
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
            shinyjs::show(id="nl_terms")
            shinyjs::hide(id="showNL")
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
            shinyjs::hide(id="nl_terms")
            shinyjs::hide(id="showNL")
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
            shinyjs::hide(id="nl_terms")
            shinyjs::hide(id="showNL")
        }
    })
    # Upload the data
    data_all <- reactive({
        req(input$mbgdata)
        dff <- input$mbgdata
        if (is.null(dff))
            return(NULL)
        if(grepl("\\.rds$", dff$name)){
            x <- as.data.frame(readRDS(dff$datapath))
            x
        }else{
            x <- as.data.frame(read_csv(dff$datapath, show_col_types=FALSE))
            # NB: do not inject helper columns here. Columns of data_all() feed the
            # variable-selection dropdowns (updateVarSelectInput), so extras would
            # show up as selectable variables. emplogit is computed in model.fit().
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
        shp_path <- paste(uploaddirectory, shpdf$name[grep(pattern="*.shp$", shpdf$name)], sep="/")
        map <- st_transform(sf::st_read(shp_path, quiet=TRUE), crs=4326)
        # map <- st_as_sf(spTransform(map, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))

        map

    })


    # Upload the grid.locations
    gridpred <- reactive({
        req(input$gridpreddata)
        dff <- input$gridpreddata
        if (is.null(dff))
            return(NULL)
        if(grepl("\\.rds$", dff$name)){
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
        if(grepl("\\.rds$", dff$name)){
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
        # updateTextInput(session, "nl_terms", value = paste(input$D, collapse = "+"))
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
        req(input$xaxis, input$yaxis, nchar(input$xaxis) > 0, nchar(input$yaxis) > 0)

        make_explore_leaflet <- function(mapdata_wgs84, fill_vals, legend_title, shp_wgs84=NULL) {
            pal <- colorNumeric("RdYlBu", fill_vals, reverse=TRUE, na.color="transparent")
            coords_ll <- st_coordinates(mapdata_wgs84)
            m <- leaflet() %>%
                addProviderTiles("CartoDB.Positron") %>%
                addCircleMarkers(
                    lng=coords_ll[,1], lat=coords_ll[,2],
                    radius=6, color=pal(fill_vals),
                    fillOpacity=0.8, stroke=FALSE,
                    popup=paste0(legend_title, ": ", round(fill_vals, 4))
                ) %>%
                addLegend("bottomright", pal=pal, values=fill_vals,
                          title=legend_title, labFormat=labelFormat(digits=3))
            if(!is.null(shp_wgs84)) {
                m <- m %>% addPolylines(data=shp_wgs84, color="black", weight=1)
            }
            m
        }

        if(input$datatype == 'continuous'){
            req(input$y, nchar(input$y) > 0)
            crs_use <- if(!is.null(input$mbgshp)) st_crs(map_all()) else st_crs(as.integer(input$crs))
            mapdata <- st_transform(st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=crs_use), 4326)
            shp_wgs84 <- if(!is.null(input$mbgshp)) st_transform(map_all(), 4326) else NULL
            make_explore_leaflet(mapdata, df[, input$y], input$y, shp_wgs84)

        } else if(input$datatype == 'prevalence'){
            req(input$p, input$m, nchar(input$p) > 0, nchar(input$m) > 0)
            crs_use <- if(!is.null(input$mbgshp)) st_crs(map_all()) else st_crs(as.integer(input$crs))
            mapdata <- st_transform(st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=crs_use), 4326)
            prev_vals <- df[, input$p] / df[, input$m]
            shp_wgs84 <- if(!is.null(input$mbgshp)) st_transform(map_all(), 4326) else NULL
            make_explore_leaflet(mapdata, prev_vals, "Empirical prevalence", shp_wgs84)

        } else {
            req(input$c, input$e, nchar(input$c) > 0, nchar(input$e) > 0)
            crs_use <- if(!is.null(input$mbgshp)) st_crs(map_all()) else st_crs(as.integer(input$crs))
            mapdata <- st_transform(st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=crs_use), 4326)
            inc_vals <- df[, input$c] / df[, input$e]
            shp_wgs84 <- if(!is.null(input$mbgshp)) st_transform(map_all(), 4326) else NULL
            make_explore_leaflet(mapdata, inc_vals, "Incidence", shp_wgs84)
        }
    })


    output$map <- renderLeaflet({
        if (is.null(explore_map_lf())) return(NULL)
        explore_map_lf()
    })

    ################### mapping without interactive map ###################################

    explore_map_st <- reactive({
        df <- data_all()
        req(input$xaxis, input$yaxis, nchar(input$xaxis) > 0, nchar(input$yaxis) > 0)

        make_explore_ggplot <- function(mapdata_sf, fill_col, legend_title, shp=NULL) {
            p <- ggplot() +
                geom_sf(data=mapdata_sf, aes(color=.data[[fill_col]]), size=2, alpha=0.8) +
                scale_color_distiller(palette="RdYlBu", direction=1, name=legend_title) +
                theme_bw() +
                theme(legend.position="right")
            if(!is.null(shp)) {
                p <- p + geom_sf(data=shp, fill=NA, color="black", linewidth=0.4)
            }
            p
        }

        if(input$datatype == 'continuous'){
            req(input$y, nchar(input$y) > 0)
            mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
            shp <- if(!is.null(input$mbgshp)) map_all() else NULL
            make_explore_ggplot(mapdata, input$y, input$y, shp)

        } else if(input$datatype == 'prevalence'){
            req(input$p, input$m, nchar(input$p) > 0, nchar(input$m) > 0)
            mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
            mapdata[["Prevalence"]] <- df[, input$p] / df[, input$m]
            shp <- if(!is.null(input$mbgshp)) map_all() else NULL
            make_explore_ggplot(mapdata, "Prevalence", "Empirical\nprevalence", shp)

        } else {
            req(input$c, input$e, nchar(input$c) > 0, nchar(input$e) > 0)
            mapdata <- st_as_sf(df, coords=c(input$xaxis, input$yaxis))
            mapdata[["incidence"]] <- df[, input$c] / df[, input$e]
            shp <- if(!is.null(input$mbgshp)) map_all() else NULL
            make_explore_ggplot(mapdata, "incidence", "Incidence", shp)
        }
    })

    output$map2 <- renderPlot({
        if (is.null(explore_map_st())) return(NULL)
        explore_map_st()
    })


    ######### Plotting the relationship between the outcomes and the covariates

    scatter_ass_plot <- reactive({
        df <- data_all()
        req(input$xaxis, nchar(input$xaxis) > 0)
        if(is.null(input$D) || length(input$D) == 0) return(NULL)

        func <- switch(input$transformcov,
                       log=log,
                       sqrt=sqrt,
                       identity)


        if(input$showNL) f0 <- create_formula(y="y", covars = strsplit(input$nl_terms, "\\+")[[1]])


        if(input$datatype=='continuous'){
            new_dat <- data.frame(df[, c(input$y, input$D), drop=FALSE])
            if (input$transformcont == "log"){
                toExclude <- names(new_dat)[1]
                new_dat[,toExclude] <- log(new_dat[,toExclude])
                new_dat2 <- pivot_longer(new_dat, cols = -all_of(toExclude), names_to = "key", values_to = "value")
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = toExclude)) +
                    geom_point() +
                    facet_wrap(facets = ~key, scales = "free_x") +
                    geom_smooth(se=F) +
                    labs(x="", y=paste0("Log-", toExclude))
                if(input$showNL){
                    if(length(attr(terms(f0),"term.labels"))>1){
                        # f0 <- create_formula(y="y", covars = strsplit(nl_terms, "\\+")[[1]])
                        ff <- drop_formula_term(the_formula = f0, var_name = all.vars(f0[[3]]))
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = ff[[x$key[1]]], method="lm", col="red"))
                    }else{
                        ff <- gsub(pattern = paste0("\\b", all.vars(f0[[3]]), "\\b"), replacement = "x",  x=attr(terms(f0),"term.labels")[1])
                        fff <- create_formula("y", ff)
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = fff, method="lm", col="red"))
                    }
                    pp <- pp + p_smooth
                }
                pp
            }else{
                toExclude <- names(new_dat)[1]
                # new_dat[,toExclude] <- func(new_dat[,toExclude])
                new_dat2 <- pivot_longer(new_dat, cols = -all_of(toExclude), names_to = "key", values_to = "value")
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = toExclude)) +
                    geom_point() +
                    facet_wrap(facets = ~key, scales = "free_x") +
                    geom_smooth(se=F) +
                    labs(x="")
                if(input$showNL){
                    if(length(attr(terms(f0),"term.labels"))>1){
                        # f0 <- create_formula(y="y", covars = strsplit(nl_terms, "\\+")[[1]])
                        ff <- drop_formula_term(the_formula = f0, var_name = all.vars(f0[[3]]))
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = ff[[x$key[1]]], method="lm", col="red"))
                    }else{
                        ff <- gsub(pattern = paste0("\\b", all.vars(f0[[3]]), "\\b"), replacement = "x",  x=attr(terms(f0),"term.labels")[1])
                        fff <- create_formula("y", ff)
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = fff, method="lm", col="red"))
                    }
                    pp <- pp + p_smooth
                }
                pp

            }
        }else if (input$datatype=='prevalence'){

            if(input$transformprev == "logit"){
                new_dat <- data.frame(df[, c(input$p, input$m, input$D), drop=FALSE])
                new_dat[,"Emplogit"] <- log((new_dat[,input$p] + 0.5)/(new_dat[, input$m] - new_dat[,input$p] + 0.5))
                new_dat2 <- pivot_longer(new_dat[, -c(1,2), drop=FALSE], cols = -Emplogit, names_to = "key", values_to = "value")
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "Emplogit")) +
                    geom_point() +
                    facet_wrap(facets = ~key, scales = "free_x") +
                    geom_smooth(se=F) +
                    labs(x="", y=paste0("Emp-logit prevalence"))
                if(input$showNL){
                    if(length(attr(terms(f0),"term.labels"))>1){
                        # f0 <- create_formula(y="y", covars = strsplit(nl_terms, "\\+")[[1]])
                        ff <- drop_formula_term(the_formula = f0, var_name = all.vars(f0[[3]]))
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = ff[[x$key[1]]], method="lm", col="red"))
                    }else{
                        ff <- gsub(pattern = paste0("\\b", all.vars(f0[[3]]), "\\b"), replacement = "x",  x=attr(terms(f0),"term.labels")[1])
                        fff <- create_formula("y", ff)
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = fff, method="lm", col="red"))
                    }
                    pp <- pp + p_smooth
                }
                pp
            }else if (input$transformprev == "log"){
                new_dat <- data.frame(df[, c(input$p, input$m, input$D), drop=FALSE])
                new_dat[,"logprev"] <- log((new_dat[,input$p])/(new_dat[, input$m]))
                new_dat2 <- pivot_longer(new_dat[, -c(1,2), drop=FALSE], cols = -logprev, names_to = "key", values_to = "value")
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "logprev")) +
                    geom_point() +
                    facet_wrap(facets = ~key, scales = "free_x") +
                    geom_smooth(se=F) +
                    labs(x="", y=paste0("Log-prevalence"))
                if(input$showNL){
                    if(length(attr(terms(f0),"term.labels"))>1){
                        # f0 <- create_formula(y="y", covars = strsplit(nl_terms, "\\+")[[1]])
                        ff <- drop_formula_term(the_formula = f0, var_name = all.vars(f0[[3]]))
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = ff[[x$key[1]]], method="lm", col="red"))
                    }else{
                        ff <- gsub(pattern = paste0("\\b", all.vars(f0[[3]]), "\\b"), replacement = "x",  x=attr(terms(f0),"term.labels")[1])
                        fff <- create_formula("y", ff)
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = fff, method="lm", col="red"))
                    }
                    pp <- pp + p_smooth
                }
                pp
            }else{
                new_dat <- data.frame(df[, c(input$p, input$m, input$D), drop=FALSE])
                new_dat[,"pprev"] <- as.numeric((new_dat[,input$p])/(new_dat[, input$m]))
                new_dat2 <- pivot_longer(new_dat[, -c(1,2), drop=FALSE], cols = -pprev, names_to = "key", values_to = "value")
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "pprev")) +
                     geom_point() +
                    facet_wrap(facets = ~key, scales = "free_x") +
                    geom_smooth(se=F) +
                    labs(x="", y=paste0("Prevalence"))
                if(input$showNL){
                    if(length(attr(terms(f0),"term.labels"))>1){
                        # f0 <- create_formula(y="y", covars = strsplit(nl_terms, "\\+")[[1]])
                        ff <- drop_formula_term(the_formula = f0, var_name = all.vars(f0[[3]]))
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = ff[[x$key[1]]], method="lm", col="red"))
                    }else{
                        ff <- gsub(pattern = paste0("\\b", all.vars(f0[[3]]), "\\b"), replacement = "x",  x=attr(terms(f0),"term.labels")[1])
                        fff <- create_formula("y", ff)
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = fff, method="lm", col="red"))
                    }
                    pp <- pp + p_smooth
                }
                pp

            }
        } else{

            if (input$transformcnt == "log"){
                new_dat <- data.frame(df[, c(input$c, input$e, input$D), drop=FALSE])
                new_dat[,"logincidence"] <- log((new_dat[,input$c])/(new_dat[, input$e]))
                new_dat2 <- pivot_longer(new_dat[, -c(1,2), drop=FALSE], cols = -logincidence, names_to = "key", values_to = "value")
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "logincidence")) +
                    geom_point() +
                    facet_wrap(facets = ~key, scales = "free_x") +
                    geom_smooth(se=F) +
                    labs(x="", y=paste0("Log-incidence"))
                if(input$showNL){
                    if(length(attr(terms(f0),"term.labels"))>1){
                        # f0 <- create_formula(y="y", covars = strsplit(nl_terms, "\\+")[[1]])
                        ff <- drop_formula_term(the_formula = f0, var_name = all.vars(f0[[3]]))
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = ff[[x$key[1]]], method="lm", col="red"))
                    }else{
                        ff <- gsub(pattern = paste0("\\b", all.vars(f0[[3]]), "\\b"), replacement = "x",  x=attr(terms(f0),"term.labels")[1])
                        fff <- create_formula("y", ff)
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = fff, method="lm", col="red"))
                    }
                    pp <- pp + p_smooth
                }
                pp
            }else{
                new_dat <- data.frame(df[, c(input$c, input$e, input$D), drop=FALSE])
                new_dat[,"iincidence"] <- (new_dat[,input$c])/(new_dat[, input$e])
                new_dat2 <- pivot_longer(new_dat[, -c(1,2), drop=FALSE], cols = -iincidence, names_to = "key", values_to = "value")
                new_dat2[,names(new_dat2)[3]] <- func(new_dat2[,names(new_dat2)[3]])
                pp <- ggplot(new_dat2, aes_string(x = names(new_dat2)[3], y = "iincidence")) +
                    geom_point() +
                    facet_wrap(facets = ~key, scales = "free_x") +
                    geom_smooth(se=F) +
                    labs(x="", y=paste0("Incidence"))
                if(input$showNL){
                    if(length(attr(terms(f0),"term.labels"))>1){
                        # f0 <- create_formula(y="y", covars = strsplit(nl_terms, "\\+")[[1]])
                        ff <- drop_formula_term(the_formula = f0, var_name = all.vars(f0[[3]]))
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = ff[[x$key[1]]], method="lm", col="red"))
                    }else{
                        ff <- gsub(pattern = paste0("\\b", all.vars(f0[[3]]), "\\b"), replacement = "x",  x=attr(terms(f0),"term.labels")[1])
                        fff <- create_formula("y", ff)
                        new_dat3 <- subset(new_dat2, key %in% all.vars(f0[[3]]))
                        p_smooth <- by(new_dat3, new_dat3$key,
                                       function(x) geom_smooth(data=x, se=F, formula = fff, method="lm", col="red"))
                    }
                    pp <- pp + p_smooth
                }
                pp

            }
        }

    })


    output$Plot <- renderPlot({
        if (is.null(scatter_ass_plot())) return(NULL)
        scatter_ass_plot() + theme_bw()
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
                # y <- input$y




                fml <- return_formula(y=input$y, covars = input$D, nl_terms = input$nl_terms)

                # fml <- as.formula(paste(paste0(input$y, " ~ ", paste(input$D, collapse= "+"))))
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




                fml <- return_formula(y=input$p, covars = input$D, nl_terms = input$nl_terms)

                # xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                m <- model.frame(fml, df)
                xmat <- model.matrix(fml, m)
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

                fml <- return_formula(y=input$c, covars = input$D, nl_terms = input$nl_terms)

                # xmat <- as.matrix(cbind(1, df[, input$D, drop=FALSE]))
                m <- model.frame(fml, df)
                xmat <- model.matrix(fml, m)
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
        var_plot_sum()$pl + theme_bw()
    })

    output$summary <- renderPrint({
        req(!is.null(var_plot_sum()), !is.null(var_plot_sum()$summ))
        var_plot_sum()$summ
    })


    model.fit <- reactive({
        withProgress(message="Fitting model...", value=0, {
            # Guard: Stan not implemented; INLA requires the package
            backend_sel <- if(!is.null(input$backend)) input$backend else "riskmap"
            if(backend_sel == "stan") {
                showNotification("Stan backend is not yet implemented. Please select RiskMap or INLA.",
                                 type="warning", duration=6)
                return(NULL)
            }
            if(backend_sel == "inla_na") {
                showNotification("INLA is not installed. Install it from r-inla.org, then restart the app.",
                                 type="error", duration=8)
                return(NULL)
            }

            df <- data_all()  # already a plain data.frame

            # UTM code for projection
            if(input$maptype == 'view') {
                coords_tmp  <- df[, c(input$xaxis, input$yaxis), drop=FALSE]
                utmcode_int <- as.integer(lonlat2UTM(as.numeric(coords_tmp[1, ])))
                input_crs   <- as.integer(input$crs)
            } else {
                utmcode_int <- NULL
                input_crs   <- NULL
            }

            nugget_val <- if(input$includenugget == 1) input$nu else 0

            build_gp_formula <- function(fml_base) {
                as.formula(paste0(
                    deparse(fml_base),
                    " + gp(", input$xaxis, ", ", input$yaxis,
                    ", kappa=", input$kappa, ", nugget=", nugget_val, ")"
                ))
            }

            # glgpm uses deparse(substitute(den)) so den must be a bare symbol.
            # Keep `data` as a symbol too: passing the literal data frame embeds
            # it in the stored call (res$call <- match.call()), which summary()
            # and print() then dump in full. Binding it in a local environment
            # keeps the printed Call as `data = data`.
            call_glgpm <- function(...) {
                args <- list(...)
                env  <- new.env(parent = parent.frame())
                env$data  <- args$data
                args$data <- as.name("data")
                cl <- as.call(c(list(as.name("glgpm")), args))
                eval(cl, envir = env)
            }

            use_inla <- has_inla && backend_sel == "inla"

            if(use_inla) {
                # ---- INLA path ----
                incProgress(0.1, message="Building INLA mesh...")

                # Coordinates in model projection
                coords_sf <- st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=as.integer(input$crs))
                if(isTRUE(input$maptype == 'view') && !is.null(utmcode_int)) {
                    coords_mat <- st_coordinates(st_transform(coords_sf, crs=utmcode_int))
                } else {
                    coords_mat <- st_coordinates(coords_sf)
                }

                phi_val <- input$phi
                mbg_mesh <- INLA::inla.mesh.2d(
                    loc      = coords_mat,
                    max.edge = c(phi_val * 0.5, phi_val * 2),
                    cutoff   = max(phi_val * 0.05, diff(range(coords_mat[, 1])) / 100)
                )
                mbg_spde <- INLA::inla.spde2.matern(mesh=mbg_mesh, alpha=2)
                A_fit    <- INLA::inla.spde.make.A(mesh=mbg_mesh, loc=coords_mat)

                # Determine response variable and build base formula
                resp_var <- switch(input$datatype,
                    continuous = input$y,
                    prevalence = if(input$fitlinear == "linearmodel") "emplogit" else input$p,
                    count      = input$c
                )
                if(input$datatype == "prevalence" && input$fitlinear == "linearmodel") {
                    df[["emplogit"]] <- log((df[[input$p]] + 0.5) / (df[[input$m]] - df[[input$p]] + 0.5))
                }
                fml_base <- if(is.null(input$D)) {
                    as.formula(paste0(resp_var, " ~ 1"))
                } else {
                    return_formula(y=resp_var, covars=input$D, nl_terms=input$nl_terms)
                }

                # Design matrix
                X_fit <- model.matrix(update(fml_base, NULL ~ .), data=df)
                colnames(X_fit)[colnames(X_fit) == "(Intercept)"] <- "Intercept"
                covar_names <- colnames(X_fit)
                X_df        <- as.data.frame(X_fit)

                idx_spde    <- INLA::inla.spde.make.index("spatial_field", n.spde=mbg_spde$n.spde)
                inla_fml_str <- paste0("y_inla ~ -1 + ", paste(covar_names, collapse=" + "),
                                       " + f(spatial_field, model=mbg_spde)")
                inla_fml <- as.formula(inla_fml_str)

                incProgress(0.3, message="Running INLA...")

                if(input$datatype == "continuous" || (input$datatype == "prevalence" && input$fitlinear == "linearmodel")) {
                    y_inla    <- df[[resp_var]]
                    stack_fit <- INLA::inla.stack(
                        data=list(y_inla=y_inla), A=list(A_fit, 1),
                        effects=list(idx_spde, X_df), tag="fit"
                    )
                    inla_fit <- INLA::inla(
                        inla_fml, family="gaussian",
                        data=INLA::inla.stack.data(stack_fit),
                        control.predictor=list(A=INLA::inla.stack.A(stack_fit), compute=FALSE),
                        control.compute=list(config=TRUE), verbose=FALSE
                    )

                } else if(input$datatype == "prevalence") {
                    y_inla    <- df[[input$p]]
                    Ntrials   <- df[[input$m]]
                    stack_fit <- INLA::inla.stack(
                        data=list(y_inla=y_inla, Ntrials=Ntrials), A=list(A_fit, 1),
                        effects=list(idx_spde, X_df), tag="fit"
                    )
                    inla_fit <- INLA::inla(
                        inla_fml, family="binomial",
                        Ntrials=INLA::inla.stack.data(stack_fit)$Ntrials,
                        data=INLA::inla.stack.data(stack_fit),
                        control.predictor=list(A=INLA::inla.stack.A(stack_fit), compute=FALSE),
                        control.compute=list(config=TRUE), verbose=FALSE
                    )

                } else {  # count / Poisson
                    y_inla    <- df[[input$c]]
                    X_df$log_offset <- log(pmax(df[[input$e]], 1e-8))
                    stack_fit <- INLA::inla.stack(
                        data=list(y_inla=y_inla), A=list(A_fit, 1),
                        effects=list(idx_spde, X_df), tag="fit"
                    )
                    inla_fml_count <- as.formula(paste0(
                        "y_inla ~ -1 + ", paste(covar_names, collapse=" + "),
                        " + offset(log_offset) + f(spatial_field, model=mbg_spde)"
                    ))
                    inla_fit <- INLA::inla(
                        inla_fml_count, family="poisson",
                        data=INLA::inla.stack.data(stack_fit),
                        control.predictor=list(A=INLA::inla.stack.A(stack_fit), compute=FALSE),
                        control.compute=list(config=TRUE), verbose=FALSE
                    )
                }

                incProgress(0.2, message="Drawing posterior samples...")
                n_post   <- min(if(!is.null(input$mcmcNsim)) input$mcmcNsim else 1000, 2000)
                post_samp <- INLA::inla.posterior.sample(n_post, inla_fit)

                inla_fit$app_backend     <- "inla"
                inla_fit$app_mesh        <- mbg_mesh
                inla_fit$app_spde        <- mbg_spde
                inla_fit$app_utmcode     <- utmcode_int
                inla_fit$app_fml         <- fml_base
                inla_fit$app_X_fit       <- X_fit
                inla_fit$app_covar_names <- covar_names
                inla_fit$app_coords_utm  <- coords_mat
                inla_fit$app_post_samp   <- post_samp
                inla_fit$app_datatype    <- input$datatype
                inla_fit$app_stack       <- stack_fit
                incProgress(0.1, message="Done.")
                inla_fit

            } else {
                # ---- RiskMap MCMC path ----
                incProgress(0.1, message="Preparing model...")

                if(input$datatype == 'continuous') {
                    fml_base <- if(is.null(input$D)) {
                        as.formula(paste0(input$y, " ~ 1"))
                    } else {
                        return_formula(y=input$y, covars=input$D, nl_terms=input$nl_terms)
                    }
                    incProgress(0.3, message="Running MCMC (may take a few minutes)...")
                    fit <- glgpm(
                        formula       = build_gp_formula(fml_base),
                        data          = df, family="gaussian",
                        crs=input_crs, convert_to_crs=utmcode_int,
                        scale_to_km=(input$maptype=='view'),
                        control_mcmc=set_control_sim(n_sim=1000, linear_model=TRUE),
                        start_pars=list(phi=input$phi), messages=FALSE
                    )
                    fit$app_backend <- "riskmap"
                    fit$app_fml     <- fml_base
                    fit$app_utmcode <- utmcode_int
                    incProgress(0.6, message="Done.")
                    fit

                } else if(input$datatype == 'prevalence') {
                    if(input$fitlinear == "binomialmodel") {
                        fml_base <- if(is.null(input$D)) {
                            as.formula(paste0(input$p, " ~ 1"))
                        } else {
                            return_formula(y=input$p, covars=input$D, nl_terms=input$nl_terms)
                        }
                        incProgress(0.3, message="Running MCMC (may take several minutes)...")
                        fit <- call_glgpm(
                            formula=build_gp_formula(fml_base), data=df, family="binomial",
                            den=as.name(input$m), crs=input_crs, convert_to_crs=utmcode_int,
                            scale_to_km=(input$maptype=='view'),
                            control_mcmc=set_control_sim(n_sim=input$mcmcNsim, burnin=input$mcmcNburn, thin=input$mcmcNthin),
                            start_pars=list(phi=input$phi), return_samples=TRUE, messages=FALSE
                        )
                        fit$app_backend <- "riskmap"
                        fit$app_fml     <- fml_base
                        fit$app_utmcode <- utmcode_int
                        incProgress(0.6, message="Done.")
                        fit

                    } else {  # linearmodel
                        df[, "emplogit"] <- log((df[, input$p] + 0.5) / (df[, input$m] - df[, input$p] + 0.5))
                        fml_base <- if(is.null(input$D)) {
                            as.formula("emplogit ~ 1")
                        } else {
                            return_formula(y="emplogit", covars=input$D, nl_terms=input$nl_terms)
                        }
                        incProgress(0.3, message="Running MCMC...")
                        fit <- glgpm(
                            formula=build_gp_formula(fml_base), data=df, family="gaussian",
                            crs=input_crs, convert_to_crs=utmcode_int,
                            scale_to_km=(input$maptype=='view'),
                            control_mcmc=set_control_sim(n_sim=1000, linear_model=TRUE),
                            start_pars=list(phi=input$phi), messages=FALSE
                        )
                        fit$app_backend <- "riskmap"
                        fit$app_fml     <- fml_base
                        fit$app_utmcode <- utmcode_int
                        incProgress(0.6, message="Done.")
                        fit
                    }

                } else {  # count / Poisson
                    fml_base <- if(is.null(input$D)) {
                        as.formula(paste0(input$c, " ~ 1"))
                    } else {
                        return_formula(y=input$c, covars=input$D, nl_terms=input$nl_terms)
                    }
                    incProgress(0.3, message="Running MCMC (may take several minutes)...")
                    fit <- call_glgpm(
                        formula=build_gp_formula(fml_base), data=df, family="poisson",
                        den=as.name(input$e), crs=input_crs, convert_to_crs=utmcode_int,
                        scale_to_km=(input$maptype=='view'),
                        control_mcmc=set_control_sim(n_sim=input$mcmcNsim, burnin=input$mcmcNburn, thin=input$mcmcNthin),
                        start_pars=list(phi=input$phi), return_samples=TRUE, messages=FALSE
                    )
                    fit$app_backend <- "riskmap"
                    fit$app_fml     <- fml_base
                    fit$app_utmcode <- utmcode_int
                    incProgress(0.6, message="Done.")
                    fit
                }
            }
        })
    }) |>
        bindCache(
            input$mbgdata$datapath,
            input$datatype, input$fitlinear,
            input$p, input$m, input$y, input$c, input$e,
            input$xaxis, input$yaxis, input$crs, input$maptype,
            input$phi, input$kappa, input$includenugget, input$nu,
            input$mcmcNsim, input$mcmcNburn, input$mcmcNthin,
            paste0(sort(input$D), collapse=","), input$nl_terms,
            if(!is.null(input$backend)) input$backend else "riskmap",
            cache = "session"
        ) |>
        bindEvent(input$ShowEst, ignoreNULL=TRUE, ignoreInit=TRUE)

    output$estsummary <- renderPrint({
        if (is.null(model.fit())) return(NULL)
        summary(model.fit())
    })

    output$tab <- renderTable({
        if (is.null(model.fit())) return(NULL)
        tab <- as.data.frame(to_table(model.fit()))
        # to_table() carries parameter names as row names; renderTable() drops
        # row names by default, so surface them as an explicit first column.
        tab <- cbind(Parameter = rownames(tab), tab)
        rownames(tab) <- NULL
        tab
    }, striped=TRUE, hover=TRUE, bordered=TRUE)


    pred.fit <- eventReactive(input$ShowPred, {
        withProgress(message="Computing spatial predictions...", value=0, {
            req(model.fit())
            fit         <- model.fit()
            df          <- data_all()
            utmcode_int <- fit$app_utmcode
            view_mode   <- (input$maptype == 'view')

            incProgress(0.1, message="Building prediction grid...")

            # ---- Build prediction grid ----
            if(is.null(input$gridpreddata)) {
                if(is.null(input$mbgshp)) {
                    if(view_mode) {
                        hull_sf  <- convex_hull_sf(st_as_sf(df, coords=c(input$xaxis, input$yaxis), crs=as.integer(input$crs)))
                        grid_utm <- create_grid(hull_sf, spat_res=input$resolution, grid_crs=utmcode_int)
                        grid_sfc <- st_geometry(st_transform(grid_utm, crs=as.integer(input$crs)))
                    } else {
                        hull_sfc <- st_convex_hull(st_union(st_as_sf(df, coords=c(input$xaxis, input$yaxis))))
                        grid_sfc <- st_geometry(st_make_grid(hull_sfc, cellsize=input$resolution, what="centers"))
                        grid_utm <- NULL
                    }
                } else {
                    shp <- map_all()
                    if(view_mode) {
                        hull_sf  <- st_transform(shp, crs=as.integer(input$crs))
                        grid_utm <- create_grid(hull_sf, spat_res=input$resolution, grid_crs=utmcode_int)
                        grid_sfc <- st_geometry(st_transform(grid_utm, crs=as.integer(input$crs)))
                    } else {
                        grid_sfc <- st_geometry(st_make_grid(shp, cellsize=input$resolution, what="centers"))
                        grid_utm <- NULL
                    }
                }
            } else {
                gridpred_df <- gridpred()
                if(view_mode) {
                    grid_sfc <- st_geometry(st_as_sf(gridpred_df, coords=c(1, 2), crs=as.integer(input$crs)))
                } else {
                    grid_sfc <- st_geometry(st_as_sf(gridpred_df, coords=c(1, 2)))
                }
                grid_utm <- NULL
            }

            # ---- Predictors ----
            pred_vars <- if(is.null(input$predictorsdata)) NULL else data.frame(predictors())

            # ---- INLA prediction path ----
            if(!is.null(fit$app_backend) && fit$app_backend == "inla") {
                incProgress(0.2, message="Projecting INLA posterior to grid...")

                mbg_mesh    <- fit$app_mesh
                post_samp   <- fit$app_post_samp
                covar_names <- fit$app_covar_names

                # Prediction coordinates in model CRS (UTM or original)
                if(isTRUE(view_mode) && !is.null(utmcode_int)) {
                    grid_sf_tmp      <- st_set_crs(st_sf(geometry=grid_sfc), as.integer(input$crs))
                    grid_pred_coords <- st_coordinates(st_transform(grid_sf_tmp, crs=utmcode_int))
                } else {
                    grid_pred_coords <- st_coordinates(grid_sfc)
                }
                n_pred <- nrow(grid_pred_coords)

                # Projection matrix from mesh to prediction locations
                A_pred <- INLA::inla.spde.make.A(mesh=mbg_mesh, loc=grid_pred_coords)

                # Design matrix at prediction locations
                if(!is.null(pred_vars)) {
                    X_pred_mat <- model.matrix(update(fit$app_fml, NULL ~ .), data=pred_vars)
                    colnames(X_pred_mat)[colnames(X_pred_mat) == "(Intercept)"] <- "Intercept"
                } else {
                    X_pred_mat <- matrix(0, nrow=n_pred, ncol=length(covar_names),
                                         dimnames=list(NULL, covar_names))
                    X_pred_mat[, "Intercept"] <- 1
                }

                # Extract spatial field and fixed effects from posterior samples
                samp_names  <- rownames(post_samp[[1]]$latent)
                spatial_idx <- grep("^spatial_field:", samp_names)
                spatial_mat <- do.call(cbind, lapply(post_samp, function(s) s$latent[spatial_idx, 1]))

                beta_rows <- match(paste0(covar_names, ":1"), samp_names)
                beta_mat  <- do.call(cbind, lapply(post_samp, function(s) s$latent[beta_rows, 1]))

                # Linear predictor at prediction grid
                lp_mat <- as.matrix(A_pred %*% spatial_mat) + X_pred_mat %*% beta_mat

                # Inverse link
                response_mat <- switch(fit$app_datatype,
                    continuous = lp_mat,
                    prevalence = plogis(lp_mat),
                    count      = exp(lp_mat)
                )

                # Coordinates for raster (use UTM if available for regular grid)
                if(view_mode && !is.null(grid_utm)) {
                    coords_ras <- st_coordinates(grid_utm)[, 1:2]
                    coords_crs <- utmcode_int
                } else {
                    coords_ras <- grid_pred_coords
                    coords_crs <- if(view_mode) as.integer(input$crs) else NA_integer_
                }
                res_df <- data.frame(coords_ras, response_mat)
                attr(res_df, "utmcode_int") <- utmcode_int
                attr(res_df, "coords_crs")  <- coords_crs
                attr(res_df, "view_mode")   <- view_mode
                incProgress(0.5, message="Done.")
                return(res_df)
            }

            # ---- RiskMap MCMC prediction path ----
            incProgress(0.2, message="Running predictive simulations...")
            is_linear <- (input$datatype == 'continuous') ||
                         (input$datatype == 'prevalence' && input$fitlinear == 'linearmodel')
            ctrl_pred <- if(is_linear) {
                set_control_sim(n_sim=1000, linear_model=TRUE)
            } else {
                set_control_sim(n_sim=input$mcmcNsim, burnin=input$mcmcNburn, thin=input$mcmcNthin)
            }

            pred_re <- pred_over_grid(
                object=fit, grid_pred=grid_sfc, predictors=pred_vars,
                control_sim=ctrl_pred, type="marginal", messages=FALSE
            )
            pred_targets <- pred_target_grid(pred_re)
            lp           <- pred_targets$lp_samples

            response_samples <- switch(input$datatype,
                continuous = lp,
                prevalence = plogis(lp),
                count      = exp(lp)
            )

            if(view_mode && !is.null(grid_utm)) {
                coords_ras <- st_coordinates(grid_utm)[, 1:2]
                coords_crs <- utmcode_int
            } else {
                coords_ras <- st_coordinates(grid_sfc)[, 1:2]
                coords_crs <- if(view_mode) as.integer(input$crs) else NA_integer_
            }

            res_df <- data.frame(coords_ras, response_samples)
            attr(res_df, "utmcode_int") <- utmcode_int
            attr(res_df, "coords_crs")  <- coords_crs
            attr(res_df, "view_mode")   <- view_mode
            incProgress(0.5, message="Done.")
            res_df
        })
    })

    # ---- Shared helpers ----
    make_pred_raster <- function(ras_dff, coords_crs) {
        r <- terra::rast(ras_dff, type="xyz")
        if(!is.null(coords_crs) && !is.na(coords_crs)) {
            terra::crs(r) <- paste0("EPSG:", coords_crs)
        }
        r
    }

    # Compute raster value column from pred.fit output
    pred_value_raster <- reactive({
        if(is.null(pred.fit())) return(NULL)
        all_df      <- pred.fit()
        samples     <- all_df[, -c(1:2), drop=FALSE]
        utmcode_int <- attr(all_df, "utmcode_int")
        coords_crs  <- attr(all_df, "coords_crs")
        view_mode   <- attr(all_df, "view_mode")

        map_choice <- switch(input$datatype,
            continuous = input$predtomapcont,
            prevalence = input$predtomapprev,
            count      = input$predtomapcount
        )

        value_vec <- if(map_choice == "meann") {
            rowMeans(samples)
        } else if(map_choice == "sdd") {
            apply(samples, 1, sd)
        } else if(map_choice == "exprob") {
            rowMeans(samples > input$threshold)
        } else {
            apply(samples, 1, quantile, probs=input$quantprob)
        }

        ras_dff <- data.frame(all_df[, 1:2], value=value_vec)
        r       <- make_pred_raster(ras_dff, coords_crs)
        list(r=r, utmcode_int=utmcode_int, coords_crs=coords_crs, view_mode=view_mode)
    })

    pred_leaflet_map <- reactive({
        rv <- pred_value_raster()
        if(is.null(rv)) return(NULL)
        r          <- rv$r
        coords_crs <- rv$coords_crs
        view_mode  <- rv$view_mode

        map_choice <- switch(input$datatype,
            continuous = input$predtomapcont,
            prevalence = input$predtomapprev,
            count      = input$predtomapcount
        )
        map_title <- switch(map_choice,
            meann  = switch(input$datatype, continuous="Mean outcome", prevalence="Prevalence", count="Incidence"),
            sdd    = "Std. error",
            exprob = paste0("Ex-prob >", input$threshold * 100, "%"),
            quant  = paste0(input$quantprob * 100, "% Quantile")
        )

        # Project to WGS84 for leaflet (only if raster has a non-WGS84 CRS)
        needs_project <- !is.null(coords_crs) && !is.na(coords_crs) && coords_crs != 4326L
        if(needs_project) {
            r_wgs84 <- terra::project(r, "EPSG:4326")
        } else if(!is.null(coords_crs) && !is.na(coords_crs)) {
            r_wgs84 <- r  # already WGS84
        } else {
            terra::crs(r) <- "EPSG:4326"
            r_wgs84 <- r
        }
        r_stars <- stars::st_as_stars(r_wgs84)
        vals    <- na.omit(terra::values(r_wgs84))
        pal     <- colorNumeric("RdYlBu", vals, na.color="transparent", reverse=TRUE)

        leaflet() %>%
            addProviderTiles("CartoDB.Positron") %>%
            leafem::addStarsImage(r_stars, colors=pal, opacity=0.75, project=FALSE) %>%
            addLegend("bottomright", pal=pal, values=vals, title=map_title,
                      labFormat=labelFormat(digits=3))
    })

    pred_ggplot_map <- reactive({
        rv <- pred_value_raster()
        if(is.null(rv)) return(NULL)
        r <- rv$r

        map_choice <- switch(input$datatype,
            continuous = input$predtomapcont,
            prevalence = input$predtomapprev,
            count      = input$predtomapcount
        )
        map_title <- switch(map_choice,
            meann  = switch(input$datatype, continuous="Mean outcome", prevalence="Prevalence", count="Incidence"),
            sdd    = "Std. error",
            exprob = paste0("Ex-prob\n>", input$threshold * 100, "%"),
            quant  = paste0(input$quantprob * 100, "%\nQuantile")
        )

        ggplot() +
            tidyterra::geom_spatraster(data=r, aes(fill=value)) +
            scale_fill_distiller(palette="RdYlBu", direction=1, name=map_title, na.value=NA) +
            theme_bw() +
            theme(legend.position="right")
    })

    output$predmap <- renderLeaflet({
        if(is.null(pred_leaflet_map())) return(NULL)
        pred_leaflet_map()
    })

    output$predmap2 <- renderPlot({
        if(is.null(pred_ggplot_map())) return(NULL)
        pred_ggplot_map()
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
            parasumm = summary(model.fit())
        }else{
            parasumm = NULL
        }

        ##set which one to map
        if(any(input$whattoshow == 'fig5')){
            pred_map = pred_ggplot_map()
        }else{
            pred_map = NULL
        }


        # LLM text: use edited textarea values if available, else fall back to generated text
        llm_data  <- if (!is.null(input$edit_llm_data)  && nchar(trimws(input$edit_llm_data))  > 0) input$edit_llm_data  else NULL
        llm_variog <- if (!is.null(input$edit_llm_variog) && nchar(trimws(input$edit_llm_variog)) > 0) input$edit_llm_variog else NULL
        llm_est   <- if (!is.null(input$edit_llm_est)   && nchar(trimws(input$edit_llm_est))   > 0) input$edit_llm_est   else NULL
        llm_pred  <- if (!is.null(input$edit_llm_pred)  && nchar(trimws(input$edit_llm_pred))  > 0) input$edit_llm_pred  else NULL

        # correlation function for equations section
        cov_model  <- if (!is.null(input$functions)) input$functions else "matern"
        has_nugget <- if (!is.null(input$includenugget)) input$includenugget == 1 else TRUE
        kappa_val  <- if (!is.null(input$kappa)) input$kappa else 0.5
        fit_linear <- if (!is.null(input$fitlinear)) input$fitlinear else "binomialmodel"

        params <- list(
            nameofanalysis = input$datatype,
            cov_model      = cov_model,
            has_nugget     = has_nugget,
            kappa_val      = kappa_val,
            fit_linear     = fit_linear,
            exploremap     = exploremap,
            scatterplot    = scatterplot,
            varplot        = varplot,
            parasumm       = parasumm,
            predmap        = pred_map,
            llm_data       = llm_data,
            llm_variog     = llm_variog,
            llm_est        = llm_est,
            llm_pred       = llm_pred
        )
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
