# Set of functions used to create / handle / analyse models
##### Libraries ##### ---------------------------------------------------------
library(broom.mixed)
library(dotwhisker)
library(ggeffects)
library(ggpmisc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(abind)
library(readr)
library(Hmsc)
library(cli)
library(sf)

library(lme4)
library(lmerTest)


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters


##### Functions ##### ---------------------------------------------------------
# A function that takes parameters and makes the path to a folder's run.
# ARGS:
#   - folder: a string.
#   - combination: a list of parameters.
#   - k_fold: a numeric.
make_run_path <- function(folder, combination, k_fold) {
    # check that there is at MUST one parameter with a list of values
    if (sum(lapply(combination, length) > 1) > 1) {
        # more than one parameter with a list of values : 
        # -> I have not setup this so... This should definetly not happen.
        stop(paste(
            "Error in combination, found more than one parameter",
            "with a vector of values of length > 1."
        ))
    } else if (sum(lapply(combination, length) > 1) == 1) {
        # exatcly one, perfect: select the associated parameter's name
        loop_on <- names(which(lapply(combination, length) > 1))
        cli_alert_info(paste(
            "Auto-selection of parameter to loop on:", loop_on))
    } else {
        # none parameter with a list of value: assign a harmless one
        loop_on <- "TRAIN_SIZES"
    }


    # Every parameter has one element, except 'loop_on'
    # for each element in 'loop_on' the path is different
    run_paths <- c()
    for (loop_element in combination[[loop_on]]) {
        # For each new iteration, reset the parameters to use.
        local_com <- list()
        local_com$R_EFFECTS <- combination$R_EFFECTS
        local_com$STRATEGIES <- combination$STRATEGIES
        local_com$NEW_SAMPLE_SIZES <- combination$NEW_SAMPLE_SIZES
        local_com$TRAIN_SIZES <- combination$TRAIN_SIZES
        local_com$HMSC_XFORMULAS <- combination$HMSC_XFORMULAS

        # Overwrite parameter with multiple elements with one element 
        # (the one from the current iteration)
        # + check: if parameter is a formula, transform it to 
        # the number of variables in the formula
        local_com[[loop_on]] <- loop_element
            if (is.list(local_com$HMSC_XFORMULAS)) {
                if (is_formula(local_com$HMSC_XFORMULAS[[1]])) {
                    local_com$HMSC_XFORMULAS <- local_com$HMSC_XFORMULAS[[1]]
                } else {
                    stop("Function was not made to handle loop_element as list when not a formula.")
                }
            }
        
        # if strategy is none, overwrite n_new_samples to avoid errors
        if ((local_com$STRATEGIES == "none") & (local_com$NEW_SAMPLE_SIZES != 0)) {
            cli_alert_warning(paste0(
                "Strategy is 'none' but n_new_samples is '", 
                local_com$NEW_SAMPLE_SIZES, 
                "'. Overwriting with n_new_samples = '0' to avoid errors."))
            local_com$NEW_SAMPLE_SIZES <- 0
        }

        # if n_new_samples is 0, no need to apply a strategy, use "none"
        if ((local_com$NEW_SAMPLE_SIZES == 0) & (local_com$STRATEGIES != "none")) {
            cli_alert_warning(paste0(
                "n_new_samples is '0' but strategy is '", local_com$STRATEGIES, 
                "'. Overwriting with strategy = 'none' to avoid errors."))
            local_com$STRATEGIES <- 'none'
        }     
        
        # Create path string and add it to the list
        run_paths <- c(
            run_paths,
            file.path(
                folder, paste0(
                    "model_random-", local_com$R_EFFECTS,
                    "_strategy-", local_com$STRATEGIES,
                    "_new-samples-", local_com$NEW_SAMPLE_SIZES,
                    "_training-size-", local_com$TRAIN_SIZES,
                    "_n-variables-", length(all.vars(local_com$HMSC_XFORMULAS)),
                    "_k", k_fold
                )
            )
        )
    }
    
    return(unique(run_paths))
}

# A function to turn lists of lists into a dataframe with combinations as rows.
# ARGS:
#   - combinations: a list of lists containing elements 
#       TRAIN_SIZES, R_EFFECTS, STRATEGIES, HMSC_XFORMULAS, NEW_SAMPLE_SIZE.
build_param_grid <- function(combinations) {
    # create a grid, each line is a combination of variables
    # that way, we loop on this table instead of having nested loops
    grids <- lapply(combinations, function(p) {
        expand.grid(
            train_size = p$TRAIN_SIZES,
            r_effect = p$R_EFFECTS,
            strategy = p$STRATEGIES,
            formulas = p$HMSC_XFORMULAS,
            n_new_samples = p$NEW_SAMPLE_SIZES,
            stringsAsFactors = FALSE
        )
    })
    dplyr::distinct(dplyr::bind_rows(grids))
}

# A function to prepare dataset for Hmsc training, outputs a model template
# ARGS:
#   - subdataset: a dataframe. Must contain columns listed in x_cols and y_cols.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - y_cols: a list of strings. The columns containing species occurrences.
#   - formula: a formula for the Hmsc model. Based on names in x_cols/y_cols.
#   - random_effect: whether to add 
#       a random effect to Hmsc ("points", "carres" or "spatial") 
#       or not ("none", default).
prepare_hmsc_training <- function(
        subdataset, x_cols, y_cols, formula, random_effect = "none"){
    cli_alert_info("Preparing training set for HMSC...")

    # Format y and x data to a format that Hmsc accepts.
    ydata <- as.matrix(subdataset[y_cols]) 
    xdata <- as.data.frame(setNames(
        lapply(x_cols, function(col) subdataset[[col]]),
        x_cols)
    )
    studyDesign <- data.frame(
        points = as.factor(subdataset$point), 
        carres = as.factor(subdataset$carre), 
        spatial = as.factor(subdataset$id_point_annee)
    )

    # Depending on the random effect, the model has to be set up differently
    # FYI: here, we use occurence data. For that, a "probit" distribution is
    # the most logical option (but others are available such as "normal",
    # "poisson", "lognormal poisson")
    if (random_effect == "none") {
        # No effect: easy, do not add anything to the model.
        cli_alert_info("Creation of Hmsc object without random effects...")
        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata, 
            XFormula = formula, 
            distr = "probit",           
            studyDesign = studyDesign)
   
        
    } else if ((random_effect == "points") | (random_effect == "carres")) {
        # Random effect as "unit": points are grouped by squares (carres) or 
        # by single point observation. This assumes independant random effect
        # for each unit. Much faster to compute than rL.spatial but it means
        # that squares not sampled have no random effected fitted to them.
        cli_alert_info("Creation of Hmsc object with 'units' effect...")
        rL.units = HmscRandomLevel(
            units = unique(studyDesign[[random_effect]])) # random unit effect 
        rL.units = setPriors(
            rL.units, nfMin =1,  nfMax = 1) # limit number of latent variables
        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata, 
            ranLevels = setNames(list(rL.units), random_effect),
            XFormula = formula, 
            distr = "probit",
            studyDesign = studyDesign)

        
    }  else if (random_effect == "spatial") {
        # Random effect as "spatial": we consider the coordinates of each point
        # sampled. The random effect is a function of the distance between the
        # points. Takes more time to compute but (usually), yields to better 
        # results in prediction.
        cli_alert_info("Creation of Hmsc object with 'spatial' effect...")
        # convert coordinates to metric
        coords_sf <- st_as_sf(
            subdataset, 
            coords = c("LON", "LAT"), 
            crs = 4326)
        coords_proj <- st_transform(coords_sf, crs = 2154)

        # associate coordinate to the name of each sampled point
        xy <- st_coordinates(coords_proj)
        rownames(xy) <- as.character(subdataset$id_point_annee)
        colnames(xy) <- c("longitude_grid_2154", "latitude_grid_2154")

        # There should be one observation per point (because i_point_annee)
        # but we must make sure or else Hmsc will fail.
        if (any(duplicated(xy))) {
            stop("Duplicate coordinates found across distinct points; check LON/LAT data.")
            # xy <- xy[!duplicated(rownames(xy)), ] # quick fix for this error
        }
        
        # There are two possibilities : map the entire grid of points or use
        # nearest neighbour approximation. After ~1000 points its better to use
        # the approximation for faster computation.
        if (nrow(xy) < 1000) {
            rL.spatial <- HmscRandomLevel(sData = xy)
        } else {
            rL.spatial <- HmscRandomLevel(
                sData = xy, sMethod = "NNGP", nNeighbours = 10)
        }
        rL.spatial = setPriors(
            rL.spatial, nfMin =1,  nfMax = 1) # limit number of latent variables

        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata, 
            ranLevels = list("spatial" = rL.spatial),
            XFormula = formula, 
            distr = "probit",
            studyDesign = studyDesign)

    } else {
        stop(paste0("Random effect must be one of 'none', 'points',", 
        " 'carres', or 'spatial'. Got ", random_effect))
    }

    cli_alert_success("Created Hmsc object!\n\n")
    return(hmsc_object)
}

# A function to fit a Hmsc model and save results.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - save_to: a string. The path where the file will be saved (end with .rds).
#   - nchains: a numeric. The number of chains to run.
#   - thin: a numeric. The number of steps between each recording of a sample.
#   - nsamples: a numeric. The number of samples to collect
#   - ntransient: a numeric. The number of steps to wait for before 
#       collecting samples.
#   - freq_verbose: a numeric (default is 100). The frequency of verbose 
#       messages, you get one every "freq_verbose" step.
#   - allow_parallel: a boolean (default is TRUE). 
#       If nchains >1, allows to compute chains in parallel for faster fitting. 
#       Removes verbose messages during fitting.
fitting_hmsc <- function(
    hM, 
    nchains,
    thin,
    nsamples,
    ntransient,
    freq_verbose = 100,
    allow_parallel = TRUE,
    save_to = NULL) {
    # Start clock to get runtime
    start <- Sys.time()
    cli_alert_info(paste0("Started fitting at: ", start))

    # When HMSC runs chains in parallel, it does not show progress bar
    # Since I was expecting some update messages during the runs at first, 
    # its probably better to print some warning messages.
    if (!is.null(freq_verbose) & (nchains > 1) & allow_parallel) {
        cli_alert_warning(paste0("Cannot display fitting progress when running",
        " chains in parallel."))
        cli_alert_warning("Set `allow_parallel` to `FALSE` to see progress.")
    }

    # Depending on allow_parallel = TRUE/FALSE, overwrite nParallel
    fitted.hmsc <-  sampleMcmc(
        hM, 
        thin = thin, 
        samples = nsamples, 
        transient = ntransient, 
        nChains = nchains, 
        nParallel = if (allow_parallel) nchains else 1, 
        updater=list(GammaEta=FALSE),
        verbose = freq_verbose)
    
    # Print runtime
    stop <- Sys.time()
    cli_alert_info(paste0("Completed fitting at: ", stop))
    cli_alert_info(paste0("Time elapsed: ", round(stop-start, 2)))

    # Save run 
    if (!is.null(save_to)) {
        cli_alert_info("Saving model...")
        saveRDS(fitted.hmsc, file = save_to)
        cli_alert_success("Model saved!\n\n")    
    }
    return(fitted.hmsc)
}

# A function to display convergence diagnostics for a Hmsc model.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - nchains: a numeric. The number of chains to run.
#   - thin: a numeric. The number of steps between each recording of a sample.
#   - save_folder: a string. The path where plotted PDFs will be saved.
convergence_hmsc <- function(hM, nchains, thin, save_folder) {
    # Hmsc object can be (and should be) converted to the commonly used Coda
    # format in Bayesian statistics.
    coda_outputs <- convertToCodaObject(hM)
    
    # Summary plots (not necessary here but might be useful one day):
    # MCMCsummary(object = coda_outputs$Beta, round = 2) 
    # MCMCplot(object = coda_outputs$Beta)

    # HMSC has many parameters, in our case we are only interested in:
    #   - Beta: fixed effects
    #   - Omega: random variation in co-occurence
    for (param in c("Beta", "Omega")) {
        # skip if name not in parameters (e.g. no omega if no random effect)
        if (!(param %in% names(coda_outputs))) {
            next
        } 
        cli_alert_info(paste0("*** [Parameter: ", param, "] ***"))

        ### Convergence diagnostics
        # Convert Omega output to match Beta format
        if (param == "Omega") {
            chains <- coda_outputs$Omega[[1]]  # 3D array: [iter, sp, sp]
        } else {
            chains <- coda_outputs[[param]]
        }
        
        ## Traceplot (Rhat and effective size)      
        cli_alert_info("Computation of traceplots, can take some time...")
        # Fast version:
        tryCatch({ # avoids debugger mode
            MCMCtrace(
                object = chains,
                pdf = TRUE,
                filename = file.path(
                    save_folder, 
                    paste0(param, "_all_traceplots.pdf")),
                ind = TRUE,
                open_pdf = FALSE,
                plot = TRUE,
                Rhat = TRUE, 
                n.eff = TRUE, 
                type = "both" # explicitly request trace + density
            )}, 
        error = function(e) {
            message("Error in MCMCtrace: ", e$message)
        })
        # # Slower but "more beautiful" version
        # traceplots <- ggplot_custom_MCMCtrace(
        #     coda_object = coda_fitted_model$Beta,
        #     show_Rhat = TRUE,
        #     show_Neff = TRUE)       
        # for (i_plot in seq_along(traceplots)) {
        #     standardised_ggplot_save(
        #         figure = (traceplots[[i_plot]]$trace + traceplots[[i_plot]]$density), 
        #         save_path = file.path(
        #             save_folder, 
        #             paste0("Beta_", i_plot, "_traceplot.pdf")))
        # }
        cli_alert_info("Traceplots saved!")


        ## Effective size
        cli_alert_info("-> Effective size:")
        eff_size <- effectiveSize(chains)
        eff_size_plot <- ggplot_bars(
                as.data.frame(eff_size), "eff_size", 
                breaks = ceiling(c(0, seq(100, max(eff_size), length.out = 19))),
                underlayers = list(
                    bad = annotate("rect", 
                        xmin = -Inf, xmax = 100, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[1], alpha = 0.5),
                    acceptable = annotate("rect", 
                        xmin = 100, xmax = 400, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[5], alpha = 0.5),
                    good = annotate("rect", 
                        xmin = 400, xmax = Inf, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[2], alpha = 0.5),
                    captions = labs(caption = "Green: good, Orange: acceptable, Red: bad."))
                )
        print(eff_size_plot)
        interpret_diagnostics(
            eff_size, bad = 100, good = 400, order = "high_better", mode = "quick")
        standardised_ggplot_save(
            figure = eff_size_plot, 
            save_path = file.path(
                save_folder, 
                paste0("hist_neff_", param, ".pdf")))

        ## Gelman-Rubin convergence diagnostic
        if (nchains > 1) {
            cli_alert_info("-> Gelman-Rubin convergence diagnostic:")
            psrf <- gelman.diag(
                chains,  multivariate = FALSE)$psrf[, "Point est."]
            psrf_plot <- ggplot_bars(
                as.data.frame(psrf), "psrf", bins = 10,                
                underlayers = list(
                    good = annotate("rect", 
                        xmin = -Inf, xmax = 1.05, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[2], alpha = 0.5),
                    acceptable = annotate("rect", 
                        xmin = 1.05, xmax = 1.1, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[5], alpha = 0.5),
                    good = annotate("rect", 
                        xmin = 1.1, xmax = Inf, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[1], alpha = 0.5),
                    captions = labs(caption = "Green: good, Orange: acceptable, Red: bad."))
                )
            print(psrf_plot)
            interpret_diagnostics(psrf, bad=1.1, good=1.05, mode = "quick")
            standardised_ggplot_save(
                figure = psrf_plot, 
                save_path = file.path(
                    save_folder, 
                    paste0("hist_psrf_", param, ".pdf")))

            # ## Geweke diagnostic
            # cli_alert_info("-> Geweke convergence diagnostic:")
            # # Rule of thumb: <2 (no proof of non-convergence), >2 (monitor convergence closely)
            # geweke_estimates_all <- geweke.diag(chains)
            # for (i in seq(nchain(chains))) {
            #     cli_alert_info(paste0("[Chain ", i, "]"))
            #     interpret_diagnostics(geweke_estimates_all[[i]][[1]], good=2, mode = "quick")
            # }

            ## Autocorrelation
            cli_alert_info("-> Autocorrelation:")
            for (lag in c(20, 50)) {
                autocorr_estimates <- autocorr.diag(
                    chains, lags=c(lag/thin))
                cli_alert_info(paste0("[Lag ", lag, "]:"))
                interpret_diagnostics(autocorr_estimates, good=0.1, mode = "quick")
            }
            
        }
    }
    cat("\n")
}

# A function to display convergence diagnostics for a Hmsc model.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - subset: a dataframe. Must contain columns listed in x_cols and y_cols.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - sp_cols: a list of strings. The columns containing species occurrences.
evaluate_hmsc_performances <- function(hM, subset, x_cols, sp_cols) {
    local_preds_list <- predict_hmsc(
        hM = hM, 
        df = subset, 
        x_variables = x_cols)
    local_preds <- abind(local_preds_list, along = 3) # model predictions
    local_Y <- as.matrix(subset[sp_cols])             # actual observations
    
    # Extract metric by comparing predictions and observed values
    evaluateModelFitCustom(hM = hM, Y = local_Y, predY = local_preds)
}

# A function to display XX and XY associations after Hmsc fitting.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - x_groups_cats: a list of integers. For each x variable, 
#       a number assigning it to a group (for variance partitioning)
#   - x_groups_names: a list of strings. 
#       A label for each number in x_groups_cat.
#   - save_folder: a string. The path where plotted PDFs will be saved.
#   - supportLevel: a numeric between 0 and 1. 
#       The minimum confidence to display results (default is 0.95)
analyses_hmsc <- function(
        hM, save_folder, x_groups_cats, x_groups_names, supportLevel = 0.95) {
    # X-Y associations
    for (param in c("Beta", "Omega")) {
        if (is.null(hM$ranLevels) & (param=="Omega")){
            next
        }
        post_association = getPostEstimate(hM, parName = param)
        XY_grid <- ggplot_custom_plotBeta(
            hM, post = post_association, supportLevel = supportLevel)
        standardised_ggplot_save(
            figure = XY_grid, 
            save_path = file.path(save_folder, paste0(param, "_XY_associations.pdf")))
    }

    if (!is.null(hM$ranLevels)) {
        rand_XX_grid <- ggplot_custom_random_corr_associations(
            hM, supportLevel = supportLevel)
        standardised_ggplot_save(
            figure = rand_XX_grid, 
            save_path = file.path(save_folder, "random_XX_associations.pdf"))
    }
    
    # Variance partitionning 
    vp = computeVariancePartitioning(
        hM, 
        group = x_groups_cats, # c(1,2,2)
        groupnames = x_groups_names) # c("habitat","climate"))
    variance_bars <- ggplot_custom_plotVariancePartitioning(hM, VP = vp)
    standardised_ggplot_save(
        figure = variance_bars, 
        save_path = file.path(save_folder, "variance_partitioning.pdf"))
}

# A function to compute the habitat suitability map for a species 
# based on Hmsc predictions.
# ARGS:
#   - df: a dataframe. Must contain columns listed in x_cols and y_cols.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - parent_folder: a string. The path containing the subfolder for each model.
#   - loop_prefix: a string. The prefix of the subfolder name.
#   - loop_element: a string. Middle element for subfolder name.
#   - loop_suffix: a string. The suffix of the subfolder name.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#       Goes after suffix.
#   - xlabel: a string. The label for the x-axis of the plot 
#       (default is "Model type").
#   - sp: a string. The name of a species (in df) to plot.
LEGACY_map_results_hmsc <- function(
        df,
        x_cols,
        parent_folder, 
        loop_prefix, 
        loop_element, 
        loop_suffix, 
        k_fold, 
        sp) {

    local_folder <- file.path(
            parent_folder,
            paste0(loop_prefix, loop_element, loop_suffix, k_fold))
    local_model <- readRDS(file.path(local_folder, "train_outputs.rds")) 

    cli_alert_info("Making predictions, can take a few minutes...")
    full_preds_list <- predict_hmsc(
        hM = local_model, df = df, x_variables = x_cols)
    cli_alert_info("Predictions are ready.")

    # Predictions are DISTRIBUTIONS: let's take the average for each point
    cli_alert_info("Converting distributions to point average...")
    full_preds_avg <- apply(
        simplify2array(full_preds_list), c(1,2), mean)
    full_preds_sd <- apply(
        simplify2array(full_preds_list), c(1,2), sd)
    df_preds <- df |>
        mutate(suitability = full_preds_avg[, sp]) |>
        mutate(uncertainty = full_preds_sd[, sp]) |>
        mutate(subset = ifelse(
            id_point_annee %in% k_fold_points$training_points[[k_fold]], 
            "train",
            ifelse(
                id_point_annee %in% k_fold_points$val_training_points[[k_fold]], 
                "val", 
                "test")))
    # To cite 10.5281/ZENODO.11067678: "The models cannot provide a direct 
    # indication of the probability of presence – they instead give an 
    # index of habitat suitability, measured on a scale from 0 to 1."
    cli_alert_info(paste0("\tSelected species: ", sp))

    # # plot training samples
    # p0.1 <- ggplot_categorical_df_on_background_map(
    #     background_map = ggplot_get_france_base_map("national"), 
    #     df = df_preds |> filter(subset == "train"), 
    #     LON = "LON",
    #     LAT = "LAT",
    #     column = NULL,
    #     legend_title = "Training location")
    # print(p0.1)
    # standardised_ggplot_save(
    #     figure = p0.1, 
    #     save_path = file.path(local_folder, "training_samples.pdf"))
    # p0.2 <- ggplot_categorical_df_on_background_map(
    #     background_map = ggplot_get_france_base_map("national"), 
    #     df = df_preds |> filter(subset == "test"), 
    #     LON = "LON",
    #     LAT = "LAT",
    #     column = NULL,
    #     legend_title = "Training location")
    # print(p0.2)
    # standardised_ggplot_save(
    #     figure = p0.2, 
    #     save_path = file.path(local_folder, "test_samples.pdf"))

    # show map of suitability per site
    p1.1 <- ggplot_quantitative_df_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            df = df_preds, 
            LON = "LON",
            LAT = "LAT",
            column = "suitability",
            unit = paste0("Estimated HSI for ", sp)) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p1.1)
    standardised_ggplot_save(
        figure = p1.1, 
        save_path = file.path(local_folder, "suitability_for_species.pdf"))
    shp1.2 <- interpolate_scattered_points_to_hexagons(
            df = df_preds,  
            column = "suitability",
            res_km = RES_KM, 
            LON = "LON",
            LAT = "LAT",
            idp = 2,           
            maxdist_m = 100)
    p1.2 <- ggplot_quantitative_shapefile_on_background_map(
        background_map = ggplot_get_france_base_map("national"),
        shapefile = shp1.2,
        layer_name = "interpolated_value",
        unit = paste0("Estimated HSI for ", sp),
        limits = NULL,
        precision_auto_limits = 1) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p1.2)
    standardised_ggplot_save(
        figure = p1.2, 
        save_path = file.path(local_folder, "suitability_for_species_hexagons.pdf"))    
    
    

    # show map of uncertainty per site (function of suitability)
    p2.1 <- ggplot_quantitative_df_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            df = df_preds, 
            LON = "LON",
            LAT = "LAT",
            column = "uncertainty",
            unit = paste0("Standard deviation of HSI for ", sp)) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p2.1)
    standardised_ggplot_save(
        figure = p2.1, 
        save_path = file.path(local_folder, "certainty_of_suitability.pdf"))
    shp2.2 <- interpolate_scattered_points_to_hexagons(
            df = df_preds,  
            column = "uncertainty",
            res_km = RES_KM, 
            LON = "LON",
            LAT = "LAT",
            idp = 2,           
            maxdist_m = 100)
    p2.2 <- ggplot_quantitative_shapefile_on_background_map(
        background_map = ggplot_get_france_base_map("national"),
        shapefile = shp2.2,
        layer_name = "interpolated_value",
        unit = paste0("Standard deviation of HSI for ", sp),
        limits = NULL,
        precision_auto_limits = 1e-5) +
        labs(caption = "HSI = 'Habitat Suitability Index'")
    print(p2.2)
    standardised_ggplot_save(
        figure = p2.2, 
        save_path = file.path(local_folder, "certainty_of_suitability_hexagons.pdf"))    

    return(df_preds)
}

# A function to automatise verbose interpretation of diagnostic vectors.
# ARGS:
#   - vector: a vector of values to analyse.
#   - bad: a numeric. The threshold under which values reveal a bad fit.
#   - good: a numeric. The threshold over which values reveal a good fit.
#   - order: a string. Indicates if lower is better ("low_better", default) 
#       or higer is bette (high_better).
#   - mode: a sting. Indicates if showing full analysis ("full", default) 
#       or a quick 1-line summary ("quick").
interpret_diagnostics <- function(
    vector, 
    good, 
    bad = NULL, 
    mode = "full",
    order = "low_better") {
    
    if (order == "low_better") {
        n_good <- sum(vector < good)    

        if (is.null(bad)) {
            n_bad <- sum(vector > good)
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: <", good," (good)"))
            }

            n_acceptable <- 0
        } else {
            n_acceptable <- sum((vector > good) & (vector < bad))
            n_bad <- sum(vector > bad)    
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: <", good," (good), ",
                    good, "-", bad, ", (acceptable), >",
                    bad, " (bad)"))
            }
        }
        
    } else if (order == "high_better") {
        n_good <- sum(vector > good)    

        if (is.null(bad)) {
            n_bad <- sum(vector < good)  
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: >", good," (good)"))
            }

            n_acceptable <- 0
        } else {
            n_acceptable <- sum((vector < good) & (vector > bad))
            n_bad <- sum(vector < bad)
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: >", good," (good), ",
                    good, "-", bad, ", (acceptable), <",
                    bad, " (bad)"))   
            }
        }
     
    }
    else {
        stop(paste0("Mode should be one of 'low_better' or 'high_better'. ",
        "Got '", mode, "' ."))
    }
    
    if (mode == "full") {
        if (n_good > 0) {
            cli_alert_success(paste0(
            "- Number of 'good' estimates: ", n_good, 
            " (", round(100*n_good/length(vector), 2), "% of given values)"
            ))
        }
        if ((n_acceptable > 0) & !(is.null(bad))) {
            cli_alert_info(paste0(
            "- Number of 'acceptable' estimates: ", n_acceptable, 
            " (", round(100*n_acceptable/length(vector), 2), "% of given values)"
            ))  
        }
        if (n_bad > 0) {
            cli_alert_info(paste0(
            "- Number of 'bad' estimates: ", n_bad, 
            " (", round(100*n_bad/length(vector), 2), "% of given values)"
            ))
        }
    } else if (mode == "quick") {
        if (is.null(bad)) {
            cli_alert_info(paste0(
                n_good, " (", round(100*n_good/length(vector), 2), "%) are good ",
                "and ", n_bad ," (", round(100*n_bad/length(vector), 2),"%) are bad."))
        } else {
            cli_alert_info(paste0(
                n_good, " (", round(100*n_good/length(vector), 2), "%) are good, ",
                n_acceptable, " (", round(100*n_acceptable/length(vector), 2), "%) are acceptable ",
                "and ", n_bad ," (", round(100*n_bad/length(vector), 2),"%) are bad."))
        }
        
    }

}

# A function that automates the process of making predictions with Hmsc models.
# ARGS:
#   - hM: a Hmsc fitted model object.
#   - df: a dataframe with columns "point", "id_point_annee" 
#       and names in x_variables.
#   - x_variables: a list of strings. The names of columns to keep in data.
predict_hmsc <- function(hM, df, x_variables) {
    # Format explanatory variables to Hmsc expected format
    XData <- as.data.frame(setNames(
        lapply(x_variables, function(col) df[[col]]),
        x_variables))

    if ("spatial" %in% names(hM$ranLevels)) {
        # Format coordinates associated to each point
        coords_sf <- sf::st_as_sf(df, coords = c("LON", "LAT"), crs = 4326)
        coords_proj <- sf::st_transform(coords_sf, crs = 2154)
        xy_new <- sf::st_coordinates(coords_proj)
        rownames(xy_new) <- as.character(df$id_point_annee)
        colnames(xy_new) <- c("longitude_grid_2154", "latitude_grid_2154")

        # Use Gradient (instead of studyDesign), for spatial gradients
        Gradient <- prepareGradient(
            hM,
            XDataNew = XData,
            sDataNew = list(spatial = xy_new)
        )

        # Make prediction on new dataset
        return(predict(hM, Gradient = Gradient, expected = TRUE))

    } else {
        # Format study design (grouping of samples together)
        studyDesign <- data.frame(
            points = as.factor(df$point),
            carres = as.factor(df$carre), 
            spatial = as.factor(df$id_point_annee)
        )

        # Make prediction on new dataset
        return(predict(
            hM, XData = XData, studyDesign = studyDesign, expected = TRUE))
    }
}

# A function that mimciks Hmsc::evaluateModelFit but can also work on 
# non-training data.
# ARGS :
#   - hM: a Hmsc fitted model object.
#   - y: a matrix of species observation ("ground truth").
#   - predY: the predictions made by the model.
evaluateModelFitCustom <- function(hM, Y, predY) {

    ns <- ncol(Y) # number of samples per observation/species
    mPredY <- apply(predY, c(1, 2), mean) # mean prediction per obs/species

    # Initialise metrics to compute
    RMSE <- rep(NA, ns)     # RMSE (the lower the better)
    AUC <- rep(NA, ns)      # AUC (the closer to 1, the better)
    TjurR2 <- rep(NA, ns)   # Tjur R² (% of variance explained)

    # For each sample
    for (j in seq_len(ns)) {
        sel <- !is.na(Y[, j])   # extract observations/species
        obs <- Y[sel, j]        # get observed value
        pred <- mPredY[sel, j]  # get predicted value

        # compute RMSE
        RMSE[j] <- sqrt(mean((obs - pred)^2))

        # compute AUC (only meaningful if both 0s and 1s present)
        if (length(unique(obs)) == 2) {
            AUC[j] <- as.numeric(pROC::auc(obs, pred, quiet = TRUE))
        }

        # compute Tjur R2: difference in mean predicted probability between
        # presences and absences
        if (length(unique(obs)) == 2) {
            TjurR2[j] <- mean(pred[obs == 1]) - mean(pred[obs == 0])
        }
    }

    names(RMSE) <- names(AUC) <- names(TjurR2) <- colnames(Y)
    return(list(RMSE = RMSE, AUC = AUC, TjurR2 = TjurR2))
}

# A function that computes the uncertainty of a Hmsc model on its predictions.
# ARGS:
#   - hM: a Hmsc fitted model object.
#   - df: a dataframe with columns "point", "id_point_annee" 
#       and names in x_variables.
#   - x_variables: a list of strings. The names of columns to keep in data.
get_uncertainty_hmsc <- function(hM, df, x_cols) {
    # Get predictions
    predicted_occurrences <- predict_hmsc(
        hM = hM, df = df, x_variables = x_cols)
    
    # predicted_occurrences contains a list of length = number of samples.
    # for each sample, we get a matrix of n_obs x n_species
    # HERE, WE DEFINE UNCERTAINTY AS
    # the average of the standard deviation accross observed species
    sd_point_sp <- apply(simplify2array(predicted_occurrences), c(1,2), sd)
    uncertainty_per_point <- as_tibble(sd_point_sp) |> 
        mutate(average_sd = rowMeans(across(everything())))
    
    return(uncertainty_per_point$average_sd)
}

# ====================== HELPERS for scores functions =========================
# A function to abbreviate hyphenated loop_element labels for the x-axis: 
# a single word is left as-is, a multi-word (e.g. "foo-bar-baz") element 
# becomes initials (e.g. "FBB").
# ARGS:
#   - x: a vector of string.
abbreviate_loop_labels <- function(x) {
    sapply(strsplit(x, "-"), function(words) {
        if (length(words) == 1) {
            words
        } else {
            paste0(toupper(substr(words, 1, 1)), collapse = "")
        }
    })
}

# A function that prints the plot, optionally saves it, and emits cli messages.
# ARGS:
#   - p: a ggplot object.
#   - save_to: a string. A path which should end with .pdf. 
#       If NULL, does not save pdf.
#   - what: a string. Type of success in message.
finalize_plot <- function(p, what = "performances", save_to = NULL) {
    print(p)
    if (!is.null(save_to)) {
        standardised_ggplot_save(figure = p, save_path = save_to)
        cli_alert_success(paste0("Plot of ", what, " saved!"))
    }
    cli_alert_success(paste0("Plot of ", what, " ready!\n\n"))
}

# A function that reads a single metric's column from one model-run folder's
# `{subset_name}_scores.csv`. Centralises the MSE/RMSE/AUC/TjurR2 handling.
# ARGS:
#   - run_path: a string. Path to parent folder of `{subset_name}_scores.csv`.
#   - subset_name: a string. Name of the subset in `{subset_name}_scores.csv`.
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2). 
#       Default is NULL (output all dataframe). 
#       If a metric is given, outputs only the column of that metric.
load_metric_scores <- function(run_path, subset_name, metric) {
    # SECURITY : auto-fetches name of the strategy and number of new samples
    # if strategy and new sample mismatch, stop functions
    # I have quick fixes ready in comments but these checks are handled sooner
    # now, and should not happen here in the new versions of the script.
    pattern_strategy <- "strategy-([^_]+)_new"
    pattern_new_samples <- "new-samples-([^_]+)_training"
    extracted_strategy <- regmatches(
        run_path, regexec(pattern_strategy, run_path))[[1]][2]
    extracted_new_samples <- regmatches(
        run_path, regexec(pattern_new_samples, run_path))[[1]][2]
    if (extracted_strategy == "none") {
        if (extracted_new_samples != "0") {
            stop(paste0(
                "Found unconsistent number of new samples '", 
                extracted_new_samples, "' for strategy '", 
                extracted_strategy , "'."))
            # # Fix (LEGACY)
            # cli_alert_warning(paste0(
            #     "Found unconsistent number of new samples '", 
            #     extracted_new_samples, "' for strategy '", 
            #     extracted_strategy , "'. Replacing by '0'."))
            # run_path <- sub(
            #     pattern_new_samples, "new-samples-0_training", 
            #     run_path)
        }
    }
    if (extracted_new_samples == "0") {
        if (extracted_strategy != "none") {
            stop(paste0(       
                "Found unconsistent strategy name '", 
                extracted_strategy, "' for number of new samples '", 
                extracted_new_samples , "'."))
            # # Fix (LEGACY)
            # cli_alert_warning(paste0(
            #     "Found unconsistent strategy name '", 
            #     extracted_strategy, "' for number of new samples '", 
            #     extracted_new_samples , "'. Replacing by 'none'."))
            # run_path <- sub(
            #     pattern_strategy, "strategy-none_new", 
            #     run_path)
        }
    }    

    # Fetch file (suppress import messages)
    file_path <- file.path(run_path, paste0(subset_name, "_scores.csv"))
    df <- read_csv(file_path, show_col_types = FALSE)

    # add MSE (RMSE is computed b default, MSE is RMSE squared)
    df$MSE <- df$RMSE^2

    # output either full dataframe or single column vector
    if (is.null(metric)) {
        df
    } else if (metric %in% c("MSE", "RMSE", "AUC", "TjurR2")) {
        df[[metric]]
    } else {
        stop(paste0("Metric '", metric, "' is not handled by this function."))
    }
}

# A function that loads reference + "other" (looped) scores for a single metric across all
# k_folds, subsets and loop_elements.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - reference_model_combination: a list of single parameters. 
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - loop_on: a string. The name of a parameter in a combination. 
#       The parameter which has several values.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - subset_names: a list string. Usually c("train", "val", "test").
load_reference_and_other_scores <- function(
        parent_folder, reference_model_combination, loop_model_combination,
        loop_on, k_fold, metric, subset_names) {
    cli_alert_info("Loading base results...")
    ref_scores <- list()
    for (subset_name in subset_names) {
        ref_scores[[subset_name]] <- lapply(seq(k_fold), function(k) {
            run_path <- make_run_path(
                parent_folder, reference_model_combination, k)
            load_metric_scores(run_path, subset_name, metric)
        })
    }

    cli_alert_info("Loading other results...")
    other_scores <- list()
    for (subset_name in subset_names) {
        other_scores[[subset_name]] <- list()
        run_paths <- make_run_path(
            parent_folder, loop_model_combination, "")
        loop_elements <- loop_model_combination[[loop_on]]

        for (i in 1:length(run_paths)) {
            loop_element <- loop_elements[i]
            if (is.list(loop_element)) {
                if (is_formula(loop_element[[1]])) {
                    loop_element <- length(all.vars(loop_element[[1]]))
                } else {
                    stop("Function was not made to handle loop_element as list when not a formula.")
                }
            }
            
            other_scores[[subset_name]][[loop_element]] <- lapply(
                seq(k_fold), function(k) {
                    local_run_path <- paste0(run_paths[i], k)
                    load_metric_scores(local_run_path, subset_name, metric)
                }
            )
        }
    }

    list(ref_scores = ref_scores, other_scores = other_scores)
}

# A function computes per-loop_element / per-subset (/ per-species) 
# differences and stats between "other" and "reference" scores.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - reference_model_combination: a list of single parameters. 
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements). 
#   - loop_on: a string. The name of a parameter in a combination. 
#       The parameter which has several values.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - group_species: whether to take the mean 
#       accross all k_folds and species (TRUE, default) or 
#       only accross k_folds (FALSE).
#   - species_names: names of species in CSV 
#       (rownames are not available from csvs).
compute_score_diffs <- function(
        parent_folder, reference_model_combination, loop_model_combination, 
        loop_on, k_fold, metric = "MSE", 
        subset_names = c("train", "val", "test"),
        group_species = TRUE, species_names = NULL) {

    # Auto load of reference and loop results
    loaded <- load_reference_and_other_scores(
        parent_folder, reference_model_combination, loop_model_combination, 
        loop_on, k_fold, metric, subset_names)
    ref_scores <- loaded$ref_scores
    other_scores <- loaded$other_scores

    # Number of species = length of a score vector: build names.
    if (!group_species) {
        n_species <- length(ref_scores[[subset_names[1]]][[1]])
        if (is.null(species_names)) {
            species_names <- paste0("species_", seq_len(n_species))
        } else if (length(species_names) != n_species) {
            stop(paste0(
                "`species_names` has length ", length(species_names),
                " but the score files have ", n_species, " rows (species)."
            ))
        }
    }

    # For each element to loop on, store its score compared to reference
    loop_elements <- loop_model_combination[[loop_on]]
    scores_df <- tibble()
    raw_diffs_df <- tibble()

    for (subset_name in subset_names) {
        for (loop_element in loop_elements) {

            if (is_formula(loop_element)) {
                loop_element <- length(all.vars(loop_element))
            }

            # (k_fold x species) matrix of raw differences before any
            # averaging, so both fold- and species-level variability remain
            # available for the CI/SD computation below.
            diff_matrix <- do.call(rbind, lapply(seq(k_fold), function(k) {
                other_scores[[subset_name]][[loop_element]][[k]] -
                    ref_scores[[subset_name]][[k]]
            }))
            # rows = folds, columns = species

            if (group_species) {
                # Average across species within each fold first, so each
                # fold contributes exactly one observation; the interval
                # then reflects fold-to-fold (across k_fold) variability.
                diff_obs_list <- list(
                    all = rowMeans(diff_matrix, na.rm = TRUE))
            } else {
                # Keep species separate: for each species, the observations
                # are that species' differences across the k folds.
                n_species <- ncol(diff_matrix)
                diff_obs_list <- setNames(
                    lapply(seq_len(n_species), function(s) diff_matrix[, s]),
                    species_names
                )
            }

            for (obs_name in names(diff_obs_list)) {
                diff_obs <- diff_obs_list[[obs_name]]
                species_val <- if (group_species) NA_character_ else obs_name

                # raw observations (used by the boxplot)
                raw_diffs_df <- bind_rows(
                    raw_diffs_df,
                    tibble(
                        diff_value = diff_obs,
                        subset = subset_name,
                        loop_element = loop_element,
                        species = species_val))

                # summary stats (used by the dot-whisker plot, and returned
                # alongside raw_diffs for the boxplot)
                n <- sum(!is.na(diff_obs))
                avg_metric <- mean(diff_obs, na.rm = TRUE)
                sd_metric <- sd(diff_obs, na.rm = TRUE)

                # t critical value is more appropriate than a normal
                # approximation when n (e.g. k_fold) is small; NA when
                # there are fewer than 2 observations.
                crit_value <- if (n > 1) qt(0.975, df = n - 1) else NA_real_
                conf_margin <- crit_value * sd_metric / sqrt(n)

                scores_df <- bind_rows(
                    scores_df,
                    tibble(
                        average_metric = avg_metric,
                        sd_lower = avg_metric - sd_metric,
                        sd_upper = avg_metric + sd_metric,
                        ci_lower = avg_metric - conf_margin,
                        ci_upper = avg_metric + conf_margin,
                        n_obs = n,
                        conf_margin = conf_margin,
                        sd_metric = sd_metric,
                        subset = subset_name,
                        loop_element = loop_element,
                        species = species_val))
            }
        }
    }

    # if loop_elements are formulas, apply transformation to get numerics
    if (all(sapply(loop_elements, is_formula))) {
        loop_elements <- sort(sapply(sapply(loop_elements, all.vars), length))
    }

    # Ensure right formatting (mainly that factors are notnumerics)
    scores_df <- scores_df |>
        mutate(loop_element = factor(loop_element, levels = loop_elements)) |>
        mutate(subset = factor(subset, levels = subset_names))
    raw_diffs_df <- raw_diffs_df |>
        mutate(loop_element = factor(loop_element, levels = loop_elements)) |>
        mutate(subset = factor(subset, levels = subset_names))

    # Add column for species if necessary
    if (!group_species) {
        scores_df <- scores_df |> 
            mutate(species = factor(species, levels = species_names))
        raw_diffs_df <- raw_diffs_df |> 
            mutate(species = factor(species, levels = species_names))
    }

    list(summary = scores_df, raw_diffs = raw_diffs_df, ref = ref_scores)
}

# A function that computes the number of species that reached a certain
# level of improvement compared to a reference, in proportion.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - reference_model_combination: a list of single parameters. 
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - loop_on: a string. The name of a parameter in a combination. 
#       The parameter which has several values.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - species_names: names of species in CSV 
#       (rownames are not available from csvs).
#   - proportion: a numeric between 0 and 1. The threshold to count a species
#       score as improved. 
compute_number_of_improvements <- function(
        parent_folder, reference_model_combination, loop_model_combination, 
        loop_on, k_fold, metric = "MSE", 
        subset_names = c("train", "val", "test"),
        species_names = NULL, proportion = 1/100) {

    # Auto load of reference and loop results
    loaded <- load_reference_and_other_scores(
        parent_folder, reference_model_combination, loop_model_combination, 
        loop_on, k_fold, metric, subset_names)
    ref_scores <- loaded$ref_scores
    other_scores <- loaded$other_scores

    # Number of species = length of a score vector: build names.
    n_species <- length(ref_scores[[subset_names[1]]][[1]])
    if (is.null(species_names)) {
        species_names <- paste0("species_", seq_len(n_species))
    } else if (length(species_names) != n_species) {
        stop(paste0(
            "`species_names` has length ", length(species_names),
            " but the score files have ", n_species, " rows (species)."
        ))
    }

    # For each element to loop on, store its score compared to reference
    loop_elements <- loop_model_combination[[loop_on]]
    prop_diffs_df <- tibble()

    for (subset_name in subset_names) {
        for (loop_element in loop_elements) {

            if (is_formula(loop_element)) {
                loop_element <- length(all.vars(loop_element))
            }

            # (k_fold x species) matrix of raw differences before any
            # averaging, so both fold- and species-level variability remain
            # available for the CI/SD computation below.
            diff_matrix <- do.call(rbind, lapply(seq(k_fold), function(k) {
                (other_scores[[subset_name]][[loop_element]][[k]] -
                    ref_scores[[subset_name]][[k]]) /
                    ref_scores[[subset_name]][[k]]
            }))      
            
            # replace eventual NA by 0
            diff_matrix[is.na(diff_matrix)] <- 0

            # Keep species separate: for each species, the observations
            # are that species' differences across the k folds.
            n_species <- ncol(diff_matrix)
            diff_obs_list <- setNames(
                lapply(seq_len(n_species), function(s) diff_matrix[, s]),
                species_names
            )

            for (k in seq(k_fold)) {
                if (grepl("MSE", metric)) {
                    prop_diffs_df <- bind_rows(
                        prop_diffs_df,
                        tibble(
                            improvements = sum(-diff_matrix[k, ] > proportion),
                            subset = subset_name,
                            loop_element = loop_element,
                            k_fold = k))
                } else {
                    prop_diffs_df <- bind_rows(
                        prop_diffs_df,
                        tibble(
                            improvements = sum(diff_matrix[k, ] > proportion),
                            subset = subset_name,
                            loop_element = loop_element,
                            k_fold = k))
                }
            }
        }
    }

    # if loop_elements are formulas, apply transformation to get numerics
    if (all(sapply(loop_elements, is_formula))) {
        loop_elements <- sort(sapply(sapply(loop_elements, all.vars), length))
    }

    # Ensure right formatting (mainly that factors are notnumerics)
    prop_diffs_df <- prop_diffs_df |>
        mutate(loop_element = factor(loop_element, levels = loop_elements)) |>
        mutate(subset = factor(subset, levels = subset_names))

    list(prop_diffs = prop_diffs_df, n_species = n_species)
}

# =========================== Scores functions ================================
# A function to compare scores between k_folds, subset and model type.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - loop_on: a string. The name of a parameter in a combination. 
#       The parameter which has several values.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#       Goes after suffix.
#   - xlabel: a string. The xlabel for the plot (default is "Model type").
#   - ylabel: a string. The ylabel for the plot (default is "Average Score").
#   - group_species: whether to take the mean 
#       accross all k_folds and species (TRUE, default) or 
#       only accross k_folds (FALSE).
#   - species_names: names of species in CSV 
#       (rownames are not available from csvs).
#   - barplot: makes a barplot (TRUE, default) or a pointplot (FALSE).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - save_to: a string. Path which should end with .pdf. 
#       If NULL, does not save pdf.
barplot_raw_scores <- function(
        parent_folder,
        loop_model_combination,
        loop_on,
        k_fold,
        xlabel = "Model type",
        ylabel = "Average Score",
        group_species = TRUE,
        species_names = NULL,
        barplot = TRUE,
        subset_names = c("train", "val", "test"),
        save_to = NULL) {

    cli_alert_info("Fetching scores...")

    # No auto import of results in this function (mainly because I didnt need
    # it at the time...)
    # The following loop import each necessary scores.csv file and parse it 
    # into a dataframe.
    run_paths <- make_run_path(
        parent_folder, loop_model_combination, "")
    loop_elements <- loop_model_combination[[loop_on]]
    scores_df <- data.frame()
    for (k in seq(k_fold)) {
        for (i in 1:length(run_paths)) {
            run_path <- paste0(run_paths[i], k)
            loop_element <- loop_elements[i]

            if (is.list(loop_element)) {
                if (is_formula(loop_element[[1]])) {
                    loop_element <- length(all.vars(loop_element[[1]]))
                } else {
                    stop("Function was not made to handle loop_element as list when not a formula.")
                }
            }

            for (subset_name in subset_names) {
                local_csv <- load_metric_scores(
                    run_path, subset_name, metric = NULL)
                subset_scores <- local_csv |>
                    mutate(species = species_names) |>
                    pivot_longer(
                        cols = c(RMSE, AUC, TjurR2),
                        names_to = "metric", values_to = "score") |>
                    mutate(
                        loop_element = loop_element, 
                        k_fold = k, 
                        dataset = subset_name)

                scores_df <- rbind(scores_df, subset_scores)
            }
        }
    }

    # Select the columns that should be kept in the final dataframe
    if (group_species) {
        # "species" column is not kept -> species are aggregated together
        group_cols <- c("metric", "loop_element", "dataset") 
    } else {
        group_cols <- c("metric", "loop_element", "dataset", "species")
        scores_df <- scores_df |> 
            mutate(species = factor(species, levels = species_names))
    }
    # aggregate rows together based on unique combinations of group_cols
    aggregated_df <- scores_df |>
        group_by(across(all_of(group_cols))) |>
        summarise(
            avg_score = mean(score, na.rm = TRUE),
            sd_score = sd(score, na.rm = TRUE),
            .groups = "drop_last") |>
        mutate(dataset = factor(dataset, levels = subset_names))

    cli_alert_info("Creating plot...")
    if (barplot) {
        aggregated_df <- aggregated_df |>
            mutate(loop_element = factor(loop_element, levels = loop_elements))
        p <- ggplot(
                aggregated_df,
                aes(y = avg_score, 
                    x = loop_element, 
                    fill = dataset,
                    ymin = avg_score - sd_score, 
                    ymax = avg_score + sd_score)) +
            geom_bar(
                stat = "identity", 
                position = position_dodge(width = 0.66), 
                width = 0.66) +
            geom_errorbar(
                position = position_dodge(width = 0.66), 
                width = 0.2)
    } else {
        p <- ggplot(
                aggregated_df,
                aes(y = avg_score, 
                    x = loop_element,
                    ymin = avg_score - sd_score, 
                    ymax = avg_score + sd_score)) +
            geom_ribbon(alpha = 0.33, aes(fill = dataset)) +
            geom_point(size = 1, aes(color = dataset)) +
            geom_line(aes(color = dataset))
    }

    if (group_species) {
        bottom_caption <- paste(
            "SD and mean computed per k_fold",
            "(species-averaged within each fold, across k_fold).")
    } else {
        bottom_caption <- paste(
            "SD and mean computed per k_fold across k_fold, per species.")
    }

    # format formula to number of variables
    temp_comb <- loop_model_combination
    temp_comb$HMSC_XFORMULAS <- lapply(lapply(
        loop_model_combination$HMSC_XFORMULAS, all.vars), length)
    
    # add model type to captions
    model_type <- ""
    for (param in names(temp_comb)) {
        if (length(temp_comb[[param]]) > 1) {
            model_type <- paste0(
                model_type, tolower(param), "=see x-axis, ")
        } else {
            model_type <- paste0(
                model_type, tolower(param), "=", temp_comb[[param]], ", ")
        }
    }
    model_type <- paste0(substr(model_type, 1, nchar(model_type)-3), ".")
    bottom_caption <- paste0(bottom_caption, "\nModel type: ", model_type)

    p <- p + 
        labs(caption = bottom_caption)

    if (group_species) {
            p <- p + facet_wrap(~ metric) 
        } else {
            p <- p + facet_grid(metric ~ species, scales = "free_x")
        }
    
    p <- my_custom_ggplot_theme(p) + 
        scale_fill_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
        scale_color_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
        xlab(xlabel) + 
        ylab(ylabel)

    if (barplot)
        p <- p + scale_x_discrete(labels = abbreviate_loop_labels)

    finalize_plot(p, save_to, what = "performances")
    return(list(plot = p))
}

# A function to compare scores between k_folds, subset and model type.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - reference_model_combination: a list of single parameters. 
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - loop_on: a string. The name of a parameter in a combination. 
#       The parameter which has several values.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - xlabel: a string. The xlabel for the plot (default is "Effect").
#   - group_species: whether to take the mean 
#       accross all k_folds and species (TRUE, default) or 
#       only accross k_folds (FALSE).
#   - species_names: names of species in CSV 
#       (rownames are not available from csvs).
#   - save_to: a string. Path which should end with .pdf. 
#       If NULL, does not save pdf.
boxplot_compare_scores <- function(
        parent_folder,
        reference_model_combination,
        loop_model_combination,
        loop_on,
        k_fold,
        metric = "MSE",
        subset_names = c("train", "val", "test"),
        xlabel = "Effect",
        group_species = TRUE,
        species_names = NULL,
        save_to = NULL) {

    # auto compute differences between reference scores and the other scores
    diffs <- compute_score_diffs(
        parent_folder, reference_model_combination, loop_model_combination,
        loop_on, k_fold, metric, subset_names, group_species, species_names)
    scores_df <- diffs$summary
    raw_diffs_df <- diffs$raw_diffs

    # load reference values for display in captions
    refs_list <- list()
    for (subset_name in subset_names) {
        for (k in 1:k_fold) {
            run_path <- make_run_path(
                parent_folder, reference_model_combination, k)
            refs <- load_metric_scores(run_path, subset_name, metric)

            refs_list[[length(refs_list) + 1]] <- data.frame(
                subset = subset_name,
                k_fold = k,
                scores = refs,
                species = species_names,
                stringsAsFactors = FALSE
            )
        }
    }
    refs_df <- bind_rows(refs_list)
    if (group_species) {
        mean_refs_df <- refs_df |>
            group_by(subset) |>
            summarise(mean_score = mean(scores, na.rm = TRUE), .groups = "drop")
    } else {
        mean_refs_df <- refs_df |>
            group_by(subset, species) |>
            summarise(mean_score = mean(scores, na.rm = TRUE), .groups = "drop")
    }

    # just to be sure
    if (!group_species) {
        scores_df <- scores_df |> 
            mutate(species = factor(species, levels = species_names))
        raw_diffs_df <- raw_diffs_df |> 
            mutate(species = factor(species, levels = species_names))
    }

    if (group_species) {
        bottom_caption <- paste(
            "Distribution of per-fold differences",
            "(species-averaged within each fold, across k_fold).")
    } else {
        bottom_caption <- paste(
            "Distribution of per-fold differences across k_fold, per species.")
    }

    # format formula to number of variables
    temp_loop <- loop_model_combination
    temp_loop$HMSC_XFORMULAS <- lapply(lapply(
        temp_loop$HMSC_XFORMULAS, all.vars), length)
    temp_ref <- reference_model_combination
    temp_ref$HMSC_XFORMULAS <- lapply(lapply(
        temp_ref$HMSC_XFORMULAS, all.vars), length)
    
    # add model type to caption
    ref_type <- ""
    model_type <- ""
    for (param in names(temp_loop)) {
        if (length(temp_loop[[param]]) > 1) {
            model_type <- paste0(
                model_type, tolower(param), "=see x-axis, ")
        } else {
            model_type <- paste0(
                model_type, tolower(param), "=", temp_loop[[param]], ", ")
        }
        ref_type <- paste0(
                ref_type, tolower(param), "=", temp_ref[[param]], ", ")
    }
    model_type <- paste0(substr(model_type, 1, nchar(model_type)-2), ".")
    ref_type <- paste0(substr(ref_type, 1, nchar(ref_type)-2), ".")
    bottom_caption <- paste0(
        bottom_caption,
        "\nReference: ", ref_type, 
        ".\nCompared with: ", model_type)
    
    p <- ggplot(
            raw_diffs_df, 
            aes(y = diff_value, x = loop_element, fill = subset)) +
        geom_boxplot(
            position = position_dodge(width = 0.75), 
            width = 0.6, 
            outlier.shape = 16) +
        labs(caption = bottom_caption, fill = "Subset") +
        ylab(paste("Delta in average", metric)) +
        xlab(xlabel) +
        geom_hline(yintercept = 0, linetype = "dashed")

    # add text for mean reference
    if (group_species) {
        ref_values <- paste0(mean_refs_df |> 
            mutate(subset = factor(subset, levels = subset_names)) |>
            arrange(subset) |>
            mutate(label = paste0(subset, ": ", round(mean_score, 3))) |>
            pull(label), collapse=", ")

        p <- p + annotation_custom(
            grob = grid::textGrob(
                paste0("Reference scores: ", ref_values, "."),
                x = unit(0, "npc"), y = unit(0, "npc"),
                hjust = -0.02, vjust = -0.75,
                gp = grid::gpar(
                    fontsize = 9, fontface = "italic", lineheight = 0.8)
            )
        )
    }

    if (!group_species) p <- p + facet_wrap(~species, nrow = 1)

    p <- my_custom_ggplot_theme(p)  +
        scale_fill_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
        scale_x_discrete(labels = abbreviate_loop_labels)

    if (grepl("MSE", metric)) {
        p <- p + scale_y_reverse()
    }

    finalize_plot(p, save_to, what = "data")
    return(list(scores = scores_df, diffs = raw_diffs_df, plot = p))
}

# A function to show the number of species with an improvement of 1% in their
# prediction scores, accross k_folds.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - reference_model_combination: a list of single parameters. 
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - loop_on: a string. The name of a parameter in a combination. 
#       The parameter which has several values.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - xlabel: a string. The xlabel for the plot (default is "Model type").
#   - ylabel: a string. The ylabel for the plot (default is "Average Score").
#   - proportion: a numeric between 0 and 1. The threshold to count a species
#       score as improved. 
#   - save_to: a string. Path which should end with .pdf. 
#       If NULL, does not save pdf.
boxplot_sp_improvements <- function(
        parent_folder,
        reference_model_combination,
        loop_model_combination,
        loop_on,
        k_fold,
        metric = "MSE",
        proportion = 1/100,
        subset_names = c("train", "val", "test"),
        xlabel = "Effect",
        save_to = NULL) {

    # auto compute differences between reference scores and the other scores
    diffs <- compute_number_of_improvements(
        parent_folder = parent_folder, 
        reference_model_combination = reference_model_combination, 
        loop_model_combination = loop_model_combination, 
        loop_on = loop_on,
        k_fold = k_fold, 
        metric = metric,
        subset_names = subset_names,
        species_names = NULL, 
        proportion = proportion)
    prop_diffs_df <- diffs$prop_diffs
    n_species <- diffs$n_species

    # When boxes have less than 5 values, they can appear squashed
    # for visualisation purposes, its a good idea to switch from fill to color 
    # under this threshold.
    if (n_species <= 5) {
        p <- ggplot(
            prop_diffs_df, 
            aes(y = improvements, x = loop_element, color = subset))
    } else {
        p <- ggplot(
                prop_diffs_df, 
                aes(y = improvements, x = loop_element, fill = subset))        
    }
 
    # format formula to number of variables
    temp_loop <- loop_model_combination
    temp_loop$HMSC_XFORMULAS <- lapply(lapply(
        temp_loop$HMSC_XFORMULAS, all.vars), length)
    temp_ref <- reference_model_combination
    temp_ref$HMSC_XFORMULAS <- lapply(lapply(
        temp_ref$HMSC_XFORMULAS, all.vars), length)

    # create captions
    ref_type <- ""
    model_type <- ""
    bottom_caption <- paste("Number of species per box:", n_species, ".")
    for (param in names(temp_loop)) {
        if (length(temp_loop[[param]]) > 1) {
            model_type <- paste0(
                model_type, tolower(param), "=see x-axis, ")
        } else {
            model_type <- paste0(
                model_type, tolower(param), "=", temp_loop[[param]], ", ")
        }
        ref_type <- paste0(
                ref_type, tolower(param), "=", temp_ref[[param]], ", ")
    }
    model_type <- paste0(substr(model_type, 1, nchar(model_type)-2), ".")
    ref_type <- paste0(substr(ref_type, 1, nchar(ref_type)-2), ".")
    bottom_caption <- paste0(
        bottom_caption,
        "\nReference: ", ref_type, 
        ".\nCompared with: ", model_type)
    
    p <- p +
        geom_boxplot(
            position = position_dodge(width = 0.75), 
            width = 0.6, outlier.shape = 16) +
        labs(
            caption = bottom_caption, 
            fill = "Subset") +
        ylab(paste0(
            "Number of species with \nimprovement >", proportion*100,
            "% in ", metric)) +
        xlab(xlabel)

    # switch from fill to color if less than 5 species (see higher in function)
    if (n_species <= 5) {
        p <- my_custom_ggplot_theme(p, LIGHT = TRUE)  +
            scale_color_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
            scale_x_discrete(labels = abbreviate_loop_labels)
    } else {
        p <- my_custom_ggplot_theme(p)  +
            scale_fill_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
            scale_x_discrete(labels = abbreviate_loop_labels)
    }

    finalize_plot(p, save_to, what = "data")
    return(list(plot = p))
}


# A function to make statistical comparison of scores, with a factorial effect.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#       If NULL, does not save pdf.
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - loop_on: a string. The name of a parameter in a combination. 
#       The parameter which has several values.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - xlabel: a string. The xlabel for the plot (default is "Model type").
#   - ylabel: a string. The ylabel for the plot (default is "Average Score").
#   - group_species: whether to take the mean 
#       accross all k_folds and species (TRUE, default) or 
#       only accross k_folds (FALSE).
#   - species_names: names of species in CSV 
#       (rownames are not available from csvs).
#   - save_to: a string. Path which should end with .pdf. 
dotwhisker_model_scores <- function(
        parent_folder,
        loop_model_combination,
        loop_on,
        k_fold,
        metric,
        subset_names = c("train", "val", "test"),
        xlabel = "Value",
        ylabel = "Effect",
        group_species = TRUE,
        species_names = NULL,
        save_to = NULL) {
            
    loop_elements <- loop_model_combination[[loop_on]]
    if (is.numeric(loop_elements)) {
        stop(paste0(
            "This function was made to handle factors in 'loop_elements',",
            " not numerics. For numerics, Try barplot_raw_scores() with ",
            "`barplot = FALSE` instead."))
    }

    # No auto import of results in this function (mainly because I didnt need
    # it at the time...)
    # The following loop import each necessary scores.csv file and parse it 
    # into a dataframe.
    run_paths <- make_run_path(
        parent_folder, loop_model_combination, "")
    cli_alert_info("Fetching scores...")
    scores_df <- data.frame()
    for (k in seq(k_fold)) {
        for (i in 1:length(run_paths)) {
            run_path <- paste0(run_paths[i], k)
            loop_element <- loop_elements[i]

            if (is.list(loop_element)) {
                if (is_formula(loop_element[[1]])) {
                    loop_element <- length(all.vars(loop_element[[1]]))
                } else {
                    stop("Function was not made to handle loop_element as list when not a formula.")
                }
            }

            for (subset_name in subset_names) {
                local_csv <- suppressMessages(load_metric_scores(
                    run_path, subset_name, metric = NULL))
                subset_scores <- local_csv |>
                    mutate(species = factor(species_names)) |>
                    pivot_longer(
                        cols = c(all_of(metric)),
                        names_to = "score_type", values_to = "score") |>
                    mutate(
                        loop_element = loop_element, 
                        k_fold = k, 
                        dataset = subset_name)

                scores_df <- rbind(scores_df, subset_scores)
            }
        }
    }

    # One pooled "group" if group_species, otherwise one group per species   
    if (!group_species) {
        scores_df <- scores_df |> 
            mutate(species = factor(species, levels = species_names))
    }

    # Ensure formatting of columns as factors
    scores_df <- scores_df |>
        mutate(k_fold = factor(k_fold)) |>
        mutate(loop_element = factor(loop_element, levels=loop_elements))

    # fit one model per data subset
    model_list <- list()
    for (subset in subset_names) {
        if (group_species) {
            model_list[[subset]] <- lmer(score ~ loop_element + (1 | k_fold) + (1 | species),
                data = subset(scores_df, dataset == subset)) 
        } else {
            for (sp in species_names) {
                key <- paste(subset, sp, sep = "_")
                model_list[[key]] <- lmer(
                    score ~ loop_element + (1 | k_fold),
                    data = subset(scores_df, dataset == subset & species == sp)
                )
            }
        }
    }

    # tidy the table if not grouping species
    if (!group_species) {
        tidy_list <- list()
        for (subset in subset_names) {
            for (sp in species_names) {
                key <- paste(subset, sp, sep = "_")
                tidy_list[[key]] <- broom.mixed::tidy(
                    model_list[[key]], effects = "fixed", conf.int = TRUE) |>
                    mutate(subset = subset, species = sp)
            }
        }
        model_list <- bind_rows(tidy_list) |>
            filter(term != "(Intercept)") |>
            mutate(model = subset)
    }

    if (group_species) {
        bottom_caption <- paste(
            "Comparison of scores per subset, k_fold mean + 95% CI.")
    } else {
        bottom_caption <- paste(
            "Comparison of scores per subset per species,",
            "k_fold mean + 95% CI.")
    }

    # format formula to number of variables
    temp_comb <- loop_model_combination
    temp_comb$HMSC_XFORMULAS <- lapply(lapply(
        loop_model_combination$HMSC_XFORMULAS, all.vars), length)

    # add model type to captions
    model_type <- ""
    for (param in names(temp_comb)) {
        if (length(temp_comb[[param]]) > 1) {
            model_type <- paste0(
                model_type, tolower(param), "=see x-axis, ")
        } else {
            model_type <- paste0(
                model_type, tolower(param), "=", temp_comb[[param]], ", ")
        }
    }
    model_type <- paste0(substr(model_type, 1, nchar(model_type)-3), ".")
    bottom_caption <- paste0(bottom_caption, "\nModel type: ", model_type)

    p <- dwplot(model_list, 
                effects = "fixed",
                model_order = rev(subset_names),
                vars_order = paste0("loop_element", rev(loop_elements[-1]))) |>
        relabel_predictors(setNames(
            rev(loop_elements[-1]),
            paste0("loop_element", rev(loop_elements[-1])))) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        coord_flip() +
        xlab(xlabel) +
        ylab(ylabel) +
        labs(caption = bottom_caption, fill = "Subset", color = "Subset") +
        scale_y_discrete(labels = abbreviate_loop_labels)

    if (!group_species) p <- p + facet_wrap(~species, nrow = 1)

    p <- my_custom_ggplot_theme(p, LIGHT = TRUE) +
        scale_color_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1]))

    if (grepl("MSE", metric)) {
        p <- p + scale_x_reverse()
    }

    finalize_plot(p, save_to, what = "dotwhiskers")
    return(list(plot=p, models=model_list))
}

# A function to make statistical comparisons of scores, with a numeric effect.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - xlabel: a string. The xlabel for the plot (default is "Model type").
#   - ylabel: a string. The ylabel for the plot (default is "Average Score").
#   - group_species: whether to take the mean 
#       accross all k_folds and species (TRUE, default) or 
#       only accross k_folds (FALSE).
#   - species_names: names of species in CSV 
#       (rownames are not available from csvs).
#   - save_to: a string. Path which should end with .pdf. 
#       If NULL, does not save pdf.
lineplot_model_scores <- function(
        parent_folder,
        loop_model_combination,
        loop_on,
        k_fold,
        metric,
        xlabel = "Model type",
        ylabel = "Average Score",
        fit_poly = 2,
        group_species = TRUE,
        species_names = NULL,
        subset_names = c("train", "val", "test"),
        save_to = NULL) {
    
    loop_elements <- loop_model_combination[[loop_on]]
    run_paths <- make_run_path(
        parent_folder, loop_model_combination, "")

    # No auto import of results in this function (mainly because I didnt need
    # it at the time...)
    # The following loop import each necessary scores.csv file and parse it 
    # into a dataframe.
    cli_alert_info("Fetching scores...")
    scores_df <- data.frame()
    for (k in seq(k_fold)) {
        for (i in 1:length(run_paths)) {
            run_path <- paste0(run_paths[i], k)
            loop_element <- loop_elements[i]

            if (is.list(loop_element)) {
                if (is_formula(loop_element[[1]])) {
                    loop_element <- length(all.vars(loop_element[[1]]))
                } else {
                    stop("Function was not made to handle loop_element as list when not a formula.")
                }
            }

            for (subset_name in subset_names) {
                local_csv <- suppressMessages(load_metric_scores(
                    run_path, subset_name, metric = NULL))
                subset_scores <- local_csv |>
                    mutate(species = factor(species_names)) |>
                    pivot_longer(
                        cols = c(all_of(metric)),
                        names_to = "score_type", values_to = "score") |>
                    mutate(
                        loop_element = loop_element, 
                        k_fold = k, 
                        dataset = subset_name)

                scores_df <- rbind(scores_df, subset_scores)
            }
        }
    }

    # One pooled "group" if group_species, otherwise one group per species   
    if (!group_species) {
        scores_df <- scores_df |> 
            mutate(species = factor(species, levels = species_names))
    }

    ### Aggregate
    model_list <- list()
    for (subset in subset_names) {
        model_list[[subset]] <- lmer(
            score ~ poly(loop_element, fit_poly) + (1 | k_fold) + (1 | species),
            data = subset(scores_df, dataset == subset)) 
        # # The model will estimate the average effect of each value of loop_element 
        # # on the metric, with a random effect of k_fold and species
        # cat(paste(toupper(subset), "\n"))
        # print(summary(model_list[[subset]])$coefficients)
        # cat("\n")
    }

    # Fitted lines
    if (group_species) {
        pred_df <- bind_rows(lapply(names(model_list), function(subset) {
            ggpredict(model_list[[subset]], terms = "loop_element [all]") |>
                as.data.frame(terms_to_colnames = TRUE) |>
                mutate(dataset = subset)
        }))
    } else {
        pred_df <- bind_rows(lapply(names(model_list), function(subset) {
            ggpredict(model_list[[subset]],
                    terms = c("loop_element", "species"),
                    type = "random") |>
                as.data.frame(terms_to_colnames = TRUE) |>
                mutate(dataset = subset)
        }))
    }

    pred_df <- pred_df |>
        mutate(dataset = factor(dataset, levels = subset_names))

    # Scatter plot of means
    if (group_species) {
        group_these_columns <- c("loop_element", "dataset")
    } else {
        group_these_columns <- c("loop_element", "dataset", "species")
    }
    aggregated_df <- scores_df |>
        group_by(across(all_of(group_these_columns))) |>
        summarise(
            avg_score = mean(score, na.rm = TRUE),
            .groups = "drop_last") |>
        mutate(dataset = factor(dataset, levels = subset_names))
    
    cli_alert_info("Creating plot...")

    p <- ggplot() +
        # Ribbon + line from the model predictions
        geom_ribbon(data = pred_df,
                    aes(x = loop_element, ymin = conf.low, ymax = conf.high, fill = dataset),
                    alpha = 0.33) +
        geom_line(data = pred_df,
                    aes(x = loop_element, y = predicted, color = dataset),
                    linewidth = 1) +
        
        # Points + line from the aggregated data
        geom_point(data = aggregated_df,
                    aes(x = loop_element, y = avg_score, color = dataset),
                    size = 1) +
        geom_line(data = aggregated_df,
                    aes(x = loop_element, y = avg_score, color = dataset)) +
        labs(x = "loop_element", y = "Predicted metric",
            color = "Subset", fill = "Subset") +
        theme_minimal()

    if (group_species) {
        bottom_caption <- paste(
            "Comparison of scores per subset, k_fold mean + 95% CI.")
    } else {
        bottom_caption <- paste(
            "Comparison of scores per subset, per species,",
            "k_fold mean + 95% CI.")
    }

    model_type <- ""
    temp_comb <- loop_model_combination
    temp_comb$HMSC_XFORMULAS <- lapply(lapply(
        loop_model_combination$HMSC_XFORMULAS, all.vars), length)
    for (param in names(temp_comb)) {
        if (length(temp_comb[[param]]) > 1) {
            model_type <- paste0(
                model_type, tolower(param), "=see x-axis, ")
        } else {
            model_type <- paste0(
                model_type, tolower(param), "=", temp_comb[[param]], ", ")
        }
    }
    model_type <- paste0(substr(model_type, 1, nchar(model_type)-3), ".")

    bottom_caption <- paste0(bottom_caption, "\nModel type: ", model_type)

    p <- p + 
        labs(caption = bottom_caption)

    if (!group_species) {
            p <- p + facet_grid( ~ species, scales = "fixed")
        }
    
    p <- my_custom_ggplot_theme(p) + 
        scale_fill_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
        scale_color_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
        xlab(xlabel) + 
        ylab(ylabel) +
        guides(fill = "none") 

    finalize_plot(p, save_to, what = "performances")
    return(list(plot= p, data_df=aggregated_df, models=model_list))
}

# A function to make statistical comparisons of scores depending on a numeric
# and a qualitative effect.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - loop_model_combinations: list of lists of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - type_loop_on: a string. A name within a combination: category displayed 
#       as different colors.
#   - x_loop_on: a string. A name within a combination: category (numeric) 
#       displayed on the x-axis.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - relative_diff: a boolean. Whether to display the results as raw scores 
#       or relative differences.
#   - fit_poly: an integer. The value passed to poly(loop_element, fit_poly).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - xlabel: a string. The xlabel for the plot (default is "Model type").
#   - ylabel: a string. The ylabel for the plot (default is "Average Score").
#   - species_names: names of species in CSV 
#       (rownames are not available from csvs).
#   - show_n_improved: a boolean. If TRUE (default) shows the number of species
#       that improved by more than `improvement_threshold`. 
#       Overwrites relative_diff
#   - improvement_threshold: a float. When show_n_improved, the improvement in 
#       percentage to count a prediction as being better than baseline.
#       Default is 0.01 (1%).
#   - save_to: a string. Path which should end with .pdf. 
#       If NULL, does not save pdf.
multi_lineplot_model_scores <- function(
        parent_folder,
        loop_model_combinations,
        type_loop_on,
        x_loop_on,
        k_fold,
        metric,
        xlabel = "Model type",
        ylabel = "Average Score",
        species_names = NULL,
        subset_names = c("train", "val", "test"),
        fit_poly = 2,
        relative_diff = FALSE,
        improvement_threshold = 0.01,
        show_n_improved = TRUE,
        save_to = NULL) {

    # "Improvement" is only meaningful for a lower-is-better metric like MSE.
    # If you need this for a higher-is-better metric, flip the sign in the
    # rel_improvement calculation below (or set show_n_improved = FALSE to
    # fall back to the original average-score plot).
    if (show_n_improved && metric != "MSE") {
        stop("show_n_improved assumes a lower-is-better metric (MSE). ",
             "Set show_n_improved = FALSE to plot other metrics as before.")
    }

    # test if we got several models (at least 2)
    if (!all(
            names(loop_model_combinations[[1]]) ==
            names(loop_model_combinations[[2]]))) {
        stop("loop_model_combinations was not initialised correctly.")
    }

    # loop on in combination to get scores
    overall_pred_df <- NULL
    type_order <- NULL
    all_lmer_models <- NULL
    all_scores <- NULL

    for (j in 1:length(loop_model_combinations)) {
        loop_elements <- loop_model_combinations[[j]][[x_loop_on]]
        run_paths <- make_run_path(
            parent_folder, loop_model_combinations[[j]], "")
        model_type <- loop_model_combinations[[j]][[type_loop_on]]
        type_order <- c(type_order, model_type)

        # No auto import of results in this function (mainly because I didnt need
        # it at the time...)
        # The following loop import each necessary scores.csv file and parse it
        # into a dataframe.
        cli_alert_info("Fetching scores...")
        scores_df <- data.frame()
        for (k in seq(k_fold)) {
            for (i in 1:length(run_paths)) {
                run_path <- paste0(run_paths[i], k)
                loop_element <- loop_elements[i]

                if (is.list(loop_element)) {
                    if (is_formula(loop_element[[1]])) {
                        loop_element <- length(all.vars(loop_element[[1]]))
                    } else {
                        stop("Function was not made to handle loop_element as list when not a formula.")
                    }
                }

                for (subset_name in subset_names) {
                    local_csv <- suppressMessages(load_metric_scores(
                        run_path, subset_name, metric = NULL))
                    subset_scores <- local_csv |>
                        mutate(species = factor(species_names)) |>
                        pivot_longer(
                            cols = c(all_of(metric)),
                            names_to = "score_type", values_to = "score") |>
                        mutate(
                            loop_element = loop_element,
                            k_fold = k,
                            dataset = subset_name,
                            type = model_type)

                    scores_df <- rbind(scores_df, subset_scores)
                }
            }
        }

        # --- Per-species relative improvement vs a baseline ---
        # Baseline = each species' own score at the smallest loop_element
        # within its (species, k_fold, dataset, type) group (same reference
        # point the old relative_diff block used). "improved" flags species
        # whose MSE dropped by more than improvement_threshold (fraction of
        # baseline), e.g. 0.01 = 1%.
        scores_df <- scores_df |>
            group_by(species, k_fold, dataset, type) |>
            mutate(
                baseline = score[which.min(loop_element)],
                rel_improvement = (baseline - score) / baseline,
                improved = rel_improvement > improvement_threshold) |>
            ungroup()

        # Old behaviour preserved: turn score into an absolute difference
        # from baseline when relative_diff = TRUE and show_n_improved = FALSE.
        if (relative_diff) {
            scores_df <- scores_df |>
                mutate(score = score - baseline)
        }
        scores_df <- scores_df |> select(-baseline)

        # Drop the baseline loop_element itself when counting/modelling
        # "improved" species: by construction rel_improvement == 0 there,
        # so it's not a real data point for this metric and would bias the
        # model fit (and the plot) toward a forced zero at the left edge.
        if (show_n_improved) {
            scores_df <- scores_df |>
                group_by(species, k_fold, dataset, type) |>
                filter(loop_element != loop_element[which.min(loop_element)]) |>
                ungroup()
        }
       
        all_scores <- bind_rows(all_scores, scores_df)

        ### Aggregate
        model_list <- list()
        if (show_n_improved) {
            # Collapse species into a per-(loop_element, k_fold, dataset)
            # count of "improved" species, then model that count across
            # loop_element. Species can no longer be a random effect since
            # it has been summed away; k_fold remains one.
            n_improved_df <- scores_df |>
                group_by(loop_element, k_fold, dataset, type) |>
                summarise(n_improved = sum(improved, na.rm = TRUE), .groups = "drop")

            for (subset in subset_names) {
                # NOTE: a Poisson/negative-binomial GLMM (e.g. glmer with
                # family = poisson) is statistically a more natural
                # choice for count data than lmer
                model_list[[subset]] <- lmer(
                    n_improved ~ poly(loop_element, fit_poly) + (1 | k_fold),
                    data = subset(n_improved_df, dataset == subset))
            }
        } else {
            for (subset in subset_names) {
                model_list[[subset]] <- lmer(
                    score ~ poly(loop_element, fit_poly) + (1 | k_fold) + (1 | species),
                    data = subset(scores_df, dataset == subset))
            }
        }
        all_lmer_models <- c(all_lmer_models, model_list)

        # Fitted lines
        pred_df <- bind_rows(lapply(names(model_list), function(subset) {
            ggpredict(model_list[[subset]], terms = "loop_element [all]") |>
                as.data.frame(terms_to_colnames = TRUE) |>
                mutate(
                    dataset = subset,
                    type = model_type)
        }))

        pred_df <- pred_df |>
            mutate(dataset = factor(dataset, levels = subset_names))

        overall_pred_df <- bind_rows(overall_pred_df, pred_df)
    }

    overall_pred_df <- overall_pred_df |>
        mutate(type = factor(type, levels = type_order))

    # Scatter plot of means
    group_these_columns <- c("loop_element", "dataset", "type")

    if (show_n_improved) {
        # Per-point value = number of species improved > threshold in that
        # fold, averaged across k_folds.
        aggregated_df <- all_scores |>
            group_by(loop_element, k_fold, dataset, type) |>
            summarise(n_improved = sum(improved, na.rm = TRUE), .groups = "drop") |>
            group_by(across(all_of(group_these_columns))) |>
            summarise(avg_score = mean(n_improved, na.rm = TRUE), .groups = "drop_last") |>
            mutate(dataset = factor(dataset, levels = subset_names)) |>
            mutate(type = factor(type, levels = type_order))
    } else {
        aggregated_df <- all_scores |>
            group_by(across(all_of(group_these_columns))) |>
            summarise(
                avg_score = mean(score, na.rm = TRUE),
                .groups = "drop_last") |>
            mutate(dataset = factor(dataset, levels = subset_names)) |>
            mutate(type = factor(type, levels = type_order))
    }

    cli_alert_info("Creating plot...")

    p <- ggplot() +
        # Ribbon model predictions
        geom_ribbon(data = overall_pred_df,
                    aes(
                        x = loop_element,
                        ymin = conf.low, ymax = conf.high,
                        fill = type),
                    alpha = 0.25) +
        # # line from model predictions (hidden)
        # geom_line(data = overall_pred_df,
        #             aes(x = loop_element, y = predicted, color = type),
        #             linewidth = 1) +

        # Points + line from the aggregated data
        geom_point(data = aggregated_df,
                    aes(x = loop_element, y = avg_score, color = type),
                    size = 1) +
        geom_line(data = aggregated_df,
                    aes(x = loop_element, y = avg_score, color = type)) +
        labs(x = "loop_element",
             y = if (show_n_improved) "N species improved" else "Predicted metric",
             color = "Model type", fill = "Model type") +
        theme_minimal()

    bottom_caption <- if (show_n_improved) {
        paste0(
            "Number of species with a MSE improvement > ",
            improvement_threshold * 100,
            "% vs baseline.\nMean over all k_folds, per subset. 95% CI computed with quadratic lmer.")
    } else {
        paste("Comparison of means over all k_folds, per subset. 95% CI computed with quadratic lmer.")
    }

    p <- p +
        labs(caption = bottom_caption) +
        facet_grid( ~ dataset, scales = "fixed")

    p <- my_custom_ggplot_theme(p, with_palette = TRUE, LIGHT = TRUE) +
        xlab(xlabel) +
        ylab(if (show_n_improved) "Number of species improved" else ylabel) +
        guides(fill = "none")

    if ((relative_diff) & (metric == "MSE") & !show_n_improved) {
        p <- p + scale_y_reverse()
    }

    finalize_plot(p, save_to, what = "performances")
    return(list(plot = p, data_df = aggregated_df, models = all_lmer_models))
}