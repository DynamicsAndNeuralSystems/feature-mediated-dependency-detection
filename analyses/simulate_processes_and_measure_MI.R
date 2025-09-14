#-----------------------------------------
#Input parameters
#-----------------------------------------

# Generative process settings
# source_type from main.R script
ts_length <- as.numeric(ts_length)

# Driving feature settings
driving_feature_timescale <- as.numeric(driving_feature_timescale)

# Capturing feature settings
capturing_feature_timescale <- driving_feature_timescale 

# Target generation 
linear_cor <- 0.4 # Coefficient of target process and source feature 
lag_window_source_target <- 1
source_target_function <- "linear"

#Estimator (Kraskov)
theiler_window_multiplier <- 1

#zcoring features
zscore <- T

#seed 
master_seed <- 42

#number of simulations
no_simulations <- 50

#result_file
if (!exists("seed")){
  result_file <- paste0(here("results","simulation_studies"),"/new_results.csv")
} else {
  result_file <- paste0(here("results","simulation_studies"),"/results_tslength", ts_length, "process",source_type, "driving_timescale", driving_feature_timescale,"seed",seed,".csv")
}
#plot folder 
plot_folder <- here("results","simulation_studies","plots")

#-----------------------------------------
# Simulations
#-----------------------------------------
if (file.exists(result_file)){
  results <- read.csv(result_file,header=T)
} else {
  results <- data.frame(seed=integer(),
                        process=character(),
                        ts_length=integer(),
                        driving_feature=character(),
                        driving_feature_timescale = integer(),
                        capturing_feature=character(),
                        capturing_feature_timescale=integer(),
                        linear_cor=numeric(),
                        inference_type=character(),
                        TE=numeric(),
                        pValue=numeric()
  )
  write.csv(results,file=result_file,row.names = F)
}


simulate_process_and_compute_TE <- function(seed, source_type, ts_length, driving_feature_timescale, 
                                            capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier){
  print(seed)
  set.seed(seed)
  source_ts <- generate_source_var(source_type, ts_length)
  if (source_type=="random"){
    driving_features <- generate_catch22_features_sliding_window(source_ts,driving_feature_timescale)
    for (driving_feature_name in names(driving_features)){
      driving_feature <- driving_features[[driving_feature_name]]
      if (sd(driving_feature)!=0){
        new_results <- simulate_from_driving_feature(seed, source_ts, source_type, ts_length, driving_feature, driving_feature_name, driving_feature_timescale, 
                                      capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier)
      }
    }
  } else if (source_type=="bimodal_spiking"){
    driving_feature_name <- "spike_count"
  } else if (source_type=="varying_skewness"){
    driving_feature_name <- "skewness"
  } else if (source_type=="FMW"){
    driving_feature_name <- "dominant_frequency"
  } else if (source_type=="AR3"){
    driving_feature_name <- "autocorrelation_lag1"
  } 
  if (source_type!="random"){
    driving_feature <- generate_feature_sliding_window(source_ts,driving_feature_timescale,feature_name=driving_feature_name)
    new_results <-  simulate_from_driving_feature(seed, source_ts, source_type, ts_length, driving_feature, driving_feature_name, driving_feature_timescale, 
                                             capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier)
    
  }
  write.table(new_results,result_file,append=T, row.names = F, col.names=F,sep=",")
  print(paste("Processed seed",seed,"linear_cor", linear_cor))
}

simulate_from_driving_feature <- function(seed, source_ts, source_type, ts_length, driving_feature, driving_feature_name, driving_feature_timescale, 
                                           capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier){
  new_results <- data.frame(seed=integer(),
                            process=character(),
                            ts_length=numeric(),
                            driving_feature=character(),
                            driving_feature_timescale=integer(),
                            capturing_feature=character(),
                            capturing_feature_timescale=integer(),
                            linear_cor=numeric(),
                            inference_type = character(),
                            TE = numeric(),
                            pValue=numeric())
  
  dest_ts <- build_target(source_target_function=source_target_function, ts=driving_feature,linear_cor = linear_cor)
  theiler_window_signal <- max(find_theiler_window(source_ts), find_theiler_window(dest_ts)) 
  sMI <- compute_TE(theiler_window_signal, tail(source_ts,length(dest_ts)),dest_ts, 0,1,capturing_feature_timescale,1,1, 
                    theiler_window_multiplier=theiler_window_multiplier, rotating_surrogates=T) #mutual information between source and target 
  new_results <- rbind.data.frame(new_results, data.frame(seed=seed,
                                                  process=source_type,
                                                  ts_length=ts_length,
                                                  driving_feature=driving_feature_name,
                                                  driving_feature_timescale=driving_feature_timescale,
                                                  capturing_feature="MI_s",
                                                  capturing_feature_timescale=capturing_feature_timescale,
                                                  linear_cor=linear_cor,
                                                  inference_type="MI_s",
                                                  TE=sMI$TE,
                                                  pValue=sMI$pValue))
  nMI <- compute_TE(theiler_window_signal, rnorm(length(dest_ts)),dest_ts,0,1,1,1,1, 
                    theiler_window_multiplier=theiler_window_multiplier) #mutual information between a random vector and target 
  new_results <- rbind.data.frame(new_results, data.frame(seed=seed,
                                                  process=source_type,
                                                  ts_length=ts_length,
                                                  driving_feature=driving_feature_name,
                                                  driving_feature_timescale=driving_feature_timescale,
                                                  capturing_feature="null",
                                                  capturing_feature_timescale=capturing_feature_timescale,
                                                  linear_cor=linear_cor,
                                                  inference_type="null",
                                                  TE=nMI$TE,
                                                  pValue=nMI$pValue))
  capturing_features <- generate_catch22_features_sliding_window(source_ts,capturing_feature_timescale)
  if (!driving_feature_name%in%names(capturing_features)){
    driving_feature_as_capturing_feature <- generate_feature_sliding_window(source_ts,capturing_feature_timescale,feature_name=driving_feature_name)
    capturing_features[[driving_feature_name]] <- driving_feature_as_capturing_feature # add the driving feature to the capturing feature list to evaluate performance when we have capturing feature match the driving feature
  }
   for (feature in names(capturing_features))
  {
    feature_ts <- capturing_features[[feature]]
    theiler_window_feature <- max(find_theiler_window(feature_ts), find_theiler_window(dest_ts)) 
    
    if (driving_feature_timescale >= capturing_feature_timescale){
      feature_ts_capped <- tail(feature_ts,length(dest_ts))
      dest_ts_capped <- dest_ts
    } else {
      dest_ts_capped <- tail(dest_ts,length(feature_ts))
      feature_ts_capped <- feature_ts
    }
    fMI <- compute_TE(theiler_window_feature, feature_ts_capped, dest_ts_capped, 0,1,1,1,1,
                      number_surrogates = 1000, theiler_window_multiplier =theiler_window_multiplier) #mutual information between source feature and target
    new_results <- rbind.data.frame(new_results, data.frame(seed=seed,
                                                    process=source_type,
                                                    ts_length=ts_length,
                                                    driving_feature=driving_feature_name,
                                                    driving_feature_timescale=driving_feature_timescale,
                                                    capturing_feature=feature,
                                                    capturing_feature_timescale=capturing_feature_timescale,
                                                    linear_cor=linear_cor,
                                                    inference_type="MI_f",
                                                    TE=fMI$TE,
                                                    pValue=fMI$pValue))
   }
  return(new_results)
}

simulate_multiple_runs <- function(source_type, ts_length, driving_feature_timescale, 
                                   capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier, no_simulations){
  set.seed(master_seed)
  simulation_seeds <- sample.int(100,no_simulations)
  for (seed in simulation_seeds){
    if (nrow(results[results$process==source_type & 
                     results$ts_length==ts_length & 
                     results$driving_feature_timescale==driving_feature_timescale & 
                     results$capturing_feature_timescale==capturing_feature_timescale &
                     results$linear_cor==linear_cor &
                     results$seed==seed,])==0){
      simulate_process_and_compute_TE(seed, source_type, ts_length, driving_feature_timescale, 
                                                 capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier)
    }
  }
}


if (!exists("seed")){
  # for (capturing_feature_timescale in seq(50,200,25))
  # {
  #   print(paste0("capturing_feature_timescale:",capturing_feature_timescale))
  #   simulate_multiple_runs(source_type, ts_length, driving_feature_timescale,
  #                          capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier, no_simulations)
  # 
  # }
  # simulate_multiple_runs(source_type, ts_length, driving_feature_timescale,
  #                        capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier, no_simulations)

  for (linear_cor in c(0.5, 0.6, 0.7, 0.8, 0.9 ))
  {
    print(paste0("linear_cor:",linear_cor))
    simulate_multiple_runs(source_type, ts_length, driving_feature_timescale,
                                      capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier, no_simulations)

  }
  
  # for (driving_feature_timescale in c(50, 150, 200))
  # {
  #   print(paste0("driving_feature_timescale:",driving_feature_timescale))
  #   linear_cor <- 0.4
  #   capturing_feature_timescale <- driving_feature_timescale
  #   simulate_multiple_runs(source_type, ts_length, driving_feature_timescale,
  #                          capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier, no_simulations)
  #   
  # }
  
} else if (nrow(results[results$process==source_type & 
                           results$ts_length==ts_length & 
                           results$driving_feature_timescale==driving_feature_timescale & 
                           results$capturing_feature_timescale==capturing_feature_timescale &
                           results$linear_cor==linear_cor &
                           results$seed==seed,])==0) { 
  simulate_process_and_compute_TE(seed, source_type, ts_length, driving_feature_timescale, 
                                  capturing_feature_timescale, source_target_function, linear_cor, theiler_window_multiplier)
} 
