rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  # Define a function named rand_vect that takes the following arguments:
  # - N: Length of the vector
  # - M: Desired sum of the vector
  # - sd: Standard deviation for random normal distribution (default is 1)
  # - pos.only: Boolean indicating if only positive values are allowed (default is TRUE)
  
  # Generate a random vector of length N from a normal distribution with mean M/N and standard deviation sd
  vec <- rnorm(N, M/N, sd)
  
  # If the absolute sum of the vector is very close to 0 (within 0.01), add 1 to every element of the vector
  if (abs(sum(vec)) < 0.01) vec <- vec + 1 
  
  # Scale the vector to ensure that it sums up to M
  vec <- round(vec / sum(vec) * M) 
  
  # Calculate the deviation from the desired sum M
  deviation <- M - sum(vec)
  
  # Adjust the vector to match the desired sum M
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  
  # If pos.only is TRUE, ensure that all elements are positive
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

find_theiler_window <- function(ts){
  autocors <- acf(ts,lag.max=length(ts),plot=F)
  #theiler_window <- min(which(autocors$acf<0),floor(length(ts)/10)) - 2
  theiler_window <- min(which(autocors$acf<0)) -2 
  return(theiler_window)
}

count_spikes <- function(ts){
  threshold <- mean(ts) + 3*sd(ts)
  spikes <- 0
  for (x in ts){
    if (x > threshold){
      spikes <- spikes + 1
    }
  }
  return(spikes)
}

count_spike_by_fixed_threshold <- function(ts,threshold){
  spikes <- 0
  for (x in ts){
    if (x > threshold){
      spikes <- spikes + 1
    }
  }
  return(spikes)
}

set_params_default_values <- function(source_type, params) {
  params$noise_coef <- ifelse(is.null(params$noise_coef), 0.1, params$noise_coef)
  if (source_type == "auto_correlated") {
    params$ac_source_coef <- ifelse(is.null(params$ac_source_coef), 0.5, params$ac_source_coef)
  } else if (source_type == "noisy_sine") {
    params$offset <- ifelse(is.null(params$offset), 0.05, params$offset)
  } else if (source_type == "gaussian_bimodel") {
    params$alpha <- ifelse(is.null(params$alpha), 0.4, params$alpha)
    params$gs_mean <- ifelse(is.null(params$gs_mean), 4, params$gs_mean)
    params$gs_sd <- ifelse(is.null(params$gs_sd), 2, params$gs_sd)
  } else if (source_type == "AR3") {
    params$number_segments <- ifelse(is.null(params$number_segments), 10, params$number_segments)
    params$phi1_start <- ifelse(is.null(params$phi1_start), 0.8, params$phi1_start)
    params$ar3_step_size <- ifelse(is.null(params$ar3_step_size), 0.0005, params$ar3_step_size)
    params$phi2 <- ifelse(is.null(params$phi2), 0.1, params$phi2)
    params$phi3 <- ifelse(is.null(params$phi3), 0.1, params$phi3)
  } else if (source_type == "non_stationary_waveform") {
    params$wave_step_size <- ifelse(is.null(params$wave_step_size), 1, params$wave_step_size)
    params$step_size <- ifelse(is.null(params$step_size), 0.005, params$step_size)
  } else if (source_type == "FMW") {
    params$wave_step_size <- ifelse(is.null(params$wave_step_size), 1, params$wave_step_size)
    params$step_size <- ifelse(is.null(params$step_size), 0.000001, params$step_size)
    params$number_segments <- ifelse(is.null(params$number_segments), 1, params$number_segments)
  } else if (source_type == "bimodal_spiking") {
    params$mode_switch_rate <- ifelse(is.null(params$mode_switch_rate), 0.01, params$mode_switch_rate)
    params$spike_rate1 <- ifelse(is.null(params$spike_rate1), 0.1, params$spike_rate1)
    params$spike_rate2 <- ifelse(is.null(params$spike_rate2), 0.05, params$spike_rate2)
    params$delta <- ifelse(is.null(params$delta), 5, params$delta)
  } 
  return(params)
}

generate_source_var <- function(source_type, ts_length, params=NULL){
  params <- set_params_default_values(source_type, params)
  
  if (source_type=="random"){
    source <- rnorm(ts_length)
  } 
  
  else if (source_type=="auto_correlated"){
    source_1 <- rnorm(1)
    source <- c(source_1)
    for (i in seq(1,ts_length-1)){
      source_t <- ac_source_coef*source[length(source)] + params$noise_coef*rnorm(1)
      source <- c(source, source_t)
    } 
  } 
  
  else if (source_type=="noisy_sine"){
    source <- (1-params$noise_coef) * sin(seq(1,ts_length)) + params$noise_coef*rnorm(ts_length) + params$offset
  } 
  
  else if (source_type=="gaussian_bimodel"){
    source <- (1-params$alpha)*rnorm(ts_length) + params$alpha*rnorm(params$gs_mean,params$gs_sd)
  } 
  
  else if (source_type=="AR3"){
    source <- rnorm(3)
    # phi1 <- c()
    # phi1_segment_start <- params$phi1_start
    # for (i in seq(1,params$number_segments,1)){
    #   segment_length <- floor(ts_length/params$number_segments)
    #   if (i %% 2 ==0){
    #     phi1_seg <- seq(phi1_segment_start, phi1_segment_start-params$ar3_step_size*(segment_length-1), -1*params$ar3_step_size)
    #     phi1_segment_start <- phi1_segment_start-params$ar3_step_size*(segment_length-1)
    #   } else{
    #     phi1_seg <- seq(phi1_segment_start, phi1_segment_start+params$ar3_step_size*(segment_length-1)/10, 1*params$ar3_step_size/10)
    #     phi1_segment_start <- phi1_segment_start+params$ar3_step_size*(segment_length-1)/10
    #   }
    #   phi1 <- c(phi1,phi1_seg)
    # }
    phi1 <- c()
    phi1_lengths <- rand_vect(params$number_segments,ts_length, sd=20)
    phi1_unique_values <- runif(min=0.3,max=0.8,params$number_segments)
    for (j in seq(1,params$number_segments,1)){
      phi1 <- c(phi1,rep(phi1_unique_values[j],phi1_lengths[j]))
    }
    for (i in seq(1,ts_length-3)){
      source_t <- phi1[i]*source[length(source)] + params$phi2*source[length(source)-1] + params$phi3*source[length(source)-2] + params$noise_coef*rnorm(1)
      source <- c(source, source_t)
    }
  } 
  
  else if (source_type=="non_stationary_waveform"){
    if(params$wave_variating_factor=="amplitude"){
      amplitude_coefs <- seq(1,1+params$step_size*(ts_length-1),params$step_size)
      source <-  amplitude_coefs * sin(seq(0,(ts_length-1)*params$wave_step_size,params$wave_step_size))+ params$noise_coef*rnorm(ts_length) 
    } 
    
  else if (wave_variating_factor=="period"){
      amplitude <- rnorm(1)
      period_coefs <- seq(1,1+params$step_size*(ts_length-1),params$step_size)
      source <-  amplitude*sin(period_coefs*seq(0,(ts_length-1)*params$wave_step_size, params$wave_step_size))+ params$noise_coef*rnorm(ts_length)
    } 
  } 
  
  else if (source_type=="FMW"){
    segment_len <- ts_length/params$number_segments
    source <- c()
    freq_coefs_all <- c()
    integrated_freqs <- c()
    for (i in 1:params$number_segments){
      if (i%%2 == 1){
        freq_coefs <- seq(0.01,0.01+params$step_size*(segment_len-1),params$step_size)
      } else {
        freq_coefs <- seq(0.01+params$step_size*(segment_len-1),0.01,-1*params$step_size)
      }
      freq_coefs_all <- c(freq_coefs_all,freq_coefs)
      for (j in seq(1,segment_len,1)){
        if (j==1){
          integrated_freqs[j] <- 0 
        } else {
          integrated_freqs[j] <- integrated_freqs[j-1]+ (freq_coefs[j]+freq_coefs[j-1])/2 
        }
      }
      segment <- cos(2*pi*(seq(0,(segment_len-1)*params$wave_step_size,params$wave_step_size)+integrated_freqs))+ params$noise_coef*rnorm(segment_len)
      source <- c(source, segment)
    }
    
  } 
  
  else if (source_type=="noisy_waveform"){
    if(wave_variating_factor=="amplitude"){
      amplitude_coefs <- rnorm(ts_length)
      source <- amplitude_coefs * sin(seq(0,(ts_length-1)*params$wave_step_size,params$wave_step_size)) + params$offset
    } 
  } 
  
  else if (source_type=="step_wise_variance"){
    segment <- c(rep(0,params$segment_len/2),rep(1,params$segment_len/2))
    source <- rnorm(ts_length)*rep(segment,ts_length/length(segment)) + params$noise_coef*rnorm(ts_length)
  } 
  
  else if (source_type=="bimodal_spiking"){
    p <- params$mode_switch_rate
    p0_t <- 0
    p0s <- c()
    for (i in seq(1,ts_length,1)){
      p0_t <- (1-rbern(1,p))*p0_t + rbern(1,p)*(1-p0_t)
      p0s <- c(p0s,p0_t)
    }
    p0s <- ifelse(p0s==1,params$spike_rate1,params$spike_rate2)
    
    v_t <- 0
    v <- c()
    source <- c()
    for (i in seq(1,ts_length,1)){
      p0 <- p0s[i]
      p1 <- 1 - p0
      if (v_t==0){
        v_t <- (1-rbern(1,p0))*v_t + rbern(1,p0)*(1-v_t)
      } else {
        v_t <- (1-rbern(1,p1))*v_t + rbern(1,p1)*(1-v_t)
      }
      x_t <- rnorm(1) + v_t*params$delta 
      source <- c(source,x_t)
      v <- c(v,v_t)
    }
    #note: after some careful re-math work, this v is just B(1,p0s). The convoluted maths is credited to 3am. 
  }
  return(source)
}

build_target <- function(source_target_function,ts, lag_window_source_target=0, linear_cor=0.8, ac_target_coef=0.4,zscore=T){
  if (source_target_function == "identity_shift") {
    target = c(rep(0,lag_window_source_target),ts[1:(length(ts)-lag_window_source_target)]) + rnorm(length(ts))
  } else if (source_target_function == "cos"){
    target = c(rep(0,lag_window_source_target),cos(ts[1:(length(ts)-lag_window_source_target)])) + rnorm(length(ts))
  } else if (source_target_function == "sin"){
    target = c(rep(0,lag_window_source_target),sin(ts[1:(length(ts)-lag_window_source_target)])) + rnorm(length(ts))
  } else if (source_target_function == "linear"){
    if (zscore==T){
      ts <- (ts-mean(ts))/sd(ts)
    }
    target = c(rep(0,lag_window_source_target), linear_cor*ts[1:(length(ts)-lag_window_source_target)]) + (1-linear_cor)*rnorm(length(ts))
  } else if (source_target_function == "linear_autocorrelated"){
    if (zscore==T){
      ts <- (ts-mean(ts))/sd(ts)
    }
    target_1 = rnorm(1)
    target = c(target_1)
    for (i in seq(1,ts_length-1)){
      if (i < lag_window_source_target) {
        target_t = ac_target_coef*target[length(target)]
      } else
        target_t = ac_target_coef*target[length(target)] + linear_cor*ts[i-lag_window_source_target+1] + (1-linear_cor)*rnorm(1)
      target = c(target, target_t)
    }
  }
  return(target)
}

generate_features_non_overlapping <- function(ts,timescale){
  features <- data.frame(matrix(ncol = 24, nrow = 0))
  segment_start = length(ts)%%timescale + 1
  while (segment_start < length(ts)){
    segment_features <- catch22_all(ts[segment_start:(segment_start+timescale-1)], catch24=T)
    segment_features <- pivot_wider(segment_features,names_from = names,values_from = values)
    segment_features <- data.frame(segment_features)
    features <- rbind(features,segment_features)
    segment_start <- segment_start + timescale
  }
  return(features[,feature_names])
}

generate_catch22_features_sliding_window <- function(ts, timescale, step_size =1 ){
  features <- data.frame(matrix(ncol = 24, nrow = 0))
  segment_start = 1
  while (segment_start <= length(ts) - timescale + 1){
    segment_features <- catch22_all(ts[segment_start:(segment_start+timescale-1)], catch24=T)
    segment_features <- pivot_wider(segment_features,names_from = names,values_from = values)
    segment_features <- data.frame(segment_features)
    features <- rbind(features,segment_features)
    segment_start <- segment_start + step_size 
  }
  return(features)
}

generate_feature_sliding_window <- function(ts, timescale, step_size =1, feature_name ){
  segment_start <- 1
  feature_ts <- c()
  while (segment_start <= length(ts) - timescale + 1){
    if (feature_name == "spike_count"){
      segment_feature <- count_spike_by_fixed_threshold(ts[segment_start:(segment_start+timescale-1)],4)
    }
    else if (feature_name == "dominant_frequency"){
      segment_feature <- dominant_frequency(ts[segment_start:(segment_start+timescale-1)])
    }
    else if (feature_name == "skewness"){
      segment_feature <- skewness(ts[segment_start:(segment_start+timescale-1)])
    }
    else if (feature_name == "autocorrelation_lag1"){
      segment_feature <- acf(ts[segment_start:(segment_start+timescale-1)],lag.max = 1,plot=F)$acf[2]
    }
    feature_ts <- c(feature_ts,segment_feature)
    segment_start <- segment_start + step_size
  } 
  return(feature_ts)
}

dominant_frequency <- function(ts){
  ts_fft <- fft(ts)
  N <- length(ts_fft)/2
  ts_fft_mod <- Mod(ts_fft[1:N])
  dominant_frequency_index <- which(ts_fft_mod==max(ts_fft_mod))[1]
  dominant_frequency <- (dominant_frequency_index-1)/length(ts)
  return(dominant_frequency)
}

find_theiler_window <- function(ts){
  autocors <- acf(ts,lag.max=length(ts),plot=F)
  theiler_window <- min(which(autocors$acf<0),floor(length(ts)/10)) - 2
  return(theiler_window)
}

compute_TE <- function(theiler_window, source, target, k, k_tau, l, l_tau, delay, reinit=T,number_surrogates=50,auto_embedding=F, theiler_window_multiplier=1, rotating_surrogates=T) {
  if(reinit) {
    .jcall(teCalc, "V", "initialise",as.integer(k), as.integer(k_tau), as.integer(l), as.integer(l_tau), as.integer(delay))
  }
  if(auto_embedding){
    .jcall(teCalc,"V","setProperty","PROP_AUTO_EMBED_METHOD","AUTO_EMBED_METHOD_MAX_CORR_AIS_AND_TE")
  } 
  .jcall(teCalc, "V", "setProperty", "DYN_CORR_EXCL", as.character(theiler_window_multiplier*theiler_window))
  .jcall(teCalc, "V", "setObservations", source, target)
  TE <- .jcall(teCalc, "D", "computeAverageLocalOfObservations")
  if (rotating_surrogates){
    numberObservations <- .jcall(teCalc,"I","getNumObservations")
    newOrderings <- matrix(0, nrow = number_surrogates, ncol = numberObservations)
    for (s in 1:number_surrogates) {
      # Need to shift more than the window, but not so much that we push values back within their own window
      thisShift <- theiler_window + sample(numberObservations - 2*theiler_window + 1, 1)
      newOrderings[s, ] <- shifter(c(1:numberObservations),thisShift) -1 
    }
    storage.mode(newOrderings) <- "integer"
    newOrderingsArray <- .jarray(newOrderings,"[I",dispatch=TRUE)
    nullDist <- .jcall(teCalc, "Linfodynamics/utils/EmpiricalMeasurementDistribution;", "computeSignificance",newOrderingsArray)
    
  } else{
    nullDist <- .jcall(teCalc, "Linfodynamics/utils/EmpiricalMeasurementDistribution;", "computeSignificance", as.integer(number_surrogates))
  }
  pValue <- nullDist$pValue
  return(data.frame(TE=TE, pValue=pValue))
}

detect_dependency_with_catch22 <- function(source,target,feature_window,target_hist, number_surrogates=1000,theiler_window_multiplier=1,rotating_surrogates=T){
  source_features <- generate_catch22_features_sliding_window(source,feature_window)
  # MEASURE DEPENDENCIES
  # TE  between target and features
  features_results <- data.frame(Feature = character(),
                           TE = numeric(),
                           pValue = numeric(),
                           stringsAsFactors = FALSE)
  
  for (feature in names(source_features)){
    feature_vector <- source_features[[feature]]
    theiler_window_fte <- max(find_theiler_window(feature_vector),find_theiler_window(target))
    TE_result <- compute_TE(theiler_window=theiler_window_fte,
                            source=feature_vector,
                            target=tail(target,length(feature_vector)),
                            k=target_hist,
                            k_tau=1,
                            l=1,
                            l_tau=1,
                            delay=1,
                            number_surrogates=number_surrogates,
                            theiler_window_multiplier=theiler_window_multiplier,
                            rotating_surrogates=rotating_surrogates)
    features_results <-rbind.data.frame(features_results,data.frame(Feature=feature, TE=TE_result$TE,pValue=TE_result$pValue,stringsAsFactors = F))
  }
  return(features_results)
}

shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}




