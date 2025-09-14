######################### Read and process data ##########################
##########################################################################
feature_name_mapping <- read.csv(here("Feature_name_mapping.csv"),header=T)
# Define color and line type mappings
color_map <- c("null" = "#808080", "MI_s" = "black",
               "mean"="#56B4E9","std"="#0B5394",
               "mode_5"="#D55E00","mode_10"="#9e4600",
               "trev"="#58D6A7",
               "acf_timescale"="#E69F00","periodicity"="#a2a848","centroid_freq"="#FFAA00",
               "autocorrelation_lag1"="#DC97BD","spike_count"="#CC79A7","dominant_frequency"="#E261A9",
               "MI_F"="#008a7d")

# process all results
simulation_results <- read.csv(paste0(here("results","simulation_studies"),"/featureBasedDependency_simulation_results.csv"),header=T)

simulation_results$seed <- as.integer(simulation_results$seed)
simulation_results$ts_length <- as.integer(simulation_results$ts_length)
simulation_results$driving_feature_timescale <- as.integer(simulation_results$driving_feature_timescale)
simulation_results$capturing_feature_timescale <- as.integer(simulation_results$capturing_feature_timescale)
simulation_results$linear_cor <- as.numeric(simulation_results$linear_cor)
simulation_results$TE <- as.numeric(simulation_results$TE)
simulation_results$pValue <- as.numeric(simulation_results$pValue)

simulation_results$driving_feature <- ifelse(simulation_results$driving_feature %in% feature_name_mapping$Feature,
                                                feature_name_mapping$Short_name[match(simulation_results$driving_feature, feature_name_mapping$Feature)],
                                                simulation_results$driving_feature)

simulation_results$capturing_feature <- ifelse(simulation_results$capturing_feature %in% feature_name_mapping$Feature,
                                                  feature_name_mapping$Short_name[match(simulation_results$capturing_feature, feature_name_mapping$Feature)],
                                                  simulation_results$capturing_feature)

# apply Holm pvalue correction 
simulation_results <- simulation_results %>%
  group_by(seed, process, ts_length, driving_feature, driving_feature_timescale,linear_cor) %>%
  mutate(
    pValue_holm_corrected = ifelse(
      inference_type == "MI_f",
      p.adjust(pValue, method = "holm"),
      pValue
    )
  ) 
simulation_results$significant_pValue_holm <- ifelse(simulation_results$pValue_holm_corrected<0.05,1,0)

bonferroni_p_value_threshold <- 0.05/25 #apply Bonferroni correction for p-value when doing multiple hypothesis testing 
simulation_results$significant_pValue_bonferroni <- ifelse(simulation_results$inference_type%in%c("null","MI_s") & simulation_results$pValue < 0.05,1,
                                                 ifelse(simulation_results$inference_type%in%c("null","MI_s") & simulation_results$pValue >= 0.05,0,
                                                        ifelse(simulation_results$pValue < bonferroni_p_value_threshold,1,0)))



################### Random noise process results #####################
######################################################################
random_results <- simulation_results[simulation_results$process=="random",]

# Matrix plots for noise process to show similar features can capture same thing
random_results_matrix <- random_results %>%
  filter(capturing_feature != "null" & ts_length==1000) %>% 
  group_by(driving_feature, capturing_feature) %>%
  summarise(avg_capture_rate=mean(significant_pValue_holm)*100)

# Order features by feature similarity (measured by the correlation of feature values across time series for each pair of features)
file_path <- paste0(here("ordered_feature_names.txt")) # this file is the result of running "order_features_by_similarity.R"
ordered_feature_names <- read.table(file = file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
feature_name_order <- ordered_feature_names$V1 # Extract the feature names as a character vector
short_feature_order <- feature_name_mapping$Short_name[match(feature_name_order, feature_name_mapping$Feature)]

# remove features with 0 standard deviation
features_to_remove <- c("dfa","rs_range","low_freq_power","entropy_pairs","ami_timescale","periodicity")

random_results_matrix <-  random_results_matrix[!random_results_matrix$capturing_feature%in%c(features_to_remove,"MI_s")
                                                & !random_results_matrix$driving_feature%in%features_to_remove,]

random_results_matrix$driving_feature <- factor(random_results_matrix$driving_feature, levels = short_feature_order)
random_results_matrix$capturing_feature <- factor(random_results_matrix$capturing_feature, levels = rev(c(short_feature_order, "MI_s")))

# Capture rate matrix (Figure 4a)
random_results_matrix_plot_path = file.path(output_dir, "RandomResultsMatrix.pdf")
random_results_matrix_plot <- ggplot(random_results_matrix,aes(x = driving_feature, y = capturing_feature)) +
  #make heatmap with geom_tile 
  geom_tile(aes(fill = avg_capture_rate)) +
  #add the total number of devices to each tile
  geom_fit_text(aes(label = round(avg_capture_rate)), size = 8) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  theme(text = element_text(size = 12))  +
  #some formatting to make things easier to read
  scale_x_discrete(position = "top") +
  scale_fill_gradient(high = "#009688", low = "#DDF1EF") +
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank()) +
  labs(x="Driving feature", y = "Capturing feature")

ggsave(filename = random_results_matrix_plot_path, plot = random_results_matrix_plot) 

# Plots of capture rates for different capturing feature vs signal one across different TS lengths
random_results_summary <- random_results %>% 
  group_by(process,ts_length,driving_feature,driving_feature_timescale,capturing_feature, capturing_feature_timescale,linear_cor) %>%
  summarise(capture_rate_holm = mean(significant_pValue_holm)*100) %>% 
  data.frame()

random_MIF_results <- random_results %>%
  filter(!capturing_feature %in% c("null","MI_s")) %>%
  group_by(process, ts_length, driving_feature, driving_feature_timescale, capturing_feature_timescale, linear_cor, seed) %>%
  summarise(MIF_detected = max(significant_pValue_holm)) %>%
  group_by(process, ts_length, driving_feature, driving_feature_timescale, capturing_feature_timescale, linear_cor) %>%
  summarise(capture_rate_holm = mean(MIF_detected) * 100) %>%
  mutate(capturing_feature = "MI_F")
random_MIF_results <- random_MIF_results[,names(random_results_summary)]
random_results_summary <- rbind.data.frame(random_results_summary,random_MIF_results)

selected_features <- c("mean",
                       "acf_timescale",
                       "mode_5",
                       "trev",
                       "MI_F")
random_results_summary$ts_length_log <- log10(random_results_summary$ts_length)

random_results_summary$line_type <- ifelse(random_results_summary$capturing_feature == "null", "dashed", "solid")
random_results_summary$line_thickness <- ifelse(random_results_summary$capturing_feature ==random_results_summary$driving_feature, "thick", "normal")

random_results_line_plot_path = file.path(output_dir, "RandomResultsLinePlots.pdf")
random_results_line_plot <- ggplot(data = random_results_summary[random_results_summary$driving_feature %in% selected_features
                                     & random_results_summary$capturing_feature %in% c(selected_features, "null", "MI_s"),], 
       aes(ts_length, capture_rate_holm, col = capturing_feature, 
           linetype=line_type,
           size = line_thickness,
           shape=capturing_feature,
           order = ifelse(capturing_feature == "MI_F", 2, 1))) +
  geom_line() + 
  geom_point(size = 2,position = position_jitter(width = 0.04)) +
  #scale_x_log10(limits = c(min(random_results_summary$ts_length_log), max(random_results_summary$ts_length_log))) +
  scale_x_log10(breaks = 10^(0:4),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(side="b") +
  labs(x = "Time series length, T (samples)", y = "Capture rate (%)", col = "Capturing feature") +
  facet_wrap(~ driving_feature, scales = "free_x", ncol = 1) +
  theme_classic() +
  theme(strip.background = element_blank(),  # Remove the box around the title
        strip.text = element_text(size = 10),  # Adjust the title text size
        axis.text.x = element_text(hjust = 1),
        panel.grid = element_blank()) +  
  scale_color_manual(values = color_map[names(color_map)%in%c(selected_features,"null","MI_s")]) +
  scale_linetype_manual(values = c("dashed"="dashed","solid"="solid")) +
  scale_size_manual(values = c("thick" = 1.5, "normal" = 0.75)) +
  guides(col = guide_legend(override.aes = list(size = 1.5))) +
  theme(legend.key.size = unit(5, "mm"))  # Adjust legend key size for better visibility
ggsave(filename = random_results_line_plot_path, plot = random_results_line_plot, width = 5.5, height = 7)  # adjust height if needed


###################### Non-stationary processes' results #####################
##############################################################################
# function for summarising capture rate by process
process_results_summary <- function(data, summary_type,ts_length_filter, driving_feature_timescale_filter, capturing_feature_timescale_filter, linear_cor_filter) {
  
  if (summary_type == "noise"){
    filtered_data <- data %>%  
      filter(ts_length == ts_length_filter, driving_feature_timescale == driving_feature_timescale_filter, capturing_feature_timescale == capturing_feature_timescale_filter)
    result_summary <- filtered_data %>%
      group_by(process, ts_length, driving_feature, capturing_feature, noise_coef) %>%
      summarise(avg_capture_rate_holm = mean(significant_pValue_holm) * 100) %>%
      data.frame()
    
    MIF_results <- filtered_data %>%
      filter(!capturing_feature %in% c("null","MI_s","spike_count","autocorrelation_lag1","dominant_frequency")) %>%
      group_by(process, ts_length, driving_feature, noise_coef, seed) %>%
      summarise(fipi_detected = max(significant_pValue_holm)) %>%
      group_by(process, ts_length, driving_feature, noise_coef) %>%
      summarise(avg_capture_rate_holm = mean(fipi_detected) * 100) %>%
      mutate(capturing_feature = "MI_F")
    
    result_summary <- rbind.data.frame(MIF_results, result_summary)
    return(result_summary)
  }
  
  if (summary_type == "driving_timescale"){
    filtered_data <- data %>%  
      filter(ts_length == ts_length_filter, linear_cor == linear_cor_filter, capturing_feature_timescale == driving_feature_timescale)
    result_summary <- filtered_data %>%
      group_by(process, ts_length, driving_feature, capturing_feature, driving_feature_timescale) %>%
      summarise(avg_capture_rate_holm = mean(significant_pValue_holm) * 100) %>%
      data.frame()
    
    MIF_results <- filtered_data %>%
      filter(!capturing_feature %in% c("null","MI_s","spike_count","autocorrelation_lag1","dominant_frequency")) %>%
      group_by(process, ts_length, driving_feature, driving_feature_timescale, seed) %>%
      summarise(fipi_detected = max(significant_pValue_holm)) %>%
      group_by(process, ts_length, driving_feature, driving_feature_timescale) %>%
      summarise(avg_capture_rate_holm = mean(fipi_detected) * 100) %>%
      mutate(capturing_feature = "MI_F")
    
    result_summary <- rbind.data.frame(MIF_results, result_summary)
    return(result_summary)
  }
  
  if (summary_type == "capturing_timescale"){
    filtered_data <- data %>%  
      filter(ts_length == ts_length_filter, linear_cor == linear_cor_filter, driving_feature_timescale == driving_feature_timescale_filter)
    result_summary <- filtered_data %>%
      group_by(process, ts_length, driving_feature, capturing_feature, capturing_feature_timescale) %>%
      summarise(avg_capture_rate_holm = mean(significant_pValue_holm) * 100) %>%
      data.frame()
    
    MIF_results <- filtered_data %>%
      filter(!capturing_feature %in% c("null","MI_s","spike_count","autocorrelation_lag1","dominant_frequency")) %>%
      group_by(process, ts_length, driving_feature, capturing_feature_timescale, seed) %>%
      summarise(fipi_detected = max(significant_pValue_holm)) %>%
      group_by(process, ts_length, driving_feature, capturing_feature_timescale) %>%
      summarise(avg_capture_rate_holm = mean(fipi_detected) * 100) %>%
      mutate(capturing_feature = "MI_F")
    
    result_summary <- rbind.data.frame(MIF_results, result_summary)
    return(result_summary)
  }
  
}

plot_process <- function(process, selected_features, output_dir) {
  # Subset and compute extra columns
  process_results <- simulation_results[simulation_results$process==process,]
  process_results$SNR <- process_results$linear_cor^2 / (1 - process_results$linear_cor)^2
  process_results$noise_coef <- 1 - process_results$linear_cor
  
  # Capture rate summaries by different factors - noise, driving timescale, capturing timescale
  summaries <- list(
    noise = list(summary_type = "noise",
                 filters = list(driving_feature_timescale_filter = 100,
                                capturing_feature_timescale_filter = 100,
                                linear_cor_filter = NA,
                                ts_length_filter = 1000),
                 xvar = "noise_coef",
                 xlabel = "Noise coefficient",
                 filename = paste0(process, "_Noise.pdf")),
    
    driving_timescale = list(summary_type = "driving_timescale",
                             filters = list(driving_feature_timescale_filter = NA,
                                            capturing_feature_timescale_filter = NA,
                                            linear_cor_filter = 0.4,
                                            ts_length_filter = 1000),
                             xvar = "driving_feature_timescale",
                             xlabel = "Interaction timescale",
                             filename = paste0(process, "_DrivingTimescale.pdf")),
    
    capturing_timescale = list(summary_type = "capturing_timescale",
                               filters = list(driving_feature_timescale_filter = 100,
                                              capturing_feature_timescale_filter = NA,
                                              linear_cor_filter = 0.4,
                                              ts_length_filter = 1000),
                               xvar = "capturing_feature_timescale",
                               xlabel = "Capturing feature timescale",
                               filename = paste0(process, "_CapturingTimescale.pdf"))
  )
  
  # Loop through each summary type
  for (s in summaries) {
    res <- do.call(process_results_summary,
                   c(list(process_results, s$summary_type), s$filters))
    
    p <- ggplot(res[res$capturing_feature %in% selected_features,],
                aes_string(s$xvar, "avg_capture_rate_holm",
                           color = "capturing_feature")) +
      geom_line(size = 1.25) +
      geom_point(size = 2) +
      xlab(s$xlabel) + ylab("Capture rate (%)") +
      scale_color_manual(values = color_map[names(color_map) %in% selected_features]) +
      theme_classic() +
      theme(legend.text = element_text(size = 6),
            strip.background = element_blank(),
            strip.text = element_text(size = 20),
            axis.text.x = element_text(hjust = 1, size = 13.2),
            axis.text.y = element_text(size = 13.2),
            panel.grid = element_blank())
    
    ggsave(filename = file.path(output_dir, s$filename),
           plot = p, width = 6, height = 4)
  }
}

# AR3
plot_process(
  process = "AR3",
  selected_features = c("null","MI_s","std","mean", "autocorrelation_lag1","acf_timescale","MI_F"),
  output_dir = output_dir
)

# Bimodal spiking
plot_process(
  process = "bimodal_spiking",
  selected_features = c("null","MI_s","std","mean", "mode_10","spike_count","MI_F"),
  output_dir = output_dir
)

