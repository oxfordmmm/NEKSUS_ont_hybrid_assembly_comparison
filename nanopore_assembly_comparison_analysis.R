#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Nanopore Assembly Comparison - R analysis script
# R Analysis script to produce summary statistics, significance tests, and plots for Nanopore assembly comparison manuscript
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Set working directory
setwd("C:/Users/dnagy/OneDrive - Nexus365/Documents/DPhil_Clin_Medicine/DPhil/NEKSUS/pilot2/for_publication")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Install & Load Packages ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# install if needed
if (!requireNamespace("tidyverse", quietly = TRUE)) {
 install.packages("tidyverse")
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
 install.packages("gridExtra")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
 install.packages("reshape2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
 install.packages("ggrepel")
}
if (!requireNamespace("ggsignif", quietly = TRUE)) {
 install.packages("ggsignif")
}
if (!requireNamespace("gtsummary", quietly = TRUE)) {
 install.packages("gtsummary")
}
if (!requireNamespace("rstatix", quietly = TRUE)) {
 install.packages("rstatix")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
 install.packages("patchwork")
}
if (!requireNamespace("ggforce", quietly = TRUE)) {
 install.packages("ggforce")
}
if (!requireNamespace("colourspace", quietly = TRUE)) {
 install.packages("colourspace")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
 install.packages("igraph")
}
if (!requireNamespace("grid", quietly = TRUE)) {
 install.packages("grid")
}


#library(tidyverse)
# or load only required packages from tidyverse
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(purrr)

#Load other packages
library(gridExtra)
library(reshape2)
library(ggrepel)
library(ggsignif)
library(gtsummary)
library(rstatix)
library(patchwork)
library(ggforce)  
library(colorspace)
library(igraph)
library(grid)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Raw Reads ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Load data ####
raw_qc_sup <- read.csv("raw_qc_sup_cleaned.csv")
#Inspect data if required
#View(raw_qc_sup)
#colnames(raw_qc_sup)
#str(raw_qc_sup)

#~~~~~~~~~~~~~~~~~~~~~#
# * Summary Stats ####
# Define summary sats function by sample type, mean, median, mode, min, max, IQR, Q1, Q3, etc.
summary_stats <- function(dataframe, grouping_var, var_to_summarise) {
  # Group by the specified grouping variable
  grouped_data <- aggregate(dataframe[[var_to_summarise]], by = list(dataframe[[grouping_var]]), FUN = function(x) {
    c(mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      #mode = as.numeric(names(sort(table(x), decreasing = TRUE)[1])),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE),
      Q1 = quantile(x, probs = 0.25, na.rm = TRUE),
      Q3 = quantile(x, probs = 0.75, na.rm = TRUE),
      iqr = IQR(x, na.rm = TRUE))
  })
  return(grouped_data)
}

#Make summary table with mean, median, mode, Q1, Q3, min, max, range and IQR of all metrics
raw_qc_summary_table <- data.frame()
for (col_index in 6:23) {
  col_name <- colnames(raw_qc_sup)[col_index]
  summary_table <- summary_stats(raw_qc_sup, "assembler", col_name)
  summary_table <- cbind(metric = col_name, summary_table)
  raw_qc_summary_table <- rbind(raw_qc_summary_table, summary_table)
}
#View(raw_qc_summary_table)
write.table(raw_qc_summary_table, file = "raw_qc_summary_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Another summary table considering only first subsampled read set:
# Filter for only the un-subsampled and first subsampled read sets
filtered_qc <- raw_qc_sup |>
  filter(assembler %in% c("Illumina (raw)", "ONT sup (raw)", "Illumina (sub)", "ONT sup (sub)"))

# Summarise key stats for each assembler and metric
qc_summary <- filtered_qc |>
  pivot_longer(cols = c(coverage, avg_len, AvgQual), names_to = "metric", values_to = "value") |>
  group_by(assembler, metric) |>
  summarise(
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(summary = sprintf("%.1f (%.1f–%.1f)", median, Q1, Q3))

qc_summary <- data.frame(qc_summary)



# Reshape to wide format with metrics as rows and raw/sub in columns
qc_summary_table <- qc_summary |>
  mutate(assembler_type = case_when(
    str_detect(assembler, "Illumina") ~ "Illumina",
    str_detect(assembler, "ONT") ~ "ONT"
  ),
  version = case_when(
    str_detect(assembler, "\\(raw\\)") ~ "Raw",
    str_detect(assembler, "\\(sub\\)") ~ "Sub"
  )) |>
  select(metric, assembler_type, version, summary) |>
  pivot_wider(names_from = version, values_from = summary) |>
  arrange(factor(metric, levels = c("coverage", "avg_len", "AvgQual")),
          factor(assembler_type, levels = c("Illumina", "ONT")))

print(qc_summary_table)
#save
write.csv(qc_summary_table, file = "raw_qc_summary_table_gtsummary.csv")


#~~~~~~~~~~~~~~~~~~#
# * Stats Test ####
#Wilcoxon signed rank test for coverage, read length (both mean (avg_len) and median (Q2)), and AvgQual between un-subsampled (raw) and subsampled (sup).
#arrange by sample to get constant order for pairwise comparison
raw_qc_sup  <- raw_qc_sup |> arrange(sample)
filtered_qc  <- filtered_qc |> arrange(sample)

#Seaprate into illumina and ONT sets for stats test
raw_qc_sup_illumina <- raw_qc_sup |> filter(grepl("Illumina", assembler))
raw_qc_sup_ont <- raw_qc_sup |> filter(grepl("ONT", assembler))

metrics <- c("coverage", "avg_len", "Q2", "AvgQual")

# Store results
wilcox_results <- list()

# Pairwise wilcoxon signed rank test to chec for no difference between un-subsampled and subsampled (and different subsamples)
for (metric in metrics) {
  message("Running Wilcoxon tests for metric: ", metric)
  
  illumina_result <- tryCatch({
    pairwise.wilcox.test(raw_qc_sup_illumina[[metric]], raw_qc_sup_illumina$assembler,
      p.adjust.method = "bonf", paired = TRUE, exact = FALSE)
  }, error = function(e) e)
  
  ont_result <- tryCatch({
    pairwise.wilcox.test(raw_qc_sup_ont[[metric]], raw_qc_sup_ont$assembler,
      p.adjust.method = "bonf", paired = TRUE,  exact = FALSE )
  }, error = function(e) e)
  
  # Save results
  wilcox_results[[metric]] <- list(
    illumina = illumina_result,
    ont = ont_result
  )
}

# Print all results cleanly
for (metric in metrics) {
  cat("\n=== Metric:", metric, "===\n")
  cat("\n-- Illumina --\n")
  print(wilcox_results[[metric]]$illumina)
  
  cat("\n-- ONT --\n")
  print(wilcox_results[[metric]]$ont)
}


#~~~~~~~~~~~~~~#
# * Plot ####
custom_colours <- c(
  "ONT hac (raw)" = "#006666",
  "ONT hac (sub)" = "#66b2b2",
  "ONT sup (raw)" = "#36439A",
  "ONT sup (sub)" = "#4a7bb7",
  "Illumina (raw)" = "#802255",
  "Illumina (sub)" = "#CC6677")

#Write function to plot QC metrics:
plot_qc_metrics <- function(data, columns, ncol = 3, colours) {
  plot_list <- list()
  
  for (col_name in columns) {
    plot <- ggplot(data, aes(x = assembler, y = .data[[col_name]], fill = assembler)) +
      geom_boxplot(width = 0.3, color = "black") +  # Thin black outline
      geom_point(aes(fill = assembler), color = "black", size = 1, stroke = 0.3,  
                 shape = 21, alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +  
      labs(title = col_name, x = "read type") +
      scale_fill_manual(values = colours, drop = FALSE) +  
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5, size = 8)
      ) 
    
    plot_list[[col_name]] <- plot
  }
  
  arranged_plots <- do.call(gridExtra::grid.arrange, c(plot_list, ncol = ncol))
  return(arranged_plots)
}

# change column names to be more complete
colnames(filtered_qc) <- c("sample","assembler","file","format","type",
                          "No. sequences","Total No. Bases","Minimum Read Length",
                          "Mean Read Length", "Maximum Read Length",
                          "First Quartile (Q1) Read Length",  "Median (Q2) Read Length", "Third Quartile (Q3) Read Length", 
                          "No. Gaps", "N50" ,"N50_num", "% Bases >Q20 Quality Score","% Bases >Q30 Quality Score", 
                          "Mean Quality Score", "% GC Content", "No. Ambiguous Bases", "Genome Size", "Coverage"   )
raw_qc_arranged_plots <- plot_qc_metrics(data = filtered_qc, columns = colnames(filtered_qc)[6:23], 
                                         ncol = 3, colours = custom_colours)
raw_qc_arranged_plots
ggsave("for_publication/raw_qc_plots.png", raw_qc_arranged_plots, width = 12, height = 15, units = "in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Assemblies- Chromosomes (+ all contigs) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Load data ####
contigs_sup <- read.csv("contigs_summary_sup_cleaned.csv")
#View(contigs_sup)

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Summary Stats ####

# Condense to 1 row per sample/assembler combo and add chromosome circularity and completenes 
## complete = >4Mb; circularity info derived from fasta headers or from gfa files where available
condense_contigs <- function(df){
  df$sample <- sub("_failed$", "", df$sample)
  contigs_per_sample <- df |>
    group_by(sample, assembler) |>
    filter(length == max(length))|>
    ungroup()
  # Add chromosome circularity and completeness info
  contigs_per_sample <- contigs_per_sample |>
    group_by(sample, assembler) |>
    mutate(complete_circular_chromosomes = sum(length >= 4000000 & circular == "true"),
           complete_noncircular_chromosomes = sum(length >= 4000000 & circular != "true"),
           incomplete_noncircular_chromosomes = sum(length <4000000 & circular != "true"),
           incomplete_circular_chromosomes = sum(length <4000000 & circular == "true"))
  return(contigs_per_sample)
}

# Check size distribution of circular contigs
#shortest circular choromosome = 4.6Mbp, longers circular plasmid ~300 kbp, so big gap, therefore reasonable to separate chsomosomes/ plasmids based on length for this data
circular_contigs <- contigs_sup |>
  filter(circular == "true")
hist(circular_contigs$length)
#View(circular_contigs)


#call function to condense contig dataframes to 1 row per sample/assembler combo
contigs_per_sample_sup <- condense_contigs(contigs_sup)


#Summarise per assembler:
summarise_by_assembler <- function(df) {
  chromosomes_summary_by_assembler <- df |>
    group_by(assembler) |>
    summarise(
      complete_circular_chromosome_count = sum(complete_circular_chromosomes ==1, na.rm = TRUE),
      complete_noncircular_chromosome_count = sum(complete_noncircular_chromosomes ==1, na.rm = TRUE),
      all_complete_chromosome_count = sum(completeness == "complete", na.rm = TRUE),
      incomplete_noncircular_chromosome_count = sum(incomplete_noncircular_chromosomes ==1, na.rm = TRUE),
      incomplete_circular_chromosome_count = sum(incomplete_circular_chromosomes ==1, na.rm = TRUE))
  return(chromosomes_summary_by_assembler)
}

# Call function to summarise chromosomes
chromosomes_summary_sup <- summarise_by_assembler(contigs_per_sample_sup)


#~~~~~~~~~~~~~~~~~~~~#
# * Stats Test ####

#set assembler order
assembler_priority <- c("autocycler", "flye", "hybracter_long", "hybracter_hybrid", "unicycler", "unicycler_bold")
contigs_per_sample_sup$assembler <- factor(contigs_per_sample_sup$assembler, levels = assembler_priority)

#2-proportion test to calculate Chi-squared test statistic. Care if failure count <5
# Create a contingency table
contingency_table_chromosomes <- contigs_per_sample_sup |>
  group_by(assembler) |>
  summarise(count = sum(complete_circular_chromosomes ==1, na.rm = TRUE)) |>
  mutate(failure = 92-count) |>
  select(count, failure)
contingency_table_chromosomes <- as.matrix(contingency_table_chromosomes)
rownames(contingency_table_chromosomes) <- assembler_priority
  
# Perform the chi-squared test (base R stats)
chisq_test <- chisq.test(contingency_table_chromosomes)
print(chisq_test)

#Pairwise Fisher's exact test
pairwise_fisher_results <- rstatix::pairwise_fisher_test(contingency_table_chromosomes)
print(pairwise_fisher_results)

#~~~~~~~~~~~~~~~~~~~~~~~#
# * Plot ####

# Order samples based on no_circular_chromosomes for each assembler in priority order
sample_order <- contigs_per_sample_sup |>
  filter(assembler %in% assembler_priority) |> 
  select(sample, assembler, no_circular_chromosomes) |>
  pivot_wider(names_from = assembler, values_from = no_circular_chromosomes) |>
  arrange(desc(autocycler), desc(flye), desc(hybracter_long),  desc(hybracter_hybrid), desc(unicycler), desc(unicycler_bold)) |>
  pull(sample)

# Apply order to main df
contigs_per_sample_sup <- contigs_per_sample_sup |>
  mutate(sample = factor(sample, levels = sample_order)) |>
  mutate(sample_id = as.integer(sample)) # replce sample IDs with numbers

# Define the color mapping for each assembler when no_circular_chromosomes == 1
assembler_colors <- c(
  "autocycler.1" = "#332288",  # Dark Blue for 1 in autocycler
  "flye.1" = "#88ccee",        # Light Blue for 1 in flye
  "hybracter_long.1" = "#44aa99", # Teal for 1 in hybracter_long
  "hybracter_hybrid.1" = "#117733", # Green for 1 in hybracter_hybrid
  "unicycler.1" = "#cc6677",   # Red for 1 in unicycler
  "unicycler_bold.1" = "#882255",  # Dark Red for 1 in unicycler_bold
  "autocycler.0" = "#ECEADA",  # Light cream for 0 in autocycler
  "flye.0" = "#ECEADA",        # Light cream for 0 in flye
  "hybracter_long.0" = "#ECEADA", # Light cream for 0 in hybracter_long
  "hybracter_hybrid.0" = "#ECEADA", # Light cream for 0 in hybracter_hybrid
  "unicycler.0" = "#ECEADA",   # Light cream for 0 in unicycler
  "unicycler_bold.0" = "#ECEADA"  # Light cream for 0 in unicycler_bold
)

# Custom assembler labels
assembler_labels <- c(
  "flye" = "Flye",
  "autocycler" = "Autocycler",
  "unicycler" = "Unicycler",
  "unicycler_bold" = "Unicycler \n(bold)",
  "hybracter_long" = "Hybracter\n(long)", 
  "hybracter_hybrid" = "Hybracter\n(hybrid)"
)

percent_samples_fully_circularised <- contigs_per_sample_sup |>
  group_by(assembler) |>
  summarise(percent = round(sum(no_circular_chromosomes == 1) / n() * 100,0),
            count = sum(no_circular_chromosomes == 1),  # % of samples with circular chromosome
            total = n(),
            failure = total - count)

# Tile plot to show circularised vs non-circularised chromosomes
chromosomes_plot_sup <- ggplot(contigs_per_sample_sup, aes(x = assembler, y = sample_id)) +
  # Set the tile color based on the interaction between assembler and no_circular_chromosomes
  geom_tile(aes(fill = interaction(assembler, no_circular_chromosomes)), color = "white", linewidth = 0.3) +
  
  # Set the colors for each assembler (for no_circular_chromosomes == 1) and 0 as light cream
  scale_fill_manual(
    values = assembler_colors,
    na.value = "#ECEADA",  # Ensure 0 gets light cream color
    name = "No. Circular Chromosomes"
  ) +
  
  scale_x_discrete(labels = assembler_labels) +  # Custom labels
  scale_y_continuous(
    breaks = contigs_per_sample_sup$sample_id,
    labels = contigs_per_sample_sup$sample_id
  ) +
  labs(x = "Assembler", y = "Sample", title = "No. Circularised Chromosomes by Assembler") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold", angle = 0, hjust = 0.5, margin = margin(t=15)),  # Larger x-axis font
        axis.text.y = element_text(size = 5),  # Smaller y-axis font
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=20, t = 20)),  # Centered title
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10), # Increase top margin
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, face = "bold")) +
  coord_cartesian(clip = "off") + 
  # Add assembler-specific percentage annotations
  annotate("text", x = percent_samples_fully_circularised$assembler, y = length(sample_order) + 3, 
           label = paste0(round(percent_samples_fully_circularised$percent, 1), "%"), size = 7, fontface = "bold") +
  
  # Add brackets for assembler groups using annotate()
  annotate("segment", x = 0.8, xend = 3.2, y = -3, yend = -3, size = 0.8) +  # Bracket line for long-read assemblers
  annotate("segment", x = 3.8, xend = 6.2, y = -3, yend = -3, size = 0.8) +  # Bracket line for hybrid assemblers
  annotate("text", x = 2, y = -5.8, label = "Long-read only", size = 5) +  
  annotate("text", x = 5.0, y = -5.8, label = "Hybrid", size = 5)

chromosomes_plot_sup
ggsave("for_publication/chromosomes_plot_sup.png", chromosomes_plot_sup, width = 6, height = 11, units = "in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3a. Plasmids- Hybracter (hybrid) reference ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Load data ####
plasmids_match_hybracter_mash <- read.csv("plasmids_match_hybracter_mash.csv")
#View(plasmids_match_hybracter_mash)

#~~~~~~~~~~~~~~~~~~~~~~#
# * Summary Stats ####

#Calculate number of unique true plasmids
num_manual_plasmids <- length(unique(plasmids_match_hybracter_mash$plasmid_id))

# Define assembler order
assembler_order <- c("autocycler", "flye", "hybracter_long", "hybracter_hybrid", "unicycler", "unicycler_bold")

# Summary Table – Simple Status
summary_simple <- plasmids_match_hybracter_mash |>
  count(assembler_qry, status_with_mash_0.025) |>
  mutate(percent = round(n / num_manual_plasmids * 100, 1),
         n_percent = paste0(n, " (", percent, "%)")) |>
  select(-n, -percent) |>
  pivot_wider(names_from = assembler_qry, values_from = n_percent, values_fill = "0 (0%)") |>
  arrange(factor(status_with_mash_0.025, levels = c("present", "misassembled", "absent"))) |>
  rename(Status = status_with_mash_0.025)
print(summary_simple)


# Summary Table – Detailed Status
summary_detailed <- plasmids_match_hybracter_mash |>
  count(assembler_qry, status_detailed_with_mash_0.025) |>
  mutate(percent = round(n / num_manual_plasmids * 100, 1),
         n_percent = paste0(n, " (", percent, "%)")) |>
  select(-n, -percent) |>
  pivot_wider(names_from = assembler_qry, values_from = n_percent, values_fill = "0 (0%)") |>
  arrange(factor(status_detailed_with_mash_0.025, levels = c(
    "present", "circ_mismatch", "size_mismatch", "mash_mismatch",
    "circ_size_mismatch", "size_mash_mismatch", "circ_mash_mismatch", "absent"
  ))) |>
  rename(Status = status_detailed_with_mash_0.025)
print(summary_detailed)


# Optional Save as CSV
#write.csv(summary_simple, "plasmids_mash_match_hybracter_ref_simple_summary_table.csv", row.names = FALSE)
#write.csv(summary_detailed, "plasmids_mash_match_hybracter_ref_detailed_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~#
# * Stats Test ####
contingency_simple <- plasmids_match_hybracter_mash |>
  group_by(assembler_qry) |>
  summarise(present = sum(status_with_mash_0.025 == "present"),
            not_present = sum(status_with_mash_0.025 != "present"))

# Ensure expected order and numeric conversion
contingency_simple$assembler_qry <- factor(contingency_simple$assembler_qry, levels = assembler_order)
contingency_simple <- contingency_simple |> mutate(across(c(present, not_present), as.numeric))

contingency_mat <- as.matrix(contingency_simple |> select(present, not_present))
rownames(contingency_mat) <- contingency_simple$assembler_qry

# Chi-squared test
chisq_present_result <- chisq.test(contingency_mat)
print(chisq_present_result)

# Pairwise Fisher's exact test
hybracter_ref_plasmids_fisher_test_results <- rstatix::pairwise_fisher_test(contingency_mat) # function defined above
print(hybracter_ref_plasmids_fisher_test_results)


# Detailed metrics chi-squared for each
contingency_detailed <- plasmids_match_hybracter_mash |>
  group_by(assembler_qry, status_detailed_with_mash_0.025) |>
  summarise(count = n()) |>
  pivot_wider(id_cols = assembler_qry, names_from = status_detailed_with_mash_0.025, values_from = count, values_fill = 0)

# Total reference plasmids
total_plasmids <- length(unique(plasmids_match_hybracter_mash_cleaned$plasmid_id))

# Define metrics for testing
metrics <- colnames(contingency_detailed)[-1]  # excluding assembler

# Run chi-squared tests
chisq_results <- lapply(metrics, function(metric) {
  counts <- contingency_detailed[[metric]]
  others <- total_plasmids - counts
  contingency_mat <- rbind(counts, others)
  chisq.test(contingency_mat)$p.value
})

chisq_df <- data.frame(metric = metrics,p_value = unlist(chisq_results))

print(chisq_df)


#~~~~~~~~~~~~~~~~~~~~~~~~#
# * Plot ####
# Complex upset-style plot
# load extra packages if needed
#library(patchwork)
#library(ggforce)  # for geom_circle
#library(colorspace)  # for lighten

# Define assembler colors
assembler_colors <- c(
  "Autocycler" = "#332288",
  "Flye" = "#88ccee",
  "Hybracter\n(long)" = "#44aa99",
  "Hybracter\n(hybrid)" = "#117733",
  "Unicycler" = "#cc6677",
  "Unicycler\n(bold)" = "#882255"
)

status_levels <- c("absent", "misassembled", "present")
status_alpha <- rev(c(1.0, 0.5, 0.1))

# Prep long format
df_long <- plasmids_match_hybracter_mash %>%
  select(plasmid_id, assembler_label, status_with_mash_0.025) %>%
  mutate(status_with_mash_0.025 = factor(status_with_mash_0.025, levels = status_levels))

# Pivot to wide, then group by unique combinations
heatmap_data_wide <- df_long %>%
  pivot_wider(names_from = assembler_label, values_from = status_with_mash_0.025, values_fill = "absent") %>%
  group_by(across(all_of(names(assembler_colors)))) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(combo_id = row_number())  # New order: most common to least

# Long format again for heatmap
heatmap_long <- heatmap_data_wide %>%
  pivot_longer(cols = names(assembler_colors), names_to = "assembler", values_to = "status") %>%
  mutate(
    assembler = factor(assembler, levels = rev(names(assembler_colors))),  # Reverse for top-down
    combo_id = factor(combo_id, levels = unique((combo_id)))  # Keep left-to-right ordering
  )

# Get fill color with transparency
status_to_alpha <- function(status) {
  status_alpha[which(status_levels == status)]
}
status_to_base_color <- function(assembler) {
  assembler_colors[[assembler]]
}
heatmap_long <- heatmap_long %>%
  mutate(
    fill_color = mapply(function(a, s) {
      adjustcolor(status_to_base_color(as.character(a)), alpha.f = status_to_alpha(s))
    }, assembler, status)
  )



# Ensure combo_id is a factor with consistent levels
heatmap_long$combo_id <- factor(heatmap_long$combo_id, levels = unique(heatmap_long$combo_id))
heatmap_data_wide$combo_id <- factor(heatmap_data_wide$combo_id, levels = levels(heatmap_long$combo_id))


# Heatmap circles
heatmap_plot <- ggplot(heatmap_long, aes(x = combo_id, y = assembler)) +
  geom_point(aes(fill = fill_color), shape = 21, size = 5, stroke = 0.3, color = "black") +
  scale_fill_identity() +
  theme_minimal(base_size = 12) +
  labs(x = "Unique status combinations", y = "Assembler") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Top bar plot: count of plasmids per combination
top_bar <- ggplot(heatmap_data_wide, aes(x = factor(combo_id, levels = combo_id), y = count)) +
  geom_col(fill = "grey40") +
  geom_text(aes(label = count), vjust = -0.2, size = 3.2) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme_minimal(base_size = 12) +
  labs(y = "Plasmid count", x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 10, r = 5, b = 0, l = 39))

# Left stacked bar: count of each status per assembler
left_bar_data <- df_long |>
  count(assembler_label, status_with_mash_0.025) |>
  mutate(
    assembler_label = factor(assembler_label, levels = rev(names(assembler_colors))),
    
    # Reverse the levels so stacking starts with "other_misassembly" at the bottom
    status_with_mash_0.025 = factor(
      status_with_mash_0.025,
      levels = c("absent", "misassembled", "present")
    ),
    
    fill_color = mapply(function(a, s) {
      adjustcolor(assembler_colors[[as.character(a)]], alpha.f = status_to_alpha(s))
    }, assembler_label, status_with_mash_0.025)
  )

# Calculate position of label: total stacked value per assembler
label_data <- left_bar_data |>
  group_by(assembler_label) |>
  mutate(cumulative = cumsum(n)) |> 
  mutate(n_percent = paste0(round(n/cumulative *100), "%")) |>
  ungroup() |>
  filter(status_with_mash_0.025 == "present")

left_bar <- ggplot(left_bar_data, aes(x = n, y = assembler_label, fill = factor(fill_color, levels = rev(fill_color)))) +
  geom_col(position = "stack") +
  geom_text(data = label_data,
            aes(x = cumulative + 1, y = assembler_label, label = n_percent),  # Add offset to the right
            inherit.aes = FALSE,
            hjust = 0, size = 3.2) +
  scale_fill_identity() +
  theme_minimal(base_size = 12) +
  labs(x = "Plasmid count", y = NULL) +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 0)) +
  coord_cartesian(clip = "off") 


# Blank spacer
blank_spacer <- plot_spacer()

# Layout
final_plot <- (blank_spacer + top_bar) / (left_bar + heatmap_plot) +
  plot_layout(widths = c(0.25, 1), heights = c(0.3, 1))

print(final_plot)

ggsave("for_publication/plasmids_upset_plot_hybracter_ref.png", final_plot, width = 10, height = 5, units = "in")

#~~~~~~~~~~~~~~~~~~~#
#Frequency Polygon plot of contig length distributions
plot_contig_lengths <- function(df, title_lab) {
  assembler_colors <- c(
    "Autocycler" = "#332288",       # Dark Blue
    "Flye" = "#88ccee",             # Light Blue
    "Hybracter (long)" = "#44aa99",  # Teal
    "Hybracter (hybrid)" = "#117733", # Green
    "Unicycler" = "#cc6677",        # Red
    "Unicycler (bold)" = "#882255"   # Dark Red
  )
  
  assembler_labels <- c(
    "autocycler" = "Autocycler",
    "flye" = "Flye",
    "hybracter_long" = "Hybracter (long)",
    "hybracter_hybrid" = "Hybracter (hybrid)",
    "unicycler" = "Unicycler",
    "unicycler_bold" = "Unicycler (bold)"
  )
  
  # Apply labels
  df$assembler_label <- factor(df$assembler_qry, levels = names(assembler_labels), labels = assembler_labels)
  
  # Create reverse factor for plotting so desired lines are drawn on top
  df$assembler_plot <- factor(df$assembler_label, levels = rev(levels(df$assembler_label)))
  
  contig_lengths <- ggplot(df, aes(x = size_qry, colour = assembler_plot)) +
    #geom_histogram(binwidth = 0.1, alpha = 0.2, position = "identity") +
    geom_freqpoly(binwidth = 0.1, size =2, alpha = 0.9) +
    scale_color_manual(values = assembler_colors, name = "Assembler") +
    labs(x= "Length (bp)", y = "No. Plasmids", title = title_lab) +
    scale_x_log10(breaks = c(3000, 10000, 30000, 100000, 300000), labels = scales::comma) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(contig_lengths)
}

#filter for only 'present' plasmids
plasmids_match_hybracter_mash_present <- plasmids_match_hybracter_mash |> filter(status == "present")

hybracter_ref_plasmid_lengths_plot <- plot_contig_lengths(plasmids_match_hybracter_mash_present, "Freq Polygon of Contig lengths of Hybracter (hybrid) ref plasmids")
ggsave("for_publication/freq_polygon_hybracter_ref_mash_match_plasmids.png", hybracter_ref_plasmid_lengths_plot, height = 5, width = 7, units = "in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3b. Plasmids- Manually-Currated reference ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Load data ####
# load mash distances
mash <- read.csv("for_publication/mash_cleaned.csv")
#View(mash)  

#load MOB-suite data
mobsuite <- read.csv(file = "for_publication/mobsuite_cleaned.csv")
#View(mobsuite)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * clean and process plasmid data ####
# remove same-assembler matches 
mash_cleaned <- mash |> filter(qry_assembler != ref_assembler)
#View(mash_cleaned)
nrow(mash_cleaned) #147856


#collapse so only 1-way pairs exist:
# Sort by increasing mash_dist
mash_sorted <- mash_cleaned |>
  arrange(mash_dist) |>
  # Create a helper column for identifying symmetric pairs
  mutate(pair_id = pmap_chr(list(sample, qry_assembler, qry_contig, ref_assembler, ref_contig), 
                            ~ paste(sort(c(..1, ..2, ..3, ..4, ..5)), collapse = "|")))

# Group by the pair_id and sample to find symmetric pairs
mash_unique <- mash_sorted |>
  group_by(pair_id, sample) |>
  slice_head(n = 1) |>  # Keep only the first (i.e. lowest mash_dist due to sorting)
  ungroup()

# For reverse match values, join with the original table (but flipped columns)
mash_rev <- mash_sorted |>
  select(sample, qry_assembler, qry_contig, ref_assembler, ref_contig,
         mash_dist, p_value, matching_hashes) |>
  rename(
    ref_assembler_rev = qry_assembler,
    ref_contig_rev = qry_contig,
    qry_assembler_rev = ref_assembler,
    qry_contig_rev = ref_contig,
    rev_mash = mash_dist,
    rev_p_value = p_value,
    rev_matching_hashes = matching_hashes
  )

# Join to add the reverse values
mash_final <- mash_unique |>
  left_join(mash_rev,
            by = c("sample" = "sample",
                   "qry_assembler" = "qry_assembler_rev",
                   "qry_contig" = "qry_contig_rev",
                   "ref_assembler" = "ref_assembler_rev",
                   "ref_contig" = "ref_contig_rev"))


#add a column for mash difference between forwards and rev match
mash_final <- mash_final |>
  mutate(mash_diff = mash_dist - rev_mash) |>
  arrange(desc(mash_diff))

#check
sum(mash_final$mash_diff != 0) # should be 0


#keep only best matching contig from each other assembler
mash_final_best <- mash_final |>
  arrange(mash_dist) |>
  group_by(sample, qry_assembler, qry_contig, ref_assembler) |>
  mutate(match_rank = row_number()) |>
  filter(match_rank == 1) |>
  ungroup()
nrow(mash_final_best) #9288

# repeat for ref assembler, so only keeping 1 best match 
mash_final_best <- mash_final_best |>
  arrange(mash_dist) |>
  group_by(sample, ref_assembler, ref_contig, qry_assembler) |>
  mutate(match_rank = row_number()) |>
  filter(match_rank == 1) |>
  ungroup()
nrow(mash_final_best) #5860
#View(mash_final_best)


#keep only required columns
mash_final_best <- mash_final_best |>
  select(sample, qry_assembler, qry_contig, ref_assembler, ref_contig, mash_dist, p_value, matching_hashes)

# merge with contig annotation data from MOB-suite (only 1000-400,000 conitgs; chromosomes removed)
# clean up mobsuite data first: # Filter only potential plasmids data- remove contigs >400,000 and <1000 and add number of contigs per sample
mobsuite_no_chromosomes <- mobsuite |>
  filter(size <= 400000 & size >= 1000 ) |>
  arrange(sample, assembler, -size) |>
  group_by(sample_id, assembler) |>
  mutate(numcontig = n()) |>
  ungroup() |>
  select(-sample_id) |>
  arrange(sample, assembler, -size) |>
  mutate(count = row_number())
#View(mobsuite_no_chromosomes)

# separate contig id field
mobsuite_no_chromosomes <- mobsuite_no_chromosomes |>
  mutate(contig_id_cleaned = str_extract(contig_id, "^[^ ]+")) |>
  select(assembler, sample, molecule_type, primary_cluster_id, secondary_cluster_id, contig_id_cleaned, 
         size, circularity_status, rep_type.s., relaxase_type.s. , mpf_type, orit_type.s., predicted_mobility, mash_nearest_neighbor)

# merge based on ref to conitgs info.
mash_mob_merged <- mash_final_best |>
  left_join(mobsuite_no_chromosomes,
            by = c("sample" = "sample",
                   "ref_assembler" = "assembler",
                   "ref_contig" = "contig_id_cleaned"))

mash_mob_merged <- mash_mob_merged |>
  left_join(mobsuite_no_chromosomes,
            by = c("sample" = "sample",
                   "qry_assembler" = "assembler",
                   "qry_contig" = "contig_id_cleaned"),
            suffix = c("_ref", "_qry"))

#remoe NA rows- these are the chromosomes or too short contigs
mash_mob_merged <- mash_mob_merged |>
  filter(!is.na(molecule_type_ref) & !is.na(molecule_type_qry))

#normalise 'novel' cluster names, and make a rounded size column for qry and ref
mash_mob_merged <- mash_mob_merged |>
  mutate(primary_cluster_id_ref = if_else(grepl("novel_", primary_cluster_id_ref), "novel", primary_cluster_id_ref)) |>
  mutate(primary_cluster_id_qry = if_else(grepl("novel_", primary_cluster_id_qry), "novel", primary_cluster_id_qry)) |>
  mutate(size_rounded_ref = signif(size_ref, 2)) |>
  mutate(size_rounded_qry = signif(size_qry, 2)) |>
  #add size difference column
  mutate(size_diff = size_ref - size_qry)

nrow(mash_mob_merged) # 4397
#View(mash_mob_merged)


# get plasmid matching (connected components) sets from network analysis:
# load igraph package for network analysis, if not already loaded
#library(igraph)

# Set threshold
mash_thresh <- 0.025

# Create graph edges from filtered MASH pairs
edges <- mash_mob_merged |>
  filter(mash_dist < mash_thresh) |>
  mutate(
    from = paste(sample, qry_assembler, qry_contig, sep = "|"),
    to   = paste(sample, ref_assembler, ref_contig, sep = "|")
  )

# Build undirected graph
g <- graph_from_data_frame(edges |> select(from, to), directed = FALSE)

# Extract connected components
components <- components(g)

# Create data frame of nodes with membership
component_df <- data.frame(
  node = names(components$membership),
  match_set_id = components$membership
) |>
  separate(node, into = c("sample", "assembler", "contig"), sep = "\\|")
#View(component_df)

# Count number of assemblers per sample + match set
mash_plasmid_sets <- component_df |>
  group_by(sample, match_set_id) |>
  mutate(num_assemblers = n_distinct(assembler)) |>
  ungroup()
#View(mash_plasmid_sets)

#check the sizes of the match sets
num_groups <- component_df |>
  group_by(match_set_id) |>
  summarise(n = n())
table(num_groups$n)
#view(num_groups)

#check which match sets have more than 1 of an assembler present
check_repeats <- component_df |>
  group_by(sample, assembler, match_set_id) |>
  summarise(count = n()) |>
  filter(count >1)
#View(check_repeats)
unique(check_repeats$match_set_id) # 2 match sets: 52 (sample pilot_AF43) 301 (pilot_AF57)

# Manually fix groups with >6 contigs present, based on match_set_id
# remove erroneous links
plasmids_mash_manually_cleaned <- mash_mob_merged |>
  # remove incorrectly assigned contigs (manually checked based on all contig/ plasmid info)
  filter(!(sample == "pilot_AJ1" & qry_assembler == "autocycler" & qry_contig == "2" & ref_assembler %in% c("unicycler", "unicycler_bold") & ref_contig == "4")) |>
  filter(!(sample == "pilot_AF43" & qry_assembler == "autocycler" & qry_contig == "3" & ref_assembler %in% c("unicycler", "unicycler_bold") & ref_contig == "2")) |>
# remove link between 2 separate plasmids
  filter(!(sample == "pilot_AF57" & qry_assembler == "autocycler" & qry_contig == "8" & ref_assembler  == "unicycler_bold" & ref_contig == "6")) |>
  filter(!(sample == "pilot_AF57" & qry_assembler == "unicycler_bold" & qry_contig == "6" & ref_assembler == "unicycler" & ref_contig == "8")) 



#get rows to replace:
rows_to_add <- mash_final |>
  filter(
           sample == "pilot_AJ1" & qry_assembler == "autocycler" & qry_contig == "2" & ref_assembler %in% c("unicycler", "unicycler_bold") & ref_contig == "3" | # this column not found
           sample == "pilot_AF43" & qry_assembler == "autocycler" & qry_contig == "2" & ref_assembler %in% c("unicycler", "unicycler_bold") & ref_contig == "2"|
           sample == "pilot_AF57" & qry_assembler == "autocycler" & qry_contig == "7" & ref_assembler  == "unicycler_bold" & ref_contig == "6"|
           sample == "pilot_AF57" & qry_assembler == "unicycler_bold" & qry_contig == "6" & ref_assembler  == "unicycler" & ref_contig == "7"
  )
#View(rows_to_add)
# merge based on ref to conitgs info.
rows_to_add <- rows_to_add |>
  left_join(mobsuite_no_chromosomes,
            by = c("sample" = "sample",
                   "ref_assembler" = "assembler",
                   "ref_contig" = "contig_id_cleaned"))
rows_to_add <- rows_to_add |>
  left_join(mobsuite_no_chromosomes,
            by = c("sample" = "sample",
                   "qry_assembler" = "assembler",
                   "qry_contig" = "contig_id_cleaned"),
            suffix = c("_ref", "_qry"))
rows_to_add <- rows_to_add |>
  select(-c(qry_id, ref_id, pair_id, rev_mash, rev_p_value, rev_matching_hashes, mash_diff))

rows_to_add <- rows_to_add |>
  mutate(size_rounded_ref = signif(size_ref, 2),
         size_rounded_qry = signif(size_qry, 2),
         size_diff = size_ref - size_qry)

plasmids_mash_manually_cleaned <- rbind(plasmids_mash_manually_cleaned, rows_to_add)
#View(plasmids_mash_manually_cleaned)

#Clustering:
# --- Set higher threshold for all possible matches
mash_thresh <- 0.2

# --- Define all node names
df_plasmids_mash_all_nodes <- plasmids_mash_manually_cleaned |>
  select(sample, qry_assembler, qry_contig) |>
  rename(assembler = qry_assembler, contig = qry_contig) |>
  bind_rows(
    plasmids_mash_manually_cleaned |>
      select(sample, ref_assembler, ref_contig) |>
      rename(assembler = ref_assembler, contig = ref_contig)
  ) |>
  distinct()

# --- Add node ID
df_plasmids_mash_all_nodes <- df_plasmids_mash_all_nodes |>
  mutate(node = paste(sample, assembler, contig, sep = "|"))

# --- Create graph edges from filtered MASH pairs
edges <- plasmids_mash_manually_cleaned |>
  filter(mash_dist < mash_thresh) |>
  mutate(
    from = paste(sample, qry_assembler, qry_contig, sep = "|"),
    to   = paste(sample, ref_assembler, ref_contig, sep = "|")
  )

g <- graph_from_data_frame(edges |> select(from, to), directed = FALSE)
components <- components(g)

# --- Match set ID per node
component_df <- data.frame(
  node = names(components$membership),
  match_set_id = components$membership
) |>
  separate(node, into = c("sample", "assembler", "contig"), sep = "\\|")

# --- Full node info with match_set_id
all_plasmids <- df_plasmids_mash_all_nodes |>
  left_join(component_df, by = c("sample", "assembler", "contig"))

# --- Get all assembler combinations to spread
assemblers <- union(
  unique(plasmids_mash_manually_cleaned$qry_assembler),
  unique(plasmids_mash_manually_cleaned$ref_assembler)
)
# --- For each assembler, add mash_dist_* and size_diff_*
for (a in assemblers) {
  all_plasmids <- all_plasmids |>
    rowwise() |>
    mutate(
      !!paste0("mash_dist_", a) := {
        row_match <- plasmids_mash_manually_cleaned |>
          filter(sample == cur_data()$sample,
                 ((qry_assembler == assembler & qry_contig == contig & ref_assembler == a) |
                    (ref_assembler == assembler & ref_contig == contig & qry_assembler == a)))
        if (nrow(row_match) == 0) NA_real_ else row_match$mash_dist[1]
      },
      !!paste0("size_diff_", a) := {
        row_match <- plasmids_mash_manually_cleaned |>
          filter(sample == cur_data()$sample,
                 ((qry_assembler == assembler & qry_contig == contig & ref_assembler == a) |
                    (ref_assembler == assembler & ref_contig == contig & qry_assembler == a)))
        if (nrow(row_match) == 0) NA_real_ else row_match$size_diff[1]
      }
    ) |>
    ungroup()
}

# --- Compute num_non_na, dist_match, and mash_match
all_plasmids <- all_plasmids |>
  rowwise() |>
  mutate(
    num_non_na = sum(!is.na(c_across(starts_with("mash_dist_")))),
    mash_match = sum(c_across(starts_with("mash_dist_")) < 0.025, na.rm = TRUE) > num_non_na / 2,
    dist_match = sum(abs(c_across(starts_with("size_diff_"))) < 10, na.rm = TRUE) > num_non_na / 2
  ) |>
  ungroup()
#View(all_plasmids)

# merge contig info
all_plasmids_contig_info <- all_plasmids |>
  left_join(mobsuite_no_chromosomes,
            by = c("sample" = "sample",
                   "assembler" = "assembler",
                   "contig" = "contig_id_cleaned"))
#View(all_plasmids_contig_info)

# add circularity status match info, and match status
all_plasmids_status_info <- all_plasmids_contig_info |>
  # Step 1: Count number of circular plasmids per match set
  group_by(match_set_id) |>
  mutate(
    num_circular = sum(circularity_status %in% TRUE, na.rm = TRUE),
    true_plasmid_status = num_circular >= 2
  ) |>
  ungroup() |>
  
  # Step 2: Count number of TRUE mash and size matches per match set
  group_by(match_set_id) |>
  mutate(
    num_mash_match = sum(mash_match %in% TRUE, na.rm = TRUE),
    size_match = rowSums(
      across(starts_with("size_diff_"), ~ abs(.) <= 0.1),  # 10% threshold
      na.rm = TRUE
    ) >= ceiling((num_non_na %||% 1) / 2),  # avoid divide by 0
    
    num_size_match = sum(size_match %in% TRUE, na.rm = TRUE)
  ) |>
  ungroup() |>
  
  # Step 3: Determine status per plasmid
  mutate(
    status = case_when(
      mash_match & dist_match & circularity_status ~ "present",
      !mash_match & dist_match & circularity_status ~ "mash_mismatch",
      mash_match & !dist_match & circularity_status ~ "size_mismatch",
      mash_match & dist_match & !circularity_status ~ "circ_mismatch",
      mash_match & !dist_match & !circularity_status ~ "circ_size_mismatch",
      !mash_match & dist_match & !circularity_status ~ "mash_circ_mismatch",
      !mash_match & !dist_match ~ "absent",
      TRUE ~ NA_character_
    )
  )

#View(all_plasmids_status_info)
nrow(all_plasmids_status_info) #1934

#optional checks
#table(all_plasmids_status_info$status, all_plasmids_status_info$assembler, useNA = "ifany")
#table(all_plasmids_status_info$true_plasmid_status, all_plasmids_status_info$assembler, useNA = "ifany")

# filter for true plasmids only (at least 2 circular for that match set)
all_plasmids_status_true <- all_plasmids_status_info |>
  filter(true_plasmid_status)

# ensure all combos of match set and plasmid exist
# Vector of all 6 assemblers you expect
all_assemblers <- c("autocycler", "flye", "hybracter_long", "hybracter_hybrid", "unicycler", "unicycler_bold")

# Create complete dataframe with all combinations of match_set_id × assembler
all_plasmids_status_true_completed <- all_plasmids_status_true |>
  complete(match_set_id, assembler = all_assemblers) |>
  mutate(status = if_else(is.na(status), "absent", status))

#add simple status
all_plasmids_status_true_completed <- all_plasmids_status_true_completed |>
  mutate(status_simple = case_when(status == "present"~ "present",
                                   status == "absent" ~ "absent",
                                   TRUE ~ "misassembled"))


#check how many present for each match_set
match_set_check <- all_plasmids_status_true_completed |>
  group_by(match_set_id) |>
  summarise(count = n()) |>
  arrange(count) |>
  filter(count != 6)
#View(match_set_check) # should be all at 6

#save compelted table
write.csv(all_plasmids_status_true_completed, "for_publication/plasmids_match_manual_mash.csv", row.names = FALSE)

# check match set NAs!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Tables of detailed status with mash distance ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Summary table  ####
# load packages if not already loaded
#library(dplyr)
#library(tidyr)
#library(gtsummary)
#library(gt)

# Step 1: Calculate number of unique true plasmids
num_manual_plasmids <- length(unique(all_plasmids_status_true_completed$match_set_id)) #303

# Define assembler order
assembler_order <- c("autocycler", "flye", "hybracter_long", "hybracter_hybrid", "unicycler", "unicycler_bold")

# Step 2: Summary Table – Simple Status
summary_simple <- all_plasmids_status_true_completed |>
  count(assembler, status_simple) |>
  mutate(percent = round(n / num_manual_plasmids * 100, 1),
         n_percent = paste0(n, " (", percent, "%)")) |>
  select(-n, -percent) |>
  pivot_wider(names_from = assembler, values_from = n_percent, values_fill = "0 (0%)") |>
  arrange(factor(status_simple, levels = c("present", "misassembled", "absent"))) |>
  rename(Status = status_simple)

# Step 3: Summary Table – Detailed Status
summary_detailed <- all_plasmids_status_true_completed |>
  count(assembler, status) |>
  mutate(percent = round(n / num_manual_plasmids * 100, 1),
         n_percent = paste0(n, " (", percent, "%)")) |>
  select(-n, -percent) |>
  pivot_wider(names_from = assembler, values_from = n_percent, values_fill = "0 (0%)") |>
  arrange(factor(status, levels = c(
    "present", "circ_mismatch", "size_mismatch", "mash_mismatch",
    "circ_size_mismatch", "size_mash_mismatch", "circ_mash_mismatch", "absent"
  ))) |>
  rename(Status = status)

# Step 4: Save as CSV (optional)
#write.csv(summary_simple, "plasmids_mash_match_manual_simple_summary_table.csv", row.names = FALSE)
#write.csv(summary_detailed, "plasmids_mash_match_manual_detailed_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~#
# * Stats tests ####

contingency_simple <- all_plasmids_status_true_completed |>
  group_by(assembler) |>
  summarise(present = sum(status == "present"),
            not_present = sum(status != "present"))

# Ensure expected order and numeric conversion
contingency_simple$assembler <- factor(contingency_simple$assembler, levels = assembler_order)
contingency_simple <- contingency_simple |> mutate(across(c(present, not_present), as.numeric))

contingency_mat <- as.matrix(contingency_simple |> select(present, not_present))
rownames(contingency_mat) <- contingency_simple$assembler

# Chi-squared test
chisq_present_result <- chisq.test(contingency_mat)
print(chisq_present_result)

# Pairwise Fisher's exact test
manual_ref_plasmids_fisher_test_results <- rstatix::pairwise_fisher_test(contingency_mat) # function defined above
print(manual_ref_plasmids_fisher_test_results)


# Detailed metrics chi-squared for each
contingency_detailed <- all_plasmids_status_true_completed |>
  group_by(assembler, status) |>
  summarise(count = n()) |>
  pivot_wider(id_cols = assembler, names_from = status, values_from = count, values_fill = 0)

# Total reference plasmids
total_plasmids <- length(unique(all_plasmids_status_true_completed$match_set_id))

# Define metrics for testing
metrics <- colnames(contingency_detailed)[-1]  # excluding assembler

# Run chi-squared tests
chisq_results <- lapply(metrics, function(metric) {
  counts <- contingency_detailed[[metric]]
  others <- total_plasmids - counts
  contingency_mat <- rbind(counts, others)
  chisq.test(contingency_mat)$p.value
})

chisq_df <- data.frame(metric = metrics,p_value = unlist(chisq_results))

print(chisq_df)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Plots, using less detailed status ####
# load packages if not already loaded
#library(tidyr)
#library(ggplot2)
#library(patchwork)
#library(ggforce)  # for geom_circle
#library(colorspace)  # for lighten

#add assembler labels if not already in workspace
assembler_labels <- c(
  "autocycler" = "Autocycler",
  "flye" = "Flye",
  "hybracter_long" = "Hybracter\n(long)",
  "hybracter_hybrid" = "Hybracter\n(hybrid)",
  "unicycler" = "Unicycler",
  "unicycler_bold" = "Unicycler\n(bold)"
)

# Define assembler colors if not already in workspace
assembler_colors <- c(
  "Autocycler" = "#332288",
  "Flye" = "#88ccee",
  "Hybracter\n(long)" = "#44aa99",
  "Hybracter\n(hybrid)" = "#117733",
  "Unicycler" = "#cc6677",
  "Unicycler\n(bold)" = "#882255"
)

status_levels <- c("absent", "misassembled", "present")
status_alpha <- rev(c(1.0, 0.5, 0.1))

# Prep long format
df_long <- all_plasmids_status_true_completed %>%
  mutate(assembler_label = recode(assembler, !!!assembler_labels)) |>
  select(match_set_id, assembler_label, status_simple) %>%
  mutate(status_simple = factor(status_simple, levels = status_levels))

# Pivot to wide, then group by unique combinations
heatmap_data_wide <- df_long %>%
  pivot_wider(names_from = assembler_label, values_from = status_simple, values_fill = "absent") %>%
  group_by(across(all_of(names(assembler_colors)))) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(combo_id = row_number())  # New order: most common to least

# Long format again for heatmap
heatmap_long <- heatmap_data_wide %>%
  pivot_longer(cols = names(assembler_colors), names_to = "assembler", values_to = "status") %>%
  mutate(
    assembler = factor(assembler, levels = rev(names(assembler_colors))),  # Reverse for top-down
    combo_id = factor(combo_id, levels = unique((combo_id)))  # Keep left-to-right ordering
  )

# Get fill color with transparency
status_to_alpha <- function(status) {
  status_alpha[which(status_levels == status)]
}
status_to_base_color <- function(assembler) {
  assembler_colors[[assembler]]
}
heatmap_long <- heatmap_long %>%
  mutate(
    fill_color = mapply(function(a, s) {
      adjustcolor(status_to_base_color(as.character(a)), alpha.f = status_to_alpha(s))
    }, assembler, status)
  )



# Ensure combo_id is a factor with consistent levels
heatmap_long$combo_id <- factor(heatmap_long$combo_id, levels = unique(heatmap_long$combo_id))
heatmap_data_wide$combo_id <- factor(heatmap_data_wide$combo_id, levels = levels(heatmap_long$combo_id))


# Heatmap circles
heatmap_plot <- ggplot(heatmap_long, aes(x = combo_id, y = assembler)) +
  geom_point(aes(fill = fill_color), shape = 21, size = 5, stroke = 0.3, color = "black") +
  scale_fill_identity() +
  theme_minimal(base_size = 12) +
  labs(x = "Unique status combinations", y = "Assembler") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Top bar plot: count of plasmids per combination
top_bar <- ggplot(heatmap_data_wide, aes(x = factor(combo_id, levels = combo_id), y = count)) +
  geom_col(fill = "grey40") +
  geom_text(aes(label = count), vjust = -0.2, size = 3.2) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme_minimal(base_size = 12) +
  labs(y = "Plasmid count", x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 10, r = 5, b = 0, l = 39))

# Left stacked bar: count of each status per assembler
left_bar_data <- df_long |>
  count(assembler_label, status_simple) |>
  mutate(
    assembler_label = factor(assembler_label, levels = rev(names(assembler_colors))),
    
    # Reverse the levels so stacking starts with "other_misassembly" at the bottom
    status_with_mash_0.025 = factor(
      status_simple,
      levels = c("absent", "misassembled", "present")
    ),
    
    fill_color = mapply(function(a, s) {
      adjustcolor(assembler_colors[[as.character(a)]], alpha.f = status_to_alpha(s))
    }, assembler_label, status_simple)
  )

# Calculate position of label: total stacked value per assembler
label_data <- left_bar_data |>
  group_by(assembler_label) |>
  mutate(cumulative = cumsum(n)) |> 
  mutate(n_percent = paste0(round(n/cumulative *100), "%")) |>
  ungroup() |>
  filter(status_simple == "present")

left_bar <- ggplot(left_bar_data, aes(x = n, y = assembler_label, fill = factor(fill_color, levels = rev(fill_color)))) +
  geom_col(position = "stack") +
  geom_text(data = label_data,
            aes(x = cumulative + 1, y = assembler_label, label = n_percent),  # Add offset to the right
            inherit.aes = FALSE,
            hjust = 0, size = 3.2) +
  scale_fill_identity() +
  theme_minimal(base_size = 12) +
  labs(x = "Plasmid count", y = NULL) +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 0)) +
  coord_cartesian(clip = "off") 


# Blank spacer
blank_spacer <- plot_spacer()

# Layout
final_plot <- (blank_spacer + top_bar) / (left_bar + heatmap_plot) +
  plot_layout(widths = c(0.25, 1), heights = c(0.3, 1))

final_plot

#optional save
#ggsave("plasmids_mash_0_025_match_manually_curated.png", final_plot, width = 12, height = 6, units = "in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Frequency Polygon plot:
#call function defined above to plot frequency polygon
all_plasmids_status_true_completed_present <- all_plasmids_status_true_completed |>
  filter(status == "present")

manual_ref_plasmid_lengths_plot <- plot_contig_lengths(all_plasmids_status_true_completed_present, "Freq Polygon of Contig lengths of Manually-Curated ref plasmids")
#optional save
#ggsave("freq_polygon_manually-curated_ref_mash_match_plasmids.png", manual_ref_plasmid_lengths_plot, height = 5, width = 7, units = "in")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Nucloetide-level Accuracy (SNPs, Indels, Average Gene Length) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Load data ####
# load snp, indel and QV data
assembly_accuracy <- read.csv("assembly_nucleotide_accuracy_cleaned.csv")
#View(assembly_accuracy)


#load CheckM2 data
checkm2 <- read.csv("checkm2_cleaned.csv")
#View(checkm2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Summary Stats by assembler ####
summarise_metrics <- function(df, column_indexes) {
  df_summary_table <- data.frame()
  for (col_index in column_indexes) {
    col_name <- colnames(df)[col_index]
    summary_table <- summary_stats(df, "assembler", col_name)
    summary_table <- cbind(metric = col_name, summary_table)
    df_summary_table <- rbind(df_summary_table, summary_table)
  }
  return(df_summary_table)
}

#Standardise SNPs/Indels per 1Mb bp:
assembly_accuracy <- assembly_accuracy |>
  mutate(snps_per_Mb = substitution_errors/assembly_size * 1000000,
         indels_per_Mb = indel_errors/assembly_size * 1000000,
         qv_rounded = round(consensus_qv_before_polishing, 0)) |>
  select(assembler, sample, substitution_errors, indel_errors, assembly_size, consensus_quality_before_polishing, consensus_qv_before_polishing, snps_per_Mb, indels_per_Mb, assembler_tick, polishing, assembler_family, qv_rounded)


assembly_accuracy_summary <- summarise_metrics(assembly_accuracy, 3:9)
#View(assembly_accuracy_summary)
#optional save summary table
#write.csv(assembly_accuracy_summary, file = "assembly_accuracy_summary.csv", row.names = FALSE)


#Repeat for CheckM2 for Coding Density (col 7) and Average gene length (col 9)
checkm2_summary <- summarise_metrics(checkm2, c(7,9))
#View(checkm2_summary)
#optional save summary table
#write.csv(checkm2_summary, file = "checkm2_summary.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Stats Test ####
#Kruskall-Walis overall test
snps_kruskall_result <- kruskal.test(snps_per_Mb ~ assembler, data = assembly_accuracy)
indels_kruskall_result <- kruskal.test(indels_per_Mb ~ assembler, data = assembly_accuracy)
qv_kruskall_result <- kruskal.test(consensus_qv_before_polishing ~ assembler, data = assembly_accuracy)

CD_kruskall_result <- kruskal.test(Coding_Density ~ assembler, data = checkm2)
agl_kruskall_result <- kruskal.test(Average_Gene_Length ~ assembler, data = checkm2)



#Paiwise Wilcoxom signed rank test
assembly_accuracy <- assembly_accuracy |> arrange(sample)
metrics_accuracy <- c("snps_per_Mb", "indels_per_Mb", "consensus_qv_before_polishing")

wilcox_results_list <- list()
for (metric in metrics_accuracy) {
  cat("\nWilcoxon Signed-Rank Test for:", metric, "\n")
  wilcox_results <- pairwise.wilcox.test(x = assembly_accuracy[[metric]], 
                                         g = assembly_accuracy$assembler, 
                                         paired = TRUE, 
                                         p.adjust.method = "bonferroni")
  print(wilcox_results)
  # Convert results to a tidy data frame
  wilcox_df <- as.data.frame(wilcox_results$p.value) 
  # Store in the list
  wilcox_results_list[[metric]] <- wilcox_df
}

# Combine all results into a single data frame and save
final_wilcox_results <- bind_rows(wilcox_results_list)
#View(final_wilcox_results)

# save wilcox results optional
# write.csv(final_wilcox_results, "assembly_accuracy_wilcoxon_results.csv", row.names = FALSE)

#Repeat for checkm2
checkm2 <- checkm2 |> arrange(Name)
metrics_checkm2 <- c("Coding_Density", "Average_Gene_Length")

wilcox_results_list <- list()
for (metric in metrics_checkm2) {
  cat("\nWilcoxon Signed-Rank Test for:", metric, "\n")
  wilcox_results <- pairwise.wilcox.test(x = checkm2_results[[metric]], 
                                         g = checkm2_results$assembler, 
                                         paired = TRUE, 
                                         p.adjust.method = "bonferroni")
  print(wilcox_results)
  # Convert results to a tidy data frame
  wilcox_df <- as.data.frame(wilcox_results$p.value) 
  # Store in the list
  wilcox_results_list[[metric]] <- wilcox_df
}
# Combine all results into a single data frame and save
final_wilcox_results <- bind_rows(wilcox_results_list)
#View(final_wilcox_results)
#optional save pairwise wilcox results
#write.csv(final_wilcox_results, "checkm2_pairwise_wilcoxon_results.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~#
# * Plots ####
# plot function for any metric
plot_metric <- function(data, metric, metric_label, title_lab, scaling_for_labels = 0.9) {
  
  
  # Define colors based on "assembler"
  assembler_colors <- c(
    "autocycler" = "#332288",             # Dark Blue for autocycler
    "autocycler_medaka" = "#5E56A6",      # Lighter Blue for autocycler_medaka
    "autocycler_medaka_full" = "#8A84C5", # Even lighter Blue for autocycler_medaka_full
    "autocycler_pypolca" = "#B1A0E3",  # Lightest Blue for autocycler_polypolish_pypolca
    "flye" = "#88ccee",                  # Light Blue for flye
    "flye_medaka" = "#B2D8E9",           # Lighter Blue for flye_medaka
    "flye_medaka_full" = "#C7E6F2",      # Even lighter Blue for flye_medaka_full
    "flye_pypolca" = "#D9F1F9",  # Lightest Blue for flye_polypolish_pypolca
    "hybracter_long" = "#44aa99",        # Teal for hybracter_long
    "hybracter_hybrid" = "#117733",      # Green for hybracter_hybrid
    "unicycler" = "#cc6677",             # Red for unicycler
    "unicycler_bold" = "#882255"         # Dark Red for unicycler_bold
  )
  
  # Define shapes based on "polishing" category
  polishing_shapes <- c(
    "unpolished" = 4,                  # open circle
    "long-read (subsampled)" = 24,      # open Triangle
    "long-read (all reads)" = 22,       # open Square
    "short-read" = 21,                  # star
    "NA (Unicycler)" = 18,  # Solid diamond
    "NA (Unicycler-bold)" = 15  #solid square
  )
  
  # Define x-axis group labels and horizontal lines
  assembler_labels <- data.frame(
    x = c(2.5, 6.5, 9.3, 12),  # Midpoints of each group
    y = min(data[[metric]], na.rm = TRUE) * scaling_for_labels,  # Position below lowest value
    label = c("Autocycler", "Flye", "Hybracter", "Unicycler")
  )
  
  assembler_lines <- data.frame(
    x_start = c(0.7, 4.7, 8.7, 10.7),
    x_end = c(4.3, 8.3, 10.3, 12.3),
    y = min(data[[metric]], na.rm = TRUE) * scaling_for_labels  # Just below x-axis
    # change to 0.9 for pypolca plots to have some more distance
  )
  
  metric_plot <- ggplot(data, aes(x = assembler, y = !!sym(metric))) +
    # geom_violin(aes(fill = assembler), color = "black") +  # Black boxplot borders
    geom_boxplot(aes(fill = assembler), width = 0.8, outlier.shape = NA, color = "black") +  # Black boxplot borders
    geom_point(aes(shape = polishing, fill = assembler), 
               size = 1.5, stroke = 0.5, color = "black", 
               alpha = 0.8, position = position_jitter(width = 0.1, height = 0)) +  
    labs(x = "", y = metric_label , title = title_lab) +
    scale_fill_manual(values = assembler_colors, guide = "none") +  # Hide fill legend
    scale_shape_manual(name = "Polishing", values = polishing_shapes) +  # Keep only shape legend
    scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +  
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank(),  
      plot.margin = margin(10, 10, 14, 10)  
    ) +
    coord_cartesian(clip = "off") +  
    geom_segment(data = assembler_lines, aes(x = x_start, xend = x_end, y = y, yend = y),
                 inherit.aes = FALSE, color = "black", size = 0.8) +  
    geom_text(data = assembler_labels, aes(x = x, y = y, label = label),  
              vjust = 3.5, size = 5, fontface = "bold", inherit.aes = FALSE) 
  
  return(metric_plot)
}

# Call plot function
#Assembly Accuracy metrics
# SNPs / Mb
plot_accuracy <- plot_metric(assembly_accuracy, "snps_per_Mb" , "SNPs/1,000,000 bp" ,  "Substitution errors corrected by Illumina short-read alignment by Assembler (sup)", 0.9)
#optional save
#ggsave("assembly_accuracy_substitutions_perMb_plot.png", plot_accuracy, width = 7.5, height = 5, units = "in")

#Indels / Mb
plot_accuracy <- plot_metric(assembly_accuracy, "indels_per_Mb" , "Indels/1,000,000 bp" , "Insertions errors corrected by Illumina short-read alignment by Assembler (sup)", 0.9)
#optional save
#ggsave("assembly_accuracy_indels_perMb_plot.png", plot_accuracy, width = 7.5 , height = 5, units = "in")

#QV
plot_accuracy <- plot_metric(assembly_accuracy, "consensus_qv_before_polishing" , "Quality Value" , "Assembly Quality Value Before Short-read Alignment Error Correction by Assembler (sup)", 0.9)
#optional save
#ggsave("assembly_accuracy_QV_plot.png", plot_accuracy, width = 7.5, height = 5, units = "in")


#CheckM2 metrics
#Avergae Gene Length
plot_accuracy <- plot_metric(checkm2, "Average_Gene_Length" , "Average Gene Length" , "Average Gene Length by Assembler- Polisher Combination", 0.999)
#optional save
#ggsave("checkm2_average_gene_length_plot.png", plot_accuracy, width = 7.5, height = 5, units = "in")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. MLST ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Load data ####
mlst <- read.csv("mlst_cleaned.csv")


# * Summary Stats ####
summarise_mlst <- function(df) {
  mlst_summary <- df |>
    pivot_wider(
      id_cols = c(sample, species),
      names_from = assembler,
      values_from = c(mlst, gene1, gene2, gene3, gene4, gene5, gene5_2, gene6)
    ) 
  return(mlst_summary)
}

#call summary function
mlst_wide<- summarise_mlst(mlst)

# save table (optional)
#write.csv(mlst_summary, file = "mlst_summary.csv")

#Summarise how many mlsts correctly identified by each assembler

#Identify the 'correct' MLST profile/ alleles from the most common values among assemblers
# Identify column groups
mlst_cols <- grep("^mlst_", names(mlst_wide), value = TRUE)
gene_cols <- grep("^gene[0-9]_.*|^gene5_2_.*", names(mlst_wide), value = TRUE)
assemblers <- c("autocycler", "autocycler_medaka","autocycler_medaka_full", "autocycler_polypolish_pypolca",
                "flye", "flye_medaka","flye_medaka_full", "flye_polypolish_pypolca",
                "hybracter_long", "hybracter_hybrid", 
                "unicycler", "unicycler_bold")


# Function to get the most common value per row
most_common <- function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
}

# Compute "real" values
mlst_real <- mlst_wide |>
  rowwise() |>
  mutate(
    real_mlst = list(most_common(c_across(all_of(mlst_cols)))),
    real_gene1 = most_common(c_across(grep("^gene1_", names(mlst_summary_sup_contam_free), value = TRUE))),
    real_gene2 = most_common(c_across(grep("^gene2_", names(mlst_summary_sup_contam_free), value = TRUE))),
    real_gene3 = most_common(c_across(grep("^gene3_", names(mlst_summary_sup_contam_free), value = TRUE))),
    real_gene4 = most_common(c_across(grep("^gene4_", names(mlst_summary_sup_contam_free), value = TRUE))),
    real_gene5 = most_common(c_across(grep("^gene5_", names(mlst_summary_sup_contam_free), value = TRUE))),
    real_gene5_2 = most_common(c_across(grep("^gene5_2_", names(mlst_summary_sup_contam_free), value = TRUE))),
    real_gene6 = most_common(c_across(grep("^gene6_", names(mlst_summary_sup_contam_free), value = TRUE))),
    
  ) |>
  ungroup() 
#View(mlst_real)

# Cross-reference assemblers to true value and summarise which are 'correct'
for (assembler in assemblers) {
  # Ensure MLST match comparison is done on character vectors
  mlst_real[[paste0("mlst_match_", assembler)]] <- 
    as.character(mlst_real[[paste0("mlst_", assembler)]]) == as.character(mlst_real$real_mlst)
  # Extract correct gene columns
  gene_match_cols <- grep(paste0("^gene[1-6]_", assembler, "$|^gene5_2_", assembler, "$"), names(mlst_real), value = TRUE)
  real_gene_cols <- paste0("real_", sub(paste0("_", assembler, "$"), "", gene_match_cols))
  mlst_real[[paste0("num_genes_match_", assembler)]] <- 
    rowSums(mlst_real[gene_match_cols] == mlst_real[real_gene_cols], na.rm = TRUE)
}

#convert lists to strings before saving
mlst_real <- data.frame(lapply(mlst_real, function(x) {
  if (is.list(x)) sapply(x, toString) else x}), stringsAsFactors = FALSE)
#View(mlst_real)
nrow(mlst_real)


#filter out rows where no MLST typing sheme available:
mlst_real <- mlst_real |>
  filter(
    real_gene1 != "" & !is.na(real_gene1),
    real_gene2 != "" & !is.na(real_gene2),
    real_gene3 != "" & !is.na(real_gene3),
    real_gene4 != "" & !is.na(real_gene4),
    real_gene5 != "" & !is.na(real_gene5),
    real_gene5_2 != "" & !is.na(real_gene5_2),
    real_gene6 != "" & !is.na(real_gene6)
  )

nrow(mlst_real) # removed reow for S. marcascens with no typing scheme available

# Summarize per assembler
mlst_summary <- mlst_real |>
  # For MLST match columns: count TRUE values for each assembler
  summarise(across(starts_with("num_genes_match"), ~mean(., na.rm = TRUE), .names = "avg_num_genes_match_{sub('num_genes_match_', '', .col)}"),
            # Adding a column to check if 7/7 genes match, which will also count as MLST match
            across(starts_with("num_genes_match"), ~sum(. == 7, na.rm = TRUE), .names = "mlst_match_7_genes_{sub('num_genes_match_', '', .col)}"))

#select only overall summary of number of isolates mathcing all 7 alleles
mlst_summary_selected <- mlst_summary |>
  select(starts_with("mlst_match_7_genes_"))

# View the summary table
#View(mlst_summary_selected)
#write.csv(mlst_summary_selected, file = "mlst_summary.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Stats Test ####
# Chi-squared
mlst_summary_transposed <- t(mlst_summary_selected)
#View(mlst_summary_transposed)
colnames(mlst_summary_transposed) <- "correct_mlst"
mlst_summary_transposed <- data.frame(mlst_summary_transposed)

mlst_summary_transposed <- mlst_summary_transposed |>
  mutate(no_mlst_match = nrow(mlst_real)- correct_mlst) 
mlst_chi2_results <- chisq.test(mlst_summary_transposed) # 0.7965


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. AMRFinder Plus ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Load data ####
amrfinder <- read.csv("for_publication/amrfinderplus_cleaned.csv")
#View(amrfinder)
# note for 'strand', negative strand "-" loaded as NA


#~~~~~~~~~~~~~~~~~~~~~~#
# * Summary Stats ####
#produce summary by contig and gene type
summarise_elements <- function(df) {
  df_element_summary <- df |>
    group_by(across(c(assembler, sample, polishing, assembler_family, Type))) |>
    summarise(no_genes = n(), .groups = 'drop') |>
    pivot_wider(names_from = Type, values_from = no_genes) |>
    mutate(
      AMR = replace_na(AMR, 0),
      STRESS = replace_na(STRESS, 0),
      VIRULENCE = replace_na(VIRULENCE, 0)
    ) |>
    mutate(allgenes = AMR + STRESS +VIRULENCE)
  return( df_element_summary)
}

#summarise AMR, STRESS and VIRULENCE genes, and overall gene count across all samples
amrfinder_elements_by_sample <- summarise_elements(amrfinder)
#View(amrfinder_elements_by_sample)

amrfinder_summary <- summarise_metrics(amrfinder_elements_by_sample, 5:8)
#View(amrfinder_summary)
#optional save summary table
#write.csv(amrfinder_summary, file = "amrfinder_summary.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~#
# * Stats Test ####
#Stats test
#Kruskall wallis test for oveall differneces beween assemblers 
kruskal.test(AMR ~ assembler, data = amrfinder_elements_by_sample) 
kruskal.test(STRESS ~ assembler, data = amrfinder_elements_by_sample) 
kruskal.test(VIRULENCE ~ assembler, data = amrfinder_elements_by_sample) 
kruskal.test(allgenes ~ assembler, data = amrfinder_elements_by_sample) 


#Pairwise Wilcoxon rank sum 
amrfinder_elements_by_sample <- amrfinder_elements_by_sample |> arrange(sample)
amrfinder_wilcox_result <- data.frame()
for (col_index in c(5:8)) {
  col_name <- colnames(amrfinder_elements_by_sample)[col_index]
  wilcox_result <- pairwise.wilcox.test(amrfinder_elements_by_sample[[col_name]], amrfinder_elements_by_sample$assembler, p.adjust.method = "bonf", paired = TRUE)
  summary_table <- cbind(metric = col_name, wilcox_result$p.value)
  amrfinder_wilcox_result <- rbind(amrfinder_wilcox_result, summary_table)
}
#View(amrfinder_wilcox_result)
#optional save wilcox results
#write.csv(amrfinder_wilcox_result, file = "amrfinder_wilcox_result.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. Bakta Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Load data ####
bakta_by_contig <- read.csv("for_publication/bakta_by_contig_cleaned.csv")
View(bakta_by_contig)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Summary Stats ####
#Summarise by sample (so 1 row per sample/assembler combo)
summarise_by_sample <- function(df) {
  df_sum <- df |>
    group_by(assembler, sample, polishing, assembler_family) |>
    summarise(across(2:15, sum, na.rm = TRUE))
  return(df_sum)
}

bakta_by_sample<- summarise_by_sample(bakta_by_contig)
#View(bakta_by_sample)
#optional save
#write.csv(bakta_by_sample, file = "bakta_by_sample.csv", row.names = FALSE)

#summarise total cds per sample
bakta_summary <- summarise_metrics(bakta_by_sample, 5)
#View(bakta_summary)
#optional save summary table
#write.csv(bakta_summary, "bakta_summary.csv", row.names = FALSE)


# * Stats Test ####
#Kruskall-Walis overall test
bakta_kruskall_result <- kruskal.test(CDS ~ assembler, data = bakta_by_sample)


#Paiwise Wilcoxom signed rank test
bakta_by_sample <- bakta_by_sample |> arrange(sample)
metrics_bakta <- c("CDS")

wilcox_results_list <- list()
for (metric in metrics_bakta) {
  cat("\nWilcoxon Signed-Rank Test for:", metric, "\n")
  wilcox_results <- pairwise.wilcox.test(x = bakta_by_sample[[metric]], 
                                         g = bakta_by_sample$assembler, 
                                         paired = TRUE, 
                                         p.adjust.method = "bonferroni")
  print(wilcox_results)
  # Convert results to a tidy data frame
  wilcox_df <- as.data.frame(wilcox_results$p.value) 
  # Store in the list
  wilcox_results_list[[metric]] <- wilcox_df
}

# Combine all results into a single data frame and save
final_wilcox_results <- bind_rows(wilcox_results_list)
#View(final_wilcox_results)

# optional save wilcox results optional
# write.csv(final_wilcox_results, "bakta_cds_wilcoxon_results.csv", row.names = FALSE)
