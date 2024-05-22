library(tidyverse)
library(readr)

# Read in the counts table
knownCanonicalTSS_readCounts_table <- read_delim("knownCanonicalTSS.readCounts.table.txt", 
                                                 delim = "\t", escape_double = FALSE,
                                                 trim_ws = TRUE)


# Reshape all un-normalized data to long format
long_data_nonorm <- knownCanonicalTSS_readCounts_table %>%
  pivot_longer(cols = contains("me3"), # Only include the data columns
               names_to = "Condition_Replicate",
               values_to = "ReadCounts") %>%
  mutate(LogReadCounts = log(ReadCounts + 1)) # Adding 1 to avoid log(0)
# by default, log() is base e.

# Create violin plots for each sample with the un-normalized data.
ggplot(long_data_nonorm, aes(x = Condition_Replicate, y = LogReadCounts)) +
  geom_violin(trim = FALSE) + # 'trim=FALSE' to include all data points in the plot shape
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate sample names for readability
  labs(title = "Violin Plots of LogReadCounts for Each Sample",
       x = "Sample",
       y = "Raw Log Read Counts")
ggsave("allPromoters.noNorm.violinPlot.pdf")

# Reshape the raw data to long format for K27me3 samples
long_data_k27me3 <- knownCanonicalTSS_readCounts_table %>%
  pivot_longer(cols = starts_with("K27me3"),
               names_to = "Condition_Replicate",
               values_to = "ReadCounts") %>%
  mutate(LogReadCounts = log(ReadCounts + 1)) # Adding 1 to avoid log(0)
# by default, log() is base e.

# Reshape the raw data to long format for K4me3 samples
long_data_k4me3 <- knownCanonicalTSS_readCounts_table %>%
  pivot_longer(cols = starts_with("K4me3"),
               names_to = "Condition_Replicate",
               values_to = "ReadCounts") %>%
  mutate(LogReadCounts = log(ReadCounts + 1)) # Adding 1 to avoid log(0)

# Create the violin plot for the K27me3 data
ggplot(long_data_k27me3, aes(x = Condition_Replicate, y = LogReadCounts, fill = Condition_Replicate)) +
  geom_violin() +
  scale_fill_manual(values = c("K27me3.Minus.R1" = "red", "K27me3.Minus.R2" = "red",
       "K27me3.Plus.R1" = "green", "K27me3.Plus.R2" = "green")) +
  labs(title = "Log-transformed Violin Plot of K27me3",
             x = "Condition and Replicate",
             y = "Log-transformed Read Counts") +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Improve label readability
ggsave("allPromoters.K27me3.noNorm.violinPlot.pdf")

# Create the violin plot for the K4me3 data
ggplot(long_data_k4me3, aes(x = Condition_Replicate, y = LogReadCounts, fill = Condition_Replicate)) +
  geom_violin() +
  scale_fill_manual(values = c("K4me3.Minus.R1" = "red", "K4me3.Minus.R2" = "red",
                               "K4me3.Plus.R1" = "green", "K4me3.Plus.R2" = "green")) +
  labs(title = "Log-transformed Violin Plot of K4me3",
       x = "Condition and Replicate",
       y = "Log-transformed Read Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Improve label readability
ggsave("allPromoters.K4me3.noNorm.violinPlot.pdf")

# Create the violin plot for all noNorm data
ggplot(long_data_nonorm, aes(x = Condition_Replicate, y = LogReadCounts, fill = Condition_Replicate)) +
  geom_violin() +
  scale_fill_manual(values = c("K4me3.Minus.R1" = "red", "K4me3.Minus.R2" = "red",
                               "K4me3.Plus.R1" = "green", "K4me3.Plus.R2" = "green",
                               "K27me3.Minus.R1" = "red", "K27me3.Minus.R2" = "red",
                               "K27me3.Plus.R1" = "green", "K27me3.Plus.R2" = "green")) +
  labs(title = "Log-transformed Violin Plot",
       x = "Condition and Replicate",
       y = "Log-transformed Read Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Improve label readability
ggsave("allPromoters.noNorm.violinPlot.pdf")

# Now lets make scatter plots with each point representing a gene.
# Calculate the log of the mean for each pair of replicates, with a 1 pseudocount added.
data_means_Minus <- knownCanonicalTSS_readCounts_table %>%
  rowwise() %>%
  mutate(LogMean_K27me3_Minus = log(mean(c(K27me3.Minus.R1+1, K27me3.Minus.R2+1))),
         LogMean_K4me3_Minus = log(mean(c(K4me3.Minus.R1+1, K4me3.Minus.R2+1)))) %>%
  ungroup()

# Add a new column called Label, that has the gene name without the transcript detail
data_means_Minus <- data_means_Minus %>%
  mutate(Label = str_extract(Name, "^[^.]+"))


# Read in the files with the names of poised promoters
minus_gene_list <- readLines("poisedPromoters.A_minus_only.Venny.txt")
print(length(minus_gene_list))
plus_gene_list <- readLines("poisedPromoters.A_plus_only.Venny.txt")
print(length(plus_gene_list))
both_gene_list <- readLines("poisedPromoters.A_minus_A_plus.Venny.txt")
print(length(both_gene_list))

# Repeat for data_means_Plus
data_means_Plus <- knownCanonicalTSS_readCounts_table %>%
  rowwise() %>%
  mutate(LogMean_K27me3_Plus = log(mean(c(K27me3.Plus.R1+1, K27me3.Plus.R2+1))),
         LogMean_K4me3_Plus = log(mean(c(K4me3.Plus.R1+1, K4me3.Plus.R2+1)))) %>%
  ungroup()
data_means_Plus <- data_means_Plus %>%
  mutate(Label = str_extract(Name, "^[^.]+"))

# Create a scatter plot for the Minus cells
ggplot(data_means_Minus, aes(x = LogMean_K27me3_Minus, y = LogMean_K4me3_Minus)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 10))+
  scale_y_continuous(limits = c(0, 10))+
  theme_minimal() +
  labs(title = "Log of Mean Replicates Comparison, Minus Cells",
       x = "Log of Minus Mean (K27me3)",
       y = "Log of Minus Mean (K4me3)") +
  geom_point(alpha = 0.5) +  # Plot points
  geom_point(data = filter(data_means_Minus, Label %in% minus_gene_list), color = "red", alpha = 0.5)

ggsave("MinusCells.MinusPoised.min0.K4me3.vs.K27me3.scatter.pdf")

# Create a scatter plot for the Plus cells
ggplot(data_means_Plus, aes(x = LogMean_K27me3_Plus, y = LogMean_K4me3_Plus)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 10))+
  scale_y_continuous(limits = c(0, 10))+
  theme_minimal() +
  labs(title = "Log of Mean Replicates Comparison, Plus Cells",
       x = "Log of Plus Mean (K27me3)",
       y = "Log of Plus Mean (K4me3)") +
  geom_point(data = filter(data_means_Plus, Label %in% minus_gene_list), color = "red", alpha = 0.5)

  ggsave("PlusCells.MinusPoised.min0.K4me3.vs.K27me3.scatter.pdf")

# Create a scatter plot for the Minus cells
ggplot(data_means_Minus, aes(x = LogMean_K27me3_Minus, y = LogMean_K4me3_Minus)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 10))+
  scale_y_continuous(limits = c(0, 10))+
  theme_minimal() +
  labs(title = "Log10 of Mean Replicates Comparison, Minus Cells",
       x = "Log of Minus Mean (K27me3)",
       y = "Log of Minus Mean (K4me3)") +
  geom_point(alpha = 0.5) +  # Plot points
  geom_point(data = filter(data_means_Minus, Label %in% plus_gene_list), color = "green", alpha = 0.5)

ggsave("MinusCells.PlusPoised.min0.K4me3.vs.K27me3.scatter.pdf")

# Create a scatter plot for the Plus cells
ggplot(data_means_Plus, aes(x = LogMean_K27me3_Plus, y = LogMean_K4me3_Plus)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 10))+
  scale_y_continuous(limits = c(0, 10))+
  theme_minimal() +
  labs(title = "Log of Mean Replicates Comparison, Plus Cells",
       x = "Log of Plus Mean (K27me3)",
       y = "Log of Plus Mean (K4me3)") +
  geom_point(data = filter(data_means_Plus, Label %in% plus_gene_list), color = "green", alpha = 0.5)
  
  ggsave("PlusCells.PlusPoised.min0.K4me3.vs.K27me3.scatter.pdf")
  
  # Create a scatter plot for the Minus cells
  ggplot(data_means_Minus, aes(x = LogMean_K27me3_Minus, y = LogMean_K4me3_Minus)) +
    geom_point(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 10))+
    scale_y_continuous(limits = c(0, 10))+
    theme_minimal() +
    labs(title = "Log10 of Mean Replicates Comparison, Minus Cells",
         x = "Log of Minus Mean (K27me3)",
         y = "Log of Minus Mean (K4me3)") +
    geom_point(alpha = 0.5) +  # Plot points
    geom_point(data = filter(data_means_Minus, Label %in% both_gene_list), color = "yellow", alpha = 0.5)
  
  ggsave("MinusCells.BothPoised.min0.K4me3.vs.K27me3.scatter.pdf")
  
  # Create a scatter plot for the Plus cells
  ggplot(data_means_Plus, aes(x = LogMean_K27me3_Plus, y = LogMean_K4me3_Plus)) +
    geom_point(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 10))+
    scale_y_continuous(limits = c(0, 10))+
    theme_minimal() +
    labs(title = "Log10 of Mean Replicates Comparison, Plus Cells",
         x = "Log of Plus Mean (K27me3)",
         y = "Log of Plus Mean (K4me3)") +
    geom_point(data = filter(data_means_Plus, Label %in% both_gene_list), color = "yellow", alpha = 0.5)
  
  ggsave("PlusCells.BothPoised.min0.K4me3.vs.K27me3.scatter.pdf")
  
  # Create a plot of the ALDH+ poised genes in just the second replicate
  # Create a scatter plot for the Minus cells
  ggplot(data_means_Minus, aes(x = log(K27me3.Minus.R2+1), y = log(K4me3.Minus.R2+1)))+
    geom_point(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 10))+
    scale_y_continuous(limits = c(0, 10))+
    theme_minimal() +
    labs(title = "Log of Replicate2 Comparison, Minus Cells",
         x = "Log of ALDH- readcount (K27me3)",
         y = "Log of ALDH- readcount (K4me3)") +
    geom_point(alpha = 0.5) +  # Plot points
    geom_point(data = filter(data_means_Minus, Label %in% plus_gene_list), color = "green", alpha = 0.5)
  
  ggsave("MinusCells.PlusPoised.rep2.K4me3.vs.K27me3.scatter.pdf")
  
  # Create a scatter plot for the Plus cells
  ggplot(data_means_Plus, aes(x = log(K27me3.Plus.R2+1), y = log(K4me3.Plus.R2+1))) +
    geom_point(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 10))+
    scale_y_continuous(limits = c(0, 10))+
    theme_minimal() +
    labs(title = "Log of Replicate2 Comparison, Plus Cells",
         x = "Log of ALDH+ readcount (K27me3)",
         y = "Log of ALDH+ readcount (K4me3)") +
    geom_point(data = filter(data_means_Plus, Label %in% plus_gene_list), color = "green", alpha = 0.5)
  
  ggsave("PlusCells.PlusPoised.rep2.K4me3.vs.K27me3.scatter.pdf")