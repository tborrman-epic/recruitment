# Generate data for modeling question for technical interview
library(ggplot2)

set.seed(123)  # For reproducibility

mean_val <- 1
sd_val <- 3
meanlog <- log(mean_val^2 / sqrt(sd_val^2 + mean_val^2))
sdlog <- sqrt(log(1 + (sd_val^2 / mean_val^2)))
gene_activation <- rlnorm(n = 500, meanlog = meanlog, sdlog = sdlog)

# Generate ATAC-seq signal to be linearly predictive of gene_activation but with heteroskedasticity
atac_seq_signal <- gene_activation + rnorm(n = 500, mean = 0, sd = 0.02 * (gene_activation**2))
atac_seq_signal <- atac_seq_signal * 10

# Generate H3K27ac, H3K9ac, and H3K4me3 signals with multicolinearity
h3k27ac_signal <- rlnorm(n = 500, meanlog = meanlog, sd = sdlog + 0.2)
h3k9ac_signal <- h3k27ac_signal + rnorm(n = 500, mean = 0, sd = 0.3)
h3k4me3_signal <- h3k27ac_signal + rnorm(n = 500, mean = 0, sd = 0.3)

# Generate other non predictive and non correlated features
h3k27me3_signal <- rlnorm(n = 500, mean = meanlog, sd = sdlog)
distance_to_tss <- runif(n = 500, min = 0, max = 2000)
free_energy <- rnorm(n=500, mean = -50, sd = 5)

# Generate guide IDs
guide_ids <- paste0("guide_", seq_len(500))

# Generate data frame
data <- data.frame(
    guide_ids = guide_ids,
    atac_seq_signal = atac_seq_signal,
    h3k27me3_signal = h3k27me3_signal,
    h3k4me3_signal = h3k4me3_signal,
    h3k27ac_signal = h3k27ac_signal,
    h3k9ac_signal = h3k9ac_signal,
    distance_to_tss = distance_to_tss,
    free_energy = free_energy,
    gene_activation = gene_activation
)

# Save the data to a CSV file
write.csv(data, "guides_data.csv", row.names = FALSE)

# Plot all scatter plots of the data
png("modeling_question_scatter_plots.png", width = 3000, height = 3000, res = 300)
pairs(
    data[2:ncol(data)],
    pch = 20,
    col = rgb(0, 0, 0, alpha = 0.2),
    lower.panel = NULL)
dev.off()

# Easier viewing of the data with log scale
png("modeling_question_log_scale.png", width = 2000, height = 2000, res = 300)
ggplot(
    data = data,
    aes(x = atac_seq_signal, y = gene_activation)
    ) +
    geom_point(alpha = 0.2) +
    theme_minimal() + 
    scale_y_log10() +
    scale_x_log10()
dev.off()

# Scale the features for linear regression
# Note: Scaling is not strictly necessary for linear regression, but it can help with interpretation
data_scaled <- as.data.frame(scale(data[2:ncol(data)]))

# Generate linear regression model from all the features of data
full_model <- lm(gene_activation ~ ., data = data_scaled)
print(summary(full_model))

# Generate linear regression model using most predictive feature 
linear_model <- lm(gene_activation ~ atac_seq_signal, data =  data)

# Plot the data and the linear regression line showing heteroskedasticiy
png("modeling_question_linear_regression.png", width = 2000, height = 2000, res = 300)
ggplot(data=data, aes(x = atac_seq_signal, y = gene_activation)) +
    geom_point(alpha = 0.2) +
    theme_minimal() + 
    geom_abline(
        slope = coef(linear_model)[2],
        intercept = coef(linear_model)[1],
        color = "red"
    )
dev.off()
