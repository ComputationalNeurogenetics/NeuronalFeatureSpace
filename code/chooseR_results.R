
# Use mclapply based loop in chooseR_server_side.Rmd to calculate results per resolution value by using variables below.
s.data <- qread("../scATAC_data/s.data.merged.for.chooseR.qs", nthreads = 20)
npcs <- 58
resolutions <- c(seq(1,20,by=1))
assay <- "peaks"
reduction <- "lsi"
results_path <- "../chooseR/results_leiden/"


source("../chooseR/R/pipeline.R")
# Create silhouette plot
# Read in scores and calculate CIs
scores <- purrr::map(
  paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
  readRDS
)
scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()

meds <- scores %>%
  dplyr::group_by(res) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters)
  ) %>%
  tidyr::unnest_wider(boot)

counts <- purrr::map(
  paste0(results_path, "cells_per_cluster__", resolutions, ".rds"),
  readRDS
)
meds$min.cell.count <- sapply(counts,function(x){summary(as.numeric(x))[1]})
meds$first.quarter.cell.count <- sapply(counts,function(x){summary(as.numeric(x))[2]})
meds$mean.cell.count <- sapply(counts,function(x){summary(as.numeric(x))[4]})
meds$median.cell.count <- sapply(counts,function(x){summary(as.numeric(x))[3]})
meds$third.quarter.cell.count <- sapply(counts,function(x){summary(as.numeric(x))[5]})
meds$max.cell.count <- sapply(counts,function(x){summary(as.numeric(x))[6]})

writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))

# Plot cluster size statistical summary values
ggplot(meds, aes(x=1:38)) + geom_line(aes(y=min.cell.count, colour="Min")) + geom_line(aes(y=max.cell.count, colour="Max")) + geom_line(aes(y=median.cell.count, colour="Median")) + geom_line(aes(y=mean.cell.count, colour="Mean")) + geom_line(aes(y=first.quarter.cell.count, colour="1st quarter")) + geom_line(aes(y=third.quarter.cell.count, colour="3rd quarter")) + ylab("Cell count") + ggtitle("Number of cells per cluster per resolution value - statistical summary values") + theme_minimal() + scale_x_continuous(breaks = seq(1, 20, by = 1)) + scale_y_continuous(breaks = seq(0, 600, by = 50)) + theme(text=element_text(size=16)) + xlab("Resolution")

# Find thresholds
threshold <- max(meds$low_med)
choice <- as.character(
  meds %>%
    dplyr::filter(med >= threshold) %>%
    dplyr::arrange(n_clusters) %>%
    tail(n = 1) %>%
    dplyr::pull(res)
)

# And plot!
ggplot(meds, aes(factor(res), med)) +
  geom_crossbar(
    aes(ymin = low_med, ymax = high_med),
    fill = "grey",
    size = 0.25
  ) +
  geom_hline(aes(yintercept = threshold), colour = "blue") +
  geom_vline(aes(xintercept = choice), colour = "red") +
  geom_jitter(
    data = scores,
    aes(factor(res), avg_sil),
    size = 0.35,
    width = 0.15
  ) +
  scale_x_discrete("Resolution") +
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

ggsave(
  filename = paste0(results_path, "silhouette_distribution_plot.png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

# Finally, a dot plot of silhouette scores to help identify less robust clusters
# The initial pipe is to order the clusters by silhouette score
scores %>%
  dplyr::filter(res == choice) %>%
  dplyr::arrange(dplyr::desc(avg_sil)) %>%
  dplyr::mutate_at("cluster", ordered, levels = .$cluster) %>%
  ggplot(aes(factor(cluster), avg_sil)) +
  geom_point() +
  scale_x_discrete("Cluster") +
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_grid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

ggsave(
  filename = paste0(results_path, "silhouette_point_plot_", choice, ".png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

scores.out <- scores %>%
  dplyr::filter(res == choice) %>%
  dplyr::arrange(dplyr::desc(avg_sil)) %>%
  dplyr::mutate_at("cluster", ordered, levels = .$cluster)

saveRDS(scores.out,file="scores.Rds")

