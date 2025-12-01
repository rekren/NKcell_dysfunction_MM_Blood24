quant_opt_pcs <- function(seurat_obj) {
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100 # Determine percent of variation associated with each PC
  cumu <- cumsum(pct) # Calculate cumulative percents for each PC
  co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # Determine the difference between variation of PC and subsequent PC
  co2 # last point where change of % of variation is more than 0.1%.
  dims <- min(co1, co2) # Minimum of the two calculation
  dims
  plot_df <- data.frame(pct = pct,
                        cumu = cumu,
                        rank = 1:length(pct)) # Create a dataframe with values
  p <-ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > dims)) +
    geom_text(size=2.5) +
    geom_vline(xintercept = 90, color = "grey") +
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw() +  labs(y = "Per dimension std deviation (%)", x = "Cumulative variation(%)") +
    ElbowPlot(seurat_obj,ndims = 30) # Elbow plot to visualize
  return(p)
}
