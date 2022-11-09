
nps1.sorted <- subset_samples(nps1, Sorted == "Postsort") 
ord.sorted <- ordinate(nps1.sorted, method = "RDA", distance  = "euclidean")

pdf("pca_plot_5n.pdf", 6.5, 3)
ggp <- plot_ordination(nps1.sorted, ordination = ord.sorted, type = "samples", color = "Disease.status", shape = "Sort") + theme_bw()
ggp <- ggp + geom_point(size = 2) + facet_grid( ~ Age) +scale_color_manual(values=c( "#00BFC4","#F8766D"))+ stat_ellipse()
plot(ggp)
dev.off()
