library(WGCNA)

gnet <- readRDS("data/genes/Networks/Network_DS2_MM20.rds")
inet <- readRDS("data/isoforms/Networks/Network_DS2_MM20.rds")
pdf("data/figures/network_dendrograms_gene_iso.pdf", width = 4, height = 3)
plotDendroAndColors(
    dendro = gnet$dendrograms[[1]], 
    colors = labels2colors(gnet$colors), 
    "Modules", 
    dendroLabels = FALSE, 
    hang = 0.03, 
    addGuide = TRUE, 
    guideHang = 0.05, 
)
plotDendroAndColors(
    dendro = inet$dendrograms[[1]], 
    colors = labels2colors(inet$colors[inet$blockGenes[[1]]]), 
    "Modules", 
    dendroLabels = FALSE, 
    hang = 0.03, 
    addGuide = TRUE, 
    guideHang = 0.05
)
plotDendroAndColors(
    dendro = inet$dendrograms[[2]], 
    colors = labels2colors(inet$colors[inet$blockGenes[[2]]]), 
    "Modules", 
    dendroLabels = FALSE, 
    hang = 0.03, 
    addGuide = TRUE, 
    guideHang = 0.05
)
plotDendroAndColors(
    dendro = inet$dendrograms[[3]], 
    colors = labels2colors(inet$colors[inet$blockGenes[[3]]]), 
    "Modules", 
    dendroLabels = FALSE, 
    hang = 0.03, 
    addGuide = TRUE, 
    guideHang = 0.05
)
dev.off()
