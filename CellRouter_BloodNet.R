setwd("~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/bloodnet/")
source("~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/utils_RNA_seq.R")
source('~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/CellRouter_Class.R')
#folder containing .jar files for CellRouter core engine and libraries
libdir <- '/home/edroaldo/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/CellRouter/'

## original data were downloaded from: http://blood.stemcells.cam.ac.uk/data/coordinates_gene_counts_flow_cytometry.txt.gz
##we provide a cleaned file to reproduce these analysis ('expData_gene_symbols_Feb182017.R')
#reviewers can go directly to the line starting with ##START HERE

##Then, we followed the steps below:
data <- read.table('coordinates_gene_counts_flow_cytometry.txt', row.names = 1, sep='\t', header=TRUE, check.names = FALSE)
rownames(data) <- gsub("(LT.)", "LT_", rownames(data))
ndata <- data[, -c(1:26)] #dataset already filtered to include variable genes
ndata <- as.data.frame(t(ndata))

##converting ensembl gene ids to gene symbols
ids <- get(load('ids.R'))
geneTab <- ids[,c('ensembl_gene_id', "external_gene_name")]
geneTab <- geneTab[!duplicated(geneTab$ensembl_gene_id),]
rownames(geneTab) <- as.vector(geneTab$ensembl_gene_id)
colnames(geneTab) <- c('id', 'symbol')
expDat <- averageIds(ndata, geneTab, 'symbol')
m <- data[,c(1:3)]
save(expDat, file='expData_gene_symbols_Feb182017.R')
save(m, file='m.R')

##START HERE
##opening processed files, containing only gene expression data and with gene symbols.
## also load diffusion components coordinates file (m.R)
ndata <- get(load('expData_gene_symbols_Feb182017.R'))
#m contains diffusion components coordinates obtained from the original publication.
m <- get(load('m.R')) 

# Gene regulatory network reconstruction using a correlation-based version of CLR
# published in our previous work with CelNet (Cahan et al, Cell 2014)
library(igraph)
tfs <- find_tfs(species = 'Mm')
grn <- globalGRN(ndata, tfs, 5)
colnames(grn)[1:2]<-c("TG", "TF");
ggrn<- ig_tabToIgraph(grn, simplify = TRUE)
save(ggrn, file='results/GRN.R')


##opening gene regulatory network object
ggrn <- get(load('results/GRN.R'))


### Subpopulation identification and gene signatures with CellRouter
cellrouter <- CellRouter(expdata=as.data.frame(ndata), annotations=colnames(ndata))
cellrouter@rdimension <- as.data.frame(m)
cellrouter <- findsubpopulations(cellrouter, 20, 'jaccard', 'results/kNN_network.gml')
cellrouter <- diffexpr(cellrouter, pvalue = 0.05)
markers <- findmarkers(cellrouter)
write.csv(markers, file='results/Supplementary_Table_3_Gene_Signatures.csv')

## generating plot in 3D
library(rgl)
library(plot3D)
matrix <- cellrouter@rdimension
matrix2 <- matrix[rownames(cellrouter@sampTab),]
colors <- cellrouter@sampTab$colors
plot(matrix2[,1:2],col=colors,pch=19, cex=0.5, xlab='tSNE_1', ylab='tSNE_2', bty ="n")

#To keep a 3D view consistent with the original publication, we rotated and saved
#the 3D coordinates such plots showing gene expression on diffusion components 
#will be the same for each gene, without having to rotate the 3D axis again.
#http://stackoverflow.com/questions/16362381/save-the-orientation-of-a-rgl-plot3d-plot
#coordinates already saved. Do not run the next three lines again.

###plot3d(matrix2[,1], matrix2[,2], matrix2[,3], col=colors, size=3)
###pp <- par3d(no.readonly=TRUE)
###dput(pp, file="results/bloodnetView.R", control = "all")
##rgl.postscript("results/persp3dd_nolabels.pdf","pdf") 

## For some reason, the first pdf plot looks incomplete. However, when opened
#in an image editor, the figure is fine. It could be some problems with handling
#the pre-saved 3D coordinates.
pp <- dget("results/bloodnetView.R")
labels <- nodeLabels(cellrouter@sampTab,'community')
plot3d(matrix2[,1], matrix2[,2], matrix2[,3], col=colors, size=2, xlab = 'DC1', ylab='DC2', zlab='DC3')
text3d(matrix2[,1], matrix2[,2], matrix2[,3], labels, font=5, cex=5)
par3d(pp)
rgl.postscript("results/BloodNet_Subpopulations.pdf","pdf") 

genelist <- c('Gata1', 'Ccl3', 'Ctsg', 'Procr')
plot3DExpression(cellrouter, genelist, FALSE, parameters=pp, filename='results/3D_expression_')

######## Trajectory Detection using CellRouter ###
#k <- findK(cellrouter, 20)
cellrouter <- createKNN(cellrouter, 10, 'jaccard', 'results/paths/kNN_network_trajectory.gml')
filename <- "results/paths/cell_edge_weighted_network.txt"
write.table(cellrouter@graph$edges, file=filename, sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE) #input network

#selecting starting subpopulation
sources <- c('SP_2')
##all remaining subpopulations are targets
targets <- setdiff(as.vector(cellrouter@sampTab$population), sources)
methods <- c("euclidean", "maximum", "manhattan","canberra","binary", 'graph') #graph for distances in KNN
cellrouter <- findpaths(cellrouter, libdir, paste(getwd(), 'results/paths', sep='/'), method="graph")
save(cellrouter, file='results/CellRouter_BloodNet.R')


### loading cellrouter object
cellrouter <- get(load('results/CellRouter_BloodNet.R'))
ranks <- c('path_cost', 'path_flow', 'rank', 'length')
genes2use <- rownames(cellrouter@ndata)
cellrouter <- processtrajectories(cellrouter, genes2use, path.rank=ranks[3], num.cells = 3, neighs = 1)

## perform the analysis using all trajectories
names <- unique(names(cellrouter@pathsinfo$distr))
clusters.show <- names
cellrouter <- correlationpseudotime(cellrouter, type='spearman')
cellrouter <- topgenes(cellrouter, 0.85, 0.15)
cellrouter <- smoothdynamics(cellrouter, names)
cellrouter <- clusterGenesPseudotime(cellrouter, 5)
save(cellrouter, file='results/CellRouter_BloodNet_Processed.R')


## Loading cellrouter object with processed trajectories
cellrouter <- get(load('results/CellRouter_BloodNet_Processed.R'))

transitions <- c('SP_2.SP_8','SP_2.SP_11','SP_2.SP_4') #transitions to be analyzed
p <- 'SP_2.SP_4' #focus on lymphoid differentiation
x <- grnscores(cellrouter, tfs, transitions, direction='both', flip=TRUE, q.up=0.8, q.down=0.2, dir.targets='up', columns=3, 7, 5, 'results/lineage_regulators_score')
m2 <- plottr(cellrouter, p, x$SP_2.SP_4$scores, cluster=FALSE, 2, 8.5, 5, paste('results/', p, 'up_down_diff_dynamics.pdf',sep=''))
plotpaths(cellrouter, transitions, c('Il12rb1', 'Satb1', 'Gata2', 'Klf1', 'Gfi1'), columns = 2, width=8, height = 5, file_prefix = 'results/selected_genes_trajectories')


ids <- get(load('ids.R'))
cc <- get(load('cell_cycle_genes_December_30_2015.R'))
#GO enrichment on selected trajectories
cellrouter <- pathwayenrichment(cellrouter, transitions, 'mouse','org.Mm.eg.db', ids)
pathwaycluster(cellrouter, cellrouter@pathwayenrichment$UP$GOBP, 10, TRUE, 5, 5, 'results/Supplementary_Table_4_GO.pdf')

##characterize lymphoid development
plotClusterHeatmap(cellrouter, transitions, 10, 1.5, 3, 'results/dynamics.pdf')
plotclusters(cellrouter, p, 2, 400, 450, 'results/cluster_dynamics.png')

### enrichment for clustered gene trends
geneList <- list()
clustering <- cellrouter@clusters$SP_2.SP_4$clustering
for(c in unique(clustering)){
  geneList[[paste('cl',c,sep='')]] <- names(clustering[which(clustering == c)])
}

geneList <- lapply(geneList, function(x){convertIDs(ids, x,from='external_gene_name', to="entrezgene")})
geneList <- lapply(geneList, names)
geneList <- lapply(geneList, function(x){x[!is.na(x)]})
ck3 <- compareCluster(geneCluster = geneList, fun = "enrichGO", ont='BP', OrgDb='org.Mm.eg.db', pvalueCutoff = 0.05, readable=T)
write.csv(ck3@compareClusterResult, file='results/Supplementary_Table_5_GO.csv')
pdf(file='results/GO_enrichment_SP_2.SP_4.pdf', width=8, height=7)
plot(ck3, showCategory=10, title='Gene Ontology Biological Processes')
dev.off()


##intermediate subpopulations
transitions <- c('SP_2.SP_9','SP_2.SP_12','SP_2.SP_13')
x <- grnscores(cellrouter, tfs, transitions, direction='up', flip=TRUE, q.up=10, dir.targets='up', columns=3, width=7, height=3, filename='results/lineage_regulators_score_intermediate_up')
p <- 'SP_2.SP_9'
m2 <- plottr(cellrouter, p, x$SP_2.SP_9$scores, cluster=FALSE, 1, 4, 4, paste('results/', p, 'up_down_diff_dynamics.pdf',sep=''))
p <- 'SP_2.SP_12'
m2 <- plottr(cellrouter, p, x$SP_2.SP_12$scores, cluster=FALSE, 1, 4, 4, paste('results/', p, 'up_down_diff_dynamics.pdf',sep=''))
p <- 'SP_2.SP_13'
m2 <- plottr(cellrouter, p, x$SP_2.SP_13$scores, cluster=FALSE, 1, 4, 4, paste('results/', p, 'up_down_diff_dynamics.pdf',sep=''))


