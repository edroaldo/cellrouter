#--> OK. Analysis reproductible!
setwd("~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/mesoderm/")
source("~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/utils_RNA_seq.R")
source('~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/CellRouter_Class.R')
libdir <- '~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/CellRouter/'

#library('scLVM')
library('Rtsne')
library('DESeq2')

#Uncompress the file counts.gz inside the folder mesoderm, using: tar -xvzf counts.gz 
## The file counts.gz can also be downloaded from http://gastrulation.stemcells.cam.ac.uk/scialdone2016
#to obtain the file 'counts.txt'
# Folow the instructions below to reproduce the analysis

## Preprocessing steps
counts <- read.table('counts.txt', sep="", check.names = FALSE)
samples <- read.table('original_metadata.txt', sep="", header=TRUE, row.names=1)
dds <- DESeqDataSetFromMatrix(countData =  counts,
                              colData = samples,
                              design = ~ embryoStage)
dds <- estimateSizeFactors(dds)
nCounts <- counts(dds, normalized=TRUE)

### map to gene symbols
ids <- get(load('ids.R'))
geneTab <- ids[,c('ensembl_gene_id', "external_gene_name")]
geneTab <- geneTab[!duplicated(geneTab$ensembl_gene_id),]
rownames(geneTab) <- as.vector(geneTab$ensembl_gene_id)
colnames(geneTab) <- c('id', 'symbol')
expDat <- averageIds(nCounts, geneTab, 'symbol')
save(expDat, file='norm_counts_geneSymbols.R')

#### convert list of high variable genes from ensembl ids to gene symbols
ids <- get(load('ids.R'))
load("high_var_genes.RData")#HVG across all cells
symbols <- convertIDs(ids, high.var.genes.all, from='ensembl_gene_id', to="external_gene_name")
symbols <- intersect(names(symbols), rownames(expDat))
genes<-symbols
save(genes, file='variable_genes.R') #we also provide this gene list

nCounts <- get(load('norm_counts_geneSymbols.R'))
#var <- apply(nCounts, 1, var)
#nCounts <- nCounts[which(var > 0),]

### tSNE dimenionality reduction using most variable genes, as performed in original paper
cor.mat<-cor(nCounts[genes,],method="spearman")
dissim<-(1-cor.mat)/2
dist.mat<-as.dist(dissim)
set.seed(0)
tsne.seed.0<-Rtsne(dist.mat, is_distance=TRUE)

#plot tSNE coordinates with cells colored by clusters identified in the orignal paper
plot(pch=20,tsne.seed.0$Y[,1],tsne.seed.0$Y[,2], col=as.character(samples[colnames(nCounts),"cluster"]))
m <- as.data.frame(tsne.seed.0$Y[,1:2])
rownames(m) <- colnames(nCounts)
colnames(m) <- c('tSNE1', 'tSNE2')
save(m, file='rdimension.R') #we provide this file, to keep consistent with our paper


###### Opening files to start the analysis ###
#nCounts <- get(load('norm_counts_geneSymbols.R'))
#m <- get(load('rdimension.R'))
nCounts <- nCounts[which(apply(nCounts, 1, var) > 0),]
nCounts <- nCounts[rowMeans(nCounts) > 1,]
pca <- prcomp(t(nCounts), scale=TRUE, center=TRUE)
loadings <- pca$rotation
num_pc <- 3
quantile <- 0.95
genes2use <- unique(as.vector(unlist(apply(loadings[,1:num_pc], 2, function(x){names(x[which(abs(x) >= quantile(x, quantile))])}))))
genes2use <- unique(c(genes2use, get(load('variable_genes.R'))))
save(genes2use, file='genes2use.R')

# This files were generated from the previous code. 
nCounts <- get(load('norm_counts_geneSymbols.R'))
genes2use <- get(load('genes2use.R'))
m <- get(load('rdimension.R'))

# Gene regulatory network reconstruction using a correlation-based version of CLR
# published in our previous work with CelNet (Cahan et al, Cell 2014)
tfs <- find_tfs(species = 'Mm')
grn <- globalGRN(nCounts[genes2use,], tfs, 5)
colnames(grn)[1:2]<-c("TG", "TF");
ggrn<- ig_tabToIgraph(grn, simplify = TRUE) #simplify changes the weights...
save(ggrn, file='results/GRN.R')
#####

ggrn <- get(load('results/GRN.R'))

### Subpopulation identification and gene signatures with CellRouter
cellrouter <- CellRouter(expdata=as.data.frame(nCounts), annotations=colnames(nCounts))
cellrouter@rdimension <- as.data.frame(m)
cellrouter <- findsubpopulations(cellrouter, 5, 'jaccard', 'results/kNN_network.gml')
cellrouter <- diffexpr(cellrouter, pvalue = 0.05)
markers <- findmarkers(cellrouter)
write.csv(markers, file='results/Supplementary_Table_5_Gene_Signatures.csv')

######## Trajectory Detection using CellRouter ###
#k <- findK(cellrouter, 20)
cellrouter <- createKNN(cellrouter, 15, 'jaccard', 'results/paths/kNN_network_trajectory.gml')
filename <- "results/paths/cell_edge_weighted_network.txt"
write.table(cellrouter@graph$edges, file=filename, sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE) #input network

#selecting starting subpopulations
sources <- c('SP_11', 'SP_6', 'SP_12') #mesoderm dynamics and endothelium to blood dynamics
#remaining subpopulations will be targets
targets <- setdiff(as.vector(cellrouter@sampTab$population), sources)
methods <- c("euclidean", "maximum", "manhattan","canberra","binary", 'graph') #graph for distances in KNN
cellrouter <- findpaths(cellrouter, libdir, paste(getwd(), 'results/paths', sep='/'), method="graph")
save(cellrouter, file='results/CellRouter_Mesoderm.R')

cellrouter <- get(load('results/CellRouter_Mesoderm.R'))
ranks <- c('path_cost', 'path_flow', 'rank', 'length')
cellrouter@ndata <- log2(cellrouter@ndata + 1)
cellrouter <- processtrajectories(cellrouter, genes2use, path.rank=ranks[3], 
                                  num.cells = 3, neighs = 1)

## perform the analysis using all trajectories
names <- unique(names(cellrouter@pathsinfo$distr))
clusters.show <- names
cellrouter <- correlationpseudotime(cellrouter, type='spearman')
cellrouter <- topgenes(cellrouter, 0.85, 0.15)
cellrouter <- smoothdynamics(cellrouter, names)
cellrouter <- clusterGenesPseudotime(cellrouter, 5)
save(cellrouter, file='results/CellRouter_Mesoderm_Processed.R')

## Loading cellrouter object with processed trajectories
cellrouter <- get(load('results/CellRouter_Mesoderm_Processed.R'))

transitions <- c('SP_6.SP_15', 'SP_11.SP_15', 'SP_12.SP_4', 'SP_12.SP_10')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=15, dir.targets='up', 
               columns=2, width=7, height=4, flip=FALSE, filename='results/lineage_regulators_score_up')
p <- 'SP_6.SP_15'
m2 <- plottr(cellrouter, p, x[[p]]$scores, cluster=TRUE, 1, 5, 5.5, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))
plottrajectories(cellrouter, p, c('Lefty2','Tbx6','Bmp4', 'Podxl'), rescale = TRUE, columns=1, width=4,height=1.5, filename='results/selected_genes_mesoderm_SP_6.SP_15.pdf')

p <- 'SP_12.SP_4'
m2 <- plottr(cellrouter, p, x[[p]]$scores, cluster=TRUE, 2, 8, 2.5, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))

plotbranch(cellrouter, 'up','SP_12.SP_4', 'SP_6.SP_15', 1, width=3, height=4, filename='results/SP_12.SP_4_up_branch_dynamics.pdf')
plotbranch(cellrouter, 'down','SP_12.SP_4', 'SP_6.SP_15', 1, width=3, height=4, filename='results/SP_12.SP_4_down_branch_dynamics.pdf')

##Enrichment analysis for genes dinamically regulated during erythroid specification
ids <- get(load('ids.R'))
geneList <- list(up=names(cellrouter@top.correlations$up$SP_12.SP_4),
                 down=names(cellrouter@top.correlations$down$SP_12.SP_4))
geneList <- lapply(geneList, function(x){convertIDs(ids, x,from='external_gene_name', to="entrezgene")})
geneList <- lapply(geneList, names)
geneList <- lapply(geneList, function(x){x[!is.na(x)]})
ck3 <- compareCluster(geneCluster = geneList, fun = "enrichGO", ont='BP', OrgDb='org.Mm.eg.db', pvalueCutoff = 0.05)
#write.csv(ck3@compareClusterResult, file='results/GO_enrichment_SP_12.SP_4.csv')
write.csv(ck3@compareClusterResult, file='results/Supplementary_Table_7_GO.csv')
pdf(file='results/GO_enrichment_SP_12.SP_4.pdf', width=6, height=5)
plot(ck3, showCategory=10, title='Gene Ontology Biological Processes')
dev.off()

### GO enrichment for genes regulated during mersoderm development
geneList <- list(up=names(cellrouter@top.correlations$up$SP_6.SP_15),
                 down=names(cellrouter@top.correlations$down$SP_6.SP_15))
geneList <- lapply(geneList, function(x){convertIDs(ids, x,from='external_gene_name', to="entrezgene")})
geneList <- lapply(geneList, names)
geneList <- lapply(geneList, function(x){x[!is.na(x)]})
ck3 <- compareCluster(geneCluster = geneList, fun = "enrichGO", ont='BP', OrgDb='org.Mm.eg.db', pvalueCutoff = 0.05)
#write.csv(ck3@compareClusterResult, file='results/GO_enrichment_SP_6.SP_15.csv')
#pdf(file='results/GO_enrichment_SP_6.SP_15.pdf', width=8, height=5)
write.csv(ck3@compareClusterResult, file='results/Supplementary_Table_6_GO.csv')
pdf(file='results/GO_enrichment_SP_6.SP_15.pdf', width=8, height=5)
plot(ck3, showCategory=10, title='Gene Ontology Biological Processes')
dev.off()