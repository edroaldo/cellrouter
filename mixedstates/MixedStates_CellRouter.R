source('~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/CellRouter_Class.R')
source("~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/utils_RNA_seq.R")
setwd("~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/mixedstates/")
libdir <- '~/Documents/Projects/CellRouter/Manuscript/Submission_1/scripts/CellRouter/'


library(scLVM)
library(tsne)

## The file 'Olsson_RSEM_SingleCellRNASeq.csv' is too big to upload to github.
## The original data can be downloaded from https://www.synapse.org/#!Synapse:syn4975057/files/
## To make easier to reproduce the analysis, we uplodaed a cleaned version of this dataset to github('Olsson_clean_dataset.R')
## Therefore, the analysis can be started from the #START HERE line, below

###guideGenes <- read.csv('ICGS-s5.txt', sep='\t', row.names = 1)
###guideGenes <- guideGenes[-1,]
###ndata <- read.csv('Olsson_RSEM_SingleCellRNASeq.csv', row.names = 1)
###ndata <- ndata[,1:382]
###ndata <- ndata[which(apply(ndata, 1, var) > 0),]
###save(ndata, file='Olsson_clean_dataset.R')
#write.csv(ndata, 'Olsson_RSEM_SingleCellRNASeq_382_single_cells.csv')

#START HERE
guideGenes <- read.csv('ICGS-s5.txt', sep='\t', row.names = 1)
guideGenes <- guideGenes[-1,]
ndata <- get(load('Olsson_clean_dataset.R'))

pca <- prcomp(t(ndata), scale=TRUE, center=TRUE)
loadings <- pca$rotation
num_pc <- 5
quantile <- 0.975
genes2use <- unique(as.vector(unlist(apply(loadings[,1:num_pc], 2, 
                      function(x){print(quantile(x, quantile))
                      names(x[which(abs(x) >= quantile(x, quantile))])}))))

di <- dist(t(ndata[rownames(guideGenes),]))
m <- tsne(di, initial_config=cmdscale(di,k=2)) #rseed=15555
colnames(m) <- c('tSNE1', 'tSNE2')
save(m, file='results/rdimension.R')

######## tSNE annotated by clusters identified by the original publication #######
tmp <- as.data.frame(t(read.csv('ICGS-s5.txt', sep='\t', row.names = 1)))
tmp <- tmp[-1,]
s <- gsub("^[^.]*.", "", rownames(tmp))
s <- gsub("(CD34..LK.)", "LK.", s)
m <- m[s,]
colnames(m) <- c('tSNE1', 'tSNE2')
scores <- as.data.frame(m)
scores$group <- tmp$`column_clusters-flat`
col.rainbow <- rainbow(length(unique(scores$group)))
pdf(file="results/tSNE_mixed_states_clusters.pdf", width=4.2, height=2.5)
p1 <- ggplot(scores,aes(x = tSNE1, y=tSNE2, colour=factor(group))) + geom_point(size=1.5) +
  ggtitle('')  + scale_color_manual(name = "", values=c('yellow', 'darkgoldenrod',
                                                        'orange', 'red','purple',
                                                        'darkblue', 'lightblue',
                                                        'darkgreen','green')) + theme_bw() +
        theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1))

plot(p1)
dev.off()


# Gene regulatory network reconstruction using a correlation-based version of CLR
# published in our previous work with CelNet (Cahan et al, Cell 2014)
tfs <- find_tfs(species = 'Mm')
grn <- globalGRN(ndata[genes2use,], tfs, 5)
colnames(grn)[1:2]<-c("TG", "TF");
ggrn<- ig_tabToIgraph(grn, simplify = TRUE) #simplify changes the weights...
save(ggrn, file='results/GRN.R')

## open GRN object
ggrn <- get(load('results/GRN.R'))

### Subpopulation identification and gene signatures with CellRouter
m <- get(load('results/rdimension.R'))
cellrouter <- CellRouter(expdata=ndata[genes2use,], annotations=colnames(ndata))
cellrouter@rdimension <- as.data.frame(m)
cellrouter <- findsubpopulations(cellrouter, 4, 'jaccard', 'results/kNN_network.gml')
cellrouter <- diffexpr(cellrouter, pvalue = 0.05)
markers <- findmarkers(cellrouter)
write.csv(markers, file='results/Supplementary_Table_3_Gene_Signatures.csv')

genelist <- c('Mecom', 'Flt3', 'Meis1', 'Gata2', 'Mpl', 'Vwf', 'Pf4', 'Itga2b',
              'Epor', 'Klf1', 'Gata1', 'Gfi1b', 'Cebpa', 'Spi1', 'Ctsg', 'Mpo',
              'Elane', 'Slpi', 'Ly86', 'Csf1r', 'Irf8', 'Klf4', 'Gfi1', 'Cebpe', 
              'S100a8', 'Camp', 'Mmp9')


######## Trajectory Detection using CellRouter ###
#k <- findK(cellrouter, 20)
cellrouter <- createKNN(cellrouter, 15, 'jaccard', 'results/paths/kNN_network_trajectory.gml')
filename <- "results/paths/cell_edge_weighted_network.txt"
write.table(cellrouter@graph$edges, file=filename, sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE) #input network

#selecting starting subpopulations
sources <- c('SP_3', 'SP_4', 'SP_2', 'SP_1')
#remaining subpopulations are targets
targets <- setdiff(as.vector(cellrouter@sampTab$population), sources)
methods <- c("euclidean", "maximum", "manhattan","canberra","binary", 'graph') #graph for distances in KNN
cellrouter <- findpaths(cellrouter, libdir, paste(getwd(), 'results/paths', sep='/'), method="graph")
save(cellrouter, file='results/CellRouter_Mixed_States.R')

#load
cellrouter <- get(load('results/CellRouter_Mixed_States.R'))
ranks <- c('path_cost', 'path_flow', 'rank', 'length')
cellrouter <- processtrajectories(cellrouter, genes2use, path.rank=ranks[3], 
                                  num.cells = 3, neighs = 1)
## perform the analysis using all trajectories
## This will take longer because there are more starting subpopulations
names <- unique(names(cellrouter@pathsinfo$distr))
clusters.show <- names
cellrouter <- correlationpseudotime(cellrouter, type='spearman')
cellrouter <- topgenes(cellrouter, 0.85, 0.15)
cellrouter <- smoothdynamics(cellrouter, names)
cellrouter <- clusterGenesPseudotime(cellrouter, 10)
save(cellrouter, file='results/CellRouter_Mixed_States_Processed.R')

##open more recent CellRouter object
cellrouter <- get(load('results/CellRouter_Mixed_States_Processed.R'))

dftfs <- plotheatmapCategory(cellrouter, threshold=6, genelist=genelist, width=9, height=5)

trajectories <- c('SP_4.SP_8', 'SP_4.SP_17', 'SP_3.SP_8', 'SP_3.SP_17',
                  'SP_4.SP_14', 'SP_4.SP_15', 'SP_3.SP_14', 'SP_3.SP_15')
plottrajectories(cellrouter, trajectories, c('Gfi1', 'Irf8'), rescale=FALSE, columns=4, width=10, height=3, filename='results/dynamics_mono_granu_sub_res.pdf')
plotDRExpression(cellrouter, c('Gfi1', 'Irf8'), TRUE, 1, 3, 3, 'results/genes_monogran.pdf')


#### GRN scores for selected transitions ####
### regulators up->targets up
transitions <- c('SP_4.SP_5','SP_4.SP_19','SP_4.SP_15', 'SP_4.SP_8')
p <- 'SP_4.SP_8'
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=15, dir.targets='up', 
               columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
scores <- x[[p]]$scores
rgrn <- induced.subgraph(ggrn, vids=unlist(neighborhood(graph=ggrn,order=1,nodes=names(scores))))
subnet <- ig_convertSmall(rgrn, exponent = 0.5)
subnet <- ig_NiceGraph(cellrouter, subnet, p, vLabels = c('Regulator'))

## this gml files should be imported into a network visualization software such as cytoscape or gephi
write.graph(subnet, file=paste('results/', p, '_network.gml', sep=''), format='gml')
m2 <- plottr(cellrouter, p, x$SP_4.SP_8$scores, cluster=TRUE, 2, 8, 2.5, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))

p <- 'SP_4.SP_5'
m2 <- plottr(cellrouter, p, x$SP_4.SP_5$scores, cluster=TRUE, 2, 8, 3, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))
p <- 'SP_4.SP_19'
m2 <- plottr(cellrouter, p, x$SP_4.SP_19$scores, cluster=TRUE, 2, 8, 3, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))
p <- 'SP_4.SP_15'
m2 <- plottr(cellrouter, p, x$SP_4.SP_15$scores, cluster=TRUE, 2, 8, 3, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))

paths <- c('SP_4.SP_8', 'SP_3.SP_8')
plotpaths(cellrouter, paths, c('Gfi1', 'S100a8', 'Cebpe'), columns = 1, width=4, height = 6, file_prefix = 'results/selected_genes_granulocytes')

#paths <- c('SP_4.SP_8', 'SP_3.SP_8', 'SP_4.SP_15')
#plottrajectories(cellrouter, paths, c('Gfi1', 'S100a8', 'Cebpe'), rescale = TRUE, columns=1, width=3, height=4.5, filename='results/dynamics_granu.pdf')
#plotDRExpression(cellrouter, c('Gfi1', 'S100a8', 'Cebpe'), TRUE, 1, 3, 5, 'results/genes_panelE_tSNE.pdf')

ids <- get(load('ids.R'))
cc <- get(load('cell_cycle_genes_December_30_2015.R')) #remove cell cyle genes

#GO enrichment on selected trajectories
paths <- c('SP_4.SP_19', 'SP_4.SP_5', 'SP_4.SP_8', 'SP_4.SP_17', 'SP_4.SP_15')
cellrouter <- pathwayenrichment(cellrouter, paths, 'mouse','org.Mm.eg.db', ids)
pathwaycluster(cellrouter, cellrouter@pathwayenrichment$UP$GOBP, 10, TRUE, 5, 5, 'results/Supplementary_Table_4_GO.pdf')
