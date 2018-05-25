##### CellRouter Class #####
#install.packages("devtools")
#devtools::install_github('hadley/devtools')
#devtools::install_github("klutometis/roxygen")
#library(roxygen2)

suppressWarnings(suppressMessages(require('reshape')))
suppressWarnings(suppressMessages(require('reshape2')))
suppressWarnings(suppressMessages(require('pheatmap')))
suppressWarnings(suppressMessages(library('scales')))
#library('geomnet')
suppressWarnings(suppressMessages(require('ggplot2')))
suppressWarnings(suppressMessages(require('grid')))

CellRouter <- setClass("CellRouter", slots=
                        c(rawdata="data.frame", ndata="data.frame", scale.data="data.frame",
                          sampTab="data.frame", rdimension="data.frame", pca="list", tsne="list", dc="list", dr.custom="list",
                          var.genes="character",graph="list",signatures="list", sources="character", targets="character",
                          directory="list", paths="data.frame", networks="list",
                          genes.trajectory="character", pathsinfo="list",
                          dynamics="list", clusters="list", correlation="list",
                          top.correlations="list", pathwayenrichment="list"))

#' Check CellRouter object
#' @param object CellRouter object
setValidity("CellRouter",
            function(object){
              msg <- NULL
              if(!is.data.frame(object@rawdata)){
                msg <- c(msg, "expression data must be a data.frame")
              }else if(sum(apply(is.na(object@rawdata), 1, sum) > 0)){
                msn <- c(msg, "expression data must not have NAs")
              }
              if(is.null(msg)){
                TRUE
              }else{
                msg
              }
            })

#' Scale and center the data. Individually regress variables provided in vars.regress using using a linear model. Other models are under development.
#' @param object CellRouter object
#' @param genes.use Vector of genes to scale/center. Default is all genes in object@@ndata
#' @param scale.max Max value in scaled data. Default id 10
#' @param vars.regress Variables to regress out
#'
#' @return CellRouter object
#' @export
setGeneric("scaleData", function(object, genes.use=NULL, scale.max=10, vars.regress=NULL) standardGeneric("scaleData"))
setMethod("scaleData",
          signature="CellRouter",
          definition=function(object, genes.use, scale.max, vars.regress){

            if(is.null(genes.use)){
              #genes.use <- object@var.genes
              genes.use <- rownames(object@ndata)
            }
            data.use <- object@ndata[genes.use,]

            if(!is.null(vars.regress)){
              print('Regression...')
              cat(vars.regress, '\n')
              #data.use <- object@ndata[genes.use, , drop = FALSE];
              #gene.expr <- as.matrix(x = data.use[genes.use, , drop = FALSE])
              gene.expr <- as.matrix(object@ndata[genes.use, , drop = FALSE])
              latent.data <- as.data.frame(object@sampTab[,vars.regress])
              rownames(latent.data) <- rownames(object@sampTab)
              colnames(latent.data) <- vars.regress

              new.data <- sapply(X=genes.use, FUN=function(x){
                regression.mat <- cbind(latent.data, gene.expr[x,])#cbind(latent.data, gene.expr[x,])
                colnames(regression.mat) <- c(colnames(latent.data), "GENE")
                fmla <- as.formula(
                  object = paste0(
                    "GENE ",
                    " ~ ",
                    paste(vars.regress, collapse = "+")
                  )
                )
                lm(formula = fmla, data = regression.mat)$residuals
              })
              data.use <- t(new.data)
            }
            #scale.data <- t(scale(t(object@ndata[genes.use,])))
            scale.data <- t(scale(t(data.use)))
            scale.data[is.na(scale.data)] <- 0
            scale.data[which(scale.data > scale.max)] <- scale.max
            object@scale.data <- as.data.frame(scale.data)

            gc(verbose = FALSE)
            object
          }
)
#' Principal component analysis
#' @param object CelLRouter object
#' @param num.pcs Number of principlal components to compute
#' @param genes.use Genes used for principlam component analysis. Default is all genes in object@@ndata
#' @param seed seed
#' @export

setGeneric("computePCA", function(object, num.pcs, genes.use=NULL, seed=1) standardGeneric("computePCA"))
setMethod("computePCA",
          signature="CellRouter",
          definition=function(object, num.pcs, genes.use, seed){
            #computePCA <- function(data, num.pcs, seed=7){
            library(irlba)
            if (!is.null(seed)) {
              set.seed(seed = seed)
            }
            if(is.null(genes.use)){
              #genes.use <- object@var.genes
              genes.use <- rownames(object@ndata)
            }
            pca <- irlba(A = t(object@scale.data[genes.use,]), nv = num.pcs)
            gene.loadings <- pca$v
            #cell.embeddings <- pca$u
            cell.embeddings <- pca$u %*% diag(pca$d)
            sdev <- pca$d/sqrt(max(1, ncol(data) - 1))
            rownames(gene.loadings) <- rownames(object@scale.data)
            colnames(gene.loadings) <- paste0('PC', 1:num.pcs)
            rownames(cell.embeddings) <- colnames(object@scale.data)
            colnames(cell.embeddings) <- colnames(gene.loadings)

            object@pca <- list(gene.loadings = gene.loadings, cell.embeddings=cell.embeddings, sdev=sdev)
            object@rdimension <- as.data.frame(object@pca$cell.embeddings)

            object
          }
)
#' Perform dimensionality reduction using t-SNE
#' @param object do something
#' @param num.pcs Number of principal components used for dimensionality reduction using t-SNE
#' @param perplexity Perplexity
#' @param max_iter, Max number of iterations
#' @param seed seed
#' @import Rtsne
#' @export
setGeneric("computeTSNE", function(object, num.pcs, perplexity=30, max_iter=2000, seed=7) standardGeneric("computeTSNE"))
setMethod("computeTSNE",
          signature="CellRouter",
          definition=function(object, num.pcs, perplexity, max_iter, seed){

            #computeTSNE <- function(pca, num.pcs, perplexity=40, max_iter=2000, seed=7){
            library(Rtsne)
            if (!is.null(seed)) {
              set.seed(seed = seed)
            }

            #tsne.done <- Rtsne(pca$cell.embeddings[,1:num.pcs], perplexity = perplexity, max_iter = max_iter) #implement this type of tsne analysis in cellrouter, using PCA....
            pca <- object@pca
            #tsne.done <- Rtsne(pca$cell.embeddings[,1:num.pcs], max_iter = max_iter) #implement this type of tsne analysis in cellrouter, using PCA....
            tsne.done <- Rtsne(pca$cell.embeddings[,1:num.pcs], perplexity=perplexity, max_iter = max_iter) #implement this type of tsne analysis in cellrouter, using PCA....
            #tsne.done <- Rtsne(pca$cell.embeddings[,1:num.pcs])
            #plot(tsne.done$Y, xlab= 't-SNE 1',ylab= 't-SNE 2', pch=20)
            m <- tsne.done$Y
            rownames(m) <- rownames(pca$cell.embeddings)
            colnames(m) <- c('tSNE 1', 'tSNE 2')
            object@tsne <- list(cell.embeddings=m)

            object
          }
)
#' Dimensionality reduction using diffusion components
#' @param object CellRouter objecr
#' @param genes.use Genes used for dimensionality reduction
#' @param k Parameter k to be used by the DiffusionMap function from destiny pakcage. default is 20
#' @param sigma Parameter sigma to be used by the DiffusionMap function from destiny pakcage. default is local
#' @param seed seed
#' @export

computeDC <- function(object, genes.use=NULL, k=20,sigma='local', seed=1){
  library(destiny) #anyoing error with DLLs all the time...
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  if(is.null(genes.use)){
    genes.use <- object@var.genes
  }

  pca <- object@pca
  diff.comp <- DiffusionMap(as.matrix(t(object@scale.data[genes.use,])), k=20, sigma='local')
  dc <- eigenvectors(diff.comp)
  rownames(dc) <- colnames(object@scale.data)
  object@dc <- list(cell.embeddings=dc)

  object
}

setGeneric("customSpace", function(object, matrix) standardGeneric("customSpace"))
setMethod("customSpace",
          signature="CellRouter",
          definition=function(object, matrix){
            object@dr.custom <- list(cell.embeddings=matrix)
            return(object)
          }
)


#' Normalize the data
#' @param object CellRouter object
#' @export
setGeneric("Normalize", function(object) standardGeneric("Normalize"))
setMethod("Normalize",
          signature="CellRouter",
          definition=function(object){

            x <- object@ndata
            x <- t(t(x)/apply(x,2,sum))
            x <- log1p(x * 10000)
            object@ndata <- as.data.frame(x)
            return(object)
          }
)

#' Identify clusters baed on graph-clustering or model based clustering
#' @param object CellRouter object
#' @param method Method: graph-based clustering or model-based clustering
#' @param k number of nearest neighbors to build a k-nearest neighbors graph
#' @param num.pcs number of principal components that will define the space from where the kNN graph is identified. For example, if num.pcs=10, the kNN graph will be created from a 10-dimensional PCA space
#' @param sim.type Updates the kNN graph to encode cell-cell similarities. Only the Jaccard similarity is implemented in the current version
#' @export

setGeneric("findClusters", function(object, method='graph.clustering', k=20, num.pcs=20, sim.type='jaccard') standardGeneric("findClusters"))
setMethod("findClusters",
          signature = "CellRouter",
          definition = function(object, method, k, num.pcs, sim.type){
            if(method=='graph.clustering'){
              cat('Graph-based clustering\n')
              cat('k: ', k, '\n')
              cat('similarity type: ', sim.type, '\n')
              cat('number of principal components: ', num.pcs, '\n')
              object <- graphClustering(object, k=k, num.pcs=num.pcs, sim.type)
            }else if(method=='model.clustering'){
              cat('Model-based clustering\n')
              cat('number of principal components: ', num.pcs)
              object <- modelClustering(object, num.pcs=num.pcs)
            }
            object
          }
)
#' Model-based clustering using the Mclust package
#' @param object CellRouter object
#' @param num.pcs number of principal components that will define the space used as input to perform model-based clustering
#' @export
setGeneric("modelClustering", function(object, num.pcs) standardGeneric("modelClustering"))
setMethod("modelClustering",
          signature = "CellRouter",
          definition = function(object, num.pcs){
              library(mclust)
              sampTab <- object@sampTab
              colname <- 'population'
              matrix <- object@pca$cell.embeddings[,1:num.pcs]
              mclust <- Mclust(matrix)
              sampTab[names(mclust$classification),colname] <- as.character(mclust$classification)
              sampTab <- sampTab[order(as.numeric(sampTab[[colname]])),]

              colors <- cRampClust(1:length(unique(sampTab[[colname]])), 8)
              names(colors) <- unique(sampTab[[colname]])

              replicate_row <- as.vector(unlist(lapply(split(sampTab, sampTab[[colname]]), nrow)))
              colors_row <- rep(colors, times=replicate_row)
              sampTab[,'colors'] <- colors_row

              object@sampTab <- sampTab
              object
            }
          )

#' Graph-based clustering
#' @param object CellRouter object
#' @param k number of nearest neighbors to build a k-nearest neighbors graph
#' @param num.pcs number of principal components that will define the space from where the kNN graph is identified. For example, if num.pcs=10, the kNN graph will be created from a 10-dimensional PCA space
#' @param sim.type Updates the kNN graph to encode cell-cell similarities. Only the Jaccard similarity is implemented in the current version
#' @param filename Save .gml file containing the kNN graph
#' @export

setGeneric("graphClustering", function(object, k=5, num.pcs, sim.type="jaccard",
                                          filename="graph_subpopulations.gml") standardGeneric("graphClustering"))
setMethod("graphClustering",
          signature = "CellRouter",
          definition = function(object, k, num.pcs, sim.type, filename){
            library('cccd')
            library('proxy') # Library of similarity/dissimilarity measures for 'dist()'
            #matrix <- object@rdimension
            sampTab <- object@sampTab
            matrix <- object@pca$cell.embeddings[,1:num.pcs]

            print('building k-nearest neighbors graph')
            dm <- as.matrix(dist(matrix))
            h <- nng(dx=dm,k=k)
            if(sim.type == 'jaccard'){
              sim <- similarity.jaccard(h, vids=V(h), loops=FALSE)
            }else if(sim.type == 'invlog'){
              sim <- similarity.invlogweighted(h, vids=V(h))
            }
            el <- get.edgelist(h)
            weights <- sim[el]

            E(h)$weight  <- weights
            V(h)$name <- rownames(matrix)
            edges <- as.data.frame(get.edgelist(h))
            rownames(edges) <- paste(edges$V1, edges$V2, sep='_')
            edges$weight <- as.numeric(E(h)$weight)


            ## Community detection to discover subpopulation structure
            print('discoverying subpopulation structure')
            comms <- multilevel.community(as.undirected(h), weights = E(h)$weight)
            V(h)$comms <- membership(comms)
            cell.comms <- commToNames(comms, '') #SP means SubPopulations
            allcells <- as.vector(unlist(cell.comms))

            ## Making sure that color mappings are correct
            sampTab <- sampTab[allcells,] #changesorder of cells in the table
            sampTab$population <- ''
            sampTab$colors <- ''
            comm.colors <- cRampClust(unique(membership(comms)), 8)
            #comm.colors <- cRampClust(unique(membership(comms)), length(unique(membership(comms))))
            names(comm.colors) <- names(cell.comms)
            for(comm in names(cell.comms)){
              sampTab[cell.comms[[comm]], 'population'] <- comm
              sampTab[cell.comms[[comm]], 'colors'] <- comm.colors[comm]
            }
            #sampTab$community <- as.vector(unlist(lapply(strsplit(sampTab$population, split="_"), "[", 2)))
            sampTab$community <- as.vector(sampTab$population)

            ## mapping information to the igraph object
            V(h)[rownames(sampTab)]$subpopulation <- sampTab$colors
            V(h)[rownames(sampTab)]$colors <- sampTab$colors
            V(h)[names(nodeLabels(sampTab,'community'))]$label <- nodeLabels(sampTab, 'community')
            V(h)$size <- 5
            E(h)$arrow.size <- 0.01
            colors <- rainbow(max(membership(comms)))

            #print("plotting graph in RStudio")
            #plot(h,vertex.color=V(h)$colors, vertex.frame.color=V(h)$colors, layout=as.matrix(matrix[,dim.plot]))
            #print('done plotting graph')

            ## Useful information about the graph
            graph <- list()
            graph[['network']] <- h
            graph[['edges']] <- edges
            graph[['similarity_matrix']] <- sim
            graph[['subpopulation']] <- cell.comms
            graph[['communities']] <- comms

            print('updating CellRouter object')
            object@graph <- graph
            if(nrow(object@rawdata) > 0){
             object@rawdata <- object@rawdata[,rownames(sampTab)]
            }
            object@ndata <- object@ndata[,rownames(sampTab)]
            object@sampTab <- sampTab

            write.graph(graph = h, file = filename, format = 'gml')

            rm(h)
            rm(edges)
            rm(sim)

            return(object)
          }
)


#https://briatte.github.io/ggnetwork/

#' Plot kNN graph`
#' @param object CellRouter object
#' @param reduction.type The reduced dimension space used to visualize the kNN graph: tsne, pca, dc or custom
#' @param column.ann column in the metadata table used to annotate the kNN graph. For example, clusters, sorted cell populations
#' @param column.color column in the metadata table corresponding to color used to annotate the kNN graph. Should correspond to the metadata in column.ann
#' @param dims.use dimensions to plot
#' @param width width of output file
#' @param height height og outpur file
#' @param filename name of pdf file generated
#' @import ggplot2
#' @export

setGeneric("plotKNN", function(object, reduction.type="tsne", column.ann, column.color,
                               dims.use=c(1,2), width=10, height=10, filename='kNN_graph.pdf') standardGeneric("plotKNN"))
setMethod("plotKNN",
          signature = "CellRouter",
          definition = function(object, reduction.type, column.ann, column.color,
                                dims.use, width, height, filename){
            library(ggnetwork)
            h <- object@graph$network
            matrix <- slot(object, reduction.type)$cell.embeddings[V(h)$name,dims.use]
            colors <- unique(object@sampTab[[column.color]])
            names(colors) <- unique(as.vector(object@sampTab[[column.ann]]))

            g <- ggnetwork(h, layout=as.matrix(matrix), na.rm=TRUE)
            pdf(file=filename, width=width, height=height)
            g2 <- ggplot(g, aes(x = x, y = y, xend = xend, yend = yend)) +
              geom_edges(alpha=0.3) +
              geom_nodes(aes(color = factor(comms)), size=1) +
              theme_blank() + scale_color_manual("", values=colors) +
              guides(col=guide_legend(direction="vertical", keywidth = 0.75, keyheight = 0.85, override.aes = list(size=3)))
            plot(g2)
            dev.off()
            plot(g2)
          }
)

#' Build kNN graph
#' @param object CellRouter object
#' @param k number of nearest neighbors to build a k-nearest neighbors graph for trajectory reconstruction
#' @param column.ann Column in the metadata table specifying the transitions to be identified. For example, if 'population' is provided, transitions will be identified between clusters previously identified. However, sorted cell populations or customized states can also be used. Check our tutorials for detailed examples.
#' @param num.pcs number of principal components that will define the space from where the kNN graph is identified. For example, if num.pcs=10, the kNN graph will be created from a 10-dimensional PCA space
#' @param sim.type Updates the kNN graph to encode cell-cell similarities. Only the Jaccard similarity is implemented in the current version
#' @param filename name of gml file containing the kNN graph
#' @export

setGeneric("buildKNN", function(object, k=5, column.ann, num.pcs=20, sim.type="jaccard", filename="graph_clusters.gml") standardGeneric("buildKNN"))
setMethod("buildKNN",
          signature = "CellRouter",
          definition = function(object, k, column.ann,  num.pcs, sim.type, filename){

            suppressWarnings(suppressMessages(library('cccd')))
            suppressWarnings(suppressMessages(library('proxy'))) # Library of similarity/dissimilarity measures for 'dist()'
            matrix <- object@pca$cell.embeddings[,1:num.pcs]
            sampTab <- object@sampTab
            smapTab <- sampTab[order(sampTab[[column.ann]]),]

            print('building k-nearest neighbors graph')
            dm <- as.matrix(dist(matrix))
            h <- nng(dx=dm,k=k)
            if(sim.type == 'jaccard'){
              sim <- similarity.jaccard(h, vids=V(h), loops=FALSE)
            }else if(sim.type == 'invlog'){
              sim <- similarity.invlogweighted(h, vids=V(h))
            }
            el <- get.edgelist(h)
            weights <- sim[el]

            E(h)$weight  <- weights
            V(h)$name <- rownames(matrix)
            edges <- as.data.frame(get.edgelist(h))
            rownames(edges) <- paste(edges$V1, edges$V2, sep='_')
            edges$weight <- as.numeric(E(h)$weight)

            V(h)[rownames(sampTab)]$comms <- as.vector(sampTab[[column.ann]]) #cluster->celltype
            cell.comms <- split(sampTab, sampTab[[column.ann]])
            cell.comms <- lapply(cell.comms, rownames)
            #allcells <- as.vector(unlist(cell.comms))
            #sampTab <- sampTab[allcells,] #change order of cells in the table
            #V(h)[names(nodeLabels(sampTab,column.ann))]$label <- nodeLabels(sampTab, column.ann)

            ## Useful information about the graph
            graph <- list()
            graph[['network']] <- h
            graph[['edges']] <- edges
            graph[['similarity_matrix']] <- sim
            graph[['subpopulation']] <- cell.comms
            #graph[['communities']] <- comms

            print('updating CellRouter object')
            object@graph <- graph
            object@rawdata <- object@rawdata[,rownames(sampTab)]
            object@ndata <- object@ndata[,rownames(sampTab)]
            object@sampTab <- sampTab
            write.graph(graph = h, file = filename, format = 'gml')

            rm(h)
            rm(edges)
            rm(sim)

            return(object)
          }
)
#' Plot heatmap with gene signatures
#' @param object CellRouter object
#' @param markers Genes preferentially expressed in each column.ann. For example, in clusters or sorted populations
#' @param column.ann Column in the metadata table used to annotate the kNN graph. For example, clusters, sorted cell populations
#' @param column.color Color
#' @param num.cells Number of cells to show in the heatmap
#' @param threshold Threshold used to center the data
#' @param genes.show #Vector of gene names to show in the heatmap
#' @param low Color for low expression
#' @param intermediate Color for intermediate expression
#' @param high Color for high expression
#' @param width width
#' @param height height
#' @param filename filename
#' @export

plotSignaturesHeatmap <- function(object, markers, column.ann, column.color, num.cells=NULL, threshold=2, genes.show=NULL,
                                  low='purple', intermediate='black', high='yellow', width, height, filename){
  
  if(is.null(num.cells)){
    cells.keep <- rownames(object@sampTab)
  }else{
    #cells.use <- object@sampTab %>% group_by_(column.ann) %>% sample_n(size = num.cells, replace = TRUE)
    cells.use <- split(object@sampTab, object@sampTab[[column.ann]])
    cells.use <- lapply(cells.use, function(x){
      if(nrow(x) < num.cells){
        cells.use.x <- x[sample(rownames(x), size = nrow(x)),]
      }else{
        cells.use.xx <- x[sample(rownames(x), size = num.cells),]
      }
    })
    cells.use.tmp <- do.call(rbind, cells.use)
    cells.keep <- as.vector(cells.use.tmp$sample_id)
  }

  #data <- object@ndata[,cells.use]
  matrix <- center_with_threshold(object@ndata[,cells.keep], threshold)

  paletteLength <- 100
  #myColor <- colorRampPalette(c("purple","black","yellow"))(paletteLength)
  myColor <- colorRampPalette(c(low, intermediate, high))(paletteLength)
  myBreaks <- c(seq(min(matrix), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(matrix)/paletteLength, max(matrix), length.out=floor(paletteLength/2)))


  library(data.table)
  markers2 <- as.data.frame(markers)
  #markers2 <- as.data.frame(markers)
  #markers2 <- as.data.table(markers2)[, .SD[which.max(fc.column)], by=gene]
  #markers2 <- as.data.frame(markers2)
  #markers2 <- as.data.frame(markers)
  #markers2 <- markers2[!duplicated(markers2$gene),] #need to make sure there is no duplicated element...
  sampTab <- object@sampTab
  sampTab <- sampTab[cells.keep,]

  if(column.ann == 'population'){
    markers2 <- markers2[order(as.numeric(markers2$population)),]
    rownames(markers2) <- as.vector(markers2$gene)
    sampTab <- sampTab[order(as.numeric(sampTab$population)),]
  }else{
    markers2 <- markers2[order(as.character(markers2$population)),]
    rownames(markers2) <- as.vector(markers2$gene)
    sampTab <- sampTab[order(as.character(sampTab[[column.ann]])),]
  }

  #clusters <- as.vector(object@sampTab$population)
  clusters <- as.vector(sampTab[[column.ann]])
  names(clusters) <- rownames(sampTab)
  #clusters <- clusters[order(clusters)]
  ann_col <- data.frame(population=as.vector(clusters), stringsAsFactors = FALSE)
  rownames(ann_col) <- names(clusters)

  ann_row <- data.frame(signature=as.vector(markers2$population), stringsAsFactors = FALSE)
  rownames(ann_row) <- as.vector(markers2$gene)

  #colors <- cRampClust(cluster.vector, 8)
  #names(colors) <- cluster.vector
  colors <- unique(sampTab[[column.color]])
  names(colors) <- unique(as.vector(sampTab[[column.ann]]))

  color_lists <- list(population=colors, signature=colors)

  #replicate_row <- as.vector(unlist(lapply(split(ann_row, ann_row$signature), nrow)))
  #colors_row <- rep(colors, times=replicate_row)
  #replicate_col <- as.vector(unlist(lapply(split(ann_col, ann_col$population), nrow)))
  #colors_col <- rep(colors, times=replicate_col)
  index <- getIndexes(ann_col, ann_row, order.columns = unique(ann_col$population), order.rows = unique(ann_row$signature))

  if(is.null(genes.show)){
    genes.show <- markers2 %>% group_by(population) %>% top_n(5, fc)
    genes.show <- as.vector(genes.show$gene)
    selected <- as.vector(markers2$gene)
    selected[!(selected %in% genes.show)] <- ""
  }else{
    selected <- as.vector(markers2$gene)
    selected[!(selected %in% genes.show)] <- ""
  }

  pheatmap(matrix[rownames(ann_row),rownames(ann_col)], cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, breaks=myBreaks,
           gaps_row = index$rowsep, gaps_col = index$colsep, annotation_col = ann_col, annotation_row = ann_row, annotation_colors = color_lists,
           labels_row = selected,
           labels_col = rep("", ncol(matrix)), width=width, height=height, filename=filename)

  pheatmap(matrix[rownames(ann_row),rownames(ann_col)], cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, breaks=myBreaks,
           gaps_row = index$rowsep, gaps_col = index$colsep, annotation_col = ann_col, annotation_row = ann_row, annotation_colors = color_lists,
           labels_row = selected,
           labels_col = rep("", ncol(matrix)))
  
  gc(verbose = FALSE)
  #pdf(file='results/heatmap_2.pdf', width=10, height=8)
  #heatmap.2(as.matrix(matrix[rownames(ann_row),rownames(ann_col)]), col=myColor,trace="none",
  #          density.info="none", scale="none",margin=c(5,5), key=TRUE, Colv=F, Rowv=F,
  #          srtCol=60, dendrogram="none", cexCol=0.75, cexRow=0.65, labRow=FALSE,
  #          labCol = FALSE, symm=T,symkey=T,
  ##          ColSideColors=colors_col,
  #          RowSideColors=colors_row,
  #          colsep=index$colsep, rowsep=index$rowsep, sepcolor = 'black')
  #dev.off()
}


#' Compute fold changes and find gene signatures
#' @param object CellRouter object
#' @param column Column in the metadata table to group cells for differential expression. For example, if 'population' is specified, population-specific gene signatures will be identified
#' @param pos.only Only uses genes upregulated
#' @param fc.threshold FOld change threshold
#' @export

setGeneric("computeFC", function(object, column='population', pos.only=TRUE, fc.threshold=0.25) standardGeneric("computeFC"))
#computeFC(object, column, pos.only, fc.threshold)
setMethod("computeFC",
          signature = "CellRouter",
          definition = function(object, column, pos.only, fc.threshold){
            print('discovering subpopulation-specific gene signatures')
            expDat <- object@ndata
            #membs <- as.vector(object@sampTab$population)
            membs <- as.vector(object@sampTab[[column]])
            diffs <- list()
            for(i in unique(membs)){
              cat('cluster ', i, '\n')
              if(sum(membs == i) == 0) next
              m <- if(sum(membs != i) > 1) apply(expDat[, membs != i], 1, mean) else expDat[, membs != i]
              n <- if(sum(membs == i) > 1) apply(expDat[, membs == i], 1, mean) else expDat[, membs == i]

              #pv <- binompval(m/sum(m),sum(n),n)
              #d <- data.frame(mean.np=m, mean.p=n, fc=n-m, pv=pv) #log scale
              d <- data.frame(mean.np=m, mean.p=n, fc=n-m) #log scale
              if(pos.only){
                genes.use <- rownames(d[which(d$fc > fc.threshold),])
              }else{
                genes.use <- rownames(d[which(abs(d$fc) > fc.threshold),])
              }
              m <- m[genes.use]
              n <- n[genes.use]
              print(length(genes.use))
              #diffs[[i]] <- d
              #d <- data.frame(mean.np=m, mean.p=n, fc=n/m, pv=pv) #linear scale
              #d <- d[!is.infinite(d$fc),]
              #d <- d[which(d$pv < 0.05),]
              #d <- d[order(d$pv, decreasing=FALSE),]
              #d <- d[order(d$mean.p, decreasing=TRUE),]
              #diffs[[i]] <- d[which(d$pv < pvalue & d$fc > foldchange),]

              #for wilcox test
              coldata <- object@sampTab
              coldata[membs == i, "group"] <- "Group1"
              coldata[membs != i, "group"] <- "Group2"

              coldata$group <- factor(x = coldata$group)
              countdata.test <- as.matrix(expDat[genes.use, rownames(x = coldata)])

              p_val <- sapply(X = 1:nrow(x = countdata.test), FUN = function(x) {
                return(wilcox.test(countdata.test[x, ] ~ coldata$group)$p.value)
              })
              #to.return <- data.frame(p_val, row.names = rownames(countdata.test))
              d2 <- data.frame(mean.np=m, mean.p=n, fc=n-m, pv=p_val, p.adj=p.adjust(p_val, method='bonferroni'))
              #diffs[[i]] <- d2[which(d2$pv < pvalue),]
              diffs[[i]] <- d2
            }
            object@signatures <- diffs
            return (object)
          }
)

#' Find gene signatures based on a template-matching approach
#' @param expDat Expression matrix
#' @param sampTab Sample annotation table
#' @param qtile qyantile to select top genes correlated with the idealized expression pattern
#' @param remove Remove overlaping genes
#' @param dLevel Groups to compare
#' @export

ranked_findSpecGenes<-function# find genes that are preferentially expressed in specified samples
(expDat, ### expression matrix
 sampTab, ### sample table
 qtile=0.95, ### quantile
 remove=FALSE,
 dLevel="population_name" #### annotation level to group on
){
  cat("Template matching...\n")
  myPatternG<-cn_sampR_to_pattern(as.vector(sampTab[,dLevel]));
  specificSets<-apply(myPatternG, 1, cn_testPattern, expDat=expDat);

  # adaptively extract the best genes per lineage
  cat("First pass identification of specific gene sets...\n")
  cvalT<-vector();
  ctGenes<-list();
  ctNames<-unique(as.vector(sampTab[,dLevel]));
  for(ctName in ctNames){
    x<-specificSets[[ctName]];
    cval<-quantile(x$cval, qtile, na.rm = TRUE);
    tmp<-rownames(x[x$cval>cval,]);
    specificSets[[ctName]] <- specificSets[[ctName]][tmp,]
    ctGenes[[ctName]]<-tmp;
    cvalT<-append(cvalT, cval);
  }
  if(remove){
    cat("Remove common genes...\n");
    # now limit to genes exclusive to each list
    specGenes<-list();
    for(ctName in ctNames){
      others<-setdiff(ctNames, ctName);
      x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
      specificSets[[ctName]] <- specificSets[[ctName]][x,]
      specificSets[[ctName]]$gene <- rownames(specificSets[[ctName]])
      specGenes[[ctName]]<-x;
    }
    result <- specGenes
  }else {
    result <- ctGenes;
  }
  specificSets <- lapply(specificSets, function(x){x[order(x$cval, decreasing = TRUE),]})
  specificSets <- lapply(specificSets, function(x){colnames(x) <- c('tm.pval', 'cval', 'tm.padj', 'gene'); x})

  specificSets
}
#' Helper function
#' @param sampR sampR
cn_sampR_to_pattern<-function# return a pattern for use in cn_testPattern (template matching)
(sampR){
  d_ids<-unique(as.vector(sampR));
  nnnc<-length(sampR);
  ans<-matrix(nrow=length(d_ids), ncol=nnnc);
  for(i in seq(length(d_ids))){
    x<-rep(0,nnnc);
    x[which(sampR==d_ids[i])]<-1;
    ans[i,]<-x;
  }
  colnames(ans)<-as.vector(sampR);
  rownames(ans)<-d_ids;
  ans;
}
#' Helper function
#' @param pattern pattern
#' @param expDat expression matrix
cn_testPattern<-function(pattern, expDat){
  pval<-vector();
  cval<-vector();
  geneids<-rownames(expDat);
  llfit<-ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE);
  xxx<-matrix( unlist(llfit$coef), ncol=8,byrow=TRUE);
  ccorr<-xxx[,6];
  cval<- sqrt(as.numeric(llfit$summary[,2])) * sign(ccorr);
  pval<-as.numeric(xxx[,8]);

  #qval<-qvalue(pval)$qval;
  holm<-p.adjust(pval, method='holm');
  #data.frame(row.names=geneids, pval=pval, cval=cval, qval=qval, holm=holm);
  data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
}

#' Identify gene signatures
#' @param object CellRouter object
#' @param column Specify the groups to compare
#' @param test.use Differential expression test to use. Default is wilcox. Alternative is based on template-matching. Possible values: wilcox or template
#' @param pos.only Use only upregulated genes
#' @param fc.threshold Fold change threshold
#' @param fc.tm Wheter to include fold change values in the template-matching differential expression analysis
#' @export

findSignatures <- function(object, column='population', test.use='wilcox', pos.only=TRUE, fc.threshold=0.25, fc.tm=FALSE){
  if(test.use == 'wilcox'){
    cat('Calculating fold changes...', '\n')
    object <- computeFC(object, column, pos.only, fc.threshold)
    markers <- findmarkers(object)
    #markers$gene <- rownames(markers)

  }else if(test.use == 'template'){
    cat('Calculating template-matchings...', '\n')
    signatures <- ranked_findSpecGenes(object@ndata, object@sampTab, qtile=0.99, remove=TRUE, dLevel = column)
    #mylist <- lapply(signatures, functionai(x){x$assignment})
    mylist <- signatures
    for(i in 1:length(mylist) ){ mylist[[i]] <- cbind(mylist[[i]], population=rep(names(mylist[i]), nrow(mylist[[i]])) ) }
    markers <- as.data.frame(do.call(rbind, mylist))
    rownames(markers) <- as.vector(markers.s$gene)
    if(fc.tm){
      object <- computeFC(object, column, pos.only, fc.threshold)
      for(signature in names(object@signatures)){
        markers[rownames(signatures[[signature]]), 'log2FC'] <- object@signatures[[signature]][rownames(signatures[[signature]]), 'fc']
        markers[rownames(signatures[[signature]]), 'log2FC_pval'] <- object@signatures[[signature]][rownames(signatures[[signature]]), 'pv']
        markers[rownames(signatures[[signature]]), 'log2FC_p.adj'] <- object@signatures[[signature]][rownames(signatures[[signature]]), 'p.adj']
      }
    }
  }
  markers
}

##Create a table of genes, fc_subpopulation, subpopulation with max expression
setGeneric("findmarkers", function(object) standardGeneric("findmarkers"))
setMethod("findmarkers",
          signature = "CellRouter",
          definition = function(object){
            print('finding subpopulation markers')
            genes <- unique(as.vector(unlist(lapply(object@signatures, rownames))))
            df <- data.frame(matrix(0, nrow=length(genes), ncol=length(object@signatures)))
            rownames(df) <- genes
            colnames(df) <- names(object@signatures)
            for(gene in rownames(df)){
              for(pop in colnames(df)){
                if(gene %in% rownames(object@signatures[[pop]])){
                  df[gene, pop] <- object@signatures[[pop]][gene, 'fc']
                }else{
                  df[gene, pop] <- 0
                }
              }
            }
            x <- apply(df, 1, function(x){names(x)[which(x == max(x))][1]})
            xx <- apply(df, 1, function(x){max(x)})
            df$population <- x
            df$fc <- xx
            #df <- df[order(df$population), c('population', 'fc')]
            df <- df[, c('population', 'fc')]
            for(pop in names(object@signatures)){
              dfx <- df[which(df$population == pop),]
              df[rownames(dfx), 'pval'] <- object@signatures[[pop]][rownames(dfx), 'pv']
              df[rownames(dfx), 'p.adj'] <- object@signatures[[pop]][rownames(dfx), 'p.adj']
            }
            df$gene <- rownames(df)
            df <- df[which(df$fc > 0),] #new line...
            #df <- df[order(nchar(df$population)),]
            return(df)
          }
)


#' Predict a gene regulatory network
#' @param object CellRouter object
#' @param species species
#' @param genes.use genes to include in the gene regulatory network
#' @param zscore zscore threshold to identify putative regulatory interactions
#' @param filename filename of GRN data
#' @export

setGeneric("buildGRN", function(object, species, genes.use=NULL, zscore=5, filename='GRN.R') standardGeneric("buildGRN"))
setMethod("buildGRN",
          signature = "CellRouter",
          definition = function(object, species, genes.use, zscore=5, filename='GRN.R'){
            if(is.null(genes.use)){
              genes.use <- rownames(object@ndata)
            }

            tfs <- find_tfs(species = species)
            grn <- globalGRN(object@ndata[genes.use,], tfs, zscore)
            colnames(grn)[1:2]<-c("TG", "TF");
            ggrn<- ig_tabToIgraph(grn, simplify = TRUE)
            x <- list(GRN=ggrn, tfs=tfs)
            save(x, file=filename)
            
            return(x)
          }
)

#' Plot violin plot
#' @param object CellRouter object
#' @param geneList gene list to plot
#' @param column column to group on
#' @param cols how many clumns in the output figure
#' @param width width
#' @param height height
#' @param filename filename
#' @import ggplot2
#' @export

plotViolin <- function(object, geneList, column, column.color, cols, width=10, height=5, filename){
  plots <- list()
  sampTab <- object@sampTab
  expDat <- object@ndata
  T0 <- expDat
  allgenes <- data.frame()
  for(g in geneList){
    #cat(time, ' ', dim(T0), '\n')
    genes <- as.data.frame(t(T0[g,]))
    genes$gene <- g
    genes$conditions <- as.vector(sampTab[,column])
    genes.m <- melt(genes, id.var=c('gene',"conditions"))
    allgenes <- rbind(allgenes, genes.m)
  }

  colors <- unique(sampTab[[column.color]])
  names(colors) <- unique(sampTab[[column]])
  order <- unique(allgenes$conditions)
  order <- order[order(order, decreasing=FALSE)]
  allgenes$conditions <- factor(allgenes$conditions, levels=order)
  p <- ggplot(allgenes, aes(x=conditions, y=value, fill=conditions)) +
    geom_violin(scale="width") + stat_summary(fun.y=mean,geom='point', size=1) + #geom_boxplot(alpha=.9) +
    theme_bw() + xlab("") + ylab("") + theme(legend.position="none") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1),
          strip.text.y = element_text(angle=180),
          panel.border=element_rect(fill = NA, colour=alpha('black', 0.75),size=1),
          strip.background = element_rect(colour="white", fill="white"),
          panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(0, "lines")) +
    scale_fill_manual("", values=colors) +
    facet_wrap(~variable, ncol = cols, strip.position = "left") #+ coord_flip()
  #facet_grid(variable ~ ., scales='free') + coord_flip()
  print(p)

  pdf(file=filename, width=width, height=height)
  #multiplot(plotlist = plots, cols=cols)
  print(p)
  dev.off();
  #multiplot(plotlist = plots, cols=cols)
  gc()
}


#' Predict cell cycle phase
#' @param object CellRouter object
#' @param columns columns to be selected from the metadata table
#' @export
predictCellCyle <- function(object, columns){
  cc.scores <- object@sampTab[, columns]
  x <- apply(cc.scores, 1, function(x){
    if(all(x < 0)){
      return('G1')
    }else{
      return(names(x)[which(x == max(x))][1])
    }
  })
  object <- addInfo(object, x, colname = 'Phase')

  object
}

#' Predict state based on scored gene lists
#' @param object CellRouter object
#' @param columns columns to select from the metadata table
#' @param col.name column name to be added to the metadata table after state prediction
#' @export

predictState <- function(object, columns, col.name){
  cc.scores <- object@sampTab[, columns]
  x <- apply(cc.scores, 1, function(x){
    if(all(x < 0)){
      return('Unknown')
    }else{
      return(names(x)[which(x == max(x))][1])
    }
  })
  object <- addInfo(object, x, colname = col.name)

  object
}

#' Score gene sets
#' @param object CellRouter object
#' @param gene.list gene lists for which the scores will be calculated
#' @param  bins #number of bins to split expression data
#' @param genes.combine default is all genes in object@@ndata
#' @export
scoreGeneSets <- function(object, genes.list, bins=25, genes.combine=NULL){
  if(is.null(genes.combine)){
    genes.combine=rownames(object@ndata)
  }
  ctrl.size <- min(unlist(lapply(genes.list, length)))
  genes.list <- lapply(genes.list, function(x){intersect(x, rownames(object@ndata))})
  cluster.length <- length(x = genes.list)
  data.avg <- Matrix::rowMeans(object@ndata[genes.combine, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(
    x = data.avg,
    m = round(length(x = data.avg) / bins)
  ))
  names(data.cut) <- names(data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)

  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(genes.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(sample(data.cut[which(data.cut == data.cut[genes.use[j]])], size = ctrl.size,replace = FALSE))
      )
    }
  }
  ctrl.use <- lapply(ctrl.use, unique)
  ctrl.scores <- matrix(data = numeric(length = 1L),nrow = length(ctrl.use),ncol = ncol(object@ndata))

  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(object@ndata[genes.use, ])
  }

  genes.scores <- matrix(data = numeric(length = 1L),nrow = cluster.length,ncol = ncol(object@ndata))
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    genes.scores[i, ] <- Matrix::colMeans(object@ndata[genes.use, , drop = FALSE])
    #data.use <- object@ndata[genes.use, , drop = FALSE]
    #genes.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  scores <- genes.scores - ctrl.scores
  rownames(scores) <- names(genes.list)#paste0(col.name, 1:cluster.length)
  scores <- as.data.frame(t(scores))
  rownames(scores) <- colnames(object@ndata)

  #object@sampTab <- cbind(object@sampTab, scores[rownames(object@sampTab),])
  #object@sampTab <- addInfo(object, colname=colnames(scores)) #cbind(object@sampTab, scores[rownames(object@sampTab),])
  for(col.name in colnames(scores)){
    object <- addInfo(object, scores, colname = col.name, metadata.column = col.name)
  }
  object
}


##########################################################

#' Initialize CellRouter object
#' @param .Object object
#' @param rawdata raw data provided as input
#' @param path path from which the raw data will be loaded
#' @param min.genes keep only cells expressing at least min.gene
#' @param min.cells Keep only genes expressed in at least min.cells
#' @param is.expr Threshold to determine whether a gene is expressed or not. By default, genes with raw counts > 0 are considered as expressed.
#' @export

setMethod("initialize",
          signature = "CellRouter",
          definition = function(.Object, rawdata=NULL, path, min.genes=0, min.cells=0, is.expr=0){
            print("Initializing CellRouter object")
            if(is.null(rawdata)){
              rawdata <- as.data.frame(get(load(path)))
            }

            num.genes <- colSums(rawdata > is.expr)
            num.mol <- colSums(rawdata)
            cells.use <- names(num.genes[which(num.genes > min.genes)])
            expdat <- rawdata[, cells.use]
            rawdata <- rawdata[, cells.use]
            genes.use <- rownames(rawdata)
            if (min.cells > 0) {
              num.cells <- rowSums(rawdata > 0)
              genes.use <- names(num.cells[which(num.cells >= min.cells)])
              rawdata <-rawdata[genes.use, ]

            }
            nGene <- num.genes[cells.use]
            nUMI <- num.mol[cells.use]

            .Object@sampTab <- data.frame(sample_id=colnames(rawdata), nGene=nGene, nUMI=nUMI, conditions=colnames(expdat))
            rownames(.Object@sampTab) <- .Object@sampTab$sample_id
            .Object@rawdata <- rawdata
            .Object@ndata <- rawdata
            validObject(.Object)
            return(.Object)
          }
)

#' Add medata information to CellRouter metadata in object@@sampTab
#' @param object CellRouter object
#' @param metadata Metadata to be added
#' @param colname column name to be added to object@@sampTab in case metadata is a vector
#' @param metadata.column column to selected from the metadata to be added and included in object@@sampTab
#' @export

addInfo <- function(object, metadata, colname, metadata.column='population'){ #uupdate to include data.frames as well...
  sampTab <- object@sampTab
  if(class(metadata) == 'data.frame'){
    sampTab[rownames(metadata), colname] <- as.vector(metadata[[metadata.column]])
  }else{
    sampTab[names(metadata), colname] <- metadata
  }

  colors <- cRampClust(1:length(unique(sampTab[[colname]])), 8) #change $ of colors more properly...
  names(colors) <- unique(sampTab[[colname]])

  replicate_row <- as.vector(unlist(lapply(split(sampTab, sampTab[[colname]]), nrow)))
  colors_row <- rep(colors, times=replicate_row)
  color.column <- paste(colname, 'color',sep='_')
  sampTab[,color.column] <- colors_row

  object@sampTab <- sampTab
  object
}

#' Quality control. Filter out cells based on variables provided in the parameter variables
#' @param object CellRouter object
#' @param variables Filter out cells based on these variables, such as number of detected genes or mitochondrial content
#' @param threshold.low Cells with values lower than the ones provided here are filtered out
#' @param threshold.high Cells with values higher than the ones provided here are filtered out
#' @export

filterCells <- function(object, variables, thresholds.low, threshold.high){
  sampTab <- object@sampTab

  for(v in 1:length(variables)){
    sampTab <- sampTab[which(as.vector(sampTab[,variables[v]]) < thresholds.high[v] & as.vector(sampTab[,variables[v]]) > thresholds.low[v]),]
  }

  object@rawdata <- object@rawdata[,rownames(sampTab)]
  object@ndata <- object@ndata[,rownames(sampTab)]
  object@sampTab <- sampTab

  object
}


#' Proportion plot.
#' @param object CellRouter object
#' @param condition Column in the metadata table specifying an annotation, such as sorted populations
#' @param population Column in the metddata table specifying another annotation
#' @param width width
#' @param height height
#' @param filename filename
#' @import ggplot2
#' @export
plotProportion <- function(object, condition, population, width, height, filename){
  samples <- object@sampTab
  data2 <- data.frame(cells=samples[[condition]], classification=samples[[population]])
  data2$classification <- factor(data2$classification, levels=unique(as.vector(data2$classification)))
  colors <- as.vector(unique(samples$colors))
  names(colors) <- unique(as.vector(samples$population))

  pdf(file=filename, width = width, height=height)
  g <- ggplot(data2,aes(x = classification, fill = cells)) +
    geom_bar(position = "fill", color='black') +
    theme_bw() + theme(legend.position='right',
                       legend.key.size = unit(0.3, "cm"),
                       legend.text=element_text(size=7)) +
    xlab("") + ylab("Proportion") +
    theme(axis.text.x = element_text(size=10, angle=45, hjust=1), axis.title.y = element_text(size = rel(1), angle = 90)) +
    #scale_fill_manual("", values=colors)
  scale_fill_brewer("", palette = 'Paired')
  print(g)
  dev.off()

}

#' Identify trajectories connecting source and target populations in the kNN graph
#' @param object CellRouter object
#' @param column Column in the metadata table specifying wheter transitions are between clusters or other annotations, such as sorted populations
#' @param libdir Path to Java libraries required
#' @param maindir Directory
#' @param method Method used to determine the source and target cells based on the source and target populations
#' @export
#'
setGeneric("findPaths", function(object, column='population',libdir, maindir, method) standardGeneric("findPaths"))
setMethod("findPaths",
          signature="CellRouter",
          definition=function(object, column, libdir, maindir, method){
            curdir <- getwd()
            dirs <- list()

            if(method %in% c("euclidean", "maximum", "manhattan","canberra","binary")){
              #bla2 <- as.data.frame(as.matrix(dist(object@rdimension, method = method)))
              tmp <- as.data.frame(as.matrix(dist(object@rdimension, method = method)))
            }else{
              g <- object@graph$network
              ##g <- simplify(g, remove.loops = TRUE)
              ##E(g)$weight <- -1*E(g)$weight
              ##bla2 <- as.data.frame(shortest.paths(g, v=V(g), to=V(g)))
              #bla2 <- as.data.frame(igraph::distances(g, v=V(g), to=V(g), weights = NULL, algorithm = "bellman-ford"))
              tmp <- as.data.frame(igraph::distances(g, v=V(g), to=V(g), weights = NULL, algorithm = "bellman-ford"))
            }
            sampTab <- object@sampTab
            b <- list()
            for(p1 in sources){
              cellsp1 <- as.vector(sampTab[which(sampTab[[column]] == p1), 'sample_id'])
              for(p2 in targets){
                if(p1 != p2){
                  cellsp2 <- as.vector(sampTab[which(sampTab[[column]] == p2), 'sample_id'])
                  #x <- bla2[cellsp1, cellsp2]
                  x <- tmp[cellsp1, cellsp2]
                  x$population1 <- p1
                  x$population2 <- p2
                  x$merge <- rownames(x)

                  #bla3 <- melt(x, id.vars=c('population1', 'population2', 'merge'))
                  #bla3 <- bla3[order(bla3$value, decreasing = TRUE),]
                  tmp2 <- melt(x, id.vars=c('population1', 'population2', 'merge'))
                  tmp2 <- tmp2[order(tmp2$value, decreasing = TRUE),]

                  #b[[p1]][[p2]] <- bla3[1,]
                  b[[p1]][[p2]] <- tmp2[1,]
                }
              }
            }
            #for(i in 1:length(sources)){
            #  for(j in 1:length(targets)){
            for(i in names(b)){
              for(j in names(b[[i]])){

                #ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
                s <- as.vector(b[[i]][[j]][,'merge'])
                t <- as.vector(b[[i]][[j]][,'variable'])
                dir <- paste(i, j, sep='.')
                cat('-------------------Transition:', dir, ' -----------------------\n')
                dir.create(file.path(maindir, dir), showWarnings = FALSE)
                system(paste('chmod 777 ', file.path(maindir, dir), sep=''))
                setwd(file.path(maindir, dir))
                dirs[[dir]] <- file.path(maindir, dir)
                write.table(s, file='cell_sources.txt', sep=',', row.names=FALSE, quote=FALSE)
                write.table(t, file='cell_sinks.txt', sep=',', row.names=FALSE, quote=FALSE)

                sink("CellRouter.sh")
                cat('#!/bin/sh\n')
                cat('# run CellRouter, run!\n')
                #cat('LIB_DIR=~/NetBeansProjects/CellRouter/\n')
                cat(paste('LIB_DIR=', libdir, sep=''), '\n')
                cat('LIB=$LIB_DIR/dist/CellRouter.jar:$LIB_DIR/lib/commons-cli-1.2.jar\n')
                cat(paste('OUT_DIR=', maindir, sep=''), '\n')

                cat(paste('NET=', maindir, '/cell_edge_weighted_network.txt', sep=''), '\n')
                cat('network=Cells_FlowNetwork\n')
                cat('source=cell_sources.txt\n')
                cat('sink=cell_sinks.txt\n')
                cat('topPaths=25\n')

                cat('echo "1) Computing flow network"\n')
                cat('java -cp $LIB "cellrouter.CellRouter" -NET $NET -source $source -sink $sink -topPaths $topPaths -out $OUT_DIR -C -f $network\n')
                sink()

                system('chmod 777 CellRouter.sh')
                system('./CellRouter.sh')

                setwd(file.path(curdir))
              }
            }
            object@directory <- dirs

            return(object)
          }
)

#' Process trajectories
#' @param object CellRouter object
#' @param genes Vector of gene names that are used for trajectory analysis
#' @param path.rank How to rank paths. See tutorials for examples
#' @param num.cells Trajectories should contain at least num.cells
#' @param neighs The size of the neighborhood in kNN graph used to smoothen kinetic profiles
#' @param column.ann Transitions between the cell states provided here are identified, such as clusters or sorted populations
#' @param column.color The colors corresponsding to the annotation in column.ann
#' @export
setGeneric("processTrajectories", function(object, genes, path.rank, num.cells, neighs,
                                           column.ann='population', column.color='colors') standardGeneric("processTrajectories"))
setMethod("processTrajectories",
          signature="CellRouter",
          definition=function(object, genes, path.rank, num.cells, neighs, column.ann, column.color){
            library(igraph)

            object@genes.trajectory <- genes

            print('parsing trajectory information')
            paths <- lapply(object@directory, function(x){read.csv(paste(x, 'Cells_FlowNetwork_all_paths.txt', sep='/'),
                                                                            sep="\t", stringsAsFactors=FALSE)})
            paths <- paths[lapply(paths, nrow) > 1]
            #paths <- lapply(paths, function(x){x[order(x$path_flow, decreasing = TRUE),]})
            if(path.rank == "path_cost"){
              paths <- lapply(paths, function(x){x[order(x[,path.rank], decreasing = FALSE),]})
            }else{
              paths <- lapply(paths, function(x){x[order(x[,path.rank], decreasing = TRUE),]})
            }
            paths <- lapply(paths, function(x){x[duplicated(x$path),]})
            paths <- do.call(rbind, lapply(paths, function(x){x[1,]}))

            #remove empty paths
            paths <- paths[complete.cases(paths),]
            #opening flow networks
            networks <- lapply(object@directory,
                               function(x){
                                 file <- list.files(x, '*.gml*');
                                 if(length(file) > 0){
                                   cat(file, '\n');
                                   read.graph(paste(x, 'Cells_FlowNetwork_all_paths_subnet.gml', sep='/'), format = 'gml')
                                 }
                               })
            networks <- networks[lapply(networks, length) > 0]

            paths$ID <- paste('path', 1:nrow(paths), sep='_')
            paths$population <- rownames(paths)
            object@paths <- paths
            object@networks <- networks

            object@pathsinfo <- pathsinfo(object, num.cells = num.cells, neighs = neighs, column.ann = column.ann, column.color = column.color)
            return (object)
          }

)

#' Helper function
#' @param object CellRouter object
#' @param num.cells num.cells
#' @param neighs neighs
#' @param column.ann column.ann
#' @param column.color column.color

setGeneric("pathsinfo", function(object, num.cells, neighs, column.ann='community', column.color='colors') standardGeneric("pathsinfo"))
setMethod("pathsinfo",
          signature="CellRouter",
          definition=function(object, num.cells, neighs, column.ann, column.color){

            paths <- object@paths
            keep <- intersect(rownames(object@ndata), object@genes.trajectory)
            expDat <- object@ndata[keep,]
            sampTab <- object@sampTab
            print(neighs)

            networks <- object@networks
            means <- list()
            o <- neighs
            for(transition in rownames(paths)){
              g <- networks[[transition]]
              path <- paths[transition, 'path']
              split_paths <- strsplit(as.vector(path), "->")[[1]]
              cells <- split_paths[c(-1,-length(split_paths))]
              mean <- list()
              for(cell in cells){
                neighs <- igraph::induced.subgraph(graph=g,vids=unlist(igraph::neighborhood(graph=g,order=o,nodes=cell)))
                neigh.names <- V(neighs)$name
                mean[[cell]] <- apply(expDat[,neigh.names], 1, mean)
                #mean <- append(mean, apply(expDat[,neigh.names], 1, mean))
              }
              mean.df <- do.call(cbind, mean)
              means[[transition]] <- mean.df
            }

            i <- 1
            path_info <- list()
            pathsDF <- data.frame()
            for(transition in rownames(paths)){
              path <- paths[transition, 'path']
              #for(path in as.vector(paths$path)){
              split_paths <- strsplit(as.vector(path), "->")[[1]]
              cells <- split_paths[c(-1,-length(split_paths))]
              if(length(cells) > num.cells){ #only include paths with more than num.cells cells
                #expression <- expDat[, cells]
                expression <- means[[transition]]
                #path_name <- paste("path_", i, sep="")
                path_name <- paths$population[i]
                print(path_name)

                path_info[['distr']][[path_name]] <- expression #rownames(expression) #
                path_info[['path']][[path_name]] <- cells
                path_info[['pathInfo']][[path_name]] <- data.frame(path=paste(as.vector(cells), collapse=','),
                                                                   source=cells[1], sink=cells[length(cells)], cost=paths$path_cost[i],
                                                                   source_population=sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),column.ann],
                                                                   target_population=sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), column.ann])
                pathsDF[path_name, 'source'] <- cells[1]
                pathsDF[path_name, 'sink'] <- cells[length(cells)]
                pathsDF[path_name, 'cost'] <- paths$path_cost[i]
                pathsDF[path_name, 'source_population'] <- as.character(sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),column.ann])
                pathsDF[path_name, 'source_color'] <- as.character(sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),column.color])
                pathsDF[path_name, 'target_population'] <- as.character(sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), column.ann])
                pathsDF[path_name, 'target_color'] <- as.character(sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), column.color])

                pseudotime <- (0:(length(cells)-1)) / length(cells) #stick to this
                names(pseudotime) <- cells
                path_info[['pseudotime']][[path_name]] <- pseudotime
                path_info$path_data <- pathsDF
                i <- i+1
              }
            }

            return (path_info)
          }
)

#' Plot kinetic trends for genes in each transcriptopnal cluster identified
#' @param object CellRouter object
#' @param p Selected trajectory
#' @param columns Number of columns in the output figure
#' @param width width
#' @param height height
#' @param filename filename
#' @import ggplot2
#' @export
setGeneric('plotClusters', function(object, p, columns, width, height, filename) standardGeneric('plotClusters'))
setMethod('plotClusters',
          signature="CellRouter",
          definition=function(object, p, columns, width, height, filename){
            # g <- names(object@top.correlations[[direction]][[p1]]) #branch of interest
            # g <- intersect(g, rownames(cellrouter@dynamics[[p2]])) #other branch
            # df <- object@dynamics[[p1]][g,]
            # c <- hclust(dist(df, method='euclidean'), method='ward.D')
            # df <- df[c$order,]
            # df <- as.data.frame(t(scale(t(df))))
            #
            plots <- list()
            clustering <- object@clusters[[p]]$clustering
            for(cl in unique(clustering)){
              print(cl)
              g <- names(clustering[which(clustering==cl)])
              df <- object@dynamics[[p]][g,]
              c <- hclust(dist(df, method='euclidean'), method='ward.D')
              df <- df[c$order,]
              #df <- as.data.frame(t(scale(t(df))))
              df <- center_with_threshold(df, 2)
              matrix <- as.data.frame(df)
              matrix$gene <- rownames(df)
              matrix.m <- melt(matrix, id.var="gene")
              matrix.m$gene <- factor(rownames(matrix), levels=rev(rownames(matrix)))
              g1 <- ggplot(matrix.m, aes(variable, gene)) + geom_tile(aes(fill = value)) +
                scale_fill_gradientn("zscore",colours=c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2")) + theme_bw() +
                #xlab("CellRouter trajectory") + ylab("") +
                xlab("") + ylab(paste(length(g), ' genes')) +
                #scale_fill_gradient2("zscore", low="navy", high="red") + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
                theme(legend.position="none", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
                      axis.title.y = element_text(size = rel(2), angle = 90), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
                      panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) #+
                #ggtitle(cl)

              file <- paste(filename, '_', cl, '.png', sep='')
              png(file=file, width=200, height=150)# units = 'in', res = 1000)
              #pdf(file=file, width=3, height=2)
              plot(g1)
              #multiplot(plotlist = plots, cols=columns)
              dev.off()

              plots[[cl]] <- g1
            }

            file <- paste(filename, 'combined_', '.png',sep='')
            #pdf(file=file, width=width, height=height)
            png(file=file, width=width, height=height)
            multiplot(plotlist = plots, cols=columns)
            dev.off()
            multiplot(plotlist = plots, cols=columns)
          }
)

#' Plot branch-specific kinetic patterns
#' @param object CellRouter object
#' @param direction Plot genes up or down-regulated along trajectories
#' @param p1 Trajectory 1
#' @param p2 Trajecotory 2
#' @param columns Number of columns in the figure generated
#' @param width width
#' @param height height
#' @param filename filename
#' #' @import ggplot2
#' @export
setGeneric('plotbranch', function(object, direction, p1, p2, columns, width, height, filename) standardGeneric('plotbranch'))
setMethod('plotbranch',
          signature="CellRouter",
          definition=function(object, direction, p1, p2, columns, width, height, filename){
            g <- names(object@top.correlations[[direction]][[p1]]) #branch of interest
            g <- intersect(g, rownames(cellrouter@dynamics[[p2]])) #other branch
            df <- object@dynamics[[p1]][g,]
            c <- hclust(dist(df, method='euclidean'), method='ward.D')
            df <- df[c$order,]
            df <- as.data.frame(t(scale(t(df))))

            plots <- list()
            matrix <- as.data.frame(df)
            matrix$gene <- rownames(df)
            matrix.m <- melt(matrix, id.var="gene")
            matrix.m$gene <- factor(rownames(matrix), levels=rev(rownames(matrix)))
            g1 <- ggplot(matrix.m, aes(variable, gene)) + geom_tile(aes(fill = value)) +
              scale_fill_gradientn("zscore",colours=c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2")) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              #scale_fill_gradient2("zscore", low="navy", high="red") + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              theme(legend.position="right", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
                    axis.title.y = element_text(size = rel(0.3), angle = 90), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
                    panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
              ggtitle(p1)

            plots[[p1]] <- g1


            df <- object@dynamics[[p2]][g,] #dynamics of genes in p1 in branch p2
            #c <- hclust(dist(df, method='euclidean'), method='ward.D')
            df <- df[c$order,]
            df <- as.data.frame(t(scale(t(df))))

            matrix <- as.data.frame(df)
            matrix$gene <- rownames(df)
            matrix.m <- melt(matrix, id.var="gene")
            matrix.m$gene <- factor(rownames(matrix), levels=rev(rownames(matrix)))
            g2 <- ggplot(matrix.m, aes(variable, gene)) + geom_tile(aes(fill = value)) +
              scale_fill_gradientn("zscore",colours=c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2")) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              #scale_fill_gradient2("zscore", low="navy", high="red") + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              theme(legend.position="right", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
                    axis.title.y = element_text(size = rel(0.3), angle = 90), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
                    panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
              ggtitle(p2)

            plots[[p2]] <- g2
            #pdf(file=filename, width=width, height=height)
            png(file=filename, width = width, height = height, units = 'in', res = 500)
            multiplot(plotlist = plots, cols=columns)
            dev.off()

            multiplot(plotlist = plots, cols=columns)
          }
)
#' Plot derivative and kinetic patterns os predicted regulators of cell-fate transitions
#' @param object CellRouter object
#' @param p Trajectory
#' @param scores Scores of transcriptional regulators
#' @param cluster Cluster kinetic patterns. Default is TRUE
#' @param columns Number of columns in the figure generated
#' @param width width
#' @param height height
#' @param filename filename
#' @export

setGeneric('plottr', function(object, p, scores, cluster=TRUE, columns, width, height, filename) standardGeneric('plottr'))
setMethod('plottr',
          signature="CellRouter",
          definition=function(object, p, scores, cluster=TRUE, columns, width, height, filename){
            colors <- c("navy","white","red")
            #derivative plot
            matrix <- object@dynamics[[p]][names(scores),]
            time <- 1:501
            matrix <- as.data.frame(t(apply(matrix, 1, function(x){diff(x)/diff(time)}))) #derivative analysis
            colnames(matrix) <- 1:500
            #newr <- c(-max(matrix), max(matrix))
            #matrix <- as.data.frame(t(apply(matrix, 1, function(x){rescale(x, newr)})))
            hc <- hclust(dist(matrix), method='ward.D')
            if(cluster==TRUE){
              matrix <- matrix[hc$order,]
            }
            matrix2 <- matrix

             paletteLength <- 100
             myColor <- colorRampPalette(c("navy","white","red"))(paletteLength)
             myBreaks <- c(seq(min(matrix), 0, length.out=ceiling(paletteLength/2) + 1),
                           seq(max(matrix)/paletteLength, max(matrix), length.out=floor(paletteLength/2)))
             #pheatmap(matrix, color = myColor, breaks = myBreaks,
            #          cluster_cols = FALSE, clustering_method = 'ward.D', main=p,
            #          show_colnames = FALSE)

            matrix$cluster <- rownames(matrix)
            matrix.m <- melt(matrix, id.var="cluster")
            matrix.m$cluster <- factor(rownames(matrix), levels=rev(rownames(matrix)))
            g2 <- ggplot(matrix.m, aes(variable, cluster)) + geom_tile(aes(fill = value)) +
              #scale_fill_gradientn("Derivative",colours=colors) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              scale_fill_gradient2("Derivative", low="navy", high="red") + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              #scale_fill_gradient("", low="navy", high) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              theme(legend.position="bottom", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
                    axis.title.y = element_text(size = rel(0.3), angle = 90), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x=element_blank(),axis.ticks=element_blank(),
                    panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
              ggtitle(p)

            #matrix <- object@dynamics[[p]][names(scores),]
            matrix <- object@dynamics[[p]][rownames(matrix),]
            matrix <- as.data.frame(t(apply(matrix, 1, function(x){rescale(x, c(0,1))})))
            matrix$cluster <- rownames(matrix)
            matrix.m <- melt(matrix, id.var="cluster")
            matrix.m$cluster <- factor(rownames(matrix), levels=rev(rownames(matrix)))
            pdf(file=filename, width=width, height=height)
            g1 <- ggplot(matrix.m, aes(variable, cluster)) + geom_tile(aes(fill = value)) +
              scale_fill_gradientn("Expression",colours=colors) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              #scale_fill_gradient("", low="navy", high) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              theme(legend.position="bottom", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
                    axis.title.y = element_text(size = rel(0.3), angle = 90), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x=element_blank(),axis.ticks=element_blank(),
                    panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
              ggtitle(p)

            plots <- list(expression=g2, derivative=g1)

            pdf(file=filename, width=width, height=height)
            #png(file=filename, width = width, height = height, units = 'in', res = 1000)
            multiplot(plotlist = plots, cols=columns)
            dev.off()

            multiplot(plotlist = plots, cols=columns)

            #matrix2
            # matrix <- object@dynamics[[p]][names(scores),]
            # time <- 1:501
            # matrix <- as.data.frame(t(apply(matrix, 1, function(x){diff(x)/diff(time)}))) #derivative analysis
            # paletteLength <- 100
            # myColor <- colorRampPalette(c("navy","white","red"))(paletteLength)
            # myBreaks <- c(seq(min(matrix), 0, length.out=ceiling(paletteLength/2) + 1),
            #               seq(max(matrix)/paletteLength, max(matrix), length.out=floor(paletteLength/2)))
            # pheatmap(matrix, color = myColor, breaks = myBreaks,
            #          cluster_cols = FALSE,  clustering_method = 'ward.D', main=p,
            #          show_colnames = FALSE, filename=filename)

          }
)


#' Compute GRN scores for transcriptional regulators
#' @param object CellRouter object
#' @param ggrn Gene regulatory network
#' @param tfs Vector of gene names of transcriptional regulators
#' @param transitions Selected transitions of interest
#' @param direction Use transcriptional regulators that are up or down-regulated or both
#' @param dir.targets The predicted targets are up or down-regulated
#' @param q.up Cutoff to select top q.up transcriptional regulators
#' @param q.down Cutoff to selec top q.down transcriptional regulators
#' @param columns Number of columns in the figure generated
#' @param width width
#' @param height height
#' @param flip Plot in the horizontal or vertical
#' @param filename filename
#' @export

setGeneric('grnscores', function(object, ggrn, tfs, transitions, direction, dir.targets='up',
                                 q.up=0.95, q.down=0.05,
                                 columns=1, width, height, flip, filename) standardGeneric('grnscores'))
setMethod('grnscores',
          signature="CellRouter",
          definition=function(object, ggrn, tfs, transitions, direction,
                              dir.targets, q.up, q.down,
                              columns, width, height, flip, filename){
            plots <- list()
            #ps <- list()
            allscores <- list()
            for(p in transitions){
              if(direction=='up'){
                tfs.transition <- intersect(tfs,names(object@top.correlations[['up']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
              }else if(direction == 'down'){
                tfs.transition <- intersect(tfs,names(object@top.correlations[['down']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
              }else{
                tfs.transition.up <- intersect(tfs,names(object@top.correlations[['up']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
                tfs.transition.down <- intersect(tfs,names(object@top.correlations[['down']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
                tfs.transition <- c(tfs.transition.up, tfs.transition.down)
              }

              tfs.transition <- intersect(V(ggrn)$name, tfs.transition)
              tfs.transition <- object@correlation[[p]][tfs.transition]
              averages <- vector()
              num.genes <- vector()
              names <- vector()
              tf.targets <- list()
              for(r in names(tfs.transition)){
                rgrn <- igraph::induced.subgraph(ggrn, vids=unlist(igraph::neighborhood(graph=ggrn,order=1,nodes=r)))
                x <- object@top.correlations[[dir.targets]][[p]]
                genes <- intersect(V(rgrn)$name, names(x)) #subnetwork active during transition p
                corrs <- object@correlation[[p]][genes]
                if(length(corrs) == 0){
                  #cat(r, 'has no targets\n')
                }else if(length(corrs) == 1 & names(corrs) == r){
                  #cat(r, 'regulates only itself\n')
                }else if(length(corrs) > 0){ #at least one target required
                  tf.targets[[r]] <- names(corrs)
                  averages <- append(averages, mean(corrs))
                  num.genes <- append(num.genes, length(corrs))
                  names <- append(names, r)
                }
                #df <- data.frame(gene=as.vector(names(bla)), corr=as.numeric(bla), TF=rep(r, length(bla)), stringsAsFactors = FALSE)
                #ps[[r]] <- df
              }
              names(averages) <- names
              names(num.genes) <- names
              aux <- averages[is.na(averages)]
              averages <- averages[!is.na(averages)]
              averages <- averages[order(averages, decreasing=TRUE)]
              num.genes <- num.genes[names(averages)]
              #rescale num.genes
              num.genes <- rescale(num.genes, newrange = c(0.01,1)) #when it rescaled, it changes the scores when only up or down-regulated genes are included...
              averages <- abs(averages)
              scores <- tfs.transition[names(averages)] * averages * num.genes #* number of genes regulated
              #scores <- scores[order(scores, decreasing=TRUE)]
              #scores <- scores[1:num.genes]
              ## if up or down, q.up or q.down are the top genes
              ##if both, q.up or q.down are quantiles
              if(direction=='up'){
                scores <- scores[order(scores, decreasing=TRUE)]
                scores <- scores[1:q.up]
                colors <- c('white', 'orange', 'red')
              }else if(direction=='down'){
                scores <- scores[order(scores, decreasing=FALSE)]
                scores <- scores[1:q.down]
                colors <- c('blue','lightblue','white')
              }else{
                scores <- scores[which(scores > quantile(scores, q.up) | scores < quantile(scores, q.down))]
                colors <- c('blue', 'white', 'red')
              }
              scores <- scores[order(scores, decreasing=TRUE)]
              targets <- tf.targets[names(scores)]
              allscores[[p]] <- list(scores=scores, targets=targets)
              df <- data.frame(gene=names(scores), score=as.numeric(scores))
              df <- df[order(df$score, decreasing = TRUE),]
              angle=30
              if(flip){
                df$gene <- factor(df$gene, levels=rev(df$gene))
              }else{
                df$gene <- factor(df$gene, levels=df$gene)
                angle=45
              }

              #pdf(file=filename, width=width, height=height)
              g <- ggplot(df, aes(x=gene, y=score, fill=score)) + #geom_violin(scale="width") + stat_summary(fun.y=mean,geom='point') + #geom_boxplot(alpha=.9) +
                geom_bar(stat = 'identity', color='black') + #geom_boxplot(alpha=.9) +
                scale_fill_gradientn("",colours=colors) +
                #ggplot(means, aes(x=p, y=mean, fill=p)) + geom_bar(stat='identity') +
                theme_bw() + xlab("") + ylab("GRN score") + theme(legend.position="none") +
                theme(axis.text.x = element_text(size=rel(1), angle=angle, hjust=1)) +
                ggtitle(paste(p)) +
                theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      plot.background=element_blank(),
                      panel.border = element_blank(),
                      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) #+
              #coord_flip()
              if(flip){
                g <- g + coord_flip()
              }
              plots[[p]] <- g
              #panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1))
              #print(g)
              #dev.off()
            }
            file <- paste(filename, '_target_direction_', dir.targets, '.pdf', sep='')
            pdf(file = file, width=width, height=height)
            multiplot(plotlist = plots,cols=columns)
            dev.off()

            multiplot(plotlist = plots,cols=columns)

            allscores
            #list(scores=scores, targets=tf.targets)
          }
)



#' Smooth dynamics along trajectories. Helpful to cluster complex transcriptopnal patterns
#' @param object CellRouter object
#' @param names Selected trajectories. Default consists of all trajectories
#' @export
setGeneric("smoothDynamics", function(object, names) standardGeneric("smoothDynamics"))
setMethod("smoothDynamics",
          signature="CellRouter",
          definition=function(object, names){

            dynamics <- list()
            #ndata <- object@ndata
            expDat <- object@pathsinfo$distr
            #geneList <- object@genes.trajectory
            for(path in names(object@pathsinfo$distr[names])){
              print(path)
              x_axis <- as.numeric(object@pathsinfo$pseudotime[[path]])
              #geneList <- intersect(geneList, rownames(expDat[[path]]))
              #geneList <- rownames(expDat[[path]][geneList,])
              geneList <- rownames(expDat[[path]]) #object@pathsinfo$distr[[path]]
              smoothDynamics <- data.frame(matrix(0, nrow=length(geneList), ncol=501))
              rownames(smoothDynamics) <- geneList
              #cells <- object@pathsinfo$path[[path]]
              #tmp.ndata <- ndata[,cells]
              for(gene_id in geneList){
                y_axis <- as.numeric(expDat[[path]][gene_id, ])
                #y_axis <- as.numeric(tmp.ndata[gene_id, ])
                lo <- loess(y_axis~x_axis)
                xl <- seq(min(x_axis),max(x_axis), (max(x_axis) - min(x_axis))/500)
                y_axis <- predict(lo,xl)
                #y_axis <- rescale(y_axis, newrange = c(0,1))

                smoothDynamics[gene_id, ] <- as.numeric(y_axis)
              }
              dynamics[[path]] <- smoothDynamics
            }
            object@dynamics <- dynamics

            return(object)
          }
)

#' Cluster kinetic profiles along each trajectory into num.clusters
#' @param object CellRouter object
#' @param num.clusters Number of clusters
#' @export

setGeneric("clusterGenesPseudotime", function(object, num.clusters) standardGeneric("clusterGenesPseudotime"))
setMethod("clusterGenesPseudotime",
          signature="CellRouter",
          definition=function(object, num.clusters){

            library(cluster)
            clusters <- list()
            trajectories <- object@dynamics

            for(trajectory in names(trajectories)){
              cat(trajectory, '\n')
              matrix <- trajectories[[trajectory]]
              x <- suppressWarnings(clustergenes(matrix, num.clusters))
              clusters[[trajectory]] <- x
            }
            object@clusters <- clusters

            return(object)
          }
)
#' Rank genes in each transcriptional cluster based on their correlation with the mean kinetic profile
#' @param object CellRouter object
#' @param num.genes Top num.genes in each transcriptional cluster
#' @export

setGeneric("rankGenesTranscriptionalClusters", function(object, num.genes) standardGeneric("rankGenesTranscriptionalClusters"))
setMethod("rankGenesTranscriptionalClusters",
          signature="CellRouter",
          definition=function(object, num.genes){

            clusters <- object@clusters
            clusterTrends <- list()
            for(trajectory in names(clusters)){
              tmp <- list()
              exprs <- clusters[[trajectory]][['exprs']]
              clustering <- clusters[[trajectory]][['clustering']]
              for(cluster in 1:max(clustering)){
                x <- clustering[which(clustering==cluster)]
                xx <- exprs[names(x),]
                xxx <- apply(xx, 2, mean)
                cors <- apply(exprs, 1, function(x){cor(as.numeric(x), as.numeric(xxx),method = 'spearman')})
                cors <- cors[order(cors, decreasing=TRUE)]
                cors <- cors[1:num.genes]
                #threshold.up <- quantile(cors, max_quantile, na.rm=TRUE)
                #threshold.down <- quantile(cors, min_quantile, na.rm=TRUE)
                #cors.pos <- cors[cors > threshold.up]
                #cors.neg <- cors[cors < threshold.down]
                #cors.pos <- cors.pos[order(cors.pos, decreasing=TRUE)]
                #cat(length(cors.pos), threshold.up, length(cors), cluster, '\n')
                #cors.neg <- cors.neg[order(cors.neg)]
                df <- data.frame(genes=names(cors), cors=cors, cluster=cluster, trend='positive')
                tmp[[paste('cl',cluster,sep='')]] <- df
                #clusterTrends[['up']][[paste('cl',cluster,sep='')]][[trajectory]] <- cors.up
                #clusterTrends[['down']][[paste('cl',cluster,sep='')]][[trajectory]] <- cors.down
              }
              clusterTrends[[trajectory]] <- do.call(rbind, tmp)

            }
            clusterTrends
          }
)

#' Plot top-ranked genes in transcriptional clusters
#' @param rankedDynamics Ranked dynamics as calculated by rankGenesTranscriptionalClusters
#' @param traj.name Selected trajectory
#' @param num.genes Number of genes to plot
#' @param columns Number of columns in the figure gernerated
#' @param width width
#' @param height height
#' @param filename filename
#' @export

setGeneric("plotGenesInClusters", function(rankedDynamics, traj.name, num.genes, columns=5, width=23, height=15, filename) standardGeneric("plotGenesInClusters"))
setMethod("plotGenesInClusters",
          signature="list",
          definition=function(rankedDynamics, traj.name, num.genes, columns, width, height, filename){

            trajectory <- rankedDynamics[[traj.name]]
            traj.plots <- list()
            for(cl.name in unique(as.vector(rankedDynamics[[traj.name]]$cluster))){
              x <- trajectory[which(trajectory$cluster == cl.name),]
              x <- as.vector(x[1:num.genes,'genes'])
              g <- plottrajectory(cellrouter, cellrouter@pathsinfo$path[[traj.name]], x, width=4, height=1.5, filename='results/test.pdf')
              g <- g + ggtitle(paste('cluster',cl.name,sep='_'))
              g <- g + guides(col=guide_legend(direction="vertical",keywidth = 0.5, keyheight = 0.85))
              g <- g + theme(legend.text = element_text(colour="black", size=7))
              traj.plots[[paste('cluster',cl.name,sep='_')]] <- g
            }
            pdf(file = filename, width=width, height=height)
            multiplot(plotlist = traj.plots,cols=columns)
            dev.off()
          }
)
#' Compute a correlation of each gene along each trajectory
#' @param object CellRouter object
#' @param type Correlation method: Spearman or Perason correlation
#' @export
setGeneric("correlationPseudotime", function(object, type) standardGeneric("correlationPseudotime"))
setMethod("correlationPseudotime",
          signature="CellRouter",
          definition=function(object, type){
            pathsInfo <- object@pathsinfo
            genelist <- object@genes.trajectory
            correlations <- list()

            #ndata <- object@ndata
            print('computing correlation with the pseudotime')
            for(path in names(pathsInfo[['distr']])){
              cat(path, '\n')
              x <- pathsInfo[['pseudotime']][[path]]
              genes <- intersect(genelist, rownames(pathsInfo[['distr']][[path]]))
              #genes <- intersect(genelist, pathsInfo[['distr']][[path]])
              #cells <- object@pathsinfo$path[[path]]
              #tmp.ndata <- ndata[,cells]
              for(gene in genes){
                y <- as.numeric(pathsInfo[['distr']][[path]][gene,])
                #y <- as.numeric(tmp.ndata[gene, ])
                if(type == 'slope'){
                  df <- data.frame(x=x, y=y)
                  df <- lm(y ~ x, data=df)
                  correlations[[path]][[gene]] <- df$coefficients[2]
                }else if(type == 'pearson'){
                  cors <- suppressWarnings(cor.test(x, y, method="pearson"))
                  correlations[[path]][[gene]] <- as.numeric(cors$estimate)
                }else{
                  cors <- suppressWarnings(cor.test(x, y, method="spearman"))
                  correlations[[path]][[gene]] <- as.numeric(cors$estimate)
                }

              }
            }
            object@correlation <- correlations

            return(object)
          }
)


#' Find top genes more highly correlated with pseudotime
#' @param object Cellrouter object
#' @param max.quantile quantile to select positively correlated genes
#' @param min.quantile quantile to select negatively correlated genes
#' @export
setGeneric("topGenes", function(object, max.quantile, min.quantile) standardGeneric("topGenes"))
setMethod("topGenes",
          signature="CellRouter",
          definition=function(object, max.quantile, min.quantile){
            slopes <- object@correlation
            geneTrends <- list()
            plots <- list()

            for(path in names(slopes)){
              threshold.up <- quantile(slopes[[path]], max.quantile, na.rm=TRUE)
              threshold.down <- quantile(slopes[[path]], min.quantile, na.rm=TRUE)

              x <- slopes[[path]]
              genesUP <- x[which(x > threshold.up)]
              genesDOWN <- x[which(x < threshold.down)]
              geneTrends[['up']][[path]] <- genesUP[order(genesUP, decreasing=TRUE)]
              geneTrends[['down']][[path]] <- genesDOWN[order(genesDOWN)]
            }
            object@top.correlations <- geneTrends
            return(object)
          }
)

#' Plot mean kinetic profile for each transcritopnal cluster
#' @param object CellRouter object
#' @param show Selected trajectories
#' @param width width
#' @param height height
#' @param columns number of columns in the figure generated
#' @param file filename
#' export
#'
setGeneric("plotClusterHeatmap", function(object, show, width, height, columns, file) standardGeneric("plotClusterHeatmap"))
setMethod("plotClusterHeatmap",
          signature="CellRouter",
          definition=function(object, show, width, height, columns, file){

            clusters <- object@clusters[show]
            plots <- list()
            for(trajectory in names(clusters)){
              exprs <- clusters[[trajectory]][['exprs']]
              clustering <- clusters[[trajectory]][['clustering']]
              df <- data.frame(matrix(0, nrow=max(clustering), ncol=501))
              for(cluster in 1:max(clustering)){
                x <- clustering[which(clustering==cluster)]
                xx <- exprs[names(x),]
                xxx <- apply(xx, 2, mean)
                df[cluster, ] <- rescale(xxx, c(0,1))
                #df[cluster, ] <- xxx
              }
              colors <- c("navy","yellow","red")
              matrix <- df
              matrix$cluster <- rownames(matrix)
              matrix.m <- melt(matrix, id.var="cluster")
              matrix.m$cluster <- factor(rownames(df), levels=rev(rownames(df)))
              g <- ggplot(matrix.m, aes(variable, cluster)) + geom_tile(aes(fill = value)) +
                scale_fill_gradientn("", colours=colors) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
                theme(legend.position="right", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
                      axis.title.y = element_text(size = rel(0.3), angle = 90), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.text.x=element_blank(),axis.ticks=element_blank(),
                      panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) + ggtitle(trajectory)
              plots[[trajectory]] <- g
            }

            pdf(file=file, width=width, height=height)
            multiplot(plotlist = plots, cols=columns)
            dev.off();
            multiplot(plotlist = plots, cols=columns)
          }
)

#Does it need threshod parameter?
#' Plot heatmap of selected genes along selected trajectories
#' @param object CellRouter object
#' @param paths Selected trajectories
#' @param genelist Selected genes
#' @param threshold threshold
#' @param width width
#' @param height height
#' @param dir directory
#' @export

setGeneric("plotPathHeatmap", function(object, paths, genelist, threshold=2, width, height, dir) standardGeneric("plotPathHeatmap"))
setMethod("plotPathHeatmap",
          signature="CellRouter",
          definition=function(object, paths, genelist, threshold, width, height, dir){

            #plotPathHeatmap <- function(corsPaths, pathsInfo, graph, num_genes, width, height, dir){
            corsPaths <- object@correlation
            pathsInfo <- object@pathsinfo
            sampTab <- object@sampTab

            for(path in paths){
              genelist2 <- intersect(genelist, rownames(pathsInfo$distr[[path]]))
              tmpexpr <- pathsInfo$distr[[path]][genelist2,]
              andf <- data.frame(sampTab[pathsInfo$path[[path]], 'community',])
              rownames(andf) <- pathsInfo$path[[path]]
              colnames(andf) <- c('subpopulation')

              target_colors <- unique(sampTab[pathsInfo$path[[path]], 'colors',])
              names(target_colors) <- unique(andf$subpopulation)
              ann_colors = list(
                subpopulation = target_colors
              )
              from <- sapply(strsplit(path, split='.', fixed=TRUE), function(x){x[1]})
              to <- sapply(strsplit(path, split='.', fixed=TRUE), function(x){x[2]})
              title <- paste('Transition ', from, ' ', to, sep='')
              file <- paste(dir, 'heatmap_top_', path, '.pdf', sep='')
              labels <- sapply(strsplit(rownames(tmpexpr), split='__', fixed=TRUE), function(x){x[1]})
              #pheatmap(center_with_threshold(tmpexpr, threshold), cluster_rows = FALSE,
              #pheatmap(tmpexpr, cluster_rows = TRUE,
              pheatmap(tmpexpr, cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       annotation_col = andf, annotation_colors=ann_colors,
                       show_colnames = FALSE, border=FALSE, main=title, filename = file,
                       width = width, height=height, labels_row = labels)

              #to show plot in the report
              pheatmap(tmpexpr, cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       annotation_col = andf, annotation_colors=ann_colors,
                       show_colnames = FALSE, border=FALSE, main=title, labels_row = labels)

            }
          }
)


#' Plot genes along a trajectory
#' @param object CellRouter object
#' @param trajectory Selected trajectories
#' @param geneList Gene list
#' @param columns Number of columns in the figure generated
#' @param width width
#' @param height height
#' @param filename filename
#' @import stats
#' @export

setGeneric("plottrajectory", function(object, trajectory, geneList, columns=5, width=23, height=15, filename) standardGeneric("plottrajectory"))
setMethod("plottrajectory",
          signature="CellRouter",
          definition=function(object, trajectory, geneList,  columns, width, height, filename){

            library(smoother)

            plots <- list()
            x_axis <- 1:length(trajectory) #trajectory
            for(gene_id in geneList){
              y_axis <- as.numeric(object@ndata[gene_id, trajectory])
              lo <- loess(y_axis~x_axis)
              xl <- seq(min(x_axis),max(x_axis), (max(x_axis) - min(x_axis))/1000)
              y_axis <- predict(lo,xl)
              y_axis <- rescale(y_axis, newrange = c(0,1))

              df <- data.frame(cells=1:length(y_axis), Expression=as.numeric(y_axis))
              df$gene <- gene_id
              df$cells <- factor(df$cells, levels=df$cells)
              num_subpops <- length(unique(df$population))
              plots[[gene_id]] <- df
            }
            tables <- do.call(rbind, plots)
            labels <- x <- sapply(strsplit(as.vector(tables$gene), split='__', fixed=TRUE), function(x){x[1]})
            tables$gene <- labels
            tables$Expression <- rescale(tables$Expression, newrange = c(0,1))
            g1 <- ggplot(tables, aes(x=cells, y=Expression, group=gene, colour=gene)) +
              theme_bw() + geom_line(size=1) + xlab('CellRouter trajectory') +
              guides(col=guide_legend(direction="vertical")) + #, nrow = 2
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    axis.text.x=element_blank(), axis.ticks=element_blank(),
                    legend.position = "right",
                    panel.border = element_blank()) +
              theme(axis.line.x = element_line(color="black", size = 0.5),
                    axis.line.y = element_line(color="black", size = 0.5))+
              scale_color_manual("", values=rainbow(length(geneList)))
              #scale_color_brewer("", palette = 'Set1')

            pdf(filename, width=width, height=height)
            print(g1)
            dev.off()
            print(g1)

            g1
          }
)
#' Plot genes along a trajectory
#' @param object CellRouter object
#' @param trajectories Selected trajectories
#' @param genelist Gene list
#' @param rescale Whether to rescale curves simultaneously or individualy
#' @param columns Number of columns in the figure generated
#' @param width width
#' @param height height
#' @param filename filename
#' @export

setGeneric("plottrajectories", function(object, trajectories, geneList, rescale, columns=5, width=23, height=15, filename) standardGeneric("plottrajectories"))
setMethod("plottrajectories",
          signature="CellRouter",
          definition=function(object, trajectories, geneList, rescale, columns, width, height, filename){

            library(smoother)

            plotlist <- list()
            alldfs <- data.frame()
            for(t in trajectories){
              plots <- list()
              trajectory <- object@pathsinfo$path[[t]]
              x_axis <- 1:length(trajectory) #trajectory
              geneList2 <- intersect(geneList, rownames(object@pathsinfo$distr[[t]]))
              for(gene_id in geneList2){
                #y_axis <- as.numeric(object@ndata[gene_id, trajectory])
                y_axis <- as.numeric(object@pathsinfo$distr[[t]][gene_id, ])
                lo <- loess(y_axis~x_axis)
                xl <- seq(min(x_axis),max(x_axis), (max(x_axis) - min(x_axis))/1000)
                y_axis <- predict(lo,xl)
                if(rescale){
                  y_axis <- rescale(y_axis, newrange = c(0,1))
                }

                df <- data.frame(cells=1:length(y_axis), Expression=as.numeric(y_axis))
                df$gene <- gene_id
                df$cells <- factor(df$cells, levels=df$cells)
                num_subpops <- length(unique(df$population))
                plots[[gene_id]] <- df
              }
              tables <- do.call(rbind, plots)
              labels <- x <- sapply(strsplit(as.vector(tables$gene), split='__', fixed=TRUE), function(x){x[1]})
              tables$gene <- labels
              if(!rescale){
                tables$Expression <- rescale(tables$Expression, newrange = c(0,1))
              }
              tables$trajectory <- t
              alldfs <- rbind(alldfs, tables)
            }
            pdf(filename, width=width, height=height)
            g1 <- ggplot(alldfs, aes(x=cells, y=Expression, group=gene, colour=gene)) +
              theme_bw() + geom_line(size=1) + xlab('CellRouter trajectory') +
              guides(col=guide_legend(direction="vertical")) + #, nrow = 2
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    axis.text.x=element_blank(), axis.ticks=element_blank(),
                    legend.position = "right",
                    panel.border = element_blank(),
                    strip.background = element_rect(colour="white", fill="white")) +
              theme(axis.line.x = element_line(color="black", size = 0.5),
                    axis.line.y = element_line(color="black", size = 0.5))+
              scale_color_manual("", values=rainbow(length(geneList))) +
              facet_wrap(~trajectory, ncol = columns, scales='free')
            #scale_color_brewer("", palette = 'Set1')
            #plotlist[[t]] <- g1
            #print(g1)
            #multiplot(plotlist = plotlist, cols = columns)
            print(g1)
            dev.off()

            plot(g1)

          }
)

#' Plot genes along selected trajectories showing each single-cell in the trajectory and curve fit
#' @param object CellRouter object
#' @param paths Selected trajectories
#' @param genelist Gene list
#' @param columns Number of columns in the figure generated
#' @param width width
#' @param height height
#' @param filename filename
#' @export

setGeneric("plotpaths", function(object, paths, genelist, columns=5, width=23, height=15, file_prefix) standardGeneric("plotpaths"))
setMethod("plotpaths",
          signature="CellRouter",
          definition=function(object, paths, genelist, columns=5, width=23, height=15, file_prefix){

            path_distr <- cellrouter@pathsinfo
            sampTab <- object@sampTab
            cellDistances <- object@graph$similarity_matrix
            #paths <- names(path_distr$distr)

            for(path in paths){
              print(path)
              plots <- list()
              for(gene_id in genelist){
                if(gene_id %in% rownames(path_distr$distr[[path]])){
                  x_axis <- 1:length(path_distr$distr[[path]][gene_id,])
                  y_axis <- path_distr$distr[[path]][gene_id,]
                  df <- data.frame(cells=x_axis, Expression=as.numeric(y_axis), population=sampTab[names(y_axis),'population'])
                  df$cells <- factor(df$cells, levels=df$cells)
                  df$population <- factor(df$population, levels=unique(df$population))
                  num_subpops <- length(unique(df$population))
                  colors <- sampTab[names(y_axis), 'colors']
                  names(colors) <- as.vector(df$population)

                  g1 <- ggplot(df, aes(x=cells, y=Expression, group=1, colour=population)) +
                    geom_point(size=5) + stat_smooth() + theme_bw() + ggtitle(paste(gene_id, path, sep='--')) +
                    scale_colour_manual(values = colors) + xlab("Trajectory") +
                    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.text.x=element_blank(),
                          panel.background=element_blank(),
                          plot.background=element_blank(),
                          panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1))

                  plots[[gene_id]] <- g1
                }else{
                  cat(gene_id, ' is not regulated through trajectory: ', path, '\n')
                }
              }

              if(!is.null(cellDistances)){
                #path similarity matrix
                matrix <- as.data.frame(cellDistances[as.vector(path_distr$path[[path]]), as.vector(path_distr$path[[path]])])
                matrix$cells <- rownames(matrix)
                matrix.m <- melt(matrix, id.var="cells")
                matrix.m$cells <- factor(as.vector(path_distr$path[[path]]), levels=as.vector(path_distr$path[[path]]))

                g2 <- ggplot(matrix.m, aes(cells, variable)) + geom_tile(aes(fill = value)) +
                  scale_fill_gradient2(low = "blue",mid="white", high = "red") + theme_bw() + xlab("Cells") + ylab("Cells") +
                  theme(legend.position="none", axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        #axis.title.y = element_text(size = rel(0.3), angle = 90),
                        axis.ticks=element_blank(),
                        #axis.title.y=element_blank(),
                        panel.background=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        plot.background=element_blank(),
                        panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
                  ggtitle(path)
                plots[[path]] <- g2
              }
              filename <- paste(file_prefix, "_path_distribution_", path, ".pdf", sep="")
              pdf(file=filename, width=width, height=height)
              multiplot(plotlist = plots, cols=columns)
              dev.off();
              multiplot(plotlist = plots, cols=columns)
            }
          }
)

#' Plot reduced dimension space
#' @param object CellRouter object
#' @param reduction.type The dimension reduction space to be used: pca, tsne, dc of diffusion components and custom
#' @param dims.use Dimensions to use
#' @param annotation Column in the metadata table to annotate the figure
#' @param annotation.color Correspondinf color column for the annotation
#' @param showlabels Show labels in the plot.
#' @param width width
#' @param height height
#' @param filename filename
#' @export

plotReducedDimension <- function(object, reduction.type='tsne', dims.use=c(1,2), annotation, annotation.color, showlabels, width, height, filename){
  matrix <- slot(object, reduction.type)$cell.embeddings[,dims.use]
  scores <- as.data.frame(matrix)
  scores <- scores[rownames(object@sampTab),]
  colnames(scores) <- c('x', 'y')
  scores$group <- factor(as.vector(object@sampTab[[annotation]]), levels=unique(as.vector(object@sampTab[[annotation]])))
  colors <- unique(object@sampTab[[annotation.color]])
  names(colors) <- unique(as.vector(object@sampTab[[annotation]]))

  g1 <- ggplot(scores,aes(x = x, y=y, colour=group)) + theme(legend.position='right') + geom_point(size=0.5)
  
  pdf(file=filename, width=width, height=height)
  if(showlabels==TRUE){
    centers <- scores %>% dplyr::group_by(group) %>% summarize(x = median(x = x), y = median(x = y))
    g1 <- g1 + geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) + geom_text(data = centers, mapping = aes(label = group), size = 5, colour='black')
  }

  if(reduction.type == 'tsne'){
    xlab <- paste('tSNE ', dims.use[1], sep=' ')
    ylab <- paste('tSNE ', dims.use[2], sep=' ')
  }else if(reduction.type == 'pca'){
    xlab <- paste('PC', dims.use[1], sep='')
    ylab <- paste('PC', dims.use[2], sep='')
  }else if(reduction.type == 'dc'){
    xlab <- paste('DC', dims.use[1], sep='')
    ylab <- paste('DC', dims.use[2], sep='')
  }else{
    xlab <- 'Dim 1'
    ylab <- 'Dim 2'
  }

  g1 <- g1 + theme_bw() +
    ggtitle('')  +  scale_color_manual("", values=colors) + #scale_color_brewer("", palette = 'Paired') +
    xlab(xlab) + ylab(ylab) +
    theme(panel.border = element_rect(fill = NA, colour = "white"),
          strip.background = element_blank()) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) +
    guides(col=guide_legend(direction="vertical", keywidth = 0.75, keyheight = 0.85, override.aes = list(size=3)))
  print(g1)
  dev.off()
  plot(g1)
}

#' Show gene expression in the space of reduced dimensionality
#' @param object CellRouter object
#' @param genelist Gene list to show
#' @param reduction.type The dimension reduction space to be used: pca, tsne, dc of diffusion components and custom
#' @param threshold Threshold to rescale gene expression
#' @param dims.use Dimensions to use
#' @param columns Number of columns in the figure generated
#' @param width width
#' @param height height
#' @param filename filename
#' @export

setGeneric("plotDRExpression", function(object, genelist, reduction.type='tsne',threshold=2, dims.use=c(1,2), columns=5, width=23, height=15, filename) standardGeneric("plotDRExpression"))
setMethod("plotDRExpression",
          signature="CellRouter",
          definition=function(object, genelist, reduction.type,threshold, dims.use, columns=5, width=23, height=15, filename){

            #plotDRExpression <- function(expDat, matrix, geneList, width=10, height=3.5, num_columns=2, filename){
            matrix <- as.data.frame(slot(object, reduction.type)$cell.embeddings[,dims.use])
            plots <- list()
            scores <- matrix
            colnames(scores) <- c('Dim_1', 'Dim_2')

            if(reduction.type == 'tsne'){
              xlab <- paste('tSNE ', dims.use[1], sep=' ')
              ylab <- paste('tSNE ', dims.use[2], sep=' ')
            }else if(reduction.type == 'pca'){
              xlab <- paste('PC', dims.use[1], sep='')
              ylab <- paste('PC', dims.use[2], sep='')
            }else if(reduction.type == 'DC'){
              xlab <- paste('DC', dims.use[1], sep='')
              ylab <- paste('DC', dims.use[2], sep='')
            }else{
              xlab <- 'Dim 1'
              ylab <- 'Dim 2'
            }

            #x <- center_with_threshold(object@ndata[genelist, rownames(matrix)], 3)
            x <- object@ndata[genelist, rownames(matrix)]
            x <- as.data.frame(t(apply(x, 1, function(y){rescale(y, c(-1*threshold,threshold))})))
            #x <- center_with_threshold(object@ndata[genelist, rownames(matrix)], 2)
            #x <- x[genelist,]
            gc()
            dfs <- data.frame()
            for(gene in genelist){
              #expr <- object@ndata[gene,rownames(matrix)]
              expr <- x[gene,]
              scores$GENE <- as.numeric(expr)
              scores$gene <- gene
              dfs <- rbind(dfs, scores)
              #plots[[gene]] <- p1
            }
            dfs$gene <- factor(dfs$gene, levels=genelist)
            pdf(file=filename, width=width, height=height)
            #multiplot(plotlist = plots, cols=columns)
            p1 <- ggplot(dfs,aes(x = Dim_1, y=Dim_2, colour=GENE)) + geom_point(size=0.1) + theme_bw() +
              #scale_colour_gradientn(name = "Expression", colours=rev(rainbow(4))) +
              scale_colour_gradientn("Relative expression", colours=c("midnightblue","dodgerblue3","white", "darkorange2", "yellow")) +
              #scale_colour_gradient2("Relative expression",low="midnightblue", high="darkorange2") +
              ylab(ylab) + xlab(xlab) +
              theme(panel.border = element_rect(fill = NA, colour = "white"),
                    strip.background = element_blank()) +
              theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(), panel.background = element_blank()) +
              theme(legend.position="bottom",
                    strip.background = element_rect(colour="white", fill="white")) +
              guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
                     size = guide_legend(title.position="top", title.hjust = 0.5)) +

              #facet_wrap(~gene, ncol = columns, scales='free')
              facet_wrap(~gene, ncol = columns)

            print(p1)
            dev.off();
            #multiplot(plotlist = plots, cols=columns)
            plot(p1)
          }
)

#' Plot a network module enriched around top transcriptional regulators of cell-fate transitions
#' @param x GRN scores as calculated by the grnscores function
#' @param ggrn Gene regulatory network
#' @param genelist Vector of gene names containing the transcriptional regulators from which the network modules will be identified
#' @param nrow Number of rows in the figure generated
#' @param ncol Number of columns in the figure generated
#' @param width width
#' @param height height
#' @param filename filename
#' @export

regulatornetwork <- function(x, ggrn, genelist, nrow,
                             ncol, width, height, filename){
  nets <- list()
  nets2 <- list()
  for(t in names(x)){
    for(n in genelist){
      genes <- x[[t]]$targets[[n]]
      if(length(genes) > 0){
        genes <- c(n, genes)
        rgrn <- igraph::induced.subgraph(ggrn,
                                         vids=unlist(igraph::neighborhood(graph=ggrn,order=0,
                                                                          nodes=genes)))
        remove <- V(rgrn)$name[igraph::degree(rgrn)==0]
        rgrn <- igraph::delete.vertices(rgrn, remove)
        V(rgrn)$color <- cellrouter@correlation[[t]][V(rgrn)$name]
        V(rgrn)$size <- log(igraph::degree(rgrn)+1)

        g <- fortify(rgrn)
        g$network <- n
        g$transition <- t
        name <- paste(t,n,sep='_')
        nets2[[name]] <- rgrn
        nets[[name]] <- g
      }
    }
  }
  l <- do.call(rbind, nets)
  l$network <- factor(l$network, levels=regulators)
  l$label[which(l$type == 'Target')] <- ""
  ##choose a threshold for each network independently.
  q <- quantile(l$color, 0.5)
  xxx <- l[which(l$type == 'Regulator'),]
  xxx <- as.vector(xxx[which(xxx$color < q), 'label'])
  xxx <- xxx[!(xxx %in% genelist)]
  l$label[l$label %in% xxx] <- ""

  set.seed(1)
  pdf(file=filename, width=width, height=height)
  g <- ggplot(data = l, aes(from_id = from, to_id = to)) +
    #geom_net(aes(colour = props, size=size, label=label), layout.alg = "kamadakawai",
    geom_net(aes(colour = color, size=size, label=label), layout.alg = "kamadakawai",
             labelon = TRUE, vjust = -0.6, ecolour = "grey60",
             directed =FALSE, fontsize = 2, ealpha = 0.1, labelcolour = 'black',
             fiteach=T, arrowsize = 1) +
    xlim(c(-0.1, 1.1)) + ylim(-0.1, 1.1) +
    theme_net() + theme(legend.position = "bottom") +
    #scale_colour_manual(values = c('lightgreen', 'red','orange')) +
    #scale_colour_brewer("", palette = 'Paired') +
    scale_colour_gradientn(colours = c('cyan', 'white', 'red')) +
    facet_grid(network~transition, scales='free') +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          strip.background = element_rect(colour="black"))
  #scale_colour_brewer("", palette = 'Set2')
  print(g)
  dev.off()

  plot(g)
  nets2
}


#' Pathway enrichment analysis on genes up or down-regulated along trajectories
#' @param object CellRouter object
#' @param names Selected trajectories
#' @param cc List of cell cycle genes or other gene lists to be removed from the genes up or down-regulated in each trajectory
#' @param species specifies: Hs for Homo Sapiers or Mm for Mus Musculus
#' @param annotation organism-specific annotations: 'mouse','org.Mm.eg.db' or 'human' 'org.Hs.eg.db'
#' @param ids table containing mappings between gene identifiers
#' @export

setGeneric("pathwayenrichment", function(object, names, cc=NULL, species, annotation, ids) standardGeneric("pathwayenrichment"))
setMethod("pathwayenrichment",
          signature="CellRouter",
          definition=function(object, names, cc=NULL, species, annotation, ids){

            ##up-regulated
            upregulated <- list()
            geneNames <- lapply(object@top.correlations$up, names)
            if(!is.null(cc)){
              cat('removing genes...\n')
              geneNames <- lapply(geneNames, function(x){setdiff(x, cc)}) #remove a selected gene set
            }
            geneList <- lapply(geneNames, function(x){convertIDs(ids, x, from='external_gene_name', to="entrezgene")})
            geneList <- lapply(geneList, names)
            geneList <- lapply(geneList, function(x){x[!is.na(x)]})
            print('pathway enrichment for up-regulated genes')
            ck1 <- compareCluster(geneCluster = geneList[names], fun = "enrichPathway", organism = species, pvalueCutoff = 0.05, readable=T)
            #print('Here1')
            #ck2 <- compareCluster(geneCluster = geneList[names], fun = "enrichKEGG", organism = 'mouse', pvalueCutoff = 0.05)
            ck3 <- compareCluster(geneCluster = geneList[names], fun = "enrichGO", ont='BP', OrgDb=annotation, pvalueCutoff = 0.05, readable=T)
            #print('Here1')
            upregulated[['REACTOME']] <- ck1
            #upregulated[['KEGG']] <- ck2
            upregulated[['GOBP']] <- ck3

            #down-regulated
            downregulated <- list()
            geneNames <- lapply(object@top.correlations$down, names)
            geneList <- lapply(geneNames, function(x){convertIDs(ids, x, from='external_gene_name', to="entrezgene")})
            geneList <- lapply(geneList, names)
            print('pathway enrichment for down-regulated genes')
            ck1 <- compareCluster(geneCluster = geneList[names], fun = "enrichPathway", organism = species, pvalueCutoff = 0.05, readable=T)
            #ck2 <- compareCluster(geneCluster = geneList[names], fun = "enrichKEGG", organism = 'mouse', pvalueCutoff = 0.05)
            ck3 <- compareCluster(geneCluster = geneList[names], fun = "enrichGO", ont='BP', OrgDb=annotation, pvalueCutoff = 0.05, readable=T)
            #downregulated[['REACTOME']] <- ck1
            #downregulated[['KEGG']] <- ck2
            downregulated[['GOBP']] <- ck3

            object@pathwayenrichment <- list("UP"=upregulated, "DOWN"=downregulated)

            return(object)
          })

#' Plot pathway enrichment analysis
#' @param object CellRouter object
#' @param pathway Which pathway database: GO, Reactome.
#' @param numpathways Number of pathways to show for each trajectory
#' @param logTransform Log-transform the p-values
#' @param width width
#' @param height height
#' @param filename filename
#' @export
setGeneric("pathwaycluster", function(object, pathway, numpathways=5, logTransform, width, height, filename) standardGeneric("pathwaycluster"))
setMethod("pathwaycluster",
          signature="CellRouter",
          definition=function(object, pathway, numpathways, logTransform, width, height, filename){
            pathsInfo <- object@pathsinfo
            #pathwaycluster <- function(pathsInfo, ck, num_pathways=5, crows, ccols, filename, width, height){
            #library(tidyr)
            #include <- object@pathwayenrichment$UP$GOBP@compareClusterResult$Cluster %in% names

            #b <- pathway@compareClusterResult[include,]
            b <- pathway@compareClusterResult
            write.csv(b, paste(filename, '.csv', sep=''))
            #b <- b[complete.cases(b$Description),]
            b2 <- lapply(split(b, as.vector(b$Cluster)),
                         function(x){x[order(as.vector(x$p.adjust), decreasing=FALSE),][1:numpathways,]})
            b <- do.call(rbind, b2)
            b <- b[complete.cases(b),]
            #b$Description <-  make.unique(as.vector(b$Description), sep='_')
            pathways <- data.frame(matrix(1, nrow=length(unique(b$Description)), ncol=length(unique(b$Cluster))))
            rownames(pathways) <- as.vector(unique(b$Description))
            colnames(pathways) <- as.vector(unique(b$Cluster))
            for(c in colnames(pathways)){
              for(d in rownames(pathways)){
                x <- b[which(as.vector(b$Cluster) == c & as.vector(b$Description) == d),]
                if(nrow(x) == 1){
                  pathways[as.character(x$Description), as.character(x$Cluster)] <- as.numeric(x$p.adjust)
                }
              }
            }
            if(logTransform){
              pathways <- -log(pathways)
            }
            pathways <- pathways[,order(colnames(pathways))]

            ann <- pathsInfo$path_data[colnames(pathways),]
            source_colors <- unique(ann$source_color)
            names(source_colors) <- unique(ann$source_population)
            target_colors <- unique(ann$target_color)
            names(target_colors) <- unique(ann$target_population)
            ann_colors = list(
              source_population = source_colors,
              target_population = target_colors
            )

            #hmcols <- colorRampPalette(c("midnightblue","dodgerblue3","white","orange", "darkred"))(100)
            pheatmap(pathways, border=TRUE, border_color = 'gray',
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     cellwidth=8, cellheight = 8, annotation_col = ann[,c(4,6)],
                     annotation_colors=ann_colors, clustering_method = 'ward.D', filename=filename)

            pheatmap(pathways, border=TRUE, border_color = 'gray',
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     cellwidth=8, cellheight = 8, annotation_col = ann[,c(4,6)],
                     annotation_colors=ann_colors, clustering_method = 'ward.D')

            pathways
          }
)

#' Convert gene identifiers
#' @param ids Table containing mappings betweej gene identifiers
#' @param ens identifier
#' @param from Gene identifier
#' @param to Gene identifier
convertID <- function(ids, ens, from, to){
  symbol <-  ids[which(ids[,from]==ens), to]
  if(length(symbol) == 0){
    ens
  }else{
    symbol[1]
  }
}
#' Convert gene lists from one identifier to another
#' @param ids Table containing mappings betweej gene identifiers
#' @param geneList Gene list to be converted
#' @param from Gene identifier
#' @param to Gene identifier

convertIDs <- function(ids, geneList, from, to){
  genes <- vector()
  for(gene in geneList){
    g <- convertID(ids, gene, from, to)
    genes <- append(genes, g)
  }
  names(geneList) <- genes
  genes <- genes[!is.na(genes)]
  geneList
}

#' Cluster genes in num.clusters kinetic trends
#' @param fits fits
#' @param num.clusters num.clusters
clustergenes <- function(fits, num.clusters){
  expr_matrix <- fits
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0,]
  expr_matrix <- t(scale(t(log10(expr_matrix))))
  expr_matrix <- expr_matrix[is.nan(rowSums(expr_matrix)) == FALSE,]
  expr_matrix[is.na(expr_matrix)] <- 0

  var <- apply(expr_matrix, 1, var)
  var <- var[var > 0]
  n <- as.dist((1 - cor(t(expr_matrix[names(var),])))/2)
  clusters<-pam(n,num.clusters)
  clusters$exprs <- expr_matrix

  clusters
}



##Create colors from numeric values
#' Helper function
#' @param x values to be converted to colors
range01 <- function(x)(x-min(x))/diff(range(x))

#' Helper function
#' @param x vaues to be converted to colors

cRamp <- function(x){
  cols <- colorRamp(rev(rainbow(4)))(range01(x))
  #cols <- colorRamp(c('deepskyblue3', 'white', 'red'))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}
#' Helper function
#' @param x values to be converted to colors
cRamp2 <- function(x){
  #cols <- colorRamp(rev(rainbow(4)))(range01(x))
  cols <- colorRamp(c('white', "goldenrod1","darkorange2"))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}

#' Helper function
#' @param x values to be converted to columns
#' @param num_colors Number of colors to interpolate

cRampClust <- function(x, num_colors){
  library(RColorBrewer)
  #cols <- colorRamp(rev(rainbow(num_colors)))(range01(x))
  cols <- colorRamp(rev(brewer.pal(num_colors, "Paired")))(range01(x))
  #cols <- colorRamp(rev(brewer.pal(num_colors, "Set1")))(range01(x))
  #cols <- colorRamp(c('deepskyblue3', 'lightblue', 'orange','darkred'))(range01(x))
  #color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  #cols <- color[sample(num_colors)]
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
  #cols
}

#' Center the data
#' @param center_data data to be processed
#' @param threhsold Max value present in center_data

center_with_threshold <- function(center_data, threshold){
  # Center data (automatically ignores zeros)
  center_data <- center_data - rowMeans(center_data, na.rm=TRUE)
  # Cap values between threshold and -threshold and recenter
  center_data[center_data > threshold] <- threshold
  center_data[center_data < (-1 * threshold)] <- -1 * threshold
  #center_data <- center_data - rowMeans(center_data, na.rm=TRUE)
  return(center_data)
}


### GRN reconstruction code: from CellNet (Cahan et al., Cell 2014) #############
#' GRN reconstruction
#' @param corrMat correlation matrix
mat_zscores<-function(corrMat){
  z_row<-scale(t(corrMat))**2;
  cat(dim(z_row),"\n");
  z_col<-scale(corrMat)**2;
  cat(dim(z_col),"\n");
  ans<-sqrt( z_row+z_col);
  ans;
}

#' Create global gene regulatory network
#' @param expr expression matrix
#' @param tfs list of transcriptional regulators
#' @param threshold zscore threshold to identify regulatory interactions

globalGRN <- function(expr, tfs, threshold){
  tfs <- intersect(tfs, rownames(expr))
  corrAll <- abs(cor(t(expr)))
  zscores <- mat_zscores(corrAll)
  zscores <- zscores[,tfs]
  grnTable <- extractRegulators(zscores, corrAll, rownames(expr), threshold)

  grnTable
}

#' Helper function to extract transcriptional regulators
#' @param zscores zscores
#' @param corrMatrix correlation matrix
#' @param genes genes
#' @param threshold threshold

extractRegulators <- function(zscores, corrMatrix, genes, threshold){
  targets <- vector()
  regulators <- vector()
  zscoresX <- vector()
  correlations <- vector()

  targets <- rep('', 1e6)
  regulators <- rep('', 1e6)
  zscoresX <- rep('', 1e6)
  correlations <- rep('', 1e6)

  i <- 1
  j <- 1
  for(target in genes){
    x <- zscores[target,]
    regs <- names(which(x > threshold))
    if(length(regs) > 0){
      zzs <- x[regs]
      ncount <- length(zzs)
      corrs <- corrMatrix[target, regs]
      j <- i + ncount - 1
      targets[i:j] <- rep(target, ncount)
      regulators[i:j] <- regs
      zscoresX[i:j] <- zzs
      correlations[i:j] <- corrs
      i <- j + 1
    }
  }
  targets <- targets[1:j]
  regulators <- regulators[1:j]
  zscoresX <- zscoresX[1:j]
  correlations <- correlations[1:j]

  data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations)
}

#' Helper function to create an igraph object for the gene regulatory network
#' @param grnTab Table containing gene regulatory relationships
#' @param simplify Simplify the graph. For example, remove loops
#' @param directed Create a directed or undirected graph
#' @param weights weights
#' @import igraph

ig_tabToIgraph<-function# return a iGraph object
(grnTab, ### table of TF, TF, maybe zscores, maybe correlations
 simplify=FALSE, # failed when iranges is loaded...
 directed=FALSE,
 weights=TRUE
){
  tmpAns<-as.matrix(grnTab[,c("TF", "TG")]);
  regs<-as.vector(unique(grnTab[,"TF"]));
  ###cat("Length TFs:", length(regs), "\n");
  targs<-setdiff( as.vector(grnTab[,"TG"]), regs);

  ###  cat("Length TGs:", length(targs), "\n");
  myRegs<-rep("Regulator", length=length(regs));
  myTargs<-rep("Target", length=length(targs));

  types<-c(myRegs, myTargs);
  verticies<-data.frame(name=c(regs,targs), label=c(regs,targs), type=types);

  iG<-graph.data.frame(tmpAns,directed=directed,v=verticies);

  if(weights){
    #E(iG)$weight<-grnTab$weight;
    #E(iG)$weight<-as.numeric(as.character(grnTab$zscore));
    #print(range(as.numeric(as.character(grnTab$corr))))
    E(iG)$weight<-as.numeric(as.character(grnTab$corr));
    #print(range(E(iG)$weight))
  }

  if(simplify){
    #iG<-igraph::simplify(iG, remove.loops = TRUE, remove.multiple = FALSE);
    iG<-igraph::simplify(iG);
    #print(range(E(iG)$weight))
  }
  V(iG)$nEnts<-1;

  iG;
}
## Assign names to subpopulations
commToNames<-function(commObj, prefix){
  ans<-list();
  comms<-communities(commObj);
  for(i in seq(length(comms))){
    #nname<-paste(prefix,"_sn_",i,sep='');
    nname<-paste(prefix, "", i,sep='');
    #ans[[nname]]<-commObj$names[comms[[i]]];
    ans[[nname]]<-comms[[i]]; #new igraph version
  }
  ans;
}
## labels for subpopulations (for plotting)
nodeLabels <- function(sampTab, column){
  x <- sampTab
  x$label <- ''
  xx <- vector()
  names <- vector()
  #for(c in unique(as.vector(x$community))){
  for(c in unique(as.vector(x[[column]]))){
    #tmp <- x[which(x$community == c), ]
    tmp <- x[which(x[[column]] == c), ]
    #label <- c(unique(tmp$community), rep('', nrow(tmp)-1))
    label <- c(unique(tmp[[column]]), rep('', nrow(tmp)-1))
    tmp$label <- sample(label)
    xx <- append(xx, as.vector(tmp$label))
    names <- append(names, rownames(tmp))
  }
  names(xx) <- names
  xx
}

############### end of GRN code #######

#' Heler function to create heatmaps with gene signatures
#' @param ann_col Column annotations
#' @param ann_row Row annotations
#' @param order.columns Order of columns
#' @param order.rows Order of rows

getIndexes <- function(ann_col, ann_row, order.columns, order.rows){
  ann_col$ID <- as.vector(1:nrow(ann_col))
  ref_groups <- split(ann_col, as.factor(ann_col$population))
  ref_groups <- lapply(ref_groups, function(x){x$ID})
  #ref_groups <- ref_groups[order(nchar(names(ref_groups)))]
  ref_groups <- ref_groups[order.columns]

  ### Figure out where to draw lines between subtypes in the heatmap
  ref_seps <- c()
  i_cur_idx <- 0
  order_idx <- c()

  for(ref_grp in ref_groups){
    i_cur_idx <- i_cur_idx + length(ref_grp)
    ref_seps <- c(ref_seps, i_cur_idx)
    order_idx <- c(order_idx, ref_grp)
  }
  #ref_seps <- ref_seps[1:(length(ref_seps) - 1)]
  #ref_seps <- c(ref_seps, c(0,nrow(ann_col)))

  ### Figure out where to draw lines between gene signatures in the heatmap
  #genesDF <- data.frame(clusters=rep(seq(1:4), times=lapply(subtypeSigs, length)), ID=1:length(comb.subtype.sigs))
  #genesDF <- data.frame(clusters=graph$sampTab$population, ID=1:length(comb.subtype.sigs))
  ann_row$ID <- as.vector(1:nrow(ann_row))
  ref_groups <- split(ann_row, as.factor(ann_row$signature))
  ref_groups <- lapply(ref_groups, function(x){x$ID})
  ref_groups <- ref_groups[order.rows]
  #ref_groups <- ref_groups[order(nchar(names(ref_groups)))]

  ref_seps_c <- c()
  i_cur_cdx <- 0
  order_cdx <- c()
  for(ref_grp in ref_groups){
    i_cur_cdx <- i_cur_cdx + length(ref_grp)
    ref_seps_c <- c(ref_seps_c, i_cur_cdx)
    order_cdx <- c(order_cdx, ref_grp)
  }
  #ref_seps_c <- ref_seps_c[1:(length(ref_seps_c) - 1)]
  #ref_seps_c <- c(ref_seps_c, c(0,nrow(ann_row)))

  list(colsep=ref_seps, rowsep=ref_seps_c)

}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

#' Multiplot function
#' @param plotlist list of plots
#' @param file filename
#' @param cols number columns
#' @param layout layout

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  #plots <- c(list(...), plotlist) #comment this to plot several plots in list
  plots <- plotlist
  numPlots = length(plots)


  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Helper function to identify transcriptional regulators based on gene ontology annotations
#' @param species species: Hs for human or Mm for mouse

find_tfs<-function# find transcript factors
(species='Hs' # species abbreviation
){

  #cat("Loading gene annotations ...\n")
  suppressWarnings(suppressMessages(require(GO.db)));

  if(species=='Hs'){
    suppressWarnings(suppressMessages(require(org.Hs.eg.db)));
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }
  else{
    suppressWarnings(suppressMessages(require(org.Mm.eg.db)));
    egSymbols<-as.list(org.Mm.egSYMBOL);
    goegs<-as.list(org.Mm.egGO2ALLEGS);
  }

  goterms<-as.list(GOTERM);
  goids<-names(goegs);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts=='BP'];
  goids<-names(unlist(bps));

  cat("matching gene symbols and annotations")
  gobpList<-list();
  for(goid in goids){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gobpList[[goterm]]<-genes;
  }

  ### newHsTRs<-gobpList[['regulation of transcription, DNA-dependent']];
  regNames<-names(gobpList)[grep("regulation of transcription", names(gobpList))];
  trs<- unique(unlist(gobpList[regNames]));
  #cat("Regulation of transcription: ", length(trs),"\n");
  #cat(regNames, '\n')

  mfs<-onts[onts=='MF'];
  goidsMF<-names(unlist(mfs));

  gomfList<-list();
  for(goid in goidsMF){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gomfList[[goterm]]<-genes;
  }
  dbs<-gomfList[['DNA binding']];
  #dbs<-gomfList[['RNA binding']];
  #cat("DNA binding: ", length(dbs),"\n");
  #cat("RNA binding: ", length(dbs),"\n");
  sort(intersect(trs, dbs));
  #sort(dbs)
}
