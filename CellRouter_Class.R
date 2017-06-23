##### CellRouter Class #####
## add require statements here
## make clustering independent in CellRouter->Cluster memberships as input
require('reshape')
require('reshape2')
require('pheatmap')
require('clusterProfiler')
require('ReactomePA')
require('plotrix')
require('tsne')
require('igraph')

CellRouter <- setClass("CellRouter", slots=
                        c(expdata="data.frame", ndata="data.frame",
                          sampTab="data.frame", rdimension="data.frame", graph="list",
                          signatures="list", sources="character", targets="character",
                          directory="list", paths="data.frame", networks="list", 
                          genes.trajectory="character", pathsinfo="list",
                          dynamics="list", clusters="list", correlation="list",
                          top.correlations="list", davidenrichment="list", pathwayenrichment="list"))

### Function to validate CellRouter object
setValidity("CellRouter",
            function(object){
              msg <- NULL
              if(!is.data.frame(object@expdata)){
                msg <- c(msg, "expression data must be a data.frame")
              }else if(nrow(object@expdata) < 2){
                msg <- c(msg, "expression data must have more than 2 columns")
              }else if(sum(apply(is.na(object@expdata), 1, sum) > 0)){
                msn <- c(msg, "expression data must not have NAs")
              }#else if(sum(apply(object@expdata, 1, min) < 0)){
              #  msg <- c(msg, "negative values are not allowed in expression data")
              #}
              if(is.null(msg)){
                TRUE
              }else{
                msg
              }
            })

### boxplot distribution
boxplotGenes <- function(object, geneList, column, cols, width=10, height=5, filename){
  plots <- list()
  expDat <- object@ndata
  samples <- object@sampTab
  T0 <- expDat
  for(g in geneList){
    genes <- as.data.frame(t(T0[g,]))
    genes$gene <- g
    genes$conditions <- as.vector(samples[,column])
    genes.m <- melt(genes, id.var=c('gene',"conditions"))
    genes.m$conditions <- factor(genes.m$conditions, levels=unique(samples$population))
    
    p <- ggplot(genes.m, aes(x=conditions, y=value, fill=conditions)) + geom_boxplot(alpha=.9) + 
      theme_bw() + xlab("") + ylab("Gene expression") + theme(legend.position="none") +
      theme(axis.text.x = element_text(size=rel(1), angle=45, hjust=1)) + ggtitle(g) +
      scale_fill_manual("", values=unique(samples$colors))
    
    plots[[g]] <- p
  }
  #pdf(file=filename, width=10, height=7)
  pdf(file=filename, width=width, height=height)
  multiplot(plotlist = plots, cols=cols)
  dev.off();
}



########################### GRN reconstruction code: from CellNet #############
mat_zscores<-function# computes sqrt(zscore_row + zscore_col) .. see JJ Faith et al 2007
(corrMat # does it matter whether TF is row or col?
){
  z_row<-scale(t(corrMat))**2;
  cat(dim(z_row),"\n");
  z_col<-scale(corrMat)**2;
  cat(dim(z_col),"\n");
  ans<-sqrt( z_row+z_col);
  ans;  
}

globalGRN <- function(expr, tfs, threshold){
  tfs <- intersect(tfs, rownames(expr))
  corrAll <- abs(cor(t(expr)))
  zscores <- mat_zscores(corrAll)
  zscores <- zscores[,tfs]
  grnTable <- extractRegulators(zscores, corrAll, rownames(expr), threshold)
  
  grnTable
}

#extract transcriptional regulators
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

ig_scaleV<-function# return a vector of scaled sizes for a vector of verticies
(vals, # values associated with verticies.  Can be number of sub-net nodes (members) or degree (number of edges))
 sf=5, # scaling factor, so sf=5 means that the maximum vertix will have cex=5)
 minVal=2
){
  vals<-vals-min(vals); # shift such that m
  vals<-vals/max(vals); # scale such that range vals == 0,1
  minVal+(vals*sf); # scale to sf
}


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
    print(range(as.numeric(as.character(grn$corr))))
    E(iG)$weight<-as.numeric(as.character(grnTab$corr));   
    print(range(E(iG)$weight))
  }
  
  if(simplify){
    #iG<-igraph::simplify(iG, remove.loops = TRUE, remove.multiple = FALSE);
    iG<-igraph::simplify(iG);
    print(range(E(iG)$weight))
  }
  V(iG)$nEnts<-1;
  
  iG;
}



ig_convertSmall<-function# change igraph attributes so that it is suitable for plotting a small GRN
(iG, ###  in which regulators and targets are already specified
 rCol="#F46D43", # default color for regulator nodes
 tCol="#66C2A5",  # default color for target nodes
 vScale=1,
 exponent=2.6
){
  
  E(iG)$color<-rgb(0,0,.5,alpha=.05);
  
  V(iG)[type=='Regulator']$label<-V(iG)[type=='Regulator']$name;
  #V(iG)[type=='Target']$label<-"";
  V(iG)$size<-ig_scaleV( igraph::degree(iG), sf=vScale, minVal=4);
  V(iG)[type=='Regulator']$shape<-"csquare";
  V(iG)[type=='Regulator']$label.color<-"black";
  V(iG)[type=='Regulator']$label.cex<-.75;
  V(iG)[type=='Target']$label.color<-"black";
  V(iG)[type=='Target']$label.cex<-.5;
  V(iG)[type=='Target']$shape<-"circle";
  #V(iG)[type=='Target']$frame.color=NA;
  #V(iG)$label.degree<- -1*pi/2;
  #V(iG)$label.dist=.25;
  nRegs<-length(V(iG)[type=='Regulator']$name);
  nTargs<-length(V(iG)[type=='Target']$name);
  #rCols<-unlist(lapply(rep(rCol, nRegs), mp_makeTransparent ) );
  rCols<-rep(rCol, nRegs);
  tCols<-rep(tCol, nTargs);
  #tCols<-unlist(lapply(rep(tCol, nTargs), mp_makeTransparent ) );
  V(iG)[type=='Regulator']$color<-rCols;
  V(iG)[type=='Target']$color<-tCols;  
  #iG$layout<-layout.lgl#layout.fruchterman.reingold(iG, niter=10000)#layout.fruchterman.reingold(iG);
  
  iG;
}
mp_assignColors<-function # takes a vector of values and returns a vecotr of colors (low=blue, high=red)
(vect){
  require(gplots);
  xcols<-rev(heat.colors(length(vect))) #bluered(length(vect));
  mycols<-rep('', length(vect));
  mycols[order(vect)]<-xcols;
  mycols;
}



ig_NiceGraph<-function# returns a pretty graph given a grnTab and expression data
(object,
 myG,
 p,# result of running induced.subgraph(graphX, v=genes), and ig_convertSmall
 vLabels='Regulator'
 # bold=NULL){
){
  #  myG<-ig_convertSmall(myG);
  # calculate correlations between TF->TG
  expDat <- object@ndata
  lEdges<-length(E(myG));
  myCors<-rep(0, length=lEdges);
  edgeColors<-rep("lightblue", length=lEdges);
  for(i in seq(lEdges)){
    edgeI<-V(myG)$name[get.edge(myG, i)];
    tf<-edgeI[1];
    tg<-edgeI[2];
    #cat(tf, "->", tg,":");
    myCors[i]<-sign(cor(as.numeric(expDat[tf,]), as.numeric(expDat[tg,])));
    if(myCors[i]>0){
      edgeColors[i]<-"red";
    }
  }
  E(myG)$color<-edgeColors;  
  E(myG)$arrow.size=.2;
  
  ## node colors
  genes<-V(myG)$name;
  genes <- intersect(names(object@correlation[[p]]), genes)
  #geneVals<-apply(expDat[genes,], 1, mean);
  geneVals<-object@correlation[[p]][genes]
  V(myG)$color<-cRamp2(geneVals)#mp_assignColors(geneVals);
  
  # node size
  outdegree<-igraph::degree(myG, mode='out');
  V(myG)$size<-ig_scaleV( outdegree, sf=10, 4);
  
    #V(myG)[V(myG)]$label<-'';
    # node labels
    if(length(vLabels)==1){
      if(vLabels=="Target"){
        V(myG)[V(myG)$type=="Regulator"]$label<-'';
      }
      else{
        if(vLabels=="Regulator"){
          V(myG)[V(myG)$type=="Target"]$label<-'';
        }
      }
    }
    else{
      others<-setdiff(V(myG)$name, vLabels);
      xi<-match(others, V(myG)$name);
      V(myG)[xi]$label<-'';
    }

  if(FALSE){
    # label face
    if(!is.null(bold)){
      xi<-match(bold, V(myG)$name);
      cat(xi,"\n");
      V(myG)[xi]$label.font=2;
    }
  }
  
  myG;
}



############### end of GRN code #######

### Function to initialize a CellRouter class with the expression data
setMethod("initialize", 
          signature = "CellRouter", 
          definition = function(.Object, expdata, annotations){
            print("Initializing CellRouter object")
            
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@sampTab <- data.frame(sample_id=colnames(expdata), conditions=annotations)
            rownames(.Object@sampTab) <- .Object@sampTab$sample_id
            validObject(.Object)
            return(.Object)
          }
)

setGeneric("trajectoryOverlap", function(object, direction, trajectory.1, trajectory.2, filename, ...) standardGeneric("trajectoryOverlap"))
setMethod("trajectoryOverlap", 
          signature="CellRouter",
          definition=function(object, direction, trajectory.1, trajectory.2, filename){
            
            library("Vennerable");
            venn <- list();
            venn[[trajectory.1]] <- names(object@top.correlations[[direction]][[trajectory.1]])
            venn[[trajectory.2]] <- names(object@top.correlations[[direction]][[trajectory.2]])
            vennD <- Venn(venn);
            filename <- paste(filename, '_',trajectory.1, '_',trajectory.2, '_', direction, '.pdf',sep='')
            pdf(file=filename, width=3, height=3);
            plot(vennD,doWeights=TRUE)
            dev.off();
          }
)

## Filter dataset, normalize and scale
## Basic filtering on cells and genes.
## Median normalization, log transformation
## Also implement a total count normalization (as in indrops paper)
## Users can directly specify a filtered, normalized dataset if they want to.
setGeneric("filterdata", function(object, min.cells=5, min.genes=2000, min.ExprSum=60, max.expr=5000, scale=TRUE, center=TRUE, ...) standardGeneric("filterdata"))
setMethod("filterdata", 
          signature="CellRouter",
          definition=function(object, min.cells, min.genes, min.expr, max.expr, scale, center, ...){
            print("Filtering data")
            print("Filtering cells and genes...")
            
            #Cell filtering
            num.genes <- colSums(object@expdata > 0);
            cells.use <- names(num.genes[which(num.genes > min.genes)])
            tmp.data <- object@expdata[,cells.use]
            
            #Gene filtering
            num.cells <- rowSums(tmp.data > 0);
            genes.use <- names(num.cells[which(num.cells > min.cells)])
            genes.use1 <- names(which(rowSums(object@expdata[,cells.use] > min.expr) > min.ExprSum))
            genes.use2 <- rownames(object@expdata)[apply(object@expdata[,cells.use], 1, max) < max.expr]
            genes.use <- Reduce(list(genes.use, genes.use1, genes.use2), intersect)
            tmp.data <- tmp.data[genes.use, ]
            
            ## Meadian normalization (Cell, Bipolar dataset)
            col.sums <- apply(tmp.data, 2, sum)
            median.expr <- median(col.sums)
            norm.data <- median.expr * scale(tmp.data, center=TRUE, scale=col.sums)
            rm(tmp.data)
            object@ndata <- as.data.frame(log(norm.data + 1))
            
            return(object)
          }
)
cGrad <- function(x, num_colors){
  cols <- colorRamp(c('blue', 'orange', 'red'))(range01(x))
  #cols <- colorRamp(c('deepskyblue3', 'white', 'red'))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}

## Implement function to reduce dimension first ####

### Discovering supopulation structure ##
setGeneric("findsubpopulations", function(object, k=5, sim.type="jaccard", filename="graph_subpopulations.gml", ...) standardGeneric("findsubpopulations"))

setMethod("findsubpopulations",
          signature = "CellRouter",
          definition = function(object, k, sim.type, filename, ...){

            library('cccd')
            library('proxy') # Library of similarity/dissimilarity measures for 'dist()'
            matrix <- object@rdimension
            sampTab <- object@sampTab
            
            print('building k-nearest neighbors graph')
            dm <- as.matrix(dist(matrix))
            h <- nng(dx=dm,k=k)
            V(h)$name <- rownames(matrix)
            if(sim.type == 'jaccard'){
              sim <- similarity.jaccard(h, vids=V(h), loops=FALSE)
            }else if(sim.type == 'invlog'){
              sim <- similarity.invlogweighted(h, vids=V(h))
            }
            rownames(sim) <- V(h)$name
            colnames(sim) <- V(h)$name
            
            edges <- as.data.frame(get.edgelist(h))
            rownames(edges) <- paste(edges$V1, edges$V2, sep='_')
            edges$weight <- 0
            for(i in 1:nrow(edges)){
              edge <- edges[i,]
              from <- as.character(edge$V1)
              to <- as.character(edge$V2)
              metric <- as.numeric(sim[from, to])
              index <- paste(from, to, sep='_')
              edges[index, 'weight'] <- metric
            }
            E(h)$weight <- edges$weight
            
            ## Community detection to discover subpopulation structure
            print('discoverying subpopulation structure')
            comms <- multilevel.community(as.undirected(h), weights = E(h)$weight)
            V(h)$comms <- membership(comms)
            cell.comms <- commToNames(comms, 'SP') #SP means SubPopulations
            allcells <- as.vector(unlist(cell.comms))
            
            ## Making sure that color mappings are correct
            sampTab <- sampTab[allcells,] #changesorder of cells in the table
            sampTab$population <- ''
            sampTab$colors <- ''
            comm.colors <- cRampClust(unique(membership(comms)), length(unique(membership(comms))))
            names(comm.colors) <- names(cell.comms)
            for(comm in names(cell.comms)){
              sampTab[cell.comms[[comm]], 'population'] <- comm
              sampTab[cell.comms[[comm]], 'colors'] <- comm.colors[comm]
            }
            sampTab$community <- as.vector(unlist(lapply(strsplit(sampTab$population, split="_"), "[", 2)))
            
            ## mapping information to the igraph object
            V(h)[rownames(sampTab)]$subpopulation <- sampTab$colors
            V(h)[rownames(sampTab)]$colors <- sampTab$colors
            V(h)[names(nodeLabels(sampTab,'community'))]$label <- nodeLabels(sampTab, 'community')
            V(h)$size <- 5
            E(h)$arrow.size <- 0.01
            colors <- rainbow(max(membership(comms)))
            
            print("plotting graph in RStudio")
            plot(h,vertex.color=V(h)$colors, vertex.frame.color=V(h)$colors, layout=as.matrix(matrix))
            print('done plotting graph')
            
            ## Useful information about the graph
            graph <- list()
            graph[['network']] <- h
            graph[['edges']] <- edges
            graph[['similarity_matrix']] <- sim
            graph[['subpopulation']] <- cell.comms
            graph[['communities']] <- comms
            
            print('updating CellRouter object')
            object@graph <- graph
            object@expdata <- object@expdata[,rownames(sampTab)]
            object@ndata <- object@ndata[,rownames(sampTab)]
            object@sampTab <- sampTab
            
            write.graph(graph = h, file = filename, format = 'gml')
            
            rm(h)
            rm(edges)
            rm(sim)
            
            return(object)
          }
)

setGeneric("findK", function(object, numK=20) standardGeneric("findK"))            
setMethod("findK",
          signature="CellRouter",
          definition=function(object, numK){
            library(cccd)
            library('proxy') # Library of similarity/dissimilarity measures for 'dist()'
            dm <- as.matrix(dist(object@rdimension))
            for(k in 2:numK){
              h <- nng(dx=dm,k=k)
              if(is.connected(h, mode="weak")){
                cat('minimum value of K that provides a connected graph: ', k, '\n')
                return(k);
              }
            }
          }
)

setGeneric("createKNN", function(object, k, sim.type, filename) standardGeneric("createKNN"))
setMethod("createKNN",
          signature="CellRouter",
          definition=function(object, k, sim.type, filename){
            
            library(cccd)
            library('proxy') # Library of similarity/dissimilarity measures for 'dist()'
            
            matrix <- object@rdimension
            sampTab <- object@sampTab
            
            dm <- as.matrix(dist(matrix))
            h <- nng(dx=dm,k=k)
            V(h)$name <- rownames(matrix)
            if(sim.type == 'jaccard'){
              sim <- similarity.jaccard(h, vids=V(h), loops=FALSE)
            }else if(sim.type == 'invlog'){
              sim <- similarity.invlogweighted(h, vids=V(h))
            }
            rownames(sim) <- V(h)$name
            colnames(sim) <- V(h)$name
            
            edges <- as.data.frame(get.edgelist(h))
            rownames(edges) <- paste(edges$V1, edges$V2, sep='_')
            edges$weight <- 0
            for(i in 1:nrow(edges)){
              edge <- edges[i,]
              from <- as.character(edge$V1)
              to <- as.character(edge$V2)
              metric <- as.numeric(sim[from, to])
              index <- paste(from, to, sep='_')
              edges[index, 'weight'] <- metric
            }
            
            E(h)$weight <- edges$weight
            comms <- cellrouter@graph$communities
            V(h)$comms <- membership(comms)
            cell.comms <- commToNames(comms, 'SP')
            allcells <- as.vector(unlist(cell.comms))
            
            sampTab <- sampTab[allcells,] #change order of cells in the table
            sampTab$population <- ''
            sampTab$colors <- ''
            comm.colors <- cRampClust(unique(comms$membership), length(unique(comms$membership)))
            names(comm.colors) <- names(cell.comms)
            for(comm in names(cell.comms)){
              sampTab[cell.comms[[comm]], 'population'] <- comm
              sampTab[cell.comms[[comm]], 'colors'] <- comm.colors[comm]
            }
            sampTab$community <- as.vector(unlist(lapply(strsplit(sampTab$population, split="_"), "[", 2)))
            V(h)[rownames(sampTab)]$subpopulation <- sampTab$colors
            V(h)[rownames(sampTab)]$colors <- sampTab$colors
            V(h)[names(nodeLabels(sampTab,'community'))]$label <- nodeLabels(sampTab, 'community') #as.vector(sampTab$community)
            #outdegree<-degree(h);
            #V(h)$size<-ig_scaleV(outdegree, sf=10, 4);
            V(h)$size <- 2
            E(h)$arrow.size <- 0.01
            colors <- rainbow(max(membership(comms)))
            plot(h,vertex.color=V(h)$colors, vertex.frame.color=V(h)$colors, layout=as.matrix(matrix))
            
            graph <- list()
            graph[['network']] <- h
            graph[['edges']] <- edges
            graph[['similarity_matrix']] <- sim
            graph[['subpopulation']] <- cell.comms
            graph[['communities']] <- comms
            
            print('updating CellRouter object')
            object@graph <- graph
            object@expdata <- object@expdata[,rownames(sampTab)]
            object@ndata <- object@ndata[,rownames(sampTab)]
            object@sampTab <- sampTab
            
            write.graph(graph = h, file = filename, format = 'gml')
            
            rm(h)
            rm(edges)
            rm(sim)
            
            return(object)
          }
)



### finding gene signatures for each subpopulation (genes more highly expression in each subpopulation)
setGeneric("diffexpr", function(object, foldchange=1, pvalue=0.01) standardGeneric("diffexpr"))

setMethod("diffexpr",
          signature = "CellRouter",
          definition = function(object, foldchange, pvalue){
            print('discovering subpopulation-specific gene signatures')
            expDat <- object@ndata
            membs <- as.vector(object@sampTab$population)
            diffs <- list()
            for(i in unique(membs)){
              if(sum(membs == i) == 0) next
              m <- if(sum(membs != i) > 1) apply(expDat[, membs != i], 1, mean) else expDat[, membs != i]
              n <- if(sum(membs == i) > 1) apply(expDat[, membs == i], 1, mean) else expDat[, membs == i]
              
              pv <- binompval(m/sum(m),sum(n),n)
              d <- data.frame(mean.np=m, mean.p=n, fc=n/m, pv=pv)
              d <- d[!is.infinite(d$fc),]
              d <- d[which(d$pv < 0.05),]
              #d <- d[order(d$pv, decreasing=FALSE),]
              #d <- d[order(d$mean.p, decreasing=TRUE),]
              #diffs[[i]] <- d[which(d$pv < pvalue & d$fc > foldchange),]
              diffs[[i]] <- d
            }
            object@signatures <- diffs
            return (object)
          }
)

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
            #df <- df[order(nchar(df$population)),]
            return(df)
          }
)

######## plot dot plot with marker genes and heatmap: do it later... ###
setGeneric("dotplot", function(object, sampTab, genes.use, thr, logtransform=TRUE,width, height, filename) standardGeneric("dotplot"))
setMethod("dotplot",
          signature = "CellRouter",
          definition = function(object, sampTab, genes.use, thr, logtransform, width, height, filename){
            perc <- data.frame(matrix(0, nrow=length(genes.use), ncol=0))
            exp <- perc
            rownames(perc) <- genes.use
            
            for(i in unique(sampTab$population)){
              cells.population <- rownames(sampTab[which(sampTab$population == i),])
              p <- apply(object@ndata[genes.use, cells.population], 1, function(x){sum(x>thr)/length(x)})
              perc <- cbind(perc, p)
              
              v <- apply(object@ndata[genes.use, cells.population], 1, mean)
              exp <- cbind(exp, v)
            }
            colnames(perc) <- unique(sampTab$population)
            colnames(exp) <- unique(sampTab$population)
            rownames(exp) <- sapply(strsplit(rownames(exp), split='__', fixed=TRUE), function(x){x[1]})
            perc$gene <- rownames(perc)
            perc = melt(perc, id.vars  = 'gene')
            exp$gene <- rownames(exp)
            exp <- melt(exp, id.vars='gene')
            exp$size <- perc$value*100
            if(logtransform){
              exp$value <- log(exp$value + 1)
            }
            
            pdf(file=filename, width=width, height=height)
            g <- ggplot(exp, aes(gene, variable)) + geom_point(aes(colour=value, size=size)) +
              theme_bw() + xlab("") + ylab("") +
              theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.text.x=element_blank(),axis.ticks=element_blank(),
                    panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
              scale_color_gradient("average",low ="blue", high = "red")
            print(g)
            dev.off()
          }
)

setGeneric("markersheatmap", function(object, markers, width, height, filename) standardGeneric("markersheatmap"))
setMethod("markersheatmap",
          signature = "CellRouter",
          definition = function(object, markers, width, height, filename){
            
            sampTab <- object@sampTab
            df <- center_with_threshold(object@ndata, 2)
            paletteLength <- 50
            myColor <- colorRampPalette(c("blue","white","red"))(paletteLength)
            myBreaks <- c(seq(min(df), 0, length.out=ceiling(paletteLength/2) + 1), 
                          seq(max(df)/paletteLength, max(df), length.out=floor(paletteLength/2)))
            
            ann_col <- data.frame(population=as.vector(sampTab[,'population']))
            rownames(ann_col) <- rownames(sampTab)
            
            ann_row <- data.frame(signature=as.vector(markers$population))
            rownames(ann_row) <- rownames(markers)
            
            index <- getIndexes(ann_col, ann_row)
            tmp <- sampTab[match(as.vector(ann_row$signature), sampTab$population),'colors']
            
            pdf(file='results/heatmap_test.pdf', width=10, height=7)
            #png(filename, width = width, height = height, units = 'in', res = 1000)
            hmcols<- colorRampPalette(c("darkblue", "black", "yellow"))
            heatmap.2(as.matrix(df[rownames(markers),]), col=hmcols,trace="none",
                      density.info="none", scale="none",margin=c(10,10), key=TRUE, Colv=F, Rowv=F,
                      srtCol=60, dendrogram="none", cexCol=0.75, cexRow=0.65, labRow=FALSE, 
                      labCol = FALSE, symm=T,symkey=T, ColSideColors=sampTab$colors,
                      RowSideColors=tmp, colsep=index$colsep, rowsep=index$rowsep, sepcolor = 'black')
            dev.off()
          }
)

#flow network from all cells in one population against all others cells in all
#other populations...
setGeneric("findpathsallcells", function(object, maindir) standardGeneric("findpathsallcells"))
setMethod("findpathsallcells",
          signature="CellRouter",
          definition=function(object, maindir){
            curdir <- getwd()
            dirs <- list()
            #for(i in 1:length(sources)){
            #  for(j in 1:length(targets)){
            sampTab <- object@sampTab
            for(source in sources){
              for(target in targets){
                #ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
                s <- rownames(sampTab[which(sampTab$population == source),])
                t <- rownames(sampTab[which(sampTab$population == target),])
                dir <- paste(source, target, sep='.')
                cat('--------------------------:', dir, '\n')
                dir.create(file.path(maindir, dir), showWarnings = FALSE)
                system(paste('chmod 777 ', file.path(maindir, dir), sep=''))
                setwd(file.path(maindir, dir))
                dirs[[dir]] <- file.path(maindir, dir)
                write.table(s, file='cell_sources.txt', sep=',', row.names=FALSE, quote=FALSE)
                write.table(t, file='cell_sinks.txt', sep=',', row.names=FALSE, quote=FALSE)
                
                sink("CellRouter.sh")
                cat('#!/bin/sh\n')
                cat('# run CellRouter, run!\n')
                cat('LIB_DIR=~/NetBeansProjects/CellRouter/\n')
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

setGeneric("findpaths", function(object, libdir, maindir, method) standardGeneric("findpaths"))
setMethod("findpaths",
          signature="CellRouter",
          definition=function(object, libdir, maindir, method){
            curdir <- getwd()
            dirs <- list()
            
            if(method %in% c("euclidean", "maximum", "manhattan","canberra","binary")){
              bla2 <- as.data.frame(as.matrix(dist(object@rdimension, method = method)))
            }else{
              g <- object@graph$network
              ##g <- simplify(g, remove.loops = TRUE)
              ##E(g)$weight <- -1*E(g)$weight
              ##bla2 <- as.data.frame(shortest.paths(g, v=V(g), to=V(g)))
              bla2 <- as.data.frame(igraph::distances(g, v=V(g), to=V(g), weights = NULL, algorithm = "bellman-ford"))
            }
            sampTab <- object@sampTab
            b <- list()
            for(p1 in sources){
              cellsp1 <- as.vector(sampTab[which(sampTab$population == p1), 'sample_id'])
              for(p2 in targets){
                if(p1 != p2){
                  cellsp2 <- as.vector(sampTab[which(sampTab$population == p2), 'sample_id'])
                  x <- bla2[cellsp1, cellsp2]
                  x$population1 <- p1
                  x$population2 <- p2
                  x$merge <- rownames(x)
                  
                  bla3 <- melt(x, id.vars=c('population1', 'population2', 'merge'))
                  bla3 <- bla3[order(bla3$value, decreasing = TRUE),]
                  b[[p1]][[p2]] <- bla3[1,]
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
                cat('--------------------------:', dir, '\n')
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



setGeneric("findpaths2", function(object, maindir) standardGeneric("findpaths2"))
setMethod("findpaths2",
          signature="CellRouter",
          definition=function(object, maindir){
            curdir <- getwd()
            dirs <- list()
            for(i in 1:length(sources)){
              for(j in 1:length(targets)){
                #ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
                s <- sources[i]
                t <- targets[j]
                dir <- paste(names(s), names(t), sep='.')
                cat('--------------------------:', dir, '\n')
                dir.create(file.path(maindir, dir), showWarnings = FALSE)
                system(paste('chmod 777 ', file.path(maindir, dir), sep=''))
                setwd(file.path(maindir, dir))
                dirs[[dir]] <- file.path(maindir, dir)
                write.table(sources[i], file='cell_sources.txt', sep=',', row.names=FALSE, quote=FALSE)
                write.table(targets[j], file='cell_sinks.txt', sep=',', row.names=FALSE, quote=FALSE)
                
                sink("CellRouter.sh")
                cat('#!/bin/sh\n')
                cat('# run CellRouter, run!\n')
                cat('LIB_DIR=~/NetBeansProjects/CellRouter/\n')
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


## remove genes with zero-variance, 
## find informative genes by PCA analysis
## parse paths
##make it optional: low cost paths, high flow paths, all paths
##for the paper, I think all paths should be included, to increase relevance!
setGeneric("processtrajectories2", function(object, genes, num.cells) standardGeneric("processtrajectories2"))
setMethod("processtrajectories2",
          signature="CellRouter",
          definition=function(object, genes, num.cells){
            library(igraph)
            
            #print('finding informative genes using principal component analysis')
            #pca <- prcomp(t(cellrouter@ndata), scale=TRUE, center=TRUE)
            #loadings <- pca$rotation
            #genes <- unique(as.vector(unlist(apply(loadings[,1:num.pc], 2, 
            #                                          function(x){names(x[which(abs(x) >= quantile(x, quantile))])}))))
            
            object@genes.trajectory <- genes
            
            print('parsing trajectory information')
            #opening the top path
            #paths <- do.call(rbind, lapply(object@directory, function(x){read.csv(paste(x, 'Cells_FlowNetwork_all_paths.txt', sep='/'),
            #                                                                sep="\t", stringsAsFactors=FALSE)[1,]}))
            
            paths <- lapply(object@directory, function(x){read.csv(paste(x, 'Cells_FlowNetwork_all_paths.txt', sep='/'),
                                                                   sep="\t", stringsAsFactors=FALSE)})
            paths <- lapply(paths, function(x){x[order(x$path_flow, decreasing = TRUE),]})
            paths <- lapply(paths, function(x){x[duplicated(x$path),]})
            for(p in names(paths)){
              paths[[p]]$population <- rep(p, nrow(paths[[p]]))
            }
            paths <- do.call(rbind, lapply(paths, function(x){x}))
            #paths <- do.call(rbind, lapply(paths, function(x){x[1,]}))
            
            #remove empty paths
            paths <- paths[complete.cases(paths),]
            #opening flow networks
            networks <- lapply(cellrouter@directory, 
                               function(x){
                                 file <- list.files(x, '*.gml*');
                                 if(length(file) > 0){
                                   cat(file, '\n');
                                   read.graph(paste(x, 'Cells_FlowNetwork_all_paths_subnet.gml', sep='/'), format = 'gml')
                                 }
                               })
            networks <- networks[lapply(networks, length) > 0]
            
            paths$ID <- paste('path', 1:nrow(paths), sep='__')
            #paths$population <- rownames(paths)
            object@paths <- paths
            object@networks <- networks
            
            object@pathsinfo <- pathsinfo2(object, num.cells)
            return (object)
          }          
          
)


## remove genes with zero-variance, 
## find informative genes by PCA analysis
## parse paths
##make it optional: low cost paths, high flow paths, all paths
##for the paper, I think all paths should be included, to increase relevance!
setGeneric("processtrajectories", function(object, genes, path.rank, num.cells, neighs) standardGeneric("processtrajectories"))
setMethod("processtrajectories",
          signature="CellRouter",
          definition=function(object, genes, path.rank, num.cells, neighs){
            library(igraph)
            
            #print('finding informative genes using principal component analysis')
            #pca <- prcomp(t(cellrouter@ndata), scale=TRUE, center=TRUE)
            #loadings <- pca$rotation
            #genes <- unique(as.vector(unlist(apply(loadings[,1:num.pc], 2, 
            #                                          function(x){names(x[which(abs(x) >= quantile(x, quantile))])}))))
            
            object@genes.trajectory <- genes
            
            print('parsing trajectory information')
            #opening the top path
            #paths <- do.call(rbind, lapply(object@directory, function(x){read.csv(paste(x, 'Cells_FlowNetwork_all_paths.txt', sep='/'),
            #                                                                sep="\t", stringsAsFactors=FALSE)[1,]}))
            
            paths <- lapply(object@directory, function(x){read.csv(paste(x, 'Cells_FlowNetwork_all_paths.txt', sep='/'),
                                                                            sep="\t", stringsAsFactors=FALSE)})
            #paths <- lapply(paths, function(x){x[order(x$path_flow, decreasing = TRUE),]})
            if(path.rank == "path_cost"){
              paths <- lapply(paths, function(x){x[order(x[,path.rank], decreasing = FALSE),]})
            }else{
              paths <- lapply(paths, function(x){x[order(x[,path.rank], decreasing = TRUE),]})
            }
            #paths <- lapply(paths, function(x){x[order(x$length, decreasing = TRUE),]})
            #paths <- lapply(paths, function(x){x[order(x$path_cost, decreasing = FALSE),]})
            paths <- lapply(paths, function(x){x[duplicated(x$path),]})
            
            #paths <- do.call(rbind, lapply(paths, function(x){x}))
            paths <- do.call(rbind, lapply(paths, function(x){x[1,]}))
            
            #remove empty paths
            paths <- paths[complete.cases(paths),]
            #opening flow networks
            networks <- lapply(cellrouter@directory, 
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
            
            object@pathsinfo <- pathsinfo(object, num.cells = num.cells, neighs = neighs)
            return (object)
          }          
          
)

##same genes stored in object@genes.trajectory
# setGeneric("processTrajectoriesSubpopulation", function(object, subpopulation,...) standardGeneric("processTrajectoriesSubpopulation"))
# setMethod("processTrajectoriesSubpopulation",
#           signature="CellRouter",
#           definition=function(object, num.pc, quantile){
#             library(igraph)
#             
#             genes <- object@genes.trajectory
#             
#             print('parsing trajectory information')
#             #opening all path
#             paths <- do.call(rbind, lapply(object@directory[subpopulation], function(x){read.csv(paste(x, 'Cells_FlowNetwork_all_paths.txt', sep='/'),
#                                                                                   sep="\t", stringsAsFactors=FALSE)}))
# 
#             paths$ID <- paste('path', 1:nrow(paths), sep='_')
#             paths$population <- rownames(paths)
#             object@subpopulation$paths <- paths
#             
#             object@subpopulation$pathsinfo <- pathsinfo(object)
#             return (object)
#           }          
#           
# )
# 
# ##write a subpopulation-specific pathsinfo...

setGeneric("pathsinfo2", function(object, num.cells) standardGeneric("pathsinfo2"))
setMethod("pathsinfo2",
          signature="CellRouter",
          definition=function(object, num.cells){
            
            paths <- object@paths
            expDat <- object@ndata[object@genes.trajectory,]
            sampTab <- object@sampTab
            
            
            path_info <- list()
            pathsDF <- data.frame()
            
            pathslist <- split(paths, paths$population)
            for(p in names(pathslist)){
              i <- 1
              xpath <- pathslist[[p]]
              for(path in as.vector(xpath$path)){
                split_paths <- strsplit(as.vector(path), "->")[[1]]
                cells <- split_paths[c(-1,-length(split_paths))]
                if(length(cells) > num.cells){ #only include paths with more than num.cells cells
                  expression <- expDat[, cells]
                  ##remove genes with zero variance along a path
                  var <- apply(expression, 1, var)
                  expression <- expression[which(var != 0),]
                  
                  #path_name <- paste("path_", i, sep="")
                  path_name <- paste(xpath$population[i], "__", i, sep='')
                  print(path_name)
                  
                  path_info[[p]][['distr']][[path_name]] <- expression
                  path_info[[p]][['path']][[path_name]] <- cells
                  path_info[[p]][['pathInfo']][[path_name]] <- data.frame(path=paste(as.vector(cells), collapse=','), 
                                                                     source=cells[1], sink=cells[length(cells)], cost=paths$path_cost[i],
                                                                     source_population=sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),'community'],
                                                                     target_population=sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), 'community'])
                  pathsDF[path_name, 'source'] <- cells[1]
                  pathsDF[path_name, 'sink'] <- cells[length(cells)]
                  pathsDF[path_name, 'cost'] <- paths$path_cost[i]
                  pathsDF[path_name, 'source_population'] <- as.character(sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),'community'])
                  pathsDF[path_name, 'source_color'] <- as.character(sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),'colors'])
                  pathsDF[path_name, 'target_population'] <- as.character(sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), 'community'])
                  pathsDF[path_name, 'target_color'] <- as.character(sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), 'colors'])
                  
                  pseudotime <- (0:(length(cells)-1)) / length(cells) #stick to this
                  names(pseudotime) <- cells
                  path_info[[p]][['pseudotime']][[path_name]] <- pseudotime
                  path_info[[p]]$path_data <- pathsDF
                  i <- i+1
                }
              }
            }
            
            return (path_info)
          }
)

setGeneric('plotgraph', function(object, ggrn, scores, transition, vLabels, filename, ...) standardGeneric('plotgraph'))
setMethod('plotgraph',
          signature="CellRouter",
          definition=function(object, ggrn, scores, vLabels, filname){
            
            rgrn <- induced.subgraph(ggrn, vids=unlist(neighborhood(graph=ggrn,order=1,nodes=names(scores))))
            subnet <- ig_convertSmall(rgrn, exponent = 0.5)
            subnet <- ig_NiceGraph(object, subnet, transition, vLabels = vLavels)
            write.graph(subnet, file='results/SP_4.SP_5_network.gml', format='gml')
          }
)

setGeneric('plotclusters', function(object, p, columns, width, height, filename,...) standardGeneric('plotclusters'))
setMethod('plotclusters',
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
                xlab("") + ylab("") +
                #scale_fill_gradient2("zscore", low="navy", high="red") + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
                theme(legend.position="right", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
                      axis.title.y = element_text(size = rel(0.3), angle = 90), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
                      panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
                ggtitle(cl)
              
              plots[[cl]] <- g1
            }
            
            #pdf(file=filename, width=width, height=height)
            png(file=filename, width=width, height=height)
            multiplot(plotlist = plots, cols=columns)
            dev.off()
          }
)


setGeneric('plotbranch', function(object, direction, p1, p2, columns, width, height, filename,...) standardGeneric('plotbranch'))
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
            pdf(file=filename, width=width, height=height)
            multiplot(plotlist = plots, cols=columns)
            dev.off()
          }
)

setGeneric('plottr', function(object, p, scores, cluster=TRUE, columns, width, height, filename,...) standardGeneric('plottr'))
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
             pheatmap(matrix, color = myColor, breaks = myBreaks, 
                      cluster_cols = FALSE, clustering_method = 'ward.D', main=p,
                      show_colnames = FALSE)
            
            matrix$cluster <- rownames(matrix)
            matrix.m <- melt(matrix, id.var="cluster")
            matrix.m$cluster <- factor(rownames(matrix), levels=rev(rownames(matrix)))
            g2 <- ggplot(matrix.m, aes(variable, cluster)) + geom_tile(aes(fill = value)) +
              #scale_fill_gradientn("Derivative",colours=colors) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              scale_fill_gradient2("Derivative", low="navy", high="red") + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              #scale_fill_gradient("", low="navy", high) + theme_bw() + xlab("CellRouter trajectory") + ylab("") +
              theme(legend.position="right", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
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
              theme(legend.position="right", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
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
            
            matrix2
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

##grn score on branch-specific genes
setGeneric('grntransition', function(object, tfs, transitions, dir.targets='up',
                                 q.up=0.95, q.down=0.05, 
                                 columns=1, width, height, filename,...) standardGeneric('grntransition'))
setMethod('grntransition',
          signature="CellRouter",
          definition=function(object, tfs, transitions,
                              dir.targets, q.up, q.down,
                              columns, width, height, filename){
            plots <- list()
            ##obtain gene sets unique to each trajectory...
            transition.genes<-list();
            for(t in transitions){
              others<-setdiff(transitions, t);
              tmp <- object@top.correlations
              x.up<-setdiff(names(tmp[['up']][[t]]), unlist(lapply(tmp[['up']][others], names)));
              x.down<-setdiff(names(tmp[['down']][[t]]), unlist(lapply(tmp[['down']][others], names)));
              transition.genes[['up']][[t]]<-object@correlation[[t]][x.up];
              transition.genes[['down']][[t]]<-object@correlation[[t]][x.down];
            }
            
            #ps <- list()
            allscores <- list()
            for(p in transitions){
              tfs.transition.up <- intersect(tfs,names(transition.genes[['up']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
              tfs.transition.down <- intersect(tfs,names(transition.genes[['down']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
              tfs.transition <- c(tfs.transition.up, tfs.transition.down)
              tfs.transition <- intersect(V(ggrn)$name, tfs.transition)
              tfs.transition <- object@correlation[[p]][tfs.transition]
              averages <- vector()
              names <- vector()
              tf.targets <- list()
              proportions <- vector()
              for(r in names(tfs.transition)){
                #all genes connected to r
                rgrn <- induced.subgraph(ggrn, vids=unlist(neighborhood(graph=ggrn,order=1,nodes=r)))
                #x <- object@top.correlations[[dir.targets]][[p]]
                x <- transition.genes[[dir.targets]][[p]]
                #all genes connected to r and regulated along trajectory p
                genes <- intersect(V(rgrn)$name, names(x))
                bla <- cellrouter@correlation[[p]][genes]
                if(length(bla) == 0){
                  cat(r, 'has no targets\n')
                }else if(length(bla) == 1 & names(bla) == r){
                  cat(r, 'only regulates only itself\n')
                }else if(length(bla) > 0){
                  tf.targets[[r]] <- names(bla)
                  averages <- append(averages, mean(bla))
                  proportion <- length(genes) / length(V(rgrn)$name)
                  proportions <- append(proportions, proportion)
                  names <- append(names, r)
                }
                #df <- data.frame(gene=as.vector(names(bla)), corr=as.numeric(bla), TF=rep(r, length(bla)), stringsAsFactors = FALSE)
                #ps[[r]] <- df
              }
              names(averages) <- names
              names(proportions) <- names
              aux <- averages[is.na(averages)]
              print(length(aux))
              averages <- averages[!is.na(averages)]
              averages <- averages[order(averages, decreasing=TRUE)]
              proportions <- proportions[names(averages)]
              scores <- tfs.transition[names(averages)] * averages * proportions
              #scores <- scores[order(scores, decreasing=TRUE)]
              #scores <- scores[1:num.genes]
              scores <- scores[which(scores > quantile(scores, q.up) | scores < quantile(scores, q.down))]
              scores <- scores[order(scores, decreasing=TRUE)]
              allscores[[p]] <- list(scores=scores, targets=tf.targets)
              df <- data.frame(gene=names(scores), score=as.numeric(scores))
              df <- df[order(df$score, decreasing = TRUE),]
              df$gene <- factor(df$gene, levels=rev(df$gene))
              colors <- c('blue', 'white', 'red')
              #pdf(file=filename, width=width, height=height)
              g <- ggplot(df, aes(x=gene, y=score, fill=score)) + #geom_violin(scale="width") + stat_summary(fun.y=mean,geom='point') + #geom_boxplot(alpha=.9) + 
                geom_bar(stat = 'identity', color='black') + #geom_boxplot(alpha=.9) +
                scale_fill_gradientn("",colours=colors) +
                #ggplot(means, aes(x=p, y=mean, fill=p)) + geom_bar(stat='identity') + 
                theme_bw() + xlab("") + ylab("GRN score") + theme(legend.position="none") +
                theme(axis.text.x = element_text(size=rel(1), angle=00, hjust=1)) + 
                ggtitle(paste(p)) +
                theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      plot.background=element_blank(),
                      panel.border = element_blank(),
                      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) +
                coord_flip()
              plots[[p]] <- g
              #panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1))
              #print(g)
              #dev.off()
            }
            file <- paste(filename, '_target_direction_', dir.targets, '.pdf', sep='')
            pdf(file = file, width=width, height=height)
            multiplot(plotlist = plots,cols=columns)
            dev.off()
            allscores
            #list(scores=scores, targets=tf.targets)
          }
)

setGeneric('grnscores', function(object, tfs, transitions, direction, dir.targets='up',
                                 q.up=0.95, q.down=0.05,
                                 columns=1, width, height, flip, filename,...) standardGeneric('grnscores'))
setMethod('grnscores',
          signature="CellRouter",
          definition=function(object, tfs, transitions, direction,
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
                rgrn <- induced.subgraph(ggrn, vids=unlist(neighborhood(graph=ggrn,order=1,nodes=r)))
                x <- object@top.correlations[[dir.targets]][[p]]
                genes <- intersect(V(rgrn)$name, names(x)) #subnetwork active during transition p
                bla <- cellrouter@correlation[[p]][genes]
                if(length(bla) == 0){
                  cat(r, 'has no targets\n')
                }else if(length(bla) == 1 & names(bla) == r){
                  cat(r, 'regulates only itself\n')
                }else if(length(bla) > 0){ #at least one target required
                  tf.targets[[r]] <- names(bla)
                  averages <- append(averages, mean(bla))
                  num.genes <- append(num.genes, length(bla))
                  names <- append(names, r)
                }
                #df <- data.frame(gene=as.vector(names(bla)), corr=as.numeric(bla), TF=rep(r, length(bla)), stringsAsFactors = FALSE)
                #ps[[r]] <- df
              }
              names(averages) <- names
              names(num.genes) <- names
              aux <- averages[is.na(averages)]
              print(length(aux))
              averages <- averages[!is.na(averages)]
              averages <- averages[order(averages, decreasing=TRUE)]
              num.genes <- num.genes[names(averages)]
              #rescale num.genes
              num.genes <- rescale(num.genes, newrange = c(0.01,1)) #when it rescaled, it changes the scores when only up or down-regulated genes are included...
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
                #x <- scores[order(scores, decreasing=TRUE)]
                #x <- x[1:q.up]
                #xx <- scores[order(scores, decreasing=FALSE)]
                #xx <- xx[1:q.down]
                #scores <- c(x, xx)
                colors <- c('blue', 'white', 'red')
              }
              scores <- scores[order(scores, decreasing=TRUE)]
              allscores[[p]] <- list(scores=scores, targets=tf.targets)
              df <- data.frame(gene=names(scores), score=as.numeric(scores))
              df <- df[order(df$score, decreasing = TRUE),]
              angle=00
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
            allscores
            #list(scores=scores, targets=tf.targets)
          }
)



# setGeneric('grnscores', function(object, tfs, transitions, direction, dir.targets='up',
#                                  q.up=0.95, q.down=0.05,
#                                  columns=1, width, height, flip, filename,...) standardGeneric('grnscores'))
# setMethod('grnscores',
#           signature="CellRouter",
#           definition=function(object, tfs, transitions, direction,
#                               dir.targets, q.up, q.down,
#                               columns, width, height, flip, filename){
#             plots <- list()
#             #ps <- list()
#             allscores <- list()
#             for(p in transitions){
#               if(direction=='up'){
#                 tfs.transition <- intersect(tfs,names(object@top.correlations[['up']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
#               }else if(direction == 'down'){
#                 tfs.transition <- intersect(tfs,names(object@top.correlations[['down']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
#               }else{
#                 tfs.transition.up <- intersect(tfs,names(object@top.correlations[['up']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
#                 tfs.transition.down <- intersect(tfs,names(object@top.correlations[['down']][[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
#                 tfs.transition <- c(tfs.transition.up, tfs.transition.down)
#               }
#               
#               tfs.transition <- intersect(V(ggrn)$name, tfs.transition)
#               tfs.transition <- object@correlation[[p]][tfs.transition]
#               averages <- vector()
#               num.genes <- vector()
#               names <- vector()
#               tf.targets <- list()
#               for(r in names(tfs.transition)){
#                 rgrn <- induced.subgraph(ggrn, vids=unlist(neighborhood(graph=ggrn,order=1,nodes=r)))
#                 x <- object@top.correlations[[dir.targets]][[p]]
#                 genes <- intersect(V(rgrn)$name, names(x)) #subnetwork active during transition p
#                 bla <- cellrouter@correlation[[p]][genes]
#                 if(length(bla) == 0){
#                   cat(r, 'has no targets\n')
#                 }else if(length(bla) == 1 & names(bla) == r){
#                   cat(r, 'regulates only itself\n')
#                 }else if(length(bla) > 0){ #at least one target required
#                   tf.targets[[r]] <- names(bla)
#                   averages <- append(averages, mean(bla))
#                   num.genes <- append(num.genes, length(bla))
#                   names <- append(names, r)
#                 }
#                 #df <- data.frame(gene=as.vector(names(bla)), corr=as.numeric(bla), TF=rep(r, length(bla)), stringsAsFactors = FALSE)
#                 #ps[[r]] <- df
#               }
#               names(averages) <- names
#               names(num.genes) <- names
#               aux <- averages[is.na(averages)]
#               print(length(aux))
#               averages <- averages[!is.na(averages)]
#               averages <- averages[order(averages, decreasing=TRUE)]
#               num.genes <- num.genes[names(averages)]
#               #rescale num.genes
#               num.genes <- rescale(num.genes, newrange = c(0.01,1))
#               scores <- tfs.transition[names(averages)] * averages * num.genes #* number of genes regulated
#               #scores <- scores[order(scores, decreasing=TRUE)]
#               #scores <- scores[1:num.genes]
#               ## if up or down, q.up or q.down are the top genes
#               ##if both, q.up or q.down are quantiles
#               if(direction=='up'){
#                 scores <- scores[order(scores, decreasing=TRUE)]
#                 scores <- scores[1:q.up]
#                 colors <- c('white', 'orange', 'red')
#               }else if(direction=='down'){
#                 scores <- scores[order(scores, decreasing=FALSE)]
#                 scores <- scores[1:q.down]
#                 colors <- c('white','lightblue','blue')
#               }else{
#                 #scores <- scores[which(scores > quantile(scores, q.up) | scores < quantile(scores, q.down))]
#                 x <- scores[order(scores, decreasing=TRUE)]
#                 x <- x[1:q.up]
#                 xx <- scores[order(scores, decreasing=FALSE)]
#                 xx <- xx[1:q.down]
#                 scores <- c(x, xx)
#                 colors <- c('blue', 'white', 'red')
#               }
#               #scores <- scores[order(scores, decreasing=TRUE)]
#               allscores[[p]] <- list(scores=scores, targets=tf.targets)
#               df <- data.frame(gene=names(scores), score=as.numeric(scores))
#               df <- df[order(df$score, decreasing = TRUE),]
#               angle=00
#               if(flip){
#                 df$gene <- factor(df$gene, levels=rev(df$gene))
#               }else{
#                 df$gene <- factor(df$gene, levels=df$gene)
#                 angle=45
#               }
#               
#               #pdf(file=filename, width=width, height=height)
#               g <- ggplot(df, aes(x=gene, y=score, fill=score)) + #geom_violin(scale="width") + stat_summary(fun.y=mean,geom='point') + #geom_boxplot(alpha=.9) +
#                 geom_bar(stat = 'identity', color='black') + #geom_boxplot(alpha=.9) +
#                 scale_fill_gradientn("",colours=colors) +
#                 #ggplot(means, aes(x=p, y=mean, fill=p)) + geom_bar(stat='identity') +
#                 theme_bw() + xlab("") + ylab("GRN score") + theme(legend.position="none") +
#                 theme(axis.text.x = element_text(size=rel(1), angle=angle, hjust=1)) +
#                 ggtitle(paste(p)) +
#                 theme(panel.background=element_blank(),
#                       panel.grid.major=element_blank(),
#                       panel.grid.minor=element_blank(),
#                       plot.background=element_blank(),
#                       panel.border = element_blank(),
#                       axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
#                       axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) #+
#                 #coord_flip()
#               if(flip){
#                 g <- g + coord_flip()
#               }
#               plots[[p]] <- g
#               #panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1))
#               #print(g)
#               #dev.off()
#             }
#             file <- paste(filename, '_target_direction_', dir.targets, '.pdf', sep='')
#             pdf(file = file, width=width, height=height)
#             multiplot(plotlist = plots,cols=columns)
#             dev.off()
#             allscores
#             #list(scores=scores, targets=tf.targets)
#           }
# )


#regulator up, targets up -> inducing expression during differentiation
#regulator down, targets up -> de-repression during differentiation
#regulator up, targets down -> repression during differentiation
#regulator down, targets down -> not inducing expression during differentiation
# setGeneric('grnscores', function(object, tfs, transitions, dir.targets='up',
#                                  q.up=0.95, q.down=0.05, 
#                                  columns=1, width, height, filename,...) standardGeneric('grnscores'))
# setMethod('grnscores',
#           signature="CellRouter",
#           definition=function(object, tfs, transitions, 
#                               dir.targets, q.up, q.down,
#                               columns, width, height, filename){
#             plots <- list()
#             #ps <- list()
#             allscores <- list()
#             for(p in transitions){
#               tfs.transition <- intersect(tfs,names(object@correlation[[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
#               tfs.transition <- intersect(V(ggrn)$name, tfs.transition)
#               tfs.transition <- object@correlation[[p]][tfs.transition]
#               averages <- vector()
#               names <- vector()
#               tf.targets <- list()
#               for(r in names(tfs.transition)){
#                 rgrn <- induced.subgraph(ggrn, vids=unlist(neighborhood(graph=ggrn,order=1,nodes=r)))
#                 x <- object@correlation[[p]]
#                 if(dir.targets == 'up'){
#                   x <- x[x > 0]
#                 }else{
#                   x <- x[x < 0]
#                 }
#                 genes <- intersect(V(rgrn)$name, names(x))
#                 bla <- cellrouter@correlation[[p]][genes]
#                 if(length(bla) == 0){
#                   cat(r, 'has no targets')
#                 }else if(length(bla) == 1 & names(bla) == r){
#                   cat(r, 'has only regulates only itself\n')
#                 }else if(length(bla) > 0){
#                   tf.targets[[r]] <- names(bla)
#                   averages <- append(averages, mean(bla))
#                   names <- append(names, r)
#                   
#                 }
#                 #df <- data.frame(gene=as.vector(names(bla)), corr=as.numeric(bla), TF=rep(r, length(bla)), stringsAsFactors = FALSE)
#                 #ps[[r]] <- df
#               }
#               names(averages) <- names
#               aux <- averages[is.na(averages)]
#               print(length(aux))
#               averages <- averages[!is.na(averages)]
#               averages <- averages[order(averages, decreasing=TRUE)]
#               scores <- tfs.transition[names(averages)] * averages
#               #scores <- scores[order(scores, decreasing=TRUE)]
#               #scores <- scores[1:num.genes]
#               scores <- scores[which(scores > quantile(scores, q.up) | scores < quantile(scores, q.down))]
#               scores <- scores[order(scores, decreasing=TRUE)]
#               allscores[[p]] <- list(scores=scores, targets=tf.targets)
#               df <- data.frame(gene=names(scores), score=as.numeric(scores))
#               df <- df[order(df$score, decreasing = TRUE),]
#               df$gene <- factor(df$gene, levels=rev(df$gene))
#               colors <- c('blue', 'white', 'red')
#               #pdf(file=filename, width=width, height=height)
#               g <- ggplot(df, aes(x=gene, y=score, fill=score)) + #geom_violin(scale="width") + stat_summary(fun.y=mean,geom='point') + #geom_boxplot(alpha=.9) + 
#                 geom_bar(stat = 'identity', color='black') + #geom_boxplot(alpha=.9) +
#                 scale_fill_gradientn("",colours=colors) +
#                 #ggplot(means, aes(x=p, y=mean, fill=p)) + geom_bar(stat='identity') + 
#                 theme_bw() + xlab("") + ylab("GRN score") + theme(legend.position="none") +
#                 theme(axis.text.x = element_text(size=rel(1), angle=00, hjust=1)) + 
#                 ggtitle(paste(p)) +
#                 theme(panel.background=element_blank(),
#                       panel.grid.major=element_blank(),
#                       panel.grid.minor=element_blank(),
#                       plot.background=element_blank(),
#                       panel.border = element_blank(),
#                       axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
#                       axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) +
#                 coord_flip()
#               plots[[p]] <- g
#               #panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1))
#               #print(g)
#               #dev.off()
#             }
#             file <- paste(filename, '_target_direction_', dir.targets, '.pdf', sep='')
#             pdf(file = file, width=width, height=height)
#             multiplot(plotlist = plots,cols=columns)
#             dev.off()
#             allscores
#             #list(scores=scores, targets=tf.targets)
#           }
# )



setGeneric('grnscore', function(object, tfs, p, num.genes, width, height, filename) standardGeneric('grnscore'))
setMethod('grnscore',
          signature="CellRouter",
          definition=function(object, tfs, p, num.genes, width, height, filename){
            #TFs and targets positively correlated with pseudotime...
            #tfs.transition <- intersect(tfs,names(cellrouter@correlation[[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
            #tfs.transition <- cellrouter@correlation[[p]][tfs.transition]
            tfs.transition <- intersect(tfs,names(object@top.correlations$up[[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
            tfs.transition <- intersect(V(ggrn)$name, tfs.transition)
            tfs.transition <- object@top.correlations$up[[p]][tfs.transition]
            #tfs.transition <- tfs.transition[which(tfs.transition > 0)]
            ps <- list()
            averages <- vector()
            names <- vector()
            tf.targets <- list()
            for(r in names(tfs.transition)){
              rgrn <- induced.subgraph(ggrn, vids=unlist(neighborhood(graph=ggrn,order=1,nodes=r)))
              genes <- intersect(V(rgrn)$name, names(cellrouter@top.correlations$up[[p]]))
              bla <- cellrouter@top.correlations$up[[p]][genes]
              #bla <- bla[bla > 0]
              tf.targets[[r]] <- names(bla)
              averages <- append(averages, mean(bla))
              names <- append(names, r)
              #df <- data.frame(gene=as.vector(names(bla)), corr=as.numeric(bla), TF=rep(r, length(bla)), stringsAsFactors = FALSE)
              #ps[[r]] <- df
            }
            names(averages) <- names
            averages <- averages[order(averages, decreasing=TRUE)]
            scores <- tfs.transition[names(averages)] * averages
            scores <- scores[order(scores, decreasing=TRUE)]
            scores <- scores[1:num.genes]
            df <- data.frame(gene=names(scores), score=as.numeric(scores))
            df <- df[order(df$score, decreasing = TRUE),]
            df$gene <- factor(df$gene, levels=df$gene)
            colors <- c('white', 'orange', 'red')
            pdf(file=filename, width=width, height=height)
            g <- ggplot(df, aes(x=gene, y=score, fill=score)) + #geom_violin(scale="width") + stat_summary(fun.y=mean,geom='point') + #geom_boxplot(alpha=.9) + 
              geom_bar(stat = 'identity', color='black') + #geom_boxplot(alpha=.9) +
              scale_fill_gradientn("",colours=colors) +
              #ggplot(means, aes(x=p, y=mean, fill=p)) + geom_bar(stat='identity') + 
              theme_bw() + xlab("") + ylab("GRN score") + theme(legend.position="none") +
              theme(axis.text.x = element_text(size=rel(1), angle=45, hjust=1)) + ggtitle(p) +
              theme(panel.background=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.background=element_blank(),
                    panel.border = element_blank(),
                    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) 
            #panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1))
            print(g)
            dev.off()
            list(scores=scores, targets=tf.targets)
          }
)


setGeneric('grndynamics', function(object, tfs, p, num.genes) standardGeneric('grndynamics'))
setMethod('grndynamics',
          signature="CellRouter",
          definition=function(object, tfs, p, num.genes){
            tfs.transition <- intersect(tfs,names(cellrouter@correlation[[p]])) #intersect(tfs,names(cellrouter@top.correlations$up$SP_3.SP_8))
            tfs.transition <- cellrouter@correlation[[p]][tfs.transition]
            tfs.transition <- tfs.transition[order(tfs.transition, decreasing = TRUE)]
            tfs.transition <- tfs.transition[which(tfs.transition > 0)]
            ps <- list()
            averages <- vector()
            names <- vector()
            for(r in names(tfs.transition)){
              rgrn <- induced.subgraph(ggrn, vids=unlist(neighborhood(graph=ggrn,order=1,nodes=r)))
              genes <- intersect(V(rgrn)$name, names(cellrouter@correlation[[p]]))
              bla <- cellrouter@correlation[[p]][genes]
              bla <- bla[bla > 0]
              averages <- append(averages, mean(bla))
              names <- append(names, r)
              df <- data.frame(gene=as.vector(names(bla)), corr=as.numeric(bla), TF=rep(r, length(bla)), stringsAsFactors = FALSE)
              ps[[r]] <- df
            }
            names(averages) <- names
            averages <- averages[order(averages, decreasing=TRUE)]
            #scores <- tfs.transition[names(averages)] * averages
            #scores <- scores[order(scores, decreasing=TRUE)]
            #scores <- scores[1:num.genes]
            averages <- averages[1:num.genes]
            ps <- ps[names(averages)]
            df <- do.call(rbind, ps)
            df$TF <- factor(as.vector(df$TF), levels=names(averages))
            #means <- data.frame(p=names(averages), mean=averages)
            #means$p <- factor(names(averages), levels=names(averages))
            
            ggplot(df, aes(x=TF, y=corr, fill=TF)) + geom_violin(scale="width") + stat_summary(fun.y=mean,geom='point') + #geom_boxplot(alpha=.9) + 
              #geom_boxplot(alpha=.9) +
            #ggplot(means, aes(x=p, y=mean, fill=p)) + geom_bar(stat='identity') + 
              theme_bw() + xlab("") + ylab("Correlation") + theme(legend.position="none") +
              theme(axis.text.x = element_text(size=rel(1), angle=45, hjust=1)) + ggtitle(p) +
              theme(axis.line=element_blank(),
                  #axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.y=element_blank(),
                  panel.background=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank(),
                  panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1))
            
          }
)


#I will do a simpler and stringent co-expression based GRN approach...
# setGeneric('grndynamics', function(object, tfs, p, cor, num.genes) standardGeneric('grndynamics'))
#   setMethod('grndynamics',
#             signature="CellRouter",
#             definition=function(object, tfs, p, cor, num.genes){
#               paths <- cellrouter@pathsinfo$path[p]
#               sampTab <- cellrouter@sampTab
#               sgrn <- list()
#                 subpopulations <- unique(sampTab[paths[[p]], 'population'])
#                 #genes <- names(cellroute@top.correlations$up[[p]])
#                 for(s in subpopulations){
#                   sc <- rownames(sampTab[sampTab$population == s,])
#                   sce <- cellrouter@ndata[,sc]
#                   sce <- sce[which(apply(sce, 1, var) > 0),]
#                   #tfs2 <- intersect(tfs, rownames(sce))
#                   print(dim(sce))
#                   net <- cor(t(sce))
#                   net2 <- apply(net, 1, function(x){x[abs(x) > cor]})
#                   net3 <- net2[lapply(net2, length) > num.genes]
#                   tr <- unique(c(names(net3), as.vector(unlist(lapply(net3, names)))))
#                   grn <- net[tr, tr]
#                   grn <- igraph::graph.adjacency(grn, mode='undirected', weighted = TRUE)
#                   grn <- igraph::simplify(grn, remove.loops = TRUE)
#                   V(grn)$type <- 'Target'
#                   V(grn)[intersect(V(grn)$name, tfs)]$type <- 'Regulator'
#                   grn <- ig_convertSmall(grn, exponent = 0.5)
#                   sgrn[[s]] <- grn
#                 }
#                 sgrn
#             }
# )


setGeneric("pathsinfo", function(object, num.cells, neighs) standardGeneric("pathsinfo"))
setMethod("pathsinfo",
          signature="CellRouter",
          definition=function(object, num.cells, neighs){
            
            paths <- object@paths
            expDat <- object@ndata[object@genes.trajectory,]
            sampTab <- object@sampTab
            print(neighs)
            
            #mean of neighboring cells
            # means <- list()
            # for(pop in unique(cellrouter@sampTab$population)){
            #   sx <- apply(cellrouter@ndata[,which(cellrouter@sampTab$population == pop)], 1, mean)
            #   #sx <- sx[which(sx > 0)]
            #   means[[pop]] <- sx#1-rescale(sx, c(0,1))
            # }
            # means <- do.call(cbind, means)
            # 
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
                neighs <- induced.subgraph(graph=g,vids=unlist(neighborhood(graph=g,order=o,nodes=cell)))
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
                #bla <- object@sampTab[cells,'population']
                #names(bla) <- cells
                #bb <- means[rownames(expression), bla]
                ####colnames(bb) <- colnames(expression)
                #expression <- (expression + bb)/2
                ####expression <- bb
                ##remove genes with zero variance along a path
                var <- apply(expression, 1, var)
                expression <- expression[which(var != 0),]
                
                #path_name <- paste("path_", i, sep="")
                path_name <- paths$population[i]
                print(path_name)
                
                path_info[['distr']][[path_name]] <- expression #rownames(expression) #
                path_info[['path']][[path_name]] <- cells
                path_info[['pathInfo']][[path_name]] <- data.frame(path=paste(as.vector(cells), collapse=','), 
                                                                   source=cells[1], sink=cells[length(cells)], cost=paths$path_cost[i],
                                                                   source_population=sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),'community'],
                                                                   target_population=sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), 'community'])
                pathsDF[path_name, 'source'] <- cells[1]
                pathsDF[path_name, 'sink'] <- cells[length(cells)]
                pathsDF[path_name, 'cost'] <- paths$path_cost[i]
                pathsDF[path_name, 'source_population'] <- as.character(sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),'community'])
                pathsDF[path_name, 'source_color'] <- as.character(sampTab[grep(paste('^',cells[1],'$',sep=''), rownames(sampTab)),'colors'])
                pathsDF[path_name, 'target_population'] <- as.character(sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), 'community'])
                pathsDF[path_name, 'target_color'] <- as.character(sampTab[grep(paste('^',cells[length(cells)],'$',sep=''), rownames(sampTab)), 'colors'])
                
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

#### smooth dynamics along trajectories to help cluster gene expression trends
setGeneric("smoothdynamics2", function(object, names) standardGeneric("smoothdynamics2"))
setMethod("smoothdynamics2",
          signature="CellRouter",
          definition=function(object, names){
            
            dynamics <- list()
            expDat <- object@pathsinfo$distr
            geneList <- cellrouter@genes.trajectory
            for(path in names(expDat[names])){
              print(path)
              x_axis <- as.numeric(object@pathsinfo$pseudotime[[path]])
              geneList <- rownames(expDat[[path]][geneList,])
              geneList <- intersect(geneList, rownames(expDat[[path]]))
              smoothDynamics <- data.frame(matrix(0, nrow=length(geneList), ncol=1001))
              rownames(smoothDynamics) <- geneList
              for(gene_id in geneList){
                y_axis <- as.numeric(expDat[[path]][gene_id, ])
                lo <- loess(y_axis~x_axis)
                xl <- seq(min(x_axis),max(x_axis), (max(x_axis) - min(x_axis))/1000)
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


#### smooth dynamics along trajectories to help cluster gene expression trends
setGeneric("smoothdynamics", function(object, names) standardGeneric("smoothdynamics"))
setMethod("smoothdynamics",
          signature="CellRouter",
          definition=function(object, names){
            
            dynamics <- list()
            #ndata <- object@ndata
            expDat <- object@pathsinfo$distr
            geneList <- object@genes.trajectory
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

###cluster genes along trajectories
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
              x <- clustergenes(matrix, num.clusters)
              clusters[[trajectory]] <- x
            }
            object@clusters <- clusters
            
            return(object)
          }
)

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

setGeneric("plotGenesInClusters", function(rankedDynamics, traj.name, num.genes, columns=5, width=23, height=15, filename, ...) standardGeneric("plotGenesInClusters"))
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

setGeneric("correlationpseudotime2", function(object, type) standardGeneric("correlationpseudotime2"))
setMethod("correlationpseudotime2",
          signature="CellRouter",
          definition=function(object, type){
            xx <- object@pathsinfo
            genelist <- cellrouter@genes.trajectory
            print('computing correlation with the pseudotime')
            values <- list()
            for(t in names(xx)){
              cat(t, '\n')
              pathsInfo <- xx[[t]]  
              correlations <- list()
              for(path in names(pathsInfo[['distr']])){
                #cat(path, '\n')
                x <- pathsInfo[['pseudotime']][[path]]
                genes <- intersect(genelist, rownames(pathsInfo[['distr']][[path]]))
                for(gene in genes){
                  y <- as.numeric(pathsInfo[['distr']][[path]][gene,])
                  if(type == 'slope'){
                    df <- data.frame(x=x, y=y)
                    df <- lm(y ~ x, data=df)
                    correlations[[path]][[gene]] <- df$coefficients[2]
                  }else if(type == 'pearson'){
                    cors <- cor.test(x, y, method="pearson")
                    correlations[[path]][[gene]] <- as.numeric(cors$estimate)
                  }else{
                    cors <- cor.test(x, y, method="spearman")
                    correlations[[path]][[gene]] <- as.numeric(cors$estimate)
                  }
                }
              }
              values[[t]] <- correlations
            }
            #object@correlation <- correlations
            object@correlation <- values
            
            return(object)
          }
)


setGeneric("correlationpseudotime", function(object, type) standardGeneric("correlationpseudotime"))
setMethod("correlationpseudotime",
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
                  cors <- cor.test(x, y, method="pearson")
                  correlations[[path]][[gene]] <- as.numeric(cors$estimate)
                }else{
                  cors <- cor.test(x, y, method="spearman")
                  correlations[[path]][[gene]] <- as.numeric(cors$estimate)
                }
                
              }
            }
            object@correlation <- correlations
            
            return(object)
          }
)

### Find top genes more highly correlated with pseudotime
setGeneric("topgenes2", function(object, max.quantile, min.quantile) standardGeneric("topgenes2"))
setMethod("topgenes2",
          signature="CellRouter",
          definition=function(object, max.quantile, min.quantile){
            slopes <- object@correlation
            geneTrends <- list()
            plots <- list()
            
            for(transition in names(slopes)){
              for(path in names(slopes[[transition]])){
                threshold.up <- quantile(slopes[[transition]][[path]], max.quantile, na.rm=TRUE)
                threshold.down <- quantile(slopes[[transition]][[path]], min.quantile, na.rm=TRUE)
                
                x <- slopes[[transition]][[path]]
                genesUP <- x[which(x > threshold.up)]
                genesDOWN <- x[which(x < threshold.down)]
                geneTrends[[transition]][['up']][[path]] <- genesUP[order(genesUP, decreasing=TRUE)]
                geneTrends[[transition]][['down']][[path]] <- genesDOWN[order(genesDOWN)]
              }
            }
            object@top.correlations <- geneTrends
            return(object)
          }
)


### Find top genes more highly correlated with pseudotime
setGeneric("topgenes", function(object, max.quantile, min.quantile) standardGeneric("topgenes"))
setMethod("topgenes",
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

#create a data.frame of genes by paths with the correlation with pseudotime
setGeneric("pathsdf", function(object, threshold=0.7) standardGeneric("pathsdf"))
setMethod("pathsdf",
          signature="CellRouter",
          definition=function(object, threshold){
            
            cors <- object@top.correlations
            all.genes <- vector()
            for(path in names(cors[['up']])){
              genes.up <- names(cors[['up']][[path]])
              genes.down <- names(cors[['down']][[path]])
              x <- c(genes.up, genes.down)
              all.genes <- append(all.genes, x)
            }
            all.genes <- unique(all.genes)
            df <- data.frame()
            for(gene in all.genes){
              for(path in names(object@correlation)){
                df[gene, path] <- object@correlation[[path]][gene]
              }
            }
            df[is.na(df)] <- 0 #replace NAs by zero
            #df[df > threshold] <- threshold
            #df[df < -threshold] <- -threshold
            
            df
          }
)

## heatmap of genes by paths with correlation with pseudotime
setGeneric("plotheatmap", function(object, paths, threshold, crows=5, ccols=1, width=1, height=3.5, filename="results/cor.path_heatmap.pdf") standardGeneric("plotheatmap"))
setMethod("plotheatmap",
          signature="CellRouter",
          definition=function(object, paths, threshold, crows, ccols, width, height, filename){
            df <- pathsdf(object, threshold)
            
            pathsInfo <- object@pathsinfo
            ann <- pathsInfo$path_data
            
            df <- df[,paths]
            ann <- ann[paths,]
            source_colors <- unique(ann$source_color)
            names(source_colors) <- unique(ann$source_population)
            target_colors <- unique(ann$target_color)
            names(target_colors) <- unique(ann$target_population)
            ann_colors = list(
              source_population = source_colors,
              target_population = target_colors
            )
            names <- paste(colnames(df), as.vector(ann$source_population), as.vector(ann$target_population), sep="_")
            paletteLength <- 100
            myColor <- colorRampPalette(c("darkblue","black","yellow"))(paletteLength)
            myBreaks <- c(seq(min(df), 0, length.out=ceiling(paletteLength/2) + 1), 
                          seq(max(df)/paletteLength, max(df), length.out=floor(paletteLength/2)))
            
            a <- pheatmap(df, color = myColor, breaks = myBreaks,  clustering_method = 'ward.D',
                          annotation_col=ann[,c(4,6)], cutree_cols = ccols, show_rownames = FALSE, show_colnames=TRUE,
                          annotation_colors=ann_colors, clustering_distance_rows="correlation", 
                          clustering_distance_cols="correlation", border_color="black", border=TRUE, fontsize=8,
                          cluster_rows = TRUE, cluster_cols = TRUE, labels_col = names, width=width, height=height, filename=filename)
            
          df
          }
)

## heatmap of genes by paths with correlation with pseudotime
setGeneric("plotheatmapCategory", function(object, threshold, genelist, crows=1, ccols=1, width=6, height=3.5, filename="results/cor.path_heatmap.pdf") standardGeneric("plotheatmapCategory"))
setMethod("plotheatmapCategory",
          signature="CellRouter",
          definition=function(object, threshold, genelist, crows, ccols, width, height, filename){
            df <- pathsdf(cellrouter, threshold)
            df <- df[intersect(rownames(df), genelist),]
            pathsInfo <- object@pathsinfo
            ann <- pathsInfo$path_data
            source_colors <- unique(ann$source_color)
            names(source_colors) <- unique(ann$source_population)
            target_colors <- unique(ann$target_color)
            names(target_colors) <- unique(ann$target_population)
            ann_colors = list(
              source_population = source_colors,
              target_population = target_colors
            )
            names <- paste(colnames(df), as.vector(ann$source_population), as.vector(ann$target_population), sep="_")
            paletteLength <- 100
            myColor <- colorRampPalette(c("lightblue","black","yellow"))(paletteLength)
            myBreaks <- c(seq(min(df), 0, length.out=ceiling(paletteLength/2) + 1), 
                          seq(max(df)/paletteLength, max(df), length.out=floor(paletteLength/2)))
            
            a <- pheatmap(df, color = myColor, breaks = myBreaks,  clustering_method = 'ward.D',
                          annotation_col=ann[,c(4,6)], cutree_cols = ccols, show_rownames = TRUE, show_colnames=TRUE,
                          annotation_colors=ann_colors, clustering_distance_rows="correlation", 
                          clustering_distance_cols="correlation", border_color="black", border=TRUE, fontsize=8,
                          cluster_rows = TRUE, cluster_cols = TRUE, labels_col = names, width=width, height=height, filename=filename)
            
            df
          }
)


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

#plot a heatmap of top genes ordered by pseudotime, data is centered and thresholded
#for visualization
setGeneric("plotPathHeatmap", function(object, num.genes, logTransform=TRUE, threshold=2, width, height, dir) standardGeneric("plotPathHeatmap"))
setMethod("plotPathHeatmap", 
          signature="CellRouter",
          definition=function(object, num.genes, logTransform, threshold, width, height, dir){
            
            #plotPathHeatmap <- function(corsPaths, pathsInfo, graph, num_genes, width, height, dir){
            corsPaths <- object@top.correlations
            pathsInfo <- object@pathsinfo
            sampTab <- object@sampTab
            
            for(path in names(corsPaths$up)){
              tmpexpr <- pathsInfo$distr[[path]][c(names(corsPaths$up[[path]])[1:num.genes],
                                                   names(corsPaths$down[[path]])[1:num.genes]),]
              if(logTransform){
                tmpexpr <- log(tmpexpr + 1)
              }
              andf <- data.frame(sampTab[pathsInfo$path[[path]], 'community',])
              rownames(andf) <- pathsInfo$path[[path]]
              colnames(andf) <- c('subpopulation')
              
              target_colors <- unique(sampTab[pathsInfo$path[[path]], 'colors',])
              names(target_colors) <- unique(andf$subpopulation)
              ann_colors = list(
                subpopulation = target_colors
              )
              from <- sapply(strsplit(path, split='_', fixed=TRUE), function(x){x[2]})
              to <- sapply(strsplit(path, split='_', fixed=TRUE), function(x){x[4]})
              title <- paste('Transition ', from, ' ', to, sep='')
              file <- paste(dir, 'heatmap_top_', num.genes,'_', path, '.pdf', sep='')
              labels <- sapply(strsplit(rownames(tmpexpr), split='__', fixed=TRUE), function(x){x[1]})
              pheatmap(center_with_threshold(tmpexpr, threshold), cluster_rows = FALSE, 
                       cluster_cols = FALSE, annotation_col = andf, annotation_colors=ann_colors, 
                       show_colnames = FALSE, border=FALSE, main=title, filename = file, 
                       width = width, height=height, labels_row = labels)
            }
          }
)

### plot genes along a trajectory
setGeneric("plottrajectory", function(object, trajectory, geneList, columns=5, width=23, height=15, filename, ...) standardGeneric("plottrajectory"))
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

setGeneric("plottrajectories", function(object, trajectories, geneList, rescale, columns=5, width=23, height=15, filename, ...) standardGeneric("plottrajectories"))
setMethod("plottrajectories",
          signature="CellRouter",
          definition=function(object, trajectories, geneList, rescale, columns, width, height, filename){
            
            library(smoother)
            
            plotlist <- list()
            for(t in trajectories){
              plots <- list()
              trajectory <- object@pathsinfo$path[[t]]
              x_axis <- 1:length(trajectory) #trajectory
              for(gene_id in geneList){
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
              g1 <- ggplot(tables, aes(x=cells, y=Expression, group=gene, colour=gene)) +
                theme_bw() + geom_line(size=1) + xlab('CellRouter trajectory') +
                guides(col=guide_legend(direction="vertical")) + #, nrow = 2
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      axis.text.x=element_blank(), axis.ticks=element_blank(),
                      legend.position = "right",
                      panel.border = element_blank()) +
                theme(axis.line.x = element_line(color="black", size = 0.5),
                      axis.line.y = element_line(color="black", size = 0.5))+
                scale_color_manual("", values=rainbow(length(geneList))) +
                ggtitle(t)
              #scale_color_brewer("", palette = 'Set1')
              plotlist[[t]] <- g1
              
            }
            pdf(filename, width=width, height=height)
            #print(g1)
            multiplot(plotlist = plotlist, cols = columns)
            dev.off()
            print(g1)
            
            
          }
)


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
                  df$population <- factor(df$population, levels=df$population)
                  num_subpops <- length(unique(df$population))
                  colors <- sampTab[names(y_axis), 'colors']
                  names(colors) <- as.vector(df$population)
                  
                  g1 <- ggplot(df, aes(x=cells, y=Expression, group=1, colour=population)) +
                    geom_point(size=5) + stat_smooth() + theme_bw() + ggtitle(paste(gene_id, path, sep='--')) +
                    scale_colour_manual(values = colors) + xlab("Trajectory") +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          axis.text.x=element_blank(),
                          panel.background=element_blank(),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
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

setGeneric("plot3DExpression", function(object, genelist, logTransform=TRUE, parameters, filename,...) standardGeneric("plot3DExpression"))
setMethod("plot3DExpression",
          signature="CellRouter",
          definition=function(object, genelist, logTransform, parameters,filename){
            matrix <- object@rdimension
            for(gene in genelist){
              expr <- object@ndata[gene,rownames(matrix)]
              cols <- cGrad(as.numeric(object@ndata[gene,rownames(matrix)]))
              if(logTransform){
                expr <- log(expr + 1)
              }
              plot3d(matrix[,1], matrix[,2], matrix[,3], col=cols, size=5, xlab = 'DC1', ylab='DC2', zlab='DC3')
              #text3d(matrix2[,1], matrix2[,2], matrix2[,3], labels, font=5, cex=5)
              par3d(parameters)
              file <- paste(filename, '_', gene, '.pdf',sep='')
              rgl.postscript(file,"pdf") 
            }
          }
)


## color Dimensionality Reduction plot by expression of selected genes
setGeneric("plotDRExpression", function(object, genelist, logTransform=TRUE, columns=5, width=23, height=15, filename) standardGeneric("plotDRExpression"))
setMethod("plotDRExpression",
          signature="CellRouter",
          definition=function(object, genelist, logTransform, columns=5, width=23, height=15, filename){
            
            #plotDRExpression <- function(expDat, matrix, geneList, width=10, height=3.5, num_columns=2, filename){
            matrix <- object@rdimension
            plots <- list()
            scores <- as.data.frame(matrix)
            for(gene in genelist){
              expr <- object@ndata[gene,rownames(matrix)]
              if(logTransform){
                expr <- log(expr + 1)
              }
              scores$GENE <- as.numeric(expr)
              p1 <- ggplot(scores,aes(x = tSNE1, y=tSNE2, colour=GENE)) + geom_point(size=0.1) + theme_bw() +
                ggtitle(gene)  + scale_colour_gradientn(name = "Expression", colours=rev(rainbow(4))) +
                theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
                theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
                theme(panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
                theme(legend.position="right")
              
              plots[[gene]] <- p1
            }
            pdf(file=filename, width=width, height=height)
            multiplot(plotlist = plots, cols=columns)
            dev.off();
            multiplot(plotlist = plots, cols=columns)
            
          }
)

setGeneric("davidclustering", function(object, names, tmpdir, currdir, email, ids) standardGeneric("davidclustering"))
setMethod("davidclustering",
          signature="CellRouter",
          definition=function(object, names, tmpdir, currdir, email, ids){
            
            library(FGNet)
            library(rJava)
            library(RDAVIDWebService)
            
            enrichment <- list()
            #API_defaults <- c(overlap=3L, initialSeed=3L, finalSeed=3L, linkage=0.5, kappa=50L)
            for(transition in names(object@clusters[names])){
              clusters <- object@clusters[[transition]]$clustering
              genes.clusters <- split(clusters, clusters)
              genes.clusters <- lapply(genes.clusters, function(x){names(x)})
              geneList <- lapply(genes.clusters, function(x){convertIDs(ids, x,from='external_gene_name', to="ensembl_gene_id")})
              
              setwd(tmpdir)
              david.results <- lapply(geneList, function(x){
                fea_david(names(x), 
                          annotations=c("GOTERM_BP_ALL",
                                        'KEGG_PATHWAY',
                                        "REACTOME_PATHWAY"),
                                         geneLabels=x,
                                         email=email,
                                         argsWS=c(overlap=3L, initialSeed=3L, finalSeed=3L, linkage=0.5, kappa=50L))
              })
              enrichment[[transition]] <- david.results
              setwd(currdir)
              #save(david.results, file=paste(dir, population, '_DAVID_clustering.R', sep=''))
            }
            object@davidenrichment <- enrichment
            return(object)
          }
          
)
###plot DAVID clustering results for the different gene clusters
setGeneric("plotenrichment", function(object, threshold=2, logTransform=2, width=23, height=15, dir) standardGeneric("plotenrichment"))
setMethod("plotenrichment",
          signature="CellRouter",
          definition=function(object, threshold=2, logTransform=2, width=23, height=15, dir){
            enrichment <- object@davidenrichment
            for(transition in names(enrichment)){
              david.results <- enrichment[[transition]]
              
              david.results.clean <- lapply(david.results, function(x){x$clusters[which(as.numeric(as.character(x$clusters$ClusterEnrichmentScore)) > threshold),]})
              #f <- paste(dir, '_DAVID_enrichment_', transition, '.csv')
              #write.csv(david.results.clean, file=f)   
              
              all.terms <- unique(as.vector(unlist(lapply(david.results.clean, function(x){as.vector(x$keyWordsTerm)}))))
              df <- data.frame(matrix(0, nrow=length(david.results), ncol=length(all.terms)))
              rownames(df) <- names(david.results.clean)
              colnames(df) <- all.terms#mgKeyTerm#
              for(p in rownames(df)){
                for(e in colnames(df)){
                  x <- david.results.clean[[p]]#david.results[[p]][['clusters']]
                  xx <- x[which(as.vector(x$keyWordsTerm) == e),]
                  if(nrow(xx) == 1){
                    df[p, e] <- as.numeric(as.character(xx$ClusterEnrichmentScore))
                  }else if(nrow(xx) > 1){
                    print(paste('Different clusters with same keyWordsTerm: ', e, sep=''))
                    xx.mean <- mean(as.numeric(as.character(xx[,'ClusterEnrichmentScore'])))
                    df[p, e] <- xx.mean
                  }
                }
              }
              if(logTransform){
                df <- log(df+1)
              }
              pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE, clustering_method = 'ward.D',
                       width=width, height=height, filename=paste(dir, transition, '_DAVID_clustering.pdf', sep=''),
                       border_color = 'black')
              
            }
          }
)

setGeneric("pathwayenrichment", function(object, names, species, annotation, ids) standardGeneric("pathwayenrichment"))
setMethod("pathwayenrichment",
          signature="CellRouter",
          definition=function(object, names, species, annotation, ids){
            
            ##up-regulated
            upregulated <- list()
            geneNames <- lapply(object@top.correlations$up, names)
            geneNames <- lapply(geneNames, function(x){setdiff(x, cc)}) #remove cell cycle genes
            geneList <- lapply(geneNames, function(x){convertIDs(ids, x, from='external_gene_name', to="entrezgene")})
            geneList <- lapply(geneList, names)
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
            
            pathways
          }
)

#CellRouter_utils_symbols.R/ cellrouter_analysis_waterfall_datasets.R
correlationdistribution <- function(object, selected.paths, method, pseudotime, direction, colors, width, height, filename){
#correlationdistribution <- function(corsPaths,selected_paths, overallPseudotime, num_paths, direction, colors, width, height, filename){
  corsPaths <- object@top.correlations
  x <- list()
  for(path in names(corsPaths[[direction]][selected.paths])){
    if(method=='Monocle'){
      x[[path]] <- data.frame(CellRouter=corsPaths[[direction]][[path]], Monocle=pseudotime[names(corsPaths[[direction]][[path]])],
                            pathid=rep(path, length(corsPaths[[direction]][[path]])))
    }else if(method=='Waterfall'){
    x[[path]] <- data.frame(CellRouter=corsPaths[[direction]][[path]], Waterfall=pseudotime[names(corsPaths[[direction]][[path]])],
                            pathid=rep(path, length(corsPaths[[direction]][[path]])))
    }
  }
  x <- do.call(rbind, x)
  df <- x

  df.m <- melt(df, id.vars=c('pathid'))
  pdf(file=filename, width=width, height=height)
  g <- ggplot(df.m, aes(x=pathid, y=value, fill=variable)) + geom_boxplot(alpha=.9) + #geom_violin(scale="width", colour="black") +
    theme_bw() + xlab("") + ylab("Correlation with pseudotime") + theme(legend.position="none") +
    #theme(axis.text.x = element_text(size=rel(1), angle=00, hjust=1)) +
    theme(axis.line=element_blank(),
          #axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.y=element_blank(),
          legend.position="top",
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
    scale_fill_manual(name="", values=c(colors)) + coord_flip()
  print(g)
  dev.off()
}


meanexpressionsubpopulations <- function(object, geneList, column, filename, cols, width=10, height=5){
  plots <- list()
  sampTab <- object@sampTab
  graph <- object@graph
  expDat <- object@ndata
  T0 <- expDat
  for(g in geneList){
    averages <- vector()
    name.averages <- vector()
    for(i in unique(sampTab$population)){
      cells.population <- rownames(sampTab[which(sampTab$population == i),])
      p <- mean(as.numeric(object@ndata[g, cells.population]))
      averages <- append(p, averages)
      name.averages <- append(i, name.averages)
    }
    names(averages) <- name.averages
    averages <- averages[order(averages, decreasing=FALSE)]
    means <- data.frame(p=names(averages), mean=averages)
    means$p <- factor(names(averages), levels=names(averages))
    
    colors <- unique(sampTab$colors)
    names(colors) <- unique(sampTab$population)
    #p <- ggplot(genes.m, aes(x=conditions, y=value, fill=conditions)) + geom_boxplot(alpha=.9, notch = TRUE) + 
    p <- ggplot(means, aes(x=p, y=mean, fill=p)) + 
      geom_bar(stat ="identity") +
      #geom_jitter(position='identity', size=1) + 
      theme_bw() + xlab("") + ylab("Gene expression") + theme(legend.position="none") +
      theme(axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x=element_blank(),axis.ticks=element_blank(),
            panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) + ggtitle(g) +
      scale_fill_manual("", values=colors)
    
    plots[[g]] <- p
  }
  #pdf(file=filename, width=10, height=7)
  pdf(file=filename, width=width, height=height)
  multiplot(plotlist = plots, cols=cols)
  dev.off();
  multiplot(plotlist = plots, cols=cols)
}


#expressionsubpopulation <- function(expDat, geneList, graph, column, filename, cols, width=10, height=5){
expressionsubpopulations <- function(object, geneList, column, filename, cols, width=10, height=5){
  plots <- list()
  sampTab <- object@sampTab
  graph <- object@graph
  expDat <- object@ndata
  T0 <- expDat
  for(g in geneList){
    #cat(time, ' ', dim(T0), '\n')
    genes <- as.data.frame(t(T0[g,]))
    genes$gene <- g
    genes$conditions <- as.vector(sampTab[,column])
    genes.m <- melt(genes, id.var=c('gene',"conditions"))
    
    averages <- vector()
    name.averages <- vector()
    for(i in unique(sampTab$population)){
      cells.population <- rownames(sampTab[which(sampTab$population == i),])
      p <- mean(as.numeric(object@ndata[g, cells.population]))
      averages <- append(p, averages)
      name.averages <- append(i, name.averages)
    }
    names(averages) <- name.averages
    averages <- averages[order(averages, decreasing=FALSE)]
    genes.m$conditions <- factor(genes.m$conditions, levels=names(averages))
    colors <- unique(sampTab$colors)
    names(colors) <- unique(sampTab$population)
    #p <- ggplot(genes.m, aes(x=conditions, y=value, fill=conditions)) + geom_boxplot(alpha=.9, notch = TRUE) + 
    p <- ggplot(genes.m, aes(x=conditions, y=value, fill=conditions)) + 
      geom_violin(scale="width") + stat_summary(fun.y=mean,geom='point') + #geom_boxplot(alpha=.9) + 
      #geom_jitter(position='identity', size=1) + 
      theme_bw() + xlab("") + ylab("Gene expression") + theme(legend.position="none") +
      theme(axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x=element_blank(),axis.ticks=element_blank(),
            panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) + ggtitle(g) +
      scale_fill_manual("", values=colors)
    
    plots[[g]] <- p
  }
  #pdf(file=filename, width=10, height=7)
  pdf(file=filename, width=width, height=height)
  multiplot(plotlist = plots, cols=cols)
  dev.off();
  multiplot(plotlist = plots, cols=cols)
}




convertID <- function(ids, ens, from, to){
  symbol <-  ids[which(ids[,from]==ens), to]
  if(length(symbol) == 0){
    ens
  }else{
    symbol[1]
  }
}

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

## Assign names to subpopulations
commToNames<-function(commObj, prefix){
  ans<-list();
  comms<-communities(commObj);
  for(i in seq(length(comms))){
    #nname<-paste(prefix,"_sn_",i,sep='');
    nname<-paste(prefix, "_", i,sep='');
    #ans[[nname]]<-commObj$names[comms[[i]]];
    ans[[nname]]<-comms[[i]]; #new igraph version
  }
  ans;
}

##Create colors from numeric values
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(rev(rainbow(4)))(range01(x))
  #cols <- colorRamp(c('deepskyblue3', 'white', 'red'))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}
cRamp2 <- function(x){
  #cols <- colorRamp(rev(rainbow(4)))(range01(x))
  cols <- colorRamp(c('white', "goldenrod1","darkorange2"))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}
cRampClust <- function(x, num_colors){
  library(RColorBrewer)
  cols <- colorRamp(rev(rainbow(num_colors)))(range01(x))
  #cols <- colorRamp(c('deepskyblue3', 'lightblue', 'orange','darkred'))(range01(x))
  #color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  #cols <- color[sample(num_colors)]
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
  #cols
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



averageIds <- function(expDat,geneTab,nameCol="symbol"){
  if(!is.matrix(expDat)){
    expDat<-as.matrix(expDat);
  }
  rownames(geneTab)<-as.character(geneTab$id);
  sameProbes<-intersect(rownames(expDat), as.character(geneTab$id));
  
  expDat<-expDat[sameProbes,];
  geneTab<-geneTab[sameProbes,];
  
  eids<-unique(as.vector(geneTab[,nameCol]));
  uSymbols<-vector(length=length(eids));
  ans<-matrix(nrow=length(eids), ncol=ncol(expDat));
  for(i in seq(length(eids))){
    eid<-eids[i];
    xi <-  which( geneTab[,nameCol]==eid );
    desProbes <- as.character(geneTab[xi,]$id);
    if(length(xi)>1){
      ans[i,]<- apply(expDat[desProbes,], 2, mean);
    }
    else{
      ans[i,]<-expDat[desProbes,];
    }
    uSymbols[i]<-as.vector(geneTab[ xi[1] ,nameCol]);
  }
  rownames(ans)<-uSymbols;
  colnames(ans)<-colnames(expDat);
  return (as.data.frame(ans))
}

binompval <- function(p,N,n){
  pval   <- pbinom(n,round(N,0),p,lower.tail=TRUE)
  pval[!is.na(pval) & pval > 0.5] <- 1-pval[!is.na(pval) & pval > 0.5]
  return(pval)
}

center_with_threshold <- function(center_data, threshold){
  # Center data (automatically ignores zeros)
  center_data <- center_data - rowMeans(center_data, na.rm=TRUE)
  # Cap values between threshold and -threshold and recenter
  center_data[center_data > threshold] <- threshold
  center_data[center_data < (-1 * threshold)] <- -1 * threshold
  #center_data <- center_data - rowMeans(center_data, na.rm=TRUE)
  return(center_data)
}

getIndexes <- function(ann_col, ann_row){
  ann_col$ID <- as.vector(1:nrow(ann_col))
  ref_groups <- split(ann_col, as.factor(ann_col$population))
  ref_groups <- lapply(ref_groups, function(x){x$ID})
  ref_groups <- ref_groups[order(nchar(names(ref_groups)))]
  
  ### Figure out where to draw lines between subtypes in the heatmap
  ref_seps <- c()
  i_cur_idx <- 0
  order_idx <- c()
  
  for(ref_grp in ref_groups){
    i_cur_idx <- i_cur_idx + length(ref_grp)
    ref_seps <- c(ref_seps, i_cur_idx)
    order_idx <- c(order_idx, ref_grp)
  }
  ref_seps <- ref_seps[1:(length(ref_seps) - 1)]
  ref_seps <- c(ref_seps, c(0,nrow(ann_col)))
  
  ### Figure out where to draw lines between gene signatures in the heatmap
  #genesDF <- data.frame(clusters=rep(seq(1:4), times=lapply(subtypeSigs, length)), ID=1:length(comb.subtype.sigs))
  #genesDF <- data.frame(clusters=graph$sampTab$population, ID=1:length(comb.subtype.sigs))
  ann_row$ID <- as.vector(1:nrow(ann_row))
  ref_groups <- split(ann_row, as.factor(ann_row$signature))
  ref_groups <- lapply(ref_groups, function(x){x$ID})
  ref_groups <- ref_groups[order(nchar(names(ref_groups)))]
  
  ref_seps_c <- c()
  i_cur_cdx <- 0
  order_cdx <- c()
  for(ref_grp in ref_groups){
    i_cur_cdx <- i_cur_cdx + length(ref_grp)
    ref_seps_c <- c(ref_seps_c, i_cur_cdx)
    order_cdx <- c(order_cdx, ref_grp)
  }
  ref_seps_c <- ref_seps_c[1:(length(ref_seps_c) - 1)]
  ref_seps_c <- c(ref_seps_c, c(0,nrow(ann_row)))
  
  list(colsep=ref_seps, rowsep=ref_seps_c)
  
}
