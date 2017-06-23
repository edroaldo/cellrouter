require('ggplot2')
require("DESeq")
require('mclust')
require('grid')
require('scde')
require('gplots')
require('genefilter')
#require("edgeR")
#require("limma")

##### Write function to perform an in deep principal component analysis - See Nature paper on lung development... ####

#### Write the code to account for cell cycle variation in single-cell data ###
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
  
  data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations, stringsAsFactors = FALSE)
}


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


#### enrichment analysis
sigOverlapXXX<-function(sigs, geneSet, gsName, universe){
  ans<-data.frame();
  
  queryGenes<-intersect(universe, geneSet);
  nGenesQuery<-length(queryGenes);
  
  queryNames<-rep(gsName, length(sigs));
  annNames<-names(sigs);
  pvals<-rep(0, length(sigs));
  #nBoths<-rep(0, length(sigs));
  #expecteds<-rep(0, length(sigs));
  #ratios<-rep(0, length(sigs));
  genes<-rep('', length(sigs));
  nGS<-rep(0, length(sigs));
  nAnn<-rep(0, length(sigs));
  
  for(i in seq(length(sigs))){
    # set up contingency table
    # only include genes that are both in the signature AND in the universe
    annGenes <- intersect( universe, sigs[[i]]);
    
    # genesBoth = genes that are both in the query genese (module) and in the annotation (from the GMT file)
    genesBoth<-intersect(annGenes,queryGenes);
    nBoth<-length(genesBoth);
    
    genesQueryOnly <- setdiff(queryGenes, genesBoth);
    genesAnnOnly   <- setdiff(annGenes,   genesBoth);
    genesNeither   <- setdiff(universe,   c(queryGenes, annGenes));
    
    nQueryOnly     <- length(genesQueryOnly);
    nAnnOnly       <- length(genesAnnOnly  );
    nNeither       <- length(genesNeither  );
    
    mt<-matrix( c( nBoth,nQueryOnly,nAnnOnly, nNeither ), nrow=2, byrow=T);
    #print(mt)
    tResult<-fisher.test(mt);
    #tResult<-chisq.test(mt);
    
    pvals[i]<-tResult$p.value;
    #nBoths[i]<-tResult$observed[1];
    #expecteds[i]<-tResult$expected[1];
    #ratios[i]<-tResult$observed[1]/tResult$expected[1];
    genes[i]<-paste(genesBoth,collapse=',');
    nGS[i]<-length(queryGenes);
    nAnn[i]<-length(annGenes);
  }
  ans<-data.frame(queryName=queryNames,
                  querySize=nGS,
                  annName=annNames,
                  annSize=nAnn,
                  pval=pvals,
                  #nBoth=nBoths,
                  #expected=expecteds,
                  #ratio=ratios,
                  genes=genes);
  
  ans<-cbind(ans, holm=p.adjust(ans$pval, method='holm'));
  
  # re-order columns
  #ans<-ans[,c("queryName", "querySize", "annName", "annSize", "holm", "ratio", "nBoth", "expected", "genes")];
  ans<-ans[,c("queryName", "querySize", "annName", "annSize", "holm", "genes")];
  ans;
}

# runs above for multiple gene sets sigOversXXX Patrick cahan 
# Say you have a 
# (1) a gene annotation list (named list of like 'dna replication'=c(genes involved in dna rep), 'apoptosis'=c(genes involved in apoptosis), etc)
# (2) a list of gene sets (i.e. your modules)
# (3) vector of all genes on the array
# then you can test whether the overlap between each annotation and each module is greater/less than expected by chance by running:
# gobpEnrich<-sigOversXXX(geneAnnList, moduleList, allGenes)

sigOversXXX<-function(sigObj, geneLists, universe){
  ans<-data.frame();
  ans<-list();
  ai<-1;
  namesGL<-names(geneLists);
  for( i in seq(length(geneLists))){
    ans[[namesGL[i]]]<- sigOverlapXXX(sigObj,geneLists[[i]],names(geneLists)[i], universe);
    ans[[namesGL[i]]] <- ans[[namesGL[i]]][which(ans[[namesGL[i]]]$holm < 0.05),]
  }
  ans;
}

WGCNA_findBeta<-function(expDat){
  # estimate scale-free network fit for a range of Beta
  #
  # Args:
  #   expDat:
  #      matrix of gene expression data
  #
  # Returns:
  #      list of $powerEstimate (automaitcally selected beta) and $fitIndices (table of all values)
  #      prints out Rsqr fits
  #        
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  pickSoftThreshold(t(expDat), powerVector = powers, verbose = 5);
}

WGCNA_findModules<-function(expDat, beta, minModuleSize=15, deepSplit=2){
  # find modules
  #
  # Args:
  #   expDat:matrix of gene expression data
  #   beta: beta
  #   minModuleSize: 
  #   deepSplit:
  #
  # Returns:
  #      list of genes in each module  
  cat("WGCNA-ing ...");
  results <- list()
  wgcnaRes<-blockwiseModules(t(expDat),
                             maxBlockSize=9000,
                             power=beta,
                             minModuleSize=minModuleSize,
                             deepSplit=deepSplit,
                             pamRespectsDendro=FALSE,
                             mergeCutHeight=0.10,
                             numericLabels=TRUE,
                             saveTOMs=F,
                             verbose=5,
                             networkType='signed',
                             TOMType='signed');
  cat(" done\n");
  results[['result']] <- wgcnaRes
  results[['modules']] <- .convertMods(wgcnaRes, rownames(expDat));
  
  results
}

.convertMods<-function(wgcnaRes, geneVect){
  ans<-list();
  uMods<-unique(wgcnaRes$colors);
  modNames<-colnames(wgcnaRes$MEs);
  for(modName in modNames){
    suff<-strsplit(modName,"ME")[[1]][2];
    ans[[modName]]<-geneVect[which(wgcnaRes$colors==suff)];
  }
  ans;
}


plotPCAExpression <- function(expression, pca, geneList, width=10, height=3.5, num_columns=2, filename){
  plots <- list()
  scores <- as.data.frame(pca$x)
  #scores <- as.data.frame(pca@scores)
  for(gene in geneList){
    scores$GENE <- as.vector(expression[gene,])
    p1 <- ggplot(scores,aes(x = PC1, y=PC2, colour=GENE)) + geom_point(size=5) + theme_bw() +
      ggtitle(gene)  + scale_colour_gradientn(name = "Expression", colours=rev(rainbow(4)))
    plots[[gene]] <- p1
  }
  pdf(file=filename, width=width, height=height)
  multiplot(plotlist = plots, cols=num_columns)
  dev.off();
}


#### Write function to perform single-cell differential expression using SCDE ###
doSCDE <- function(){
  data(es.mef.small)
  sg <- factor(gsub("(MEF|ESC).*","\\1",colnames(es.mef.small)), levels=c("ESC", "MEF"))
  names(sg) <- colnames(es.mef.small)
  table(sg)
  
  #clean up dataset. 
  cd <- es.mef.small
  cd <- cd[rowSums(cd) > 0,]   #Omit genes never detected
  cd <- cd[colSums(cd) > 1e4,] #Omit cells with very poor coverage
  
  #Fitting error models
  n.cores <- 5
  o.ifm <- scde.error.models(counts=cd, groups=sg, n.cores=n.cores, threshold.segmentation = T, 
                             save.crossfit.plot=F, save.model.plots=F, verbose=1)
  valid.cells <- o.ifm$corr.a > 0
  table(valid.cells)
  o.ifm <- o.ifm[valid.cells,]
  o.prior <- scde.expression.prior(models=o.ifm, counts=cd, length.out = 400, show.plot=F)
  
  groups <- factor(gsub("(MEF|ESC).*","\\1", rownames(o.ifm), levels(c("ESC", "MEF"))))
  names(groups) <- rownames(o.ifm)
  
  #rundifferential expression tests on all genes
  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups=groups, n.randomizations=100, n.cores = n.cores, verbose=1)
  head(ediff[order(ediff$Z, decreasing=T),])
  write.table(ediff[order(abs(ediff$Z),decreasing=T),],file="results.txt",row.names=T,col.names=T,sep="\t",quote=F)
  
  
}
######PCA analysis using prcomp ######
doPCA <- function(stDat,
                  expDat, #typically, normalized expression values (log2(FPKM+1) in this context)
                  conditions,
                  colors,
                  symbols,
                  pcaX,
                  pcaY,
                  filename){
  
  cat("counts dimension: ", dim(expDat),"\n")
  if(length(conditions) == 2){ #include a subset of samples
    stDat <- stDat[which(stDat$conditions == conditions[1] | stDat$conditions == conditions[2]),]
  }
  expDat <- expDat[,rownames(stDat)]
  #expDat <- expDat[rowSums(expDat) > 0,]
  
  cat("expDat dimension: ", dim(expDat),"\n")
  cat("samples dimension: ", dim(stDat),"\n")
  
  #genes <- apply(expDat, 1, var)
  #ord <- order(genes, decreasing=TRUE)[1:500]
  #sc.pca <- prcomp(t(expDat[ord,]), center = TRUE, scale = FALSE)
  sc.pca <- prcomp(t(expDat), center = TRUE, scale = FALSE)
  var <- sc.pca$sdev^2
  var.pca <- var / sum(var) * 100
  
  xlab=paste(pcaX, " st component (",round(var.pca[pcaX])," % of variance)",sep="")
  ylab=paste(pcaY, " nd component (",round(var.pca[pcaY])," % of variance)",sep="")
  
  pdf(file=filename, width=5, height=5);
  plot(sc.pca$x[,c(pcaX,pcaY)], pch=symbols[stDat$conditions], col=colors[as.character(stDat$conditions)], 
       xlab=xlab, ylab=ylab, cex = 2)
  legend(100, 0.8, as.vector(unique(samples$conditions)), pch=unique(symbols[samples$conditions]),
  col=unique(colors[as.character(samples$conditions)]))
  dev.off();
  
  return (sc.pca)
}

######Intra- and inter-condition distances ######
calcIntraConditionDistances <- function(data, 
                                        description,
                                        condition){
  stDat <- data[["samples"]]
  expDat <- data[["log2.norm"]]
  stDat <- stDat[which(stDat[[description]] == condition), ]
  expDat <- expDat[,rownames(stDat)]
  dist <- data.frame(dist=as.vector(1-cor(expDat, method="pearson")))
  
  return(dist)
}

calcInterConditionDistances <- function(data, 
                                        description,
                                        condition1,
                                        condition2){
  stDat <- data[["samples"]]
  expDat <- data[["log2.norm"]]
  stDatCondition1 <- stDat[which(stDat[[description]] == condition1), ]
  expDatCondition1 <- expDat[,rownames(stDatCondition1)]
  
  stDatCondition2 <- stDat[which(stDat[[description]] == condition2), ]
  expDatCondition2 <- expDat[,rownames(stDatCondition2)]
  
  dist <- data.frame(dist=as.vector(1-cor(expDatCondition1, expDatCondition2, method="pearson")))
  
  return(dist)
}

plotDistances <- function(condition1, 
                          condition2,
                          label,
                          filename){
  
    pdf(file=filename, width=3, height=2);
    p <- ggplot() + geom_density(alpha=.5, aes(x=dist), fill="red", colour="red", data=condition1) + 
      geom_density(alpha=.5, aes(x=dist), fill="blue", colour="blue", data=condition2) + theme_bw() +
      theme(legend.position=c(0.7, 0.7), legend.text=element_text(size=10)) + xlab("Distance") + ylab("Density") + 
      theme(axis.text.x = element_text(size=10, angle=00, hjust=1), axis.title.y = element_text(size = rel(0.9), angle = 90), 
      panel.border = element_blank(), axis.line = element_line()) + ggtitle(label)
    
    print(p)
    dev.off();
  
}

##### single-cell state classification #####
stateClassification <- function(expDat,
                                myData,
                                colors.class,
                                symbols.class,
                                filename){
  myClust <- Mclust(myData) #pay attention to unclassified cells. Probably column 'z' in myClust object.
  #plot(myClust, what = c("BIC", "classification"))
  #table(myClust$classification)
  
  data <- data.frame(cells=expDat[["samples"]][["conditions"]], classification=as.character(myClust$classification))
  file <- paste(filename, "_proportion.pdf",sep="")
  pdf(file=file, width=3, height=3);
  p <- ggplot(data,aes(x = cells, fill = classification)) + 
    geom_bar(position = "fill") + theme_bw() + theme(legend.position='right', legend.key.size = unit(0.3, "cm"), legend.text=element_text(size=7)) +
    xlab("") + ylab("Proportion") + 
    theme(axis.text.x = element_text(size=10, angle=45, hjust=1), axis.title.y = element_text(size = rel(1), angle = 90)) +
    scale_fill_manual("", values=c(colors.class))
    print(p)
  dev.off()
  
  #### PCA result coloring by most likely state classification ######
  file <- paste(filename, "_PCA.pdf",sep="")
  pdf(file=file, width=5, height=5);
  plot(sc.pca$x[,c(1,2)], pch=symbols.class[data$classification], col=colors.class[as.character(data$classification)], cex = 2)
  dev.off();
  
  return (myClust)
}

###### Plot DESeq-genearted principal component analysis ########
doDESeq_PCA <- function(stDat,
                        counts,
                        conditions,
                        filename){
  cat("counts dimension: ", dim(counts),"\n")
  if(length(conditions) == 2){ #include a subset of samples
    stDat <- stDat[which(stDat$conditions == conditions[1] | stDat$conditions == conditions[2]),]
  }
  
  aux.counts <- counts
  aux.counts <- aux.counts[,rownames(stDat)]
  aux.counts <- aux.counts[rowSums(aux.counts) > 0, ]
  aux.counts <- aux.counts[,colSums(aux.counts) > 1e4]
  stDat <- stDat[colnames(aux.counts), ]
  
  cat("aux.counts dimension: ", dim(aux.counts),"\n")
  cat("samples dimension: ", dim(stDat),"\n")
  
  cds <- newCountDataSet(aux.counts, stDat$conditions)
  cds <- estimateSizeFactors(cds)
  cdsB <- estimateDispersions(cds, method="blind")
  vsd <- varianceStabilizingTransformation(cdsB)
  
  pdf(file=filename, width=5, height=5)
  print(plotPCA(vsd))
  dev.off()
}

doDESeq <- function(stDat,
                    counts,
                    conditions){
  
  cat("Comparing ", conditions[2], "to", conditions[1], '\n')
  cat('sample dimensions: ', dim(stDat), '\n')
  cat('counts dimensions: ', dim(counts), '\n')
  
  aux.counts <- counts[,colSums(counts) > 3e4]
  stDat <- stDat[colnames(aux.counts), ]
  cds <- newCountDataSet(aux.counts, factor(as.vector(stDat$conditions)))
  #use <- (rowSums(counts(cds)) > quantile(rowSums(counts(cds)), probs=0.40))
  #cds <- cds[use,]
  
  cds <- estimateSizeFactors(cds)
  #cds<-estimateDispersions(cds, method="per-condition", fitType="local", sharingMode="maximum")
  cds<-estimateDispersions(cds) #for Cynthia
  results <- nbinomTest(cds, conditions[1], conditions[2])
  rownames(results) <- results$id
  sig.results<-results[results$padj < 0.05 & !is.na(results$padj) & !is.infinite(results$log2FoldChange),]
  sig.results <- sig.results[order(sig.results$log2FoldChange, decreasing=TRUE), ]
  data <- list()
  data[["raw"]] <- results
  data[["clean"]] <- sig.results
  data[['clean_results']] <- counts(cds)
  
  return (data)
}

doDESeqBlind <- function(stDat,
                         counts,
                         column,
                         conditions){
  
  cat("Comparing ", conditions[2], "to", conditions[1], '\n')
  cat('sample dimensions: ', dim(stDat), '\n')
  cat('counts dimensions: ', dim(counts), '\n')
  
  aux.counts <- counts[,colSums(counts) > 3e4]
  stDat <- stDat[colnames(aux.counts), ]
  cds <- newCountDataSet(aux.counts, factor(as.vector(stDat[,column])))
  #use <- (rowSums(counts(cds)) > quantile(rowSums(counts(cds)), probs=0.40))
  #cds <- cds[use,]
  
  cds <- estimateSizeFactors(cds)
  cds<-estimateDispersions(cds, method="blind", fitType="local", sharingMode="fit-only")
  results <- nbinomTest(cds, conditions[1], conditions[2])
  rownames(results) <- results$id
  sig.results<-results[results$padj < 0.05 & !is.na(results$padj) & !is.infinite(results$log2FoldChange),]
  sig.results <- sig.results[order(sig.results$log2FoldChange, decreasing=TRUE), ]
  data <- list()
  data[["raw"]] <- results
  data[["clean"]] <- sig.results
  data[['clean_results']] <- counts(cds)
  
  return (data)
}


doDESeqBlind <- function(stDat,
                    counts,
                    conditions){
  
  cat("Comparing ", conditions[2], "to", conditions[1], '\n')
  cat('sample dimensions: ', dim(stDat), '\n')
  cat('counts dimensions: ', dim(counts), '\n')
  
  aux.counts <- counts[,colSums(counts) > 3e4]
  stDat <- stDat[colnames(aux.counts), ]
  cds <- newCountDataSet(aux.counts, factor(as.vector(stDat$conditions)))
  #use <- (rowSums(counts(cds)) > quantile(rowSums(counts(cds)), probs=0.40))
  #cds <- cds[use,]
  
  cds <- estimateSizeFactors(cds)
  cds<-estimateDispersions(cds, method="blind", fitType="local", sharingMode="fit-only")
  results <- nbinomTest(cds, conditions[1], conditions[2])
  rownames(results) <- results$id
  sig.results<-results[results$padj < 0.05 & !is.na(results$padj) & !is.infinite(results$log2FoldChange),]
  sig.results <- sig.results[order(sig.results$log2FoldChange, decreasing=TRUE), ]
  data <- list()
  data[["raw"]] <- results
  data[["clean"]] <- sig.results
  data[['clean_results']] <- counts(cds)
  
  return (data)
}


norm.DESeq <- function(stDat,
                       counts){
  
  aux <- list()
  aux.counts <- counts[,colSums(counts) > 3e4]
  #aux.counts <- counts[,colSums(counts) > 1e3]
  #aux.counts <- counts
  stDat <- stDat[colnames(aux.counts), ]
  cds <- newCountDataSet(aux.counts, factor(as.vector(stDat$conditions)))
  
  cds <- estimateSizeFactors(cds)
  aux[['samples']] <- stDat #sample table after filtering samples with overall ounts < 3e4
  aux[['norm']] <- counts(cds, normalized=TRUE) #normalized counts
  aux[['log2.norm']] <- log2(counts(cds, normalized=TRUE) + 1)  #log2 of normalized counts
  
  cdsB <- estimateDispersions(cds, method="blind", fitType="local")
  vsd <- getVarianceStabilizedData(cdsB)
  aux[['var.norm']] <- vsd
  
  return (aux)
}

detectGenes <- function(expr, min_expr=0.1){
  detected <- list()  
  cells_detected <- do.call(rbind, apply(expr, 2, function(x){return (data.frame(num_genes_detected=sum(x > min_expr)))}))
  genes_detected <- do.call(rbind, apply(expr, 1, function(x){return (data.frame(num_cells_detected=sum(x > min_expr)))}))
  detected[['genes_detected']] <- genes_detected
  detected[['cells_detected']] <- cells_detected
  
  detected
}


cn_findSpecGenes<-function# find genes that are preferentially expressed in specified samples
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
    ctGenes[[ctName]]<-tmp;
    cvalT<-append(cvalT, cval);
  }
  #specificSets[['306G']][1:5,]
  if(remove){
    cat("Prune common genes...\n");
    # now limit to genes exclusive to each list
    specGenes<-list();
    for(ctName in ctNames){
      others<-setdiff(ctNames, ctName);
      x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
      specGenes[[ctName]]<-x;
    }
    result <- specGenes
  }else {
    result <- ctGenes;
  }
  result
}

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
    ctGenes[[ctName]]<-tmp;
    cvalT<-append(cvalT, cval);
  }
  if(remove){
    cat("Prune common genes...\n");
    # now limit to genes exclusive to each list
    specGenes<-list();
    for(ctName in ctNames){
      others<-setdiff(ctNames, ctName);
      x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
      specGenes[[ctName]]<-x;
    }
    result <- specGenes
  }else {
    result <- ctGenes;
  }
  results <- list()
  results[['genes']] <- result
  results[['table']] <- specificSets
  
  results
}

# cn_findSpecGenes<-function# find genes that are preferentially expressed in specified samples
# (expDat, ### expression matrix
#  sampTab, ### sample table
#  qtile=0.95, ### quantile
#  dLevel="population_name" #### annotation level to group on
# ){
#   
#   cat("Template matching...\n")
#   myPatternG<-cn_sampR_to_pattern(as.vector(sampTab[,dLevel]));
#   specificSets<-apply(myPatternG, 1, cn_testPattern, expDat=expDat);
#   
#   # adaptively extract the best genes per lineage
#   cat("First pass identification of specific gene sets...\n")
#   cvalT<-vector();
#   ctGenes<-list();
#   ctNames<-unique(as.vector(sampTab[,dLevel]));
#   for(ctName in ctNames){
#     x<-specificSets[[ctName]];
#     cval<-quantile(x$cval, qtile);
#     tmp<-rownames(x[x$cval>cval,]);
#     ctGenes[[ctName]]<-tmp;
#     cvalT<-append(cvalT, cval);
#   }
#   
#   cat("Prune common genes...\n");
#   # now limit to genes exclusive to each list
#   specGenes<-list();
#   for(ctName in ctNames){
#     others<-setdiff(ctNames, ctName);
#     x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
#     specGenes[[ctName]]<-x;
#   }
#   specGenes;
# }
# 
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


###quick test
# library("pasilla")
# datafile = system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
# datafile
# 
# ## Read in the data making the row names the first column
# counttable <- read.table(datafile, header=T, row.names=1)
# head(counttable) 
# meta <- data.frame(
#   row.names=colnames(counttable),
#   condition=c("untreated", "untreated", "untreated", "untreated", "treated", "treated", "treated"),
#   libType=c("single", "single", "paired", "paired", "single", "paired", "paired"))
# meta$conditions <- relevel(meta$condition, ref="untreated")
# meta 
# 
# res <- doDESeq(meta, counttable, c("untreated", "treated"))
# 
# d <- newCountDataSet(counttable, meta$condition)
# ## Estimate library size and dispersion
# d <- estimateSizeFactors(d)
# d <- estimateDispersions(d) 
# plotDispEsts(d, main="DESeq: Per-gene dispersion estimates")
# results <- nbinomTest(d, "untreated", "treated")
##############


######### MDS plot using edgeR ##########3
doEdgeR_MDS <- function(stDat,
                        counts,
                        conditions,
                        colors,
                        symbols,
                        filename){
  cat("counts dimension: ", dim(counts),"\n")
  if(length(conditions) == 2){ #include a subset of samples
    stDat <- stDat[which(stDat$conditions == conditions[1] | stDat$conditions == conditions[2]),]
  }
  
  aux.counts <- counts
  aux.counts <- aux.counts[,rownames(stDat)]
  aux.counts <- aux.counts[rowSums(aux.counts) > 0, ]
  aux.counts <- aux.counts[,colSums(aux.counts) > 1e4]
  stDat <- stDat[colnames(aux.counts), ]
  
  cat("aux.counts dimension: ", dim(aux.counts),"\n")
  cat("samples dimension: ", dim(stDat),"\n")
  
  d <- DGEList(counts=aux.counts, group=stDat$conditions)
  d <- calcNormFactors(d)
  
  pdf(file=filename, width=5, height=5)
  plotMDS(d, labels=NULL, pch=symbols[stDat$conditions], col=colors[as.character(stDat$conditions)],
               xlab="", ylab="")
  dev.off()
}

doEdgeR <- function(stDat,
                    counts,
                    conditions,
                    file){
  
  cat("counts dimension: ", dim(counts),"\n")
  if(length(conditions) == 2){ #include a subset of samples
    stDat <- stDat[which(stDat$conditions == conditions[1] | stDat$conditions == conditions[2]),]
  }
  
  aux.counts <- counts
  aux.counts <- aux.counts[,rownames(stDat)]
  aux.counts <- aux.counts[rowSums(aux.counts) > 0, ]
  aux.counts <- aux.counts[,colSums(aux.counts) > 1e4]
  stDat <- stDat[colnames(aux.counts), ]
  
  cat("aux.counts dimension: ", dim(aux.counts),"\n")
  cat("samples dimension: ", dim(stDat),"\n")
  
  d <- DGEList(counts=aux.counts, group=stDat$conditions)
  d <- calcNormFactors(d)
  d <- estimateCommonDisp(d)
  d <- estimateTagwiseDisp(d)
  
  de <- exactTest(d, pair=conditions)
  #de <- exactTest(d, pair=c("306N", "306G"))
  results <- topTags(de, n=nrow(d))
  
  markers <- c('Cdkn2a', 'Casp8', 'Cdkn1a', 'Serpine1', 'Igfbp5', 'Igfbp2', 'Mmp3', 'Mmp13')
  markersFC <- data.frame(results[markers, c('logFC', 'PValue')])
  markersFC$genes <- rownames(markersFC)
  markersFC$genes <- factor(markersFC$genes, levels=markersFC$genes[order(markersFC$logFC)])
  
  
  #file <- "results/edgeR_marker_genes_306G_vs_306N"
  #file <- "results/edgeR_marker_genes_310G_vs_310N"
  file.pdf <- paste(file, 'pdf',sep=".")
  file.csv <- paste(file, 'csv',sep=".")
  full.csv <- paste(file, '_full_results.csv',sep="")
  pdf(file=file.pdf, width=2.5, height=2.5);
  ggplot(data=markersFC, aes(x=genes, y=logFC))+ 
    p <- geom_bar(stat="identity", width=.9) +
    theme_bw() + theme(legend.position=c(0.7, 0.7), legend.text=element_text(size=7)) +
    xlab("") + ylab("Log2 Fold Change") + 
    theme(axis.text.x = element_text(size=8, angle=30, hjust=1), axis.title.y = element_text(size = rel(0.7), angle = 90))
    print(p)
  dev.off();
  cat(file.csv, "\n", full.csv, "\n")
  write.csv(markersFC, file=file.csv)
  write.csv(results, file=full.csv)
}


########### Plotting marker distributions across all cells ########
plotDistributions <- function(stDat,
                              expDat, #genes by cells expression matrix (log2(FPKM), TPM, whatever)
                              condition, #annotations specifying a subset of samples
                              geneList, #list for each expression distributions will be plotted
                              width=6,
                              height=3.5,
                              num_columns=2,
                              filename){ 
  
  if(condition != "all"){
    stDat <- stDat[which(stDat$conditions == condition), ]
  }
  expDat <- expDat[,rownames(stDat)]
  expDat <- t(expDat[geneList,])
  
  plots <- list();
  
  for(i in seq(1:length(geneList))){
    df <- data.frame(gene=expDat[,geneList[i]])
    #p <- ggplot(df, aes(x=gene)) + geom_histogram(
    p <- ggplot(df, aes(x=gene)) + geom_histogram(aes(y=..density..),
      binwidth=1, colour="black", fill="green") + 
      #binwidth=.5, colour="black", fill="green") + 
      geom_density(alpha=.2, fill="#FF6666") +
      xlab("log2[normalized counts]") + ylab("Density") + 
      ggtitle(as.character(geneList[i])) +
      theme_bw() + theme(axis.title=element_text(size=10)) +
      theme(axis.text.x = element_text(size=8, angle=0, hjust=1),
            axis.text.y = element_text(size=8, hjust=1)) #+
      #geom_vline(aes(xintercept=mean(gene, na.rm=T)),   # Ignore NA values for mean
      #           color="red", linetype="dashed", size=1)
      
    plots[[i]] <- p
  }
  
  #filename <- paste("results/markers_distributions_log2FPKM_", condition, ".pdf", sep="")
  #filename <- paste("results/markers_distributions_DESeq_", condition, ".pdf", sep="")
  pdf(file=filename, width=width, height=height)
  multiplot(plotlist = plots, cols=num_columns)
  dev.off();
  
  return (plots)
}

plotDistributionsOverlaid <- function(stDat,
                              expDat, #genes by cells expression matrix (log2(FPKM), TPM, whatever)
                              condition, #annotations specifying a subset of samples
                              condition2,
                              geneList, #list for each expression distributions will be plotted
                              filename){ 
  
  if(condition != "all"){
    stDat1 <- stDat[which(stDat$conditions == condition), ]
    stDat2 <- stDat[which(stDat$conditions == condition2), ]
  }
  expDat1 <- expDat[,rownames(stDat1)]
  expDat1 <- t(expDat1[geneList,])
  
  expDat2 <- expDat[, rownames(stDat2)]
  expDat2 <- t(expDat2[geneList,])
  
  plots <- list();
  
  for(i in seq(1:length(geneList))){
    df1 <- data.frame(gene=expDat1[,geneList[i]])
    df2 <- data.frame(gene=expDat2[,geneList[i]])
    df1$condition <- condition
    df2$condition <- condition2
    df <- rbind(df1, df2)
    #p <- ggplot(df, aes(x=gene)) + geom_histogram(
    p <- ggplot(df, aes(x=gene, fill=condition)) + geom_histogram(alpha=0.7, show_guide=FALSE, position='dodge', aes(y=..density..),
                                                  binwidth=1, colour="black")  +
      #binwidth=.5, colour="black", fill="green") + 
      geom_density(alpha=.2, fill="#FF6666") +
      xlab("log2[norm counts]") + ylab("density") + 
      ggtitle(as.character(geneList[i])) +
      theme_bw() + theme(axis.title=element_text(size=10)) +
      theme(axis.text.x = element_text(size=8, angle=0, hjust=1),
            axis.text.y = element_text(size=8, hjust=1)) +
      #geom_vline(aes(xintercept=mean(gene, na.rm=T)),   # Ignore NA values for mean
      #           color="red", linetype="dashed", size=1) +
      scale_fill_manual("", values=c('green', 'red'))
      
    plots[[i]] <- p
  }
  
  #filename <- paste("results/markers_distributions_log2FPKM_", condition, ".pdf", sep="")
  #filename <- paste("results/markers_distributions_DESeq_", condition, ".pdf", sep="")
  pdf(file=filename, width=30, height=30)
  multiplot(plotlist = plots, cols=10)
  dev.off();
  
  return (plots)
}

#####
find_genes_go<-function# find transcript factors
(species='Hs', # species abbreviation
 process,
 category
){
  
  cat("Loading gene annotations ...\n")
  require(GO.db);
  
  if(species=='Hs'){
    require(org.Hs.eg.db);
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }
  else{
    require(org.Mm.eg.db);
    egSymbols<-as.list(org.Mm.egSYMBOL);
    goegs<-as.list(org.Mm.egGO2ALLEGS);
  }
  
  goterms<-as.list(GOTERM);
  goids<-names(goegs);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts==category];
  #bps<-onts[onts=='BP'];
  #bps<-onts[onts=='CC'];
  #bps<-onts[onts=='MF'];
  goids<-names(unlist(bps));
  
  cat("matching gene symbols and annotations")
  gobpList<-list();
  for(goid in goids){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gobpList[[goterm]]<-genes;
  }
  #save(gobpList, file="GO_MF_database_September_14_2015.R")
  regNames<-names(gobpList)[grep(process, names(gobpList))];
  genes<- unique(unlist(gobpList[regNames]));
  #save(genes, file='cell_adhesion_GO0007155_August_18_2015.R')
  
  return (list(regNames=regNames, genes=genes)) 
}



######## Obtain genes under specific GO terms #####
find_genes_go2<-function# find transcript factors
(species='Hs' # species abbreviation
){
  
  cat("Loading gene annotations ...\n")
  require(GO.db);
  
  if(species=='Hs'){
    require(org.Hs.eg.db);
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }
  else{
    require(org.Mm.eg.db);
    egSymbols<-as.list(org.Mm.egSYMBOL);
    goegs<-as.list(org.Mm.egGO2ALLEGS);
  }
  
  goterms<-as.list(GOTERM);
  goids<-names(goegs);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts=='BP'];
  #bps<-onts[onts=='CC'];
  goids<-names(unlist(bps));
  
  cat("matching gene symbols and annotations")
  gobpList<-list();
  for(goid in goids){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gobpList[[goterm]]<-genes;
  }
  regNames<-names(gobpList)[grep("miRNA", names(gobpList))];
  #regNames<-names(gobpList)[grep("^cell cycle$", names(gobpList))];
  
  #regNames<-names(gobpList)[grep("immune response", names(gobpList))];
  #regNames<-names(gobpList)[grep("glucocorticoid receptor signaling pathway", names(gobpList))];
  #regNames<-names(gobpList)[grep("hypoxia", names(gobpList))];
  #regNames<-names(gobpList)[grep("transforming growth factor beta receptor", names(gobpList))];
  #regNames<-names(gobpList)[grep("cell migration", names(gobpList))];
  #regNames<-names(gobpList)[grep("^cell adhesion", names(gobpList))];
  #genes<- unique(unlist(gobpList[regNames[1]]));
  genes<- unique(unlist(gobpList[regNames]));
  save(genes, file='cell_cycle_genes_December_30_2015.R')
  #save(genes, file='immune_related_genes_December_21_2015.R')
  #write.csv(regNames, file='immune_related_processes_December_21_2015.csv')
  #save(genes, file='cell_adhesion_GO0007155_August_18_2015.R')
  #genes<- unique(unlist(gobpList[regNames[30]]));
  #save(genes, file='cell_migration_GO0030334_August_18_2015.R')
  #save(genes, file='cell_migration_GO0030334_August_04_2015.R')
  #save(genes, file='results_August/Tgfb_related_pathways_August_04_2015.R')
  #save(genes, file='results_August/hypoxia_related_pathways_August_04_2015.R')
  #save(genes, file='glucocorticoid_receptor_signaling pathway_July172015.R')
  #regNames<-names(gobpList)[grep("glucocorticoid receptor signaling pathway", names(gobpList))];
  #regNames<-names(gobpList)[grep("regulation of transcription", names(gobpList))];
  #regNames<-names(gobpList)[grep("nuclear receptor", names(gobpList))];
  #genes<- unique(unlist(gobpList[regNames]));
  #regNames<-names(gobpList)[grep("cell surface", names(gobpList))];
  #regNames<-names(gobpList)[grep("extracellular region", names(gobpList))];
  #genes<- unique(unlist(gobpList[regNames[[1]]]));
  #save(genes, file='cell_surface_GO0009986.R')
  #save(genes, file='extracellular region_GO0005576_human.R')
  #save(genes, file='PPAR_signaling_pathway_mouse.R')
  return (genes) 
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

###### SPIA #####
spia2 <- function (de = NULL, all = NULL, organism = "hsa", data.dir = NULL, 
          pathids = NULL, nB = 2000, plots = FALSE, verbose = TRUE, 
          beta = NULL, combine = "fisher") 
{
  if (is.null(de) | is.null(all)) {
    stop("de and all arguments can not be NULL!")
  }
  rel <- c("activation", "compound", "binding/association", 
           "expression", "inhibition", "activation_phosphorylation", 
           "phosphorylation", "inhibition_phosphorylation", "inhibition_dephosphorylation", 
           "dissociation", "dephosphorylation", "activation_dephosphorylation", 
           "state change", "activation_indirect effect", "inhibition_ubiquination", 
           "ubiquination", "expression_indirect effect", "inhibition_indirect effect", 
           "repression", "dissociation_phosphorylation", "indirect effect_phosphorylation", 
           "activation_binding/association", "indirect effect", 
           "activation_compound", "activation_ubiquination")
  if (is.null(beta)) {
    beta = c(1, 0, 0, 1, -1, 1, 0, -1, -1, 0, 0, 1, 0, 1, 
             -1, 0, 1, -1, -1, 0, 0, 1, 0, 1, 1)
    names(beta) <- rel
  }
  else {
    if (!all(names(beta) %in% rel) | length(names(beta)) != 
          length(rel)) {
      stop(paste("beta must be a numeric vector of length", 
                 length(rel), "with the following names:", "\n", 
                 paste(rel, collapse = ",")))
    }
  }
  .myDataEnv <- new.env(parent = emptyenv())
  datload <- paste(organism, "SPIA", sep = "")
  if (is.null(data.dir)) {
    if (!paste(datload, ".RData", sep = "") %in% dir(system.file("extdata", 
                                                                 package = "SPIA"))) {
      cat("The KEGG pathway data for your organism is not present in the extdata folder of the SPIA package!!!")
      cat("\n")
      cat("Please generate one first using makeSPIAdata and specify its location using data.dir argument or copy it in the extdata folder of the SPIA package!")
    }
    else {
      load(file = paste(system.file("extdata", package = "SPIA"), 
                        paste("/", organism, "SPIA", sep = ""), ".RData", 
                        sep = ""), envir = .myDataEnv)
    }
  }
  if (!is.null(data.dir)) {
    if (!paste(datload, ".RData", sep = "") %in% dir(data.dir)) {
      cat(paste(data.dir, " does not contin a file called ", 
                paste(datload, ".RData", sep = "")))
    }
    else {
      load(file = paste(data.dir, paste(datload, ".RData", 
                                        sep = ""), sep = ""), envir = .myDataEnv)
    }
  }
  datpT = .myDataEnv[["path.info"]]
  if (!is.null(pathids)) {
    if (all(pathids %in% names(datpT))) {
      datpT = datpT[pathids]
    }
    else {
      stop(paste("pathids must be a subset of these pathway ids: ", 
                 paste(names(datpT), collapse = " "), sep = " "))
    }
  }
  datp <- list()
  path.names <- NULL
  hasR <- NULL
  for (jj in 1:length(datpT)) {
    sizem <- dim(datpT[[jj]]$activation)[1]
    s <- 0
    con <- 0
    for (bb in 1:length(rel)) {
      con = con + datpT[[jj]][[rel[bb]]] * abs(sign(beta[rel[bb]]))
      s = s + datpT[[jj]][[rel[bb]]] * beta[rel[bb]]
    }
    z = matrix(rep(apply(con, 2, sum), dim(con)[1]), dim(con)[1], 
               dim(con)[1], byrow = TRUE)
    z[z == 0] <- 1
    datp[[jj]] <- s/z
    path.names <- c(path.names, datpT[[jj]]$title)
    hasR <- c(hasR, datpT[[jj]]$NumberOfReactions >= 1)
  }
  names(datp) <- names(datpT)
  names(path.names) <- names(datpT)
  tor <- lapply(datp, function(d) {
    sum(abs(d))
  }) == 0 | hasR | is.na(path.names)
  datp <- datp[!tor]
  path.names <- path.names[!tor]
  IDsNotP <- names(de)[!names(de) %in% all]
  if (length(IDsNotP)/length(de) > 0.01) {
    stop("More than 1% of your de genes have IDs are not present in the reference array!. Are you sure you use the right reference array?")
  }
  if (!length(IDsNotP) == 0) {
    cat("The following IDs are missing from all vector...:\n")
    cat(paste(IDsNotP, collapse = ","))
    cat("\nThey were added to your universe...")
    all <- c(all, IDsNotP)
  }
  if (length(intersect(names(de), all)) != length(de)) {
    stop("de must be a vector of log2 fold changes. The names of de should be included in the refference array!")
  }
  ph <- pb <- pcomb <- nGP <- pSize <- smPFS <- tA <- tAraw <- KEGGLINK <- NULL
  set.seed(1)
  if (plots) {
    pdf("SPIAPerturbationPlots.pdf")
  }
  for (i in 1:length(names(datp))) {
    path <- names(datp)[i]
    M <- datp[[path]]
    diag(M) <- diag(M) - 1
    X <- de[rownames(M)]
    noMy <- sum(!is.na(X))
    nGP[i] <- noMy
    okg <- intersect(rownames(M), all)
    ok <- rownames(M) %in% all
    pSize[i] <- length(okg)
    if ((noMy) > 0 & (abs(det(M)) > 1e-07)) {
      gnns <- paste(names(X)[!is.na(X)], collapse = "+")
      KEGGLINK[i] <- paste("http://www.genome.jp/dbget-bin/show_pathway?", 
                           organism, names(datp)[i], "+", gnns, sep = "")
      X[is.na(X)] <- 0
      pfs <- solve(M, -X)
      smPFS[i] <- sum(pfs - X)
      tAraw[i] <- smPFS[i]
      if (plots) {
        par(mfrow = c(1, 2))
        plot(X, pfs - X, main = paste("pathway ID=", 
                                      names(datp)[i], sep = ""), xlab = "Log2 FC", 
             ylab = "Perturbation accumulation (Acc)", cex.main = 0.8, 
             cex.lab = 1.2)
        abline(h = 0, lwd = 2, col = "darkgrey")
        abline(v = 0, lwd = 2, col = "darkgrey")
        points(X[abs(X) > 0 & X == pfs], pfs[abs(X) > 
                                               0 & X == pfs] - X[abs(X) > 0 & X == pfs], col = "blue", 
               pch = 19, cex = 1.4)
        points(X[abs(X) > 0 & X != pfs], pfs[abs(X) > 
                                               0 & X != pfs] - X[abs(X) > 0 & X != pfs], col = "red", 
               pch = 19, cex = 1.4)
        points(X[abs(X) == 0 & X == pfs], pfs[abs(X) == 
                                                0 & X == pfs] - X[abs(X) == 0 & X == pfs], 
               col = "black", pch = 19, cex = 1.4)
        points(X[abs(X) == 0 & X != pfs], pfs[abs(X) == 
                                                0 & X != pfs] - X[abs(X) == 0 & X != pfs], 
               col = "green", pch = 19, cex = 1.4)
      }
      ph[i] <- phyper(q = noMy - 1, m = pSize[i], n = length(all) - 
                        pSize[i], k = length(de), lower.tail = FALSE)
      pfstmp <- NULL
      for (k in 1:nB) {
        x <- rep(0, length(X))
        names(x) <- rownames(M)
        x[ok][sample(1:sum(ok), noMy)] <- as.vector(sample(de, 
                                                           noMy))
        tt <- solve(M, -x)
        pfstmp <- c(pfstmp, sum(tt - x))
      }
      mnn <- median(pfstmp)
      pfstmp <- pfstmp - mnn
      ob <- smPFS[i] - mnn
      tA[i] <- ob
      if (ob > 0) {
        pb[i] <- sum(pfstmp >= ob)/length(pfstmp) * 2
        if (pb[i] <= 0) {
          pb[i] <- 1/nB/100
        }
        if (pb[i] > 1) {
          pb[i] <- 1
        }
      }
      if (ob < 0) {
        pb[i] <- sum(pfstmp <= ob)/length(pfstmp) * 2
        if (pb[i] <= 0) {
          pb[i] <- 1/nB/100
        }
        if (pb[i] > 1) {
          pb[i] <- 1
        }
      }
      if (ob == 0) {
        if (all(pfstmp == 0)) {
          pb[i] <- NA
        }
        else {
          pb[i] <- 1
        }
      }
      if (plots && sd(pfstmp)!=0) {
        bwidth = sd(pfstmp)/4
        if (bwidth > 0) {
          plot(density(pfstmp, bw = bwidth), cex.lab = 1.2, 
               col = "black", lwd = 2, main = paste("pathway ID=", 
                                                    names(datp)[i], "  P PERT=", round(pb[i], 
                                                                                       5), sep = ""), xlim = c(min(c(tA[i] - 
                                                                                                                       0.5, pfstmp)), max(c(tA[i] + 0.5, pfstmp))), 
               cex.main = 0.8, xlab = "Total Perturbation Accumulation (TA)")
        }
        else {
          pfsTab = table(pfstmp)
          plot(as.numeric(names(pfsTab)), as.numeric(pfsTab), 
               cex.lab = 1.2, col = "black", main = paste("pathway ID=", 
                                                          names(datp)[i], "  P PERT=", round(pb[i], 
                                                                                             5), sep = ""), xlim = c(min(c(tA[i] - 
                                                                                                                             0.5, pfstmp)), max(c(tA[i] + 0.5, pfstmp))), 
               cex.main = 0.8, xlab = "Total Perturbation Accumulation (TA)", 
               ylab = "frequency")
        }
        abline(v = 0, col = "grey", lwd = 2)
        abline(v = tA[i], col = "red", lwd = 3)
      }
      pcomb[i] <- combfunc(pb[i], ph[i], combine)
    }
    else {
      pb[i] <- ph[i] <- smPFS[i] <- pcomb[i] <- tAraw[i] <- tA[i] <- KEGGLINK[i] <- NA
    }
    if (verbose) {
      cat("\n")
      cat(paste("Done pathway ", i, " : ", substr(path.names[names(datp)[i]], 
                                                  1, 30), "..", sep = ""))
    }
  }
  if (plots) {
    par(mfrow = c(1, 1))
    dev.off()
  }
  pcombFDR = p.adjust(pcomb, "fdr")
  phFdr = p.adjust(ph, "fdr")
  pcombfwer = p.adjust(pcomb, "bonferroni")
  Name = path.names[names(datp)]
  Status = ifelse(tA > 0, "Activated", "Inhibited")
  res <- data.frame(Name, ID = names(datp), pSize, NDE = nGP, 
                    pNDE = ph, tA, pPERT = pb, pG = pcomb, pGFdr = pcombFDR, 
                    pGFWER = pcombfwer, Status, KEGGLINK, stringsAsFactors = FALSE)
  res <- res[!is.na(res$pNDE), ]
  res <- res[order(res$pG), ]
  rownames(res) <- NULL
  res
}


WGCNA_findBeta<-function(expDat){
  # estimate scale-free network fit for a range of Beta
  #
  # Args:
  #   expDat:
  #      matrix of gene expression data
  #
  # Returns:
  #      list of $powerEstimate (automaitcally selected beta) and $fitIndices (table of all values)
  #      prints out Rsqr fits
  #        
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  pickSoftThreshold(t(expDat), powerVector = powers, verbose = 5);
}

WGCNA_findModules<-function(expDat, beta, minModuleSize=15, deepSplit=2){
  # find modules
  #
  # Args:
  #   expDat:matrix of gene expression data
  #   beta: beta
  #   minModuleSize: 
  #   deepSplit:
  #
  # Returns:
  #      list of genes in each module  
  cat("WGCNA-ing ...");
  results <- list()
  wgcnaRes<-blockwiseModules(t(expDat),
                             maxBlockSize=9000,
                             power=beta,
                             minModuleSize=minModuleSize,
                             deepSplit=deepSplit,
                             pamRespectsDendro=FALSE,
                             mergeCutHeight=0.10,
                             numericLabels=TRUE,
                             saveTOMs=F,
                             verbose=5,
                             networkType='signed',
                             TOMType='signed');
  cat(" done\n");
  results[['result']] <- wgcnaRes
  results[['modules']] <- .convertMods(wgcnaRes, rownames(expDat));
  
  results
}

.convertMods<-function(wgcnaRes, geneVect){
  ans<-list();
  uMods<-unique(wgcnaRes$colors);
  modNames<-colnames(wgcnaRes$MEs);
  for(modName in modNames){
    suff<-strsplit(modName,"ME")[[1]][2];
    ans[[modName]]<-geneVect[which(wgcnaRes$colors==suff)];
  }
  ans;
}

find_tfs<-function# find transcript factors
(species='Hs' # species abbreviation
){
  
  cat("Loading gene annotations ...\n")
  require(GO.db);
  
  if(species=='Hs'){
    require(org.Hs.eg.db);
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }
  else{
    require(org.Mm.eg.db);
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
  cat("Regulation of transcription: ", length(trs),"\n");
  cat(regNames, '\n')
  
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
  cat("DNA binding: ", length(dbs),"\n");
  #cat("RNA binding: ", length(dbs),"\n");
  sort(intersect(trs, dbs));
  #sort(dbs)
}
