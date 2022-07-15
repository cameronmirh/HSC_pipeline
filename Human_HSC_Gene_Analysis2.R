###############################################################################################
# This script is to identify metabolic genes that correlate with different haematopoietic
# trajectories using microarray data https://www.ncbi.nlm.nih.gov/pubmed/?term=21241896.
#
# @author: Cameron Mirhossaini (cameronmirhossaini@gmail.com) ; Jason Cosgrove (jason.cosgrove@curie.fr)
# @date:   28/08/2019
###############################################################################################

#load in the required libraries
library(affy)
library(simpleaffy)
library(affyPLM)
library(limma)
library(dplyr)
library(GEOquery)
library(container)
library("hgu133a.db")
library("annotate")
library("hgu95av2.db") 
library(gplots)
library(ggplot2)
library(xgboost)
library(pROC)
library(magrittr)
library(pheatmap)
library(dichromat)


#-------------------------------------------
# Global Variables
#-------------------------------------------

#sample data used for analysis
load(file = "/Users/cameronmirhossaini/Documents/2018-19/4 Summer/HumanMetabolism/dmapExpression.rdata")

#characteristics of each sample
load(file = "/Users/cameronmirhossaini/Documents/2018-19/4 Summer/HumanMetabolism/samples.rdata")


#compilation of genes important in metabolism from KEGG database
load(file = "/Users/cameronmirhossaini/Documents/2018-19/4 Summer/HumanMetabolism/keggGenes.rdata")

#dictionary that maps probe ID to genes
load(file="/Users/cameronmirhossaini/Documents/2018-19/4 Summer/MLScript/probeToGene.Rdata")

#these dictionaries are global variables that will be used throughout the script
linToType <-Dict$new()
linToTissue <- Dict$new()
linToMat <- Dict$new()

typeToDescr <-Dict$new()
sampleToLin <-Dict$new()

#characterizations of cell types into lineages
lymphoid <- c("B Cell Lineage", "T Cell Lineage" , "Innate Lymphocyte Lineage", "NK Cell Lineage")
myeloid <- c("Dendritic Cell Lineage" , "Neutrophil Lineage", "Mast Cell Lineage",  "Macrophage Lineage",             
             "Basophil Lineage",  "Eosinophil Lineage")


#*******************************************
#-------------------------------------------
# Helper Functions
#-------------------------------------------
#*******************************************

#only run once per script
#loads all dictionaries with data from samples in order to optimize functions
#  used to change data visualization in PCA 
loadDictionairies <- function() {
  mpp <- 0
  rpp <- 0
  lym <- 0 
  mye <- 0
  ery <- 0
  meg <- 0
  
  for (i in 1:nrow(samples)) {
    lin <- toString(samples[i,5])
    type <- ""
    
    #convert each cell type into more general cell lineage
    if(lin %in% lymphoid) {
      lin <- paste("Lymphoid.",lym, sep="")
      type <- paste(samples[i,2], ".", lym, sep="")
      lym <- 1 + lym
    }
    
    if(lin %in% myeloid) {
      lin <- paste("Myeloid.", mye,sep ="")
      type <- paste(samples[i,2], ".", mye, sep="")
      mye <- 1 + mye
    }
    
    if(lin ==  "Erythrocyte Lineage") {
      lin <- paste("Erythroid.", ery,sep ="")
      type <- paste(samples[i,2], ".", ery, sep="")
      ery <- 1 + ery
    }
    
    if(lin ==  "Megakaryocyte Lineage") {
      lin <- paste("Meg.", meg, sep="")
      type <- paste(samples[i,2], ".", meg, sep="")
      meg <- 1 + meg
    }
    
    if(lin ==  "Multi Potential Progenitor") {
      lin <- paste("MPP.", mpp, sep = "")
      type <- paste(samples[i,2], ".", mpp, sep="")
      mpp <- mpp + 1
    }
    
    if(lin ==  "Restricted Potential Progenitor")  {
      lin <- paste("RPP.", rpp,  sep = "")
      type <- paste(samples[i,2], ".", rpp, sep="")
      rpp <- rpp + 1
    }
    
    sampleToLin$add(toString(samples[i,1]), lin)
    linToType$add(lin, type)
    typeToDescr$add(type, paste(samples[i, 4]))
    linToMat$add(lin, paste(samples[i,6]))
    linToTissue$add(lin, paste(samples[i,3]))
  }
}


#Returns matrix with filtered cell lineages/cell types/samples using drops matrix
#rows = cell lineage
#column = genes
#filtering does not include RPP/MPP filtering and must occur when columns are cell lineages
formatData <- function(isMetabolic, drops) {
  dmap <- dmapExpression
  
  
  #convert column names from sample ID to cell lineages
  for(i in 1:ncol(dmap)){
    colnames(dmap)[i] <- sampleToLin$peek(colnames(dmap)[i])
  }

  #convert row names from probeID to gene symbols
  for (i in 1:nrow(dmap)) {
    probeId <- rownames(dmap)[i]
    geneSym <- probeToGene$peek(rownames(dmap)[i])
    #if probeID has valid hash pairing, proceed
    if(!is.null(geneSym) && !is.na(geneSym)) {
      #if gene is duplicate, then concat probeId to it
      if(geneSym %in% rownames(dmap)) {
        rownames(dmap)[i] <- paste(geneSym, "_", probeId, sep = "")
      } else {
        rownames(dmap)[i] <- geneSym
      }
    } 
  }
  
  #find which specific columns to drop
  cols <- c()
  if(length(drops) != 0) {
    for(i in 1:length(drops)) {
      if (substr(drops[i], nchar(drops[i]), nchar(drops[i])) == ".") {
        cols <- c(cols, grep(drops[i], colnames(dmap), value=T))
      } else {
        cols <- c(cols, drops[i])
      }
    }
  }
  
  #filter out non-desired samples (columns)
  if(length(cols) != 0) {
    dmap <- dmap[,!(colnames(dmap) %in% cols) ]
  }
  
  #find intersect of genes from our sample and all metabolic genes
  interestingGenes <- intersect(keggGenes, rownames(dmap))
  
  #filter on metabolic genes
  if(isMetabolic) dmap <- dmap[interestingGenes,]
  
  return(t(data.frame(dmap)))
  
  
}

#for training on Mature Cells
#returns matrix of samples predicted incorrectly
#prints out Confusion Matrix & AUC
#***Make parameters for xgb
trainGradientBoost <-function(matrixWID, matrixNoID, nfolds, numToLin) {
  
  labels <- as.factor(rownames(matrixNoID))
  weights <- getWeights(labels, matrixNoID, 2)
  
  
  nfolds <<- 10   # CV folds
  nrounds = 16  # boosting iterationsm; was 20
  # convert the data matrix to an xgb.matrix object
  xgb_data = xgb.DMatrix(as.matrix(matrixNoID), label=as.numeric(labels)-1)
  #must run matrix without number tags
  mdXGB <- xgb.cv(data = xgb_data, label = as.numeric(labels)-1, nfold = nfolds,
              params = list(subsample=0.8, colsample_bytree=0.5,colsample_bynode = 0.5,weight = weights,   lambda = 1),
              nrounds = nrounds,
              prediction = T, verbose = T,
              metrics = 'merror', objective = "multi:softprob", num_class = length(unique(labels)),
              callbacks = list(cb.cv.predict(save_models = TRUE)))
  
  #print Confusion Matrix
  xgb.preds = apply(mdXGB$pred, 1, which.max) 
  print(table(xgb.preds, labels))
  
  #create a dictionary that links lineage to its numeric entry
  predMat <- table(xgb.preds, labels)
  
  for(i in 1:ncol(predMat)) {
    numToLin$add(rownames(predMat)[i],colnames(predMat)[i])
  }
  
  #print AUC
  roc_obj<- multiclass.roc(labels, xgb.preds)
  print(auc(roc_obj))
  
 return(mdXGB)
  
}

#uses global mdXGB to predict the cell lineages of progenitors
predictProgenitors <- function(nfolds, numToLin, mdXGB) {
  progenitorMatrix <- formatData(T, c("Myeloid.","Lymphoid.", "Erythroid.", "Meg."))
  
  # for each cell sample, lets generate a lineage prediction from each of our 10 models
  # and store it in a matrix
  # we have 10 models because we did 10-fold cross validation
  output.mat <- matrix(0, nrow = length(rownames(progenitorMatrix)), ncol = nfolds)
  output.mat.numeric <- matrix(0, nrow = length(rownames(progenitorMatrix)), ncol = nfolds)
  rownames(output.mat) <- rownames(progenitorMatrix)
  for(i in 1:length(rownames(progenitorMatrix))){
    for(j in 1:nfolds){
      index <- which.max(predict(mdXGB$models[[j]],t(as.matrix(progenitorMatrix[i,]))))
      output.mat[i,j] <- numToLin$peek(index)
      output.mat.numeric[i,j] <- index
    }
  }
  
  #change row names of output matrix from cell lineage to cell types
  for(i in 1:nrow(output.mat)){
    rownames(output.mat)[i] <- paste(linToType$peek(rownames(output.mat)[i]))
  }
  
  consensus.prediction <- data.frame(apply(output.mat, 1,Consensus))
  print(consensus.prediction)
  return(consensus.prediction)
  
}



#Identifies most important Genes used in Gradient Boost Algortithm 
# ***update would to add customization on how to filter: Via Gain, Frequency, or Cover
identifyImportantGenes <- function(data, mdXGB) {
  
  all_features = lapply(mdXGB$models, function(x) {xgb.importance(model = x, feature_names = colnames(data))}$Feature)
  all_features = unique(unlist(all_features))
  
  
  all_variable_importance = data.frame()
  for(fold in 1:nfolds){
    # xgb object has 1 model per fold of the CV. This loop combines variable importances for all folds
    fold_variable_importance = data.frame(Fold=rep(fold,length(all_features)), Feature = all_features, 
                                          Gain = rep(0,length(all_features)),
                                          Cover = rep(0,length(all_features)),
                                          Frequency = rep(0,length(all_features)))
    importance_matrix <- xgb.importance(model = mdXGB$models[[fold]], feature_names = colnames(data))
    for(i in 1:length(all_features)){
      row = importance_matrix[importance_matrix$Feature==all_features[i]]
      if(nrow(row) > 0){
        fold_variable_importance[i,] = list(fold, all_features[i], row$Gain, row$Cover, row$Frequency)
      }
    }
    all_variable_importance = rbind(all_variable_importance, fold_variable_importance) 
  }
  
  
  genes.ordered <-  all_variable_importance %>%  group_by(Feature) %>% 
    summarise(Gain = mean(Gain), Cover = mean(Cover), Frequency = mean(Frequency)) %>%
    dplyr::mutate(Feature = reorder(Feature, Gain)) %>% dplyr::filter(Gain > 0.0005)  
  
  #these genes should be used for k-means clustering and heat map
  genes.ordered <- genes.ordered[order(genes.ordered$Gain, decreasing = T),]
  
  
  #Graph Top Genes
  print(all_variable_importance %>%  group_by(Feature) %>% 
          summarise(Gain = mean(Gain), Cover = mean(Cover), Frequency = mean(Frequency)) %>%
          dplyr::mutate(Feature = reorder(Feature, Gain)) %>% dplyr::filter(Gain > 0.002)  %>% 
          ggplot(aes(x=Feature, y=Gain )) + 
          geom_bar(stat="identity") +
          coord_flip() +
          theme_classic()
  )
  print(paste(genes.ordered$Feature))
  return(paste(genes.ordered$Feature))
  
}

#Runs K-mean clustering
kmeansClustering <-function(genes, matrixNoID, clusters) {
  geneCluster <- kmeans(t(matrixNoID)[genes,],clusters)
  str(geneCluster)
  
  cluster1 <- c()
  cluster2 <- c()
  cluster3 <- c()
  cluster4 <- c()
  for(i in 1:length(geneCluster$cluster)) {
    if (geneCluster$cluster[i] == 1) cluster1 <- c(cluster1, toString(genes[i]))
    if (geneCluster$cluster[i] == 2) cluster2 <- c(cluster2, toString(genes[i]))
    if (geneCluster$cluster[i] == 3) cluster3 <- c(cluster3, toString(genes[i]))
    if (geneCluster$cluster[i] == 4) cluster4 <- c(cluster4, toString(genes[i]))
  }
  print(cluster1)
  print(cluster2)
  print(cluster3)
  print(cluster4)
  return(c(cluster1,cluster2,cluster3,cluster4))
}

#create heatmap using data
heatmapSimulator <- function(genes, matrixWID) {
  bk = c(seq(-2,2,by = .01))
  redblue<- colorRampPalette(c("dark blue","white","red")) (n = length(bk)-1 )
  
  
  
  pheatmap(t(matrixWID)[genes,], scale = "row", 
           breaks = bk,
           color = redblue,  
           fontsize_col = 3,
           fontsize_row = 4
  )
  
  
  
}

#generates as many Box Plots of distributions per lineage for each gene
#parameters is just a file name/location
generateBoxPlots <- function(genes, tmatrixNoID, fileNameBP) {
  pdf(paste("/Users/cameronmirhossaini/Documents/2018-19/4 Summer/MLScript/Box_Plots/", fileNameBP, ".pdf", sep=""))
  for(i in 1:length(genes)) {
    invisible(print(boxplot(tmatrixNoID[genes[i], grepl("Lymphoid",colnames(tmatrixNoID))],
                            tmatrixNoID[genes[i], grepl("Myeloid",colnames(tmatrixNoID))], 
                            tmatrixNoID[genes[i], grepl("Erythroid",colnames(tmatrixNoID))],
                            tmatrixNoID[genes[i], grepl("Meg",colnames(tmatrixNoID))],
                  names = c("Lymphoid","Myeloid","Erythroid", "Meg"),
                  main = paste(genes[i], "Gene Distribution"),
                  xlab = "Lineage"
          
                  
                  
    )))
  }
  
  dev.off()

}


#colorbyWhat: 1 == Lineage
#      : 2 == Tissue
#      : 3 == Cell Type
#      : 4 == maturity
# pcaMax: how many charts you want
generatePCPlots <- function(matrixWID, incorrectSamples, colorbyWhat, pcaMax, isText) {
  
  pca <- prcomp(matrixWID, scale=TRUE)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  plot_list = list()
  for (i in 2:pcaMax) {
    sp <-plotPCA(pca,i-1,i,incorrectSamples, isText, colorbyWhat, pca.var, pca.var.per)
    plot_list[[i-1]] <- sp
    print(plot_list[[i-1]])
  }
  #dev.off()
  
}

#takes two PC's and plots them
#must input how to sort samples and if you want text or not
plotPCA <- function(pca, pc1, pc2, incorrectSamples, isText,colorbyWhat, pca.var, pca.var.per) {
  
  pca.data <- data.frame(Sample=rownames(pca$x),
                         X=pca$x[,pc1],
                         Y=pca$x[,pc2])
  
  pca.dataCopy <- pca.data
  pca.data <- resortPCAData(pca.data, colorbyWhat)
  
  top.genes1 <- findTopGenes(pca, pc1, 3)
  top.genes2 <- findTopGenes(pca, pc2, 3)
  
  # Remove the decimal dot and the numbers behind it

  if(!isText) pca.data$Sample <- strtrim(pca.data$Sample, 3)
  
  sp <- ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
    xlab(paste("PC", pc1, " - ", pca.var.per[pc1], "% ", top.genes1[1], " ", top.genes1[2], " ", top.genes1[3], " ", sep="")) +
    ylab(paste("PC", pc2, " - ", pca.var.per[pc2], "% ", top.genes2[2], " ", top.genes1[2], " ", top.genes1[3], " ", sep="")) +
    theme_bw() +
    ggtitle(paste("PC", pc1, " vs. ", "PC", pc2, sep=""))
  
  #modify points on plot to either be sorted by text or color
  if(isText) {
    sp <- sp+geom_text(size=.5, check_overlap = TRUE)
  } else {
    sp <- sp + geom_point() + aes(colour = factor(Sample))
  } 
  
  return(sp +  geom_text(data = pca.dataCopy %>% dplyr::filter(Sample %in% incorrectSamples), colour="black", size=2))
  
  
}



#This function takes in the following parameters:
#isMetabolic -- boolean; T if you want to keep only metabolic genes throughout this analysis
#cellsToDrop -- vector; something like c("RPP.","MPP.", "Meg.", "Myeloid.6","Myeloid.7","Myeloid.8","Myeloid.9","Myeloid.10")
#           if you want to remove entire cell lineage, write with a period ex Lymphoid. removes all Lymphoid while Lymphoid.0
#           only removes Lymphoid.0
#           to identify which cells correspond to which number, use linToType$values()
#kclust -- integer; how many clusters you want kclusters to make. Should be as many lineages that you have
#bpFileName -- string; file that will store all boxplots of gene distibution for each cell type
#colorByWhat -- int; when doing PCA, decides how you want the points to be colored
#             #colorbyWhat: 1 == Lineage
#                         : 2 == Tissue
#                         : 3 == Cell Type
#                         : 4 == maturity
#nPCAs -- int; how many PCA charts do you want
#               recommended is 7
runAnalysis <-function(isMetabolic, cellsToDrop, kclust, bpFileName, colorByWhat, nPCAs) {
  
  print("Formatting data...")
  matrixWID <- formatData(isMetabolic, cellsToDrop)
  matrixNoID <<- t(removeCellID(t(matrixWID)))
  
  numToLin <<- Dict$new()
  nfolds <- 10
  
  print("Preparing Gradient Boost...")
  mdXGB <- trainGradientBoost(matrixWID, matrixNoID, nfolds, numToLin)
  
  print("Predicting Progenitors...")
  progenitorResults <- predictProgenitors(nfolds, numToLin, mdXGB)
  
  incorrectPred <- findIncorrectSorting(matrixWID, data.frame(mdXGB$pred))
  genes <- identifyImportantGenes(matrixNoID, mdXGB)
  print("printing kclusters...")
  ordered.clusters <- kmeansClustering(genes, matrixNoID, kclust)
  print("formatting heat map...")
  heatmapSimulator(genes, matrixWID)
  print("generating box plots.previous info will be lost in console.")
  cat ("Press [enter] to continue")
  line <- readline()
  generateBoxPlots(genes, t(matrixNoID), bpFileName)
  generatePCPlots(matrixWID, incorrectPred[,"trueLin"], colorByWhat, nPCAs, F)
  dev.off()
  
}

#-------------------------------------------
# Minor Helper Functions
#-------------------------------------------

#for each cell type find out what the most popular lineage classification is and store it in a vector
Consensus <- function(x){return(names(which.max(table(x))))}

## get the name of the top numGene measurements (genes) that contribute most to pcX.
findTopGenes <- function(pca, x, numGenes) {
  loading_scores <- pca$rotation[,x]
  gene_scores <- abs(loading_scores) ## get the magnitudes
  gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
  top.genes <- names(gene_score_ranked[1:numGenes])
  return(top.genes)
}

#filters PCA data by user's input, either Tissue, Cell Type, Maturity, or Lineage
resortPCAData <- function(pca.data, byWhat) {
  keyToValue <- Dict$new()
  if(byWhat == 1) return(pca.data)
  if(byWhat == 2) keyToValue <- linToTissue
  if(byWhat == 3) keyToValue <- linToType
  if(byWhat == 4) keyToValue <-linToMat
  
  pca.data$Sample <- paste(pca.data$Sample)
  for(i in 1:length(pca.data$Sample)) {
    pca.data$Sample[i] <- keyToValue$peek(pca.data$Sample[i])
  }
  
  return(pca.data)
}


#creates matrix of all true lineage, predicted lineage, and true cell type 
findIncorrectSorting <- function(matrixWID, results) {
  colNum <- ncol(results)
  results$lineage_predicted <- findPredicted(results)
  results$lineage_true <- rownames(matrixWID)
  results$celltype <- findCellType(results, colNum)
  return(condenseIncorrectSorting(results, colNum))
}

#condenses Incorrectly sorted matrix into legible and clean format
condenseIncorrectSorting <- function(results, colNum) {
  trueLin <- c()
  predLin <- c()
  cellType <- c()
  for (i in 1:nrow(results)) {
    if(results[i,colNum+3] != "") {
      predLin <- c(predLin, results[i,colNum+1])
      trueLin <- c(trueLin, results[i,colNum+2])
      cellType <- c(cellType, results[i,colNum+3])
    }
    
  }
  mat <- cbind(trueLin,predLin, cellType)
  print(mat)
  return(mat)
  
}

#find cell type of incorrectly predicted samples
findCellType <- function(results, colNum) {
  pred <- 1
  true <- 2
  cellTypes <- c()
  for (i in 1:nrow(results)) {
    if(substr(results[i,colNum+pred],1,3) != substr(results[i,colNum+true],1,3)) {
      cellTypes <- c(cellTypes, linToType$peek(results[i,colNum+true]))
    } else {
      cellTypes <- c(cellTypes,"")
    }
  }
  return(cellTypes)
  
}


#gets max of each row and converts column into lineage
findPredicted <- function(results) {
  predictions <- c()
  for (r in 1:nrow(results)) {
    predictions <- c(predictions, numToLin$peek(which.max(results[r,])))
  }
  return(predictions)
  
}


#variable sc is scale for Megs if you so choose to keep Megs in the analysis
getWeights <- function(labels, matrixNoID, sc) {
  
  #weights are calculated by dividing the number of cells from the minimum lineage 
  # and dividing it by the cell count from each lineage
  minCount <- min(table(rownames(matrixNoID)))
  eryWeight <- minCount / table(rownames(matrixNoID))["Erythroid"]
  lymWeight <- minCount / table(rownames(matrixNoID))["Lymphoid"]
  megWeight <- (minCount / table(rownames(matrixNoID))["Meg"]) * sc #scale bc very few meg samples
  myeWeight <- minCount / table(rownames(matrixNoID))["Myeloid"] 
  
  #we have unequal sample sizes for each lineage so we need to add a weight variable that accounts for this
  weights <- c()
  for(i in 1:length(labels)){
    if(labels[i] == "Lymphoid") weights <- c(weights, lymWeight)
    if(labels[i] == "Myeloid") weights <- c(weights, myeWeight)
    if(labels[i] == "Erythroid") weights <- c(weights, eryWeight)
    if(labels[i] == "Meg") weights <- c(weights, megWeight)
  }
  
  return(weights)
  
}


#only applicable when lineages are columns
#removes period and ID from lineages for XGBoost Algorithm and PCA
removeCellID <- function(mymatrix) {
  for (i in 1:ncol(mymatrix)) {
    if(substr(colnames(mymatrix)[i],1,3) == "Lym") colnames(mymatrix)[i] <- "Lymphoid" 
    if(substr(colnames(mymatrix)[i],1,3) == "Mye") colnames(mymatrix)[i] <- "Myeloid" 
    if(substr(colnames(mymatrix)[i],1,3) == "Ery") colnames(mymatrix)[i] <- "Erythroid" 
    if(substr(colnames(mymatrix)[i],1,3) == "Meg") colnames(mymatrix)[i] <- "Meg" 
  }
  return(mymatrix)
  
}



#*******************************************
#-------------------------------------------
# Main
#-------------------------------------------
#*******************************************

loadDictionairies()
runAnalysis(T,c("RPP.","MPP.", "Meg.", "Myeloid.6","Myeloid.7","Myeloid.8","Myeloid.9","Myeloid.10"), 3, "Random_File", 3, 7)
















