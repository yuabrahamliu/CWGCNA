#Clustering####

signadj <- function(line){
  
  if (sign(line[1]) == -1) {
    line <- line * -1
  }
  return(line)
  
}

baseCCA <- function(dat1, 
                    dat2, 
                    numcc = 30, 
                    seednum = 2022){
  
  if(is.null(dat2)){
    dat2 <- dat1
  }
  
  set.seed(seednum)
  
  numcc <- min(c(numcc, nrow(dat1), nrow(dat2)))
  
  sampleids <- intersect(colnames(dat1), colnames(dat2))
  
  dat1 <- dat1[, sampleids, drop = FALSE]
  dat2 <- dat2[, sampleids, drop = FALSE]
  sampleids1 <- paste0(sampleids, '_dat1')
  sampleids2 <- paste0(sampleids, '_dat2')
  
  gc()
  mat <- crossprod(x = t(dat1), y = t(dat2))
  
  gc()
  mat <- mat/(ncol(dat1) - 1)
  
  suppressMessages(
    
    cca.svd <- tryCatch({
      
      irlba::irlba(A = mat, nv = numcc)
      
    }, error = function(err){
      
      svd(x = mat, nu = numcc, nv = numcc)
      
    }, warning = function(war){
      
      svd(x = mat, nu = numcc, nv = numcc)
      
    })
    
  )
  
  gc()
  
  u <- apply(X = cca.svd$u, MARGIN = 2, FUN = signadj)
  v <- apply(X = cca.svd$v, MARGIN = 2, FUN = signadj)
  
  cca.samples1 <- crossprod(x = dat1, y = u)
  cca.samples2 <- crossprod(x = dat2, y = v)
  cca.samples <- rbind(cca.samples1, cca.samples2)
  
  colnames(cca.samples) <- paste0("CC", 1:numcc)
  rownames(cca.samples) <- c(sampleids1, sampleids2)
  
  res <- list()
  res$cca.samples <- cca.samples
  res$u <- u
  res$v <- v
  res$sigular <- cca.svd$d
  
  return(res)
  
}

calccasamples <- function(dat1, 
                          dat2, 
                          numcc = 30, 
                          seednum = 2022){
  
  #CCA
  cca.res <- baseCCA(dat1 = dat1, 
                     dat2 = dat2, 
                     numcc = numcc, 
                     seednum = seednum)
  
  #weighting by variance
  
  if(!is.null(dat2)){
    ucov <- cov(t(dat1))
    
    gc()
    uvars <- crossprod(cca.res$u, ucov)
    uvars <- crossprod(t(uvars), cca.res$u)
    uvars <- diag(uvars)
    
    vcov <- cov(t(dat2))
    
    gc()
    vvars <- crossprod(cca.res$v, vcov)
    vvars <- crossprod(t(vvars), cca.res$v)
    vvars <- diag(vvars)
    
    uweight <- min(sum(uvars)/nrow(ucov), 1)
    vweight <- min(sum(vvars)/nrow(vcov), 1)
    
    uweight <- uweight/(uweight + vweight)
    vweight <- 1 - uweight
    
    cca.samples1 <- cca.res$cca.samples[1:(0.5*nrow(cca.res$cca.samples)), , 
                                        drop = FALSE]
    cca.samples2 <- cca.res$cca.samples[(0.5*nrow(cca.res$cca.samples) + 1): 
                                          nrow(cca.res$cca.samples), , drop = FALSE]
    
    cca.samples <- uweight*cca.samples1 + vweight*cca.samples2
    row.names(cca.samples) <- gsub(pattern = '_dat1', replacement = '', 
                                   x = row.names(cca.samples))
  }else{
    
    cca.samples <- cca.res$cca.samples[1:(0.5*nrow(cca.res$cca.samples)), , 
                                       drop = FALSE]
  }
  
  row.names(cca.samples) <- gsub(pattern = '_dat1', replacement = '', 
                                 x = row.names(cca.samples))
  
  return(cca.samples)
  
}

singlecluster <- function(i = 1, 
                          dats, 
                          numcc = 20, 
                          seednum = 2022, 
                          k = c(2), 
                          consensuscluster = FALSE){
  
  if(consensuscluster == TRUE){
    
    for(j in 1:length(dats)){
      
      set.seed(seednum + i)
      dats[[j]] <- dats[[j]][sample(x = 1:nrow(dats[[j]]), 
                                    size = nrow(dats[[j]]), 
                                    replace = TRUE), , drop = FALSE]
      
    }
    
  }
  
  if(length(dats) == 1){
    
    cca.samples <- calccasamples(dat1 = dats[[1]], 
                                 dat2 = NULL, 
                                 numcc = numcc, 
                                 seednum = seednum)
    
  }else{
    
    cca.samples <- t(dats[[1]])
    
    for(j in 2:length(dats)){
      
      cca.samples <- calccasamples(dat1 = t(cca.samples), 
                                   dat2 = dats[[j]], 
                                   numcc = numcc, 
                                   seednum = seednum)
      
      
    }
    
  }
  
  k <- k[k <= nrow(cca.samples)]
  
  for(i in 1:length(k)){
    
    if(i == 1){
      preds <- list()
      ivis <- list()
      
    }
    
    #clustering
    pred <- kmeans(x = cca.samples, centers = k[i])
    
    preds[[i]] <- pred
    names(preds)[i] <- as.character(paste0('k = ', k[i]))
    
    res <- list(preds = preds)
    
    if(consensuscluster == FALSE){
      #IVI
      ivi <- unlist(clusterCrit::intCriteria(traj = cca.samples,
                                             part = pred$cluster,
                                             crit = c('Silhouette',
                                                      'Calinski_Harabasz',
                                                      'Davies_Bouldin')))
      
      
      
      
      ivis[[i]] <- ivi
      names(ivis)[i] <- as.character(paste0('k = ', k[i]))
      
      res <- list(preds = preds, ivis = ivis)
      
    }
    
  }
  
  if(consensuscluster == FALSE){
    
    res$cca.samples <- cca.samples
    
  }
  
  return(res)
  
}

orgconsensus <- function(consensusreslist, 
                         textsize = 13, 
                         titleprefix = NULL, 
                         plot = TRUE){
  
  if(!is.null(titleprefix)){
    titleprefix <- paste0(titleprefix, ' consensus clustering (cluster = ')
  }else{
    titleprefix <- 'Consensus clustering (cluster = '
  }
  
  kres <- list()
  
  i <- 1
  for(i in 1:length(consensusreslist)){
    
    consensusres <- consensusreslist[[i]]$preds
    
    j <- 1
    for(j in 1:length(consensusres)){
      
      singleconsensusres <- consensusres[[j]]$cluster
      
      if(i == 1){
        kres[[j]] <- list()
      }
      
      kres[[j]][[i]] <- singleconsensusres
      
      names(kres[[j]])[i] <- paste0('consensus = ', i)
      
      names(kres)[j] <- names(consensusres)[j]
      
    }
    
  }
  
  kres <- lapply(X = kres, FUN = do.call, what = rbind)
  
  hammings <- list()
  consensusclusterslist <- list()
  k <- 1
  for(k in 1:length(kres)){
    
    kresdat <- kres[[k]]
    
    ncluster <- as.integer(sub(pattern = 'k = ', replacement = '', 
                               x = names(kres)[k]))
    
    hamming <- apply(X = kresdat, MARGIN = 2, 
                     FUN = function(x, dat = kresdat){colSums(x == dat)})
    hamming <- as.matrix(hamming)
    hamming <- hamming/nrow(kresdat)
    hammings[[k]] <- hamming
    names(hammings)[k] <- names(kres)[k]
    
    #heatmap
    
    consensusclusters <- cutree(hclust(dist(hamming)), k = ncluster)
    annos <- data.frame(clusters = consensusclusters, stringsAsFactors = FALSE)
    consensusclusterslist[[k]] <- consensusclusters
    names(consensusclusterslist)[k] <- names(kres)[k]
    
    if(plot == TRUE){
      
      print(
        
        pheatmap::pheatmap(hamming, 
                           color = rev(colorRampPalette(heat.colors(100))(100)), 
                           show_colnames = FALSE, show_rownames = FALSE, 
                           
                           annotation_col = annos, 
                           annotation_names_col = FALSE, 
                           annotation_colors = list(clusters = scales::hue_pal()(ncluster)), 
                           annotation_legend = FALSE, 
                           
                           annotation_row = annos, 
                           annotation_names_row = FALSE, 
                           
                           cluster_cols = TRUE, 
                           cluster_rows = TRUE, 
                           border_color = NA, 
                           
                           main = paste0(titleprefix, 
                                         ncluster, ', consensus = ', nrow(kresdat), ')'), 
                           fontsize = textsize)
        
      )
      
    }
    
  }
  
  
  res <- list(kreses = kres, consensusres = consensusclusterslist)
  
  return(res)
  
}

orgsingle <- function(singlereslist){
  
  kres <- list()
  ivis <- list()
  i <- 1
  for(i in 1:length(singlereslist$preds)){
    
    singleres <- singlereslist$preds[[i]]$cluster
    singleivi <- singlereslist$ivis[[i]]
    
    kres[[i]] <- singleres
    ivis[[i]] <- singleivi
    
    names(kres)[i] <- names(ivis)[i] <- names(singlereslist$preds)[i]
    
  }
  
  res <- list(kreses = kres, 
              ivis = ivis)
  
}

pairedCCA <- function(dat1, 
                      dat2, 
                      varfeatures = NULL, 
                      scale = TRUE, 
                      numcc = 30, 
                      seednum = 2022, 
                      perplexity = 4, 
                      k, 
                      consensus = 1, 
                      threads = 1, 
                      plot = FALSE, 
                      titlefix = NULL, 
                      textsize = 13, 
                      titlesize = 15, 
                      face = 'bold'){
  
  if(!is.null(titlefix)){
    titlesuffix <- paste0(' ', titlefix)
  }else{
    titlesuffix <- NULL
  }
  
  if(!is.null(dat2)){
    
    sharedsamples <- intersect(colnames(dat1), colnames(dat2))
    if(length(sharedsamples) > 0){
      dat1 <- dat1[, sharedsamples, drop = FALSE]
      dat2 <- dat2[, sharedsamples, drop = FALSE]
    }
    
  }
  
  #Variable feature selection
  if(!is.null(varfeatures)){
    
    dat1features <- featuresampling(betas = dat1, 
                                    topfeatures = varfeatures, 
                                    pddat = NULL, 
                                    variancetype = 'sd', 
                                    anova = FALSE)
    
    dat1 <- dat1features$betas[dat1features$varicanceres[row.names(dat1features$betas),] != 0, 
                               , drop = FALSE]
    
    if(!is.null(dat2)){
      
      
      dat2features <- featuresampling(betas = dat2, 
                                      topfeatures = varfeatures, 
                                      pddat = NULL, 
                                      variancetype = 'sd', 
                                      anova = FALSE)
      
      
      dat2 <- dat2features$betas[dat2features$varicanceres[row.names(dat2features$betas),] != 0, 
                                 , drop = FALSE]
    }
    
  }
  
  #Scaling
  if(scale == TRUE){
    
    dat1samples <- colnames(dat1)
    dat1 <- apply(X = dat1, MARGIN = 1, FUN = scale)
    row.names(dat1) <- dat1samples
    dat1 <- t(dat1)
    
    if(!is.null(dat2)){
      
      dat2samples <- colnames(dat2)
      dat2 <- apply(X = dat2, MARGIN = 1, FUN = scale)
      row.names(dat2) <- dat2samples
      dat2 <- t(dat2)
      
    }
    
  }
  
  if(consensus <= 1){
    
    clusterreslist <- singlecluster(i = 1, 
                                    dat1 = dat1, 
                                    dat2 = dat2, 
                                    numcc = numcc, 
                                    seednum = seednum, 
                                    k = k, 
                                    consensuscluster = FALSE)
    
    
    
    singleres <- orgsingle(singlereslist = clusterreslist)
    
    cca.samples <- clusterreslist$cca.samples
    
    res <- singleres
    
  }else{
    
    
    
    signadj <- function(line){
      
      if (sign(line[1]) == -1) {
        line <- line * -1
      }
      return(line)
      
    }
    
    baseCCA <- function(dat1, 
                        dat2, 
                        numcc = 30, 
                        seednum = 2022){
      
      if(is.null(dat2)){
        dat2 <- dat1
      }
      
      set.seed(seednum)
      
      numcc <- min(c(numcc, nrow(dat1), nrow(dat2)))
      
      sampleids <- intersect(colnames(dat1), colnames(dat2))
      
      dat1 <- dat1[, sampleids, drop = FALSE]
      dat2 <- dat2[, sampleids, drop = FALSE]
      sampleids1 <- paste0(sampleids, '_dat1')
      sampleids2 <- paste0(sampleids, '_dat2')
      
      gc()
      mat <- crossprod(x = t(dat1), y = t(dat2))
      
      gc()
      mat <- mat/(ncol(dat1) - 1)
      
      cca.svd <- irlba::irlba(A = mat, nv = numcc)
      gc()
      
      u <- apply(X = cca.svd$u, MARGIN = 2, FUN = signadj)
      v <- apply(X = cca.svd$v, MARGIN = 2, FUN = signadj)
      
      cca.samples1 <- crossprod(x = dat1, y = u)
      cca.samples2 <- crossprod(x = dat2, y = v)
      cca.samples <- rbind(cca.samples1, cca.samples2)
      
      colnames(cca.samples) <- paste0("CC", 1:numcc)
      rownames(cca.samples) <- c(sampleids1, sampleids2)
      
      res <- list()
      res$cca.samples <- cca.samples
      res$u <- u
      res$v <- v
      res$sigular <- cca.svd$d
      
      return(res)
      
    }
    
    calccasamples <- function(dat1, 
                              dat2, 
                              numcc = 20, 
                              seednum = 2022){
      
      #CCA
      cca.res <- baseCCA(dat1 = dat1, 
                         dat2 = dat2, 
                         numcc = numcc, 
                         seednum = seednum)
      
      #weighting by variance
      
      if(!is.null(dat2)){
        ucov <- cov(t(dat1))
        
        gc()
        uvars <- crossprod(cca.res$u, ucov)
        uvars <- crossprod(t(uvars), cca.res$u)
        uvars <- diag(uvars)
        
        vcov <- cov(t(dat2))
        
        gc()
        vvars <- crossprod(cca.res$v, vcov)
        vvars <- crossprod(t(vvars), cca.res$v)
        vvars <- diag(vvars)
        
        uweight <- min(sum(uvars)/nrow(ucov), 1)
        vweight <- min(sum(vvars)/nrow(vcov), 1)
        
        uweight <- uweight/(uweight + vweight)
        vweight <- 1 - uweight
        
        cca.samples1 <- cca.res$cca.samples[1:(0.5*nrow(cca.res$cca.samples)), , 
                                            drop = FALSE]
        cca.samples2 <- cca.res$cca.samples[(0.5*nrow(cca.res$cca.samples) + 1): 
                                              nrow(cca.res$cca.samples), , drop = FALSE]
        
        cca.samples <- uweight*cca.samples1 + vweight*cca.samples2
        row.names(cca.samples) <- gsub(pattern = '_dat1', replacement = '', 
                                       x = row.names(cca.samples))
      }else{
        
        cca.samples <- cca.res$cca.samples[1:(0.5*nrow(cca.res$cca.samples)), , 
                                           drop = FALSE]
      }
      
      row.names(cca.samples) <- gsub(pattern = '_dat1', replacement = '', 
                                     x = row.names(cca.samples))
      
      return(cca.samples)
      
    }
    
    
    singlecluster <- function(i = 1, 
                              dat1, 
                              dat2, 
                              numcc = 20, 
                              seednum = 2022, 
                              k = 2, 
                              consensuscluster = FALSE){
      
      if(consensuscluster == TRUE){
        
        set.seed(seednum + i)
        dat1 <- dat1[sample(x = 1:nrow(dat1), size = nrow(dat1), replace = TRUE), , 
                     drop = FALSE]
        
        if(!is.null(dat2)){
          
          set.seed(seednum + i)
          dat2 <- dat2[sample(x = 1:nrow(dat2), size = nrow(dat2), replace = TRUE), , 
                       drop = FALSE]
          
        }else{
          dat2 <- NULL
        }
        
      }
      
      cca.samples <- calccasamples(dat1 = dat1, 
                                   dat2 = dat2, 
                                   numcc = numcc, 
                                   seednum = seednum)
      
      
      for(i in 1:length(k)){
        
        if(i == 1){
          preds <- list()
          ivis <- list()
          
        }
        
        #clustering
        pred <- kmeans(x = cca.samples, centers = k[i])
        
        preds[[i]] <- pred
        names(preds)[i] <- as.character(paste0('k = ', k[i]))
        
        res <- list(preds = preds)
        
        if(consensuscluster == FALSE){
          #IVI
          ivi <- unlist(clusterCrit::intCriteria(traj = cca.samples,
                                                 part = pred$cluster,
                                                 crit = c('Silhouette',
                                                          'Calinski_Harabasz',
                                                          'Davies_Bouldin')))
          
          
          
          
          ivis[[i]] <- ivi
          names(ivis)[i] <- as.character(paste0('k = ', k[i]))
          
          res <- list(preds = preds, ivis = ivis)
          
        }
        
      }
      
      return(res)
      
    }
    
    if(threads == 1){
      
      clusterreslist <- list()
      i <- 1
      for(i in seq(1, consensus, 1)){
        
        clusterres <- singlecluster(i = i, 
                                    dat1 = dat1, 
                                    dat2 = dat2, 
                                    numcc = numcc, 
                                    seednum = seednum, 
                                    k = k, 
                                    consensuscluster = TRUE)
        clusterreslist[[i]] <- clusterres
        
      }
      
    }else{
      
      #library(doParallel)
      
      cores <- parallel::detectCores()
      cl <- parallel::makeCluster(min(threads, cores))
      
      doParallel::registerDoParallel(cl)
      
      #date()
      `%dopar%` <- foreach::`%dopar%`
      clusterreslist <- foreach::foreach(iseqs = seq(1, consensus, 1), 
                                         .export = NULL) %dopar% {
                                           singlecluster(i = iseqs, 
                                                         dat1 = dat1, 
                                                         dat2 = dat2, 
                                                         numcc = numcc, 
                                                         seednum = seednum, 
                                                         k = k, 
                                                         consensuscluster = TRUE)}
      
      parallel::stopCluster(cl)
      
      unregister_dopar()
      
    }
    
    
    consensusres <- orgconsensus(consensusreslist = clusterreslist, 
                                 textsize = textsize, 
                                 titleprefix = titlefix, 
                                 plot = plot)
    
    
    cca.samples <- calccasamples(dat1 = dat1, 
                                 dat2 = dat2, 
                                 numcc = numcc, 
                                 seednum = seednum)
    
    for(j in 1:length(consensusres$consensusres)){
      
      if(j == 1){
        ivis <- list()
      }
      
      ivi <- unlist(clusterCrit::intCriteria(traj = cca.samples,
                                             part = consensusres$consensusres[[j]],
                                             crit = c('Silhouette',
                                                      'Calinski_Harabasz',
                                                      'Davies_Bouldin')))
      
      
      
      
      ivis[[j]] <- ivi
      names(ivis)[j] <- names(consensusres$consensusres)[j]
      
    }
    
    consensusres$ivis <- ivis
    
    res <- consensusres
    
  }
  
  
  samplenum <- ncol(dat1)
  
  if(plot == TRUE){
    
    j <- 1
    for(j in 1:length(res$kreses)){
      
      ncluster <- as.integer(sub(pattern = 'k = ', replacement = '', 
                                 x = names(res$kreses)[j]))
      subtitle <- paste0(samplenum, ' samples; ', 
                         ncluster, 
                         ' clusters')
      colors <- scales::hue_pal()(ncluster)
      
      if(consensus == 1){
        response <- paste0('Cluster', res$kreses[[j]])
      }else{
        response <- paste0('Cluster', res$consensusres[[j]])
      }
      
      
      pcadat <- tsneplot(tsnematrix = cca.samples, 
                         perplexity = perplexity, 
                         colors = colors, 
                         response = response, 
                         titlesuffix = titlesuffix, 
                         subtitle = subtitle, 
                         responsevarname = 'Cluster', 
                         titlesize = titlesize, 
                         face = face, 
                         textsize = textsize, 
                         labelpoints = FALSE, 
                         embedding = 'pca', 
                         seednum = seednum)
      
      tsnedat <- tsneplot(tsnematrix = cca.samples, 
                          perplexity = perplexity, 
                          colors = colors, 
                          response = response, 
                          titlesuffix = titlesuffix, 
                          subtitle = subtitle, 
                          responsevarname = 'Cluster', 
                          titlesize = titlesize, 
                          face = face, 
                          textsize = textsize, 
                          labelpoints = FALSE, 
                          embedding = 'tsne', 
                          seednum = seednum)
      
      umapdat <- tsneplot(tsnematrix = cca.samples, 
                          perplexity = perplexity, 
                          colors = colors, 
                          response = response, 
                          titlesuffix = titlesuffix, 
                          subtitle = subtitle, 
                          responsevarname = 'Cluster', 
                          titlesize = titlesize, 
                          face = face, 
                          textsize = textsize, 
                          labelpoints = FALSE, 
                          embedding = 'umap', 
                          seednum = seednum)
      
    }
    
  }
  
  
  res$cca.samples <- cca.samples
  
  return(res)
  
}

#'Perform PCA/CCA-based omics clustering
#'
#'Cluster single-omic data with PCA + k-means or cluster multi-omics data with 
#'  CCA + k-means. 
#'
#'@param dats A list with each element as each omic data of the samples to be 
#'  clustered. Each omic data should be a matrix with columns as samples and 
#'  rows as features. The column names are the sample names and the row names 
#'  are the feature names.If only one omic data need to be clustered, it need 
#'  to be a list with one element as the omic data matrix. If multi-omics data 
#'  need to be clustered, the elements (omics) order in the list can influence 
#'  the multi-CCA result, and it is suggested to order the elements (omics) 
#'  according to their feature numbers, i.e., the omic data matrix with the 
#'  largest feature number should be the first element, and the one with the 
#'  second largest feature number should be the second element, and so on.
#'@param varfeatures For each omic data in \code{dats}, the standard deviation 
#'  of the features can be calculated and the top variable features will be 
#'  selected and used for clustering. This parameter defines the number of the 
#'  top variable features will be selected. It should be an integer, such as 
#'  10000, so that the top 10000 most variable features of each omic will be 
#'  selected, and if the omic data have a feature number less than 10000, all 
#'  the features will be used. Can also be NULL, so that this variable feature 
#'  selection step will be skipped and all the omics will use all the features 
#'  for the clustering. Default is NULL.
#'@param scale Whether the features of the omics need to be standardized, so 
#'  that all of them have a mean across all the samples as 0 and a standard 
#'  deviation of 1, and the clustering will be performed after that. Default 
#'  is TRUE.
#'@param numcc If only one omic is provided to \code{dats}, this function will 
#'  perform PCA on the data and use the PCs (principle components) with the 
#'  largest eigen values for the following k-means clustering to get the final 
#'  sample clusters. If mulit-omics data are provided to \code{dats}, CCA will 
#'  be used to compress the multi-omics data, and the function will use the 
#'  top CCs (canonical components) with the largest singular values for the 
#'  following k-means clustering. This parameter \code{numcc} defines how many 
#'  top PCs or CCs will be used for the k-means. Default is 30.
#'@param seednum Random seed number to be used for the CCA/PCA step fulfilled 
#'  with the IRLBA algorithm (implicitly restarted Lanczos bidiagonalization). 
#'  It is also used as the random seed for the feature bootstrapping step of 
#'  the consensus clustering mode, which can be activated by setting another 
#'  parameter \code{consensus} as an integer greater than 1. In addition, it 
#'  is also the random seed for the tSNE and UMAP plots to show the clustering 
#'  result. Default is 2022.
#'@param perplexity The perplexity parameter for the tSNE plot generated to 
#'  show the clustering result. Default is 4. 
#'@param k The CCA/PCA step compresses the multi-omics or single-omic data. 
#'  Then, k-means clustring will be used on the compressed data to get the 
#'  final clustering result. This parameter is used to transfer the cluster 
#'  number to the k-means step. To evaluate the clustering performance, the 
#'  function will calculate several internal validation indices, including 
#'  Silhouette index, Calinski index, and Davies-Bouldin index. If the optimal 
#'  cluster number for the k-means step is not known, can transfer a vector 
#'  with the candiate cluster numbers as elements, so k-means clustering will 
#'  use each of them as the cluster number and evaluate the clustering results 
#'  by calculating the interanl validation indices. The optimal number should 
#'  have the largest Silhouette and Calinski indices, but the smallest Davies- 
#'  Bouldin index. If the best cluster number is known, it can be transferred 
#'  this parameter as a single integer. Then, only it will be used for the 
#'  k-means clustering. Default is an integer 2. 
#'@param consensus If this parameter is set as the integer 1, the normal mode 
#'  will be used, i.e., the original data will be used for the CCA/PCA and the 
#'  k-means clustering steps. If this parameter is set as an integer > 1, the 
#'  consensus mode will be used, and the features in each omic will undergo a 
#'  bootstrapping process for several times. Each time, a bootstrapped dataset 
#'  will be generated and the CCA/PCA and k-means steps will be used on this 
#'  one. This parameter defines the number of bootstrapping. Then, all the 
#'  sample clustering results from the bootstrapping data will be aggregated 
#'  to get the consensus result. The consensus mode can be very time-consuming 
#'  because the CCA/PCA step need to be used on each bootstrapped data. The 
#'  default value of this parameter is 1, so that the clustering will be under 
#'  the normal mode. 
#'@param threads Number of threads need for parallelization. Default is 1.
#'@param plot Whether need to generate the PCA, tSNE, and UMAP plots to show 
#'  the clustering results, and whether to generate the heatmap plots to show 
#'  the consensus mode results. Default is FALSE. 
#'@param titlefix The character string need to be shown in the titles of the 
#'  PCA, tSNE, UMAP, and heatmap plots. It can be set as any character string 
#'  need to be shown in the titles. Default is NULL.
#'@param textsize The font size of the plot texts. Default is 13. 
#'@param titlesize The font size of the PCA, tSNE, and UMAP plots titles. Its 
#'  default is 15.
#'@param face The font face of the PCA, tSNE, and UMAP plots texts. Default is 
#'  "bold".
#'@return A list with several slots will be returned. The slot named "kreses" 
#'  contains the final sample cluster labels for each candidate cluster number 
#'  defined by the parameter \code{k}. The slot "ivis" contains the internal 
#'  validation indices for each cluster number. The "cca.samples" slot is the 
#'  CCA-compressed data. If the consensus clustering mode is used, the slot 
#'  "kreses" will be the sample cluster labels for each bootstrapped dataset, 
#'  and another slot "consensusres" will be the final consensus cluster labels 
#'  for the samples.
#'  
#'  
#'  
#'@examples
#'library(CWGCNA)
#'
#'betas <- system.file("extdata", "placentabetas.rds", package = "CWGCNA")
#'betas <- readRDS(betas)
#'
#'pds <- system.file("extdata", "placentapds.rds", package = "CWGCNA")
#'pds <- readRDS(pds)
#'
#'#Extract the DNAm probe data for the 101 preeclampsia samples
#'prepds <- subset(pds, Group == "Preeclampsia")
#'row.names(prepds) <- 1:nrow(prepds)
#'
#'prebetas <- betas[, prepds$sampleid, drop = FALSE]
#'
#'#Clustering
#'presubtyperes <- multiCCA(dats = list(prebetas), 
#'                          
#'  k = 2, consensus = 1, 
#'                          
#'  seednum = 2022, threads = 6, plot = TRUE, titlefix = "Preeclampsia", 
#'  titlesize = 18, textsize = 16)
#'
#'
#'
#'@export
multiCCA <- function(dats, 
                     varfeatures = NULL, 
                     scale = TRUE, 
                     numcc = 30, 
                     seednum = 2022, 
                     perplexity = 4, 
                     k = c(2), 
                     consensus = 1, 
                     threads = 1, 
                     plot = FALSE, 
                     titlefix = NULL, 
                     textsize = 13, 
                     titlesize = 15, 
                     face = 'bold'){
  
  if(!is.null(titlefix)){
    titlesuffix <- paste0(' ', titlefix)
  }else{
    titlesuffix <- NULL
  }
  
  if(length(dats) > 1){
    
    sharedsamples <- colnames(dats[[1]])
    
    for(i in 2:length(dats)){
      
      sharedsamples <- intersect(sharedsamples, colnames(dats[[i]]))
      
    }
    
    if(length(sharedsamples) > 0){
      
      for(i in 1:length(dats)){
        
        dats[[i]] <- dats[[i]][, sharedsamples, drop = FALSE]
        
      }
      
    }
    
  }
  
  
  
  #Variable feature selection
  if(!is.null(varfeatures)){
    
    for(i in 1:length(dats)){
      
      datfeatures <- featuresampling(betas = dats[[i]], 
                                     topfeatures = varfeatures, 
                                     pddat = NULL, 
                                     variancetype = 'sd', 
                                     anova = FALSE)
      
      dats[[i]] <- datfeatures$betas[datfeatures$varicanceres[row.names(datfeatures$betas),] != 0, 
                                     , drop = FALSE]
      
    }
    
    
  }
  
  #Scaling
  if(scale == TRUE){
    
    for(i in 1:length(dats)){
      
      datsamples <- colnames(dats[[i]])
      dats[[i]] <- apply(X = dats[[i]], MARGIN = 1, FUN = scale)
      row.names(dats[[i]]) <- datsamples
      dats[[i]] <- t(dats[[i]])
      
    }
    
  }
  
  if(consensus <= 1){
    
    clusterreslist <- singlecluster(i = 1, 
                                    dats = dats, 
                                    numcc = numcc, 
                                    seednum = seednum, 
                                    k = k, 
                                    consensuscluster = FALSE)
    
    singleres <- orgsingle(singlereslist = clusterreslist)
    
    cca.samples <- clusterreslist$cca.samples
    
    res <- singleres
    
  }else{
    
    signadj <- function(line){
      
      if (sign(line[1]) == -1) {
        line <- line * -1
      }
      return(line)
      
    }
    
    baseCCA <- function(dat1, 
                        dat2, 
                        numcc = 30, 
                        seednum = 2022){
      
      if(is.null(dat2)){
        dat2 <- dat1
      }
      
      set.seed(seednum)
      
      numcc <- min(c(numcc, nrow(dat1), nrow(dat2)))
      
      sampleids <- intersect(colnames(dat1), colnames(dat2))
      
      dat1 <- dat1[, sampleids, drop = FALSE]
      dat2 <- dat2[, sampleids, drop = FALSE]
      sampleids1 <- paste0(sampleids, '_dat1')
      sampleids2 <- paste0(sampleids, '_dat2')
      
      gc()
      mat <- crossprod(x = t(dat1), y = t(dat2))
      
      gc()
      mat <- mat/(ncol(dat1) - 1)
      
      suppressMessages(
        
        cca.svd <- tryCatch({
          
          irlba::irlba(A = mat, nv = numcc)
          
        }, error = function(err){
          
          svd(x = mat, nu = numcc, nv = numcc)
          
        }, warning = function(war){
          
          svd(x = mat, nu = numcc, nv = numcc)
          
        })
        
      )
      
      gc()
      
      u <- apply(X = cca.svd$u, MARGIN = 2, FUN = signadj)
      v <- apply(X = cca.svd$v, MARGIN = 2, FUN = signadj)
      
      cca.samples1 <- crossprod(x = dat1, y = u)
      cca.samples2 <- crossprod(x = dat2, y = v)
      cca.samples <- rbind(cca.samples1, cca.samples2)
      
      colnames(cca.samples) <- paste0("CC", 1:numcc)
      rownames(cca.samples) <- c(sampleids1, sampleids2)
      
      res <- list()
      res$cca.samples <- cca.samples
      res$u <- u
      res$v <- v
      res$sigular <- cca.svd$d
      
      return(res)
      
    }
    
    calccasamples <- function(dat1, 
                              dat2, 
                              numcc = 20, 
                              seednum = 2022){
      
      #CCA
      cca.res <- baseCCA(dat1 = dat1, 
                         dat2 = dat2, 
                         numcc = numcc, 
                         seednum = seednum)
      
      #weighting by variance
      
      if(!is.null(dat2)){
        ucov <- cov(t(dat1))
        
        gc()
        uvars <- crossprod(cca.res$u, ucov)
        uvars <- crossprod(t(uvars), cca.res$u)
        uvars <- diag(uvars)
        
        vcov <- cov(t(dat2))
        
        gc()
        vvars <- crossprod(cca.res$v, vcov)
        vvars <- crossprod(t(vvars), cca.res$v)
        vvars <- diag(vvars)
        
        uweight <- min(sum(uvars)/nrow(ucov), 1)
        vweight <- min(sum(vvars)/nrow(vcov), 1)
        
        uweight <- uweight/(uweight + vweight)
        vweight <- 1 - uweight
        
        cca.samples1 <- cca.res$cca.samples[1:(0.5*nrow(cca.res$cca.samples)), , 
                                            drop = FALSE]
        cca.samples2 <- cca.res$cca.samples[(0.5*nrow(cca.res$cca.samples) + 1): 
                                              nrow(cca.res$cca.samples), , drop = FALSE]
        
        cca.samples <- uweight*cca.samples1 + vweight*cca.samples2
        row.names(cca.samples) <- gsub(pattern = '_dat1', replacement = '', 
                                       x = row.names(cca.samples))
      }else{
        
        cca.samples <- cca.res$cca.samples[1:(0.5*nrow(cca.res$cca.samples)), , 
                                           drop = FALSE]
      }
      
      row.names(cca.samples) <- gsub(pattern = '_dat1', replacement = '', 
                                     x = row.names(cca.samples))
      
      return(cca.samples)
      
    }
    
    
    singlecluster <- function(i = 1, 
                              dats, 
                              numcc = 20, 
                              seednum = 2022, 
                              k = 2, 
                              consensuscluster = FALSE){
      
      if(consensuscluster == TRUE){
        
        for(j in 1:length(dats)){
          
          set.seed(seednum + i)
          dats[[j]] <- dats[[j]][sample(x = 1:nrow(dats[[j]]), 
                                        size = nrow(dats[[j]]), 
                                        replace = TRUE), , drop = FALSE]
          
        }
        
      }
      
      if(length(dats) == 1){
        
        cca.samples <- calccasamples(dat1 = dats[[1]], 
                                     dat2 = NULL, 
                                     numcc = numcc, 
                                     seednum = seednum)
        
      }else{
        
        cca.samples <- t(dats[[1]])
        
        for(j in 2:length(dats)){
          
          cca.samples <- calccasamples(dat1 = t(cca.samples), 
                                       dat2 = dats[[j]], 
                                       numcc = numcc, 
                                       seednum = seednum)
          
          
        }
        
      }
      
      k <- k[k <= nrow(cca.samples)]
      
      for(i in 1:length(k)){
        
        if(i == 1){
          preds <- list()
          ivis <- list()
          
        }
        
        #clustering
        pred <- kmeans(x = cca.samples, centers = k[i])
        
        preds[[i]] <- pred
        names(preds)[i] <- as.character(paste0('k = ', k[i]))
        
        res <- list(preds = preds)
        
        if(consensuscluster == FALSE){
          #IVI
          ivi <- unlist(clusterCrit::intCriteria(traj = cca.samples,
                                                 part = pred$cluster,
                                                 crit = c('Silhouette',
                                                          'Calinski_Harabasz',
                                                          'Davies_Bouldin')))
          
          
          
          
          ivis[[i]] <- ivi
          names(ivis)[i] <- as.character(paste0('k = ', k[i]))
          
          res <- list(preds = preds, ivis = ivis)
          
        }
        
      }
      
      if(consensuscluster == FALSE){
        
        res$cca.samples <- cca.samples
        
      }
      
      return(res)
      
    }
    
    if(threads == 1){
      
      clusterreslist <- list()
      i <- 1
      for(i in seq(1, consensus, 1)){
        
        clusterres <- singlecluster(i = i, 
                                    dats = dats, 
                                    numcc = numcc, 
                                    seednum = seednum, 
                                    k = k, 
                                    consensuscluster = TRUE)
        clusterreslist[[i]] <- clusterres
        
      }
      
    }else{
      
      #library(doParallel)
      
      cores <- parallel::detectCores()
      cl <- parallel::makeCluster(min(threads, cores))
      
      doParallel::registerDoParallel(cl)
      
      #date()
      `%dopar%` <- foreach::`%dopar%`
      clusterreslist <- foreach::foreach(iseqs = seq(1, consensus, 1), 
                                         .export = NULL) %dopar% {
                                           singlecluster(i = iseqs, 
                                                         dats = dats, 
                                                         numcc = numcc, 
                                                         seednum = seednum, 
                                                         k = k, 
                                                         consensuscluster = TRUE)}
      
      parallel::stopCluster(cl)
      
      unregister_dopar()
      
    }
    
    
    consensusres <- orgconsensus(consensusreslist = clusterreslist, 
                                 textsize = textsize, 
                                 titleprefix = titlefix, 
                                 plot = plot)
    
    
    if(length(dats) == 1){
      
      cca.samples <- calccasamples(dat1 = dats[[1]], 
                                   dat2 = NULL, 
                                   numcc = numcc, 
                                   seednum = seednum)
      
    }else{
      
      cca.samples <- t(dats[[1]])
      
      for(j in 2:length(dats)){
        
        cca.samples <- calccasamples(dat1 = t(cca.samples), 
                                     dat2 = dats[[j]], 
                                     numcc = numcc, 
                                     seednum = seednum)
        
        
      }
      
    }
    
    
    for(j in 1:length(consensusres$consensusres)){
      
      if(j == 1){
        ivis <- list()
      }
      
      ivi <- unlist(clusterCrit::intCriteria(traj = cca.samples,
                                             part = consensusres$consensusres[[j]],
                                             crit = c('Silhouette',
                                                      'Calinski_Harabasz',
                                                      'Davies_Bouldin')))
      
      
      
      
      ivis[[j]] <- ivi
      names(ivis)[j] <- names(consensusres$consensusres)[j]
      
    }
    
    consensusres$ivis <- ivis
    
    res <- consensusres
    
  }
  
  
  samplenum <- ncol(dats[[1]])
  
  if(plot == TRUE){
    
    j <- 1
    for(j in 1:length(res$kreses)){
      
      ncluster <- as.integer(sub(pattern = 'k = ', replacement = '', 
                                 x = names(res$kreses)[j]))
      subtitle <- paste0(samplenum, ' samples; ', 
                         ncluster, 
                         ' clusters')
      colors <- scales::hue_pal()(ncluster)
      
      if(consensus == 1){
        response <- paste0('Cluster', res$kreses[[j]])
      }else{
        response <- paste0('Cluster', res$consensusres[[j]])
      }
      
      
      pcadat <- tsneplot(tsnematrix = cca.samples, 
                         perplexity = perplexity, 
                         colors = colors, 
                         response = response, 
                         titlesuffix = titlesuffix, 
                         subtitle = subtitle, 
                         responsevarname = 'Cluster', 
                         titlesize = titlesize, 
                         face = face, 
                         textsize = textsize, 
                         labelpoints = FALSE, 
                         embedding = 'pca', 
                         seednum = seednum)
      
      tsnedat <- tsneplot(tsnematrix = cca.samples, 
                          perplexity = perplexity, 
                          colors = colors, 
                          response = response, 
                          titlesuffix = titlesuffix, 
                          subtitle = subtitle, 
                          responsevarname = 'Cluster', 
                          titlesize = titlesize, 
                          face = face, 
                          textsize = textsize, 
                          labelpoints = FALSE, 
                          embedding = 'tsne', 
                          seednum = seednum)
      
      umapdat <- tsneplot(tsnematrix = cca.samples, 
                          perplexity = perplexity, 
                          colors = colors, 
                          response = response, 
                          titlesuffix = titlesuffix, 
                          subtitle = subtitle, 
                          responsevarname = 'Cluster', 
                          titlesize = titlesize, 
                          face = face, 
                          textsize = textsize, 
                          labelpoints = FALSE, 
                          embedding = 'umap', 
                          seednum = seednum)
      
    }
    
  }
  
  
  res$cca.samples <- cca.samples
  
  return(res)
  
}

calwgcnasamples <- function(dat, 
                            seednum = 2022, 
                            
                            sftpowers = NULL, 
                            powers = seq(1, 20, 1), 
                            rsqcutline = 0.8, 
                            minclustersize = 50, 
                            mergecutheight = 0.2, 
                            maxblocksize = 5000, 
                            
                            plot = FALSE, 
                            titleprefix, 
                            threads = 1){
  
  if(!is.null(sftpowers)){
    
    if(length(unique(sftpowers)) == 1){
      
      wgcnares <- wgcnamain(dat = dat, 
                            sftpower = unique(sftpowers), 
                            minclustersize = minclustersize, 
                            mergecutheight = mergecutheight, 
                            returntom = FALSE, 
                            plot = plot, 
                            featuretype = 'feature', 
                            plottitleprefix = titleprefix, 
                            threads = threads, 
                            maxblocksize = maxblocksize, 
                            seed = seednum)
      
      
    }else{
      
      sftpowers <- NULL
      
    }
    
  }
  
  if(is.null(sftpowers)){
    
    sftpowers <- choosepower(datlist = list(dat), 
                             i = 1, 
                             powers = powers, 
                             rsqcutline = rsqcutline, 
                             plot = FALSE)
    
    if(is.null(sftpowers)){
      return(NULL)
    }
    
    wgcnares <- wgcnamain(dat = dat, 
                          sftpower = as.vector(unlist(sftpowers$powerEstimate)), 
                          minclustersize = minclustersize, 
                          mergecutheight = mergecutheight, 
                          returntom = FALSE, 
                          plot = plot, 
                          featuretype = 'feature', 
                          plottitleprefix = titleprefix, 
                          threads = threads, 
                          maxblocksize = maxblocksize, 
                          seed = seednum)
    
  }
  
  gc()
  
  wgcna.samples <- wgcnares$mergedmelist
  
  return(wgcna.samples)
  
}

singlewgcnacluster <- function(i = 1, 
                               dats, 
                               seednum = 2022, 
                               k = c(2), 
                               consensuscluster = FALSE, 
                               
                               sftpowers = NULL, 
                               powers = seq(1, 20, 1), 
                               rsqcutline = 0.8, 
                               minclustersize = 50, 
                               mergecutheight = 0.2, 
                               maxblocksize = 5000, 
                               
                               plot = FALSE, 
                               titlefix, 
                               threads = 1){
  
  if(consensuscluster == TRUE){
    
    plot <- FALSE
    
    for(j in 1:length(dats)){
      
      set.seed(seednum + i)
      dats[[j]] <- dats[[j]][sample(x = 1:nrow(dats[[j]]), 
                                    size = nrow(dats[[j]]), 
                                    replace = TRUE), , drop = FALSE]
      
    }
    
  }
  
  if(length(dats) == 1){
    
    wgcna.samples <- calwgcnasamples(dat = dats[[1]], 
                                     seednum = seednum, 
                                     
                                     sftpowers = sftpowers, 
                                     powers = powers, 
                                     rsqcutline = rsqcutline, 
                                     minclustersize = minclustersize, 
                                     mergecutheight = mergecutheight, 
                                     maxblocksize = maxblocksize, 
                                     
                                     plot = plot, 
                                     titleprefix = paste0(titlefix, ' dataset'), 
                                     threads = threads)
    
    if(is.null(wgcna.samples)){
      return(NULL)
    }
    
  }else{
    
    wgcna.samples.list <- list()
    
    for(j in 1:length(dats)){
      
      wgcna.samples <- calwgcnasamples(dat = dats[[j]], 
                                       seednum = seednum, 
                                       
                                       sftpowers = sftpowers, 
                                       powers = powers, 
                                       rsqcutline = rsqcutline, 
                                       minclustersize = minclustersize, 
                                       mergecutheight = mergecutheight, 
                                       maxblocksize = maxblocksize, 
                                       
                                       plot = plot, 
                                       titleprefix = paste0(titlefix, ' dataset', j), 
                                       threads = threads)
      
      if(!is.null(wgcna.samples)){
        
        colnames(wgcna.samples) <- paste0(colnames(wgcna.samples), 
                                          '_dataset', j)
        
      }
      
      wgcna.samples.list[[j]] <- wgcna.samples
      
    }
    
    
    wgcna.samples <- do.call(cbind, wgcna.samples.list)
    
    if(is.null(wgcna.samples)){
      return(NULL)
    }
    
  }
  
  ME0idx <- grepl(pattern = '^ME0', x = colnames(wgcna.samples))
  wgcna.samples <- wgcna.samples[, !ME0idx, drop = FALSE]
  
  if(is.null(wgcna.samples)){
    return(NULL)
  }
  
  if(nrow(wgcna.samples) == 0){
    return(NULL)
  }
  
  wgcna.samples <- as.matrix(wgcna.samples)
  
  k <- k[k <= nrow(wgcna.samples)]
  
  for(i in 1:length(k)){
    
    if(i == 1){
      preds <- list()
      ivis <- list()
      
    }
    
    #clustering
    pred <- kmeans(x = wgcna.samples, centers = k[i])
    
    preds[[i]] <- pred
    names(preds)[i] <- as.character(paste0('k = ', k[i]))
    
    res <- list(preds = preds)
    
    if(consensuscluster == FALSE){
      #IVI
      ivi <- unlist(clusterCrit::intCriteria(traj = wgcna.samples,
                                             part = pred$cluster,
                                             crit = c('Silhouette',
                                                      'Calinski_Harabasz',
                                                      'Davies_Bouldin')))
      
      
      
      
      ivis[[i]] <- ivi
      names(ivis)[i] <- as.character(paste0('k = ', k[i]))
      
      res <- list(preds = preds, ivis = ivis)
      
    }
    
  }
  
  if(consensuscluster == FALSE){
    
    res$wgcna.samples <- wgcna.samples
    
  }
  
  return(res)
  
}

#'Perform WGCNA-based omics clustering
#'
#'Cluster single-omic or multi-omics data with WGCNA + k-means. 
#'
#'@param dats A list with each element as each omic data of the samples to be 
#'  clustered. Each omic data should be a matrix with columns as samples and 
#'  rows as features. The column names are the sample names and the row names 
#'  are the feature names.If only one omic data need to be clustered, it need 
#'  to be a list with one element as the omic data matrix. 
#'@param varfeatures For each omic data in \code{dats}, the standard deviation 
#'  of the features can be calculated and the top variable features will be 
#'  selected and used for clustering. This parameter defines the number of the 
#'  top variable features will be selected. It should be an integer, such as 
#'  10000, so that the top 10000 most variable features of each omic will be 
#'  selected, and if the omic data have a feature number less than 10000, all 
#'  the features will be used. Can also be NULL, so that this variable feature 
#'  selection step will be skipped and all the omics will use all the features 
#'  for the clustering. Default is 5000, so that for each omic, the top 5000 
#'  most variable features will be used for WGCNA module calling, and then, 
#'  the WGCNA module eigengenes from all the omics will be combined together 
#'  into one compressed data followed by k-means clustering.
#'@param scale Whether the features of the omics need to be standardized, so 
#'  that all of them have a mean across all the samples as 0 and a standard 
#'  deviation of 1, and the clustering will be performed after that. Default 
#'  is TRUE.
#'@param sftpowers The soft-thresholding power to be used for WGCNA network 
#'  construction. Should be an integer or NULL. Default is NULL, and in this 
#'  case, for each omic, the function will respectively calculate and pick up 
#'  an appropriate power for it, from the values provided to another parameter 
#'  \code{powers}. If \code{sftpowers} is an integer, all the omic data will 
#'  use this same soft-thresholding power for the WGCNA network construction. 
#'@param powers A vector used to provide candidate soft-thresholding power for 
#'  WGCNA network construction. Should be a vector with integers as elements. 
#'  If the parameter \code{sftpowers} is NULL, for each omic, the function 
#'  will calculate the scale free topology fitting index for each element of 
#'  this vector and use the one with the largest index as the final power. Its 
#'  default value is the integer sequence from 1 to 20.
#'@param rsqcutline A float defining the desired minimum scale free topology 
#'  fitting index to select the soft-thresholding power from the candidates 
#'  included in the parameter \code{powers}, in the case that the parameter 
#'  \code{sftpowers} is NULL and the soft-thresholding powers for the omics 
#'  need to be calculated by the function. Default is 0.8.
#'@param minclustersize The minimum module size for module detection during 
#'  WGCNA network construction. Default is 50. 
#'@param mergecutheight The dendrogram cut height for WGCNA module merging. It 
#'  default is 0.2.
#'@param maxblocksize Integer giving maximum block size for module detection 
#'  during WGCNA network construction. Default is 50. 
#'@param seednum Random seed number to be used for the WGCNA module calling 
#'  step. It is also used as the random seed for the feature bootstrapping 
#'  step of the consensus clustering mode, which can be activated by setting 
#'  another parameter \code{consensus} as an integer greater than 1. It is 
#'  also the random seed for the tSNE and UMAP plots to show the clustering 
#'  result. Default is 2022.
#'@param perplexity The perplexity parameter for the tSNE plot generated to 
#'  show the clustering result. Default is 4. 
#'@param k The WGCNA step compresses the multi-omics or single-omic data by 
#'  combining the omics WGCNA eigengenes. Then, k-means clustring will be used 
#'  on the compressed data to get the final clustering result. This parameter 
#'  is used to transfer the cluster number to the k-means step. To evaluate 
#'  the clustering performance, the function will calculate several internal 
#'  validation indices, including Silhouette index, Calinski-Harabasz index, 
#'  and Davies-Bouldin index. If the optimal cluster number for the k-means 
#'  step is not known, can transfer a vector with the candiate cluster numbers 
#'  as elements, so k-means clustering will use each of them as the cluster 
#'  number and evaluate their clustering results by calculating the interanl 
#'  validation indices. The optimal number should have the largest Silhouette 
#'  and Calinski-Harabasz indices, but the smallest Davies-Bouldin index. If 
#'  the best cluster number is known, it can be transferred to the parameter a 
#'  single integer so only it will be used for the k-means clustering. Default 
#'  is an integer 2. 
#'@param consensus If this parameter is set as the integer 1, the normal mode 
#'  will be used, i.e., the original data will be used for the WGCNA and the 
#'  k-means clustering steps. If this parameter is set as an integer > 1, the 
#'  consensus mode will be used, and the features in each omic will undergo a 
#'  bootstrapping process for several times. Each time, a bootstrapped dataset 
#'  will be generated and the WGCNA and k-means steps will be used on this 
#'  one. This parameter defines the number of bootstrapping. Then, all the 
#'  sample clustering results from the bootstrapping data will be aggregated 
#'  to get the consensus result. The consensus mode can be very time-consuming 
#'  because the WGCNA step need to be performed on each bootstrapped data. The 
#'  default value of this parameter is 1, so that the clustering will be under 
#'  the normal mode. 
#'@param threads Number of threads need for parallelization. Default is 1.
#'@param plot Whether need to generate the WGCNA dendrogram plot during the 
#'  WGCNA module calling, whether to generate the PCA, tSNE, and UMAP plots to 
#'  show the clustering results, and whether to generate the heatmap plots to 
#'  show the consensus mode results. Default is FALSE.
#'@param titlefix The character string need to be shown in the titles of the 
#'  WGCNA dendrogram, PCA, tSNE, UMAP, and heatmap plots. It can be set as any 
#'  character string need to be shown in the titles. Default is NULL.
#'@param textsize The font size of the plot texts. Default is 13. 
#'@param titlesize The font size of the PCA, tSNE, and UMAP plots titles. Its 
#'  default is 15.
#'@param face The font face of the PCA, tSNE, and UMAP plots texts. Default is 
#'  "bold".
#'@return A list with several slots will be returned. The slot named "kreses" 
#'  contains the final sample cluster labels for each candidate cluster number 
#'  defined by the parameter \code{k}. The slot "ivis" contains the internal 
#'  validation indices for each cluster number. The slot "wgcna.samples" has 
#'  the WGCNA-compressed data. If the consensus clustering mode is used, the 
#'  slot "kreses" will contain the sample cluster labels for each bootstrapped 
#'  dataset, and another slot "consensusres" will contain the final consensus 
#'  cluster labels for the samples.
#'@export
wgcnacluster <- function(dats, 
                         varfeatures = 5000, 
                         scale = TRUE, 
                         
                         sftpowers = NULL, 
                         powers = seq(1, 20, 1), 
                         rsqcutline = 0.8, 
                         minclustersize = 50, 
                         mergecutheight = 0.2, 
                         maxblocksize = 5000, 
                         
                         seednum = 2022, 
                         perplexity = 4, 
                         k = c(2), 
                         consensus = 1, 
                         threads = 1, 
                         plot = TRUE, 
                         titlefix = NULL, 
                         textsize = 13, 
                         titlesize = 15, 
                         face = 'bold'){
  
  if(!is.null(titlefix)){
    titlesuffix <- paste0(' ', titlefix)
  }else{
    titlesuffix <- NULL
  }
  
  if(length(dats) > 1){
    
    sharedsamples <- colnames(dats[[1]])
    
    for(i in 2:length(dats)){
      
      sharedsamples <- intersect(sharedsamples, colnames(dats[[i]]))
      
    }
    
    if(length(sharedsamples) > 0){
      
      for(i in 1:length(dats)){
        
        dats[[i]] <- dats[[i]][, sharedsamples, drop = FALSE]
        
      }
      
    }
    
  }
  
  #Variable feature selection
  if(!is.null(varfeatures)){
    
    for(i in 1:length(dats)){
      
      datfeatures <- featuresampling(betas = dats[[i]], 
                                     topfeatures = varfeatures, 
                                     pddat = NULL, 
                                     variancetype = 'sd', 
                                     anova = FALSE)
      
      dats[[i]] <- datfeatures$betas[datfeatures$varicanceres[row.names(datfeatures$betas),] != 0, 
                                     , drop = FALSE]
      
    }
    
    
  }
  
  #Scaling
  if(scale == TRUE){
    
    for(i in 1:length(dats)){
      
      datsamples <- colnames(dats[[i]])
      dats[[i]] <- apply(X = dats[[i]], MARGIN = 1, FUN = scale)
      row.names(dats[[i]]) <- datsamples
      dats[[i]] <- t(dats[[i]])
      
    }
    
  }
  
  if(consensus <= 1){
    
    clusterreslist <- singlewgcnacluster(i = 1, 
                                         dats = dats, 
                                         seednum = seednum, 
                                         k = k, 
                                         consensuscluster = FALSE, 
                                         
                                         sftpowers = sftpowers, 
                                         powers = powers, 
                                         rsqcutline = rsqcutline, 
                                         minclustersize = minclustersize, 
                                         mergecutheight = mergecutheight, 
                                         maxblocksize = maxblocksize, 
                                         
                                         plot = plot, 
                                         titlefix = titlefix, 
                                         threads = threads)
    
    if(is.null(clusterreslist)){
      return(NULL)
    }
    
    singleres <- orgsingle(singlereslist = clusterreslist)
    
    wgcna.samples <- clusterreslist$wgcna.samples
    
    res <- singleres
    
  }else{
    
    clusterreslist <- list()
    i <- 1
    for(i in seq(1, consensus, 1)){
      
      clusterres <- singlewgcnacluster(i = i, 
                                       dats = dats, 
                                       seednum = seednum, 
                                       k = k, 
                                       consensuscluster = TRUE, 
                                       
                                       sftpowers = sftpowers, 
                                       powers = powers, 
                                       rsqcutline = rsqcutline, 
                                       minclustersize = minclustersize, 
                                       mergecutheight = mergecutheight, 
                                       maxblocksize = maxblocksize, 
                                       
                                       plot = FALSE, 
                                       titlefix = titlefix, 
                                       threads = threads)
      
      
      clusterreslist[[i]] <- clusterres
      
    }
    
    if(length(clusterreslist) == 0){
      return(NULL)
    }
    
    
    consensusres <- orgconsensus(consensusreslist = clusterreslist, 
                                 textsize = textsize, 
                                 titleprefix = titlefix, 
                                 plot = plot)
    
    if(length(dats) == 1){
      
      wgcna.samples <- calwgcnasamples(dat = dats[[1]], 
                                       seednum = seednum, 
                                       
                                       sftpowers = sftpowers, 
                                       powers = powers, 
                                       rsqcutline = rsqcutline, 
                                       minclustersize = minclustersize, 
                                       mergecutheight = mergecutheight, 
                                       maxblocksize = maxblocksize, 
                                       
                                       plot = plot, 
                                       titleprefix = paste0(titlefix, ' dataset'), 
                                       threads = threads)
      
    }else{
      
      wgcna.samples.list <- list()
      
      for(j in 1:length(dats)){
        
        wgcna.samples <- calwgcnasamples(dat = dats[[j]], 
                                         seednum = seednum, 
                                         
                                         sftpowers = sftpowers, 
                                         powers = powers, 
                                         rsqcutline = rsqcutline, 
                                         minclustersize = minclustersize, 
                                         mergecutheight = mergecutheight, 
                                         maxblocksize = maxblocksize, 
                                         
                                         plot = plot, 
                                         titleprefix = paste0(titlefix, ' dataset', j), 
                                         threads = threads)
        
        if(!is.null(wgcna.samples)){
          
          colnames(wgcna.samples) <- paste0(colnames(wgcna.samples), 
                                            '_dataset', j)
          
        }
        
        wgcna.samples.list[[j]] <- wgcna.samples
        
      }
      
      
      wgcna.samples <- do.call(cbind, wgcna.samples.list)
      
    }
    
    if(!is.null(wgcna.samples)){
      
      ME0idx <- grepl(pattern = '^ME0', x = colnames(wgcna.samples))
      wgcna.samples <- wgcna.samples[, !ME0idx, drop = FALSE]
      
      wgcna.samples <- as.matrix(wgcna.samples)
      
      for(j in 1:length(consensusres$consensusres)){
        
        if(j == 1){
          ivis <- list()
        }
        
        ivi <- unlist(clusterCrit::intCriteria(traj = wgcna.samples,
                                               part = consensusres$consensusres[[j]],
                                               crit = c('Silhouette',
                                                        'Calinski_Harabasz',
                                                        'Davies_Bouldin')))
        
        
        
        
        ivis[[j]] <- ivi
        names(ivis)[j] <- names(consensusres$consensusres)[j]
        
      }
      
      consensusres$ivis <- ivis
      
    }
    
    res <- consensusres
    
  }
  
  
  samplenum <- ncol(dats[[1]])
  
  if(plot == TRUE){
    
    j <- 1
    for(j in 1:length(res$kreses)){
      
      ncluster <- as.integer(sub(pattern = 'k = ', replacement = '', 
                                 x = names(res$kreses)[j]))
      subtitle <- paste0(samplenum, ' samples; ', 
                         ncluster, 
                         ' clusters')
      colors <- scales::hue_pal()(ncluster)
      
      if(consensus == 1){
        response <- paste0('Cluster', res$kreses[[j]])
      }else{
        response <- paste0('Cluster', res$consensusres[[j]])
      }
      
      
      pcadat <- tsneplot(tsnematrix = wgcna.samples, 
                         perplexity = perplexity, 
                         colors = colors, 
                         response = response, 
                         titlesuffix = titlesuffix, 
                         subtitle = subtitle, 
                         responsevarname = 'Cluster', 
                         titlesize = titlesize, 
                         face = face, 
                         textsize = textsize, 
                         labelpoints = FALSE, 
                         embedding = 'pca', 
                         seednum = seednum)
      
      tsnedat <- tsneplot(tsnematrix = wgcna.samples, 
                          perplexity = perplexity, 
                          colors = colors, 
                          response = response, 
                          titlesuffix = titlesuffix, 
                          subtitle = subtitle, 
                          responsevarname = 'Cluster', 
                          titlesize = titlesize, 
                          face = face, 
                          textsize = textsize, 
                          labelpoints = FALSE, 
                          embedding = 'tsne', 
                          seednum = seednum)
      
      umapdat <- tsneplot(tsnematrix = wgcna.samples, 
                          perplexity = perplexity, 
                          colors = colors, 
                          response = response, 
                          titlesuffix = titlesuffix, 
                          subtitle = subtitle, 
                          responsevarname = 'Cluster', 
                          titlesize = titlesize, 
                          face = face, 
                          textsize = textsize, 
                          labelpoints = FALSE, 
                          embedding = 'umap', 
                          seednum = seednum)
      
    }
    
  }
  
  
  res$wgcna.samples <- wgcna.samples
  
  return(res)
  
}


#Parallelization####

#Some parallel computing going on in the background that is not getting cleaned up
#fully between runs can cause Error in summary.connection(connection) : invalid
#connection. The following function is needed to be called to fix this error.
unregister_dopar <- function(){
  
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  
}

mclapply.win <- function(X, FUN, ..., mc.cores){
  
  cl <- parallel::makeCluster(getOption("cl.cores", mc.cores))
  
  res <- parallel::parLapply(cl = cl, X = X, fun = FUN, ...)
  
  parallel::stopCluster(cl)
  
  unregister_dopar()
  
  return(res)
  
}

mclapply <- function(X, FUN, ..., mc.cores){
  
  cores <- parallel::detectCores()
  
  cl <- min(mc.cores, cores - 1)
  
  doParallel::registerDoParallel(cl)
  
  if(Sys.info()['sysname'] == 'Windows'){
    
    res <- mclapply.win(X = X, mc.cores = cl, FUN = FUN, ...)
    
  }else{
    
    res <- parallel::mclapply(X = X, mc.cores = cl, FUN = FUN, ...)
    
  }
  
  return(res)
  
}


#SVM####

subfunc_svm_e1071_linear_train_tuner_mc <- function(data.xTrain,
                                                    target.yTrain,
                                                    mod.type = "C-classification",
                                                    kernel. = "linear",
                                                    scale. = T,
                                                    C.base = 10,
                                                    C.min = -3,
                                                    C.max = 3,
                                                    n.CV = 5,
                                                    verbose = T,
                                                    seed = 1234,
                                                    parallel = T,
                                                    mc.cores = 4L, 
                                                    weighted = FALSE){
  
  
  Cost.l <- as.list(C.base^(C.min:C.max))
  
  cvfit.e1071.linear.C.tuner <- list()
  
  basesvm <- function(i,
                      seed,
                      data.xTrain,
                      target.yTrain,
                      scale.,
                      mod.type,
                      kernel.,
                      Cost.l,
                      n.CV, 
                      weighted = FALSE){
    
    if(weighted == TRUE){
      
      classweights <- 100/table(target.yTrain)
      
      set.seed(seed + 1, kind ="default")
      
      res <- e1071::svm(x = data.xTrain,
                        y = target.yTrain,
                        scale = scale.,
                        type = mod.type,
                        kernel = kernel.,
                        cross = n.CV,
                        probability = TRUE,
                        fitted = TRUE, 
                        
                        class.weights = classweights)
      
    }else{
      
      classweights <- NULL
      
      set.seed(seed + 1, kind ="default")
      
      res <- e1071::svm(x = data.xTrain,
                        y = target.yTrain,
                        scale = scale.,
                        type = mod.type,
                        kernel = kernel.,
                        cost = Cost.l[[i]],
                        cross = n.CV,
                        probability = TRUE,
                        fitted = TRUE, 
                        
                        class.weights = classweights)
      
    }
    
    return(res)
    
  }
  
  # Parallel
  if(parallel){
    
    cvfit.e1071.linear.C.tuner <- mclapply(X = seq_along(Cost.l),
                                           FUN = basesvm,
                                           seed = seed,
                                           data.xTrain = data.xTrain,
                                           target.yTrain = target.yTrain,
                                           scale. = scale.,
                                           mod.type = mod.type,
                                           kernel. = kernel.,
                                           Cost.l = Cost.l,
                                           n.CV = n.CV, 
                                           weighted = weighted, 
                                           #mc.preschedule = T,
                                           #mc.set.seed = T,
                                           mc.cores = mc.cores)
    
    unregister_dopar()
    
    return(cvfit.e1071.linear.C.tuner)
    
    # Sequential
  } else {
    cvfit.e1071.linear.C.tuner <- lapply(seq_along(Cost.l),
                                         FUN = basesvm,
                                         seed = seed,
                                         data.xTrain = data.xTrain,
                                         target.yTrain = target.yTrain,
                                         scale. = scale.,
                                         mod.type = mod.type,
                                         kernel. = kernel.,
                                         Cost.l = Cost.l,
                                         n.CV = n.CV, 
                                         weighted = weighted)
    
    return(cvfit.e1071.linear.C.tuner)
  }
}

subfunc_svm_e1071_linear_C_selector <- function(results.cvfit.e1071.linear.C.tuner,
                                                C.base = 10, 
                                                C.min = -3, 
                                                C.max = -2,
                                                n.CV = 5, 
                                                verbose = TRUE){
  
  Costs.l <- as.list(C.base^(C.min:C.max))
  
  res.cvfit.svm.accuracies.nCV <- 
    sapply(seq_along(C.base^(C.min:C.max)), 
           function(i){simplify2array(results.cvfit.e1071.linear.C.tuner[[i]]$accuracies)})
  
  colnames(res.cvfit.svm.accuracies.nCV) <- paste0("Cost_", Costs.l)
  
  rownames(res.cvfit.svm.accuracies.nCV) <- paste0("nCV", seq(1, n.CV, 1))
  
  if(verbose){
    message("\nMatrix of all CV accuracies:")
    print(res.cvfit.svm.accuracies.nCV)
  }
  
  # Average accuracy
  res.cvfit.svm.accuracies.mean <- sapply(seq_along(C.base^(C.min:C.max)), function(i){
    simplify2array(results.cvfit.e1071.linear.C.tuner[[i]]$tot.accuracy)})
  names(res.cvfit.svm.accuracies.mean) <- paste0("Cost_", Costs.l)
  # Same as: res.cvfit.svm.accuracies.mean <- apply(res.cvfit.svm.accuracies.nCV, 2, mean)
  # Print list of average CV accuracies/ $tot.accuracy
  if(verbose){
    message("\nMean CV accuracies:")
    print(res.cvfit.svm.accuracies.mean)
  }
  
  # Selection
  # Chooses the smallest C with highest 5-fold cross-validated accuracy among possible choices
  # => if C is large enough anyway (see Appendix-N.5.) doesnt make a difference
  # => saves also computation time if C is smaller # => error-margin/Nr of supp.vecs.
  C.selected <- Costs.l[[which.max(res.cvfit.svm.accuracies.mean)]]
  message("\nCost parameter with highest ", n.CV, "-fold CV accuracy : C = ", C.selected, " ; ",
          "\n Note: If more than one maximal accuracy exists, C returns the smallest cost parameter with highest accuracy.",
          "\n Once C is large than a certain value, the obtained models have similar performances",
          " (for theoretical proof see Theorem 3 of Keerthi and Lin, 2003)")
  res <- list(C.selected = C.selected,
              mtx.accuracies.nCV = res.cvfit.svm.accuracies.nCV,
              mtx.accuracies.mean = res.cvfit.svm.accuracies.mean)
  return(res)
  
  # Literature:
  # Important # A practical guide to LIBLINEAR - Fan, Chang, Hsiesh, Wang and Lin 2008
  # <https://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf>
  # Appendix-N.5. Parameter Selection:
  # 1. Solvers in LIBLINEAR are not very sensitive to C. Once C is large than a certain value,
  # the obtained models have similar performances.
  # Theoretical proof: Theorem 3 of Keerthi and Lin (2003)
}

subfunc_svm_e1071_linear_modfit_train <- function(C.tuned,
                                                  data.xTrain, target.yTrain,
                                                  results.cvfit.e1071.linear.C.tuner,
                                                  C.selector.accuracy.mean,
                                                  use.fitted = TRUE){  #res.svm.C.tuner.l
  
  message("\n\nRe-fitting training data ... ", Sys.time())
  #Costs.l <- as.list(C.base^(C.min:C.max))
  i <- which.max(C.selector.accuracy.mean)
  
  # Ver.1 - Use predict function to refit training data
  # Note If the training set was scaled by svm (done by default),
  # the new data is scaled accordingly using scale and center of the training data.
  modfit.train.svm.lin.pred <- tryCatch({
    
    predict(object = results.cvfit.e1071.linear.C.tuner[[i]],
            newdata =  data.xTrain,
            decision.values = T,
            # decision values of all binary classif. in multiclass setting are returned.
            probability = T)
    
  }, error = function(err){
    
    e1071:::predict.svm(object = results.cvfit.e1071.linear.C.tuner[[i]],
                        newdata =  data.xTrain,
                        decision.values = T,
                        # decision values of all binary classif. in multiclass setting are returned.
                        probability = T)
    
  })
  
  # Ver.2 - Use fitted() - see ??svm - Examples
  if(use.fitted){modfit.train.svm.lin.fitted <- fitted(results.cvfit.e1071.linear.C.tuner[[i]])}
  
  # Output
  res <- list(svm.e1071.model.object = results.cvfit.e1071.linear.C.tuner[[i]],
              trainfit.svm.lin1 = modfit.train.svm.lin.pred,
              trainfit.svm.lin2 = modfit.train.svm.lin.fitted) # => output file is rel. large / lot of large mtx.
  return(res)
}


train_SVM_e1071_LK <- function(y, 
                               betas.Train,
                               seed,
                               mc.cores,
                               nfolds = 5,
                               C.base = 10, 
                               C.min = -3, 
                               C.max = -2,
                               scale.internally.by.e1071.svm = TRUE,
                               mod.type = "C-classification", 
                               weighted = FALSE){
  
  
  set.seed(seed, kind = "default")
  
  cvfit.svm.e1071.linear.C.tuner <- subfunc_svm_e1071_linear_train_tuner_mc(data.xTrain = betas.Train,
                                                                            target.yTrain = y,
                                                                            mod.type = mod.type,
                                                                            kernel. = "linear",
                                                                            scale. = scale.internally.by.e1071.svm,
                                                                            C.base = C.base,
                                                                            C.min = C.min, 
                                                                            C.max = C.max,
                                                                            n.CV = nfolds,
                                                                            verbose = TRUE,
                                                                            seed = seed,
                                                                            parallel = TRUE,
                                                                            mc.cores = mc.cores, 
                                                                            weighted = weighted)
  
  
  # Extract optimal C or smallest C with highest accuracy
  C.tuned.cv <-  subfunc_svm_e1071_linear_C_selector(results.cvfit.e1071.linear.C.tuner = cvfit.svm.e1071.linear.C.tuner,
                                                     C.base = C.base,
                                                     C.min = C.min,
                                                     C.max = C.max,
                                                     n.CV = nfolds,
                                                     verbose = TRUE)
  
  
  # C.tuned.cv = list of 3: $C.selected $mtx.accuracies.nCV $mtx.accuracies.mean
  # Provide message with value
  message(paste0("Optimal cost (C) parameter: ", C.tuned.cv$C.selected))
  
  modfit.svm.linear.train <- subfunc_svm_e1071_linear_modfit_train(C.tuned = C.tuned.cv$C.selected,
                                                                   data.xTrain = betas.Train,
                                                                   target.yTrain = y,
                                                                   results.cvfit.e1071.linear.C.tuner = cvfit.svm.e1071.linear.C.tuner,
                                                                   C.selector.accuracy.mean = C.tuned.cv$mtx.accuracies.mean,
                                                                   use.fitted = TRUE)
  # uses predict supposed to scale data.xTrain / betas.Train automatically
  
  
  # CAVE conames order is not the same as in levels(y) !!!
  pred.scores.trainfit.svm.lin1 <- attr(modfit.svm.linear.train$trainfit.svm.lin1, "probabilities")
  
  # Results
  res <- list(modfit.svm.linear.train$svm.e1071.model.object,
              modfit.svm.linear.train$trainfit.svm.lin1,
              pred.scores.trainfit.svm.lin1,
              modfit.svm.linear.train$trainfit.svm.lin2,
              cvfit.svm.e1071.linear.C.tuner,
              C.tuned.cv)
  return(res)
}

#Classification####

singleclassifiernet <- function(trainvar, 
                                trainres, 
                                threads = 1, 
                                alphas = c(0.5), 
                                seednum = 2022, 
                                foldnum = 5, 
                                errortype = 'min', 
                                verbose = TRUE){
  
  a <- alphas
  #Tune the value of alpha through a line search with the parallelism
  
  for(k in 1:length(a)){
    l <- a[k]
    
    if(threads > 1){
      
      parallelval <- TRUE
      cores <- parallel::detectCores()
      cl <- parallel::makeCluster(min(threads, cores))
      
      doParallel::registerDoParallel(cl)
      
    }else{
      parallelval <- FALSE
    }
    
    set.seed(seednum)
    cv <- glmnet::cv.glmnet(x = trainvar, y = trainres$Response, family = "multinomial", 
                            nfold = foldnum, type.measure = "deviance", 
                            paralle = parallelval, alpha = l)
    
    if(threads > 1){
      parallel::stopCluster(cl)
      unregister_dopar()
    }
    
    
    if(k == 1){
      cvms.1ses <- c(cv$cvm[cv$lambda == cv$lambda.1se])
      cvms.mins <- c(cv$cvm[cv$lambda == cv$lambda.min])
      lambda.1ses <- c(cv$lambda.1se)
      lambda.mins <- c(cv$lambda.min)
      alphas <- c(l)
    }else{
      cvms.1ses <- c(cvms.1ses, cv$cvm[cv$lambda == cv$lambda.1se])
      cvms.mins <- c(cvms.mins, cv$cvm[cv$lambda == cv$lambda.min])
      lambda.1ses <- c(lambda.1ses, cv$lambda.1se)
      lambda.mins <- c(lambda.mins, cv$lambda.min)
      alphas <- c(alphas, l)
    }
    
  }
  
  search <- data.frame(cvms.1ses = cvms.1ses, cvms.mins = cvms.mins, 
                       lambda.1ses = lambda.1ses, lambda.mins = lambda.mins, 
                       alphas = alphas)
  
  if(errortype == 'min'){
    parameters <- search[search$cvms.mins == min(search$cvms.mins),]
    
    if(verbose == TRUE){
      
      cat(paste0('Choose the elastic net mixing parameter (alpha) as ', 
                 parameters$alphas, '\n'))
      cat(paste0('Choose the regularization constant (lambda) as ', 
                 signif(parameters$lambda.mins, 3), '\n'))
      
    }
    
    
    
    Alpha <- parameters$alphas
    Lambda <- parameters$lambda.mins
    
    parameters.backup <- search[search$cvms.1ses == min(search$cvms.1ses),]
    
    Alpha.backup <- parameters.backup$alphas
    Lambda.backup <- parameters.backup$lambda.1ses
    
  }else{
    parameters <- search[search$cvms.1ses == min(search$cvms.1ses),]
    
    if(verbose == TRUE){
      
      cat(paste0('Choose the elastic net mixing parameter (alpha) as ', 
                 parameters$alphas, '\n'))
      cat(paste0('Choose the regularization constant (lambda) as ', 
                 signif(parameters$lambda.1ses, 3), '\n'))
      
    }
    
    Alpha <- parameters$alphas
    Lambda <- parameters$lambda.1ses
    
    Alpha.backup <- Alpha
    Lambda.backup <- Lambda
    
  }
  
  #Elastic net##########
  
  
  
  elasticnetmodel <- tryCatch({
    
    glmnet::glmnet(x = trainvar, 
                   y = trainres$Response, 
                   family = "multinomial", 
                   lambda = Lambda, 
                   alpha = Alpha)
    
  }, error = function(err){
    
    NULL
    
  }, warning = function(war){
    
    NULL
    
  })
  
  if(is.null(elasticnetmodel)){
    
    cat(paste0('Rechoose the elastic net mixing parameter (alpha) as ', 
               Alpha.backup, '\n'))
    cat(paste0('Rechoose the regularization constant (lambda) as ', 
               signif(Lambda.backup, 3), '\n'))
    
    
    elasticnetmodel <- tryCatch({
      
      glmnet::glmnet(x = trainvar, 
                     y = trainres$Response, 
                     family = "multinomial", 
                     lambda = Lambda.backup, 
                     alpha = Alpha.backup)
      
    }, error = function(err){
      
      NULL
      
    }, warning = function(war){
      
      NULL
      
    })
    
    if(is.null(elasticnetmodel)){
      
      return(NULL)
      
    }
    
    
  }
  
  
  modelcoefs <- coef(elasticnetmodel)
  
  modelcoef <- do.call(cbind, modelcoefs)
  modelcoef <- as.matrix(modelcoef)
  
  #modelcoef <- modelcoef[rowSums(modelcoef) != 0, , drop = FALSE]
  
  colnames(modelcoef) <- paste0(colnames(modelcoef), '_', names(modelcoefs))
  
  reslist <- list(elasticnetmodel = elasticnetmodel, modelcoefs = modelcoef)
  
  return(reslist)
  
}

coefensemble <- function(modelcoeflist){
  
  multifinalfeatures <- list()
  i <- 1
  for(i in 1:ncol(modelcoeflist[[1]])){
    
    modelcoefs <- do.call(cbind, lapply(X = modelcoeflist, 
                                        FUN = function(x){x[, i, drop = FALSE]}))
    
    coefmean <- rowMeans(modelcoefs)
    coefplus <- rowSums(modelcoefs > 0)
    coefminus <- rowSums(modelcoefs < 0)
    coefsign <- coefplus - coefminus
    coefscore <- coefsign*coefmean
    coefscoresd <- sd(coefscore)
    coefscoremean <- mean(coefscore)
    coefscoreupper <- coefscoremean + 1.5*coefscoresd
    coefscorelower <- coefscoremean - 1.5*coefscoresd
    
    finalfeatures <- coefscore[coefscore < coefscorelower | coefscore > coefscoreupper]
    
    finalfeatures <- as.matrix(finalfeatures)
    finalfeatures <- finalfeatures[-1,]
    
    finalfeatures <- finalfeatures[order(-finalfeatures)]
    
    
    
    if(length(finalfeatures) == 0){
      
      finalfeatures <- coefscore[-1]
      
    }
    
    
    
    multifinalfeatures[[i]] <- finalfeatures
    
  }
  
  return(multifinalfeatures)
  
}



singleclassifier <- function(finalfeatures, 
                             resamplevar, 
                             resampleresponse, 
                             seednum = 2022, 
                             
                             method = 'SVM', 
                             C.base = 10, 
                             C.min = -3, 
                             C.max = -2){
  
  #Model training on resampled data
  subresamplevar <- resamplevar[,finalfeatures, drop = FALSE]
  
  subresamplevar <- cbind(subresamplevar, resampleresponse$Response)
  colnames(subresamplevar)[ncol(subresamplevar)] <- 'Response'
  subresamplevar <- as.data.frame(subresamplevar)
  
  
  
  if(method == 'SVM'){
    
    foldnum <- max(min(floor(nrow(subresamplevar)/10), 10), 3)
    
    suppressMessages(
      
      svm.linearcv <- train_SVM_e1071_LK(y = subresamplevar$Response,
                                         betas.Train = subresamplevar[, 1:ncol(subresamplevar)-1, drop = FALSE],
                                         seed = seednum,
                                         nfolds = foldnum,
                                         mc.cores = 1,
                                         C.base = C.base,
                                         C.min = C.min,
                                         C.max = C.max, 
                                         weighted = FALSE)
      
    )
    
    fit <- svm.linearcv[[1]]
    
  }else{
    
    #fit <- nnet::multinom(formula = Response~., data = subresamplevar, trace = FALSE)
    
    suppressMessages(
      
      fit <- tryCatch({
        
        glmnet::glmnet(x = subresamplevar[, 1:ncol(subresamplevar)-1, drop = FALSE], 
                       y = subresamplevar$Response, family = "multinomial", 
                       lambda = 0)
        
      }, error = function(err){
        
        NULL
        
      }, warning = function(war){
        
        NULL
        
      })
      
    )
    
    if(!is.null(fit)){
      
      if(is.na(fit$dev.ratio)){
        fit <- NULL
      }
      
      if(sum(fit$dfmat != fit$df) > 0){
        
        fit <- NULL
        
      }
    }
    
    if(is.null(fit)){
      
      foldnum <- max(min(floor(nrow(subresamplevar)/10), 10), 3)
      
      set.seed(seednum)
      cv <- glmnet::cv.glmnet(x = as.matrix(subresamplevar[, 1:ncol(subresamplevar)-1, 
                                                           drop = FALSE]), 
                              y = subresamplevar$Response, 
                              family = "multinomial", 
                              nfold = foldnum, 
                              type.measure = "deviance", 
                              paralle = FALSE, 
                              alpha = 0)
      
      cvms.1ses <- c(cv$cvm[cv$lambda == cv$lambda.1se])
      cvms.mins <- c(cv$cvm[cv$lambda == cv$lambda.min])
      lambda.1ses <- c(cv$lambda.1se)
      lambda.mins <- c(cv$lambda.min)
      alphas <- 0
      
      search <- data.frame(cvms.1ses = cvms.1ses, cvms.mins = cvms.mins, 
                           lambda.1ses = lambda.1ses, lambda.mins = lambda.mins, 
                           alphas = alphas)
      
      parameters <- search[search$cvms.mins == min(search$cvms.mins),]
      
      Alpha <- parameters$alphas
      Lambda <- parameters$lambda.mins
      
      
      fit <- glmnet::glmnet(x = subresamplevar[, 1:ncol(subresamplevar)-1, drop = FALSE], 
                            y = subresamplevar$Response, 
                            family = "multinomial", 
                            alpha = Alpha, 
                            lambda = Lambda)
      
    }
    
  }
  
  return(fit)
  
}



modensemble <- function(orivar, 
                        orires, 
                        modellist){
  
  renameresponse <- FALSE
  
  method <- 'glmnet'
  
  if(!is.factor(orires) & !is.numeric(orires)){
    
    responselevels <- orires
    responselevels <- unique(responselevels)
    responselevels <- responselevels[order(responselevels)]
    
    orires <- factor(orires, 
                     levels = responselevels, 
                     ordered = TRUE)
    
  }
  
  if(!is.numeric(orires)){
    
    responselevels <- levels(orires)
    orires <- as.numeric(orires)
    renameresponse <- TRUE
    
  }
  
  
  
  featurenames <- row.names(modellist[[1]]$beta[[1]])
  
  if(length(featurenames) == 0){
    
    method <- 'SVM'
    
    featurenames <- colnames(modellist[[1]]$SV)
    
  }
  
  featurenames <- gsub(pattern = "`", replacement = '', x = featurenames)
  
  subtrain <- orivar[,featurenames, drop = FALSE]
  
  subtrain <- as.data.frame(subtrain, stringsAsFactors = FALSE)
  
  trainproblist <- list()
  trainprelist <- list()
  m <- 1
  for(m in 1:length(modellist)){
    subfit <- modellist[[m]]
    
    if(method == 'SVM'){
      
      svmres <- e1071:::predict.svm(object = subfit,
                                    newdata =  as.matrix(subtrain),
                                    decision.values = TRUE,
                                    probability = TRUE)
      
      trainpre <- as.character(svmres)
      trainpre <- matrix(trainpre, ncol = 1)
      
      trainprob <- attr(svmres, 'probabilities')
      trainprob <- trainprob[, levels(svmres), drop = FALSE]
      
      trainproblist[[m]] <- trainprob
      
      trainprelist[[m]] <- trainpre
      
      
    }else{
      
      trainprob <- predict(subfit, as.matrix(subtrain), type = 'response')
      trainproblist[[m]] <- trainprob
      
      trainpre <- predict(subfit, as.matrix(subtrain), type = 'class')
      trainprelist[[m]] <- trainpre
      
    }
    
  }
  
  
  if(length(dim(trainproblist[[1]])) == 3){
    traintrueprobs <- trainproblist[[1]][,,1]
  }else{
    traintrueprobs <- trainproblist[[1]]
  }
  
  traintrueprobs[traintrueprobs != 0] <- 0
  idx <- match(orires, seq(1, ncol(traintrueprobs), 1))
  traintrueprobs[cbind(1:nrow(traintrueprobs), idx)] <- 1
  
  modeleuclosses <- c()
  
  for(n in 1:length(trainproblist)){
    
    if(length(dim(trainproblist[[n]])) == 3){
      eucloss <- lapply(X = seq(1, nrow(traintrueprobs), 1), FUN = function(x){
        dist(rbind(traintrueprobs[x,], trainproblist[[n]][x,,1]), 
             method = 'euclidean')^2
      })
      
      eucloss <- mean(unlist(eucloss))
      
      #eucloss <- mean(diag(proxy::dist(traintrueprobs, trainproblist[[n]][,,1], 
      #                                 method = 'euclidean'))^2)
      
      
    }else{
      eucloss <- lapply(X = seq(1, nrow(traintrueprobs), 1), FUN = function(x){
        dist(rbind(traintrueprobs[x,], trainproblist[[n]][x,]), 
             method = 'euclidean')^2
      })
      
      eucloss <- mean(unlist(eucloss))
      
      #eucloss <- mean(diag(proxy::dist(traintrueprobs, trainproblist[[n]], 
      #                                 method = 'euclidean'))^2)
      
      
    }
    
    modeleuclosses[n] <- eucloss
    
    
  }
  
  modeleuclosses <- as.matrix(modeleuclosses)
  colnames(modeleuclosses) <- c('training')
  row.names(modeleuclosses) <- paste0('model_', 1:nrow(modeleuclosses))
  
  
  
  
  trainpreres <- do.call(cbind, trainprelist)
  
  colnames(trainpreres) <- paste0('model_', 1:ncol(trainpreres))
  
  modelaccs <- colSums(trainpreres == orires)/length(orires)
  modelaccs <- as.matrix(modelaccs)
  colnames(modelaccs) <- c('training')
  
  #Only keep base learners with an R^2 > 0.5 in training data
  acc <- 0.5
  
  if(sum(modelaccs[,1] > acc) == 0){
    return(NULL)
  }
  
  trainproblist <- trainproblist[modelaccs[,1] > acc]
  
  
  modellist <- modellist[modelaccs[,1] > acc]
  trainpreres <- trainpreres[, modelaccs[,1] > acc, drop = FALSE]
  
  modelaccs <- modelaccs[modelaccs[,1] > acc, , drop = FALSE]
  modeleuclosses <- modeleuclosses[modelaccs[,1] > acc, , drop = FALSE]
  
  
  if(!is.matrix(trainpreres)){
    
    trainpreres <- as.matrix(trainpreres)
    colnames(trainpreres) <- 'model_1'
    
    modelaccs <- t(as.matrix(modelaccs))
    rownames(modelaccs) <- 'model_1'
    
    modeleuclosses <- t(as.matrix(modeleuclosses))
    rownames(modeleuclosses) <- 'model_1'
    
    
  }
  
  #weights <- 1/modeleuclosses[,1]
  
  #Calcuate the weight for each base learner only from training data
  weights <- 0.5 * log(modelaccs[,1]/(1 - modelaccs[,1]))
  
  if(any(is.infinite(weights))){
    if(any(is.finite(weights))){
      
      eucweights <- 1/modeleuclosses[,1]
      
      finiteaccmax <- max(weights[is.finite(weights)])
      
      finiteeucmax <- min(eucweights[weights == finiteaccmax])
      
      infiniteeucmax <- max(eucweights[is.infinite(weights)])
      
      fac <- infiniteeucmax/finiteeucmax
      
      fac <- max(length(unique(orires)) - 1, fac)
      
      weights[is.infinite(weights)] <- fac*max(weights[is.finite(weights)])
      
    }else{
      
      weights <- rep(1/length(weights), length(weights))
      names(weights) <- row.names(modelaccs)
      
    }
  }
  
  weights <- weights/sum(weights)
  weights <- as.matrix(weights)
  weights <- t(weights)
  
  
  #Ensemble the prediction from all base learners
  
  trainprobres <- Map(`*`, trainproblist, t(weights))
  trainprobres <- Reduce(`+`, trainprobres)
  trainprobres <- trainprobres/rowSums(trainprobres)
  
  trainpreres <- colnames(trainprobres)[apply(trainprobres, 1, which.max)]
  
  if(length(trainpreres) == 0){
    return(NULL)
  }
  
  rownames(weights) <- 'weights'
  
  if(is.null(colnames(weights))){
    colnames(weights) <- 'model_1'
  }
  
  trainacc <- sum(trainpreres == orires)/length(orires)
  
  
  traincomp <- cbind(trainpreres, orires)
  colnames(traincomp) <- c('Prediction', 'True')
  row.names(traincomp) <- row.names(subtrain)
  traincomp <- as.data.frame(traincomp, stringsAsFactors = FALSE)
  
  if(renameresponse == TRUE){
    
    colnames(trainprobres) <- responselevels
    
    traincomp[,1] <- responselevels[as.integer(traincomp[,1])]
    traincomp[,2] <- responselevels[as.integer(traincomp[,2])]
    
    traincomp[,1] <- factor(traincomp[,1], levels = responselevels)
    traincomp[,2] <- factor(traincomp[,2], levels = responselevels)
    
  }
  
  if(method != 'SVM'){
    
    trainprobres <- trainprobres[, , 1]
    
  }
  
  reslist <- list(baselearners = modellist, 
                  baseweights = weights, 
                  traincomp = traincomp, 
                  trainprobs = trainprobres)
  
  if(renameresponse == TRUE){
    reslist$responselevels <- responselevels
  }
  
  
  return(reslist)
  
}



singlebalancemod <- function(trainingdat, 
                             trainingpd, 
                             balanceadj = 1, 
                             samplingmethod = 'updn', 
                             nround = 10, 
                             k = 5, 
                             threads = 1, 
                             
                             method = 'SVM', 
                             C.base = 10, 
                             C.min = -3, 
                             C.max = -2, 
                             
                             alphas = c(0.5), 
                             seednum = 2022, 
                             errortype = 'min', 
                             
                             plot = FALSE, 
                             prefix, 
                             textsize = 13, 
                             verbose = TRUE){
  
  if(balanceadj == 1){
    
    balanceddats <- balancesampling(dat = trainingdat, 
                                    pddat = trainingpd, 
                                    responsename = 'Response', 
                                    nround = nround, 
                                    seed = seednum, 
                                    k = k, 
                                    threads = threads, 
                                    samplingmethod = samplingmethod)
    
    
  }else if(balanceadj == 2){
    
    balanceddats <- bootsampling(dat = trainingdat, 
                                 pddat = trainingpd, 
                                 responsename = 'Response', 
                                 nround = nround, 
                                 seed = seednum, 
                                 threads = threads)
    
  }else{
    
    balanceddats <- list(dats = list(round_1 = trainingdat), 
                         pds = list(round_1 = trainingpd))
    
    
    
  }
  
  
  probelist <- list()
  
  nullidx <- c()
  
  j <- 1
  
  for(j in 1:length(balanceddats$dats)){
    
    cat(paste0('j = ', j, '\n'))
    
    foldnum <- max(min(floor(ncol(balanceddats$dats[[j]])/10), 10), 3)
    
    elasticres <- singleclassifiernet(trainvar = t(balanceddats$dats[[j]]), 
                                      trainres = balanceddats$pds[[j]], 
                                      threads = threads, 
                                      alphas = alphas, 
                                      seednum = seednum + j, 
                                      foldnum = foldnum, 
                                      errortype = errortype, 
                                      verbose = verbose)
    
    if(is.null(elasticres)){
      probelist[[j]] <- NULL
      nullidx <- c(nullidx, j)
      
    }else{
      probelist[[j]] <- elasticres$modelcoefs
    }
    
    
  }
  
  if(length(nullidx) > 0){
    for(p in nullidx){
      balanceddats$dats[[p]] <- NULL
      balanceddats$pds[[p]] <- NULL
    }
  }
  
  probelist <- probelist[lengths(probelist) != 0]
  balanceddats$dats <- balanceddats$dats[lengths(balanceddats$dats) != 0]
  balanceddats$pds <- balanceddats$pds[lengths(balanceddats$pds) != 0]
  
  
  #Part4, feature integration
  multifinalprobes <- coefensemble(modelcoeflist = probelist)
  finalprobes <- unique(unlist(lapply(X = multifinalprobes, FUN = function(x){names(x)})))
  
  #Part5, multinomial base learner training
  modellist <- list()
  m <- 1
  for(m in 1:length(balanceddats$dats)){
    
    cat(paste0('m = ', m, '\n'))
    
    balancedsamples <- list(resamplevars = t(balanceddats$dats[[m]]), 
                            resampleresponse = balanceddats$pds[[m]])
    
    
    classfit <- singleclassifier(finalfeatures = finalprobes, 
                                 resamplevar = balancedsamples$resamplevars, 
                                 resampleresponse = balancedsamples$resampleresponse, 
                                 
                                 seednum = seednum, 
                                 
                                 method = method, 
                                 C.base = C.base, 
                                 C.min = C.min, 
                                 C.max = C.max)
    
    modellist[[m]] <- classfit
    
  }
  
  #Part6, final ensemble model
  resampleensemble <- modensemble(orivar = t(trainingdat), 
                                  orires = trainingpd$Response, 
                                  modellist = modellist)
  
  if(is.null(resampleensemble)){
    return(NULL)
  }
  
  
  modelscores <- lapply(X = multifinalprobes, 
                        FUN = function(x){x <- as.matrix(x); colnames(x)[1] <- 'score'; return(x)})
  
  
  names(modelscores) <- paste0('model_', seq(1, length(modelscores), 1))
  
  resampleensemble$modelscores <- modelscores
  
  if(plot == TRUE){
    
    heatmapplotting(orivar = t(trainingdat), 
                    comptab = resampleensemble$traincomp, 
                    featurenames = row.names(resampleensemble$modelscores[[1]]), 
                    textsize = textsize, 
                    
                    title = paste0(prefix, 
                                   c('Training from', 
                                     ' training from')[c(is.null(prefix), !is.null(prefix))])
    )
    
    
  }
  
  return(resampleensemble)
  
}



pairedbalancemod <- function(trainingdats, 
                             trainingpd, 
                             balanceadj = 1, 
                             samplingmethod = 'updn', 
                             nround = 10, 
                             k = 5, 
                             threads = 1, 
                             
                             method = 'SVM', 
                             C.base = 10, 
                             C.min = -3, 
                             C.max = -2, 
                             
                             alphas = c(0.5), 
                             seednum = 2022, 
                             errortype = 'min', 
                             
                             plot = FALSE, 
                             prefixes, 
                             textsize = 13, 
                             verbose = TRUE){
  
  datsensemble <- list()
  
  cat(paste0('Dat = 1\n'))
  datsensemble[[1]] <- singlebalancemod(trainingdat = trainingdats[[1]], 
                                        trainingpd = trainingpd, 
                                        balanceadj = balanceadj, 
                                        samplingmethod = samplingmethod, 
                                        nround = nround, 
                                        k = k, 
                                        threads = threads, 
                                        
                                        method = method, 
                                        C.base = C.base, 
                                        C.min = C.min, 
                                        C.max = C.max, 
                                        
                                        alphas = alphas, 
                                        seednum = seednum, 
                                        errortype = errortype, 
                                        
                                        plot = plot, 
                                        prefix = prefixes[1], 
                                        textsize = textsize, 
                                        verbose = verbose)
  
  
  
  if(length(trainingdats) == 1){
    
    weights <- 1
    weights <- as.matrix(weights)
    weights <- t(weights)
    colnames(weights) <- 'model_1'
    row.names(weights) <- 'weights'
    
    traincomp = datsensemble[[1]]$traincomp
    trainprobres = datsensemble[[1]]$trainprobs
    
    
    if('responselevels' %in% names(datsensemble[[1]])){
      
      if(!is.null(datsensemble[[1]]$responselevels)){
        
        traincomp$Prediction <- factor(traincomp$Prediction, 
                                       levels = datsensemble[[1]]$responselevels)
        
      }else{
        traincomp$Prediction <- as.numeric(traincomp$Prediction)
      }
      
    }else{
      traincomp$Prediction <- as.numeric(traincomp$Prediction)
    }
    
  }else{
    
    datseucloss <- c()
    
    traintrueprobs <- datsensemble[[1]]$trainprobs
    traintrueprobs[traintrueprobs != 0] <- 0
    idx <- match(datsensemble[[1]]$traincomp[,2], colnames(traintrueprobs))
    traintrueprobs[cbind(1:nrow(traintrueprobs), idx)] <- 1
    
    eucloss <- lapply(X = seq(1, nrow(traintrueprobs), 1), FUN = function(x){
      dist(rbind(traintrueprobs[x,], datsensemble[[1]]$trainprobs[x,]), 
           method = 'euclidean')^2
    })
    
    eucloss <- mean(unlist(eucloss))
    
    #eucloss <- mean(diag(proxy::dist(traintrueprobs, datsensemble[[1]]$trainprobs, 
    #                                 method = 'euclidean'))^2)
    
    datseucloss[1] <- eucloss
    
    
    
    datsacc <- c()
    
    datsacc[1] <- sum(datsensemble[[1]]$traincomp[,1] == datsensemble[[1]]$traincomp[,2])/
      nrow(datsensemble[[1]]$traincomp)
    
    trainprobses <- list()
    
    trainprobses[[1]] <- datsensemble[[1]]$trainprobs
    
    for(i in 2:length(trainingdats)){
      
      cat(paste0('Dat = ', i, '\n'))
      
      singlebalancemodres <- singlebalancemod(trainingdat = trainingdats[[i]], 
                                              trainingpd = trainingpd, 
                                              balanceadj = balanceadj, 
                                              samplingmethod = samplingmethod, 
                                              nround = nround, 
                                              k = k, 
                                              threads = threads, 
                                              
                                              
                                              method = method, 
                                              C.base = C.base, 
                                              C.min = C.min, 
                                              C.max = C.max, 
                                              
                                              
                                              alphas = alphas, 
                                              seednum = seednum, 
                                              errortype = errortype, 
                                              
                                              plot = plot, 
                                              prefix = prefixes[i], 
                                              textsize = textsize, 
                                              verbose = verbose)
      
      datsensemble[[i]] <- singlebalancemodres
      
      if((is.null(singlebalancemodres)) & (i != length(trainingdats))){
        next()
      }
      
      if(!is.null(singlebalancemodres)){
        
        traintrueprobs <- datsensemble[[i]]$trainprobs
        traintrueprobs[traintrueprobs != 0] <- 0
        idx <- match(datsensemble[[i]]$traincomp[,2], colnames(traintrueprobs))
        traintrueprobs[cbind(1:nrow(traintrueprobs), idx)] <- 1
        
        eucloss <- lapply(X = seq(1, nrow(traintrueprobs), 1), FUN = function(x){
          dist(rbind(traintrueprobs[x,], datsensemble[[i]]$trainprobs[x,]), 
               method = 'euclidean')^2
        })
        
        eucloss <- mean(unlist(eucloss))
        
        #eucloss <- mean(diag(proxy::dist(traintrueprobs, datsensemble[[i]]$trainprobs, 
        #                                 method = 'euclidean'))^2)
        
        datseucloss[i] <- eucloss
        
        
        datsacc[i] <- sum(datsensemble[[i]]$traincomp[,1] == datsensemble[[i]]$traincomp[,2])/
          nrow(datsensemble[[i]]$traincomp)
        
        trainprobses[[i]] <- datsensemble[[i]]$trainprobs
        
      }
      
      if(i == length(trainingdats)){
        
        modeleuclosses <- matrix(datseucloss, ncol = 1)
        row.names(modeleuclosses) <- paste0('model_', seq(1, length(trainingdats), 1))
        
        modelaccs <- matrix(datsacc, ncol = 1)
        row.names(modelaccs) <- paste0('model_', seq(1, length(trainingdats), 1))
        
        #weights <- 1/modeleuclosses[,1]
        
        #Calcuate the weight for each base learner only from training data
        weights <- 0.5 * log(modelaccs[,1]/(1 - modelaccs[,1]))
        
        if(any(is.infinite(weights))){
          if(any(is.finite(weights))){
            
            eucweights <- 1/modeleuclosses[,1]
            
            finiteaccmax <- max(weights[is.finite(weights)])
            
            finiteeucmax <- min(eucweights[weights == finiteaccmax])
            
            infiniteeucmax <- max(eucweights[is.infinite(weights)])
            
            fac <- infiniteeucmax/finiteeucmax
            
            fac <- max(length(unique(trainingpd$Response)) - 1, fac)
            
            weights[is.infinite(weights)] <- fac*max(weights[is.finite(weights)])
            
          }else{
            
            weights <- rep(1/length(weights), length(weights))
            names(weights) <- row.names(modelaccs)
            
          }
        }
        
        weights <- weights/sum(weights)
        weights <- as.matrix(weights)
        weights <- t(weights)
        
        #Ensemble the prediction from all base learners
        
        trainprobres <- Map(`*`, trainprobses, t(weights))
        
        trainprobres <- Reduce(`+`, trainprobres)
        trainprobres <- trainprobres/rowSums(trainprobres)
        
        trainpreres <- colnames(trainprobres)[apply(trainprobres, 1, which.max)]
        
        rownames(weights) <- 'weights'
        
        
        trainacc <- sum(trainpreres == trainingpd$Response)/length(trainingpd$Response)
        
        traincomp <- data.frame(Prediction = trainpreres, 
                                True = trainingpd$Response, 
                                stringsAsFactors = FALSE)
        row.names(traincomp) <- trainingpd$sampleid
        
        if('responselevels' %in% names(datsensemble[[1]])){
          
          if(!is.null(datsensemble[[1]]$responselevels)){
            
            traincomp$Prediction <- factor(traincomp$Prediction, 
                                           levels = datsensemble[[1]]$responselevels)
            
          }else{
            traincomp$Prediction <- as.numeric(traincomp$Prediction)
          }
          
        }else{
          traincomp$Prediction <- as.numeric(traincomp$Prediction)
        }
        
      }
      
      
    }
    
  }
  
  finalensemble <- list(baselearners = datsensemble, 
                        baseweights = weights, 
                        traincomp = traincomp, 
                        trainprobs = trainprobres, 
                        responselevels = datsensemble[[1]]$responselevels)
  
  return(finalensemble)
  
}



singleensemblepredict <- function(balanceres, 
                                  vardat){
  
  method <- 'glmnet'
  
  if('SV' %in% names(balanceres$baselearners[[1]])){
    
    method <- 'SVM'
    
  }
  
  num <- length(balanceres$baselearners)
  
  baselearnerprobs <- list()
  o <- 1
  for(o in 1:num){
    
    baselearner <- balanceres$baselearners[[o]]
    
    if(method == 'SVM'){
      
      featurenames <- colnames(baselearner$SV)
      
    }else{
      
      featurenames <- row.names(baselearner$beta[[1]])
      
    }
    
    vardat <- vardat[, featurenames, drop = FALSE]
    
    if(method == 'SVM'){
      
      baselearnerprob <- e1071:::predict.svm(object = baselearner,
                                              newdata =  as.matrix(vardat),
                                              decision.values = FALSE,
                                              probability = TRUE)
      
      baselearnerprob <- attr(baselearnerprob, 
                              'probabilities')[, levels(baselearnerprob), 
                                               drop = FALSE]
      
      
    }else{
      
      baselearnerprob <- predict(baselearner, vardat, type = 'response')[, , 1]
      
    }
    
    
    
    
    
    if('responselevels' %in% names(balanceres)){
      
      if(!is.null(balanceres$responselevels)){
        
        colnames(baselearnerprob) <- balanceres$responselevels
        
      }
      
      
      
    }
    
    
    baselearnerprobs[[o]] <- baselearnerprob*balanceres$baseweights[1, o]
    
  }
  
  
  ensembleprobs <- Reduce(`+`, baselearnerprobs)
  
  ensemblepres <- colnames(ensembleprobs)[apply(X = ensembleprobs, MARGIN = 1, 
                                                FUN = which.max)]
  ensemblepres <- data.frame(Prediction = ensemblepres, stringsAsFactors = FALSE)
  
  if('responselevels' %in% names(balanceres)){
    
    if(!is.null(balanceres$responselevels)){
      
      ensemblepres$Prediction <- factor(as.character(ensemblepres$Prediction), 
                                        levels = balanceres$responselevels)
      
    }else{
      ensemblepres$Prediction <- as.numeric(ensemblepres$Prediction)
    }
    
    
  }else{
    ensemblepres$Prediction <- as.numeric(ensemblepres$Prediction)
  }
  
  
  row.names(ensemblepres) <- row.names(baselearnerprob)
  
  return(list(pres = ensemblepres, probs = ensembleprobs))
  
}

#'Predict single-omic or multi-omics sample class
#'
#'Predict single-omic or multi-omics sample class with the model trained from 
#'  the function \code{omicsclassifier}. 
#'  
#'@param balanceres The model trained by the function \code{omicsclassifier}, 
#'  which is the slot "mod" of the returned list result of the function.
#'@param dats A list with each element as each omic data of the samples to be 
#'  classified. Each omic data should be a matrix with columns as features and 
#'  rows as samples. The column names are the feature names and the row names 
#'  are the sample names. The number of omics here should be the same as the 
#'  data used to train the model with the function \code{omicsclassifier}. If 
#'  there is only one omic, it need to be a list with one element as the omic 
#'  data matrix. 
#'@return A list with 2 slots. One named "pres" is a single-column data frame 
#'  with the predicted class labels of the samples. The other named "probs" is 
#'  a matrix of the predicted class probability distribution of the samples. 
#'@export
pairedensemblepredict <- function(balanceres, 
                                  vardats){
  
  ensemblereses <- list()
  
  ensemblereses[[1]] <- singleensemblepredict(balanceres = balanceres$baselearners[[1]], 
                                              vardat = vardats[[1]])
  
  if(length(vardats) == 1){
    
    ensembleprobs <- ensemblereses[[1]]$probs
    
    ensemblepres <- ensemblereses[[1]]$pres
    
    if('responselevels' %in% names(balanceres)){
      
      if(!is.null(balanceres$responselevels)){
        ensemblepres$Prediction <- factor(ensemblepres$Prediction, 
                                          levels = balanceres$responselevels)
      }else{
        ensemblepres$Prediction <- as.numeric(ensemblepres$Prediction)
      }
      
    }else{
      
      ensemblepres$Prediction <- as.numeric(ensemblepres$Prediction)
      
    }
    
    
    
  }else{
    
    ensembleprobs <- ensemblereses[[1]]$probs * balanceres$baseweights[1,1]
    
    for(i in 2:length(balanceres$baselearners)){
      
      ensemblereses[[i]] <- singleensemblepredict(balanceres = balanceres$baselearners[[i]], 
                                                  vardat = vardats[[i]])
      
      
      ensembleprobs <- ensembleprobs + ensemblereses[[i]]$probs * balanceres$baseweights[1,i]
      
      if(i == length(balanceres$baselearners)){
        
        ensemblepres <- colnames(ensembleprobs)[apply(X = ensembleprobs, 
                                                      MARGIN = 1, 
                                                      FUN = which.max)]
        
        ensemblepres <- data.frame(Prediction = ensemblepres, stringsAsFactors = FALSE)
        
        if('responselevels' %in% names(balanceres)){
          
          if(!is.null(balanceres$responselevels)){
            ensemblepres$Prediction <- factor(ensemblepres$Prediction, 
                                              levels = balanceres$responselevels)
          }else{
            ensemblepres$Prediction <- as.numeric(ensemblepres$Prediction)
          }
          
        }else{
          
          ensemblepres$Prediction <- as.numeric(ensemblepres$Prediction)
          
        }
        
      }
      
    }
    
  }
  
  colnames(ensemblepres) <- 'Prediction'
  row.names(ensemblepres) <- row.names(ensemblereses[[1]]$pres)
  
  res <- list()
  
  res$pres <- ensemblepres
  res$probs <- ensembleprobs
  
  return(res)
  
}

#'Perform single-omic or multi-omics classification
#'
#'Train single-omic or multi-omics classifiers via several steps, including 
#'  a data distribution adjustment step (optional), a elastic net feature 
#'  selection step, and an SVM or MLR model training step.
#'
#'@param dats A list with each element as each omic data of the samples to be 
#'  classified. Each omic data should be a matrix with columns as samples and 
#'  rows as features. The column names are the sample names and the row names 
#'  are the feature names.If only one omic data need to be classified, it need 
#'  to be a list with one element as the omic data matrix. 
#'@param truelabels A vector with the sample class labels as elements. It can 
#'  use the sample names as the element names of the vector. If it does not 
#'  have element names, its element order should be consistent with the sample 
#'  order in all the omics data of \code{dats}.
#'@param varfeatures For each omic data in \code{dats}, the standard deviation 
#'  of the features can be calculated and the top variable features will be 
#'  selected and used for clustering. This parameter defines the number of the 
#'  top variable features will be selected. It should be an integer, such as 
#'  10000, so that the top 10000 most variable features of each omic will be 
#'  selected, and if the omic data have a feature number less than 10000, all 
#'  the features will be used. Can also be NULL, so that this variable feature 
#'  selection step will be skipped and all the omics will use all the features 
#'  for the classification. Default is NULL.
#'@param scale Whether the features of the omics need to be standardized, so 
#'  that all of them have a mean across all the samples as 0 and a standard 
#'  deviation of 1, and the classification will be performed after that. Its 
#'  default is TRUE.
#'@param balanceadj Before the elastic net feature selection and the SVM/MLR 
#'  model training steps, whether to construct an ensemble framework or not. 
#'  When this parameter is set as 1, a bagging-SMOTE (synthetic minority over- 
#'  sampling technique) ensemble framework will be constructed. In this case, 
#'  the sample distribution will be adjusted with bagging and SMOTE sampling 
#'  so that all the sample classes will have the same size, which prevents the 
#'  downstream model training from being influenced by any class sample size 
#'  imbalance. The model will have better accuracy in rare sample prediction. 
#'  When this parameter is set as 2, a normal bagging framework, rather than 
#'  the bagging-SMOTE one, will be constructed. The bagging framework is also 
#'  an ensemble one but it will not adjust the sample distribution. When this 
#'  parameter is set as 3, no ensemble framework will be constructed. In this 
#'  case, the ensemble step will be skipped and only the downstream elastic 
#'  net feature selection and SVM/MLR model training steps will be performed. 
#'  If the dataset is a single-omic one, the final prediction will be returned 
#'  after these steps. However, if the dataset is a multi-omics one, these 
#'  steps will be performed on each omic, respectively, and after all of them 
#'  return the prediction results for all the omics, a further ensemble step 
#'  will be used to aggregate them into the final one. The default value of 
#'  this paramter is 1.
#'@param nround When \code{balanceadj} is 1 or 2, this parameter will be used 
#'  to determine the number of the base learner sets in the ensemble-based 
#'  mode. Default is 10.
#'@param k When \code{balanceadj} is set as 1, this parameter will be used for 
#'  the SMOTE up-sampling step of the bagging-SMOTE framework, which works on 
#'  small sample classes and synthesizes new samples to increase their sample 
#'  number to the level of original total sample number/group number. It does 
#'  this from the closet neighbors of a randomly selected sample in the same 
#'  group, this parameter is used to define how many closet neighbors will be 
#'  needed. Default is 5, meaning the top 5 closet neighbors will be included 
#'  in synthesizing the new pseudo-sample.
#'@param nfold This parameter can be an integer and if it is greater than 1, 
#'  a cross validation will be performed and if the model achieves an accuracy 
#'  greater than the cutoff defined by the parameter \code{crosstestacccut}, 
#'  the final model will be constructed from all the data transferred to the 
#'  parameter \code{dats}. If this parameter is the integer 1 or NULL, the 
#'  cross validation will be skipped and a model will be directly constructed 
#'  from the whole data. Default is 5, so a 5-fold cross validation will be 
#'  performed. 
#'@param seednum The random seed number used for the elastic net regression to 
#'  define a 10-fold cross-validation, so an optimal regularization constant 
#'  lambda can be selected for the elastic net. It will also be used for the 
#'  balanced base learner sampling process and the cross validation dataset 
#'  generation process. Default is 2022.
#'@param threads Threads number for parallelization, default is 1.
#'@param prescreen When this parameter is set as TRUE, before starting the 
#'  formal training, a preliminary screen on the features will be performed 
#'  first with ANOVA so that only the features with an ANOVA p-val less than 
#'  0.05 on the response classes will be selected for the model training. Its 
#'  default value is TRUE.
#'@param method In the SVM/MLR model training step, which algorithm should be 
#'  used. Its value can be "SVM" or "MLR". Default is "SVM".
#'@param C.base A parameter special for the SVM model training step to set its 
#'  regularization constant. This constant will be calculated by the function 
#'  as base^index, and \code{C.base} here serves as the base number. Combined 
#'  with other 2 parameters \code{C.min} and \code{C.max} serving as indexes, 
#'  it will define a regularization constant series. The start value of this 
#'  series is \code{C.base}^\code{C.min}, and the end value of this series is 
#'  \code{C.base}^\code{C.max}, while the near elements of the series have a 
#'  difference of \code{C.base} fold. If the 2 indexes are set as the same, 
#'  the series will become 1 regularization constant. The default value of the 
#'  parameter \code{C.base} here is 10. 
#'@param C.min As mentioned in the \code{C.base} part, this parameter is used
#'  as the index of the small regularization constant number to set a series
#'  for SVM. Default is -3.
#'@param C.max As mentioned in the \code{C.base} part, this parameter is used
#'  as the index of the large regularization constant number to set a series
#'  for SVM. Default is -2.
#'@param alpha The alpha parameter for elastic net feature selection, and it 
#'  controls the balance of L1 and L2 penalties. Default is 0.5.
#'@param errortype The method to choose the regularization constant lambda for 
#'  the elastic net feature selection step. If it is set as "min", the lambda 
#'  that combines with \code{alpha} and gives the minimum mean 10-fold cross- 
#'  validation error will be used, and if it is set as "1ses", the lambda that 
#'  combines with \code{alpha} and gives the cross-validation error within one 
#'  standard error of the minimum will be used. The 10-fold cross-validation 
#'  is defined by the random seed \code{seednum}. For \code{errortype}, its 
#'  default value is "min".
#'@param crosstestacccut If a cross validation is performed as defined by the 
#'  parameter \code{nfold}, the cross validation accuracy on testing data need 
#'  to be greater than this \code{crossacccut}. Default value is 0.7. If it is 
#'  not reached, the function will return NULL.
#'@param plot Whether generate heatmap plots to show the sample classification 
#'  results.
#'@param textsize The font size of the plot texts. Default is 13. 
#'@param prefixes The character string need to be shown in the titles of the 
#'  heatmap plots. It can be set as any character string need to be shown. 
#'  Default is NULL.
#'@return A list with a slot named "mod", which contains the final model from 
#'  the whole dataset. If cross validation is performed, an additional slot 
#'  named "cvtestcomps" will also be returned, which records the true labels 
#'  of the samples and their predicted labels when they are in testing dataset 
#'  during the cross validation.
#'
#'
#'
#'@examples
#'library(CWGCNA)
#'
#'betas <- system.file("extdata", "placentabetas.rds", package = "CWGCNA")
#'betas <- readRDS(betas)
#'
#'pds <- system.file("extdata", "placentapds.rds", package = "CWGCNA")
#'pds <- readRDS(pds)
#'
#'#Extract the DNAm probe data for the 101 preeclampsia samples
#'prepds <- subset(pds, Group == "Preeclampsia")
#'row.names(prepds) <- 1:nrow(prepds)
#'
#'prebetas <- betas[, prepds$sampleid, drop = FALSE]
#'
#'#Clustering
#'presubtyperes <- multiCCA(dats = list(prebetas), 
#'                          
#'  k = 2, consensus = 1, 
#'                          
#'  seednum = 2022, threads = 6, plot = TRUE, titlefix = "Preeclampsia", 
#'  titlesize = 18, textsize = 16)
#'
#'#Make the sample labels from the clustering result
#'subtypes <- paste0("Subtype", presubtyperes$kreses$`k = 2`)
#'
#'#Classification
#'\dontrun{
#'presubtypeclassifierres <- omicsclassifier(dats = list(prebetas), 
#'  truelabels = subtypes, 
#'                                           
#'  balanceadj = 1, 
#'                                           
#'  method = "SVM", 
#'                                           
#'  alphas = c(0.5), nfold = 5, 
#'                                           
#'  seednum = 2022, threads = 6, 
#'                                           
#'  plot = TRUE, prefixes = c("Preeclampsia (SVM-balance)"))
#'}
#'
#'
#'
#'@export
omicsclassifier <- function(dats, 
                            truelabels, 
                            varfeatures = NULL, 
                            scale = TRUE, 
                            balanceadj = 1, 
                            #samplingmethod = 'updn', 
                            nround = 10, 
                            k = 5, 
                            nfold = 5, 
                            seednum = 2022, 
                            threads = 1, 
                            prescreen = TRUE, 
                            
                            method = 'SVM', 
                            C.base = 10,
                            C.min = -3,
                            C.max = -2,
                            
                            alphas = c(0.5), 
                            errortype = 'min', 
                            crosstestacccut = 0.7, 
                            
                            plot = TRUE, 
                            textsize = 13, 
                            prefixes = NULL){
  
  
  samplingmethod <- 'updn'
  
  
  for(i in 1:length(dats)){
    
    row.names(dats[[i]]) <- make.names(row.names(dats[[i]]))
    
  }
  
  
  if(!is.factor(truelabels) & !is.numeric(truelabels)){
    
    responselevels <- names(table(truelabels))[order(-table(truelabels))]
    
    truelabels <- factor(truelabels, 
                         levels = responselevels,
                         ordered = TRUE)
    
  }
  
  if(!is.null(names(truelabels))){
    
    sharedsamples <- intersect(colnames(dats[[1]]), names(truelabels))
    
    if(length(sharedsamples) > 0){
      dats[[1]] <- dats[[1]][, sharedsamples, drop = FALSE]
      truelabels <- truelabels[sharedsamples]
    }
    
  }
  
  if(length(dats) > 1){
    
    sharedsamples <- colnames(dats[[1]])
    
    for(i in 2:length(dats)){
      sharedsamples <- intersect(sharedsamples, colnames(dats[[i]]))
    }
    
    if(!is.null(names(truelabels))){
      sharedsamples <- intersect(sharedsamples, names(truelabels))
      
    }
    
    if(length(sharedsamples) > 0){
      
      for(i in 1:length(dats)){
        
        dats[[i]] <- dats[[i]][, sharedsamples, drop = FALSE]
        
      }
      
      if(!is.null(names(truelabels))){
        truelabels <- truelabels[sharedsamples]
        
      }
    }
    
  }
  
  #Variable feature selection
  if(!is.null(varfeatures)){
    
    for(i in 1:length(dats)){
      
      datfeatures <- featuresampling(betas = dats[[i]], 
                                     topfeatures = varfeatures, 
                                     pddat = NULL, 
                                     variancetype = 'sd', 
                                     anova = FALSE)
      
      dats[[i]] <- datfeatures$betas[datfeatures$varicanceres[row.names(datfeatures$betas),] != 0, 
                                     , drop = FALSE]
      
    }
    
  }
  
  #Scaling
  if(scale == TRUE){
    
    for(i in 1:length(dats)){
      
      datsamples <- colnames(dats[[i]])
      dats[[i]] <- apply(X = dats[[i]], MARGIN = 1, FUN = scale)
      row.names(dats[[i]]) <- datsamples
      dats[[i]] <- t(dats[[i]])
      
    }
  }
  
  pddat <- data.frame(sampleid = colnames(dats[[1]]), 
                      #Response = as.vector(truelabels), 
                      stringsAsFactors = FALSE)
  
  pddat$Response <- truelabels
  row.names(pddat) <- 1:nrow(pddat)
  
  
  if(prescreen == TRUE){
    
    for(i in 1:length(dats)){
      
      datfeatures <- featuresampling(betas = dats[[i]], 
                                     topfeatures = nrow(dats[[i]]), 
                                     pddat = pddat, 
                                     anova = TRUE, 
                                     responsevarname = 'Response', 
                                     threads = threads, 
                                     plotannovares = FALSE)
      
      datselectedfeatures <- row.names(datfeatures$varicanceres)[
        datfeatures$varicanceres$Response_pval < 0.05]
      
      if(length(datselectedfeatures) < 10){
        datselectedfeatures <- row.names(dats[[i]])
      }
      
      dats[[i]] <- dats[[i]][datselectedfeatures, , drop = FALSE]
      
    }
    
  }
  
  testcomps <- NULL
  
  if(!is.null(nfold)){
    
    sampledivides <- list()
    
    for(i in 1:length(dats)){
      
      #Part1, divide the original dataset into train and test sets
      sampledivides[[i]] <- dividecv(i = seednum, nfold = nfold, totalpd = pddat, 
                                     totalvar = t(dats[[i]]), factortonum = FALSE)
      
      
    }
    
    
    if(!is.null(sampledivides[[1]]$testlist)){
      
      cvfeatureses <- list()
      
      for(i in 1:length(sampledivides)){
        cvfeatureses[[i]] <- list()
      }
      
      innertrainvars <- list()
      innertestvars <- list()
      
      k <- 1
      for(k in 1:length(sampledivides[[1]]$trainlist)){
        
        cat(paste0('CV = ', k, '\n'))
        
        for(i in 1:length(sampledivides)){
          
          innertrainvars[[i]] <- t(sampledivides[[i]]$trainlist[[k]]$innertrainvar)
          
          innertestvars[[i]] <- sampledivides[[i]]$testlist[[k]]$innertestvar
          
        }
        
        #Part2, elastic net on train dataset
        
        pairedensemblecv <- pairedbalancemod(trainingdats = innertrainvars, 
                                             trainingpd = sampledivides[[1]]$trainlist[[k]]$innertrainpd, 
                                             balanceadj = balanceadj, 
                                             samplingmethod = samplingmethod, 
                                             nround = nround, 
                                             k = k, 
                                             threads = threads, 
                                             
                                             method = method, 
                                             C.base = C.base, 
                                             C.min = C.min, 
                                             C.max = C.max, 
                                             
                                             alphas = alphas, 
                                             seednum = seednum, 
                                             errortype = errortype, 
                                             
                                             plot = FALSE, 
                                             prefixes = prefixes, 
                                             textsize = textsize, 
                                             verbose = FALSE)
        
        for(i in 1:length(sampledivides)){
          
          cvfeatureses[[i]][[k]] <- row.names(pairedensemblecv$baselearners[[i]]$modelscores[[1]])
          
        }
        
        #Part3, summarize prediction on test dataset
        
        traincomp <- pairedensemblecv$traincomp
        
        testprecv <- pairedensemblepredict(balanceres = pairedensemblecv, 
                                           vardats = innertestvars)
        
        
        testcomp <- cbind(testprecv$pres, sampledivides[[1]]$testlist[[k]]$innertestpd$Response)
        colnames(testcomp) <- c('Prediction', 'True')
        
        if(k == 1){
          traincomps <- traincomp
          testcomps <- testcomp
        }else{
          traincomps <- rbind(traincomps, traincomp)
          testcomps <- rbind(testcomps, testcomp)
        }
        
      }
      
      
      crosstestacc <- sum(as.character(testcomps[,1]) == as.character(testcomps[,2]))/nrow(testcomps)
      
      
      if(plot == TRUE){
        
        
        heatmapplotting(orivar = t(dats[[1]]), 
                        comptab = apply(X = testcomps, MARGIN = 2, 
                                        FUN = function(x){if(is.numeric(x)){return(factor(x))}else{return(x)}}), 
                        featurenames = unique(unlist(cvfeatureses[[1]])), 
                        title = paste0(prefixes[1], 
                                       c('Test set with features from all crosses', 
                                         ' test set with features from all crosses\n')[c(is.null(prefixes[1]), 
                                                                                         !is.null(prefixes[1]))]), 
                        textsize = textsize)
        
        if(length(dats) > 1){
          
          for(i in 2:length(dats)){
            
            heatmapplotting(orivar = t(dats[[i]]), 
                            comptab = apply(X = testcomps, MARGIN = 2, 
                                            FUN = function(x){if(is.numeric(x)){return(factor(x))}else{return(x)}}), 
                            featurenames = unique(unlist(cvfeatureses[[i]])), 
                            title = paste0(prefixes[i], 
                                           c('Test set with features from all crosses', 
                                             ' test set with features from all crosses\n')[c(is.null(prefixes[i]), 
                                                                                             !is.null(prefixes[i]))]), 
                            textsize = textsize)
            
          }
          
        }
        
      }
      
      if(is.na(crosstestacc)){
        
        cat('The ACC of the test from cross validation is NA, so drop the fitting\n')
        
        return(NULL)
        
      }else if(crosstestacc < crosstestacccut){
        
        cat(paste0('The ACC of the test from cross validation is ', 
                   signif(crosstestacc, 3), ' < ', signif(crosstestacccut, 3), 
                   ' (the ACC cutoff), so drop the fitting\n'))
        
        return(NULL)
        
      }
      
    }else{
      testcomps <- NULL
    }
    
  }
  
  pairedensemble <- pairedbalancemod(trainingdats = dats, 
                                     trainingpd = pddat, 
                                     balanceadj = balanceadj, 
                                     samplingmethod = samplingmethod, 
                                     nround = nround, 
                                     k = k, 
                                     threads = threads, 
                                     
                                     method = method, 
                                     C.base = C.base, 
                                     C.min = C.min, 
                                     C.max = C.max, 
                                     
                                     alphas = alphas, 
                                     seednum = seednum, 
                                     errortype = errortype, 
                                     
                                     plot = plot, 
                                     prefixes = prefixes, 
                                     textsize = textsize)
  
  
  res <- list(mod = pairedensemble, 
              cvtestcomps = testcomps)
  
  return(res)
  
}





