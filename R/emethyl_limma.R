#ACAT####

ACAT<-function(Pvals,weights=NULL,is.check=TRUE){
  Pvals<-as.matrix(Pvals)
  if (is.check){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
      stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
      stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(colSums(Pvals==0)>=1)
    is.one<-(colSums(Pvals==1)>=1)
    if (sum((is.zero+is.one)==2)>0){
      stop("Cannot have both 0 and 1 p-values in the same column!")
    }
    
    if (sum(is.zero)>0){
      warning("There are p-values that are exactly 0!")
    }
    if (sum(is.one)>0){
      warning("There are p-values that are exactly 1!")
    }
    
  }
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(weights)){
    is.weights.null<-TRUE
  }else{
    is.weights.null<-FALSE
    weights<-as.matrix(weights)
    if (sum(dim(weights)!=dim(Pvals))>0){
      stop("The dimensions of weights and Pvals must be the same!")
    }else if (is.check & (sum(weights<0)>0)){
      stop("All the weights must be nonnegative!")
    }else{
      w.sum<-colSums(weights)
      if (sum(w.sum<=0)>0){
        stop("At least one weight should be positive in each column!")
      }else{
        for (j in 1:ncol(weights)){
          weights[,j]<-weights[,j]/w.sum[j]
        }
      }
    }
    
  }
  
  #### check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small<-(Pvals<1e-15)
  if (is.weights.null){
    Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-1/Pvals[is.small]/pi
    cct.stat<-colMeans(Pvals)
  }else{
    Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
    cct.stat<-colSums(Pvals)
  }
  #### return the ACAT p value(s).
  pval<-pcauchy(cct.stat,lower.tail = F)
  return(pval)
}

ACAT_V<-function(G,obj,weights.beta=c(1,25),weights=NULL,mac.thresh=10){
  ### check weights
  if (!is.null(weights) && length(weights)!=ncol(G)){
    stop("The length of weights must equal to the number of variants!")
  }
  
  mac<-Matrix::colSums(G)
  ### remove SNPs with mac=0
  if (sum(mac==0)>0){
    G<-G[,mac>0,drop=FALSE]
    weights<-weights[mac>0]
    mac<-mac[mac>0]
    if (length(mac)==0){
      stop("The genotype matrix do not have non-zero element!")
    }
  }
  ### p and n
  p<-length(mac)
  n<-nrow(G)
  ###
  
  if (sum(mac>mac.thresh)==0){  ## only Burden
    pval<-Burden(G,obj, weights.beta = weights.beta, weights = weights)
  }else if (sum(mac<=mac.thresh)==0){ ## only cauchy method
    if (is.null(weights)){
      MAF<-mac/(2*n)
      W<-(dbeta(MAF,weights.beta[1],weights.beta[2])/dbeta(MAF,0.5,0.5))^2
    }else{
      W<-weights
    }
    
    Mpvals<-Get.marginal.pval(G,obj)
    pval<-ACAT(Mpvals,W)
  }else{  ## Burden + Cauchy method
    is.very.rare<-mac<=mac.thresh
    weights.sparse<-weights[is.very.rare]
    weights.dense<-weights[!is.very.rare]
    pval.dense<-Burden(G[,is.very.rare,drop=FALSE],obj, weights.beta = weights.beta, weights = weights.sparse)
    
    Mpvals<-Get.marginal.pval(G[,!is.very.rare,drop=FALSE],obj)
    
    Mpvals<-c(Mpvals,pval.dense)
    if (is.null(weights)){
      MAF<-mac/(2*n)
      mafs<-c(MAF[!is.very.rare],mean(MAF[is.very.rare])) ## maf for p-values
      W<-(dbeta(mafs,weights.beta[1],weights.beta[2])/dbeta(mafs,0.5,0.5))^2
    }else{
      W<-c(weights.dense,mean(weights.sparse))
    }
    
    
    is.keep<-rep(T,length(Mpvals))
    is.keep[which(Mpvals==1)]<-F  ## remove p-values of 1.
    pval<-ACAT(Mpvals[is.keep],W[is.keep])
  }
  return(pval)
}

NULL_Model<-function(Y,Z=NULL){
  n<-length(Y)
  #### check the type of Y
  if ((sum(Y==0)+sum(Y==1))==n){
    out_type<-"D"
  }else{
    out_type<-"C"
  }
  #### Add intercept
  Z.tilde<-cbind(rep(1,length(Y)),Z)
  if (out_type=="C"){
    #### estimate of sigma square
    Z.med<-Z.tilde%*%solve(chol(t(Z.tilde)%*%Z.tilde))   ## Z.med%*%t(Z.med) is the projection matrix of Z.tilde
    Y.res<-as.vector(Y-(Y%*%Z.med)%*%t(Z.med))
    sigma2<-sum(Y.res^2)/(n-ncol(Z.med))
    #### output
    res<-list()
    res[["out_type"]]<-out_type
    res[["Z.med"]]<-Z.med
    res[["Y.res"]]<-Y.res
    res[["sigma2"]]<-sigma2
  }else if (out_type=="D"){
    #### fit null model
    g<-glm(Y~0+Z.tilde,family = "binomial")
    prob.est<-g[["fitted.values"]]
    #### unstandarized residuals
    Y.res<-(Y-prob.est)
    ### Sigma when rho=0
    sigma2.Y<-prob.est*(1-prob.est)  ### variance of each Y_i
    ### output
    res<-list()
    res[["out_type"]]<-out_type
    res[["Z.tilde"]]<-Z.tilde
    res[["Y.res"]]<-Y.res
    res[["sigma2.Y"]]<-sigma2.Y
  }
  return(res)
}

Get.marginal.pval<-function(G,obj){
  ### check obj
  if (names(obj)[1]!="out_type"){
    stop("obj is not calculated from MOAT_NULL_MODEL!")
  }else{
    out_type<-obj[["out_type"]]
    if (out_type=="C"){
      if (!all.equal(names(obj)[2:length(obj)],c("Z.med","Y.res","sigma2"))){
        stop("obj is not calculated from MOAT_NULL_MODEL!")
      }else{
        Z.med<-obj[["Z.med"]]
        Y.res<-obj[["Y.res"]]
        n<-length(Y.res)
        SST<-obj[["sigma2"]]*(n-ncol(Z.med))
      }
    }else if (out_type=="D"){
      if (!all.equal(names(obj)[2:length(obj)],c("Z.tilde","Y.res","sigma2.Y"))){
        stop("obj is not calculated from MOAT_NULL_MODEL!")
      }else{
        Z.tilde<-obj[["Z.tilde"]]
        Y.res<-obj[["Y.res"]]
        sigma2.Y<-obj[["sigma2.Y"]]
        n<-length(Y.res)
      }
    }
  }
  
  if (class(G)!="matrix" && class(G)!="dgCMatrix"){
    stop("The class of G must be matrix or dgCMatrix!")
  }
  
  if (out_type=="C"){
    G_tX.med<-as.matrix(Matrix::crossprod(Z.med,G))
    ### Sigma^2 of G
    Sigma2.G<-Matrix::colSums(G^2)-Matrix::colSums(G_tX.med^2)
    SSR<-as.vector((Y.res%*%G)^2/Sigma2.G)
    SSR[Sigma2.G<=0]<-0
    df.2<-n-1-ncol(Z.med)
    t.stat<-suppressWarnings(sqrt(SSR/((SST-SSR)/df.2)))
    marginal.pvals<-2*pt(t.stat,(n-1-ncol(Z.med)),lower.tail = FALSE)
  }else if (out_type=="D"){
    Z.stat0<-as.vector(Y.res%*%G)
    ### Sigma when rho=0
    tG_X.tilde_sigma2<-as.matrix(Matrix::crossprod(G,Z.tilde*sigma2.Y))
    Sigma2.G<-Matrix::colSums(G^2*sigma2.Y)-diag(tG_X.tilde_sigma2%*%solve(t(Z.tilde)%*%(Z.tilde*sigma2.Y))%*%t(tG_X.tilde_sigma2))
    marginal.pvals<-2*pnorm(abs(Z.stat0)/sqrt(Sigma2.G),lower.tail = FALSE)
  }
  
  return(marginal.pvals)
}

Burden<-function(G,obj,kernel="linear.weighted",weights.beta=c(1,25),weights=NULL){
  ### check obj
  if (names(obj)[1]!="out_type"){
    stop("obj is not calculated from NULL_MODEL!")
  }else{
    out_type<-obj[["out_type"]]
    if (out_type=="C"){
      if (!all.equal(names(obj)[2:length(obj)],c("Z.med","Y.res","sigma2"))){
        stop("obj is not calculated from NULL_MODEL!")
      }else{
        Z.med<-obj[["Z.med"]]
        Y.res<-obj[["Y.res"]]/sqrt(obj[["sigma2"]])  ## rescaled residules
        n<-length(Y.res)
      }
    }else if (out_type=="D"){
      if (!all.equal(names(obj)[2:length(obj)],c("Z.tilde","Y.res","sigma2.Y"))){
        stop("obj is not calculated from NULL_MODEL!")
      }else{
        Z.tilde<-obj[["Z.tilde"]]
        Y.res<-obj[["Y.res"]]
        sigma2.Y<-obj[["sigma2.Y"]]
        n<-length(Y.res)
      }
    }
  }
  ### MAF
  MAF<-Matrix::colSums(G)/(2*dim(G)[1])
  p<-length(MAF)
  #### weights
  if (kernel=="linear.weighted"){
    if (is.null(weights)){
      W<-dbeta(MAF,weights.beta[1],weights.beta[2])
    }else{
      if (length(weights)==p){
        W<-weights
      }else{
        stop("The length of weights must equal to the number of variants!")
      }
    }
    
  }else if (kernel=="linear"){
    W<-rep(1,p)
  }else{
    stop("The kernel name is not valid!")
  }
  
  ###### if G is sparse or not
  if (class(G)=="matrix" || class(G)=="dgCMatrix"){
    if (out_type=="C"){
      Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
      Gw<-G%*%W
      sigma.z<-sqrt(sum(Gw^2)-sum((t(Z.med)%*%(Gw))^2))
    }else if (out_type=="D"){
      Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
      Gw<-as.vector(G%*%W)
      sigma.z<-sum(Gw^2*sigma2.Y)-((Gw*sigma2.Y)%*%Z.tilde)%*%solve(t(Z.tilde)%*%(Z.tilde*sigma2.Y))%*%t((Gw*sigma2.Y)%*%Z.tilde)
      sigma.z<-as.vector(sqrt(sigma.z))
    }
  }else{
    stop("The class of G must be matrix or dgCMatrix!")
  }
  
  V<-Z.stat.sum/sigma.z
  Q<-V^2   ## Q test statistic
  pval<-1-pchisq(Q,df=1)
  return(pval)
}

#Pd functions####

betatoM <- function(bval){
  
  minval <- min(bval[bval > 0])
  maxval <- max(bval[bval < 1])
  
  bval[bval <= 0] <- minval
  bval[bval >= 1] <- maxval
  
  mval <- log2(bval/(1 - bval))
  
  return(mval)
}

Mtobeta <- function(mval){
  
  bval <- 2^mval/(1 + 2^mval)
  
  minval <- min(bval[bval > 0])
  maxval <- max(bval[bval < 1])
  
  bval[bval <= 0] <- minval
  bval[bval > 1] <- maxval
  
  return(bval)
  
}

pathtrailing <- function(sym = '/', path = getwd()){
  pathelements <- strsplit(x = path, split = '', fixed = TRUE)
  lastsym <- pathelements[[1]][[length(pathelements[[1]])]]
  if(lastsym == sym){
    newpath <- pathelements[[1]][-length(pathelements[[1]])]
    newpath <- paste0(newpath, collapse = '')
    pathtrailing(sym = sym, path = newpath)
  }else{
    return(path)
  }
}

countna <- function(colval){
  nacount <- sum(is.na(colval))
  return(nacount)
  
}

#'Impute missing values in meta data
#'
#'Impute missing values in meta data using the MICE (multivariate imputations 
#'by chained equations) method.
#'
#'@param pddat Meta data data frame. The first column should be the sample 
#'  names, and the remaining ones should be all the confounding factors in the 
#'  data, including both the ones with missing values and the complete ones.
#'@param percentcutoff If a confounding factor has a missing value proportion 
#'  greater than this \code{percentcutoff} value, it will not be imputed and 
#'  will not appear in the final imputed meta data frame returned. Default is 
#'  1/4.
#'@param seednum A random seed number used by the MICE algorithm to impute the 
#'  missing values. Default is 2022.
#'@return A data frame with the original complete confounding factors and the 
#'  ones imputed as complete.
#'@export
imputemeta <- function(pddat, 
                       percentcutoff = 1/4, 
                       seednum = 2022){
  
  #library(mice)
  
  sampleidvarname <- colnames(pddat)[1]
  row.names(pddat) <- pddat[,1]
  pddat <- pddat[-1]
  
  charvars <- c()
  i <- 1
  for(i in 1:ncol(pddat)){
    
    if(is.character(pddat[,i])){
      charvar <- colnames(pddat)[i]
      charvars <- c(charvars, charvar)
      pddat[,i] <- as.factor(pddat[,i])
      
    }
  }
  
  
  pddatnacount <- apply(X = pddat, MARGIN = 2, FUN = countna)
  pddatnacount <- pddatnacount/nrow(pddat)
  reserved <- names(pddatnacount)[pddatnacount <= percentcutoff]
  
  #mice::md.pattern(newpddat)
  
  pddatimp <- mice::mice(pddat, maxit = 25, m = 5, seed = seednum)
  
  #mice::stripplot(pddatimp, Age, pch = 19, xlab = 'Imputation number')
  
  impvar <- mice::complete(pddatimp)
  
  impvar <- impvar[reserved]
  
  charvars <- intersect(charvars, reserved)
  
  i <- 1
  for(i in 1:length(charvars)){
    charvar <- charvars[i]
    impvar[,charvar] <- as.character(impvar[,charvar])
    
  }
  
  impvar$sampleid <- row.names(impvar)
  impvar <- impvar[c(sampleidvarname, colnames(impvar)[1:ncol(impvar)-1])]
  row.names(impvar) <- 1:nrow(impvar)
  
  return(impvar)
  
}

#Some parallel computing going on in the background that is not getting cleaned up
#fully between runs can cause Error in summary.connection(connection) : invalid
#connection. The following function is needed to be called to fix this error.
unregister_dopar <- function(){
  
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  
}

anovaplot <- function(restab, 
                      plottype = 'pval', 
                      titlesufix = NULL, 
                      samplesize = NULL, 
                      titlesize = 15, 
                      textsize = 13, 
                      face = 'bold'){
  
  #library(ggplot2)
  
  if(!is.null(titlesufix)){
    titlesufix <- paste0(' ', titlesufix)
  }
  
  resmean <- colMeans(restab)
  resmean <- as.matrix(resmean)
  resmean <- t(resmean)
  row.names(resmean) <- 1:nrow(resmean)
  
  if(plottype == 'MSS'){
    
    resmean <- resmean[, grepl(pattern = '_MSS', colnames(resmean)), drop = FALSE]
    colnames(resmean) <- gsub(pattern = '_MSS', replacement = '', x = colnames(resmean))
    ytitle <- 'MSS'
    
  }else if(plottype == 'pval'){
    
    resmean <- resmean[, grepl(pattern = '_pval', colnames(resmean)), drop = FALSE]
    colnames(resmean) <- gsub(pattern = '_pval', replacement = '', x = colnames(resmean))
    resmean <- -log10(resmean)
    ytitle <- '-log10(pval)'
    
  }else{
    
    resmean <- resmean[, grepl(pattern = '_F', colnames(resmean)), drop = FALSE]
    colnames(resmean) <- gsub(pattern = '_F', replacement = '', x = colnames(resmean))
    ytitle <- 'F statistic'
    
  }
  
  resmean <- t(resmean)
  colnames(resmean) <- 'Statistic'
  resmean <- as.data.frame(resmean, stringsAsFactors = FALSE)
  
  resmean$Factor <- factor(row.names(resmean), levels = row.names(resmean), 
                           ordered = TRUE)
  row.names(resmean) <- 1:nrow(resmean)
  
  if(!is.null(samplesize)){
    subtitle <- paste0('(sample size = ', samplesize, ')')
  }
  
  p <- ggplot2::ggplot(data = resmean, 
                       mapping = ggplot2::aes(x = Factor, 
                                              y = Statistic, fill = Factor))
  p <- p + ggplot2::geom_bar(stat = 'identity') + 
    ggplot2::ggtitle(paste0('Type 3 ANOVA', titlesufix, ' ', ytitle), subtitle = subtitle) + 
    ggplot2::ylab(paste0('Averaged feature ', ytitle)) + 
    ggplot2::xlab('') + 
    ggplot2::scale_fill_discrete(guide = 'none') + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = titlesize, face = face),
                   axis.text.x = ggplot2::element_text(face = face, size = textsize,
                                                       angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(face = face, size = textsize - 3),
                   plot.title = ggplot2::element_text(size = titlesize, face = face), 
                   plot.subtitle = ggplot2::element_text(size = textsize, face = face))
  
  if(plottype == 'pval'){
    p <- p + ggplot2::geom_hline(yintercept = -log10(0.05), col = 'red')
  }
  
  print(p)
  
}

#'Select data features with standard deviation (SD), median absolute deviation 
#'(MAD), or type-III ANOVA.
#'
#'Select features with the largest SD, largest MAD, or smallest ANOVA p-val. 
#'
#'@param betas Data matrix for feature selection, with features as rows and 
#'  samples as columns.
#'@param topfeatures The number of features need to select. Default is 10000. 
#'@param pddat Meta data frame. The first column should be the sample names, 
#'  and the remaining should be the phenotypic variables need to be analyzed 
#'  with ANOVA. Only needed for ANOVA-based feature selection. Default value 
#'  is NULL. 
#'@param variancetype The method need to calculate feature variance. Default 
#'  is 'sd', meaning standard deviation. Can also be 'mad', meaning median 
#'  absolute deviation.
#'@param anova Whether need to perform ANOVA-based feature selection. Default 
#'  is FALSE. If it is TRUE, and \code{pddat} is also provided, the function 
#'  will perform ANOVA for each feature in the data, so its variance explained 
#'  by each phenotypic variable in \code{pddat} can be obtained. Then, the 
#'  features will be ordered by its ANOVA p-val from the response variable in 
#'  \code{pddat} (defined by the parameter \code{responsevarname}), and those 
#'  with the smallest p-val will be selected and returned. If \code{anova} is 
#'  FALSE, or \code{pddat} is not provided, the features with the largest SD 
#'  or MAD variance will be selected, as defined by \code{variancetype}.
#'@param responsevarname The column name of the response variable in the data 
#'  frame provided to \code{pddat}. It is the variable of interest, and other 
#'  variables in \code{pddat} will be the confounding variables. Default is 
#'  NULL, which means the second column of \code{pddat} will be the response. 
#'@param filternovar Whether features with no variance across the samples need 
#'  to be excluded from the feature selection. Default is TRUE. 
#'@param threads Number of threads need for parallelization. Default is 1.
#'@param plotannovares Whether need to plot the ANOVA result. Default is TRUE. 
#'  In this case, the ANOVA results from all the features will be averaged, 
#'  including their F statistics, MSS (mean sum of the square), and p-vals, 
#'  and the final will be that for the dataset. Then, a bar plot will be made 
#'  for this dataset result.
#'@param titlesize The font size for the plot title of \code{plotannovares}. 
#'  Default is 15.
#'@param textsize The font size for the plot texts of \code{plotannovares}. 
#'  Default is 13.
#'@param face The font face for the plot texts of \code{plotannovares}. The 
#'  default value is "bold".
#'@param featuretype The type of the features in the data. Will appear in the 
#'  title of the ANOVA plot. Can be "gene", "probe", and so on. Default value 
#'  is NULL, so that this information will not appear.
#'@param plottitlesuffx The prefix of the ANOVA plot title. It can be set as 
#'  any character string need to be shown in the title. Default is NULL.
#'@return A list with 2 elements. The first is the final data matrix only with 
#'  the selected features. The second is the selection information of all the 
#'  features in the original data, such as their SD (if using SD for feature 
#'  selection), MAD (if using MAD), or their ANOVA F statistic, MSS, and p-val 
#'  (if using ANOVA feature selection). If \code{plotannovares} were TRUE, the 
#'  ANOVA result for the dataset will also be plotted.
#'@examples
#'library(CWGCNA)
#'
#'betas <- system.file("extdata", "placentabetas.rds", package = "CWGCNA")
#'betas <- readRDS(betas)
#'
#'pds <- system.file("extdata", "placentapds.rds", package = "CWGCNA")
#'pds <- readRDS(pds)
#'
#'top10k <- featuresampling(betas = betas, topfeatures = 10000, 
#'  variancetype = "sd", threads = 6)
#'  
#'anovares <- featuresampling(betas = top10k$betas, pddat = pds, anova = TRUE, 
#'  
#'  plotannovares = TRUE, featuretype = "probe", plottitlesuffix = "placenta", 
#'  titlesize = 18, textsize = 16, 
#'  
#'  threads = 6)
#'@export
featuresampling <- function(betas, 
                            topfeatures = 10000, 
                            pddat = NULL, 
                            variancetype = 'sd', 
                            anova = FALSE, 
                            responsevarname = NULL, 
                            filternovar = TRUE, 
                            threads = 1, 
                            plotannovares = TRUE, 
                            titlesize = 15, 
                            textsize = 13, 
                            face = 'bold', 
                            featuretype = NULL, 
                            plottitlesuffix = NULL){
  
  if(topfeatures <= 0){
    topfeatures <- 1000
  }
  
  if(topfeatures > nrow(betas)){
    topfeatures <- nrow(betas)
  }
  
  if(variancetype == 'mad'){
    
    calvar <- function(line){
      
      linemedian <- median(line, na.rm = TRUE)
      res <- median(abs(line - linemedian), na.rm = TRUE)
      names(res) <- 'MAD'
      
      return(res)
      
    }
    
  }else{
    
    calvar <- function(line){
      
      res <- sd(line, na.rm = TRUE)
      names(res) <- 'SD'
      
      return(res)
      
    }
    
    
  }
  
  if(anova == TRUE & (!is.null(pddat))){
    
    betas <- betas[,pddat[,1], drop = FALSE]
    
    if(filternovar == TRUE){
      
      filteridx <- apply(X = betas, MARGIN = 1, FUN = sd, na.rm = TRUE)
      
      betas <- betas[filteridx != 0, , drop = FALSE]
      
    }
    
    
    
    calvar <- function(line, pds = pddat){
      
      #library(car)
      
      newdata <- pds[, -1, drop = FALSE]
      pdnames <- colnames(newdata)
      newdata$beta <- line
      
      formstr <- paste0(pdnames, collapse = ' + ')
      formstr <- paste0('beta ~ ', formstr)
      formstr <- as.formula(formstr)
      
      fit <- lm(formstr, data = newdata)
      
      aovfit <- car::Anova(fit, type = 3, singular.ok = TRUE)
      
      F <- aovfit$`F value`
      F <- F[2:(length(F)-1)]
      names(F) <- pdnames
      F <- as.data.frame(F, stringsAsFactors = FALSE)
      F <- as.data.frame(t(F))
      row.names(F) <- 1
      colnames(F) <- paste0(colnames(F), '_F')
      
      pvals <- aovfit$`Pr(>F)`
      pvals <- pvals[2:(length(pvals)-1)]
      names(pvals) <- pdnames
      pvals <- as.data.frame(pvals, stringsAsFactors = FALSE)
      pvals <- as.data.frame(t(pvals))
      row.names(pvals) <- 1
      colnames(pvals) <- paste0(colnames(pvals), '_pval')
      
      SS <- aovfit$`Sum Sq`
      DF <- aovfit$Df
      MSS <- SS/DF
      
      MSS <- MSS[2:(length(MSS)-1)]
      names(MSS) <- pdnames
      MSS <- as.data.frame(MSS, stringsAsFactors = FALSE)
      MSS <- as.data.frame(t(MSS))
      row.names(MSS) <- 1
      colnames(MSS) <- paste0(colnames(MSS), '_MSS')
      
      lineres <- cbind(F, pvals, MSS)
      row.names(lineres) <- row.names(line)
      
      return(lineres)
      
    }
    
  }
  
  if(threads == 1){
    
    varlist <- list()
    for(i in 1:nrow(betas)){
      varlist[[i]] <- calvar(line = betas[i,])
    }
    
  }else{
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    varlist <- foreach::foreach(i = 1:nrow(betas), 
                                .export = NULL) %dopar% {
                                  calvar(line = betas[i,])
                                }
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
    
  }
  
  varres <- do.call(rbind, varlist)
  row.names(varres) <- row.names(betas)
  
  if(anova == TRUE){
    
    if(is.null(responsevarname)){
      if(!is.null(pddat)){
        responsevarname <- names(pddat)[2]
      }
      
    }
    
    if(!is.null(responsevarname)){
      varres <- varres[order(varres[,paste0(responsevarname, '_pval')]), , drop = FALSE]
    }
    
  }else{
    if(variancetype == 'mad'){
      varres <- varres[order(-varres[,'MAD']), , drop = FALSE]
    }else{
      varres <- varres[order(-varres[,'SD']), , drop = FALSE]
    }
  }
  
  
  topvarres <- varres[1:topfeatures, , drop = FALSE]
  
  
  topbetas <- betas[(row.names(betas) %in% row.names(topvarres)), , drop = FALSE]
  
  #topbetas <- betas[1:topfeatures,]
  
  res <- list(betas = topbetas, 
              varicanceres = varres)
  
  if(anova == TRUE & (!is.null(pddat))){
    
    if(is.null(featuretype)){
      featuretype <- 'feature'
    }
    
    if(plotannovares == TRUE){
      
      if(!is.null(plottitlesuffix)){
        plottitlesuffix <- paste0('for ', plottitlesuffix, ' top ')
      }else{
        plottitlesuffix <- 'for top '
      }
      
      anovaplot(restab = topvarres, plottype = 'pval', 
                titlesufix = paste0(plottitlesuffix, topfeatures, ' ', featuretype, 's'), 
                titlesize = titlesize, textsize = textsize, face = face, 
                samplesize = ncol(betas))
      anovaplot(restab = topvarres, plottype = 'MSS', 
                titlesufix = paste0(plottitlesuffix, topfeatures, ' ', featuretype, 's'), 
                titlesize = titlesize, textsize = textsize, face = face, 
                samplesize = ncol(betas))
      anovaplot(restab = topvarres, plottype = 'F', 
                titlesufix = paste0(plottitlesuffix, topfeatures, ' ', featuretype, 's'), 
                titlesize = titlesize, textsize = textsize, face = face, 
                samplesize = ncol(betas))
      
    }
    
    
    
  }
  
  return(res)
  
  
}

tsneplot <- function(tsnematrix, 
                     perplexity = 4, 
                     colors, 
                     response = NULL, 
                     titlesuffix = NULL, 
                     subtitle = NULL, 
                     responsevarname = 'Cluster', 
                     titlesize = 15, 
                     face = 'bold', 
                     textsize = 13, 
                     labelpoints = TRUE, 
                     embedding = 'tsne', 
                     seednum = 2022){
  
  if(is.null(tsnematrix)){
    return(NULL)
  }
  
  if(ncol(tsnematrix) < 2){
    return(NULL)
  }
  
  if(embedding == 'tsne'){
    set.seed(seednum) # Set a seed if you want reproducible results
    tsne_out <- Rtsne::Rtsne(tsnematrix, perplexity=perplexity, max_iter = 1000) #Run tSNE
    
    tsnedat <- data.frame(tSNE1 = tsne_out$Y[,1], tSNE2 = tsne_out$Y[,2])
    
    xlabel <- 'tSNE1'
    ylabel <- 'tSNE2'
    plotlabel <- 'tSNE'
    
  }else if(embedding == 'umap'){
    #library(umap)
    set.seed(seednum)
    umap_out <- umap::umap(tsnematrix)
    
    tsnedat <- data.frame(tSNE1 = umap_out$layout[,1], tSNE2 = umap_out$layout[,2])
    
    xlabel <- 'UMAP1'
    ylabel <- 'UMAP2'
    plotlabel <- 'UMAP'
    
  }else{
    
    #PCA for variable probes
    pcares <- prcomp(tsnematrix)
    
    tsnedat <- pcares$x[,c(1, 2)]
    tsnedat <- as.data.frame(tsnedat, stringsAsFactors = FALSE)
    colnames(tsnedat) <- c('tSNE1', 'tSNE2')
    
    xlabel <- 'PC1'
    ylabel <- 'PC2'
    plotlabel <- 'PCA'
    
  }
  
  
  
  if(is.null(subtitle)){
    featurenum <- ncol(tsnematrix)
    samplenum <- nrow(tsnematrix)
    subtitle <- paste0(samplenum, ' samples; ', featurenum, ' features')
  }
  
  if(!is.null(response)){
    
    tsnedat$Response <- response
    
    if(labelpoints == TRUE){
      
      p <- ggplot2::ggplot(data = tsnedat, mapping = ggplot2::aes(x = tSNE1, 
                                                                  y = tSNE2, 
                                                                  color = Response, 
                                                                  label = Response))
      
      print(
        
        p + ggplot2::geom_point(size = 2) + 
          ggplot2::xlab(xlabel) + 
          ggplot2::ylab(ylabel) + 
          ggplot2::ggtitle(paste0(plotlabel, titlesuffix), 
                           subtitle = paste0('(', subtitle, ')')) + 
          
          ggplot2::geom_text(hjust = 0, vjust = 0, size = annotextsize) + 
          ggplot2::scale_color_manual(values = colors) + 
          ggplot2::theme_bw() + 
          ggplot2::theme(axis.title.y = ggplot2::element_text(size = titlesize, face = face),
                         axis.title.x = ggplot2::element_text(size = titlesize, face = face), 
                         axis.text.x = ggplot2::element_text(face = face, size = textsize),
                         axis.text.y = ggplot2::element_text(face = face, size = textsize),
                         plot.title = ggplot2::element_text(size = titlesize, face = face), 
                         plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                         legend.position = 'none')
        
      )
      
    }else{
      
      p <- ggplot2::ggplot(data = tsnedat, mapping = ggplot2::aes(x = tSNE1, 
                                                                  y = tSNE2, 
                                                                  color = Response))
      
      print(
        
        p + ggplot2::geom_point(size = 2) + 
          ggplot2::xlab(xlabel) + 
          ggplot2::ylab(ylabel) + 
          ggplot2::ggtitle(paste0(plotlabel, titlesuffix), 
                           subtitle = paste0('(', subtitle, ')')) + 
          
          ggplot2::scale_color_manual(responsevarname, values = colors) + 
          ggplot2::theme_bw() + 
          ggplot2::theme(axis.title.y = ggplot2::element_text(size = titlesize, face = face),
                         axis.title.x = ggplot2::element_text(size = titlesize, face = face), 
                         axis.text.x = ggplot2::element_text(face = face, size = textsize),
                         axis.text.y = ggplot2::element_text(face = face, size = textsize),
                         plot.title = ggplot2::element_text(size = titlesize, face = face), 
                         plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                         legend.text = ggplot2::element_text(size = textsize, face = face), 
                         legend.title = ggplot2::element_text(size = textsize, face = face))
        
      )
      
    }
    
    
    
  }else{
    
    colors <- c('grey')
    
    p <- ggplot2::ggplot(data = tsnedat, mapping = ggplot2::aes(x = tSNE1, 
                                                                y = tSNE2))
    
    print(
      
      p + ggplot2::geom_point(size = 2, col = colors) + 
        ggplot2::xlab(xlabel) + 
        ggplot2::ylab(ylabel) + 
        ggplot2::ggtitle(paste0(plotlabel, titlesuffix), 
                         subtitle = paste0('(', subtitle, ')')) + 
        
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = titlesize, face = face),
                       axis.title.x = ggplot2::element_text(size = titlesize, face = face), 
                       axis.text.x = ggplot2::element_text(face = face, size = textsize),
                       axis.text.y = ggplot2::element_text(face = face, size = textsize),
                       plot.title = ggplot2::element_text(size = titlesize, face = face), 
                       plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                       legend.position = 'none')
      
    )
    
  }
  
  colnames(tsnedat)[c(1, 2)] <- c(xlabel, ylabel)
  
  return(tsnedat)
  
  
}

#betas is a matrix using features as rows and samples as columns
embeddingplots <- function(betas, 
                           pddat = NULL, 
                           responsevarname = NULL, 
                           responselevels = NULL, 
                           seednum = 2022, 
                           perplexity = 4, 
                           titlesize = 15, 
                           textsize = 13, 
                           annotextsize = 4, 
                           face = 'bold', 
                           titlesuffix = NULL, 
                           labelpoints = TRUE){
  
  colors <- c('grey')
  
  if(!is.null(pddat)){
    
    if(is.null(responsevarname) & ncol(pddat) >= 2){
      responsevarname <- names(pddat)[2]
    }
    
    if(!is.null(responsevarname)){
      
      if(!is.factor(pddat[[responsevarname]])){
        
        if(is.null(responselevels)){
          responselevels <- pddat[[responsevarname]]
          responselevels <- unique(responselevels)
          responselevels <- responselevels[order(responselevels)]
        }
        pddat$Response <- factor(pddat[[responsevarname]], 
                                 levels = responselevels, 
                                 ordered = TRUE)
      }
      
      
      if(!is.null(responselevels)){
        colors <- scales::hue_pal()(length(unique(responselevels)))
      }else{
        colors <- scales::hue_pal()(length(unique(pddat$Response)))
      }
      
    }
    
    betas <- betas[,pddat[,1]]
    
  }
  
  
  featurenum <- nrow(betas)
  samplenum <- ncol(betas)
  
  if(!is.null(titlesuffix)){
    
    titlesuffix <- paste0(' ', titlesuffix)
    
  }
  
  subtitle <- paste0(samplenum, ' samples; ', featurenum, ' features')
  
  tsnematrix <- as.matrix(t(betas))
  
  plotdat <- tsneplot(tsnematrix = tsnematrix, 
                      perplexity = perplexity, 
                      colors = colors, 
                      response = pddat$Response, 
                      titlesuffix = titlesuffix, 
                      subtitle = subtitle, 
                      responsevarname = responsevarname, 
                      titlesize = titlesize, 
                      face = face, 
                      textsize = textsize, 
                      labelpoints = labelpoints, 
                      embedding = 'pca', 
                      seednum = seednum)
  
  tsnedat <- tsneplot(tsnematrix = tsnematrix, 
                      perplexity = perplexity, 
                      colors = colors, 
                      response = pddat$Response, 
                      titlesuffix = titlesuffix, 
                      subtitle = subtitle, 
                      responsevarname = responsevarname, 
                      titlesize = titlesize, 
                      face = face, 
                      textsize = textsize, 
                      labelpoints = labelpoints, 
                      embedding = 'tsne', 
                      seednum = seednum)
  
  umapdat <- tsneplot(tsnematrix = tsnematrix, 
                      perplexity = perplexity, 
                      colors = colors, 
                      response = pddat$Response, 
                      titlesuffix = titlesuffix, 
                      subtitle = subtitle, 
                      responsevarname = responsevarname, 
                      titlesize = titlesize, 
                      face = face, 
                      textsize = textsize, 
                      labelpoints = labelpoints, 
                      embedding = 'umap', 
                      seednum = seednum)
  
  
  
  row.names(plotdat) <- row.names(tsnedat) <- row.names(umapdat) <- 
    pddat[,1]
  
  reslist <- list()
  reslist$PCA <- plotdat
  reslist$tSNE <- tsnedat
  reslist$UMAP <- umapdat
  
  return(reslist)
  
}

#Probe/Gene/DMR annotation####

#'Annotate Illumina methylation probes
#'
#'Annotate Illumina methylation probes
#'
#'@param platform The platform of the probes. Can be 27 (for 27k platform), 
#'  450 (for 450k platform), or 850 (for EPIC platform). Default is 450.
#'@param probes The names of the probes need to be annotated.
#'@return A data.frame with the probe annotation information, including the 
#'  genomic loci of the probes and their related CpG islands and genes.
#'@export
probeanno <- function(platform = 450, 
                      probes){
  
  if(platform == 27){
    annotation <- 'ilmn12.hg19'
    array <- 'IlluminaHumanMethylation27k'
  }else if(platform == 450){
    annotation <- 'ilmn12.hg19'
    array <- 'IlluminaHumanMethylation450k'
  }else if(platform == 850){
    annotation <- 'ilm10b4.hg19'
    array <- 'IlluminaHumanMethylationEPIC'
  }else{
    cat('The parameter `platform` should be provided a value from 27, 450, and 850\n')
    return(NULL)
  }
  
  annopackage <- paste0(array, 'anno.', annotation)
  
  if(!(annopackage %in% installed.packages()[,'Package'])){
    cat(paste0('Package ', annopackage, ' is needed to run this function\n'))
    return(NULL)
  }
  
  if(!('AnnotationDbi' %in% installed.packages()[,'Package'])){
    cat('Package AnnotationDbi is needed to run this function\n')
    return(NULL)
  }
  
  if(!('org.Hs.eg.db' %in% installed.packages()[,'Package'])){
    cat(paste0('Package org.Hs.eg.db is needed to run this function\n'))
    return(NULL)
  }
  
  
  if(platform == 27){
    
    probeinfo <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Other
    islandsinfo <- probeinfo[c("CPG_ISLAND_LOCATIONS", "CPG_ISLAND")]
    locinfo <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Locations
    
    selectedcol <- c("Symbol", "Distance_to_TSS")
    genecol <- 'Symbol'
    featurecol <- 'Distance_to_TSS'
    islandcol <- 'CPG_ISLAND'
    
  }else if(platform == 450){
    
    probeinfo <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
    islandsinfo <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC
    locinfo <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
    
    selectedcol <- c("UCSC_RefGene_Name", "UCSC_RefGene_Group")
    genecol <- 'UCSC_RefGene_Name'
    featurecol <- 'UCSC_RefGene_Group'
    islandcol <- 'Relation_to_Island'
  }else if(platform == 850){
    
    probeinfo <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other
    islandsinfo <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC
    locinfo <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
    
    selectedcol <- c("UCSC_RefGene_Name", "UCSC_RefGene_Group")
    genecol <- 'UCSC_RefGene_Name'
    featurecol <- 'UCSC_RefGene_Group'
    islandcol <- 'Relation_to_Island'
  }
  
  probes <- probes[probes %in% row.names(probeinfo)]
  if(length(probes) == 0){
    return(NULL)
  }
  
  
  probeinfo <- probeinfo[probes,]
  probeinfodata <- as.data.frame(probeinfo)
  
  probeinfodata <- probeinfodata[selectedcol]
  
  islandsinfo <- islandsinfo[probes,]
  islandsinfodata <- as.data.frame(islandsinfo)
  
  locinfo <- locinfo[probes,]
  locinfo <- as.data.frame(locinfo)
  
  probeinfodata <- cbind(locinfo, probeinfodata, islandsinfodata)
  probeinfodata$Probe <- row.names(probeinfodata)
  row.names(probeinfodata) <- 1:nrow(probeinfodata)
  
  if(platform == 27){
    probeinfodata$ENTREZID <- probeinfo$Gene_ID
  }
  
  
  orgnizeprobeinfo <- function(colvec){
    parseelement <- function(element){
      elementlist <- unlist(strsplit(x = element, split = ';', fixed = TRUE))
      
      if(length(elementlist) == 0){
        elementlist <- ''
      }
      
      return(elementlist)
    }
    collist <- lapply(colvec, parseelement)
    return(collist)
  }
  
  genenamelist <- orgnizeprobeinfo(colvec = probeinfodata[,genecol])
  featurelist <- orgnizeprobeinfo(colvec = probeinfodata[,featurecol])
  poslist <- orgnizeprobeinfo(colvec = probeinfodata[,islandcol])
  
  genenamelistlens <- unlist(lapply(X = genenamelist, FUN = length))
  featurelistlens <- unlist(lapply(X = featurelist, FUN = length))
  poslistlens <- unlist(lapply(X = poslist, FUN = length))
  
  if(sum(genenamelistlens != 1) == 0 & platform == 27){
    probeinfodata$ENTREZID <- gsub(pattern = 'GeneID:', replacement = '', 
                                   x = probeinfodata$ENTREZID, fixed = TRUE)
  }
  
  if(sum(genenamelistlens == 0) == 0){
    
    datnames <- colnames(probeinfodata)
    
    chrvec <- rep(probeinfodata[,1], times = genenamelistlens)
    locvec <- rep(probeinfodata[,2], times = genenamelistlens)
    strandvec <- rep(probeinfodata[,3], times = genenamelistlens)
    
    genenamevec <- unlist(genenamelist)
    featurevec <- unlist(featurelist)
    
    islandvec <- rep(probeinfodata[,6], times = genenamelistlens)
    posvec <- rep(probeinfodata[,7], times = genenamelistlens)
    probevec <- rep(probeinfodata[,8], times = genenamelistlens)
    
    if(platform == 27){
      
      geneidlist <- orgnizeprobeinfo(colvec = probeinfodata$ENTREZID)
      geneidvec <- unlist(geneidlist)
      geneidvec <- gsub(pattern = 'GeneID:', replacement = '', 
                        x = geneidvec, fixed = TRUE)
      
      probeinfodata <- tryCatch({
        data.frame(chrvec, locvec, strandvec, 
                   genenamevec, 
                   featurevec, 
                   islandvec, posvec, probevec, 
                   geneidvec, 
                   stringsAsFactors = FALSE)
      }, error = function(err){
        probeinfodata
      })
      
      names(probeinfodata) <- datnames[1:ncol(probeinfodata)]
      
      probeinfodata <- unique(probeinfodata)
      
    }else{
      
      suppressMessages(unigeneidvec <- 
                         AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, 
                                               keys = unique(genenamevec), 
                                               columns = 'ENTREZID', 
                                               keytype = 'SYMBOL'))
      
      names(unigeneidvec)[1] <- selectedcol[1]
      
      probeinfodata <- tryCatch({
        data.frame(chrvec, locvec, strandvec, 
                   genenamevec, 
                   featurevec, 
                   islandvec, posvec, probevec, 
                   stringsAsFactors = FALSE)
      }, error = function(err){
        probeinfodata
      })
      
      names(probeinfodata) <- datnames[1:ncol(probeinfodata)]
      probeinfodata <- unique(probeinfodata)
      
      if(sum(!(unique(probeinfodata[,selectedcol[1]]) %in% 
               unigeneidvec[,selectedcol[1]])) == 0){
        
        probeinfodata <- merge(probeinfodata, unigeneidvec, 
                               by = selectedcol[1], 
                               sort = FALSE)
        
      }
      
      probeinfodata <- unique(probeinfodata)
      
    }
    
    if('Distance_to_TSS' %in% datnames){
      
      probeinfodata$Distance_to_TSS <- as.integer(probeinfodata$Distance_to_TSS)
      
    }
    
    if('CPG_ISLAND' %in% datnames){
      
      probeinfodata$CPG_ISLAND <- as.logical(probeinfodata$CPG_ISLAND)
      
    }
    
  }
  
  if(platform == 27){
    
    resnames <- c('Probe', 'chr', 'pos', 'strand', 
                  'CPG_ISLAND_LOCATIONS', 'CPG_ISLAND', 
                  'Symbol', 'ENTREZID', 'Distance_to_TSS')
    
  }else{
    
    resnames <- c('Probe', 'chr', 'pos', 'strand', 
                  'Islands_Name', 'Relation_to_Island', 
                  'UCSC_RefGene_Name', 'ENTREZID', 'UCSC_RefGene_Group')
    
  }
  
  probeinfodata <- probeinfodata[,resnames[unique(c(1:7, 
                                                    (ncol(probeinfodata) - 1), 9))]]
  
  row.names(probeinfodata) <- 1:nrow(probeinfodata)
  
  
  return(probeinfodata)
  
}

summaryfeature <- function(dat, featurecolidx){
  
  if(!('plyr' %in% installed.packages()[,'Package'])){
    cat('Package plyr is needed to run this function\n')
    return(NULL)
  }
  
  calmean <- function(block){
    gene <- unique(block$gene)
    subblock <- block[-1]
    submean <- colMeans(subblock)
    submatrix <- data.frame(submean)
    submatrix <- t(submatrix)
    row.names(submatrix) <- gene
    submatrix <- as.data.frame(submatrix)
    return(submatrix)
    
  }
  
  features <- dat[,featurecolidx]
  featurefreqs <- table(features)
  unifeatures <- names(featurefreqs[featurefreqs == 1])
  mulfeatures <- names(featurefreqs[featurefreqs > 1])
  
  unipart <- dat[dat[,featurecolidx] %in% unifeatures,]
  mulpart <- dat[dat[,featurecolidx] %in% mulfeatures,]
  
  mulpart <- plyr::ddply(.data = mulpart, 
                         .variables = c(names(dat)[featurecolidx]), 
                         .fun = calmean)
  
  row.names(mulpart) <- mulpart[,featurecolidx]
  mulpart <- mulpart[-featurecolidx]
  
  row.names(unipart) <- unipart[,featurecolidx]
  unipart <- unipart[-featurecolidx]
  
  finaldat <- rbind(unipart, mulpart)
  finaldat <- as.matrix(finaldat)
  
  return(finaldat)
  
}

#'Summarize the methylation beta values of probes to genes
#'
#'Summarize the methylation beta values of probes to genes by averaging the 
#'  probes located closely to the TSS of a gene.
#'
#'@param betadat A matrix recording the beta values of methylation probes for 
#'  samples. Each column represents one sample and each row represents one 
#'  probe. The row names are the probe names while the column names should be 
#'  sample IDs.
#'@param platform The platform of the probes. Can be 27 (for 27k platform), 
#'  450 (for 450k platform), or 850 (for EPIC platform). Default is 450.
#'@param range27k A positive number or a vector with two positive numbers. If 
#'  the data is based on 27k platform, it is needed to define which probes 
#'  could be considered as related to a specific gene, and only the ones with 
#'  a distance to the TSS of a gene less than the maximum value and greater 
#'  than the minimum value of \code{range27k} will be considered as related to 
#'  the gene, and the beta values of these probes will be averaged to get the 
#'  gene beta value. If it is a single number, the probes with a distance less 
#'  than this number and greater than 0 will be attributed to a gene. Default 
#'  is 200.
#'@param group450k850k A vector or single string. If the data is based on 450k 
#'  or EPIC platform, this parameter is needed to define which probes could be 
#'  considered as related to a specific gene. Only the ones located in the 
#'  gene regions included in this parameter will be considered as belong to 
#'  the gene. Its value can be selected from "TSS200", "TSS1500", "1stExon", 
#'  "5'UTR", '3'UTR", and "Body". Default is the vector c("TSS200", "TSS1500", 
#'  "1stExon"), which means probes within these 3 regions of a gene will be 
#'  attributed to the gene and their beta values will be averaged to get the 
#'  gene beta value.
#'@param includemultimatch Several probes can be attributed to more than one 
#'  gene. If this parameter is TRUE, these probes will be involved into the 
#'  beta value calculation for all their related genes. Otherwise, they will 
#'  be discarded, so that the beta values of all the genes are averaged only 
#'  from their uniquely related probes. Default is FALSE.
#'@param checkbetaval When this parameter is set as TRUE, the function will 
#'  check if the data values are between 0 and 1 first, and if not, it will 
#'  drop the downstream computation and return NULL. Default is TRUE.
#'@return A matrix recording the summarized gene beta values for samples.
#'@examples
#'library(CWGCNA)
#'
#'betas <- system.file("extdata", "placentabetas.rds", package = "CWGCNA")
#'betas <- readRDS(betas)
#'
#'betasgene <- probestogenes(betadat = betas, 
#'  group450k850k = c("TSS200", "TSS1500", "1stExon"))
#'@export
probestogenes <- function(betadat, 
                          platform = 450, 
                          range27k = 200, 
                          group450k850k = c('TSS200', 'TSS1500', '1stExon'), 
                          includemultimatch = FALSE, 
                          checkbetaval = TRUE){
  
  if(checkbetaval == TRUE){
    if(min(betadat) < 0 | max(betadat) > 1){
      
      cat('Need methylation beta value to run this function and values < 0 or > 1 cannot be used\n')
      return(NULL)
      
    }
  }
  
  
  beforeprobes <- row.names(betadat)
  
  tssinfo <- probeanno(platform = platform, probes = beforeprobes)
  
  if(is.null(tssinfo)){
    return(NULL)
  }
  
  if(platform == 27){
    if(length(range27k) == 1){
      range27k <- c(0, range27k)
    }else{
      range27k <- c(min(range27k), max(range27k))
    }
    
    tssinfor <- tssinfo[!is.na(tssinfo$Distance_to_TSS),]
    tssinfor <- tssinfor[(tssinfor$Distance_to_TSS >= as.numeric(min(range27k)) & 
                            tssinfor$Distance_to_TSS < as.numeric(max(range27k))),]
    tssinfor <- tssinfor[,c('Probe', 
                            'Symbol', 'ENTREZID')]
  }else{
    tssinfor <- tssinfo[tssinfo$UCSC_RefGene_Group %in% group450k850k,]
    tssinfor <- tssinfor[,c('Probe', 
                            'UCSC_RefGene_Name', 'ENTREZID')]
  }
  
  tssinfor <- unique(tssinfor)
  
  probes <- tssinfor$Probe
  
  if(includemultimatch == FALSE){
    
    probes <- probes[probes %in% names(table(probes)[table(probes) == 1])]
    
  }
  
  tssinfor <- subset(tssinfor, Probe %in% probes)
  
  tssinfor <- tssinfor[order(tssinfor[,2], tssinfor[,3]),]
  row.names(tssinfor) <- 1:nrow(tssinfor)
  
  
  tssinfor$genename <- paste0(tssinfor[,2], '::', tssinfor[,3])
  tssinfor$genename <- gsub(pattern = '^::', replacement = '', 
                            x = tssinfor$genename)
  tssinfor$genename <- gsub(pattern = '::$', replacement = '', 
                            x = tssinfor$genename)
  tssinfor$genename <- gsub(pattern = '^NA::', replacement = '', 
                            x = tssinfor$genename)
  tssinfor$genename <- gsub(pattern = '::NA$', replacement = '', 
                            x = tssinfor$genename)
  
  
  betadat <- betadat[tssinfor$Probe, , drop = FALSE]
  betadat <- as.data.frame(betadat, stringsAsFactors = FALSE)
  betadat$genename <- tssinfor$genename
  betadat <- betadat[,c(ncol(betadat), 1:(ncol(betadat) - 1))]
  
  betadat <- summaryfeature(dat = betadat, featurecolidx = 1)
  
  if(is.null(betadat)){
    return(NULL)
  }
  
  genenames <- row.names(betadat)
  geneorders <- order(genenames)
  betadat <- betadat[geneorders, , drop = FALSE]
  
  return(betadat)
  
}

#'Annotate genes
#'
#'Annotate gene coordinates and functions
#'
#'@param genesymbols A vector containing the symbols of genes to be annotated. 
#'  Default is NULL.
#'@param geneentrezs A vector with the ENTREZ IDs of genes to be annotated. 
#'  If \code{genesymbols} is NULL, this value will be used. Default is NULL. 
#'@return A data frame with the gene annotation information.
#'@export
geneanno <- function(genesymbols = NULL, 
                     geneentrezs = NULL){
  
  
  summaryannotation <- get('summaryannotation')
  #summaryannotation <- 
  #  readRDS('C:/Users/yuabr/Desktop/Transfer/codetransfer/emethyl/files/summaryannotation.rds')
  
  
  if(!is.null(genesymbols)){
    
    genesymbols <- genesymbols[!is.na(genesymbols)]
    genesymbols <- toupper(genesymbols)
    
    totalprobegeneannotation <- subset(summaryannotation, 
                                       SYMBOL %in% genesymbols | 
                                         preferred_name %in% genesymbols)
    
    part1 <- subset(totalprobegeneannotation, 
                    SYMBOL %in% genesymbols)
    part1$input <- part1$SYMBOL
    
    part2 <- subset(totalprobegeneannotation, 
                    !(SYMBOL %in% genesymbols))
    part2$input <- part2$preferred_name
    
    totalprobegeneannotation <- rbind(part1, part2)
    totalprobegeneannotation <- unique(totalprobegeneannotation)
    
    
    others <- genesymbols[!(genesymbols %in% 
                              c(totalprobegeneannotation$SYMBOL, 
                                totalprobegeneannotation$preferred_name))]
    
    
  }else if(!is.null(geneentrezs)){
    
    geneentrezs <- as.character(geneentrezs)
    geneentrezs <- geneentrezs[!is.na(geneentrezs)]
    totalprobegeneannotation <- subset(summaryannotation, 
                                       ENTREZID %in% geneentrezs)
    totalprobegeneannotation$input <- totalprobegeneannotation$ENTREZID
    totalprobegeneannotation <- unique(totalprobegeneannotation)
    
    others <- geneentrezs[!(geneentrezs %in% totalprobegeneannotation$ENTREZID)]
    
  }else{
    
    return(NULL)
  }
  
  if(length(others) > 0){
    others <- data.frame(ENTREZID = NA, 
                         SYMBOL = NA, 
                         chr = NA, 
                         start = NA, 
                         end = NA, 
                         strand = NA, 
                         preferred_name = NA, 
                         UniProt = NA, 
                         UniProt.name = NA, 
                         protein.size = NA, 
                         annotation = NA, 
                         input = others, 
                         stringsAsFactors = FALSE)
    totalprobegeneannotation <- rbind(totalprobegeneannotation, others)
    
  }
  
  
  
  totalprobegeneannotation <- unique(totalprobegeneannotation)
  totalprobegeneannotation <- totalprobegeneannotation[
    c('input', 'ENTREZID', 'SYMBOL', 'chr', 'start', 'end', 'strand', 
      'preferred_name', 'UniProt', 'UniProt.name', 'protein.size', 'annotation')
  ]
  
  row.names(totalprobegeneannotation) <- 1:nrow(totalprobegeneannotation)
  
  
  return(totalprobegeneannotation)
  
  
}

coverredTSS <- function(fragments, tssradius = NULL, 
                        ignorestrand = TRUE){
  
  if(!('GenomicFeatures' %in% installed.packages()[,'Package'])){
    cat('Package GenomicFeatures is needed to run this function\n')
    return(NULL)
  }
  
  if(!('TxDb.Hsapiens.UCSC.hg19.knownGene' %in% installed.packages()[,'Package'])){
    cat('Package TxDb.Hsapiens.UCSC.hg19.knownGene is needed to run this function\n')
    return(NULL)
  }
  
  if(!('IRanges' %in% installed.packages()[,'Package'])){
    cat('Package IRanges is needed to run this function\n')
    return(NULL)
  }
  
  if(!('GenomicRanges' %in% installed.packages()[,'Package'])){
    cat('Package GenomicRanges is needed to run this function\n')
    return(NULL)
  }
  
  if(!('AnnotationDbi' %in% installed.packages()[,'Package'])){
    cat('Package AnnotationDbi is needed to run this function\n')
    return(NULL)
  }
  
  if(!('org.Hs.eg.db' %in% installed.packages()[,'Package'])){
    cat('Package org.Hs.eg.db is needed to run this function\n')
    return(NULL)
  }
  
  genecoords <- suppressMessages(
    GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene))
  generanges <- genecoords@ranges
  
  geneids <- genecoords$gene_id
  geneseqs <- genecoords@seqnames
  genestrands <- genecoords@strand
  
  genetssp <- generanges@start[as.vector(genestrands == '+')]
  genetssm <- generanges@start[as.vector(genestrands == '-')] + 
    generanges@width[as.vector(genestrands == '-')] - 1
  
  tss <- rep(0, length(genecoords))
  tss[as.vector(genestrands == '+')] <- genetssp
  tss[as.vector(genestrands == '-')] <- genetssm
  
  tssranges <- IRanges::IRanges(start = tss, width = 1)
  tsscoords <- GenomicRanges::GRanges(seqnames = geneseqs, 
                                      ranges = tssranges, 
                                      strand = genestrands, 
                                      gene_id = geneids)
  
  
  fragchrs <- gsub(pattern = ':.*$', replacement = '', x = fragments)
  fragcoords <- gsub(pattern = '^chr.*:', replacement = '', x = fragments)
  fragstarts <- gsub(pattern = '-.*$', replacement = '', x = fragcoords)
  fragends <- gsub(pattern = '^.*-', replacement = '', x = fragcoords)
  fragstarts <- as.numeric(fragstarts)
  fragends <-as.numeric(fragends)
  
  fragranges <- IRanges::IRanges(start = fragstarts, end = fragends)
  fragranges <- GenomicRanges::GRanges(seqnames = fragchrs, 
                                       ranges = fragranges, strand = '*', 
                                       fragmentname = fragments)
  
  dis <- GenomicRanges::distanceToNearest(x = tsscoords, 
                                          subject = fragranges, 
                                          ignore.strand = ignorestrand)
  disvec <- dis@elementMetadata$distance
  
  rangeinfo <- fragranges[dis@to,]
  geneinfo <- genecoords[dis@from,]
  rangeinfo <- as.data.frame(rangeinfo)
  names(geneinfo) <- 1:length(geneinfo)
  geneinfo <- as.data.frame(geneinfo)
  names(rangeinfo) <- paste('range', names(rangeinfo), sep = '_')
  names(geneinfo) <- paste('gene', names(geneinfo), sep = '_')
  overlap <- cbind(rangeinfo, geneinfo)
  overlap$tssdistance <- disvec
  overlap <- unique(overlap)
  
  
  suppressMessages(genesyms <- 
                     AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, 
                                           keys = overlap$gene_gene_id, 
                                           columns = 'SYMBOL', 
                                           keytype = 'ENTREZID'))
  overlap <- cbind(overlap, genesyms)
  
  generes <- data.frame(frag = overlap$range_fragmentname, 
                        seqnames = overlap$gene_seqnames, 
                        start = overlap$gene_start, 
                        end = overlap$gene_end, 
                        width = overlap$gene_width, 
                        strand = overlap$gene_strand, 
                        geneid = overlap$ENTREZID, 
                        genename = overlap$SYMBOL, 
                        TSSdistance = overlap$tssdistance, 
                        stringsAsFactors = FALSE)
  
  if(!is.null(tssradius)){
    generes <- subset(generes, TSSdistance <= tssradius)
  }
  
  return(generes)
}

#Limma functions####

#betas is a matrix using features as rows and samples as columns
orgnas <- function(betas){
  
  if(sum(is.infinite(betas)) > 0){
    
    betas[is.infinite(betas)] <- NA
    
  }
  
  if(sum(!complete.cases(betas)) > 0){
    
    probeids <- row.names(betas)
    vals <- impute::impute.knn(betas)
    betas <- vals$data
    
    betas <- betas[complete.cases(betas),]
    betas <- betas[is.finite(rowSums(betas)),]
    
  }
  
  return(betas)
  
}

#'Find DNA methylation regions (DMRs) relevant to specific phenotypic variable 
#'
#'Cluster the probes into DNA methylation regions (DMRs) and identify the ones 
#'  relevant to specific phenotypic variable.
#'
#'@param dat A matrix with the beta values of methylation probes for samples. 
#'  Each column represents one sample and each row represents one probe. The 
#'  row names are the probe names and the column names should be sample IDs. 
#'@param pddat Meta data frame. The first column should be the sample names 
#'  with a column name as "sampleid", and the remaining should be phenotypic 
#'  variables.
#'@param responsevarname The column name of the response variable in the data 
#'  frame provided to \code{pddat}. It is the variable of interest. Default 
#'  is NULL, which means the second column of \code{pddat} is the response. 
#'  This response variable should be a binary or continuous variable.
#'@param responselevels If the response variable is binary, it will be treated 
#'  as a factor, and this parameter is needed to define the factor level, so 
#'  it should be a vector defining this level. Only needed if the response in 
#'  \code{pddat} is a character vector but not a factor. In this case, it will 
#'  be converted to a factor by the function, and \code{responselevels} will 
#'  be used to define the factor level. Default is NULL, meaning this level 
#'  will be set by the function following the character order of the elements 
#'  in the response variable.
#'@param confoundings The column name of the confounding factors in the data 
#'  frame provided to \code{pddat}. Default is NULL, which means there is no 
#'  confounding factor in the data. 
#'@param platform The platform of the methylation data provided to \code{dat}. 
#'  Can be set as 450 (for 450k platform), or 850 (for EPIC platform). Default 
#'  is 450.
#'@param maxgap An integer indicating the cutoff of probe-probe distance when 
#'  clustering the DNA methylation probes in \code{dat} to DMRs. In the DNA 
#'  metylation data, if the distance between 2 neighbor probes is less than 
#'  \code{maxgap}, they will be clustered into the same DMR. Then, this large 
#'  DMR will be further divided into smaller DMRs according to the values of 
#'  its probes. For a set of consecutive probes in the large DMR, if all the 
#'  probe values were greater than a cutoff value, which is defined by the 
#'  function from a bootstrapping distribution, this set will be a small DMR. 
#'  On the other hand, if all the probes in a set is less than the minus value 
#'  of the cutoff, it will be defined as another small DMR. These small DMRs 
#'  will be returned as DMRs significantly related to the response variable.
#'@param Bvalue An integer denoting the number of sampling when computing the 
#'  bootstrapping distributions. This defaults to 100.
#'@param randomseed Random seed for the bootstrapping sampling. Default value 
#'  is 2022.
#'@param threads Number of threads need for parallelization. Default is 1.
#'@return A data frame recording the DMRs related to the response variable, 
#'  with their information such as genomic coordinates, bootstrappoing p-val, 
#'  etc.
#'@export
sigdmrs <- function(dat, 
                    pddat, 
                    responsevarname = NULL, 
                    responselevels = NULL, 
                    confoundings = NULL, 
                    platform = 450, 
                    maxgap = 300, 
                    Bvalue = 100, 
                    randomseed = 2022, 
                    threads = 1){
  
  
  if(is.null(responsevarname) & ncol(pddat) >= 2){
    responsevarname <- names(pddat)[2]
  }
  
  if(is.null(responsevarname)){
    return(NULL)
  }
  
  #library(limma)
  
  simpddat <- pddat[c('sampleid', responsevarname, confoundings)]
  simpddat <- simpddat[complete.cases(simpddat),]
  
  if(!is.factor(pddat[[responsevarname]]) & !is.numeric(pddat[[responsevarname]])){
    
    if(is.null(responselevels)){
      responselevels <- simpddat[[responsevarname]]
      responselevels <- unique(responselevels)
      responselevels <- responselevels[order(responselevels)]
    }
    simpddat$Response <- factor(simpddat[[responsevarname]], 
                                levels = responselevels, 
                                ordered = TRUE)
  }else{
    
    simpddat$Response <- simpddat[[responsevarname]]
    
  }
  
  
  vars <- paste(c('Response', confoundings), collapse = ' + ')
  
  design <- model.matrix(as.formula(paste0('~ ', vars)), 
                         data = simpddat)
  
  ###########################
  ## following parameters are specifically for Bumphunter method.
  #Convert probe value to gene value manually
  .default.450k.annotation <- "ilmn12.hg19"
  .default.epic.annotation <- "ilm10b4.hg19"
  
  if(platform == 450){
    array <- "IlluminaHumanMethylation450k"
    annotation <- .default.450k.annotation
    loci <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
  }else{
    array <- "IlluminaHumanMethylationEPIC"
    annotation <- .default.epic.annotation
    loci <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
  }
  
  annopackage <- paste0(array, 'anno.', annotation)
  
  locus <- intersect(row.names(dat), row.names(loci))
  dat <- dat[locus,]
  locus <- loci[locus,]
  
  cl <- bumphunter::clusterMaker(locus$chr, locus$pos, maxGap = maxgap)
  
  if(threads == 1){
    
    set.seed(randomseed)
    tab <- bumphunter::bumphunter(dat, 
                                  design, 
                                  locus$chr, 
                                  locus$pos, 
                                  cl, 
                                  coef = 2, 
                                  nullMethod = 'bootstrap', 
                                  cutoff = NULL, 
                                  pickCutoff = TRUE, 
                                  smooth = TRUE, smoothFunction = bumphunter::loessByCluster,
                                  useWeights = FALSE, permutations = NULL,
                                  B = Bvalue)
    #The function `regionFinder` is used in the final steps of `bumphunter`. It 
    #simply finds regions that are above a certain threshold.
    #The clusters were first defined by `cl`, which is the result of `clusterMaker`. 
    #Then, within each cluster, `regionFinder` further divides this cluster INTO 
    #SMALLER REGIONS (DMRs) according to their values (> cutoff or < -cutoff).
    #For example, within a `cl` cluster, if all the probes in a sub-region (a 
    #series of consecutive probes) > cutoff, while both the up and downstream sub-
    #regions are not > cutoff (are < -cutoff or almost = 0). Then, this sub-region 
    #and the up and downstream sub-regions will be separated into 3 SMALLER REGIONS. 
    #The middle one (the current region) is > cutoff, while the up and downstream 
    #regions will be < -cutoff or almost = 0. The final DMRs defined by `bumphunter` 
    #are these SMALLER REGIONS after further division, not the original clusters 
    #defined by `cl`.
    #If a cluster is divided by `regionFinder` into several smaller regions, each 
    #smaller region will account for one row of the result table of `regionFinder`, 
    #and the column `cluster` in it records the original cluster IDs of the smaller 
    #regions, so if 2 smaller regions are from the same `cl` cluster, their values 
    #in the column `cluster` will be the same, although they are 2 different 
    #regions for 2 different rows in the table.
    #The smaller regions > cutoff or < -cutoff will be returned by `bumphunter` 
    #as significant DMRs, other regions not passing the cutoff will not be returned. 
    #The cutoff can be set by the parameter `cutoff` of `bumphunter` definitely, 
    #or by setting the parameter `pickCutoff` as TRUE and then letting `bumphunter` 
    #itself define a cutoff with the permutation distribution (This is the case here).
    
    
    #Estimate regions for which a genomic profile deviates from its baseline value. 
    #Originally implemented to detect differentially methylated genomic regions 
    #between two populations.
    #y is the normalized methy matrix, one row is one methy locus, one column is 
    #one sample
    #X is the pddat matrix, the first column is the interception, the second column is 
    #the variable case-control
    #coef	An integer denoting the column of the design matrix containing the 
    #     covariate of interest. The hunt for bumps will be only be done for the 
    #     estimate of this coefficient.
    #nullMethod	Method used to generate null candidate regions, must be one of 
    #           'bootstrap' or 'permutation' (defaults to 'permutation'). 
    #           However, if covariates in addition to the outcome of interest are 
    #           included in the design matrix (ncol(design)>2), the 'permutation' 
    #           approach is not recommended.
    #B	An integer denoting the number of resamples to use when computing null 
    #   distributions. This defaults to 0.
    
    #The main output is a table of candidate regions with permutation or 
    #bootstrap-based family-wide error rates (FWER) and p-values assigned.
    
    
  }else{
    
    cores <- parallel::detectCores()
    cores <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cores)
    
    set.seed(randomseed)
    tab <- bumphunter::bumphunter(dat, 
                                  design, 
                                  locus$chr, 
                                  locus$pos, 
                                  cl, 
                                  coef = 2, 
                                  nullMethod = 'bootstrap', 
                                  cutoff = NULL, pickCutoff = TRUE, 
                                  smooth = TRUE, smoothFunction = bumphunter::loessByCluster,
                                  useWeights = FALSE, permutations = NULL,
                                  B = Bvalue)
    
    parallel::stopCluster(cores)
    
    unregister_dopar()
    
    
  }
  
  
  bumptable <- tab$table
  
  if(is.na(bumptable) | is.null(bumptable)){
    return(NULL)
  }
  
  names(bumptable)[1] <- 'seqnames'
  bmp <- bumptable
  
  bmp$width <- abs(bmp$end - bmp$start) + 1
  bmp$strand <- factor('*')
  bmpnames <- c('seqnames', 'start', 'end', 'width', 'strand', 'value', 'area', 'cluster', 'indexStart', 
                'indexEnd', 'L', 'clusterL', 'p.value', 'fwer', 'p.valueArea', 'fwerArea')
  bmp <- bmp[bmpnames]
  bmp <- bmp[order(bmp$p.value),]
  rowindex <- paste0('DMR_', 1:nrow(bmp))
  row.names(bmp) <- rowindex
  
  return(bmp)
  
  
}

#dat is a matrix using features as rows and samples as columns
#pds data.frame must contain a column named "sampleid"
diffsites <- function(dats, 
                      pddats, 
                      i = 1, 
                      responsevarname = NULL, 
                      responselevels = NULL, 
                      confoundings = NULL){
  
  
  dat <- dats[[i]]
  pddat <- pddats[[i]]
  
  if(is.null(responsevarname) & ncol(pddat) >= 2){
    responsevarname <- names(pddat)[2]
  }
  
  if(is.null(responsevarname)){
    return(NULL)
  }
  
  #library(limma)
  
  simpddat <- pddat[c('sampleid', responsevarname, confoundings)]
  simpddat <- simpddat[complete.cases(simpddat),]
  
  if(!is.factor(pddat[[responsevarname]]) & !is.numeric(pddat[[responsevarname]])){
    
    if(is.null(responselevels)){
      responselevels <- simpddat[[responsevarname]]
      responselevels <- unique(responselevels)
      responselevels <- responselevels[order(responselevels)]
    }
    simpddat$Response <- factor(simpddat[[responsevarname]], 
                                levels = responselevels, 
                                ordered = TRUE)
  }else{
    
    simpddat$Response <- simpddat[[responsevarname]]
    
  }
  
  
  vars <- paste(c('Response', confoundings), collapse = ' + ')
  
  design <- model.matrix(as.formula(paste0('~ ', vars)), 
                         data = simpddat)
  
  betasub <- dat[,simpddat$sampleid, drop = FALSE]
  
  
  
  fit1 <- limma::lmFit(betasub, design)
  fit1 <- limma::eBayes(fit1)
  
  allg.limma <- limma::topTable(fit1, coef = 2, n = dim(fit1)[1])
  
  if((nrow(allg.limma) == 1) & (nrow(betasub) == 1)){
    row.names(allg.limma) <- row.names(betasub)
  }
  
  return(allg.limma)
  
}

overlappedprobes <- function(totalprobes, 
                             removereg, 
                             platform = 450, 
                             ignorestrand = TRUE){
  
  if(platform == 850){
    chrinfo <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
  }else if(platform == 450){
    chrinfo <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
  }else if(platform == 27){
    chrinfo <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Locations
  }
  
  
  if(class(removereg)[1] == 'GRanges'){
    
    removedmrgr <- removereg
    
  }else if(class(removereg)[1] == 'data.frame'){
    
    removedmrir <- IRanges::IRanges(start = removereg$start, end = removereg$end)
    removedmrgr <- GenomicRanges::GRanges(seqnames = removereg$seqnames, 
                                          ranges = removedmrir, 
                                          strand = removereg$strand)
    
  }
  
  removedmrgr <- GenomicRanges::reduce(removedmrgr)
  
  
  
  totalprobes <- intersect(totalprobes, row.names(chrinfo))
  totalprobes <- chrinfo[totalprobes,]
  
  totalprobesir <- IRanges::IRanges(start = totalprobes$pos, width = 1)
  totalprobesgr <- GenomicRanges::GRanges(seqnames = totalprobes$chr, 
                                          ranges = totalprobesir, 
                                          strand = totalprobes$strand)
  
  
  intprobes <- GenomicRanges::findOverlaps(query = totalprobesgr, 
                                           subject = removedmrgr, 
                                           type = c('within'), 
                                           ignore.strand = ignorestrand)
  
  overprobes <- totalprobes[intprobes@from,]
  overdmrs <- removedmrgr[intprobes@to,]
  
  overlapmapping <- overprobes
  
  if(nrow(overlapmapping) > 0){
    
    overdmrs <- paste0(overdmrs@seqnames, ':', overdmrs@ranges@start, '-', 
                       overdmrs@ranges@start + overdmrs@ranges@width - 1)
    
    overlapmapping$frag <- overdmrs
    
    overlapmapping <- as.data.frame(overlapmapping, stringsAsFactors = FALSE)
    
    res <- overlapmapping
  }else{
    res <- NULL
  }
  
  return(res)
  
}

convertgenenames <- function(totalgenes){
  
  formergenenames <- gsub(pattern = '::.*$', replacement = '', x = totalgenes)
  latergenenames <- gsub(pattern = '^.*::', replacement = '', x = totalgenes)
  
  totalgenenames <- unique(c(formergenenames, latergenenames))
  totalgenenames <- totalgenenames[totalgenenames != '']
  
  suppressWarnings(geneids <- as.numeric(totalgenenames))
  genesyms <- totalgenenames[is.na(geneids)]
  geneids <- geneids[!is.na(geneids)]
  
  if(length(geneids) == 0){
    geneids <- NULL
  }
  
  if(length(genesyms) > 0){
    
    suppressMessages(
      
      convertedsyms <- tryCatch({
        
        AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, 
                              keys = genesyms, 
                              columns = 'SYMBOL', 
                              keytype = 'ENSEMBL')
        
      }, error = function(err){
        
        NULL
        
      })
      
    )
    
    genesyms <- unique(c(genesyms, convertedsyms))
    
  }else{
    genesyms <- NULL
  }
  
  genenames <- list(genesyms = genesyms, 
                    geneids = geneids)
  
  return(genenames)
  
}

ensemblidx <- function(syms, 
                       genenames, 
                       genesym){
  
  syms <- syms[complete.cases(syms),]
  syms <- unique(syms)
  dupsyms <- unique(syms$ENSEMBL[duplicated(syms$ENSEMBL)])
  syms <- subset(syms, !(ENSEMBL %in% dupsyms))
  dupsyms <- unique(syms$SYMBOL[duplicated(syms$SYMBOL)])
  syms <- subset(syms, !(SYMBOL %in% dupsyms))
  
  syms <- unique(syms)
  
  symsidx <- match(x = syms$ENSEMBL, table = genenames)
  
  idx <- rep(FALSE, length(genenames))
  idx[symsidx] <- syms$SYMBOL %in% genesym
  
  return(idx)
  
}

idxgenenames <- function(totalgenes = row.names(limmares), 
                         genesym = removedfeatures$genesym, 
                         geneid = removedfeatures$geneid){
  
  formergenenames <- gsub(pattern = '::.*$', replacement = '', x = totalgenes)
  latergenenames <- gsub(pattern = '^.*::', replacement = '', x = totalgenes)
  
  #match genesyms
  idx1 <- formergenenames %in% genesym
  idx2 <- latergenenames %in% genesym
  
  #match geneids
  idx3 <- formergenenames %in% geneid
  idx4 <- latergenenames %in% geneid
  
  #match ensemble
  
  suppressMessages(
    
    formersyms <- tryCatch({
      
      AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, 
                            keys = formergenenames, 
                            columns = 'SYMBOL', 
                            keytype = 'ENSEMBL')
      
    }, error = function(err){
      
      NULL
      
    })
    
  )
  
  suppressMessages(
    
    latersyms <- tryCatch({
      
      AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, 
                            keys = latergenenames, 
                            columns = 'SYMBOL', 
                            keytype = 'ENSEMBL')
      
    }, error = function(err){
      
      NULL
      
    })
    
  )
  
  
  if(!is.null(formersyms)){
    
    idx5 <- ensemblidx(syms = formersyms, 
                       genenames = formergenenames, 
                       genesym = genesym)
    
  }else{
    idx5 <- rep(FALSE, length(formergenenames))
  }
  
  if(!is.null(latersyms)){
    
    idx6 <- ensemblidx(syms = latersyms, 
                       genenames = latergenenames, 
                       genesym = genesym)
    
  }else{
    idx6 <- rep(FALSE, length(latergenenames))
  }
  
  idx <- idx1 | idx2 | idx3 | idx4 | idx5 | idx6
  
  return(idx)
  
}

overlappedgenes <- function(totalgenes, 
                            removereg, 
                            tssradius = 1500, 
                            ignorestrand = TRUE){
  
  
  if(class(removereg)[1] == 'GRanges'){
    
    removedmrgr <- removereg
    
  }else if(class(removereg)[1] == 'data.frame'){
    
    removedmrir <- IRanges::IRanges(start = removereg$start, end = removereg$end)
    removedmrgr <- GenomicRanges::GRanges(seqnames = removereg$seqnames, 
                                          ranges = removedmrir, 
                                          strand = removereg$strand)
  }
  
  removedmrgr <- GenomicRanges::reduce(removedmrgr)
  
  
  
  dmrfrags <- paste0(removedmrgr@seqnames, ':', 
                     removedmrgr@ranges@start, '-', 
                     removedmrgr@ranges@start + removedmrgr@ranges@width - 1)
  
  genemapping <- coverredTSS(fragments = dmrfrags, tssradius = tssradius, 
                             ignorestrand = ignorestrand)
  
  if(is.null(genemapping)){
    return(NULL)
  }
  
  genenames <- convertgenenames(totalgenes = totalgenes)
  
  if(('genesyms' %in% names(genenames)) & ('geneids' %in% names(genenames))){
    
    submappings <- subset(genemapping, (geneid %in% genenames$geneids) | 
                            (genename %in% genenames$genesyms))
    
  }else if('genesyms' %in% names(genenames)){
    
    submappings <- subset(genemapping, genename %in% genenames$genesyms)
    
  }else if('geneids' %in% names(genenames)){
    
    submappings <- subset(genemapping, genename %in% genenames$geneids)
    
  }else{
    submappings <- NULL
  }
  
  res <- NULL
  
  if(!is.null(submappings)){
    
    if(nrow(submappings) > 0){
      
      submappings <- submappings[,c('genename', 'geneid', 'seqnames', 
                                    'start', 'end', 'width', 'strand', 
                                    'frag', 'TSSdistance')]
      colnames(submappings)[c(1, 2, 3, 
                              4, 5, 6, 7)] <- c('genesym', 'geneid', 'chr', 
                                                'gene_start', 'gene_end', 'gene_width', 
                                                'gene_strand')
      row.names(submappings) <- 1:nrow(submappings)
      
      res <- submappings
    }
  }
  
  return(res)
  
}

overlappeddmrs <- function(totaldmrs, 
                           removereg, 
                           ignorestrand = TRUE){
  
  if(class(removereg)[1] == 'GRanges'){
    
    removedmrgr <- removereg
    
  }else if(class(removereg)[1] == 'data.frame'){
    
    removedmrir <- IRanges::IRanges(start = removereg$start, end = removereg$end)
    removedmrgr <- GenomicRanges::GRanges(seqnames = removereg$seqnames, 
                                          ranges = removedmrir, 
                                          strand = removereg$strand)
  }
  
  removedmrgr <- GenomicRanges::reduce(removedmrgr)
  
  dmrstart <- totaldmrs[,grepl(pattern = 'start', x = names(totaldmrs))]
  dmrend <- totaldmrs[,grepl(pattern = 'end', x = names(totaldmrs))]
  dmrchr <- totaldmrs[,grepl(pattern = 'chr', x = names(totaldmrs))]
  dmrir <- IRanges::IRanges(start = dmrstart, end = dmrend)
  dmrgr <- GenomicRanges::GRanges(seqnames = dmrchr, 
                                  ranges = dmrir)
  
  
  intdmrs <- GenomicRanges::findOverlaps(query = dmrgr, subject = removedmrgr, type = c('any'), 
                                         ignore.strand = ignorestrand)
  
  overtotaldmrs <- totaldmrs[intdmrs@from,]
  overremoveddmrs <- removedmrgr[intdmrs@to,]
  
  overlapmapping <- overtotaldmrs
  
  if(nrow(overlapmapping) > 0){
    
    overremoveddmrs <- paste0(overremoveddmrs@seqnames, ':', overremoveddmrs@ranges@start, '-', 
                              overremoveddmrs@ranges@start + overremoveddmrs@ranges@width - 1)
    
    overlapmapping$frag <- overremoveddmrs
    
    overlapmapping <- as.data.frame(overlapmapping, stringsAsFactors = FALSE)
    
    res <- overlapmapping
  }else{
    res <- NULL
  }
  
  return(res)
  
}

#dat is a matrix using features as rows and samples as columns
#pddat data.frame must contain a column named "sampleid"
balancesampling <- function(dat,
                            pddat, 
                            responsename, 
                            nround = 10,
                            seed = 2022, 
                            k = 5, 
                            threads = 1, 
                            samplingmethod = 'updn'){
  
  dat <- t(dat)
  
  if(is.null(responsename) & ncol(pddat) >= 2){
    responsename <- names(pddat)[2]
  }
  
  if(is.null(responsename)){
    return(NULL)
  }
  
  
  
  sampleidmapping <- pddat[,c('sampleid', responsename)]
  names(sampleidmapping) <- c('sampleid', 'responsetype')
  
  dat <- dat[sampleidmapping$sampleid, , drop = FALSE]
  
  responsetypes <- unique(sampleidmapping$responsetype)
  
  
  
  samplingpds <- function(samplingsamples, 
                          pddat){
    
    samples <- gsub(pattern = '\\.[0-9]*$', replacement = '', x = samplingsamples)
    samplingpds <- pddat[match(samples, pddat$sampleid), , drop = FALSE]
    samplingpds$sampleid <- samplingsamples
    row.names(samplingpds) <- 1:nrow(samplingpds)
    
    return(samplingpds)
    
  }
  
  createnames <- function(rawnames){
    
    reversenames <- FALSE
    if(sum(grepl(pattern = '-', x = rawnames)) > 0){
      
      rawnames <- gsub(pattern = '-', replacement = '___', x = rawnames)
      reversenames <- TRUE
    }
    
    
    newnames <- make.names(names = rawnames, unique = TRUE)
    suffix <- gsub(pattern = '^.*\\.', replacement = '', x = newnames)
    suffix[suffix %in% rawnames] <- ''
    suffix <- paste0('.', suffix)
    suffix[suffix == '.'] <- ''
    newnames <- paste0(rawnames, suffix)
    
    if(reversenames == TRUE){
      newnames <- gsub(pattern = '___', replacement = '-', x = newnames)
    }
    
    return(newnames)
    
  }
  
  singleroundsampling <- function(responsetypes,
                                  sampleidmapping,
                                  i = 2022,
                                  dat, 
                                  k = 5, 
                                  samplingmethod = 'updn'){
    
    if(samplingmethod == 'up'){
      
      targetnum <- max(table(sampleidmapping$responsetype))
      
    }else if(samplingmethod == 'dn'){
      
      targetnum <- min(table(sampleidmapping$responsetype))
      
    }else{
      
      targetnum <- ceiling(nrow(sampleidmapping)/length(responsetypes))
      
    }
    
    
    
    j <- 1
    for(j in 1:length(responsetypes)){
      
      ResponseType <- responsetypes[j]
      sub <- subset(sampleidmapping, responsetype == ResponseType)
      
      subnum <- nrow(sub)
      
      if(subnum >= targetnum & subnum >= 1){
        
        sampleids <- sub$sampleid
        
        set.seed(i)
        samples <- sample(x = sampleids,
                          size = targetnum,
                          replace = FALSE)
        #Set replace to FALSE (random sampling), not TRUE (bootstrapping), 
        #because if using bootstrapping, the downstream limma diff feature 
        #analysis will return many false positive differential features 
        #due to the repeated features introduced by this bootstrapping step
        
        datsub <- dat[samples, , drop = FALSE]
        
        row.names(datsub) <- createnames(rawnames = row.names(datsub))
        
        
      }else if(subnum < targetnum & subnum > 1){
        
        subvar <- dat[sub$sampleid, , drop = FALSE]
        
        knnres <- dbscan::kNN(x = subvar, k = min(k, nrow(subvar) - 1))
        knnres <- knnres$id
        knnresidx <- seq(1, nrow(knnres)*ncol(knnres), by = 1)
        knnresidxmat <- matrix(data = knnresidx, 
                               nrow = nrow(knnres), 
                               byrow = FALSE)
        row.names(knnresidxmat) <- 1:nrow(knnresidxmat)
        colnames(knnresidxmat) <- 1:ncol(knnresidxmat)
        
        set.seed(i)
        sampledidx <- sample(x = knnresidx, size = targetnum - subnum, replace = TRUE)
        
        smoteidx <- do.call(rbind, 
                            lapply(X = sampledidx, 
                                   FUN = function(x){which(knnresidxmat == x, arr.ind = TRUE)}))
        row.names(smoteidx) <- 1:nrow(smoteidx)
        
        set.seed(i)
        zetas <- runif(nrow(smoteidx))
        
        
        synthesizeddat <- lapply(X = seq(1:nrow(smoteidx)), 
                                 FUN = function(x){subvar[row.names(knnres)[smoteidx[x, 1]],] + 
                                     zetas[x]*(subvar[knnres[smoteidx[x, 1], smoteidx[x, 2]],] - 
                                                 subvar[row.names(knnres)[smoteidx[x, 1]],])})
        names(synthesizeddat) <- row.names(knnres)[smoteidx[,1]]
        
        synthesizeddat <- do.call(rbind, synthesizeddat)
        
        
        datsub <- rbind(subvar, synthesizeddat)
        
        row.names(datsub) <- createnames(rawnames = row.names(datsub))
        
      }else if(subnum == 1){
        
        subvar <- dat
        
        knnres <- dbscan::kNN(x = subvar, k = min(k, nrow(subvar) - 1))
        
        mindist <- min(knnres$dist[knnres$dist > 0])
        mindistidx <- which(knnres$dist == mindist, arr.ind = TRUE)
        
        difference <- subvar[knnres$id[mindistidx[1, 1], mindistidx[1, 2]],] - 
          subvar[row.names(knnres$id)[mindistidx[1, 1]],]
        
        set.seed(i)
        zetas <- runif(targetnum - subnum)
        
        synthesizeddat <- lapply(X = 1:length(zetas), 
                                 FUN = function(x){subvar[sub$sampleid[1],] + zetas[x]*difference})
        
        names(synthesizeddat) <- rep(sub$sampleid[1], length(synthesizeddat))
        
        synthesizeddat <- do.call(rbind, synthesizeddat)
        
        datsub <- rbind(subvar[sub$sampleid[1], , drop = FALSE], synthesizeddat)
        
        row.names(datsub) <- createnames(rawnames = row.names(datsub))
        
      }
      
      if(j == 1){
        datsubs <- datsub
      }else{
        datsubs <- rbind(datsubs, datsub)
      }
      
    }
    
    datsubs <- t(datsubs)
    
    return(datsubs)
    
  }
  
  
  if(threads == 1){
    
    sampleddats <- list()
    
    for(i in 1:nround){
      
      tmpdat <- singleroundsampling(responsetypes = responsetypes,
                                    sampleidmapping = sampleidmapping,
                                    i = seed + i - 1,
                                    dat = dat, 
                                    k = k, 
                                    samplingmethod = samplingmethod)
      
      tmppd <- samplingpds(samplingsamples = colnames(tmpdat), 
                           pddat = pddat)
      
      sampleddats[[i]] <- tmpdat
      
      
    }
    
    
  }else{
    
    iseqs <- 1:nround
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    sampleddats <- foreach::foreach(i = iseqs,
                                    #.export = ls(name = globalenv())) %dopar% {
                                    .export = NULL) %dopar% {
                                      singleroundsampling(responsetypes = responsetypes,
                                                          sampleidmapping = sampleidmapping,
                                                          i = seed +i - 1,
                                                          dat = dat, 
                                                          k = k, 
                                                          samplingmethod = samplingmethod)
                                    }
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
  }
  
  sampledpds <- list()
  for(i in 1:nround){
    
    sampledpds[[i]] <- samplingpds(samplingsamples = colnames(sampleddats[[i]]), 
                                   pddat = pddat)
    
    
  }
  
  
  names(sampleddats) <- names(sampledpds) <- paste0('round_', 1:nround)
  
  res <- list(dats = sampleddats, 
              pds = sampledpds)
  
  return(res)
  
}




bootsampling <- function(dat,
                         pddat, 
                         responsename, 
                         nround = 10,
                         seed = 2022, 
                         threads = 1){
  
  dat <- t(dat)
  
  if(is.null(responsename) & ncol(pddat) >= 2){
    responsename <- names(pddat)[2]
  }
  
  if(is.null(responsename)){
    return(NULL)
  }
  
  
  
  sampleidmapping <- pddat[,c('sampleid', responsename)]
  names(sampleidmapping) <- c('sampleid', 'responsetype')
  
  dat <- dat[sampleidmapping$sampleid, , drop = FALSE]
  
  responsetypes <- unique(sampleidmapping$responsetype)
  
  
  
  samplingpds <- function(samplingsamples, 
                          pddat){
    
    samples <- gsub(pattern = '\\.[0-9]*$', replacement = '', x = samplingsamples)
    samplingpds <- pddat[match(samples, pddat$sampleid), , drop = FALSE]
    samplingpds$sampleid <- samplingsamples
    row.names(samplingpds) <- 1:nrow(samplingpds)
    
    return(samplingpds)
    
  }
  
  createnames <- function(rawnames){
    
    reversenames <- FALSE
    if(sum(grepl(pattern = '-', x = rawnames)) > 0){
      
      rawnames <- gsub(pattern = '-', replacement = '___', x = rawnames)
      reversenames <- TRUE
    }
    
    
    newnames <- make.names(names = rawnames, unique = TRUE)
    suffix <- gsub(pattern = '^.*\\.', replacement = '', x = newnames)
    suffix[suffix %in% rawnames] <- ''
    suffix <- paste0('.', suffix)
    suffix[suffix == '.'] <- ''
    newnames <- paste0(rawnames, suffix)
    
    if(reversenames == TRUE){
      newnames <- gsub(pattern = '___', replacement = '-', x = newnames)
    }
    
    return(newnames)
    
  }
  
  singleroundsampling <- function(responsetypes,
                                  sampleidmapping,
                                  i = 2022,
                                  dat){
    
    
    j <- 1
    for(j in 1:length(responsetypes)){
      
      ResponseType <- responsetypes[j]
      sub <- subset(sampleidmapping, responsetype == ResponseType)
      
      subnum <- nrow(sub)
      
      targetnum <- subnum
      
      
        
      sampleids <- sub$sampleid
      
      set.seed(i)
      samples <- sample(x = sampleids,
                        size = targetnum,
                        replace = TRUE)
      
      datsub <- dat[samples, , drop = FALSE]
      
      row.names(datsub) <- createnames(rawnames = row.names(datsub))
      
      
      if(j == 1){
        datsubs <- datsub
      }else{
        datsubs <- rbind(datsubs, datsub)
      }
      
    }
    
    datsubs <- t(datsubs)
    
    return(datsubs)
    
  }
  
  
  if(threads == 1){
    
    sampleddats <- list()
    
    for(i in 1:nround){
      
      tmpdat <- singleroundsampling(responsetypes = responsetypes,
                                    sampleidmapping = sampleidmapping,
                                    i = seed + i - 1,
                                    dat = dat)
      
      tmppd <- samplingpds(samplingsamples = colnames(tmpdat), 
                           pddat = pddat)
      
      sampleddats[[i]] <- tmpdat
      
      
    }
    
    
  }else{
    
    iseqs <- 1:nround
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    sampleddats <- foreach::foreach(i = iseqs,
                                    #.export = ls(name = globalenv())) %dopar% {
                                    .export = NULL) %dopar% {
                                      singleroundsampling(responsetypes = responsetypes,
                                                          sampleidmapping = sampleidmapping,
                                                          i = seed +i - 1,
                                                          dat = dat)
                                    }
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
  }
  
  sampledpds <- list()
  for(i in 1:nround){
    
    sampledpds[[i]] <- samplingpds(samplingsamples = colnames(sampleddats[[i]]), 
                                   pddat = pddat)
    
    
  }
  
  
  names(sampleddats) <- names(sampledpds) <- paste0('round_', 1:nround)
  
  res <- list(dats = sampleddats, 
              pds = sampledpds)
  
  return(res)
  
}










combinepval <- function(limmareslist, 
                        method = 'Cauchy'){
  
  probes <- row.names(limmareslist[[1]])
  pvallist <- list()
  logfclist <- list()
  i <- 1
  for(i in 1:length(limmareslist)){
    
    limmareslist[[i]] <- limmareslist[[i]][probes,]
    
    pvals <- limmareslist[[i]]$P.Value
    logfcs <- limmareslist[[i]]$logFC
    
    names(pvals) <- probes
    names(logfcs) <- probes
    
    pvallist[[i]] <- pvals
    logfclist[[i]] <- logfcs
    
  }
  
  names(pvallist) <- names(logfclist) <- names(limmareslist)
  pvaldat <- do.call(rbind, pvallist)
  logfcdat <- do.call(rbind, logfclist)
  
  if(method == 'Fisher'){
    
    pcombs <- apply(X = pvaldat, MARGIN = 2, FUN = function(x){metap::sumlog(x)$p})
    padjs <- p.adjust(pcombs, method = 'BH')
    
    
  }else if(method == 'harmonicmean'){
    
    pcombs <- apply(X = pvaldat, MARGIN = 2, 
                    FUN = harmonicmeanp::p.hmp, 
                    w = rep(1/nrow(pvaldat), nrow(pvaldat)), 
                    L = nrow(pvaldat))
    padjs <- p.adjust(pcombs, method = 'BH')
    
  }else if(method == 'Cauchy'){
    
    #pcombs <- apply(X = pvaldat, MARGIN = 2, FUN = ACAT::ACAT)
    pcombs <- apply(X = pvaldat, MARGIN = 2, FUN = ACAT)
    padjs <- p.adjust(pcombs, method = 'BH')
    
  }
  
  logfcsems <- apply(X = logfcdat, MARGIN = 1, FUN = function(x){sqrt(var(x)/length(x))})
  ivws <- 1/logfcsems/sum(1/logfcsems)
  logfccombs <- (t(logfcdat) %*% ivws)[,1]
  
  #suppressWarnings(
    #if(is.na(logfccombs)){
    #  logfccombs <- colMeans(logfcdat)
    #}
  #)
  
  res <- list(pcombs = pcombs, 
              padjs = padjs, 
              logfccombs = logfccombs)
  
  return(res)
  
}

#Get volcano plot for limma comparison result
singlevolcano <- function(limmares, 
                          betaval = TRUE, 
                          absxcutoff = 0, 
                          pvalcut = 0.05, 
                          pvalcolname = 'adj.P.Val', 
                          confoundings = NULL, 
                          responsevarname, 
                          titleprefix = NULL, 
                          featuretype = 'probe', 
                          titlesize = 15, 
                          textsize = 13, 
                          face = 'bold', 
                          labelnum = 5, 
                          annotextsize = 4, 
                          
                          responsecause = NULL, 
                          responseres = NULL){
  
  if(absxcutoff < 0){
    absxcutoff <- -absxcutoff
  }
  
  if(betaval == TRUE){
    xaxis <- 'Difference'
    
    if(abs(absxcutoff) > 1){
      absxcutoff <- 0
    }
    
  }else{
    xaxis <- 'log2FC'
  }
  
  vars <- paste(c(responsevarname, confoundings), collapse = ' + ')
  
  
  if('remove' %in% names(limmares)){
    
    sigg.limma <- limmares[(limmares[,pvalcolname] < pvalcut) & 
                             (abs(limmares$logFC) > absxcutoff) & 
                             (limmares$remove == FALSE),]
    nonsigg.limma <- limmares[(limmares[,pvalcolname] >= pvalcut) | 
                                (abs(limmares$logFC) <= absxcutoff) | 
                                (limmares$remove == TRUE),]
    
  }else{
    
    sigg.limma <- limmares[(limmares[,pvalcolname] < pvalcut) & 
                             (abs(limmares$logFC) > absxcutoff),]
    nonsigg.limma <- limmares[(limmares[,pvalcolname] >= pvalcut) | 
                                (abs(limmares$logFC) <= absxcutoff),]
    
  }
  
  if(('color' %in% names(limmares)) & 
     ('Type' %in% names(limmares))){
    
    both_pos <- limmares
    names(both_pos)[names(both_pos) == 'color'] <- 'denscolor'
    both_pos$probe <- row.names(both_pos)
    types <- unique(both_pos$Type)
    subtitlename <- paste0(types, ' = ', as.vector(table(limmares$Type)[types]))
    subtitlename <- rev(subtitlename)
    subtitlename <- paste(subtitlename, collapse = ', ')
    subtitlename <- paste0('(', subtitlename, ')')
    
    both_pos$Sample_Group <- subtitlename
    row.names(both_pos) <- 1:nrow(both_pos)
    
    legendmapping <- both_pos[c('Type', 'denscolor')]
    legendmapping <- unique(legendmapping)
    legendmapping <- legendmapping[match(rev(types), legendmapping$Type),]
    row.names(legendmapping) <- 1:nrow(legendmapping)
    
    
  }else{
    
    if(nrow(sigg.limma) == 0){
      
      both_pos <- limmares[c('logFC', pvalcolname)]
      both_pos$Type <- 'NonSig'
      both_pos$denscolor <- '#C0C0C0'
      both_pos$red <- 192
      both_pos$green <- 192
      both_pos$blue <- 192
      both_pos$probe <- row.names(both_pos)
      
      subtitlename <- paste0('(NC = ', nrow(both_pos), ')')
      
      both_pos$Sample_Group <- subtitlename
      row.names(both_pos) <- 1:nrow(both_pos)
      
    }else{
      
      #Draw probe valcano
      both_pos <- sigg.limma[c('logFC', pvalcolname)]
      
      both_pos <- as.data.frame(both_pos, stringsAsFactors = FALSE)
      
      both_pos$Type <- 'NonSig'
      both_pos$Type[both_pos$logFC > absxcutoff] <- 'UP'
      both_pos$Type[both_pos$logFC < -absxcutoff] <- 'DN'
      
      hypernum <- nrow(subset(both_pos, Type == 'UP'))
      hyponum <- nrow(subset(both_pos, Type == 'DN'))
      
      both_pos_nonsig <- subset(both_pos, Type == 'NonSig')
      both_pos <- subset(both_pos, Type != 'NonSig')
      
      #library(ggplot2)
      #library(RColorBrewer)
      
      if(nrow(both_pos) > 50){
        myColor <- densCols(both_pos$logFC, -log10(both_pos[,pvalcolname]), 
                            colramp = colorRampPalette(rev(rainbow(10, end=4/6))))
      }else{
        myColor <- rep('blue', nrow(both_pos))
      }
      
      
      both_pos$denscolor <- myColor
      rgbmat <- t(col2rgb(myColor))
      rgbmat <- as.data.frame(rgbmat, stringsAsFactors = FALSE)
      both_pos <- cbind(both_pos, rgbmat)
      both_pos <- both_pos[order(-both_pos$blue, both_pos$red, both_pos$green),]
      both_pos1 <- subset(both_pos, blue >= red)
      both_pos2 <- subset(both_pos, blue < red)
      both_pos2 <- both_pos2[order(-both_pos2$blue, both_pos2$red, -both_pos2$green),]
      both_pos <- rbind(both_pos1, both_pos2)
      
      both_nonsig <- nonsigg.limma[c('logFC', pvalcolname)]
      
      if(nrow(both_nonsig) > 0){
        both_nonsig$Type <- 'NonSig'
        both_nonsig <- rbind(both_nonsig, both_pos_nonsig)
        both_nonsig$denscolor <- '#C0C0C0'
        both_nonsig$red <- 192
        both_nonsig$green <- 192
        both_nonsig$blue <- 192
        
        both_pos <- rbind(both_nonsig, both_pos)
        
      }
      
      nonsignum <- nrow(both_nonsig)
      
      
      
      both_pos$probe <- row.names(both_pos)
      
      ncnum <- nrow(both_pos) - hypernum - hyponum
      
      
      pieces <- c(paste0('UP = ', hypernum), paste0('DN = ', hyponum), 
                  paste0('NC = ', ncnum))
      
      pieces <- pieces[c(hypernum, hyponum, ncnum) != 0]
      subtitlename <- paste(pieces, collapse = ', ')
      subtitlename <- paste0('(', subtitlename, ')')
      
      
      
      both_pos$Sample_Group <- subtitlename
      row.names(both_pos) <- 1:nrow(both_pos)
      
    }
    
    
  }
  
  
  
  names(both_pos)[names(both_pos) == pvalcolname] <- 'pval'
  
  if(pvalcolname == 'P.Value'){
    pvalaxisname <- 'p-value'
  }else{
    pvalaxisname <- 'adjusted p-value'
  }
  
  if(is.null(titleprefix)){
    maintitle <- paste0('Differential ', featuretype, 's (', vars, ')')
  }else{
    maintitle <- paste0(titleprefix, ' differential ', featuretype, 's (', vars, ')')
  }
  
  if(nchar(maintitle) > 80){
    
    if(is.null(titleprefix)){
      maintitle <- paste0('Differential ', featuretype, 
                          's(', responsevarname, ' + confoundings)')
    }else{
      maintitle <- paste0(titleprefix, ' differential ', featuretype, 
                          's(', responsevarname, ' + confoundings)')
    }
    
  }
  
  
  if(!('label' %in% names(both_pos))){
    
    both_pos$label <- FALSE
    
    if(!is.null(labelnum)){
      
      labelnum <- as.numeric(labelnum)
      
      if(is.na(labelnum)){
        labelnum <- NULL
      }
    }
    
    if(!is.null(labelnum)){
      
      subup <- subset(both_pos, Type == 'UP')
      subdn <- subset(both_pos, Type == 'DN')
      
      if(nrow(subup) > 0){
        labelupnum <- min(labelnum, nrow(subup))
        
        labelupprobes <- subup[order(subup$pval, -abs(subup$logFC)),]$probe[seq(1, labelupnum, 1)]
        
      }else{
        labelupnum <- 0
        
        labelupprobes <- NULL
      }
      
      if(nrow(subdn) > 0){
        labeldnnum <- min(labelnum, nrow(subdn))
        
        labeldnprobes <- subdn[order(subdn$pval, -abs(subdn$logFC)),]$probe[seq(1, labeldnnum, 1)]
        
      }else{
        labeldnnum <- 0
        
        labeldnprobes <- NULL
      }
      
      both_pos$label[both_pos$probe %in% c(labelupprobes, labeldnprobes)] <- TRUE
      
      
    }
    
  }
  
  
  if(('color' %in% names(limmares)) & 
     ('Type' %in% names(limmares))){
    
    p <- ggplot2::ggplot(both_pos, ggplot2::aes(x = logFC, y = -log10(pval), 
                                                color = Type))
    
    print(
      
      p + ggplot2::geom_point(position='jitter') + 
        ggplot2::xlab(xaxis) + ggplot2::ylab(paste0('-log10(', pvalaxisname, ')')) + 
        ggplot2::ggtitle(maintitle, 
                         subtitle = subtitlename) + 
        ggplot2::geom_hline(yintercept = -log10(pvalcut), color = 'red') + 
        ggplot2::geom_vline(xintercept = absxcutoff, color = 'blue') + 
        ggplot2::geom_vline(xintercept = -absxcutoff, color = 'blue') + 
        ggplot2::geom_hline(yintercept = 0, color = 'gray', size = 0.5) + 
        ggplot2::geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
        ggplot2::theme_bw() + 
        ggplot2::scale_color_manual(name = 'Type',
                                    breaks = legendmapping$Type,
                                    values = legendmapping$denscolor, 
                                    guide = 'legend') + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                       plot.title = ggplot2::element_text(size = titlesize, face = face),
                       plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                       axis.text = ggplot2::element_text(size = textsize, face = face), 
                       axis.title = ggplot2::element_text(size = textsize, face = face), 
                       plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"), 
                       legend.title = ggplot2::element_text(size = textsize, face = face), 
                       legend.text = ggplot2::element_text(size = textsize, face = face)) + 
        
        ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(label == TRUE, as.character(probe),'')), 
                                 size = annotextsize, fontface = face, max.overlaps = 1000, 
                                 color = 'black')
      
      
    )
    
  }else{
    
    p <- ggplot2::ggplot(both_pos, ggplot2::aes(x = logFC, y = -log10(pval)))
    
    print(
      
      p + ggplot2::geom_point(color = both_pos$denscolor, position='jitter') + 
        ggplot2::xlab(xaxis) + ggplot2::ylab(paste0('-log10(', pvalaxisname, ')')) + 
        ggplot2::ggtitle(maintitle, 
                         subtitle = subtitlename) + 
        ggplot2::geom_hline(yintercept = -log10(pvalcut), color = 'red') + 
        ggplot2::geom_vline(xintercept = absxcutoff, color = 'blue') + 
        ggplot2::geom_vline(xintercept = -absxcutoff, color = 'blue') + 
        ggplot2::geom_hline(yintercept = 0, color = 'gray', size = 0.5) + 
        ggplot2::geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
        ggplot2::theme_bw() + 
        ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                       plot.title = ggplot2::element_text(size = titlesize, face = face),
                       plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                       axis.text = ggplot2::element_text(size = textsize, face = face), 
                       axis.title = ggplot2::element_text(size = textsize, face = face), 
                       plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")) + 
        
        ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(label == TRUE, as.character(probe),'')), 
                                 size = annotextsize, fontface = face, max.overlaps = 1000)
      
      
    )
    
  }
  
  row.names(both_pos) <- both_pos$probe
  both_pos <- both_pos[,c('pval', 'logFC', 'Type', 'denscolor', 'label')]
  names(both_pos)[1] <- pvalcolname
  names(both_pos)[4] <- 'color'
  names(both_pos)[5] <- 'labeled'
  
  return(both_pos)
  
}

#'Perform ensemble-based limma or normal limma 
#'
#'Perform ensemble-based limma or normal limma to find features relevant to 
#'  specific phenotypic variable. 
#'  
#'@param dat A matrix with the feature values for samples. Each column is one 
#'  sample and each row represents one feature. The row names are the feature 
#'  names and the column names should be sample IDs. 
#'@param pddat Meta data frame. The first column should be the sample names 
#'  with a column name as "sampleid", and the remaining should be phenotypic 
#'  variables.
#'@param responsevarname The column name of the response variable in the data 
#'  frame provided to \code{pddat}. It is the variable of interest. Default 
#'  is NULL, which means the second column of \code{pddat} is the response. 
#'  This response variable should be a binary or continuous variable.
#'@param responselevels If the response variable is binary, it will be treated 
#'  as a factor, and this parameter is needed to define the factor element 
#'  level, so it should be a vector defining this level. Only needed if the 
#'  response in \code{pddat} is a character vector but not a factor. In this 
#'  case, it will be converted to a factor by the function, and this paramter 
#'  \code{responselevels} will be used to define the factor element level. Its 
#'  default is NULL, meaning that this element level will be set automatically 
#'  following the character order of the elements in the response variable.
#'@param confoundings The column name of the confounding factors in the data 
#'  frame provided to \code{pddat}. Default is NULL, which means there is no 
#'  confounding factor in the data. 
#'@param balanceadj Whether to perform ensemble-based limma or not. When it 
#'  is set as TRUE. The function will perform ensemble-based limma. It means 
#'  a sampling process will be used on the samples to generate several base 
#'  learner datasets, and the sample groups of each base learner set will be 
#'  adjusted during the sampling so that they will have the same size. Then, 
#'  limma will be used on each base learner to call the differential features, 
#'  and their results will finally be ensembled together to get the result for 
#'  the whole data. If it were set as FALSE, normal limma will be performed. 
#'  Default is FALSE.
#'@param samplingmethod When \code{balanceadj} is TRUE, this parameter will 
#'  be needed to determine how to sample the data to get the base learner sets 
#'  for the ensemeble mode. If it is set as "updn", to make the sample groups 
#'  have the same sample size in a base leaner dataset, a large group will be 
#'  down-sampled randomly, and a small group will be up-sampled with SMOTE 
#'  (synthetic minority over-sampling technique), so that the final group size 
#'  will be the total sample number/group number. If it is set as "up", the 
#'  samples in a small group will be up-sampled, so its sample size will be 
#'  the large group size. If this parameter is set as "dn", the samples in the 
#'  large group will be down-sampled so that its sample size will be the small 
#'  sample size. Default is "updn".
#'@param nround When \code{balanceadj} is TRUE, this parameter will be used 
#'  to determine the number of the base learner sets in the ensemble-based 
#'  mode. Default is 10.
#'@param seed Random seed for the sampling process to generate base learner 
#'  sets in ensemble-based limma.
#'@param k When \code{balanceadj} is TRUE, and \code{samplingmethod} is "updn" 
#'  or "up", this parameter will be used for up-sampling. Because up-sampling 
#'  is based on SMOTE, which needs to synthesize new samples from the closet 
#'  neighbors of a randomly selected sample in the same group, this parameter 
#'  is used to define how many closet neighbors will be needed. Default is 5, 
#'  meaning the top 5 closet neighbors will be used for the synthesis.
#'@param threads Number of threads need for parallelization. Default is 1.
#'@param method When \code{balanceadj} is TRUE, this parameter will be needed 
#'  for aggregating the results from all the base learners in ensemble-based 
#'  limma. Each base learner will return a set of p-vals for all the features 
#'  in the dataset, so for each feature, it will finally get several p-vals 
#'  from different base learners, and to combine them into one value, some 
#'  p-val combination methods can be used. This parameter is used to choose 
#'  the method, can be "Cauchy", "harmonicmean", or "Fisher". Default value is 
#'  "Cauchy", which uses ACAT (aggregated Cauchy association test) for the 
#'  combination. ACAT is suitable to combine strongly dependent p-vals, which 
#'  is the case of ensemble-based limma because its base learner sets are from 
#'  the same original dataset.
#'@param removereg A GRanges object or a data frame that records the genomic 
#'  coordinates of any regions that need to be excluded from the analysis, can 
#'  be the genomic regions related to any confounding factors. If a feature 
#'  called by normal limma or ensemble-based limma is located there, it will be 
#'  removed. For regions related to confounding factors, this region removal 
#'  can help avoid influence from the confounding factors. However, because 
#'  both normal limma and ensemble-based limma can adjust confoundings during 
#'  their regression step, a further region removal process will have very 
#'  limited effects. Hence, this step is unnecessary, and can be skipped by 
#'  setting this parameter \code{removereg} as NULL. If need to transfer a 
#'  data frame, it should have columns named as "seqnames", "start", "end", 
#'  and "strand", which help to define the genomic coordinates of the regions. 
#'  Each row of the data frame should be a region. The default value is NULL.
#'@param featuretype The type of the features in the data. It will be used 
#'  when \code{removereg} is not NULL and judge whether any features located 
#'  in the genomic regions to be removed. Can be "gene" or "probe" (for probes 
#'  in the DNA methylation data). If another parameter \code{plot} is TRUE, 
#'  it will also appear in the title of the volcano plot generated by this 
#'  function. Default value is "probe".
#'@param platform If the data to be analyzed is methylation probe data, this 
#'  parameter will be used to indicate the platform of the data. Can be set as 
#'  450 (for 450k platform), 850 (for EPIC), or 27 (27K platform). Its default 
#'  value is 450.
#'@param ignorestrand When the parameter \code{removereg} is not NULL and some 
#'  genomic regions need to be excluded from the analysis, this parameter will 
#'  be used to judge whether any features located in the genomic regions. If 
#'  it is FALSE, the strand information of the features should match that of 
#'  the regions if they are judged as within the regions, and if it is TRUE, 
#'  the strand information will not be considered. Default is TRUE. 
#'@param  tssradius When \code{removereg} is not NULL and \code{featuretype} 
#'  is "gene". This parameter will be used to judge whether a gene located in 
#'  any of the genomic regions. If the TSS of a gene is located within the 
#'  regions, or not within, but the distance between its TSS and the end of 
#'  any regions is less than the value of this \code{tssradius}, this gene 
#'  will be judged as located within the regions (if the strand information 
#'  required by \code{ignorestrand} is also fulfilled). Then, the gene will 
#'  be excluded from the analysis. Default is 1500 (for 1500 bp).
#'@param plot Whether need to plot a volcano plot to show the differential 
#'  features. Default is FALSE. 
#'@param isbetaval If the feature values in the data are log2 transformed, 
#'  set it as FALSE, if they are not log2 transformed, i.e., they are gene 
#'  counts or methylation beta values, set this parameter as TRUE. Will only 
#'  be needed when \code{plot} is TRUE.
#'@param absxcutoff The cutoff on log2FC (for log2 transformed features) or 
#'  inter-group difference (for non-log2 transformed feature counts or beta 
#'  values) to judge whether a feature is significantly different between the 
#'  sample groups. Default is 0, meaning any log2FC or difference value will 
#'  not influence the judgment on feature significance.
#'@param pvalcutoff The p-val (or adjusted p-val) cutoff to judge whether a 
#'  feature is significantly different between the groups or not. Default is 
#'  0.05.
#'@param pvalcolname Which p-val should be used to judge the significance of 
#'  the features. Can be "P.Val" or "adj.P.Val". Default is "adj.P.Val".
#'@param titleprefix The prefix of the volcano plot title. It can be set as 
#'  any character string need to be shown in the title. Default is NULL.
#'@param titlesize The font size for the plot title of \code{plot}. Default 
#'  is 15.
#'@param textsize The font size for the plot texts of \code{plot}. Default is 
#'  13.
#'@param face The font face for the plot texts of \code{plot}. Default value 
#'  is "bold".
#'@param labelnum In the volcano plot, the names of the top up-regulated and 
#'  top down-regulated features can be labeled. This number indicates how many 
#'  top features in the groups should be labeled. If it is NULL, no feature 
#'  name will be labeled. Can also be any number. Default is NULL.
#'@param annotextsize If \code{labelnum} is not NULL, this parameter will be 
#'  needed to set the font size of the gene names to be labeled in the volcano 
#'  plot. Default is 4.
#'@return If \code{plot} were FALSE, will only return a data frame recording 
#'  the normal limma or ensemble-based limma results for all the features in 
#'  the data, including the adjusted p-val, log2FC, etc. If \code{plot} were 
#'  TRUE, will return a list containing not only the limma results, but also 
#'  the plot data, and a volcano plot will be generated.
#'@export
difffeatures <- function(dat, 
                         pddat, 
                         responsevarname = NULL, 
                         responselevels = NULL, 
                         confoundings = NULL, 
                         balanceadj = FALSE, 
                         samplingmethod = 'updn', 
                         nround = 10,
                         seed = 2022, 
                         k = 5, 
                         threads = 1, 
                         method = 'Cauchy', 
                         removereg = NULL, 
                         featuretype = 'probe', 
                         platform = 450, 
                         ignorestrand = TRUE, 
                         tssradius = 1500, 
                         plot = FALSE, 
                         isbetaval = TRUE, 
                         absxcutoff = 0, 
                         pvalcutoff = 0.05, 
                         pvalcolname = 'adj.P.Val', 
                         titleprefix = NULL, 
                         titlesize = 15, 
                         textsize = 13, 
                         face = 'bold', 
                         labelnum = NULL, 
                         annotextsize = 4){
  
  if(is.null(responsevarname) & ncol(pddat) >= 2){
    responsevarname <- names(pddat)[2]
  }
  
  if(is.null(responsevarname)){
    return(NULL)
  }
  
  
  
  if(!is.factor(pddat[[responsevarname]]) & !is.numeric(pddat[[responsevarname]])){
    
    if(is.null(responselevels)){
      responselevels <- pddat[[responsevarname]]
      responselevels <- unique(responselevels)
      responselevels <- responselevels[order(responselevels)]
    }
    
    Response <- factor(pddat[[responsevarname]], 
                       levels = responselevels, 
                       ordered = TRUE)
  }else{
    
    Response <- pddat[[responsevarname]]
    
  }
  
  Responsetype <- judgevartype(varvals = Response)
  
  if((Responsetype != 'binary') & (Responsetype != 'continuous')){
    
    cat('This function only supports binary or continuous variables but the current response is not.\n')
    
    return(NULL)
    
  }
  
  if((Responsetype != 'binary') & (balanceadj == TRUE)){
    
    cat('The ensemble mode only supports binary variable, but the current response is not.\nHence, the mode has been changed to normal mode.\n')
    
    balanceadj <- FALSE
    
  }
  
  
  if(balanceadj == TRUE){
    
    balanceddats <- balancesampling(dat = dat, 
                                    pddat = pddat, 
                                    responsename = responsevarname, 
                                    nround = nround, 
                                    seed = seed, 
                                    k = k, 
                                    threads = threads, 
                                    samplingmethod = samplingmethod)
    
    if(threads == 1){
      
      limmareslist <- list()
      
      iseqs <- seq(1, nround, 1)
      
      for(i in iseqs){
        
        limmares <- diffsites(dats = balanceddats$dats, 
                              pddats = balanceddats$pds, 
                              i = i, 
                              responsevarname = responsevarname, 
                              responselevels = responselevels, 
                              confoundings = confoundings)
        
        limmareslist[[i]] <- limmares
        
        
      }
      
    }else{
      
      diffsites <- function(dats, 
                            pddats, 
                            i = 1, 
                            responsevarname = NULL, 
                            responselevels = NULL, 
                            confoundings = NULL){
        
        
        dat <- dats[[i]]
        pddat <- pddats[[i]]
        
        if(is.null(responsevarname) & ncol(pddat) >= 2){
          responsevarname <- names(pddat)[2]
        }
        
        if(is.null(responsevarname)){
          return(NULL)
        }
        
        #library(limma)
        
        simpddat <- pddat[c('sampleid', responsevarname, confoundings)]
        simpddat <- simpddat[complete.cases(simpddat),]
        
        if(!is.factor(pddat[[responsevarname]]) & !is.numeric(pddat[[responsevarname]])){
          
          if(is.null(responselevels)){
            responselevels <- simpddat[[responsevarname]]
            responselevels <- unique(responselevels)
            responselevels <- responselevels[order(responselevels)]
          }
          simpddat$Response <- factor(simpddat[[responsevarname]], 
                                      levels = responselevels, 
                                      ordered = TRUE)
        }else{
          
          simpddat$Response <- simpddat[[responsevarname]]
          
        }
        
        
        vars <- paste(c('Response', confoundings), collapse = ' + ')
        
        design <- model.matrix(as.formula(paste0('~ ', vars)), 
                               data = simpddat)
        
        betasub <- dat[,simpddat$sampleid, drop = FALSE]
        
        
        
        fit1 <- limma::lmFit(betasub, design)
        fit1 <- limma::eBayes(fit1)
        
        allg.limma <- limma::topTable(fit1, coef = 2, n = dim(fit1)[1])
        
        if((nrow(allg.limma) == 1) & (nrow(betasub) == 1)){
          row.names(allg.limma) <- row.names(betasub)
        }
        
        return(allg.limma)
        
      }
      
      
      iseqs <- seq(1, nround, 1)
      
      #library(doParallel)
      
      cores <- parallel::detectCores()
      cl <- parallel::makeCluster(min(threads, cores))
      
      doParallel::registerDoParallel(cl)
      
      #date()
      `%dopar%` <- foreach::`%dopar%`
      limmareslist <- foreach::foreach(i = iseqs,
                                       #.export = ls(name = globalenv())) %dopar% {
                                       .export = NULL) %dopar% {
                                         diffsites(dats = balanceddats$dats, 
                                                   pddats = balanceddats$pds, 
                                                   i, 
                                                   responsevarname = responsevarname, 
                                                   responselevels = responselevels, 
                                                   confoundings = confoundings)
                                       }
      
      parallel::stopCluster(cl)
      
      unregister_dopar()
      
      
    }
    
    names(limmareslist) <- names(balanceddats$dats)
    
    combineres <- combinepval(limmareslist = limmareslist, method = method)
    
    metalimmares <- data.frame(logFC = combineres$logfccombs, 
                               P.Value = combineres$pcombs, 
                               adj.P.Val = combineres$padjs, 
                               stringsAsFactors = FALSE)
    
    limmares <- metalimmares
    
    
  }else{
    
    balanceddats <- list(dats = list(round_1 = dat), 
                         pds = list(round_1 = pddat))
    
    limmares <- diffsites(dats = balanceddats$dats, 
                          pddats = balanceddats$pds, 
                          i = 1, 
                          responsevarname = responsevarname, 
                          responselevels = responselevels, 
                          confoundings = confoundings)
    
  }
  
  
  if(!is.null(removereg)){
    
    if(featuretype == 'probe'){
      removedfeatures <- overlappedprobes(totalprobes = row.names(limmares), 
                                          removereg = removereg, 
                                          platform = platform, 
                                          ignorestrand = ignorestrand)
    }else if(featuretype == 'gene'){
      removedfeatures <- overlappedgenes(totalgenes = row.names(limmares), 
                                         removereg = removereg, 
                                         tssradius = tssradius, 
                                         ignorestrand = ignorestrand)
    }else{
      removedfeatures <- NULL
    }
    
    
    if(!is.null(removedfeatures)){
      
      if(featuretype == 'probe'){
        removedfeatures <- unique(row.names(removedfeatures))
        
        limmares$remove <- FALSE
        limmares$remove[row.names(limmares) %in% removedfeatures] <- TRUE
        
      }else if(featuretype == 'gene'){
        
        removegenes <- idxgenenames(totalgenes = row.names(limmares), 
                                    genesym = removedfeatures$genesym, 
                                    geneid = removedfeatures$geneid)
        limmares$remove <- removegenes
        
      }else{
        limmares$remove <- FALSE
      }
      
    }else{
      limmares$remove <- FALSE
    }
    
  }
  
  
  
  if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val', 
                                 'pvalue', 'p-value', 'p.value')){
    pvalcolname <- 'P.Value'
  }else{
    pvalcolname <- 'adj.P.Val'
  }
  
  #otherpvalcolname <- setdiff(c('P.Value', 'adj.P.Val'), pvalcolname)
  
  limmares <- limmares[order(limmares[,pvalcolname], 
                             #limmares[otherpvalcolname], 
                             -abs(limmares$logFC)), , drop = FALSE]
  
  
  if(plot == TRUE){
    
    if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val', 
                                   'pvalue', 'p-value', 'p.value')){
      pvalcolname <- 'P.Value'
    }else{
      pvalcolname <- 'adj.P.Val'
    } 
    
    plotdat <- singlevolcano(limmares = limmares, 
                             betaval = isbetaval, 
                             absxcutoff = absxcutoff, 
                             pvalcut = pvalcutoff, 
                             pvalcolname = pvalcolname, 
                             confoundings = confoundings, 
                             responsevarname = responsevarname, 
                             titleprefix = titleprefix, 
                             featuretype = featuretype, 
                             titlesize = titlesize, 
                             textsize = textsize, 
                             face = face, 
                             labelnum = labelnum, 
                             annotextsize = annotextsize)
    
  }
  
  if(plot == FALSE){
    
    res <- limmares
    
  }else{
    
    res <- list(limmares = limmares, plotdat = plotdat)
    
  }
  
  return(res)
  
}



