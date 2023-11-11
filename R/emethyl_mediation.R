#Mediation analysis#####

judgevartype <- function(varvals){
  
  if(length(unique(varvals)) < 2){
    stop("The source variable must have at least two unique values.")
  }else if(length(unique(varvals)) == 2){
    vartype <- "binary"
  }else if(is.factor(varvals) | is.character(varvals)){
    vartype <- "multinomial"
  }else{
    vartype <- "continuous"
  }
  
  return(vartype)
  
}

singlemediation <- function(j, 
                            mediatordat, 
                            predictdat, 
                            target, 
                            source, 
                            confoundings, 
                            interaction = FALSE, 
                            IPW, 
                            responselevels, 
                            mod11form, 
                            mod1form, 
                            mod2form){
  
  mediator <- colnames(mediatordat)[j]
  mediatorval <- mediatordat[, j, drop = FALSE]
  names(mediatorval) <- 'mediator'
  trainingdat <- cbind(predictdat, mediatorval)
  trainingdat <- trainingdat[complete.cases(trainingdat),]
  
  if(IPW == TRUE & length(confoundings) > 0){
    
    trainingdat <- trainingdat[c(target, source, 'mediator', confoundings, 'invwt')]
    
  }else{
    
    trainingdat <- trainingdat[c(target, source, 'mediator', confoundings)]
    
  }
  
  if(!is.factor(trainingdat[[target]]) & !is.numeric(trainingdat[[target]])){
    
    if(is.null(responselevels)){
      responselevels <- trainingdat[[target]]
      responselevels <- unique(responselevels)
      responselevels <- responselevels[order(responselevels)]
    }
    trainingdat[[target]] <- factor(trainingdat[[target]], 
                                    levels = responselevels, 
                                    ordered = TRUE)
    #trainingdat[[target]] <- as.numeric(trainingdat[[target]]) - 1
  }
  
  targettype <- judgevartype(varvals = trainingdat[[target]])
  
  if(targettype == 'binary'){
    
    #no interaction mod
    mod11 <- glm(mod11form, family = binomial(link = "logit"), data = trainingdat)
    #interaction mod
    mod1 <- glm(mod1form, family = binomial(link = "logit"), data = trainingdat)
    
  }else if(targettype == 'continuous'){
    
    #no interaction mod
    mod11 <- lm(mod11form, data = trainingdat)
    #interaction mod
    mod1 <- lm(mod1form, data = trainingdat)
    
  }
  
  if(sum(is.na(mod11$coefficients)) > 0){
    return(NULL)
  }
  
  if(sum(is.na(mod1$coefficients)) > 0){
    return(NULL)
  }
  
  #mediator prediction mod
  if(IPW == TRUE & length(confoundings) > 0){
    mod2 <- lm(mod2form, data = trainingdat, weights = invwt)
  }else{
    mod2 <- lm(mod2form, data = trainingdat)
  }
  
  if(sum(is.na(mod2$coefficients)) > 0){
    return(NULL)
  }
  
  #coef for source (interaction mod)
  theta1 <- mod1$coefficients[2]
  sd.theta1 <- summary(mod1)$coefficients[2, 2]
  pval.theta1 <- summary(mod1)$coefficients[2, 4]
  
  #coef for mediator (interaction mod)
  theta2 <- mod1$coefficients[3]
  sd.theta2 <- summary(mod1)$coefficients[3, 2]
  pval.theta2 <- summary(mod1)$coefficients[3, 4]
  
  #coef for interaction (interaction mod)
  interidx <- length(mod1$coefficients)
  theta3 <- mod1$coefficients[interidx]
  sd.theta3 <- summary(mod1)$coefficients[interidx, 2]
  pval.theta3 <- summary(mod1)$coefficients[interidx, 4]
  
  #coef for source (no interactio mod)
  theta11 <- mod11$coefficients[2]
  sd.theta11 <- summary(mod11)$coefficients[2, 2]
  pval.theta11 <- summary(mod11)$coefficients[2, 4]
  
  #coef for mediator (no interaction mod)
  theta21 <- mod11$coefficients[3]
  sd.theta21 <- summary(mod11)$coefficients[3, 2]
  pval.theta21 <- summary(mod11)$coefficients[3, 4]
  
  #intercept (mediator prediction mod)
  beta0 <- mod2$coefficients[1]
  
  #coef for source (mediator prediction mod)
  beta1 <- mod2$coefficients[2]
  sd.beta1 <- summary(mod2)$coefficients[2, 2]
  pval.beta1 <- summary(mod2)$coefficients[2, 4]
  
  #coef for other confoundings (mediator prediction mod)
  conftailidx <- length(mod2$coefficients)
  if(conftailidx >= 3){
    beta2 <- mod2$coefficients[3:conftailidx]
  }else{
    beta2 <- NULL
  }
  
  #no interaction mod coefs
  mod11_spon_coef <- as.numeric(c(theta11, sd.theta11, pval.theta11, 
                                  theta21, sd.theta21, pval.theta21, 
                                  0, 0, NA, 
                                  beta1, sd.beta1, pval.beta1))
  
  #interaction mod coefs
  mod1_spon_coef <- as.numeric(c(theta1, sd.theta1, pval.theta1, 
                                 theta2, sd.theta2, pval.theta2, 
                                 theta3, sd.theta3, pval.theta3, 
                                 beta1, sd.beta1, pval.beta1))
  
  names(mod1_spon_coef) <- 
    names(mod11_spon_coef) <- 
    c('theta1', 'sd.theta1', 'pval.theta1', 
      'theta2', 'sd.theta2', 'pval.theta2', 
      'theta3', 'sd.theta3', 'pval.theta3', 
      'beta1', 'sd.beta1', 'pval.beta1')
  
  #(Residual standard error)^2 for mediator prediction model
  sigma2 <- (summary(mod2)$sigma)^2
  
  
  
  #Set all source value to 0 (a0) or 1 (a)
  a <- 1
  a0 <- 0
  
  
  
  continueconfs <- intersect(colnames(trainingdat), 
                             gsub(pattern = '`', replacement = '', 
                                  names(beta2)))
  
  discreteconfs <- setdiff(names(beta2), continueconfs)
  
  #Set discrete confoundings to 0 and continuous confoundings to their medians
  if(length(continueconfs) > 0){
    continueconfs <- trainingdat[,continueconfs, drop = FALSE]
    continueconfs <- apply(X = continueconfs, MARGIN = 2, 
                           FUN = function(x){x <- rep(median(x), length(x))})
    row.names(continueconfs) <- row.names(trainingdat)
  }else{
    continueconfs <- NULL
  }
  
  if(length(discreteconfs) > 0){
    
    discreteconfs <- matrix(data = 0, nrow = nrow(trainingdat), 
                            ncol = length(discreteconfs))
    colnames(discreteconfs) <- setdiff(names(beta2), 
                                       intersect(colnames(trainingdat), 
                                                 gsub(pattern = '`', 
                                                      replacement = '', 
                                                      names(beta2))))
    row.names(discreteconfs) <- row.names(trainingdat)
    
  }else{
    discreteconfs <- NULL
  }
  
  confs <- cbind(discreteconfs, continueconfs)
  confs <- confs[,names(beta2)]
  
  c.all <- confs
  
  #Set mediator value to its median value
  med <- median(trainingdat$mediator)
  
  #For interaction mod, set mod_spon_coef as mod1_spon_coef, 
  #for no interaction mod, set mod_spon_coef as mod11_spon_coef
  
  if(interaction == TRUE){
    mod_spon_coef <- mod1_spon_coef
  }else{
    mod_spon_coef <- mod11_spon_coef
  }
  
  #Controlled Direct Effect (CDE) (interaction mod, 
  #need to consider interaction term, so theta3 is included)
  cde <- as.numeric((mod_spon_coef[['theta1']] + 
                       mod_spon_coef[['theta3']]*med)*(a - a0))
  
  #Natural Indirect Effect (NIE)
  #The change in the odds of target in assocaition with a defined 
  #change in source (a0 to a), while holding the meidator at the 
  #level it would have naturally been at when source is at a specific 
  #level (a) (interaction mod, need to consider interaction term, 
  #so theta3 is included)
  nie <- as.numeric((mod_spon_coef[['theta2']]*mod_spon_coef[['beta1']] + 
                       mod_spon_coef[['theta3']]*mod_spon_coef[['beta1']]*a)*(a - a0))
  
  
  #Natural Direct Effect (NDE)
  #The change in the odds of target in association with a defined 
  #change in source (a0 to a), while holding the mediator at the 
  #level it would have naturally been at when source is at the 
  #original level (a0) (interaction mod, need to consider interaction 
  #term, so theta3 is included)
  
  confterm <- tryCatch({
    c.all %*% beta2
  }, error = function(err){
    0
  })
  
  
  if(targettype == 'continuous'){
    
    nde <- as.numeric(
      (mod_spon_coef[['theta1']] + 
         mod_spon_coef[['theta3']]*(beta0 + 
                                      mod_spon_coef[['beta1']]*a0 + 
                                      #mod_spon_coef[['beta1']]*a + 
                                      confterm))*(a - a0)
    )
    
    if(unique(nie*nde) > 0){
      pct.med <- nie/(nde + nie)
    }else{
      pct.med <- abs(nie/nde)
    }
    
    
    
  }else if(targettype == 'binary'){
    
    nde <- as.numeric(
      (mod_spon_coef[['theta1']] + 
         mod_spon_coef[['theta3']]*(beta0 + 
                                      mod_spon_coef[['beta1']]*a0 + 
                                      #mod_spon_coef[['beta1']]*a + 
                                      confterm + 
                                      mod_spon_coef[['theta2']]*sigma2))*(a - a0) + 
        0.5*mod_spon_coef[['theta3']]^2*sigma2*(a^2 - a0^2)
    )
    
    #pct.med <- exp(nde)*(exp(nie) - 1)/(exp(nde)*exp(nie) - 1)
    if(unique(nie*nde) > 0){
      pct.med <- nie/(nde + nie)
    }else{
      pct.med <- abs(nie/nde)
    }
    
  }
  
  
  ## bootstrapping to estimate the standard errors
  nboot <- bootnum
  theta11.boot <- theta21.boot <- 
    theta1.boot <- theta2.boot <- theta3.boot <- 
    beta0.boot <- beta1.boot <- 
    sigma2.boot <- 
    a0.boot <- a.boot <- med.boot <- cde.boot <- nie.boot <- rep(0, nboot)
  if(!is.null(confs)){
    beta2.boot <- mat.or.vec(nboot, ncol(confs))
  }else{
    beta2.boot <- NULL
  }
  
  #mat.or.vec create an nr by nc zero matrix if nc is greater than 1, 
  #and a zero vector of legnth nr if nc equals 1.
  nde.boot <- pct.med.boot <- mat.or.vec(nboot, nrow(trainingdat))
  #n is sample number
  
  B <- 1
  for(B in 1:nboot){
    set.seed(B)
    trainingdat.boot <- trainingdat[sample(1:nrow(trainingdat), 
                                           nrow(trainingdat), 
                                           replace = TRUE), ]
    
    if(IPW == TRUE){
      fac.boot <- sum(trainingdat.boot$invwt)/sum(trainingdat$invwt)
      trainingdat.boot$invwt <- trainingdat.boot$invwt/fac.boot
      
    }
    
    if(targettype == 'binary'){
      
      #no interaction mod
      mod11.boot <- glm(mod11form, family = binomial(link = "logit"), data = trainingdat.boot)
      #interaction mod
      mod1.boot <- glm(mod1form, family = binomial(link = "logit"), data = trainingdat.boot)
      
    }else if(targettype == 'continuous'){
      
      #no interaction mod
      mod11.boot <- lm(mod11form, data = trainingdat.boot)
      #interaction mod
      mod1.boot <- lm(mod1form, data = trainingdat.boot)
      
    }
    
    #mediator prediction mod
    if(IPW == TRUE & length(confoundings) > 0){
      mod2.boot <- lm(mod2form, data = trainingdat.boot, weights = invwt)
    }else{
      mod2.boot <- lm(mod2form, data = trainingdat.boot)
    }
    
    #coef for source (interaction mod)
    theta1.boot[B] <- mod1.boot$coefficients[2]
    #coef for mediator (interaction mod)
    theta2.boot[B] <- mod1.boot$coefficients[3]
    #coef for interaction (interaction mod)
    interidx <- length(mod1.boot$coefficients)
    theta3.boot[B] <- mod1$coefficients[interidx]
    
    #coef for source (no interactio mod)
    theta11.boot[B] <- mod11.boot$coefficients[2]
    #coef for mediator (no interaction mod)
    theta21.boot[B] <- mod11.boot$coefficients[3]
    
    #intercept (mediator prediction mod)
    beta0.boot[B] <- mod2.boot$coefficients[1]
    #coef for source (mediator prediction mod)
    beta1.boot[B] <- mod2.boot$coefficients[2]
    
    #coef for other confoundings (mediator prediction mod)
    conftailidx <- length(mod2.boot$coefficients)
    
    if(conftailidx >= 3){
      
      beta2boot <- mod2.boot$coefficients[3:conftailidx]
      remains <- setdiff(names(beta2), names(beta2boot))
      completes <- c(names(beta2boot), remains)
      beta2boot <- c(beta2boot, rep(0, length(remains)))
      names(beta2boot) <- completes
      beta2boot <- beta2boot[names(beta2)]
      beta2.boot[B,] <- beta2boot
      
    }
    
    #(Residual standard error)^2 for mediator prediction model
    sigma2.boot[B] <- (summary(mod2.boot)$sigma)^2
    
    
    
    a.boot[B] <- 1
    a0.boot[B] <- 0
    
    
    
    #Set discrete confoundings to 0 and continuous confoundings to their medians
    
    if(!is.null(c.all)){
      if(!is.null(continueconfs)){
        continueconfs.boot <- trainingdat.boot[,colnames(continueconfs), 
                                               drop = FALSE]
        continueconfs.boot <- apply(X = continueconfs.boot, MARGIN = 2, 
                                    FUN = function(x){x <- rep(median(x), length(x))})
        row.names(continueconfs.boot) <- row.names(trainingdat.boot)
      }else{
        continueconfs.boot <- NULL
      }
      
      if(!is.null(discreteconfs)){
        discreteconfs.boot <- matrix(data = 0, nrow = nrow(trainingdat.boot), 
                                     ncol = ncol(discreteconfs))
        colnames(discreteconfs.boot) <- colnames(discreteconfs)
        row.names(discreteconfs.boot) <- row.names(trainingdat.boot)
        
      }else{
        discreteconfs.boot <- NULL
      }
      
      c.all.boot <- cbind(discreteconfs.boot, continueconfs.boot)
      c.all.boot <- c.all.boot[,colnames(c.all), drop = FALSE]
    }else{
      c.all.boot <- NULL
    }
    
    med.boot[B] <- median(trainingdat.boot$mediator)
    
    mod11_spon_coef.boot <- as.numeric(c(theta11.boot[B], 
                                         theta21.boot[B], 
                                         0, 
                                         beta1.boot[B]))
    mod1_spon_coef.boot <- as.numeric(c(theta1.boot[B], 
                                        theta2.boot[B], 
                                        theta3.boot[B], 
                                        beta1.boot[B]))
    
    names(mod11_spon_coef.boot) <- names(mod1_spon_coef.boot) <- 
      c('theta1.boot', 'theta2.boot', 'theta3.boot', 'beta1.boot')
    
    #For interaction mod, set mod_spon_coef as mod1_spon_coef, 
    #for no interaction mod, set mod_spon_coef as mo11_spon_coef
    
    if(interaction == TRUE){
      mod_spon_coef.boot <- mod1_spon_coef.boot
    }else{
      mod_spon_coef.boot <- mod11_spon_coef.boot
    }
    
    #Controlled Direct Effect (CDE) (interaction mod, need to consider interaction term, 
    #so theta3 is included)
    cde.boot[B] <- as.numeric((mod_spon_coef.boot[['theta1.boot']] + 
                                 mod_spon_coef.boot[['theta3.boot']]*med.boot[B])*
                                (a.boot[B] - a0.boot[B]))
    
    #Natural Indirect Effect (NIE)
    nie.boot[B] <- as.numeric((mod_spon_coef.boot[['theta2.boot']]*mod_spon_coef.boot[['beta1.boot']] + 
                                 mod_spon_coef.boot[['theta3.boot']]*mod_spon_coef.boot[['beta1.boot']]*a.boot[B])*
                                (a.boot[B] - a0.boot[B]))
    
    
    #Natural Direct Effect (NDE)
    if(!is.null(confterm)){
      confterm.boot <- c.all.boot %*% beta2.boot[B,]
    }else{
      confterm.boot <- 0
    }
    
    
    
    
    
    
    if(targettype == 'continuous'){
      
      nde.boot[B,] <- as.numeric(
        (mod_spon_coef.boot[['theta1.boot']] + 
           mod_spon_coef.boot[['theta3.boot']]*(beta0.boot[B] + 
                                                  mod_spon_coef.boot[['beta1.boot']]*a0.boot[B] + 
                                                  confterm.boot))*(a.boot[B] - a0.boot[B])
      )
      
      
      if(unique(nie.boot[B]*nde.boot[B,]) > 0){
        pct.med.boot[B,] <- nie.boot[B]/(nde.boot[B,] + nie.boot[B])
      }else{
        pct.med.boot[B,] <- abs(nie.boot[B]/nde.boot[B,])
      }
      
      
    }else if(targettype == 'binary'){
      
      nde.boot[B,] <- as.numeric(
        (mod_spon_coef.boot[['theta1.boot']] + 
           mod_spon_coef.boot[['theta3.boot']]*(beta0.boot[B] + 
                                                  mod_spon_coef.boot[['beta1.boot']]*a0.boot[B] + 
                                                  confterm.boot + 
                                                  mod_spon_coef.boot[['theta2.boot']]*sigma2.boot[B]))*
          (a.boot[B] - a0.boot[B]) + 
          0.5*mod_spon_coef.boot[['theta3.boot']]^2*sigma2.boot[B]*(a.boot[B]^2 - a0.boot[B]^2)
      )
      
      
      
      #pct.med.boot[B,] <- exp(nde.boot[B,])*(exp(nie.boot[B]) - 1)/(exp(nde.boot[B,])*exp(nie.boot[B]) - 1)
      
      if(unique(nie.boot[B]*nde.boot[B,]) > 0){
        pct.med.boot[B,] <- nie.boot[B]/(nde.boot[B,] + nie.boot[B])
      }else{
        pct.med.boot[B,] <- abs(nie.boot[B]/nde.boot[B,])
      }
      
    }
    
    
    
  }
  
  summarycols <- c(setdiff(colnames(trainingdat)[4:max(4, 
                                                       ncol(trainingdat))], 
                           c('invwt')), 
                   #'CDE', 'CDE_lower', 'CDE_upper', 
                   'NDE', 'NDE_lower', 'NDE_upper', 
                   'NIE', 'NIE_lower', 'NIE_upper', 
                   'Prop', 'Prop_lower', 'Prop_upper')
  summarycols <- summarycols[!is.na(summarycols)]
  summary.tab <- mat.or.vec(nrow(trainingdat), length(summarycols))
  colnames(summary.tab) <- summarycols
  rownames(summary.tab) <- row.names(confs)
  
  if(!is.null(c.all)){
    for(B in 1:ncol(c.all)){
      
      summary.tab[,B] <- c.all[,B]
      
    }
  }
  
  #summary.tab[,'CDE'] <- cde
  #summary.tab[,'CDE_lower'] <- cde.boot[order(cde.boot)][nboot*0.025]
  #summary.tab[,'CDE_upper'] <- cde.boot[order(cde.boot)][nboot*0.975]
  
  
  summary.tab[,'NDE'] <- nde
  summary.tab[,'NDE_lower'] <- apply(nde.boot, 2, function(x){x[order(x)][nboot*0.025]})
  summary.tab[,'NDE_upper'] <- apply(nde.boot, 2, function(x){x[order(x)][nboot*0.975]})
  
  
  summary.tab[,'NIE'] <- nie
  summary.tab[,'NIE_lower'] <- nie.boot[order(nie.boot)][nboot*0.025]
  summary.tab[,'NIE_upper'] <- nie.boot[order(nie.boot)][nboot*0.975]
  
  
  
  summary.tab[,'Prop'] <- pct.med
  
  summary.tab[,'Prop_lower'] <- apply(pct.med.boot, 2, function(x){x[order(x)][nboot*0.025]})
  summary.tab[,'Prop_upper'] <- apply(pct.med.boot, 2, function(x){x[order(x)][nboot*0.975]})
  
  
  
  summary.tab <- as.data.frame(summary.tab, stringsAsFactors = FALSE)
  summary.tab$mediator <- mediator
  summary.tab <- unique(summary.tab)
  row.names(summary.tab) <- 1:nrow(summary.tab)
  
  cat(paste0('mediator number = ', j, '\n'))
  
  return(summary.tab)
  
}

mediationmod <- function(mediatordat, 
                         predictdat, 
                         source, 
                         target, 
                         interaction = FALSE, 
                         bootnum = 100, 
                         IPW = TRUE, 
                         responselevels = NULL, 
                         sourcelevels = NULL, 
                         threads = 1){
  
  colnames(predictdat)[colnames(predictdat) == source] <- 'source'
  
  colnames(predictdat)[colnames(predictdat) == target] <- 'target'
  
  confoundings <- colnames(predictdat)
  confoundings <- confoundings[-match('source', confoundings)]
  confoundings <- confoundings[-match('target', confoundings)]
  if(length(confoundings) == 0){
    confoundingstr <- NULL
  }else{
    confoundingstr <- paste(confoundings, collapse = '` + `')
    confoundingstr <- paste0('`', confoundings, '`')
    confoundingstr <- paste(confoundingstr, collapse = ' + ')
  }
  
  
  if(!is.factor(predictdat[['source']]) & !is.numeric(predictdat[['source']])){
    
    if(is.null(sourcelevels)){
      sourcelevels <- predictdat[['source']]
      sourcelevels <- unique(sourcelevels)
      sourcelevels <- sourcelevels[order(sourcelevels)]
    }
    predictdat[['source']] <- factor(predictdat[['source']], 
                                     levels = sourcelevels, 
                                     ordered = TRUE)
    #predictdat[['source']] <- as.numeric(predictdat[['source']]) - 1
  }
  
  
  mod11str <- paste0('`', target, '` ~ `', source, '` + mediator', 
                     c(' + ', '')[c(!is.null(confoundingstr), is.null(confoundingstr))], 
                     confoundingstr)
  mod1str <- paste0('`', target, '` ~ `', source, '`*mediator', 
                    c(' + ', '')[c(!is.null(confoundingstr), is.null(confoundingstr))], 
                    confoundingstr)
  mod2str <- paste0('mediator ~ `', source, '`', 
                    c(' + ', '')[c(!is.null(confoundingstr), is.null(confoundingstr))], 
                    confoundingstr)
  
  mod11form <- as.formula(mod11str)
  mod1form <- as.formula(mod1str)
  mod2form <- as.formula(mod2str)
  
  #IPW model###
  
  if(IPW == TRUE & !is.null(confoundingstr)){
    
    ipwstr <- paste0('source ~ ', confoundingstr)
    ipwform <- as.formula(ipwstr)
    
    samplenames <- row.names(predictdat)
    `%>%` <- magrittr::`%>%`
    
    vartype <- judgevartype(varvals = predictdat$source)
    
    
    
    if(vartype == 'continuous'){
      
      #invwt_continuous_weightit <- WeightIt::weightit(ipwform, 
      #                                                data = predictdat, 
      #                                                stabilize = TRUE)
      
      
      suppressMessages(
        
        invwt_continuous_weightit <- tryCatch({
          
          WeightIt::weightit(ipwform, 
                             data = predictdat, 
                             stabilize = TRUE)
          
        }, warning = function(war){
          
          NULL
          
        })
        
      )
      
      
        
      if(is.null(invwt_continuous_weightit)){
        
        invwt_continuous_weightit <- WeightIt::weightit(ipwform, 
                                                        data = predictdat, 
                                                        stabilize = TRUE)
        
        invwt_continuous_weightit <- WeightIt::trim(invwt_continuous_weightit, 
                                                    at = 0.9)
        
        cat('Extreme weights have been trimed to avoid effective sample size reduction\n')
        
      }
      
      
      
      #summary(invwt_continuous_weightit)
      
      predictdat <- predictdat %>%
        plyr::mutate(invwt = invwt_continuous_weightit$weights)
      
    }else{
      
      #invwt_discrete_weightit <- WeightIt::weightit(ipwform, 
      #                                              data = predictdat, 
      #                                              estimand = "ATE", 
      #                                              method = "ps")
      
      
      suppressMessages(
        
        invwt_discrete_weightit <- tryCatch({
          
          WeightIt::weightit(ipwform, 
                             data = predictdat, 
                             estimand = "ATE", 
                             method = "ps")
          
        }, warning = function(war){
          
          NULL
          
        })
        
      )
      
      
        
      if(is.null(invwt_discrete_weightit)){
        
        invwt_discrete_weightit <- WeightIt::weightit(ipwform, 
                                                      data = predictdat, 
                                                      estimand = "ATE", 
                                                      method = "ps")
        
        invwt_discrete_weightit <- WeightIt::trim(invwt_discrete_weightit, 
                                                  at = 0.9)
        
        cat('Extreme weights have been trimed to avoid effective sample size reduction\n')
        
      }
      
      
      
      #summary(invwt_discrete_weightit)
      
      predictdat <- predictdat %>%
        plyr::mutate(invwt = invwt_discrete_weightit$weights)
    }
    
    
    
    row.names(predictdat) <- samplenames
    predictdat <- predictdat[complete.cases(predictdat),]
    predictdat <- predictdat[predictdat[,ncol(predictdat)] >= 0,]
    
  }
  
  colnames(predictdat)[colnames(predictdat) == 'source'] <- source
  colnames(predictdat)[colnames(predictdat) == 'target'] <- target
  
  ##Organize data
  mediatordat <- as.data.frame(mediatordat)
  mediatordat <- mediatordat[row.names(predictdat), , drop = FALSE]
  
  
  
  jseqs <- 1:ncol(mediatordat)
  
  if(threads == 1){
    
    reslist <- list()
    j <- 1
    
    for(j in jseqs){
      
      singleres <- singlemediation(j = j, 
                                   mediatordat = mediatordat, 
                                   predictdat = predictdat, 
                                   target = target, 
                                   source = source, 
                                   confoundings = confoundings, 
                                   interaction = interaction, 
                                   IPW = IPW, 
                                   responselevels = responselevels, 
                                   mod11form = mod11form, 
                                   mod1form = mod1form, 
                                   mod2form = mod2form)
      
      reslist[[j]] <- singleres
      
    }
    
  }else{
    
    judgevartype <- function(varvals){
      
      if(length(unique(varvals)) < 2){
        stop("The source variable must have at least two unique values.")
      }else if(length(unique(varvals)) == 2){
        vartype <- "binary"
      }else if(is.factor(varvals) | is.character(varvals)){
        vartype <- "multinomial"
      }else{
        vartype <- "continuous"
      }
      
      return(vartype)
      
    }
    
    
    singlemediation <- function(j, 
                                mediatordat, 
                                predictdat, 
                                target, 
                                source, 
                                confoundings, 
                                interaction = FALSE, 
                                IPW, 
                                responselevels, 
                                mod11form, 
                                mod1form, 
                                mod2form){
      
      mediator <- colnames(mediatordat)[j]
      mediatorval <- mediatordat[, j, drop = FALSE]
      names(mediatorval) <- 'mediator'
      trainingdat <- cbind(predictdat, mediatorval)
      trainingdat <- trainingdat[complete.cases(trainingdat),]
      
      if(IPW == TRUE & length(confoundings) > 0){
        
        trainingdat <- trainingdat[c(target, source, 'mediator', confoundings, 'invwt')]
        
      }else{
        
        trainingdat <- trainingdat[c(target, source, 'mediator', confoundings)]
        
      }
      
      if(!is.factor(trainingdat[[target]]) & !is.numeric(trainingdat[[target]])){
        
        if(is.null(responselevels)){
          responselevels <- trainingdat[[target]]
          responselevels <- unique(responselevels)
          responselevels <- responselevels[order(responselevels)]
        }
        trainingdat[[target]] <- factor(trainingdat[[target]], 
                                        levels = responselevels, 
                                        ordered = TRUE)
        #trainingdat[[target]] <- as.numeric(trainingdat[[target]]) - 1
      }
      
      targettype <- judgevartype(varvals = trainingdat[[target]])
      
      if(targettype == 'binary'){
        
        #no interaction mod
        mod11 <- glm(mod11form, family = binomial(link = "logit"), data = trainingdat)
        #interaction mod
        mod1 <- glm(mod1form, family = binomial(link = "logit"), data = trainingdat)
        
      }else if(targettype == 'continuous'){
        
        #no interaction mod
        mod11 <- lm(mod11form, data = trainingdat)
        #interaction mod
        mod1 <- lm(mod1form, data = trainingdat)
        
      }
      
      if(sum(is.na(mod11$coefficients)) > 0){
        return(NULL)
      }
      
      if(sum(is.na(mod1$coefficients)) > 0){
        return(NULL)
      }
      
      #mediator prediction mod
      if(IPW == TRUE & length(confoundings) > 0){
        mod2 <- lm(mod2form, data = trainingdat, weights = invwt)
      }else{
        mod2 <- lm(mod2form, data = trainingdat)
      }
      
      if(sum(is.na(mod2$coefficients)) > 0){
        return(NULL)
      }
      
      #coef for source (interaction mod)
      theta1 <- mod1$coefficients[2]
      sd.theta1 <- summary(mod1)$coefficients[2, 2]
      pval.theta1 <- summary(mod1)$coefficients[2, 4]
      
      #coef for mediator (interaction mod)
      theta2 <- mod1$coefficients[3]
      sd.theta2 <- summary(mod1)$coefficients[3, 2]
      pval.theta2 <- summary(mod1)$coefficients[3, 4]
      
      #coef for interaction (interaction mod)
      interidx <- length(mod1$coefficients)
      theta3 <- mod1$coefficients[interidx]
      sd.theta3 <- summary(mod1)$coefficients[interidx, 2]
      pval.theta3 <- summary(mod1)$coefficients[interidx, 4]
      
      #coef for source (no interactio mod)
      theta11 <- mod11$coefficients[2]
      sd.theta11 <- summary(mod11)$coefficients[2, 2]
      pval.theta11 <- summary(mod11)$coefficients[2, 4]
      
      #coef for mediator (no interaction mod)
      theta21 <- mod11$coefficients[3]
      sd.theta21 <- summary(mod11)$coefficients[3, 2]
      pval.theta21 <- summary(mod11)$coefficients[3, 4]
      
      #intercept (mediator prediction mod)
      beta0 <- mod2$coefficients[1]
      
      #coef for source (mediator prediction mod)
      beta1 <- mod2$coefficients[2]
      sd.beta1 <- summary(mod2)$coefficients[2, 2]
      pval.beta1 <- summary(mod2)$coefficients[2, 4]
      
      #coef for other confoundings (mediator prediction mod)
      conftailidx <- length(mod2$coefficients)
      if(conftailidx >= 3){
        beta2 <- mod2$coefficients[3:conftailidx]
      }else{
        beta2 <- NULL
      }
      
      #no interaction mod coefs
      mod11_spon_coef <- as.numeric(c(theta11, sd.theta11, pval.theta11, 
                                      theta21, sd.theta21, pval.theta21, 
                                      0, 0, NA, 
                                      beta1, sd.beta1, pval.beta1))
      
      #interaction mod coefs
      mod1_spon_coef <- as.numeric(c(theta1, sd.theta1, pval.theta1, 
                                     theta2, sd.theta2, pval.theta2, 
                                     theta3, sd.theta3, pval.theta3, 
                                     beta1, sd.beta1, pval.beta1))
      
      names(mod1_spon_coef) <- 
        names(mod11_spon_coef) <- 
        c('theta1', 'sd.theta1', 'pval.theta1', 
          'theta2', 'sd.theta2', 'pval.theta2', 
          'theta3', 'sd.theta3', 'pval.theta3', 
          'beta1', 'sd.beta1', 'pval.beta1')
      
      #(Residual standard error)^2 for mediator prediction model
      sigma2 <- (summary(mod2)$sigma)^2
      
      
      
      #Set all source value to 0 (a0) or 1 (a)
      a <- 1
      a0 <- 0
      
      
      
      continueconfs <- intersect(colnames(trainingdat), 
                                 gsub(pattern = '`', replacement = '', 
                                      names(beta2)))
      
      discreteconfs <- setdiff(names(beta2), continueconfs)
      
      #Set discrete confoundings to 0 and continuous confoundings to their medians
      if(length(continueconfs) > 0){
        continueconfs <- trainingdat[,continueconfs, drop = FALSE]
        continueconfs <- apply(X = continueconfs, MARGIN = 2, 
                               FUN = function(x){x <- rep(median(x), length(x))})
        row.names(continueconfs) <- row.names(trainingdat)
      }else{
        continueconfs <- NULL
      }
      
      if(length(discreteconfs) > 0){
        
        discreteconfs <- matrix(data = 0, nrow = nrow(trainingdat), 
                                ncol = length(discreteconfs))
        colnames(discreteconfs) <- setdiff(names(beta2), 
                                           intersect(colnames(trainingdat), 
                                                     gsub(pattern = '`', 
                                                          replacement = '', 
                                                          names(beta2))))
        row.names(discreteconfs) <- row.names(trainingdat)
        
      }else{
        discreteconfs <- NULL
      }
      
      confs <- cbind(discreteconfs, continueconfs)
      confs <- confs[,names(beta2)]
      
      c.all <- confs
      
      #Set mediator value to its median value
      med <- median(trainingdat$mediator)
      
      #For interaction mod, set mod_spon_coef as mod1_spon_coef, 
      #for no interaction mod, set mod_spon_coef as mod11_spon_coef
      
      if(interaction == TRUE){
        mod_spon_coef <- mod1_spon_coef
      }else{
        mod_spon_coef <- mod11_spon_coef
      }
      
      #Controlled Direct Effect (CDE) (interaction mod, 
      #need to consider interaction term, so theta3 is included)
      cde <- as.numeric((mod_spon_coef[['theta1']] + 
                           mod_spon_coef[['theta3']]*med)*(a - a0))
      
      #Natural Indirect Effect (NIE)
      #The change in the odds of target in assocaition with a defined 
      #change in source (a0 to a), while holding the meidator at the 
      #level it would have naturally been at when source is at a specific 
      #level (a) (interaction mod, need to consider interaction term, 
      #so theta3 is included)
      nie <- as.numeric((mod_spon_coef[['theta2']]*mod_spon_coef[['beta1']] + 
                           mod_spon_coef[['theta3']]*mod_spon_coef[['beta1']]*a)*(a - a0))
      
      
      #Natural Direct Effect (NDE)
      #The change in the odds of target in association with a defined 
      #change in source (a0 to a), while holding the mediator at the 
      #level it would have naturally been at when source is at the 
      #original level (a0) (interaction mod, need to consider interaction 
      #term, so theta3 is included)
      
      confterm <- tryCatch({
        c.all %*% beta2
      }, error = function(err){
        0
      })
      
      
      if(targettype == 'continuous'){
        
        nde <- as.numeric(
          (mod_spon_coef[['theta1']] + 
             mod_spon_coef[['theta3']]*(beta0 + 
                                          mod_spon_coef[['beta1']]*a0 + 
                                          #mod_spon_coef[['beta1']]*a + 
                                          confterm))*(a - a0)
        )
        
        if(unique(nie*nde) > 0){
          pct.med <- nie/(nde + nie)
        }else{
          pct.med <- abs(nie/nde)
        }
        
        
        
      }else if(targettype == 'binary'){
        
        nde <- as.numeric(
          (mod_spon_coef[['theta1']] + 
             mod_spon_coef[['theta3']]*(beta0 + 
                                          mod_spon_coef[['beta1']]*a0 + 
                                          #mod_spon_coef[['beta1']]*a + 
                                          confterm + 
                                          mod_spon_coef[['theta2']]*sigma2))*(a - a0) + 
            0.5*mod_spon_coef[['theta3']]^2*sigma2*(a^2 - a0^2)
        )
        
        #pct.med <- exp(nde)*(exp(nie) - 1)/(exp(nde)*exp(nie) - 1)
        if(unique(nie*nde) > 0){
          pct.med <- nie/(nde + nie)
        }else{
          pct.med <- abs(nie/nde)
        }
        
      }
      
      
      ## bootstrapping to estimate the standard errors
      nboot <- bootnum
      theta11.boot <- theta21.boot <- 
        theta1.boot <- theta2.boot <- theta3.boot <- 
        beta0.boot <- beta1.boot <- 
        sigma2.boot <- 
        a0.boot <- a.boot <- med.boot <- cde.boot <- nie.boot <- rep(0, nboot)
      if(!is.null(confs)){
        beta2.boot <- mat.or.vec(nboot, ncol(confs))
      }else{
        beta2.boot <- NULL
      }
      
      #mat.or.vec create an nr by nc zero matrix if nc is greater than 1, 
      #and a zero vector of legnth nr if nc equals 1.
      nde.boot <- pct.med.boot <- mat.or.vec(nboot, nrow(trainingdat))
      #n is sample number
      
      B <- 1
      for(B in 1:nboot){
        set.seed(B)
        trainingdat.boot <- trainingdat[sample(1:nrow(trainingdat), 
                                               nrow(trainingdat), 
                                               replace = TRUE), ]
        
        if(IPW == TRUE){
          fac.boot <- sum(trainingdat.boot$invwt)/sum(trainingdat$invwt)
          trainingdat.boot$invwt <- trainingdat.boot$invwt/fac.boot
          
        }
        
        if(targettype == 'binary'){
          
          #no interaction mod
          mod11.boot <- glm(mod11form, family = binomial(link = "logit"), data = trainingdat.boot)
          #interaction mod
          mod1.boot <- glm(mod1form, family = binomial(link = "logit"), data = trainingdat.boot)
          
        }else if(targettype == 'continuous'){
          
          #no interaction mod
          mod11.boot <- lm(mod11form, data = trainingdat.boot)
          #interaction mod
          mod1.boot <- lm(mod1form, data = trainingdat.boot)
          
        }
        
        #mediator prediction mod
        if(IPW == TRUE & length(confoundings) > 0){
          mod2.boot <- lm(mod2form, data = trainingdat.boot, weights = invwt)
        }else{
          mod2.boot <- lm(mod2form, data = trainingdat.boot)
        }
        
        #coef for source (interaction mod)
        theta1.boot[B] <- mod1.boot$coefficients[2]
        #coef for mediator (interaction mod)
        theta2.boot[B] <- mod1.boot$coefficients[3]
        #coef for interaction (interaction mod)
        interidx <- length(mod1.boot$coefficients)
        theta3.boot[B] <- mod1$coefficients[interidx]
        
        #coef for source (no interactio mod)
        theta11.boot[B] <- mod11.boot$coefficients[2]
        #coef for mediator (no interaction mod)
        theta21.boot[B] <- mod11.boot$coefficients[3]
        
        #intercept (mediator prediction mod)
        beta0.boot[B] <- mod2.boot$coefficients[1]
        #coef for source (mediator prediction mod)
        beta1.boot[B] <- mod2.boot$coefficients[2]
        
        #coef for other confoundings (mediator prediction mod)
        conftailidx <- length(mod2.boot$coefficients)
        
        if(conftailidx >= 3){
          
          beta2boot <- mod2.boot$coefficients[3:conftailidx]
          remains <- setdiff(names(beta2), names(beta2boot))
          completes <- c(names(beta2boot), remains)
          beta2boot <- c(beta2boot, rep(0, length(remains)))
          names(beta2boot) <- completes
          beta2boot <- beta2boot[names(beta2)]
          beta2.boot[B,] <- beta2boot
          
        }
        
        #(Residual standard error)^2 for mediator prediction model
        sigma2.boot[B] <- (summary(mod2.boot)$sigma)^2
        
        
        
        a.boot[B] <- 1
        a0.boot[B] <- 0
        
        
        
        #Set discrete confoundings to 0 and continuous confoundings to their medians
        
        if(!is.null(c.all)){
          if(!is.null(continueconfs)){
            continueconfs.boot <- trainingdat.boot[,colnames(continueconfs), 
                                                   drop = FALSE]
            continueconfs.boot <- apply(X = continueconfs.boot, MARGIN = 2, 
                                        FUN = function(x){x <- rep(median(x), length(x))})
            row.names(continueconfs.boot) <- row.names(trainingdat.boot)
          }else{
            continueconfs.boot <- NULL
          }
          
          if(!is.null(discreteconfs)){
            discreteconfs.boot <- matrix(data = 0, nrow = nrow(trainingdat.boot), 
                                         ncol = ncol(discreteconfs))
            colnames(discreteconfs.boot) <- colnames(discreteconfs)
            row.names(discreteconfs.boot) <- row.names(trainingdat.boot)
            
          }else{
            discreteconfs.boot <- NULL
          }
          
          c.all.boot <- cbind(discreteconfs.boot, continueconfs.boot)
          c.all.boot <- c.all.boot[,colnames(c.all), drop = FALSE]
        }else{
          c.all.boot <- NULL
        }
        
        med.boot[B] <- median(trainingdat.boot$mediator)
        
        mod11_spon_coef.boot <- as.numeric(c(theta11.boot[B], 
                                             theta21.boot[B], 
                                             0, 
                                             beta1.boot[B]))
        mod1_spon_coef.boot <- as.numeric(c(theta1.boot[B], 
                                            theta2.boot[B], 
                                            theta3.boot[B], 
                                            beta1.boot[B]))
        
        names(mod11_spon_coef.boot) <- names(mod1_spon_coef.boot) <- 
          c('theta1.boot', 'theta2.boot', 'theta3.boot', 'beta1.boot')
        
        #For interaction mod, set mod_spon_coef as mod1_spon_coef, 
        #for no interaction mod, set mod_spon_coef as mo11_spon_coef
        
        if(interaction == TRUE){
          mod_spon_coef.boot <- mod1_spon_coef.boot
        }else{
          mod_spon_coef.boot <- mod11_spon_coef.boot
        }
        
        #Controlled Direct Effect (CDE) (interaction mod, need to consider interaction term, 
        #so theta3 is included)
        cde.boot[B] <- as.numeric((mod_spon_coef.boot[['theta1.boot']] + 
                                     mod_spon_coef.boot[['theta3.boot']]*med.boot[B])*
                                    (a.boot[B] - a0.boot[B]))
        
        #Natural Indirect Effect (NIE)
        nie.boot[B] <- as.numeric((mod_spon_coef.boot[['theta2.boot']]*mod_spon_coef.boot[['beta1.boot']] + 
                                     mod_spon_coef.boot[['theta3.boot']]*mod_spon_coef.boot[['beta1.boot']]*a.boot[B])*
                                    (a.boot[B] - a0.boot[B]))
        
        
        #Natural Direct Effect (NDE)
        if(!is.null(confterm)){
          confterm.boot <- c.all.boot %*% beta2.boot[B,]
        }else{
          confterm.boot <- 0
        }
        
        
        
        
        
        
        if(targettype == 'continuous'){
          
          nde.boot[B,] <- as.numeric(
            (mod_spon_coef.boot[['theta1.boot']] + 
               mod_spon_coef.boot[['theta3.boot']]*(beta0.boot[B] + 
                                                      mod_spon_coef.boot[['beta1.boot']]*a0.boot[B] + 
                                                      confterm.boot))*(a.boot[B] - a0.boot[B])
          )
          
          
          if(unique(nie.boot[B]*nde.boot[B,]) > 0){
            pct.med.boot[B,] <- nie.boot[B]/(nde.boot[B,] + nie.boot[B])
          }else{
            pct.med.boot[B,] <- abs(nie.boot[B]/nde.boot[B,])
          }
          
          
        }else if(targettype == 'binary'){
          
          nde.boot[B,] <- as.numeric(
            (mod_spon_coef.boot[['theta1.boot']] + 
               mod_spon_coef.boot[['theta3.boot']]*(beta0.boot[B] + 
                                                      mod_spon_coef.boot[['beta1.boot']]*a0.boot[B] + 
                                                      confterm.boot + 
                                                      mod_spon_coef.boot[['theta2.boot']]*sigma2.boot[B]))*
              (a.boot[B] - a0.boot[B]) + 
              0.5*mod_spon_coef.boot[['theta3.boot']]^2*sigma2.boot[B]*(a.boot[B]^2 - a0.boot[B]^2)
          )
          
          
          
          #pct.med.boot[B,] <- exp(nde.boot[B,])*(exp(nie.boot[B]) - 1)/(exp(nde.boot[B,])*exp(nie.boot[B]) - 1)
          
          if(unique(nie.boot[B]*nde.boot[B,]) > 0){
            pct.med.boot[B,] <- nie.boot[B]/(nde.boot[B,] + nie.boot[B])
          }else{
            pct.med.boot[B,] <- abs(nie.boot[B]/nde.boot[B,])
          }
          
        }
        
        
        
      }
      
      summarycols <- c(setdiff(colnames(trainingdat)[4:max(4, 
                                                           ncol(trainingdat))], 
                               c('invwt')), 
                       #'CDE', 'CDE_lower', 'CDE_upper', 
                       'NDE', 'NDE_lower', 'NDE_upper', 
                       'NIE', 'NIE_lower', 'NIE_upper', 
                       'Prop', 'Prop_lower', 'Prop_upper')
      summarycols <- summarycols[!is.na(summarycols)]
      summary.tab <- mat.or.vec(nrow(trainingdat), length(summarycols))
      colnames(summary.tab) <- summarycols
      rownames(summary.tab) <- row.names(confs)
      
      if(!is.null(c.all)){
        for(B in 1:ncol(c.all)){
          
          summary.tab[,B] <- c.all[,B]
          
        }
      }
      
      #summary.tab[,'CDE'] <- cde
      #summary.tab[,'CDE_lower'] <- cde.boot[order(cde.boot)][nboot*0.025]
      #summary.tab[,'CDE_upper'] <- cde.boot[order(cde.boot)][nboot*0.975]
      
      
      summary.tab[,'NDE'] <- nde
      summary.tab[,'NDE_lower'] <- apply(nde.boot, 2, function(x){x[order(x)][nboot*0.025]})
      summary.tab[,'NDE_upper'] <- apply(nde.boot, 2, function(x){x[order(x)][nboot*0.975]})
      
      
      summary.tab[,'NIE'] <- nie
      summary.tab[,'NIE_lower'] <- nie.boot[order(nie.boot)][nboot*0.025]
      summary.tab[,'NIE_upper'] <- nie.boot[order(nie.boot)][nboot*0.975]
      
      
      
      summary.tab[,'Prop'] <- pct.med
      
      summary.tab[,'Prop_lower'] <- apply(pct.med.boot, 2, function(x){x[order(x)][nboot*0.025]})
      summary.tab[,'Prop_upper'] <- apply(pct.med.boot, 2, function(x){x[order(x)][nboot*0.975]})
      
      
      
      summary.tab <- as.data.frame(summary.tab, stringsAsFactors = FALSE)
      summary.tab$mediator <- mediator
      summary.tab <- unique(summary.tab)
      row.names(summary.tab) <- 1:nrow(summary.tab)
      
      cat(paste0('mediator number = ', j, '\n'))
      
      return(summary.tab)
      
    }
    
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    reslist <- foreach::foreach(j = jseqs,
                                #.export = ls(name = globalenv())) %dopar% {
                                .export = NULL) %dopar% {
                                  
                                  singlemediation(j = j, 
                                                  mediatordat = mediatordat, 
                                                  predictdat = predictdat, 
                                                  target = target, 
                                                  source = source, 
                                                  confoundings = confoundings, 
                                                  interaction = interaction, 
                                                  IPW = IPW, 
                                                  responselevels = responselevels, 
                                                  mod11form = mod11form, 
                                                  mod1form = mod1form, 
                                                  mod2form = mod2form)
                                }
    
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
  }
  
  if(length(reslist) == ncol(mediatordat)){
    names(reslist) <- colnames(mediatordat)
  }else if(length(reslist) > 0){
    names(reslist) <- colnames(mediatordat)[1:length(reslist)]
  }
  
  return(reslist)
  
  
}

mediationanno <- function(reslines, propna = TRUE){
  
  reslines$mediation <- FALSE
  reslines$complete <- FALSE
  reslines$inconsistent <- FALSE
  
  reslines$mediation[reslines$NIE_lower*reslines$NIE_upper > 0] <- TRUE
  reslines$complete[(reslines$NDE_lower*reslines$NDE_upper < 0) & 
                      reslines$mediation] <- TRUE
  reslines$inconsistent[(reslines$NDE_low*reslines$NDE_upper > 0) & 
                          (reslines$NIE_lower*reslines$NIE_upper > 0) & 
                          (reslines$NDE*reslines$NIE < 0)] <- TRUE
  
  reslines <- reslines[order(-reslines$mediation, 
                             -reslines$complete, 
                             -reslines$inconsistent, 
                             -abs(reslines$Prop)),]
  
  if(propna == TRUE){
    
    reslines$Prop[(reslines$mediation == FALSE) | 
                    (reslines$mediation == TRUE & reslines$complete == TRUE) | 
                    (reslines$mediation == TRUE & reslines$inconsistent == TRUE)] <- 
      NA
    reslines$Prop_lower[is.na(reslines$Prop)] <- NA
    reslines$Prop_upper[is.na(reslines$Prop)] <- NA
    
  }
  
  row.names(reslines) <- 1:nrow(reslines)
  
  return(reslines)
  
}

orgreslist <- function(reslist, 
                       targetname, 
                       sourcename, 
                       anno = TRUE, 
                       propna = TRUE){
  
  i <- 1
  for(i in 1:length(reslist)){
    
    resline <- reslist[[i]]
    
    resline$source <- sourcename
    resline$target <- targetname
    
    resline <- resline[,c('source', 'target', 'mediator', 
                          colnames(resline)[1:(ncol(resline) - 3)])]
    
    if(i == 1){
      reslines <- resline
    }else{
      reslines <- rbind(reslines, resline)
    }
    
    
  }
  
  row.names(reslines) <- 1:nrow(reslines)
  
  if(anno == TRUE){
    
    reslines <- mediationanno(reslines = reslines, propna = propna)
    
    qualifiedsubs <- subset(reslines, mediation == TRUE)
    
    if(nrow(qualifiedsubs) == 0){
      res <- list(completeres = reslines, 
                  subres = NULL)
    }else{
      
      row.names(qualifiedsubs) <- 1:nrow(qualifiedsubs)
      
      res <- list(completeres = reslines, 
                  subres = qualifiedsubs)
    }
    
  }else{
    
    res <- list(completeres = reslines, 
                subres = NULL)
  }
  
  return(res)
  
}

singleci <- function(sub){
  
  n <- nrow(sub)
  
  #sub$CDE <- mean(sub$CDE)
  sub$NDE <- mean(sub$NDE)
  sub$NIE <- mean(sub$NIE)
  #sub$TT <- mean(sub$TT)
  sub$Prop <- mean(sub$Prop)
  
  #sub$CDE_lower <- sqrt(mean(sub$CDE_lower^2))
  #sub$CDE_upper <- sqrt(mean(sub$CDE_upper^2))
  sub$NDE_lower <- sqrt(mean(sub$NDE_lower^2))
  sub$NDE_upper <- sqrt(mean(sub$NDE_upper^2))
  sub$NIE_lower <- sqrt(mean(sub$NIE_lower^2))
  sub$NIE_upper <- sqrt(mean(sub$NIE_upper^2))
  #sub$TT_lower <- sqrt(mean(sub$TT_lower^2))
  #sub$TT_upper <- sqrt(mean(sub$TT_upper^2))
  sub$Prop_lower <- sqrt(mean(sub$Prop_lower^2))
  sub$Prop_upper <- sqrt(mean(sub$Prop_upper^2))
  
  #sub$CDE_lower <- sub$CDE - sub$CDE_lower
  #sub$CDE_upper <- sub$CDE + sub$CDE_upper
  sub$NDE_lower <- sub$NDE - sub$NDE_lower
  sub$NDE_upper <- sub$NDE + sub$NDE_upper
  sub$NIE_lower <- sub$NIE - sub$NIE_lower
  sub$NIE_upper <- sub$NIE + sub$NIE_upper
  #sub$TT_lower <- sub$TT - sub$TT_lower
  #sub$TT_upper <- sub$TT + sub$TT_upper
  sub$Prop_lower <- sub$Prop - sub$Prop_lower
  sub$Prop_upper <- sub$Prop + sub$Prop_upper
  
  sub <- sub[, c('source', 'target', 'mediator', 
                 'NDE', 'NDE_lower', 'NDE_upper', 
                 'NIE', 'NIE_lower', 'NIE_upper', 
                 'Prop', 'Prop_lower', 'Prop_upper'), drop = FALSE]
  
  sub <- unique(sub)
  
  return(sub)
  
  
}

combineci <- function(cires){
  
  #cires$CDE_lower <- cires$CDE - cires$CDE_lower
  #cires$CDE_upper <- cires$CDE_upper - cires$CDE
  
  cires$NDE_lower <- cires$NDE - cires$NDE_lower
  cires$NDE_upper <- cires$NDE_upper - cires$NDE
  
  cires$NIE_lower <- cires$NIE - cires$NIE_lower
  cires$NIE_upper <- cires$NIE_upper - cires$NIE
  
  #cires$TT_lower <- cires$TT - cires$TT_lower
  #cires$TT_upper <- cires$TT_upper - cires$TT
  
  cires$Prop_lower <- cires$Prop - cires$Prop_lower
  cires$Prop_upper <- cires$Prop_upper - cires$Prop
  
  cires <- cires[order(cires$mediator, cires$source, cires$target),]
  
  combineres <- plyr::ddply(.data = cires, 
                            .variables = c('mediator', 'source', 'target'), 
                            .fun = singleci)
  
  row.names(combineres) <- 1:nrow(combineres)
  
  return(combineres)
  
}

getgeneanno <- function(inputs){
  
  annoentrezs <- annosymbols <- unique(inputs)
  
  suppressWarnings(annoentrezs <- annoentrezs[!is.na(as.numeric(annoentrezs))])
  suppressWarnings(annosymbols <- annosymbols[is.na(as.numeric(annosymbols))])
  
  if(length(annoentrezs) > 0){
    entrezannos <- geneanno(geneentrezs = annoentrezs)
    if(nrow(entrezannos) > 0){
      
      entrezannos <- subset(entrezannos, !is.na(ENTREZID) | !is.na(SYMBOL))
      
      entrezannos <- entrezannos[match(annoentrezs, 
                                       entrezannos$input)[!is.na(match(annoentrezs, 
                                                                       entrezannos$input))], , 
                                 drop = FALSE]
      entrezannos <- entrezannos[, 2:ncol(entrezannos), drop = FALSE]
      
      if(!is.null(entrezannos)){
        
        if(nrow(entrezannos) > 0){
          row.names(entrezannos) <- 1:nrow(entrezannos)
        }else{
          entrezannos <- NULL
        }
        
      }
      
    }
    
  }else{
    entrezannos <- NULL
  }
  
  if(length(annosymbols) > 0){
    symbolannos <- geneanno(genesymbols = annosymbols)
    if(nrow(symbolannos) > 0){
      
      symbolannos <- subset(symbolannos, !is.na(ENTREZID) | !is.na(SYMBOL))
      
      symbolannos <- symbolannos[match(annosymbols, 
                                       symbolannos$input)[!is.na(match(annosymbols, 
                                                                       symbolannos$input))], , 
                                 drop = FALSE]
      
      symbolannos <- symbolannos[, 2:ncol(symbolannos), drop = FALSE]
      
      if(!is.null(symbolannos)){
        
        if(nrow(symbolannos) > 0){
          row.names(symbolannos) <- 1:nrow(symbolannos)
        }else{
          symbolannos <- NULL
        }
        
      }
      
    }
  }else{
    symbolannos <- NULL
  }
  
  inputsannos <- unique(rbind(entrezannos, symbolannos))
  
  if(!is.null(inputsannos)){
    inputsannos <- inputsannos[(!is.na(inputsannos$ENTREZID)) | 
                                 (!is.na(inputsannos$SYMBOL)),]
    if(nrow(inputsannos) == 0){
      inputsannos <- NULL
    }else{
      row.names(inputsannos) <- 1:nrow(inputsannos)
    }
  }
  
  return(inputsannos)
  
}

singleorglassodat <- function(i = 1, 
                              responseset, 
                              pddat){
  
  if(is.vector(responseset)){
    responsevals <- responseset
  }else{
    
    responseset <- responseset[pddat$sampleid, , drop = FALSE]
    responsevals <- responseset[, i]
  }
  
  simpddat <- pddat
  simpddat$Response <- responsevals
  simpddat <- simpddat[complete.cases(simpddat),]
  
  return(simpddat)
  
}

dividecv <- function(i = 1, nfold = 5, totalpd, totalvar, factortonum = TRUE){
  
  if(is.factor(totalpd$Response) & factortonum == TRUE){
    totalpd$Response <- as.numeric(totalpd$Response)
  }
  
  trainlist  <- list()
  testlist <- list()
  
  innertestsize <- floor(nrow(totalpd)/nfold)
  innertestsizes <- rep(innertestsize, nfold - 1)
  innertestsizes <- c(innertestsizes, nrow(totalpd) - sum(innertestsizes))
  
  if(innertestsize <= 0 | nfold == 1 | nfold <= 0){
    
    reslist <- list(trainlist = trainlist, 
                    testlist = NULL)
    
    return(reslist)
    
  }
  
  innertestidxes <- rep(seq(1, nfold, 1), innertestsizes)
  set.seed(i)
  innertest_idxes <- sample(x = 1:nrow(totalpd), size = nrow(totalpd), replace = FALSE)
  innertest_idxes <- split(x = innertest_idxes, innertestidxes)
  
  j <- 1
  for(j in 1:length(innertest_idxes)){
    
    innertest_idx <- innertest_idxes[[j]]
    
    
    innertestpd <- totalpd[innertest_idx,]
    innertrainpd <- totalpd[-innertest_idx,]
    
    innertrainpd <- innertrainpd[order(innertrainpd$sampleid, innertrainpd$Response),]
    innertestpd <- innertestpd[order(innertestpd$sampleid, innertestpd$Response),]
    row.names(innertrainpd) <- 1:nrow(innertrainpd)
    row.names(innertestpd) <- 1:nrow(innertestpd)
    
    innertrainvar <- totalvar[innertrainpd$sampleid,]
    innertestvar <- totalvar[innertestpd$sampleid,]
    
    
    trainlist[[j]] <- list(innertrainvar = innertrainvar, 
                           innertrainpd = innertrainpd)
    
    testlist[[j]] <- list(innertestvar = innertestvar, 
                          innertestpd = innertestpd)
  }
  
  reslist <- list(trainlist = trainlist, 
                  testlist = testlist)
  
  return(reslist)
  
}

singleelasticnet <- function(trainvar, 
                             trainres, 
                             threads = 1, 
                             alpha = 0.5, 
                             seednum = 1, 
                             foldnum = 10, 
                             errortype = 'min', 
                             verbose = TRUE){
  
  #Elastic net on bootstrapped/SMOTE train dataset
  #library(glmnet)
  #library(doParallel)
  #library(foreach)
  
  a <- alpha
  #Tune the value of alpha through a line search with the parallelism
  l <- a
  
  if(threads > 1){
    
    parallelval <- TRUE
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
  }else{
    parallelval <- FALSE
  }
  
  set.seed(seednum)
  cv <- glmnet::cv.glmnet(x = trainvar, y = trainres, family = "gaussian", 
                          nfold = foldnum, type.measure = "deviance", 
                          paralle = parallelval, alpha = l)
  
  if(threads > 1){
    parallel::stopCluster(cl)
  }
  
  gc()
  
  
  cvms.1ses <- c(cv$cvm[cv$lambda == cv$lambda.1se])
  cvms.mins <- c(cv$cvm[cv$lambda == cv$lambda.min])
  lambda.1ses <- c(cv$lambda.1se)
  lambda.mins <- c(cv$lambda.min)
  alphas <- c(l)
  
  
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
  }
  
  #Elastic net##########
  
  elasticnetmodel <- glmnet::glmnet(x = trainvar, y = trainres, 
                                    family = "gaussian", 
                                    lambda = Lambda, 
                                    alpha = Alpha)
  modelcoefs <- coef(elasticnetmodel)
  modelcoefs <- as.matrix(modelcoefs)
  
  reslist <- list(elasticnetmodel = elasticnetmodel, modelcoefs = modelcoefs)
  
  return(reslist)
  
}

prescreenprobes <- function(oripd, orivar){
  
  cor.pval <- function(x, y){
    
    pvalue <- cor.test(x = x, y = y)$p.value
    
    return(pvalue)
  }
  
  if(ncol(oripd) > 1){
    orivar <- orivar[oripd$sampleid, , drop = FALSE]
    
    pccpvals <- apply(X = orivar, MARGIN = 2, 
                      FUN = cor.pval, y = oripd$Response)
  }else{
    orivar <- orivar[row.names(oripd), , drop = FALSE]
    pccpvals <- apply(X = orivar, MARGIN = 2, 
                      FUN = cor.pval, y = oripd[,1])
  }
  
  
  pccpvals <- pccpvals[!is.na(pccpvals)]
  
  if(length(pccpvals) == 0){
    return(orivar)
  }
  
  screenfeatures <- pccpvals[pccpvals < 0.05]
  
  if(length(screenfeatures) <= 1){
    return(orivar)
  }
  
  if(ncol(oripd) > 0){
    screenvar <- orivar[, names(screenfeatures), drop = FALSE]
  }
  
  
  
  return(screenvar)
}

scatterplotting <- function(comptab, 
                            featurenum = NULL, 
                            title, 
                            colorful = FALSE, 
                            titlesize = 15, 
                            textsize = 13, 
                            face = 'bold'){
  
  samplenum <- nrow(comptab)
  
  if(!is.null(featurenum)){
    subtitle <- paste0('(', featurenum, ' features on ', samplenum, ' samples)')
  }else{
    subtitle <- paste0('(', samplenum, ' samples)')
    
  }
  
  comptab <- as.data.frame(comptab, stringsAsFactors = FALSE)
  pcc <- cor(comptab$Prediction, comptab$True)
  pcc <- signif(pcc, 3)
  
  lmfit <- lm(formula = Prediction~True, data = comptab)
  inter <- as.vector(coefficients(lmfit))[1]
  slope <- as.vector(coefficients(lmfit))[2]
  
  if(inter > 0){
    equalsign <- ' + '
    absinter <- signif(inter, 3)
  }else if(inter < 0){
    equalsign <- ' - '
    absinter <- signif(abs(inter), 3)
  }else{
    equalsign <- ''
    absinter <- ''
  }
  
  form <- paste0('y = ', signif(slope, 3), 'x', equalsign, absinter)
  
  lm_eqn <- function(comptab){
    m <- lm(Prediction ~ True, comptab);
    eq <- substitute(~~bolditalic(R^2~"="~r2), 
                     list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }
  
  xrange <- range(comptab$True)
  yrange <- range(comptab$Prediction)
  xlen <- abs(xrange[2] - xrange[1])
  ylen <- abs(yrange[2] - yrange[1])
  x1 <- xrange[1] + xlen*1/4
  y1 <- yrange[1] + ylen*3/4
  x2 <- xrange[1] + xlen*3/4
  y2 <- yrange[1] + ylen*1/4
  
  
  myColor <- rep('blue', nrow(comptab))
  
  if(colorful == TRUE & nrow(comptab) > 50){
    
    myColor <- densCols(x = comptab$True, y = comptab$Prediction, 
                        colramp = colorRampPalette(rev(rainbow(10, end=4/6))))
    
  }
  
  comptab$denscolor <- myColor
  rgbmat <- t(col2rgb(myColor))
  rgbmat <- as.data.frame(rgbmat, stringsAsFactors = FALSE)
  comptab <- cbind(comptab, rgbmat)
  comptab <- comptab[order(-comptab$blue, comptab$red, comptab$green),]
  comptab1 <- subset(comptab, blue >= red)
  comptab2 <- subset(comptab, blue < red)
  comptab2 <- comptab2[order(-comptab2$blue, comptab2$red, -comptab2$green),]
  comptab <- rbind(comptab1, comptab2)
  
  
  p <- ggplot2::ggplot(comptab, ggplot2::aes(x=True, y=Prediction))
  
  print(
    p + ggplot2::geom_point(color=comptab$denscolor, position='jitter', size = 2.5) + 
      ggplot2::xlab('True values') + 
      ggplot2::ylab('Predicted values') + 
      ggplot2::ggtitle(title, subtitle = subtitle) + 
      ggplot2::geom_abline(slope = slope, intercept = inter, color = 'red', size = 1) + 
      ggplot2::geom_text(x = x1, y = y1, label = lm_eqn(comptab), parse = TRUE, 
                         size = 5, color = 'red') + 
      ggplot2::annotate('text', label = paste0('bolditalic("', form, '")'), 
                        x = x2, y = y2, size = 5, color = 'red', parse = TRUE) + 
      ggplot2::theme_bw() + 
      
      ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face),
                     plot.subtitle = ggplot2::element_text(size = textsize, face = face), 
                     axis.text = ggplot2::element_text(size = textsize, face = face), 
                     axis.title = ggplot2::element_text(size = textsize, face = face))
  )
  
  
}

heatmapplotting <- function(orivar, 
                            comptab, 
                            featurenames = NULL, 
                            title, 
                            textsize = 13, 
                            show_rownames = FALSE){
  
  if(!is.null(featurenames)){
    heatmapmatrix <- orivar[, featurenames, drop = FALSE]
  }else{
    heatmapmatrix <- orivar
  }
  
  heatmapmatrix <- t(heatmapmatrix)
  
  comptab <- as.data.frame(comptab, stringsAsFactors = FALSE)
  
  
  heatmappd <- comptab[order(comptab$Prediction),]
  heatmapmatrix <- heatmapmatrix[, row.names(heatmappd), drop = FALSE]
  colnames(heatmappd) <- c('Predicted Response', 'True Response')
  
  showrownames <- FALSE
  fontsizerow <- textsize - 1
  
  clusterrows <- TRUE
  features <- 'features'
  if(nrow(heatmapmatrix) == 1){
    
    clusterrows <- FALSE
    features <- 'feature'
    
  }
  
  
  
  if(nrow(heatmapmatrix) <= 20){
    showrownames <- TRUE
  }
  
  if(show_rownames == FALSE){
    showrownames <- FALSE
  }
  
  
  
  
  suppressMessages(
    
    tryCatch({
      
      print(
        
        pheatmap::pheatmap(heatmapmatrix, scale = 'row', 
                           color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
                           show_colnames = FALSE, show_rownames = showrownames, 
                           main = paste0(title, ' ', 
                                         nrow(heatmapmatrix), ' ', features, ' on ', 
                                         ncol(heatmapmatrix), ' sampels'), 
                           annotation_col = heatmappd, 
                           cluster_cols = TRUE, 
                           cluster_rows = clusterrows, 
                           border_color = NA, 
                           fontsize = textsize, 
                           fontsize_row = fontsizerow)
        
      )
      
    }, error = function(err){
      
      print(
        
        pheatmap::pheatmap(heatmapmatrix[apply(heatmapmatrix, 1, var) != 0, , 
                                         drop = FALSE], 
                           scale = 'row', 
                           color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
                           show_colnames = FALSE, show_rownames = showrownames, 
                           main = paste0(title, ' ', 
                                         nrow(heatmapmatrix[apply(heatmapmatrix, 1, var) != 0, , 
                                                            drop = FALSE]), ' ', features, ' on ', 
                                         ncol(heatmapmatrix), ' sampels'), 
                           annotation_col = heatmappd, 
                           cluster_cols = TRUE, 
                           cluster_rows = clusterrows, 
                           border_color = NA, 
                           fontsize = textsize, 
                           fontsize_row = fontsizerow)
        
      )
      
    })
    
  )
  
  
}

singleprobeselection <- function(orivar, 
                                 oripd, 
                                 nfold = 5, 
                                 seednum = 1, 
                                 cores = 1, 
                                 alpha = 0.5, 
                                 errortype = 'min', 
                                 prescreen = TRUE, 
                                 crosstestrcut = 0.25, 
                                 plotting = FALSE, 
                                 plotprefix = NULL, 
                                 titlesize = 15, 
                                 textsize = 13, 
                                 face = 'bold', 
                                 verbose = TRUE){
  
  if(prescreen == TRUE){
    
    orivar <- prescreenprobes(oripd = oripd, orivar = orivar)
    
  }
  
  #Part1, divide the original dataset into train and test sets
  sampledivide <- dividecv(i = seednum, nfold = nfold, totalpd = oripd, totalvar = orivar)
  
  if(!is.null(sampledivide$testlist)){
    
    k <- 1
    for(k in 1:length(sampledivide$trainlist)){
      
      #Part2, elastic net on train dataset
      foldnum <- max(min(floor(nrow(sampledivide$trainlist[[k]]$innertrainvar)/10), 10), 3)
      
      elasticres <- singleelasticnet(trainvar = sampledivide$trainlist[[k]]$innertrainvar, 
                                     trainres = sampledivide$trainlist[[k]]$innertrainpd$Response, 
                                     threads = cores, 
                                     alpha = alpha, 
                                     seednum = seednum, 
                                     foldnum = foldnum, 
                                     errortype = errortype, 
                                     verbose = verbose)
      
      #Part3, summarize prediction on train and test datasets
      trainpre <-predict(elasticres$elasticnetmodel, 
                         s = elasticres$elasticnetmodel$lambda, 
                         newx = sampledivide$trainlist[[k]]$innertrainvar)
      
      testpre <- predict(elasticres$elasticnetmodel, 
                         s = elasticres$elasticnetmodel$lambda, 
                         newx = sampledivide$testlist[[k]]$innertestvar)
      
      traincomp <- cbind(trainpre, sampledivide$trainlist[[k]]$innertrainpd$Response)
      colnames(traincomp) <- c('Prediction', 'True')
      
      testcomp <- cbind(testpre, sampledivide$testlist[[k]]$innertestpd$Response)
      colnames(testcomp) <- c('Prediction', 'True')
      
      if(k == 1){
        traincomps <- traincomp
        testcomps <- testcomp
      }else{
        traincomps <- rbind(traincomps, traincomp)
        testcomps <- rbind(testcomps, testcomp)
      }
      
    }
    
    crosstestr <- cor(testcomps)[1, 2]^2
    
    if(is.na(crosstestr)){
      
      cat('The R square of the test from cross validation is NA, so drop the fitting\n')
      
      return(NULL)
      
    }else if(crosstestr < crosstestrcut){
      
      cat(paste0('The R square of the test from cross validation is ', 
                 signif(crosstestr, 3), ' < ', signif(crosstestrcut, 3), 
                 ' (the R square cutoff), so drop the fitting\n'))
      
      return(NULL)
      
    }else{
      
      if(verbose == TRUE){
        cat(paste0('The R square of the test from cross validation is ', 
                   signif(crosstestr, 3), ' >= ', signif(crosstestrcut, 3), 
                   ' (the R square cutoff), so continue the fitting\n'))
      }
      
    }
    
  }else{
    testcomps <- NULL
  }
  
  
  #Part2, elastic net on train dataset
  foldnum <- max(min(floor(nrow(orivar)/10), 10), 3)
  
  elasticres <- singleelasticnet(trainvar = orivar, 
                                 trainres = oripd$Response, 
                                 threads = cores, 
                                 alpha = alpha, 
                                 seednum = seednum, 
                                 foldnum = foldnum, 
                                 errortype = errortype, 
                                 verbose = verbose)
  
  #Part3, summarize prediction on train and test datasets
  oripre <- predict(elasticres$elasticnetmodel, 
                    s = elasticres$elasticnetmodel$lambda, 
                    newx = orivar)
  
  oricomp <- cbind(oripre, oripd$Response)
  colnames(oricomp) <- c('Prediction', 'True')
  
  
  
  if(plotting == TRUE){
    
    featurenum <- sum(elasticres$modelcoefs[-1,] != 0)
    
    if(is.null(plotprefix)){
      completeprefix <- 'Complete data'
      
      if(!is.null(testcomps)){
        testprefix <- paste0('Testing data from ', nfold, ' fold cross validation')
      }
      
    }else{
      plotprefix <- paste(toupper(substr(plotprefix, 1, 1)), 
                          substr(plotprefix, 2, nchar(plotprefix)), sep = '')
      completeprefix <- paste0(plotprefix, ' complete data')
      
      if(!is.null(testcomps)){
        testprefix <- paste0(plotprefix, ' testing data from ', nfold, 
                             ' fold cross validation')
      }
      
    }
    
    scatterplotting(comptab = oricomp, 
                    featurenum = featurenum, 
                    title = completeprefix, colorful = TRUE, 
                    titlesize = titlesize, 
                    textsize = textsize, 
                    face = face)
    
    
    if(!is.null(testcomps)){
      scatterplotting(comptab = testcomps, 
                      featurenum = NULL, 
                      title = testprefix, colorful = TRUE)
    }
    
    features <- elasticres$modelcoefs[-1,1]
    features <- features[features != 0]
    features <- features[order(features)]
    
    heatmapplotting(orivar = orivar, 
                    comptab = oricomp, 
                    featurenames = names(features), 
                    title = completeprefix, 
                    textsize = textsize)
    
    
  }
  
  
  res <- list(elasticmodel = elasticres$elasticnetmodel, 
              modelcoeffs = elasticres$modelcoefs, 
              completecomp = oricomp, 
              testcomp = testcomps)
  
  return(res)
  
  
}

extractlassoprobes <- function(res){
  
  if(!('plyr' %in% installed.packages()[,'Package'])){
    cat('Package plyr is needed to run this function\n')
    return(NULL)
  }
  
  if('modelcoeffs' %in% names(res)){
    
    probesub <- res$modelcoeffs
    probesub  <- probesub[complete.cases(probesub),]
    probesub <- as.matrix(probesub)
    probesub <- probesub[rowSums(probesub) != 0,]
    probesub <- as.matrix(probesub)
    probesub <- probesub[row.names(probesub) != '(Intercept)',1:ncol(probesub)]
    probesub <- as.matrix(probesub)
    
    probesub <- rowMeans(probesub)
    probesub <- as.matrix(probesub)
    Probes <- row.names(probesub)
    
    probesub <- as.data.frame(probesub, stringsAsFactors = FALSE)
    probesub$Probe <- Probes
    row.names(probesub) <- Probes
    
  }else if('modellist' %in% names(res)){
    
    if(length(res$modellist) >= 1){
      
      for(i in 1:length(res$modellist)){
        
        probes <- res$modellist[[i]]$modelcoeffs
        probesub  <- probesub[complete.cases(probesub),]
        probesub <- as.matrix(probesub)
        probes <- probes[rowSums(probes) != 0,]
        probes <- as.matrix(probes)
        probes <- probes[row.names(probes) != '(Intercept)',1:ncol(probes)]
        probes <- as.matrix(probes)
        
        probes <- rowMeans(probes)
        probes <- as.matrix(probes)
        Probes <- row.names(probesub)
        
        probes <- as.data.frame(probes, stringsAsFactors = FALSE)
        probesub$Probe <- Probes
        row.names(probesub) <- Probes
        
        if(i == 1){
          probesub <- probes
        }else{
          probesub <- rbind(probesub, probes)
        }
        
      }
      
    }else{
      
      return(NULL)
      
    }
    
  }else{
    return(NULL)
  }
  
  row.names(probesub) <- 1:nrow(probesub)
  
  calsum <- function(block){
    
    probe <- unique(block$probe)
    subblock <- block[-grep('Probe', colnames(probesub))]
    subsum <- colSums(subblock)
    submatrix <- data.frame(subsum)
    submatrix <- t(submatrix)
    row.names(submatrix) <- probe
    submatrix <- as.data.frame(submatrix)
    return(submatrix)
    
  }
  
  probesub <- plyr::ddply(.data = probesub, .variables = c('Probe'), 
                          .fun = calsum)
  names(probesub)[2] <- 'Coeff'
  probesub <- probesub[order(-probesub$Coeff),]
  row.names(probesub) <- 1:nrow(probesub)
  
  return(probesub)
  
}

mediatorannotation <- function(mediators, 
                               platform = 450){
  
  mediatornames <- c(sub(pattern = '::.*$', replacement = '', 
                         mediators))
  mediatornames <- c(mediatornames, 
                     c(sub(pattern = '^.*::', replacement = '', 
                           mediators)))
  mediatornames <- unique(mediatornames)
  
  probeannores <- probeanno(platform = platform, probes = mediatornames)
  
  mediatorgeneannos <- getgeneanno(inputs = mediatornames)
  
  if(!is.null(probeannores)){
    
    probeannores <- probeannores[match(mediatornames, 
                                       probeannores$Probe)[!is.na(
                                         match(mediatornames, 
                                               probeannores$Probe))], , 
                                 drop = FALSE]
    
    probeannores <- unique(probeannores)
    
    row.names(probeannores) <- 1:nrow(probeannores)
    
    probegeneannos <- getgeneanno(inputs = c(probeannores$UCSC_RefGene_Name, 
                                             probeannores$ENTREZID))
    
  }else{
    probegeneannos <- NULL
  }
  
  geneannores <- unique(rbind(mediatorgeneannos, probegeneannos))
  
  res <- list()
  
  res$geneannores <- geneannores
  
  if(!is.null(probeannores)){
    res$probeannores <- probeannores
  }
  
  
  return(res)
  
}

#'Perform elastic net regression and/or mediation analysis 
#'
#'Perform elastic net regression and/or mediation analysis between different 
#'  omics. 
#'
#'@param i This function performs elastic net regression using one feature in 
#'  the data \code{candrnadats} as response, or performs mediation test using  
#'  it as the potential mediator. The data \code{candrnadats} is a matrix with 
#'  rows as samples and columns as features. The parameter \code{i} here is an 
#'  column index to select one feature from the multiple features of the data. 
#'  Default is 1.
#'@param candrnadats The data of the features. It is a matrix with each row as 
#'  a sample and each column as a feature. The row names are sample names and 
#'  the column names are feature names. The parameter \code{i} will serve as 
#'  the column index to select one column from this data as the response of 
#'  the elastic net regression or as the potential mediator of the mediation 
#'  test. Typically an RNA data matrix.
#'@param oripddats The phenotypic data of the samples. It is a data frame with 
#'  the first column named as "sampleid" and recording the sample names same 
#'  as the row names of the matrix \code{candrnadats} and that of the matrix 
#'  \code{candmethyldats}. It should also has another column as the response 
#'  variable in the mediation analysis, so that two mediation relationships 
#'  will be tested as "response (contained in this \code{oripddats}) -> RNA 
#'  feature (the feature in \code{candrnadats} selected by \code{i}) -> result 
#'  feature (from \code{candmethyldats})", and its opposite as "result feature 
#'  -> RNA feature -> response". The column name of the response variable need 
#'  to be transferred to another parameter \code{responsevarname}, so that it 
#'  could be identified as the response, but if it is not transferred, the 
#'  second column of this \code{oripddats} data frame will be the response. 
#'  For the mediation analysis, this response variable should be a continuous 
#'  or a binary variable. 
#'@param candmethyldats The function uses elastic net regression to fit one 
#'  feature in \code{candrnadats}, and the predictor features to fit it are 
#'  from this parameter \code{candmethyldats}. It is also a matrix with rows 
#'  representing samples and columns representing features. All its features 
#'  will be used in the elastic net model and then the model will select the 
#'  ones able to fit the one feature in \code{candrnadats}. In the mediation 
#'  models, the one feature in \code{candrnadats} will serve as the mediator 
#'  and it will pair with each feature selected by the elastic net to build 
#'  several mediation models. If the elastic net regression step is skipped by 
#'  the parameters \code{onlyelasticnet} and \code{onlymediation}, all the 
#'  features in \code{candmethyldats} will be used to pair with the one in 
#'  \code{candrnadats} to construct several mediation models. The data of this 
#'  \code{candmethyldats} is typically a DNA methylation matrix.
#'@param seednum The random seed number used for the elastic net regression to 
#'  define a 10-fold cross-validation, so an optimal regularization constant 
#'  lambda can be selected for the elastic net. Default is 2022. 
#'@param threads Threads number for parallelization, default is 1.
#'@param alpha The alpha parameter for elastic net regression, controlling the 
#'  balance of L1 and L2 penalties. Default is 1, meaning the model will only 
#'  use L1 penalty, so it is actually a LASSO model.
#'@param errortype The method to choose the regularization constant lambda for 
#'  the elastic net regression model. If it is set as "min", the lambda that 
#'  combines with \code{alpha} and gives the minimum 10-fold cross-validation 
#'  error will be used, and if it is set as "1ses", the lambda that combines 
#'  with \code{alpha} and gives the cross-validation error within one standard 
#'  error of the minimum will be used. The 10-fold cross-validation is defined 
#'  by the random seed \code{seednum}. The default value of \code{errortype} 
#'  is "min".
#'@param prescreen If this parameter is set as TRUE, before the elastic net 
#'  training, a preliminary screen on the features of \code{candmethyldats} 
#'  will be performed first by calculating Pearson correlation coefficients 
#'  to the response variable (the feature from \code{candrnadats}), and the 
#'  ones with such a p-value less than 0.05 will be selected for the elastic 
#'  net model training. Default is TRUE.
#'@param crossrcut If the elastic net regression step is followed by mediation 
#'  analysis by the parameters \code{onlyelasticnet} and \code{onlymediation}, 
#'  the elastic net model predicted response should have an R square greater 
#'  than this \code{crossrcut}, i.e., its Pearson correlation coefficient with 
#'  the original feature value from \code{candrnadats} should be greater than 
#'  the square root of this value. Default value is 0.5, so the corresponding 
#'  Pearson correlation coefficient should be greater than 0.707. If it is not 
#'  reached, the function will return NULL.
#'@param plotting If this parameter is TRUE, a scatter plot will be generated 
#'  to compare the elastic net predicted response value and the true response 
#'  value (the feature from \code{candrnadats}), and a heatmap will be plotted 
#'  to show the elastic net selected features for the fitting (the features in 
#'  \code{candmethyldats}). Default is FALSE.
#'@param plotprefix The characters that will appear in the title of the plots 
#'  when \code{plotting} is TRUE. Default is NULL. 
#'@param responsevarname Defines which column of \code{oripddats} serves as 
#'  the response variable in the mediation analysis. It is the column name. 
#'  If it is NULL, the second column of \code{oripddats} will be the response. 
#'  For the mediation analysis, this response variable should be a continuous 
#'  or a binary one. Default is NULL. 
#'@param confoundings Defines which columns of \code{oripddats} serve as the 
#'  confounding factors in the mediation analysis. It is a vector containing 
#'  the confounding column names. Its default is NULL, meaning no confoundings 
#'  need to be included in the mediation analysis.
#'@param responselevels If the response variable in the mediation analysis is 
#'  binary, it will be converted to a factor by the function. This parameter 
#'  is used to define the element level of the factor. It should be a vector 
#'  with the elements ordered following the factor level. It can also be NULL, 
#'  so that the level will follow the alpha-beta order of the factor elements. 
#'  Default is NULL. 
#'@param methylplatform If the matrix transferred to \code{candmethyldats} is 
#'  DNA methylation microarray data, this parameter is used to define the data 
#'  platform. Can be 27 (for 27k platform), 450 (for 450k platform), or 850 
#'  (for EPIC platform). Default is 450.
#'@param onlyelasticnet Whether the function need to perform the elastic net 
#'  model training only and skip the mediation analysis. Default is FALSE. 
#'@param onlymediation Whether the function need to perform the mediation test 
#'  only and skip the elastic net model training. Default is FALSE. 
#'@param returncompleteres Whether all the mediation test results need to be 
#'  returned, including the insignificant ones. If it is FALSE, only the final 
#'  significant results will be returned. Default is FALSE. 
#'@param propna For the insignificant mediation test relatinships, whether the 
#'  insignificant mediation proportion should be returned as it is or directly 
#'  labeled as NA because they are insignficant. Default is TRUE, meaning they 
#'  will be directly labeled as NA. 
#'@return By the parameters \code{onlyelasticnet} and \code{onlymediation}, if 
#'  only the elastic net regression is performed, a list will be returned with 
#'  slots named as "elasticnetcoefs" and "completecomp". The former records 
#'  the regression coefficients of the features from \code{candmethyldats} to 
#'  fit the feature in \code{candrnadats}, and the latter is the predicted and 
#'  ture values of that feature in \code{candrnadats}. If only the mediation 
#'  analysis is performed, a list with a slot named as "mediationres" will be 
#'  returned, which records the significant mediation relationships from the 
#'  mediation analysis. In addition, other slots named as "probeannores" and 
#'  "geneannores" may also included in this result list if the features in the 
#'  significant mediation relationships have definite probe and gene function 
#'  information. Besides, if the parameter \code{returncompleteres} is set as 
#'  TRUE, all the tested mediaiton relationship results will be returned in a 
#'  slot named "completeres", no matter they are significant or not. If both 
#'  the elastic net regresson and the mediation analysis are performed by the 
#'  parameters \code{onlyelasticnet} and \code{onlymediation}, all the slots 
#'  above will be returned, including "elasticnetcoefs", "mediaitonres", etc. 
#'  However, if no siginficant mediation relationship can be detected by the 
#'  mediation analysis step, or the elastic net regression does not pass the 
#'  R square cutoff defined by the parameter \code{crossrcut}, the final list 
#'  will not have the slots "mediationres", "probeannores", and "geneannores". 
#'@export
singlelassomediation <- function(i = 1, 
                                 candrnadats, 
                                 oripddats, 
                                 candmethyldats, 
                                 
                                 seednum = 2022, 
                                 threads = 1, 
                                 alpha = 1, 
                                 errortype = 'min', 
                                 prescreen = TRUE, 
                                 crossrcut = 0.5, 
                                 plotting = FALSE, 
                                 plotprefix = NULL, 
                                 
                                 responsevarname = NULL, 
                                 confoundings = NULL, 
                                 responselevels = NULL, 
                                 
                                 methylplatform = 450, 
                                 onlyelasticnet = FALSE, 
                                 onlymediation = FALSE, 
                                 returncompleteres = FALSE, 
                                 propna = TRUE){
  
  #returncompleteres <- FALSE
  
  if(is.null(responsevarname) & ncol(oripddats) >= 2){
    responsevarname <- names(oripddats)[2]
  }
  
  if(is.null(responsevarname)){
    return(NULL)
  }
  
  mediatordat <- candrnadats[, i, drop = FALSE]
  mediatordat <- mediatordat[oripddats$sampleid, , drop = FALSE]
  candidatedats <- candmethyldats[oripddats$sampleid,]
  singlegenefeatures <- data.frame(Probe = colnames(candmethyldats), Coeff = NA, 
                                   stringsAsFactors = FALSE)
  
  
  
  if(!is.factor(oripddats[[responsevarname]]) & !is.numeric(oripddats[[responsevarname]])){
    
    if(is.null(responselevels)){
      responselevels <- oripddats[[responsevarname]]
      responselevels <- unique(responselevels)
      responselevels <- responselevels[order(responselevels)]
    }
    
    Response <- factor(oripddats[[responsevarname]], 
                       levels = responselevels, 
                       ordered = TRUE)
  }else{
    
    Response <- oripddats[[responsevarname]]
    
  }
  
  Responsetype <- judgevartype(varvals = Response)
  
  if((Responsetype != 'binary') & (Responsetype != 'continuous')){
    
    cat('The mediation analysis only supports binary or continuous variables 
        but the current response is not, so the mediation analysis cannot be performed\n')
    
    if(onlymediation == TRUE & onlyelasticnet == FALSE){
      
      return(NULL)
      
    }else{
      
      onlyelasticnet <- TRUE
      onlymediation <- FALSE
      
    }
    
  }
  
  
  
  
  
  if(onlymediation == TRUE & onlyelasticnet == TRUE){
    onlymediation <- FALSE
    onlyelasticnet <- FALSE
  }
  
  if(onlymediation == FALSE){
    
    if(onlyelasticnet == TRUE){
      nfold <- 1
    }
    
    singlegenepd <- singleorglassodat(i = i, 
                                      responseset = candrnadats, 
                                      pddat = oripddats)
    
    singlegenelassores <- singleprobeselection(orivar = candmethyldats, 
                                               oripd = singlegenepd, 
                                               nfold = 1, 
                                               seednum = seednum, 
                                               cores = threads, 
                                               alpha = alpha, 
                                               errortype = errortype, 
                                               prescreen = prescreen, 
                                               plotting = plotting, 
                                               plotprefix = plotprefix, 
                                               verbose = FALSE)
    
    if(is.null(singlegenelassores)){
      return(NULL)
    }
    
    if(onlyelasticnet == TRUE){
      
      res <- list()
      res$elasticnetcoefs <- singlegenelassores$modelcoeffs
      res$completecomp <- singlegenelassores$completecomp
      
      return(res)
    }
    
    singlegenefeatures <- extractlassoprobes(res = singlegenelassores)
    
    candidatedats <- candidatedats[, singlegenefeatures$Probe, drop = FALSE]
    
    
    crosstestr <- cor(singlegenelassores$completecomp)[1, 2]^2
    
    if(is.na(crosstestr)){
      
      return(NULL)
      
    }else if(crosstestr < crossrcut){
      
      return(NULL)
      
    }
    
    
    
  }
  
  
  
  
  completeres <- list()
  completerevres <- list()
  j <- 1
  for(j in 1:ncol(candidatedats)){
    
    #cat(paste0('j = ', j, '\n'))
    
    candidatedat <- candidatedats[, j, drop = FALSE]
    candidatename <- colnames(candidatedats)[j]
    
    otherdat <- oripddats[, c(responsevarname, confoundings), drop = FALSE]
    otherdat <- cbind(otherdat, candidatedat)
    row.names(otherdat) <- row.names(mediatordat)
    
    
    mediationreslist <- mediationmod(mediatordat = mediatordat, 
                                     predictdat = otherdat, 
                                     source = candidatename, 
                                     target = responsevarname, 
                                     interaction = FALSE, 
                                     bootnum = 100, 
                                     IPW = TRUE, 
                                     responselevels = responselevels, 
                                     threads = threads)
    
    
    if(length(mediationreslist) == 0){
      
      mediationrevres <- list()
      
      subrevres <- NULL
      
    }else{
      
      mediationreslist <- mediationreslist[!unlist(lapply(mediationreslist, is.null))]
      
      if(length(mediationreslist) > 0){
        
        mediationres <- orgreslist(reslist = mediationreslist, 
                                   targetname = responsevarname, 
                                   sourcename = candidatename, 
                                   propna = propna)
        
        subres <- mediationres$subres
        
      }else{
        
        mediationres <- list()
        
        subres <- NULL
        
      }
      
    }
    
    
      
    
    mediationrevreslist <- mediationmod(mediatordat = mediatordat, 
                                        predictdat = otherdat, 
                                        source = responsevarname, 
                                        target = candidatename, 
                                        interaction = FALSE, 
                                        bootnum = 100, 
                                        IPW = TRUE, 
                                        sourcelevels = responselevels, 
                                        threads = threads)
    
    if(length(mediationrevreslist) == 0){
      
      mediationrevres <- list()
      
      subrevres <- NULL
      
    }else{
      
      mediationrevreslist <- mediationrevreslist[!unlist(lapply(mediationrevreslist, is.null))]
      
      
      if(length(mediationrevreslist) > 0){
        
        mediationrevres <- orgreslist(reslist = mediationrevreslist, 
                                      targetname = candidatename, 
                                      sourcename = responsevarname, 
                                      propna = propna)
        
        subrevres <- mediationrevres$subres
        
      }else{
        
        mediationrevres <- list()
        
        subrevres <- NULL
      }
      
      
    }
    
    
    
    
    if(!is.null(subres) & !is.null(subrevres)){
      sharedmediators <- intersect(subres$mediator, subrevres$mediator)
      subres <- subset(subres, !(mediator %in% sharedmediators))
      subrevres <- subset(subrevres, !(mediator %in% sharedmediators))
    }
    
    qualifies <- rbind(subres, subrevres)
    
    if(!is.null(qualifies)){
      if(nrow(qualifies) > 0){
        row.names(qualifies) <- 1:nrow(qualifies)
      }
    }
    
    if(j == 1){
      qualires <- qualifies
    }else{
      qualires <- rbind(qualires, qualifies)
    }
    
    if(returncompleteres == TRUE){
      completeres[[j]] <- mediationres$completeres
      completerevres[[j]] <- mediationrevres$completeres
    }
    
  }
  
  if(returncompleteres == TRUE){
    
    completeres <- do.call(rbind, completeres)
    completerevres <- do.call(rbind, completerevres)
    completeres <- unique(rbind(completeres, completerevres))
    
  }
  
  if(is.null(qualires)){
    
    if(returncompleteres == FALSE){
      return(NULL)
    }else{
      res <- list()
      res$completeres <- completeres
      return(res)
    }
    
  }
  
  if(!is.null(qualires)){
    if(nrow(qualires) > 0){
      row.names(qualires) <- 1:nrow(qualires)
      
      mediatorannores <-  mediatorannotation(mediators = qualires$mediator, 
                                             platform = methylplatform)
      
      targetannores <- mediatorannotation(mediators = qualires$target, 
                                          platform = methylplatform)
      
      sourceannores <- mediatorannotation(mediators = qualires$source, 
                                          platform = methylplatform)
      
      subannores <- unique(rbind(mediatorannores$geneannores, 
                                 targetannores$geneannores, 
                                 sourceannores$geneannores))
      subprobeannores <- unique(rbind(mediatorannores$probeannores, 
                                      targetannores$probeannores, 
                                      sourceannores$probeannores))
      
      finalannores <- subannores
      finalprobeannores <- subprobeannores
      
      
      if(!is.null(finalannores)){
        if(nrow(finalannores) > 0){
          row.names(finalannores) <- 1:nrow(finalannores)
        }
      }
      
      if(!is.null(finalprobeannores)){
        if(nrow(finalprobeannores) > 0){
          row.names(finalprobeannores) <- 1:nrow(finalprobeannores)
        }
      }
      
      res <- list(mediationres = qualires)
      
      res$elasticnetcoefs <- NULL
      res$rsquare <- NULL
      
      if(onlymediation == FALSE){
        
        res$elasticnetcoefs <- singlegenelassores$modelcoeffs
        res$rsquare <- crosstestr
        
      }
      
      if(!is.null(finalprobeannores)){
        
        if(nrow(finalprobeannores)){
          res$probeannores <- finalprobeannores
        }
        
      }
      
      if(!is.null(finalannores)){
        if(nrow(finalannores) > 0){
          res$geneannores <- finalannores
        }
      }
      
    }else{
      
      if(returncompleteres == FALSE){
        return(NULL)
      }else{
        res <- list()
        res$completeres <- completeres
        return(res)
      }
      
    }
    
  }else{
    
    if(returncompleteres == FALSE){
      return(NULL)
    }else{
      res <- list()
      res$completeres <- completeres
      return(res)
    }
    
  }
  
  if(returncompleteres == TRUE){
    
    if(!is.null(res)){
      res$completeres <- completeres
    }else{
      res <- list()
      res$completeres <- completeres
    }
    
  }
  
  return(res)
  
}

probeselection <- function(modelcoeflist){
  
  modelcoefs <- do.call(cbind, modelcoeflist)
  coefmean <- rowMeans(modelcoefs)
  coefplus <- rowSums(modelcoefs > 0)
  coefminus <- rowSums(modelcoefs < 0)
  coefsign <- coefplus - coefminus
  coefscore <- coefsign*coefmean
  coefscoresd <- sd(coefscore)
  coefscoremean <- mean(coefscore)
  coefscoreupper <- coefscoremean + 2*coefscoresd
  coefscorelower <- coefscoremean - 2*coefscoresd
  
  finalfeatures <- coefscore[coefscore < coefscorelower | coefscore > coefscoreupper]
  
  finalfeatures <- as.matrix(finalfeatures)
  finalfeatures <- finalfeatures[-1,]
  
  finalfeatures <- finalfeatures[order(-finalfeatures)]
  
  return(finalfeatures)
  
}

singlelinear <- function(finalfeatures, 
                         resamplevar, 
                         resampleresponse){
  
  #Model training on resampled data
  subresamplevar <- resamplevar[,finalfeatures, drop = FALSE]
  
  subresamplevar <- cbind(subresamplevar, resampleresponse$Response)
  colnames(subresamplevar)[ncol(subresamplevar)] <- 'Response'
  subresamplevar <- as.data.frame(subresamplevar)
  
  fit <- lm(Response ~ ., data = subresamplevar)
  
  return(fit)
  
}

weightlinears <- function(innertrainvar, 
                          innertrainres, 
                          modellist, 
                          rsquarecut = 0.5){
  
  featurenames <- names(modellist[[1]]$coefficients)[-1]
  
  if(length(featurenames) == 0){
    return(NULL)
  }
  
  featurenames <- gsub(pattern = "`", replacement = '', x = featurenames)
  subtrain <- innertrainvar[,featurenames, drop = FALSE]
  
  subtrain <- as.data.frame(subtrain, stringsAsFactors = FALSE)
  
  trainprelist <- list()

  for(m in 1:length(modellist)){
    subfit <- modellist[[m]]
    
    trainpre <- predict(subfit, subtrain)
    
    trainprelist[[m]] <- trainpre
    
  }
  
  trainpreres <- do.call(cbind, trainprelist)
  
  colnames(trainpreres) <- paste0('model_', 1:ncol(trainpreres))
  
  trainpccs <- cor(trainpreres, innertrainres)
  modelpccs <- trainpccs
  colnames(modelpccs) <- 'training'
  
  #Only keep base learners with an R^2 > 0.5 in training data
  Rsquare <- rsquarecut
  
  if(sum(modelpccs[,1] > sqrt(Rsquare)) == 0){
    return(NULL)
  }
  
  modellist <- modellist[modelpccs[,1] > sqrt(Rsquare)]
  trainpreres <- trainpreres[,modelpccs[,1] > sqrt(Rsquare), drop = FALSE]
  
  modelpccs <- modelpccs[modelpccs[,1] > sqrt(Rsquare), , drop = FALSE]
  
  if(!is.matrix(trainpreres)){
    
    trainpreres <- as.matrix(trainpreres)
    colnames(trainpreres) <- 'model_1'
    
    modelpccs <- t(as.matrix(modelpccs))
    rownames(modelpccs) <- 'model_1'
  }
  
  #Calcuate the weight for each base learner only from training data
  weights <- 0.5 * log(modelpccs[,1]^2/(1 - modelpccs[,1]^2))
  weights <- weights/sum(weights)
  weights <- as.matrix(weights)
  weights <- t(weights)
  
  
  #Ensemble the prediction from all base learners
  trainpreres <- trainpreres %*% t(weights)
  rownames(weights) <- 'weights'
  
  if(is.null(colnames(weights))){
    colnames(weights) <- 'model_1'
  }
  
  colnames(trainpreres) <- 'Prediction'
  
  traincomp <- cbind(trainpreres, innertrainres)
  colnames(traincomp)[2] <- 'True'
  
  reslist <- list(baselearners = modellist, 
                  baseweights = weights, 
                  traincomp = traincomp)
  
  return(reslist)
  
}

filtermediator <- function(mediationres){
  
  tags <- apply(X = mediationres[, c('source', 'target'), drop = FALSE], MARGIN = 1, 
                FUN = function(x){return(paste0(c(x[1], x[2])[order(c(x[1], x[2]))], 
                                                collapse = '___'))})
  tags <- as.vector(tags)
  duptags <- unique(tags[duplicated(tags)])
  
  filterredmediatonres <- mediationres[!(tags %in% duptags), , drop = FALSE]
  
  if(nrow(filterredmediatonres) == 0){
    return(NULL)
  }else{
    row.names(filterredmediatonres) <- 1:nrow(filterredmediatonres)
    return(filterredmediatonres)
  }
  
}

mediationheatmap <- function(orivar, 
                             comptab, 
                             annotab, 
                             featurenames, 
                             title, 
                             textsize = 13, 
                             show_rownames = TRUE, 
                             quantilecuts = NULL, 
                             cellwidth = NA, 
                             cellheight = NA){
  
  if(!is.null(featurenames)){
    heatmapmatrix <- orivar[,featurenames, drop = FALSE]
  }else{
    heatmapmatrix <- orivar
  }
  
  heatmapmatrix <- t(heatmapmatrix)
  
  heatmappd <- as.data.frame(comptab, stringsAsFactors = FALSE)
  
  heatmapmatrix <- heatmapmatrix[,row.names(heatmappd), drop = FALSE]
  
  showrownames <- show_rownames
  fontsizerow <- textsize - 1
  
  clusterrows <- TRUE
  features <- 'features'
  if(nrow(heatmapmatrix) == 1){
    
    clusterrows <- FALSE
    features <- 'feature'
    
  }
  
  if(nrow(heatmapmatrix) > 40){
    showrownames <- FALSE
  }
  
  
  if(!is.null(quantilecuts)){
    
    points <- quantile(heatmapmatrix, c(0.25, 0.75))
    
    midLength <- length(heatmapmatrix[heatmapmatrix >= points[1] & 
                                        heatmapmatrix <= points[2]])
    
    midColor <- colorRampPalette(colors = c('blue', 'white', 'red'))(midLength)
    
    smallLength <- length(heatmapmatrix[heatmapmatrix < points[1]])
    smallColor <- rep(midColor[1], smallLength)
    
    greatLegnth <- length(heatmapmatrix[heatmapmatrix > points[2]])
    greatColor <- rep(midColor[length(midColor)], greatLegnth)
    
    colorList <- c(smallColor, midColor, greatColor)
    
  }else{
    colorList <- colorRampPalette(colors = c('blue', 'white', 'red'))(100)
  }
  
  
  suppressMessages(
    
    tryCatch({
      
      print(
        
        pheatmap::pheatmap(heatmapmatrix, scale = 'row', 
                           color = colorList, 
                           show_colnames = FALSE, show_rownames = showrownames, 
                           main = paste0(title, ' mediation on ', 
                                         nrow(heatmapmatrix), ' ', features, ' and ', 
                                         ncol(heatmapmatrix), ' sampels'), 
                           annotation_col = heatmappd, 
                           annotation_row = annotab, 
                           cluster_cols = TRUE, 
                           cluster_rows = clusterrows, 
                           border_color = NA, 
                           fontsize = textsize, 
                           fontsize_row = fontsizerow, 
                           cellwidth = cellwidth, 
                           cellheight = cellheight)
        
      )
      
    }, error = function(err){
      
      print(
        
        pheatmap::pheatmap(heatmapmatrix[apply(heatmapmatrix, 1, var) != 0, , 
                                         drop = FALSE], 
                           scale = 'row', 
                           color = colorList, 
                           show_colnames = FALSE, show_rownames = showrownames, 
                           main = paste0(title, ' mediation on ', 
                                         nrow(heatmapmatrix[apply(heatmapmatrix, 1, var) != 0, , 
                                                            drop = FALSE]), ' ', features, ' and ', 
                                         ncol(heatmapmatrix), ' sampels'), 
                           annotation_col = heatmappd, 
                           annotation_row = annotab, 
                           cluster_cols = TRUE, 
                           cluster_rows = clusterrows, 
                           border_color = NA, 
                           fontsize = textsize, 
                           fontsize_row = fontsizerow, 
                           cellwidth = cellwidth, 
                           cellheight = cellheight)
        
      )
      
    })
    
  )
  
  heatmappd <- heatmappd[order(-heatmappd[,1], heatmappd[,2]), , drop = FALSE]
  
  suppressMessages(
    
    tryCatch({
      
      print(
        
        pheatmap::pheatmap(heatmapmatrix[,row.names(heatmappd), drop = FALSE], scale = 'row', 
                           color = colorList, 
                           show_colnames = FALSE, show_rownames = showrownames, 
                           main = paste0(title, ' mediation on ', 
                                         nrow(heatmapmatrix), ' ', features, ' and ', 
                                         ncol(heatmapmatrix), ' sampels'), 
                           annotation_col = heatmappd, 
                           annotation_row = annotab, 
                           cluster_cols = FALSE, 
                           cluster_rows = clusterrows, 
                           border_color = NA, 
                           fontsize = textsize, 
                           fontsize_row = fontsizerow, 
                           cellwidth = cellwidth, 
                           cellheight = cellheight)
        
      )
      
    }, error = function(err){
      
      print(
        
        pheatmap::pheatmap(heatmapmatrix[apply(heatmapmatrix, 1, var) != 0, 
                                         row.names(heatmappd), 
                                         drop = FALSE], 
                           scale = 'row', 
                           color = colorList, 
                           show_colnames = FALSE, show_rownames = showrownames, 
                           main = paste0(title, ' mediation on ', 
                                         nrow(heatmapmatrix[apply(heatmapmatrix, 1, var) != 0, , 
                                                            drop = FALSE]), ' ', features, ' and ', 
                                         ncol(heatmapmatrix), ' sampels'), 
                           annotation_col = heatmappd, 
                           annotation_row = annotab, 
                           cluster_cols = FALSE, 
                           cluster_rows = clusterrows, 
                           border_color = NA, 
                           fontsize = textsize, 
                           fontsize_row = fontsizerow, 
                           cellwidth = cellwidth, 
                           cellheight = cellheight)
        
      )
      
    })
    
  )
  
  
  
  
  
}

#'Perform multi-omics regulatory analysis
#'
#'Perform multi-omics regulatory analysis between different omics. 
#'
#'@param methyldat This function performs LASSO plus mediation analysis with 
#'  the features in the data \code{rnadat} and \code{methyldat}. The parameter 
#'  \code{methyldat} here accepts a matrix with rows as samples and columns as 
#'  features, where the row names are sample names and the column names are 
#'  feature names. For a feature from \code{rnadat}, the function will use the 
#'  LASSO method to select features from \code{methyldat} so that they can fit 
#'  the feature from \code{rnadat}. Then, in mediation analysis, that feature 
#'  will become the mediator and pair with the LASSO selected \code{methyldat} 
#'  features to build mediation models. The data of \code{methyldat} can be a 
#'  DNA methylation matrix. 
#'@param rnadat It is a matrix with each row as a sample and each column as a 
#'  feature. The row names are sample names and the column names are feature 
#'  names. For a feature from \code{rnadat}, the features in \code{methyldat} 
#'  will be used to fit it with LASSO regression. Then, in mediation analysis, 
#'  that feature will become the mediator and pair with the LASSO selected 
#'  \code{methyldat} features to build mediation models. Typically, the data 
#'  of \code{rnadat} is an RNA matrix.
#'@param pddats The phenotypic data of the samples. It is a data frame with 
#'  the first column named as "sampleid" and recording the sample names same 
#'  as the row names of the matrices \code{rnadat} and \code{methyldat}. It 
#'  should also have another column as the response variable in the mediation 
#'  analysis, so that two mediation relationships will be tested as "response 
#'  (contained in \code{pddat}) -> RNA feature (the feature in \code{rnadat}) 
#'  -> result feature (from \code{methyldat})", and its opposite as "result 
#'  feature -> RNA feature -> response". The column name of the response need 
#'  to be transferred to another parameter \code{responsevarname}, so that it 
#'  could be identified as the response, but if it is not transferred, the 
#'  second column of this \code{pddat} data frame will be the response. For 
#'  the mediation analysis, this response variable should be a continuous or a 
#'  binary variable. 
#'@param responsevarname Defines which column of \code{pddat} serves as the 
#'  response variable in the mediation analysis. It is the column name. If it 
#'  is NULL, the second column of \code{pddat} will be the response. For the 
#'  mediation analysis, this response variable should be a continuous or a 
#'  binary one. Default is NULL. 
#'@param responselevels If the response variable in the mediation analysis is 
#'  binary, it will be converted to a factor by the function. This parameter 
#'  is used to define the element level of the factor. It should be a vector 
#'  with the elements ordered following the factor level. It can also be NULL, 
#'  so that the level will follow the alpha-beta order of the factor elements. 
#'  Default is NULL. 
#'@param confoundings Defines which columns of the \code{pddat} data frame are 
#'  the confounding factors in the mediation analysis. It should be a vector 
#'  containing the confounding column names. Its default is NULL, meaning no 
#'  confoundings need to be included in the mediation analysis.
#'@param calldifffeatures The features in \code{methyldat} and \code{rnadat} 
#'  can be directly used for the LASSO plus mediation analysis, or they can be 
#'  filtered first with limma, and then only the differential features between 
#'  sample groups will be used for the LASSO and mediation analysis. If this 
#'  parameter is set as TRUE, the differential feature calling step will be 
#'  performed, and if it is set as FALSE, it will not be performed. Default is 
#'  FALSE.
#'@param topdiffrnafeaturenum If \code{calldifffeatures} is set as TRUE, this 
#'  parameter will be used to define how many top differential features in the 
#'  \code{rnadat} will be used for the LASSO plus mediation analysis. For each 
#'  of them, it will be used in the LASSO regression so that the differential 
#'  features in \code{methyldat} can fit it. Then, in the mediation analysis, 
#'  it will serve as the mediator, and will pair with each \code{methyldat} 
#'  feature selected by LASSO. For such a feature pair, 2 mediation directions 
#'  will be tested as "response (contained in \code{pddat}) -> RNA feature 
#'  (the feature in \code{rnadat}) -> result feature (from \code{methyldat})", 
#'  and its opposite as "result feature -> RNA feature -> response". Default 
#'  is 10.
#'@param balanceadj Whether to use the ensemble-based mode. When it is set 
#'  as TRUE, a sampling process will be performed on the samples to generate 
#'  several base learner datasets, and the sample groups of each base learner 
#'  set will be adjusted during the sampling so that they will have the same 
#'  size. Then, for the differential feature selection, limma will be used on 
#'  each base learner to call the differential features, and their results 
#'  will finally be ensembled together to get the result for the whole data. 
#'  If it is set as FALSE, normal limma will be performed. Similarly, for the 
#'  LASSO and mediation analysis steps, they can also be performed on the base 
#'  learner datasets with an ensemble mode, or directly on the whole data with 
#'  a normal mode, depending on this parameter. Default is FALSE.
#'@param samplingmethod When \code{balanceadj} is TRUE, this parameter will 
#'  be needed to determine how to sample the data to get the base learner sets 
#'  for ensemeble mode. If it is set as "updn", to make the sample groups have 
#'  the same sample size in a base leaner set, a large group will be randomly 
#'  down-sampled randomly, and a small group will be up-sampled with SMOTE 
#'  (synthetic minority over-sampling technique), so that the final group size 
#'  will be the total sample number/group number. If it is set as "up", the 
#'  samples in a small group will be up-sampled, so its sample size will be 
#'  the large group size. If this parameter is set as "dn", the samples in the 
#'  large group will be down-sampled so that its sample size will be the small 
#'  sample size. Default is "updn".
#'@param nround When \code{balanceadj} is TRUE, this parameter will be used 
#'  to determine the number of the base learner sets in ensemble-based mode. 
#'  Default is 10.
#'@param seed Random seed for the sampling process to generate base learner 
#'  sets in the ensemble-based mode. For the LASSO regression step, it also 
#'  defines the 10-fold cross-validation to choose the optimal regularization 
#'  constant lambda. The default value is 2022.
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
#'@param methylremovereg A GRanges object or a data frame with the genomic 
#'  coordinates of any regions that need to be excluded from the differential 
#'  feature selection step for \code{methyldat}, can be the genomic regions 
#'  related to any confounding factors. If a feature called by normal limma or 
#'  ensemble-based limma is located there, this feature will be removed. For 
#'  regions related to confounding factors, this region removal can help avoid 
#'  influence from the confoundings. However, because both normal limma and 
#'  ensemble-based limma can adjust confounding factors by themselves, this 
#'  further region removal process will have very limited effects. Hence, this 
#'  step is unnecessary, and can be skipped by setting this parameter as NULL. 
#'  If need to transfer a GRanges or data frame, the data frame should contain 
#'  columns named as "seqnames", "start", "end", and "strand", which help to 
#'  indicate the genomic coordinates of the regions. Each row of the data frame 
#'  should be a region. The default value is NULL.
#'@param rnaremovereg Similar to \code{methylremovereg}, but this parameter is 
#'  used for the genomic region removal for the limma result of \code{rnadat}, 
#'  not \code{methyldat}. Default is NULL.
#'@param ignorestrand When \code{methylremovereg} or \code{rnaremovereg} is 
#'  not NULL and some genomic regions need to be excluded from the analysis, 
#'  this parameter will be used to judge whether any features located in the 
#'  genomic regions. If it is FALSE, the strand information of the features 
#'  should match the regions' strand if they are judged as within the regions, 
#'  and if it is TRUE, the strand information will not be considered. Default 
#'  is TRUE. 
#'@param methylfeaturetype The type of the features in the methylation data. 
#'  It will be used when \code{methylremovereg} is not NULL and judge whether 
#'  any features located in the genomic regions to be removed. Can be "gene" 
#'  or "probe" (for probes in the DNA methylation data). If another parameter 
#'  \code{plot} is TRUE, it will also appear in the title of the volcano plot 
#'  generated by this function. Default value is "probe".
#'@param tssradius When \code{rnaremovereg} or \code{methylremovereg} is not 
#'  NULL and need to judge whether a gene is in the regions to be removed, 
#'  this parameter will be used. If the TSS of a gene is located within the 
#'  regions, or not within, but the distance between its TSS and the end of 
#'  any regions is less than the value of this \code{tssradius}, this gene 
#'  will be judged as located within the regions (if the strand information 
#'  required by \code{ignorestrand} is also fulfilled). Then, the gene will 
#'  be excluded from the analysis. Default is 1500 (for 1500 bp).
#'@param plot Whether need to plot a volcano plot to show the differential 
#'  features if the parameter \code{calldifffeatures} is TRUE. It can also 
#'  control whether need to plot a heatmap to show the significant mediation 
#'  relationships identified. 
#'@param titleprefix The prefix of the volcano plot title. It can be set as 
#'  any character string need to be shown in the title. Default is NULL.
#'@param labelnum In the volcano plot, the names of the top up-regulated and 
#'  top down-regulated features can be labeled. This number indicates how many 
#'  top features in the groups should be labeled. If it is NULL, no feature 
#'  name will be labeled. Can also be any number. Default is NULL.
#'@param isbetaval_methyl If the feature values in \code{methyldat} are log2 
#'  transformed, set it as FALSE, if they are not log2 transformed, i.e., they 
#'  are gene counts or methylation beta values, set this parameter as TRUE. It 
#'  will only be needed when \code{plot} is TRUE.
#'@param isbetaval_rna Similar to \code{isbetaval_methyl}, but it is used to 
#'  indicate whether the feature values in \code{rnadat} are log2 transformed 
#'  or gene counts, not for \code{methyldat}.
#'@param absxcutoff The cutoff on log2FC (for log2 transformed features) or 
#'  inter-group difference (for non-log2 transformed feature counts or beta 
#'  values) to judge whether a feature is significantly different between the 
#'  sample groups. Default is 0, meaning any log2FC or difference value will 
#'  not influence the judgment on feature significance, when the parameter 
#'  \code{calldifffeatures} is set as TRUE.
#'@param pvalcutoff When \code{calldifffeatures} is TRUE, this parameter can 
#'  define the p-val (or adjusted p-val) cutoff to judge whether a feature is 
#'  significantly different between the groups or not. Default is 0.05.
#'@param pvalcolname When \code{calldifffeatures} is TRUE, this parameter can 
#'  define which p-val should be used to judge the statistical significance of 
#'  the features. Can be "P.Val" or "adj.P.Val". Default is "adj.P.Val".
#'@param titlesize The font size of the volcano plot title. Default is 15.
#'@param textsize The font size of the volcano plot and heatmap texts. Default 
#'  is 13.
#'@param face The font face of the volcano plot texts. Default is "bold".  
#'@param annotextsize If \code{labelnum} is not NULL, this parameter will be 
#'  needed to set the font size of the gene names to be labeled in the volcano 
#'  plot. Default is 4.
#'@param errortype The method to choose the regularization constant lambda for 
#'  the LASSO regression. If it is "min", the lambda that gives the minimum 
#'  mean 10-fold cross-validation error will be used , and if it is "1ses", 
#'  the lambda that gives the cross-validation error within one standard error 
#'  of the minimum will be used. The 10-fold cross-validation is defined by 
#'  the random seed \code{seed}. The default value of this \code{errortype} is 
#'  "min".
#'@param prescreen If this parameter is TRUE, before the LASSO training, a 
#'  preliminary screen on the features of \code{methyldat} will be performed 
#'  by calculating Pearson correlation coefficients to the response variable 
#'  (the feature from \code{rnadat}), and the ones with such a p-value less 
#'  than 0.05 will be selected for the LASSO model training. Default is TRUE.
#'@param crossrcut Before using a feature in \code{rnadat} as the mediator of 
#'  mediation analysis, the LASSO model predicted response should have an R 
#'  square greater than this \code{crossrcut}, i.e., its Pearson correlation 
#'  coefficient with the original \code{rnadat} feature value should be larger 
#'  than the square root of this value. Default is 0.5, so the corresponding 
#'  Pearson correlation coefficient should be greater than 0.707. If it is not 
#'  reached, the mediation analysis for this \code{rnadat} feature will not be 
#'  performed.
#'@param methylplatform If the data of \code{methyldat} is methylation probe 
#'  data, this parameter will be used to indicate the platform of the data. 
#'  Can be set as 450 (for 450k platform), 850 (for EPIC platform), or 27 (27K 
#'  platform). Its default value is 450.
#'@return A list with several slots will be returned. One named "mediationres" 
#'  records the significant mediation relationships from mediation analysis. 
#'  Another named "elasticnetcoefs" records the regression coefficients of the 
#'  LASSO regression models. The "elasticnetrsquares" slot is the R square of 
#'  the LASSO regression models. The "annolist" slot contains the biological 
#'  information of the genes or methylation probes with significant mediation 
#'  relationships. If the analysis are performed in the ensemble mode, another 
#'  slot named "elasticnetweights" will also be returned, which shows the base 
#'  learner weights in the ensemble-based LASSO models. In addition, if the 
#'  parameter \code{calldifffeatures} is set as TRUE, the differential feature 
#'  results will also be returned, with the slots named as "sigg.rnalimma", 
#'  "sigg.methyllimma", "nonsigg.rnalimma", and "nonsigg.methyllimma", which 
#'  record the significantly and insignificantly differential features in the 
#'  data of \code{rnadat} and \code{methyldat}.
#'@export
lassomediation <- function(methyldat, 
                           rnadat, 
                           pddat, 
                           
                           responsevarname = NULL, 
                           responselevels = NULL, 
                           confoundings = NULL, 
                           
                           calldifffeatures = FALSE, 
                           topdiffrnafeaturenum = 10, 
                           balanceadj = TRUE, 
                           samplingmethod = 'updn', 
                           nround = 10,
                           seed = 2022, 
                           k = 5, 
                           threads = 1, 
                           method = 'Cauchy', 
                           
                           methylremovereg = NULL, 
                           rnaremovereg = NULL, 
                           ignorestrand = TRUE, 
                           methylfeaturetype = 'probe', 
                           tssradius = 1500, 
                           
                           plot = FALSE, 
                           titleprefix = NULL, 
                           labelnum = 5, 
                           isbetaval_methyl = TRUE, 
                           isbetaval_rna = FALSE, 
                           absxcutoff = 0, 
                           pvalcutoff = 0.05, 
                           pvalcolname = 'adj.P.Val', 
                           titlesize = 15, 
                           textsize = 13, 
                           face = 'bold', 
                           annotextsize = 4, 
                           
                           errortype = 'min', 
                           prescreen = TRUE, 
                           crossrcut = 0.5, 
                           
                           methylplatform = 450){
  
  if(nrow(methyldat) == 1){
    
    cat('The methylation data only contains one feature, so the function `singlelassomediation` 
        should be used to perform the mediation analysis instead (Set its parameter 
        `onlyelasticnet` as FALSE and `onlymediation` as TRUE)')
    
    return(NULL)
    
  }
  
  pddat <- subset(pddat, sampleid %in% unique(c(colnames(methyldat), 
                                                colnames(rnadat))))
  methyldat <- methyldat[, pddat$sampleid, drop = FALSE]
  rnadat <- rnadat[, pddat$sampleid, drop = FALSE]
  
  
  
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
    
    cat('The mediation analysis only supports binary or continuous variables 
        but the current response is not, so the analysis cannot continue\n')
    
    return(NULL)
    
  }
  
  if((Responsetype != 'binary') & (balanceadj == TRUE)){
    
    cat('The ensemble mode only supports binary variable, but the current response is not.\nHence, the mode has been changed to normal mode.\n')
    
    balanceadj <- FALSE
    
  }
  
  
  
  if(balanceadj == TRUE){
    
    methylbalanceddats <- balancesampling(dat = methyldat, 
                                          pddat = pddat, 
                                          responsename = responsevarname, 
                                          nround = nround, 
                                          seed = seed, 
                                          k = k, 
                                          threads = threads, 
                                          samplingmethod = samplingmethod)
    
    rnabalanceddats <- balancesampling(dat = rnadat, 
                                       pddat = pddat, 
                                       responsename = responsevarname, 
                                       nround = nround, 
                                       seed = seed, 
                                       k = k, 
                                       threads = threads, 
                                       samplingmethod = samplingmethod)
    
    
  }
  
  #limma
  if(calldifffeatures == TRUE){
    
    if(balanceadj == TRUE){
      
      methyllimmareslist <- list()
      
      iseqs <- seq(1, nround, 1)
      
      for(i in iseqs){
        
        methyllimmares <- diffsites(dats = methylbalanceddats$dats, 
                                    pddats = methylbalanceddats$pds, 
                                    i = i, 
                                    responsevarname = responsevarname, 
                                    responselevels = responselevels, 
                                    confoundings = confoundings)
        
        methyllimmareslist[[i]] <- methyllimmares
        
        
      }
      
      
      names(methyllimmareslist) <- names(methylbalanceddats$dats)
      
      methylcombineres <- combinepval(limmareslist = methyllimmareslist, 
                                      method = method)
      
      metamethyllimmares <- data.frame(logFC = methylcombineres$logfccombs, 
                                       P.Value = methylcombineres$pcombs, 
                                       adj.P.Val = methylcombineres$padjs, 
                                       stringsAsFactors = FALSE)
      
      methyllimmares <- metamethyllimmares
      
      
      
      rnalimmareslist <- list()
      
      for(i in iseqs){
        
        rnalimmares <- diffsites(dats = rnabalanceddats$dats, 
                                 pddats = rnabalanceddats$pds, 
                                 i = i, 
                                 responsevarname = responsevarname, 
                                 responselevels = responselevels, 
                                 confoundings = confoundings)
        
        rnalimmareslist[[i]] <- rnalimmares
        
        
      }
      
      
      names(rnalimmareslist) <- names(rnabalanceddats$dats)
      
      rnacombineres <- combinepval(limmareslist = rnalimmareslist, 
                                   method = method)
      
      metarnalimmares <- data.frame(logFC = rnacombineres$logfccombs, 
                                    P.Value = rnacombineres$pcombs, 
                                    adj.P.Val = rnacombineres$padjs, 
                                    stringsAsFactors = FALSE)
      
      rnalimmares <- metarnalimmares
      
      
    }else{
      
      methylbalanceddats <- list(dats = list(round_1 = methyldat), 
                                 pds = list(round_1 = pddat))
      
      methyllimmares <- diffsites(dats = methylbalanceddats$dats, 
                                  pddats = methylbalanceddats$pds, 
                                  i = 1, 
                                  responsevarname = responsevarname, 
                                  responselevels = responselevels, 
                                  confoundings = confoundings)
      
      
      rnabalanceddats <- list(dats = list(round_1 = rnadat), 
                              pds = list(round_1 = pddat))
      
      rnalimmares <- diffsites(dats = rnabalanceddats$dats, 
                               pddats = rnabalanceddats$pds, 
                               i = 1, 
                               responsevarname = responsevarname, 
                               responselevels = responselevels, 
                               confoundings = confoundings)
      
    }
    
    if(!is.null(methylremovereg)){
      
      if(methylfeaturetype == 'probe'){
        methylremovedfeatures <- overlappedprobes(totalprobes = row.names(methyllimmares), 
                                                  removereg = methylremovereg, 
                                                  platform = methylplatform, 
                                                  ignorestrand = ignorestrand)
      }else if(methylfeaturetype == 'gene'){
        methylremovedfeatures <- overlappedgenes(totalgenes = row.names(methyllimmares), 
                                                 removereg = methylremovereg, 
                                                 tssradius = tssradius, 
                                                 ignorestrand = ignorestrand)
      }else{
        methylremovedfeatures <- NULL
      }
      
      
      if(!is.null(methylremovedfeatures)){
        
        if(methylfeaturetype == 'probe'){
          methylremovedfeatures <- unique(row.names(methylremovedfeatures))
          
          methyllimmares$remove <- FALSE
          methyllimmares$remove[row.names(methyllimmares) %in% methylremovedfeatures] <- TRUE
          
        }else if(methylfeaturetype == 'gene'){
          
          methylremovegenes <- idxgenenames(totalgenes = row.names(methyllimmares), 
                                            genesym = methylremovedfeatures$genesym, 
                                            geneid = methylremovedfeatures$geneid)
          
          methyllimmares$remove <- methylremovegenes
          
        }else{
          
          methyllimmares$remove <- FALSE
        }
        
      }else{
        methyllimmares$remove <- FALSE
      }
      
    }
    
    if(!is.null(rnaremovereg)){
      
      rnaremovedfeatures <- overlappedgenes(totalgenes = row.names(rnalimmares), 
                                            removereg = rnaremovereg, 
                                            tssradius = tssradius, 
                                            ignorestrand = ignorestrand)
      
      
      if(!is.null(rnaremovedfeatures)){
        
        rnaremovegenes <- idxgenenames(totalgenes = row.names(rnalimmares), 
                                       genesym = rnaremovedfeatures$genesym, 
                                       geneid = rnaremovedfeatures$geneid)
        
        rnalimmares$remove <- rnaremovegenes
        
        
        
      }else{
        rnalimmares$remove <- FALSE
      }
      
    }
    
    if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val', 
                                   'pvalue', 'p-value', 'p.value')){
      pvalcolname <- 'P.Value'
    }else{
      pvalcolname <- 'adj.P.Val'
    }
    
    #otherpvalcolname <- setdiff(c('P.Value', 'adj.P.Val'), pvalcolname)
    
    methyllimmares <- methyllimmares[order(methyllimmares[,pvalcolname], 
                                           #methyllimmares[otherpvalcolname], 
                                           -abs(methyllimmares$logFC)), , drop = FALSE]
    
    rnalimmares <- rnalimmares[order(rnalimmares[,pvalcolname], 
                                     #rnalimmares[otherpvalcolname], 
                                     -abs(rnalimmares$logFC)), , drop = FALSE]
    
    
    
    if('remove' %in% names(methyllimmares)){
      
      sigg.methyllimma <- methyllimmares[(methyllimmares[,pvalcolname] < pvalcutoff) & 
                                           (abs(methyllimmares$logFC) > absxcutoff) & 
                                           (methyllimmares$remove == FALSE),]
      nonsigg.methyllimma <- methyllimmares[(methyllimmares[,pvalcolname] >= pvalcutoff) | 
                                              (abs(methyllimmares$logFC) <= absxcutoff) | 
                                              (methyllimmares$remove == TRUE),]
      
    }else{
      
      sigg.methyllimma <- methyllimmares[(methyllimmares[,pvalcolname] < pvalcutoff) & 
                                           (abs(methyllimmares$logFC) > absxcutoff),]
      nonsigg.methyllimma <- methyllimmares[(methyllimmares[,pvalcolname] >= pvalcutoff) | 
                                              (abs(methyllimmares$logFC) <= absxcutoff),]
      
    }
    
    if('remove' %in% names(rnalimmares)){
      
      sigg.rnalimma <- rnalimmares[(rnalimmares[,pvalcolname] < pvalcutoff) & 
                                     (abs(rnalimmares$logFC) > absxcutoff) & 
                                     (rnalimmares$remove == FALSE),]
      nonsigg.rnalimma <- rnalimmares[(rnalimmares[,pvalcolname] >= pvalcutoff) | 
                                        (abs(rnalimmares$logFC) <= absxcutoff) | 
                                        (rnalimmares$remove == TRUE),]
      
    }else{
      
      sigg.rnalimma <- rnalimmares[(rnalimmares[,pvalcolname] < pvalcutoff) & 
                                     (abs(rnalimmares$logFC) > absxcutoff),]
      nonsigg.rnalimma <- rnalimmares[(rnalimmares[,pvalcolname] >= pvalcutoff) | 
                                        (abs(rnalimmares$logFC) <= absxcutoff),]
      
    }
    
    if(plot == TRUE){
      
      if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val', 
                                     'pvalue', 'p-value', 'p.value')){
        pvalcolname <- 'P.Value'
      }else{
        pvalcolname <- 'adj.P.Val'
      } 
      
      methylplotdat <- singlevolcano(limmares = methyllimmares, 
                                     betaval = isbetaval_methyl, 
                                     absxcutoff = absxcutoff, 
                                     pvalcut = pvalcutoff, 
                                     pvalcolname = pvalcolname, 
                                     confoundings = confoundings, 
                                     responsevarname = responsevarname, 
                                     titleprefix = paste0(titleprefix, ' methylation'), 
                                     featuretype = methylfeaturetype, 
                                     titlesize = titlesize, 
                                     textsize = textsize, 
                                     face = face, 
                                     labelnum = labelnum, 
                                     annotextsize = annotextsize)
      
      rnaplotdat <- singlevolcano(limmares = rnalimmares, 
                                  betaval = isbetaval_rna, 
                                  absxcutoff = absxcutoff, 
                                  pvalcut = pvalcutoff, 
                                  pvalcolname = pvalcolname, 
                                  confoundings = confoundings, 
                                  responsevarname = responsevarname, 
                                  titleprefix = paste0(titleprefix, ' RNA'), 
                                  featuretype = 'gene', 
                                  titlesize = titlesize, 
                                  textsize = textsize, 
                                  face = face, 
                                  labelnum = labelnum, 
                                  annotextsize = annotextsize)
      
      
    }
    
    if(nrow(sigg.methyllimma) > 1){
      candmethylfeatures <- row.names(sigg.methyllimma)
      
    }else if(nrow(sigg.methyllimma) == 1){
      
      candmethylfeatures <- row.names(methyldat)
      
      cat('Only one significant methylation diff feature can be called, will use 
          the original methylation features in the following analysis\n')
      
    }else{
      cat('No significant methylation diff features can be called, the analysis cannot continue\n')
      return(NULL)
    }
    
    
    if(nrow(sigg.rnalimma) > 0){
      candrnafeatures <- row.names(sigg.rnalimma)
      
      if(!is.null(topdiffrnafeaturenum)){
        candrnafeatures <- candrnafeatures[seq(1, 
                                               min(topdiffrnafeaturenum, 
                                                   length(candrnafeatures)), 1)]
      }
      
    }else{
      cat('No significant RNA diff features can be called, the analysis cannot continue\n')
      return(NULL)
    }
    
  }
  
  #mediation analysis
  
  elasticnetcoefs <- list()
  elasticnetrsquares <- list()
  annolist <- list()
  
  elasticnetweights <- list()
  
  if(balanceadj == TRUE){
    
    if(calldifffeatures == FALSE){
      candmethylfeatures <- row.names(methyldat)
      candrnafeatures <- row.names(rnadat)
    }
    
    qualifieses <- list()
    j <- 1
    for(j in 1:length(candrnafeatures)){
      
      cat(paste0('j = ', j, '\n'))
      
      candrnafeature <- candrnafeatures[j]
      mediatorname <- candrnafeature
      
      if(prescreen == TRUE){
        
        candmethylfeatures <- colnames(prescreenprobes(oripd = t(rnadat[candrnafeature, , 
                                                                        drop = FALSE]), 
                                                       orivar = t(methyldat[candmethylfeatures, , 
                                                                            drop = FALSE])))
        
      }
      
      
      probelist <- list()
      n <- 1
      for(n in 1:length(methylbalanceddats$dats)){
        
        cat(paste0('n = ', n, '\n'))
        
        candmethyldats <- methylbalanceddats$dats[[n]][candmethylfeatures, , drop = FALSE]
        candmethyldats <- t(candmethyldats)
        
        candrnadats <- rnabalanceddats$dats[[n]][candrnafeature, , drop = FALSE]
        candrnadats <- t(candrnadats)
        
        
        elasticres <- singlelassomediation(i = 1, 
                                           candrnadats = candrnadats, 
                                           oripddats = methylbalanceddats$pds[[n]], 
                                           candmethyldats = candmethyldats, 
                                           
                                           seednum = seed, 
                                           threads = threads, 
                                           alpha = 1, 
                                           errortype = errortype, 
                                           prescreen = FALSE, 
                                           crossrcut = 0, 
                                           plotting = FALSE, 
                                           
                                           responsevarname = responsevarname, 
                                           confoundings = confoundings, 
                                           responselevels = responselevels, 
                                           
                                           methylplatform = methylplatform, 
                                           onlyelasticnet = TRUE, 
                                           onlymediation = FALSE)
        
        gc()
        
        probelist[[n]] <- elasticres$elasticnetcoefs
        
      }
      
      
      #feature integration
      finalprobes <- probeselection(modelcoeflist = probelist)
      
      if(is.null(finalprobes)){
        next()
      }
      
      if(length(finalprobes) == 0){
        next()
      }
      
      
      
      #linear base learner training
      modellist <- list()
      p <- 1
      for(p in 1:length(methylbalanceddats$dats)){
        
        resampleresponse <- rnabalanceddats$pds[[p]]
        resampleresponse$Response <- rnabalanceddats$dats[[p]][candrnafeature,]
        
        balancedsamples <- list(resamplevars = t(methylbalanceddats$dats[[p]]), 
                                resampleresponse = resampleresponse)
        linearfit <- singlelinear(finalfeatures = names(finalprobes), 
                                  resamplevar = balancedsamples$resamplevars, 
                                  resampleresponse = balancedsamples$resampleresponse)
        
        modellist[[p]] <- linearfit
        
      }
      
      #final ensemble model
      resampleensemble <- weightlinears(innertrainvar = t(methyldat), 
                                        innertrainres = rnadat[candrnafeature,], 
                                        modellist = modellist, 
                                        rsquarecut = 0.5)
      if(is.null(resampleensemble)){
        next()
      }
      
      if(cor(resampleensemble$traincomp)[1, 2] < sqrt(crossrcut)){
        next()
      }
      
      
      modelscores <- as.matrix(finalprobes)
      colnames(modelscores)[1] <- 'score'
      if(length(resampleensemble$baselearners) >= 1){
        
        for(n in 1:length(resampleensemble$baselearners)){
          
          modelcoeff <- resampleensemble$baselearners[[n]]$coefficients
          names(modelcoeff) <- gsub(pattern = "`", replacement = '', 
                                    x = names(modelcoeff))
          modelcoeff <- modelcoeff[c('(Intercept)', row.names(modelscores))]
          modelcoeff <- as.matrix(modelcoeff)
          if(n == 1){
            modelcoeffs <- modelcoeff
          }else{
            modelcoeffs <- cbind(modelcoeffs, modelcoeff)
          }
          
        }
        
        colnames(modelcoeffs) <- paste0('model_', seq(1, ncol(modelcoeffs), 1))
      }
      
      
      modelcoeffs <- modelcoeffs[complete.cases(modelcoeffs), , drop = FALSE]
      modelcoeffs <- as.matrix(modelcoeffs)
      colnames(modelcoeffs) <- paste0('model_', seq(1, ncol(modelcoeffs), 1))
      modelscores <- modelscores[row.names(modelcoeffs)[row.names(modelcoeffs) != 
                                                          '(Intercept)'],]
      modelscores <- as.matrix(modelscores)
      colnames(modelscores)[1] <- 'score'
      
      
      mediationreses <- list()
      m <- 1
      for(m in 1:length(methylbalanceddats$dats)){
        
        cat(paste0('m = ', m, '\n'))
        
        candmethyldats <- methylbalanceddats$dats[[m]][names(finalprobes), , drop = FALSE]
        candmethyldats <- t(candmethyldats)
        
        candrnadats <- rnabalanceddats$dats[[m]][candrnafeature, , drop = FALSE]
        candrnadats <- t(candrnadats)
        
        
        mediationres <- singlelassomediation(i = 1, 
                                             candrnadats = candrnadats, 
                                             oripddats = methylbalanceddats$pds[[m]], 
                                             candmethyldats = candmethyldats, 
                                             
                                             seednum = seed, 
                                             threads = threads, 
                                             alpha = 1, 
                                             errortype = errortype, 
                                             prescreen = prescreen, 
                                             crossrcut = crossrcut, 
                                             plotting = FALSE, 
                                             
                                             responsevarname = responsevarname, 
                                             confoundings = confoundings, 
                                             responselevels = responselevels, 
                                             
                                             methylplatform = methylplatform, 
                                             onlyelasticnet = FALSE, 
                                             onlymediation = TRUE, 
                                             returncompleteres = TRUE, 
                                             propna = FALSE)
        
        gc()
        
        mediationreses[[names(methylbalanceddats$dats)[m]]] <- mediationres$completeres
        
        
      }
      
      
      mediationreses <- do.call(rbind, mediationreses)
      row.names(mediationreses) <- 1:nrow(mediationreses)
      
      completemediationres <- combineci(cires = mediationreses)
      
      completemediationres <- mediationanno(reslines = completemediationres, 
                                            propna = TRUE)
      
      qualifies <- subset(completemediationres, mediation == TRUE)
      qualifieses[[length(qualifieses) + 1]] <- qualifies
      
      if(!is.null(qualifies)){
        if(nrow(qualifies) > 0){
          row.names(qualifies) <- 1:nrow(qualifies)
          
          
          mediatorannores <-  mediatorannotation(mediators = qualifies$mediator, 
                                                 platform = methylplatform)
          
          targetannores <- mediatorannotation(mediators = qualifies$target, 
                                              platform = methylplatform)
          
          sourceannores <- mediatorannotation(mediators = qualifies$source, 
                                              platform = methylplatform)
          
          subannores <- unique(rbind(mediatorannores$geneannores, 
                                     targetannores$geneannores, 
                                     sourceannores$geneannores))
          subprobeannores <- unique(rbind(mediatorannores$probeannores, 
                                          targetannores$probeannores, 
                                          sourceannores$probeannores))
          
          finalannores <- subannores
          finalprobeannores <- subprobeannores
          
          elasticnetcoefs[[mediatorname]] <- data.frame(Probe = names(finalprobes), 
                                                        stringsAsFactors = FALSE)
          
          elasticnetcoefs[[mediatorname]] <- cbind(elasticnetcoefs[[mediatorname]], 
                                                   modelcoeffs[names(finalprobes),])
          
          elasticnetrsquares[[mediatorname]] <- cor(resampleensemble$traincomp)[1, 2]^2
          
          elasticnetweights[[mediatorname]] <- resampleensemble$baseweights
          
          
          if(!is.null(finalannores)){
            if(nrow(finalannores) > 0){
              row.names(finalannores) <- 1:nrow(finalannores)
            }
          }
          
          if(!is.null(finalprobeannores)){
            if(nrow(finalprobeannores) > 0){
              row.names(finalprobeannores) <- 1:nrow(finalprobeannores)
            }
          }
          
          if(!is.null(finalprobeannores) | !is.null(finalannores)){
            
            annolist[[mediatorname]] <- list()
            
            if(!is.null(finalprobeannores)){
              annolist[[mediatorname]]$probeannores <- finalprobeannores
            }
            
            if(nrow(finalannores) > 0){
              annolist[[mediatorname]]$geneannores <- finalannores
            }
            
          }
          
          cat(paste0('Mediator ', mediatorname, ' has been completed\n'))
          
        }
      }
      
    }
    
    
    qualifies <- do.call(rbind, qualifieses)
    qualifies <- unique(qualifies)
    
    if(!is.null(qualifies)){
      
      if(nrow(qualifies) > 0){
        row.names(qualifies) <- 1:nrow(qualifies)
      }else{
        return(NULL)
      }
      
    }else{
      return(NULL)
    }
    
    
    
    qualifies <- filtermediator(mediationres = qualifies)
    
    if(is.null(qualifies)){
      return(NULL)
    }
    
    for(i in 1:length(annolist)){
      
      if(names(annolist)[i] %in% qualifies$mediator){
        
        subqualifies <- subset(qualifies, mediator == names(annolist)[i])
        row.names(subqualifies) <- 1:nrow(subqualifies)
        
        mediatorannores <-  mediatorannotation(mediators = subqualifies$mediator, 
                                               platform = methylplatform)
        
        targetannores <- mediatorannotation(mediators = subqualifies$target, 
                                            platform = methylplatform)
        
        sourceannores <- mediatorannotation(mediators = subqualifies$source, 
                                            platform = methylplatform)
        
        subannores <- unique(rbind(mediatorannores$geneannores, 
                                   targetannores$geneannores, 
                                   sourceannores$geneannores))
        subprobeannores <- unique(rbind(mediatorannores$probeannores, 
                                        targetannores$probeannores, 
                                        sourceannores$probeannores))
        
        finalannores <- subannores
        finalprobeannores <- subprobeannores
        
        
        if(!is.null(finalannores)){
          if(nrow(finalannores) > 0){
            row.names(finalannores) <- 1:nrow(finalannores)
          }
        }
        
        if(!is.null(finalprobeannores)){
          if(nrow(finalprobeannores) > 0){
            row.names(finalprobeannores) <- 1:nrow(finalprobeannores)
          }
        }
        
        if(!is.null(finalprobeannores)){
          annolist[[i]]$probeannores <- finalprobeannores
        }
        
        if(nrow(finalannores) > 0){
          annolist[[i]]$geneannores <- finalannores
        }
        
      }
      
    }
    
    for(annolistname in names(annolist)){
      
      if(!(annolistname %in% qualifies$mediator)){
        
        annolist[[annolistname]] <- NULL
        
        elasticnetcoefs[[annolistname]] <- NULL
        
        elasticnetrsquares[[annolistname]] <- NULL
        
        elasticnetweights[[annolistname]] <- NULL
        
      }
      
    }
    
    
    res <- list(mediationres = qualifies, 
                elasticnetcoefs = elasticnetcoefs, 
                elasticnetrsquares = elasticnetrsquares, 
                elasticnetweights = elasticnetweights, 
                annolist = annolist)
    
    if(calldifffeatures == TRUE){
      res$sigg.rnalimma <- sigg.rnalimma
      res$sigg.methyllimma <- sigg.methyllimma
      res$nonsigg.rnalimma <- nonsigg.rnalimma
      res$nonsigg.methyllimma <- nonsigg.methyllimma
      
    }
    
    if(plot == TRUE){
      
      #res <- readRDS('luad.multimediationres.balance.rds')
      
      unimediators <- unique(res$mediationres$mediator)
      
      i <- 1
      for(i in 1:length(unimediators)){
        
        unimediator <- unimediators[i]
        
        sub <- subset(res$mediationres, mediator == unimediator)
        
        selectedfeatures <- setdiff(c(sub$source, sub$target), responsevarname)
        
        annotab <- data.frame(Dir = rep('Fwd', length(selectedfeatures)), 
                              stringsAsFactors = FALSE)
        row.names(annotab) <- selectedfeatures
        annotab$Dir[row.names(annotab) %in% res$mediationres$target] <- 'Rev'
        annotab$Dir <- factor(annotab$Dir, levels = c('Fwd', 'Rev'), 
                              ordered = TRUE)
        
        
        unimediatordat <- rnadat[unimediator,]
        
        responsedat <- pddat[responsevarname][,1]
        
        comptab <- data.frame(mediator = unimediatordat, 
                              response = responsedat, 
                              stringsAsFactors = FALSE)
        colnames(comptab) <- c(unimediator, responsevarname)
        
        mediationheatmap(orivar = t(methyldat), 
                         comptab = comptab, 
                         annotab = annotab, 
                         featurenames = selectedfeatures, 
                         title = unimediator, 
                         textsize = textsize, 
                         show_rownames = TRUE, 
                         quantilecuts = NULL)
        
      }
      
    }
    
    return(res)
    
  }else{
    
    if(calldifffeatures == TRUE){
      
      candmethyldats <- t(methyldat[candmethylfeatures, , drop = FALSE])
      candrnadats <- t(rnadat[candrnafeatures, , drop = FALSE])
      
    }else{
      
      candmethyldats <- t(methyldat)
      candrnadats <- t(rnadat)
      
    }
    
    jseqs <- seq(1, ncol(candrnadats), 1)
    
    qualifieses <- list()
    j <- 1
    for(j in jseqs){
      
      cat(paste0('j = ', j, '\n'))
      
      mediationres <- singlelassomediation(i = j, 
                                           candrnadats = candrnadats, 
                                           oripddats = pddat, 
                                           candmethyldats = candmethyldats, 
                                           
                                           seednum = seed, 
                                           threads = threads, 
                                           alpha = 1, 
                                           errortype = errortype, 
                                           prescreen = prescreen, 
                                           crossrcut = crossrcut, 
                                           plotting = FALSE, 
                                           
                                           responsevarname = responsevarname, 
                                           confoundings = confoundings, 
                                           responselevels = responselevels, 
                                           
                                           methylplatform = methylplatform, 
                                           propna = TRUE)
      
      gc()
      
      if(!is.null(mediationres)){
        
        mediatorname <- colnames(candrnadats)[j]
        
        qualifieses[[length(qualifieses) + 1]] <- mediationres$mediationres
        
        mediationres$mediationres <- NULL
        
        elasticnetcoefs[[mediatorname]] <- mediationres$elasticnetcoefs
        
        mediationres$elasticnetcoefs <- NULL
        
        elasticnetrsquares[[mediatorname]] <- mediationres$rsquare
        
        mediationres$rsquare <- NULL
        
        annolist[[mediatorname]] <- mediationres
        
        cat(paste0('Mediator ', mediatorname, ' has been completed\n'))
        
      }
      
    }
    
    qualifies <- do.call(rbind, qualifieses)
    qualifies <- unique(qualifies)
    
    if(!is.null(qualifies)){
      
      if(nrow(qualifies) > 0){
        row.names(qualifies) <- 1:nrow(qualifies)
      }else{
        return(NULL)
      }
      
    }else{
      return(NULL)
    }
    
    
    
    qualifies <- filtermediator(mediationres = qualifies)
    
    if(is.null(qualifies)){
      return(NULL)
    }
    
    for(i in 1:length(annolist)){
      
      if(names(annolist)[i] %in% qualifies$mediator){
        
        subqualifies <- subset(qualifies, mediator == names(annolist)[i])
        row.names(subqualifies) <- 1:nrow(subqualifies)
        
        mediatorannores <-  mediatorannotation(mediators = subqualifies$mediator, 
                                               platform = methylplatform)
        
        targetannores <- mediatorannotation(mediators = subqualifies$target, 
                                            platform = methylplatform)
        
        sourceannores <- mediatorannotation(mediators = subqualifies$source, 
                                            platform = methylplatform)
        
        subannores <- unique(rbind(mediatorannores$geneannores, 
                                   targetannores$geneannores, 
                                   sourceannores$geneannores))
        subprobeannores <- unique(rbind(mediatorannores$probeannores, 
                                        targetannores$probeannores, 
                                        sourceannores$probeannores))
        
        finalannores <- subannores
        finalprobeannores <- subprobeannores
        
        
        if(!is.null(finalannores)){
          if(nrow(finalannores) > 0){
            row.names(finalannores) <- 1:nrow(finalannores)
          }
        }
        
        if(!is.null(finalprobeannores)){
          if(nrow(finalprobeannores) > 0){
            row.names(finalprobeannores) <- 1:nrow(finalprobeannores)
          }
        }
        
        if(!is.null(finalprobeannores)){
          annolist[[i]]$probeannores <- finalprobeannores
        }
        
        if(nrow(finalannores) > 0){
          annolist[[i]]$geneannores <- finalannores
        }
        
      }
      
    }
    
    for(annolistname in names(annolist)){
      
      if(!(annolistname %in% qualifies$mediator)){
        
        annolist[[annolistname]] <- NULL
        
        elasticnetcoefs[[annolistname]] <- NULL
        
        elasticnetrsquares[[annolistname]] <- NULL
        
        #elasticnetweights[[annolistname]] <- NULL
        
      }
      
    }
    
    
    
    res <- list(mediationres = qualifies, 
                elasticnetcoefs = elasticnetcoefs, 
                elasticnetrsquares = elasticnetrsquares, 
                annolist = annolist)
    
    if(calldifffeatures == TRUE){
      res$sigg.rnalimma <- sigg.rnalimma
      res$sigg.methyllimma <- sigg.methyllimma
      res$nonsigg.rnalimma <- nonsigg.rnalimma
      res$nonsigg.methyllimma <- nonsigg.methyllimma
      
    }
    
    if(plot == TRUE){
      
      #res <- readRDS('luad.multimediationres.balance.rds')
      
      unimediators <- unique(res$mediationres$mediator)
      
      i <- 1
      for(i in 1:length(unimediators)){
        
        unimediator <- unimediators[i]
        
        sub <- subset(res$mediationres, mediator == unimediator)
        
        selectedfeatures <- setdiff(c(sub$source, sub$target), responsevarname)
        
        annotab <- data.frame(Dir = rep('Fwd', length(selectedfeatures)), 
                              stringsAsFactors = FALSE)
        row.names(annotab) <- selectedfeatures
        annotab$Dir[row.names(annotab) %in% res$mediationres$target] <- 'Rev'
        annotab$Dir <- factor(annotab$Dir, levels = c('Fwd', 'Rev'), 
                              ordered = TRUE)
        
        
        unimediatordat <- rnadat[unimediator,]
        
        responsedat <- pddat[responsevarname][,1]
        
        comptab <- data.frame(mediator = unimediatordat, 
                              response = responsedat, 
                              stringsAsFactors = FALSE)
        colnames(comptab) <- c(unimediator, responsevarname)
        
        mediationheatmap(orivar = t(methyldat), 
                         comptab = comptab, 
                         annotab = annotab, 
                         featurenames = selectedfeatures, 
                         title = unimediator, 
                         textsize = textsize, 
                         show_rownames = TRUE, 
                         quantilecuts = NULL)
        
      }
      
    }
    
    return(res)
    
  }
  
  
}






