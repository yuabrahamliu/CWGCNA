#Enrichment functions####

getcorrgenes <- function(probes,
                         pairedRNA,
                         pairedmethyl,
                         cutoff = 0.5,
                         generegions = c('TSS200'),
                         platform = 450, 
                         maptogenes = TRUE){
  
  probes <- intersect(probes, row.names(pairedmethyl))
  pairedprobes <- pairedmethyl[probes, , drop = FALSE]
  
  if(maptogenes == TRUE){
    
    pairedprobes <- probestogenes(betadat = pairedprobes,
                                  platform = platform,
                                  group450k850k = generegions,
                                  includemultimatch = FALSE)
    
  }
  
  pairedprobes <- t(pairedprobes)
  pairedgenes <- t(pairedRNA)
  
  posgenelist <- c()
  neggenelist <- c()
  
  genecors <- t(cor(pairedprobes, pairedgenes))
  
  posgenelist <- c(posgenelist,
                   row.names(which(x = genecors > cutoff, arr.ind = TRUE)))
  neggenelist <- c(neggenelist,
                   row.names(which(x = genecors < -cutoff, arr.ind = TRUE)))
  posgenelist <- unique(posgenelist)
  neggenelist <- unique(neggenelist)
  
  sharedgeneslist <- intersect(posgenelist, neggenelist)
  posgenelist <- setdiff(posgenelist, sharedgeneslist)
  neggenelist <- setdiff(neggenelist, sharedgeneslist)
  
  res <- list(posgenes = posgenelist,
              neggenes = neggenelist)
  
  return(res)
  
  
}

getrnagenes <- function(upprobes = NULL, 
                        dnprobes = NULL, 
                        pairedRNA,
                        pairedmethyl,
                        abscut = 0.5,
                        generegions = c('TSS200'),
                        platform = 450, 
                        maptogenes = TRUE){
  
  if(length(upprobes) > 0){
    
    upprobescorgenes <- getcorrgenes(probes = upprobes,
                                     pairedRNA = pairedRNA,
                                     pairedmethyl = pairedmethyl,
                                     cutoff = abscut,
                                     generegions = generegions,
                                     platform = platform, 
                                     maptogenes = maptogenes)
    
    upprobesposgenes <- upprobescorgenes$posgenes
    upprobesneggenes <- upprobescorgenes$neggenes
    
  }else{
    
    upprobescorgenes <- NULL
    
    upprobesposgenes <- NULL
    upprobesneggenes <- NULL
    
  }
  
  
  if(length(dnprobes) > 0){
    
    dnprobescorgenes <- getcorrgenes(probes = dnprobes,
                                     pairedRNA = pairedRNA,
                                     pairedmethyl = pairedmethyl,
                                     cutoff = abscut,
                                     generegions = generegions,
                                     platform = platform, 
                                     maptogenes = maptogenes)
    
    dnprobesposgenes <- dnprobescorgenes$posgenes
    dnprobesneggenes <- dnprobescorgenes$neggenes
    
  }else{
    
    dnprobescorgenes <- NULL
    
    dnprobesposgenes <- NULL
    dnprobesneggenes <- NULL
    
  }
  
  enhancegenes <- unique(c(upprobesposgenes, dnprobesneggenes))
  inhibitgenes <- unique(c(upprobesneggenes, dnprobesposgenes))
  
  sharedgenes <- intersect(enhancegenes, inhibitgenes)
  enhancegenes <- setdiff(enhancegenes, sharedgenes)
  inhibitgenes <- setdiff(inhibitgenes, sharedgenes)
  
  res <- list(upprobescorgenes = upprobescorgenes,
              dnprobescorgenes = dnprobescorgenes,
              enhancegenes = enhancegenes,
              inhibitgenes = inhibitgenes)
  
  return(res)
  
}

getprobegenes <- function(upprobes = NULL, 
                          dnprobes = NULL, 
                          platform = 450,
                          generegions = c('TSS200')){
  
  if(length(upprobes) > 0){
    
    upgenes <- probeanno(platform = platform,
                         probes = upprobes)
    
    if(!is.null(upgenes)){
      upgenes <- subset(upgenes, UCSC_RefGene_Group %in% generegions)
      
      upgenes <- unique(upgenes$UCSC_RefGene_Name)
    }else{
      upgenes <- NULL
    }
    
  }else{
    
    upgenes <- NULL
    
  }
  
  
  if(length(dnprobes) > 0){
    
    dngenes <- probeanno(platform = platform,
                         probes = dnprobes)
    
    if(!is.null(dngenes) > 0){
      dngenes <- subset(dngenes, UCSC_RefGene_Group %in% generegions)
      
      dngenes <- unique(dngenes$UCSC_RefGene_Name)
    }else{
      dngenes <- NULL
    }
    
  }else {
    
    dngenes <- NULL
    
  }
  
  res <- list(upgenes = upgenes, dngenes = dngenes)
  
  return(res)
  
}

getenrichrres <- function(genes,
                          write = FALSE,
                          prefix,
                          dbs = c('GO_Biological_Process_2018',
                                  'GO_Molecular_Function_2018',
                                  'BioPlanet_2019', 'WikiPathways_2019_Human',
                                  'KEGG_2019_Human', 'BioCarta_2016',
                                  'Reactome_2016', 'NCI-Nature_2016',
                                  'Panther_2016')){
  genes <- genes[!is.na(genes)]
  genes <- unique(genes)
  
  if(length(genes) == 0){
    return(NULL)
  }
  
  if(sum(genes != '') == 0){
    return(NULL)
  }
  
  library(enrichR)
  
  enriched <- enrichr(genes = genes, databases = dbs)
  
  enrichedsum <- do.call(rbind, enriched)
  
  if(nrow(enrichedsum) == 0){
    return(NULL)
  }
  
  enrichedsum$Database <- row.names(enrichedsum)
  enrichedsum$Database <- gsub(pattern = '\\..*$', replacement = '',
                               x = enrichedsum$Database)
  enrichedsum$Term <- gsub(pattern = '_.*$', replacement = '',
                           x = enrichedsum$Term)
  
  enrichedsum <- enrichedsum[c('Term', 'Database', 'Adjusted.P.value', 'P.value',
                               'Odds.Ratio', 'Combined.Score',
                               'Overlap', 'Genes')]
  
  if(sum(enrichedsum$Adjusted.P.value < 0.05) >= 10){
    enrichedsum <- subset(enrichedsum, Adjusted.P.value < 0.05)
  }else{
    enrichedsum <- subset(enrichedsum, P.value < 0.01)
    if(nrow(enrichedsum) == 0){
      return(NULL)
    }
  }
  
  enrichedsum <- enrichedsum[order(enrichedsum$Adjusted.P.value,
                                   enrichedsum$P.value,
                                   -enrichedsum$Combined.Score,
                                   -enrichedsum$Odds.Ratio),]
  row.names(enrichedsum) <- 1:nrow(enrichedsum)
  enrichedsum$Annotation <- prefix
  
  
  if(write == TRUE){
    
    stamp <- Sys.time()
    stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
    stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
    stamp <- paste0('.', stamp)
    
    write.table(enrichedsum,
                paste0(prefix, '_enrichr', stamp, '.txt'),
                sep = '\t',
                row.names = FALSE,
                quote = FALSE)
    
  }
  
  return(enrichedsum)
  
}

Enrich <- function(enhancegenes,
                   inhibitgenes,
                   pairedRNA = NULL,
                   pairedmethyl = NULL,
                   abscut = 0.5,
                   dbs = c('GO_Biological_Process_2018',
                           'BioPlanet_2019',
                           'Reactome_2016'),
                   write = FALSE){
  
  
  #library(eClock)
  
  if(length(inhibitgenes) > 0){
    upanno <- geneanno(genesymbols = inhibitgenes)
  }else{
    upanno <- NULL
  }
  
  if(length(enhancegenes) > 0){
    dnanno <- geneanno(genesymbols = enhancegenes)
  }else{
    dnanno <- NULL
  }
  
  
  
  if(is.null(pairedRNA) | is.null(pairedmethyl)){
    
    upenrich <- getenrichrres(genes = inhibitgenes,
                              write = write,
                              prefix = paste0('hyper'),
                              dbs = dbs)
    
    dnenrich <- getenrichrres(genes = enhancegenes,
                              write = write,
                              prefix = paste0('hypo'),
                              dbs = dbs)
    
    
    res <- list(hypergeneanno = upanno,
                hypogeneanno = dnanno,
                hypergeneenrich = upenrich,
                hypogeneenrich = dnenrich)
    
  }else{
    
    upenrich <- getenrichrres(genes = inhibitgenes,
                              write = write,
                              prefix = paste0('inhibit'),
                              dbs = dbs)
    
    dnenrich <- getenrichrres(genes = enhancegenes,
                              write = write,
                              prefix = paste0('enhance'),
                              dbs = dbs)
    
    
    res <- list(enhancegeneanno = dnanno,
                inhibitgeneanno = upanno,
                enhancegeneenrich = dnenrich,
                inhibitgeneenrich = upenrich)
    
  }
  
  
  return(res)
  
}

getsinglename <- function(singlegenename, entrez = TRUE){
  
  genenames <- c(sub(pattern = '::.*$', replacement = '', singlegenename))
  genenames <- c(genenames, c(sub(pattern = '^.*::', replacement = '', singlegenename)))
  genenames <- unique(genenames)
  
  geneentrezs <- genesymbols <- genenames
  
  if(entrez == TRUE){
    suppressWarnings(genenameres <- geneentrezs[!is.na(as.numeric(geneentrezs))])
    
    if(length(genenameres) == 0){
      
      library(org.Hs.eg.db)
      
      suppressMessages(
        
        suppressWarnings(
          
          genenameres <- tryCatch({
            AnnotationDbi::select(x = org.Hs.eg.db, 
                                  keys = geneentrezs[is.na(as.numeric(geneentrezs))], 
                                  columns = 'ENTREZID', 
                                  keytype = 'SYMBOL')
          }, error = function(err){
            NA
          })
          
        )
        
      )
      
      suppressWarnings(
        
        if(!is.na(genenameres)){
          genenameres <- unique(genenameres$ENTREZID)
          if(length(genenameres) > 1){
            genenameres <- NA
          }
        }
        
      )
      
      
    }
    
  }else{
    
    suppressWarnings(genenameres <- genesymbols[is.na(as.numeric(geneentrezs))])
    
    if(length(genenameres) == 0){
      
      library(org.Hs.eg.db)
      
      suppressMessages(
        
        suppressWarnings(
          
          genenameres <- tryCatch({
            AnnotationDbi::select(x = org.Hs.eg.db, 
                                  keys = genesymbols[!is.na(as.numeric(geneentrezs))], 
                                  columns = 'SYMBOL', 
                                  keytype = 'ENTREZID')
          }, error = function(err){
            NA
          })
          
        )
        
      )
      
      suppressWarnings(
        
        if(!is.na(genenameres)){
          genenameres <- unique(genenameres$SYMBOL)
          if(length(genenameres) > 1){
            genenameres <- NA
          }
        }
        
      )
      
    }
    
  }
  
  if(length(genenameres) == 0){
    genenameres <- NA
  }
  
  return(genenameres)
  
}

#'Perform correlation-based gene functional enrichment or EnrichR enrichment
#'
#'Perform correlation-based gene functional enrichment or EnrichR enrichment. 
#'
#'@param hyperprobes The vector with the hypermethylated methylation probes or 
#'  DOWN-regulated RNA genes to be analyzed. Default is NULL.
#'@param hypoprobes The vector with the hypomethylated methylation probes or 
#'  UP-regulated RNA genes to be analyzed. Default value is NULL. Can transfer 
#'  2 vectors to the 2 parameters \code{hyperprobes} and \code{hypoprobes}, or 
#'  transfer 1 vector to 1 of these 2 parameters and leave the other as NULL. 
#'@param pairedRNA The RNA part of the paired RNA-DNA methylation dataset. If 
#'  both this parameter and \code{pairedmethyl} are not NULL, a correlation-
#'  based method will be used to find the genes whose RNA expression level in 
#'  the RNA data significantly correlated with any of the \code{hyperprobes} 
#'  and \code{hypoprobes} in the paired methylation data. The probe values in 
#'  the methylation data will be converted to gene values first, and then the 
#'  correlation analysis will be performed. In addition to conducting EnrichR 
#'  enrichment for the correlated genes finally identified, each of the genes 
#'  will also get its own function annotated. Should transfer a matrix to this 
#'  parameter, with columns representing samples and rows representing genes. 
#'  Row names are gene symbols and column names are sample IDs. Default value 
#'  is NULL, and in this case, the correlation-based method will not be used, 
#'  and instead, function annotation and enrichment analysis will be directly 
#'  conducted on the hyper and hypomethylated methylation probes transferred 
#'  to \code{hyperprobes} and \code{hypoprobes}, with these probes mapped to 
#'  corresponding genes first. Or, if \code{hyperprobes} and \code{hypoprobes} 
#'  have already been gene symbols, the mapping step will be skipped.
#'@param pairedmethyl The methylation part of the paired RNA-DNA methylation 
#'  dataset. Similar to \code{pairedRNA}, if both these 2 parameters are not 
#'  NULL, the correlation-based method will be used, but if any one of them is 
#'  NULL, EnrichR enrichment will be directly conducted on probes or genes in 
#'  \code{hyperprobes} and \code{hypoprobes}. \code{pairedmethyl} should be a 
#'  beta value matrix with columns representing samples and rows for PROBEs. 
#'  Row names are PROBE names. Column names are sample names. The sample names 
#'  should be the same as the ones in \code{pairedRNA}, because they are data 
#'  for paired samples. Default value of \code{pairedmethyl} is NULL.
#'@param abscut When the correlation method is used to select the genes with a 
#'  significant correlation to the probes, the ones with an absolute Pearson 
#'  correlation coefficient value greater than this parameter value will be 
#'  selected and used for function analysis. Default is 0.7.
#'@param generegions When \code{hyperprobes} and \code{hypoprobes} contain the 
#'  names of methylation probes, this parameter will be used when converting 
#'  them into genes. That is, if a probe located in the regions defined by it, 
#'  the corresponding gene will be used for the function analysis. Default is 
#'  "TSS200", so probes within TSS200 regions will be used for gene mapping. 
#'  It can also be a vector, such as \code{c("TSS200", "TSS1500", "1stExon")}, 
#'  so that probes within any of these 3 regions will be used for the mapping. 
#'  This gene mapping process will be used when any of \code{pairedRNA} and 
#'  \code{pairedmethyl} is NULL, so the EnrichR enrichment will be performed 
#'  directly on the genes mapped from the parameters \code{hyperprobes} and 
#'  \code{hypoprobes}. When both \code{pairedRNA} and \code{pairedmethyl} are 
#'  not NULL, the correlation-based analysis will be used, and this mapping 
#'  will also be used because before screening for the correlated genes, the 
#'  the methylation probes will be merged into gene methylation values first, 
#'  and then used for calculating the correlation with the gene RNA values in 
#'  the paired RNA data. The gene methylation values are summarized according 
#'  to the regions in this parameter \code{generegions}.
#'@param platform The platform of the probes. Can be 27 (for 27k platform), 
#'  450 (for 450k platform), or 850 (for EPIC platform). Default is 450. 
#'@param dbs The databases for gene functional enrichment analysis. Default is 
#'  \code{c("GO_Biological_Process_2018", "BioPlanet_2019", "Reactome_2016")}.
#'@param write A logical value indicating whether the gene function results
#'  need to be written into txt files in the working directory. The default
#'  value is FALSE.
#'@return A list with four slots recording the gene annotation and functional 
#'  enrichment results for the enhanced and inhibited genes (when correlation- 
#'  based method is used) or for the genes mapped from the hypomethylated and 
#'  hypermethylated probes (when correlation-based method is not used). When 
#'  \code{hyperprobes} and \code{hypoprobes} contain gene names rather than 
#'  probes, the slot \code{hypergeneenrich} records the result of the genes in 
#'  \code{hyperprobes}, which are DOWN-regulated genes. In constrast, the slot 
#'  \code{hypogeneenrich} recordes the results for \code{hypoprobes}, which 
#'  UP-regulated genes. If \code{write} is TRUE, txt files will be generated 
#'  to save these results.
#'@export
corenrich <- function(hyperprobes = NULL, 
                      hypoprobes = NULL, 
                      pairedRNA = NULL,
                      pairedmethyl = NULL,
                      abscut = 0.7,
                      generegions = c('TSS200'),
                      platform = 450,
                      dbs = c('GO_Biological_Process_2018',
                              'BioPlanet_2019',
                              'Reactome_2016'),
                      write = FALSE){
  
  hyperprobes <- unique(hyperprobes)
  hypoprobes <- unique(hypoprobes)
  
  if(is.null(pairedRNA) | is.null(pairedmethyl)){
    
    siggenes <- getprobegenes(upprobes = hyperprobes, 
                              dnprobes = hypoprobes, 
                              platform = platform, 
                              generegions = generegions)
    
    enhancegenes <- siggenes$dngenes
    inhibitgenes <- siggenes$upgenes
    
    if(is.null(enhancegenes) & is.null(inhibitgenes)){
      enhancegenes <- getsinglename(singlegenename = hypoprobes, entrez = FALSE)
      inhibitgenes <- getsinglename(singlegenename = hyperprobes, entrez = FALSE)
    }
    
    
  }else{
    
    siggenes <- getrnagenes(upprobes = hyperprobes, 
                            dnprobes = hypoprobes, 
                            pairedRNA = pairedRNA,
                            pairedmethyl = pairedmethyl,
                            abscut = abscut, 
                            maptogenes = TRUE)
    enhancegenes <- siggenes$enhancegenes
    inhibitgenes <- siggenes$inhibitgenes
    
  }
  
  sharedgenes <- intersect(enhancegenes, inhibitgenes)
  enhancegenes <- setdiff(enhancegenes, sharedgenes)
  inhibitgenes <- setdiff(inhibitgenes, sharedgenes)
  
  
  
  res <- Enrich(enhancegenes = enhancegenes,
                inhibitgenes = inhibitgenes,
                pairedRNA = pairedRNA,
                pairedmethyl = pairedmethyl,
                abscut = abscut,
                dbs = dbs,
                write = write)
  
  return(res)
  
}

#TopoGO####

getmultinamedat <- function(multigenenames, entrez = TRUE){
  
  mappings <- NULL
  
  genenames <- c(sub(pattern = '::.*$', replacement = '', multigenenames))
  genenames <- c(genenames, c(sub(pattern = '^.*::', replacement = '', multigenenames)))
  genenames <- unique(genenames)
  
  geneentrezs <- genesymbols <- genenames
  
  if(entrez == TRUE){
    
    suppressWarnings(genenameres <- geneentrezs[!is.na(as.numeric(geneentrezs))])
    
    if(length(genenameres) == 0){
      
      library(org.Hs.eg.db)
      
      suppressMessages(
        
        suppressWarnings(
          
          genenameres <- tryCatch({
            AnnotationDbi::select(x = org.Hs.eg.db, 
                                  keys = geneentrezs[is.na(as.numeric(geneentrezs))], 
                                  columns = 'ENTREZID', 
                                  keytype = 'SYMBOL')
          }, error = function(err){
            NULL
          })
          
        )
        
      )
      
      suppressWarnings(
        
        if(!is.null(genenameres)){
          
          if(nrow(genenameres) == 0){
            
            mappings <- NULL
            
          }else{
            
            genenameres <- genenameres[complete.cases(genenameres), , drop = FALSE]
            
            if(nrow(genenameres) == 0){
              
              mappings <- NULL
              
            }else{
              
              genenameres <- unique(genenameres)
              
              dupsymbols <- unique(genenameres$SYMBOL[duplicated(genenameres$SYMBOL)])
              dupentrezs <- unique(genenameres$ENTREZID[duplicated(genenameres$ENTREZID)])
              
              genenameres <- subset(genenameres, !(SYMBOL %in% dupsymbols) & 
                                      !(ENTREZID %in% dupentrezs))
              
              if(nrow(genenameres) == 0){
                
                mappings <- NULL
                
              }else{
                
                mappings <- genenameres
                mappings$genenames <- paste0(genenameres$SYMBOL, '::', 
                                             genenameres$ENTREZID)
                
                mappings <- mappings[c('SYMBOL', 'ENTREZID', 'genenames')]
                mappings <- unique(mappings)
                row.names(mappings) <- 1:nrow(mappings)
                
                
              }
              
            }
            
          }
          
        }else{
          mappings <- NULL
        }
        
      )
      
      
    }
    
  }else{
    
    suppressWarnings(genenameres <- genesymbols[is.na(as.numeric(geneentrezs))])
    
    if(length(genenameres) == 0){
      
      library(org.Hs.eg.db)
      
      suppressMessages(
        
        suppressWarnings(
          
          genenameres <- tryCatch({
            AnnotationDbi::select(x = org.Hs.eg.db, 
                                  keys = genesymbols[!is.na(as.numeric(geneentrezs))], 
                                  columns = 'SYMBOL', 
                                  keytype = 'ENTREZID')
          }, error = function(err){
            NULL
          })
          
        )
        
      )
      
      suppressWarnings(
        
        if(!is.null(genenameres)){
          
          if(nrow(genenameres) == 0){
            
            mappings <- NULL
            
          }else{
            
            genenameres <- genenameres[complete.cases(genenameres), , drop = FALSE]
            
            if(nrow(genenameres) == 0){
              
              mappings <- NULL
              
            }else{
              
              genenameres <- unique(genenameres)
              
              dupsymbols <- unique(genenameres$SYMBOL[duplicated(genenameres$SYMBOL)])
              dupentrezs <- unique(genenameres$ENTREZID[duplicated(genenameres$ENTREZID)])
              
              genenameres <- subset(genenameres, !(SYMBOL %in% dupsymbols) & 
                                      !(ENTREZID %in% dupentrezs))
              
              if(nrow(genenameres) == 0){
                
                mappings <- NULL
                
              }else{
                
                mappings <- genenameres
                mappings$genenames <- paste0(genenameres$SYMBOL, '::', 
                                             genenameres$ENTREZID)
                
                mappings <- mappings[c('SYMBOL', 'ENTREZID', 'genenames')]
                mappings <- unique(mappings)
                row.names(mappings) <- 1:nrow(mappings)
                
                
              }
              
            }
            
            
          }
          
          
        }else{
          mappings <- NULL
        }
        
      )
      
    }
    
  }
  
  return(mappings)
  
}

#make whole genome node background

makemappingtable <- function(mappinglist){
  
  mappingtable <- data.frame(GOID = rep(names(mappinglist), 
                                        as.vector(lapply(X = mappinglist, 
                                                         FUN = length))), 
                             ENTREZID = as.vector(unlist(mappinglist)), 
                             stringsAsFactors = FALSE)
  mappingtable <- unique(mappingtable)
  row.names(mappingtable) <- 1:nrow(mappingtable)
  
  return(mappingtable)
  
}

makenodebackground <- function(species = 'human'){
  
  #keggnames <- readRDS('C:/Users/yuabr/Desktop/Transfer/codetransfer/emethyl/files/keggnames.rds')
  keggnames <- get('keggnames')
  
  if(species == 'human'){
    
    library(org.Hs.eg.db)
    
    gg <- as.list(org.Hs.egGO2EG)
    gg <- makemappingtable(mappinglist = gg)
    
    kk <- as.list(org.Hs.egPATH2EG)
    kk <- makemappingtable(mappinglist = kk)
    names(kk)[1] <- 'KEGGID'
    
    library(reactome.db)
    rr <- as.list(reactomePATHID2EXTID)
    rr <- makemappingtable(mappinglist = rr)
    names(rr)[1] <- 'REACTOMEID'
    rr <- rr[grepl(pattern = '^R\\-HSA\\-[0-9]*$', x = rr$REACTOMEID),]
    rr <- unique(rr)
    row.names(rr) <- 1:nrow(rr)
    
    gg <- gg[c('ENTREZID', 'GOID')]
    rr <- rr[c('ENTREZID', 'REACTOMEID')]
    kk <- kk[c('ENTREZID', 'KEGGID')]
    
    res <- list(gos = gg, reactomes = rr, keggs = kk)
    
    
  }else if(species == 'mouse'){
    
    library(org.Mm.eg.db)
    mgg <- as.list(org.Mm.egGO2EG)
    mgg <- makemappingtable(mappinglist = mgg)
    
    mkk <- as.list(org.Mm.egPATH2EG)
    mkk <- makemappingtable(mappinglist = mkk)
    names(mkk)[1] <- 'KEGGID'
    
    library(reactome.db)
    mrr <- as.list(reactomePATHID2EXTID)
    mrr <- makemappingtable(mappinglist = mrr)
    names(mrr)[1] <- 'REACTOMEID'
    mrr <- mrr[grepl(pattern = '^R\\-MMU\\-[0-9]*$', x = mrr$REACTOMEID),]
    mrr <- unique(mrr)
    row.names(mrr) <- 1:nrow(mrr)
    
    mgg <- mgg[c('ENTREZID', 'GOID')]
    mrr <- mrr[c('ENTREZID', 'REACTOMEID')]
    mkk <- mkk[c('ENTREZID', 'KEGGID')]
    
    res <- list(gos = mgg, reactomes = mrr, keggs = mkk)
    
  }
  
  gos <- res$gos
  library(GO.db)
  suppressMessages(gonames <- AnnotationDbi::select(x = GO.db, 
                                                    keys = unique(gos$GOID), 
                                                    columns = 'TERM', 
                                                    keytype = 'GOID'))
  suppressMessages(goonts <- AnnotationDbi::select(x = GO.db, 
                                                   keys = unique(gos$GOID), 
                                                   columns = 'ONTOLOGY', 
                                                   keytype = 'GOID'))
  gonames <- subset(gonames, GOID %in% gos$GOID)
  goonts <- subset(goonts, GOID %in% gos$GOID)
  
  gos <- merge(x = gos, y = gonames, by = 'GOID', all = TRUE)
  gos <- merge(x = gos, y = goonts, by = 'GOID', all = TRUE)
  gos <- gos[c('ENTREZID', 'GOID', 'TERM', 'ONTOLOGY')]
  row.names(gos) <- 1:nrow(gos)
  
  gobps <- subset(gos, ONTOLOGY == 'BP')
  goccs <- subset(gos, ONTOLOGY == 'CC')
  gomfs <- subset(gos, ONTOLOGY == 'MF')
  
  gobps <- gobps[c('ENTREZID', 'GOID', 'TERM')]
  goccs <- goccs[c('ENTREZID', 'GOID', 'TERM')]
  gomfs <- gomfs[c('ENTREZID', 'GOID', 'TERM')]
  names(gobps)[2] <- 'BPID'
  names(goccs)[2] <- 'CCID'
  names(gomfs)[2] <- 'MFID'
  
  row.names(gobps) <- 1:nrow(gobps)
  row.names(goccs) <- 1:nrow(goccs)
  row.names(gomfs) <- 1:nrow(gomfs)
  
  res$gos <- gos
  res$gobps <- gobps
  res$goccs <- goccs
  res$gomfs <- gomfs
  
  
  keggs <- res$keggs
  
  keggnames <- subset(keggnames, KEGGID %in% keggs$KEGGID)
  reskeggnames <- merge(keggs, keggnames, by = c('KEGGID'), all = TRUE)
  reskeggnames <- unique(reskeggnames)
  row.names(reskeggnames) <- 1:nrow(reskeggnames)
  reskeggnames <- reskeggnames[c('ENTREZID', 'KEGGID', 'KEGG')]
  res$keggs <- reskeggnames
  
  
  reactomes <- res$reactomes
  library(reactome.db)
  suppressMessages(reactomenames <- 
                     AnnotationDbi::select(x = reactome.db, 
                                           keys = unique(reactomes$REACTOMEID), 
                                           columns = 'PATHNAME', keytype = 'PATHID'))
  names(reactomenames) <- c('REACTOMEID', 'REACTOME')
  reactomenames <- subset(reactomenames, REACTOMEID %in% reactomenames$REACTOMEID)
  
  reactomes <- merge(x = reactomes, y = reactomenames, by = 'REACTOMEID', all = TRUE)
  reactomes <- reactomes[c('ENTREZID', 'REACTOMEID', 'REACTOME')]
  row.names(reactomes) <- 1:nrow(reactomes)
  
  reactomes$REACTOME <- gsub(pattern = '^Homo sapiens\\: ', replacement = '', 
                             reactomes$REACTOME)
  reactomes$REACTOME <- gsub(pattern = '^Mus musculus\\: ', replacement = '', 
                             reactomes$REACTOME)
  
  res$reactomes <- reactomes
  
  return(res)
  
}

#humanbacks <- makenodebackground(species = 'human')
#mousebacks <- makenodebackground(species = 'mouse')

#saveRDS(humanbacks, 'humanbacks.rds')
#saveRDS(mousebacks, 'mousebacks.rds')

searchgos <- function(backgos, 
                      entrezdats){
  
  gores <- merge(entrezdats, backgos, by = 'ENTREZID')
  
  if(!is.null(gores)){
    if(nrow(gores) > 0){
      row.names(gores) <- 1:nrow(gores)
    }
  }
  
  return(gores)
}

searchtable <- function(allgenes = NULL, 
                        allentrez = NULL, 
                        allsymbol = NULL, 
                        species = 'human'){
  
  if(species == 'mouse'){
    #backs <- 
    #  readRDS('C:/Users/yuabr/Desktop/Transfer/codetransfer/emethyl/files/mousebacks.rds')
    backs <- get('mousebacks')
  }else{
    #backs <- 
    #  readRDS('C:/Users/yuabr/Desktop/Transfer/codetransfer/emethyl/files/humanbacks.rds')
    backs <- get('humanbacks')
  }
  
  if(is.null(allgenes) & is.null(allentrez) & is.null(allsymbol)){
    
    stop('Should provide the gene symbols or Entrez IDs\n')
    
  }
  
  library(org.Hs.eg.db)
  
  if(!is.null(allgenes)){
    
    allentrez <- gsub(pattern = '^.*::', replacement = '', x = allgenes)
    suppressWarnings(allentrez <- as.numeric(allentrez))
    allentrez <- as.character(allentrez)
    
    allsymbol <- gsub(pattern = '::.*$', replacement = '', x = allgenes)
    suppressWarnings(allsymbol[!is.na(as.numeric(allsymbol))] <- NA)
    
    
    mapping <- data.frame(GENE = allgenes, SYMBOL = allsymbol, ENTREZID = allentrez, 
                          stringsAsFactors = FALSE)
    mapping <- unique(mapping)
    mapping <- mapping[!is.na(mapping$ENTREZID), , drop = FALSE]
    
  }else if(!is.null(allentrez)){
    
    suppressWarnings(allentrez <- as.numeric(allentrez))
    allentrez <- as.character(allentrez)
    allentrez <- allentrez[!is.na(allentrez)]
    allentrez <- unique(allentrez)
    allentrez <- as.character(allentrez)
    
    suppressMessages(mapping <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                                      keys = allentrez, 
                                                      columns = 'SYMBOL', 
                                                      keytype = 'ENTREZID'))
    
    mapping <- unique(mapping)
    mapping <- mapping[c('SYMBOL', 'ENTREZID')]
    mapping$GENE <- mapping$ENTREZID
    mapping <- mapping[c('GENE', 'SYMBOL', 'ENTREZID')]
    mapping <- unique(mapping)
    mapping <- mapping[!is.na(mapping$ENTREZID), , drop = FALSE]
    
  }else if(!is.null(allsymbol)){
    
    suppressWarnings(allsymbol <- allsymbol[is.na(as.numeric(allsymbol))])
    
    suppressMessages(mapping <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                                      keys = allsymbol, 
                                                      columns = 'ENTREZID', 
                                                      keytype = 'SYMBOL'))
    mapping <- unique(mapping)
    mapping <- mapping[c('SYMBOL', 'ENTREZID')]
    mapping$GENE <- mapping$SYMBOL
    mapping <- mapping[c('GENE', 'SYMBOL', 'ENTREZID')]
    mapping <- unique(mapping)
    mapping <- mapping[!is.na(mapping$ENTREZID), , drop = FALSE]
    
  }
  
  allbps <- searchgos(backgos = backs$gobps, 
                      entrezdats = mapping)
  allccs <- searchgos(backgos = backs$goccs, 
                      entrezdats = mapping)
  allmfs <- searchgos(backgos = backs$gomfs, 
                      entrezdats = mapping)
  allkeggs <- searchgos(backgos = backs$keggs, 
                        entrezdats = mapping)
  allpaths <- searchgos(backgos = backs$reactomes, 
                        entrezdats = mapping)
  
  colnames(allbps)[colnames(allbps) == 'BPID'] <- 
    colnames(allccs)[colnames(allccs) == 'CCID'] <- 
    colnames(allmfs)[colnames(allmfs) == 'MFID'] <- 
    colnames(allkeggs)[colnames(allkeggs) == 'KEGGID'] <- 
    colnames(allpaths)[colnames(allpaths) == 'REACTOMEID'] <- 'GOID'
  
  colnames(allkeggs)[colnames(allkeggs) == 'KEGG'] <- 'TERM'
  colnames(allpaths)[colnames(allpaths) == 'REACTOME'] <- 'TERM'
  
  
  res <- list(mapping = mapping, 
              allbps = allbps, 
              allpaths = allpaths, 
              allccs = allccs, 
              allmfs = allmfs, 
              allkeggs = allkeggs)
  
  return(res)
  
}

addzeroedges <- function(shufflegenepairweights, 
                         count = FALSE){
  
  uniquenodenum <- length(unique(c(shufflegenepairweights$fromNode, 
                                   shufflegenepairweights$toNode)))
  
  completeshufflegenepairweights <- matrix(data = 0, 
                                           nrow = uniquenodenum, 
                                           ncol = uniquenodenum)
  row.names(completeshufflegenepairweights) <- 
    colnames(completeshufflegenepairweights) <- 
    unique(c(shufflegenepairweights$fromNode, 
             shufflegenepairweights$toNode))
  fromnodeidx <- match(shufflegenepairweights$fromNode, 
                       colnames(completeshufflegenepairweights))
  tonodeidx <- match(shufflegenepairweights$toNode, 
                     colnames(completeshufflegenepairweights))
  
  if(count == TRUE){
    shufflegenepairweights$weight[shufflegenepairweights$weight != 0] <- 1
  }
  
  
  idx <- matrix(c(fromnodeidx, tonodeidx), byrow = FALSE, ncol = 2)
  completeshufflegenepairweights[idx] <- shufflegenepairweights$weight
  idx <- matrix(c(tonodeidx, fromnodeidx), byrow = FALSE, ncol = 2)
  completeshufflegenepairweights[idx] <- shufflegenepairweights$weight
  
  inf <- min(completeshufflegenepairweights) - 1
  
  completecyt <- 
    WGCNA::exportNetworkToCytoscape(completeshufflegenepairweights, 
                                    edgeFile = NULL, 
                                    nodeFile = NULL, 
                                    weighted = TRUE, 
                                    threshold = inf, 
                                    nodeNames = row.names(
                                      completeshufflegenepairweights))
  completeshufflegenepairweights <- completecyt$edgeData
  rm(completecyt)
  completeshufflegenepairweights$ME <- unique(shufflegenepairweights$ME)
  shufflegenepairweights <- completeshufflegenepairweights[,colnames(
    shufflegenepairweights)[colnames(shufflegenepairweights) %in% 
                              colnames(completeshufflegenepairweights)], 
    drop = FALSE]
  rm(completeshufflegenepairweights)
  
  return(shufflegenepairweights)
  
}

calweights <- function(block){
  
  goname <- unique(block$GOID)
  termname <- unique(block$TERM)
  numblock <- block$weight
  numsum <- sum(numblock)
  resblock <- data.frame(GOID = goname, TERM = termname, weight = numsum, 
                         stringsAsFactors = FALSE)
  resblock <- unique(resblock)
  row.names(resblock) <- 1:nrow(resblock)
  return(resblock)
  
}

singlefisher <- function(line){
  
  mat <- as.numeric(c(line[3], line[4], line[5], line[6]))
  mat <- matrix(mat, byrow = TRUE, ncol = 2)
  
  fisherres <- fisher.test(x = mat, alternative = 'greater')
  #fisherres <- fisher.test(x = mat)
  
  pval <- fisherres$p.value
  
  return(pval)
  
}

weightnorm <- function(weightsvector){
  
  weightsmax <- max(weightsvector)
  weightsmin <- min(weightsvector)
  
  if(weightsmax == 1 & weightsmin == 1){
    normweights <- weightsvector
  }else{
    #preweights <- (weightsmax - weightsvector)/(weightsmax - weightsmin)
    preweights <- (weightsvector - weightsmin)/(weightsmax - weightsmin)
    normweights <- preweights - mean(preweights) + 1
  }
  
  return(normweights)
  
}

lgammacombination <- function(N, n){
  
  lres <- lgamma(N + 1) - lgamma(n + 1) - lgamma(N - n + 1)
  
  return(lres)
  
}

integratehyper <- function(total, totalblack, sub, subblack){
  
  singlegammahyperfunction <- function(Subblack, 
                                       Total = total, 
                                       Totalblack = totalblack, 
                                       Sub = sub){
    Totalwhite <- Total - Totalblack
    Subwhite <- Sub - Subblack
    
    lres <- lgammacombination(N = Totalblack, n = Subblack) + 
      lgammacombination(N = Totalwhite, n = Subwhite) - 
      lgammacombination(N = Total, n = Sub)
    
    res <- exp(lres)
    
    return(res)
    
    
  }
  
  if(as.integer(total) == total & 
     as.integer(totalblack) == totalblack & 
     as.integer(sub) == sub & 
     as.integer(subblack) == subblack){
    
    res <- 0
    
    if(min(totalblack, sub) >= subblack + 1){
      
      for(i in seq(subblack + 1, min(totalblack, sub), 1)){
        res <- res + singlegammahyperfunction(Subblack = i, 
                                              Total = total, 
                                              Totalblack = totalblack, 
                                              Sub = sub)
      }
      
    }
    
  }else{
    res <- integrate(f = singlegammahyperfunction, lower = subblack, 
                     upper = min(totalblack, sub))
    res <- res$value
    
  }
  
  
  return(res)
  
}

singleenrich <- function(line, hypertest = TRUE){
  
  if(hypertest == TRUE){
    
    pval <- integratehyper(total = as.numeric(line[6]), 
                           totalblack = as.numeric(line[5]), 
                           sub = as.numeric(line[4]), 
                           subblack = as.numeric(line[3]))
    
  }else{
    
    pval <- singlefisher(line = line)
    
  }
  
  
  odds <- (as.numeric(line[3])/as.numeric(line[4]))/
    (as.numeric(line[5])/as.numeric(line[6]))
  
  res <- c(pval, odds)
  names(res) <- c('pval', 'odds')
  
  return(res)
  
}

genego <- function(geneweights, 
                   backweights = NULL, 
                   backgos, 
                   hypertest = TRUE, 
                   genecount = FALSE){
  
  if(is.null(backweights)){
    genecount <- TRUE
  }
  
  if(genecount == TRUE){
    
    genecounts <- table(geneweights$GENE)
    genecounts <- genecounts[geneweights$GENE]
    geneweights$weight <- genecounts
    
    if(!is.null(backweights)){
      backcounts <- table(backweights$GENE)
      backcounts <- backcounts[backweights$GENE]
      backweights$weight <- backcounts
    }
    
  }
  
  
  
  if(!is.null(backweights)){
    
    if(hypertest == TRUE){
      
      backweights$weight <- weightnorm(weightsvector = backweights$weight)
      
    }
    
    geneweights <- backweights[match(geneweights$GENE, backweights$GENE),]
    
    geneweightscols <- c('GENE', setdiff(colnames(geneweights), 
                                         colnames(backgos)))
    
    geneweights <- merge(backgos, geneweights[, geneweightscols, drop = FALSE], 
                         by = c('GENE'))
    
    geneweights <- geneweights[order(geneweights$GENE, geneweights$GOID),]
    
    row.names(geneweights) <- 1:nrow(geneweights)
    
    
    generes <- plyr::ddply(.data = geneweights, .variables = c('GOID'), 
                           .fun = calweights)
    
    backweightscols <- c('GENE', setdiff(colnames(backweights), 
                                         colnames(backgos)))
    
    backweights <- merge(backgos, backweights[, backweightscols, drop = FALSE], 
                         by = c('GENE'))
    
    backweights <- backweights[order(backweights$GENE, backweights$GOID),]
    
    row.names(backweights) <- 1:nrow(backweights)
    
    backres <- plyr::ddply(.data = backweights, .variables = c('GOID'), 
                           .fun = calweights)
    
  }else{
    
    geneweightscols <- c('GENE', setdiff(colnames(geneweights), 
                                         colnames(backgos)))
    
    geneweights <- merge(backgos, geneweights[, geneweightscols, drop = FALSE], 
                         by = c('GENE'))
    
    geneweights <- geneweights[order(geneweights$GENE, geneweights$GOID),]
    
    row.names(geneweights) <- 1:nrow(geneweights)
    
    
    generes <- plyr::ddply(.data = geneweights, .variables = c('GOID'), 
                           .fun = calweights)
    
    
    
    backgos$weight <- 1
    backres <- plyr::ddply(.data = backgos, .variables = c('GOID'), 
                           .fun = calweights)
    
  }
  
  weightres <- merge(generes, backres, by = c('GOID', 'TERM'))
  colnames(weightres) <- c('GOID', 'TERM', 'geneweight', 'backweight')
  backsum <- sum(backres$weight)
  genesum <- sum(generes$weight)
  weightres$genesum <- genesum
  weightres$backsum <- backsum
  weightres <- weightres[,c('GOID', 'TERM', 
                            'geneweight', 'genesum', 'backweight', 'backsum'), 
                         drop = FALSE]
  
  
  
  
  #P-val calculation start
  
  res <- list()
  i <- 1
  
  for(i in 1:nrow(weightres)){
    
    res[[i]] <- singleenrich(line = weightres[i,], hypertest = hypertest)
    
  }
  
  weightpvals <- do.call(rbind, res)
  row.names(weightpvals) <- 1:nrow(weightpvals)
  
  #P-val calculation complete
  
  
  weightpvals <- as.data.frame(weightpvals, stringsAsFactors = FALSE)
  weightpvals$GOID <- weightres$GOID
  weightpvals$TERM <- weightres$TERM
  weightpvals$padj <- p.adjust(p = weightpvals$pval, method = 'BH')
  weightpvals <- weightpvals[c('GOID', 'TERM', 'pval', 'padj', 'odds')]
  weightpvals <- weightpvals[order(weightpvals$padj, weightpvals$pval, 
                                   -weightpvals$odds),]
  row.names(weightpvals) <- 1:nrow(weightpvals)
  
  geneweights <- geneweights[order(geneweights$GOID, 
                                   geneweights$GENE, 
                                   geneweights$weight),]
  
  geneweights <- geneweights[c('GOID', 'TERM', 'GENE', 'weight', 
                               setdiff(names(geneweights), 
                                       c('GOID', 'TERM', 'GENE', 'weight', 'ME')))]
  
  
  if(nrow(geneweights) > 0){
    row.names(geneweights) <- 1:nrow(geneweights)
  }
  
  if(!is.null(backweights)){
    
    backweights <- backweights[order(backweights$GOID, 
                                     backweights$GENE, 
                                     backweights$weight),]
    
    backweights <- backweights[c('GOID', 'TERM', 'GENE', 'weight', 
                                 setdiff(names(backweights), 
                                         c('GOID', 'TERM', 'GENE', 'weight')))]
    
    row.names(backweights) <- 1:nrow(backweights)
    
  }else{
    backweights <- NULL
  }
  
  
  
  res <- list(enrichpvals = weightpvals, 
              geneweights = geneweights, 
              backweights = backweights)
  
  return(res)
  
}

subsharedgo <- function(i = 1, 
                        goids, 
                        frombackgos, 
                        tobackgos){
  
  goid <- goids[i]
  
  frombackgo <- subset(frombackgos, GOID == goid)
  tobackgo <- subset(tobackgos, GOID == goid)
  
  edgebackmat <- merge(frombackgo, tobackgo, by = c('GOID', 'TERM'))
  names(edgebackmat) <- c('GOID', 'TERM', 'fromNode', 'toNode')
  edgebackmat <- edgebackmat[c('fromNode', 'toNode', 'GOID', 'TERM')]
  edgebackmat <- subset(edgebackmat, fromNode != toNode)
  
  if(nrow(edgebackmat) == 0){
    return(edgebackmat)
  }else{
    edgebackmat <- unique(edgebackmat)
    row.names(edgebackmat) <- 1:nrow(edgebackmat)
    return(edgebackmat)
  }
  
}

sharedgo <- function(edgedat, 
                     backgos, 
                     threads = 1){
  
  backgos <- unique(backgos)
  
  frombackgos <- subset(backgos, GENE %in% edgedat$fromNode, 
                        select = c('GENE', 'GOID', 'TERM'))
  tobackgos <- subset(backgos, GENE %in% edgedat$toNode, 
                      select = c('GENE', 'GOID', 'TERM'))
  
  goids <- intersect(frombackgos$GOID, tobackgos$GOID)
  
  if(threads == 1){
    
    edgegolist <- list()
    i <- 1
    for(i in 1:length(goids)){
      
      edgegolist[[i]] <- subsharedgo(i = i, 
                                     goids = goids, 
                                     frombackgos = frombackgos, 
                                     tobackgos = tobackgos)
      
    }
    
  }else{
    
    subsharedgo <- function(i = 1, 
                            goids, 
                            frombackgos, 
                            tobackgos){
      
      goid <- goids[i]
      
      frombackgo <- subset(frombackgos, GOID == goid)
      tobackgo <- subset(tobackgos, GOID == goid)
      
      edgebackmat <- merge(frombackgo, tobackgo, by = c('GOID', 'TERM'))
      names(edgebackmat) <- c('GOID', 'TERM', 'fromNode', 'toNode')
      edgebackmat <- edgebackmat[c('fromNode', 'toNode', 'GOID', 'TERM')]
      edgebackmat <- subset(edgebackmat, fromNode != toNode)
      
      if(nrow(edgebackmat) == 0){
        return(edgebackmat)
      }else{
        edgebackmat <- unique(edgebackmat)
        row.names(edgebackmat) <- 1:nrow(edgebackmat)
        return(edgebackmat)
      }
      
    }
    
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    edgegolist <- foreach::foreach(i = 1:length(goids),
                                   #.export = ls(name = globalenv())) %dopar% {
                                   .export = NULL) %dopar% {
                                     
                                     subsharedgo(i = i, 
                                                 goids = goids, 
                                                 frombackgos = frombackgos, 
                                                 tobackgos = tobackgos)
                                   }
    
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
    
  }
  
  edgego <- do.call(rbind, edgegolist)
  
  if(nrow(edgego) > 0){
    edgego <- unique(edgego)
    row.names(edgego) <- 1:nrow(edgego)
  }
  
  return(edgego)
  
}

extractpairgos <- function(pairdat, 
                           backpairgos){
  
  sharedgenes <- intersect(pairdat$GENE, backpairgos$GENE)
  pairdat <- subset(pairdat, GENE %in% sharedgenes)
  pairdat$GOID <- backpairgos$GOID[match(pairdat$GENE, backpairgos$GENE)]
  
  pairdat$TERM <- backpairgos$TERM[match(pairdat$GENE, backpairgos$GENE)]
  
  row.names(pairdat) <- 1:nrow(pairdat)
  return(pairdat)
  
}

getsinglepairgoid <- function(line, backgenegos = backgenegos){
  
  fromgoid <- subset(backgenegos, GENE == line$fromNode)
  togoid <- subset(backgenegos, GENE == line$toNode)
  
  sharedgoids <- intersect(fromgoid$GOID, togoid$GOID)
  sharedgoids <- subset(fromgoid, GOID %in% sharedgoids, 
                        select = c('GOID', 'TERM'))
  sharedgoids <- unique(sharedgoids)
  
  sharedgoidres <- line[rep(seq(nrow(line)), each = nrow(sharedgoids)),]
  
  if(nrow(sharedgoidres) > 0){
    row.names(sharedgoidres) <- 1:nrow(sharedgoidres)
  }
  
  sharedgoidres$GOID <- sharedgoids$GOID
  sharedgoidres$TERM <- sharedgoids$TERM
  
  return(sharedgoidres)
  
}

extractpairgosfromgenegos <- function(pairdat, 
                                      backgenegos){
  
  pairdat <- subset(pairdat, fromNode %in% backgenegos$GENE & 
                      toNode %in% backgenegos$GENE)
  
  pairgodat <- plyr::ddply(.data = pairdat, .variables = c('fromNode', 'toNode'), 
                           .fun = getsinglepairgoid)
  
  row.names(pairgodat) <- 1:nrow(pairgodat)
  return(pairgodat)
  
}

calbackpairweight <- function(sub){
  
  n <- nrow(sub)
  backpairweight <- n*(n - 1)/2
  res <- data.frame(GOID = unique(sub$GOID), 
                    TERM = unique(sub$TERM), 
                    backpairweight = backpairweight, 
                    stringsAsFactors = FALSE)
  
  return(res)
  
}

pairgo <- function(genepairweights, 
                   backpairweights = NULL, 
                   backpairgos = NULL, 
                   backgos = NULL, 
                   hypertest = TRUE, 
                   paircount = FALSE, 
                   threads = 1){
  
  if(is.null(backpairgos) & is.null(backgos)){
    return(NULL)
  }
  
  if(is.null(backpairweights)){
    paircount <- TRUE
  }
  
  genepairweights$GENE <- paste0(genepairweights$fromNode, '_', 
                                 genepairweights$toNode)
  
  if(!is.null(backpairweights)){
    
    backpairweights$GENE <- paste0(backpairweights$fromNode, '_', 
                                   backpairweights$toNode)
    
    if(is.null(backpairgos)){
      #pairgo
      backpairgos <- sharedgo(edgedat = backpairweights, 
                              backgos = backgos, 
                              threads = threads)
    }
    
  }
  
  
  if(paircount == TRUE){
    
    genepaircounts <- table(genepairweights$GENE)
    genepaircounts <- genepaircounts[genepairweights$GENE]
    genepairweights$weight <- genepaircounts
    
    if(!is.null(backpairweights)){
      backpaircounts <- table(backpairweights$GENE)
      backpaircounts <- backpaircounts[backpairweights$GENE]
      backpairweights$weight <- backpaircounts
    }
    
  }
  
  if(!is.null(backpairgos)){
    backpairgosrev <- backpairgos
    backpairgosrev$fromNode <- backpairgos$toNode
    backpairgosrev$toNode <- backpairgos$fromNode
    backpairgos <- rbind(backpairgos, backpairgosrev)
    
    backpairgos$GENE <- paste0(backpairgos$fromNode, '_', 
                               backpairgos$toNode)
    
    backpairgos <- unique(backpairgos)
    
  }
  
  
  if(!is.null(backpairweights)){
    
    if(hypertest == TRUE){
      backpairweights$weight <- weightnorm(weightsvector = backpairweights$weight)
    }
    
    genepairweights <- backpairweights[match(genepairweights$GENE, backpairweights$GENE),]
    
    
    suppressMessages(
      
      genepairweights <- tryCatch({
        
        extractpairgos(pairdat = genepairweights, 
                       backpairgos = backpairgos)
        
      }, error = function(err){
        
        NULL
        
      })
      
    )
    
    if(is.null(genepairweights)){
      return(NULL)
    }
    
    backpairweights <- extractpairgos(pairdat = backpairweights, 
                                      backpairgos = backpairgos)
    
    
    
    backpairweights <- backpairweights[order(backpairweights$GENE, 
                                             backpairweights$GOID),]
    genepairweights <- genepairweights[order(genepairweights$GENE, 
                                             genepairweights$GOID),]
    row.names(backpairweights) <- 1:nrow(backpairweights)
    row.names(genepairweights) <- 1:nrow(genepairweights)
    
    
    
    
    genepairres <- plyr::ddply(.data = genepairweights, .variables = c('GOID'), 
                               .fun = calweights)
    backpairres <- plyr::ddply(.data = backpairweights, .variables = c('GOID'), 
                               .fun = calweights)
    backsum <- sum(backpairres$weight)
    
  }else{
    
    if(!is.null(backpairgos)){
      
      suppressMessages(
        
        genepairweights <- tryCatch({
          
          extractpairgos(pairdat = genepairweights, 
                         backpairgos = backpairgos)
          
        }, error = function(err){
          
          NULL
          
        })
        
      )
      
      if(is.null(genepairweights)){
        return(NULL)
      }
      
      
      
      genepairweights <- genepairweights[order(genepairweights$GENE, 
                                               genepairweights$GOID),]
      row.names(genepairweights) <- 1:nrow(genepairweights)
      
      genepairres <- plyr::ddply(.data = genepairweights, .variables = c('GOID'), 
                                 .fun = calweights)
      
      
      #backpairgos$weight <- 1
      backpairgos$weight <- 0.5
      
      backpairres <- plyr::ddply(.data = backpairgos, .variables = c('GOID'), 
                                 .fun = calweights)
      
      
      backsum <- sum(backpairres$weight)
      
    }else{
      
      suppressMessages(
        
        genepairweights <- tryCatch({
          
          extractpairgosfromgenegos(pairdat = genepairweights, 
                                    backgenegos = backgos)
          
        }, error = function(err){
          
          NULL
          
        })
        
      )
      
      if(is.null(genepairweights)){
        return(NULL)
      }
      
      
      genepairweights <- genepairweights[order(genepairweights$GENE, 
                                               genepairweights$GOID),]
      row.names(genepairweights) <- 1:nrow(genepairweights)
      
      genepairres <- plyr::ddply(.data = genepairweights, .variables = c('GOID'), 
                                 .fun = calweights)
      
      backgosub <- subset(backgos, GOID %in% unique(c(genepairweights$GOID)))
      backgosub <- unique(backgosub)
      
      backpairres <- plyr::ddply(.data = backgosub, .variables = c('GOID'), 
                                 .fun = calbackpairweight)
      backsum <- nrow(backgos)*(nrow(backgos) - 1)/2
      
    }
    
  }
  
  
  weightres <- merge(genepairres, backpairres, by = c('GOID', 'TERM'))
  
  colnames(weightres) <- c('GOID', 'TERM', 'genepairweight', 'backpairweight')
  
  genesum <- sum(genepairres$weight)
  
  weightres$genepairsum <- genesum
  weightres$backpairsum <- backsum
  
  
  weightres <- weightres[,c('GOID', 'TERM', 'genepairweight', 'genepairsum', 
                            'backpairweight', 'backpairsum'), 
                         drop = FALSE]
  
  #P-val calculation start
  
  res <- list()
  i <- 1
  
  for(i in 1:nrow(weightres)){
    
    res[[i]] <- singleenrich(line = weightres[i,], hypertest = hypertest)
    
  }
  
  weightpvals <- do.call(rbind, res)
  row.names(weightpvals) <- 1:nrow(weightpvals)
  
  #P-val calculation complete
  
  weightpvals <- as.data.frame(x = weightpvals, stringsAsFactors = FALSE)
  weightpvals$GOID <- weightres$GOID
  weightpvals$TERM <- weightres$TERM
  weightpvals$padj <- p.adjust(p = weightpvals$pval, method = 'BH')
  weightpvals <- weightpvals[c('GOID', 'TERM', 'pval', 'padj', 'odds')]
  
  weightpvals <- weightpvals[order(weightpvals$padj, weightpvals$pval, 
                                   -weightpvals$odds),]
  row.names(weightpvals) <- 1:nrow(weightpvals)
  
  genepairweights <- genepairweights[order(genepairweights$GOID, 
                                           genepairweights$GENE, 
                                           genepairweights$weight),]
  
  
  genepairweights <- genepairweights[c('GOID', 'TERM', 'fromNode', 'toNode', 'weight', 
                                       setdiff(names(genepairweights), 
                                               c('GOID', 'TERM', 'fromNode', 
                                                 'toNode', 'weight', 'GENE', 'ME')))]
  
  
  
  if(nrow(genepairweights) > 0){
    row.names(genepairweights) <- 1:nrow(genepairweights)
  }
  
  if(!is.null(backpairweights)){
    backpairweights <- backpairweights[order(backpairweights$GOID, 
                                             backpairweights$GENE, 
                                             backpairweights$weight),]
    
    backpairweights <- backpairweights[c('GOID', 'TERM', 'fromNode', 'toNode', 'weight', 
                                         setdiff(names(backpairweights), 
                                                 c('GOID', 'TERM', 'fromNode', 
                                                   'toNode', 'weight', 'GENE')))]
    
    
    
    row.names(backpairweights) <- 1:nrow(backpairweights)
    
  }else{
    backpairweights <- NULL
  }
  
  
  res <- list(enrichpvals = weightpvals, 
              genepairweights = genepairweights, 
              backpairweights = backpairweights)
  
  return(res)
  
}

calshufflepval <- function(i = 1, 
                           generes, 
                           shufflenum, 
                           weights, 
                           seed){
  
  counts <- unique(generes$count)
  counts <- counts[order(counts)]
  
  subcounts <- subset(generes, count == counts[i])
  
  
  sampledweights <- c()
  
  for(j in 1:shufflenum){
    set.seed(seed + j)
    sampledweight <- sample(x = weights, 
                            size = counts[i], 
                            replace = FALSE)
    sampledweights <- c(sampledweights, sum(sampledweight))
    
  }
  
  shufflepvals <- c()
  for(k in 1:nrow(subcounts)){
    
    shufflepval <- ecdf(sampledweights)(subcounts$weight[k])
    shufflepvals <- c(shufflepvals, shufflepval)
    
  }
  
  names(shufflepvals) <- subcounts$GOID
  
  return(shufflepvals)
  
}

adjshuffleres <- function(shuffleres){
  
  pvals1 <- 1 - shuffleres$quantile
  pvals2 <- shuffleres$quantile
  
  pvals <- pvals1
  pvals[pvals > 0.5] <- pvals2[pvals > 0.5]
  
  pvals <- 2*pvals
  
  shuffleres$pval <- pvals
  shuffleres$padj <- p.adjust(p = pvals, method = 'BH')
  
  shuffleres <- shuffleres[order(-shuffleres$quantile, 
                                 shuffleres$padj, shuffleres$pval, 
                                 -shuffleres$weight, -shuffleres$count),]
  
  if(nrow(shuffleres) > 0){
    row.names(shuffleres) <- 1:nrow(shuffleres)
  }
  
  return (shuffleres)
  
}

geneshufflego <- function(geneweights, 
                          genegos, 
                          shufflenum = 1000, 
                          seed = 2022, 
                          threads = 1, 
                          genecount = FALSE, 
                          shuffletwotail = FALSE){
  
  
  if(genecount == TRUE){
    
    genecounts <- table(geneweights$GENE)
    genecounts <- genecounts[geneweights$GENE]
    geneweights$weight <- genecounts
    
  }
  
  
  #weights <- geneweights$weight
  
  geneweightscols <- c('GENE', setdiff(colnames(geneweights), 
                                       colnames(genegos)))
  
  geneweights <- merge(genegos, geneweights[, geneweightscols, drop = FALSE], 
                       by = c('GENE'))
  
  
  
  weights <- geneweights$weight
  
  
  
  generes <- plyr::ddply(.data = geneweights, .variables = c('GOID'), .fun = calweights)
  
  gocounts <- table(geneweights$GOID)
  gocounts <- gocounts[generes$GOID]
  generes$count <- as.vector(gocounts)
  
  iseq <- seq(1, length(unique(generes$count)), 1)
  
  if(threads == 1){
    
    shufflepvals <- list()
    i = 1
    for(i in iseq){
      
      shufflepvals[[i]] <- calshufflepval(i = i, 
                                          generes = generes, 
                                          shufflenum = shufflenum, 
                                          weights = weights, 
                                          seed = seed)
      
    }
    
    
  }else{
    
    calshufflepval <- function(i = 1, 
                               generes, 
                               shufflenum, 
                               weights, 
                               seed){
      
      counts <- unique(generes$count)
      counts <- counts[order(counts)]
      
      subcounts <- subset(generes, count == counts[i])
      
      
      sampledweights <- c()
      
      for(j in 1:shufflenum){
        set.seed(seed + j)
        sampledweight <- sample(x = weights, 
                                size = counts[i], 
                                replace = FALSE)
        sampledweights <- c(sampledweights, sum(sampledweight))
        
      }
      
      shufflepvals <- c()
      for(k in 1:nrow(subcounts)){
        
        shufflepval <- ecdf(sampledweights)(subcounts$weight[k])
        shufflepvals <- c(shufflepvals, shufflepval)
        
      }
      
      names(shufflepvals) <- subcounts$GOID
      
      return(shufflepvals)
      
    }
    
    
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    shufflepvals <- foreach::foreach(i = iseq, 
                                     #.export = ls(name = globalenv())) %dopar% {
                                     .export = NULL) %dopar% {
                                       
                                       calshufflepval(i = i, 
                                                      generes = generes, 
                                                      shufflenum = shufflenum, 
                                                      weights = weights, 
                                                      seed = seed)
                                     }
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
    
    
  }
  
  shufflepvals <- unlist(shufflepvals)
  shufflepvals <- shufflepvals[generes$GOID]
  
  generes$pval <- 1 - shufflepvals
  generes$padj <- p.adjust(p = generes$pval, method = 'BH')
  generes$quantile <- shufflepvals
  generes <- generes[c('GOID', 'TERM', 'pval', 'padj', 'quantile', 
                       'weight', 'count')]
  generes <- generes[order(generes$padj, generes$pval, 
                           -generes$weight, -generes$count),]
  row.names(generes) <- 1:nrow(generes)
  
  if('ME' %in% colnames(geneweights)){
    ME <- 'ME'
  }else{
    ME <- NULL
  }
  
  geneweights <- geneweights[c('GENE', 'GOID', 'weight', ME, 
                               setdiff(colnames(geneweights), 
                                       c('GENE', 'GOID', 'weight', 'ME')))]
  
  if(shuffletwotail == TRUE){
    if(!is.null(generes)){
      if(nrow(generes) > 0){
        
        generes <- adjshuffleres(shuffleres = generes)
        
      }
    }
    
  }
  
  
  
  res <- list(shuffleenrichpvals = generes, 
              geneweights = geneweights)
  
  
  return(res)
  
}

pairshufflego <- function(genepairweights, 
                          genepairgos = NULL, 
                          genegos = NULL,
                          shufflenum = 1000, 
                          seed = 2022, 
                          threads = 1, 
                          paircount = FALSE, 
                          addzeroedgesforshuffle = TRUE, 
                          shuffletwotail = FALSE){
  
  if(is.null(genepairgos) & is.null(genegos)){
    return(NULL)
  }
  
  genepairweights$GENE <- paste0(genepairweights$fromNode, '_', 
                                 genepairweights$toNode)
  
  if(is.null(genepairgos)){
    #pairgo
    genepairgos <- sharedgo(edgedat = genepairweights, 
                            backgos = genegos, 
                            threads = threads)
  }
  
  genepairgos$GENE <- paste0(genepairgos$fromNode, '_', 
                             genepairgos$toNode)
  
  
  if(paircount == TRUE){
    
    genepairweights$weight[genepairweights$weight != 0] <- 1
    uniquenodenum <- length(unique(c(genepairweights$fromNode, 
                                     genepairweights$toNode)))
    if(uniquenodenum*(uniquenodenum - 1)/2 != nrow(genepairweights)){
      
      genepairweights <- addzeroedges(genepairweights = genepairweights, 
                                      count = TRUE)
      
    }
    
  }
  
  if(paircount == FALSE & addzeroedgesforshuffle == TRUE){
    
    uniquenodenum <- length(unique(c(genepairweights$fromNode, 
                                     genepairweights$toNode)))
    if(uniquenodenum*(uniquenodenum - 1)/2 != nrow(genepairweights)){
      
      genepairweights <- addzeroedges(genepairweights = shufflegenepairweights, 
                                      count = FALSE)
      
      
    }
    
  }
  
  genepairweights$GENE <- paste0(genepairweights$fromNode, '_', 
                                 genepairweights$toNode)
  
  #weights <- genepairweights$weight
  
  genepairweightscols <- c('GENE', setdiff(colnames(genepairweights), 
                                           colnames(genepairgos)))
  
  genepairweights <- merge(genepairgos, genepairweights[, genepairweightscols, 
                                                        drop = FALSE], 
                           by = c('GENE'))
  
  
  genepairweights <- genepairweights[c('GENE', 'fromNode', 'toNode', 
                                       'GOID', 'TERM', 'weight', 'ME')]
  names(genepairweights) <- c('GENE', 'fromNode', 'toNode', 
                              'GOID', 'TERM', 'weight', 'ME')
  
  
  
  weights <- genepairweights$weight
  
  
  
  genepairres <- plyr::ddply(.data = genepairweights, 
                             .variables = c('GOID'), 
                             .fun = calweights)
  
  gopaircounts <- table(genepairweights$GOID)
  gopaircounts <- gopaircounts[genepairres$GOID]
  genepairres$count <- as.vector(gopaircounts)
  
  iseq <- seq(1, length(unique(genepairres$count)), 1)
  
  if(threads == 1){
    
    shufflepvals <- list()
    i = 1
    for(i in iseq){
      
      shufflepvals[[i]] <- calshufflepval(i = i, 
                                          generes = genepairres, 
                                          shufflenum = shufflenum, 
                                          weights = weights, 
                                          seed = seed)
      
    }
    
    
  }else{
    
    calshufflepval <- function(i = 1, 
                               generes, 
                               shufflenum, 
                               weights, 
                               seed){
      
      counts <- unique(generes$count)
      counts <- counts[order(counts)]
      
      subcounts <- subset(generes, count == counts[i])
      
      
      sampledweights <- c()
      
      for(j in 1:shufflenum){
        set.seed(seed + j)
        sampledweight <- sample(x = weights, 
                                size = counts[i], 
                                replace = FALSE)
        sampledweights <- c(sampledweights, sum(sampledweight))
        
      }
      
      shufflepvals <- c()
      for(k in 1:nrow(subcounts)){
        
        shufflepval <- ecdf(sampledweights)(subcounts$weight[k])
        shufflepvals <- c(shufflepvals, shufflepval)
        
      }
      
      names(shufflepvals) <- subcounts$GOID
      
      return(shufflepvals)
      
    }
    
    
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    shufflepvals <- foreach::foreach(i = iseq, 
                                     #.export = ls(name = globalenv())) %dopar% {
                                     .export = NULL) %dopar% {
                                       calshufflepval(i = i, 
                                                      generes = genepairres, 
                                                      shufflenum = shufflenum, 
                                                      weights = weights, 
                                                      seed = seed)
                                     }
    
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
    
    
  }
  
  shufflepvals <- unlist(shufflepvals)
  shufflepvals <- shufflepvals[genepairres$GOID]
  
  genepairres$pval <- 1 - shufflepvals
  genepairres$padj <- p.adjust(p = genepairres$pval, method = 'BH')
  genepairres$quantile <- shufflepvals
  genepairres <- genepairres[c('GOID', 'TERM', 'pval', 'padj', 'quantile', 
                               'weight', 'count')]
  genepairres <- genepairres[order(genepairres$padj, genepairres$pval, 
                                   -genepairres$weight, -genepairres$count),]
  row.names(genepairres) <- 1:nrow(genepairres)
  
  if('ME' %in% colnames(genepairweights)){
    ME <- 'ME'
  }else{
    ME <- NULL
  }
  
  genepairweights <- genepairweights[c('fromNode', 'toNode', 'GOID', 
                                       'weight', ME, 
                                       setdiff(colnames(genepairweights), 
                                               c('GENE', 'fromNode', 'toNode', 
                                                 'GOID', 'weight', 'ME')))]
  
  if(shuffletwotail == TRUE){
    if(!is.null(genepairres)){
      if(nrow(genepairres) > 0){
        
        genepairres <- adjshuffleres(shuffleres = genepairres)
        
      }
    }
    
  }
  
  
  
  res <- list(shuffleenrichpvals = genepairres, 
              genepairweights = genepairweights)
  
  return(res)
  
}

maingo <- function(geneweights,
                   backweights, 
                   backgos, 
                   backpairweights, 
                   genegos, 
                   genepairweights, 
                   shufflegenepairweights, 
                   genecount = FALSE, 
                   paircount = FALSE, 
                   addzeroedgesforshuffle = TRUE, 
                   threads = 1, 
                   enrichtype = c('gene', 'pair'), 
                   seed = 2022, 
                   shufflenum = 1000, 
                   shuffletwotail = FALSE, 
                   pairenrichgo = TRUE){
  
  if('gene' %in% enrichtype){
    #genego
    genegores <- genego(geneweights = geneweights, 
                        backweights = backweights, 
                        backgos = backgos, 
                        hypertest = TRUE, 
                        genecount = genecount)
    
  }else{
    genegores <- NULL
  }
  
  if('pair' %in% enrichtype & pairenrichgo == TRUE){
    
    #pairgo
    backpairgos <- sharedgo(edgedat = backpairweights, 
                            backgos = backgos, 
                            threads = threads)
    
    #The edges for pairgo is the edges in the top 1% edges set of the whole 
    #genome, not the complete as pairshufflego
    pairgores <- tryCatch({
      
      pairgo(genepairweights = genepairweights, 
             backpairweights = backpairweights, 
             backpairgos = backpairgos, 
             hypertest = TRUE, 
             paircount = paircount, 
             threads = threads)
      
    }, error = function(err){
      
      NULL
      
    })
    
    
  }else{
    
    pairgores <- NULL
    
  }
  
  if(('gene' %in% enrichtype) & !is.null(shufflenum)){
    
    #geneshufflego
    geneshufflegores <- geneshufflego(geneweights = geneweights, 
                                      genegos = genegos, 
                                      shufflenum = shufflenum, 
                                      seed = seed, 
                                      threads = threads, 
                                      genecount = genecount, 
                                      shuffletwotail = shuffletwotail)
    
  }else{
    geneshufflegores <- NULL
  }
  
  if(('pair' %in% enrichtype) & !is.null(shufflenum)){
    
    #pairshufflego
    megenepairgos <- sharedgo(edgedat = shufflegenepairweights, 
                              backgos = backgos, 
                              threads = threads)
    
    #The edges for pairshufflego is the complete edges, not the top 1% as 
    #pairgo
    pairshufflegores <- pairshufflego(genepairweights = shufflegenepairweights, 
                                      genepairgos = megenepairgos, 
                                      shufflenum = shufflenum, 
                                      seed = seed, 
                                      threads = threads, 
                                      paircount = paircount, 
                                      addzeroedgesforshuffle = addzeroedgesforshuffle, 
                                      shuffletwotail = shuffletwotail)
    
  }else{
    pairshufflegores <- NULL
  }
  
  
  res <- list(generes = genegores$enrichpvals, 
              pairres = pairgores$enrichpvals, 
              geneshuffleres = geneshufflegores$shuffleenrichpvals, 
              pairshuffleres = pairshufflegores$shuffleenrichpvals, 
              
              geneweights = genegores$geneweights, 
              pairweights = pairgores$genepairweights, 
              pairshuffleweights = pairshufflegores$genepairweights)
  
  
  return(res)
  
}

#'Perform structure-based gene functional enrichment or normal enrichment
#'
#'Perform structure-based gene functional enrichment or normal enrichment. 
#'
#'@param geneweights If need to perform enrichment analysis on the gene level, 
#'  this parameter can be used. It accepts a data frame with columns as "GENE" 
#'  and "weight". The former records gene symbols or ENTREZ IDs, or mixed gene 
#'  names as sybmol::ENTREZ, such as "CACNA1S::779", where "CACNA1S" is the 
#'  gene symbol and "779" is the ENTREZ ID. The other column "weight" records 
#'  the gene weights, so that a gene weight-based enrichment analysis can be 
#'  performed from it. The traditional gene count-based analysis is a special 
#'  case of this gene weight-based one, with the weights of all the genes as 
#'  1. This parameter can also accept a vector of gene symbols or ENTREZ IDs, 
#'  or the mixed gene names, and all the genes in this vector will be assigned 
#'  a weight as 1, so the gene count-based analysis will be performed. Default 
#'  is NULL, which means the enrichment will not be performed on gene level, 
#'  because this function can also perform the enrichment on gene-gene pair 
#'  level.
#'@param backweights The background for the gene level enrichment. Similar to 
#'  \code{geneweights}, it can be a data frame with columns named as "GENE" 
#'  and "weight", or a vector recording the gene names. The difference is that 
#'  it is used to transfer the gene background to perform the hypergeometric 
#'  enrichment, not to transfer the genes to be analyzed. Default is NULL, and 
#'  in this case, all the human genomic genes will be used as the background. 
#'  In this case, the weights of all the genes in the background will be 1 and 
#'  the enrichment will be the gene count-based one, not the weight-based one. 
#'  The background is only needed for the hypergeometric enrichment analysis. 
#'  The function can also perform a gene weight shuffling enrichment, and in 
#'  this case, the background is not needed.
#'@param genepairweights If need to perform enrichment analysis on the gene- 
#'  gene pair level, this parameter can be used. It accepts a data frame with 
#'  columns as "fromNode", "toNode" and "weight". The "fromNode" and "toNode" 
#'  columns are similar to the "GENE" column in \code{geneweights}, recording 
#'  gene symbols or ENTREZ IDs, or mixed gene names. Another column "weight" 
#'  records the gene pair weights, so that a gene pair weight-based enrichment 
#'  analysis can be performed from it. It can also perform a gene pair count- 
#'  based enrichment if all the values in the column "weight" are 1, or there 
#'  is no column named "weight" in the data frame. The gene pair count-based 
#'  enrichment is a special case of the weight-based one. The default value of 
#'  this parameter is NULL, which means the enrichment will not be performed 
#'  on gene-gene pair level, because this function can also perform the a gene 
#'  level enrichment analysis.
#'@param backpairweights The background for the gene pair level enrichment. 
#'  Similar to \code{genepairweights}, it should be a data frame with columns 
#'  named as "fromNode", "toNode" and "weight". The difference is that it is 
#'  used to transfer the gene pair background to perform the hypergeometric 
#'  enrichment, not to transfer the gene pairs to be analyzed. Its default is 
#'  NULL, and in this case, the gene pair weight-based hypergeometric analysis 
#'  will not be performed, because this gene pair background is only necessary 
#'  for hypergeometric enrichment, and the function will also perform a gene 
#'  pair weight shuffling enrichment, where the background is not needed.
#'@param genecount For the enrichment on gene level, if this parameter is set 
#'  as TRUE, all the gene weights will be set as 1 and the enrichment will be 
#'  performed on the gene count level. Default is FALSE.
#'@param paircount For the enrichment on gene pair level, if this parameter is 
#'  set as TRUE, all the pair weights will be set as 1 and the enrichment will 
#'  be performed on the gene pair count level. Default is FALSE. 
#'@param addzeroedgesforshuffle For the gene pair weight shuffling enrichment, 
#'  if 2 genes do not constitute a pair, e.g., in a gene network, there is no 
#'  edge between 2 genes, but this parameter is set as TRUE, the gene pair can 
#'  still participate in the gene pair weight shuffling process, just with a 
#'  weight of 0. If it is set as FALSE, this gene pair will not participate in 
#'  the shuffling process. Default is TRUE. 
#'@param threads Number of threads to be used for parallelization, default is 
#'  1. 
#'@param datasests The gene functional datasets to be used for the gene and 
#'  gene pair level enrichment, can be a vector including any of the 5 dataset 
#'  names: "bp", "cc", "mf", "kegg", and "reactome", representing the GOBP, 
#'  GOCC, GOMF, KEGG, and REACTOME datasets. Default is with all these 5. 
#'@param combinedatasetreses If transferring multiple functional datasets to 
#'  the parameter \code{datasets}, and this parameter is set as FALSE, each 
#'  dataset will get a separate results set for the enrichment from it. On the 
#'  other hand, if this parameter is TRUE, all the results will be combined 
#'  together and then returned. Default is TRUE. 
#'@param seed Random seed for the gene or gene pair weight shuffling method. 
#'  Default is 2022.
#'@param shufflenum The shuffling number for the weight shuffling method. Its 
#'  default value is 1000.
#'@param shuffletwotail Whether the enrichment p-val returned by the shuffling 
#'  method should be two tailed or single tailed. Default is FALSE, meaning it 
#'  should be a single tailed one, representing the probability of obtaining a 
#'  greater weight sum for the genes or pairs belonging to a functional term. 
#'@return A list with several slots. If gene level enrichment is performed, it 
#'  will contain a slot named "geneenrichres", which includes the gene level 
#'  hypergeometric enrichment, weight shuffling enrichment results, and the 
#'  functional terms of each gene. If gene pair level enrichment is performed, 
#'  it will contain a slot named "pairenrichres", which includes the gene pair 
#'  level hypergeometric enrichment, weight shuffling enrichment results, and 
#'  the functional terms of each gene pair, which are the intersection terms 
#'  of its 2 genes' functional terms. These results can also be the results 
#'  based on gene counts and pair counts, instead of weights. If the parameter 
#'  \code{combinedatasetreses} is FALSE, the results can be further split into 
#'  different dataset results according to the ones assigned by the parameter 
#'  \code{datasets}.
#'@export
topoenrich <- function(geneweights = NULL, 
                       backweights = NULL, 
                       
                       genepairweights = NULL, 
                       backpairweights = NULL, 
                       
                       genecount = FALSE, 
                       paircount = FALSE, 
                       addzeroedgesforshuffle = TRUE, 
                       #species = 'human', 
                       threads = 1, 
                       datasets = c('bp', 'cc', 'mf', 'kegg', 'reactome'), 
                       combinedatasetreses = TRUE, 
                       seed = 2022, 
                       shufflenum = 1000, 
                       shuffletwotail = FALSE){
  
  datasetnames <- paste0('all', datasets, 's')
  datasetnames[datasetnames == 'allreactomes'] <- 'allpaths'
  
  geneenrichres <- NULL
  pairenrichres <- NULL
  
  if(!is.null(geneweights)){
    
    if(is.vector(geneweights)){
      genecount <- TRUE
      genes <- geneweights
      geneentrezs <- unlist(lapply(genes, getsinglename, entrez = TRUE))
      geneweights <- data.frame(GENE = geneweights, ENTREZID = geneentrezs, 
                                stringsAsFactors = FALSE)
      geneweights <- subset(geneweights, !is.na(ENTREZID))
      geneentrezs <- geneweights$ENTREZID
      
      if(nrow(geneweights) > 0){
        row.names(geneweights) <- 1:nrow(geneweights)
      }
      
    }else{
      genes <- geneweights$GENE
      geneentrezs <- unlist(lapply(genes, getsinglename, entrez = TRUE))
      geneweights$ENTREZID <- geneentrezs
      geneweights <- subset(geneweights, !is.na(ENTREZID))
      geneentrezs <- geneweights$ENTREZID
      
      if(nrow(geneweights) > 0){
        row.names(geneweights) <- 1:nrow(geneweights)
      }
    }
    
    genemapping <- geneweights[c('GENE', 'ENTREZID')]
    geneweights$GENE <- geneweights$ENTREZID
    
    
    if(!is.null(backweights)){
      
      if(is.vector(backweights)){
        
        genecount <- TRUE
        allgenes <- backweights
        allgeneentrezs <- unlist(lapply(allgenes, getsinglename, entrez = TRUE))
        
        backweights <- data.frame(GENE = backweights, ENTREZID = allgeneentrezs, 
                                  stringsAsFactors = FALSE)
        backweights <- subset(backweights, is.na(ENTREZID))
        allgeneentrezs <- backweights$ENTREZID
        
        if(nrow(backweights) > 0){
          row.names(backweights) <- 1:nrow(backweights)
        }
        
      }else{
        
        allgenes <- backweights$GENE
        allgeneentrezs <- unlist(lapply(allgenes, getsinglename, entrez = TRUE))
        backweights$ENTREZID <- allgeneentrezs
        backweights <- subset(backweights, !is.na(ENTREZID))
        allgeneentrezs <- backweights$ENTREZID
        
        if(nrow(backweights) > 0){
          row.names(backweights) <- 1:nrow(backweights)
        }
      }
      
      backmapping <- backweights[c('GENE', 'ENTREZID')]
      backweights$GENE <- backweights$ENTREZID
      
    }else{
      
      genecount <- TRUE
      
      if(species == 'mouse'){
        #backs <- 
        #  readRDS('C:/Users/yuabr/Desktop/Transfer/codetransfer/emethyl/files/mousebacks.rds')
        backs <- get('mousebacks')
      }else{
        #backs <- 
        #  readRDS('C:/Users/yuabr/Desktop/Transfer/codetransfer/emethyl/files/humanbacks.rds')
        backs <- get('humanbacks')
      }
      
      allgeneentrezs <- unique(c(backs$gos$ENTREZID, 
                                 backs$reactomes$ENTREZID, 
                                 backs$keggs$ENTREZID))
      
    }
    
    genomebackground <- searchtable(allentrez = allgeneentrezs, 
                                    species = 'human')
    
    backgos <- genomebackground
    backgos$mapping <- NULL
    
    genegos <- list()
    
    geneenrichres <- list()
    
    for(i in 1:length(backgos)){
      
      slotname <- names(genomebackground)[i + 1]
      
      if(!(slotname %in% datasetnames)){
        next()
      }
      
      slotdat <- genomebackground[[slotname]]
      
      genegos[[slotname]] <- subset(slotdat, ENTREZID %in% geneentrezs)
      
      
      slotenrichres <- maingo(geneweights = geneweights,
                              backweights = backweights, 
                              backgos = slotdat, 
                              genegos = genegos[[slotname]], 
                              genecount = genecount, 
                              paircount = paircount, 
                              addzeroedgesforshuffle = addzeroedgesforshuffle, 
                              threads = threads, 
                              enrichtype = c('gene'), 
                              seed = seed, 
                              shufflenum = shufflenum, 
                              shuffletwotail = shuffletwotail)
      
      slotgeneres <- slotenrichres$generes
      slotgeneshuffleres <- slotenrichres$geneshuffleres
      slotgenegos <- slotenrichres$geneweights
      
      slotgenegos$GENE <- genemapping$GENE[match(slotgenegos$ENTREZID, 
                                                 genemapping$ENTREZID)]
      
      dataset <- datasets[match(slotname, datasetnames)]
      
      geneenrichres[[dataset]] <- list(generes = slotgeneres, 
                                       geneshuffleres = slotgeneshuffleres, 
                                       genegos = slotgenegos)
      
    }
    
    if(combinedatasetreses == TRUE){
      
      for(j in 1:length(geneenrichres)){
        
        dataset <- names(geneenrichres)[j]
        
        generes <- geneenrichres[[j]]$generes
        geneshuffleres <- geneenrichres[[j]]$geneshuffleres
        genegos <- geneenrichres[[j]]$genegos
        
        if(!is.null(generes)){
          generes$dataset <- dataset
        }
        
        if(!is.null(geneshuffleres)){
          geneshuffleres$dataset <- dataset
        }
        
        if(!is.null(genegos)){
          genegos$dataset <- dataset
        }
        
        
        if(j == 1){
          genereses <- generes
          geneshufflereses <- geneshuffleres
          genegoses <- genegos
        }else{
          genereses <- rbind(genereses, generes)
          geneshufflereses <- rbind(geneshufflereses, geneshuffleres)
          genegoses <- rbind(genegoses, genegos)
        }
        
      }
      
      if(!is.null(genereses)){
        if(nrow(genereses) > 0){
          
          genereses <- unique(genereses)
          genereses <- genereses[order(genereses$padj, genereses$pval, 
                                       -genereses$odds),]
          row.names(genereses) <- 1:nrow(genereses)
        }
      }
      
      if(!is.null(geneshufflereses)){
        
        if(nrow(geneshufflereses) > 0){
          
          geneshufflereses <- unique(geneshufflereses)
          geneshufflereses <- geneshufflereses[order(geneshufflereses$padj, 
                                                     geneshufflereses$pval, 
                                                     -geneshufflereses$weight),]
          row.names(geneshufflereses) <- 1:nrow(geneshufflereses)
        }
        
      }
      
      if(!is.null(genegoses)){
        if(nrow(genegoses) > 0){
          
          genegoses <- unique(genegoses)
          genegoses <- genegoses[order(genegoses$GOID, genegoses$TERM),]
          row.names(genegoses) <- 1:nrow(genegoses)
        }
      }
      
      
      geneenrichres <- list()
      
      geneenrichres$generes <- genereses
      geneenrichres$geneshuffleres <- geneshufflereses
      geneenrichres$genegos <- genegoses
      
    }
    
  }
  
  
  if(!is.null(genepairweights)){
    
    if(!('weight' %in% colnames(genepairweights))){
      paircount <- TRUE
    }else if(mean(genepairweights$weight) == 1 & 
             sd(genepairweights$weight) == 0){
      paircount <- TRUE
    }
    
    fromnodes <- genepairweights$fromNode
    tonodes <- genepairweights$toNode
    
    fromnodeentrezs <- unlist(lapply(fromnodes, getsinglename, entrez = TRUE))
    tonodeentrezs <- unlist(lapply(tonodes, getsinglename, entrez = TRUE))
    
    genepairweights$fromNodeENTREZ <- fromnodeentrezs
    genepairweights$toNodeENTREZ <- tonodeentrezs
    genepairweights <- subset(genepairweights, !is.na(fromNodeENTREZ) & 
                                !is.na(toNodeENTREZ))
    
    fromnodeentrezs <- genepairweights$fromNodeENTREZ
    tonodeentrezs <- genepairweights$toNodeENTREZ
    
    geneentrezs <- unique(c(fromnodeentrezs, tonodeentrezs))
    
    if(!is.null(genepairweights)){
      if(nrow(genepairweights) > 0){
        row.names(genepairweights) <- 1:nrow(genepairweights)
      }
      
    }
    
    
    genepairmapping <- genepairweights[c('fromNode', 'toNode', 
                                         'fromNodeENTREZ', 'toNodeENTREZ')]
    genepairweights$fromNode <- genepairweights$fromNodeENTREZ
    genepairweights$toNode <- genepairweights$toNodeENTREZ
    
    pairenrichgo <- TRUE
    if(!is.null(backpairweights)){
      
      if(!('weight' %in% colnames(backpairweights))){
        paircount <- TRUE
      }else if(mean(backpairweights$weight) == 1 & 
               sd(backpairweights$weight) == 0){
        paircount <- TRUE
      }
      
      allfromnodes <- backpairweights$fromNode
      alltonodes <- backpairweights$toNode
      
      allfromnodeentrezs <- unlist(lapply(allfromnodes, getsinglename, entrez = TRUE))
      alltonodeentrezs <- unlist(lapply(alltonodes, getsinglename, entrez = TRUE))
      
      backpairweights$fromNodeENTREZ <- allfromnodeentrezs
      backpairweights$toNodeENTREZ <- alltonodeentrezs
      
      backpairweights <- subset(backpairweights, !is.na(fromNodeENTREZ) & 
                                  !is.na(toNodeENTREZ))
      
      allfromnodeentrezs <- backpairweights$fromNodeENTREZ
      alltonodeentrezs <- backpairweights$toNodeENTREZ
      
      allgeneentrezs <- unique(c(allfromnodeentrezs, alltonodeentrezs))
      
      if(nrow(backpairweights) > 0){
        row.names(backpairweights) <- 1:nrow(backpairweights)
      }
      
      backpairmapping <- backpairweights[c('fromNode', 'toNode', 
                                           'fromNodeENTREZ', 'toNodeENTREZ')]
      backpairweights$fromNode <- backpairweights$fromNodeENTREZ
      backpairweights$toNode <- backpairweights$toNodeENTREZ
      
    }else{
      
      paircount <- TRUE
      pairenrichgo <- FALSE
      
      if(species == 'mouse'){
        #backs <- 
        #  readRDS('C:/Users/yuabr/Desktop/Transfer/codetransfer/emethyl/files/mousebacks.rds')
        backs <- get('mousebacks')
      }else{
        #backs <- 
        #  readRDS('C:/Users/yuabr/Desktop/Transfer/codetransfer/emethyl/files/humanbacks.rds')
        backs <- get('humanbacks')
      }
      
      allgeneentrezs <- unique(c(backs$gos$ENTREZID, 
                                 backs$reactomes$ENTREZID, 
                                 backs$keggs$ENTREZID))
      
    }
    
    genomebackground <- searchtable(allentrez = allgeneentrezs, 
                                    species = 'human')
    
    backgos <- genomebackground
    backgos$mapping <- NULL
    
    genegos <- list()
    pairenrichres <- list()
    for(i in 1:length(backgos)){
      
      slotname <- names(genomebackground)[i + 1]
      
      if(!(slotname %in% datasetnames)){
        next()
      }
      
      slotdat <- genomebackground[[slotname]]
      
      genegos[[slotname]] <- subset(slotdat, ENTREZID %in% geneentrezs)
      
      slotenrichres <- maingo(backgos = slotdat, 
                              backpairweights = backpairweights, 
                              genegos = genegos[[slotname]], 
                              genepairweights = genepairweights, 
                              shufflegenepairweights = genepairweights, 
                              genecount = genecount, 
                              paircount = paircount, 
                              addzeroedgesforshuffle = addzeroedgesforshuffle, 
                              threads = threads, 
                              enrichtype = c('pair'), 
                              seed = seed, 
                              shufflenum = shufflenum, 
                              shuffletwotail = shuffletwotail, 
                              pairenrichgo = pairenrichgo)
      
      
      slotpairres <- slotenrichres$pairres
      slotpairshuffleres <- slotenrichres$pairshuffleres
      slotpairgos <- slotenrichres$pairweights
      
      slotpairgos$fromNode <- genepairmapping$fromNode[match(slotpairgos$fromNodeENTREZ, 
                                                             genepairmapping$fromNodeENTREZ)]
      
      slotpairgos$toNode <- genepairmapping$toNode[match(slotpairgos$toNodeENTREZ, 
                                                         genepairmapping$toNodeENTREZ)]
      
      
      dataset <- datasets[match(slotname, datasetnames)]
      
      pairenrichres[[dataset]] <- list(pairres = slotpairres, 
                                       pairshuffleres = slotpairshuffleres, 
                                       pairgos = slotpairgos)
      
      
    }
    
    if(combinedatasetreses == TRUE){
      
      for(j in 1:length(pairenrichres)){
        
        dataset <- names(pairenrichres)[j]
        
        pairres <- pairenrichres[[j]]$pairres
        pairshuffleres <- pairenrichres[[j]]$pairshuffleres
        pairgos <- pairenrichres[[j]]$pairgos
        
        if(!is.null(pairres)){
          pairres$dataset <- dataset
        }
        
        if(!is.null(pairshuffleres)){
          pairshuffleres$dataset <- dataset
        }
        
        if(!is.null(pairgos)){
          pairgos$dataset <- dataset
        }
        
        
        if(j == 1){
          pairreses <- pairres
          pairshufflereses <- pairshuffleres
          pairgoses <- pairgos
        }else{
          pairreses <- rbind(pairreses, pairres)
          pairshufflereses <- rbind(pairshufflereses, pairshuffleres)
          pairgoses <- rbind(pairgoses, pairgos)
        }
        
      }
      
      if(!is.null(pairreses)){
        if(nrow(pairreses) > 0){
          
          pairreses <- unique(pairreses)
          pairreses <- pairreses[order(pairreses$padj, pairreses$pval, 
                                       -pairreses$odds),]
          row.names(pairreses) <- 1:nrow(pairreses)
        }
      }
      
      if(!is.null(pairshufflereses)){
        
        if(nrow(pairshufflereses) > 0){
          
          pairshufflereses <- unique(pairshufflereses)
          pairshufflereses <- pairshufflereses[order(pairshufflereses$padj, 
                                                     pairshufflereses$pval, 
                                                     -pairshufflereses$weight),]
          row.names(pairshufflereses) <- 1:nrow(pairshufflereses)
        }
        
      }
      
      if(!is.null(pairgoses)){
        if(nrow(pairgoses) > 0){
          
          pairgoses <- unique(pairgoses)
          pairgoses <- pairgoses[order(pairgoses$GOID, pairgoses$TERM),]
          row.names(pairgoses) <- 1:nrow(pairgoses)
        }
      }
      
      
      pairenrichres <- list()
      
      pairenrichres$pairres <- pairreses
      pairenrichres$pairshuffleres <- pairshufflereses
      pairenrichres$pairgos <- pairgoses
      
    }
    
  }
  
  finalres <- list()
  
  if(length(geneenrichres) > 0 & length(pairenrichres) > 0){
    
    finalres$geneenrichres <- geneenrichres
    finalres$pairenrichres <- pairenrichres
    
  }else if(length(geneenrichres) > 0 & length(pairenrichres) == 0){
    finalres$geneenrichres <- geneenrichres
  }else if(length(geneenrichres) == 0 & length(pairenrichres) > 0){
    finalres$pairenrichres <- pairenrichres
  }
  
  return(finalres)
  
}

#Loci distribution####

probesiteenrich <- function(sub,
                            probenames,
                            dataname,
                            types){
  
  subsub <- subset(sub, Probe %in% probenames)
  subsub <- subsub[subsub[,ncol(subsub)] %in% names(types),]
  subsubtypes <- table(subsub[,ncol(subsub)])
  
  for(i in 1:length(types)){
    singletype <- names(types)[i]
    othertypes <- names(types)[-i]
    
    a11 <- sum(subsubtypes[singletype])
    if(is.na(a11)){
      a11 <- 0
    }
    a12 <- sum(types[singletype])
    
    othertypes1 <- intersect(names(subsubtypes), othertypes)
    a21 <- sum(subsubtypes[othertypes1])
    othertypes2 <- intersect(names(types), othertypes)
    a22 <- sum(types[othertypes2])
    
    fishermat <- matrix(c(a11, a12, a21, a22), byrow = TRUE, nrow = 2)
    fisherp <- fisher.test(fishermat)$p.val
    fisherp <- signif(fisherp, 3)
    
    odds <- (a11/a21)/(a12/a22)
    
    percent <- a11/(a11 + a21)
    if(i == 1){
      fisherps <- fisherp
      oddses <- odds
      percents <- percent
    }else{
      fisherps <- c(fisherps, fisherp)
      oddses <- c(oddses, odds)
      percents <- c(percents, percent)
    }
    
    names(fisherps)[length(fisherps)] <- names(percents)[length(percents)] <-
      names(oddses)[length(oddses)] <- singletype
  }
  
  bardata <- data.frame(regionname = names(percents), percents = percents,
                        fisherps = fisherps, oddses = oddses, stringsAsFactors = FALSE)
  bardata$dataname <- dataname
  
  return(bardata)
}

plotdatorganize <- function(oridat,
                            alldataname){
  
  oridat$dir <- 'NC'
  oridat$dir[oridat$fisherps < 0.05 & oridat$oddses > 1] <- 'UP'
  oridat$dir[oridat$fisherps < 0.05 & oridat$oddses < 1] <- 'DN'
  
  oridat$SS <- -log10(oridat$fisherps)
  oridat$SS[oridat$oddses < 1] <- -oridat$SS[oridat$oddses < 1]
  
  row.names(oridat) <- 1:nrow(oridat)
  
  oridat$color <- 'gray'
  oridat$color[oridat$dir == 'UP'] <- 'blue'
  oridat$color[oridat$dir == 'DN'] <- 'red'
  
  oridat$dir[oridat$dataname == alldataname] <- ''
  
  return(oridat)
  
}

#'Attribute methylation probes to different genomic regions
#'
#'Attribute methylation probes to different genomic regions and show their
#'distribution patterns
#'
#'@param platform The platform of the probes need to be annotated. Can be set
#'  as 450 or 850. Default is 450.
#'@param allprobes The background probes to perform the Fisher's exact test to
#'  check the enrichment of probes in different genomic regions. Should be a
#'  vector with the names of the background probes as elements. Default value
#'  is NULL, and in this case, all the probes provided by another parameter
#'  \code{targetprobelist} will be used as the background probes together.
#'@param targetprobelist A list with each slot recording the probe names of
#'  a specific probe group need to check the genomic region enrichment status
#'  and the names of the slots are the group names. For example, it can be a
#'  list containing 2 slots with one named "Hypermethylated" and the other as
#'  "Hypomethylated", and the probes in each group will be mapped to different
#'  genomic regions and then for each region, its enrichment for the group
#'  probes compared with the background ones will be checked using Fisher's
#'  exact test. For each probe group, a barplot will be generated to show the
#'  probe proportion of each genomic region within that group, and whether the
#'  probes in each genomic region is significantly up or down enriched will be
#'  labeled.
#'@param plotbackground If this parameter is TRUE, a barplot showing the probe
#'  proportion of each genomic region in the background will be generated, but
#'  the enrichment significance cannot be labeled. Default is FALSE.
#'@param removedup Some methylation probes can be mapped to multiple genomic
#'  island regions or gene regions and for these special ones, whether remove
#'  them from the analysis or not. Defult is TRUE, meaning to remove them.
#'@param titlesuffix A suffix string that can be added to the plot title. Its 
#'  default value is NULL.
#'@param titlesize Font size of the plot title. Default is 20.
#'@param textsize Font size of the legend title, legend text, axis label, ect. 
#'  Default is 15.
#'@param annotextsize Font size of the annotation text in the plot. Default is 
#'  4.
#'@param face Font face of the plot.
#'@return Will return a list containing 2 slots. One indicates the statistics
#'  after attributing the probes into different island regions, and the other
#'  indicates after attributing them into different gene regions, including
#'  the percentage of the probes in each region, their Fisher's exact p-value,
#'  etc. Also, the corresponding barplots will be generated.
#'@export
locidistribution <- function(platform = 450,
                             allprobes = NULL,
                             targetprobelist,
                             plotbackground = FALSE,
                             removedup = TRUE,
                             titlesuffix = NULL, 
                             titlesize = 20,
                             textsize = 15,
                             annotextsize = 4,
                             face = NULL){
  
  if(!is.null(titlesuffix)){
    titlesuffix <- paste0(' ', titlesuffix)
  }
  
  allprobeanno <- probeanno(platform = platform, probes = allprobes)
  
  if(removedup == TRUE){
    
    dups <- unique(allprobeanno$Probe[duplicated(allprobeanno$Probe)])
    allprobeanno <- subset(allprobeanno, !(Probe %in% dups))
  }
  
  islandspos <- allprobeanno[c('Probe', 'chr', 'pos', 'strand',
                               'Islands_Name', 'Relation_to_Island')]
  tssinfo <- allprobeanno[c('Probe', 'chr', 'pos', 'strand',
                            'UCSC_RefGene_Name', 'ENTREZID', 'UCSC_RefGene_Group')]
  tssinfo$UCSC_RefGene_Group[tssinfo$UCSC_RefGene_Group == ''] <- 'Intergenic'
  
  islandspostype <- table(islandspos$Relation_to_Island)
  tssinfotype <- table(tssinfo$UCSC_RefGene_Group)
  
  alldataname <- paste0('All ', length(allprobes), ' sites')
  
  
  
  islandall <- probesiteenrich(sub = islandspos,
                               probenames = allprobes,
                               dataname = alldataname,
                               types = islandspostype)
  
  tssinfoall <- probesiteenrich(sub = tssinfo,
                                probenames = allprobes,
                                dataname = alldataname,
                                types = tssinfotype)
  
  if(is.null(names(targetprobelist))){
    
    if(length(targetprobelist) == 1){
      names(targetprobelist) <- 'Targeted sites'
    }else{
      names(targetprobelist) <- paste0('Targeted sites group ', i)
    }
    
  }
  
  i <- 1
  for(i in 1:length(targetprobelist)){
    
    targetdataname <- names(targetprobelist)[i]
    
    targetprobes <- targetprobelist[[i]]
    
    islandtop <- probesiteenrich(sub = islandspos,
                                 probenames = targetprobes,
                                 dataname = targetdataname,
                                 types = islandspostype)
    
    if(i == 1){
      islandbar <- islandtop
    }else{
      islandbar <- rbind(islandbar, islandtop)
    }
    
    
    tssinfotop <- probesiteenrich(sub = tssinfo,
                                  probenames = targetprobes,
                                  dataname = targetdataname,
                                  types = tssinfotype)
    
    
    if(i == 1){
      tssinfobar <- tssinfotop
    }else{
      tssinfobar <- rbind(tssinfobar, tssinfotop)
    }
    
    #print(i)
    
    
  }
  
  if(plotbackground == TRUE){
    
    islandbar <- rbind(islandall, islandbar)
    tssinfobar <- rbind(tssinfoall, tssinfobar)
  }
  
  
  islandbar <- plotdatorganize(oridat = islandbar,
                               alldataname = alldataname)
  tssinfobar <- plotdatorganize(oridat = tssinfobar,
                                alldataname = alldataname)
  
  
  if(plotbackground == TRUE){
    islandbar$dataname <- factor(islandbar$dataname,
                                 levels = c(alldataname, names(targetprobelist)),
                                 ordered = TRUE)
    tssinfobar$dataname <- factor(tssinfobar$dataname,
                                  levels = c(alldataname, names(targetprobelist)),
                                  ordered = TRUE)
    
  }else{
    islandbar$dataname <- factor(islandbar$dataname,
                                 levels = names(targetprobelist),
                                 ordered = TRUE)
    tssinfobar$dataname <- factor(tssinfobar$dataname,
                                  levels = names(targetprobelist),
                                  ordered = TRUE)
  }
  
  #library(ggplot2)
  
  xcoord6 <- c(0.6, 0.75, 0.9, 1.1, 1.25, 1.4)
  xcoord7 <- c(0.6, 0.75, 0.9, 1, 1.1, 1.25, 1.4)
  xcoord8 <- c(0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35)
  
  if(length(unique(islandbar$regionname)) == 6){
    xcoord <- xcoord6
  }else if(length(unique(islandbar$regionname)) == 7){
    xcoord <- xcoord7
  }else if(length(unique(islandbar$regionname)) == 8){
    xcoord <- xcoord8
  }
  
  p <- ggplot2::ggplot(islandbar, ggplot2::aes(x = 1,
                                               y = percents,
                                               fill = regionname))
  
  print(
    
    p + ggplot2::geom_bar(stat = 'identity', position = 'dodge') +
      ggplot2::xlab(NULL) + ggplot2::ylab('Percentage of Sites') +
      ggplot2::ggtitle(paste0('Island region', titlesuffix)) +
      ggplot2::scale_fill_discrete(name = 'Island region') +
      
      ggplot2::facet_wrap(ggplot2::vars(dataname)) +
      
      ggplot2::geom_text(data = islandbar,
                         
                         x = rep(xcoord,
                                 length(unique(islandbar$dataname))),
                         y = max(islandbar$percents)*0.9,
                         ggplot2::aes(label = dir),
                         
                         size = annotextsize,
                         color = islandbar$color,
                         angle = 90,
                         fontface = 'italic') +
      
      ggplot2::theme_bw() +
      
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = textsize, face = face), 
                     axis.title.y = ggplot2::element_text(size = textsize, face = face),
                     legend.title = ggplot2::element_text(size = textsize, face = face),
                     legend.text = ggplot2::element_text(size = textsize, face = face),
                     plot.title = ggplot2::element_text(size = titlesize, face = face),
                     plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                     strip.text.x = ggplot2::element_text(size = textsize, face = face))
    
    
    
  )
  
  if(length(unique(tssinfobar$regionname)) == 6){
    xcoord <- xcoord6
  }else if(length(unique(tssinfobar$regionname)) == 7){
    xcoord <- xcoord7
  }else if(length(unique(tssinfobar$regionname)) == 8){
    xcoord <- xcoord8
  }
  
  p <- ggplot2::ggplot(tssinfobar, ggplot2::aes(x = 1,
                                                y = percents,
                                                fill = regionname))
  
  
  print(
    
    p + ggplot2::geom_bar(stat = 'identity', position = 'dodge') +
      ggplot2::xlab(NULL) + ggplot2::ylab('Percentage of Sites') +
      ggplot2::ggtitle(paste0('Gene region', titlesuffix)) +
      ggplot2::scale_fill_discrete(name = 'Gene region') +
      
      ggplot2::facet_wrap(ggplot2::vars(dataname)) +
      
      ggplot2::geom_text(data = tssinfobar,
                         
                         x = rep(xcoord,
                                 length(unique(tssinfobar$dataname))),
                         y = max(tssinfobar$percents)*0.9,
                         ggplot2::aes(label = dir),
                         
                         size = annotextsize,
                         color = tssinfobar$color,
                         angle = 90,
                         fontface = 'italic') +
      
      ggplot2::theme_bw() +
      
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(size = textsize, face = face), 
                     axis.title.y = ggplot2::element_text(size = textsize, face = face),
                     legend.title = ggplot2::element_text(size = textsize, face = face),
                     legend.text = ggplot2::element_text(size = textsize, face = face),
                     plot.title = ggplot2::element_text(size = titlesize, face = face),
                     plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                     strip.text.x = ggplot2::element_text(size = textsize, face = face))
    
    
    
  )
  
  
  res <- list(islandinfo = islandbar[,1:(ncol(islandbar) - 1)],
              tssinfo = tssinfobar[,1:(ncol(tssinfobar) - 1)])
  
  return(res)
  
}




