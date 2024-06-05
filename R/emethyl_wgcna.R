#blockwiseconsensusmodules####

blockwiseConsensusModules = function(
    multiExpr, 
    
    # Data checking options
    checkMissingData = TRUE,
    
    # Blocking options
    blocks = NULL, 
    maxBlockSize = 5000, 
    blockSizePenaltyPower = 5,
    nPreclusteringCenters = NULL,
    randomSeed = 54321,
    
    # individual TOM information
    
    individualTOMInfo = NULL,
    useIndivTOMSubset = NULL,
    
    # Network construction arguments: correlation options
    corType = "pearson",
    maxPOutliers = 1,
    quickCor = 0,
    pearsonFallback = "individual", 
    cosineCorrelation = FALSE,
    
    # Adjacency function options
    power = 6, 
    networkType = "unsigned", 
    checkPower = TRUE,
    replaceMissingAdjacencies = FALSE,
    
    # Topological overlap options
    
    TOMType = "unsigned",           
    TOMDenom = "min",
    suppressNegativeTOM = FALSE,
    
    # Save individual TOMs?
    
    saveIndividualTOMs = TRUE,
    individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
    
    # Consensus calculation options: network calibration
    
    networkCalibration = c("single quantile", "full quantile", "none"),
    
    ## Save scaled TOMs? <-- leave this option for users willing to run consensusTOM on its own
    #saveScaledIndividualTOMs = FALSE,
    #scaledIndividualTOMFilePattern = "scaledIndividualTOM-Set%s-Block%b.RData",
    
    # Simple quantile calibration options
    
    calibrationQuantile = 0.95,
    sampleForCalibration = TRUE, sampleForCalibrationFactor = 1000,
    getNetworkCalibrationSamples = FALSE, 
    
    # Consensus definition
    
    consensusQuantile = 0,
    useMean = FALSE,
    setWeights = NULL,
    
    # Saving the consensus TOM
    
    saveConsensusTOMs = FALSE, 
    consensusTOMFilePattern = "consensusTOM-block.%b.RData",
    
    # Internal handling of TOMs
    
    useDiskCache = TRUE, chunkSize = NULL,
    cacheBase = ".blockConsModsCache",
    cacheDir = ".",
    
    # Alternative consensus TOM input from a previous calculation 
    
    consensusTOMInfo = NULL,
    
    # Basic tree cut options 
    
    deepSplit = 2, 
    detectCutHeight = 0.995, minModuleSize = 20,
    checkMinModuleSize = TRUE,
    
    # Advanced tree cut opyions
    
    maxCoreScatter = NULL, minGap = NULL,
    maxAbsCoreScatter = NULL, minAbsGap = NULL,
    minSplitHeight = NULL, minAbsSplitHeight = NULL,
    
    useBranchEigennodeDissim = FALSE,
    minBranchEigennodeDissim = mergeCutHeight,
    
    stabilityLabels = NULL,
    minStabilityDissim = NULL,
    
    pamStage = TRUE,  pamRespectsDendro = TRUE,
    
    # Gene joining and removal from a module, and module "significance" criteria
    
    reassignThresholdPS = 1e-4,
    trimmingConsensusQuantile = consensusQuantile,
    # minKMEtoJoin =0.7, 
    minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
    minKMEtoStay = 0.2,
    
    # Module eigengene calculation options
    
    impute = TRUE,
    trapErrors = FALSE,
    
    # Module merging options
    
    equalizeQuantilesForModuleMerging = FALSE,
    quantileSummaryForModuleMerging = "mean",
    mergeCutHeight = 0.15, 
    mergeConsensusQuantile = consensusQuantile,
    
    
    # Output options
    
    numericLabels = FALSE,
    
    # General options
    nThreads = 0,
    verbose = 2, indent = 0,
    ...){
  
  library(WGCNA)
  
  spaces = indentSpaces(indent);
  
  dataSize = checkSets(multiExpr);
  nSets = dataSize$nSets;
  nGenes = dataSize$nGenes;
  
  if (length(power)!=1)
  {
    if (length(power)!=nSets)
      stop("Invalid arguments: Length of 'power' must equal number of sets given in 'multiExpr'.");
  } else {
    power = rep(power, nSets);
  }
  
  seedSaved = FALSE;
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
      seedSaved = TRUE;
      savedSeed = .Random.seed
    } 
    set.seed(randomSeed);
  }
  
  if ( (consensusQuantile < 0) | (consensusQuantile > 1) ) 
    stop("'consensusQuantile' must be between 0 and 1.");
  
  if (checkMinModuleSize & (minModuleSize > nGenes/2))
  {
    minModuleSize = nGenes/2;
    warning("blockwiseConsensusModules: minModuleSize appeared too large and was lowered to", 
            minModuleSize,
            ". If you insist on keeping the original minModuleSize, set checkMinModuleSize = FALSE.");
  }
  
  if (verbose>0) 
    printFlush(paste(spaces, "Calculating consensus modules and module eigengenes", 
                     "block-wise from all genes"));
  
  originalGeneNames = mtd.colnames(multiExpr);
  if (is.null(originalGeneNames)) originalGeneNames = spaste("Column.", 1:nGenes)
  
  originalSampleNames = mtd.apply(multiExpr, function(x)
  {
    out = rownames(x);
    if (is.null(out)) out = spaste("Row.", 1:nrow(x));
    out;
  });
  
  branchSplitFnc = NULL;
  minBranchDissimilarities = numeric(0);
  externalSplitFncNeedsDistance = logical(0);
  if (useBranchEigennodeDissim)
  {
    branchSplitFnc = "mtd.branchEigengeneDissim";
    minBranchDissimilarities = minBranchEigennodeDissim;
    externalSplitFncNeedsDistance = FALSE;
  } 
  
  if (!is.null(stabilityLabels))
  {
    branchSplitFnc = c(branchSplitFnc, "branchSplitFromStabilityLabels");
    minBranchDissimilarities = c(minBranchDissimilarities, minStabilityDissim);
    externalSplitFncNeedsDistance = c(externalSplitFncNeedsDistance, FALSE);
  }
  
  # Basic checks on consensusTOMInfo
  
  if (!is.null(consensusTOMInfo))
  {
    
    .checkComponents(consensusTOMInfo, c("saveConsensusTOMs", "individualTOMInfo", "goodSamplesAndGenes"));
    
    if (length(consensusTOMInfo$individualTOMInfo$blocks)!=nGenes)
      stop("Inconsistent number of genes in 'consensusTOMInfo$individualTOMInfo$blocks'.");
    
    if (!is.null(consensusTOMInfo$consensusQuantile) &&
        (consensusQuantile!=consensusTOMInfo$consensusQuantile) )
      warning(immediate. = TRUE,
              "blockwiseConsensusModules: given (possibly default) 'consensusQuantile' is different\n",
              "from the value recorded in 'consensusTOMInfo'. This is normally undesirable and may\n",
              "indicate a mistake in the function call.");
  }
  
  # Handle "other arguments"
  
  args = list(...);
  if (is.null(args$reproduceBranchEigennodeQuantileError))
  {
    reproduceBranchEigennodeQuantileError = FALSE;
  } else reproduceBranchEigennodeQuantileError = args$reproduceBranchEigennodeQuantileError;
  
  # If topological overlaps weren't calculated yet, calculate them.
  
  removeIndividualTOMsOnExit = FALSE;
  nBlocks.0 = length(unique(blocks));
  
  if (is.null(individualTOMInfo))
  {
    if (is.null(consensusTOMInfo))
    {
      individualTOMInfo = blockwiseIndividualTOMs(multiExpr = multiExpr, 
                                                  checkMissingData = checkMissingData,
                                                  blocks = blocks,
                                                  maxBlockSize = maxBlockSize,
                                                  blockSizePenaltyPower = blockSizePenaltyPower,
                                                  nPreclusteringCenters = nPreclusteringCenters,
                                                  randomSeed = randomSeed,
                                                  corType = corType,
                                                  maxPOutliers = maxPOutliers,
                                                  quickCor = quickCor,
                                                  pearsonFallback = pearsonFallback,
                                                  cosineCorrelation = cosineCorrelation,
                                                  power = power,
                                                  networkType = networkType, 
                                                  replaceMissingAdjacencies= replaceMissingAdjacencies,
                                                  TOMType = TOMType,
                                                  TOMDenom = TOMDenom,
                                                  suppressNegativeTOM = suppressNegativeTOM,
                                                  saveTOMs = useDiskCache | nBlocks.0>1,
                                                  individualTOMFileNames = individualTOMFileNames,
                                                  nThreads = nThreads,
                                                  verbose = verbose, indent = indent);
      removeIndividualTOMsOnExit = TRUE;
    } else
      individualTOMInfo = consensusTOMInfo$individualTOMInfo;
  } 
  
  if (is.null(useIndivTOMSubset))
  {  
    if (individualTOMInfo$nSets != nSets)
      stop(paste("Number of sets in individualTOMInfo and in multiExpr do not agree.\n",
                 "  To use a subset of individualTOMInfo, set useIndivTOMSubset appropriately."));
    
    useIndivTOMSubset = c(1:nSets);
  }
  
  if (length(useIndivTOMSubset)!=nSets)
    stop("Length of 'useIndivTOMSubset' must equal the number of sets in 'multiExpr'");
  
  if (length(unique(useIndivTOMSubset))!=nSets)
    stop("Entries of 'useIndivTOMSubset' must be unique");
  
  if (any(useIndivTOMSubset<1) | any(useIndivTOMSubset>individualTOMInfo$nSets))
    stop("All entries of 'useIndivTOMSubset' must be between 1 and the number of sets in individualTOMInfo");
  
  # if ( (minKMEtoJoin >1) | (minKMEtoJoin  <0) ) stop("minKMEtoJoin  must be between 0 and 1.");
  
  intNetworkType = individualTOMInfo$intNetworkType;
  intCorType = individualTOMInfo$intCorType;
  
  corFnc = match.fun(WGCNA:::.corFnc[intCorType]);
  corOptions = list(use = 'p');
  
  fallback = pmatch(pearsonFallback, WGCNA:::.pearsonFallbacks)
  
  nSamples = dataSize$nSamples;
  
  allLabels = rep(0, nGenes);
  allLabelIndex = NULL;
  
  # Restrict data to goodSamples and goodGenes
  
  gsg = individualTOMInfo$goodSamplesAndGenes;
  
  # Restrict gsg to used sets
  
  gsg$goodSamples = gsg$goodSamples[useIndivTOMSubset];
  if (!gsg$allOK)
    multiExpr = mtd.subset(multiExpr, gsg$goodSamples, gsg$goodGenes);
  
  # prepare scaled and imputed multiExpr.
  multiExpr.scaled = mtd.apply(multiExpr, scale);
  hasMissing = unlist(multiData2list(mtd.apply(multiExpr, function(x) { any(is.na(x)) })));
  # Impute those that have missing data
  multiExpr.scaled.imputed = mtd.mapply(function(x, doImpute) 
  { if (doImpute) t(impute.knn(t(x))$data) else x },
  multiExpr.scaled, hasMissing);
  
  nGGenes = sum(gsg$goodGenes);
  nGSamples = rep(0, nSets);
  for (set in 1:nSets) nGSamples[set] = sum(gsg$goodSamples[[ set ]]);
  
  blocks = individualTOMInfo$blocks;
  gBlocks = individualTOMInfo$gBlocks;
  
  blockLevels = sort(unique(gBlocks));
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels);
  
  if (is.null(chunkSize)) chunkSize = as.integer(WGCNA:::.largestBlockSize/nSets)
  
  reassignThreshold = reassignThresholdPS^nSets;
  
  consMEs = vector(mode = "list", length = nSets);
  dendros = list();
  
  maxUsedLabel = 0;
  collectGarbage();
  # Here's where the analysis starts
  
  removeConsensusTOMOnExit = FALSE;
  if (is.null(consensusTOMInfo) && (nBlocks==1 || saveConsensusTOMs || getNetworkCalibrationSamples))
  {
    consensusTOMInfo = consensusTOM(
      individualTOMInfo = individualTOMInfo,
      useIndivTOMSubset = useIndivTOMSubset,
      
      networkCalibration = networkCalibration,
      saveCalibratedIndividualTOMs = FALSE,
      
      calibrationQuantile = calibrationQuantile,
      sampleForCalibration = sampleForCalibration,
      sampleForCalibrationFactor = sampleForCalibrationFactor,
      getNetworkCalibrationSamples = getNetworkCalibrationSamples,
      
      consensusQuantile = consensusQuantile,
      useMean = useMean,
      setWeights = setWeights,
      
      # Return options
      saveConsensusTOMs = saveConsensusTOMs,
      consensusTOMFilePattern = consensusTOMFilePattern,
      returnTOMs = nBlocks==1,
      
      # Internal handling of TOMs
      useDiskCache = useDiskCache, 
      chunkSize = chunkSize,
      cacheBase = cacheBase,
      cacheDir = cacheDir,
      verbose = verbose, indent = indent);
    removeConsensusTOMOnExit = !saveConsensusTOMs;
  }
  
  blockwiseConsensusCalculation = is.null(consensusTOMInfo);
  
  for (blockNo in 1:nBlocks)
  {
    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."));
    # Select block genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]];
    nBlockGenes = length(block);
    
    selExpr = mtd.subset(multiExpr, , block);
    errorOccurred = FALSE;
    
    if (blockwiseConsensusCalculation)
    {
      # This code is only reached if input saveConsensusTOMs is FALSE and there are at least 2 blocks.
      consensusTOMInfo = consensusTOM(
        individualTOMInfo = individualTOMInfo,
        useIndivTOMSubset = useIndivTOMSubset,
        
        useBlocks = blockNo,
        
        networkCalibration = networkCalibration,
        saveCalibratedIndividualTOMs = FALSE,
        
        calibrationQuantile = calibrationQuantile,
        sampleForCalibration = sampleForCalibration,
        sampleForCalibrationFactor = sampleForCalibrationFactor,
        getNetworkCalibrationSamples = FALSE,
        
        consensusQuantile = consensusQuantile,
        useMean = useMean,
        setWeights = setWeights,
        
        saveConsensusTOMs = FALSE,
        returnTOMs = TRUE,
        
        useDiskCache = useDiskCache, 
        chunkSize = chunkSize,
        cacheBase = cacheBase,
        cacheDir = cacheDir);
      
      consTomDS = consensusTOMInfo$consensusTOM[[1]];
      # Remove the consensus TOM from the structure.
      consensusTOMInfo$consensusTOM[[1]] = NULL;
      consensusTOMInfo$consensusTOM = NULL;
    } else {
      if (consensusTOMInfo$saveConsensusTOMs)
      {
        consTomDS = WGCNA:::.loadObject(file = consensusTOMInfo$TOMFiles[blockNo],
                                size = nBlockGenes * (nBlockGenes-1)/2);
      } else 
        consTomDS = consensusTOMInfo$consensusTOM[[blockNo]];
    }
    
    # Temporary "cast" so fastcluster::hclust doesn't complain about non-integer size.
    
    attr(consTomDS, "Size") = as.integer(attr(consTomDS, "Size"));
    
    consTomDS = 1-consTomDS;
    
    collectGarbage();
    
    if (verbose>2) printFlush(paste(spaces, "....clustering and detecting modules.."));
    errorOccured = FALSE;
    dendros[[blockNo]] = fastcluster::hclust(consTomDS, method = "average");
    if (verbose > 8)
    {
      if (interactive())
        plot(dendros[[blockNo]], labels = FALSE, main = paste("Block", blockNo));
    }
    
    externalSplitOptions = list();
    e.index = 1;
    if (useBranchEigennodeDissim)
    {
      externalSplitOptions[[e.index]] = list(multiExpr = mtd.subset(multiExpr.scaled.imputed,, block),
                                             corFnc = corFnc, corOptions = corOptions,
                                             consensusQuantile = consensusQuantile,
                                             signed = networkType %in% c("signed", "signed hybrid"),
                                             reproduceQuantileError = reproduceBranchEigennodeQuantileError);
      e.index = e.index +1;
    }
    if (!is.null(stabilityLabels))
    {
      externalSplitOptions[[e.index]] = list(stabilityLabels = stabilityLabels);
      e.index = e.index + 1;
    }
    collectGarbage();
    #blockLabels = try(cutreeDynamic(dendro = dendros[[blockNo]], 
    blockLabels = cutreeDynamic(dendro = dendros[[blockNo]], 
                                distM = as.matrix(consTomDS), 
                                deepSplit = deepSplit,
                                cutHeight = detectCutHeight, minClusterSize = minModuleSize, 
                                method ="hybrid", 
                                maxCoreScatter = maxCoreScatter, minGap = minGap,
                                maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                                minSplitHeight = minSplitHeight, minAbsSplitHeight = minAbsSplitHeight,
                                
                                externalBranchSplitFnc = branchSplitFnc, 
                                minExternalSplit = minBranchDissimilarities,
                                externalSplitOptions = externalSplitOptions,
                                externalSplitFncNeedsDistance = externalSplitFncNeedsDistance,
                                assumeSimpleExternalSpecification = FALSE,
                                
                                pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
                                verbose = verbose-3, indent = indent + 2)
    #verbose = verbose-3, indent = indent + 2), silent = TRUE);
    if (verbose > 8)
    {
      print(table(blockLabels));
      if (interactive())
        plotDendroAndColors(dendros[[blockNo]], labels2colors(blockLabels), dendroLabels = FALSE, 
                            main = paste("Block", blockNo));
    }
    if (inherits(blockLabels, 'try-error'))
    {
      (if (verbose>0) printFlush else warning)
      (paste(spaces, "blockwiseConsensusModules: cutreeDynamic failed:\n    ", spaces, 
             blockLabels, "\n", spaces, "    Error occured in block", blockNo, "\n",
             spaces, "   Continuing with next block. "));
      next;
    } else {
      blockLabels[blockLabels>0] = blockLabels[blockLabels>0] + maxUsedLabel;
      maxUsedLabel = max(blockLabels);
    }
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1) 
      {
        printFlush(paste(spaces, "No modules detected in block", blockNo,
                         "--> continuing with next block."))
      }
      next;
    }
    
    # Calculate eigengenes for this batch
    
    if (verbose>2) printFlush(paste(spaces, "....calculating eigengenes.."));
    blockAssigned = c(1:nBlockGenes)[blockLabels!=0];
    blockLabelIndex = as.numeric(levels(as.factor(blockLabels[blockAssigned])));
    blockConsMEs = try(multiSetMEs(selExpr, universalColors = blockLabels,
                                   excludeGrey = TRUE, grey = 0, impute = impute,
                                   # trapErrors = TRUE, returnValidOnly = TRUE, 
                                   verbose = verbose-4, indent = indent + 3), silent = TRUE);
    if (inherits(blockConsMEs, 'try-error'))
    {
      if (verbose>0)
      {
        printFlush(paste(spaces, "*** multiSetMEs failed with the message:"));
        printFlush(paste(spaces, "     ", blockConsMEs));
        printFlush(paste(spaces, "*** --> Ending module detection here"));
      } else warning(paste("blocwiseConsensusModules: multiSetMEs failed with the message: \n",
                           "      ", blockConsMEs, "\n--> continuing with next block."));
      next;
    }
    
    deleteModules = NULL;
    changedModules = NULL;
    
    # find genes whose closest module eigengene has cor higher than minKMEtoJoin and assign them 
    # Removed - should not change blocks before clustering them
    #unassGenes = c(c(1:nGGenes)[-block][allLabels[-block]==0], block[blockLabels==0]);
    #if (length(unassGenes) > 0)
    #{
    #blockKME = array(0, dim = c(length(unassGenes), ncol(blockConsMEs[[1]]$data), nSets));
    #corEval = parse(text = paste(.corFnc[intCorType], 
    #"( multiExpr[[set]]$data[, unassGenes], blockConsMEs[[set]]$data,", 
    #.corOptions[intCorType], ")"))
    #for (set in 1:nSets) blockKME[, , set] = eval(corEval);
    #if (intNetworkType==1) blockKME = abs(blockKME);
    #consKME = as.matrix(apply(blockKME, c(1,2), min));
    #consKMEmax = apply(consKME, 1, max);
    #closestModule = blockLabelIndex[apply(consKME, 1, which.max)];
    #assign = (consKMEmax >= minKMEtoJoin );
    #if (sum(assign>0))
    #{
    #allLabels[unassGenes[assign]] = closestModule[assign]; 
    #changedModules = union(changedModules, closestModule[assign]);
    #}
    #rm(blockKME, consKME, consKMEmax);
    #}
    
    collectGarbage();
    
    # Check modules: make sure that of the genes present in the module, at least a minimum number
    # have a correlation with the eigengene higher than a given cutoff, and that all member genes have
    # the required minimum consensus KME
    
    if (verbose>2) 
      printFlush(paste(spaces, "....checking consensus modules for statistical meaningfulness.."));
    
    for (mod in 1:ncol(blockConsMEs[[1]]$data))
    {
      modGenes = (blockLabels==blockLabelIndex[mod]);
      KME = matrix(0, nrow = sum(modGenes), ncol = nSets);
      corEval = parse(text = paste(WGCNA:::.corFnc[intCorType], 
                                   "( selExpr[[set]]$data[, modGenes], blockConsMEs[[set]]$data[, mod]", 
                                   prepComma(WGCNA:::.corOptions[intCorType]), ")"))
      for (set in 1:nSets) KME[, set] = as.vector(eval(corEval));
      if (intNetworkType==1) KME = abs(KME);
      consKME = apply(KME, 1, quantile, probs = trimmingConsensusQuantile, names = FALSE, na.rm = TRUE);
      if (sum(consKME>minCoreKME) < minCoreKMESize) 
      {
        blockLabels[modGenes] = 0;
        deleteModules = union(deleteModules, mod);
        if (verbose>3) 
          printFlush(paste(spaces, "    ..deleting module ",blockLabelIndex[mod], 
                           ": of ", sum(modGenes), 
                           " total genes in the module only ",  sum(consKME>minCoreKME), 
                           " have the requisite high correlation with the eigengene in all sets.", sep=""));
      } else if (sum(consKME<minKMEtoStay)>0)
      {
        if (verbose > 3) 
          printFlush(paste(spaces, "    ..removing", sum(consKME<minKMEtoStay),
                           "genes from module", blockLabelIndex[mod], "because their KME is too low."));
        blockLabels[modGenes][consKME < minKMEtoStay] = 0;
        if (sum(blockLabels[modGenes]>0) < minModuleSize) 
        {
          deleteModules = union(deleteModules, mod);
          blockLabels[modGenes] = 0;
          if (verbose>3) 
            printFlush(paste(spaces, "    ..deleting module ",blockLabelIndex[mod], 
                             ": not enough genes in the module after removal of low KME genes.", sep=""));
        } else {
          changedModules = union(changedModules, blockLabelIndex[mod]);
        }
      }
    }
    
    # Remove marked modules
    
    if (!is.null(deleteModules)) 
    {
      for (set in 1:nSets) blockConsMEs[[set]]$data = blockConsMEs[[set]]$data[, -deleteModules];
      modGenes = is.finite(match(blockLabels, blockLabelIndex[deleteModules]));
      blockLabels[modGenes] = 0;
      modAllGenes = is.finite(match(allLabels, blockLabelIndex[deleteModules]));
      allLabels[modAllGenes] = 0;
      blockLabelIndex = blockLabelIndex[-deleteModules];
    }
    
    # Check whether there's anything left
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1) 
      {
        printFlush(paste(spaces, "  ..No significant modules detected in block", blockNo))
        printFlush(paste(spaces, "  ..continuing with next block."));
      }
      next;
    }
    
    # Update module eigengenes
    
    for (set in 1:nSets) 
      if (is.null(dim(blockConsMEs[[set]]$data))) 
        dim(blockConsMEs[[set]]$data) = c(length(blockConsMEs[[set]]$data), 1);
    
    if (is.null(consMEs[[1]]))
    {
      for (set in 1:nSets) consMEs[[set]] = list(data = blockConsMEs[[set]]$data);
    } else for (set in 1:nSets)
      consMEs[[set]]$data = cbind(consMEs[[set]]$data, blockConsMEs[[set]]$data);
    
    # Update allLabels
    
    allLabelIndex = c(allLabelIndex, blockLabelIndex);
    allLabels[gsg$goodGenes][block[blockAssigned]] = blockLabels[blockAssigned];
    
    collectGarbage();
    
  }
  
  # Check whether any of the already assigned genes (in this or previous blocks) should be re-assigned
  
  if (verbose>2) 
    printFlush(paste(spaces, "....checking for genes that should be reassigned.."));
  
  deleteModules = NULL;
  goodLabels = allLabels[gsg$goodGenes];
  if (sum(goodLabels!=0) > 0)
  {
    propLabels = goodLabels[goodLabels!=0];
    assGenes = (c(1:nGenes)[gsg$goodGenes])[goodLabels!=0];
    corEval = parse(text = paste(WGCNA:::.corFnc[intCorType], 
                                 "(multiExpr[[set]]$data[, goodLabels!=0], consMEs[[set]]$data",
                                 prepComma(WGCNA:::.corOptions[intCorType]), ")"));
    nMods = ncol(consMEs[[1]]$data);
    lpValues = array(0, dim = c(length(propLabels), nMods, nSets));
    sumSign = array(0, dim = c(length(propLabels), nMods));
    if (verbose>3) 
      printFlush(paste(spaces, "......module membership p-values.."));
    for (set in 1:nSets) 
    {
      KME = eval(corEval);
      if (intNetworkType == 1) KME = abs(KME)
      lpValues[,,set] = -2*log(corPvalueFisher(KME, nGSamples[set], twoSided = FALSE));
      sumSign = sumSign + sign(KME);
    }
    if (verbose>3) 
      printFlush(paste(spaces, "......module membership scores.."));
    scoreAll = as.matrix(apply(lpValues, c(1,2), sum)) * (nSets + sumSign)/(2*nSets);
    scoreAll[!is.finite(scoreAll)] = 0.001 # This low should be enough
    bestScore = apply(scoreAll, 1, max);
    if (intNetworkType==1) sumSign = abs(sumSign);
    if (verbose>3) 
    {
      cat(paste(spaces, "......individual modules.."));
      pind = initProgInd();
    }
    for (mod in 1:nMods)
    {
      modGenes = c(1:length(propLabels))[propLabels==allLabelIndex[mod]];
      scoreMod = scoreAll[modGenes, mod];
      candidates = (bestScore[modGenes] > scoreMod);
      candidates[!is.finite(candidates)] = FALSE;
      if (sum(candidates) > 0)
      {
        pModule = pchisq(scoreMod[candidates], nSets, log.p = TRUE)
        whichBest = apply(scoreAll[modGenes[candidates], , drop = FALSE], 1, which.max);
        pBest = pchisq(bestScore[modGenes[candidates]], nSets, log.p = TRUE);
        reassign =  ifelse(is.finite(pBest - pModule), 
                           ( (pBest - pModule) < log(reassignThreshold) ), 
                           FALSE);
        if (sum(reassign)>0)
        {
          allLabels[assGenes[modGenes[candidates][reassign]]] = whichBest[reassign];
          changedModules = union(changedModules, whichBest[reassign]);
          if (length(modGenes)-sum(reassign) < minModuleSize)
          {
            deleteModules = union(deleteModules, mod);
          } else
            changedModules = union(changedModules, mod);
        }
      }
      if (verbose > 3) pind = updateProgInd(mod/nMods, pind);
    }
    rm(lpValues, sumSign, scoreAll);
    if (verbose > 3) printFlush("");
  }
  
  # Remove marked modules
  
  if (!is.null(deleteModules)) 
  {
    # for (set in 1:nSets) consMEs[[set]]$data = consMEs[[set]]$data[, -deleteModules];
    modGenes = is.finite(match(allLabels, allLabelIndex[deleteModules]));
    allLabels[modGenes] = 0;
    # allLabelIndex = allLabelIndex[-deleteModules];
  }
  
  if (verbose>1) printFlush(paste(spaces, "..merging consensus modules that are too close.."));
  #print(table(allLabels));
  #print(is.numeric(allLabels))
  if (numericLabels) {
    colors = allLabels
  } else {
    colors = labels2colors(allLabels)
  }
  mergedColors = colors;
  mergedMods = try(mergeCloseModules(multiExpr, colors[gsg$goodGenes], 
                                     equalizeQuantiles = equalizeQuantilesForModuleMerging,
                                     quantileSummary = quantileSummaryForModuleMerging,
                                     consensusQuantile = mergeConsensusQuantile,
                                     cutHeight = mergeCutHeight, 
                                     relabel = TRUE, impute = impute,
                                     verbose = verbose-2, indent = indent + 2), silent = TRUE);
  if (inherits(mergedMods, 'try-error'))
  {
    if (verbose>0) 
    {
      printFlush(paste(spaces, 'blockwiseConsensusModules: mergeCloseModule failed with this message:\n',
                       spaces, '    ', mergedMods, spaces,
                       '---> returning unmerged consensus modules'));
    } else warning(paste('blockwiseConsensusModules: mergeCloseModule failed with this message:\n     ',
                         mergedMods, '---> returning unmerged consensus modules'));
    MEs = try(multiSetMEs(multiExpr, universalColors = colors[gsg$goodGenes]
                          # trapErrors = TRUE, returnValidOnly = TRUE
    ), silent = TRUE);
    if (inherits(MEs, 'try-error'))
    {
      warning(paste('blockwiseConsensusModules: ME calculation failed with this message:\n     ',
                    MEs, '---> returning empty module eigengenes'));
      allSampleMEs = NULL;
    } else {
      if (!MEs[[1]]$allOK) mergedColors[gsg$goodGenes] = MEs[[1]]$validColors;
      allSampleMEs = vector(mode = "list", length = nSets);
      names(allSampleMEs) = names(multiExpr);
      for (set in 1:nSets)
      {
        allSampleMEs[[set]] =
          list(data = as.data.frame(matrix(NA, nrow = nSamples[set], ncol = ncol(MEs[[set]]$data))));
        allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = MEs[[set]]$data[,];
        names(allSampleMEs[[set]]$data) = names(MEs[[set]]$data);
        rownames(allSampleMEs[[set]]$data) = make.unique(originalSampleNames[[set]]$data)
      }
    }
  } else {
    mergedColors[gsg$goodGenes] = mergedMods$colors;
    allSampleMEs = vector(mode = "list", length = nSets);
    names(allSampleMEs) = names(multiExpr);
    for (set in 1:nSets)
    {
      allSampleMEs[[set]] = 
        list(data = as.data.frame(matrix(NA, nrow = nSamples[set], 
                                         ncol = ncol(mergedMods$newMEs[[1]]$data))));
      allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = mergedMods$newMEs[[set]]$data[,];
      names(allSampleMEs[[set]]$data) = names(mergedMods$newMEs[[set]]$data);
      rownames(allSampleMEs[[set]]$data) = make.unique(originalSampleNames[[set]]$data);
    }
  }
  
  if (seedSaved) .Random.seed <<- savedSeed;
  
  if (removeConsensusTOMOnExit) 
  {
    WGCNA:::.checkAndDelete(consensusTOMInfo$TOMFiles);
    consensusTOMInfo$TOMFiles = NULL;
  }
  
  if (removeIndividualTOMsOnExit)
  {
    WGCNA:::.checkAndDelete(individualTOMInfo$actualTOMFileNames);
    individualTOMInfo$actualTOMFileNames = NULL;
  }
  
  # Under no circumstances return consensus TOM or individual TOM similarities within the returned list.
  consensusTOMInfo$consensusTOM = NULL;
  individualTOMInfo$TOMSimilarities = NULL
  
  
  
  detach('package:WGCNA', unload = TRUE)
  detach('package:fastcluster', unload = TRUE)
  detach('package:dynamicTreeCut', unload = TRUE)
  
  
  
  names(mergedColors) = names(colors) = originalGeneNames;
  list(colors = mergedColors, 
       unmergedColors = colors,
       multiMEs = allSampleMEs, 
       goodSamples = gsg$goodSamples, 
       goodGenes = gsg$goodGenes, 
       dendrograms = dendros,
       TOMFiles = consensusTOMInfo$TOMFiles,
       blockGenes = individualTOMInfo$blockGenes,
       blocks = blocks,
       originCount = consensusTOMInfo$originCount,
       networkCalibrationSamples = consensusTOMInfo$networkCalibrationSamples, 
       individualTOMInfo = individualTOMInfo,
       consensusTOMInfo = if (saveConsensusTOMs) consensusTOMInfo else NULL,
       consensusQuantile = consensusQuantile
  )
}

#WGCNA#####

covarpcs <- function(covardat){
  
  pcs <- prcomp(covardat, scale. = TRUE)
  
  cumvar <- cumsum(pcs$sdev^2/sum(pcs$sdev^2))
  pcnum <- which(x = cumvar > 0.8)[1]
  pcnum <- min(pcnum, 5)
  
  toppcs <- pcs$x[, c(1:pcnum), drop = FALSE]
  
  #Note, The signs of the columns of the rotation matrix are
  #arbitrary, and so may differ between different programs for
  #PCA, and even between different builds of R.
  #Hence, need to check the correlation between pc1 and the
  #row means of expsub, if it is less than 0, it means the
  #sign of pc1 is reversed, and need to reverse it back
  pc1cor <- cor(rowMeans(apply(covardat, 2, scale)), toppcs[,1])
  if(pc1cor < 0){
    toppcs <- -toppcs
  }
  
  return(toppcs)
  
}

#dat is a matrix using features as rows and samples as columns
covadj <- function(dat,
                   pddat, 
                   confoundings = NULL){
  
  if(is.null(confoundings)){
    return(dat)
  }
  
  confoundings <- colnames(pddat)[colnames(pddat) %in% confoundings]
  
  simpddat <- pddat[, confoundings, drop = FALSE]
  row.names(simpddat) <- pddat$sampleid
  
  for(i in 1:ncol(simpddat)){
    
    if(!is.numeric(simpddat[,i])){
      
      if(!is.factor(simpddat[,i])){
        
        simpddat[,i] <- factor(simpddat[,i])
        
      }
      
      simpddat[,i] <- as.numeric(simpddat[,i])
      
      if(length(unique(simpddat[,i])) > 1){
        simpddat[,i] <- simpddat[,i] - 1
      }
    }
    
  }
  
  simpddat <- simpddat[complete.cases(simpddat),]
  
  pcdat <- covarpcs(covardat = simpddat)
  pcnames <- colnames(pcdat)
  row.names(pcdat) <- row.names(simpddat)
  pcdat <- as.data.frame(pcdat, stringsAsFactors = FALSE)
  
  vars <- paste(pcnames, collapse = ' + ')
  
  design <- model.matrix(as.formula(paste0('~ ', vars)), 
                         data = pcdat)
  
  betasub <- dat[,row.names(simpddat)]
  
  fit1 <- limma::lmFit(betasub, design)
  fit1 <- limma::eBayes(fit1)
  
  res <- limma::residuals.MArrayLM(object = fit1, y = betasub)
  res <- apply(X = res, MARGIN = 1, FUN = scale)
  res <- t(res)
  colnames(res) <- colnames(betasub)
  
  return(res)
  
}

#data in datlist is a matrix using features as rows and sampels as columns
choosepower <- function(datlist, 
                        i = 1,
                        powers = seq(1, 20, 1), 
                        rsqcutline = 0.8, 
                        plot = FALSE, 
                        titlesize, 
                        textsize, 
                        face, 
                        signtype = 'unsigned'){
  
  dat <- datlist[[i]]
  
  dat <- t(dat)
  
  sft <- WGCNA::pickSoftThreshold(data = dat, 
                                  powerVector = powers, 
                                  RsquaredCut = rsqcutline, 
                                  verbose = 5, 
                                  networkType = signtype)
  
  if(is.na(sft$powerEstimate)){
    
    if(max(sft$fitIndices$SFT.R.sq) > 0.5){
      sft$powerEstimate <- sft$fitIndices$Power[match(max(sft$fitIndices$SFT.R.sq), 
                                                      sft$fitIndices$SFT.R.sq)]
      
      warning(paste0('No power with a soft R square > ', rsqcutline, 
                     ', use soft R square = ', max(sft$fitIndices$SFT.R.sq), ' instead\n'))
    }else{
      return(NULL)
    }
    
  }
  
  if(plot == TRUE & length(datlist) == 1){
    
    sftdat <- data.frame(powers = rep(sft$fitIndices[,1], 2), 
                         vals = c(-sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
                                  sft$fitIndices[,5]), 
                         metrics = rep(c('Signed R Square', 'Mean Connectivity'), 
                                       each = nrow(sft$fitIndices)), 
                         stringsAsFactors = FALSE)
    
    sftdat$metrics <- factor(x = sftdat$metrics, levels = c('Signed R Square', 
                                                            'Mean Connectivity'), 
                             ordered = TRUE)
    
    p <- ggplot2::ggplot(sftdat, ggplot2::aes(x = powers, y = vals))
    
    print(
      
      p + ggplot2::geom_point(size = 2, color = scales::hue_pal()(1)) + 
        ggplot2::geom_line(size = 1, color = scales::hue_pal()(1)) + 
        ggplot2::xlab('Soft Threshold (Power)') + ggplot2::ylab('Value') + 
        ggplot2::ggtitle('Scale free topology for soft-thresholding') + 
        ggplot2::geom_vline(xintercept = sft$powerEstimate, color = 'blue', size = 1.5) + 
        ggplot2::facet_wrap(ggplot2::vars(metrics), scales = 'free_y', ncol = 2) + 
        ggplot2::theme_bw() + 
        ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face), 
                       legend.title = ggplot2::element_text(size = textsize, face = face), 
                       legend.text = ggplot2::element_text(size = textsize, face = face), 
                       axis.text = ggplot2::element_text(size = textsize, face = face), 
                       axis.title = ggplot2::element_text(size = textsize, face = face), 
                       strip.text = ggplot2::element_text(size = textsize, face = face))
      
    )
    
  }
  
  return(sft)
  
}

choosecommonsft <- function(sftreslist, 
                            rsqcutline = 0.8, 
                            mergesft = TRUE){
  
  powerestimates <- lapply(X = sftreslist, FUN = function(x){x[[1]]})
  fitidiceses <- lapply(X = sftreslist, FUN = function(x){x[[2]]})
  
  powerestimates <- unlist(powerestimates)
  
  if(length(powerestimates) < length(fitidiceses)){
    return(NULL)
  }
  
  if(mergesft == FALSE){
    res <- list(power = powerestimates)
    return(res)
  }
  
  unipowerestimates <- unique(powerestimates)
  
  if(length(unipowerestimates) == 1){
    res <- list(power = unipowerestimates)
    return(res)
  }
  
  qualipowerestimates <- lapply(X = fitidiceses, FUN = function(x){x$Power[
    x$SFT.R.sq > rsqcutline]})
  qualipowerestimates <- unlist(qualipowerestimates)
  qualipowerestimates <- as.numeric(names(table(qualipowerestimates))[
    table(qualipowerestimates) == length(fitidiceses)])
  
  if(length(qualipowerestimates) >= 1){
    res <- list(power = min(qualipowerestimates))
    return(res)
  }else{
    res <- list(power = powerestimates)
    return(res)
  }
  
}

esftpower <- function(dats, 
                      powers = c(1, 20, 1), 
                      rsqcutline = 0.8, 
                      threads = 1, 
                      mergesft = TRUE, 
                      signtype = 'unsigned'){
  
  choosepower <- function(datlist, 
                          i = 1,
                          powers = seq(1, 20, 1), 
                          rsqcutline = 0.8, 
                          plot = FALSE, 
                          titlesize, 
                          textsize, 
                          face, 
                          signtype = 'unsigned'){
    
    dat <- datlist[[i]]
    
    dat <- t(dat)
    
    sft <- WGCNA::pickSoftThreshold(data = dat, 
                                    powerVector = powers, 
                                    RsquaredCut = rsqcutline, 
                                    verbose = 5, 
                                    networkType = signtype)
    
    if(is.na(sft$powerEstimate)){
      
      if(max(sft$fitIndices$SFT.R.sq) > 0.5){
        sft$powerEstimate <- sft$fitIndices$Power[match(max(sft$fitIndices$SFT.R.sq), 
                                                        sft$fitIndices$SFT.R.sq)]
        
        warning(paste0('No power with a soft R square > ', rsqcutline, 
                       ', use soft R square = ', max(sft$fitIndices$SFT.R.sq), ' instead\n'))
      }else{
        return(NULL)
      }
      
    }
    
    if(plot == TRUE){
      
      sftdat <- data.frame(powers = rep(sft$fitIndices[,1], 2), 
                           vals = c(-sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
                                    sft$fitIndices[,5]), 
                           metrics = rep(c('Signed R Square', 'Mean Connectivity'), 
                                         each = nrow(sft$fitIndices)), 
                           stringsAsFactors = FALSE)
      
      sftdat$metrics <- factor(x = sftdat$metrics, levels = c('Signed R Square', 
                                                              'Mean Connectivity'), 
                               ordered = TRUE)
      
      p <- ggplot2::ggplot(sftdat, ggplot2::aes(x = powers, y = vals))
      
      print(
        
        p + ggplot2::geom_point(size = 2, color = scales::hue_pal()(1)) + 
          ggplot2::geom_line(size = 1, color = scales::hue_pal()(1)) + 
          ggplot2::xlab('Soft Threshold (Power)') + ggplot2::ylab('Value') + 
          ggplot2::ggtitle('Scale free topology for soft-thresholding') + 
          ggplot2::geom_vline(xintercept = sft$powerEstimate, color = 'blue', size = 1.5) + 
          ggplot2::facet_wrap(ggplot2::vars(metrics), scales = 'free_y', ncol = 2) + 
          ggplot2::theme_bw() + 
          ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face), 
                         legend.title = ggplot2::element_text(size = textsize, face = face), 
                         legend.text = ggplot2::element_text(size = textsize, face = face), 
                         axis.text = ggplot2::element_text(size = textsize, face = face), 
                         axis.title = ggplot2::element_text(size = textsize, face = face), 
                         strip.text = ggplot2::element_text(size = textsize, face = face))
        
      )
      
    }
    
    return(sft)
    
  }
  
  if(threads == 1){
    
    sftreslist <- list()
    
    for(i in 1:length(dats)){
      
      sftres <- choosepower(datlist = dats, 
                            i = i, 
                            powers = powers, 
                            rsqcutline = rsqcutline, 
                            signtype = signtype)
      
      sftreslist[[i]] <- sftres
      
    }
    
  }else{
    
    
    iseqs <- seq(1, length(dats), 1)
    
    #library(doParallel)
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    #date()
    `%dopar%` <- foreach::`%dopar%`
    sftreslist <- foreach::foreach(i = iseqs,
                                   #.export = ls(name = globalenv())) %dopar% {
                                   .export = NULL) %dopar% {
                                     choosepower(datlist = dats, 
                                                 i, 
                                                 powers = powers, 
                                                 rsqcutline = rsqcutline, 
                                                 signtype = signtype)
                                   }
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
    
    
  }
  
  commonsft <- choosecommonsft(sftreslist = sftreslist, rsqcutline = rsqcutline, 
                               mergesft = mergesft)
  
  commonsft$power <- commonsft
  commonsft$sftpowers <- sftreslist
  
  return(commonsft)
  
}

maphuecolors <- function(numericvec){
  
  uninumericvec <- unique(numericvec)
  uninumericvec <- uninumericvec[order(uninumericvec)]
  uninumericvec <- unique(c(0, uninumericvec))
  unicolorvec <- c('#C0C0C0', scales::hue_pal()(length(uninumericvec) - 1))
  names(unicolorvec) <- as.character(uninumericvec)
  colorvec <- unicolorvec[as.character(numericvec)]
  colorvec <- as.vector(colorvec)
  return(colorvec)
  
}

wgcnamain <- function(dat, 
                      sftpower, 
                      minclustersize = 50, 
                      mergecutheight = 0.2, 
                      maxblocksize = 5000, 
                      returntom = FALSE, 
                      plot = FALSE, 
                      featuretype = 'gene', 
                      plottitleprefix = NULL, 
                      threads = 0, 
                      seed = 2022, 
                      
                      signtype = 'unsigned'){
  
  if(!is.numeric(mergecutheight)){
    mergecutheight <- WGCNA::dynamicMergeCut(n = ncol(dats[[1]]), mergeCor = 0.8)
  }
  
  cor <- WGCNA::cor
  
  stamp <- Sys.time()
  stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
  stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
  stamp <- paste0('blockwiseTOM.', stamp)
  
  wgcnares <-  WGCNA::blockwiseModules(datExpr = t(dat), 
                                       power = as.vector(unlist(sftpower)),
                                       TOMType = signtype, 
                                       networkType = signtype, 
                                       minModuleSize = minclustersize,
                                       deepSplit = TRUE, 
                                       reassignThreshold = 0, 
                                       mergeCutHeight = mergecutheight,
                                       numericLabels = TRUE, 
                                       pamRespectsDendro = FALSE,
                                       saveTOMs = returntom, 
                                       saveTOMFileBase = stamp, 
                                       verbose = 3, 
                                       nThreads = threads, 
                                       maxBlockSize = maxblocksize, 
                                       randomSeed = seed)
  
  if(returntom == TRUE){
    
    if(length(wgcnares$TOMFiles) == 1){
      tomfile <- wgcnares$TOMFiles
      env2load <- environment()
      loadtom <- load(file = tomfile, envir = env2load)
      tom <- get(x = loadtom, envir = env2load)
      
      tom <- as.matrix(tom)
      colnames(tom) <- row.names(tom) <- 
        names(wgcnares$colors)[wgcnares$blockGenes[[1]]]
      
      unlink(tomfile)
      
    }else{
      tom <- list()
      i <- 1
      for(i in 1:length(wgcnares$TOMFiles)){
        tomfile <- wgcnares$TOMFiles[i]
        env2load <- environment()
        loadtom <- load(file = tomfile, envir = env2load)
        tom[[i]] <- get(x = loadtom, envir = env2load)
        
        tom[[i]] <- as.matrix(tom[[i]])
        colnames(tom[[i]]) <- row.names(tom[[i]]) <- 
          names(wgcnares$colors)[wgcnares$blockGenes[[i]]]
        
        unlink(tomfile)
      }
    }
    
  }
  
  
  MEs <- wgcnares$MEs
  numericlabels <- wgcnares$colors
  
  if(plot == TRUE){
    
    
    #colorlabels <- WGCNA::labels2colors(labels = numericlabels)
    colorlabels <- maphuecolors(numericvec = numericlabels)
    names(colorlabels) <- names(numericlabels)
    
    if(is.null(plottitleprefix)){
      
      captialstr <- paste0(toupper(substr(featuretype, 1, 1)), 
                           substr(featuretype, 2, nchar(featuretype)))
      
      plottitleprefix <- paste0(captialstr, ' dendrogram and module colors')
    }else{
      plottitleprefix <- paste0(plottitleprefix, ' ', featuretype, 
                                ' dendrogram and module colors')
    }
    
    
    if(length(wgcnares$dendrograms) == 1){
      
      print(
        
        WGCNA::plotDendroAndColors(wgcnares$dendrograms[[1]], 
                                   colorlabels, 
                                   groupLabels = 'Module colors', 
                                   dendroLabels = FALSE, 
                                   hang = 0.03,
                                   addGuide = TRUE, 
                                   guideHang = 0.05, 
                                   main = plottitleprefix)
        
      )
    }else{
      
      # Plot the dendrogram and the module colors underneath for each block
      for(j in 1:length(wgcnares$dendrograms)){
        
        print(
          
          WGCNA::plotDendroAndColors(wgcnares$dendrograms[[j]], 
                                     colorlabels[wgcnares$blockGenes[[j]]], 
                                     groupLabels = 'Module colors', 
                                     dendroLabels = FALSE, 
                                     hang = 0.03,
                                     addGuide = TRUE, 
                                     guideHang = 0.05,
                                     main = paste0(plottitleprefix, ' in block ', j))
          
        )
      }
    }
    
    
  }
  
  res <- list(mergedmelist = MEs, numericlabels = numericlabels)
  if(returntom == TRUE){
    res$tom <- tom
  }
  if(plot == TRUE){
    res$colorlabels <- colorlabels
  }
  
  return(res)
  
}

createmultiexp <- function(dats){
  
  multiexpr <- vector(mode = "list", length = length(dats))
  
  i <- 1
  for(i in 1:length(dats)){
    multiexpr[[i]] <- list(data = as.data.frame(t(dats[[i]])))
    names(multiexpr[[i]]$data) <- row.names(dats[[i]])
    rownames(multiexpr[[i]]$data) <- colnames(dats[[i]])
    
  }
  
  return(multiexpr)
  
}

ewgcna <- function(dats = balanceddats$dats, 
                   powers = sftpowers$power, 
                   minclustersize = 50, 
                   mergecutheight = 0.2, 
                   maxblocksize = 5000, 
                   returnconsensustom = FALSE, 
                   plot = FALSE, 
                   featuretype = 'gene', 
                   plottitleprefix = NULL, 
                   threads = 0, 
                   seed = 2022, 
                   
                   signtype = 'unsigned'){
  
  if(length(as.vector(unlist(powers))) == 1){
    powers <- rep(powers, length(dats))
  }
  
  if(!is.numeric(mergecutheight)){
    mergecutheight <- WGCNA::dynamicMergeCut(n = ncol(dats[[1]]), mergeCor = 0.8)
  }
  
  cor <- WGCNA::cor
  
  stamp <- Sys.time()
  stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
  stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
  consensusstamp <- paste0("consensusTOM-block.%b.", stamp, ".RData")
  
  multidats <- createmultiexp(dats = dats)
  
  
  wgcnareslist <- blockwiseConsensusModules(multiExpr = multidats, 
                                            power = powers, 
                                            minModuleSize = minclustersize, 
                                            deepSplit = TRUE, 
                                            pamRespectsDendro = FALSE, 
                                            mergeCutHeight = mergecutheight, 
                                            numericLabels = TRUE, 
                                            minKMEtoStay = 0, 
                                            saveConsensusTOMs = returnconsensustom, 
                                            consensusTOMFilePattern = consensusstamp, 
                                            verbose = 5, 
                                            nThreads = threads, 
                                            maxBlockSize = maxblocksize, 
                                            randomSeed = seed, 
                                            
                                            networkType = signtype, 
                                            TOMType = signtype)
  
  #wgcnareslist <- readRDS('wgcnareslist.rds')
  
  if(returnconsensustom == TRUE){
    
    if(length(wgcnareslist$consensusTOMInfo$TOMFiles) == 1){
      consensustomfile <- wgcnareslist$consensusTOMInfo$TOMFiles
      env2load <- environment()
      loadconsensustom <- load(file = consensustomfile, envir = env2load)
      consensustom <- get(x = loadconsensustom, envir = env2load)
      
      consensustom <- as.matrix(consensustom)
      colnames(consensustom) <- row.names(consensustom) <- 
        names(wgcnareslist$colors)[wgcnareslist$blockGenes[[1]]]
      
      unlink(consensustomfile)
      
    }else{
      consensustom <- list()
      i <- 1
      for(i in 1:length(wgcnareslist$consensusTOMInfo$TOMFiles)){
        consensustomfile <- wgcnareslist$consensusTOMInfo$TOMFiles[i]
        env2load <- environment()
        loadconsensustom <- load(file = consensustomfile, envir = env2load)
        consensustom[[i]] <- get(x = loadconsensustom, envir = env2load)
        
        consensustom[[i]] <- as.matrix(consensustom[[i]])
        colnames(consensustom[[i]]) <- row.names(consensustom[[i]]) <- 
          names(wgcnareslist$colors)[wgcnareslist$blockGenes[[i]]]
        
        unlink(consensustomfile)
      }
    }
    
  }
  
  
  MElist <- list()
  for(k in 1:length(wgcnareslist$multiMEs)){
    MElist[[k]] <- wgcnareslist$multiMEs[[k]]$data
    names(MElist)[k] <- paste0('round_', k)
  }
  numericlabels <- wgcnareslist$colors
  
  if(plot == TRUE){
    
    #colorlabels <- WGCNA::labels2colors(labels = numericlabels)
    colorlabels <- maphuecolors(numericvec = numericlabels)
    names(colorlabels) <- names(numericlabels)
    
    if(is.null(plottitleprefix)){
      plottitleprefix <- paste0('Consensus ', featuretype, ' dendrogram and module colors')
    }else{
      plottitleprefix <- paste0(plottitleprefix, ' consensus ', featuretype, 
                                ' dendrogram and module colors')
    }
    
    
    if(length(wgcnareslist$dendrograms) == 1){
      
      print(
        
        WGCNA::plotDendroAndColors(wgcnareslist$dendrograms[[1]], 
                                   colorlabels[wgcnareslist$blockGenes[[1]]], 
                                   groupLabels = 'Module colors', 
                                   dendroLabels = FALSE, 
                                   hang = 0.03,
                                   addGuide = TRUE, 
                                   guideHang = 0.05, 
                                   main = plottitleprefix)
        
      )
    }else{
      
      # Plot the dendrogram and the module colors underneath for each block
      for(j in 1:length(wgcnareslist$dendrograms)){
        
        print(
          
          WGCNA::plotDendroAndColors(wgcnareslist$dendrograms[[j]], 
                                     colorlabels[wgcnareslist$blockGenes[[j]]], 
                                     groupLabels = 'Module colors', 
                                     dendroLabels = FALSE, 
                                     hang = 0.03,
                                     addGuide = TRUE, 
                                     guideHang = 0.05,
                                     main = paste0(plottitleprefix, ' in block ', j))
          
        )
      }
    }
    
    
  }
  
  res <- list(mergedmelist = MElist, numericlabels = numericlabels)
  
  if(returnconsensustom == TRUE){
    res$consensustom <- consensustom
  }
  
  if(plot == TRUE){
    res$colorlabels <- colorlabels
  }
  
  return(res)
  
}

limmabarplots <- function(limmares, 
                          betaval = TRUE, 
                          diffcutoff = 0, 
                          pvalcut = 0.05, 
                          pvalcolname = 'adj.P.Val', 
                          confoundings = NULL, 
                          responsevarname = NULL, 
                          titleprefix = NULL, 
                          titlesize = 15, 
                          textsize = 13, 
                          face = 'bold', 
                          mesizes = NULL, 
                          xtextangle = 0){
  
  if(betaval == TRUE){
    facetname <- 'Difference'
  }else{
    facetname <- 'log2FC'
  }
  
  vars <- paste(c(responsevarname, confoundings), collapse = ' + ')
  
  sigg.limma <- limmares[(limmares[,pvalcolname] < pvalcut) & 
                           (abs(limmares$logFC) > diffcutoff),]
  nonsigg.limma <- limmares[(limmares[,pvalcolname] >= pvalcut) | 
                              (abs(limmares$logFC) <= diffcutoff),]
  
  if(nrow(sigg.limma) == 0){
    
    both_pos <- limmares[c('logFC', pvalcolname)]
    both_pos$Type <- 'NonSig'
    both_pos$ME <- row.names(both_pos)
    
    row.names(both_pos) <- 1:nrow(both_pos)
    
  }else{
    
    #Draw probe valcano
    both_pos <- sigg.limma[c('logFC', pvalcolname)]
    
    both_pos <- as.data.frame(both_pos, stringsAsFactors = FALSE)
    
    both_pos$Type <- 'NonSig'
    both_pos$Type[both_pos$logFC > diffcutoff] <- 'UP'
    both_pos$Type[both_pos$logFC < -diffcutoff] <- 'DN'
    
    hypernum <- nrow(subset(both_pos, Type == 'UP'))
    hyponum <- nrow(subset(both_pos, Type == 'DN'))
    
    both_pos_nonsig <- subset(both_pos, Type == 'NonSig')
    both_pos <- subset(both_pos, Type != 'NonSig')
    
    both_nonsig <- nonsigg.limma[c('logFC', pvalcolname)]
    both_nonsig$Type <- 'NonSig'
    both_nonsig <- rbind(both_nonsig, both_pos_nonsig)
    
    nonsignum <- nrow(both_nonsig)
    
    both_pos <- rbind(both_nonsig, both_pos)
    
    both_pos$ME <- row.names(both_pos)
    
    ncnum <- nrow(both_pos) - hypernum - hyponum
    
    row.names(both_pos) <- 1:nrow(both_pos)
    
  }
  
  names(both_pos)[names(both_pos) == pvalcolname] <- 'pval'
  
  if(pvalcolname == 'P.Value'){
    pvalfacetname <- 'p-value'
  }else{
    pvalfacetname <- 'adjusted p-value'
  }
  
  if(is.null(titleprefix)){
    maintitle <- paste0('Differential WGCNA modules (', vars, ')')
  }else{
    maintitle <- paste0(titleprefix, ' differential WGCNA modules (', vars, ')')
  }
  
  subtitle <- NULL
  
  if(!is.null(mesizes)){
    if(length(mesizes) <= 5){
      
      subtitle <- paste0('ME', seq(0, length(mesizes) - 1, 1), 
                         ' = ', as.vector(mesizes))
      subtitle <- paste0('(', paste(subtitle, collapse = '; '), ')')
      
    }
    
  }
  
  
  
  barplotdat <- reshape::melt(both_pos, c('Type', 'ME'))
  barplotdat$variable <- as.character(barplotdat$variable)
  barplotdat$variable[barplotdat$variable == 'pval'] <- pvalfacetname
  
  mes <- unique(as.numeric(gsub(pattern = 'ME', replacement = '', x = barplotdat$ME)))
  mes <- mes[order(mes)]
  mes <- paste0('ME', mes)
  barplotdat$ME <- factor(x = barplotdat$ME, levels = mes, ordered = TRUE)
  barplotdat$value[barplotdat$variable == pvalfacetname] <- 
    -log10(barplotdat$value[barplotdat$variable == pvalfacetname])
  barplotdat$variable <- factor(x = barplotdat$variable, 
                                levels = c(pvalfacetname, 'logFC'), 
                                ordered = TRUE)
  
  
  
  
  
  data_hline <- data.frame(variable = c(setdiff(barplotdat$variable, 'logFC'), 
                                        rep('logFC', 3)), 
                           hline = NA, 
                           stringsAsFactors = FALSE)
  data_hline$hline[data_hline$variable == pvalfacetname] <- -log10(pvalcut)
  data_hline$hline[data_hline$variable == 'logFC'] <- c(diffcutoff, -diffcutoff, 0)
  data_hline$color <- c('red', 'blue', 'blue', 'gray')
  
  
  facetnames <- list(pvalfacetname = paste0('-log10(', pvalfacetname, ')'), 
                     logFC = facetname)
  names(facetnames) <- c(pvalfacetname, 'logFC')
  
  facetlabeller <- function(variable, value){
    return(facetnames[value])
  }
  
  
  p <- ggplot2::ggplot(barplotdat, ggplot2::aes(x = ME, y = value))
  
  print(
    p + ggplot2::geom_bar(ggplot2::aes(fill = ME),
                          stat = 'identity') + 
      ggplot2::xlab(label = '') +
      ggplot2::ylab(label = 'Value') +
      ggplot2::ggtitle(maintitle, subtitle = subtitle) + 
      ggplot2::scale_fill_manual(values = c('#C0C0C0', 
                                            scales::hue_pal()(length(unique(barplotdat$ME)) - 1))) + 
      ggplot2::facet_wrap(ggplot2::vars(variable), scales = 'free', ncol = 2, 
                          labeller = facetlabeller) + 
      
      ggplot2::geom_hline(data = data_hline, ggplot2::aes(yintercept = hline), 
                          color = data_hline$color) + 
      ggplot2::theme_bw() +
      
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = face, size = textsize),
                     axis.title.y = ggplot2::element_text(face = face, size = textsize), 
                     axis.text.x = ggplot2::element_text(face = face, size = textsize, 
                                                         angle = xtextangle),
                     plot.title = ggplot2::element_text(size = titlesize, face = face),
                     plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                     strip.text = ggplot2::element_text(size = textsize - 1, face = face), 
                     legend.position = 'none')
  )
  
  
  row.names(both_pos) <- both_pos$ME
  both_pos <- both_pos[,c('pval', 'logFC', 'Type')]
  names(both_pos)[1] <- pvalcolname
  
  return(both_pos)
  
}

calconn <- function(edgedat){
  
  dat1 <- edgedat[c('fromNode', 'toNode', 'weight')]
  dat2 <- edgedat[c('toNode', 'fromNode', 'weight')]
  names(dat1) <- names(dat2) <- c('node1', 'node2', 'weight')
  dat <- rbind(dat1, dat2)
  dat <- unique(dat)
  
  if(!is.null(dat)){
    if(nrow(dat) > 0){
      row.names(dat) <- 1:nrow(dat)
    }
  }
  
  
  dat <- dat[c('node1', 'weight')]
  names(dat)[1] <- 'node'
  
  #library(plyr)
  
  calsum <- function(block){
    
    nodename <- unique(block[,1])
    numblock <- block[c(2)]
    numsum <- colSums(numblock)
    resblock <- data.frame(node = nodename, conn = numsum, stringsAsFactors = FALSE)
    row.names(resblock) <- 1:nrow(resblock)
    return(resblock)
    
  }
  
  dat <- dat[order(dat$node),]
  
  if(!is.null(dat)){
    if(nrow(dat) > 0){
      row.names(dat) <- 1:nrow(dat)
    }
  }
  
  conres <- plyr::ddply(.data = dat, .variables = c('node'), .fun = calsum)
  
  return(conres)
}

orgcyt <- function(cyt){
  
  edgedata <- cyt$edgeData
  nodedata <- cyt$nodeData
  edgedata <- edgedata[c(1, 2, 3)]
  nodedata <- nodedata[c(1)]
  
  res <- list(edgedata = edgedata, 
              nodedata = nodedata)
  
  return(res)
}

calcompleteconn <- function(tom){
  
  nodeconn <- rowSums(tom)
  nodeconn <- data.frame(node = names(nodeconn), 
                         conn = as.vector(nodeconn), 
                         stringsAsFactors = FALSE)
  nodeconn <- nodeconn[order(nodeconn$node),]
  row.names(nodeconn) <- 1:nrow(nodeconn)
  
  return(nodeconn)
  
}

prepareenrich <- function(wgcnares, 
                          plot = TRUE, 
                          seed = 2022, 
                          titleprefix = NULL){
  
  if(!is.null(dim(wgcnares$mergedmelist))){
    menames <- colnames(wgcnares$mergedmelist)
  }else{
    menames <- colnames(wgcnares$mergedmelist[[1]])
    names(wgcnares)[names(wgcnares) == 'consensustom'] <- 'tom'
  }
  
  menames <- setdiff(menames, 'ME0')
  
  if(mode(wgcnares$tom) == 'list'){
    return(NULL)
  }
  
  weightcutoff <- quantile(wgcnares$tom, 0.99)
  
  modulenodedata <- data.frame(nodeName = character(0), 
                               ME = character(0), 
                               stringsAsFactors = FALSE)
  moduleedgedata <- data.frame(fromNode = character(0), 
                               toNode = character(0), 
                               weight = numeric(0), 
                               ME = character(0), 
                               stringsAsFactors = FALSE)
  modulenodeconn <- data.frame(node = character(0), 
                               conn = numeric(0), 
                               ME = character(0), 
                               stringsAsFactors = FALSE)
  
  completemodulenodedata <- data.frame(nodeName = character(0), 
                                       ME = character(0), 
                                       stringsAsFactors = FALSE)
  completemoduleedgedata <- data.frame(fromNode = character(0), 
                                       toNode = character(0), 
                                       weight = numeric(0), 
                                       ME = character(0), 
                                       stringsAsFactors = FALSE)
  completemodulenodeconn <- data.frame(node = character(0), 
                                       conn = numeric(0), 
                                       ME = character(0), 
                                       stringsAsFactors = FALSE)
  
  p <- 1
  for(p in 1:length(menames)){
    
    mename <- menames[p]
    
    targetme <- gsub(pattern = 'ME', replacement = '', x = mename)
    targetme <- as.numeric(targetme)
    targetmefeatures <- names(wgcnares$numericlabels[wgcnares$numericlabels == targetme])
    
    subTOM <- wgcnares$tom[targetmefeatures, targetmefeatures]
    
    subcyt <- WGCNA::exportNetworkToCytoscape(subTOM, 
                                              edgeFile = NULL, 
                                              nodeFile = NULL, 
                                              weighted = TRUE, 
                                              threshold = weightcutoff, 
                                              nodeNames = row.names(subTOM))
    
    #if(nrow(subcyt$edgeData) == 0 | nrow(subcyt$nodeData) == 0){
      #next()
    #}
    
    subcytres <- orgcyt(cyt = subcyt)
    subnodedata <- subcytres$nodedata
    subedgedata <- subcytres$edgedata
    
    
    
    if(!is.null(subnodedata)){
      if(nrow(subnodedata) > 0){
        subnodedata$ME <- mename
      }
    }
    
    if(!is.null(subedgedata)){
      if(nrow(subedgedata) > 0){
        subedgedata$ME <- mename
      }
    }
    
    subnodeconn <- calconn(edgedat = subedgedata)
    
    if(!is.null(subnodeconn)){
      if(nrow(subnodeconn) > 0){
        subnodeconn$ME <- mename
      }
    }
    
    modulenodedata <- rbind(modulenodedata, subnodedata)
    moduleedgedata <- rbind(moduleedgedata, subedgedata)
    modulenodeconn <- rbind(modulenodeconn, subnodeconn)
    
    
    subcyt <- WGCNA::exportNetworkToCytoscape(subTOM, 
                                              edgeFile = NULL, 
                                              nodeFile = NULL, 
                                              weighted = TRUE, 
                                              threshold = 0, 
                                              nodeNames = row.names(subTOM))
    
    subcytres <- orgcyt(cyt = subcyt)
    
    subnodedata <- subcytres$nodedata
    subedgedata <- subcytres$edgedata
    subnodedata$ME <- mename
    subedgedata$ME <- mename
    subnodeconn <- calcompleteconn(tom = subTOM)
    subnodeconn$ME <- mename
    
    completemodulenodedata <- rbind(completemodulenodedata, subnodedata)
    completemoduleedgedata <- rbind(completemoduleedgedata, subedgedata)
    completemodulenodeconn <- rbind(completemodulenodeconn, subnodeconn)
    
  }
  
  
  cyt <- WGCNA::exportNetworkToCytoscape(wgcnares$tom, 
                                         edgeFile = NULL, 
                                         nodeFile = NULL, 
                                         weighted = TRUE, 
                                         threshold = weightcutoff, 
                                         nodeNames = row.names(wgcnares$tom))
  cytres <- orgcyt(cyt = cyt)
  
  nodedata <- cytres$nodedata
  edgedata <- cytres$edgedata
  nodeconn <- calconn(edgedat = edgedata)
  
  #library(igraph)
  
  inet <- igraph::graph_from_data_frame(d = edgedata, directed = FALSE, vertices = nodedata)
  
  vertexattr <- igraph::vertex.attributes(inet)
  edgeattr <- igraph::edge.attributes(inet)
  
  
  labels <- 
    scales::hue_pal()(max(length(unique(wgcnares$numericlabels[nodedata$nodeName])) - 
                            1, 1))
  
  if(length(labels) == 1){
    labels <- '#C0C0C0'
  }else{
    labels <- c('#C0C0C0', labels)
  }
  
  labels <- labels[as.vector(wgcnares$numericlabels[nodedata$nodeName] + 1)]
  vsize <- log10(nodeconn$conn/max(nodeconn$conn))
  vsize <- scales::rescale(vsize, to = c(1, 10))
  esize <- log10(edgeattr$weight)
  esize <- scales::rescale(esize, to = c(1, 10))
  
  vertexattr$nodeColors <- labels
  igraph::V(inet)$size <- vsize
  igraph::E(inet)$width <- esize
  
  
  cyt <- WGCNA::exportNetworkToCytoscape(wgcnares$tom, 
                                         edgeFile = NULL, 
                                         nodeFile = NULL, 
                                         weighted = TRUE, 
                                         threshold = 0, 
                                         nodeNames = row.names(wgcnares$tom))
  cytres <- orgcyt(cyt = cyt)
  
  completenodedata <- cytres$nodedata
  completeedgedata <- cytres$edgedata
  completenodeconn <- calcompleteconn(tom = wgcnares$tom)
  
  res <- list(backnodedata = nodeconn, 
              backedgedata = edgedata, 
              nodedata = modulenodeconn, 
              edgedata = moduleedgedata, 
              inet = inet, 
              completebacknodedata = completenodeconn, 
              completebackedgedata = completeedgedata, 
              completenodedata = completemodulenodeconn, 
              completeedgedata = completemoduleedgedata)
  
  
  
  if(plot == TRUE){
    
    set.seed(seed)
    
    netstyle <- igraph::layout_with_fr(inet, coords = NULL, 
                                       dim = 2, niter = 1000, 
                                       grid = c("nogrid"))
    
    stamp <- Sys.time()
    stamp <- gsub(pattern = ' ', replacement = '_', x = stamp)
    stamp <- gsub(pattern = ':', replacement = '-', x = stamp)
    stamp <- paste0('.', stamp)
    
    #pdf(paste0('wgcnaplot', stamp, '.pdf'), height=10, width=10, useDingbats=FALSE)
    
    if(!is.null(titleprefix)){
      main <- paste0(titleprefix, ' WGCNA modules on top 1% topological overlap edges')
      if(nchar(main) > 60){
        main <- paste0(titleprefix, ' WGCNA modules on top 1%\ntopological overlap edges')
      }
      
    }else{
      main <- 'WGCNA modules on top 1% topological overlap edges'
    }
    
    
    print(
      
      plot(inet, vertex.label = NA, vertex.color = vertexattr$nodeColors, 
           edge.color = '#C0C0C0', 
           layout = netstyle, 
           main = main)
      
      
    )
    #dev.off()
    
  }
  
  
  return(res)
  
}

orgwgcnagores <- function(gores){
  
  mename <- unique(gores$pairshuffleweights$ME)
  
  if(!is.null(gores$generes)){
    if(nrow(gores$generes) > 0){
      gores$generes$ME <- mename
    }
  }
  
  if(!is.null(gores$pairres)){
    if(nrow(gores$pairres) > 0){
      gores$pairres$ME <- mename
    }
    
  }
  
  if(!is.null(gores$geneshuffleres)){
    if(nrow(gores$geneshuffleres) > 0){
      gores$geneshuffleres$ME <- mename
    }
    
  }
  
  if(!is.null(gores$pairshuffleres)){
    if(nrow(gores$pairshuffleres) > 0){
      gores$pairshuffleres$ME <- mename
    }
  }
  
  if(!is.null(gores$geneweights)){
    if(nrow(gores$geneweights) > 0){
      gores$geneweights$ME <- mename
    }
  }
  
  if(!is.null(gores$pairweights)){
    if(nrow(gores$pairweights) > 0){
      gores$pairweights$ME <- mename
    }
  }
  
  
  return(gores)
  
}

#'Perform causal WGCNA analysis
#'
#'Perform WGCNA analysis followed by network-based functional enrichment and 
#'  causal inference (based on mediation analysis) for WGCNA modules.
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
#'@param responselevels If the variable of \code{responsevarname} is discrete, 
#'  it will be treated as a factor, and this parameter is needed to define the 
#'  factor element level, so it should be a vector defining this level. Only 
#'  needed if the response in \code{pddat} is a character vector but not a 
#'  factor. In this case, it will be converted to a factor by the function, 
#'  and this paramter \code{responselevels} will be used to define the factor 
#'  element level. Its default is NULL, meaning that this element level will 
#'  be set automatically following the character order of the elements in the 
#'  response variable.
#'@param confoundings The column name of the confounding factors in the data 
#'  frame provided to \code{pddat}. Default is NULL, which means there is no 
#'  confounding factor in the data. 
#'@param topvariancetype The method need to calculate feature variance. Its 
#'  default is 'sd', meaning standard deviation. Can also be 'mad', meaning 
#'  median absolute deviation.
#'@param topvaricance The number of top variance features need to be included 
#'  in the WGCNA analysis. Default is 5000. 
#'@param balanceadj Whether to use the ensemble-based mode. When it is set 
#'  as TRUE, a sampling process will be performed on the samples to generate 
#'  several base learner datasets, and the sample groups of each base learner 
#'  set will be adjusted during the sampling so that they will have the same 
#'  size. After calling the WGCNA modules, a differential feature selection 
#'  step will be performed to find the modules with significantly different 
#'  eigengene values between the sample groups, and also find the features 
#'  with significantly different values between the sample groups and within 
#'  each module. At this stage, limma will be used, and if the parameter here 
#'  is set as TRUE, limma will be used on each base learner first to call the 
#'  differential features. After that, their results will finally be ensembled 
#'  together to get the result for the whole data. If this parameter is set as 
#'  FALSE, normal limma will be used. Similarly, for the mediation analysis to 
#'  find the module features mediating the causal relationships between the 
#'  module and the response variable, it can also be performed on the base 
#'  learner datasets with an ensemble mode, or directly on the whole data with 
#'  a normal mode, depending on this parameter. Default is TRUE. 
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
#'@param seed Random seed for some random processes, such as base learner sets 
#'  generation, shuffling-based functional enrichment, etc.
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
#'@param sftpowers The soft-thresholding power to be used for WGCNA network 
#'  construction. Should be an integer or NULL. Default is NULL, and in this 
#'  case, the function will calculate and pick up an appropriate power from 
#'  the values provided to another parameter \code{powers}.
#'@param powers A vector used to provide candidate soft-thresholding power for 
#'  WGCNA network construction. Should be a vector with integers as elements. 
#'  If the parameter \code{sftpowers} is NULL, the function will calculate the 
#'  scale free topology fitting index for each element of this vector and use 
#'  the one with the largest index as the final power. Default is the integer 
#'  sequence from 1 to 20.
#'@param rsqcutline A float defining the desired minimum scale free topology 
#'  fitting index to select the soft-thresholding power from the candidates 
#'  included in the parameter \code{powers}, in the case that the parameter 
#'  \code{sftpowers} is NULL and the soft-thresholding powers for the data 
#'  need to be calculated by the function. Default is 0.8.
#'@param minclustersize The minimum module size for module detection during 
#'  WGCNA network construction. Default is 50. 
#'@param mergecutheight The dendrogram cut height for WGCNA module merging. It 
#'  default is 0.2.
#'@param maxblocksize Integer giving maximum block size for module detection 
#'  during WGCNA network construction. Default is 50. 
#'@param removereg A GRanges object or a data frame that records the genomic 
#'  coordinates of any regions that need to be excluded from the analysis, can 
#'  be the genomic regions related to any confounding factors. If a feature 
#'  is located there, it will be removed at the beginning. For regions related 
#'  to confounding factors, this region removal can help avoid influence from 
#'  the confounding factors. However, for the limma and mediation analyses, 
#'  because they can adjust confoundings, this region removal will have very 
#'  limited effects. Hence, this step is unnecessary, and can be skipped by 
#'  setting this parameter \code{removereg} as NULL. If need to transfer a 
#'  data frame, it should have columns named as "seqnames", "start", "end", 
#'  and "strand", which help to define the genomic coordinates of the regions. 
#'  Each row of the data frame should be a region. The default value is NULL.
#'@param featuretype The type of the features in the data. It will be used 
#'  when \code{removereg} is not NULL and judge whether any features located 
#'  in the genomic regions to be removed. Can be "gene" or "probe" (for probes 
#'  in the DNA methylation data). If another parameter \code{plot} is TRUE, 
#'  it will also appear in the title of the dendrogram plot and volcano plot 
#'  generated by this function. In addition, only when this \code{featuretype} 
#'  is "gene", the functional enrichment analysis will be performed after the 
#'  WGCNA network construction. Default value is "gene".
#'@param topoenrichment After the WGCNA module detection, this function can 
#'  perform network-based gene functional enrichment on each module, i.e., 
#'  the enrichment analysis is performed not only on the gene level, but also 
#'  on the gene-gene pair level. If this parameter is TRUE, the enrichment 
#'  will be performed. If it is FALSE, this step will be skipped. Default is 
#'  FALSE. 
#'@param shufflenum For the functional enrichment analysis, the function can 
#'  perform it with both a hypergeometric enrichment method, and a gene pair 
#'  weight or gene weight shuffling enrichment method. This parameter is used 
#'  to define the shuffling number for the weight shuffling method. Default is 
#'  1000.
#'@param mediation After the WGCNA module detection and limma analysis, this 
#'  function can perform causal inference (mediation analysis) on each module 
#'  with a significant eigengene difference between the sample groups defined 
#'  by the parameter \code{responsevarname} in the data frame \code{pddat}. In 
#'  each of these modules, each module feature with a significant difference 
#'  between the sample groups will be used as the potential mediator of this 
#'  mediation analysis (causal inference) and for each feature in each module, 
#'  2 kinds of mediation relationships will be tested. The first is "response 
#'  (in \code{pddat}) -> module feature -> module eigengene", and the other is 
#'  the opposite direction "module eigengene -> module feature -> response". 
#'  For the mediation analysis (causal inference), the response variable need 
#'  to be either a continuous or a binary variable. If the parameter here is 
#'  TRUE, the mediation analysis will be performed to test the above 2 kinds 
#'  of mediation relationships and will return the significant ones in the 
#'  final result. If this parameter is FALSE, the mediation analysis will not 
#'  be performed. Default is FALSE.
#'@param topn If the parameter \code{mediation} is set as TRUE, this parameter 
#'  \code{topn} will be used to define the number of top differential modules 
#'  needed to be analyzed by the mediation analysis. If it is set as 1, the 
#'  most differential module between the sample groups will be analyzed by the 
#'  mediation model. If it is set as 2, the top 2 most differential modules 
#'  will go through the analysis, respectively. If it is set as 3, the top 3 
#'  most differential modules will be analyzed, etc. The default value is NULL, 
#'  meaning that all the differential modules will go through the mediation 
#'  analysis.
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
#'@param plot Whether need to generate several plots during the analysis, i.e. 
#'  the dendrogram plot and the network plot to show the WGCNA modules, the 
#'  bar plot to show the significantly differential WGCNA module eigengenes 
#'  between the sample groups, the volcano plots to show the significantly 
#'  differential module features between the sample groups and within each 
#'  module, the volcano plots to show the significant mediation relationships 
#'  within each module. Default is FALSE so that no plots will be generated.
#'@param titleprefix The prefix of the titles of the dendrogram plot, network 
#'  plot, barplot and volcano plot. It can be set as any character string need 
#'  to be shown in the titles. Default is NULL.
#'@param isbetaval If the feature values in the data are log2 transformed, 
#'  set it as FALSE, if they are not log2 transformed, i.e., they are gene 
#'  counts or methylation beta values, set this parameter as TRUE. Will only 
#'  be needed when \code{plot} is TRUE.
#'@param diffcutoff The cutoff on the inter-group difference to judge whether 
#'  the WGCNA module eigengenes are significantly different between the sample 
#'  groups in the limma analysis. Default is 0, meaning any difference value 
#'  will not influence the judgment of the module eigengene difference, and 
#'  the p-val is the only criterion to judge the significance of eigengene 
#'  difference.
#'@param pvalcutoff The p-val (or adjusted p-val) cutoff to judge whether a 
#'  feature is significantly different between the groups or not in the limma 
#'  analysis. Default is 0.05.
#'@param pvalcolname Which p-val should be used to judge the significance of 
#'  the features in the limma analysis. Can be "P.Val" or "adj.P.Val". Default 
#'  is "adj.P.Val".
#'@param absxcutoff The cutoff on log2FC (for log2 transformed features) or 
#'  inter-group difference (for non-log2 transformed feature counts or beta 
#'  values) to judge whether a feature is significantly different between the 
#'  sample groups in the limma analysis. Default is 0, meaning any log2FC or 
#'  difference value will not influence the judgment on feature significance.
#'@param titlesize The font size of the barplot and volcano plot titles. Its 
#'  default is 15.
#'@param textsize The font size of the barplot and volcano plot texts. Default 
#'  is 13.
#'@param face The font face of the barplot and volcano plot texts. Default is 
#'  "bold".
#''@param labelnum In the volcano plot, the names of the top up-regulated and 
#'  top down-regulated features can be labeled. This number indicates how many 
#'  top features in the groups should be labeled. If it is NULL, no feature 
#'  name will be labeled. Can also be any number. Default is NULL.
#'@param annotextsize If \code{labelnum} is not NULL, this parameter will be 
#'  needed to set the font size of the gene names to be labeled in the volcano 
#'  plot. Default is 4.
#'@param savetom Whether the topological overlap matrix (TOM) of the WGCNA 
#'  analysis need to be returned or not. Default is TRUE.
#'@param signtype To construct signed or unsigned WGCNA network, can be set 
#'  as "signed" or "unsigned". Default is "unsigned".
#'@return A list with several slots will be returned. One named "wgcnares" 
#'  records the WGCNA module detection results, including the module labels 
#'  of the features, the soft-thresholding power used for the module calling, 
#'  and the module eigengenes of the samples. If the paramter \code{savetom} 
#'  is TRUE, this slot will also contain the TOM matrix of the WGCNA analysis. 
#'  Another slot named "limmares" records the limma analysis results on the 
#'  module eigengene difference between sample groups. The "melimmares" slot 
#'  is the limma results on the feature difference between sample groups and 
#'  within each module. There are also several slots recording the functional 
#'  enrichment results for the modules. The slots named as "genepaths" and 
#'  "pairpaths" are the hypergeometric enrichment results on the gene level 
#'  and gene-gene pair level, respectively, and the functional terms of the 
#'  genes or gene pairs are based on the Reactome database. The gene terms are 
#'  directly the ones in this database, and the terms of a gene pair is the 
#'  intersection of its 2 genes'. Similarly, the slots named "genebps" and 
#'  "pairbps" are the hypergeometric enrichment results on gene and gene pair 
#'  levels, but the functional term database is the GOBP database, rather than 
#'  Reactome. The result list also has other slots named "geneshufflepaths", 
#'  "pairshufflepaths", "geneshufflebps", and "pairshufflebps", which are the 
#'  corresponding functional results generated with the gene and gene pair 
#'  weight shuffling method, not the hypergeometric test method. It should be 
#'  noted that, for the results on gene pairs, because the WGCNA network is 
#'  a fully-connected network and the gene-gene combination can generate a 
#'  huge number of gene pairs, which will require huge computational cost to 
#'  get their functional terms via the intersection computation. Hence, the 
#'  hypergeometric results on the gene pairs only restricted to the gene pairs 
#'  with the top 1% weight of the whole WGCNA network. The detailed functional 
#'  terms of the gene and gene pairs are in the slots of "geneweightspaths", 
#'  "pairweightspaths", and "pairshuffleweightspaths", etc, corresponding to 
#'  the results in the slots "genepaths", "pairpaths", and "pairshufflepaths", 
#'  etc. For the mediation analysis results, they are included in another slot 
#'  named "mediationres", a data frame showing the significant medation casual 
#'  relationships for the WGCNA modules. If the mediators are genes, the slots 
#'  "mediatorgenedat" contains the biological information of the genes, and if 
#'  the mediators are DNA methylation probes, the slot "mediatorprobedat" are 
#'  the information of these probes.
#'@export
diffwgcna <- function(dat, 
                      pddat, 
                      responsevarname = NULL, 
                      responselevels = NULL, 
                      confoundings = NULL,  
                      
                      topvaricancetype = 'sd', 
                      topvaricance = 5000, 
                      
                      balanceadj = TRUE, 
                      samplingmethod = 'updn', 
                      nround = 10,
                      seed = 2022, 
                      k = 5, 
                      threads = 1, 
                      method = 'Cauchy', 
                      
                      sftpowers = NULL, 
                      powers = seq(1, 20, 1), 
                      rsqcutline = 0.8, 
                      
                      minclustersize = 50, 
                      mergecutheight = 0.2, 
                      maxblocksize = 5000, 
                      
                      removereg = NULL, 
                      featuretype = 'gene', 
                      topoenrichment = FALSE, 
                      #species = 'human', 
                      shufflenum = 1000, 
                      mediation = FALSE, 
                      
                      topn = NULL, 
                      
                      platform = 450, 
                      ignorestrand = TRUE, 
                      tssradius = 1500, 
                      plot = FALSE, 
                      titleprefix = NULL, 
                      isbetaval = TRUE, 
                      diffcutoff = 0, 
                      pvalcutoff = 0.05, 
                      pvalcolname = 'adj.P.Val', 
                      absxcutoff = 0, 
                      titlesize = 15, 
                      textsize = 13, 
                      face = 'bold', 
                      labelnum = NULL, 
                      annotextsize = 4, 
                      savetom = TRUE, 
                      
                      signtype = 'unsigned'){
  
  if(featuretype == 'gene' & topoenrichment == TRUE){
    
    #cat('Note: the gene name will be converted to the form of "symbol::entrez"\n')
    
    mappings1 <- getmultinamedat(multigenenames = row.names(dat), entrez = TRUE)
    mappings2 <- getmultinamedat(multigenenames = row.names(dat), entrez = FALSE)
    
    mappings <- rbind(mappings1, mappings2)
    mappings <- unique(mappings)
    
    if(!is.null(mappings)){
      
      sub1 <- mappings[c('SYMBOL', 'genenames')]
      sub2 <- mappings[c('ENTREZID', 'genenames')]
      sub3 <- mappings[c('genenames', 'genenames')]
      colnames(sub1) <- colnames(sub2) <- colnames(sub3) <- c('orinames', 'genenames')
      sub <- rbind(sub1, sub2, sub3)
      sub <- unique(sub)
      row.names(sub) <- 1:nrow(sub)
      
      idx <- match(row.names(dat), sub[,c('orinames')])
      others <- dat[is.na(idx), , drop = FALSE]
      dat <- dat[!is.na(idx), , drop = FALSE]
      idx <- idx[!is.na(idx)]
      
      row.names(dat) <- sub[,c('genenames')][idx]
      dat <- rbind(dat, others)
      
    }
    
  }
  
  if(!is.null(removereg)){
    
    if(featuretype == 'probe'){
      removedfeatures <- overlappedprobes(totalprobes = row.names(dat), 
                                          removereg = removereg, 
                                          platform = platform, 
                                          ignorestrand = ignorestrand)
    }else if(featuretype == 'gene'){
      removedfeatures <- overlappedgenes(totalgenes = row.names(dat), 
                                         removereg = removereg, 
                                         tssradius = tssradius, 
                                         ignorestrand = ignorestrand)
    }else{
      removedfeatures <- NULL
    }
    
    
    if(!is.null(removedfeatures)){
      
      if(featuretype == 'probe'){
        removedfeatures <- unique(row.names(removedfeatures))
        
        dat <- dat[!(row.names(dat) %in% removedfeatures), , drop = FALSE]
        
      }else if(featuretype == 'gene'){
        
        removegenes <- idxgenenames(totalgenes = row.names(dat), 
                                    genesym = removedfeatures$genesym, 
                                    geneid = removedfeatures$geneid)
        
        dat <- dat[!removegenes, , drop = FALSE]
        
      }else{
        
        dat <- dat
      }
      
    }else{
      dat <- dat
    }
    
  }
  
  
  if(is.null(responsevarname) & ncol(pddat) >= 2){
    responsevarname <- names(pddat)[2]
  }
  
  if(is.null(responsevarname)){
    return(NULL)
  }
  
  
  returnearly <- FALSE
  noenrichcount <- 0
  
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
    
    returnearly <- TRUE
    
    if(mediation == TRUE){
      cat('The mediation analysis will be skipped because it only supports binary 
        or continuous variables but the current response is not.\n')
    }
    
    mediation <- FALSE
    
  }
  
  if((Responsetype != 'binary') & (balanceadj == TRUE)){
    
    cat('The ensemble mode only supports binary variable, but the current response is not.\nHence, the mode has been changed to normal mode.\n')
    
    balanceadj <- FALSE
    
  }
  
  
  
  if(is.numeric(topvaricance)){
    
    topvaricance <- min(topvaricance, nrow(dat))
    
    dat <- featuresampling(betas = dat, 
                           topfeatures = topvaricance, 
                           pddat = NULL, 
                           variancetype = topvaricancetype, 
                           anova = FALSE, 
                           responsevarname = NULL, 
                           threads = threads)
    
    dat <- dat$betas
    
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
  }
  
  if(FALSE){
    
    balanceddats <- balancesampling(dat = dat, 
                                    pddat = pddat, 
                                    responsename = responsevarname, 
                                    nround = nround, 
                                    seed = seed, 
                                    k = k, 
                                    threads = threads, 
                                    samplingmethod = samplingmethod)
    
    if(is.null(sftpowers)){
      
      sftpowers <- esftpower(dats = balanceddats$dats, 
                             powers = powers, 
                             rsqcutline = rsqcutline, 
                             threads = threads, 
                             mergesft = TRUE, 
                             signtype = signtype)
      
      
    }else{
      
      if(length(sftpowers) != length(balanceddats$dats)){
        
        sftpowers <- esftpower(dats = balanceddats$dats, 
                               powers = powers, 
                               rsqcutline = rsqcutline, 
                               threads = threads, 
                               mergesft = TRUE, 
                               signtype = signtype)
        
      }else{
        
        sftpowers <- list(power = sftpowers)
        
      }
      
    }
    
    wgcnares <- ewgcna(dats = balanceddats$dats, 
                       powers = as.vector(unlist(sftpowers$power)), 
                       minclustersize = minclustersize, 
                       mergecutheight = mergecutheight, 
                       returnconsensustom = TRUE, 
                       plot = plot, 
                       featuretype = featuretype, 
                       plottitleprefix = titleprefix, 
                       threads = threads, 
                       maxblocksize = maxblocksize, 
                       seed = seed, 
                       
                       signtype = signtype)
    
  }else{
    
    if(!is.null(sftpowers)){
      
      if(length(unique(sftpowers)) == 1){
        
        wgcnares <- wgcnamain(dat = dat, 
                              sftpower = unique(sftpowers), 
                              minclustersize = minclustersize, 
                              mergecutheight = mergecutheight, 
                              returntom = TRUE, 
                              plot = plot, 
                              featuretype = featuretype, 
                              plottitleprefix = titleprefix, 
                              threads = threads, 
                              maxblocksize = maxblocksize, 
                              seed = seed, 
                              
                              signtype = signtype)
        
        
      }else{
        
        sftpowers <- NULL
        
      }
      
    }
    
    if(is.null(sftpowers)){
      
      sftpowers <- choosepower(datlist = list(dat), 
                               i = 1, 
                               powers = powers, 
                               rsqcutline = rsqcutline, 
                               plot = FALSE, 
                               signtype = signtype)
      
      if(is.null(sftpowers)){
        return(NULL)
      }
      
      wgcnares <- wgcnamain(dat = dat, 
                            sftpower = as.vector(unlist(sftpowers$powerEstimate)), 
                            minclustersize = minclustersize, 
                            mergecutheight = mergecutheight, 
                            returntom = TRUE, 
                            plot = plot, 
                            featuretype = featuretype, 
                            plottitleprefix = titleprefix, 
                            threads = threads, 
                            maxblocksize = maxblocksize, 
                            seed = seed, 
                            
                            signtype = signtype)
      
    }
    
  }
  
  gc()
  
  if(returnearly == TRUE & topoenrichment == FALSE){
    
    res <- list(wgcnares = wgcnares, sftpowers = sftpowers)
    
    if(savetom == FALSE & balanceadj == FALSE){
      res$wgcnares$tom <- NULL
    }else if(savetom == FALSE & balanceadj == TRUE){
      res$wgcnares$consensustom <- NULL
    }
    
    
    return(res)
    
  }
  
  
  #limma on module level
  
  if(returnearly == FALSE){
    
    cat('Start limma analysis on WGCNA modules\n')
    
    if(balanceadj == TRUE){
      
      #wgcnares <- readRDS('wgcna.balance.rds')
      
      
      
      
      MElist <- list()
      
      for(k in 1:length(balanceddats$dats)){
        
        MElist[[k]] <- WGCNA::moduleEigengenes(t(balanceddats$dats[[k]]), 
                                               colors = wgcnares$numericlabels)
        MElist[[k]] <- MElist[[k]]$eigengenes
        
        names(MElist)[k] <- paste0('round_', k)
        
      }
      
      wgcnares$mergedmelist <- MElist
      
      
      
      
      
      mergedmelist <- list()
      for(i in 1:length(wgcnares$mergedmelist)){
        mergedmelist[[i]] <- t(wgcnares$mergedmelist[[i]])
      }
      names(mergedmelist) <- names(wgcnares$mergedmelist)
      
      
      if(threads == 1){
        
        limmareslist <- list()
        
        iseqs <- seq(1, nround, 1)
        
        for(i in iseqs){
          
          limmares <- diffsites(dats = mergedmelist, 
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
                                           diffsites(dats = mergedmelist, 
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
      
      #wgcnares <- readRDS('wgcna.nobalance.rds')
      
      mergedmelist <- list(t(wgcnares$mergedmelist))
      
      balanceddats <- list(dats = list(dat), 
                           pds = list(pddat))
      
      limmares <- diffsites(dats = mergedmelist, 
                            pddats = balanceddats$pds, 
                            i = 1, 
                            responsevarname = responsevarname, 
                            responselevels = responselevels, 
                            confoundings = confoundings)
      
    }
    
    if(plot == TRUE){
      
      if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val', 
                                     'pvalue', 'p-value', 'p.value')){
        pvalcolname <- 'P.Value'
      }else{
        pvalcolname <- 'adj.P.Val'
      } 
      
      limmabarplotsangel <- 0
      
      if(nrow(limmares) > 10){
        limmabarplotsangel <- 45
      }
      
      plotdat <- limmabarplots(limmares = limmares, 
                               betaval = isbetaval, 
                               diffcutoff = diffcutoff, 
                               pvalcut = pvalcutoff, 
                               pvalcolname = pvalcolname, 
                               confoundings = confoundings, 
                               responsevarname = responsevarname, 
                               titleprefix = titleprefix, 
                               titlesize = titlesize, 
                               textsize = textsize, 
                               face = face, 
                               mesizes = table(wgcnares$numericlabels), 
                               xtextangle = limmabarplotsangel)
      
      
      
    }
    
    sigmes <- limmares[(abs(limmares$logFC) > diffcutoff) & 
                         (limmares[,pvalcolname] < pvalcutoff), , drop = FALSE]
    
    if(nrow(sigmes) > 0){
      sigmes <- row.names(sigmes)
    }else{
      sigmes <- NULL
    }
    
    gc()
    
  }
  
  
  #GO analysis
  if(balanceadj == TRUE){
    mergedmenames <- colnames(wgcnares$mergedmelist[[1]])
  }else{
    mergedmenames <- colnames(wgcnares$mergedmelist)
  }
  
  if(length(mergedmenames) > 1){
    
    #saveRDS(wgcnares, 'wgcnares.rds')
    
    enrichdats <- prepareenrich(wgcnares = wgcnares, plot = plot, seed = seed, 
                                titleprefix = titleprefix)
    
    if(is.null(enrichdats)){
      topoenrichment <- FALSE
    }
    
    ###modified!!!!
    if(featuretype == 'gene' & topoenrichment == TRUE & !is.null(enrichdats)){
      
      cat('Start gene function enrichment on WGCNA modules\n')
      
      genomebackground <- searchtable(allgenes = enrichdats$completebacknodedata$node, 
                                      species = 'human')
      
      q <- 1
      startq <- 1
      
      for(q in 1:(length(mergedmenames) - 1)){
        
        mename <- mergedmenames[q]
        
        if(!(mename %in% c(enrichdats$completenodedata$ME, 
                           enrichdats$completeedgedata$ME, 
                           enrichdats$edgedata$ME, 
                           enrichdats$nodedata$ME))){
          
          if(q == startq){
            startq <- startq + 1
          }
          
          noenrichcount <- noenrichcount + 1
          
          next()
          
        }
        
        
        geneweights <- subset(enrichdats$completenodedata, ME == mename)
        names(geneweights)[c(1, 2)] <- c('GENE', 'weight')
        row.names(geneweights) <- 1:nrow(geneweights)
        
        backweights <- enrichdats$completebacknodedata
        names(backweights)[c(1, 2)] <- c('GENE', 'weight')
        
        genepairweights <- subset(enrichdats$edgedata, ME == mename)
        
        if(!is.null(genepairweights)){
          if(nrow(genepairweights) > 0){
            row.names(genepairweights) <- 1:nrow(genepairweights)
          }
        }
        
        backpairweights <- enrichdats$backedgedata
        
        shufflegenepairweights <- subset(enrichdats$completeedgedata, ME == mename)
        row.names(shufflegenepairweights) <- 1:nrow(shufflegenepairweights)
        
        mebackground <- searchtable(allgenes = geneweights$GENE)
        
        pathres <- maingo(geneweights = geneweights, 
                          backweights = backweights, 
                          backgos = genomebackground$allpaths, 
                          backpairweights = backpairweights, 
                          genegos = mebackground$allpaths, 
                          genepairweights = genepairweights, 
                          shufflegenepairweights = shufflegenepairweights, 
                          genecount = FALSE, 
                          paircount = FALSE, 
                          addzeroedgesforshuffle = TRUE, 
                          threads = threads, 
                          seed = seed, 
                          shufflenum = shufflenum)
        
        bpres <- maingo(geneweights = geneweights, 
                        backweights = backweights, 
                        backgos = genomebackground$allbps, 
                        backpairweights = backpairweights, 
                        genegos = mebackground$allbps, 
                        genepairweights = genepairweights, 
                        shufflegenepairweights = shufflegenepairweights, 
                        genecount = FALSE, 
                        paircount = FALSE, 
                        addzeroedgesforshuffle = TRUE, 
                        threads = threads, 
                        seed = seed, 
                        shufflenum = shufflenum)
        
        pathres <- orgwgcnagores(gores = pathres)
        bpres <- orgwgcnagores(gores = bpres)
        
        
        genepath <- pathres$generes
        pairpath <- pathres$pairres
        geneshufflepath <- pathres$geneshuffleres
        pairshufflepath <- pathres$pairshuffleres
        
        geneweightspath <- pathres$geneweights
        pairweightspath <- pathres$pairweights
        pairshuffleweightspath <- pathres$pairshuffleweights
        
        
        genebp <- bpres$generes
        pairbp <- bpres$pairres
        geneshufflebp <- bpres$geneshuffleres
        pairshufflebp <- bpres$pairshuffleres
        
        geneweightsbp <- bpres$geneweights
        pairweightsbp <- bpres$pairweights
        pairshuffleweightsbp <- bpres$pairshuffleweights
        
        
        if(q == startq){
          genepaths <- genepath
          pairpaths <- pairpath
          geneshufflepaths <- geneshufflepath
          pairshufflepaths <- pairshufflepath
          
          geneweightspaths <- geneweightspath
          pairweightspaths <- pairweightspath
          pairshuffleweightspaths <- pairshuffleweightspath
          
          genebps <- genebp
          pairbps <- pairbp
          geneshufflebps <- geneshufflebp
          pairshufflebps <- pairshufflebp
          
          geneweightsbps <- geneweightsbp
          pairweightsbps <- pairweightsbp
          pairshuffleweightsbps <- pairshuffleweightsbp
          
        }else{
          genepaths <- rbind(genepaths, genepath)
          pairpaths <- rbind(pairpaths, pairpath)
          geneshufflepaths <- rbind(geneshufflepaths, geneshufflepath)
          pairshufflepaths <- rbind(pairshufflepaths, pairshufflepath)
          
          geneweightspaths <- rbind(geneweightspaths, geneweightspath)
          pairweightspaths <- rbind(pairweightspaths, pairweightspath)
          pairshuffleweightspaths <- rbind(pairshuffleweightspaths, pairshuffleweightspath)
          
          genebps <- rbind(genebps, genebp)
          pairbps <- rbind(pairbps, pairbp)
          geneshufflebps <- rbind(geneshufflebps, geneshufflebp)
          pairshufflebps <- rbind(pairshufflebps, pairshufflebp)
          
          geneweightsbps <- rbind(geneweightsbps, geneweightsbp)
          pairweightsbps <- rbind(pairweightsbps, pairweightsbps)
          pairshuffleweightsbps <- rbind(pairshuffleweightsbps, pairshuffleweightsbp)
          
        }
        
        #cat(paste0('q = ', q, '\n'))
      }
      
    }
    
  }
  
  gc()
  
  #limma within module and mediation analysis
  
  melimmares <- NULL
  
  sigg.featurelimma <- NULL
  
  mediationdat <- NULL
  
  if(returnearly == FALSE){
    
    if(length(sigmes) > 0){
      
      if(mediation == FALSE){
        cat('Start limma analysis within WGCNA modules\n')
      }else{
        cat('Start limma analysis within WGCNA modules and mediation analysis\n')
      }
      
      
      melimmares <- list()
      meplotdat <- list()
      mediationdat <- list()
      mediatorprobedat <- list()
      mediatorgenedat <- list()
      
      if(is.null(topn)){
        
        topn <- length(sigmes)
        
      }else{
        
        topn <- max(min(length(sigmes), round(topn)), 0)
        
      }
      
      l <- 1
      for(l in 1:length(sigmes)){
        
        sigmename <- sigmes[l]
        sigme <- gsub(pattern = 'ME', replacement = '', x = sigmename)
        sigme <- as.numeric(sigme)
        sigmefeatures <- names(wgcnares$numericlabels[wgcnares$numericlabels == sigme])
        
        
        if(balanceadj == TRUE){
          
          featurebalanceddats <- balanceddats
          
          for(m in 1:length(balanceddats$dats)){
            featurebalanceddats$dats[[m]] <- 
              featurebalanceddats$dats[[m]][sigmefeatures, , drop = FALSE]
          }
          
          
          featurelimmareslist <- list()
          
          iseqs <- seq(1, nround, 1)
          
          for(i in iseqs){
            
            featurelimmares <- diffsites(dats = featurebalanceddats$dats, 
                                         pddats = balanceddats$pds, 
                                         i = i, 
                                         responsevarname = responsevarname, 
                                         responselevels = responselevels, 
                                         confoundings = confoundings)
            
            featurelimmareslist[[i]] <- featurelimmares
            
            
          }
          
          
          names(featurelimmareslist) <- names(balanceddats$dats)
          
          combineres <- combinepval(limmareslist = featurelimmareslist, method = method)
          
          metafeaturelimmares <- data.frame(logFC = combineres$logfccombs, 
                                            P.Value = combineres$pcombs, 
                                            adj.P.Val = combineres$padjs, 
                                            stringsAsFactors = FALSE)
          
          featurelimmares <- metafeaturelimmares
          
          
        }else{
          
          
          featurebalanceddats <- list(dats = list(round_1 = dat[sigmefeatures, , drop = FALSE]), 
                                      pds = list(round_1 = pddat))
          
          featurelimmares <- diffsites(dats = featurebalanceddats$dats, 
                                       pddats = balanceddats$pds, 
                                       i = 1, 
                                       responsevarname = responsevarname, 
                                       responselevels = responselevels, 
                                       confoundings = confoundings)
          
        }
        
        
        if(!is.null(removereg)){
          
          if(featuretype == 'probe'){
            removedfeatures <- overlappedprobes(totalprobes = row.names(featurelimmares), 
                                                removereg = removereg, 
                                                platform = platform, 
                                                ignorestrand = ignorestrand)
          }else if(featuretype == 'gene'){
            removedfeatures <- overlappedgenes(totalgenes = row.names(featurelimmares), 
                                               removereg = removereg, 
                                               tssradius = tssradius, 
                                               ignorestrand = ignorestrand)
          }else{
            removedfeatures <- NULL
          }
          
          
          if(!is.null(removedfeatures)){
            
            if(featuretype == 'probe'){
              removedfeatures <- unique(row.names(removedfeatures))
              
              featurelimmares$remove <- FALSE
              featurelimmares$remove[row.names(featurelimmares) %in% removedfeatures] <- TRUE
              
            }else if(featuretype == 'gene'){
              
              removegenes <- idxgenenames(totalgenes = row.names(featurelimmares), 
                                          genesym = removedfeatures$genesym, 
                                          geneid = removedfeatures$geneid)
              featurelimmares$remove <- removegenes
              
            }else{
              featurelimmares$remove <- FALSE
            }
            
          }else{
            featurelimmares$remove <- FALSE
          }
          
        }
        
        
        if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val', 
                                       'pvalue', 'p-value', 'p.value')){
          pvalcolname <- 'P.Value'
        }else{
          pvalcolname <- 'adj.P.Val'
        }
        
        #otherpvalcolname <- setdiff(c('P.Value', 'adj.P.Val'), pvalcolname)
        
        featurelimmares <- featurelimmares[order(featurelimmares[,pvalcolname], 
                                                 #featurelimmares[otherpvalcolname], 
                                                 -abs(featurelimmares$logFC)), , drop = FALSE]
        
        
        if('remove' %in% names(featurelimmares)){
          
          sigg.featurelimma <- featurelimmares[(featurelimmares[,pvalcolname] < 
                                                  pvalcutoff) & 
                                                 (abs(featurelimmares$logFC) > 
                                                    absxcutoff) & 
                                                 (featurelimmares$remove == FALSE),]
          nonsigg.featurelimma <- featurelimmares[(featurelimmares[,pvalcolname] >= 
                                                     pvalcutoff) | 
                                                    (abs(featurelimmares$logFC) <= 
                                                       absxcutoff) | 
                                                    (featurelimmares$remove == TRUE),]
          
        }else{
          
          sigg.featurelimma <- featurelimmares[(featurelimmares[,pvalcolname] < 
                                                  pvalcutoff) & 
                                                 (abs(featurelimmares$logFC) > 
                                                    absxcutoff),]
          nonsigg.featurelimma <- featurelimmares[(featurelimmares[,pvalcolname] >= 
                                                     pvalcutoff) | 
                                                    (abs(featurelimmares$logFC) <= 
                                                       absxcutoff),]
          
        }
        
        if(!is.null(nrow(sigg.featurelimma))){
          
          if(nrow(sigg.featurelimma) > 0 & mediation == TRUE & l <= topn){
            
            mediators <- row.names(sigg.featurelimma)
            
            probeannores <- NULL
            geneannores <- NULL
            
            if(balanceadj == TRUE){
              n <- 1
              for(n in 1:length(balanceddats$dats)){
                
                mediatordat <- balanceddats$dats[[n]][mediators, , drop = FALSE]
                mediatordat <- t(mediatordat)
                
                otherdat <- balanceddats$pds[[n]][, c(responsevarname, confoundings), 
                                                  drop = FALSE]
                otherdat <- cbind(otherdat, wgcnares$mergedmelist[[n]][, sigmename, 
                                                                       drop = FALSE])
                row.names(otherdat) <- row.names(mediatordat)
                
                
                mediationreslist <- mediationmod(mediatordat = mediatordat, 
                                                 predictdat = otherdat, 
                                                 source = sigmename, 
                                                 target = responsevarname, 
                                                 interaction = FALSE, 
                                                 bootnum = 100, 
                                                 IPW = TRUE, 
                                                 responselevels = responselevels, 
                                                 threads = threads)
                
                
                mediationres <- orgreslist(reslist = mediationreslist, 
                                           targetname = responsevarname, 
                                           sourcename = sigmename, 
                                           anno = FALSE, 
                                           propna = FALSE)
                
                if(n == 1){
                  
                  completeres <- mediationres$completeres
                  
                }else{
                  
                  completeres <- rbind(completeres, mediationres$completeres)
                  
                }
                
                
                mediationrevreslist <- mediationmod(mediatordat = mediatordat, 
                                                    predictdat = otherdat, 
                                                    source = responsevarname, 
                                                    target = sigmename, 
                                                    interaction = FALSE, 
                                                    bootnum = 100, 
                                                    IPW = TRUE, 
                                                    sourcelevels = responselevels, 
                                                    threads = threads)
                
                
                mediationrevres <- orgreslist(reslist = mediationrevreslist, 
                                              targetname = sigmename, 
                                              sourcename = responsevarname, 
                                              anno = FALSE, 
                                              propna = FALSE)
                
                
                if(n == 1){
                  
                  completerevres <- mediationrevres$completeres
                  
                }else{
                  
                  completerevres <- rbind(completerevres, mediationrevres$completeres)
                  
                }
                
              }
              
              completeres <- combineci(cires = completeres)
              
              completeres <- mediationanno(reslines = completeres, propna = TRUE)
              
              subres <- subset(completeres, mediation == TRUE)
              
              
              completerevres <- combineci(cires = completerevres)
              
              completerevres <- mediationanno(reslines = completerevres, propna = TRUE)
              
              subrevres <- subset(completerevres, mediation == TRUE)
              
              
              if(nrow(subres) > 0 & nrow(subrevres) > 0){
                sharedmediators <- intersect(subres$mediator, subrevres$mediator)
                subres <- subset(subres, !(mediator %in% sharedmediators))
                subrevres <- subset(subrevres, !(mediator %in% sharedmediators))
              }
              
              qualifies <- rbind(subres, subrevres)
              
              
              if(!is.null(qualifies)){
                if(nrow(qualifies) > 0){
                  row.names(qualifies) <- 1:nrow(qualifies)
                  
                  mediatorannores <-  mediatorannotation(mediators = qualifies$mediator, 
                                                         platform = platform)
                  
                  
                  geneannores <- mediatorannores$geneannores
                  
                  probeannores <- mediatorannores$probeannores
                  
                }
              }
              
            }else{
              
              mediatordat <- balanceddats$dats[[1]][mediators, , drop = FALSE]
              mediatordat <- t(mediatordat)
              
              otherdat <- balanceddats$pds[[1]][, c(responsevarname, confoundings), drop = FALSE]
              otherdat <- cbind(otherdat, wgcnares$mergedmelist[, sigmename, drop = FALSE])
              row.names(otherdat) <- row.names(mediatordat)
              
              mediationreslist <- mediationmod(mediatordat = mediatordat, 
                                               predictdat = otherdat, 
                                               source = sigmename, 
                                               target = responsevarname, 
                                               interaction = FALSE, 
                                               bootnum = 100, 
                                               IPW = TRUE, 
                                               responselevels = responselevels, 
                                               threads = threads)
              
              mediationres <- orgreslist(reslist = mediationreslist, 
                                         targetname = responsevarname, 
                                         sourcename = sigmename, 
                                         propna = TRUE)
              
              subres <- mediationres$subres
              
              
              mediationrevreslist <- mediationmod(mediatordat = mediatordat, 
                                                  predictdat = otherdat, 
                                                  source = responsevarname, 
                                                  target = sigmename, 
                                                  interaction = FALSE, 
                                                  bootnum = 100, 
                                                  IPW = TRUE, 
                                                  sourcelevels = responselevels, 
                                                  threads = threads)
              
              
              mediationrevres <- orgreslist(reslist = mediationrevreslist, 
                                            targetname = sigmename, 
                                            sourcename = responsevarname, 
                                            propna = TRUE)
              
              subrevres <- mediationrevres$subres
              
              if(!is.null(subres) & !is.null(subrevres)){
                sharedmediators <- intersect(subres$mediator, subrevres$mediator)
                subres <- subset(subres, !(mediator %in% sharedmediators))
                subrevres <- subset(subrevres, !(mediator %in% sharedmediators))
              }
              
              qualifies <- rbind(subres, subrevres)
              
              if(!is.null(qualifies)){
                if(nrow(qualifies) > 0){
                  row.names(qualifies) <- 1:nrow(qualifies)
                  
                  mediatorannores <-  mediatorannotation(mediators = qualifies$mediator, 
                                                         platform = platform)
                  
                  geneannores <- mediatorannores$geneannores
                  
                  probeannores <- mediatorannores$probeannores
                  
                }
              }
              
            }
            
            mediationdat[[l]] <- qualifies
            
            if(!is.null(probeannores)){
              probeannores$ME <- sigmename
              probeannores <- unique(probeannores)
              mediatorprobedat[[l]] <- probeannores
            }
            
            if(!is.null(geneannores)){
              if(nrow(geneannores) > 0){
                
                geneannores$ME <- sigmename
                geneannores <- unique(geneannores)
                mediatorgenedat[[l]] <- geneannores
                
              }
            }
            
          }
        }
        
        
        
        if(plot == TRUE){
          
          if(tolower(pvalcolname) %in% c('pval', 'p-val', 'p.val', 
                                         'pvalue', 'p-value', 'p.value')){
            pvalcolname <- 'P.Value'
          }else{
            pvalcolname <- 'adj.P.Val'
          }
          
          if(!is.null(titleprefix)){
            featureplottitle <- paste0(titleprefix, ' ', sigmename)
          }else{
            featureplottitle <- sigmename
          }
          
          featureplotdat <- singlevolcano(limmares = featurelimmares, 
                                          betaval = isbetaval, 
                                          absxcutoff = absxcutoff, 
                                          pvalcut = pvalcutoff, 
                                          pvalcolname = pvalcolname, 
                                          confoundings = confoundings, 
                                          responsevarname = responsevarname, 
                                          titleprefix = featureplottitle, 
                                          featuretype = featuretype, 
                                          titlesize = titlesize, 
                                          textsize = textsize, 
                                          face = face, 
                                          labelnum = labelnum, 
                                          annotextsize = annotextsize)
          
          
          featureplotdat$ME <- sigmename
          meplotdat[[l]] <- featureplotdat
          #names(meplotdat)[l] <- sigmename
          
          if(!is.null(nrow(sigg.featurelimma))){
            if(nrow(sigg.featurelimma) > 0 & mediation == TRUE & l <= topn){
              if(!is.null(qualifies)){
                if(nrow(qualifies) > 0){
                  
                  qualifiesplotdat <- featurelimmares[qualifies$mediator, , drop = FALSE]
                  othersplot <- featurelimmares[setdiff(row.names(featurelimmares), 
                                                        row.names(qualifiesplotdat)), , 
                                                drop = FALSE]
                  othersplot$color <- '#C0C0C0'
                  qualifiesplotdat$color <- '#C0C0C0'
                  qualifiesplotdat$color[qualifies$target == sigmename] <- '#0000FF'
                  qualifiesplotdat$color[qualifies$source == sigmename] <- '#FF0000'
                  qualifiesplotdat <- qualifiesplotdat[rev(seq(1, nrow(qualifiesplotdat))), , 
                                                       drop = FALSE]
                  qualifiesplotdat <- rbind(othersplot, qualifiesplotdat)
                  qualifiesplotdat$Type <- 'NS'
                  qualifiesplotdat$Type[qualifiesplotdat$color == '#0000FF'] <- 'Rev'
                  qualifiesplotdat$Type[qualifiesplotdat$color == '#FF0000'] <- 'Fwd'
                  
                  qualifiesplotdat$label <- FALSE
                  
                  if(!is.null(labelnum)){
                    
                    labelnum <- as.numeric(labelnum)
                    
                    if(is.na(labelnum)){
                      labelnum <- NULL
                    }
                  }
                  
                  if(!is.null(labelnum)){
                    qualifiesplotdat$label[row.names(qualifiesplotdat) %in% 
                                             row.names(tail(subset(qualifiesplotdat, 
                                                                   Type == 'Rev'), labelnum))] <- 
                      TRUE
                    qualifiesplotdat$label[row.names(qualifiesplotdat) %in% 
                                             row.names(tail(subset(qualifiesplotdat, 
                                                                   Type == 'Fwd'), labelnum))] <- 
                      TRUE
                    
                  }
                  
                  if(!is.null(titleprefix)){
                    mediationplottitle <- paste0(titleprefix, ' mediation for ', sigmename)
                  }else{
                    mediationplottitle <- paste0('Mediation for ', sigmename)
                  }
                  
                  mediationplotdat <- singlevolcano(limmares = qualifiesplotdat, 
                                                    betaval = isbetaval, 
                                                    absxcutoff = absxcutoff, 
                                                    pvalcut = pvalcutoff, 
                                                    pvalcolname = pvalcolname, 
                                                    confoundings = confoundings, 
                                                    responsevarname = responsevarname, 
                                                    titleprefix = mediationplottitle, 
                                                    featuretype = featuretype, 
                                                    titlesize = titlesize, 
                                                    textsize = textsize, 
                                                    face = face, 
                                                    labelnum = labelnum, 
                                                    annotextsize = annotextsize)
                  
                }
              }
              
            }
          }
          
          
          
          
          
        }
        
        featurelimmares$ME <- sigmename
        melimmares[[l]] <- featurelimmares
        #names(melimmares)[l] <- sigmename
        
        
      }
      
      melimmares <- do.call(rbind, melimmares)
      mediationdat <- do.call(rbind, mediationdat)
      
      if(length(mediatorprobedat) > 0){
        mediatorprobedat <- do.call(rbind, mediatorprobedat)
      }
      
      if(length(mediatorgenedat) > 0){
        mediatorgenedat <- do.call(rbind, mediatorgenedat)
      }
      
      if(plot == TRUE){
        meplotdat <- do.call(rbind, meplotdat)
      }
      
    }
    
  }
  
  
  if(returnearly == TRUE & length(mergedmenames) > 1 & 
     noenrichcount != (length(mergedmenames) - 1)){
    
    res <- list(wgcnares = wgcnares, sftpowers = sftpowers)
    
    res$genepaths <- genepaths
    res$pairpaths <- pairpaths
    res$geneshufflepaths <- geneshufflepaths
    res$pairshufflepaths <- pairshufflepaths
    res$genebps <- genebps
    res$pairbps <- pairbps
    res$geneshufflebps <- geneshufflebps
    res$pairshufflebps <- pairshufflebps
    
    res$geneweightspaths <- geneweightspaths
    res$pairweightspaths <- pairweightspaths
    res$pairshuffleweightspaths <- pairshuffleweightspaths
    res$geneweightsbps <- geneweightsbps
    res$pairweightsbps <- pairweightsbps
    res$pairshuffleweightsbps <- pairshuffleweightsbps
    
    if(savetom == FALSE & balanceadj == FALSE){
      res$wgcnares$tom <- NULL
    }else if(savetom == FALSE & balanceadj == TRUE){
      res$wgcnares$consensustom <- NULL
    }
    
    
    return(res)
    
  }
  
  
  if(plot == FALSE){
    
    res <- list(limmares = limmares, wgcnares = wgcnares, sftpowers = sftpowers, 
                melimmares = melimmares)
    
  }else{
    
    res <- list(limmares = plotdat, wgcnares = wgcnares, sftpowers = sftpowers, 
                melimmares = meplotdat)
    
  }
  
  
  
  if(!is.null(mediationdat)){
    
    if(!is.null(nrow(mediationdat))){
      
      if(nrow(mediationdat) > 0 & mediation == TRUE){
        
        res$mediationres <- mediationdat
        
        if(length(mediatorprobedat) > 0){
          res$mediatorprobedat <- mediatorprobedat
        }
        if(length(mediatorgenedat) > 0){
          res$mediatorgenedat <- mediatorgenedat
        }
        
      }
    }
    
  }
  
  
  if(featuretype == 'gene' & topoenrichment == TRUE & length(mergedmenames) > 1 & 
     noenrichcount != (length(mergedmenames) - 1)){
    
    res$genepaths <- genepaths
    res$pairpaths <- pairpaths
    res$geneshufflepaths <- geneshufflepaths
    res$pairshufflepaths <- pairshufflepaths
    res$genebps <- genebps
    res$pairbps <- pairbps
    res$geneshufflebps <- geneshufflebps
    res$pairshufflebps <- pairshufflebps
    
    res$geneweightspaths <- geneweightspaths
    res$pairweightspaths <- pairweightspaths
    res$pairshuffleweightspaths <- pairshuffleweightspaths
    res$geneweightsbps <- geneweightsbps
    res$pairweightsbps <- pairweightsbps
    res$pairshuffleweightsbps <- pairshuffleweightsbps
    
  }
  
  
  if(savetom == FALSE & balanceadj == FALSE){
    res$wgcnares$tom <- NULL
  }else if(savetom == FALSE & balanceadj == TRUE){
    res$wgcnares$consensustom <- NULL
  }
  
  
  return(res)
  
  
}

