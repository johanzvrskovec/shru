#Based on the amazing work by Grotzinger, A. D. et al. Nat. Hum. Behav. 3, 513–525 (2019) and Bulik-Sullivan, B. K. et al. Nat. Genet. 47, 291–295 (2015).
#LDSC++ is branched of the multivariate ldsc function in GenomicSEM: https://github.com/GenomicSEM/GenomicSEM
#LDSC++ is a new take on ld score regression containing exprerimental implementations (variable block-count sampling, association statistics Chi-square weighting, GWAS sumstats imputation weighting, correlation bias correction, etc)
#Forked and modified by Johan Zvrskovec (2020)

ldscpp <- function(
                     traits, 
                     sample.prev,
                     population.prev, 
                     ld = NULL,
                     wld = NULL,
                     filepathLD = NULL,
                     filepathWLD = NULL,
                     df_ld = NULL, #ldscores as a ready dataframe
                     filepathVariantsAdditional = NULL,
                     readyMergedTraits = NULL,
                     referenceAncestrySetting = "MIX",
                     ancestrySetting=NULL, #ancestry setting, list with entries for each dataset
                     trait.names = NULL,
                     cap.NEFF = NULL, #list per dataset, T/F, cap NEFF with maximum N or not (not recommended for 'real' NEFF (from sub-cohorts)), DEFAULT = TRUE
                     sep_weights = FALSE,
                     chr = 23,
                     n.blocks = NA,
                     blocksizeM = 20000,
                     blocksizeCM = NA, #we have tested with 10 
                     minBlocksizeCM = NA, #5 as a good setting if anything
                     minBlocksizeM = 5000, #from ~ LDSC common blocksize at 200 blocks and HM3 reference
                     ldsc.log = NULL,
                     select=FALSE,
                     GC=NULL, #vector with correction denominators
                     filter.chisq.max = NA,
                     filter.chisq.min = NA, #0.01579077, #qchisq(p = 0.1,df = 1)
                     filter.zz.min = NA, #0.01579077, #qchisq(p = 0.1,df = 1) #this is not recommended anymore
                     filter.info = NA, #this is the common more than (mt) filter, .6, #better to use the INFO -weighting below
                     filter.info.lt = NA, #this is a less than (lt) filter for producing results stratified by INFO. genotyped variants are considered to have INFO=1, and will be filered away with this filter !=NA.
                     filter.simpute.fraction = NA,
                     filter.sinfo.imputed = 0.9, #filter for supermunge LD-IMP imputation quality score, for fully imputed pairwise variant combinations
                     filter.sinfo = 0.9, #filter for supermunge LD-IMP imputation quality score, for partially imputed pairwise variant combinations
                     filter.ldimp.k =NA, #filter for supermunge LD-IMP imputation contributing variants
                     resetImputedN=F, #reset N for imputed variants if the N has been set to reflect GWAS sumstats imputation quality (eg. SSimp)
                     filter.maf=0.01,
                     filter.mhc=NA_integer_, #can be either 37 or 38 for filtering the MHC region according to either grch37 or grch38
                     filter.region.df=NULL,
                     filter.gwasssimputedonly=F,
                     N=NULL,
                     M.mode="maf", #maf/real - maf: the size of the LD-score library with MAF>0.05, real: use the actual number of variants included in the trait/LD-score overlap
                     forceN=FALSE,
                     leave.chr=NULL, #vector with chromosome numbers to leave out of the analysis
                     limited=FALSE, #only run for the first outer loop over traits, for pairwise comparisons with the first specified trait
                     nThreads=6,
                     preweight=T, #perform pre-weighting - correct for 1.) chi-square statistics being correlated with the LD scores and 2.) chi-square statistics heteroskedasticity
                     preweight.alternativeCorrelationCorrection=T, #Use an extra 1 + in the denominator for the correlation correction weight.
                     preweight.ChiSquare=F, #perform pre-weighting of association statistics in ChiSquare format - needs preweight=T
                     preweight.INFO=T, #perform pre-weighting of INFO; additionally genotype imputed effects as part of pre-weighting - needs preweight=T
                     preweight.SINFO=T, #perform pre-weighting of SINFO; additionally weight GWAS sumstat imputed effects as part of pre-weighting - needs preweight=T
                     correctAttenuationBias = F, #perform attenuation bias correction using the product moment correlation coefficient method
                     attenuationBiasCorrectionFactorExponent = 1/4, #Exponent of the attenuation bias correction factor to control how much correction to apply. 1 = full correction. 1/2 is the square root of the factor etc.
                     doubleRegressionRoutine = F, #This feature is controversial and needs more testing. Removed T as default. Use with care! Update: Original LDSC does not seem to perform double regressions. See irwls.py and class method wls().
                     reinflateImputedGWAS=T,
                     reinflateImputedGWAS.inflationFactorExponent = 1/2, #Exponent of the inflation factor before applying it for re-inflating LD-IMP imputed variants
                     resamplingMethod="vbcs", #vbcs/jn/bs(not implemented yet), vbcs=variable block-count sampling, jn=block jackknife
                     k.folds=0, #folds to prepare for later cross validation, set to > 0 to run additional folds using the train ratio to set the size of the train subsets
                     trainRatio=0.8, #train to total ratio for the fold subsets
                     exportReadyMergedTraits=F, #export merged datasets as a list, for multiple runs
                     force.M = NULL, #set the M value (number of variants in the LD score library)
                     covariance_std_value_lower_limit = 0.05, #Mod addition - lower limit of genetic covariance to allow for the denominator in standardising covG S.E.
                     cov.SE.p.liab.test.cDivS=0.000485594, #calculated variant 1 median of full LDSC++ S.SE^2
                     verbose = F #T -> print cell-wise results
                     
) {
  
  
  # #defaults - for testing
  # ld = NULL
  # wld = NULL
  # filepathLD = NULL
  # filepathWLD = NULL
  # df_ld = NULL
  # filepathVariantsAdditional = NULL
  # readyMergedTraits = NULL
  # referenceAncestrySetting = "MIX"
  # ancestrySetting = NULL
  # trait.names = NULL
  # cap.NEFF = NULL
  # sep_weights = FALSE
  # chr = 23
  # n.blocks = NA
  # blocksizeM = 20000
  # blocksizeCM = NA
  # minBlocksizeCM = NA
  # minBlocksizeM = 5000
  # ldsc.log = NULL
  # select=FALSE
  # GC=NULL #vector with correction denominators
  # filter.chisq.max = NA
  # filter.chisq.min = NA
  # filter.zz.min = NA #0.01579077
  # filter.info = NA
  # filter.info.lt = NA
  # filter.simpute.fraction = NA
  # filter.sinfo.imputed = 0.9 #filter for supermunge LD-IMP imputation quality score, for fully imputed pairwise variant combinations
  # filter.sinfo = 0.9 #filter for supermunge LD-IMP imputation quality score, for partially imputed pairwise variant combinations
  # filter.ldimp.k =NA #filter for supermunge LD-IMP imputation contributing variants
  # resetImputedN = F
  # ldimp.adjustN = F
  # filter.maf=0.01
  # filter.mhc=NA_integer_ #can be either 37 or 38 for filtering the MHC region according to either grch37 or grch38
  # filter.region.df=NULL
  # filter.gwasssimputedonly=F
  # N=NULL
  # M.mode="maf"
  # forceN=FALSE
  # leave.chr=NULL #vector with chromosome numbers to leave out of the analysis
  # limited=FALSE #only run for the first outer loop over traits, for pairwise comparisons with the first specified trait
  # nThreads=6
  # preweight=T
  # preweight.alternativeCorrelationCorrection=T
  # preweight.ChiSquare=F
  # preweight.INFO=T
  # preweight.SINFO=T
  # correctAttenuationBias = F
  # attenuationBiasCorrectionFactorExponent = 1/3
  # doubleRegressionRoutine = F
  # reinflateImputedGWAS=T
  # reinflateImputedGWAS.inflationFactorExponent = 1/2 #Exponent of the inflation factor before applying it for re-inflating LD-IMP imputed variants
  # resamplingMethod="cv"
  # k.folds=0 #folds to prepare for later cross validation
  # trainRatio=0.8 #train to total ratio for the fold subsets
  # exportReadyMergedTraits=F
  # force.M = NULL
  # covariance_std_value_lower_limit = 0.05
  # cov.SE.p.liab.test.cDivS=0.000485594
  # verbose = F
  
  #test for dr dev
  # traits = gwas.meta.sel$filePath
  # sample.prev = gwas.meta.sel$samplePrevalence
  # population.prev = gwas.meta.sel$populationPrevalence
  # ld = ld_path
  # trait.names = gwas.meta.sel$code
  # filter.maf = 0.01
  # resamplingMethod = "vbcs"
  # doubleRegressionRoutine = T
  
  # traits = gwas.meta.sel$filePath
  # sample.prev = gwas.meta.sel$samplePrevalence
  # population.prev = gwas.meta.sel$populationPrevalence
  # ld = file.path(folderpath.data,"ld_scores","w_ld.GLAD_EDGI_NBR.keep.b38.gcta")
  # trait.names = gwas.meta.sel$code
  # #no info filter!
  # filter.maf = 0.01
  # resamplingMethod = "vbcs"
  # doubleRegressionRoutine = T

  # #small test
  # traits = p$sumstats[c("SMRV01","SMRV02"),]$mungedpath.supermunge.1kg.orig.unfiltered
  # sample.prev =  p$sumstats[c("SMRV01","SMRV02"),]$samplePrevalence
  # population.prev = p$sumstats[c("SMRV01","SMRV02"),]$populationPrevalence
  # trait.names = p$sumstats[c("SMRV01","SMRV02"),]$code
  # #n.blocks = 200
  # ld = p$folderpath.data.mvLDSC.ld.GLADp
  # ldsc.log = p$setup.code.date
  # preweight.alternativeCorrelationCorrection = F
  # preweight.ChiSquare = F
  # correctAttenuationBias = F
  # doubleRegressionRoutine = F
  # preweight.INFO=F
  # resamplingMethod="vbcs"
  # blocksizeCM = NA #inactivate CM window def
  # filter.info = 0.6
  # filter.maf = 0.01
  # filter.zz.min = 0
  
  # #large/small test
  # traits = p$sumstats[c("BIPO02","BODY14"),]$mungedpath.supermunge.1kg.orig.unfiltered
  # sample.prev =  p$sumstats[c("BIPO02","BODY14"),]$samplePrevalence
  # population.prev = p$sumstats[c("BIPO02","BODY14"),]$populationPrevalence
  # trait.names = p$sumstats[c("BIPO02","BODY14"),]$code
  # #ld = p$folderpath.data.mvLDSC.ld.hm3
  # filepathLD = p$filepath.SNPReference.1kg
  # #filepathVariantsAdditional = p$filepath.SNPReference.1kg
  # #filter.region.df = p$highld_b38
  # N = p$sumstats[c("BIPO02","BODY14"),]$n_total
  # ldsc.log = p$setup.code.date
  # #ancestrySetting = "EUR"
  # #n.blocks = 200
  # filter.zz.min = 0
  # preweight.ChiSquare = F
  # correctAttenuationBias = F
  # doubleRegressionRoutine = F
  # preweight.INFO=F
  # resamplingMethod="vbcs"
  # blocksizeCM = NA
  # filter.info = 0.6
  # filter.info.lt = 0.8
  # filter.maf = 0.01
  # force.M = 7184778
  
  #sim test
  # traits = p$sumstats.sel.sim$mungedpath.ldsc.hm3
  # #traits = p$sumstats.sel.sim$mungedpath.supermunge.hm3.unfiltered, #do we use the ldsc=munged or the unfiltered supermunged?
  # sample.prev =  p$sumstats.sel.sim$samplePrevalence.balanced
  # population.prev = p$sumstats.sel.sim$populationPrevalence
  # trait.names = p$sumstats.sel.sim$code
  # ld = p$folderpath.data.mvLDSC.ld.hm3
  # n.blocks = 200
  # ldsc.log = p$setup.code.date
  # preweight.alternativeCorrelationCorrection = F
  # preweight.ChiSquare = F
  # correctAttenuationBias = F
  # doubleRegressionRoutine = F
  # preweight.INFO=F
  # resamplingMethod="jn"
  # blocksizeCM = NA #inactivate CM window def
  # filter.info = 0.6
  # filter.maf = 0.01
  # filter.zz.min = 0
  
  cat("\n\n\nLDSC++\t\tSHRU package version 1.4.1\n") #UPDATE DISPLAYED VERSION HERE!!!!
  #this is not written in the log btw
  
  LOG <- function(..., print = TRUE) {
    msg <- paste0(...)
    if (print) print(msg)
    cat(msg, file = log.file, sep = "\n", append = TRUE)
  }
  
  time <- proc.time()
  
  begin.time <- Sys.time()
  
  # Dimensions
  if(length(traits)>0) n.traits <- length(traits)
  if(length(readyMergedTraits)>0) n.traits <- length(readyMergedTraits$all_y)
  
  n.V <- n.traits * (n.traits + 1) / 2
  
  #mod addition - set trait names to always have a value
  if(is.null(trait.names)) trait.names <- paste0("V", 1:n.traits)
  
  #mod addition - set NEFF cap setting to always have a value
  if(is.null(cap.NEFF)) cap.NEFF <- rep(T,n.traits)
  
  if(is.null(ancestrySetting)) {ancestrySetting<-list("MIX")}
  if(length(ancestrySetting)<n.traits) ancestrySetting <- shru::padListRight(ancestrySetting,padding = ancestrySetting[1],targetLength = n.traits)
  
  #for backwards compatibility
  if(resamplingMethod=='cv') resamplingMethod<-'vbcs'
  
  if(is.null(ldsc.log)){
    logtraits<-gsub(".*/","",traits)
    log2<-paste(logtraits,collapse="_")
    if(object.size(log2) > 200){
      log2<-substr(log2,1,100)
    }
    log.file <- file(paste0(log2, "_ldsc.log"),open="wt")
  }else{log.file<-file(paste0(ldsc.log, "_ldsc.log"),open="wt")}
  
  LOG("Multivariate ld-score regression of ", length(traits), " traits ",
      "(", paste(traits, collapse = " "), ")", " began at: ", begin.time)
  
  LOG("resamplingMethod=",resamplingMethod)
  LOG("filter.maf=",filter.maf)
  LOG("filter.info=",filter.info)
  LOG("filter.info.lt=",filter.info.lt)
  LOG("filter.chisq.max=",filter.chisq.max)
  LOG("filter.chisq.min=",filter.chisq.min)
  LOG("filter.simpute.fraction=",filter.simpute.fraction)
  LOG("filter.sinfo.imputed=",filter.sinfo.imputed)
  LOG("filter.sinfo=",filter.sinfo)
  LOG("filter.ldimp.k=",filter.ldimp.k)
  LOG("doubleRegressionRoutine=",doubleRegressionRoutine)
  LOG("preweight.ChiSquare=",preweight.ChiSquare)
  LOG("preweight.INFO=",preweight.INFO)
  LOG("preweight.SINFO=",preweight.SINFO)
  LOG("correctAttenuationBias=",correctAttenuationBias)
  LOG("attenuationBiasCorrectionFactorExponent=",attenuationBiasCorrectionFactorExponent)
  LOG("reinflateImputedGWAS=",reinflateImputedGWAS)
  LOG("reinflateImputedGWAS.inflationFactorExponent=",reinflateImputedGWAS.inflationFactorExponent)

  
  if(select == "ODD" | select == "EVEN"){
    odd<-seq(1,chr,2)
    even<-seq(2,chr,2)
  }
  
  
  
  if(!(is.null(trait.names))){
    check_names<-stringr::str_detect(trait.names, "-")
    if(any(check_names))
      warning("Your trait names specified include mathematical arguments (e.g., + or -) that will be misread by lavaan. Please rename the traits using the trait.names argument.")
  }
  
  if(length(traits)==1)
    warning("Our version of ldsc requires 2 or more traits. Please include an additional trait.")
  
  
  # Storage:
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  cov.pos <- matrix(NA,nrow=n.traits,ncol=n.traits) #this is the scaled covariance
  cov.pos.p <- matrix(NA,nrow=n.traits,ncol=n.traits)
  cov.neg <- matrix(NA,nrow=n.traits,ncol=n.traits) #this is the scaled covariance
  cov.neg.p <- matrix(NA,nrow=n.traits,ncol=n.traits)
  cov <- matrix(NA,nrow=n.traits,ncol=n.traits) #this is the scaled covariance
  cov.p <- matrix(NA,nrow=n.traits,ncol=n.traits)
  cov.signed <- matrix(NA,nrow=n.traits,ncol=n.traits)
  cov.p.signed <- matrix(NA,nrow=n.traits,ncol=n.traits)
  cov.unsigned <- matrix(NA,nrow=n.traits,ncol=n.traits)
  cov.p.unsigned <- matrix(NA,nrow=n.traits,ncol=n.traits)
  cov.blocks <- matrix(NA,nrow=n.traits,ncol=n.traits)
  #V.hold <- matrix(NA,nrow=n.blocks,ncol=n.V)
  lV <- c() #new storage for results of varying blocksize (unscaled!!)
  lV.pos <- c()
  lV.neg <- c()
  lV.signed <- c()
  lV.unsigned <- c()
  Liab.S <- rep(1, n.traits)
  #Liab.S.unsigned <- rep(1, n.traits)
  I.pos <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.neg <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.signed <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.unsigned <- matrix(NA,nrow=n.traits,ncol=n.traits)
  #I.unsigned <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.pos.SE <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.neg.SE <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.SE <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.signed.SE <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.unsigned.SE <- matrix(NA,nrow=n.traits,ncol=n.traits)
  #mod addition - dataset statistics
  impstats <- as.data.frame(matrix(data = NA,nrow = 0,ncol = 0))
  #cellstats matrices
  mCellstats.imputed.m <- matrix(NA,nrow=n.traits,ncol=n.traits)
  mCellstats.imputed.full.m <- matrix(NA,nrow=n.traits,ncol=n.traits)
  mCellstats.imputed.partial.m <- matrix(NA,nrow=n.traits,ncol=n.traits)
  mCellstats.nonimputed.m <- matrix(NA,nrow=n.traits,ncol=n.traits)
  LD.SS<-c()
  nVar.per.block<-c()
  
  
  
  if(length(readyMergedTraits)>0){
    all_y <- readyMergedTraits$all_y
    impstats <- readyMergedTraits$impstats
    median_all_y <- readyMergedTraits$median_all_y
    mgc <- readyMergedTraits$mgc
    M.tot <- readyMergedTraits$M.tot
    rm(readyMergedTraits)
    
  } else {
    #########  READ LD SCORES:
    LOG("Reading in LD scores")
    x<-NULL
    if(!is.null(df_ld)){
      x<-df_ld
      rm(df_ld)
    } else if(!is.null(filepathLD)){
      x <- suppressMessages(fread(file = filepathLD, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T))
      #x <- suppressMessages(fread(file = filepathLD, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T, nrows = 20000000))
    } else {
      #mod addition - read all available ld-scores in folder - use data.table
      x <- do.call("rbind", lapply(list.files(path = ld, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
        suppressMessages(fread(file = file.path(ld, i), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T))
      }))
    }
    
    # if(select == FALSE){
    # x <- do.call("rbind", lapply(1:chr, function(i) {
    #   suppressMessages(read_delim(
    #     file.path(ld, paste0(i, ".l2.ldscore.gz")),
    #     delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    # }))
    # }
    # 
    # if(select == "ODD"){
    #   x <- do.call("rbind", lapply(odd, function(i) {
    #     suppressMessages(read_delim(
    #       file.path(ld, paste0(i, ".l2.ldscore.gz")),
    #       delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    #   }))
    # }
    # 
    # if(select == "EVEN"){
    #   x <- do.call("rbind", lapply(even, function(i) {
    #     suppressMessages(read_delim(
    #       file.path(ld, paste0(i, ".l2.ldscore.gz")),
    #       delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    #   }))
    # }
    # 
    # if(is.numeric(select)){
    #   x <- do.call("rbind", lapply(select, function(i) {
    #     suppressMessages(read_delim(
    #       file.path(ld, paste0(i, ".l2.ldscore.gz")),
    #       delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    #   }))
    # }
    # 
    #mod addition, correct SNP to lower case in case it is upper case
    x$SNP<-tolower(x$SNP)
    #x$CM <- NULL
    #x$MAF <- NULL
    
    
    #mod addition - select standard LD score if multiple - from supermunge
    #old setting
    # if(ancestrySetting!="ANY" & !any(colnames(x)=="L2")){
    #   if(any(colnames(x)==paste0("L2.",ancestrySetting))){
    #     sAncestryL2<-c(paste0("L2.",ancestrySetting))
    #     x$L2<-x[,..sAncestryL2]
    #   }
    # }
    
    #mod addition - keep valid SNP ID's, LD-score values
    x<-x[!is.na(SNP),] #& is.finite(L2) & L2!=0
    
    setkeyv(x, cols = c("SNP"))
    
    print("Head of LD scores")
    print(x)
    
    ######### READ weights:
    w<-NULL  
    if(!is.null(filepathWLD)){
      w <- suppressMessages(fread(file = filepathWLD, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T))
    } else if(!is.null(wld)){
      w <- do.call("rbind", lapply(list.files(path = wld, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
        suppressMessages(fread(file = file.path(wld, i), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T))
      }))
    }
      
    if(!is.null(w)){
      w<-w[!is.na(w$SNP),]
      
      setkeyv(w, cols = c("SNP"))
      
      # if(select == FALSE){
      # w <- do.call("rbind", lapply(1:chr, function(i) {
      #   suppressMessages(read_delim(
      #     file.path(wld, paste0(i, ".l2.ldscore.gz")),
      #     delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      # }))
      # }
      # if(select == "EVEN"){
      #   w <- do.call("rbind", lapply(even, function(i) {
      #     suppressMessages(read_delim(
      #       file.path(wld, paste0(i, ".l2.ldscore.gz")),
      #       delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      #   }))
      # }
      # if(select == "ODD"){
      #   w <- do.call("rbind", lapply(even, function(i) {
      #     suppressMessages(read_delim(
      #       file.path(wld, paste0(i, ".l2.ldscore.gz")),
      #       delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      #   }))
      # }
      #   if(is.numeric(select)){
      #     w <- do.call("rbind", lapply(select, function(i) {
      #       suppressMessages(read_delim(
      #         file.path(wld, paste0(i, ".l2.ldscore.gz")),
      #         delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      #     }))
      # }
      
      #mod addition, correct SNP to lower case in case it is upper case
      w$SNP<-tolower(w$SNP)
      #w$CM <- NULL
      #w$MAF <- NULL
      
      colnames(w)[colnames(w)=="L2"] <- "wLD"
      
    }
    
    
    #merge in additional columns - rely on the SNP column
    if(!is.null(filepathVariantsAdditional)){
      x.additional <- suppressMessages(fread(file = filepathVariantsAdditional, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T))
      setkeyv(x, cols = c("SNP"))
      setkeyv(x.additional, cols = c("SNP"))
      x[x.additional,on=c('SNP'), c("CHR","BP","MAF","CM"):=list(i.CHR,i.BP,i.MAF,i.CM)]
      #merged<-merge(x = x, y = x.additional, by = "SNP", all.x = T)
      cat("Updated LD scores with CHR,BP,MAF,CM from specified file by SNP:",nrow(x[is.finite(MAF),]),",",nrow(x[is.finite(CM),]),"/",nrow(x))
      rm(x.additional)
    }
    
    
    x.colnames <- shru::stdGwasColumnNames(columnNames = colnames(x),missingEssentialColumnsStop = NULL, ancestrySetting = referenceAncestrySetting, warnings=F) #needs the shru-package
    x.copy <- x
    colnames(x.copy) <- x.colnames$std
    M.tot<-0
    if(!is.null(force.M)){
      M.tot <- force.M
      LOG("Total M (forced from argument): ",M.tot)
    } else if (!is.null(ld)){
      #fallback read from file
      m <- do.call("rbind", lapply(list.files(path = ld, pattern = "\\.l2\\.M_5_50$"), function(i) {
        suppressMessages(fread(file = file.path(ld, i), na.strings =c(".",NA,"NA",""), encoding = "UTF-8", fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, header = F))
      }))
      M.tot <- sum(m)
      rm(m)
      LOG("Total M (read from ld score files): ",M.tot)
    }
    if(any(colnames(x.copy)=="FRQ") & M.tot==0) {
      M.tot<-nrow(x.copy[FRQ>0.05 & FRQ<0.95,]) #may be pushed downward from 0.05 as we have larger samples now compared to when original ldsc was published?
      LOG("Total M (from LD reference, MAF>0.05): ",M.tot," using ancestry setting: ",referenceAncestrySetting)
    }
    if(M.tot==0){
      stop("Please provide a value for the (total) M!")
    }
    rm(x.copy)
    
    ### READ ALL CHI2 + MERGE WITH LDSC FILES
  
    iTrait <- 1
    all_y <- c()
    for(chi1 in traits) {
      #chi1<-traits[iTrait]
      ## READ chi2
      ## mod change - use fread - using data.table!
      y1 <- fread(file = chi1, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = F)
      LOG("Read in summary statistics [", iTrait, "/", n.traits, "] from: ", chi1)
      LOG("Columns: ")
      LOG(paste(colnames(y1)," "))
      LOG("Dataset specific ancestry setting: ",ancestrySetting[iTrait])
      LOG("Cap NEFF: ",cap.NEFF[[iTrait]])
      
      
      if(any(colnames(y1)=="INFO")) impstats[iTrait,c("hasinfo")] <- c(T)
      if(any(colnames(y1)=="SINFO")) impstats[iTrait,c("hassinfo")] <- c(T)
      
      
      if(any(colnames(y1)=="LD_IMP")) y1[,SINFO:=LD_IMP][,LD_IMP:=NULL] #TEMPORARY - rename old LD_IMP columns
      
      #mod addition - statistics before filters
      impstats[iTrait,c("m","medianChiSquare.tot")] <- c(nrow(y1), median((y1[is.finite(Z^2),]$Z^2),na.rm=T))
      
      ##mod addition, check existence of FRQ, otherwise, use MAF
      if(!any(colnames(y1)=="FRQ") & any(colnames(y1)=="MAF")) y1$FRQ<-y1$MAF
      if(!any(colnames(y1)=="NEFF") & any(colnames(y1)=="NEF")) y1[,NEFF:=NEF][,NEF:=NULL]
      
      
      #mod addition
      #added filters to accommodate strangely parsed sumstats
      if(any(colnames(y1)=="BETA") & any(colnames(y1)=="SE")) y1<-y1[!is.na(SNP) & is.finite(BETA) & is.finite(SE),]
      if(any(colnames(y1)=="Z")) y1<-y1[is.finite(Z),]
      
      #mod addition, harmonise character case
      y1$SNP<-tolower(as.character(y1$SNP))
      if(any(colnames(y1)=="A1")) y1$A1<-toupper(as.character(y1$A1))
      if(any(colnames(y1)=="A2")) y1$A2<-toupper(as.character(y1$A2))
      y1$Z<-as.numeric(y1$Z)
      #mod addition - use FRQ - frequency of A1 rather than MAF. use MAF as FRQ in case there is no FRQ.
      if('FRQ' %in% names(y1)) y1$FRQ<-as.numeric(y1$FRQ)
      if('MAF' %in% names(y1)) y1$MAF<-as.numeric(y1$MAF)
      if('INFO' %in% names(y1)) y1$INFO<-as.numeric(y1$INFO)
      if('SINFO' %in% names(y1)) y1$SINFO<-as.numeric(y1$SINFO)
      
      #mod addition - use NEF as N
      if(!any(colnames(y1)=="N") & any(colnames(y1)=="NEFF")) y1$N<-y1$NEFF
      
      #mod addition - checks, use explicit new N value
      LOG("N - median, min, max: ",median(y1$N, na.rm = T),", ",min(y1$N, na.rm = T),", ", max(y1$N, na.rm = T))
      if(!is.null(N) & length(N)>=length(traits)) {
        LOG("Specified N: ",N[[iTrait]])
        if(forceN) { 
          #LOG("Set explicit N=",N[[iTrait]])
          if(
            abs((N[[iTrait]] - median(y1$N, na.rm = T))/median(y1$N, na.rm = T))>0.05
            |
            (N[[iTrait]]-min(y1$N, na.rm = T))/N[[iTrait]]>0.05
            |
            (N[[iTrait]]-max(y1$N, na.rm = T))/N[[iTrait]]>0.05
          ) {
            warning("Large (>5%) N discrepancies found between provided and existing N!")
            LOG("Large (>5%) N discrepancies found between provided and existing N!")
          }
          y1$N<-N[[iTrait]]
          
          
        } else if(!is.na(N[[iTrait]])){
          cond<-is.na(y1$N) | N[[iTrait]] < y1$N
          y1$N<-ifelse(cond,N[[iTrait]],y1$N)
          if(sum(cond)>0) LOG("Added explicit N=",N[[iTrait]], " to ",sum(cond)," SNPs where it was missing or larger than the specified value.")
        }
        
        LOG("New N - median, min, max: ",median(y1$N, na.rm = T),", ",min(y1$N, na.rm = T),", ", max(y1$N, na.rm = T))
      }
      
      #use NEFF_CAPPED rather than N for binary traits if available
      if(cap.NEFF[[iTrait]] & !is.na(sample.prev[[iTrait]]) & !is.na(population.prev[[iTrait]]) & any(colnames(y1)=="NEFF")){
        if(any(colnames(y1)=="N")){
          LOG("Attention: Using CAPPED NEFF for N and setting the sample prevalence to a balanced design (0.5).")
          maxN<-max(y1[is.finite(N),]$N,na.rm = T)
          y1[,NEFF_CAPPED:=shru::clipValues(NEFF,max = 1.1*eval(maxN), min = 0.5*eval(maxN))]
          y1[,N:=NEFF_CAPPED]
        } else {
          LOG("Attention: Using NEFF (NOT CAPPED) for N and setting the sample prevalence to a balanced design (0.5).")
          y1[,N:=NEFF]
        }
        if(any(colnames(y1)=="NEFF")) y1[,NEFF:=NULL]
        if(any(colnames(y1)=="NEFF_CAPPED")) y1[,NEFF_CAPPED:=NULL]
        if(any(colnames(y1)=="N_CAS")) y1[,N_CAS:=NULL]
        if(any(colnames(y1)=="N_CON")) y1[,N_CON:=NULL]
        sample.prev[[iTrait]]<-0.5
        LOG("New N - median, min, max: ",median(y1$N, na.rm = T),", ",min(y1$N, na.rm = T),", ", max(y1$N, na.rm = T))
      }
      
      #LD-IMP actions - temporary
      #compute mean chi-square for genotyped variants
      if(any(colnames(y1)=="SINFO")) {
          if(any(colnames(y1)=="BETA.I")) {
            meanChi2.genotyped <- mean(y1[!(is.finite(BETA.I) & BETA.I==BETA) & is.finite(Z^2),]$Z^2,na.rm=T) 
          } else{ 
              meanChi2.genotyped <- mean(y1[!is.finite(SINFO) & is.finite(Z^2),]$Z^2,na.rm=T)
              }
        } else {
          meanChi2.genotyped <- mean(y1[is.finite(Z^2),]$Z^2,na.rm=T) 
      }
      if(any(colnames(y1)=="SINFO") & any(colnames(y1)=="BETA.I")){
        
        #fix imputation quality scores == 1 (most probable not imputed data - for SSimp)
        if(nrow(y1[SINFO==1,])>0.05*nrow(y1[is.finite(SINFO),])) y1[SINFO==1,SINFO:=NA_real_]
        
        if(filter.gwasssimputedonly){
          #filter out imputed only variants, even if they come from validation runs
          LOG("From ",nrow(y1), " variants, excluding ",nrow(y1[!is.finite(SINFO),]), " that are non-imputed.")
          y1<-y1[is.finite(SINFO),]
          y1[,EFFECT:=BETA.I][,SE:=SE.I][,Z:=EFFECT/SE]
          y1[,P:=2*pnorm(q = abs(Z),mean = 0, sd = 1, lower.tail = F)]
          y1[,IMPUTED:=T]
        } else if(any(colnames(y1)=="BETA.I")) {
          y1[,IMPUTED:=(is.finite(BETA.I) & BETA.I==BETA)]
        } else {
          y1[,IMPUTED:=(is.finite(SINFO))]
        }
        
        # re-caclulate SINFO
        # if(any(colnames(y1)=="LDIMP.L2.SUM") & any(colnames(y1)=="LDIMP.K")) {
        #   y1[,SINFO:=NA_real_]
        #   y1[IMPUTED==T & is.finite(L2_REF),SINFO:=pnorm(q = LDIMP.K*L2_REF/(LDIMP.L2.SUM), mean = 1, sd = 1)]
        # }
        
        #remove validation SINFO values
        y1[IMPUTED==F,SINFO:=NA_real_]
        
        LOG("SINFO stats")
        if(any(colnames(y1)=="LDIMP.K")) LOG("K quartiles:", quantile(y1$LDIMP.K, na.rm = T))
        LOG("SINFO quality quartiles:", quantile(y1$SINFO, na.rm = T))
        LOG("SINFO quality mean: ",mean(y1$SINFO, na.rm = T))
      }
      
      #imputation statistics
      if(any(colnames(y1)=="SINFO")){
        impstats[iTrait,c("m.imputed","medianChiSquare.nonimputed","medianChiSquare.imputed")] <- c(nrow(y1[is.finite(SINFO),]), median((y1[is.na(SINFO) & is.finite(Z^2),]$Z^2),na.rm=T), median((y1[is.finite(SINFO) & is.finite(Z^2),]$Z^2),na.rm=T))
      } else {
        impstats[iTrait,c("m.imputed","medianChiSquare.nonimputed","medianChiSquare.imputed")] <- c(0,NA_real_,NA_real_)
      }
      
      
      
      
      ##mod addition - remove MHC region based on position
      #references
      #https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
      #https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
      if(!is.na(filter.mhc)){
        if(any(colnames(y1)=="CHR") & any(colnames(y1)=="BP")){
          if(filter.mhc==37) {
            rm <- (!is.na(y1$CHR) & !is.na(y1$BP) & y1$CHR=="6" & y1$BP>=28477797 & y1$BP<=33448354)
            y1 <- y1[!rm, ]
            LOG("Removing ", sum(rm), " SNPs in the GRCh",filter.mhc," MHC; ", nrow(y1), " remain")
          } else if (filter.mhc==38) {
            rm <- (!is.na(y1$CHR) & !is.na(y1$BP) & y1$CHR=="6" & y1$BP>=28510120 & y1$BP<=33480577)
            y1 <- y1[!rm, ]
            LOG("Removing ", sum(rm), " SNPs in the GRCh",filter.mhc," MHC; ", nrow(y1), " remain")
          } else {
            LOG("Warning: Invalid assembly version provided - no filtering of the MHC was done!")
          }
        } else {
          LOG("Warning: No chromosome or base-pair position information available - no filtering of the MHC was done!")
        }
      }
      
      #mod addition
      #remove custom regions according to specified dataframe
      if(!is.null(filter.region.df)){
        if(any(colnames(y1)=="CHR") & any(colnames(y1)=="BP")){
          if(!(any(colnames(filter.region.df)=="CHR") & any(colnames(filter.region.df)=="BP" & any(colnames(filter.region.df)=="BP2")))) stop("Dataframe containing regions to be removed must contain the columns CHR,BP, and BP2!")
          #setDT(filter.region.df)
          setkeyv(filter.region.df, cols = c("CHR","BP","BP2"))
          y1.nSNP<-nrow(y1)
          setkeyv(y1, cols = c("CHR","BP"))
          for(isegment in 1:nrow(filter.region.df)){
            #isegment<-1
            y1 <- y1[!is.na(CHR) & !is.na(BP) & !(CHR==filter.region.df$CHR[isegment] & BP>=filter.region.df$BP[isegment] & BP<=filter.region.df$BP2[isegment]),]
          }
          LOG("Removing ",as.character(y1.nSNP-nrow(y1))," custom region variants, ", nrow(y1), " remain")
        } else {
          LOG("Warning: No chromosome or base-pair position information available - no filtering of custom provided regions was done!")
        }
      }
      
      ##mod addition
      ## REMOVE SNPs INFO >1, <0 etc.
      if("INFO" %in% names(y1)){
        #>1
        nrm<-nrow(y1[INFO>1.0, ])
        if(nrm>0){
          y1[INFO>1.0, INFO:=1.0]
          LOG("WARNING: Setting ", nrm, " SNPs with INFO >1 to 1")
        }
        #<0
        nrm<-nrow(y1[INFO<0, ])
        if(nrm>0){
          y1[INFO<0, INFO:=0]
          LOG("WARNING: Setting ", nrm, " SNPs with INFO <0 to 0")
        }
        LOG(paste0("INFO deciles:",
                   paste(round(quantile(y1$INFO,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), na.rm=T),3), collapse =" ")
                   ))
      }
      if("SINFO" %in% names(y1)){
        #>1
        nrm<-nrow(y1[SINFO>1.0, ])
        if(nrm>0){
          y1[SINFO>1.0, SINFO:=1.0]
          LOG("WARNING: Setting ", nrm, " SNPs with SINFO >1 to 1")
        }
        #<0
        nrm<-nrow(y1[SINFO<0, ])
        if(nrm>0){
          y1[SINFO<0, SINFO:=0]
          LOG("WARNING: Setting ", nrm, " SNPs with SINFO <0 to 0")
        }
        LOG(paste0("SINFO deciles:",
                   paste(round(quantile(y1$SINFO,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), na.rm=T),3), collapse =" ")
        ))
      }
      
      
      ##mod addition
      ## REMOVE SNPs MAF<filter.maf and INFO<filter.info
      if(!is.na(filter.maf)){
        if("FRQ" %in% names(y1)){
          rm <- (!is.na(y1$FRQ) & ((y1$FRQ<filter.maf & y1$FRQ<0.5) | (1-y1$FRQ)<filter.maf))
          y1 <- y1[!rm, ]
          LOG("Removing ", sum(rm), " SNPs with MAF <", filter.maf, "; ", nrow(y1), " remain")
        } else {
          LOG("Warning: The dataset does not contain a FRQ or MAF column to apply the specified filter on.")
        }
      }
      ##mod addition
      if(!is.na(filter.info)){
        if("INFO" %in% names(y1)){
          rm <- (!is.na(y1$INFO) & y1$INFO<filter.info)
          y1 <- y1[!rm, ]
          LOG("Removing ", sum(rm), " SNPs with INFO <", filter.info, "; ", nrow(y1), " remain")
        } else {
          LOG("Warning: The dataset does not contain an INFO column to apply the specified filter on.")
        }
      }

      if(!is.na(filter.info.lt)){
        if("INFO" %in% names(y1)){
         #considers NA or INFO ~1 being genotyped rather than imputed
          rm <- (!is.na(y1$INFO) & (y1$INFO>filter.info.lt & y1$INFO<0.9999))
          y1 <- y1[!rm, ]
          LOG("Removing ", sum(rm), " SNPs (imputed only) with INFO >", filter.info.lt, "; ", nrow(y1), " remain")
         
        } else {
          LOG("Warning: The dataset does not contain an INFO column to apply the specified filter on.")
        }
      }
      
      ##mod addition
      if(!is.na(filter.ldimp.k)){
        if(any(colnames(y1)=="LDIMP.K")){
          rm <- (!is.na(y1$LDIMP.K) & y1$LDIMP.K<filter.ldimp.k)
          y1 <- y1[!rm, ]
          LOG("Removing ", sum(rm), " SNPs with LDIMP.K <", filter.ldimp.k, "; ", nrow(y1), " remain")
        } else {
          LOG("Warning: The dataset does not contain an LDIMP.K column to apply the specified filter on.")
        }
      }
      
      ##mod addition - remove variants in specified chromosomes
      if(length(leave.chr)>0){
        if("CHR" %in% names(y1)){
          for(chrToLeave in leave.chr){
            rm <- (!is.na(y1$CHR) & y1$CHR==chrToLeave)
            y1 <- y1[!rm, ]
            LOG("Removing ", sum(rm), " SNPs in chromosome with number ", chrToLeave, "; ", nrow(y1), " remain")
          }
        } else {
          LOG("Warning: The dataset does not contain a CHR column to apply the specified filter on.")
        }
      }
      
     
      
      ## Merge files
      
      #mod change - merge using data.table
      setkeyv(y1, cols = c("SNP"))
      
      x.columnNames<-stdGwasColumnNames(colnames(x),missingEssentialColumnsStop = NULL,ancestrySetting = ancestrySetting[iTrait], warnings=F) #requires the SHRU package!!
      merged<-x
      colnames(merged)<-x.columnNames$std
      
      xcols <- c("SNP")
      if(any(colnames(merged)=="CHR")) {
        merged[,CHR_REF:=CHR]
        xcols <- c(xcols,"CHR_REF")
      }
      if(any(colnames(merged)=="BP")) {
        merged[,BP_REF:=BP]
        xcols <- c(xcols,"BP_REF")
      }
      if(any(colnames(merged)=="A1")) {
        merged[,A1_REF:=A1]
        xcols <- c(xcols,"A1_REF")
      }
      if(any(colnames(merged)=="A2")) {
        merged[,A2_REF:=A2]
        xcols <- c(xcols,"A2_REF")
      }
      if(any(colnames(merged)=="FRQ")) {
        merged[,FRQ_REF:=FRQ]
        xcols <- c(xcols,"FRQ_REF")
      }
      if(any(colnames(merged)=="FRQFB")) {
        xcols <- c(xcols,"FRQFB")
      }
      if(any(colnames(merged)=="L2")) xcols <- c(xcols,"L2")
      if(any(colnames(merged)=="L2FB")) xcols <- c(xcols,"L2FB")
      if(any(colnames(merged)=="CM")) xcols <- c(xcols,"CM")
      merged<-merged[,..xcols]
      
      #merge with current dataset
      merged<-y1[merged, on=c("SNP"), nomatch=0]
      
      #Force CHR from ref if exists
      if(any(colnames(merged)=="CHR_REF")) merged[,CHR:=CHR_REF]
      #Force BP from ref if exists
      if(any(colnames(merged)=="BP_REF")) merged[,BP:=BP_REF]
      
      #FRQ and fallback FRQ
      if(any(colnames(merged)=="FRQ_REF")){
        if(any(colnames(merged)=="FRQ")){
          if(nrow(merged[is.na(FRQ),])>0){
            LOG("Missing FRQ found, attempting to infer from reference FRQ: ",as.character(nrow(merged[is.na(FRQ),])))
            merged[is.na(FRQ) & is.finite(FRQ_REF),FRQ:=as.numeric(FRQ_REF)]
          }
        } else {
          LOG("Missing FRQ found, attempting to infer from reference FRQ")
            merged[,FRQ:=as.numeric(FRQ_REF)]
        }
      }
      
      if(any(colnames(merged)=="FRQFB")){
        if(any(colnames(merged)=="FRQ")){
          if(nrow(merged[is.na(FRQ),])>0){
            LOG("FRQ still missing, attempting to infer from reference fallback FRQ: ",as.character(nrow(merged[is.na(FRQ),])))
            merged[is.na(FRQ) & is.finite(FRQFB),FRQ:=as.numeric(FRQFB)]
          }
        } else {
          LOG("Missing FRQ found, attempting to infer from reference fallback FRQ")
          merged[,FRQ:=as.numeric(FRQFB)]
        }
      }
      
      #fallback L2
      if(any(colnames(merged)=="L2FB")){
        if(any(colnames(merged)=="L2")){
          if(nrow(merged[is.na(L2) | L2==0,])>0){
            LOG("L2 still missing, attempting to infer from reference fallback L2: ",as.character(nrow(merged[(is.na(L2) | L2==0) & is.finite(L2FB),])))
            merged[(is.na(L2) | L2==0) & is.finite(L2FB),L2:=as.numeric(L2FB)]
          }
        } else {
          LOG("Missing L2 found, attempting to infer from reference fallback L2")
          merged[,L2:=as.numeric(L2FB)]
        }
      }
      
      #mod change - LD reference QC moved here!
      if(any(colnames(merged)=="L2")) merged <- merged[is.finite(L2) & L2!=0,] #& is.finite(L2) & L2!=0
      
      
      #mod addition - extra sync with reference alleles if available columns - from supermunge conditions
      if(any(colnames(merged)=="A1") & any(colnames(merged)=="A2") & any(colnames(merged)=="A1_REF") & any(colnames(merged)=="A2_REF")){
        merged[,cond.invertedAlleleOrder:=((A2!=A2_REF & A1==A2_REF) | (A1!=A1_REF & A2 ==A1_REF))]
        if(any(colnames(merged)=="cond.invertedAlleleOrder")) {
          ## Invert alleles
          merged[,A1:=ifelse(cond.invertedAlleleOrder, A2, A1)]
          merged[,A2:=ifelse(cond.invertedAlleleOrder, A1, A2)]
          ##invert FRQ
          if(any(colnames(merged)=="FRQ")) merged[cond.invertedAlleleOrder==T,FRQ:=1-FRQ]
          ##invert effect
          merged[,Z:=ifelse(cond.invertedAlleleOrder,(Z*-1),Z)]
          
          LOG("ATTENTION: Corrected discrepancies between reference allele configuration and dataset: ",nrow(merged[cond.invertedAlleleOrder==T,]))
        }
        
      }
      
      #add weights - after inferring fallback L2!
      if(!is.null(w)){
        merged<-merged[w[,c("SNP","wLD")], on=c("SNP"), nomatch=0]
      } else {
        merged[,wLD:=L2]
      }
      y1.nrow<-nrow(y1)
      rm(y1)
      
      #old
      # merged <- merge(y1[, y1.columns], w[, c("SNP", "wLD")], by = "SNP", sort = FALSE, incomparables=NA)
      # merged <- merge(merged, x, by = "SNP", sort = FALSE, incomparables=NA)
      # merged <- merged[with(merged, order(CHR, BP)), incomparables=NA]
      
      LOG("Length of unique SNPs: ",length(unique(merged$SNP)), " vs total no. SNPs: ", nrow(merged))
      
      #remove duplicates across the SNP column - same algorithm as in supermunge(?)
      if(length(unique(merged$SNP)) < nrow(merged) & any(colnames(merged)=="FRQ")){
        LOG("Removing residual SNP id duplicates.")
        merged$order<-1:nrow(merged)
        if(any(colnames(merged)=="L2")){
          merged<-merged[order(-FRQ,-L2),]
        } else {
          merged<-merged[order(-FRQ),]
        }
        m.unique<-merged[, .(SNP = head(SNP,1),order = head(order,1)), by = c("SNP")]
        setkeyv(m.unique, cols = c("order"))
        setkeyv(merged, cols = c("order"))
        merged<-merged[m.unique, on=c("order"), nomatch=0][,order:=NULL]
      }
      
      LOG("Out of ", y1.nrow, " SNPs, ", nrow(merged), " remain after merging with LD-score files")
      LOG("Columns after merge: ", paste(colnames(merged),collapse = " ")) ##mod addition
      if(all(is.na(merged$N))) LOG("Warning: The data has no N values!!") ##mod addition
      
      
      #emergency fix to adjust N for imputed variants, no matter what was entered for them previously
      # This is mainly for SSimp as it sets the N value in relation to the imputation quality.
      if(any(colnames(merged)=="SINFO")){
        #highQN <- quantile(merged$N, c(.75))[[1]]
        if(resetImputedN){
          meanN <- mean(merged$N,na.rm=T)
          if(nrow(merged[is.finite(SINFO) & N< 0.95 * eval(meanN),])/nrow(merged[is.finite(SINFO),]) > 0.05){
            merged[!is.na(SINFO),]$N <- meanN
            LOG("Fix: Adjust N to mean N for all imputed variants")
          }
        }
        LOG("Imputed variants only N quartiles:")
        LOG(quantile(merged[!is.na(SINFO),]$N, na.rm = T))
      }
      
      #mod addition - adjust for genomic inflation using either the inflation factor or ldsc intercept
      if(!is.null(GC)){
        #setDT(merged) #to make computations faster
        LOG("\n",trait.names[iTrait],", GC=", GC[iTrait])
        merged[,Z:=Z/GC[iTrait]]
        merged[,EFFECT:=SE*Z]
        merged[,P:=2*pnorm(q = abs(Z),mean = 0, sd = 1, lower.tail = F)]
        #merged<-as.data.frame(merged) #restore merged
      }
      
      
      #re-inflate GWAS summary statistics imputed effects by using the median Chi-square association as a guidemark
      #median version
      # if(any(colnames(merged)=="SINFO")){ #assume that all SINFO variants are actually imputed variants
      #   #assume that all SINFO variants are actually imputed variants by this point
      #   medianChi2.genotyped <- median(merged[is.na(SINFO),]$Z^2,na.rm=T)
      #   LOG("Median Chi^2 (genotyped):",medianChi2.genotyped)
      #   medianChi2.imputed <- median(merged[is.finite(SINFO) & SINFO > 0,]$Z^2,na.rm=T)
      #   LOG("Median Chi^2 (GWAS sumstats imputed):",medianChi2.imputed)
      #   l_reinflation<-medianChi2.imputed/medianChi2.genotyped
      #   merged[is.finite(SINFO) & SINFO > 0,Z:=Z/eval(l_reinflation^reinflateImputedGWAS.inflationFactorExponent)]
      #   merged[is.finite(SINFO) & SINFO > 0,EFFECT:=SE*Z]
      #   merged[is.finite(SINFO) & SINFO > 0,P:=2*pnorm(q = abs(Z),mean = 0, sd = 1, lower.tail = F)]
      #   medianChi2.imputed <- median(merged[is.finite(SINFO) & SINFO > 0,]$Z^2,na.rm=T)
      #   LOG("Median Chi^2 (GWAS sumstats imputed), after re-inflation:",medianChi2.imputed)
      # }
      
      
      ## REMOVE SNPS with excess chi-square:
      
      if(is.na(filter.chisq.max)){
        filter.chisq.max <- max(0.001 * max(merged$N), 80)
      }
      if(!is.na(filter.chisq.max)){
        rm <- (merged$Z^2 > filter.chisq.max)
        merged <- merged[!rm, ]
        
        LOG("Removing ", sum(rm), " SNPs with Chi^2 > ", filter.chisq.max, "; ", nrow(merged), " remain")
      }
      
      ##mod addition - remove SNPs with insufficient chi-square
      if(!is.na(filter.chisq.min)){
        rm <- (merged$Z^2 < filter.chisq.min)
        merged <- merged[!rm, ]
        LOG("Removing ", sum(rm), " SNPs with Chi^2 < ", filter.chisq.min, "; ", nrow(merged), " remain")
      }
      
      #re-inflate GWAS summary statistics imputed effects by using the median Chi-square association as a guidemark
      #mean version - benefits from coming after the removal of extreme associations
      if(any(colnames(merged)=="SINFO")){ #assume that all SINFO variants are actually imputed variants
        #assume that all SINFO variants are actually imputed variants by this point
        #Uses the meanChiSquare calculated earlier - before any genotyoped variants have been filtered out
        LOG("Mean Chi^2 (genotyped):",meanChi2.genotyped)
        meanChi2.imputed <- mean(merged[is.finite(SINFO) & SINFO > 0,]$Z^2,na.rm=T)
        LOG("Mean Chi^2 (GWAS sumstats imputed):",meanChi2.imputed)
        l_reinflation<-meanChi2.imputed/meanChi2.genotyped
        merged[is.finite(SINFO) & SINFO > 0,Z:=Z/eval(l_reinflation^reinflateImputedGWAS.inflationFactorExponent)]
        merged[is.finite(SINFO) & SINFO > 0,EFFECT:=SE*Z]
        merged[is.finite(SINFO) & SINFO > 0,P:=2*pnorm(q = abs(Z),mean = 0, sd = 1, lower.tail = F)]
        meanChi2.imputed <- median(merged[is.finite(SINFO) & SINFO > 0,]$Z^2,na.rm=T)
        LOG("Mean  Chi^2 (GWAS sumstats imputed), after re-inflation:",meanChi2.imputed)
        if(!is.na(filter.chisq.max)){
          rm <- (merged$Z^2 > filter.chisq.max)
          merged <- merged[!rm, ]
          LOG("(again after re-inflation) Removing ", sum(rm), " SNPs with Chi^2 > ", filter.chisq.max, "; ", nrow(merged), " remain")
        }
      }
      
      #mod addition - do strict selection of columns before the NA filter below
      #colToKeep<-c("SNP","CHR","BP","A1","A2","N","Z","L2","wLD","P")
      colToKeep<-c("SNP","A1","A2","Z","L2","wLD")
      if(any(colnames(merged)=="CHR")) colToKeep<- c(colToKeep,"CHR")
      if(any(colnames(merged)=="BP")) colToKeep<- c(colToKeep,"BP")
      if(any(colnames(merged)=="N")) colToKeep<- c(colToKeep,"N")
      if(any(colnames(merged)=="P")) colToKeep<- c(colToKeep,"P")
      if(any(colnames(merged)=="CM")) colToKeep<- c(colToKeep,"CM")
      if(any(colnames(merged)=="INFO")) colToKeep<- c(colToKeep,"INFO")
      if(any(colnames(merged)=="SINFO")) colToKeep<- c(colToKeep,"SINFO")
      merged <- merged[,..colToKeep]
      
      ##mod addition - remove any rows with na values as was originally performed earlier before merge
      if(any(is.na(merged))){
        nBefore <- nrow(merged)
        
        notNecessaryColumns <- c("INFO","SINFO")
        noNAColumns <- colToKeep[!(colToKeep %in% notNecessaryColumns)]
        merged <- na.omit(merged,cols = noNAColumns)
        LOG("Removing ", nBefore-nrow(merged), " SNPs with NA-values present; ", nrow(merged), " remain")
      }
      
      #mod addition - set index key - order will be finally set later before block construction, after merge with second dataset
      setkeyv(merged,cols = c("SNP"))
      
      impstats[iTrait,c("m.merged","meanChiSquare.merged")] <- c(nrow(merged), mean((merged$Z^2),na.rm=T))
      if(any(colnames(merged)=="SINFO")) {
        if(any(colnames(merged)=="BETA.I") & any(colnames(merged)=="SE.I")) {
          merged[,d.z:=(BETA.I/SE.I)-(BETA/SE)]
        } else {
          merged[,d.z:=NA_real_]
        }
        impstats[iTrait,c("m.merged.imputed","meanChiSquare.merged.imputed","meanChiSquare.merged.nonimputed","rmse_z")] <- c(nrow(merged[is.finite(SINFO),]),mean((merged[is.finite(SINFO) & is.finite(Z),]$Z^2),na.rm=T),mean((merged[is.na(SINFO) & is.finite(Z),]$Z^2),na.rm=T), sqrt(mean(merged[is.finite(d.z),]$d.z^2, na.rm=T)))
      } else {
        impstats[iTrait,c("m.merged.imputed","meanChiSquare.merged.imputed","meanChiSquare.merged.nonimputed","rmse_z")] <- c(NA_real_,NA_real_,NA_real_,NA_real_)
      }
      all_y[iTrait]<-list(merged)
      rm(merged)
      iTrait<-iTrait+1
    }
      
    #mod addition - record the median number of variants
    n_all_y <- as.numeric(lapply(all_y,function(x)(nrow(x))))
    LOG("N per dataset")
    median_all_y <- median(n_all_y,na.rm = T)
    LOG("The total number of variants (non unique) in the analysis is ",sum(n_all_y,na.rm = T))
    
    #mod addition - record the number of overlapping variants
    snplist<-x$SNP
    for(j in 1:n.traits){
      #j<-1
      all_y[[j]]$SNP
      snplist<-snplist[snplist %in% all_y[[j]]$SNP]
    }
    mgc <- length(unique(snplist))
    LOG("The overall genetic variant intersect across datasets is ",round(mgc/median_all_y,digits = 3)," of the median dataset (",mgc," vs. ",median_all_y, ")")
    if(mgc/median_all_y<0.5) LOG("Warning: The overall genetic variant intersect across datasets is less than 50% of the median dataset size!")
    
    #some memory cleanup
    rm(x)
    if(!is.null(w)) rm(w)
    #mod addition
    gc() #do garbage collect if this can help with out of memory issues.
  
  } # read files block
  
  impstats$trait.name <- trait.names
  impstats$file <- traits
  rownames(impstats) <- trait.names
  
  if(exportReadyMergedTraits){
    LOG("Returning merged traits and ending program!")
    return(list(all_y=all_y,impstats=impstats,median_all_y=median_all_y, mgc=mgc, M.tot=M.tot))
  }
  
  if(resamplingMethod=="jn"){
    LOG("Using jackknife resampling")
  } else { #if(resamplingMethod="vbcs")
    LOG("Using variable block-count sampling")
  }
  
  if(preweight.alternativeCorrelationCorrection){
    LOG("Using alternate correlation correction")
  }
  
  LOG("M = ", M.tot)
  
  #fold loop
  n.folds<-2*k.folds
  foldResults<-list()
  if(k.folds>0){
    for(iKfold in 1:k.folds){
      foldResults[[iKfold]]<-list(train=c(),test=c())
    }
  }
  
  for(iFold in 0:n.folds){
    #iFold<-0
      
    # count the total nummer of runs, both loops
    s <- 1
    #s <- 2
    
    for(j in 1:n.traits){
      #j<-1
      chi1 <- traits[j]
      
      y1 <- all_y[[j]]
      y1$chi1 <- y1$Z^2
      
      #moved here as they are the same for j
      samp.prev <- sample.prev[j]
      pop.prev <- population.prev[j]
      
     
      
      for(k in j:length(traits)){
        #k<-2
        
        ##### GENETIC COVARIANCE code
        
        #if(j != k){
        
        LOG("     ", print = FALSE)
        
        chi2 <- traits[k]
        LOG("Calculating genetic covariance [", s, "/", n.V, "] for traits: ", trait.names[j], " and ", trait.names[k])
        
        #mod change - data.table merge
        merged <- y1[all_y[[k]], on=c("SNP"), nomatch=0]
        n.snps <- nrow(merged)
        
        #filter on fold
        if(k.folds>0){
          foldW<-floor(n.snps/k.folds)
          trainW<-floor(trainRatio*n.snps)
          iKfold <- floor((iFold+1)/2)
          
          #training set limits
          foldILow<-(iKfold-1)*foldW + 1
          foldIHigh<-foldILow + trainW - 1
          
          if(foldIHigh<=n.snps){
            foldRowSelectionTrain<-foldILow:foldIHigh
          } else {
            foldIHigh<-foldIHigh-n.snps
            foldRowSelectionTrain<-c(1:foldIHigh,foldILow:n.snps)
          }
          
          if(iFold %% 2 == 0){
            #test
            merged <- merged[!foldRowSelectionTrain,]
            LOG("Fold ",iKfold, " (test) with train selection as [",foldILow,",",foldIHigh,"]/",n.snps)
          } else {
            #train
            merged <- merged[foldRowSelectionTrain,]
            LOG("Fold ",iKfold, " (train) with train selection as [",foldILow,",",foldIHigh,"]/",n.snps)
          }
          
          n.snps <- nrow(merged)
        }
        
        LOG("Trait pair has an overlap of ",n.snps," variants, ",round(n.snps/median_all_y,digits = 3)," of the median merged dataset.")
        
        #filter based on SINFO score and completeness of data
        if(any(colnames(merged)=="SINFO") ){ #& any(colnames(merged)=="i.SINFO") #in case only one of the datasets has SINFO
          if(!any(colnames(merged)=="i.SINFO")) merged[,i.SINFO:=NA_real_]
          if(!is.na(filter.sinfo) | !is.na(filter.sinfo.imputed)){
            merged <- merged[(is.na(SINFO) & is.na(i.SINFO)) 
                        | (is.na(SINFO) & (is.na(filter.sinfo) | i.SINFO>eval(filter.sinfo)))
                        | (is.na(i.SINFO) & (is.na(filter.sinfo) | SINFO>eval(filter.sinfo)))
                        | ((SINFO>eval(filter.sinfo.imputed) & i.SINFO>eval(filter.sinfo.imputed)) | (is.na(filter.sinfo.imputed) & (!is.na(SINFO) & !is.na(i.SINFO)))),]
          }
          
          m.nonimputed <- nrow(merged[is.na(SINFO) & is.na(i.SINFO),])
          
          if(!is.na(filter.simpute.fraction)){
            merged.imp<-head(merged[is.finite(SINFO) | is.finite(i.SINFO),][order(-(is.na(SINFO) | is.na(i.SINFO)),-(SINFO+i.SINFO),-(L2+i.L2)),c("SNP","SINFO","i.SINFO")],filter.simpute.fraction*m.nonimputed) #use partial imputations first
            merged[merged.imp,on=c("SNP"), keep:=T]
            rm(merged.imp)
            merged <- merged[(is.na(SINFO) & is.na(i.SINFO)) | keep,]
          }
          
          LOG(n.snps - nrow(merged), " variants excluded based on SINFO")
          
        }
        
        #LOG("Colnames of merged y: ",colnames(y))
        #old
        #y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")], by = "SNP", sort = FALSE, incomparables=NA)
        
        #mod change - improved flipped effect direction condition - this should not be much of an issue if reference harmonisation is working properly from munge or from the new harmonisation step here.
        merged[,Z:=ifelse((A2!=i.A2 & A1==i.A2) | (A1!=i.A1 & A2 ==i.A1), -Z, Z)]
        #y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x) #orig
        #y <- as.tibble(y) #to conform with any non-data.table syntax below #testing to remove this for the below to work with data table
        
        merged$ZZ <- merged$i.Z * merged$Z
        merged$chi2 <- merged$i.Z^2
        #merged$ZZ.var<-2
        #merged <- na.omit(y) #removed this to avoid deleting rows with missing imputation info - watch out for errors later!!
        if(is.finite(filter.zz.min)){
          n.snps <- nrow(merged)
          #NEW! ZZ filter - filter on ZZ rather than individual trait Chi-square stats
          merged <- merged[abs(ZZ)>filter.zz.min, ]
          nzzexcluded <- n.snps - nrow(merged)
          if(nzzexcluded>0) LOG(nzzexcluded, " variants excluded based on ZZ")
        }
        
        
        #cellstats - after final filter
        
        if(any(colnames(merged)=="SINFO") ){
          mCellstats.imputed.m[k,j]<-mCellstats.imputed.m[j,k]<-nrow(merged[is.finite(SINFO) | is.finite(i.SINFO),])
          mCellstats.imputed.full.m[k,j]<-mCellstats.imputed.full.m[j,k]<-nrow(merged[is.finite(SINFO) & is.finite(i.SINFO),])
          mCellstats.imputed.partial.m[k,j]<-mCellstats.imputed.partial.m[j,k]<-nrow(merged[(is.finite(SINFO) & is.na(i.SINFO)) | (is.na(SINFO) & is.finite(i.SINFO)),])
          mCellstats.nonimputed.m[k,j]<-mCellstats.nonimputed.m[j,k]<-m.nonimputed
        } else {
          mCellstats.imputed.m[k,j]<-mCellstats.imputed.m[j,k]<-0
          mCellstats.imputed.full.m[k,j]<-mCellstats.imputed.full.m[j,k]<-0
          mCellstats.imputed.partial.m[k,j]<-mCellstats.imputed.partial.m[j,k]<-0
          mCellstats.nonimputed.m[k,j]<-mCellstats.nonimputed.m[j,k]<-nrow(merged)
        }
        
        
        n.snps <- nrow(merged)
        
        LOG(n.snps, " SNPs remain after merging (and merge filters: SINFO, ZZ)")
        
        #mod addition - compute M from the actual number of variants in the merged dataset
        #M_5_50 total sum was 7925897 for the new reference panel LD scores
        
        
        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        
        
        #### MAKE WEIGHTS:
        
        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        if(preweight.alternativeCorrelationCorrection){
          #LOG("Using alternate correlation correction")
          merged$ld <- pmax(merged$L2, 0)
          merged$w.ld <- pmax(merged$wLD, 0)
          merged$oc.w <- (1/(1+merged$w.ld)) #mod change - added 1+ to the denominator as this corresponds more to the original LDSC how it is described in the publication (i.e. the l(S) weight is actually 1 + the ld-score), and makes more sense to change the weight from 1 to less in proportion to the ld-score, rather than also upweighting the variant in case of LD-score<1. HOWEVER!, the ld scores used are restricted to a minimum of 1, so not clear if needed anymore, or if the pmax limit should be set to 0 in combination with this weighting correction.
        } else {
          merged$ld <- pmax(merged$L2, 1)
          merged$w.ld <- pmax(merged$wLD, 1)
          merged$oc.w <- 1/merged$w.ld
        }
        merged$c <- tot.agg*merged$N/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        
        merged[,chiSquare.w:=1]
        if(preweight.ChiSquare) merged[,chiSquare.w:=abs(Z)^2]
        merged[,gimp.w:=1]
        if(any(colnames(merged)=="INFO") & preweight.INFO) merged[is.finite(INFO),gimp.w:=INFO]
        merged[,imp.w:=1]
        if(any(colnames(merged)=="SINFO") & preweight.SINFO) merged[is.finite(SINFO),imp.w:=SINFO]
        merged$w <- sqrt(merged$het.w*merged$oc.w*merged$chiSquare.w*merged$gimp.w*merged$imp.w)
        # merged$w <- merged$het.w*merged$oc.w
        # merged$initial.w <- sqrt(merged$w)
        
        tot.agg2 <- (M.tot*(mean(merged$chi2)-1))/mean(merged$L2*merged$i.N)
        tot.agg2 <- max(tot.agg2,0)
        tot.agg2 <- min(tot.agg2,1)
        if(preweight.alternativeCorrelationCorrection){
          #LOG("Using alternate correlation correction")
          merged$ld2 <- pmax(merged$i.L2, 0)
          merged$w.ld2 <- pmax(merged$wLD, 0)
          merged$oc.w2 <- (1/(1+merged$w.ld2)) #mod change - see above
        } else {
          merged$ld2 <- pmax(merged$i.L2, 1)
          merged$w.ld2 <- pmax(merged$wLD, 1)
          merged$oc.w2 <- 1/merged$w.ld2
        }
        merged$c2 <- tot.agg2*merged$i.N/M.tot
        merged$het.w2 <- 1/(2*(1+(merged$c2*merged$ld2))^2)
        
        merged[,chiSquare.w2:=1]
        if(preweight.ChiSquare) merged[,chiSquare.w2:=abs(i.Z)^2]
        merged[,gimp.w2:=1]
        if(any(colnames(merged)=="i.INFO") & preweight.INFO) merged[is.finite(i.INFO),gimp.w2:=i.INFO]
        merged[,imp.w2:=1]
        if(any(colnames(merged)=="i.SINFO") & preweight.SINFO) merged[is.finite(i.SINFO),imp.w2:=i.SINFO]
        merged$w2 <- sqrt(merged$het.w2*merged$oc.w2*merged$chiSquare.w2*merged$gimp.w2*merged$imp.w2)
        # merged$w2 <- merged$het.w2*merged$oc.w2
        # merged$initial.w2 <- sqrt(merged$w2)
        
        #merged$weights_cov <- (merged$w + merged$w2)/2 #TEST!!! keep the original scale of the weights
        merged$weights_cov <- (merged$w + merged$w2)/sum(merged$w + merged$w2, na.rm = T) #mod addition - remove NAs here
        #merged$weights_cov <- (merged$initial.w + merged$initial.w2)/sum(merged$initial.w + merged$initial.w2 ) #original
        
        merged[,L2:=(L2+i.L2)/2][,i.L2:=NULL]
        
        #https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
        N.bar <- sqrt(mean(merged$N)*mean(merged$i.N))
        merged[,N.bar:=exp((log(N)+log(i.N))/2)] #mod addition, for the per-variant N computations
        
        ## preweight LD and chi:
        merged[,c("weighted.LD","weighted.intercept","weighted.chi"):=list(weights_cov*L2,weights_cov*intercept,weights_cov*ZZ)]
        
        if(preweight){
          merged[,c("m.LD","m.intercept","m.Chi"):=list(weighted.LD,weighted.intercept,weighted.chi)]
        } else {
          merged[,c("m.LD","m.intercept","m.Chi"):=list(L2,intercept,ZZ)]
        }
        
        ## Perfrom analysis:
        n.annot <- 1
        
        #set strict unified ordering before block construction
        if(any(colnames(merged)=="CHR") & any(colnames(merged)=="BP")){
          setorder(merged,CHR,BP,SNP)
        } else if(any(colnames(merged)=="CHR")){
          setorder(merged,CHR,SNP)
        } else {
          setorder(merged,SNP)
        }
        
        #block construction follows
        if(!is.finite(n.blocks) & (is.finite(blocksizeCM) | is.finite(minBlocksizeCM))){
          #determine size and loop over each block immediately
          
          if(!any(colnames(merged)=="CM")){
            stop("There is no CM column available to use for cM aware blocks. Please use an LD-score reference containing the CM column, or turn off the cM blocksize by setting blocksizeCM=NA and minBlocksizeCM=NA.")
          }
          
          #setkeyv(merged, cols = c("CHR","BP","CM")) #old sort setting
          chrs <- unique(merged$CHR)
          xty.block.values <- c()
          xtx.block.values <- c()
          xty.block.values.pos <- xty.block.values
          xtx.block.values.pos <- xtx.block.values
          xty.block.values.neg <- xty.block.values
          xtx.block.values.neg <- xtx.block.values
          LDvsChiCorrelation.block.values<-c()
          LDvsChiCorrelation.block.values.pos<-c()
          LDvsChiCorrelation.block.values.neg<-c()
          LD.block.SS<-c()
          nVars.block<-c()
          
          mPos<-c()
          mNeg<-c()
          
          iBlock <- 1
          for(cChr in chrs){
            #cChr <- 1
            merged.chr<-merged[CHR==eval(cChr),][,m:=.I]
            setkeyv(merged.chr, cols = c("m","CM"))
            m.max<-nrow(merged.chr)
            CM.max<-max(merged.chr$CM)
            m.from <- 0
            CM.from <- min(merged.chr$CM) - 0.01
            
            while(T){
              
              if(is.finite(minBlocksizeCM) & is.finite(blocksizeM)) {
                merged.block <- merged.chr[m>eval(m.from) & CM>eval(CM.from) & (m <= (eval(m.from) + blocksizeM) | CM <= (eval(CM.from) + minBlocksizeCM)),]
              } else { #is.finite(blocksizeCM)
                merged.block <- merged.chr[m>eval(m.from) & CM>eval(CM.from) & CM <= (eval(CM.from) + blocksizeCM),]
              }
              
              
              if(nrow(merged.block)>=minBlocksizeM){
                
                #original - using mean N later
                mLD<-as.matrix(cbind(merged.block$m.LD,merged.block$m.intercept))
                mChi<-as.matrix(cbind(merged.block$m.Chi))
                
                if(doubleRegressionRoutine){
                  mPos[iBlock]<-sum(mChi>0)
                  mNeg[iBlock]<-sum(mChi<0)
                  #if(sum(mChi>0)>1){
                  xty.block.values.pos[iBlock] <- list(as.data.frame(t(t(mLD[which(mChi>0),]) %*% mChi[which(mChi>0),])))
                  xtx.block.values.pos[iBlock] <- list(as.data.frame(t(t(mLD[which(mChi>0),]) %*% mLD[which(mChi>0),])))
                  LDvsChiCorrelation.block.values.pos[iBlock]<-cor(x = mLD[which(mChi>0),1],y = mChi[which(mChi>0)])
                  #}
                  if(mNeg[iBlock]>1){
                    xty.block.values.neg[iBlock] <- list(as.data.frame(t(t(mLD[which(mChi<0),]) %*% mChi[which(mChi<0),])))
                    xtx.block.values.neg[iBlock] <- list(as.data.frame(t(t(mLD[which(mChi<0),]) %*% mLD[which(mChi<0),])))
                    LDvsChiCorrelation.block.values.neg[iBlock]<-cor(x = mLD[which(mChi<0),1],y = mChi[which(mChi<0)])
                  }
                } 
                #else { #run this always for testing
                  xty.block.values[iBlock] <- list(as.data.frame(t(t(mLD) %*% mChi)))
                  xtx.block.values[iBlock] <- list(as.data.frame(t(t(mLD) %*% mLD)))
                  LDvsChiCorrelation.block.values[iBlock]<-cor(x = mLD,y = mChi)
                #}
                LD.block.SS[iBlock]<-sum((merged.block$m.LD - mean(merged.block$m.LD, na.rm=T))^2, na.rm = T)
                nVars.block[iBlock]<-nrow(merged.block)
              
                iBlock <- iBlock + 1
              }
              
              if(is.na(blocksizeCM)){
                m.from <- max(c(max(merged.block$m) ,m.from + blocksizeM))
              } else {
                m.from <- max(max(merged.block$m), m.from + 1)
              }
              
              if(is.na(blocksizeCM)){
                CM.from <- quantile(merged.block$CM, c(.99))[[1]] #to avoid extreme outlier cM values
              } else {
                CM.from <- CM.from + blocksizeCM
              }
              
              if(m.from>=m.max | CM.from >CM.max) break
              
            }
          }
          
          xty.block.values <- as.matrix(rbindlist(xty.block.values))
          xtx.block.values <- as.matrix(rbindlist(xtx.block.values))
          xty.block.values.pos <- as.matrix(rbindlist(xty.block.values.pos))
          xtx.block.values.pos <- as.matrix(rbindlist(xtx.block.values.pos))
          xty.block.values.neg <- as.matrix(rbindlist(xty.block.values.neg))
          xtx.block.values.neg <- as.matrix(rbindlist(xtx.block.values.neg))
          
          n.blocksToUse <- nrow(xty.block.values)
          
          #fallback for later checks
          if(nrow(xty.block.values.pos)<1) xty.block.values.pos <- matrix(data=NA,nrow=n.blocksToUse,ncol =(n.annot+1))
          if(nrow(xtx.block.values.pos)<1) xtx.block.values.pos <- matrix(data=NA,nrow =((n.annot+1)* n.blocksToUse),ncol =(n.annot+1))
          if(nrow(xty.block.values.neg)<1) xty.block.values.neg <- matrix(data=NA,nrow=n.blocksToUse,ncol =(n.annot+1))
          if(nrow(xtx.block.values.neg)<1) xtx.block.values.neg <- matrix(data=NA,nrow =((n.annot+1)* n.blocksToUse),ncol =(n.annot+1))
          
          LOG(n.blocksToUse, " blocks to be used")
          
          
        } else {
          
          #original implementation with initial set blocksize 
          if(!is.finite(n.blocks) & is.finite(blocksizeM)){
            n.blocksToUse <- floor(nrow(merged)/blocksizeM)
          } else {
            n.blocksToUse <- n.blocks
          }
          
          LOG(n.blocksToUse, " blocks to be used")
    
          xty.block.values <- matrix(data=NA,nrow=n.blocksToUse,ncol =(n.annot+1))
          xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocksToUse),ncol =(n.annot+1))
          xty.block.values.pos <- xty.block.values
          xtx.block.values.pos <- xtx.block.values
          xty.block.values.neg <- xty.block.values
          xtx.block.values.neg <- xtx.block.values
          LDvsChiCorrelation.block.values<-c()
          LDvsChiCorrelation.block.values.pos<-c()
          LDvsChiCorrelation.block.values.neg<-c()
          LD.block.SS<-c()
          nVars.block<-c()
          
          mPos<-c()
          mNeg<-c()
          
          #original - using mean N later
          mLD<-as.matrix(cbind(merged$m.LD,merged$m.intercept))
          mChi<-as.matrix(cbind(merged$m.Chi))
          
          select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocksToUse+1)))
          select.to <- c(select.from[2:n.blocksToUse]-1,n.snps)
          replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
          replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
          
          for(i in 1:n.blocksToUse){
            #i<-1
            cMChi <- mChi[select.from[i]:select.to[i],]
            cMLD <- mLD[select.from[i]:select.to[i],]
            
            if(doubleRegressionRoutine){
              mPos[i]<-sum(cMChi>0)
              mNeg[i]<-sum(cMChi<0)
              #if(sum(cMChi>0)>1){
                xty.block.values.pos[i,] <- t(t(cMLD[which(cMChi>0),])%*% cMChi[which(cMChi>0)])
                xtx.block.values.pos[replace.from[i]:replace.to[i],] <- as.matrix(t(cMLD[which(cMChi>0),])%*% cMLD[which(cMChi>0),])
                LDvsChiCorrelation.block.values.pos[i]<-cor(x = cMLD[which(cMChi>0),1],y = cMChi[which(cMChi>0)])
              #}
              if(mNeg[i]>1){
                xty.block.values.neg[i,] <- t(t(cMLD[which(cMChi<0),])%*% cMChi[which(cMChi<0)])
                xtx.block.values.neg[replace.from[i]:replace.to[i],] <- as.matrix(t(cMLD[which(cMChi<0),])%*% cMLD[which(cMChi<0),])
                LDvsChiCorrelation.block.values.neg[i]<-cor(x = cMLD[which(cMChi<0),1],y = cMChi[which(cMChi<0)])
              }
            }
            #else { #run this always for testing
              xty.block.values[i,] <- t(t(cMLD)%*% cMChi)
              xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(cMLD)%*% cMLD)
              LDvsChiCorrelation.block.values[i]<-cor(x = cMLD[,1],y = cMChi)
           # }
              
            LD.block.SS[i]<-sum((cMLD[,1] - mean(cMLD[,1], na.rm=T))^2, na.rm = T) #we don't really use this - YET!
            nVars.block[i]<-nrow(cMLD)
            
          }
          
        } #end of block values construction
        
        hasPositiveResultsCheck<-hasNegativeResultsCheck<-F #init
        mPosTotAvg<-mNegTotAvg<-0
        if(doubleRegressionRoutine){
          
          mPosTotAvg<-mean(mPos,na.rm = T)
          mNegTotAvg<-mean(mNeg,na.rm = T)
          
          xty.pos <- as.matrix(colSums(xty.block.values.pos,na.rm = T))
          xtx.pos <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
          for(i in 1:nrow(xtx.pos)){xtx.pos[i,] <- t(colSums(xtx.block.values.pos[seq(from=i,to=nrow(xtx.block.values.pos),by=ncol(mLD)),],na.rm = T))}
          delete.values.pos <- matrix(data=NA,nrow=nrow(xty.block.values.pos),ncol =(n.annot+1))
          hasPositiveResultsCheck <- (!all(is.na(xtx.pos)) && !all(xtx.pos==0) && !all(is.na(xty.pos)) && !all(xty.pos==0))
          
          xty.neg <- as.matrix(colSums(xty.block.values.neg,na.rm = T))
          xtx.neg <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
          for(i in 1:nrow(xtx.neg)){xtx.neg[i,] <- t(colSums(xtx.block.values.neg[seq(from=i,to=nrow(xtx.block.values.neg),by=ncol(mLD)),],na.rm = T))}
          delete.values.neg <- matrix(data=NA,nrow=nrow(xty.block.values.neg),ncol =(n.annot+1))
          hasNegativeResultsCheck <- (!all(is.na(xtx.neg)) && !all(xtx.neg==0) && !all(is.na(xty.neg)) && !all(xty.neg==0))
        }
        #else {
          xty <- as.matrix(colSums(xty.block.values))
          xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
          for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(mLD)),]))}
          delete.values <- matrix(data=NA,nrow=nrow(xty.block.values),ncol =(n.annot+1))
          #delete.values.unsigned <- matrix(data=NA,nrow=nrow(xty.block.values),ncol =(n.annot+1)) #old implementation
        #}
          
        #this should be made compatible with dr-only runs
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        
        if(resamplingMethod=="jn"){
          
          if(doubleRegressionRoutine){
            for(i in 1:n.blocksToUse){
              #i<-1
              xty.delete.pos <- xty.pos-xty.block.values.pos[i,]
              xtx.delete.pos <- xtx.pos-xtx.block.values.pos[delete.from[i]:delete.to[i],]
              delete.values.pos[i,] <- solve(xtx.delete.pos, xty.delete.pos, tol=1e-30)
              
              xty.delete.neg <- xty.neg-xty.block.values.neg[i,]
              xtx.delete.neg <- xtx.neg-xtx.block.values.neg[delete.from[i]:delete.to[i],]
              delete.values.neg[i,] <- solve(xtx.delete.neg, xty.delete.neg, tol=1e-30)
            }
          }
          #else {
          for(i in 1:n.blocksToUse){
            #i<-1
            xty.delete <- xty-xty.block.values[i,]
            xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
            delete.values[i,] <- solve(xtx.delete, xty.delete, tol=1e-30)
          }
          #}
          
          if(doubleRegressionRoutine){
            reg.pos<-reg.neg<-matrix(NA,1,nrow(xty.pos))
            if(hasPositiveResultsCheck) reg.pos <- solve(xtx.pos, xty.pos)
            if(hasNegativeResultsCheck) reg.neg <- solve(xtx.neg, xty.neg)
          }
          reg <- solve(xtx, xty)
          
          #fallback for non-dr runs
          reg.pos<-reg.neg<-matrix(NA,1,nrow(xty)) 
          xty.delete.pos<-xty.delete.neg<-matrix(NA,nrow(xty.delete),ncol(xty.delete))
          xtx.delete.pos<-xtx.delete.neg<-matrix(NA,nrow(xtx.delete),ncol(xtx.delete))
          delete.values.pos<-delete.values.neg<-matrix(NA,nrow(delete.values),ncol(delete.values))
          
          #convert coeficients to heritability
          intercept <- reg[2]
          intercept.pos <- reg.pos[2]
          intercept.neg <- reg.neg[2]
          
          #coefs <- reg[1]/N.bar
          reg.tot <- (reg[1]/N.bar)*M.tot #replaced m with M.tot, they should be the same
          reg.tot.pos <- (reg.pos[1]/N.bar)*M.tot
          reg.tot.neg <- (reg.neg[1]/N.bar)*M.tot
          
          
          tot.delete.values.pos <- delete.values.pos[,1:n.annot]
          pseudo.values.pos <- matrix(data=NA,nrow=n.blocksToUse,ncol=length(reg.pos))
          for(i in 1:n.blocksToUse){pseudo.values.pos[i,] <- (n.blocksToUse*reg.pos)-((n.blocksToUse-1)* delete.values.pos[i,])}
          tot.delete.values.neg <- delete.values.neg[,1:n.annot]
          pseudo.values.neg <- matrix(data=NA,nrow=n.blocksToUse,ncol=length(reg.neg))
          for(i in 1:n.blocksToUse){pseudo.values.neg[i,] <- (n.blocksToUse*reg.neg)-((n.blocksToUse-1)* delete.values.neg[i,])}
          tot.delete.values <- delete.values[,1:n.annot]
          pseudo.values <- matrix(data=NA,nrow=n.blocksToUse,ncol=length(reg))
          #colnames(pseudo.values)<- colnames(weighted.LD)
          for(i in 1:n.blocksToUse){pseudo.values[i,] <- (n.blocksToUse*reg)-((n.blocksToUse-1)* delete.values[i,])}
          
          #pos
          jackknife.cov.pos <- cov(pseudo.values.pos)/n.blocksToUse #should be n.blocksToUse^2 ??? we do not correct this
          jackknife.se.pos <- sqrt(diag(jackknife.cov.pos))
          intercept.se.pos <- jackknife.se.pos[length(jackknife.se.pos)]
          coef.cov.pos <- jackknife.cov.pos[1:n.annot,1:n.annot]/(N.bar^2)
          cat.cov.pos <- coef.cov.pos*(M.tot %*% t(M.tot))
          tot.cov.pos <- sum(cat.cov.pos)
          tot.se.pos <- sqrt(tot.cov.pos)
          lV.pos[s]<-list(pseudo.values.pos[,1])
          
          #neg
          jackknife.cov.neg <- cov(pseudo.values.neg)/n.blocksToUse #should be n.blocksToUse^2 ??? we do not correct this
          jackknife.se.neg <- sqrt(diag(jackknife.cov.neg))
          intercept.se.neg <- jackknife.se.neg[length(jackknife.se.neg)]
          coef.cov.neg <- jackknife.cov.neg[1:n.annot,1:n.annot]/(N.bar^2)
          cat.cov.neg <- coef.cov.neg*(M.tot %*% t(M.tot))
          tot.cov.neg <- sum(cat.cov.neg)
          tot.se.neg <- sqrt(tot.cov.neg)
          lV.neg[s]<-ifelse(
            j == k,
            list(pseudo.values.pos[,1]),
            list(pseudo.values.neg[,1])) #set diagonal elements to equal the results from the positive regressions

          #signed consensus
          intercept.signed <- (mPosTotAvg*replaceNa(intercept.pos)+mNegTotAvg*replaceNa(intercept.neg))/(mPosTotAvg+mNegTotAvg)
          intercept.se.signed <- (mPosTotAvg*replaceNa(intercept.se.pos)+mNegTotAvg*replaceNa(intercept.se.neg))/(mPosTotAvg+mNegTotAvg) #do we want to do statistical inference here rather than weighted average? do we know the correlations between the variables?
          reg.tot.signed <- (mPosTotAvg*replaceNa(reg.tot.pos)+mNegTotAvg*replaceNa(reg.tot.neg))/(mPosTotAvg+mNegTotAvg)
          tot.se.signed <- (mPosTotAvg*replaceNa(tot.se.pos)+mNegTotAvg*replaceNa(tot.se.neg))/(mPosTotAvg+mNegTotAvg)
          if(hasNegativeResultsCheck){
            lV.signed[s]<-list(
              ((mPosTotAvg*pseudo.values.pos+mNegTotAvg*pseudo.values.neg)/(mPosTotAvg+mNegTotAvg))[,1]
              )
          } else {
            lV.signed[s]<-list(pseudo.values.pos[,1]) #we assume there to always be some positive results
          }
          
          #unsigned consensus
          intercept.unsigned <- (mPosTotAvg*abs(replaceNa(intercept.pos))+mNegTotAvg*abs(replaceNa(intercept.neg)))/(mPosTotAvg+mNegTotAvg)
          intercept.se.unsigned <- (mPosTotAvg*replaceNa(intercept.se.pos)+mNegTotAvg*replaceNa(intercept.se.neg))/(mPosTotAvg+mNegTotAvg) #should equal the signed se
          reg.tot.unsigned <- (mPosTotAvg*abs(replaceNa(reg.tot.pos))+mNegTotAvg*abs(replaceNa(reg.tot.neg)))/(mPosTotAvg+mNegTotAvg)
          tot.se.unsigned <- (mPosTotAvg*replaceNa(tot.se.pos)+mNegTotAvg*replaceNa(tot.se.neg))/(mPosTotAvg+mNegTotAvg) #should equal the signed se
          if(hasNegativeResultsCheck){
          lV.unsigned[s]<-list(
            ((mPosTotAvg*abs(pseudo.values.pos)+mNegTotAvg*abs(pseudo.values.neg))/(mPosTotAvg+mNegTotAvg))[,1]
          )
          } else {
            lV.unsigned[s]<-list(pseudo.values.pos[,1]) #we assume there to always be some positive results
          }
          
          
          #old
          jackknife.cov <- cov(pseudo.values)/n.blocksToUse #should be n.blocksToUse^2 ??? we do not correct this
          jackknife.se <- sqrt(diag(jackknife.cov))
          intercept.se <- jackknife.se[length(jackknife.se)]
          coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
          cat.cov <- coef.cov*(M.tot %*% t(M.tot))
          tot.cov <- sum(cat.cov)
          tot.se <- sqrt(tot.cov)
          lV[s]<-list(pseudo.values[,1])
          
          
          N.vec[1, s] <- N.bar
        
        } else { #if(resamplingMethod="vbcs")
          #LOG("Using variable block-count sampling")
          attenuationCorrectionFactor <-c()
          for(i in 1:n.blocksToUse){
            #i<-1
            #i<-6
            #i<-22
            #i<-209
            #i<-141
            xty.delete <- xty.block.values[i,]
            xtx.delete <- xtx.block.values[delete.from[i]:delete.to[i],]
            cDelete.values<-try(solve(xtx.delete, xty.delete, tol=1e-30))
            if(inherits(cDelete.values,"try-error")){
              warning(paste0("Error when solving general block ",i))
            } else {
              delete.values[i,] <- cDelete.values
            }
            
            if(doubleRegressionRoutine & all(!is.na(xty.block.values.pos[i,])) & all(!is.na(xtx.block.values.pos[delete.from[i]:delete.to[i],]))){
              xty.delete <- xty.block.values.pos[i,]
              xtx.delete <- xtx.block.values.pos[delete.from[i]:delete.to[i],]
              cDelete.values.pos <- try(solve(xtx.delete, xty.delete, tol=1e-30))
              if(inherits(cDelete.values.pos,"try-error")){
                warning(paste0("Error when solving positive block ",i))
              } else {
                # #unsigned first since we re-use delete.values for the final results
                # delete.values.unsigned[i,1]<-(mPos[i] * delete.values[i,1] + mNeg[i] * (-1) * cDelete.values.neg[1])/(mPos[i]+mNeg[i])
                # delete.values.unsigned[i,2]<-(mPos[i] * delete.values[i,2] + mNeg[i] * (-1) * cDelete.values.neg[2])/(mPos[i]+mNeg[i])
                # delete.values.pos[i,1] <- (mPos[i] * delete.values[i,1] + mNeg[i] * cDelete.values.neg[1])/(mPos[i]+mNeg[i])
                # delete.values.pos[i,2] <- (mPos[i] * delete.values[i,2] + mNeg[i] * cDelete.values.neg[2])/(mPos[i]+mNeg[i])
                delete.values.pos[i,] <- cDelete.values.pos
              }
              #LDvsChiCorrelation.block.values[i] <- (mPos[i] * LDvsChiCorrelation.block.values[i] + mNeg[i] * (-1) * LDvsChiCorrelation.block.values.neg[i])/(mPos[i]+mNeg[i])
              
            }
            
            if(doubleRegressionRoutine & all(!is.na(xty.block.values.neg[i,])) & all(!is.na(xtx.block.values.neg[delete.from[i]:delete.to[i],]))){
              xty.delete <- xty.block.values.neg[i,]
              xtx.delete <- xtx.block.values.neg[delete.from[i]:delete.to[i],]
              cDelete.values.neg <- try(solve(xtx.delete, xty.delete, tol=1e-30))
              if(inherits(cDelete.values.neg,"try-error")){
                warning(paste0("Error when solving negative block ",i))
              } else {
                #unsigned first since we re-use delete.values for the final results
                # delete.values.unsigned[i,1]<-(mPos[i] * delete.values[i,1] + mNeg[i] * (-1) * cDelete.values.neg[1])/(mPos[i]+mNeg[i])
                # delete.values.unsigned[i,2]<-(mPos[i] * delete.values[i,2] + mNeg[i] * (-1) * cDelete.values.neg[2])/(mPos[i]+mNeg[i])
                # delete.values[i,1] <- (mPos[i] * delete.values[i,1] + mNeg[i] * cDelete.values.neg[1])/(mPos[i]+mNeg[i])
                # delete.values[i,2] <- (mPos[i] * delete.values[i,2] + mNeg[i] * cDelete.values.neg[2])/(mPos[i]+mNeg[i])
                delete.values.neg[i,] <- cDelete.values.neg
              }
              #LDvsChiCorrelation.block.values[i] <- (mPos[i] * LDvsChiCorrelation.block.values[i] + mNeg[i] * (-1) * LDvsChiCorrelation.block.values.neg[i])/(mPos[i]+mNeg[i])
              
            }
           
            # bMeanChi2<-mean(abs(merged[select.from[i]:select.to[i],]$ZZ))
            # bMedianChi2<-median(abs(merged[select.from[i]:select.to[i],]$ZZ))
            
            #TEST!! change from LDvsChiCorrelation.block.values2 to LDvsChiCorrelation.block.values (actual corrected variables)
            attenuationCorrectionFactor[i] <- (abs(1/LDvsChiCorrelation.block.values[i]))^(attenuationBiasCorrectionFactorExponent) #careful interpretation of the attenuation factor
            #sign(LDvsChiCorrelation.block.values.2[i])* #should the factor be signed? Probably not - see the condition in Frost et al.
            # cat("\n",i,attenuationCorrectionFactor[i],bMeanChi2,bMedianChi2,delete.values[i,])
            # if(i %in% c(6,10,141,65,22,209)){
            #   cMerged<-merged[select.from[i]:select.to[i],]
            #   m.Chi.max<-mean(cMerged$weights_cov,na.rm=T)*1
            #   if(!preweight.ChiSquare) {
            #     png(file.path(p$folderpath.plots,paste0(trait.names[j],"_",trait.names[k],"_block",i,"LDChi2.png")), width = 8, height = 5, units = 'in', res = 300)
            #     plot(x = mLD[select.from[i]:select.to[i],1], y =  mChi[select.from[i]:select.to[i],], main = paste0("",i))
            #     abline(a = delete.values[i,2], b = delete.values[i,1])
            #     abline(h = m.Chi.max)
            #     abline(h = -(m.Chi.max))
            #     dev.off()
            #   }
            #   png(file.path(p$folderpath.plots,paste0(trait.names[j],"_",trait.names[k],"_block",i,"LDChi2Corr.png")), width = 8, height = 5, units = 'in', res = 300)
            #   if(preweight.ChiSquare) {
            #     plot(x = mLD[select.from[i]:select.to[i],1], y =  mChi[select.from[i]:select.to[i],], main = paste0("",i))
            #     abline(a = delete.values[i,2], b = delete.values[i,1])
            #     abline(h = m.Chi.max)
            #     abline(h = -(m.Chi.max))
            #   } else {
            #     plot(x = mLD[select.from[i]:select.to[i],1]*abs(mChi[select.from[i]:select.to[i],]), y = mChi[select.from[i]:select.to[i],]*abs(mChi[select.from[i]:select.to[i],]), main = paste0("",i))
            #   }
            #   dev.off()
            # }
            #plot(x = mLD[select.from[i]:select.to[i],1]*mChi[select.from[i]:select.to[i],], y = mChi[select.from[i]:select.to[i],], main = paste0("",i))
            #plot(x = mLD[select.from[i]:select.to[i],1]*mChi[select.from[i]:select.to[i],], y = mChi[select.from[i]:select.to[i],]*mChi[select.from[i]:select.to[i],], main = paste0("",i))
            #plot(x = log(mLD[select.from[i]:select.to[i],1]*mChi[select.from[i]:select.to[i],]), y =  log(mChi[select.from[i]:select.to[i],]), main = paste0("",i))
            
          }
          
          overallAttenuationCorrectionFactor <- mean(attenuationCorrectionFactor,na.rm = T)
          
          #TEST - reverse intercept to Chi-square uncorrected intercept scale
          if(preweight.ChiSquare) delete.values[,2] <- delete.values[,2] / mean(merged[is.finite(ZZ),]$ZZ,na.rm=T)
          
          if(correctAttenuationBias == T){
            #multiply with the correction factor
            delete.values[,1] <- delete.values[,1] * attenuationCorrectionFactor
            #delete.values[,1] <- delete.values[,1] * ifelse(delete.values[,1]>0,attenuationCorrectionFactor,1)
            #delete.values[,1] <- delete.values[,1] * ifelse(delete.values[,1]>0, overallAttenuationCorrectionFactor, 1)
            #delete.values[,1] <- delete.values[,1] * overallAttenuationCorrectionFactor
            #TODO - correct the intercept to correspond to the new slope - is this needed?
          }
          
          #fallback for non-dr runs
          delete.values.pos<-delete.values.neg<-matrix(NA,nrow(delete.values),ncol(delete.values))
          
          LOG("Correlation bias correction factor (mean) =",round(overallAttenuationCorrectionFactor,digits = 3))
          # #TEST!!! restore Chi2 weighting
          # if(chi2Weight){
          #   delete.values[,1] <- delete.values[,1]
          # }
          
          
          #TEST of using weighted mean and median as central tendency measure in covariance
          #weighted mean
          
          #jackknife.cov <- cov.wt(delete.values,wt = blockWeight)/n.blocksToUse
          
          # jackknife.cov<-matrix(NA,nrow = 2, ncol = 2)
          # jackknife.cov[1,1]<-(sum((delete.values[,1] - weighted.mean(delete.values[,1],blockWeight))^2)/nrow(delete.values))/n.blocksToUse
          # jackknife.cov[2,2]<-(sum((delete.values[,2] - weighted.mean(delete.values[,2],blockWeight))^2)/nrow(delete.values))/n.blocksToUse
          # jackknife.cov[1,2]<-jackknife.cov[2,1]<-(sum((delete.values[,1] - weighted.mean(delete.values[,1],blockWeight))*(delete.values[,2] - weighted.mean(delete.values[,2],blockWeight)))/nrow(delete.values))/n.blocksToUse
          
          #median
          # jackknife.cov<-matrix(NA,nrow = 2, ncol = 2)
          # jackknife.cov[1,1]<-(sum((delete.values[,1] - median(delete.values[,1]))^2)/nrow(delete.values))/n.blocksToUse
          # jackknife.cov[2,2]<-(sum((delete.values[,2] - median(delete.values[,2]))^2)/nrow(delete.values))/n.blocksToUse
          # jackknife.cov[1,2]<-jackknife.cov[2,1]<-(sum((delete.values[,1] - median(delete.values[,1]))*(delete.values[,2] - median(delete.values[,2])))/nrow(delete.values))/n.blocksToUse
          
          N.vec[1, s] <- N.bar
          
          #pos
          jackknife.cov.pos <- cov(delete.values.pos)/n.blocksToUse #Should be n.blocks -1 ? #should be n.blocksToUse^2 ??? we do not correct this
          jackknife.se.pos <- sqrt(diag(jackknife.cov.pos))
          intercept.se.pos <- jackknife.se.pos[length(jackknife.se.pos)]
          coef.cov.pos <- jackknife.cov.pos[1:n.annot,1:n.annot]/(N.bar^2)
          cat.cov.pos <- coef.cov.pos*(M.tot %*% t(M.tot))
          tot.cov.pos <- sum(cat.cov.pos)
          tot.se.pos <- sqrt(tot.cov.pos)
          lV.pos[s]<-list(delete.values.pos[,1])
          
          #neg
          jackknife.cov.neg <- cov(delete.values.neg)/n.blocksToUse #Should be n.blocks -1 ? #should be n.blocksToUse^2 ??? we do not correct this
          jackknife.se.neg <- sqrt(diag(jackknife.cov.neg))
          intercept.se.neg <- jackknife.se.neg[length(jackknife.se.neg)]
          coef.cov.neg <- jackknife.cov.neg[1:n.annot,1:n.annot]/(N.bar^2)
          cat.cov.neg <- coef.cov.neg*(M.tot %*% t(M.tot))
          tot.cov.neg <- sum(cat.cov.neg)
          tot.se.neg <- sqrt(tot.cov.neg)
          lV.neg[s]<-ifelse(
            j == k,
            list(delete.values.pos[,1]),
            list(delete.values.neg[,1])) #set diagonal elements to equal the results from the positive regressions
          
          #old
          jackknife.cov <- cov(delete.values)/n.blocksToUse #Should be n.blocks -1 ?
          jackknife.se <- sqrt(diag(jackknife.cov))
          intercept.se <- jackknife.se[length(jackknife.se)]
          coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
          cat.cov <- coef.cov*(M.tot %*% t(M.tot))
          tot.cov <- sum(cat.cov)
          tot.se <- sqrt(tot.cov)
          lV[s]<-list(delete.values[,1])
          
  
          #TODO -fix so it works with n.annot
          reg.pos <- t(as.matrix(data.frame(v1=mean(delete.values.pos[is.finite(delete.values.pos[,1]),1]),v2=mean(delete.values.pos[is.finite(delete.values.pos[,2]),2]))))
          reg.neg <- t(as.matrix(data.frame(v1=mean(delete.values.neg[is.finite(delete.values.neg[,1]),1]),v2=mean(delete.values.neg[is.finite(delete.values.neg[,2]),2]))))
          reg <- t(as.matrix(data.frame(v1=mean(delete.values[is.finite(delete.values[,1]),1]),v2=mean(delete.values[is.finite(delete.values[,2]),2]))))
          #reg <- t(as.matrix(data.frame(v1=median(delete.values[is.finite(delete.values[,1]),1]),v2=median(delete.values[is.finite(delete.values[,2]),2])))) #test with median
          
          #convert coeficients to heritability
          intercept <- reg[[2]]
          intercept.pos <- reg.pos[[2]]
          intercept.neg <- reg.neg[[2]]
          
          #coefs <- reg[[1]]/N.bar
          reg.tot <- (reg[[1]]/N.bar)*M.tot #replaced m with M.tot, they should be the same
          reg.tot.pos <- (reg.pos[[1]]/N.bar)*M.tot
          reg.tot.neg <- (reg.neg[[1]]/N.bar)*M.tot
          
          
          #signed consensus
          intercept.signed <- (mPosTotAvg*replaceNa(intercept.pos)+mNegTotAvg*replaceNa(intercept.neg))/(mPosTotAvg+mNegTotAvg)
          intercept.se.signed <- (mPosTotAvg*replaceNa(intercept.se.pos)+mNegTotAvg*replaceNa(intercept.se.neg))/(mPosTotAvg+mNegTotAvg) #do we want to do statistical inference here rather than weighted average? do we know the correlations between the variables?
          reg.tot.signed <- (mPosTotAvg*replaceNa(reg.tot.pos)+mNegTotAvg*replaceNa(reg.tot.neg))/(mPosTotAvg+mNegTotAvg)
          tot.se.signed <- (mPosTotAvg*replaceNa(tot.se.pos)+mNegTotAvg*replaceNa(tot.se.neg))/(mPosTotAvg+mNegTotAvg)
          if(hasNegativeResultsCheck){
            lV.signed[s]<-list(
              ((mPosTotAvg*delete.values.pos+mNegTotAvg*delete.values.neg)/(mPosTotAvg+mNegTotAvg))[,1]
            )
          } else {
            lV.signed[s]<-list(delete.values.pos[,1]) #we assume there to always be some positive results
          }
          
          #unsigned consensus
          intercept.unsigned <- (mPosTotAvg*abs(replaceNa(intercept.pos))+mNegTotAvg*abs(replaceNa(intercept.neg)))/(mPosTotAvg+mNegTotAvg)
          intercept.se.unsigned <- (mPosTotAvg*replaceNa(intercept.se.pos)+mNegTotAvg*replaceNa(intercept.se.neg))/(mPosTotAvg+mNegTotAvg) #should equal the signed se
          reg.tot.unsigned <- (mPosTotAvg*abs(replaceNa(reg.tot.pos))+mNegTotAvg*abs(replaceNa(reg.tot.neg)))/(mPosTotAvg+mNegTotAvg)
          tot.se.unsigned <- (mPosTotAvg*replaceNa(tot.se.pos)+mNegTotAvg*replaceNa(tot.se.neg))/(mPosTotAvg+mNegTotAvg) #should equal the signed se
          if(hasNegativeResultsCheck){
            lV.unsigned[s]<-list(
              ((mPosTotAvg*abs(delete.values.pos)+mNegTotAvg*abs(delete.values.neg))/(mPosTotAvg+mNegTotAvg))[,1]
            )
          } else {
            lV.unsigned[s]<-list(delete.values.pos[,1]) #we assume there to always be some positive results
          }
          
          
          #old
          # if(any(!is.na(delete.values.unsigned))){
          #   jackknife.cov.unsigned <- cov(delete.values.unsigned)/n.blocksToUse #Original
          #   jackknife.se.unsigned <- sqrt(diag(jackknife.cov.unsigned))
          #   intercept.se.unsigned <- jackknife.se.unsigned[length(jackknife.se.unsigned)]
          #   coef.cov.unsigned <- jackknife.cov.unsigned[1:n.annot,1:n.annot]/(N.bar^2)
          #   cat.cov.unsigned <- coef.cov.unsigned*(M.tot %*% t(M.tot))
          #   tot.cov.unsigned <- sum(cat.cov.unsigned)
          #   tot.se.unsigned <- sqrt(tot.cov.unsigned)
          #   lV.unsigned[s]<-list(delete.values.unsigned[,1])
          #   
          #   
          #   #TODO -fix so it works with n.annot
          #   reg.unsigned <- t(as.matrix(data.frame(v1=mean(delete.values.unsigned[is.finite(delete.values.unsigned[,1]),1]),v2=mean(delete.values.unsigned[is.finite(delete.values.unsigned[,2]),2]))))
          #   
          #   #convert coeficients to heritability
          #   intercept <- reg[[2]]
          #   intercept.unsigned <- reg.unsigned[[2]]
          #   coefs <- reg[[1]]/N.bar
          #   coefs.unsigned <- reg.unsigned[[1]]/N.bar
          #   reg.tot <- coefs*M.tot #replaced m with M.tot, they should be the same
          #   reg.tot.unsigned <- coefs.unsigned*M.tot
          # }
  
        }
        
        if(j==k & is.na(pop.prev)==F & is.na(samp.prev)==F){
          conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
          Liab.S[j] <- conversion.factor

          if(verbose){
            LOG("     ", print = FALSE)
            LOG("Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait ",
                chi1, " when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.")
            LOG("     ", print = FALSE)
          }
        }
        
        cov.pos[k, j] <- cov.pos[j, k] <- reg.tot.pos
        cov.neg[k, j] <- cov.neg[j, k] <- reg.tot.neg
        cov[k, j] <- cov[j, k] <- reg.tot
        cov.signed[k, j] <- cov.signed[j, k] <- reg.tot.signed
        cov.unsigned[k, j] <- cov.unsigned[j, k] <- reg.tot.unsigned
        I.pos[k, j] <- I.pos[j, k] <- intercept.pos
        I.neg[k, j] <- I.neg[j, k] <- intercept.neg
        I[k, j] <- I[j, k] <- intercept
        I.signed[k, j] <- I.signed[j, k] <- intercept.signed
        I.unsigned[k, j] <- I.unsigned[j, k] <- intercept.unsigned
        I.pos.SE[k, j] <- I.pos.SE[j, k] <- intercept.se.pos
        I.neg.SE[k, j] <- I.neg.SE[j, k] <- intercept.se.neg
        I.SE[k, j] <- I.SE[j, k] <- intercept.se
        I.signed.SE[k, j] <- I.signed.SE[j, k] <- intercept.se.signed
        I.unsigned.SE[k, j] <- I.unsigned.SE[j, k] <- intercept.se.unsigned
        cov.pos.p[k, j] <- cov.pos.p[j, k] <- 2 * pnorm(abs(reg.tot.pos / tot.se.pos), lower.tail = FALSE) #this is the same as when computed for the liability scale
        cov.neg.p[k, j] <- cov.neg.p[j, k] <- 2 * pnorm(abs(reg.tot.neg / tot.se.neg), lower.tail = FALSE) #this is the same as when computed for the liability scale
        cov.p[k, j] <- cov.p[j, k] <- 2 * pnorm(abs(reg.tot / tot.se), lower.tail = FALSE) #this is the same as when computed for the liability scale
        cov.p.signed[k, j] <- cov.p.signed[j, k] <- 2 * pnorm(abs(reg.tot.signed / tot.se.signed), lower.tail = FALSE) #this is the same as when computed for the liability scale
        cov.p.unsigned[k, j] <- cov.p.unsigned[j, k] <- 2 * pnorm(abs(reg.tot.unsigned / tot.se.unsigned), lower.tail = FALSE) #this is the same as when computed for the liability scale
        cov.blocks[k, j] <- cov.blocks[j, k] <- n.blocksToUse
        if(verbose){
          if(any(!is.na(delete.values.pos)) | any(!is.na(delete.values.neg))){
            
            LOG("Results following")
            LOG("Mean abs(Z*Z): ", round(mean(abs(merged$ZZ)), 4), ", Median abs(Z*Z): ", round(median(abs(merged$ZZ)), 4))
            LOG("Cross trait Intercept, aligned/opposing: ", round(intercept.pos, 4), " (", round(intercept.se.pos, 4), ") / ",round(intercept.neg, 4), " (", round(intercept.se.neg, 4), ")")
            LOG("Total Observed Scale Genetic Covariance (g_cov), aligned/opposing: ", round(reg.tot.pos, 4), " (", round(tot.se.pos, 4), ") / ", round(reg.tot.neg, 4), " (", round(tot.se.neg, 4), ")")
            LOG("g_cov Z, aligned/opposing: ", round(reg.tot.pos / tot.se.pos,3), " / ", round(reg.tot.neg / tot.se.neg,3))
            LOG("g_cov P-value, aligned/opposing: ", round(cov.pos.p[k, j], 6), " / ", round(cov.neg.p[k, j], 6))
            
          } else {
            
            LOG("Results following")
            LOG("Mean abs(Z*Z): ", round(mean(abs(merged$ZZ)), 4), ", Median abs(Z*Z): ", round(median(abs(merged$ZZ)), 4))
            LOG("Cross trait Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")")
            LOG("Total Observed Scale Genetic Covariance (g_cov): ", round(reg.tot, 4), " (", round(tot.se, 4), ")")
            LOG("g_cov Z: ", round(reg.tot / tot.se,3))
            LOG("g_cov P-value: ", round(cov.p[k, j], 6))
            
          }
        }
        
        #LD square sums per block
        LD.SS[s]<-list(LD.block.SS)
        #number of vars (m) per block
        nVar.per.block[s]<-list(nVars.block)
        
      
        #}
        
        ### Total count
        s <- s + 1
      }
      
      #mod addition
      gc() #do garbage collect if this can help with out of memory issues.
      
      if(limited) break
    } #cell loop
    
    
  
    lV.lengths <- c()
    #THESE SHOULD BE THE SAME LENGTHS
    for(i in 1:ncol(N.vec)){
      if(!is.null(lV.pos)) lV.lengths[i]<-length(lV.pos[[i]]) else if(!is.null(lV.neg)) lV.lengths[i]<-length(lV.neg[[i]]) else if(!is.null(lV)) lV.lengths[i]<-length(lV[[i]])
    }
    
    mV.pos <- matrix(data=NA,nrow = max(lV.lengths,na.rm = T), ncol = ncol(N.vec))
    mV.neg <- matrix(data=NA,nrow = max(lV.lengths,na.rm = T), ncol = ncol(N.vec))
    mV <- matrix(data=NA,nrow = max(lV.lengths,na.rm = T), ncol = ncol(N.vec))
    mV.signed <- matrix(data=NA,nrow = max(lV.lengths,na.rm = T), ncol = ncol(N.vec))
    mV.unsigned <- matrix(data=NA,nrow = max(lV.lengths,na.rm = T), ncol = ncol(N.vec))
    for(i in 1:ncol(N.vec)){
      #i<-1
      mV.pos[,i]<-padList(lV.pos[[i]],padding = NA_real_, targetLength = nrow(mV.pos))
      mV.neg[,i]<-padList(lV.neg[[i]],padding = NA_real_, targetLength = nrow(mV.neg))
      mV[,i]<-padList(lV[[i]],padding = NA_real_, targetLength = nrow(mV))
      mV.signed[,i]<-padList(lV.signed[[i]],padding = NA_real_, targetLength = nrow(mV))
      mV.unsigned[,i]<-padList(lV.unsigned[[i]],padding = NA_real_, targetLength = nrow(mV))
    }
    
    # *****S specific
    ## cov is already scaled to Mldsr, Nbar, and B(n blocks)
    ### Scale S to liability:
    ratio.liability <- tcrossprod(sqrt(Liab.S))
    S <- cov * ratio.liability
    S.pos <- cov.pos * ratio.liability
    S.neg <- cov.neg * ratio.liability
    S.signed <- cov.signed * ratio.liability
    S.unsigned <- cov.unsigned * ratio.liability
    
    diag(S.neg)<-diag(S.pos) #special treatment for the negative results as there can't be negative diagonal covariance
    diag(cov.neg.p)<-diag(cov.pos.p)
    
    
    #name traits according to trait.names argument
    #use general format of V1-VX if no names provided
    colnames(S) <- rownames(S) <- if (is.null(trait.names)) paste0("V", 1:ncol(S)) else trait.names
    colnames(S.pos) <- rownames(S.pos) <- if (is.null(trait.names)) paste0("V", 1:ncol(S.pos)) else trait.names
    colnames(S.neg) <- rownames(S.neg) <- if (is.null(trait.names)) paste0("V", 1:ncol(S.neg)) else trait.names
    colnames(S.signed) <- rownames(S.signed) <- colnames(S.pos)
    colnames(S.unsigned) <- rownames(S.unsigned) <- colnames(S.pos)
    
    #more names
    colnames(I)<-rownames(I)<-colnames(S)
    colnames(I.SE)<-rownames(I.SE)<-colnames(S)
    colnames(cov.p)<-rownames(cov.p)<-colnames(S)
    colnames(I.pos)<-rownames(I.pos)<-colnames(S.pos)
    colnames(I.pos.SE)<-rownames(I.pos.SE)<-colnames(S.pos)
    colnames(cov.pos.p)<-rownames(cov.pos.p)<-colnames(S.pos)
    colnames(I.neg)<-rownames(I.neg)<-colnames(S.neg)
    colnames(I.neg.SE)<-rownames(I.neg.SE)<-colnames(S.neg)
    colnames(cov.neg.p)<-rownames(cov.neg.p)<-colnames(S.neg)
    colnames(I.signed)<-rownames(I.signed)<-colnames(S.signed)
    colnames(I.signed.SE)<-rownames(I.signed.SE)<-colnames(S.signed)
    colnames(cov.p.signed)<-rownames(cov.p.signed)<-colnames(S.signed)
    colnames(I.unsigned)<-rownames(I.unsigned)<-colnames(S.unsigned)
    colnames(I.unsigned.SE)<-rownames(I.unsigned.SE)<-colnames(S.unsigned)
    colnames(cov.p.unsigned)<-rownames(cov.p.unsigned)<-colnames(S.unsigned)

    colnames(cov.blocks)<-rownames(cov.blocks)<-colnames(S)
    
    #calculate the ratio of the rescaled and original S matrices
    scale.liability.lt <- gdata::lowerTriangle(ratio.liability, diag = TRUE)
    # scale.liability.pos.lt <- gdata::lowerTriangle(ratio.liability, diag = TRUE)
    # scale.liability.neg.lt <- gdata::lowerTriangle(ratio.liability, diag = TRUE)
    # scale.liability.signed.lt <- gdata::lowerTriangle(ratio.liability, diag = TRUE)
    # scale.liability.unsigned.lt <- gdata::lowerTriangle(ratio.liability, diag = TRUE)
    
    #mod additions - initialise S_Stand, set all NA element in S to 0 so not to have NA values
    S_Stand<-matrix(NA, nrow(S), nrow(S))
    S[is.na(S)]<-0 #maybe add an option to not do this?
    S.pos_Stand<-matrix(NA, nrow(S.pos), nrow(S.pos))
    S.pos[is.na(S.pos)]<-0
    S.neg_Stand<-matrix(NA, nrow(S.neg), nrow(S.neg))
    S.neg[is.na(S.neg)]<-0
    S.signed_Stand<-matrix(NA, nrow(S.signed), nrow(S.signed))
    S.signed[is.na(S.signed)]<-0
    S.unsigned_Stand<-matrix(NA, nrow(S.unsigned), nrow(S.unsigned))
    S.unsigned[is.na(S.unsigned)]<-0
    
    ##calculate standardized results to print genetic correlations to log and screen
    #mod change - will run with negative heritabilities
    ratio.standardisation <- tcrossprod(1 / sqrt(abs(diag(S)))) #mod addition  - take the absolute of the diagonal rather than the actual values in case there are negative heritabilities
    if(doubleRegressionRoutine) ratio.standardisation <- tcrossprod(1 / sqrt(abs(diag(S.pos))))
    # ratio.standardisation.neg <- tcrossprod(1 / sqrt(abs(diag(S.neg)))) #this diagonal equals the S.pos
    # ratio.standardisation.signed <- tcrossprod(1 / sqrt(abs(diag(S.signed)))) #this diagonal equals the S.pos
    # ratio.standardisation.unsigned <- tcrossprod(1 / sqrt(abs(diag(S.unsigned)))) #this diagonal equals the S.pos
    S_Stand <- S * ratio.standardisation
    S.pos_Stand <- S.pos * ratio.standardisation
    S.neg_Stand <- S.neg * ratio.standardisation
    S.signed_Stand <- S.signed * ratio.standardisation
    S.unsigned_Stand <- S.unsigned * ratio.standardisation

    colnames(S_Stand)<-rownames(S_Stand)<-colnames(S)
    colnames(S.pos_Stand)<-rownames(S.pos_Stand)<-colnames(S.pos)
    colnames(S.neg_Stand)<-rownames(S.neg_Stand)<-colnames(S.neg)
    colnames(S.signed_Stand)<-rownames(S.signed_Stand)<-colnames(S.signed)
    colnames(S.unsigned_Stand)<-rownames(S.unsigned_Stand)<-colnames(S.unsigned)
    
    #calculate the ratio of the rescaled and original S matrices
    scale.standardisation.lt <- gdata::lowerTriangle(ratio.standardisation, diag = TRUE)
    
    #deprecated
    # ### We need capped absolute covG estimates for the later standardised S.E.s
    # S.abs.capped<-abs(S)
    # S.abs.capped[!is.na(abs(S_Stand)) & abs(S_Stand) < covariance_std_value_lower_limit]<-covariance_std_value_lower_limit*abs(S)[!is.na(abs(S_Stand)) & abs(S_Stand) < covariance_std_value_lower_limit] #these are low capped to avoid close to zero denominator inflation of the test variables
    # 
    # S.abs.capped.unsigned<-abs(S.unsigned)
    # S.abs.capped.unsigned[!is.na(abs(S_Stand.unsigned)) & abs(S_Stand.unsigned) < covariance_std_value_lower_limit]<-covariance_std_value_lower_limit*abs(S.unsigned)[!is.na(abs(S_Stand.unsigned)) & abs(S_Stand.unsigned) < covariance_std_value_lower_limit] #these are low capped to avoid close to zero denominator inflation of the test variables
    
    # *****V specific
    ## Scale V to N per study (assume m constant)
    ### cov.blocks is NOT previously scaled to Mldsr, Nbar, and B(n blocks)
    # /!\ crossprod instead of tcrossprod because N.vec is a one-row matrix
    #mod addition/change - use cov.blocks matrix rather than n.blocks (constant) as this is trait specific
    #!!!we do not use the corrected cov.blocks -1
    #These are on the covariance scale rather than the linear/mean scale (?)
    #original
    scale.nPerStudy <- 1/crossprod(N.vec * (sqrt(gdata::lowerTriangle(cov.blocks, diag = TRUE)) / M.tot))
    scale.nPerStudy.noblocks <- 1/crossprod(N.vec / M.tot)
    #scale.nPerStudy <- as.matrix(M.tot / (sqrt(N.vec) * sqrt(gdata::lowerTriangle(cov.blocks, diag = TRUE))))
    #scale.nPerStudy.noblocks <- as.matrix(M.tot / sqrt(N.vec))
    
    #original
    v.out <- cov(mV,use="pairwise.complete.obs") * scale.nPerStudy
    v.out.pos <- cov(mV.pos,use="pairwise.complete.obs") * scale.nPerStudy
    v.out.neg <- cov(mV.neg,use="pairwise.complete.obs") * scale.nPerStudy
    v.out.signed <- cov(mV.signed,use="pairwise.complete.obs") * scale.nPerStudy
    v.out.unsigned <- cov(mV.unsigned,use="pairwise.complete.obs") * scale.nPerStudy
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    V <- v.out * tcrossprod(scale.liability.lt)
    V.pos <- v.out.pos * tcrossprod(scale.liability.lt)
    V.neg <- v.out.neg * tcrossprod(scale.liability.lt)
    V.signed <- v.out.signed * tcrossprod(scale.liability.lt)
    V.unsigned <- v.out.unsigned * tcrossprod(scale.liability.lt)
    
    #S.E. of S
    r<-nrow(S)
    
    S.SE<-matrix(0, r, r)
    S.SE[lower.tri(S.SE,diag=TRUE)] <-sqrt(diag(V))
    S.SE[upper.tri(S.SE)]<-t(S.SE)[upper.tri(S.SE)]
    colnames(S.SE)<-rownames(S.SE)<-colnames(S)
    
    S.pos.SE<-matrix(0, r, r)
    S.pos.SE[lower.tri(S.pos.SE,diag=TRUE)] <-sqrt(diag(V.pos))
    S.pos.SE[upper.tri(S.pos.SE)]<-t(S.pos.SE)[upper.tri(S.pos.SE)]
    colnames(S.pos.SE)<-rownames(S.pos.SE)<-colnames(S.pos)
    
    S.neg.SE<-matrix(0, r, r)
    S.neg.SE[lower.tri(S.neg.SE,diag=TRUE)] <-sqrt(diag(V.neg))
    S.neg.SE[upper.tri(S.neg.SE)]<-t(S.neg.SE)[upper.tri(S.neg.SE)]
    colnames(S.neg.SE)<-rownames(S.neg.SE)<-colnames(S.neg)
    
    S.signed.SE<-matrix(0, r, r)
    S.signed.SE[lower.tri(S.signed.SE,diag=TRUE)] <-sqrt(diag(V.signed))
    S.signed.SE[upper.tri(S.signed.SE)]<-t(S.signed.SE)[upper.tri(S.signed.SE)]
    colnames(S.signed.SE)<-rownames(S.signed.SE)<-colnames(S.signed)
    
    S.unsigned.SE<-matrix(0, r, r)
    S.unsigned.SE[lower.tri(S.unsigned.SE,diag=TRUE)] <-sqrt(diag(V.unsigned))
    S.unsigned.SE[upper.tri(S.unsigned.SE)]<-t(S.unsigned.SE)[upper.tri(S.unsigned.SE)]
    colnames(S.unsigned.SE)<-rownames(S.unsigned.SE)<-colnames(S.unsigned)
    
    
    #liability scale p-values - NEW
    cov.p.liab<-matrix(0, r, r)
    cov.p.liab[lower.tri(cov.p.liab,diag=TRUE)]<-2 * pnorm(abs(S[lower.tri(S,diag=TRUE)] / S.SE[lower.tri(S.SE,diag=TRUE)]), lower.tail = FALSE)
    cov.p.liab[upper.tri(cov.p.liab)]<-t(cov.p.liab)[upper.tri(cov.p.liab)]
    colnames(cov.p.liab)<-rownames(cov.p.liab)<-colnames(S)
    
    cov.pos.p.liab<-matrix(0, r, r)
    cov.pos.p.liab[lower.tri(cov.pos.p.liab,diag=TRUE)]<-2 * pnorm(abs(S.pos[lower.tri(S.pos,diag=TRUE)] / S.pos.SE[lower.tri(S.pos.SE,diag=TRUE)]), lower.tail = FALSE)
    cov.pos.p.liab[upper.tri(cov.pos.p.liab)]<-t(cov.pos.p.liab)[upper.tri(cov.pos.p.liab)]
    colnames(cov.pos.p.liab)<-rownames(cov.pos.p.liab)<-colnames(S.pos)
    
    cov.neg.p.liab<-matrix(0, r, r)
    cov.neg.p.liab[lower.tri(cov.neg.p.liab,diag=TRUE)]<-2 * pnorm(abs(S.neg[lower.tri(S.neg,diag=TRUE)] / S.neg.SE[lower.tri(S.neg.SE,diag=TRUE)]), lower.tail = FALSE)
    cov.neg.p.liab[upper.tri(cov.neg.p.liab)]<-t(cov.neg.p.liab)[upper.tri(cov.neg.p.liab)]
    colnames(cov.neg.p.liab)<-rownames(cov.neg.p.liab)<-colnames(S.neg)
    
    cov.signed.p.liab<-matrix(0, r, r)
    cov.signed.p.liab[lower.tri(cov.signed.p.liab,diag=TRUE)]<-2 * pnorm(abs(S.signed[lower.tri(S.signed,diag=TRUE)] / S.signed.SE[lower.tri(S.signed.SE,diag=TRUE)]), lower.tail = FALSE)
    cov.signed.p.liab[upper.tri(cov.signed.p.liab)]<-t(cov.signed.p.liab)[upper.tri(cov.signed.p.liab)]
    colnames(cov.signed.p.liab)<-rownames(cov.signed.p.liab)<-colnames(S.signed)
    
    cov.unsigned.p.liab<-matrix(0, r, r)
    cov.unsigned.p.liab[lower.tri(cov.unsigned.p.liab,diag=TRUE)]<-2 * pnorm(abs(S.unsigned[lower.tri(S.unsigned,diag=TRUE)] / S.unsigned.SE[lower.tri(S.unsigned.SE,diag=TRUE)]), lower.tail = FALSE)
    cov.unsigned.p.liab[upper.tri(cov.unsigned.p.liab)]<-t(cov.unsigned.p.liab)[upper.tri(cov.unsigned.p.liab)]
    colnames(cov.unsigned.p.liab)<-rownames(cov.unsigned.p.liab)<-colnames(S.unsigned)
    
    #liability scale s.e./variance p-values - NEW - NOT FINISHED! - we may use the same test as for the difference in variance.
    
    #calculation of reference constants for tests (c)
    # var1 <- abs((p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S.SE^2)/p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S)
    # var2 <- abs((p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S.SE^2)/(p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S*p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$cov.blocks))
    # median(var1[lower.tri(var1,diag = T)])
    # median(var2[lower.tri(var2,diag = T)])
    
    #one-sided only for now!!!
    cov.SE.p.liab<-matrix(0, r, r)
    cov.SE.p.liab[lower.tri(cov.p.liab,diag=TRUE)]<- pchisq(
      q = ((S.SE[lower.tri(S.SE,diag=TRUE)])^2)/((cov.blocks-1)[lower.tri(cov.blocks,diag=TRUE)])-(cov.SE.p.liab.test.cDivS*S)[lower.tri(S,diag=TRUE)],
      df = (cov.blocks-1)[lower.tri(cov.blocks,diag=TRUE)], lower.tail = FALSE)
    
    cov.SE.p.liab[upper.tri(cov.SE.p.liab)]<-t(cov.SE.p.liab)[upper.tri(cov.SE.p.liab)]
    colnames(cov.SE.p.liab)<-colnames(S)
    rownames(cov.SE.p.liab)<-colnames(S)
    cov.SE.p.liab.unsigned<-cov.SE.p.liab #temporary
    
    ###mod addition - also compute standardised S.E's - have better chi-squared properties (Pearson's formula/test) - NOT USED!!!!
    #ONLY FOR THE V DIAGONAL AS OF NOW
    #V.std <- matrix(data = NA, nrow = ncol(mV), ncol=ncol(mV))
    cov.out.std.vec <- matrix(data = NA, nrow = 1, ncol=ncol(mV))
    cov.vec <- gdata::lowerTriangle(cov, diag = TRUE)
    #S.abs.capped.vec <- gdata::lowerTriangle(S.abs.capped, diag = TRUE) #we don't need this now!!!
    
    scale.nPerStudy.noblocks.vec <- diag(scale.nPerStudy.noblocks)
    #scale.liability.lt.vec<-diag(tcrossprod(scale.liability.lt))
    
    #this is on the linear/mean scale
    #the per-study scaling is on the covariance scale, so we need the square root of this
    mV.scaled <- sqrt(scale.nPerStudy.noblocks.vec) * scale.liability.lt * mV
    
    cov.out.std.vec <- colSums((mV.scaled - cov.vec)^2, na.rm = T)/(diag(V)*gdata::lowerTriangle(cov.blocks, diag = TRUE))

    #old
    # for(j in 1:ncol(mV)){
    #   #j<-1
    #   cov.out.std.vec[j]<-
    #     scale.liability.lt.vec[j] * sum(((mV[,j]-cov.vec[j])^2) ,na.rm = T)/(S.abs.capped.vec[j]) #do not use the scale.nPerStudy as it contains the block-count
    # }
    
    cov.out.std_normal.vec<-(cov.out.std.vec-lV.lengths+1)/sqrt(2*(lV.lengths-1))
    cov.out.std_normal.vec.p<-2*pnorm(abs(cov.out.std_normal.vec), lower.tail = FALSE) #two sided!!
    
    #Pearsons's standardised S.E. of S - for chi-square distributed tests
    r<-nrow(S)
    S.SE.std<-matrix(0, r, r)
    #S.SE.std.unsigned<-matrix(0, r, r) #you can only get these of signed values as of now
    S.SE.std[lower.tri(S.SE.std,diag=TRUE)] <- S.SE.std[upper.tri(S.SE.std,diag=TRUE)] <- sqrt(cov.out.std.vec)
    #S.SE.std.unsigned[lower.tri(S.SE.std.unsigned,diag=TRUE)] <- sqrt(cov.out.std.unsigned.vec)
    #S.SE.std[upper.tri(S.SE.std)]<-t(S.SE.std)[upper.tri(S.SE.std)] #old
    colnames(S.SE.std)<-colnames(S)
    rownames(S.SE.std)<-colnames(S)
    
    #Pearsons's standardised variances. of S - for chi-square distributed tests - as std normal
    r<-nrow(S)
    S.VAR.std_normal<-matrix(0, r, r)
    S.VAR.std_normal[lower.tri(S.VAR.std_normal,diag=TRUE)] <- S.VAR.std_normal[upper.tri(S.VAR.std_normal,diag=TRUE)] <- cov.out.std_normal.vec
    colnames(S.VAR.std_normal)<-colnames(S)
    rownames(S.VAR.std_normal)<-colnames(S)
    
    r<-nrow(S)
    S.VAR.std_normal.p<-matrix(0, r, r)
    S.VAR.std_normal.p[lower.tri(S.VAR.std_normal.p,diag=TRUE)] <- S.VAR.std_normal.p[upper.tri(S.VAR.std_normal.p,diag=TRUE)] <- cov.out.std_normal.vec.p
    colnames(S.VAR.std_normal.p)<-colnames(S)
    rownames(S.VAR.std_normal.p)<-colnames(S)

    #mod additions - initialise V_Stand, set all NA element in V to 0 so not to have NA values
    V_Stand<-matrix(NA, nrow(V), nrow(V))
    V[is.na(V)]<-0
    #rescale the sampling correlation matrix by the appropriate diagonals
    V_Stand <- V * tcrossprod(scale.standardisation.lt)
    
    V.pos_Stand<-matrix(NA, nrow(V.pos), nrow(V.pos))
    V.pos[is.na(V.pos)]<-0
    #rescale the sampling correlation matrix by the appropriate diagonals
    V.pos_Stand <- V.pos * tcrossprod(scale.standardisation.lt)
    
    V.neg_Stand<-matrix(NA, nrow(V.neg), nrow(V.neg))
    V.neg[is.na(V.neg)]<-0
    #rescale the sampling correlation matrix by the appropriate diagonals
    V.neg_Stand <- V.neg * tcrossprod(scale.standardisation.lt)
    
    V.signed_Stand<-matrix(NA, nrow(V.signed), nrow(V.signed))
    V.signed[is.na(V.signed)]<-0
    #rescale the sampling correlation matrix by the appropriate diagonals
    V.signed_Stand <- V.signed * tcrossprod(scale.standardisation.lt)
    
    V.unsigned_Stand<-matrix(NA, nrow(V.unsigned), nrow(V.unsigned))
    V.unsigned[is.na(V.unsigned)]<-0
    #rescale the sampling correlation matrix by the appropriate diagonals
    V.unsigned_Stand <- V.unsigned * tcrossprod(scale.standardisation.lt)
    
    #enter SEs from diagonal of standardized V
    r<-nrow(S)
    S_Stand.SE<-matrix(0, r, r)
    S_Stand.SE[lower.tri(S_Stand.SE,diag=TRUE)] <- S_Stand.SE[upper.tri(S_Stand.SE,diag=TRUE)] <- sqrt(diag(V_Stand))
    colnames(S_Stand.SE)<-rownames(S_Stand.SE)<-colnames(S)
    
    r<-nrow(S.pos)
    S.pos_Stand.SE<-matrix(0, r, r)
    S.pos_Stand.SE[lower.tri(S.pos_Stand.SE,diag=TRUE)] <- S.pos_Stand.SE[upper.tri(S.pos_Stand.SE,diag=TRUE)] <- sqrt(diag(V.pos_Stand))
    colnames(S.pos_Stand.SE)<-rownames(S.pos_Stand.SE)<-colnames(S.pos)
    

    S.neg_Stand.SE<-matrix(0, r, r)
    S.neg_Stand.SE[lower.tri(S.neg_Stand.SE,diag=TRUE)] <- S.neg_Stand.SE[upper.tri(S.neg_Stand.SE,diag=TRUE)] <- sqrt(diag(V.neg_Stand))
    colnames(S.neg_Stand.SE)<-rownames(S.neg_Stand.SE)<-colnames(S.neg)
    

    S.signed_Stand.SE<-matrix(0, r, r)
    S.signed_Stand.SE[lower.tri(S.signed_Stand.SE,diag=TRUE)] <- S.signed_Stand.SE[upper.tri(S.signed_Stand.SE,diag=TRUE)] <- sqrt(diag(V.signed_Stand))
    colnames(S.signed_Stand.SE)<-rownames(S.signed_Stand.SE)<-colnames(S.signed)
    

    S.unsigned_Stand.SE<-matrix(0, r, r)
    S.unsigned_Stand.SE[lower.tri(S.unsigned_Stand.SE,diag=TRUE)] <- S.unsigned_Stand.SE[upper.tri(S.unsigned_Stand.SE,diag=TRUE)] <- sqrt(diag(V.unsigned_Stand))
    colnames(S.unsigned_Stand.SE)<-rownames(S.unsigned_Stand.SE)<-colnames(S.unsigned)
    
    
    if(doubleRegressionRoutine){
      
      LOG(c("     ", "     "), print = FALSE)
      LOG("(Liability) Scale Results")
      
      for(j in 1:n.traits){
        if(is.null(trait.names)){
          chi1<-traits[j]
        }else{chi1 <- trait.names[j]}
        for(k in j:length(traits)){
          if(j == k){
            LOG("     ", print = FALSE)
            LOG("Total Liability Scale h2 for: ", chi1,": ", round(S.pos[j, j], 3), " (", round(S.pos.SE[j, j], 4), ")"," [", round(S.pos.SE[k, j]/S.pos[k, j], 3),"]")
          } #round(sign(S.VAR.std_normal[k, j])*sqrt(abs(S.VAR.std_normal[k, j]))
          
          if(j != k){
            if(is.null(trait.names)){
              chi2<-traits[k]
            }else{chi2 <- trait.names[k]}
            
            LOG("Total Liability Scale covG for ", chi1, " and ",chi2,", aligned/opposing/consensus(signed): ", round(S.pos[k, j], 3), " (", round(S.pos.SE[k, j], 4), ")"," [", round(S.pos.SE[k, j]/S.pos[k, j], 3),"]"," / ", round(S.neg[k, j], 3), " (", round(S.neg.SE[k, j], 4), ")"," [", round(S.neg.SE[k, j]/S.neg[k, j], 3),"]"," / ", round(S.signed[k, j], 3), " (", round(S.signed.SE[k, j], 4), ")"," [", round(S.signed.SE[k, j]/S.signed[k, j], 3),"]")
            
            LOG("     ", print = FALSE)
          }
        }
      }
      
      LOG(c("     ", "     "), print = FALSE)
      LOG("Genetic Correlation Results")
      
      for(j in 1:n.traits){
        if(is.null(trait.names)){
          chi1<-traits[j]
        }else{chi1 <- trait.names[j]}
        for(k in j:length(traits)){
          if(j != k){
            if(is.null(trait.names)){
              chi2<-traits[k]
            }else{chi2 <- trait.names[k]}
              
            LOG("Genetic Correlation between ", chi1, " and ", chi2)
            
            # LOG("Signed/unsigned: ",
            #       round(S.signed_Stand[k, j], 4), " (", round(S.signed_Stand.SE[k, j], 4), ") / ",round(S.unsigned_Stand[k, j], 4), " (", round(S.unsigned_Stand.SE[k, j], 4), ")")
            LOG("Aligned/opposing/consensus(signed): ",
                round(S.pos_Stand[k, j], 4), " (", round(S.pos_Stand.SE[k, j], 4), ") / ",round(S.neg_Stand[k, j], 4), " (", round(S.neg_Stand.SE[k, j], 4), ") / ",round(S.signed_Stand[k, j], 4), " (", round(S.signed_Stand.SE[k, j], 4), ")")
            LOG("     ", print = FALSE)
          }
        }
      }
      
      
      
      if(!all(diag(S.pos) > 0)) {
        warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
        LOG("Your genetic covariance matrix includes traits estimated to have a negative heritability.", print = FALSE)
        #LOG("Genetic correlation results could not be computed due to negative heritability estimates.")
      }
      
      
      
    } else {
      #old print routine
      LOG(c("     ", "     "), print = FALSE)
      LOG("(Liability) Scale Results")
      
      for(j in 1:n.traits){
        if(is.null(trait.names)){
          chi1<-traits[j]
        }else{chi1 <- trait.names[j]}
        for(k in j:length(traits)){
          if(j == k){
            LOG("     ", print = FALSE)
            LOG("Total Liability Scale h2 for: ", chi1,": ", round(S[j, j], 3), " (", round(S.SE[j, j], 4), ")"," [", round(S.SE[k, j]/S[k, j], 3),"]"," {", round(sign(S.VAR.std_normal[k, j])*sqrt(abs(S.VAR.std_normal[k, j])), 3),"}")
          }
          
          if(j != k){
            if(is.null(trait.names)){
              chi2<-traits[k]
            }else{chi2 <- trait.names[k]}
            
            LOG("Total Liability Scale covG for ", chi1, " and ",chi2,": ", round(S[k, j], 3), " (", round(S.SE[k, j], 4), ")"," [", round(S.SE[k, j]/S[k, j], 3),"]")
           #round(sign(S.VAR.std_normal[k, j])*sqrt(abs(S.VAR.std_normal[k, j])), 3)
            LOG("     ", print = FALSE)
          }
        }
      }
      
      LOG(c("     ", "     "), print = FALSE)
      LOG("Genetic Correlation Results")
      
      for(j in 1:n.traits){
        if(is.null(trait.names)){
          chi1<-traits[j]
        }else{chi1 <- trait.names[j]}
        for(k in j:length(traits)){
          if(j != k){
            if(is.null(trait.names)){
              chi2<-traits[k]
            }else{chi2 <- trait.names[k]}
            
            
              
            LOG("Genetic Correlation between ", chi1, " and ", chi2, ": ",
                round(S_Stand[k, j], 4), " (", round(S_Stand.SE[k, j], 4), ")")
              
            LOG("     ", print = FALSE)
          }
        }
      }
      
      
      
      if(!all(diag(S) > 0)) {
        warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
        LOG("Your genetic covariance matrix includes traits estimated to have a negative heritability.", print = FALSE)
        #LOG("Genetic correlation results could not be computed due to negative heritability estimates.")
      }
      
    }
    
    end.time <- Sys.time()
    
    total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
    mins <- floor(floor(total.time)/60)
    secs <- floor(total.time-mins*60)
    
    LOG("     ", print = FALSE)
    LOG("LDSC finished running at ", end.time)
    LOG("Running LDSC for all files took ", mins, " minutes and ", secs, " seconds")
    LOG("     ", print = FALSE)
    
    #print stats
    #print(impstats)
    
    
    for(j in 1:n.traits){
      #j<-1
      impstats[j,c("m.cell.mean.imputed")] <- c(mean(mCellstats.imputed.m[j,],na.rm=T))
      impstats[j,c("m.cell.mean.imputed.full")] <- c(mean(mCellstats.imputed.full.m[j,],na.rm=T))
      impstats[j,c("m.cell.mean.imputed.partial")] <- c(mean(mCellstats.imputed.partial.m[j,],na.rm=T))
      impstats[j,c("m.cell.mean.nonimputed")] <- c(mean(mCellstats.nonimputed.m[j,],na.rm=T))
    }
    
    cResults <- list(V=V, S=S, I=I, N=N.vec, m=M.tot, S.SE=S.SE, mgc=mgc, 
                     V_Stand=V_Stand, S_Stand=S_Stand, S_Stand.SE=S_Stand.SE,
                     #cov.p=cov.p, #this is actually always the same as when computed on the liability scale
                     I.SE=I.SE,
                     
                     cov.p.liab=cov.p.liab,
                     cov.SE.p.liab=cov.SE.p.liab,
                     
                     pos=list(
                       V=V.pos,
                       S=S.pos,
                       I=I.pos,
                       S.SE=S.pos.SE,
                       I.SE=I.pos.SE,
                       cov.p.liab=cov.pos.p.liab,
                       V_Stand=V.pos_Stand,
                       S_Stand=S.pos_Stand,
                       S_Stand.SE=S.pos_Stand.SE
                       ),
                     
                     neg=list(
                       V=V.neg,
                       S=S.neg,
                       I=I.neg,
                       S.SE=S.neg.SE,
                       I.SE=I.neg.SE,
                       cov.p.liab=cov.neg.p.liab,
                       V_Stand=V.neg_Stand,
                       S_Stand=S.neg_Stand,
                       S_Stand.SE=S.neg_Stand.SE
                     ),
                     
                     signed=list(
                       V=V.signed,
                       S=S.signed,
                       I=I.signed,
                       S.SE=S.signed.SE,
                       I.SE=I.signed.SE,
                       cov.p.liab=cov.unsigned.p.liab,
                       V_Stand=V.signed_Stand,
                       S_Stand=S.signed_Stand,
                       S_Stand.SE=S.signed_Stand.SE
                     ),
                     
                     unsigned=list(
                       V=V.unsigned,
                       S=S.unsigned,
                       I=I.unsigned,
                       S.SE=S.unsigned.SE,
                       I.SE=I.unsigned.SE,
                       cov.p.liab=cov.unsigned.p.liab,
                       V_Stand=V.unsigned_Stand,
                       S_Stand=S.unsigned_Stand,
                       S_Stand.SE=S.unsigned_Stand.SE
                     ),
                    
                     
                     cov.blocks=cov.blocks,
                     
                     S.SE.std=S.SE.std,S.VAR.std_normal=S.VAR.std_normal,S.VAR.std_normal.p=S.VAR.std_normal.p,
                     blockValues.LDSR_beta=lV,
                     LD.SS=LD.SS,
                     nVar.per.block=nVar.per.block,
                     impstats=impstats, m.imputed=mCellstats.imputed.m, m.imputed.full=mCellstats.imputed.full.m, m.imputed.partial=mCellstats.imputed.partial.m, m.nonimputed=mCellstats.nonimputed.m
                     )
    
    
    
    if(iFold==0){
      primeResults<-cResults
    }
    
    
    if(k.folds>0){
      if(iFold %% 2 == 0) foldResults[[iKfold]]$test<-cResults else foldResults[[iKfold]]$train<-cResults
    }
  
  } #fold loop
  
  if(k.folds>0) primeResults$fold<-foldResults
  
  
  flush(log.file)
  close(log.file)
  
  return(primeResults)
  
}
