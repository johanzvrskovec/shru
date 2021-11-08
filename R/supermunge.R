#WORK IN PROGRESS
#Based on the amazing work by Grotzinger, A. D. et al. Nat. Hum. Behav. 3, 513–525 (2019) and Bulik-Sullivan, B. K. et al. Nat. Genet. 47, 291–295 (2015).



stdGwasColumnNames <- function(columnNames, stopOnMissingEssential=T,
                               c.SNP = c("SNP","PREDICTOR","SNPID","MARKERNAME","MARKER_NAME","SNPTESTID","ID_DBSNP49","RSID","ID","RS_NUMBER","MARKER", "RS", "RSNUMBER", "RS_NUMBERS", "SNP.NAME","SNP ID", "SNP_ID","LOCATIONALID","ASSAY_NAME"),
                               c.A1 = c("A1","ALLELE1","ALLELE_1","INC_ALLELE","EA","A1_EFFECT","REF","EFFECT_ALLELE","RISK_ALLELE","EFFECTALLELE","EFFECT_ALL","REFERENCE_ALLELE","REF_ALLELE","REFERENCEALLELE","EA","ALLELE_1","INC_ALLELE","ALLELE1","A","A_1","CODED_ALLELE","TESTED_ALLELE"),
                               c.A2 = c("A2","ALLELE2","ALLELE_2","OTHER_ALLELE","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA","ALT","A2_OTHER","NONREF_ALLELE","NEFFECT_ALLELE","NEFFECTALLELE","NONEFFECT_ALLELE","OTHER_ALL","OTHERALLELE","NONEFFECTALLELE","ALLELE0","ALLELE_0","ALT_ALLELE","A_0","NONCODED_ALLELE"),
                               #c.EFFECT = c("EFFECT","OR","B","BETA","LOG_ODDS","EFFECTS","SIGNED_SUMSTAT","EST"),
                               c.BETA = c("BETA","B","EFFECT_BETA","EFFECT","EFFECTS","SIGNED_SUMSTAT","EST","GWAS_BETA","EFFECT_A1","EFFECTA1","EFFECT_NW"),
                               c.OR = c("OR","LOG_ODDS","OR","ODDS-RATIO","ODDS_RATIO","ODDSRATIO","OR(MINALLELE)","OR.LOGISTIC","OR_RAN","OR(A1)"),
                               c.SE = c("SE","STDER","STDERR","STD","STANDARD_ERROR","OR_SE","STANDARDERROR", "STDERR_NW","META.SE","SE_DGC","SE.2GC"),
                               c.Z = c("Z","ZSCORE","Z-SCORE","ZSTAT","ZSTATISTIC","GC_ZSCORE","BETAZSCALE"),
                               c.INFO = c("INFO","IMPINFO","IMPQUALITY", "INFO.PLINK", "INFO_UKBB","INFO_UKB"),
                               c.P = c("P","PVALUE","PVAL","P_VALUE","GC_PVALUE","WALD_P","P.VAL","GWAS_P","P-VALUE","P-VAL","FREQUENTIST_ADD_PVALUE","P.VALUE","P_VAL","SCAN-P","P.LMM","META.PVAL","P_RAN","P.ADD","P_BOLT_LMM"),
                               c.N = c("N","WEIGHT","NCOMPLETESAMPLES","TOTALSAMPLESIZE","TOTALN","TOTAL_N","N_COMPLETE_SAMPLES","N_TOTAL","N_SAMPLES","N_ANALYZED","NSAMPLES","SAMPLESIZE","SAMPLE_SIZE","TOTAL_SAMPLE_SIZE","TOTALSAMPLESIZE"),
                               c.N_CAS = c("N_CAS","NCASE","N_CASE","N_CASES","NCAS","NCA","NCASES","CASES","CASES_N","FRQ_A"),
                               c.N_CON = c("N_CON","NCONTROL","N_CONTROL","N_CONTROLS","NCON","NCO","N_CON","NCONTROLS","CONTROLS","CONTROLS_N","FRQ_U"),
                               c.NEF = c("NEF","NEFF","NEFFECTIVE","NE"),
                               #include FRQ_A?
                               c.FRQ = c("FRQ","MAF","AF","CEUAF","FREQ","FREQ1","EAF","FREQ1.HAPMAP","FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU","EFFECT_ALLELE_FREQ","FREQ.A1","F_A","F_U","FREQ_A","FREQ_U","MA_FREQ","MAF_NW","FREQ_A1","A1FREQ","CODED_ALLELE_FREQUENCY","FREQ_TESTED_ALLELE_IN_HRS","EAF_HRC","EAF_UKB"),
                               c.CHR = c("CHR","CH","CHROMOSOME","CHROM","CHR_BUILD38","CHR_BUILD37","CHR_BUILD36","CHR_B38","CHR_B37","CHR_B36","CHR_ID","SCAFFOLD","HG19CHR","CHR.HG19","CHR_HG19","HG18CHR","CHR.HG18","CHR_HG18","CHR_BP_HG19B37","HG19CHRC"),
                               c.BP = c("BP","ORIGBP","POS","POSITION","LOCATION","PHYSPOS","GENPOS","CHR_POSITION","POS_B38","POS_BUILD38","POS_B37","POS_BUILD37","BP_HG19B37","POS_B36","POS_BUILD36","POS.HG19","POS.HG18","POS_HG19","POS_HG18","BP_HG19","BP_HG18","BP.GRCH38","BP.GRCH37","POSITION(HG19)","POSITION(HG18)","POS(B38)","POS(B37)")
){
  #test
  #columnNames<-cSumstats.names
  
  columnNames.upper<-toupper(columnNames)
  #names(columnNames)<-columnNames
  columnNames.orig<-columnNames
  
  columnNames[columnNames.upper %in% c.SNP] <- c.SNP[1]
  columnNames[columnNames.upper %in% c.A1] <- c.A1[1]
  columnNames[columnNames.upper %in% c.A2] <- c.A2[1]
  #columnNames[columnNames.upper %in% c.EFFECT] <- c.EFFECT[1]
  #if(any(columnNames==c.EFFECT[1])) columnNames[columnNames.upper %in% c.Z] <- c.Z[1] else columnNames[columnNames.upper %in% c.Z] <- c.EFFECT[1]
  columnNames[columnNames.upper %in% c.BETA] <- c.BETA[1]
  columnNames[columnNames.upper %in% c.OR] <- c.OR[1] 
  columnNames[columnNames.upper %in% c.Z] <- c.Z[1] 
  columnNames[columnNames.upper %in% c.SE] <- c.SE[1]
  columnNames[columnNames.upper %in% c.INFO] <- c.INFO[1]
  columnNames[columnNames.upper %in% c.P] <- c.P[1]
  columnNames[columnNames.upper %in% c.N] <- c.N[1]
  columnNames[columnNames.upper %in% c.N_CAS] <- c.N_CAS[1]
  columnNames[columnNames.upper %in% c.N_CON] <- c.N_CON[1]
  columnNames[columnNames.upper %in% c.NEF] <- c.NEF[1]
  columnNames[columnNames.upper %in% c.FRQ] <- c.FRQ[1]
  columnNames[columnNames.upper %in% c.CHR] <- c.CHR[1]
  columnNames[columnNames.upper %in% c.BP] <- c.BP[1]
  
  if(stopOnMissingEssential){
    # Stop if any of these columns are not found
    if(!any(columnNames=="SNP")) stop("\nCould not find the 'SNP' column.\n")
    if(!any(columnNames=="A1")) stop("\nCould not find the 'A1' column.\n")
    if(!any(columnNames=="A2")) stop("\nCould not find the 'A2' column.\n")
  }
  
  if(!any(columnNames=="P")) warning("\nCould not find the P-value column. Standard is 'P'.\n")
  if(!any(columnNames=="BETA") & !any(columnNames=="OR" & !any(columnNames=="Z"))) warning("Could not find any effect column.\n")
  if(!any(columnNames=="SNP")) warning("\nCould not find the 'SNP' column.\n")
  if(!any(columnNames=="A1")) warning("\nCould not find the 'A1' column.\n")
  if(!any(columnNames=="A2")) warning("\nCould not find the 'A2' column.\n")
  if(!any(columnNames=="FRQ")) warning("\nCould not find the 'FRQ' column.\n")
  
  # Warn if multiple of these columns are found
  if(sum(columnNames=="SNP")>1) warning("\nMultiple 'SNP' columns found!\n")
  if(sum(columnNames=="P")>1) warning("\nMultiple 'P' columns found!\n")
  if(sum(columnNames=="A1")>1) warning("\nMultiple 'A1' columns found!\n")
  if(sum(columnNames=="A2")>1) warning("\nMultiple 'A2' columns found!\n")
  if(sum(columnNames=="BETA")>1) warning("\nMultiple 'BETA' columns found!\n")
  if(sum(columnNames=="OR")>1) warning("\nMultiple 'OR' columns found!\n")
  if(sum(columnNames=="Z")>1) warning("\nMultiple 'Z' columns found!\n")
  if(sum(columnNames=="FRQ")>1) warning("\nMultiple 'FRQ' columns found!\n")
  
  return(data.frame(std=columnNames,orig=columnNames.orig))
}

#ref, plink chromosome numbering: https://zzz.bwh.harvard.edu/plink/data.shtml
parseSNPColumnAsRSNumber <- function(text){
  #decide if BGENIE SNP format using top 100,000 SNPs
  #TODO this condition may be improved to not rely on the number of variants being >100,000
  #test
  #text<-files[[i]]$SNP
  if(sum(grepl(pattern = "^\\d+:\\w+_\\w+_\\w+", x= head(x = text, n=100000)))>90000){
    #extract and format rs-no
    indexesLengths<-regexec(pattern = "^\\d+:(\\w+)_\\w+_\\w+", text=text)
    matches<-regmatches(text,indexesLengths)
    return(lapply(X = matches, FUN = function(x)paste0("rs",x[2])))
  }
  
  text<-sub(pattern = "^XY:",replacement = "25:",x = text)
  text<-sub(pattern = "^X:",replacement = "23:",x = text)
  text<-sub(pattern = "^Y:",replacement = "24:",x = text)
  text<-sub(pattern = "^MT:",replacement = "26:",x = text)
  text<-sub(pattern = "^chr",replacement = "",x = text)
  text<-sub(pattern = "_",replacement = ":",x = text)
  
  return(text)
}

parseCHRColumn <- function(text){
  text<-sub(pattern = "^XY",replacement = "25",x = text)
  text<-sub(pattern = "^X",replacement = "23",x = text)
  text<-sub(pattern = "^Y",replacement = "24",x = text)
  text<-sub(pattern = "^MT",replacement = "26",x = text)
  text<-sub(pattern = "^chr",replacement = "",x = text)
  return(text)
}

#test
# filePaths = p$munge$filesToUse,
# refFilePath = p$filepath.SNPReference.1kg,
# #refFilePath = p$filepath.SNPReference.hm3,
# #mask = mask,
# traitNames = p$munge$traitNamesToUse,
# setChangeEffectDirectionOnAlleleFlip = T, #T=same behaviour as genomic SEM
# produceCompositeTable = F,
# N = p$munge$NToUse,
# OLS=p$munge$OLSToUse,
# linprob=p$munge$linprobToUse,
# se.logit = p$munge$se.logitToUse,
# prop=p$munge$propToUse,
# pathDirOutput = p$folderpath.data.sumstats.munged,
# invertEffectDirectionOn = c("ANXI04") #ANXI4 has an inverted effect direction for some reason


#test2
# list_df = lfGwasList
# ref_df = ref
# traitNames = paste0(cModel$code,".F",1:length(lfGwasList))
# setChangeEffectDirectionOnAlleleFlip = T #T=same behaviour as genomic SEM
# #N - precomputed for each SNP in the earlier processing step
# pathDirOutput = project$folderpath.data.sumstats.munged


#test3
# filePaths = project$sumstats.sel$cleanedpath[23]
# refFilePath=project$filepath.SNPReference.1kg
# traitNames=project$sumstats.sel$code[23]
# setChangeEffectDirectionOnAlleleFlip=T #set to TRUE to emulate genomic sem munge
# produceCompositeTable=T
# N=project$sumstats.sel$n_total[23]
# prop=(project$sumstats.sel$n_case/project$sumstats.sel$n_total)[23]
# OLS=project$sumstats.sel$dependent_variable.OLS[23]
# linprob=project$sumstats.sel$dependent_variable.linprob[23]
# se.logit=project$sumstats.sel$se.logit[23]

supermunge <- function(
  list_df=NULL,
  filePaths=NULL,
  ref_df=NULL,
  refFilePath=NULL,
  traitNames=NULL,
  setChangeEffectDirectionOnAlleleFlip=T, #set to TRUE to emulate genomic sem munge
  produceCompositeTable=F,
  N=NULL,
  forceN=T,
  prop=NULL,
  OLS=NULL,
  linprob=NULL,
  se.logit=NULL,
  pathDirOutput=".",
  keepIndel=T,
  harmoniseAllelesToReference=F,
  doChrSplit=F,
  doStatistics=F,
  mask=NULL,
  stopOnMissingEssential=T,
  maxSNPDistanceBpPadding=0,
  invertEffectDirectionOn=NULL,
  process=T,
  standardiseEffectsToExposure=F,
  writeOutput=T,
  info.filter=NULL,
  frq.filter=NULL,
  mhc.filter=NULL #can be either 37 or 38 for filtering the MHC region according to either grch37 or grch38
){
  
  timeStart <- Sys.time()
  
  if(length(list_df)>0){
    if(is.null(traitNames)){
      traitNames<-names(list_df)
    }
    
    if(!is.null(mask)){
      list_df<-list_df[mask]
      if(!is.null(traitNames)) traitNames<-traitNames[mask]
    }
    
    ds.length <- length(list_df)
  } else {
    if(is.null(traitNames)){
      traitNames<-basename(filePaths)
    }
    
    if(!is.null(mask)){
      filePaths<-filePaths[mask]
      if(!is.null(traitNames)) traitNames<-traitNames[mask]
    }
    
    ds.length <- length(filePaths)
  }
  
  
  #settings similar to GenomicSEM sumstats function
  if(is.null(OLS)){
    OLS<-rep(FALSE,ds.length)
  }
  
  if(is.null(linprob)){
    linprob<-rep(FALSE,ds.length)
  }
  
  if(is.null(se.logit)){
    se.logit<-rep(FALSE,ds.length)
  }
  
  cat("\n\n\nS U P E R ★ M U N G E\n")
  
  cat("\n--------------------------------\nSettings:")
  cat("\nkeepIndel=",keepIndel)
  cat("\nharmoniseAllelesToReference=",harmoniseAllelesToReference)
  cat("\nmaxSNPDistanceBpPadding=",maxSNPDistanceBpPadding)
  cat("\nchangeEffectDirectionOnAlleleFlip=",setChangeEffectDirectionOnAlleleFlip)
  if(length(invertEffectDirectionOn)>0){
    cat("\ninvertEffectDirectionOn=", paste(invertEffectDirectionOn,sep = ","))
  }
  cat("\n--------------------------------\n")
  
  ref<-NULL
  if(!is.null(ref_df)){
    ref<-ref_df
    cat(paste0("\nUsing reference from provided dataframe.\n"))
  } else if(!is.null(refFilePath)){
    cat(paste0("\nReading reference file...\n"))
    ref <- read.table(refFilePath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
    cat(paste0("\nRead reference file:\n",refFilePath))
  } 
  
  variantTable<-NA
  if(!is.null(ref)){
    # Transform to data table
    ref <- setDT(ref)
    # Column harmonisation
    ref.keys<-c('SNP')
    ref$SNP <- tolower(as.character(ref$SNP))
    ref$A1 <- toupper(as.character(ref$A1))
    ref$A2 <- toupper(as.character(ref$A2))
    if('CHR' %in% names(ref)) {
      ref$CHR <- toupper(as.character(ref$CHR))
      ref.keys<-c(ref.keys,'CHR') 
    }
    if('BP' %in% names(ref)) {
      ref$BP <- as.integer(ref$BP)
      ref.keys<-c(ref.keys,'BP')
    }
    
    variantTable<-NA
    if(produceCompositeTable){
      #create composite variant table
      variantTable<-ref
      variantTable <- setDT(variantTable)
      setkeyv(variantTable, cols = ref.keys)
      #check keys with key(variantTable)
    } else variantTable<-c()
    
    #rename reference columns as to distinguish them from the dataset columns
    names(ref)<-paste0(names(ref),"_REF")
    setkeyv(ref, cols = paste0(ref.keys,"_REF"))
    #check keys with key(ref)
    
  } else {
    warning("\nRunning without reference.\n")
  }
  
  sumstats.meta<-data.table(name=traitNames,file_path=filePaths,n_snp_raw=NA_integer_,n_snp_res=NA_integer_)

  for(iFile in 1:ds.length){
    #for testing!
    #iFile=1
    
    #set changeEffectDirectionOnAlleleFlip -this has to be reset for each file/dataset
    changeEffectDirectionOnAlleleFlip<-NULL
    if(!is.null(setChangeEffectDirectionOnAlleleFlip)){
      changeEffectDirectionOnAlleleFlip<-setChangeEffectDirectionOnAlleleFlip
    }
    
    timeStart.ds <- Sys.time()
    
    
    if(!is.null(list_df)){
      cat(paste("\n\nSupermunging\t",traitNames[iFile],"\n @ dataset", iFile,"\n"))
      cSumstats <- list_df[[iFile]]
    } else {
      cFilePath<-filePaths[iFile]
      cat(paste("\n\nSupermunging\t",traitNames[iFile],"\nFile:", cFilePath,"\n"))
      cSumstats <- read.table(cFilePath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
    }
    cat("\nReading.")
    
    cSumstats.meta<-data.table(message=NA_character_,display=NA_character_)
    cSumstats.warnings<-list()
    cSumstats.nSNP.raw<-nrow(cSumstats)
    sumstats.meta[iFile,c("n_snp_raw")]<-cSumstats.nSNP.raw
    cSumstats.meta<-rbind(cSumstats.meta,list("# Input variant rows",as.character(cSumstats.nSNP.raw)))
    if(!is.null(ref)) {
      cSumstats.meta<-rbind(cSumstats.meta,list("# Reference variants",as.character(nrow(ref))))
    }
    cat(".")
    
    # Give sumstats new standardised column names
    cSumstats.names <- stdGwasColumnNames(columnNames = names(cSumstats), stopOnMissingEssential = stopOnMissingEssential)
    cSumstats.names.string <-""
    #apply(cSumstats.names, MARGIN = 1, FUN = function(c){cat(c[2],"\t-> ",c[1],"\n")})
    for(iName in 1:nrow(cSumstats.names)){
      cSumstats.names.string<-paste(paste(cSumstats.names$orig,"\t->",cSumstats.names$std), collapse = '\n')
    }
    colnames(cSumstats) <- cSumstats.names$std
    cat(".")
    
    
    #deal with duplicate columns - use the first occurrence
    ##SNP
    iDup<-grep(pattern = "^SNP$",colnames(cSumstats))
    if(length(iDup)>1){
      iDup<-iDup[2:length(iDup)]
      colnames(cSumstats)[iDup]<-"XSNP"
    }
    
    ##BP
    if('BP' %in% names(cSumstats)){
      iDup<-grep(pattern = "^BP$",colnames(cSumstats))
      if(length(iDup)>1){
        iDup<-iDup[2:length(iDup)]
        colnames(cSumstats)[iDup]<-"XBP"
      }
    }
    
    ##FRQ
    if('FRQ' %in% names(cSumstats)){
      iDup<-grep(pattern = "^FRQ$",colnames(cSumstats))
      if(length(iDup)>1){
        iDup<-iDup[2:length(iDup)]
        colnames(cSumstats)[iDup]<-"XFRQ"
      }
    }
    
    ##P
    if('P' %in% names(cSumstats)){
      iDup<-grep(pattern = "^P$",colnames(cSumstats))
      if(length(iDup)>1){
        iDup<-iDup[2:length(iDup)]
        colnames(cSumstats)[iDup]<-"XP"
      }
    }
    
    # Transform to data table
    cSumstats <- setDT(cSumstats)
    # Column harmonisation
    cSumstats.keys<-c('SNP')
    cSumstats$SNP <- as.character(cSumstats$SNP)
    cSumstats$A1 <- toupper(as.character(cSumstats$A1))
    cSumstats$A2 <- toupper(as.character(cSumstats$A2))
    if('Z' %in% names(cSumstats)) cSumstats$Z<-as.numeric(cSumstats$Z)
    if('FRQ' %in% names(cSumstats)) cSumstats$FRQ<-as.numeric(cSumstats$FRQ)
    if('MAF' %in% names(cSumstats)) cSumstats$MAF<-as.numeric(cSumstats$MAF)
    if('INFO' %in% names(cSumstats)) cSumstats$INFO<-as.numeric(cSumstats$INFO)
    if('EFFECT' %in% names(cSumstats)) cSumstats$EFFECT<-as.numeric(cSumstats$EFFECT)
    if('SE' %in% names(cSumstats)) cSumstats$SE<-as.numeric(cSumstats$SE)
    if('BETA' %in% names(cSumstats)) cSumstats$BETA<-as.numeric(cSumstats$BETA)
    if('OR' %in% names(cSumstats)) cSumstats$OR<-as.numeric(cSumstats$OR)
    cat(".")
    
    #parse SNP if needed
    cSumstats$SNP<-tolower(parseSNPColumnAsRSNumber(cSumstats$SNP))
    cat(".")
    
    if('CHR' %in% names(cSumstats)) {
      cSumstats$CHR <- toupper(parseCHRColumn(as.character(cSumstats$CHR)))
      cSumstats.keys<-c(cSumstats.keys,'CHR') 
    }
    cat(".")
    if('BP' %in% names(cSumstats)) {
      cSumstats$BP <- as.integer(cSumstats$BP)
      cSumstats.keys<-c(cSumstats.keys,'BP')
    }
    cat(".")
    
    #set keys when columns stable
    setkeyv(cSumstats,cols = cSumstats.keys)
    cat(".")
    
    #Set effect to column standard - EFFECT
    #TODO Fix according to new column standard
    if(!("EFFECT" %in% colnames(cSumstats))){
      if("BETA" %in% colnames(cSumstats)) cSumstats$EFFECT<-cSumstats$BETA else
        if("OR" %in% colnames(cSumstats)) cSumstats$EFFECT<-cSumstats$OR
    }
    cat(".")
    
    #Remove MHC region based on position
    #references
    #https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
    #https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
    if(!is.null(mhc.filter)){
      if(any(colnames(cSumstats)=="CHR") & any(colnames(cSumstats)=="BP")){
        cSumstats.nSNP<-nrow(cSumstats)
        if(mhc.filter==37) cSumstats <- cSumstats[!is.na(CHR) & !is.na(BP) & CHR=="6" & BP>=28477797 & BP<=33448354, ] else if (mhc.filter==38) cSumstats <- cSumstats[!is.na(CHR) & !is.na(BP) & CHR=="6" & BP>=28510120 & BP<=33480577, ] else cSumstats.warnings<-c(cSumstats.warnings,"Invalid assembly version provided - no filtering of the MHC was done!")
        cSumstats.meta<-rbind(cSumstats.meta,list(paste("Removed variants; GRCh",mhc.filter,"MHC"),as.character(cSumstats.nSNP-nrow(cSumstats))))
      } else {
        cSumstats.warnings<-c(cSumstats.warnings,"No chromosome or base-pair position information available - no filtering of the MHC was done!")
      }
    }
    cat(".")
    
    #Filter variants MAF<frq.filter
    if(!is.null(frq.filter)){
      if("FRQ" %in% names(cSumstats)){
        rm <- (!is.na(cSumstats$FRQ) & ((cSumstats$FRQ<frq.filter & cSumstats$FRQ<0.5) | (1-cSumstats$FRQ)<frq.filter))
        cSumstats <- cSumstats[!rm, ]
        cSumstats.meta<-rbind(cSumstats.meta,list(paste("Removed variants; MAF <",frq.filter),as.character(sum(rm))))
      } else {
        cSumstats.warnings<-c(cSumstats.warnings,"The dataset does not contain a FRQ or MAF column to apply the specified filter on.")
      }
    }
    cat(".")
    
    #Filter variants INFO<info.filter
    if(!is.null(info.filter)){
      if("INFO" %in% names(cSumstats)){
        rm <- (!is.na(cSumstats$INFO) & cSumstats$INFO<info.filter)
        cSumstats <- cSumstats[!rm, ]
        cSumstats.meta<-rbind(cSumstats.meta,list(paste("Removed variants; INFO <",info.filter),as.character(sum(rm))))
      } else {
        cSumstats.warnings<-c(cSumstats.warnings,"The dataset does not contain an INFO column to apply the specified filter on.")
      }
    }
    cat(".")
    
    if(process){
      cat("\nProcessing.")
      # QC, and data management before merge with reference
      
      ## Remove SNPs with missing P
      if("P" %in% colnames(cSumstats)) {
        cSumstats.n<-nrow(cSumstats)
        cSumstats<-cSumstats[which(!is.na(cSumstats$P)),]
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; missing P",as.character(cSumstats.n-nrow(cSumstats))))
      }
      cat(".")
      
      ## Remove SNPs with missing effects
      if("EFFECT" %in% colnames(cSumstats)) {
        cSumstats.n<-nrow(cSumstats)
        cSumstats<-cSumstats[which(!is.na(cSumstats$EFFECT)),]
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; missing EFFECT",as.character(cSumstats.n-nrow(cSumstats))))
      }
      cat(".")
      
      ##Alleles, deal with indels
      if(keepIndel == T){
        cSumstats$A1 <- as.character(toupper(cSumstats$A1))
        cSumstats$A2 <- as.character(toupper(cSumstats$A2))
      } else if(keepIndel == F){
        cSumstats$A1 <- as.character(toupper(cSumstats$A1), c("A", "C", "G", "T"))
        cSumstats$A2 <- as.character(toupper(cSumstats$A2), c("A", "C", "G", "T"))
        cSumstats.meta<-rbind(cSumstats.meta,list("Discarded indels (A1)",as.character(count(is.na(cSumstats$A1)))))
        cSumstats.meta<-rbind(cSumstats.meta,list("Discarded indels (A2)",as.character(count(is.na(cSumstats$A2)))))
      }
      cat(".")
      
      ## Remove duplicated variants across SNP, A1 and A2
      if(any(colnames(cSumstats)=="SNP") & any(colnames(cSumstats)=="A1") & any(colnames(cSumstats)=="A2")) {
        cSumstats.n <- nrow(cSumstats)
        cSumstats <- unique(cSumstats,by = c("SNP","A1","A2"))
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed SNPs; duplicate SNP,A1,A2",as.character(cSumstats.n-nrow(cSumstats))))
      }
      cat(".")
      
      # Merge with reference
      if(!is.null(ref)){
        #Aligning and validating with reference file
        cSumstats.n<-nrow(cSumstats)
        #cat("\nValidating dataset \tnSNP =",cSumstats.n,"\nwith reference \t\tnSNP =", nrow(ref))
        
        
        #Join with reference on SNP rsID, only keeping SNPs with rsIDs part of the reference
        #https://stackoverflow.com/questions/34644707/left-outer-join-with-data-table-with-different-names-for-key-variables/34645997#34645997
        
        cSumstats.merged.snp<-ref[cSumstats, on=c(SNP_REF='SNP'), nomatch=0]
        #cSumstats[ref, on=c(SNP='SNP_REF'), CHR_REF:=i.CHR_REF] #experimental
        #cSumstats <- merge(ref,cSumstats,by="SNP",all.x=F,all.y=F) #old
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed SNPs; rsID not in ref",as.character(cSumstats.n-nrow(cSumstats.merged.snp))))
        #cat("\nRemoved SNPs with rsIDs not present in the reference:\t\t",cSumstats.n-nrow(cSumstats.merged.snp))
        cat("..")
        
        if('CHR' %in% names(cSumstats.merged.snp) && 'CHR_REF' %in% names(cSumstats.merged.snp))
        {
          cSumstats.merged.snp.n<-nrow(cSumstats.merged.snp)
          cSumstats.merged.snp<-cSumstats.merged.snp[CHR==CHR_REF]
          #cat("\nRemoved SNPs not matching the reference chromosome:\t\t",cSumstats.merged.snp.n-nrow(cSumstats.merged.snp))
          cSumstats.meta<-rbind(cSumstats.meta,list("Removed SNPs; CHR not matching ref",as.character(cSumstats.merged.snp.n-nrow(cSumstats.merged.snp))))
        }
        cat(".")
        
        if('BP' %in% names(cSumstats.merged.snp) && 'BP_REF' %in% names(cSumstats.merged.snp)){
          cSumstats.merged.snp.maxAlleleLength<-max(nchar(cSumstats.merged.snp$A1),nchar(cSumstats.merged.snp$A2),nchar(cSumstats.merged.snp$A1_REF),nchar(cSumstats.merged.snp$A2_REF))
          cSumstats.merged.snp.n<-nrow(cSumstats.merged.snp)
          cSumstats.merged.snp<-cSumstats.merged.snp[BP < BP_REF + cSumstats.merged.snp.maxAlleleLength + maxSNPDistanceBpPadding & BP > BP_REF - cSumstats.merged.snp.maxAlleleLength - maxSNPDistanceBpPadding]
          #cat("\nRemoved SNPs outside the specified bp window +-bp",(cSumstats.merged.snp.maxAlleleLength + maxSNPDistanceBpPadding -1),"from the reference position:\t\t",cSumstats.merged.snp.n-nrow(cSumstats.merged.snp))
          cSumstats.meta<-rbind(cSumstats.meta,list(paste("Removed SNPs; outside bp window +-bp",as.character(cSumstats.merged.snp.maxAlleleLength + maxSNPDistanceBpPadding -1)),as.character(cSumstats.merged.snp.n-nrow(cSumstats.merged.snp))))
        }
        cat(".")
        
        #replace missing columns
        cSumstats.merged.snp$SNP<-cSumstats.merged.snp$SNP_REF
        
        cSumstats.merged.pos<-NULL
        if('CHR' %in% names(cSumstats) && 'BP' %in% names(cSumstats) && 'CHR_REF' %in% names(ref) && 'BP_REF' %in% names(ref)) {
          #Join with reference on position rather than rsID
          #cSumstats.merged.pos<-ref[cSumstats, on=c(CHR_REF='CHR' , 'BP' < BP_REF + cSumstats.merged.snp.maxAlleleLength + maxSNPDistanceBpPadding & 'BP' > BP_REF - cSumstats.merged.snp.maxAlleleLength - maxSNPDistanceBpPadding)]
          
          cSumstats.merged.pos<-ref[cSumstats, on=c(CHR_REF='CHR' , BP_REF='BP'), nomatch=0]
          #cSumstats.merged.pos<-cSumstats.merged.pos[!(cSumstats.merged.pos$SNP_REF %in% cSumstats.merged.snp$SNP_REF),]
          
          #replace missing columns
          cSumstats.merged.pos$CHR<-cSumstats.merged.pos$CHR_REF
          cSumstats.merged.pos$BP<-cSumstats.merged.pos$BP_REF
        }
        cat("..")
        
        #using cSumstats to store the plain or merged result
        if(is.null(cSumstats.merged.pos)){
          cSumstats<-cSumstats.merged.snp
        } else {
          #merge merged datasets
          cSumstats.merged.pos.salvaged<-cSumstats.merged.pos[!(cSumstats.merged.pos$SNP_REF %in% cSumstats.merged.snp$SNP_REF),]
          cSumstats<-rbindlist(list(cSumstats.merged.snp,cSumstats.merged.pos.salvaged), use.names=T)
          #cat("\nSalvaged SNPs by merging on SNP position rather than rsID:",nrow(cSumstats.merged.pos.salvaged))
          cSumstats.meta<-rbind(cSumstats.meta,list("Salvaged SNPs by ref locus",as.character(nrow(cSumstats.merged.pos.salvaged))))
          cSumstats.merged.pos<-NULL
          #TODO Make more memory friendly version using https://www.biostars.org/p/432389/
        }
        cSumstats.merged.snp<-NULL
      }
      cat(".")
      
      # More QC and data management, after merge with reference
      
      #store original allele order and frequency info
      cSumstats$A1_ORIG<-cSumstats$A1
      cSumstats$A2_ORIG<-cSumstats$A2
      if(any(colnames(cSumstats)=="FRQ")) cSumstats$FRQ_ORIG<-cSumstats$FRQ
      
      
      if(!is.null(ref)){
        
        ##Synchronise SNP with reference SNP
        cSumstats$SNP<-cSumstats$SNP_REF
        
        ## Add in chr and bp from ref if not present in datasets
        if(!any(colnames(cSumstats)=="CHR") & any(colnames(cSumstats)=="CHR_REF")) {
            cSumstats$CHR<-cSumstats$CHR_REF
            cSumstats.keys<-c(cSumstats.keys,'CHR')
        }
        if(!any(colnames(cSumstats)=="BP") & any(colnames(cSumstats)=="BP_REF")){
            cSumstats$BP<-cSumstats$BP_REF
            cSumstats.keys<-c(cSumstats.keys,'BP') 
        }
        if(!any(colnames(cSumstats)=="FRQ") & any(colnames(cSumstats)=="FRQ_REF")){
          cSumstats$FRQ<-cSumstats$FRQ_REF
        }
      }
      cat(".")
      
      #restore sumstats data table keys
      setkeyv(cSumstats,cols = cSumstats.keys)
      cat(".")
      
      if(!is.null(ref)){
        ## Remove SNPs where alleles are not matching at least one of the reference alleles
        cSumstats.n<-nrow(cSumstats)
        cond.removeNonmatching<-(cSumstats$A1 != (cSumstats$A1_REF) & cSumstats$A1 != (cSumstats$A2_REF)) & (cSumstats$A2 != (cSumstats$A1_REF)  & cSumstats$A2 != (cSumstats$A2_REF))
        cSumstats<-cSumstats[which(!cond.removeNonmatching), ]
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed SNPs; A1 or A2 not matching any ref allele",as.character(sum(cond.removeNonmatching))))
        sumstats.meta[iFile,c("Removed, nonmatching ref alleles")]<-sum(cond.removeNonmatching)
      }
      cat(".")
      
      ## N
      if(any(colnames(cSumstats)=="N_CAS") && any(colnames(cSumstats)=="N_CON")) {
        ### Calculate total N from number of cases and number of controls if they are present. Overwrite any specific total N.
        cSumstats.meta<-rbind(cSumstats.meta,list("N",paste("N_CAS + N_CON")))
        cSumstats$N <- cSumstats$N_CAS + cSumstats$N_CON
      } else if(!forceN & !is.null(N) & length(N)>=iFile & any(colnames(cSumstats)=="N")){
        num<-sum(is.na(cSumstats$N))
        if(num>0){
        cSumstats[which(is.na(cSumstats$N)),]$N<-N[iFile]
        cSumstats.meta<-rbind(cSumstats.meta,list("N",paste("Set to",N[iFile],"for",num," vars missing N")))
        }
      } else if(!is.null(N) & length(N)>=iFile) {
        #cat("\nUsing the explicitly specified N for the whole dataset:",N[iFile])
        cSumstats.meta<-rbind(cSumstats.meta,list("N",paste("Set to",N[iFile])))
        cSumstats$N<-N[iFile]
      } else if(!("N" %in% colnames(cSumstats))) {
        cSumstats.warnings<-c(cSumstats.warnings,"\nNo N column detected!")
        cSumstats.meta<-rbind(cSumstats.meta,list("N","Warning: Not detected!"))
        cSumstats$N<-NA_integer_
      }
      cat(".")
      
      ## Establish allele order from the reference
      cond.invertedAlleleOrder<-NULL
      if(!is.null(ref)){
        #cond.invertedAlleleOrder<-(cSumstats$A1 != cSumstats$A1_REF & cSumstats$A2 == cSumstats$A1_REF) #the same condition as in GenomicSEM munge.
        #cond.invertedAlleleOrder<-(cSumstats$A2 != cSumstats$A2_REF & cSumstats$A1 == cSumstats$A2_REF)
        cond.invertedAlleleOrder<-((cSumstats$A2 != cSumstats$A2_REF & cSumstats$A1 == cSumstats$A2_REF) | (cSumstats$A1 != cSumstats$A1_REF & cSumstats$A2 == cSumstats$A1_REF)) #experimental - seems to work similar to the GenomicSEM implementation
      }
      cat(".")
      
      ## Invert alleles or harmonise (including their FRQ) to reference
      if(!is.null(ref) & harmoniseAllelesToReference){
        # Fix A1 and A2 to reflect the reference alleles
        cSumstats$A1<-cSumstats$A1_REF
        cSumstats$A2<-cSumstats$A2_REF
        if(any(colnames(cSumstats)=="MAF_REF")) cSumstats$FRQ<-cSumstats$MAF_REF
      } else if(!is.null(cond.invertedAlleleOrder)) { ## or Invert alleles  -FRQ is dealt with below
        cSumstats$A1<-ifelse(cond.invertedAlleleOrder, cSumstats$A2_ORIG, cSumstats$A1)
        cSumstats$A2<-ifelse(cond.invertedAlleleOrder, cSumstats$A1_ORIG, cSumstats$A2)
      }
      cat(".")
      
      # FRQ
      #cond.invertedFRQ<-NULL
      if(any(colnames(cSumstats)=="FRQ")) {
        ### Has FRQ
        #### Check if value is within limits [0,1]
        if(any(cSumstats$FRQ>1) || any(cSumstats$FRQ<0)) {
          stop(paste0('\nThere are FRQ values larger than 1 (',sum(cSumstats$FRQ>1),') or less than 0 (',sum(cSumstats$FRQ<0),') which is outside of the possible FRQ range.'))
        }
        
        ### Invert FRQ based on the previous reference matching
        if(!is.null(cond.invertedAlleleOrder) & !harmoniseAllelesToReference) {
          alleleFRQ <- ifelse(cond.invertedAlleleOrder, (1-cSumstats$FRQ), cSumstats$FRQ)
          if(mean(alleleFRQ)<mean(cSumstats$FRQ)){
            cSumstats$FRQ<-alleleFRQ
            cSumstats.meta<-rbind(cSumstats.meta,list("FRQ","Fitted (flipped) according to reference allele order"))
            sumstats.meta[iFile,c("FRQ.flipped")]<-T
          } else {
            cond.invertedMAF<-cSumstats$FRQ > .5
            cSumstats$FRQ<-ifelse(cond.invertedMAF,1-cSumstats$FRQ,cSumstats$FRQ)
            cSumstats.meta<-rbind(cSumstats.meta,list("FRQ","Set to MAF"))
            sumstats.meta[iFile,c("FRQ.flipped")]<-F
          }
        }
        ### Compute MAF
        cond.invertedMAF<-cSumstats$FRQ > .5
        cSumstats$MAF<-ifelse(cond.invertedMAF,1-cSumstats$FRQ,cSumstats$FRQ)
      } else {
        ### Does not have FRQ
        #### Note that FRQ is not present
        sumstats.meta[iFile,c("no_FRQ")]<-T
        cSumstats.warnings<-c(cSumstats.warnings,"No FRQ column present!")
        if(!is.null(ref)){
          #Set FRQ from ref if not present
          cSumstats$FRQ<-cSumstats$MAF_REF
          cSumstats.warnings<-c(cSumstats.warnings,"Inferring FRQ from reference!")
        } else {
          #### Add empty FRQ here for consistency
          cSumstats$FRQ<-NA_real_
        }
      }
      cat(".")
      
      ## Compute variance of individual variant effects according to 2pq
      cSumstats$VSNP<-2*cSumstats$FRQ*(1-cSumstats$FRQ)
      cat(".")
      
      
      # P validity checks before the Z-score calculation
      ### Check the values of the P-column
      if("P" %in% colnames(cSumstats)) {
        if((sum(cSumstats$P > 1) + sum(cSumstats$P < 0)) > 100){
          cSumstats.warnings<-c(cSumstats.warnings,"\nThe P column contains numerous values outside of the expected bounds [0,1]. This can indicate that the column is misinterpreted.")
        }
      }
      cat(".")
      
      # EFFECT and SE standardisations and corrections according to effect type (OLS, log OR or non log OR)
      ## adapted from the GenomicSEM sumstats-function
      
      #produce the unstandardised regression beta from Z if no EFFECT present
      if(any(colnames(cSumstats)=="Z") && !any(colnames(cSumstats)=="EFFECT")) {
        ## Compute BETA/EFFECT from Z if present
        cSumstats$EFFECT <- cSumstats$Z/sqrt(cSumstats$N * cSumstats$VSNP) 
        cSumstats$SE <- cSumstats$EFFECT/cSumstats$Z
        cSumstats[is.na(cSumstats$SE),]$SE<-1 #explicitly set NA SE to 1
        cSumstats.meta<-rbind(cSumstats.meta,list("BETA","Calculated from Z"))
      }
      cat(".")
      
      if(any(colnames(cSumstats)=="EFFECT")) {
        
        ## Determine effect type, and set effect to log(EFFECT) if odds ratio
        if(round(median(cSumstats$EFFECT,na.rm=T)) == 1) {
          ###is odds ratio
          cSumstats$EFFECT<-log(cSumstats$EFFECT)
          sumstats.meta[iFile,c("effect_type")]<-"OR"
          cSumstats.meta<-rbind(cSumstats.meta,list("EFFECT","OR  =>ln(OR)"))
          if(any(colnames(cSumstats)=="BETA")) {
            cSumstats.warnings<-c(cSumstats.warnings,"The effect format being an ODDS RATIO may not be compatible with the original variable naming scheme!")
            sumstats.meta[iFile,c("effect_type_warning")]<-T
          }
        } else {
          ###is NOT odds ratio
          sumstats.meta[iFile,c("effect_type")]<-"non_OR"
          cSumstats.meta<-rbind(cSumstats.meta,list("EFFECT","NON OR"))
          if(any(colnames(cSumstats)=="OR")) {
            cSumstats.warnings<-c(cSumstats.warnings,"The effect format NOT being an ODDS RATIO may not compatible with the original variable naming scheme!")
            sumstats.meta[iFile,c("effect_type_warning")]<-T
          }
        }
        cat(".")
        
        ## Compute Z score (standardised beta) - used for effect corrections further and is later corrected accordingly
        if(any(colnames(cSumstats)=="P")) {
          if(any(colnames(cSumstats)=="Z")) cSumstats$Z_ORIG<-cSumstats$Z #save original Z-score
          cSumstats$Z <- sign(cSumstats$EFFECT) * sqrt(qchisq(cSumstats$P,1,lower=F))
          cSumstats.meta<-rbind(cSumstats.meta,list("Z","Calculated from P and sign(EFFECT)"))
        } else if(any(colnames(cSumstats)=="SE")){
          if(any(colnames(cSumstats)=="Z")) cSumstats$Z_ORIG<-cSumstats$Z #save original Z-score
          cSumstats$Z <- cSumstats$EFFECT/cSumstats$SE #is this less reliable as we cannot know the scale of SE?
          cSumstats.meta<-rbind(cSumstats.meta,list("Z","Calculated from EFFECT and SE"))
        } else {
          cSumstats.meta<-rbind(cSumstats.meta,list("Z","NOT calculated since no P or SE"))
        }
        cat(".")
        
        ## Inspect new and old Z-values
        if(any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="Z_ORIG")){
          if(mean(cSumstats$Z)-mean(cSumstats$Z_ORIG)>1) cSumstats.warnings<-c(cSumstats.warnings,"New Z differ from old by more than 1sd!")
        }
        
        
          
      } else {
        ## Does not have EFFECT
        ### Note that EFFECT is not present
        cSumstats.warnings<-c(cSumstats.warnings,"No EFFECT column present.")
        sumstats.meta[iFile,c("no_EFFECT")]<-T
      }
      cat(".")
      
      
      
      
      ##invert overall effect if specified
      if(any(colnames(cSumstats)=="EFFECT") & !is.null(invertEffectDirectionOn)){
        if(any(invertEffectDirectionOn==traitNames[iFile])){
          cSumstats$EFFECT<-cSumstats$EFFECT*-1
          cSumstats.meta<-rbind(cSumstats.meta,list("Inverted overall effect",as.character(length(cSumstats$EFFECT))))
        }
      }
      cat(".")
      
      #compute minimum variance for later calculations
      if(any(colnames(cSumstats)=="SE")){
        minv<-min(cSumstats$SE, na.rm = T)^2
      }
      
      #compare hypothesised inverted allele effects with non-inverted allele effects for validation
      if(is.null(changeEffectDirectionOnAlleleFlip)) changeEffectDirectionOnAlleleFlip<-T
      
      if(any(colnames(cSumstats)=="EFFECT") & !is.null(cond.invertedAlleleOrder)){
        if(any(cond.invertedAlleleOrder)){
          sumstats.meta[iFile,c("Inverted allele order variants")]<-sum(cond.invertedAlleleOrder)
          if(any(colnames(cSumstats)=="SE")){
            cSumstats.meta<-rbind(cSumstats.meta,list("Mean effect","ivw"))
            meffects.reference<-weighted.mean(cSumstats$EFFECT[!cond.invertedAlleleOrder], w = 1/(minv + cSumstats$SE[!cond.invertedAlleleOrder]^2), na.rm = T)
            meffects.candidate<-weighted.mean(cSumstats$EFFECT[cond.invertedAlleleOrder], w = 1/(minv + cSumstats$SE[cond.invertedAlleleOrder]^2), na.rm = T)
          } else {
            cSumstats.meta<-rbind(cSumstats.meta,list("Mean effect","plain"))
            meffects.reference<-mean(cSumstats$EFFECT[!cond.invertedAlleleOrder],na.rm = T)
            meffects.candidate<-mean(cSumstats$EFFECT[cond.invertedAlleleOrder],na.rm = T)
          }
          
          meffects.candidate.inverted<-meffects.candidate*-1
          sdeffects.reference<-sd(cSumstats$EFFECT[!cond.invertedAlleleOrder],na.rm = T)
          sdeffects.candidate<-sd(cSumstats$EFFECT[cond.invertedAlleleOrder],na.rm = T)
          cSumstats.meta<-rbind(cSumstats.meta,list("Number variants, reference, candidate:",paste0(as.character(length(cSumstats$EFFECT[!cond.invertedAlleleOrder])),",",as.character(length(cSumstats$EFFECT[cond.invertedAlleleOrder])))))
          cSumstats.meta<-rbind(cSumstats.meta,list("Mean reference effect (sd)",paste0(as.character(round(meffects.reference,digits = 5))," (",round(sdeffects.reference,digits = 5),")")))
          cSumstats.meta<-rbind(cSumstats.meta,list("Mean candidate effect (sd)",paste0(as.character(round(meffects.candidate,digits = 5))," (",round(sdeffects.candidate,digits = 5),")")))
          #cSumstats.meta<-rbind(cSumstats.meta,list("Mean candidate effect, inverted",as.character(round(abs(meffects.candidate.inverted),digits = 5))))
          #cSumstats.meta<-rbind(cSumstats.meta,list("Delta effect, plain",as.character(round(abs(meffects.reference-meffects.candidate),digits = 5))))
          #cSumstats.meta<-rbind(cSumstats.meta,list("Delta effect, inverted",as.character(round(abs(meffects.reference-meffects.candidate.inverted),digits = 5))))
          
          sumstats.meta[iFile,c("Reference variants")]<-length(cSumstats$EFFECT[!cond.invertedAlleleOrder])
          sumstats.meta[iFile,c("Candidate variants")]<-length(cSumstats$EFFECT[cond.invertedAlleleOrder])
          sumstats.meta[iFile,c("Mean reference effect")]<-round(meffects.reference,digits = 5)
          sumstats.meta[iFile,c("Mean reference effect sd")]<-round(sdeffects.reference,digits = 5)
          sumstats.meta[iFile,c("Mean candidate effect")]<-round(meffects.candidate,digits = 5)
          sumstats.meta[iFile,c("Mean candidate effect sd")]<-round(sdeffects.candidate,digits = 5)
          
          if(abs(meffects.reference-meffects.candidate.inverted)>(abs(meffects.reference-meffects.candidate)+1*sdeffects.reference)) {
            cSumstats.meta<-rbind(cSumstats.meta,list("Delta effect","inverted > plain+1sd"))
            sumstats.meta[iFile,c("Delta effect, inverted > plain+1sd")]<-T
            if(is.null(changeEffectDirectionOnAlleleFlip)){
              changeEffectDirectionOnAlleleFlip<-F #inactivate the correction of effect direction on allele flip because of less credible new effect mean.
            } else if(changeEffectDirectionOnAlleleFlip){
              cSumstats.warnings<-c(cSumstats.warnings,"\nChange effect direction on allele flip specified, but the mean effect difference between plain and inverted variants is much larger than between plain and non-inverted, indicating that inverted effects may be invalid.")
            }
            
          } else {
            cSumstats.meta<-rbind(cSumstats.meta,list("Delta effect, relationship","inverted < plain+1sd"))
            sumstats.meta[iFile,c("Delta effect, inverted > plain+1sd")]<-F
            if(!is.null(changeEffectDirectionOnAlleleFlip)){
              if(!changeEffectDirectionOnAlleleFlip) cSumstats.warnings<-c(cSumstats.warnings,"\nChange effect direction on allele flip inactivated, but the mean effect difference between untouched and inverted variants is much smaller than between untouched and non-inverted, indicating that inverted effects for these variants may still be valid.")
            }
          }
        }
      }
      cat(".")
      
      ## EFFECT direction
      ##invert effect if inverted allele order
      if(any(colnames(cSumstats)=="EFFECT")& !is.null(cond.invertedAlleleOrder) & changeEffectDirectionOnAlleleFlip) {
        if(any(cond.invertedAlleleOrder)) cSumstats$EFFECT<-ifelse(cond.invertedAlleleOrder,(cSumstats$EFFECT*-1),cSumstats$EFFECT)
        if(any(colnames(cSumstats)=="SE")){
          meffects.new<-weighted.mean(cSumstats$EFFECT, w = 1/(minv + cSumstats$SE^2), na.rm = T)
          cSumstats.meta<-rbind(cSumstats.meta,list("New effect mean",as.character(meffects.new)))
        } else {
          meffects.new<-mean(cSumstats$EFFECT,na.rm = T)
          cSumstats.meta<-rbind(cSumstats.meta,list("New effect mean",as.character(meffects.new)))
        }
      }
      sumstats.meta[iFile,c("changeEffectDirectionOnAlleleFlip")]<-changeEffectDirectionOnAlleleFlip
      cat(".")
      
      
      if(!is.null(cond.invertedAlleleOrder)) cSumstats.meta<-rbind(cSumstats.meta,list(
        paste0(
          "Modified SNPs; inverted allele order [",ifelse(any(colnames(cSumstats)=="FRQ"),"FRQ",""),",",
          ifelse(changeEffectDirectionOnAlleleFlip,"EFFECT",""),"]"),
        as.character(sum(cond.invertedAlleleOrder))))
      
      
      ## Compute Z score (standardised beta) and P or update it according to any corrections just done
      if(any(colnames(cSumstats)=="SE")) {
        cSumstats$Z <- cSumstats$EFFECT/cSumstats$SE
        cSumstats$P <- 2*pnorm(q = abs(cSumstats$Z),mean = 0, sd = 1, lower.tail = F)
      }
      cat(".")
      
    } #end of process stuff
    
    if(standardiseEffectsToExposure & any(colnames(cSumstats)=="EFFECT")) {
      ##standardise outcome standardised BETA to exposure as well (fully standardised)
      ##correct binary trait regressions with a 'unit-variance-liability scaling'
      if(OLS[iFile]){
        #Has OLS unstandardised beta
        cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","OLS"))
        if(any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="N")) {
          cSumstats$EFFECT <- cSumstats$Z/sqrt(cSumstats$N * cSumstats$VSNP) #standardisation
          cSumstats$SE <- abs(cSumstats$EFFECT/cSumstats$Z) #standardisation as this is derived form the standardised EFFECT
          cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Z, N, UVL std => BETA,SE"))
        } else stop("\nCould not compute BETA,SE because of missing Z or N!\n")
      
        
      } else if(linprob[iFile]){
        #Has effect based on a linear estimator for a binary outcome (rather than a logistic model)
        cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Binary, linear"))
        if(any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="N")){
          if(is.null(prop[iFile]) | is.na(prop[iFile])) stop("\nCould not perform correction of linear BETA,SE to liability scale because of missing or invalid prop argument!\n")
          cSumstats$EFFECT <- cSumstats$Z/sqrt(prop[iFile]*(1-prop[iFile]) * cSumstats$N * cSumstats$VSNP) #standardisation
          cSumstats$SE <- 1/sqrt(prop[iFile]*(1-prop[iFile]) * cSumstats$N * cSumstats$VSNP) #standardisation
          cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Z, N,propCaCo, UVL std => BETA,SE"))
        } else stop("\nCould not compute BETA,SE because of missing Z or N!\n")
      
        correctionTerm <- sqrt((cSumstats$EFFECT^2) + (pi^2)/3)
        cSumstats$EFFECT <- cSumstats$EFFECT/correctionTerm #residual variance correction
        cSumstats$SE <- cSumstats$SE/correctionTerm #residual variance correction
       
      } else {
        #Has effect based on a logistic estimator for a binary outcome, OR or logistic beta
        cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Binary, logistic"))
        if(any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="N")) {
          cSumstats$EFFECT <- cSumstats$Z/sqrt(cSumstats$N * cSumstats$VSNP) #standardisation
          cSumstats$SE <- abs(cSumstats$EFFECT/cSumstats$Z) #standardisation as this is derived form the standardised EFFECT
          cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Z, N, UVL std=> BETA,SE"))
        } else stop("\nCould not compute BETA,SE because of missing Z or N!\n")
      
        correctionTerm <- sqrt((cSumstats$EFFECT^2) + (pi^2)/3)
        cSumstats$EFFECT <- cSumstats$EFFECT/correctionTerm #residual variance correction
        #cSumstats$SE <- cSumstats$SE/correctionTerm #residual variance correction
        if(se.logit[iFile]){
          cSumstats$SE <- cSumstats$SE/correctionTerm #residual variance correction
          cSumstats.meta <- rbind(cSumstats.meta,list("SE","logit"))
        } else {
          #transform to logit SE 
          cSumstats$SE <- (cSumstats$SE/exp(cSumstats$EFFECT))/correctionTerm #UV correction
          cSumstats.meta <- rbind(cSumstats.meta,list("SE","SE(OR) => logit"))
        }
        
      }
      cat(".")
      
      ## Compute Z,P again to update it according to any corrections just done
      if(any(colnames(cSumstats)=="SE")) {
        cSumstats$Z <- cSumstats$EFFECT/cSumstats$SE
        cSumstats$P <- 2*pnorm(q = abs(cSumstats$Z),mean = 0, sd = 1, lower.tail = F)
      }
      cat(".")
      
    }
    
    ## Check effect value credibility
    if(any(colnames(cSumstats)=="SE")) {
      if(any(colnames(cSumstats)=="MAF")) {
        # only for common + uncommon (non-rare) SNPs if MAF available
        if(mean(abs(cSumstats[MAF>0.001,EFFECT/SE]), na.rm = T) > 5) cSumstats.warnings<-c(cSumstats.warnings,"\nCommon variant EFFECT/SE ratio >5 which could be a cause of misspecified/misinterpreted arguments!\n")
      } else {
        if(mean(abs(cSumstats[,EFFECT/SE]), na.rm = T) > 5) cSumstats.warnings<-c(cSumstats.warnings,"\nOverall EFFECT/SE ratio >5 which could be a cause of misspecified/misinterpreted arguments!\n")
      }
    }
    
    #NA values check
    if(any(is.na(cSumstats))) cSumstats.warnings<-c(cSumstats.warnings,"\nNull values detected among results!\n")
    
    # output handling below
    
    #cSumstats<-cSumstats[,-which(names(cSumstats) %in% c("A1_REF","A2_REF"))]
    
    #output.colnames<- c("SNP","N","Z","A1","A2")
    output.colnames<- c("SNP")
    output.colnames<- c(output.colnames,c("A1","A2"))
    #output.colnames<- c(output.colnames,c("A1","A2","A1_ORIG","A2_ORIG"))
    if("CHR" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"CHR")
    if("BP" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"BP")
    if("FRQ" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"FRQ")
    #if("MAF" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"MAF")
    output.colnames<- c(output.colnames,c("P"))
    if("EFFECT" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"EFFECT")
    if("SE" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"SE")
    if("Z" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"Z")
    if("P" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"P")
    if("N" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N")
    if("N_CAS" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N_CAS")
    if("N_CON" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N_CON")
    
    output.colnames.more<-colnames(cSumstats)[!(colnames(cSumstats) %in% output.colnames)]
    output.colnames.all<-c(output.colnames,output.colnames.more)
    cSumstats<-subset(cSumstats,select = output.colnames) #only output standardised columns
    cSumstats.meta<-rbind(cSumstats.meta,list("SNPs after supermunge",as.character(nrow(cSumstats))))
    
    
    #merge with variantTable
    if(produceCompositeTable){
      cName.beta <- paste0("BETA.",traitNames[iFile])
      cName.se <- paste0("SE.",traitNames[iFile])
      toJoin <- cSumstats[,c("SNP","EFFECT","SE")]
      #colnames(toJoin) <- c("SNP",cName.beta,cName.se)
      #update variantTable by reference
      #ref https://stackoverflow.com/questions/44433451/r-data-table-update-join
      #ref about data.table improvements https://stackoverflow.com/questions/42537520/data-table-replace-data-using-values-from-another-data-table-conditionally/42539526#42539526
      variantTable[toJoin, on=c(SNP='SNP'), c(cName.beta,cName.se):=.(i.EFFECT,i.SE)]
    }
    
    cat("\nProcessing done!\n\n")
    
    cat("\nDataset columns interpreted as:\n",cSumstats.names.string)
    cat("\nData processing results:\n")
    #print(readr::format_delim(as.data.frame(cSumstats.meta),delim = "\t",col_names = F,quote_escape = F, eol = "\n"))
    apply(cSumstats.meta, MARGIN = 1, FUN = function(x){
      cat(as.character(x[1]),"\t\t\t",as.character(x[2]),"\n")
    })
    
    if(length(cSumstats.warnings)>0){
      cat("\nShowing warnings below.\n")
      lapply(X = cSumstats.warnings, FUN = cat)
    } else {
      cat("\nNo warnings detected.\n")
    }
    
    nfilepath<-file.path(pathDirOutput,traitNames[iFile])
    if(!doChrSplit & writeOutput){
      cat("\nSaving supermunged dataset...\n\n")
      write.table(x = cSumstats,file = nfilepath,sep="\t", quote = FALSE, row.names = F, append = F)
      nfilepath.gzip<-gzip(nfilepath)
      cat(paste("\nSupermunged dataset saved as", nfilepath.gzip, "in the specified output directory."))
    }
    
    #addition: producing per-chromosome files in a folder, as RAISS columns
    if(doChrSplit & writeOutput) {
      cat("\nSaving supermunged dataset...\n\n")
      if("CHR" %in% colnames(cSumstats)){
        dir.create(paste0(nfilepath,".chr"), showWarnings = FALSE)
        #TODO Adapt to new numeric chromosome numbering
        validChromosomes<-c(1:22,"X","Y","XY","MT") #as per Plink standard
        for(chr in validChromosomes){
          output.chr<-output[which(output$CHR==chr),c("SNP","ORIGBP","A1","A2","Z")]
          colnames(output.chr)<-c("rsID","pos","A0","A1","Z")
          write.table(x = output.chr,file = file.path(paste0(nfilepath,".chr"), paste0("z_",traitNames[iFile],"_",chr,".txt")),sep="\t", quote = FALSE, row.names = F, append = F)
        }
      } else stop("\nSplit by chromosome specified, but the dataset does not have a CHR column.")
      
      cat(paste("\nOne file per chromosome have been saved under", paste0(nfilepath,".chr"), "in the specified output directory."))
      
    }
    
    timeStop.ds <- Sys.time()
    timeDiff <- difftime(time1=timeStop.ds,time2=timeStart.ds,units="sec")
    timeDiff.minutes <- floor(floor(timeDiff)/60)
    timeDiff.seconds <- timeDiff-timeDiff.minutes*60
    
    cat("\nSupermunge of ",traitNames[iFile]," was done in",timeDiff.minutes, "minutes and",timeDiff.seconds," seconds.\n")
    gc() #do garbage collect if this can help with out of memory issues.
  }
  
  timeStop <- Sys.time()
  timeDiff <- difftime(time1=timeStop,time2=timeStart,units="sec")
  timeDiff.minutes <- floor(floor(timeDiff)/60)
  timeDiff.seconds <- timeDiff-timeDiff.minutes*60
  
  cat("\nSupermunge of all datasets was done in",timeDiff.minutes, "minutes and",timeDiff.seconds,"seconds.")
  
  return(list(
    meta=as.data.frame(sumstats.meta),
    last=as.data.frame(cSumstats),
    composite=as.data.frame(variantTable)
  )
    )
  
}
