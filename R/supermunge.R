#Johan Zvrskovec, 2021 
#Based on the fantastic work by Grotzinger, A. D. et al. Nat. Hum. Behav. 3, 513–525 (2019) and Bulik-Sullivan, B. K. et al. Nat. Genet. 47, 291–295 (2015).



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
                               c.P = c("P","PVALUE","PVAL","P_VALUE","GC_PVALUE","WALD_P","P.VAL","GWAS_P","P-VALUE","P-VAL","FREQUENTIST_ADD_PVALUE","P.VALUE","P_VAL","SCAN-P","P.LMM","META.PVAL","P_RAN","P.ADD","P_BOLT_LMM","PVAL_ESTIMATE"),
                               c.N = c("N","WEIGHT","NCOMPLETESAMPLES","TOTALSAMPLESIZE","TOTALN","TOTAL_N","N_COMPLETE_SAMPLES","N_TOTAL","N_SAMPLES","N_ANALYZED","NSAMPLES","SAMPLESIZE","SAMPLE_SIZE","TOTAL_SAMPLE_SIZE","TOTALSAMPLESIZE"),
                               c.N_CAS = c("N_CAS","NCASE","N_CASE","N_CASES","NCAS","NCA","NCASES","CASES","CASES_N","FRQ_A"),
                               c.N_CON = c("N_CON","NCONTROL","N_CONTROL","N_CONTROLS","NCON","NCO","N_CON","NCONTROLS","CONTROLS","CONTROLS_N","FRQ_U"),
                               c.NEF = c("NEF","NEFF","NEFFECTIVE","NE"),
                               #include FRQ_A?
                               c.FRQ = c("FRQ","MAF","AF","CEUAF","FREQ","FREQ1","EAF","FREQ1.HAPMAP","FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU","EFFECT_ALLELE_FREQ","FREQ.A1","F_A","F_U","FREQ_A","FREQ_U","MA_FREQ","MAF_NW","FREQ_A1","A1FREQ","CODED_ALLELE_FREQUENCY","FREQ_TESTED_ALLELE","FREQ_TESTED_ALLELE_IN_HRS","EAF_HRC","EAF_UKB","EAF_EUR_UKB","FREQ_TESTED_ALLELE"),
                               c.CHR = c("CHR","CH","CHROMOSOME","CHROM","CHR_BUILD38","CHR_BUILD37","CHR_BUILD36","CHR_B38","CHR_B37","CHR_B36","CHR_ID","SCAFFOLD","HG19CHR","CHR.HG19","CHR_HG19","HG18CHR","CHR.HG18","CHR_HG18","CHR_BP_HG19B37","HG19CHRC","#CHROM"),
                               c.BP = c("BP","BP1","ORIGBP","POS","POSITION","LOCATION","PHYSPOS","GENPOS","CHR_POSITION","POS_B38","POS_BUILD38","POS_B37","POS_BUILD37","BP_HG19B37","POS_B36","POS_BUILD36","POS.HG19","POS.HG18","POS_HG19","POS_HG18","BP_HG19","BP_HG18","BP.GRCH38","BP.GRCH37","POSITION(HG19)","POSITION(HG18)","POS(B38)","POS(B37)"),
                               c.BP2 =c("BP2"),
                               c.DF = c("DF","CHISQ_DF")
){
  #test
  #columnNames<-cSumstats.names
  columnNames<-as.character(columnNames)
  columnNames.upper<-toupper(columnNames)
  #names(columnNames)<-columnNames
  columnNames.orig<-columnNames
  
  columnNames[columnNames.upper %in% c.SNP] <- c.SNP[1]
  columnNames[columnNames.upper %in% c.A1] <- c.A1[1]
  columnNames[columnNames.upper %in% c.A2] <- c.A2[1]
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
  columnNames[columnNames.upper %in% c.BP2] <- c.BP2[1]
  columnNames[columnNames.upper %in% c.DF] <- c.DF[1]
  
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
  
  return(data.frame(std=as.character(columnNames),orig=as.character(columnNames.orig)))
  
}

#ref, plink chromosome numbering: https://zzz.bwh.harvard.edu/plink/data.shtml
#chromosome Un: https://genome.ucsc.edu/FAQ/FAQdownloads.html#download11
parseSNPColumnAsRSNumber <- function(text){
  #decide if BGENIE SNP format using top 100,000 SNPs
  #TODO this condition may be improved to not rely on the number of variants being >100,000
  #test
  #text<-files[[i]]$SNP
  if(sum(grepl(pattern = "^\\d+:\\w+_\\w+_\\w+", x= head(x = text, n=100000)))>90000){
    #extract and format rs-no
    indexesLengths<-regexec(pattern = "^\\d+:(\\w+)_\\w+_\\w+", text=text)
    matches<-regmatches(text,indexesLengths)
    return(unlist(lapply(X = matches, FUN = function(x)paste0("rs",x[2]))))
  }
  
  #rsXXXX:A1:A2 - format
  if(sum(grepl(pattern = "^\\w+:\\w+:\\w+", x= head(x = text, n=100000)))>90000){
    indexesLengths<-regexec(pattern = "^(\\w+):\\w+", text=text)
    matches<-regmatches(text,indexesLengths)
    return(unlist(lapply(X = matches, FUN = function(x) ifelse(is.na(x[2]),x[1],x[2]))))
  }
  
  text<-sub(pattern = "^chr",replacement = "",x = text, ignore.case = T)
  text<-sub(pattern = "^XY:",replacement = "25:",x = text, ignore.case = T)
  text<-sub(pattern = "^X:",replacement = "23:",x = text, ignore.case = T)
  text<-sub(pattern = "^Y:",replacement = "24:",x = text, ignore.case = T)
  text<-sub(pattern = "^MT:",replacement = "26:",x = text, ignore.case = T)
  text<-sub(pattern = "^M:",replacement = "26:",x = text, ignore.case = T)
  text<-sub(pattern = "^Un:",replacement = "0:",x = text, ignore.case = T)
  text<-sub(pattern = "_",replacement = ":",x = text, ignore.case = T)
  
  return(text)
}

#ref, plink chromosome numbering: https://zzz.bwh.harvard.edu/plink/data.shtml
#chromosome Un: https://genome.ucsc.edu/FAQ/FAQdownloads.html#download11
#this does make the output numeric
parseCHRColumn <- function(text){
  text<-trimws(text)
  text<-sub(pattern = "^chr",replacement = "",x = text, ignore.case = T)
  text<-sub(pattern = "^XY",replacement = "25",x = text, ignore.case = T)
  text<-sub(pattern = "^X",replacement = "23",x = text, ignore.case = T)
  text<-sub(pattern = "^Y",replacement = "24",x = text, ignore.case = T)
  text<-sub(pattern = "^MT",replacement = "26",x = text, ignore.case = T)
  text<-sub(pattern = "^M",replacement = "26",x = text, ignore.case = T)
  text<-sub(pattern = "^Un",replacement = "0",x = text, ignore.case = T)
  
  return(as.integer(text))
}

readFile <- function(filePath,nThreads=5){
  return(data.table::fread(file = filePath, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = F))
}

#test
# list_df=list(highld=p$highld_b37)
# chainFilePath = "../data/alignment_chains/hg19ToHg38.over.chain.gz"

#test with hard coded values
# filePaths = "../data/gwas_sumstats/raw/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"
# refFilePath = "../data/variant_lists/hc1kgp3.b38.eur.l2.jz2022.gz"
# rsSynonymsFilePath = "../data/variant_lists/dbsnp151.synonyms.gz"
# traitNames = "BODY11"
# N = 681275
# pathDirOutput = "../data/gwas_sumstats/munged_1kg_eur_supermunge"


# 
# #test with settings from analysis script
# filePaths = p$munge$filesToUse
# refFilePath = p$filepath.SNPReference.1kg
# #rsSynonymsFilePath = p$filepath.rsSynonyms.dbSNP151
# traitNames = p$munge$traitNamesToUse
# imputeFromLD=T
# ##process=F
# #region.imputation.filter_df=p$highld #provide the high-ld regions to not use for imputation
# #produceVariantTable = T
# N = p$munge$NToUse
# pathDirOutput = p$folderpath.data.sumstats.imputed
# #chainFilePath = file.path(p$folderpath.data,"alignment_chains","hg19ToHg38.over.chain.gz")
# test=T
# 

#set default params for test
# list_df=NULL
# filePaths=NULL
# ref_df=NULL
# refFilePath=NULL
# rsSynonymsFilePath=NULL
# ldDirPath=NULL
# chainFilePath = NULL
# traitNames=NULL
# ancestrySetting=c("ANY")
# setChangeEffectDirectionOnAlleleFlip=T #set to TRUE to emulate genomic sem munge
# produceCompositeTable=F
# imputeFromLD=F
# imputeAdjustN=T
# imputeFrameLenBp=500000 #500000 for comparison with SSIMP and ImpG
# imputeFrameLenCM=0.5 #frame size in cM, will override the bp frame length - set to NULL if you want to use the bp-window argument
# N=NULL
# forceN=F
# prop=NULL
# OLS=NULL
# linprob=NULL
# se.logit=NULL
# liftover=NULL
# pathDirOutput="."
# keepIndel=T
# harmoniseAllelesToReference=F
# harmoniseBPToReference=T
# doChrSplit=F
# doStatistics=F
# mask=NULL
# stopOnMissingEssential=F
# maxSNPDistanceBpPadding=0
# invertEffectDirectionOn=NULL
# parse=T
# process=T
# standardiseEffectsToExposure=F
# writeOutput=T
# info.filter=NULL
# maf.filter=NULL
# mhc.filter=NULL #can be either 37 or 38 for filtering the MHC region according to either grch37 or grch38
# region.filter_df=NULL #dataframe with columns CHR,BP1,BP2 specifying regions to be removed, due to high LD for example.
# region.imputation.filter_df=NULL #dataframe with columns CHR,BP1,BP2 specifying regions to be excluded from acting as support for imputation, due to high LD for example.
# GC="none" #"reinflate"
# nThreads = 5
# lossless = F
# test = F


supermunge <- function(
  list_df=NULL,
  filePaths=NULL,
  ref_df=NULL,
  refFilePath=NULL,
  rsSynonymsFilePath=NULL,
  ldDirPath=NULL,
  chainFilePath = NULL, #chain file for lift-over
  traitNames=NULL,
  ancestrySetting=c("ANY"), #ancestry setting list per dataset
  setChangeEffectDirectionOnAlleleFlip=T, #set to TRUE to emulate genomic sem munge
  produceCompositeTable=F,
  imputeFromLD=F,
  imputeAdjustN=F,
  imputeFrameLenBp=500000, #500000 for comparison with SSIMP and ImpG
  imputeFrameLenCM=0.5, #frame size in cM, will override the bp frame length - set to NULL if you want to use the bp-window argument
  N=NULL,
  forceN=F,
  prop=NULL,
  OLS=NULL,
  linprob=NULL,
  se.logit=NULL,
  liftover=NULL,
  pathDirOutput=".",
  keepIndel=T,
  harmoniseAllelesToReference=F,
  harmoniseBPToReference=T,
  doChrSplit=F,
  doStatistics=F,
  mask=NULL,
  stopOnMissingEssential=F,
  maxSNPDistanceBpPadding=5,
  invertEffectDirectionOn=NULL,
  parse=T,
  process=T,
  standardiseEffectsToExposure=F,
  writeOutput=T,
  info.filter=NULL,
  maf.filter=NULL,
  mhc.filter=NULL, #can be either 37 or 38 for filtering the MHC region according to either grch37 or grch38
  region.filter_df=NULL, #dataframe with columns CHR,BP1,BP2 specifying regions to be removed, due to high LD for example.
  region.imputation.filter_df=NULL, #dataframe with columns CHR,BP1,BP2 specifying regions to be excluded from acting as support for imputation, due to high LD for example.
  GC="none", #"reinflate",
  nThreads = 5,
  lossless = F, #If true, include all original and additional columns, otherwise restrict output to standard column set (default)
  test = F #set test mode - for quickly testing the function
){
  
  timeStart <- Sys.time()
  
  if(length(list_df)>0){
    if(is.null(traitNames)){
      if(is.data.frame(list_df)) list_df<-list(trait1=list_df)
      traitNames<-names(list_df)
    }
    
    if(!is.null(mask)){
      list_df<-list_df[mask]
      if(!is.null(traitNames)) traitNames<-traitNames[mask]
    }
    
    nDatasets <- length(list_df)
  } else {
    if(is.null(traitNames)){
      traitNames<-basename(filePaths)
    }
    
    if(!is.null(mask)){
      filePaths<-filePaths[mask]
      if(!is.null(traitNames)) traitNames<-traitNames[mask]
    }
    
    nDatasets <- length(filePaths)
  }
  
  if(length(ancestrySetting)<nDatasets){
    ancestrySetting<-unlist(rep(ancestrySetting,nDatasets))
  }
  
  
  #settings similar to GenomicSEM sumstats function
  ## Considers everything as OLS datasets if nothing specified however
  if(is.null(OLS)){
    OLS<-rep(TRUE,nDatasets)
  }
  
  if(is.null(linprob)){
    linprob<-rep(FALSE,nDatasets)
  }
  
  if(is.null(se.logit)){
    se.logit<-rep(FALSE,nDatasets)
  }
  
  if(is.null(liftover)){
    liftover<-rep(!is.null(chainFilePath),nDatasets)
  }
  
  cat("\n\n\nS U P E R ★ M U N G E\n")
  cat("\n",nDatasets,"dataset(s) provided")
  cat("\n--------------------------------\nSettings:")
  
  cat("\nkeepIndel=",keepIndel)
  cat("\nharmoniseAllelesToReference=",harmoniseAllelesToReference)
  cat("\nchangeEffectDirectionOnAlleleFlip=",setChangeEffectDirectionOnAlleleFlip)
  cat("\nprocess=",process)
  cat("\nimputeFromLD=",imputeFromLD)
  if(imputeFromLD) cat("\nimputeFrameLenBp=",imputeFrameLenBp)
  if(imputeFromLD) cat("\nimputeFrameLenCM=",imputeFrameLenCM)
  cat("\nproduceCompositeTable=",produceCompositeTable)
  if(length(invertEffectDirectionOn)>0) cat("\ninvertEffectDirectionOn=", paste(invertEffectDirectionOn,sep = ","))
  cat("\npathDirOutput=",pathDirOutput)
  cat("\n--------------------------------\n")
  
  #read reference variant list
  ref<-NULL
  if(!is.null(ref_df)){
    ref<-ref_df
    cat("\nUsing reference variants from provided dataframe.\n")
  } else if(!is.null(refFilePath)){
    cat(paste0("\nReading reference variant file..."))
    ref<-fread(file = refFilePath, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = F)
    #ref <- read.table(refFilePath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
    cat(paste0("\nRead reference variant file:\n",refFilePath))
    if(test){
      ref<-ref[1:(nrow(ref)/4)]
    }
  }
  
  #read SNP id synonyms
  idSynonyms<-NA
  if(!is.null(rsSynonymsFilePath)){
    cat("\nReading variant ID synonyms...")
    
    idSynonyms <- fread(file = rsSynonymsFilePath, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",header = F, fill=T,  blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = F, skip = 2, sep = "")
    cat(".")
    #idSynonyms.map<-as.data.frame(matrix(data = NA, nrow = 0, ncol = 0))
    idSynonyms$parts<-strsplit(x = idSynonyms$V1, split = " ")
    cat(".")
    idSynonyms[,V1:=NULL]
    #idSynonyms[,parts2:=paste0("rs",parts)]
    #idSynonyms$parts2<-lapply(idSynonyms$parts,FUN = function(x){paste0("rs",x)})
    idSynonyms$parts.n<-lapply(idSynonyms$parts,FUN = function(x){length(x)})
    idSynonyms[,parts.n:=as.integer(parts.n)]
    idSynonyms$first.rs<-lapply(idSynonyms$parts,FUN = function(x){x[[1]]}) #current
    idSynonyms[,first.rs:=paste0("rs",first.rs)]
    idSynonyms[,first.rs:=as.character(first.rs)]
    setkeyv(idSynonyms, cols = c("first.rs","parts.n"))
    cat(".")
    
    uniqueSynonyms<-data.table(SNP=paste0("rs",unique(unlist(idSynonyms$parts))))
    setkeyv(uniqueSynonyms, cols = c("SNP"))
    cat(".")
    
    # idSynonyms2<-read.table(file = rsSynonymsFilePath,header = F,na.strings = c(".",NA,"NA",""), blank.lines.skip = T, encoding = "UTF-8", fill = T)
    # rm(idSynonyms2)
    cat("Done!\n")
  }
  
  
  if(!is.null(ldDirPath) & is.null(ref)) stop("You must have a specified reference to append LD scores to!")
  
  variantTable<-NA
  if(!is.null(ref)){
    # ref should be a data.table at this point
    
    # Column harmonisation
    ref.keys<-c("SNP","A1","A2")
    #let's assume the ref is properly formatted and interpreted here with PLINK numeric chromosome numbers
    # ref[,SNP:=tolower(as.character(SNP))]
    # ref[,A1:=toupper(as.character(A1))]
    # ref[,A2:=toupper(as.character(A2))]
    if('CHR' %in% names(ref)) {
      ref.keys<-c(ref.keys,'CHR') 
    }
    if('BP' %in% names(ref)) {
      #ref[,BP:=as.integer(BP)]
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
    
    
    #read and merge with ld scores from directory
    #TODO Fix so it can work with any type of ld-score files in a folder rather than a set of chromosomes
    if(!is.null(ldDirPath)){
      cat("\nReading LD-scores from specified directory...")
      ldscores<-c()
      for(iChr in 1:22){
        #iChr<-1
        ldscores[[iChr]]<-suppressMessages(read_delim(
          file.path(ldDirPath, paste0(iChr, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
        ldscM<-suppressMessages(read_csv(file.path(ldDirPath, paste0(iChr, ".l2.M_5_50")), col_names = FALSE))
        ldscores[[iChr]]$M<-ldscM[[1]]
      }
      ldscores<-rbindlist(ldscores)
      setkeyv(ldscores, cols = c("SNP"))
      #ref.allsnps<-data.table(SNP=unique(ref[,c("SNP")]))
      #colnames(ref.allsnps)<-c("SNP") #needed for some strange error when SNP column is renamed here
      #setkeyv(ref.allsnps, cols = "SNP")
      #ref[ldscores[ref.allsnps, on='SNP'], on='SNP', c('L2','M') := list(i.L2,i.M)]
      ref[ldscores, on='SNP', c('L2','M') := list(i.L2,i.M)]
      cat("Done!\n")
    }
    
    #rename reference columns as to distinguish them from the dataset columns
    colnames(ref)<-paste0(names(ref),"_REF")
    setkeyv(ref, cols = paste0(ref.keys,"_REF"))
    #check keys with key(ref)
    
    
  } else {
    warning("\nRunning without reference.\n")
  }
  
  
  sumstats.meta<-data.table(name=traitNames,file_path=ifelse(is.null(filePaths),NA_character_,filePaths),n_snp_raw=NA_integer_,n_snp_res=NA_integer_)

  for(iFile in 1:nDatasets){
    #for testing!
    #iFile=3
    timeStart.ds <- Sys.time()
    
    #temporary variables which has to be reset for each file/dataset
    hasN<-F
    hasNEF<-F
    
    #set changeEffectDirectionOnAlleleFlip -this has to be reset for each file/dataset
    changeEffectDirectionOnAlleleFlip<-NULL
    if(!is.null(setChangeEffectDirectionOnAlleleFlip)){
      changeEffectDirectionOnAlleleFlip<-setChangeEffectDirectionOnAlleleFlip
    }
    
    #print per-dataset info and setting
    cat(paste("\n\nSupermunging\t",traitNames[iFile],"\n @ dataset", iFile,"\n"))
    cat("\nN=",N[iFile])
    cat("\nOLS=",OLS[iFile])
    cat("\nlinprob=",linprob[iFile])
    cat("\nse.logit=",se.logit[iFile])
    cat("\nprop=",prop[iFile])
    if(!is.null(list_df)){
      cSumstats <- list_df[[iFile]]
    } else {
      cFilePath<-filePaths[iFile]
      cat(paste("\nFile:", cFilePath,"\n"))
      cSumstats<-fread(file = cFilePath, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = F)
      #cSumstats <- read.table(cFilePath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
    }
    cat("\nReading.")
    
    if(test){
      cSumstats<-cSumstats[1:1000000,]
    }
    
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
    cSumstats.names <- stdGwasColumnNames(columnNames = colnames(cSumstats), stopOnMissingEssential = stopOnMissingEssential)
    cSumstats.names.string <-""
    for(iName in 1:nrow(cSumstats.names)){
      cSumstats.names.string<-paste(paste(cSumstats.names$orig,"\t->",cSumstats.names$std), collapse = '\n')
    }
    # print(cSumstats.names)
    # print(typeof(cSumstats.names$std))
    
    colnames(cSumstats) <- as.character(cSumstats.names$std)
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
    
    
    #Add in backstop SNP column
    if(!any(colnames(cSumstats)=="SNP")){
      cSumstats$SNP<-paste0("sm",1:nrow(cSumstats))
    }
    
    
    # Transform to data table
    #cSumstats <- setDT(cSumstats) #redundant now
    # Column harmonisation
    cSumstats.keys<-c('SNP')
    cSumstats[,SNP:=as.character(SNP)]
    if(any(colnames(cSumstats)=="A1")) {
      cSumstats[,A1:=as.character(A1)]
      cSumstats.keys<-c(cSumstats.keys,'A1')
    }
    if(any(colnames(cSumstats)=="A2")) {
      cSumstats[,A2:=as.character(A2)]
      cSumstats.keys<-c(cSumstats.keys,'A2')
    }

    if(any(colnames(cSumstats)=="Z")) cSumstats[,Z:=as.numeric(Z)]
    if(any(colnames(cSumstats)=="FRQ")) cSumstats[,FRQ:=as.numeric(FRQ)]
    #if(any(colnames(cSumstats)=="MAF")) cSumstats[,MAF:=as.numeric(MAF)]
    if(any(colnames(cSumstats)=="INFO")) cSumstats[,INFO:=as.numeric(INFO)]
    #if(any(colnames(cSumstats)=="EFFECT")) cSumstats[,EFFECT:=as.numeric(EFFECT)]
    if(any(colnames(cSumstats)=="SE")) cSumstats[,SE:=as.numeric(SE)]
    if(any(colnames(cSumstats)=="BETA")) cSumstats[,BETA:=as.numeric(BETA)]
    if(any(colnames(cSumstats)=="OR")) cSumstats[,OR:=as.numeric(OR)]
    if(any(colnames(cSumstats)=="N")) {
      cSumstats[,N:=as.numeric(N)]
      hasN<-T
      }
    if(any(colnames(cSumstats)=="NEF")) {
      cSumstats[,NEF:=as.numeric(NEF)]
      hasNEF<-T
      }
    cat(".")
    
    #parse SNP if needed
    if(parse){
      cSumstats[,SNP:=parseSNPColumnAsRSNumber(SNP)]
      cat(".")
    }
    
    if(any(colnames(cSumstats)=="CHR")) {
      if(parse) { 
        cSumstats[,CHR:=parseCHRColumn(CHR)] 
        } else {
        cSumstats[,CHR:=as.integer(CHR)]
        }
      cSumstats.keys<-c(cSumstats.keys,'CHR') 
    }
    cat(".")
    
    if('BP' %in% names(cSumstats)) {
      cSumstats[,BP:=as.integer(BP)]
      cSumstats.keys<-c(cSumstats.keys,'BP')
    }
    cat(".")
    
    if('BP2' %in% names(cSumstats)) {
      cSumstats[,BP2:=as.integer(BP2)]
      cSumstats.keys<-c(cSumstats.keys,'BP2')
    }
    cat(".")
    
    #set keys when columns stable
    setkeyv(cSumstats,cols = cSumstats.keys)
    cat(".")
    
    #Set effect to column standard - EFFECT
    #TODO Fix according to new column standard
    if(!("EFFECT" %in% colnames(cSumstats))){
      if("BETA" %in% colnames(cSumstats)) cSumstats[,EFFECT:=BETA] else
        if("OR" %in% colnames(cSumstats)) cSumstats[,EFFECT:=OR]
    }
    cat(".")
    
    
    #Filter variants MAF<maf.filter
    if(!is.null(maf.filter)){
      if("FRQ" %in% names(cSumstats)){
        rm <- (!is.na(cSumstats$FRQ) & ((cSumstats$FRQ<maf.filter & cSumstats$FRQ<0.5) | (1-cSumstats$FRQ)<maf.filter))
        cSumstats <- cSumstats[!rm, ]
        cSumstats.meta<-rbind(cSumstats.meta,list(paste("Removed variants; MAF <",maf.filter),as.character(sum(rm))))
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
    
    
    #lift-over to new coordinates before using coordinates
    if(!is.null(chainFilePath) & liftover[iFile] & any(colnames(cSumstats)=="CHR") & any(colnames(cSumstats)=="BP")){
      #chain file format reference: http://genome.ucsc.edu/goldenPath/help/chain.html
      
      #check if the file is the same build as the reference if present
      if(!is.null(ref)){
        cSumstatsBuildCheck<-cSumstats
        cSumstatsBuildCheck[ref, on=c(CHR='CHR_REF' , BP='BP_REF'), c('buildcheck') :=list(T)]
        if(nrow(cSumstatsBuildCheck[buildcheck==T,])>0.75*nrow(cSumstats)){
          cSumstats.warnings<-c(cSumstats.warnings,paste0("Dataset has more than a 75% overlap in genetic coordinates with the used reference variants (",nrow(cSumstatsBuildCheck[buildcheck==T,])," rows). Liftover may be unnecessary for this dataset."))
        }
        rm(cSumstatsBuildCheck)
      }
      cat(".")
      
      chain <- fread(file = chainFilePath, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = F, check.names = T, fill = T, blank.lines.skip = T, data.table = T,showProgress = F, nThread=nThreads)
      
      chain$row <- 1:nrow(chain)
      
      chains.dt <- chain[V1=="chain",]
      chains.dt[,V1:=NULL]
      colnames(chains.dt) <- c("score","tName","tSize","tStrand","tStart","tEnd","qName","qSize","qStrand","qStart","qEnd","id","row")
      
      
      #tName
      chains.dt$tName<-parseCHRColumn(chains.dt$tName)
      chains.dt$strange_tName<-grepl(pattern = "^[^_]+_", chains.dt$tName)
      indexesLengths<-regexec(pattern = "^([^_]+)_", text=chains.dt$tName)
      matches<-regmatches(chains.dt$tName,indexesLengths)
      chains.dt$tName[chains.dt$strange_tName] <- unlist(lapply(matches[chains.dt$strange_tName],FUN = function(x)x[2]))
      
      #qName
      chains.dt$qName<-parseCHRColumn(chains.dt$qName)
      chains.dt$strange_qName<-grepl(pattern = "^[^_]+_", chains.dt$qName)
      indexesLengths<-regexec(pattern = "^([^_]+)_", text=chains.dt$qName)
      matches<-regmatches(chains.dt$qName,indexesLengths)
      chains.dt$qName[chains.dt$strange_qName] <- unlist(lapply(matches[chains.dt$strange_qName],FUN = function(x)x[2]))
      
      #chromosome Un
      #https://genome.ucsc.edu/FAQ/FAQdownloads.html#download11
      
      chains.dt[,score:=as.numeric(score)]
      chains.dt[,tName:=as.integer(tName)]
      chains.dt[,qName:=as.integer(qName)]
      
      setkeyv(chains.dt,cols = c("tStart","tEnd","qStart","qEnd","id","row"))
      #chains.dt<-chains.dt[order("row")]
      
      segments.dt <- chain[V1!="chain",]
      segments.dt[,c("size","dt","dq"):=tstrsplit(V1,split = "\t",fixed = T)] #tstrsplit!
      segments.dt[,colnames(segments.dt)[!colnames(segments.dt) %in% c("row","size","dt","dq")]:=NULL]
      setkeyv(segments.dt,cols = c("row","size","dt","dq"))
      chains.dt<-chains.dt[order(chains.dt$row),]
      
      #update segments with chain id
      for(i in 1:nrow(chains.dt)){
        #i<-1
        cRow<-chains.dt[i,c("row")][[1]]
        cId<-chains.dt[i,c("id")][[1]]
        segments.dt[row>cRow,chain:=cId]
      }
      
      #segments.dt[, cumsize := cumsum(size), by=list(chain)]
      
      rm(chain) #we don't need chain anymore
      cat(".")
      
      #liftover
      
      #update cSumstats with chain id - sort on chain score!! =BLAT score?, -> overwrite lower score assignments later in the loop
      chains.dt<-chains.dt[order(chains.dt$score),]
      for(i in 1:nrow(chains.dt)){
        #i<-1
        cCHR<-chains.dt[i,c("tName")][[1]]
        cStartBP<-chains.dt[i,c("tStart")][[1]]
        cEndBP<-chains.dt[i,c("tEnd")][[1]]
        cNCHR<-chains.dt[i,c("qName")][[1]]
        cNStartBP<-chains.dt[i,c("qStart")][[1]]
        cId<-chains.dt[i,c("id")][[1]]
        cStrange_tName<-chains.dt[i,c("strange_tName")][[1]]
        cStrange_qName<-chains.dt[i,c("strange_qName")][[1]]
        # hits<-mTot[CHR==eval(cCHR) & BP>=eval(cStartBP) & BP<=eval(cEndBP),]
        # num<-nrow(hits)
        cSumstats[CHR==eval(cCHR) & BP>=eval(cStartBP) & BP<=eval(cEndBP),c('chain','strange_tName','strange_qName','NCHR','NBP') :=list(eval(cId),eval(cStrange_tName),eval(cStrange_qName),eval(cNCHR),(BP-eval(cStartBP)+eval(cNStartBP)))]
        if(any(colnames(cSumstats)=="BP2")){
          cSumstats[CHR==eval(cCHR) & BP2>=eval(cStartBP) & BP2<=eval(cEndBP),c('chain2','strange_tName2','strange_qName2','NCHR2','NBP2') :=list(eval(cId),eval(cStrange_tName),eval(cStrange_qName),eval(cNCHR),(BP2-eval(cStartBP)+eval(cNStartBP)))]
        }
        
      }
      
      # nrow(cSumstats)
      # 
      # remapped <- cSumstats[CHR!=NCHR | BP!= NBP,]
      # nrow(remapped)
      # 
      # unmapped <- cSumstats[is.na(NCHR) | is.na(NBP),]
      # nrow(unmapped)
      
      cSumstats<-cSumstats[!(is.na(NCHR) | is.na(NBP)),][,c("CHR_ORIG","BP_ORIG"):=list(CHR,BP)][,c("CHR","BP"):=list(NCHR,NBP)]
      #nrow(cSumstats[CHR!=NCHR | BP!= NBP,]) #check - should be 0
      
      if(any(colnames(cSumstats)=="BP2")){
        cSumstats<-cSumstats[!(is.na(NCHR2) | is.na(NBP2)),][,c("CHR2_ORIG","BP2_ORIG"):=list(CHR,BP2)][,BP2:=NBP2][,segmentChrMismatch:=(NCHR!=NCHR2)]
        #nrow(cSumstats[CHR!=NCHR2 | BP!= NBP2,]) #check - should be 0
        
      }
      
      cSumstats[,NCHR:=NULL][,NBP:=NULL]
      if(any(colnames(cSumstats)=="BP2")) cSumstats[,NCHR2:=NULL][,NBP2:=NULL]
      
      cat("🏋")
    }
    
    
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
    
    #remove custom regions according to specified dataframe
    if(!is.null(region.filter_df)){
      if(any(colnames(cSumstats)=="CHR") & any(colnames(cSumstats)=="BP")){
        if(!(any(colnames(region.filter_df)=="CHR") & any(colnames(region.filter_df)=="BP" & any(colnames(region.filter_df)=="BP2")))) stop("Dataframe containing regions to be removed must contain the columns CHR, BP, and BP2!")
        setDT(region.filter_df)
        setkeyv(region.filter_df, cols = c("CHR","BP","BP2"))
        cSumstats.nSNP<-nrow(cSumstats)
        for(isegment in 1:nrow(region.filter_df)){
          #isegment<-1
          cSumstats <- cSumstats[!(CHR==region.filter_df$CHR[isegment] & BP>=region.filter_df$BP[isegment] & BP<=region.filter_df$BP2[isegment]),]
        }
        cSumstats.meta<-rbind(cSumstats.meta,list(paste("Removed variants; custom regions"),as.character(cSumstats.nSNP-nrow(cSumstats))))
      } else {
        cSumstats.warnings<-c(cSumstats.warnings,"No chromosome or base-pair position information available - no filtering of custom provided regions was done!")
      }
    }
    cat(".")
    
    #update variant ID's from synonym list
    if(!is.na(idSynonyms)){
      # idSynonymsSplitMaxlength<-0
      # crap <- lapply(idSynonyms$parts,FUN = function(x){if(length(x)>idSynonymsSplitMaxlength) idSynonymsSplitMaxlength<<- length(x)})
      # rm(crap)
      
      for(iSynonym in 2:50){ #max 50 synonyms considered
        #iSynonym<-2
        cSynonym<-idSynonyms[parts.n>=iSynonym,]
        if(nrow(cSynonym)>0){
          cSynonym$SNP<-lapply(cSynonym$parts,FUN = function(x){ifelse(iSynonym<=length(x),x[[iSynonym]],NA)})
          cSynonym[,SNP:=as.character(SNP)][,SNP:=paste0("rs",SNP)][,c("SNP","first.rs")]
          setkeyv(cSynonym, cols = c("SNP"))
          #temp<-cSumstats[cSynonym,on=c(SNP="SNP"), nomatch=0]
          cSumstats[cSynonym,on=c(SNP="SNP"), SNP:=i.first.rs]
        } else break
      }
      cat("📛")
    }
    
    
    if(process){
      cat("\nProcessing.")
      # QC, and data management before merge with reference
      
      ## Remove SNPs with missing P
      if(any(colnames(cSumstats)=="P")) {
        cSumstats.n<-nrow(cSumstats)
        cSumstats<-cSumstats[!is.na(P),]
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; missing P",as.character(cSumstats.n-nrow(cSumstats))))
      }
      cat(".")
      
      ## Remove SNPs with missing effects
      if(any(colnames(cSumstats)=="EFFECT")) {
        cSumstats.n<-nrow(cSumstats)
        cSumstats<-cSumstats[!is.na(EFFECT),]
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; missing EFFECT",as.character(cSumstats.n-nrow(cSumstats))))
      }
      cat(".")
      
      ##Alleles, deal with indels
      if(keepIndel == T){
        #we already formatted these earlier
        # cSumstats$A1 <- as.character(toupper(cSumstats$A1))
        # cSumstats$A2 <- as.character(toupper(cSumstats$A2))
      } else if(keepIndel == F){
        cSumstats$A1 <- as.character(cSumstats$A1, c("A", "C", "G", "T"))
        cSumstats$A2 <- as.character(cSumstats$A2, c("A", "C", "G", "T"))
        cSumstats.meta<-rbind(cSumstats.meta,list("Discarded indels (A1)",as.character(count(is.na(cSumstats$A1)))))
        cSumstats.meta<-rbind(cSumstats.meta,list("Discarded indels (A2)",as.character(count(is.na(cSumstats$A2)))))
      }
      cat(".")
      
      ## Remove duplicated variants across SNP, A1 and A2 - DEPRECATED! - MOVED TO LATER
      # if(any(colnames(cSumstats)=="SNP") & any(colnames(cSumstats)=="A1") & any(colnames(cSumstats)=="A2")) {
      #   cSumstats.n <- nrow(cSumstats)
      #   cSumstats <- unique(cSumstats,by = c("SNP","A1","A2"))
      #   #TODO replace with a solution according to this syntax:
      #   #ul2<-l2[, .(L2 = head(L2,1)), by = c("CHR","SNP","BP")]
      #   cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; duplicate SNP,A1,A2",as.character(cSumstats.n-nrow(cSumstats))))
      # }
      # cat(".")
      
      # Merge with reference
      if(!is.null(ref)){
        #Aligning and validating with reference file
        cSumstats.n<-nrow(cSumstats)
        
        
        #Join with reference on SNP rsID, only keeping SNPs with rsIDs part of the reference
        #https://stackoverflow.com/questions/34644707/left-outer-join-with-data-table-with-different-names-for-key-variables/34645997#34645997
        cSumstats.merged.snp<-ref[cSumstats, on=c(SNP_REF="SNP"), nomatch=0]
        #replace missing columns
        cSumstats.merged.snp[,SNP:=SNP_REF]
        #cSumstats.merged.snp[,DBPORIG:=BP-BP_ORIG][,DBPREF:=BP-BP_REF]
        
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; rsID not in ref",as.character(cSumstats.n-nrow(cSumstats.merged.snp))))
        cat(".")
        
        if(any(colnames(cSumstats.merged.snp)=="CHR") && any(colnames(cSumstats.merged.snp)=="CHR_REF"))
        {
          cSumstats.merged.snp.n<-nrow(cSumstats.merged.snp)
          cSumstats.merged.snp<-cSumstats.merged.snp[CHR==CHR_REF]
          cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; CHR not matching ref",as.character(cSumstats.merged.snp.n-nrow(cSumstats.merged.snp))))
        }
        cat(".")
        
        #Join with reference on genetic coordinates
        cSumstats.merged.pos<-NULL
        if(any(colnames(cSumstats)=="CHR") && any(colnames(cSumstats)=="BP") && any(colnames(ref)=="CHR_REF") && any(colnames(ref)=="BP_REF") && any(colnames(cSumstats)=="A1") && any(colnames(ref)=="A1_REF") && any(colnames(cSumstats)=="A2") && any(colnames(ref)=="A2_REF")) {
          #Join with reference on position rather than rsID
          cSumstats.merged.pos<-ref[cSumstats, on=c(CHR_REF='CHR' , BP_REF='BP'), nomatch=0]
          
          #replace missing columns
          cSumstats.merged.pos[,CHR:=CHR_REF]
          cSumstats.merged.pos[,BP:=BP_REF]
        }
        cat(".")
        
        if(is.null(cSumstats.merged.pos)){
          cSumstats<-cSumstats.merged.snp
        } else {
          
          #assure proper allele configuration
          if(any(colnames(cSumstats.merged.pos)=="A1") && any(colnames(cSumstats.merged.pos)=="A1_REF") && any(colnames(cSumstats.merged.pos)=="A2") && any(colnames(cSumstats.merged.pos)=="A2_REF")) cSumstats.merged.pos<-cSumstats.merged.pos[A1==A1_REF & A2==A2_REF,]
          
          #merge merged datasets
          cSumstats.merged.pos.salvaged<-cSumstats.merged.pos[!(cSumstats.merged.pos$SNP_REF %in% cSumstats.merged.snp$SNP_REF),]
          cSumstats.merged.pos.salvaged[,c("SNP","inferredFromCoord"):=list(SNP_REF,T)] #make sure that these variants are interpreted as the one inferred from the reference from now on
          cSumstats<-rbindlist(list(cSumstats.merged.snp,cSumstats.merged.pos.salvaged), use.names=T, fill = T)
          #cat("\nSalvaged SNPs by merging on SNP position rather than rsID:",nrow(cSumstats.merged.pos.salvaged))
          cSumstats.meta<-rbind(cSumstats.meta,list("Salvaged SNPs by ref locus",as.character(nrow(cSumstats.merged.pos.salvaged))))
          cSumstats.merged.pos<-NULL
          cSumstats.merged.pos.salvaged<-NULL
          #TODO Make more memory friendly version using https://www.biostars.org/p/432389/
        }
        cSumstats.merged.snp<-NULL
        cat(".")
        
        #compare coordinates with reference
        if(any(colnames(cSumstats)=="BP") && any(colnames(cSumstats)=="BP_REF")){
          cSumstats.maxAlleleLength<-max(nchar(cSumstats$A1),nchar(cSumstats$A2),nchar(cSumstats$A1_REF),nchar(cSumstats$A2_REF))
          cSumstats.n<-nrow(cSumstats)
          cSumstats[,c("DBP","ADBP"):=list((BP-BP_REF),abs(BP-BP_REF))]
          #cSumstats.merged.snp[,DBPREF_ORIG:=abs(BP_ORIG-BP_REF)]
          
          cSumstats.medianDBPperCHR<-cSumstats[, .(medDBP = median(DBP)), by = c("CHR")]
          cat(".")
          
          #calculate adjusted BP to account for systematic differences between the dataset and reference
          cSumstats[cSumstats.medianDBPperCHR,on=(CHR="CHR"),c("BPADJ","ADBPADJ"):=list((BP-i.medDBP),abs((BP-i.medDBP)-BP_REF))]
          cat(".")
        }
        
      }
      
      # More QC and data management, after merge with reference
      
      #store original allele order and frequency info
      cSumstats[,A1_ORIG:=A1]
      cSumstats[,A2_ORIG:=A2]
      if(any(colnames(cSumstats)=="FRQ")) cSumstats[,FRQ_ORIG:=FRQ]
      
      
      if(!is.null(ref)){
        
        ##Synchronise SNP,BP with reference
        cSumstats[,SNP:=SNP_REF]
        if(harmoniseBPToReference) cSumstats[,BP:=BP_REF]
        
        ## Add in chr and bp from ref if not present in datasets
        if(!any(colnames(cSumstats)=="CHR")){
          cSumstats.warnings<-c(cSumstats.warnings,"No CHR column present!")
          if(any(colnames(cSumstats)=="CHR_REF")){
            cSumstats[,CHR:=CHR_REF]
            cSumstats.keys<-c(cSumstats.keys,'CHR')
            cSumstats.warnings<-c(cSumstats.warnings,"Inferring CHR from reference!")
          }
        }
        
        if(!any(colnames(cSumstats)=="FRQ")){
          sumstats.meta[iFile,c("no_FRQ")]<-T
          cSumstats.warnings<-c(cSumstats.warnings,"No FRQ column present!")
          if(any(colnames(cSumstats)=="MAF_REF")){
            cSumstats[,FRQ:=MAF_REF]
            cSumstats.warnings<-c(cSumstats.warnings,"Inferring FRQ from reference!")
          } else {
            cSumstats[,FRQ:=NA_real_]
          }
          
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
        cSumstats<-cSumstats[!cond.removeNonmatching, ]
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; A1 or A2 not matching any ref allele",as.character(sum(cond.removeNonmatching))))
        sumstats.meta[iFile,c("Removed, nonmatching ref alleles")]<-sum(cond.removeNonmatching)
      }
      cat(".")
      
      ## N
      if(any(colnames(cSumstats)=="N_CAS") && any(colnames(cSumstats)=="N_CON")) {
        ### Calculate total N from number of cases and number of controls if they are present. Overwrite any specific total N.
        cSumstats.meta<-rbind(cSumstats.meta,list("N",paste("N_CAS + N_CON")))
        cSumstats[,N:=N_CAS + N_CON]
      }
      cat(".")
      
      # if(any(colnames(cSumstats)=="N")){
      #   cSumstats.meta<-rbind(cSumstats.meta,list("N (median, min, max)",paste(median(cSumstats$N, na.rm = T),", ",min(cSumstats$N, na.rm = T),", ", max(cSumstats$N, na.rm = T))))
      # }
      
      if(!is.null(N) & length(N)>=iFile) {
        hasN<-T
        if(forceN && any(colnames(cSumstats)=="N")) {
          if(
            abs((median(cSumstats[is.finite(N),]$N, na.rm = T)-N[iFile])/median(cSumstats[is.finite(N),]$N, na.rm = T))>0.05
            |
            (N[iFile]-min(cSumstats[is.finite(N),]$N, na.rm = T))/N[iFile]>0.05
            |
            (N[iFile]-max(cSumstats[is.finite(N),]$N, na.rm = T))/N[iFile]>0.05
          ) cSumstats.warnings<-c(cSumstats.warnings,"Large (>5%) N discrepancies found between provided and existing N!")
          cSumstats$N<-N[iFile]
          cSumstats.meta<-rbind(cSumstats.meta,list("N",paste("Set to",N[iFile])))
        } else if(any(colnames(cSumstats)=="N")){
          cond<-is.na(cSumstats$N) | N[iFile] < cSumstats$N
          if(sum(cond)>0) {
            cSumstats$N<-ifelse(cond,N[iFile],cSumstats$N)
            cSumstats.meta<-rbind(cSumstats.meta,list("N",paste("Set to",N[iFile],"for",sum(cond)," NA or > specified.")))
          }
        } else cSumstats[,N:=eval(N)]
        
        cSumstats.meta<-rbind(cSumstats.meta,list("N (median, min, max)",paste(median(cSumstats[is.finite(N),]$N, na.rm = T),", ",min(cSumstats[is.finite(N),]$N, na.rm = T),", ", max(cSumstats[is.finite(N),]$N, na.rm = T))))
      } else if(!(any(colnames(cSumstats)=="N"))) {
        hasN<-F
        if(any(colnames(cSumstats)=="NEF")){
          cSumstats$N<-cSumstats$NEF
          cSumstats.meta<-rbind(cSumstats.meta,list("N","<= NEF"))
        } else {
          cSumstats.warnings<-c(cSumstats.warnings,"\nNo N column detected!")
          cSumstats.meta<-rbind(cSumstats.meta,list("N","Warning: Not detected!"))
          #cSumstats$N<-NA_integer_
        }
      }
      cat(".")
      
      ## Establish allele order from the reference
      if(!is.null(ref)){
        #cond.invertedAlleleOrder<-(cSumstats$A1 != cSumstats$A1_REF & cSumstats$A2 == cSumstats$A1_REF) #the same condition as in GenomicSEM munge.
        #cond.invertedAlleleOrder<-(cSumstats$A2 != cSumstats$A2_REF & cSumstats$A1 == cSumstats$A2_REF)
        #cond.invertedAlleleOrder<-((cSumstats$A2 != cSumstats$A2_REF & cSumstats$A1 == cSumstats$A2_REF) | (cSumstats$A1 != cSumstats$A1_REF & cSumstats$A2 == cSumstats$A1_REF)) #experimental - seems to work similar to the GenomicSEM implementation
        cSumstats[,cond.invertedAlleleOrder:=((A2!=A2_REF & A1==A2_REF) | (A1!=A1_REF & A2 ==A1_REF))]
        
      }
      cat(".")
      
      ## Invert alleles or harmonise (including their FRQ) to reference
      if(!is.null(ref) & harmoniseAllelesToReference){
        # Fix A1 and A2 to reflect the reference alleles
        cSumstats[,A1:=A1_REF]
        cSumstats[,A2:=A2_REF]
        if(any(colnames(cSumstats)=="MAF_REF")) cSumstats[,FRQ:=MAF_REF]
      } else if(any(colnames(cSumstats)=="cond.invertedAlleleOrder")) { ## Invert alleles  -FRQ is dealt with below
        cSumstats[,A1:=ifelse(cond.invertedAlleleOrder, A2_ORIG, A1)]
        cSumstats[,A2:=ifelse(cond.invertedAlleleOrder, A1_ORIG, A2)]
      }
      cat(".")
      
      # FRQ
      #cond.invertedFRQ<-NULL
      if(any(colnames(cSumstats)=="FRQ")) {
        ### Has FRQ
        #cSumstats.meta<-rbind(cSumstats.meta,list("FRQ (median, min(abs), max(abs))",paste(median(cSumstats$FRQ, na.rm = T),", ",min(abs(cSumstats$FRQ), na.rm = T),", ", max(abs(cSumstats$FRQ), na.rm = T))))
        #### Check if value is within limits [0,1]
        if(any(cSumstats$FRQ>1) || any(cSumstats$FRQ<0)) {
          stop(paste0('\nThere are FRQ values larger than 1 (',sum(cSumstats$FRQ>1),') or less than 0 (',sum(cSumstats$FRQ<0),') which is outside of the possible FRQ range.'))
        }
        
        ### Invert FRQ based on the previous reference matching
        if(any(colnames(cSumstats)=="cond.invertedAlleleOrder") & !harmoniseAllelesToReference) {
          alleleFRQ <- ifelse(cSumstats$cond.invertedAlleleOrder, (1-cSumstats$FRQ), cSumstats$FRQ)
          if(mean(alleleFRQ) < mean(cSumstats$FRQ) * 1.1){
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
      }
      cat(".")
      
      ## Compute variance of individual variant effects according to 2pq
      cSumstats[,VSNP:=2*FRQ*(1-FRQ)]
      cat(".")
      
      
      #remove duplicate ID variants, ordered by ADBPADJ, -MAF
      if(any(colnames(cSumstats)=="MAF")){
        cSumstats.n<-nrow(cSumstats)
        if(any(colnames(cSumstats)=="ADBPADJ")){
          cSumstats<-cSumstats[order(ADBPADJ, -MAF),]
        } else {
          cSumstats<-cSumstats[order(-MAF),]
        }
        cSumstats$ADBPADJ_ID<-1:nrow(cSumstats)
        cSumstats.unique<-cSumstats[, .(ADBPADJ_ID = head(ADBPADJ_ID,1)), by = c("SNP")]
        cSumstats<-cSumstats[cSumstats.unique, on=c(ADBPADJ_ID=c("ADBPADJ_ID"))]
        cSumstats[,i.SNP:=NULL]
        cSumstats.meta<-rbind(cSumstats.meta,list("Removed variants; duplicate variant ID",as.character(cSumstats.n-nrow(cSumstats))))
        cat(".")
      }
      
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
        cSumstats[,EFFECT:= Z/sqrt(N*VSNP)][,SE:=EFFECT/Z]
        cSumstats[is.na(SE),SE:=1] #explicitly set NA SE to 1
        cSumstats.meta<-rbind(cSumstats.meta,list("BETA","Calculated from Z"))
      }
      cat(".")
      
      if(any(colnames(cSumstats)=="EFFECT")) {
        
        ## Determine effect type, and set effect to log(EFFECT) if odds ratio
        if(round(median(cSumstats[is.finite(EFFECT),]$EFFECT,na.rm=T),digits = 1) == 1) {
          ###is odds ratio
          cSumstats[,EFFECT:=log(EFFECT)]
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
        if(any(colnames(cSumstats)=="Z")) cSumstats[,Z_ORIG:=Z] #save original Z-score
        if(any(colnames(cSumstats)=="P")) {
          
          cSumstats[,Z:=sign(EFFECT) * sqrt(qchisq(P,1,lower=F))]
          cSumstats.meta<-rbind(cSumstats.meta,list("Z","Calculated from P and sign(EFFECT)"))
        } else if(any(colnames(cSumstats)=="SE")){
          cSumstats[,Z:=EFFECT/SE] #is this less reliable as we cannot know the scale of SE?
          cSumstats.meta<-rbind(cSumstats.meta,list("Z","Calculated from EFFECT and SE"))
        } else {
          cSumstats.meta<-rbind(cSumstats.meta,list("Z","NOT calculated since no P or SE"))
        }
        cat(".")
        
        ## Inspect new and old Z-values
        if(any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="Z_ORIG")){
          if(mean(cSumstats[is.finite(Z),]$Z)-mean(cSumstats[is.finite(Z_ORIG),]$Z_ORIG,na.rm=T)>1) cSumstats.warnings<-c(cSumstats.warnings,"New Z differ from old by more than 1sd!")
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
          cSumstats[,EFFECT:=EFFECT*-1]
          cSumstats.meta<-rbind(cSumstats.meta,list("Inverted overall effect",as.character(length(cSumstats$EFFECT))))
        }
      }
      cat(".")
      
      
      #add missing SE
      if(!any(colnames(cSumstats)=="SE") & any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="EFFECT")){
        cSumstats[,SE:=EFFECT/Z]
      }
      
      #compute minimum variance for later calculations
      if(any(colnames(cSumstats)=="SE")){
        minv<-min(cSumstats[is.finite(SE),]$SE, na.rm = T)^2
      }
      
      #compare hypothesised inverted allele effects with non-inverted allele effects for validation
      if(is.null(changeEffectDirectionOnAlleleFlip)) changeEffectDirectionOnAlleleFlip<-T
      
      if(any(colnames(cSumstats)=="EFFECT") & any(colnames(cSumstats)=="cond.invertedAlleleOrder")){
        if(any(cSumstats$cond.invertedAlleleOrder)){
          sumstats.meta[iFile,c("Inverted allele order variants")]<-sum(cSumstats$cond.invertedAlleleOrder)
          
          cSumstats.meta<-rbind(cSumstats.meta,list("Mean effect","plain"))
          meffects.reference<-mean(cSumstats[is.finite(EFFECT) & !cond.invertedAlleleOrder,]$EFFECT,na.rm = T)
          meffects.candidate<-mean(cSumstats[is.finite(EFFECT) & cond.invertedAlleleOrder,]$EFFECT,na.rm = T)
          
          
          meffects.candidate.inverted<-meffects.candidate*-1
          sdeffects.reference<-sd(cSumstats[is.finite(EFFECT) & !cond.invertedAlleleOrder,]$EFFECT,na.rm = T)
          sdeffects.candidate<-sd(cSumstats[is.finite(EFFECT) & cond.invertedAlleleOrder,]$EFFECT,na.rm = T)
          cSumstats.meta<-rbind(cSumstats.meta,list("Number variants, reference, candidate:",paste0(as.character(length(cSumstats$EFFECT[!cSumstats$cond.invertedAlleleOrder])),",",as.character(length(cSumstats$EFFECT[cSumstats$cond.invertedAlleleOrder])))))
          cSumstats.meta<-rbind(cSumstats.meta,list("Mean reference effect (sd)",paste0(as.character(round(meffects.reference,digits = 5))," (",round(sdeffects.reference,digits = 5),")")))
          cSumstats.meta<-rbind(cSumstats.meta,list("Mean candidate effect (sd)",paste0(as.character(round(meffects.candidate,digits = 5))," (",round(sdeffects.candidate,digits = 5),")")))
          #cSumstats.meta<-rbind(cSumstats.meta,list("Mean candidate effect, inverted",as.character(round(abs(meffects.candidate.inverted),digits = 5))))
          #cSumstats.meta<-rbind(cSumstats.meta,list("Delta effect, plain",as.character(round(abs(meffects.reference-meffects.candidate),digits = 5))))
          #cSumstats.meta<-rbind(cSumstats.meta,list("Delta effect, inverted",as.character(round(abs(meffects.reference-meffects.candidate.inverted),digits = 5))))
          
          sumstats.meta[iFile,c("Reference variants")]<-length(cSumstats$EFFECT[!cSumstats$cond.invertedAlleleOrder])
          sumstats.meta[iFile,c("Candidate variants")]<-length(cSumstats$EFFECT[cSumstats$cond.invertedAlleleOrder])
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
      if(any(colnames(cSumstats)=="EFFECT") & any(colnames(cSumstats)=="cond.invertedAlleleOrder") & changeEffectDirectionOnAlleleFlip) {
        if(any(cSumstats$cond.invertedAlleleOrder)) cSumstats[,EFFECT:=ifelse(cond.invertedAlleleOrder,(EFFECT*-1),EFFECT)]
        
        meffects.new<-mean(cSumstats[is.finite(EFFECT),]$EFFECT,na.rm = T)
        
        cSumstats.meta<-rbind(cSumstats.meta,list("New effect mean",as.character(round(meffects.new,digits = 5))))
      }
      sumstats.meta[iFile,c("changeEffectDirectionOnAlleleFlip")]<-changeEffectDirectionOnAlleleFlip
      cat(".")
      
      
      if(any(colnames(cSumstats)=="cond.invertedAlleleOrder")) cSumstats.meta<-rbind(cSumstats.meta,list(
        paste0(
          "Modified SNPs; inverted allele order [",ifelse(any(colnames(cSumstats)=="FRQ"),"FRQ",""),",",
          ifelse(changeEffectDirectionOnAlleleFlip,"EFFECT",""),"]"),
        as.character(sum(cSumstats$cond.invertedAlleleOrder))))
      
      
      ## Compute Z score (standardised beta) and P or update it according to any corrections just done
      if(any(colnames(cSumstats)=="SE")) cSumstats[,Z:=EFFECT/SE]
      if(any(colnames(cSumstats)=="Z")) cSumstats[,P:=2*pnorm(q = abs(Z),mean = 0, sd = 1, lower.tail = F)]
      cat(".")
      
     
      cat("Processing done!")
    } #end of process stuff
    
    
    if(standardiseEffectsToExposure & any(colnames(cSumstats)=="EFFECT")) {
      cat("\nStandardising effect to exposure.")
      ##standardise outcome standardised BETA to exposure as well (fully standardised)
      ##correct binary trait regressions with a 'unit-variance-liability scaling'
      
      #quickfix backstop for missing N
      if(any(colnames(cSumstats)=="N")){
        if(any(is.na(cSumstats$N))){
          mN<-mean(cSumstats[is.finite(N),]$N,na.rm=T)
          cSumstats[is.na(get("N")), N:=mN]
        }
      }
      
      ## (Re-)Compute variance of individual variant effects according to 2pq, in case not already present (from processing step)
      if(!any(colnames(cSumstats)=="VSNP")) cSumstats[,VSNP:=2*FRQ*(1-FRQ)]
      
      if(OLS[iFile]){
        #Has OLS unstandardised beta
        cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","OLS"))
        if(any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="N")) {
          cSumstats[,EFFECT:=Z/sqrt(N*VSNP)] #standardisation
          cSumstats[,SE:=abs(EFFECT/Z)] #standardisation as this is derived form the standardised EFFECT
          cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Z, N, UVL std => BETA,SE"))
        } else stop("\nCould not compute BETA,SE because of missing Z or N!\n")
      
        
      } else if(linprob[iFile]){
        #Has effect based on a linear estimator for a binary outcome (rather than a logistic model)
        cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Binary, linear"))
        if(any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="N")){
          if(is.null(prop[iFile]) | is.na(prop[iFile])) stop("\nCould not perform correction of linear BETA,SE to liability scale because of missing or invalid prop argument!\n")
          cSumstats[,EFFECT:=Z/sqrt(prop[iFile]*(1-prop[iFile]) * N * VSNP)] #standardisation
          cSumstats[,SE:=1/sqrt(prop[iFile]*(1-prop[iFile]) * N * VSNP)] #standardisation
          cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Z, N,propCaCo, UVL std => BETA,SE"))
        } else stop("\nCould not compute BETA,SE because of missing Z or N!\n")
      
        cSumstats[,residualVarCorrectionTerm:=sqrt((EFFECT^2) + (pi^2)/3)]
        cSumstats[,EFFECT:=EFFECT/residualVarCorrectionTerm] #residual variance correction
        cSumstats[,SE:=cSumstats$SE/residualVarCorrectionTerm] #residual variance correction
       
      } else {
        #Has effect based on a logistic estimator for a binary outcome, OR or logistic beta
        cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Binary, logistic"))
        if(any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="N")) {
          cSumstats[,EFFECT:=Z/sqrt(N*VSNP)] #standardisation
          cSumstats[,SE:=abs(EFFECT/Z)] #standardisation as this is derived form the standardised EFFECT
          cSumstats.meta <- rbind(cSumstats.meta,list("EFFECT,SE","Z, N, UVL std=> BETA,SE"))
        } else stop("\nCould not compute BETA,SE because of missing Z or N!\n")
      
        cSumstats[,residualVarCorrectionTerm:=sqrt((EFFECT^2) + (pi^2)/3)]
        cSumstats[,EFFECT:=EFFECT/residualVarCorrectionTerm] #residual variance correction
        #cSumstats$SE <- cSumstats$SE/correctionTerm #residual variance correction
        if(se.logit[iFile]){
          cSumstats[,SE:=cSumstats$SE/residualVarCorrectionTerm] #residual variance correction
          cSumstats.meta <- rbind(cSumstats.meta,list("SE","logit"))
        } else {
          #transform to logit SE 
          cSumstats[,SE:=(SE/exp(EFFECT))/residualVarCorrectionTerm] #UV correction
          cSumstats.meta <- rbind(cSumstats.meta,list("SE","SE(OR) => logit"))
        }
        
      }
      cat(".")
      
      ## Compute Z,P again to update it according to any corrections just done
      if(any(colnames(cSumstats)=="SE")) {
        cSumstats[,Z:=EFFECT/SE]
        cSumstats[,P:=2*pnorm(q = abs(Z),mean = 0, sd = 1, lower.tail = F)]
      }
      cat("Done!\n")
      
    }
    
    
    #impute effects and standard errors; LD-IMP - highly experimental
    if(imputeFromLD){
      #impute betas using LD
      if(!(any(colnames(cSumstats)=="EFFECT") & any(colnames(cSumstats)=="SE") & any(colnames(cSumstats)=="CHR") & ((!is.null(imputeFrameLenBp) & any(colnames(cSumstats)=="BP")) | (!is.null(imputeFrameLenCM) & any(colnames(cSumstats)=="CM"))))) stop("LD imputation is not possible without the columns EFFECT,SE,CHR,(BP or CM)!")
      
      do.ldimp.cm<-!is.null(imputeFrameLenCM)
      
      frameLen<-ifelse(do.ldimp.cm,imputeFrameLenCM,imputeFrameLenBp)
      frameLenHalf<-frameLen/2
      cSumstats.merged.snp<-ref
      setkeyv(cSumstats,cols = cSumstats.keys)
      setkeyv(cSumstats.merged.snp, cols = paste0(ref.keys,"_REF"))
      
      #pick ancestry specific LD scores
      if(ancestrySetting[iFile]!="ANY" & !any(colnames(cSumstats.merged.snp)=="L2_REF")){
        if(any(colnames(cSumstats.merged.snp)==paste0("L2.",ancestrySetting[iFile],"_REF"))){
          sAncestryL2<-c(paste0("L2.",ancestrySetting[iFile],"_REF"))
          cSumstats.merged.snp$L2_REF<-cSumstats.merged.snp[,..sAncestryL2]
        }
      }
      
      if(!any(colnames(cSumstats.merged.snp)=="L2_REF") & any(grepl(pattern = "L2\\.+",x = colnames(cSumstats.merged.snp)))){
        sAncestryL2<-colnames(cSumstats.merged.snp)[grepl(pattern = "L2\\.+",x = colnames(cSumstats.merged.snp))]
        cSumstats.merged.snp$L2_REF<-cSumstats.merged.snp[,..sAncestryL2]
      }
      
      if(any(colnames(cSumstats)=="N")) cSumstats.merged.snp[cSumstats, on=c(SNP_REF='SNP'),c('BETA','SE','N') :=list(i.EFFECT,i.SE,i.N)] else cSumstats.merged.snp[cSumstats, on=c(SNP_REF='SNP'),c('BETA','SE') :=list(i.EFFECT,i.SE)]
      
      #filtering and selecting the subset to impute
      cSumstats.merged.snp.toimpute<-cSumstats.merged.snp[is.na(BETA) & MAF_REF>0.001,] #L2_REF>0
      cSumstats.merged.snp<-cSumstats.merged.snp[!is.na(BETA) & L2_REF>0 & MAF_REF>0.01,]
      
      #remove non-trustworthy variants according to specified regions in the df
      if(!is.null(region.imputation.filter_df)){
        if(!(any(colnames(region.imputation.filter_df)=="CHR") & any(colnames(region.imputation.filter_df)=="BP" & any(colnames(region.imputation.filter_df)=="BP2")))) stop("Dataframe containing regions to be excluded as imputation support must contain the columns CHR, BP, and BP2!")
        cSumstats.merged.snp.nSNP<-nrow(cSumstats.merged.snp)
        setDT(region.imputation.filter_df)
        setkeyv(region.imputation.filter_df, cols = c("CHR","BP","BP2"))
        for(isegment in 1:nrow(region.imputation.filter_df)){
          #isegment<-1
          cSumstats.merged.snp <- cSumstats.merged.snp[!(CHR_REF==region.imputation.filter_df$CHR[isegment] & BP_REF>=region.imputation.filter_df$BP[isegment] & BP_REF<=region.imputation.filter_df$BP2[isegment]),]
        }
        cSumstats.meta<-rbind(cSumstats.meta,list(paste("Ignored variants (imputation); custom regions"),as.character(cSumstats.merged.snp.nSNP-nrow(cSumstats.merged.snp))))
      }
      
      #truncate the imputation job when testing
      if(test){
        cSumstats.merged.snp.toimpute<-cSumstats.merged.snp.toimpute[1:10000,]
      }
      
      #cSumstats.merged.snp.toimpute<-cSumstats.merged.snp[!is.na(BETA),] #for validation
      setkeyv(cSumstats.merged.snp.toimpute,cols = "CHR_REF")
      chrsToImpute<-unique(cSumstats.merged.snp.toimpute$CHR_REF)
      cat(paste0("\nATTENTION!: Imputing ",nrow(cSumstats.merged.snp.toimpute)," variants!\n"))
      cat("I")
      #read intermediate results
      if(file.exists(file.path(pathDirOutput,paste0(traitNames[iFile],".LD-IMP.TEMP.Rds")))){
        intermediateResults<-readRDS(file=file.path(pathDirOutput,paste0(traitNames[iFile],".LD-IMP.TEMP.Rds")))
        cSumstats<-intermediateResults$cSumstats
        previousCHR<-intermediateResults$cCHR
        nChrIndex<-match(x = previousCHR, table = chrsToImpute)[1]+1
        if(nChrIndex<=length(chrsToImpute)) {
          chrsToImpute<-chrsToImpute[nChrIndex:length(chrsToImpute)]
        } else {
          chrsToImpute<-c()
        }
      }
      
      if(length(chrsToImpute)>0){
        for(cCHR in chrsToImpute){
          #cCHR<-10
          cSS<-cSumstats.merged.snp[CHR_REF==eval(cCHR),.(SNP=SNP_REF,BP=BP_REF,BETA,SE,L2=L2_REF,VAR=SE^2,CM=CM_REF)] #Z=EFFECT/SE
          
          if(do.ldimp.cm){
            setkeyv(cSS, cols = c("SNP","BP","CM")) #chromosome is fixed per chromosome loop
          } else {
            setkeyv(cSS, cols = c("SNP","BP")) #chromosome is fixed per chromosome loop
          }
          
          cI<-cSumstats.merged.snp.toimpute[CHR_REF==eval(cCHR),.(SNP=SNP_REF,BP=BP_REF,A1_REF,A2_REF,MAF_REF,L2=L2_REF,BETA,SE,VAR=SE^2,CM=CM_REF)]
          #cI<-cI[1:100,] #FOR TEST ONLY
          
          if(do.ldimp.cm){
            setkeyv(cI, cols = c("SNP","BP","CM"))
          } else {
            setkeyv(cI, cols = c("SNP","BP"))
          }
          
          if(nrow(cI)<1 || nrow(cSS)<1) next;
          #1 cM equates to ~1M bp
          for(i in 1L:nrow(cI)){
            #i<-1L
            cBP<-cI[i,BP]
            if(do.ldimp.cm){
              cCM<-cI[i,CM]
              frame<-cSS[BP!=cBP & CM< cCM+frameLenHalf & CM > cCM-frameLenHalf, .(BETA,SE,VAR,L2,W=L2/VAR)][,.(BETA,SE,VAR,L2,W,WBETA=W*BETA,WSE=W*SE)]
            } else {
              frame<-cSS[BP!=cBP & BP< cBP+frameLenHalf & BP > cBP-frameLenHalf, .(BETA,SE,VAR,L2,W=L2/VAR)][,.(BETA,SE,VAR,L2,W,WBETA=W*BETA,WSE=W*SE)]
            }
            k<-nrow(frame[!is.na(BETA),])
            W.sum<-sum(frame$W,na.rm = T)
            frame[,IMP.TERM:=WBETA/W.sum][,IMPSE.TERM:=WSE/W.sum]
            if(k>2){
              set(x = cI,i = i,j = "BETA.I",
                  value = sum(frame$IMP.TERM,na.rm = T)
              )
              set(x = cI,i = i,j = "SE.I",
                  value = sum(frame$IMPSE.TERM,na.rm = T)
              )
              set(x = cI,i = i,j = "K",
                  value = k
              )
              if(is.na(cI[i,L2]) | is.na(cI[i,L2])==0){
                set(x = cI,i = i,j = "LD",
                    value = median(frame$L2,na.rm = T)
                )
              }
              set(x = cI,i = i,j = "W.SUM",
                  value = W.sum
              )
              # set(x = cI,i = i,j = "INFO",
              #     value = k*cI[i,L2]/sqrt(W.sum)
              # )
            }
          } #for
          #validation
          # cI[,Z:=BETA/SE][,Z.I:=BETA.I/SE.I][,ZDIFF2:=(Z.I-Z)^2]
          # rmse<-sqrt(mean(cI$ZDIFF2,na.rm=T))
          # rmse
          # rmse2<-sqrt(median(cI$ZDIFF2,na.rm=T))
          # rmse2
          #View(cI)
          
          cI[,CHR:=eval(cCHR)][W.SUM>0,INFO:=k*L2/sqrt(W.SUM)]
          
          #add imputed variants
          if(any(colnames(cI)=="BETA.I") && any(colnames(cI)=="SE.I") && any(colnames(cI)=="K") && any(colnames(cI)=="INFO")){
            if(any(colnames(cSumstats)=="N")) {
              cI[,N:=round(mean(cSumstats.merged.snp[is.finite(N),]$N,na.rm=T))]
              cSumstats<-rbind(cSumstats,cI[,.(SNP,BP,CHR,A1=A1_REF,A2=A2_REF,FRQ=MAF_REF,N,EFFECT=BETA.I,SE=SE.I,LD_IMP.K=K,LD_IMP.SUM=W.SUM,LD_IMP=INFO)],fill=T)
            } else {
              cSumstats<-rbind(cSumstats,cI[,.(SNP,BP,CHR,A1=A1_REF,A2=A2_REF,FRQ=MAF_REF,EFFECT=BETA.I,SE=SE.I,LD_IMP.K=K,LD_IMP.SUM=W.SUM,LD_IMP=INFO)],fill=T)
            }
          }
          
          #write intermediate results
          #if(!(test & match(x = cCHR, table = chrsToImpute)[1]>length(chrsToImpute)/2)){
            saveRDS(object = list(cCHR=cCHR,cSumstats=cSumstats), file = file.path(pathDirOutput,paste0(traitNames[iFile],".LD-IMP.TEMP.Rds")))
          #}
          cat("I")
        }
      }
      
      ## Remove failed imputations
      cSumstats<-cSumstats[!is.na(EFFECT) && !is.na(SE),]
      
      ## Compute Z,P,VSNP again after imputation
      cSumstats$Z <- cSumstats$EFFECT/cSumstats$SE
      cSumstats$P <- 2*pnorm(q = abs(cSumstats$Z),mean = 0, sd = 1, lower.tail = F)
      setkeyv(cSumstats,cols = cSumstats.keys)
      cSumstats[,VSNP:=2*FRQ*(1-FRQ)]
      cat(".")
      
    }
    
    #adjust N to reflect uncertainty of imputed variants
    if((imputeFromLD | imputeAdjustN) & any(colnames(cSumstats)=="LD_IMP")){
      LD_IMP.median<-median(cSumstats[is.finite(LD_IMP),]$LD_IMP, na.rm = T)
      cSumstats[!is.na(LD_IMP),]$N <- round(N[iFile]*shru::clipValues(cSumstats[!is.na(LD_IMP),]$LD_IMP/LD_IMP.median,0,1))
      cSumstats.meta<-rbind(cSumstats.meta,list("N","Adjusted N to reflect uncertainty of imputed variants (LD-IMP)"))
    }
    
    #calculate genomic inflation factor
    medianChisq<-median(cSumstats[is.finite(Z),]$Z^2, na.rm = T)
    genomicInflationFactor<-medianChisq/qchisq(0.5,1)
    sumstats.meta[iFile,c("genomicInflationFactor")]<-genomicInflationFactor
    cSumstats.meta<-rbind(cSumstats.meta,list("Genomic inflation factor",as.character(round(genomicInflationFactor,digits = 4))))
    
    #basic re-inflation of deflated factor GWAS (typical for latent factor GWAS), using careful interpretation of the inflation (sqrt)
    if(GC=="reinflate" & genomicInflationFactor<0.9){
      meanChisq<-mean(cSumstats$Z^2) #using mean instead of median because it seems to be more stable than the median of deflated associations
      genomicInflationFactor<-meanChisq/qchisq(0.5,1)
      if(genomicInflationFactor<1){
        cSumstats[,Z:=Z/sqrt(genomicInflationFactor)]
        cSumstats[,EFFECT:=SE*Z] #attribute all of the re-inflation to the EFFECT
        #cSumstats[,SE:=EFFECT/Z] #attribute all of the re-inflation to the SE
        cSumstats[,P:=2*pnorm(q = abs(Z),mean = 0, sd = 1, lower.tail = F)]
        sumstats.meta[iFile,c("GC_reinflate")]<-T
        cSumstats.meta<-rbind(cSumstats.meta,list("GC","re-inflated"))
        medianChisq<-median(cSumstats$Z^2)
        genomicInflationFactor<-medianChisq/qchisq(0.5,1)
        cSumstats.meta<-rbind(cSumstats.meta,list("Genomic inflation factor 2",as.character(round(genomicInflationFactor,digits = 4))))
      }
    }
    
    #Calculate Effective Sample Size as advised from from the Genomic SEM Wiki
    ##citation: https://www.biorxiv.org/content/10.1101/603134v3
    if(any(colnames(cSumstats)=="EFFECT") & any(colnames(cSumstats)=="Z") & any(colnames(cSumstats)=="FRQ")){
      hasNEF <- any(colnames(cSumstats)=="NEF")
      cSumstats[,NEF:=round(((Z/EFFECT)^2)/VSNP,digits = 0)]
      cSumstats.meta <- rbind(
        cSumstats.meta,
        list("NEF (mean total, for MAF<.4, >.1 if available)",paste0(
          round(
            mean(
              ifelse(any(colnames(cSumstats)=="MAF"),
                     cSumstats[MAF<0.4&MAF>0.1&is.finite(NEF),NEF],
                     cSumstats[is.finite(NEF),]$NEF),
            na.rm=T),
          digits = 0))
          )
        )
      
      if(hasNEF & !hasN) cSumstats[,N:=NEF]
    }
    cat(".")
    
    ## Check effect value credibility
    if(any(colnames(cSumstats)=="EFFECT") & any(colnames(cSumstats)=="SE")) {
      if(any(colnames(cSumstats)=="MAF")) {
        # only for common + uncommon (non-rare) SNPs if MAF available
        if(mean(abs(cSumstats[MAF>0.001 & is.finite(EFFECT) & is.finite(SE),EFFECT/SE]), na.rm = T) > 5) cSumstats.warnings<-c(cSumstats.warnings,"\nNon-rare variant EFFECT/SE ratio >5 which could be a cause of misspecified/misinterpreted arguments!\n")
      } else {
        if(mean(abs(cSumstats[is.finite(EFFECT) & is.finite(SE),EFFECT/SE]), na.rm = T) > 5) cSumstats.warnings<-c(cSumstats.warnings,"\nOverall EFFECT/SE ratio >5 which could be a cause of misspecified/misinterpreted arguments!\n")
      }
    }
    cat(".")
    
    #rename the ambiguous EFFECT column to BETA, as it should be a regression beta at this point
    if(any(colnames(cSumstats)=="EFFECT")){
      cSumstats[,BETA:=EFFECT][,EFFECT:=NULL]
    }
    
    #NA values check
    if(any(is.na(cSumstats))) cSumstats.warnings<-c(cSumstats.warnings,"\nNA values detected among results!\n")
    
    # output columns
    output.colnames<- c("SNP")
    if(any(colnames(cSumstats)=="A1")) output.colnames<- c(output.colnames,"A1")
    if(any(colnames(cSumstats)=="A2")) output.colnames<- c(output.colnames,"A2")
    if("CHR" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"CHR")
    if("BP" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"BP")
    if("BP2" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"BP2")
    if("FRQ" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"FRQ")
    #if("MAF" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"MAF")
    if("BETA" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"BETA")
    if("SE" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"SE")
    if("Z" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"Z")
    if("P" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"P")
    if("N" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N")
    if("N_CAS" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N_CAS")
    if("N_CON" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N_CON")
    if("NEF" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"NEF")
    if("DF" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"DF")
    if("LD_IMP.K" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"LD_IMP.K")
    if("LD_IMP.SUM" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"LD_IMP.SUM")
    if("LD_IMP" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"LD_IMP")
    
    output.colnames.more<-colnames(cSumstats)[!(colnames(cSumstats) %in% output.colnames)]
    output.colnames.all<-c(output.colnames,output.colnames.more)
    if(lossless){
      cSumstats<-cSumstats[,..output.colnames.all]
    } else {
      cSumstats<-cSumstats[,..output.colnames]
    }
    
    cSumstats.meta<-rbind(cSumstats.meta,list("Variants after supermunge",as.character(nrow(cSumstats))))
    
    #merge with variantTable
    if(produceCompositeTable){
      cat("\nProducing composite variant table.\n")
      cNames.toJoin<-c("SNP","EFFECT","SE")
      if(any(colnames(cSumstats)=="FRQ")) cNames.toJoin <- c(cNames.toJoin,"FRQ")
      if(any(colnames(cSumstats)=="LD_IMP")) cNames.toJoin <- c(cNames.toJoin,"LD_IMP")
      if(any(colnames(cSumstats)=="LD_IMP.K")) cNames.toJoin <- c(cNames.toJoin,"LD_IMP.K")
      toJoin <- cSumstats[,..cNames.toJoin]
      setkeyv(toJoin,cols = 'SNP')
      #update variantTable by reference
      #ref https://stackoverflow.com/questions/44433451/r-data-table-update-join
      #ref about data.table improvements https://stackoverflow.com/questions/42537520/data-table-replace-data-using-values-from-another-data-table-conditionally/42539526#42539526
      cName.beta <- paste0("BETA.",traitNames[iFile])
      cName.se <- paste0("SE.",traitNames[iFile])
      cName.frq <- paste0("FRQ.",traitNames[iFile])
      cName.infolimp <- paste0("LD_IMP",traitNames[iFile])
      cName.k <- paste0("K.",traitNames[iFile])
      if(!any(colnames(toJoin)=="FRQ")) toJoin[,FRQ=NA_real_]
      if(!any(colnames(toJoin)=="LD_IMP")) toJoin[,LD_IMP=NA_real_]
      if(!any(colnames(toJoin)=="LD_IMP.K")) toJoin[,LD_IMP.K=NA_integer_]
      variantTable[toJoin, on='SNP', c(cName.beta,cName.se,cName.frq,cName.LD_IMP,cName.LD_IMP.K):=.(i.EFFECT,i.SE,i.FRQ,i.LD_IMP,i.LD_IMP.K)]
    }
    
    
    
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
    if(writeOutput){
        cat("\nSaving supermunged dataset...\n\n")
        if((!any(colnames(cSumstats)=="CHR") & doChrSplit)) warning("\nSplit by chromosome specified, but the dataset does not have a CHR column.")
        if(!doChrSplit | ((!any(colnames(cSumstats)=="CHR") & doChrSplit))){
          #write.table(x = cSumstats,file = nfilepath,sep="\t", quote = FALSE, row.names = F, append = F)
          #nfilepath.gzip<-gzip(nfilepath)
          fwrite(x = cSumstats,file = paste0(nfilepath,".gz"),append = F,quote = F,sep = "\t",col.names = T)
          cat(paste("\nSupermunged dataset saved as", paste0(nfilepath,".gz")))
        } else {
          dir.create(paste0(nfilepath,".chr"), showWarnings = FALSE)
          chromosomes<-unique(cSumstats$CHR)
          for(cChr in chromosomes){
            fwrite(x = cSumstats[CHR==cChr,],file = file.path(paste0(nfilepath,".chr"),paste0(traitNames[iFile],"_",cChr,".gz")),append = F,quote = F,sep = "\t",col.names = T)
          }
          cat(paste("\nSupermunged dataset saved as one file per chromosome under", paste0(nfilepath,".chr")))
      }
    }
    
    #remove intermediate results
    if(file.exists(file.path(pathDirOutput,paste0(traitNames[iFile],".LD-IMP.TEMP.Rds")))) file.remove(file.path(pathDirOutput,paste0(traitNames[iFile],".LD-IMP.TEMP.Rds")))
    
    timeStop.ds <- Sys.time()
    timeDiff <- difftime(time1=timeStop.ds,time2=timeStart.ds,units="sec")
    timeDiff.minutes <- floor(floor(timeDiff)/60)
    timeDiff.seconds <- timeDiff-timeDiff.minutes*60
    
    cat("\nSupermunge of ",traitNames[iFile]," was done in",timeDiff.minutes, "minutes and",timeDiff.seconds," seconds.\n")
    gc() #do garbage collect if this can help with out of memory issues.
  }
  
  #process the variant table further
  # if(!is.null(variantTable)){
  #   # sum of K
  #   colK<-colnames(variantTable)[grep("^K\\.", ignore.case = TRUE,colnames(variantTable))]
  #   if(length(colK)>0) variantTable$K.SUM<-rowSums(abs(variantTable[,..colK]),na.rm = T)
  #   colINFO.LIMP<-colnames(variantTable)[grep("^INFO\\.LIMP\\.", ignore.case = TRUE,colnames(variantTable))]
  #   if(length(colINFO.LIMP)>0) variantTable$INFO.LIMP.SUM<-rowSums(abs(variantTable[,..colK]),na.rm = T)
  # }
  
  timeStop <- Sys.time()
  timeDiff <- difftime(time1=timeStop,time2=timeStart,units="sec")
  timeDiff.minutes <- floor(floor(timeDiff)/60)
  timeDiff.seconds <- timeDiff-timeDiff.minutes*60
  
  cat("\nSupermunge of all datasets was done in",timeDiff.minutes, "minutes and",timeDiff.seconds,"seconds.")
  
  return(list(
    meta=as.data.frame(sumstats.meta),
    last=as.data.frame(cSumstats),
    composite=as.data.frame(variantTable),
    ref=as.data.frame(ref)
  )
    )
  
}
