#WORK IN PROGRESS
#Based on the amazing work by Grotzinger, A. D. et al. Nat. Hum. Behav. 3, 513–525 (2019) and Bulik-Sullivan, B. K. et al. Nat. Genet. 47, 291–295 (2015).

stdGwasColumnNames <- function(columnNames, stopOnMissingEssential=T, warnOnMultiple=T,
     c.SNP = c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS","MARKERNAME","ID","PREDICTOR","RS"),
     c.A1 = c("A1", "ALLELE1","ALLELE_1","EFFECT_ALLELE","INC_ALLELE","EA"),
     c.A2 = c("A2","ALLELE2","ALLELE_2","ALLELE0","OTHER_ALLELE","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA","ALT","REFERENCE_ALLELE","REF"),
     c.EFFECT = c("EFFECT","OR","B","BETA","LOG_ODDS","EFFECTS","SIGNED_SUMSTAT","Z","ZSCORE","Z-SCORE","EST","ZSTAT","ZSTATISTIC","GC_ZSCORE"),
     c.SE = c("SE","STDER"),
     c.Z = c("Z","ZSCORE","Z-SCORE","ZSTAT","ZSTATISTIC"),
     c.INFO = c("INFO","IMPINFO"),
     c.P = c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P"),
     c.N = c("N","WEIGHT","NCOMPLETESAMPLES","TOTALSAMPLESIZE","TOTALN","TOTAL_N","N_COMPLETE_SAMPLES"),
     c.N_CAS = c("N_CAS","NCASE","N_CASE","N_CASES","NCAS","NCA"),
     c.N_CON = c("N_CON","NCONTROL","N_CONTROL","N_CONTROLS","N_CON","CONTROLS_N","NCON","NCO"),
     c.FRQ = c("FRQ","MAF","AF","CEUAF","FREQ","FREQ1","EAF", "FREQ1.HAPMAP", "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU", "EFFECT_ALLELE_FREQ","FREQ.A1","FRQ_U","F_U"),
     c.CHR = c("CHR", "CH", "CHROMOSOME", "CHROM"),
     c.BP = c("BP", "ORIGBP", "POS")
                                       ){
  #test
  #columnNames<-cSumstats.names
  
  columnNames.upper<-toupper(columnNames)
  #names(columnNames)<-columnNames
  columnNames.orig<-columnNames
  
  columnNames[columnNames.upper %in% c.SNP] <- c.SNP[1]
  columnNames[columnNames.upper %in% c.A1] <- c.A1[1]
  columnNames[columnNames.upper %in% c.A2] <- c.A2[1]
  columnNames[columnNames.upper %in% c.EFFECT] <- c.EFFECT[1]
  columnNames[columnNames.upper %in% c.SE] <- c.SE[1]
  columnNames[columnNames.upper %in% c.Z] <- c.Z[1]
  columnNames[columnNames.upper %in% c.INFO] <- c.INFO[1]
  columnNames[columnNames.upper %in% c.P] <- c.P[1]
  columnNames[columnNames.upper %in% c.N] <- c.N[1]
  columnNames[columnNames.upper %in% c.N_CAS] <- c.N_CAS[1]
  columnNames[columnNames.upper %in% c.N_CON] <- c.N_CON[1]
  columnNames[columnNames.upper %in% c.FRQ] <- c.FRQ[1]
  columnNames[columnNames.upper %in% c.CHR] <- c.CHR[1]
  columnNames[columnNames.upper %in% c.BP] <- c.BP[1]
  
  if(stopOnMissingEssential){
    # Stop if any of these columns are not found
    if(!any(columnNames=="SNP")) stop("Could not find the 'SNP' column.")
    if(!any(columnNames=="P")) warning("Could not find the P-value column. Standard is 'P'.")
    if(!any(columnNames=="A1")) stop("Could not find the 'A1' column.")
    if(!any(columnNames=="A2")) stop("Could not find the 'A2' column.")
    if(!any(columnNames=="EFFECT")) warning("Could not find the 'EFFECT' column.")
    if(!any(columnNames=="FRQ")) warning("Could not find the 'FRQ' column.")
  }
  
  if(warnOnMultiple){
    # Warn if multiple of these columns are found
    if(sum(columnNames=="SNP")>1) warning("Multiple 'SNP' columns found!")
    if(sum(columnNames=="P")>1) warning("Multiple 'P' columns found!")
    if(sum(columnNames=="A1")>1) warning("Multiple 'A1' columns found!")
    if(sum(columnNames=="A2")>1) warning("Multiple 'A2' columns found!")
    if(sum(columnNames=="EFFECT")>1) warning("Multiple 'EFFECT' columns found!")
    if(sum(columnNames=="FRQ")>1) warning("Multiple 'FRQ' columns found!")
  }
  
  return(data.frame(std=columnNames,orig=columnNames.orig))
}


parseSNPColumnAsRSNumber <- function(text){
  #decide if BGENIE SNP format using top 100,000 SNPs
  #TODO this condition may be improved to not rely on the number of variants being >100,000
  #test
  #text<-files[[i]]$SNP
  if(sum(grepl(pattern = "^\\d+:\\w+_\\w+_\\w+", x= head(x = text, n=100000)))>90000){
    #extract and format rs-no
    indexesLengths<-regexec(pattern = "^\\d+:(\\w+)_\\w+_\\w+", text=text)
    matches<-regmatches(text,indexesLengths)
    return(lapply(X = matches, FUN = function(x)paste0("RS",x[2])))
  }
  
  return(text)
}

# parseCHRColumn <- function(text){
#   text<-cSumstats$CHR
#   #decide if formatting needed
#   if(sum(grepl(pattern = "0\\w|[\\W]", x= head(x = text, n=1000)))>995){
#     noleadingzero<-gsub(pattern = "0(\\w)", replacement = "\\2_" text=text)
#     return(noleadingzero)
#   }
#   
#   return(text)
# }

#test
# filePaths = project$munge$filesToUse
# refFilePath = project$filepath.SNPReference.1kg
# traitNames = project$munge$traitNamesToUse
# N = project$munge$NToUse
# OLS=NULL
# linprob=NULL
# maxSNPDistanceBpPadding=0
# pathDirOutput = project$folderpath.data.sumstats.munged
# mask<-c(F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)

supermunge <- function(filePaths,
                       refFilePath=NULL,
                       traitNames=NULL,
                       setChangeEffectDirectionOnAlleleFlip=T, #set to TRUE to emulate genomic sem munge
                       N=NULL,
                       OLS=NULL,
                       linprob=NULL,
                       pathDirOutput=".",
                       keepIndel=T,
                       harmoniseAllelesToReference=F,
                       doChrSplit=F,
                       doStatistics=F,
                       mask=NULL,
                       stopOnMissingEssential=T,
                       maxSNPDistanceBpPadding=0,
                       invertEffectDirectionOn=NULL
                       ){
  
  timeStart <- Sys.time()
  
  if(is.null(traitNames)){
    traitNames<-basename(filePaths)
  }
  
  if(!is.null(mask)){
    filePaths<-filePaths[mask]
    if(!is.null(traitNames)) traitNames<-traitNames[mask]
  }
  
  filePaths.length <- length(filePaths)
  
  if(is.null(OLS)){
    OLS<-rep(FALSE,filePaths.length)
  }
  
  if(is.null(linprob)){
    linprob<-rep(FALSE,filePaths.length)
  }
  
  if(is.null(traitNames)){
    names.beta <- paste0("beta.",1:filePaths.length)
    names.se <- paste0("se.",1:filePaths.length)
  }else{
    names.beta <- paste0("beta.",traitNames)
    names.se <- paste0("se.",traitNames)
  }
  
  cat("\n\n\nS U P E R ★ M U N G E\n")
  
  cat("\n--------------------------------\nSettings:")
  cat("\nkeepIndel=",keepIndel)
  cat("\nharmoniseAllelesToReference=",harmoniseAllelesToReference)
  cat("\nmaxSNPDistanceBpPadding=",maxSNPDistanceBpPadding)
  cat("\nchangeEffectDirectionOnAlleleFlip=",setChangeEffectDirectionOnAlleleFlip)
  cat("\n--------------------------------\n")
  
  ref<-NULL
  if(!is.null(refFilePath)){
    cat(paste0("\nReading reference file...\n"))
    ref <- read.table(refFilePath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))

    #rename reference columns as to distinguish them from the dataset columns
    names(ref)<-paste0(names(ref),"_REF")
    # Transform to data table
    ref <- setDT(ref)
    # Column harmonisation
    ref.keys<-c('SNP_REF')
    ref$SNP_REF <- as.character(toupper(ref$SNP_REF))
    ref$A1_REF <- as.character(toupper(ref$A1_REF))
    ref$A2_REF <- as.character(toupper(ref$A2_REF))
    if('CHR_REF' %in% names(ref)) {
      ref$CHR_REF <- as.character(toupper(ref$CHR_REF))
      ref.keys<-c(ref.keys,'CHR_REF') 
    }
    if('BP_REF' %in% names(ref)) {
      ref$BP_REF <- as.integer(ref$BP_REF)
      ref.keys<-c(ref.keys,'BP_REF')
      }
    setkeyv(ref,cols = ref.keys)
    cat(paste0("\nRead reference file:\n",refFilePath))
  } else {
    warning("\nRunning without reference.\n")
  }
  
  sumstats.meta<-data.table(name=traitNames,file_path=filePaths,n_snp_raw=NA_integer_,n_snp_res=NA_integer_)
  
  sumstats<-c()
  for(iFile in 1:filePaths.length){
    
    #set changeEffectDirectionOnAlleleFlip -this has to be reset for each file
    changeEffectDirectionOnAlleleFlip<-NULL
    if(!is.null(setChangeEffectDirectionOnAlleleFlip)){
      changeEffectDirectionOnAlleleFlip<-setChangeEffectDirectionOnAlleleFlip
    }
    
    #for testing!
    #iFile=1
    timeStart.ds <- Sys.time()
    cFilePath<-filePaths[iFile]
    cat(paste("\n\nSupermunging\t",traitNames[iFile],"\nFile:", cFilePath,"\n"))
    cat("\nProcessing...")
    cSumstats.meta<-data.table(message=NA_character_,display=NA_character_)
    cSumstats.warnings<-list()
    cSumstats <- read.table(cFilePath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
    cSumstats.nSNP.raw<-nrow(cSumstats)
    sumstats.meta[iFile,c("n_snp_raw")]<-cSumstats.nSNP.raw
    cSumstats.meta<-rbind(cSumstats.meta,list("Input SNPs",as.character(cSumstats.nSNP.raw)))
    if(!is.null(ref)) {
      cSumstats.meta<-rbind(cSumstats.meta,list("Reference SNPs",as.character(nrow(ref))))
    }
    cat("...")
    
    # select supporting info for the current dataset
    cOLS<-OLS[iFile]
    
    # Give sumstats new standardised column names
    
    cSumstats.names <- stdGwasColumnNames(columnNames = names(cSumstats), stopOnMissingEssential = stopOnMissingEssential)
    cSumstats.names.string <-""
    #apply(cSumstats.names, MARGIN = 1, FUN = function(c){cat(c[2],"\t-> ",c[1],"\n")})
    for(iName in 1:nrow(cSumstats.names)){
      cSumstats.names.string<-paste(paste(cSumstats.names$orig,"\t->",cSumstats.names$std), collapse = '\n')
    }
    names(cSumstats) <- cSumstats.names$std
    cat(".")
    
    # Transform to data table
    cSumstats <- setDT(cSumstats)
    # Column harmonisation
    cSumstats.keys<-c('SNP')
    cSumstats$SNP <- as.character(toupper(cSumstats$SNP))
    cSumstats$A1 <- as.character(toupper(cSumstats$A1))
    cSumstats$A2 <- as.character(toupper(cSumstats$A2))
    
    #parse SNP if needed
    cSumstats$SNP<-parseSNPColumnAsRSNumber(cSumstats$SNP)
    
    if('CHR' %in% names(cSumstats)) {
      cSumstats$CHR <- as.character(toupper(cSumstats$CHR))
      cSumstats.keys<-c(cSumstats.keys,'CHR') 
    }
    if('BP' %in% names(cSumstats)) {
      cSumstats$BP <- as.integer(cSumstats$BP)
      cSumstats.keys<-c(cSumstats.keys,'BP')
    }
    setkeyv(cSumstats,cols = cSumstats.keys)
    cat("..")
    
    
    # QC, and data management before merge with reference
    
    cSumstats.meta<-rbind(cSumstats.meta,list("SNPs before supermunge",as.character(nrow(cSumstats))))
    
    ## Add effect for datasets containing only Z
    if(!("EFFECT" %in% colnames(cSumstats)) && "Z" %in% colnames(cSumstats)){
      cSumstats$EFFECT<-cSumstats$Z
    }
    cat(".")
    
    ## Save original Z
    if("Z" %in% colnames(cSumstats)){
      cSumstats$Z_ORIG<-cSumstats$Z
    }
    cat(".")
    
    ## Remove SNPs with missing P
    if("P" %in% colnames(cSumstats)) {
      cSumstats.n<-nrow(cSumstats)
      cSumstats<-cSumstats[which(!is.na(cSumstats$P)),]
      cSumstats.meta<-rbind(cSumstats.meta,list("Removed SNPs; missing P",as.character(cSumstats.n-nrow(cSumstats))))
    }
    cat(".")
    
    ## Remove SNPs with missing effects
    if("EFFECT" %in% colnames(cSumstats)) {
      cSumstats.n<-nrow(cSumstats)
      cSumstats<-cSumstats[which(!is.na(cSumstats$EFFECT)),]
      cSumstats.meta<-rbind(list("Removed SNPs; missing EFFECT",as.character(cSumstats.n-nrow(cSumstats))))
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
    cat("..")
    
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
    
    if(!is.null(ref)){
      ## Add in chr and bp from ref if not present in datasets
      if(!any(colnames(cSumstats)=="CHR") & any(colnames(cSumstats)=="CHR_REF"))cSumstats$CHR<-cSumstats$CHR_REF
      if(!any(colnames(cSumstats)=="BP") & any(colnames(cSumstats)=="BP_REF"))cSumstats$BP<-cSumstats$BP_REF
    }
    
    if(!is.null(ref)){
      ## Remove SNPs where alleles are not matching at least one of the reference alleles
      cSumstats.n<-nrow(cSumstats)
      cond.removeNonmatching<-(cSumstats$A1 != (cSumstats$A1_REF) & cSumstats$A1 != (cSumstats$A2_REF)) & (cSumstats$A2 != (cSumstats$A1_REF)  & cSumstats$A2 != (cSumstats$A2_REF))
      cSumstats<-cSumstats[which(!cond.removeNonmatching), ]
      cSumstats.meta<-rbind(cSumstats.meta,list("Removed SNPs; A1 or A2 not matching any ref allele",as.character(sum(cond.removeNonmatching))))
      sumstats.meta[iFile,c("Removed, nonmatching ref alleles")]<-sum(cond.removeNonmatching)
    }
    cat(".")
    
    ## Establish allele order from the reference
    #store original allele order and frequency info
    cSumstats$A1_ORIG<-cSumstats$A1
    cSumstats$A2_ORIG<-cSumstats$A2
    if(any(colnames(cSumstats)=="FRQ")) cSumstats$FRQ_ORIG<-cSumstats$FRQ
    
    cond.invertedAlleleOrder<-NULL
    if(!is.null(ref)){
      #cond.invertedAlleleOrder<-(cSumstats$A1 != cSumstats$A1_REF & cSumstats$A2 == cSumstats$A1_REF) #the same condition as in GenomicSEM munge.
      #cond.invertedAlleleOrder<-(cSumstats$A2 != cSumstats$A2_REF & cSumstats$A1 == cSumstats$A2_REF)
      cond.invertedAlleleOrder<-((cSumstats$A2 != cSumstats$A2_REF & cSumstats$A1 == cSumstats$A2_REF) | (cSumstats$A1 != cSumstats$A1_REF & cSumstats$A2 == cSumstats$A1_REF)) #experimental - seems not to work
    } else if(any(colnames(cSumstats)=="FRQ")){
      cond.invertedAlleleOrder<- cSumstats$FRQ > .55
    }
    cat(".")
    
    #compare hypothesised inverted allele effects with non-inverted allele effects for validation
    if(!is.null(cond.invertedAlleleOrder)){
      if(any(cond.invertedAlleleOrder)){
        sumstats.meta[iFile,c("Inverted allele order variants")]<-sum(cond.invertedAlleleOrder)
        if(any(colnames(cSumstats)=="SE")){
          cSumstats.meta<-rbind(cSumstats.meta,list("Mean effect","ivw"))
          minv<-min(cSumstats$SE)^2
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
    
    ## Invert alleles
    if(!is.null(cond.invertedAlleleOrder)) cSumstats$A1<-ifelse(cond.invertedAlleleOrder, cSumstats$A2_ORIG, cSumstats$A1)
    if(!is.null(cond.invertedAlleleOrder)) cSumstats$A2<-ifelse(cond.invertedAlleleOrder, cSumstats$A1_ORIG, cSumstats$A2)
    cat(".")
    
    ## EFFECT
    if(is.null(changeEffectDirectionOnAlleleFlip)) changeEffectDirectionOnAlleleFlip<-T
    if("EFFECT" %in% colnames(cSumstats)) {
      
      ## Determine effect type, and set effect to log(EFFECT) if odds ratio
      if(round(median(cSumstats$EFFECT,na.rm=T)) == 1) {
        #is odds ratio
        cSumstats$EFFECT<-log(cSumstats$EFFECT)
        sumstats.meta[iFile,c("effect_type")]<-"OR"
        cSumstats.meta<-rbind(cSumstats.meta,list("EFFECT","OR  =>ln(OR)"))
        if(!any(c("OR") %in% toupper(cSumstats.names$orig)) | any(c("B","BETA","LOG_ODDS","Z","ZSCORE","EST","ZSTAT","ZSTATISTIC") %in% toupper(cSumstats.names$orig))) {
          cSumstats.warnings<-c(cSumstats.warnings,"\nThe effect format being an ODDS RATIO is not compatible with the original variable naming scheme!")
          sumstats.meta[iFile,c("effect_type_warning")]<-T
        }
      } else {
        #is NOT odds ratio
        sumstats.meta[iFile,c("effect_type")]<-"non_OR"
        cSumstats.meta<-rbind(cSumstats.meta,list("EFFECT","NON OR"))
        if(any(c("OR") %in% toupper(cSumstats.names$orig)) | !any(c("B","BETA","LOG_ODDS","Z","ZSCORE","EST","ZSTAT","ZSTATISTIC") %in% toupper(cSumstats.names$orig))) {
          cSumstats.warnings<-c(cSumstats.warnings,"\nThe effect format is not compatible with the original variable naming scheme!")
          sumstats.meta[iFile,c("effect_type_warning")]<-T
        }
      }
      
      ##invert overall effect if specified
      if(!is.null(invertEffectDirectionOn)){
        if(any(invertEffectDirectionOn==traitNames[iFile])){
          cSumstats$EFFECT<-cSumstats$EFFECT*-1
          cSumstats.meta<-rbind(cSumstats.meta,list("Inverted overall effect",as.character(length(cSumstats$EFFECT))))
        }
      }
      
      ##invert effect if inverted allele order
      if(!is.null(cond.invertedAlleleOrder) & changeEffectDirectionOnAlleleFlip) {
        if(any(cond.invertedAlleleOrder)) cSumstats$EFFECT<-ifelse(cond.invertedAlleleOrder,(cSumstats$EFFECT*-1),cSumstats$EFFECT)
          if(any(colnames(cSumstats)=="SE")){
            meffects.new<-weighted.mean(cSumstats$EFFECT, w = 1/(minv + cSumstats$SE^2), na.rm = T)
            cSumstats.meta<-rbind(cSumstats.meta,list("New effect mean",as.character(meffects.new)))
          } else {
            meffects.new<-mean(cSumstats$EFFECT,na.rm = T)
            cSumstats.meta<-rbind(cSumstats.meta,list("New effect mean",as.character(meffects.new)))
          }
      }
      
    } else {
      ### Does not have EFFECT
      #### Note that EFFECT is not present
      cSumstats.warnings<-c(cSumstats.warnings,"No EFFECT column present.")
      sumstats.meta[iFile,c("no_EFFECT")]<-T
      #### Add empty EFFECT here for consistency
      cSumstats$EFFECT<-NA_real_
    }
    cat("..")
    
    sumstats.meta[iFile,c("changeEffectDirectionOnAlleleFlip")]<-changeEffectDirectionOnAlleleFlip
    
    ## FRQ
    #cond.invertedFRQ<-NULL
    if(any(colnames(cSumstats)=="FRQ")) {
      ### Has FRQ
      #### Check if value is within limits [0,1]
      if(any(cSumstats$FRQ>1) || any(cSumstats$FRQ<0)) {
        stop(paste0('\nThere are FRQ values larger than 1 (',sum(cSumstats$FRQ>1),') or less than 0 (',sum(cSumstats$FRQ<0),') which is outside of the possible FRQ range.'))
      }
      
      ### Invert FRQ based on the previous reference matching
      if(!is.null(cond.invertedAlleleOrder)) cSumstats$FRQ<-ifelse(cond.invertedAlleleOrder, (1-cSumstats$FRQ), cSumstats$FRQ)
      
      ### Compute MAF
      cond.invertedMAF<-cSumstats$FRQ > .5
      cSumstats$MAF<-ifelse(cond.invertedMAF,1-cSumstats$FRQ,cSumstats$FRQ)
      
    } else {
      ### Does not have FRQ
      #### Note that FRQ is not present
      cSumstats.warnings<-c(cSumstats.warnings,"No FRQ column present.")
      sumstats.meta[iFile,c("no_FRQ")]<-T
      #### Add empty FRQ here for consistency
      cSumstats$FRQ<-NA_real_
    }
    cat(".")
    
    if(!is.null(cond.invertedAlleleOrder)) cSumstats.meta<-rbind(cSumstats.meta,list(
      paste0(
        "Modified SNPs; inverted allele order [",ifelse(any(colnames(cSumstats)=="FRQ"),"FRQ",""),",",
        ifelse(changeEffectDirectionOnAlleleFlip,"EFFECT",""),"]"),
      as.character(sum(cond.invertedAlleleOrder))))
    
    if(!is.null(ref) & harmoniseAllelesToReference){
      # Fix A1 and A2 to reflect the reference alleles
      cSumstats$A1<-cSumstats$A1_REF
      cSumstats$A2<-cSumstats$A2_REF
    }
    
    
    
    ## N
    if(any(colnames(cSumstats)=="N_CAS") && any(colnames(cSumstats)=="N_CON")) {
      ### Calculate total N from number of cases and number of controls if they are present. Overwrite any specific total N.
      cSumstats.meta<-rbind(cSumstats.meta,list("N",paste("N_CAS + N_CON")))
      cSumstats$N <- cSumstats$N_CAS + cSumstats$N_CON
    } else if(!is.null(N) & length(N)>=iFile) {
      #cat("\nUsing the explicitly specified N for the whole dataset:",N[iFile])
      cSumstats.meta<-rbind(cSumstats.meta,list("N",paste("Explicitly specified as ",N[iFile])))
      cSumstats$N<-N[iFile]
    } else if(!("N" %in% colnames(cSumstats))) {
      cSumstats.warnings<-c(cSumstats.warnings,"\nNo N column detected!")
      cSumstats.meta<-rbind(cSumstats.meta,list("N","Warning: Not detected!"))
      cSumstats$N<-NA_integer_
    }
    cat(".")
    
    # Some validity checks before the Z-score calculation
    ## Check the values of the P-column
    if("P" %in% colnames(cSumstats)) {
      if((sum(cSumstats$P > 1) + sum(cSumstats$P < 0)) > 100){
        cSumstats.warnings<-c(cSumstats.warnings,"\nThe P column contains numerous values outside of the expected bounds [0,1]. This can indicate that the column is misinterpreted.")
      }
      cat(".")
    }
    
    # Compute Z score (standardised beta)
    if("EFFECT" %in% colnames(cSumstats) && "P" %in% colnames(cSumstats)) {
      cSumstats$Z <- sign(cSumstats$EFFECT) * sqrt(qchisq(cSumstats$P,1,lower=F))
      cSumstats.meta<-rbind(cSumstats.meta,list("Z","Calculated from P and sign(EFFECT)"))
    } else {
      cSumstats.meta<-rbind(cSumstats.meta,list("Z","NOT calculated since no P or EFFECT"))
    }
    cat("..")
    
    # Set SNP to lower case since this seems to be more common
    cSumstats$SNP<-tolower(cSumstats$SNP)
    cat(".")
    
    #TODO Adapt summary statistics to be used for latent factor GWAS - NOT DONE!!!!!!!!
    
    
    # output handling below
    
    #cSumstats<-cSumstats[,-which(names(cSumstats) %in% c("A1_REF","A2_REF"))]
    
    #output.colnames<- c("SNP","N","Z","A1","A2")
    output.colnames<- c("SNP")
    if("N" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N")
    if("Z" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"Z")
    output.colnames<- c(output.colnames,c("A1","A2"))
    #output.colnames<- c(output.colnames,c("A1","A2","A1_ORIG","A2_ORIG"))
    if("CHR" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"CHR")
    if("BP" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"BP")
    if("FRQ" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"FRQ")
    if("MAF" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"MAF")
    output.colnames<- c(output.colnames,c("P"))
    if("EFFECT" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"EFFECT")
    if("SE" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"SE")
    if("N_CAS" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N_CAS")
    if("N_CON" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N_CON")
    
    output.colnames.more<-colnames(cSumstats)[!(colnames(cSumstats) %in% output.colnames)]
    output.colnames.all<-c(output.colnames,output.colnames.more)
    cSumstats<-subset(cSumstats,select = output.colnames.all)
    cSumstats.meta<-rbind(cSumstats.meta,list("SNPs after supermunge",as.character(nrow(cSumstats))))
    
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
    if(!doChrSplit){
      cat("\nSaving supermunged dataset...\n\n")
      write.table(x = cSumstats,file = nfilepath,sep="\t", quote = FALSE, row.names = F)
      nfilepath.gzip<-gzip(nfilepath)
      cat(paste("\nSupermunged dataset saved as", nfilepath.gzip, "in the specified output directory."))
    }
    
    #addition: producing per-chromosome files in a folder, as RAISS columns
    if(doChrSplit) {
      cat("\nSaving supermunged dataset...\n\n")
      if("CHR" %in% colnames(cSumstats)){
        dir.create(paste0(nfilepath,".chr"), showWarnings = FALSE)
        validChromosomes<-c(1:22,"X","Y","XY","MT") #as per Plink standard
        for(chr in validChromosomes){
          output.chr<-output[which(output$CHR==chr),c("SNP","ORIGBP","A1","A2","Z")]
          colnames(output.chr)<-c("rsID","pos","A0","A1","Z")
          write.table(x = output.chr,file = file.path(paste0(nfilepath,".chr"), paste0("z_",traitNames[iFile],"_",chr,".txt")),sep="\t", quote = FALSE, row.names = F)
        }
      } else stop("\nSplit by chromosome specified, but the dataset does not have a CHR column.")
      
      cat(paste("\nOne file per chromosome have been saved under", paste0(nfilepath,".chr"), "in the specified output directory."))
      
    }
    
    sumstats[[iFile]]<-cSumstats
    
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
  
  return(sumstats.meta)
  
}
