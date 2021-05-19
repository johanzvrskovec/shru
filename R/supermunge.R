#WORK IN PROGRESS
#Based on the amazing work by Grotzinger, A. D. et al. Nat. Hum. Behav. 3, 513–525 (2019) and Bulik-Sullivan, B. K. et al. Nat. Genet. 47, 291–295 (2015).

stdGwasColumnNames <- function(columnNames, stopOnMissing=T, warnOnMultiple=T,
     c.SNP = c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS","MARKERNAME","ID","PREDICTOR","RS"),
     c.A1 = c("A1", "ALLELE1","ALLELE_1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA","REF"),
     c.A2 = c("A2","ALLELE2","ALLELE_2","ALLELE0","OTHER_ALLELE","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA","ALT"),
     c.EFFECT = c("EFFECT","OR","B","BETA","LOG_ODDS","EFFECTS","SIGNED_SUMSTAT","Z","ZSCORE","Z-SCORE","EST","ZSTAT","ZSTATISTIC","GC_ZSCORE"),
     c.INFO = c("INFO"),
     c.P = c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P"),
     c.N = c("N","WEIGHT","NCOMPLETESAMPLES","TOTALSAMPLESIZE","TOTALN","TOTAL_N","N_COMPLETE_SAMPLES"),
     c.N_CAS = c("N_CAS","NCASE","N_CASE","N_CASES","NCAS","NCA"),
     c.N_CON = c("N_CON","NCONTROL","N_CONTROL","N_CONTROLS","N_CON","CONTROLS_N","NCON","NCO"),
     c.MAF = c("MAF","CEUAF","FREQ","FREQ1","EAF", "FREQ1.HAPMAP", "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU", "EFFECT_ALLELE_FREQ","FREQ.A1","FRQ_U","F_U"),
     c.CHR = c("CHR", "CH", "CHROMOSOME", "CHROM"),
     c.BP = c("BP", "ORIGBP")
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
  columnNames[columnNames.upper %in% c.INFO] <- c.INFO[1]
  columnNames[columnNames.upper %in% c.P] <- c.P[1]
  columnNames[columnNames.upper %in% c.N] <- c.N[1]
  columnNames[columnNames.upper %in% c.N_CAS] <- c.N_CAS[1]
  columnNames[columnNames.upper %in% c.N_CON] <- c.N_CON[1]
  columnNames[columnNames.upper %in% c.MAF] <- c.MAF[1]
  columnNames[columnNames.upper %in% c.CHR] <- c.CHR[1]
  columnNames[columnNames.upper %in% c.BP] <- c.BP[1]
  
  if(stopOnMissing){
    # Stop if any of these columns are not found
    if(!any(columnNames=="SNP")) stop("Could not find the 'SNP' column.")
    if(!any(columnNames=="P")) stop("Could not find the P-value column. Standard is 'P'.")
    if(!any(columnNames=="A1")) stop("Could not find the 'A1' column.")
    if(!any(columnNames=="A2")) stop("Could not find the 'A2' column.")
    if(!any(columnNames=="EFFECT")) stop("Could not find the 'EFFECT' column.")
    #if(!any(columnNames=="MAF")) stop("Could not find the 'MAF' column.")
  }
  
  if(warnOnMultiple){
    # Warn if multiple of these columns are found
    if(sum(columnNames=="SNP")>1) warning("Multiple 'SNP' columns found!")
    if(sum(columnNames=="P")>1) warning("Multiple 'P' columns found!")
    if(sum(columnNames=="A1")>1) warning("Multiple 'A1' columns found!")
    if(sum(columnNames=="A2")>1) warning("Multiple 'A2' columns found!")
    if(sum(columnNames=="EFFECT")>1) warning("Multiple 'EFFECT' columns found!")
    if(sum(columnNames=="MAF")>1) warning("Multiple 'MAF' columns found!")
  }
  
  return(data.frame(std=columnNames,orig=columnNames.orig))
}

supermunge <- function(filePaths,refFilePath=NULL,trait.names=NULL,N=NULL,path.dir.output=".",keep.indel=T,doChrSplit=F,doStatistics=F, num=NULL, stopOnMissing=T, maxSNPDistanceBpPadding=0){
  
  
  # for testing
  #library(data.table)
  #library(R.utils) #used for gzip
  #filePaths<-c("/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/cleaned/ALCD03.gz")
  #refFilePath<-"/Users/jakz/Documents/local_db/JZ_GED_PHD_C1/data/combined.hm3_1kg.snplist.vanilla.jz2020.txt"
  #trait.names=c("ALCD03")
  #path.dir.output="/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/munged"
  
  timeStart <- Sys.time()
  
  if(is.null(trait.names)){
    trait.names<-basename(filePaths)
  }
  
  if(!is.null(num)){
    filePaths<-filePaths[1:num]
  }
  
  cat("\n\n\nS U P E R * M U N G E\n")
  
  ref<-NULL
  if(!is.null(refFilePath)){
    cat("\nReading reference file:\n",refFilePath)
    ref <- read.table(refFilePath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
    #rename reference columns as to distinguish them from the dataset columns
    names(ref)<-paste0(names(ref),"_REF")
    ref <- setDT(ref)
    ref$A1_REF <- as.character(toupper(ref$A1_REF))
    ref$A2_REF <- as.character(toupper(ref$A2_REF))
    ref.keys<-c('SNP_REF')
    if('CHR_REF' %in% names(ref)) ref.keys<-c(ref.keys,'CHR_REF')
    if('BP_REF' %in% names(ref)) ref.keys<-c(ref.keys,'BP_REF')
    setkeyv(ref,cols = ref.keys)
    cat("\nReference file read.")
  } else {
    warning("\nRunning without reference to validate with.")
  }
  
  sumstats.metadata<-data.table(name=trait.names,file_path=filePaths,n_snp_raw=NA_integer_,n_snp_res=NA_integer_)
  
  sumstats<-c()
  filePaths.length <- length(filePaths)
  for(iFile in 1:filePaths.length){
    #for testing!
    #iFile=1
    timeStart.ds <- Sys.time()
    cFilePath<-filePaths[iFile]
    cat("\n\nSupermunging",trait.names[iFile],"GWAS summary statistics file:\n", cFilePath)
    cSumstats <- read.table(cFilePath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
    cSumstats.nSNP.raw<-nrow(cSumstats)
    sumstats.metadata[iFile,c("n_snp_raw")]<-cSumstats.nSNP.raw
    cat("\n",trait.names[iFile], "dataset read with",cSumstats.nSNP.raw, "counted SNPs")
    
    # Rename current sumstats with new standardised column names
    cat("\nRead column names:\n",names(cSumstats))
    cSumstats.names <- stdGwasColumnNames(columnNames = names(cSumstats), stopOnMissing = stopOnMissing)
    cat("\nInterpreting column names as:\n")
    apply(cSumstats.names, MARGIN = 1, FUN = function(c){cat(c[2],"\t-> ",c[1],"\n")})
    names(cSumstats) <- cSumstats.names$std
    
    # Transform to data table
    cSumstats <- setDT(cSumstats)
    cSumstats.keys<-c('SNP')
    if('CHR' %in% names(cSumstats)) cSumstats.keys<-c(cSumstats.keys,'CHR')
    if('BP' %in% names(cSumstats)) cSumstats.keys<-c(cSumstats.keys,'BP')
    setkeyv(cSumstats,cols = cSumstats.keys)
    
    # QC, and data management before merge with reference
    
    ## Remove SNPs with missing P
    if("P" %in% colnames(cSumstats)) {
      cSumstats.n<-nrow(cSumstats)
      cSumstats<-cSumstats[which(!is.na(cSumstats$P)),]
      cat("\nRemoved SNPs due to missing P:",(cSumstats.n-nrow(cSumstats)))
    }
    
    ## Remove SNPs with missing effects
    if("EFFECT" %in% colnames(cSumstats)) {
      cSumstats.n<-nrow(cSumstats)
      cSumstats<-cSumstats[which(!is.na(cSumstats$EFFECT)),]
      cat("\nRemoved SNPs due to missing EFFECT:",(cSumstats.n-nrow(cSumstats)))
    }
    
    ## MAF
    if(any(colnames(cSumstats)=="MAF")) {
      ### Has MAF
      #### Check if value is within limits [0,0.5] [0,1]
      if(any(cSumstats$MAF>1) || any(cSumstats$MAF<0)) {
        stop(paste0('\nThere are MAF values larger than 1 (',sum(cSumstats$MAF>1),') or less than 0 (',sum(cSumstats$MAF<0),') which is outside of the possible MAF range.'))
        }
      if(any(cSumstats$MAF>0.5)){
        convertedMAF<-sum(cSumstats$MAF > .5)
        warning(paste0('\nThere are ',convertedMAF,' MAF values larger than 0.5. These will now be converted to 1-VALUE.'))
        sumstats.metadata[iFile,c("converted_MAF")]<-convertedMAF
        cSumstats$MAF<-ifelse(cSumstats$MAF > .5, (1-cSumstats$MAF), cSumstats$MAF)
      }
    } else {
      ### Does not have MAF
      #### Note that MAF is not present
      warning("\nNo MAF column present.")
      sumstats.metadata[iFile,c("no_MAF")]<-T
      #### Add empty MAF here for consistency
      cSumstats$MAF<-NA_real_
    }
    
    ##Alleles, deal with indels
    if(keep.indel == T){
      cSumstats$A1 <- as.character(toupper(cSumstats$A1))
      cSumstats$A2 <- as.character(toupper(cSumstats$A2))
      cat("\nKeeping indel type variants. This may cause problems when aligning variants across traits and references.")
    } else if(keep.indel == F){
      cSumstats$A1 <- as.character(toupper(cSumstats$A1), c("A", "C", "G", "T"))
      cSumstats$A2 <- as.character(toupper(cSumstats$A2), c("A", "C", "G", "T"))
      cat("\nDiscarding indel type variants.")
    }
    
    if(!is.null(ref)){
      #Aligning and validating with reference file
      cSumstats.n<-nrow(cSumstats)
      cat("\nValidating dataset \tnSNP =",cSumstats.n,"\nwith reference \t\tnSNP =", nrow(ref))
      
      #Join with reference on SNP rsID, only keeping SNPs with rsIDs part of the reference
      #https://stackoverflow.com/questions/34644707/left-outer-join-with-data-table-with-different-names-for-key-variables/34645997#34645997
      
      cSumstats.merged.snp<-ref[cSumstats, on=c(SNP_REF='SNP'), nomatch=0]
      #cSumstats[ref, on=c(SNP='SNP_REF'), CHR_REF:=i.CHR_REF] #experimental
      #cSumstats <- merge(ref,cSumstats,by="SNP",all.x=F,all.y=F) #old
      cat("\nRemoved SNPs with rsIDs not present in the reference:\t\t",cSumstats.n-nrow(cSumstats.merged.snp))
      
      if('CHR' %in% names(cSumstats.merged.snp) && 'CHR_REF' %in% names(cSumstats.merged.snp))
      {
        cSumstats.merged.snp.n<-nrow(cSumstats.merged.snp)
        cSumstats.merged.snp<-cSumstats.merged.snp[CHR==CHR_REF]
        cat("\nRemoved SNPs not matching the reference chromosome:\t\t",cSumstats.merged.snp.n-nrow(cSumstats.merged.snp))
      }
      
      if('BP' %in% names(cSumstats.merged.snp) && 'BP_REF' %in% names(cSumstats.merged.snp)){
        cSumstats.merged.snp.maxAlleleLength<-max(nchar(cSumstats.merged.snp$A1),nchar(cSumstats.merged.snp$A2),nchar(cSumstats.merged.snp$A1_REF),nchar(cSumstats.merged.snp$A2_REF))
        cSumstats.merged.snp.n<-nrow(cSumstats.merged.snp)
        cSumstats.merged.snp<-cSumstats.merged.snp[BP < BP_REF + cSumstats.merged.snp.maxAlleleLength + maxSNPDistanceBpPadding & BP > BP_REF - cSumstats.merged.snp.maxAlleleLength - maxSNPDistanceBpPadding]
        cat("\nRemoved SNPs outside the specified bp window +-bp",(cSumstats.merged.snp.maxAlleleLength + maxSNPDistanceBpPadding -1),"from the reference position:\t\t",cSumstats.merged.snp.n-nrow(cSumstats.merged.snp))
      }
      
      #replace missing columns
      cSumstats.merged.snp$SNP<-cSumstats.merged.snp$SNP_REF
      
      if('CHR' %in% names(cSumstats) && 'BP' %in% names(cSumstats) && 'CHR_REF' %in% names(ref) && 'BP_REF' %in% names(ref)) {
        #Join with reference on position rather than rsID
        #cSumstats.merged.pos<-ref[cSumstats, on=c(CHR_REF='CHR' , 'BP' < BP_REF + cSumstats.merged.snp.maxAlleleLength + maxSNPDistanceBpPadding & 'BP' > BP_REF - cSumstats.merged.snp.maxAlleleLength - maxSNPDistanceBpPadding)]
        
        cSumstats.merged.pos<-ref[cSumstats, on=c(CHR_REF='CHR' , BP_REF='BP'), nomatch=0]
        #cSumstats.merged.pos<-cSumstats.merged.pos[!(cSumstats.merged.pos$SNP_REF %in% cSumstats.merged.snp$SNP_REF),]
        
        #replace missing columns
        cSumstats.merged.pos$CHR<-cSumstats.merged.pos$CHR_REF
        cSumstats.merged.pos$BP<-cSumstats.merged.pos$BP_REF
      }
      
      #using cSumstats to store the plain or merged result
      if(is.null(cSumstats.merged.pos)){
        cSumstats<-cSumstats.merged.snp
      } else {
        #merge merged datasets
        cSumstats.merged.pos.salvaged<-cSumstats.merged.pos[!(cSumstats.merged.pos$SNP_REF %in% cSumstats.merged.snp$SNP_REF),]
        cSumstats<-rbindlist(list(cSumstats.merged.snp,cSumstats.merged.pos.salvaged), use.names=T)
        cat("\nSalvaged SNPs by merging on SNP position rather than rsID:",nrow(cSumstats.merged.pos.salvaged))
        cSumstats.merged.pos<-NULL
      }
      cSumstats.merged.snp<-NULL
    }
    
    # More QC and data management, after merge with reference
    
    ## Remove SNPs where alleles are not matching at least one of the reference alleles
    cSumstats.n<-nrow(cSumstats)
    cSumstats<-cSumstats[which(!(cSumstats$A1 != (cSumstats$A1_REF)  & cSumstats$A1 != (cSumstats$A2_REF))), ]
    cat("\nRemoved SNPs where the effect allele (A1) did not match any reference allele:\t\t",cSumstats.n-nrow(cSumstats))
    
    cSumstats.n<-nrow(cSumstats)
    cSumstats<-cSumstats[which(!(cSumstats$A2 != (cSumstats$A1_REF)  & cSumstats$A2 != (cSumstats$A2_REF))), ]
    cat("\nRemoved SNPs where the other allele (A2) did not match any reference allele:\t\t",cSumstats.n-nrow(cSumstats))
    
    ## Invert effect to match the allele order in the reference
    cSumstats$EFFECT<-ifelse(
      (cSumstats$A1 != (cSumstats$A1_REF) & cSumstats$A1 == (cSumstats$A2_REF)) | (cSumstats$A2 != (cSumstats$A2_REF) & cSumstats$A2 == (cSumstats$A1_REF)),
      cSumstats$EFFECT*-1,cSumstats$EFFECT)
    
    ## Determine effect type, and set effect to log(EFFECT) if odds ratio
    if("EFFECT" %in% colnames(cSumstats)) {
      #effect1<-cSumstats$EFFECT[[1]]
      if(round(median(cSumstats$EFFECT,na.rm=T)) == 1) {
        #is odds ratio
        cSumstats$EFFECT<-log(cSumstats$EFFECT)
        sumstats.metadata[iFile,c("effect_type")]<-"OR"
        cat("\nThe EFFECT is determined to be an odds ratio. It has been transformed to a natural logarithmic scale.")
        if(!any(c("OR") %in% toupper(cSumstats.names$orig)) | any(c("B","BETA","LOG_ODDS","Z","ZSCORE","EST","ZSTAT","ZSTATISTIC") %in% toupper(cSumstats.names$orig))) {
          warning("\nThe effect format is not compatible with the original variable naming scheme!")
          sumstats.metadata[iFile,c("effect_type_warning")]<-T
        }
      } else {
        #is NOT odds ratio
        sumstats.metadata[iFile,c("effect_type")]<-"non_OR"
        cat("\nThe EFFECT is determined to NOT be an odds ratio.")
        if(any(c("OR") %in% toupper(cSumstats.names$orig)) | !any(c("B","BETA","LOG_ODDS","Z","ZSCORE","EST","ZSTAT","ZSTATISTIC") %in% toupper(cSumstats.names$orig))) {
          warning("\nThe effect format is not compatible with the original variable naming scheme!")
          sumstats.metadata[iFile,c("effect_type_warning")]<-T
        }
      }
    }
    
    ## N
    if(!is.null(N) & length(N)>=iFile) {
      cat("\nUsing the explicitly specified N for the whole dataset:",N[iFile])
      cSumstats$N<-N[iFile]
    } else if(any(colnames(cSumstats)=="N_CAS") && any(colnames(cSumstats)=="N_CON")) {
      ### Calculate total N from number of cases and number of controls if they are present. Overwrite any specific total N.
      cat("\nDetected separate columns for number of cases and number of controls. Using these to estimate total number of participants for each SNP, overwriting any previous specific values.")
      cSumstats$N <- cSumstats$N_CAS + cSumstats$N_CON
    } else if(!("N" %in% colnames(cSumstats))){
      cSumstats$N<-NA_integer_
    }
    
    # Some validity checks before the Z-score calculation
    
    ## Check the values of the P-column
    if((sum(cSumstats$P > 1) + sum(cSumstats$P < 0)) > 100){
      warning("\nThe P column contains numerous values outside of the expected bounds [0,1]. This can indicate that the column is misinterpreted.")
    }
    
    # Compute Z score (standardised beta)
    cat("\nComputing Z score (standardised beta)")
    if("EFFECT" %in% colnames(cSumstats)) {
      cSumstats$Z <- sign(cSumstats$EFFECT) * sqrt(qchisq(cSumstats$P,1,lower=F))
    }
    
    #output.colnames<- c("SNP","N","Z","A1","A2")
    output.colnames<- c("SNP")
    if("N" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N")
    if("Z" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"Z")
    output.colnames<- c(output.colnames,c("A1","A2","A1_REF","A2_REF"))
    if("CHR" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"CHR")
    if("BP" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"BP")
    if("MAF" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"MAF")
    output.colnames<- c(output.colnames,c("P"))
    if("EFFECT" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"EFFECT")
    if("SE" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"SE")
    if("N_CAS" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N_CAS")
    if("N_CON" %in% colnames(cSumstats)) output.colnames<- c(output.colnames,"N_CON")
    
    output.colnames.more<-colnames(cSumstats)[!(colnames(cSumstats) %in% output.colnames)]
    output.colnames.all<-c(output.colnames,output.colnames.more)
    cSumstats<-subset(cSumstats,select = output.colnames.all)
    cat("\nSNPs left after munge:\t\t",nrow(cSumstats))
    
    nfilepath<-file.path(path.dir.output,trait.names[iFile])
    if(!doChrSplit){
      write.table(x = cSumstats,file = nfilepath,sep="\t", quote = FALSE, row.names = F)
      nfilepath.gzip<-gzip(nfilepath)
      cat("\nThe file is saved as", nfilepath.gzip, "in the specified output directory.")
    }
    
    #addition: producing per-chromosome files in a folder, as RAISS columns
    if(doChrSplit) {
      if("CHR" %in% colnames(cSumstats)){
        dir.create(paste0(nfilepath,".chr"), showWarnings = FALSE)
        validChromosomes<-c(1:22,"X","Y","XY","MT") #as per Plink standard
        for(chr in validChromosomes){
          output.chr<-output[which(output$CHR==chr),c("SNP","ORIGBP","A1","A2","Z")]
          colnames(output.chr)<-c("rsID","pos","A0","A1","Z")
          write.table(x = output.chr,file = file.path(paste0(nfilepath,".chr"), paste0("z_",trait.names[iFile],"_",chr,".txt")),sep="\t", quote = FALSE, row.names = F)
        }
      } else stop("\nSplit by chromosome specified, but the dataset does not have a CHR column.")
      cat("\nOne file per chromosome have been saved under", paste0(nfilepath,".chr"), "in the specified output directory.")
    }
    
    sumstats[[iFile]]<-cSumstats
    
    timeStop.ds <- Sys.time()
    timeDiff <- difftime(time1=timeStop.ds,time2=timeStart.ds,units="sec")
    timeDiff.minutes <- floor(floor(timeDiff)/60)
    timeDiff.seconds <- timeDiff-timeDiff.minutes*60
    
    cat("\nSupermunge of ",trait.names[iFile]," is done in",timeDiff.minutes, "minutes and",timeDiff.seconds," seconds.")
    
  }
  
  timeStop <- Sys.time()
  timeDiff <- difftime(time1=timeStop,time2=timeStart,units="sec")
  timeDiff.minutes <- floor(floor(timeDiff)/60)
  timeDiff.seconds <- timeDiff-timeDiff.minutes*60
  
  cat("\nSupermunge of all datasets is done in",timeDiff.minutes, "minutes and",timeDiff.seconds,"seconds.")
  
  return(sumstats.metadata)
  
}
