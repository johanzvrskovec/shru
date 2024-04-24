# #test
# columnNames = cSumstats.names$orig
# missingEssentialColumnsStop = missingEssentialColumnsStop
# ancestrySetting = ancestrySetting[iFile]

stdGwasColumnNames <- function(
    columnNames,
    missingEssentialColumnsStop=c("SNP","A1","A2"),
    ancestrySetting=NA, #EUR, #ancestry setting string for the current dataset
    warnings=T
){
  
  c.SNP = c("SNP","PREDICTOR","SNPID","MARKERNAME","MARKER_NAME","SNPTESTID","ID_DBSNP49","ID","MARKER","SNP.NAME","SNP ID", "SNP_ID","LOCATIONALID","ASSAY_NAME","VARIANT_ID","VARIANT")
  c.RSID = c("RSID","RS_NUMBER","RS","RSNUMBER","RS_NUMBERS","RSID_UKB")
  c.A1 = c("A1","ALLELE1","ALLELE_1","A_1","A","ALLELE.1")
  c.A2 = c("A2","ALLELE2","ALLELE_2","A_2","ALLELE.2")
  c.A0 = c("A0","ALLELE0","ALLELE_0","A_0")
  c.AEFFECT = c("INC_ALLELE","EA","REF","A1_EFFECT","EFFECT_ALLELE","RISK_ALLELE","EFFECTALLELE","EFFECT_ALL","REFERENCE_ALLELE","REF_ALLELE","REFERENCEALLELE","EA","INC_ALLELE","CODED_ALLELE","TESTED_ALLELE","EFFECT.ALLELE..EA.")
  c.ANOEFFECT = c("OTHER_ALLELE","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA","ALT","A1_OTHER","A2_OTHER","NONREF_ALLELE","NEFFECT_ALLELE","NEFFECTALLELE","NONEFFECT_ALLELE","OTHER_ALL","OTHERALLELE","NONEFFECTALLELE","ALT_ALLELE","NONCODED_ALLELE")
  c.BETA = c("BETA","B","EFFECT_BETA","EFFECT","EFFECTS","SIGNED_SUMSTAT","EST","GWAS_BETA","EFFECT_A1","EFFECTA1","EFFECT_NW","STDBETA")
  c.OR = c("OR","LOG_ODDS","ODDS-RATIO","ODDS_RATIO","ODDSRATIO","OR(MINALLELE)","OR.LOGISTIC","OR_RAN","OR(A1)","LOGOR")
  c.SE = c("SE","STDER","STDERR","STD","STANDARD_ERROR","OR_SE","STANDARDERROR", "STDERR_NW","META.SE","SE_DGC","SE.2GC","STDERRLOGOR")
  c.Z = c("Z","ZSCORE","Z-SCORE","ZSTAT","ZSTATISTIC","GC_ZSCORE","BETAZSCALE","TEST.STATISTIC","Z_ESTIMATE")
  c.INFO = c("INFO","IMPINFO","IMPQUALITY", "INFO.PLINK", "INFO_UKBB","INFO_UKB","MININFO")
  c.SINFO = c("SINFO")
  c.P = c("P","PVALUE","PVAL","P_VALUE","GC_PVALUE","P.2GC","P.VAL","GWAS_P","P-VALUE","P-VAL","FREQUENTIST_ADD_PVALUE","P.VALUE","P_VAL","SCAN-P","P.LMM","META.PVAL","P_RAN","P.ADD","P_BOLT_LMM","PVAL_ESTIMATE","WALD_P","P_WALD")
  c.NL10P = c("NL10P","NEGLOG10_PVAL") # -log10 p-value, used in the Pan UKB format
  c.N = c("N","NCOMPLETESAMPLES","TOTALSAMPLESIZE","TOTALN","TOTAL_N","N_COMPLETE_SAMPLES","N_TOTAL","N_SAMPLES","N_ANALYZED","NSAMPLES","SAMPLESIZE","SAMPLE_SIZE","TOTAL_SAMPLE_SIZE","TOTALSAMPLESIZE","SAMPLE.SIZE","N_SUM","SIZE")
  c.N_CAS = c("N_CAS","NCASE","N_CASE","N_CASES","NCAS","NCA","NCASES","CASES","CASES_N","FRQ_A")
  c.N_CON = c("N_CON","NCONTROL","N_CONTROL","N_CONTROLS","NCON","NCO","N_CON","NCONTROLS","CONTROLS","CONTROLS_N","FRQ_U")
  c.NEFF = c("NEFF","NEF","NEFFECTIVE","NE","EFFECTIVE_N","WEIGHT")
  c.NEFF_HALF = c("NEFF_HALF","NEFFDIV2")
  #include FRQ_A?
  c.FRQ = c("FRQ","MAF","AF","CEUAF","FREQ","FREQ1","EAF","FREQ1.HAPMAP","FREQALLELE1HAPMAPCEU","FREQ.HAPMAP.CEU","FREQ.ALLELE1.HAPMAPCEU","EFFECT_ALLELE_FREQ","FREQ.A1","F_A","F_U","FREQ_A","FREQ_U","MA_FREQ","MAF_NW","FREQ_A1","A1FREQ","CODED_ALLELE_FREQUENCY","FREQ_TESTED_ALLELE","FREQ_TESTED_ALLELE_IN_HRS","EAF_HRC","EAF_UKB","EAF_EUR_UKB","FREQ_TESTED_ALLELE","EFFECT_ALLELE_FREQUENCY","EFFECT.ALLELE.FREQUENCY..EAF.")
  c.FRQ_CAS = c("FRQ_CAS","FCAS")
  c.FRQ_CON = c("FRQ_CON","FCON")
  c.CHR = c("CHR","CH","CHROMOSOME","CHROM","CHR_BUILD38","CHR_BUILD37","CHR_BUILD36","CHR_B38","CHR_B37","CHR_B36","CHR_ID","SCAFFOLD","HG19CHR","CHR.HG19","CHR_HG19","HG18CHR","CHR.HG18","CHR_HG18","CHR_BP_HG19B37","HG19CHRC","#CHROM","X.CHROM")
  c.BP = c("BP","BP1","ORIGBP","POS","POSITION","LOCATION","PHYSPOS","GENPOS","CHR_POSITION","POS_B38","POS_BUILD38","POS_B37","POS_BUILD37","BP_HG19B37","POS_B36","POS_BUILD36","POS.HG19","POS.HG18","POS_HG19","POS_HG18","BP_HG19","BP_HG18","BP.GRCH38","BP.GRCH37","POSITION(HG19)","POSITION(HG18)","POS(B38)","POS(B37)","BASE_PAIR_LOCATION")
  c.BP2 =c("BP2")
  c.L2 =c("L2","LD")
  c.DF = c("DF","CHISQ_DF")
  
  #force NA for ancestrySetting="ANY"
  if(toupper(ancestrySetting)=="ANY") ancestrySetting<-NA_character_
  
  columnNames<-as.character(columnNames)
  columnNames.upper<-toupper(columnNames)
  names(columnNames)<-columnNames.upper
  columnNames.orig<-columnNames
  
  if(any(columnNames.upper %in% c.RSID)){
    columnNames[columnNames.upper %in% c.SNP] <- "SNPALT"
    columnNames[columnNames.upper %in% c.RSID] <- c.SNP[1]
  } else {
    columnNames[columnNames.upper %in% c.SNP] <- c.SNP[1]
  }
  
  if(any(columnNames.upper %in% c.A0)){
    columnNames[columnNames.upper %in% c.A0] <- c.A2[1]
    columnNames[columnNames.upper %in% c.A1] <- c.A1[1]
  } else {
    columnNames[columnNames.upper %in% c.A1] <- c.A1[1]
    columnNames[columnNames.upper %in% c.A2] <- c.A2[1]
  }
  
  columnNames[columnNames.upper %in% c.AEFFECT] <- c.A1[1]
  columnNames[columnNames.upper %in% c.ANOEFFECT] <- c.A2[1]
  columnNames[columnNames.upper %in% c.BETA] <- c.BETA[1]
  columnNames[columnNames.upper %in% c.OR] <- c.OR[1] 
  columnNames[columnNames.upper %in% c.Z] <- c.Z[1] 
  columnNames[columnNames.upper %in% c.SE] <- c.SE[1]
  columnNames[columnNames.upper %in% c.INFO] <- c.INFO[1]
  columnNames[columnNames.upper %in% c.SINFO] <- c.SINFO[1]
  columnNames[columnNames.upper %in% c.P] <- c.P[1]
  columnNames[columnNames.upper %in% c.N] <- c.N[1]
  columnNames[columnNames.upper %in% c.N_CAS] <- c.N_CAS[1]
  columnNames[columnNames.upper %in% c.N_CON] <- c.N_CON[1]
  columnNames[columnNames.upper %in% c.NEFF] <- c.NEFF[1]
  columnNames[columnNames.upper %in% c.NEFF_HALF] <- c.NEFF_HALF[1]
  
  columnNames[columnNames.upper %in% c.FRQ] <- c.FRQ[1]
  columnNames[columnNames.upper %in% c.FRQ_CAS] <- c.FRQ_CAS[1]
  columnNames[columnNames.upper %in% c.FRQ_CON] <- c.FRQ_CON[1]
  #daner FRQ
  danerNcas<-NA_integer_
  danerNcon<-NA_integer_
  if(!any(columnNames=="FRQ") & any(startsWith(columnNames,prefix = "FRQ_A_")) & any(startsWith(columnNames,prefix = "FRQ_U_")))
  {
    danerNcasS <- columnNames[
      startsWith(columnNames,
                 prefix = "FRQ_A_")
    ][1]
    danerNcas <- as.integer(strsplit(danerNcasS,split = "_",fixed = T)[[1]][3])
    
    danerNconS <- columnNames[
      startsWith(columnNames,
                 prefix = "FRQ_U_")
    ][1]
    danerNcon <- as.integer(strsplit(danerNconS,split = "_",fixed = T)[[1]][3])
    
    if(is.finite(danerNcas) & is.finite(danerNcon))
    {
      columnNames[columnNames.upper %in% danerNcasS] <- "FRQ_CAS"
      columnNames[columnNames.upper %in% danerNconS] <- "FRQ_CON"
    }
  }
  
  #ancestry specific FRQ_CAS - Pan UKB compatibility, TODO!!! 
  
  
  #ancestry specific FRQ
  if(!is.na(ancestrySetting)){
    iMatch <- unlist(lapply(
      X = c.FRQ,
      FUN = function(x){grep(pattern = paste0("^",x,"[\\._]",toupper(ancestrySetting)),columnNames)}
    ))
    if(length(iMatch)>0){
      columnNames[iMatch] <- c.FRQ[1]
      iDup<-unlist(lapply(
        X = c.FRQ,
        FUN = function(x){grep(pattern = paste0("^",x,"[\\._].+"),columnNames)}
      ))
      if(length(iDup)>0){
        columnNames[iDup[1]]<-paste0(c.FRQ[1],"FB")
        if(length(iDup)>1){
          iDup<-iDup[2:length(iDup)]
          columnNames[iDup]<-paste0("X",c.FRQ[1],"FB")
        }
      }
      
    } else {
      if(warnings) warning("\nCould not find the ancestry specific 'FRQ' column.\n")
    }
  } else if(!any(columnNames == c.FRQ[1])) {
    #fallback
    iDup<- unlist(lapply(
      X = c.FRQ,
      FUN = function(x){grep(pattern = paste0("^",x,"\\..+"),columnNames)}
    ))
    if(length(iDup)>0) columnNames[iDup[1]] <- c.FRQ[1]
  }
  
  columnNames[columnNames.upper %in% c.CHR] <- c.CHR[1]
  columnNames[columnNames.upper %in% c.BP] <- c.BP[1]
  columnNames[columnNames.upper %in% c.BP2] <- c.BP2[1]
  columnNames[columnNames.upper %in% c.DF] <- c.DF[1]
  
  
  columnNames[columnNames.upper %in% c.L2] <- c.L2[1]
  #ancestry specific L2
  if(!is.na(ancestrySetting)){
    iMatch <- unlist(lapply(
      X = c.L2,
      FUN = function(x){grep(pattern = paste0("^",x,"[\\._]",toupper(ancestrySetting)),columnNames)}
    ))
    if(length(iMatch)>0){
      columnNames[iMatch] <- c.L2[1]
      iDup<-unlist(lapply(
        X = c.L2,
        FUN = function(x){grep(pattern = paste0("^",x,"[\\._].+"),columnNames)}
      ))
      if(length(iDup)>0){
        columnNames[iDup[1]]<-paste0(c.L2[1],"FB")
        if(length(iDup)>1){
          iDup<-iDup[2:length(iDup)]
          columnNames[iDup]<-paste0("X",c.L2[1],"FB")
        }
      }
    } else {
      if(warnings) warning("\nCould not find the ancestry specific 'L2' column.\n")
    }
  } else if(!any(columnNames == c.L2[1])) {
    #fallback
    iDup<- unlist(lapply(
      X = c.L2,
      FUN = function(x){grep(pattern = paste0("^",x,"\\..+"),columnNames)}
    ))
    if(length(iDup)>0) columnNames[iDup[1]] <- c.L2[1]
  }
  
  if(length(missingEssentialColumnsStop)>0){
    # Stop if any of these columns are not found
    if(!all(missingEssentialColumnsStop %in% columnNames)) stop(paste0("\nNot all essential columns found: ",missingEssentialColumnsStop[!missingEssentialColumnsStop %in% columnNames],"\n"))
  }
  
  if(!any(columnNames=="P") & warnings) warning("\nCould not find the P-value column. Standard is 'P'.\n")
  if(!any(columnNames=="BETA") & !any(columnNames=="OR" & !any(columnNames=="Z")) & warnings) warning("Could not find any effect column.\n")
  if(!any(columnNames=="SNP") & warnings) warning("\nCould not find the 'SNP' column.\n")
  if(!any(columnNames=="A1") & warnings) warning("\nCould not find the 'A1' column.\n")
  if(!any(columnNames=="A2") & warnings) warning("\nCould not find the 'A2' column.\n")
  if(!any(columnNames=="FRQ") & warnings) warning("\nCould not find the 'FRQ' column.\n")
  
  # Warn if multiple of these columns are found
  if(sum(columnNames=="SNP")>1 & warnings) warning("\nMultiple 'SNP' columns found!\n")
  if(sum(columnNames=="P")>1 & warnings) warning("\nMultiple 'P' columns found!\n")
  if(sum(columnNames=="A1")>1 & warnings) warning("\nMultiple 'A1' columns found!\n")
  if(sum(columnNames=="A2")>1 & warnings) warning("\nMultiple 'A2' columns found!\n")
  if(sum(columnNames=="BETA")>1 & warnings) warning("\nMultiple 'BETA' columns found!\n")
  if(sum(columnNames=="OR")>1 & warnings) warning("\nMultiple 'OR' columns found!\n")
  if(sum(columnNames=="Z")>1 & warnings) warning("\nMultiple 'Z' columns found!\n")
  if(sum(columnNames=="FRQ")>1 & warnings) warning("\nMultiple 'FRQ' columns found!\n")
  
  
  
  #rename duplicte columns
  
  if(c.SNP[1] %in% columnNames){
    columnNames.sorted<-columnNames[c.SNP]
    columnNames.sorted<-columnNames.sorted[!is.na(columnNames.sorted)]
    if(length(columnNames.sorted)>1){
      columnNames.sorted<-columnNames.sorted[2:length(columnNames.sorted)]
      columnNames[names(columnNames.sorted)]<-paste0("X",c.SNP[1])
    }
  }
  
  if(c.P[1] %in% columnNames){
    columnNames.sorted<-columnNames[c.P]
    columnNames.sorted<-columnNames.sorted[!is.na(columnNames.sorted)]
    if(length(columnNames.sorted)>1){
      columnNames.sorted<-columnNames.sorted[2:length(columnNames.sorted)]
      columnNames[names(columnNames.sorted)]<-paste0("X",c.P[1])
    }
  }
  
  if(c.A1[1] %in% columnNames){
    columnNames.sorted<-columnNames[c.A1]
    columnNames.sorted<-columnNames.sorted[!is.na(columnNames.sorted)]
    if(length(columnNames.sorted)>1){
      columnNames.sorted<-columnNames.sorted[2:length(columnNames.sorted)]
      columnNames[names(columnNames.sorted)]<-paste0("X",c.A1[1])
    }
  }
  
  if(c.A2[1] %in% columnNames){
    columnNames.sorted<-columnNames[c.A2]
    columnNames.sorted<-columnNames.sorted[!is.na(columnNames.sorted)]
    if(length(columnNames.sorted)>1){
      columnNames.sorted<-columnNames.sorted[2:length(columnNames.sorted)]
      columnNames[names(columnNames.sorted)]<-paste0("X",c.A2[1])
    }
  }
  
  if(c.BETA[1] %in% columnNames){
    columnNames.sorted<-columnNames[c.BETA]
    columnNames.sorted<-columnNames.sorted[!is.na(columnNames.sorted)]
    if(length(columnNames.sorted)>1){
      columnNames.sorted<-columnNames.sorted[2:length(columnNames.sorted)]
      columnNames[names(columnNames.sorted)]<-paste0("X",c.BETA[1])
    }
  }
  
  if(c.OR[1] %in% columnNames){
    columnNames.sorted<-columnNames[c.OR]
    columnNames.sorted<-columnNames.sorted[!is.na(columnNames.sorted)]
    if(length(columnNames.sorted)>1){
      columnNames.sorted<-columnNames.sorted[2:length(columnNames.sorted)]
      columnNames[names(columnNames.sorted)]<-paste0("X",c.OR[1])
    }
  }
  
  if(c.Z[1] %in% columnNames){
    columnNames.sorted<-columnNames[c.Z]
    columnNames.sorted<-columnNames.sorted[!is.na(columnNames.sorted)]
    if(length(columnNames.sorted)>1){
      columnNames.sorted<-columnNames.sorted[2:length(columnNames.sorted)]
      columnNames[names(columnNames.sorted)]<-paste0("X",c.Z[1])
    }
  }
  
  if(c.FRQ[1] %in% columnNames){
    columnNames.sorted<-columnNames[c.FRQ]
    columnNames.sorted<-columnNames.sorted[!is.na(columnNames.sorted)]
    if(length(columnNames.sorted)>1){
      columnNames.sorted<-columnNames.sorted[2:length(columnNames.sorted)]
      columnNames[names(columnNames.sorted)]<-paste0("X",c.FRQ[1])
    }
  }
  
  return(list(std=as.character(columnNames),orig=as.character(columnNames.orig),danerNcas=danerNcas,danerNcon=danerNcon))
  
}