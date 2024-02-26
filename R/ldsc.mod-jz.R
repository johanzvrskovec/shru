#Based on the amazing work by Grotzinger, A. D. et al. Nat. Hum. Behav. 3, 513–525 (2019) and Bulik-Sullivan, B. K. et al. Nat. Genet. 47, 291–295 (2015).
#ldsc.mod is a modified copy of the multivariate ldsc function in GenomicSEM: https://github.com/GenomicSEM/GenomicSEM
#Modified by Johan Zvrskovec
# THE ONLY DIFFERENCE HERE IS THE WAY LDSC READS LD SCORES - SHOULD OTHERWISE BE IDENTICAL TO GSEM LDSC AS OF THE TIME OF SYNC WITH THE MAIN REPO Fri Aug 18 13:45:00 2023, commit 5f5bc83c929058fc9c95e2c17e0f55c4454c368f

.LOG <- function(..., file, print = TRUE) {
  msg <- paste0(..., "\n")
  if (print) cat(msg)
  cat(msg, file = file, append = TRUE)
}

.get_renamed_colnames <- function(hold_names, userprovided, checkforsingle=c(), filename, N_provided, log.file,
                                  warnz=FALSE, warn_for_missing=c(), stop_on_missing=c(), utilfuncs=NULL) {
  interpreted_names <- list(
    SNP=c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS", "MARKERNAME", "ID","PREDICTOR","SNP_ID", "VARIANTID", "VARIANT_ID", "RSIDS"),
    A1=c("A1", "ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA","REF"),
    A2=c("A2","ALLELE2","ALLELE0","OTHER_ALLELE","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA", "ALT", "A0"),
    effect=c("OR","B","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT","EST", "BETA1", "LOGOR"),
    INFO=c("INFO", "IMPINFO"),
    P=c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P"),
    N=c("N","WEIGHT","NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N","N_COMPLETE_SAMPLES", "SAMPLESIZE", "NEFF", "N_EFF", "N_EFFECTIVE", "SUMNEFF"),
    MAF=c("MAF", "CEUAF", "FREQ1", "EAF", "FREQ1.HAPMAP", "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU", "EFFECT_ALLELE_FREQ", "FREQ.A1", "A1FREQ", "ALLELEFREQ"),
    Z=c("Z", "ZSCORE", "Z-SCORE", "ZSTATISTIC", "ZSTAT", "Z-STATISTIC"),
    SE=c("STDERR", "SE", "STDERRLOGOR", "SEBETA", "STANDARDERROR"),
    DIRECTION=c("DIRECTION", "DIREC", "DIRE", "SIGN")
  )
  full_names <- list(
    P="P-value",
    A1="effect allele",
    A2="other allele",
    effect="beta or effect",
    SNP="rs-id",
    SE="standard error",
    DIRECTION="direction"
  )
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
      assign(j, utilfuncs[[j]], envir=environment())
    }
  }
  if (all(c("ALT", "REF") %in% hold_names)) {
    .LOG(paste0("Found REF and ALT columns in the summary statistic file ", filename, ". Please note that REF will be interpreted as A1 (effect allele) and ALT as A2 (other allele)"), print=TRUE, file=log.file)
  }
  if (N_provided) {
    interpreted_names[["N"]] <- NULL
  } else {
    if ("NEFF" %in% hold_names | "N_EFF" %in% hold_names | "N_EFFECTIVE" %in% hold_names | "SUMNEFF" %in% hold_names) {
      .LOG("Found an NEFF column for sample size. \n
Please note that this is likely effective sample size and should only be used for liability h^2 conversion for binary traits and that it should reflect the sum of effective sample sizes across cohorts.\n
Be aware that some NEFF columns reflect half of the effective sample size; the function will automatically double the column names if recognized [check above in .log file to determine if this is the case].
If the Neff value is halved in the summary stats, but not recognized by the munge function, this should be manually doubled prior to running munge.", file=log.file)
    }
  }
  for (col in names(interpreted_names)) {
    if (col %in% names(userprovided)) {
      .LOG("Interpreting the ",userprovided[[col]]," column as the ",col, " column, as requested",file=log.file)
      hold_names[ hold_names == toupper(userprovided[[col]]) ] <- col
    } else if (col %in% hold_names) {
      .LOG("Interpreting the ",col," column as the ",col, " column.",file=log.file)
    } else if (any(interpreted_names[[col]] %in% hold_names)) {
      .LOG("Interpreting the ", hold_names[ hold_names %in% interpreted_names[[col]] ], " column as the ",col," column.",file=log.file)
      hold_names[ hold_names %in% interpreted_names[[col]] ] <- col
    } else if ((col == "effect")){
      if (any(interpreted_names[["Z"]] %in% hold_names)) {
        if (!warnz) {
          .LOG("Interpreting the ", hold_names[hold_names %in% interpreted_names[["Z"]] ] , " column as the ",col," column.",file=log.file)
          hold_names[hold_names %in% interpreted_names[["Z"]] ] <- col
        } else {
          .LOG("There appears to be a Z-statistic column in the summary statistic file ", filename, ". Please set linprob to TRUE for binary traits or OLS to true for continuous traits in order to back out the betas or if betas are already available remove this column.", print=FALSE, file=log.file)
          warning(paste0("There appears to be a Z-statistic column in the summary statistic file ", filename, ". Please set linprob to TRUE for binary traits or OLS to true for continuous traits in order to back out the betas or if betas are already available remove this column."))
        }
      }
    } else {
      if (col %in% warn_for_missing) {
        .LOG('Cannot find ', col, ' column, try renaming it to ', col, ' in the summary statistics file for:',filename,file=log.file)
      } else if (col %in% stop_on_missing) {
        stop(paste0('Cannot find ', col, ' column, try renaming it to ', col, ' in the summary statistics file for:',filename))
      }
    }
  }
  # Print log and throw warning messages if multiple or no columns were found for those specified in checkforsingle
  if (length(checkforsingle) > 0) {
    for (col in checkforsingle) {
      if(sum(hold_names == col) == 0) {
        .LOG('Cannot find ',full_names[[col]],' column, try renaming it ', col, ' in the summary statistics file for:',filename,file=log.file)
        warning(paste0('Cannot find ',full_names[[col]],' column, try renaming it ', col, ' in the summary statistics file for:', filename))
      }
      if(sum(hold_names == col) > 1) {
        .LOG('Multiple columns are being interpreted as the ',full_names[[col]],' column, try renaming the column you dont want interpreted to ', col, '2 in the summary statistics file for:',filename,file=log.file)
        warning(paste0('Multiple columns are being interpreted as the ',full_names[[col]],' column, try renaming the column you dont want interpreted to ', col, '2 in the summary statistics file for:', filename))
      }
    }
  }
  return(hold_names)
}

#function to rearrange the sampling covariance matrix from original order to lavaan's order:
#'k' is the number of variables in the model
#'fit' is the fit function of the regression model
#'names' is a vector of variable names in the order you used
.rearrange <- function (k, fit, names) {
  order1 <- names
  order2 <- rownames(inspect(fit)[[1]]) #order of variables
  kst <- k*(k+1)/2
  covA <- matrix(NA, k, k)
  covA[lower.tri(covA, diag = TRUE)] <- 1:kst
  covA <- t(covA)
  covA[lower.tri(covA, diag = TRUE)] <- 1:kst
  colnames(covA) <- rownames(covA) <- order1 #give A actual variable order from lavaan output
  #reorder A by order2
  covA <- covA[order2, order2] #rearrange rows/columns
  vec2 <- lav_matrix_vech(covA) #grab new vectorized order
  return(vec2)
}

##modification of trycatch that allows the results of a failed run to still be saved
.tryCatch.W.E <- function(expr) {
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler), warning = W)
}

.get_V_full <- function(k, V_LD, varSNPSE2, V_SNP) {
  ##create shell of full sampling covariance matrix
  V_Full<-diag(((k+1)*(k+2))/2)
  
  ##input the ld-score regression region of sampling covariance from ld-score regression SEs
  V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD
  
  ##add in SE of SNP variance as first observation in sampling covariance matrix
  V_Full[1,1]<-varSNPSE2
  
  ##add in SNP region of sampling covariance matrix
  V_Full[2:(k+1),2:(k+1)]<-V_SNP
  return(V_Full)
}

.get_V_SNP <- function(SE_SNP, I_LD, varSNP, GC, coords, k, i) {
  V_SNP<-diag(k)
  #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
  if(GC == "conserv"){
    for (p in 1:nrow(coords)) {
      x<-coords[p,1]
      y<-coords[p,2]
      if (x != y) {
        V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i]^2)}
      if (x == y) {
        V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i])^2
      }
    }
  }
  
  if(GC == "standard"){
    for (p in 1:nrow(coords)) {
      x<-coords[p,1]
      y<-coords[p,2]
      if (x != y) {
        V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varSNP[i]^2)}
      if (x == y) {
        V_SNP[x,x]<-(SE_SNP[i,x]*sqrt(I_LD[x,x])*varSNP[i])^2
      }
    }
  }
  
  if(GC == "none"){
    for (p in 1:nrow(coords)) {
      x<-coords[p,1]
      y<-coords[p,2]
      if (x != y) {
        V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*varSNP[i]^2)}
      if (x == y) {
        V_SNP[x,x]<-(SE_SNP[i,x]*varSNP[i])^2
      }
    }
  }
  return(V_SNP)
}

.get_Z_pre <- function(i, beta_SNP, SE_SNP, I_LD, GC) {
  if(GC == "conserv"){
    Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*diag(I_LD))
  }
  if(GC=="standard"){
    Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*sqrt(diag(I_LD)))
  }
  if(GC=="none"){
    Z_pre<-beta_SNP[i,]/SE_SNP[i,]
  }
  return(Z_pre)
}

ldsc.mod <- function(traits, sample.prev, population.prev, ld, wld,
                trait.names = NULL, sep_weights = FALSE, chr = 22,
                n.blocks = 200, ldsc.log = NULL, stand = FALSE,select=FALSE,chisq.max = NA) {
  
  time <- proc.time()
  
  begin.time <- Sys.time()
  
  if(is.null(ldsc.log)){
    logtraits<-gsub(".*/","",traits)
    log2<-paste(logtraits,collapse="_")
    if(object.size(log2) > 200){
      log2<-substr(log2,1,100)
    }
    log.file <- file(paste0(log2, "_ldsc.log"),open="wt")
  }else{log.file<-file(paste0(ldsc.log, "_ldsc.log"),open="wt")}
  
  .LOG("Multivariate ld-score regression of ", length(traits), " traits ", "(", paste(traits, collapse = " "), ")", " began at: ", begin.time, file=log.file)
  
  if(select == "ODD" | select == "EVEN"){
    odd<-seq(1,chr,2)
    even<-seq(2,chr,2)
  }
  
  # Dimensions
  n.traits <- length(traits)
  n.V <- n.traits * (n.traits + 1) / 2
  
  if(n.traits > 18){
    n.blocks<-(((n.traits+1)*(n.traits+2))/2)+1
    .LOG("     ", file=log.file, print = FALSE)
    .LOG("Setting the number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) to ", n.blocks, file=log.file)
    .LOG("This reflects the need to estimate V using at least one more block then their are nonredundant elements in the genetic covariance matrix that includes individual SNPs.", file=log.file)
    .LOG("If the n.blocks is > 1000 you should carefully inspect output for any strange results, such as extremely significant Q_SNP estimates.", file=log.file)
    .LOG("     ", file=log.file, print = FALSE)
    if(n.blocks > 1000){
      warning("The number of blocks needed to estimate V is > 1000, which may result in sampling dependencies across the blocks used to estimate standard errors and can bias results.")
    }
  }

  if(!(is.null(trait.names))){
    check_names<-stringr::str_detect(trait.names, "-")
    if(any(check_names))
      warning("Your trait names specified include mathematical arguments (e.g., + or -) that will be misread by lavaan. Please rename the traits using the trait.names argument.")
  }

  if(length(traits)==1)
    warning("Our version of ldsc requires 2 or more traits. Please include an additional trait.")


  # Storage:
  cov <- matrix(NA,nrow=n.traits,ncol=n.traits)
  V.hold <- matrix(NA,nrow=n.blocks,ncol=n.V)
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  Liab.S <- rep(1, n.traits)
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)


  #########  READ LD SCORES:
  .LOG("Reading in LD scores", file=log.file)

  if(select == FALSE){
    #mod change - read all ld scores in folder, do not use data.table format
    x <- do.call("rbind", lapply(list.files(path = ld, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
        suppressMessages(fread(file = file.path(ld, i), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, data.table=F))
      }))
    
    # x <- do.call("rbind", lapply(1:chr, function(i) {
    #   suppressMessages(read_delim(
    #     file.path(ld, paste0(i, ".l2.ldscore.gz")),
    #     delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    # }))
  }

  if(select == "ODD"){
    x <- do.call("rbind", lapply(odd, function(i) {
      suppressMessages(read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
  }

  if(select == "EVEN"){
    x <- do.call("rbind", lapply(even, function(i) {
      suppressMessages(read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
  }

  if(is.numeric(select)){
    x <- do.call("rbind", lapply(select, function(i) {
      suppressMessages(read_delim(
        file.path(ld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
  }


  x$CM <- NULL
  x$MAF <- NULL


  ######### READ weights:
  if(sep_weights){
    if(select == FALSE){
    w <- do.call("rbind", lapply(1:chr, function(i) {
      suppressMessages(read_delim(
        file.path(wld, paste0(i, ".l2.ldscore.gz")),
        delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
    }))
    }
    if(select == "EVEN"){
      w <- do.call("rbind", lapply(even, function(i) {
        suppressMessages(read_delim(
          file.path(wld, paste0(i, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      }))
    }
    if(select == "ODD"){
      w <- do.call("rbind", lapply(odd, function(i) {
        suppressMessages(read_delim(
          file.path(wld, paste0(i, ".l2.ldscore.gz")),
          delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
      }))
    }
      if(is.numeric(select)){
        w <- do.call("rbind", lapply(select, function(i) {
          suppressMessages(read_delim(
            file.path(wld, paste0(i, ".l2.ldscore.gz")),
            delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
        }))
    }
  }else{w<-x}

  w$CM <- NULL
  w$MAF <- NULL

  colnames(w)[ncol(w)] <- "wLD"

  ### READ M

  if(select == FALSE){
    #mod change - read all m-files in folder, do not use data.table format
    m <- do.call("rbind", lapply(list.files(path = ld, pattern = "\\.l2\\.M_5_50$"), function(i) {
      suppressMessages(fread(file = file.path(ld, i), na.strings =c(".",NA,"NA",""), encoding = "UTF-8", fill = T, blank.lines.skip = T, showProgress = F, nThread = 3, header = F, data.table=F))
    }))
    # m <- do.call("rbind", lapply(1:chr, function(i) {
    #   suppressMessages(read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    # }))
  }

  if(select == "EVEN"){
    m <- do.call("rbind", lapply(even, function(i) {
      suppressMessages(read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
  }

  if(select == "ODD"){
    m <- do.call("rbind", lapply(odd, function(i) {
      suppressMessages(read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
  }

  if(is.numeric(select)){
    m <- do.call("rbind", lapply(select, function(i) {
      suppressMessages(read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
    }))
  }

  M.tot <- sum(m)
  m <- M.tot

  ### READ ALL CHI2 + MERGE WITH LDSC FILES
  s <- 0

  all_y <- lapply(traits, function(chi1) {

    ## READ chi2
    y1 <- suppressMessages(na.omit(read_delim(
      chi1, delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)))

    .LOG("Read in summary statistics [", s <<- s + 1, "/", n.traits, "] from: ", chi1, file=log.file)

    ## Merge files
    merged <- merge(y1[, c("SNP", "N", "Z", "A1")], w[, c("SNP", "wLD")], by = "SNP", sort = FALSE)
    merged <- merge(merged, x, by = "SNP", sort = FALSE)
    merged <- merged[with(merged, order(CHR, BP)), ]

    .LOG("Out of ", nrow(y1), " SNPs, ", nrow(merged), " remain after merging with LD-score files", file=log.file)

    ## REMOVE SNPS with excess chi-square:

    if(is.na(chisq.max)){
    chisq.max <- max(0.001 * max(merged$N), 80)
    }
    rm <- (merged$Z^2 > chisq.max)
    merged <- merged[!rm, ]

    .LOG("Removing ", sum(rm), " SNPs with Chi^2 > ", chisq.max, "; ", nrow(merged), " remain", file=log.file)

    merged
  })

  # count the total nummer of runs, both loops
  s <- 1

  for(j in 1:n.traits){

    chi1 <- traits[j]

    y1 <- all_y[[j]]
    y1$chi1 <- y1$Z^2

    for(k in j:length(traits)){

      ##### HERITABILITY code

      if(j == k){

        .LOG("     ", "     ", file=log.file, print = FALSE)
        .LOG("Estimating heritability [", s, "/", n.V, "] for: ", chi1, file=log.file)

        samp.prev <- sample.prev[j]
        pop.prev <- population.prev[j]

        merged <- y1
        n.snps <- nrow(merged)

        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1


        #### MAKE WEIGHTS:

        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- pmax(merged$L2, 1)
        merged$w.ld <- pmax(merged$wLD, 1)
        merged$c <- tot.agg*merged$N/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        merged$weights <- merged$initial.w/sum(merged$initial.w)

        N.bar <- mean(merged$N)


        ## preweight LD and chi:

        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$chi1*merged$weights)


        ## Perfrom analysis:

        n.annot <- 1


        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)

        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocks),ncol =(n.annot+1))
        colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
        replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
        colnames(xtx)<- colnames(weighted.LD)
        for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))}

        reg <- solve(xtx, xty)
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m

        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete, xty.delete)
        }

        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}

        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)

        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)

        V.hold[,s] <- pseudo.values[,1]
        N.vec[1,s] <- N.bar

        if(is.na(pop.prev)==F & is.na(samp.prev)==F){
          conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
          Liab.S[j] <- conversion.factor
          .LOG("     ", file=log.file, print = FALSE)
          .LOG("Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait ",
              chi1, " when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.", file=log.file)
          .LOG("     ", file=log.file, print = FALSE)
        }

        cov[j,j] <- reg.tot
        I[j,j] <- intercept

        lambda.gc <- median(merged$chi1) / qchisq(0.5, df = 1)
        mean.Chi <- mean(merged$chi1)
        ratio <- (intercept - 1) / (mean.Chi - 1)
        ratio.se <- intercept.se / (mean.Chi - 1)

        .LOG("Heritability Results for trait: ", chi1, file=log.file)
        .LOG("Mean Chi^2 across remaining SNPs: ", round(mean.Chi, 4), file=log.file)
        .LOG("Lambda GC: ", round(lambda.gc, 4), file=log.file)
        .LOG("Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")", file=log.file)
        .LOG("Ratio: ", round(ratio, 4), " (", round(ratio.se, 4), ")", file=log.file)
        .LOG("Total Observed Scale h2: ", round(reg.tot, 4), " (", round(tot.se, 4), ")", file=log.file)
        .LOG("h2 Z: ", format(reg.tot / tot.se, digits = 3), file=log.file)
      }


      ##### GENETIC COVARIANCE code

      if(j != k){

        .LOG("     ", file=log.file, print = FALSE)

        chi2 <- traits[k]
        .LOG("Calculating genetic covariance [", s, "/", n.V, "] for traits: ", chi1, " and ", chi2, file=log.file)

        # Reuse the data read in for heritability
        y2 <- all_y[[k]]
        y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")], by = "SNP", sort = FALSE)

        y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x)
        y$ZZ <- y$Z.y * y$Z.x
        y$chi2 <- y$Z.y^2
        merged <- na.omit(y)
        n.snps <- nrow(merged)

        .LOG(n.snps, " SNPs remain after merging ", chi1, " and ", chi2, " summary statistics", file=log.file)

        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1


        #### MAKE WEIGHTS:

        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N.x)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- pmax(merged$L2, 1)
        merged$w.ld <- pmax(merged$wLD, 1)
        merged$c <- tot.agg*merged$N.x/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)

        tot.agg2 <- (M.tot*(mean(merged$chi2)-1))/mean(merged$L2*merged$N.y)
        tot.agg2 <- max(tot.agg2,0)
        tot.agg2 <- min(tot.agg2,1)
        merged$ld2 <- pmax(merged$L2, 1)
        merged$w.ld2 <- pmax(merged$wLD, 1)
        merged$c2 <- tot.agg2*merged$N.y/M.tot
        merged$het.w2 <- 1/(2*(1+(merged$c2*merged$ld))^2)
        merged$oc.w2 <- 1/merged$w.ld2
        merged$w2 <- merged$het.w2*merged$oc.w2
        merged$initial.w2 <- sqrt(merged$w2)


        merged$weights_cov <- (merged$initial.w + merged$initial.w2)/sum(merged$initial.w + merged$initial.w2 )

        N.bar <- sqrt(mean(merged$N.x)*mean(merged$N.y))

        ## preweight LD and chi:

        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$ZZ *merged$weights_cov)

        ## Perfrom analysis:


        n.annot <- 1


        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)

        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocks),ncol =(n.annot+1))
        colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
        replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
        colnames(xtx)<- colnames(weighted.LD)
        for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))}

        reg <- solve(xtx, xty)
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m

        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete, xty.delete)
        }

        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}

        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)

        V.hold[, s] <- pseudo.values[, 1]
        N.vec[1, s] <- N.bar

        cov[k, j] <- cov[j, k] <- reg.tot
        I[k, j] <- I[j, k] <- intercept

        .LOG("Results for genetic covariance between: ", chi1, " and ", chi2, file=log.file)
        .LOG("Mean Z*Z: ", round(mean(merged$ZZ), 4), file=log.file)
        .LOG("Cross trait Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")", file=log.file)
        .LOG("Total Observed Scale Genetic Covariance (g_cov): ", round(reg.tot, 4), " (", round(tot.se, 4), ")", file=log.file)
        .LOG("g_cov Z: ", format(reg.tot / tot.se, digits = 3), file=log.file)
        .LOG("g_cov P-value: ", format(2 * pnorm(abs(reg.tot / tot.se), lower.tail = FALSE), digits = 5), file=log.file)
      }

      ### Total count
      s <- s + 1
    }
  }


  ## Scale V to N per study (assume m constant)
  # /!\ crossprod instead of tcrossprod because N.vec is a one-row matrix
  v.out <- cov(V.hold) / crossprod(N.vec * (sqrt(n.blocks) / m))

  ### Scale S and V to liability:
  ratio <- tcrossprod(sqrt(Liab.S))
  S <- cov * ratio

  #calculate the ratio of the rescaled and original S matrices
  scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)

  #rescale the sampling correlation matrix by the appropriate diagonals
  V <- v.out * tcrossprod(scaleO)


  #name traits according to trait.names argument
  #use general format of V1-VX if no names provided
  colnames(S) <- if (is.null(trait.names)) paste0("V", 1:ncol(S)) else trait.names

  if(mean(Liab.S)!=1){
    r<-nrow(S)
    SE<-matrix(0, r, r)
    SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(V))

    .LOG(c("     ", "     "), file=log.file, print = FALSE)
    .LOG("Liability Scale Results", file=log.file)

    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j == k){
          .LOG("     ", file=log.file, print = FALSE)
          .LOG("Liability scale results for: ", chi1, file=log.file)
          .LOG("Total Liability Scale h2: ", round(S[j, j], 4), " (", round(SE[j, j], 4), ")", file=log.file)
        }

        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          .LOG("Total Liability Scale Genetic Covariance between ", chi1, " and ",
              chi2, ": ", round(S[k, j], 4), " (", round(SE[k, j], 4), ")", file=log.file)
          .LOG("     ", file=log.file, print = FALSE)
        }
      }
    }
  }


  if(all(diag(S) > 0)){

    ##calculate standardized results to print genetic correlations to log and screen
    ratio <- tcrossprod(1 / sqrt(diag(S)))
    S_Stand <- S * ratio

    #calculate the ratio of the rescaled and original S matrices
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)

    ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in
    # -> not possible because of 'all(diag(S) > 0)'
    # scaleO[is.nan(scaleO)] <- 0

    #rescale the sampling correlation matrix by the appropriate diagonals
    V_Stand <- V * tcrossprod(scaleO)

    #enter SEs from diagonal of standardized V
    r<-nrow(S)
    SE_Stand<-matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand,diag=TRUE)] <-sqrt(diag(V_Stand))


    .LOG(c("     ", "     "), file=log.file, print = FALSE)
    .LOG("Genetic Correlation Results", file=log.file)

    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          .LOG("Genetic Correlation between ", chi1, " and ", chi2, ": ",
              round(S_Stand[k, j], 4), " (", round(SE_Stand[k, j], 4), ")", file=log.file)
          .LOG("     ", file=log.file, print = FALSE)
        }
      }
    }
  }else{
    warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
    .LOG("Your genetic covariance matrix includes traits estimated to have a negative heritability.", file=log.file, print = FALSE)
    .LOG("Genetic correlation results could not be computed due to negative heritability estimates.", file=log.file)
  }

  end.time <- Sys.time()

  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- floor(total.time-mins*60)

  .LOG("     ", file=log.file, print = FALSE)
  .LOG("LDSC finished running at ", end.time, file=log.file)
  .LOG("Running LDSC for all files took ", mins, " minutes and ", secs, " seconds", file=log.file)
  .LOG("     ", file=log.file, print = FALSE)

  flush(log.file)
  close(log.file)

  if(stand){
    list(V=V,S=S,I=I,N=N.vec,m=m,V_Stand=V_Stand,S_Stand=S_Stand)
  } else {
    list(V=V,S=S,I=I,N=N.vec,m=m)
  }
}
