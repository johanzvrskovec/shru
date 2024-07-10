#Based on the amazing work by Grotzinger, A. D. et al. Nat. Hum. Behav. 3, 513–525 (2019) and Bulik-Sullivan, B. K. et al. Nat. Genet. 47, 291–295 (2015).
#ldsc.mod is a modified copy of the multivariate ldsc function in GenomicSEM: https://github.com/GenomicSEM/GenomicSEM
#Modified by Johan Zvrskovec
# THE ONLY DIFFERENCE HERE IS THE WAY LDSC READS LD SCORES - SHOULD OTHERWISE BE IDENTICAL TO GSEM LDSC AS OF THE TIME OF SYNC WITH THE MAIN REPO Fri Aug 18 13:45:00 2023, commit 5f5bc83c929058fc9c95e2c17e0f55c4454c368f

##modification of trycatch that allows the results of a failed run to still be saved
.tryCatch.W.E_shru <- function(expr) {
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler), warning = W)
}

.LOG_shru <- function(..., file, print = TRUE) {
  msg <- paste0(..., "\n")
  if (print) cat(msg)
  cat(msg, file = file, append = TRUE)
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
  
  .LOG_shru("Multivariate ld-score regression of ", length(traits), " traits ", "(", paste(traits, collapse = " "), ")", " began at: ", begin.time, file=log.file)
  
  if(select == "ODD" | select == "EVEN"){
    odd<-seq(1,chr,2)
    even<-seq(2,chr,2)
  }
  
  # Dimensions
  n.traits <- length(traits)
  n.V <- n.traits * (n.traits + 1) / 2
  
  if(n.traits > 18){
    n.blocks<-(((n.traits+1)*(n.traits+2))/2)+1
    .LOG_shru("     ", file=log.file, print = FALSE)
    .LOG_shru("Setting the number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) to ", n.blocks, file=log.file)
    .LOG_shru("This reflects the need to estimate V using at least one more block then their are nonredundant elements in the genetic covariance matrix that includes individual SNPs.", file=log.file)
    .LOG_shru("If the n.blocks is > 1000 you should carefully inspect output for any strange results, such as extremely significant Q_SNP estimates.", file=log.file)
    .LOG_shru("     ", file=log.file, print = FALSE)
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
  .LOG_shru("Reading in LD scores", file=log.file)

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

    .LOG_shru("Read in summary statistics [", s <<- s + 1, "/", n.traits, "] from: ", chi1, file=log.file)

    ## Merge files
    merged <- merge(y1[, c("SNP", "N", "Z", "A1")], w[, c("SNP", "wLD")], by = "SNP", sort = FALSE)
    merged <- merge(merged, x, by = "SNP", sort = FALSE)
    merged <- merged[with(merged, order(CHR, BP)), ]

    .LOG_shru("Out of ", nrow(y1), " SNPs, ", nrow(merged), " remain after merging with LD-score files", file=log.file)

    ## REMOVE SNPS with excess chi-square:

    if(is.na(chisq.max)){
    chisq.max <- max(0.001 * max(merged$N), 80)
    }
    rm <- (merged$Z^2 > chisq.max)
    merged <- merged[!rm, ]

    .LOG_shru("Removing ", sum(rm), " SNPs with Chi^2 > ", chisq.max, "; ", nrow(merged), " remain", file=log.file)

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

        .LOG_shru("     ", "     ", file=log.file, print = FALSE)
        .LOG_shru("Estimating heritability [", s, "/", n.V, "] for: ", chi1, file=log.file)

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
          .LOG_shru("     ", file=log.file, print = FALSE)
          .LOG_shru("Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait ",
              chi1, " when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.", file=log.file)
          .LOG_shru("     ", file=log.file, print = FALSE)
        }

        cov[j,j] <- reg.tot
        I[j,j] <- intercept

        lambda.gc <- median(merged$chi1) / qchisq(0.5, df = 1)
        mean.Chi <- mean(merged$chi1)
        ratio <- (intercept - 1) / (mean.Chi - 1)
        ratio.se <- intercept.se / (mean.Chi - 1)

        .LOG_shru("Heritability Results for trait: ", chi1, file=log.file)
        .LOG_shru("Mean Chi^2 across remaining SNPs: ", round(mean.Chi, 4), file=log.file)
        .LOG_shru("Lambda GC: ", round(lambda.gc, 4), file=log.file)
        .LOG_shru("Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")", file=log.file)
        .LOG_shru("Ratio: ", round(ratio, 4), " (", round(ratio.se, 4), ")", file=log.file)
        .LOG_shru("Total Observed Scale h2: ", round(reg.tot, 4), " (", round(tot.se, 4), ")", file=log.file)
        .LOG_shru("h2 Z: ", format(reg.tot / tot.se, digits = 3), file=log.file)
      }


      ##### GENETIC COVARIANCE code

      if(j != k){

        .LOG_shru("     ", file=log.file, print = FALSE)

        chi2 <- traits[k]
        .LOG_shru("Calculating genetic covariance [", s, "/", n.V, "] for traits: ", chi1, " and ", chi2, file=log.file)

        # Reuse the data read in for heritability
        y2 <- all_y[[k]]
        y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")], by = "SNP", sort = FALSE)

        y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x)
        y$ZZ <- y$Z.y * y$Z.x
        y$chi2 <- y$Z.y^2
        merged <- na.omit(y)
        n.snps <- nrow(merged)

        .LOG_shru(n.snps, " SNPs remain after merging ", chi1, " and ", chi2, " summary statistics", file=log.file)

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

        .LOG_shru("Results for genetic covariance between: ", chi1, " and ", chi2, file=log.file)
        .LOG_shru("Mean Z*Z: ", round(mean(merged$ZZ), 4), file=log.file)
        .LOG_shru("Cross trait Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")", file=log.file)
        .LOG_shru("Total Observed Scale Genetic Covariance (g_cov): ", round(reg.tot, 4), " (", round(tot.se, 4), ")", file=log.file)
        .LOG_shru("g_cov Z: ", format(reg.tot / tot.se, digits = 3), file=log.file)
        .LOG_shru("g_cov P-value: ", format(2 * pnorm(abs(reg.tot / tot.se), lower.tail = FALSE), digits = 5), file=log.file)
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

    .LOG_shru(c("     ", "     "), file=log.file, print = FALSE)
    .LOG_shru("Liability Scale Results", file=log.file)

    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j == k){
          .LOG_shru("     ", file=log.file, print = FALSE)
          .LOG_shru("Liability scale results for: ", chi1, file=log.file)
          .LOG_shru("Total Liability Scale h2: ", round(S[j, j], 4), " (", round(SE[j, j], 4), ")", file=log.file)
        }

        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          .LOG_shru("Total Liability Scale Genetic Covariance between ", chi1, " and ",
              chi2, ": ", round(S[k, j], 4), " (", round(SE[k, j], 4), ")", file=log.file)
          .LOG_shru("     ", file=log.file, print = FALSE)
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


    .LOG_shru(c("     ", "     "), file=log.file, print = FALSE)
    .LOG_shru("Genetic Correlation Results", file=log.file)

    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          .LOG_shru("Genetic Correlation between ", chi1, " and ", chi2, ": ",
              round(S_Stand[k, j], 4), " (", round(SE_Stand[k, j], 4), ")", file=log.file)
          .LOG_shru("     ", file=log.file, print = FALSE)
        }
      }
    }
  }else{
    warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
    .LOG_shru("Your genetic covariance matrix includes traits estimated to have a negative heritability.", file=log.file, print = FALSE)
    .LOG_shru("Genetic correlation results could not be computed due to negative heritability estimates.", file=log.file)
  }

  end.time <- Sys.time()

  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- floor(total.time-mins*60)

  .LOG_shru("     ", file=log.file, print = FALSE)
  .LOG_shru("LDSC finished running at ", end.time, file=log.file)
  .LOG_shru("Running LDSC for all files took ", mins, " minutes and ", secs, " seconds", file=log.file)
  .LOG_shru("     ", file=log.file, print = FALSE)

  flush(log.file)
  close(log.file)

  if(stand){
    list(V=V,S=S,I=I,N=N.vec,m=m,V_Stand=V_Stand,S_Stand=S_Stand)
  } else {
    list(V=V,S=S,I=I,N=N.vec,m=m)
  }
}
