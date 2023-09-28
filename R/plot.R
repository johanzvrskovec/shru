#custom pattern loading plot 1
plot.patternLoadings_old<-function(df,bar_aes,color1,color2){
ggplot(
  data = df,
  aes(x=indicator, y=pattern)) + 
  geom_bar(stat = "identity", aes(fill=bar_aes), width = 0.5) +
  scale_fill_manual(name="Indicator pattern loading", 
                    labels = c("Negative", "Positive"), 
                    values = c("-1"=color1, "1"=color2)) +
  geom_text(aes(label=paste(indicator,round(pattern,digits = 3),sep = "\n")), size=4.5) +
  #geom_text(aes(label=paste(indicator,round(pattern,digits = 3),sep = "\n")), vjust=-0.3, size=5) +
  labs(title=paste0("Factor",nFactor),
       y="Std. indicator loading",
       x="Indicator"
  ) +
  theme_minimal()
}

#stolen from https://rpubs.com/danmirman/plotting_factor_analysis
plot.patternLoadings<-function(){
  ggplot(loadings.m, aes(Test, abs(Loading), fill=Loading)) + 
    facet_wrap(~ Factor, nrow=1) + #place the factors in separate facets
    geom_bar(stat="identity") + #make the bars
    coord_flip() + #flip the axes so the test names can be horizontal  
    #define the fill color gradient: blue=positive, red=negative
    scale_fill_gradient2(name = "Loading", 
                         high = "blue", mid = "white", low = "red", 
                         midpoint=0, guide=F) +
    ylab("Loading Strength") + #improve y-axis label
    theme_bw(base_size=10) #use a black-and0white theme with set font size
}


#stolen from https://www.r-graph-gallery.com/101_Manhattan_plot.html
plot.manhattan.custom<-function(
    df,
    maxNLogP=24,
    pointColorValuesVector = rep(c("grey","#66CCCC"),23),
    y_limits=NULL,
    var="P", #or "BETA" or "SE"
    theme.color=list(contrastLight1="#66CCCC",contrastLight2="#FFCC66",contrastDark1="#2D2D2D",contrastDark2="#CC99CC")
    ){
  #df<-cSumstats
  
  df$P<- shru::clipValues(df$P,min = 10^(-maxNLogP), max = NULL)
  df<-df[!is.na(df$P),]
  df<-df[,c("SNP","BP","CHR","P","BETA","SE")]
  #df$CHR<-as.integer(df$CHR)
  
  df<-df[df$P<0.07,]
  
  df <- df %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(BPmax.chr=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(BP.tot.chr=cumsum(as.numeric(BPmax.chr))-BPmax.chr) %>%
  select(-BPmax.chr) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("CHR"="CHR"))

  df <- df %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BP.tot=BP+BP.tot.chr)

  axisdefinition = df %>% group_by(CHR) %>% summarize(center=( max(BP.tot) + min(BP.tot) ) / 2 )
  
  if(var=="BETA"){
    df$y <- df$BETA
  } else if(var=="SE"){
    df$y <- df$SE
  } else { #fallback to P
    df$y <- (-log10(df$P))
  }

  plot <- ggplot(df, aes(x=BP.tot, y=y)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.1) +
  scale_color_manual(values = rep(c("grey", theme.color$contrastLight1), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdefinition$CHR, breaks= axisdefinition$center ) +
  scale_y_continuous(expand = c(0, 0), limits=y_limits ) +     # remove space between plot area and x axis
  
  #geom_label_repel(data=head(setorder(subset(df, y>7),y),n=20), aes(label=SNP), size=3) +
    
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + 
  
  xlab(expression("Chromosome number"))
  
  if(var=="BETA"){
    plot <- plot +
    ylab("BETA")
  } else if(var=="SE"){
    plot <- plot +
    ylab("SE")
  } else { #fallback to P
    
    plot <- plot +
    
    geom_hline(yintercept = 5, color=theme.color$contrastDark1, linetype="dashed") +
      
    geom_hline(yintercept = 8, color=theme.color$contrastDark2, linetype="dashed") +
    
    ylab(expression(paste("-log"[10], plain(P))))
  }
  

  return(plot)
}


#stolen from https://gist.github.com/slowkow/9041570
plot.qq.custom <- function(ps, ci = 0.95) {
  #ps<-head(as.data.frame(cSumstats),n=1000000)$P
  #ps<-cSumstats$P
  n0  <- length(ps)
  ps <- ps[ps<0.07]
  nlogps <- (-log10(sort(ps)))
  #nlogps <- unique(round(-log10(sort(ps)),6))
  
  n  <- length(nlogps)
  df <- data.frame(
    #observed = -log10(sort(ps)),
    observed = nlogps,
    expected = head(-log10(ppoints(n0)),n),
    clower   = head(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n0, shape2 = n0:1)),n),
    cupper   = head(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n0, shape2 = n0:1)),n)
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed),
               #shape = 1,
               size = 1.1, alpha=0.8) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po)
}

#plot for genetic correlations and similar, in matrix form
#requires the corrplot package - NOT INCLUDED IN SHRU
plot.corr <- function(corr, pmat=NULL, SE=NULL, filename, addrect = NULL, is.corr = T, number.cex = 1.0, number.digits=2, newnames=NULL, title=NULL, sig.level=c(0.05)){
  # corr = p$mvLD$S.SE.deltaimpnoimp.500
  # pmat = p$mvLD$covstruct.mvLDSC.1kg.500$cov.p
  # filename = file.path(p$folderpath.plots,"covse.deltaimpnoimp.500.png")
  # is.corr = F
  # number.digits = 4
  # newnames = p$sumstats.sel$code.nice
  
  
  #palette<-colorRampPalette(c(theme.color$contrastDark3,"#FFFFFF",theme.color$contrastLight3))
  
  if(!is.null(newnames)){
    rownames(corr)<-newnames
    colnames(corr)<-newnames
    if(!is.null(pmat)){
      rownames(pmat)<-newnames
      colnames(pmat)<-newnames
    }
  }
  
  corr<-corr[order(colnames(corr)),order(colnames(corr))]
  if(!is.null(pmat)) pmat<-pmat[order(colnames(pmat)),order(colnames(pmat))]
  
  # if(is.null(SE)){
  #   corr.uppCI <- NULL
  #   corr.lowCI <- NULL
  # } else {
  #   corr.uppCI<-clipValues(corr + 1.96 * SE, -1,1)
  #   corr.lowCI<-clipValues(corr - 1.96 * SE, -1,1)
  # }
  
  par.addCoef.col <- NULL
  if(!is.corr) par.addCoef.col <- 'black'
  if(is.corr) par.col.lim <- c(-1,1)
  if(!is.corr) par.col.lim <- c(min(corr), max(corr))
  if(!is.corr & all(corr>=0)) par.col.lim <- c(-max(corr),max(corr))
  if(!is.corr & all(corr<=0)) par.col.lim <- c(min(corr),-min(corr))
  if(all(corr==0)) par.col.lim <- c(-1,1) #fallback if all corr are 0
  
  # if(!is.corr & all(corr>=0)) par.col.lim <- c(-round(max(corr),digits = 0),round(max(corr),digits = 0))
  # if(!is.corr & all(corr<=0)) par.col.lim <- c(round(min(corr)-1,digits = 0),-1*(round(min(corr)-1)))
  
  png(filename = filename, width = 9, height = 9, units = 'in', res = 300, family = "Helvetica")
  par(xpd=TRUE) #keep labels inside margins
  pal<-colorRampPalette(c("#FF6666","#EEEEEE","#6666FF"))
  corrplot(corr = corr, method = "circle", insig='blank', order = "original", p.mat = pmat, sig.level = sig.level, pch.cex = 1.5, full_col=!is.corr, na.label = "square", na.label.col = "grey30", diag=T, addrect = addrect, is.corr = is.corr, number.cex = number.cex, number.digits = number.digits, addCoef.col = par.addCoef.col, col = pal(200), col.lim = par.col.lim, mar=c(0,0,4,0), title = title)$corrPos -> p1  #COL2('RdBu')
  
  #p1$corr <- ifelse(0==-p1$y + ncol(corr)-p1$x + 1, NA_real_, p1$corr)
  if(is.corr) text(p1$x, p1$y, round(p1$corr, number.digits)) #only full coeficients for correlation plots
  
  dev.off()
  
}

# mvldsc = p$ptest$maf0_05_hc1kg
# folderpath.plots = p$folderpath.plots
# code = "maf0_05_hc1kg"
# titleTemplate = "~GSEM LDSC,"
# titleAddition = "\nHC1kG reference, MAF > 0.05"
# mvldscComparison = p$ptest$maf0_05_hm3
# titleAdditionComparison = ",\ncomparing with HM3 reference, MAF > 0.05"

#routine to test and generate plots for Genomic SEM multivariate LDSC results
plotAndTestBatteryForMVLDSC <- function(
    mvldsc,
    code,
    folderpath.plots,
    newnames=NULL,
    titleAddition="",
    titleAdditionComparison="",
    titleTemplate="",
    mvldscComparison=NULL,
    testOnlyTraitNameCodes=NULL){
  
  traitCodes<-colnames(mvldsc$S)
  
  cv.S_StandLimit <- 0.05 #0.05, or 0.1, 0.2 were suggested
  
  pcaRes <- eigen(mvldsc$V_Stand,symmetric = T)
  eigenSum<-sum(abs(pcaRes$values))
  eigenSumLimit<-0.95
  for(iEV in 1:length(pcaRes$values)){
    if(sum(abs(pcaRes$values[1:iEV]))/eigenSum > eigenSumLimit) break
  }
  
  mvldsc$cov.p.fdr2<-matrix(NA,nrow = nrow(mvldsc$cov.p),ncol = ncol(mvldsc$cov.p))
  mvldsc$cov.p.fdr2[lower.tri(mvldsc$cov.p,diag = T)]<-p.adjust2(mvldsc$cov.p[lower.tri(mvldsc$cov.p,diag = T)], method = "fdr",n = iEV)
  mvldsc$cov.p.fdr2[upper.tri(mvldsc$cov.p.fdr2,diag = F)]<-t( mvldsc$cov.p.fdr2)[upper.tri(mvldsc$cov.p.fdr2,diag = F)]
  
  if(!is.null(mvldscComparison)){
    mvldscComparison$cov.p.fdr2<-matrix(NA,nrow = nrow(mvldscComparison$cov.p),ncol = ncol(mvldscComparison$cov.p))
    mvldscComparison$cov.p.fdr2[lower.tri(mvldscComparison$cov.p,diag = T)]<-p.adjust2(mvldscComparison$cov.p[lower.tri(mvldscComparison$cov.p,diag = T)], method = "fdr",n = iEV)
    mvldscComparison$cov.p.fdr2[upper.tri(mvldscComparison$cov.p.fdr2,diag = F)]<-t( mvldscComparison$cov.p.fdr2)[upper.tri(mvldscComparison$cov.p.fdr2,diag = F)]
  }
  
  #fix to set colnames for I, pmat
  if(is.null(colnames(mvldsc$I))){
    colnames(mvldsc$I)<-colnames(mvldsc$S)
    rownames(mvldsc$I)<-rownames(mvldsc$S)
  }
  
  if(is.null(colnames(mvldsc$cov.p.fdr2))){
    colnames(mvldsc$cov.p.fdr2)<-colnames(mvldsc$S)
    rownames(mvldsc$cov.p.fdr2)<-rownames(mvldsc$S)
  }
  
  plot.corr(
    corr = clipValues(mvldsc$I,-3,3),
    pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("i.",code,".png")),
    is.corr = F,
    number.digits = 3,
    newnames = newnames,
    title = paste0("LD score regression intercepts, ",titleTemplate,titleAddition)
  )
  
  plot.corr(
    corr = mvldsc$S,
    pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("covg.",code,".png")),
    is.corr = F,
    number.digits = 3,
    newnames = newnames,
    title = paste0("Genetic covariances, ",titleTemplate,titleAddition)
  )
  
  plot.corr(
    corr = mvldsc$S.SE,
    #pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("covgse.",code,".png")),
    is.corr = F,
    number.digits = 4,
    newnames = newnames,
    title = paste0("Genetic covariances S.E., ",titleTemplate,titleAddition)
  )
  
  S_ForCV<-abs(mvldsc$S)
  g_se <- sqrt(diag(S_ForCV))
  S_se2 <- g_se %*% t(g_se)
  S_ForCV[S_ForCV < S_se2 * cv.S_StandLimit] <- S_se2 * cv.S_StandLimit #cap at min rg 0.05 to avoid extreme values
  plot.corr(
    corr = (mvldsc$S.SE/S_ForCV),
    #pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("covgse.std.",code,".png")),
    is.corr = F,
    number.digits = 4,
    newnames = newnames,
    title = paste0("Genetic covariances (covg) Coefficient of Variation (CV), ",titleTemplate,titleAddition)
  )
  
  #this is just the same as the values from the rg
  # S_StandForCV<-abs(mvldsc$S)
  # S_StandForCV[S_StandForCV<cv.S_StandLimit]<-cv.S_StandLimit
  # plot.corr(
  #   corr = (mvldsc$S.SE/S_StandForCV),
  #   #pmat = mvldsc$cov.p.fdr2,
  #   filename = file.path(folderpath.plots,paste0("covse.std.",code,".png")),
  #   is.corr = F,
  #   number.digits = 4,
  #   newnames = newnames,
  #   title = paste0("Genetic covariances Coefficient of Variation (CV), ",titleTemplate,titleAddition)
  #   )
  
  plot.corr(
    corr = clipValues(mvldsc$S_Stand,-1,1),
    pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("rg.",code,".png")),
    is.corr = T,
    number.digits = 3,
    newnames = newnames,
    title = paste0("Genetic correlations (rg), ",titleTemplate,titleAddition)
  )
  
  plot.corr(
    corr = mvldsc$S_Stand.SE,
    #pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("rgse.",code,".png")),
    is.corr = F,
    number.digits = 4,
    newnames = newnames,
    title = paste0("Genetic correlations (rg) S.E., ",titleTemplate,titleAddition)
  )
  
  
  S_StandForCV<-abs(clipValues(mvldsc$S_Stand,-1,1))
  S_StandForCV[S_StandForCV<cv.S_StandLimit]<-cv.S_StandLimit #cap at min 0.05 to avoid extreme values
  plot.corr(
    corr = (mvldsc$S_Stand.SE/S_StandForCV),
    #pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("rgse.std.",code,".png")),
    is.corr = F,
    number.digits = 4,
    newnames = newnames,
    title = paste0("Genetic correlations (rg) Coefficient of Variation (CV), ",titleTemplate,titleAddition)
  )
  
  # mvldsc<-p$mvLD$covstruct.mvLDSC.1kg.vbcs.drc.winfo
  # mvldscComparison<-p$mvLD$covstruct.mvLDSC.1kg.original
  # testOnlyTraitNameCodes<-infoTraits
  
  if(is.null(testOnlyTraitNameCodes)){
    testOnlyTraitNameCodes<-colnames(mvldsc$S)
  }
  
  if(!is.null(mvldscComparison)){
    cat("\n***Tests for ",code," ***\n")
    
    #Fisher transformation of correlations to make them approximate normal.
    SCorrected<-atanh(mvldsc$S_Stand)
    SCorrected.comparison<-atanh(mvldscComparison$S_Stand)
    
    #corrections for >=1 correlations
    cond<-!is.finite(SCorrected) & mvldsc$S_Stand>=1
    SCorrected[cond]<-sign(mvldsc$S_Stand[cond])*atanh(0.9999999999999999) #set to some large std normal #atanh(0.9999999999999999)
    cond<-!is.finite(SCorrected.comparison) & mvldscComparison$S_Stand>=1
    SCorrected.comparison[cond]<-sign(mvldscComparison$S_Stand[cond])*atanh(0.9999999999999999) #set to some large std normal #atanh(0.9999999999999999)
    
    #test correlation difference
    mTest<-SCorrected-SCorrected.comparison
    varTest<-2*abs(mvldsc$S_Stand.SE^2-mvldscComparison$S_Stand.SE^2) #the intersect
    pTest<-pnorm(q = mTest,sd = sqrt(varTest), lower.tail = F)
    
    pTest.fdr2<-matrix(NA,nrow = nrow(pTest),ncol = ncol(pTest))
    pTest.fdr2[lower.tri(pTest.fdr2,diag = T)]<-p.adjust2(pTest[lower.tri(pTest,diag = T)], method = "fdr", n = iEV)
    pTest.fdr2[upper.tri(pTest.fdr2,diag = F)]<-t(pTest.fdr2)[upper.tri(pTest.fdr2,diag = F)]
    colnames(pTest.fdr2)<-colnames(mvldsc$S)
    rownames(pTest.fdr2)<-colnames(mvldsc$S)
    
    S.SECorrected<-mvldsc$cov.blocks*(mvldsc$S_Stand.SE^2) #chi-squared with df=nblocks, scaled with the mean of the sample (mvldsc$S_Stand) - they are already centered!
    S.SECorrected.comparison<-mvldscComparison$cov.blocks*(mvldscComparison$S_Stand.SE^2)
    
    #try using the #blocks as df, transform to normal - EXPERIMENTAL - NOT DONE!!
    S.SECorrected2<-(S.SECorrected-mvldsc$cov.blocks*abs(mvldsc$S_Stand))/sqrt(2*mvldsc$cov.blocks*abs(mvldsc$S_Stand))
    S.SECorrected.comparison2<-(S.SECorrected.comparison-mvldscComparison$cov.blocks)/sqrt(2*mvldscComparison$cov.blocks)
    
    #test variance difference
    mTestSE<-S.SECorrected2-S.SECorrected.comparison2
    #varTestSE<-mvldsc$cov.blocks+mvldscComparison$cov.blocks #assume these are not covarying
    #varTest<-2*abs(2*mvldsc$cov.blocks-2*mvldscComparison$cov.blocks)# this is very unpractical as they are all 0
    pTestSE<-pnorm(q = mTestSE, sd = 1, lower.tail = T)
    
    
    dS_Stand<-mvldsc$S_Stand[testOnlyTraitNameCodes,testOnlyTraitNameCodes]-mvldscComparison$S_Stand[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    dS.S_Stand.SE<-mvldsc$S_Stand.SE[testOnlyTraitNameCodes,testOnlyTraitNameCodes]-mvldscComparison$S_Stand.SE[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    bartlett.test.labels<-c(rep("mvldsc",length(mvldsc$S[testOnlyTraitNameCodes,testOnlyTraitNameCodes])),rep("comparison",length(mvldscComparison$S[testOnlyTraitNameCodes,testOnlyTraitNameCodes])))
    
    cat("Mean abs(dS_Stand)=",mean(abs(dS_Stand)),"\n")
    #pS<-pchisq(abs(dS_Stand),df = 1, lower.tail = F) #this is pretending that the dist is chi-square
    # result<-bartlett.test(x = c(unlist(mvldsc$S_Stand[testOnlyTraitNameCodes,testOnlyTraitNameCodes]),unlist(mvldscComparison$S_Stand[testOnlyTraitNameCodes,testOnlyTraitNameCodes])), g =bartlett.test.labels)
    # print(result)
    
    # cat("Mean dS.S_Stand.SE=",mean(dS.S_Stand.SE),"\n")
    # 
    # cat("\nBartlett's test of the obtained S.E.\n")
    # result<-bartlett.test(x = c(unlist(mvldsc$S_Stand.SE[testOnlyTraitNameCodes,testOnlyTraitNameCodes]),unlist(mvldscComparison$S_Stand.SE[testOnlyTraitNameCodes,testOnlyTraitNameCodes])), g = bartlett.test.labels)
    # print(result)
    # 
    # cat("\nBartlett's test of the obtained S.E. (standardised)\n")
    # result<-bartlett.test(x = c(unlist((mvldsc$S_Stand.SE/abs(mvldsc$S_Stand))[testOnlyTraitNameCodes,testOnlyTraitNameCodes]),unlist((mvldscComparison$S_Stand.SE/abs(mvldscComparison$S_Stand))[testOnlyTraitNameCodes,testOnlyTraitNameCodes])), g = bartlett.test.labels)
    # print(result)
    
    
    
    #custom pdf difference test for the difference between standard errors or CV
    # mvldsc<-p$mvLD$covstruct.mvLDSC.1kg
    # mvldscComparison<-p$mvLD$covstruct.mvLDSC.hm3
    #variance-to-mean-ratio (relative variance, index of dispersion etc)
    
    S_StandForCV<-abs(mvldsc$S_Stand)
    S_StandForCV[S_StandForCV<cv.S_StandLimit]<-cv.S_StandLimit #cap at min 0.05 to avoid extreme values
    mvldsc$CV <- mvldsc$S_Stand.SE/S_StandForCV #square of this is X2 distr, if times (n-1) and scaled by the population variance
    mvldsc$CV.n <- mvldsc$cov.blocks
    mvldsc$testVar<-(mvldsc$CV^2)*mvldsc$CV.n
    mvldsc$covse.p <- pchisq(q = mvldsc$testVar, df = (mvldsc$CV.n-1),lower.tail = F)
    
    #add in FDR correction
    mvldsc$covse.p.fdr2<-matrix(NA,nrow = nrow(mvldsc$covse.p),ncol = ncol(mvldsc$covse.p))
    mvldsc$covse.p.fdr2[lower.tri(mvldsc$covse.p.fdr2,diag = T)]<-p.adjust2(mvldsc$covse.p[lower.tri(mvldsc$covse.p,diag = T)], method = "fdr", n = iEV)
    mvldsc$covse.p.fdr2[upper.tri(mvldsc$covse.p.fdr2,diag = F)]<-t(mvldsc$covse.p.fdr2)[upper.tri(mvldsc$covse.p.fdr2,diag = F)]
    
    S_StandForCV<-abs(mvldscComparison$S_Stand)
    S_StandForCV[S_StandForCV<cv.S_StandLimit]<-cv.S_StandLimit #cap at min 0.05 to avoid extreme values
    mvldscComparison$CV <- mvldscComparison$S_Stand.SE/S_StandForCV #square of this is X2 distr, if times (n-1) and scaled by the population variance
    mvldscComparison$CV.n <- mvldscComparison$cov.blocks
    mvldscComparison$testVar<-(mvldscComparison$CV^2)*mvldscComparison$CV.n
    mvldscComparison$covse.p <- pchisq(q = mvldscComparison$testVar, df = (mvldscComparison$CV.n-1),lower.tail = F)
    
    #add in FDR correction
    mvldscComparison$covse.p.fdr2<-matrix(NA,nrow = nrow(mvldscComparison$covse.p),ncol = ncol(mvldscComparison$covse.p))
    mvldscComparison$covse.p.fdr2[lower.tri(mvldscComparison$covse.p.fdr2,diag = T)]<-p.adjust2(mvldscComparison$covse.p[lower.tri(mvldscComparison$covse.p,diag = T)], method = "fdr", n = iEV)
    mvldscComparison$covse.p.fdr2[upper.tri(mvldscComparison$covse.p.fdr2,diag = F)]<-t(mvldscComparison$covse.p.fdr2)[upper.tri(mvldscComparison$covse.p.fdr2,diag = F)]
    
    
    cat("Mean dS.S_Stand.SE=",mean(dS.S_Stand.SE),"\n")
    
    covse.testVar.diff<-abs(mvldsc$testVar-mvldscComparison$testVar)
    covse.testVar.diff.n <- abs((mvldsc$CV.n -1) - (mvldscComparison$CV.n-1))
    covse.p.diff <- pchisq(q = covse.testVar.diff, df = covse.testVar.diff.n,lower.tail = F)
    
    covse.p.diff.fdr2<-matrix(NA,nrow = nrow(covse.p.diff),ncol = ncol(covse.p.diff))
    covse.p.diff.fdr2[lower.tri(covse.p.diff.fdr2,diag = T)]<-p.adjust2(covse.p.diff[lower.tri(covse.p.diff,diag = T)], method = "fdr", n = iEV) #use the same number of tests as for the correlation difference - just an approximation as we don't assess the standard error correlations
    covse.p.diff.fdr2[upper.tri(covse.p.diff.fdr2,diag = F)]<-t( covse.p.diff.fdr2)[upper.tri(covse.p.diff.fdr2,diag = F)]
    colnames(covse.p.diff.fdr2)<-colnames(covse.p.diff)
    rownames(covse.p.diff.fdr2)<-colnames(covse.p.diff)
    
    # tdata<-data.frame(x=rnorm(100,0.7,2),y=rnorm(100,-0.2,3))
    # tdata.full<-data.frame(var=c(tdata$x,tdata$y))
    # tdata.full[1:100,c("group")]<-"a"
    # tdata.full[101:200,c("group")]<-"b"
    # tdata.full$group<-as.factor(tdata.full$group)
    
    #levene test test
    #https://www.itl.nist.gov/div898/handbook/eda/section3/eda35a.htm
    #tres <- car::leveneTest(var ~ group, data=tdata.full) #HERE!!! TRY X2 DISTRIBUTIONS INSTEAD
    #mTot<-(sd(tdata$x)+sd(tdata$y))/2
    
    #((100 + 100 - 2)/(2-1)) * (100*(sd(tdata$x)-mTot)^2+100*(sd(tdata$x)-mTot)^2) / (100*var(tdata$x)  + 100*var(tdata$y)) #NO
    
    #Comparison plots
    dI<-mvldsc$I-mvldscComparison$I #re-compute as we need the full matrices here
    dI.totest<-dI[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    plot.corr(
      corr = dI,
      #pmat = pTest.fdr2,
      filename = file.path(folderpath.plots,paste0("di.",code,".png")),
      is.corr = F,
      number.digits = 3,
      newnames = newnames,
      title = paste0("Differences in LD score regression intercepts, ",titleTemplate,titleAddition,titleAdditionComparison," [mean=",round(mean(dI.totest[lower.tri(dI.totest)]),digits = 2),"]")
    )
    
    
    dS_Stand<-clipValues(mvldsc$S_Stand,-1,1)-clipValues(mvldscComparison$S_Stand,-1,1) #re-compute as we need the full matrices here
    dS_Stand.totest<-dS_Stand[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    plot.corr(
      corr = dS_Stand,
      pmat = pTest.fdr2,
      filename = file.path(folderpath.plots,paste0("drg.",code,".png")),
      is.corr = F,
      number.digits = 3,
      newnames = newnames,
      title = paste0("Differences in genetic correlations, ",titleTemplate,titleAddition,titleAdditionComparison," [abs(mean)=",round(mean(abs(dS_Stand.totest[lower.tri(dS_Stand.totest)])),digits = 2),"]")
    )
    
    #relative, in percent
    dS.S_Stand.SE<-100*(mvldsc$S_Stand.SE-mvldscComparison$S_Stand.SE)/abs(mvldscComparison$S_Stand.SE)
    dS.S_Stand.SE.totest<-dS.S_Stand.SE[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    plot.corr(
      corr = dS.S_Stand.SE,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("drgse.rel.",code,".png")),
      is.corr = F,
      number.digits = 1,
      newnames = newnames,
      title = paste0("Relative (%) differences in genetic corrleations S.E., ",titleTemplate,titleAddition,titleAdditionComparison," [mean=",round(mean(dS.S_Stand.SE.totest[lower.tri(dS.S_Stand.SE.totest)]),digits = 1),"]")
    )
    
    #relative, in percent
    dS_rel<-100*(mvldsc$S-mvldscComparison$S)/abs(mvldscComparison$S) #re-compute as we need the full matrices here
    dS_rel.totest<-dS_rel[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    plot.corr(
      corr = dS_rel,
      #pmat = pTest.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovg.rel.",code,".png")),
      is.corr = F,
      number.digits = 1,
      newnames = newnames,
      title = paste0("Relative (%) differences in genetic covariances, ",titleTemplate,titleAddition,titleAdditionComparison," [mean=",round(mean(dS_rel.totest[lower.tri(dS_rel.totest)]),digits = 1),"]")
    )
    
    dS.SE<-mvldsc$S.SE-mvldscComparison$S.SE
    dS.SE.totest<-dS.SE[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    plot.corr(
      corr = dS.SE,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovgse.",code,".png")),
      is.corr = F,
      number.digits = 4,
      newnames = newnames,
      title = paste0("Differences in genetic covariances S.E., ",titleTemplate,titleAddition,titleAdditionComparison," [mean=",round(mean(dS.SE.totest[lower.tri(dS.SE.totest)]),digits = 4),"]")
    )
    
    #relative, in percent
    dS.SE_rel<-100*(mvldsc$S.SE-mvldscComparison$S.SE)/abs(mvldscComparison$S.SE)
    dS.SE_rel.totest<-dS.SE_rel[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    plot.corr(
      corr = dS.SE_rel,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovgse.rel.",code,".png")),
      is.corr = F,
      number.digits = 1,
      newnames = newnames,
      title = paste0("Relative (%) differences in genetic covariances S.E., ",titleTemplate,titleAddition,titleAdditionComparison," [mean=",round(mean(dS.SE_rel.totest[lower.tri(dS.SE_rel.totest)]),digits = 1),"]")
    )
    
    #cov CV, standardised
    comparison_S_ForCV<-abs(mvldscComparison$S)
    comparison_g_se <- sqrt(diag(comparison_S_ForCV))
    comparison_S_se2 <- comparison_g_se %*% t(comparison_g_se)
    comparison_S_ForCV[comparison_S_ForCV < comparison_S_se2 * cv.S_StandLimit] <- comparison_S_se2 * cv.S_StandLimit #cap at min rg 0.05 to avoid extreme values
    dCV_cov <- ((mvldsc$S.SE/S_ForCV) - (mvldscComparison$S.SE/comparison_S_ForCV))
    dCV_cov.totest <- dCV_cov[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    plot.corr(
      corr = dCV_cov,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovse.std.",code,".png")),
      is.corr = F,
      number.digits = 3,
      newnames = newnames,
      title = paste0("Differences in genetic covariances Coefficient of Variation (CV), ",titleTemplate,titleAddition,titleAdditionComparison," [mean=",round(mean(dCV_cov.totest[lower.tri(dCV_cov.totest)]),digits = 3),"]")
    )
    
    #p-values for the difference
    plot.corr(
      corr = covse.p.diff.fdr2,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("drgse.p.",code,".png")),
      is.corr = F,
      number.digits = 3,
      newnames = newnames,
      title = paste0("rg/covg dCV significance (FDR corrected, ",iEV," effective tests), ",titleTemplate,titleAddition,titleAdditionComparison)
    )
    
    #this does not work bc of wrong p
    # plot.corr(
    #   corr = dS.S_Stand.SE,
    #   pmat = pTestSE,
    #   filename = file.path(folderpath.plots,paste0("sdrgse.",code,".png")),
    #   is.corr = F,
    #   number.digits = 4,
    #   newnames = newnames,
    #   title = paste0("Difference in LDSC++ genetic correlations S.E., ",titleTemplate,titleAdditionComparison)
    # )
    
  }
  
  
}




