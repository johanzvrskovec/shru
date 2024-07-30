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
    maxP=0.07,
    pointColorValuesVector = rep(c("grey","#66CCCC"),23),
    y_limits=NULL,
    var="P", #or "BETA" or "SE"
    theme.color=list(contrastLight1="#66CCCC",contrastLight2="#FFCC66",contrastDark1="#2D2D2D",contrastDark2="#CC99CC"),
    title=NULL
    ){
  #df<-cSumstats
  
  df<-as.data.frame(df) #in case of data.table or similar
  if(any(colnames(df)=="P")) df<-df[!is.na(df$P),] #filter NA P
  if(any(colnames(df)=="P") & is.finite(maxP)) df<-df[df$P<maxP,] #filter away low powered associations
  if(any(colnames(df)=="P") & is.finite(maxNLogP)) df$P<- shru::clipValues(df$P,min = 10^(-maxNLogP), max = NULL)
  
  cols<-c("SNP","BP","CHR")
  if(any(colnames(df)=="P")) cols<-c(cols,"P")
  if(any(colnames(df)=="BETA")) cols<-c(cols,"BETA")
  if(any(colnames(df)=="SE")) cols<-c(cols,"SE")
  
  df<-df[,cols]
  #df$CHR<-as.integer(df$CHR)
  
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
  
  if(!is.null(title)){
    plot<- plot + 
      ggtitle(title)
  }

  return(plot)
}


#stolen from https://gist.github.com/slowkow/9041570
plot.qq.custom <- function(ps,
                           ci = 0.95,
                           maxP=0.07,
                           title=NULL
                           ) {
  #ps<-head(as.data.frame(cSumstats),n=1000000)$P
  #ps<-cSumstats$P
  n0  <- length(ps)
  ps <- ps[ps<maxP]
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
  plot <- ggplot2::ggplot(df) +
    ggplot2::geom_point(aes(expected, observed),
                        #shape = 1,
                        size = 1.1, alpha=0.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    ggplot2::geom_line(aes(expected, cupper), linetype = 2) +
    ggplot2::geom_line(aes(expected, clower), linetype = 2) +
    ggplot2::xlab(log10Pe) +
    ggplot2::ylab(log10Po) +
    theme_bw()
    # theme( 
    #   legend.position="none",
    #   panel.border = element_blank(),
    #   panel.grid.major.x = element_blank(),
    #   panel.grid.minor.x = element_blank()
    # )
  
  if(!is.null(title)){
    plot<- plot + 
      ggtitle(title)
  }
  
  return(plot)
}

#plot for genetic correlations and similar, in matrix form
#requires the corrplot package - NOT INCLUDED IN SHRU
plot.corr <- function(
    corr, 
    pmat=NULL, 
    SE=NULL, 
    filename,
    addrect = NULL,
    is.corr = T,
    is.full = T,
    tl.cex = 1.0,
    number.cex = 1.0,
    number.digits=2,
    newnames=NULL,
    title=NULL,
    sig.level=c(0.05),
    par.col.lim.arg=NULL
    ){
  # corr=corrForPlot
  # pmat=pForPlot
  # SE=NULL
  # filename = file.path(p$folderpath.plots,"covG.large.png")
  # addrect = NULL
  # is.corr = T
  # is.full = T
  # tl.cex = 1.0
  # number.cex = 1.0
  # number.digits=2
  # newnames=NULL
  # title=NULL
  # sig.level=c(0.05)
  # par.col.lim.arg=NULL
  # is.corr = F
  # is.full = F
  # newnames = p$sumstats.sel.ssimp$code.nice
  # par.col.lim.arg = c(-0.9,0.9)
  
  
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
  m.type<-"full"
  if(!is.full) m.type <- "upper"
  if(!is.corr) par.addCoef.col <- 'black'
  if(is.null(par.col.lim.arg)){
    if(is.corr) par.col.lim <- c(-1,1)
    if(!is.corr) par.col.lim <- NULL #c(min(corr), max(corr))
    #if(!is.corr & all(corr>=0)) par.col.lim <- c(-max(corr),max(corr))
    #if(!is.corr & all(corr<=0)) par.col.lim <- c(min(corr),-min(corr))
    if(all(corr==0)) par.col.lim <- c(-1,1) #fallback if all corr are 0
  } else {
    par.col.lim<-par.col.lim.arg
    if(min(corr) < par.col.lim[1]) par.col.lim[1]<-min(corr)
    if(max(corr) > par.col.lim[2]) par.col.lim[2]<-max(corr)
    #if(all(corr>=0)) par.col.lim <- c(-max(par.col.lim),max(par.col.lim))
    #if(all(corr<=0)) par.col.lim <- c(min(par.col.lim),-min(par.col.lim))
  }
    
  #print(par.col.lim)
  
  pal<-colorRampPalette(c("#FF6666","#EEEEEE","#6666FF"))
  if(all(corr>=0)){
    pal<-colorRampPalette(c("#EEEEEE","#6666FF"))
  } else if(all(corr<=0)){
    pal<-colorRampPalette(c("#FF6666","#EEEEEE"))
  }
  
  # if(!is.corr & all(corr>=0)) par.col.lim <- c(-round(max(corr),digits = 0),round(max(corr),digits = 0))
  # if(!is.corr & all(corr<=0)) par.col.lim <- c(round(min(corr)-1,digits = 0),-1*(round(min(corr)-1)))
  
  png(filename = filename, width = 9, height = 9, units = 'in', res = 300, family = "Helvetica")
  par(xpd=TRUE) #keep labels inside margins
  
  corrplot(corr = corr, type = m.type, method = "circle", insig='blank', order = "original", p.mat = pmat, sig.level = sig.level, pch.cex = 1.5, full_col=TRUE, na.label = "square", na.label.col = "grey30", diag=T, addrect = addrect, is.corr = is.corr, tl.cex = tl.cex, number.cex = number.cex, number.digits = number.digits, addCoef.col = par.addCoef.col, full_col=(!is.corr), col = pal(200), col.lim = par.col.lim, mar=c(0,0,4,0), title = title)$corrPos -> p1  #COL2('RdBu')
  
  #full_col=!is.corr
  
  #p1$corr <- ifelse(0==-p1$y + ncol(corr)-p1$x + 1, NA_real_, p1$corr)
  if(is.corr) text(p1$x, p1$y, round(p1$corr, number.digits)) #only full coeficients for correlation plots
  
  dev.off()
  
}

# #TESTS
# mvldsc = pretendGCTAcovstruct
# folderpath.plots = p$folderpath.plots
# code = "gladp.gcta.reml"
# titleTemplate = "GCTA REML, GLAD+ reference"
# titleAddition = ""
#mvldscComparison =semplate$amendGsemLDSC(  p$mvLD$covstruct.GSEMmvLDSC.gladp.sim)
#titleAdditionComparison = " vs original Genomic SEM LDSC"
#df.summary=batRes$df.summary
#newnames = newnames

# mvldscComparison=NULL
# testOnlyTraitNameCodes=NULL
# doPlotting=TRUE
# code="covgTest"
# folderpath.plots=""
# newnames=NULL
# titleAddition=""
# titleAdditionComparison=""
# titleTemplate=""
# df.summary=NULL

#routine to test and generate plots for LDSC++ or Genomic SEM multivariate LDSC results
plotAndTestBatteryForMVLDSC <- function(
    mvldsc,
    mvldscComparison=NULL,
    testOnlyTraitNameCodes=NULL,
    doPlotting=TRUE,
    code="covgTest",
    folderpath.plots="",
    newnames=NULL,
    titleAddition="",
    titleAdditionComparison="",
    titleTemplate="",
    df.summary=NULL
    ){
  
  if(is.null(df.summary)) df.summary <- data.frame() #we can append rows to this
  
  traitCodes<-colnames(mvldsc$S)
  
  if(is.null(testOnlyTraitNameCodes)) testOnlyTraitNameCodes<-traitCodes
  
  cv.S_StandLimit <- 0.05 #0.05, or 0.1, 0.2 were suggested
  
  par.col.lim.arg.sym = c(-0.8,0.8)
  par.col.lim.arg.cv = c(0,0.8)
  
  effectiveNumberOfTests<-shru::getEffectiveNumberOfTests(covarianceMatrix = mvldsc$V_Stand)
  
  #using explicitly liability scale p-values
  mvldsc$cov.p.fdr2<-matrix(NA,nrow = nrow(mvldsc$cov.p.liab),ncol = ncol(mvldsc$cov.p.liab))
  mvldsc$cov.p.fdr2[lower.tri(mvldsc$cov.p.liab,diag = T)]<-p.adjust2(mvldsc$cov.p.liab[lower.tri(mvldsc$cov.p.liab,diag = T)], method = "fdr",n = effectiveNumberOfTests)
  mvldsc$cov.p.fdr2[upper.tri(mvldsc$cov.p.fdr2,diag = F)]<-t( mvldsc$cov.p.fdr2)[upper.tri(mvldsc$cov.p.fdr2,diag = F)]
  
  if(!is.null(mvldscComparison)){
    mvldscComparison$cov.p.fdr2<-matrix(NA,nrow = nrow(mvldscComparison$cov.p.liab),ncol = ncol(mvldscComparison$cov.p.liab))
    mvldscComparison$cov.p.fdr2[lower.tri(mvldscComparison$cov.p.liab,diag = T)]<-p.adjust2(mvldscComparison$cov.p.liab[lower.tri(mvldscComparison$cov.p.liab,diag = T)], method = "fdr",n = effectiveNumberOfTests)
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
  
  if(doPlotting) plot.corr(
    corr = clipValues(mvldsc$I,-3,3),
    pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("i.",code,".png")),
    is.corr = F,
    is.full = F,
    par.col.lim.arg = par.col.lim.arg.sym,
    number.digits = 3,
    newnames = newnames,
    title = paste0("LD score regression intercepts, ",titleTemplate,titleAddition)
  )
  
  df.summary[code,"mAbsCovG"]<- mean(abs(mvldsc$S[lower.tri(mvldsc$S, diag = T)]),na.rm=T)
  df.summary[code,"mAbsCovG.h2"]<- mean(abs(diag(mvldsc$S)),na.rm=T)
  if(doPlotting) plot.corr(
    corr = mvldsc$S,
    pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("covg.",code,".png")),
    is.corr = F,
    is.full = F,
    par.col.lim.arg = par.col.lim.arg.sym,
    number.digits = 3,
    newnames = newnames,
    title = paste0("Genetic covariances (covG), ",titleTemplate,titleAddition)
  )
  
  if(doPlotting) plot.corr(
    corr = mvldsc$S.SE,
    #pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("covgse.",code,".png")),
    is.corr = F,
    is.full = F,
    par.col.lim.arg = par.col.lim.arg.cv,
    number.digits = 4,
    newnames = newnames,
    title = paste0("Genetic covariances (covG) S.E., ",titleTemplate,titleAddition)
  )
  
  #use covG but apply rG limit settings - turns out brings same results as using the rG straight away
  #TODO - add to ldscpp
  S_ForCV<-abs(mvldsc$S)
  g_se <- sqrt(diag(S_ForCV))
  S_se2 <- g_se %*% t(g_se) #this is the diagonal based genetic covariance
  S_ForCV[S_ForCV < S_se2 * cv.S_StandLimit] <- (S_se2 * cv.S_StandLimit)[S_ForCV < S_se2 * cv.S_StandLimit] #cap at min rg 0.05 to avoid extreme values
  S.CV<-(mvldsc$S.SE/S_ForCV)
  df.summary[code,"mAbsCovGCV"]<- mean(abs(S.CV[lower.tri(S.CV, diag = T)]),na.rm=T)
  df.summary[code,"mAbsCovGCV.h2"]<- mean(abs(diag(S.CV)),na.rm=T)
  if(doPlotting) plot.corr(
    corr = S.CV,
    #pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("covgse.std.",code,".png")),
    is.corr = F,
    is.full = F,
    par.col.lim.arg = par.col.lim.arg.cv,
    number.digits = 4,
    newnames = newnames,
    title = paste0("Genetic covariance (covG) Coefficients of Variation (CV), ",titleTemplate,titleAddition)
  )
  
  if(doPlotting) plot.corr(
    corr = clipValues(mvldsc$S_Stand,-1,1),
    pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("rg.",code,".png")),
    is.corr = T,
    is.full = F,
    #par.col.lim.arg = par.col.lim.arg.sym,
    number.digits = 3,
    newnames = newnames,
    title = paste0("Genetic correlations (rG), ",titleTemplate,titleAddition)
  )
  
  if(doPlotting) plot.corr(
    corr = mvldsc$S_Stand.SE,
    #pmat = mvldsc$cov.p.fdr2,
    filename = file.path(folderpath.plots,paste0("rgse.",code,".png")),
    is.corr = F,
    is.full = F,
    par.col.lim.arg = par.col.lim.arg.cv,
    number.digits = 4,
    newnames = newnames,
    title = paste0("Genetic correlations (rG) S.E., ",titleTemplate,titleAddition)
  )
  
  # mvldsc<-p$mvLD$covstruct.mvLDSC.1kg.vbcs.drc.winfo
  # mvldscComparison<-p$mvLD$covstruct.mvLDSC.1kg.original
  # testOnlyTraitNameCodes<-infoTraits
  
  if(is.null(testOnlyTraitNameCodes)){
    testOnlyTraitNameCodes<-colnames(mvldsc$S)
  }
  
  
  comparisonTest <- NULL
  if(!is.null(mvldscComparison)){
    cat("\n***Tests for ",code," ***\n")
    
    #Our tests, per trait combination
    comparisonTest <- difftest.matrix(
      mValues1 = mvldsc$S,
      mStandard_errors1 = mvldsc$S.SE,
      df1 = mvldsc$cov.blocks,
      #mStandard_errors1.std = mvldsc$S.SE.std,
      mValues2 = mvldscComparison$S,
      mStandard_errors2 = mvldscComparison$S.SE,
      #mStandard_errors2.std = mvldscComparison$S.SE.std,
      df2 = mvldscComparison$cov.blocks,
      symmetric = T,
      effectiveNumberOfTests = effectiveNumberOfTests
      )
    
    
    
    #Supporting non-parametric test: Paired Samples Wilcoxon signed-rank test - performed across all trait combinations
    comparisonTest.npar.values = wilcox.test(mvldsc$S[lower.tri(mvldsc$S,diag = T)], mvldscComparison$S[lower.tri(mvldscComparison$S,diag = T)], paired = TRUE)
    comparisonTest.npar.standard_errors = wilcox.test((mvldsc$S.SE^2)[lower.tri(mvldsc$S.SE,diag=T)], (mvldscComparison$S.SE^2)[lower.tri(mvldscComparison$S.SE,diag=T)], paired = TRUE)
    
    #Comparison plots
    dI<-mvldsc$I-mvldscComparison$I #re-compute as we need the full matrices here
    dI.totest<-dI[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    df.summary[code,"mDeltaI"]<- mean(dI.totest[lower.tri(dI.totest, diag = T)],na.rm=T)
    df.summary[code,"mDeltaI.h2"]<- mean(diag(dI.totest),na.rm=T)
    if(doPlotting) plot.corr(
      corr = dI,
      #pmat = pTest.fdr2,
      filename = file.path(folderpath.plots,paste0("di.",code,".png")),
      is.corr = F,
      is.full = F,
      par.col.lim.arg = par.col.lim.arg.sym,
      number.digits = 3,
      newnames = newnames,
      title = paste0("Diff. in LD score regression intercepts, ",titleTemplate,titleAddition,titleAdditionComparison,"\n[mean=",round(df.summary[code,"mDeltaI"],digits = 3),"]"," [mean(h2 only)=",round(df.summary[code,"mDeltaI.h2"],digits = 3),"]")
    )
    
    dP<-comparisonTest$pTest.values.adj
    dP[is.na(dP)]<-1
    dS_Stand<-clipValues(mvldsc$S_Stand,-1,1)-clipValues(mvldscComparison$S_Stand,-1,1) #re-compute as we need the full matrices here
    dS_Stand.totest<-dS_Stand[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    if(doPlotting) plot.corr(
      corr = dS_Stand,
      #pmat = dP,
      filename = file.path(folderpath.plots,paste0("drg.",code,".png")),
      is.corr = F,
      is.full = F,
      par.col.lim.arg = par.col.lim.arg.sym,
      number.digits = 3,
      newnames = newnames,
      title = paste0("Diff. in rG, ",titleTemplate,titleAddition,titleAdditionComparison,"\n[abs(mean)=",round(mean(abs(dS_Stand.totest[lower.tri(dS_Stand.totest, diag = T)])),digits = 2),"] p(W.T)=",round(comparisonTest.npar.values$p.value,digits = 3))
    )
    
    #diff. plain covG - only for summary tables
    dS<-mvldsc$S-mvldscComparison$S #re-compute as we need the full matrices here
    dS[is.na(dS)]<-0 #in case NA
    dS.totest<-dS[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    df.summary[code,"mDeltaCovG"]<- mean(dS.totest[lower.tri(dS.totest, diag = T)],na.rm=T)
    df.summary[code,"mDeltaCovG.h2"]<- mean(diag(dS.totest),na.rm=T)
    # if(doPlotting) plot.corr(
    #   corr = dS_rel,
    #   #pmat = pTest.fdr2,
    #   filename = file.path(folderpath.plots,paste0("dcovg.",code,".png")),
    #   is.corr = F,
    #   is.full = F,
    #   par.col.lim.arg = c(-80,80),
    #   number.digits = 1,
    #   newnames = newnames,
    #   title = paste0("Diff. in covG, ",titleTemplate,titleAddition,titleAdditionComparison,"\n[mean=",round(df.summary[code,"mDeltaCovG"],digits = 1),"]", " [mean(h2 only)=",round(df.summary[code,"mDeltaCovG.h2"],digits = 1),"] p(W.T)=",round(comparisonTest.npar.values$p.value,digits = 3))
    # )
    
    
    #relative, in percent
    dS_rel<-100*(mvldsc$S-mvldscComparison$S)/abs(mvldscComparison$S) #re-compute as we need the full matrices here
    dS_rel[is.na(dS_rel)]<-0 #in case NA
    dS_rel.totest<-dS_rel[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    df.summary[code,"mRelDeltaCovG"]<- mean(dS_rel.totest[lower.tri(dS_rel.totest, diag = T)],na.rm=T)
    df.summary[code,"mRelDeltaCovG.h2"]<- mean(diag(dS_rel.totest),na.rm=T)
    df.summary[code,"maxRelSigDeltaCovG"] <- max(dS_rel.totest[comparisonTest$pTest.values.adj<=0.05 & is.finite(dS_rel.totest)],na.rm=T)
    df.summary[code,"minRelSigDeltaCovG"] <- min(dS_rel.totest[comparisonTest$pTest.values.adj<=0.05 & is.finite(dS_rel.totest)],na.rm=T)
    if(doPlotting) plot.corr(
      corr = dS_rel,
      #pmat = pTest.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovg.rel.",code,".png")),
      is.corr = F,
      is.full = F,
      par.col.lim.arg = c(-80,80),
      number.digits = 1,
      newnames = newnames,
      title = paste0("Rel. (%) diff. in covG, ",titleTemplate,titleAddition,titleAdditionComparison,"\n[mean=",round(df.summary[code,"mRelDeltaCovG"],digits = 1),"]", " [mean(h2 only)=",round(df.summary[code,"mRelDeltaCovG.h2"],digits = 1),"] p(W.T)=",round(comparisonTest.npar.values$p.value,digits = 3))
    )
    
    #covG CV, standardised
    comparison_S_ForCV<-abs(mvldscComparison$S)
    comparison_S_ForCV[is.na(comparison_S_ForCV)]<-0 #in case NA
    comparison_g_se <- sqrt(diag(comparison_S_ForCV))
    comparison_S_se2 <- comparison_g_se %*% t(comparison_g_se)
    comparison_S_ForCV[comparison_S_ForCV < comparison_S_se2 * cv.S_StandLimit] <- (comparison_S_se2 * cv.S_StandLimit)[comparison_S_ForCV < comparison_S_se2 * cv.S_StandLimit] #cap at min rg 0.05 to avoid extreme values
    dCV_cov <- ((mvldsc$S.SE/S_ForCV) - (mvldscComparison$S.SE/comparison_S_ForCV))
    dCV_cov[is.na(dCV_cov)]<-0 #in case NA
    dCV_cov.totest <- dCV_cov[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    df.summary[code,"mDeltaCovGCV"]<- mean(dCV_cov.totest[lower.tri(dCV_cov.totest, diag = T)],na.rm=T)
    df.summary[code,"mDeltaCovGCV.h2"]<- mean(diag(dCV_cov.totest),na.rm=T)
    if(doPlotting) plot.corr(
      corr = dCV_cov,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovgse.std.",code,".png")),
      is.corr = F,
      is.full = F,
      par.col.lim.arg = par.col.lim.arg.sym,
      number.digits = 3,
      newnames = newnames,
      title = paste0("Diff. in covG Coefficient of Variation (CV), ",titleTemplate,titleAddition,titleAdditionComparison,"\n[mean=",round(df.summary[code,"mDeltaCovGCV"],digits = 3),"]"," [mean(h2 only)=",round(df.summary[code,"mDeltaCovGCV.h2"],digits = 3),"] p(W.T)=",round(comparisonTest.npar.standard_errors$p.value,digits = 3))
    )
    
    #cov CV, standardised, relative
    nUniqueNegative <- sum(dCV_cov[lower.tri(dCV_cov,diag = T)]<0,na.rm = T)
    dCV_cov_rel <- 100*(dCV_cov/(mvldscComparison$S.SE/comparison_S_ForCV))
    dCV_cov_rel[is.na(dCV_cov_rel)]<-0 #in case NA
    dCV_cov.rel.totest <- dCV_cov_rel[testOnlyTraitNameCodes,testOnlyTraitNameCodes]
    df.summary[code,"mRelDeltaCovGCV"]<- mean(dCV_cov.rel.totest[lower.tri(dCV_cov.rel.totest, diag = T)],na.rm=T)
    df.summary[code,"mRelDeltaCovGCV.h2"]<- mean(diag(dCV_cov.rel.totest),na.rm=T)
    df.summary[code,"maxRelSigDeltaCovGCV"] <- max(dCV_cov.rel.totest[comparisonTest$pTest.standard_errors.adj<=0.05 & is.finite(dCV_cov.rel.totest)],na.rm=T)
    df.summary[code,"minRelSigDeltaCovGCV"] <- min(dCV_cov.rel.totest[comparisonTest$pTest.standard_errors.adj<=0.05 & is.finite(dCV_cov.rel.totest)],na.rm=T)
    if(doPlotting) plot.corr(
      corr = dCV_cov_rel,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovgse.std.rel.",code,".png")),
      is.corr = F,
      is.full = F,
      par.col.lim.arg = c(-80,80),
      number.digits = 2,
      newnames = newnames,
      title = paste0("Rel. (%) diff. in covG Coefficient of Variation (CV), ",titleTemplate,titleAddition,titleAdditionComparison,"\n[mean=",round(mean(dCV_cov.rel.totest[lower.tri(dCV_cov.rel.totest, diag = T)]),digits = 2),"]"," [mean(h2 only)=",round(df.summary[code,"mRelDeltaCovGCV.h2"],digits = 2),"] [#neg=",nUniqueNegative,"] p(W.T)=",round(comparisonTest.npar.standard_errors$p.value,digits = 3))
    )
    
   
    #p-values for the covG difference
    comparisonTest$pTest.values.adj[is.na(comparisonTest$pTest.values.adj)]<-1 #in case NA
    nUniqueSigDeltaCovG <- sum(comparisonTest$pTest.values.adj[lower.tri(comparisonTest$pTest.values.adj,diag = T)]<0.05,na.rm = T)
    if(doPlotting) plot.corr(
      corr = comparisonTest$pTest.values.adj,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovg.p.",code,".png")),
      is.corr = F,
      is.full = F,
      par.col.lim.arg = par.col.lim.arg.cv,
      number.digits = 3,
      newnames = newnames,
      title = paste0("Sig. level of diff. in covG (FDR corrected, ",effectiveNumberOfTests," effective tests), p(W.T)=",round(comparisonTest.npar.values$p.value,digits = 3),"\n",titleTemplate,titleAddition,titleAdditionComparison)
    )
    
    #p-values for the covG S.E. difference
    comparisonTest$pTest.standard_errors.adj[is.na(comparisonTest$pTest.standard_errors.adj)]<-1 #in case NA
    nUniqueSigDeltaCovGCV <- sum(comparisonTest$pTest.standard_errors.adj[lower.tri(comparisonTest$pTest.standard_errors.adj,diag = T)]<0.05,na.rm = T)
    if(doPlotting) plot.corr(
      corr = comparisonTest$pTest.standard_errors.adj,
      #pmat = covse.p.diff.fdr2,
      filename = file.path(folderpath.plots,paste0("dcovgse.p.",code,".png")),
      is.corr = F,
      is.full = F,
      par.col.lim.arg = par.col.lim.arg.cv,
      number.digits = 3,
      newnames = newnames,
      title = paste0("Sig. level of diff. in Var(covG) significance (FDR corrected, ",effectiveNumberOfTests," effective tests), p(W.T)=",round(comparisonTest.npar.standard_errors$p.value,digits = 3),"\n",titleTemplate,titleAddition,titleAdditionComparison)
    )
    
    #additional summary
    df.summary[code,"numNegDeltaCovGCV"]<-nUniqueNegative
    df.summary[code,"numSigDeltaCovG"]<-nUniqueSigDeltaCovG
    df.summary[code,"numSigDeltaCovGCV"]<-nUniqueSigDeltaCovGCV
    df.summary[code,"pWTCovG"]<-comparisonTest.npar.values$p.value
    df.summary[code,"pWTCovGVar"]<-comparisonTest.npar.standard_errors$p.value
    
  }
  
  
  return(list(
    effectiveNumberOfTests=effectiveNumberOfTests,
    cov.p.fdr2=mvldsc$cov.p.fdr2,
    S.CV=S.CV,
    comparisonTest=comparisonTest,
    df.summary=df.summary
    ))
}



