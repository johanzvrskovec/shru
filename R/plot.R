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