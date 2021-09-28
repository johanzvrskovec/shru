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
plot.manhattan.custom<-function(df){
  #df<-project$lfGWAS$gwas.for.display
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

plot <- ggplot(df, aes(x=BP.tot, y=-log10(Pval_Estimate))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", theme.color$contrastLight1), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdefinition$CHR, breaks= axisdefinition$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  geom_label_repel(data=subset(df, -log10(Pval_Estimate)>2.5), aes(label=SNP), size=3) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

  return(plot)
}