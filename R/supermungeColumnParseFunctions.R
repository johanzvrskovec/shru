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
  if(sum(grepl(pattern = "^rs\\w+:\\w+:\\w+", x= head(x = text, n=100000)))>90000){
    indexesLengths<-regexec(pattern = "^(rs\\w+):\\w+", text=text)
    matches<-regmatches(text,indexesLengths)
    return(unlist(lapply(X = matches, FUN = function(x) ifelse(is.na(x[2]),x[1],x[2]))))
  }
  
  #CHR:BP:A1:A2_rsXXXX;rsYYYY - format
  if(sum(grepl(pattern = "^\\w+:\\d+:\\w+:\\w+_rs", x= head(x = text, n=100000)))>90000){
    indexesLengths<-regexec(pattern = "^\\w+:\\d+:\\w+:\\w+_(.+)", text=text)
    matches<-regmatches(text,indexesLengths)
    prs <- unlist(lapply(X = matches, FUN = function(x) ifelse(is.na(x[2]),x[1],x[2])))
    prsSplit <- strsplit(prs, split = ";", fixed = T)
    return(unlist(lapply(X = prsSplit, FUN = function(x) x[1])))
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
