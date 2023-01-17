#LD Score Clustering - WORK IN PROGRESS!!
#Johan Zvrskovec, 2022

ldsccClass <- setRefClass("pgDatabaseUtility",
                                      fields = list(
                                        ds = "ANY",
                                        samplePrevalence = "ANY",
                                        populationPrevalence = "ANY",
                                        ancestrySetting = "ANY", #ancestry setting, for all datasets!
                                        ds.names = "ANY",
                                        numClusters = "numeric",
                                        cluster.M = "numeric",
                                        GC= "ANY", #vector with correction denominators
                                        filter.chisq.max = "numeric",
                                        filter.chisq.min = "numeric",
                                        filter.info = "numeric",
                                        filter.sinfo.imputed = "numeric", #filter for supermunge LD-IMP imputation quality score, for fully imputed pairwise variant combinations
                                        filter.sinfo = "numeric", #filter for supermunge LD-IMP imputation quality score, for partially imputed pairwise variant combinations
                                        filter.maf = "numeric",
                                        filter.mhc = "numeric", #can be either 37 or 38 for filtering the MHC region according to either grch37 or grch38
                                        filter.region.df = "ANY",
                                        M.mode = "character", #maf/real - maf: the size of the LD-score library with MAF>0.02, real: use the actual number of variants included in the trait/LD-score overlap
                                        N.mode = "character", #mean or variable
                                        leave.chr = "ANY", #vector with chromosome numbers to leave out of the analysis
                                        nThreads = "numeric",
                                        preweight = "logical", #perform pre-weighting (or not)
                                        resamplingMethod = "character" #jn/cv/bs(not implemented yet)
                                      ),
                                      methods = list
                                      (
                                        #this is the constructor as per convention
                                        initialize=function(
                                        ds,
                                        samplePrevalence,
                                        populationPrevalence,
                                        ancestrySetting = "ANY", #ancestry setting, for all datasets!
                                        ds.names = NULL,
                                        numClusters = NA_integer_,
                                        M = 20000,
                                        GC=NULL, #vector with correction denominators
                                        filter.chisq.max = NA_real_,
                                        filter.chisq.min = NA_real_,
                                        filter.info = .6,
                                        filter.sinfo.imputed = NA_real_, #filter for supermunge LD-IMP imputation quality score, for fully imputed pairwise variant combinations
                                        filter.sinfo = NA_real_, #filter for supermunge LD-IMP imputation quality score, for partially imputed pairwise variant combinations
                                        filter.maf=0.01,
                                        filter.mhc=NA_integer_, #can be either 37 or 38 for filtering the MHC region according to either grch37 or grch38
                                        filter.region.df=NULL,
                                        M.mode="maf", #maf/real - maf: the size of the LD-score library with MAF>0.02, real: use the actual number of variants included in the trait/LD-score overlap
                                        N.mode="variable", #mean or variable
                                        leave.chr=NULL, #vector with chromosome numbers to leave out of the analysis
                                        nThreads=5,
                                        preweight=T, #perform pre-weighting (or not)
                                        resamplingMethod="cv" #jn/cv/bs(not implemented yet)
                                      ) {
                                          
                                          ds <<- ds
                                          samplePrevalence <<- samplePrevalence
                                          populationPrevalence <<- populationPrevalence
                                          ancestrySetting <<- ancestrySetting
                                          ds.names <<- ds.names
                                          numClusters <<- numClusters
                                          cluster.M <<- cluster.M
                                          GC <<- GC
                                          filter.chisq.max <<- filter.chisq.max
                                          filter.chisq.min <<- filter.chisq.min
                                          filter.info <<- filter.info
                                          filter.sinfo.imputed <<- filter.sinfo.imputed
                                          filter.sinfo <<- filter.sinfo
                                          filter.maf <<- filter.maf
                                          filter.mhc <<- filter.mhc
                                          filter.region.df <<- filter.region.df
                                          M.mode <<- M.mode
                                          N.mode <<- N.mode
                                          leave.chr <<- leave.chr
                                          nThreads <<- nThreads
                                          preweight <<- preweight
                                          resamplingMethod <<- resamplingMethod
                                        }
                                      )
)

# we can add more methods after creating the ref class (but not more fields!)

ldsccClass$methods(
  setDs=function(ds,
                 samplePrevalence,
                 populationPrevalence
                 ){
    ds <<- ds
    samplePrevalence <<- samplePrevalence
    populationPrevalence <<- populationPrevalence
  }
)

#test
#l <- ldsccClass(ds = NULL, samplePrevalence = NULL, populationPrevalence = NULL)
