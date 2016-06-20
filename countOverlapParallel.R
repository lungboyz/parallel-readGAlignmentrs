#make the countOverlap run in parallel
library(parallel)
library(GenomicAlignments)



# the function to generate a parral version of readGAlignments()
#this accepts a TxDB that has been processed by eByg<-exonsBy(txdb, by="gene")
mkWorker <- function(eByg) {
  # make sure each of the three values we need passed 
  # are available in this environment
  library(GenomicAlignments)
  
  force(eByg)
  # define any and every function our worker function 
  # needs in this environment
  countOneBam <- function(eByg,aBam) {
    aligns <- GenomicAlignments::readGAlignments(aBam, index=aBam, use.names=T) # Substitute next two lines with this
    #print(aBam)
    GenomicRanges::countOverlaps(eByg, aligns)
  }
  # Finally: define and return our worker function.
  # The function worker's "lexical closure" 
  # (where it looks for unbound variables)
  # is mkWorker's activation/execution environment 
  # and not the usual Global environment.
  # The parallel library is willing to transport 
  # this environment (which it does not
  # do for the Global environment).
  worker <- function(aBam) {
    countOneBam(eByg,aBam)
  }
  return(worker)
}
df<-as.data.frame(matrix(nrow=56623))
#do in chunks of BAM files
b<-length(mybam)
#the first record
a<-1
#the last record
b<-length(mybam)
while(a<b){
#open the parrallel cluster  
  parallelCluster <- parallel::makeCluster(parallel::detectCores())
  print(parallelCluster)
  print(a)
#first bam record to read in chunk
  s<-a
#Last bam record to read in chunk  
  e<- a+6
  if(e>b){
    e<-b
  }
  #set a for the next start elelment
  a<-e+1
    
  countedBams <- parallel::parLapply(parallelCluster,mybam[s:e], mkWorker(eByg))

  df<-data.frame(df,(do.call(cbind, lapply(countedBams, data.frame, stringsAsFactors=FALSE))))
#close the cluster to release the memory allocation as
  #the memmory seems to become slowly consumed
    stopCluster(parallelCluster)
}

colnames(df)<-names(countedBams)
head(df)


