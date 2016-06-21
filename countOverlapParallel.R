#make the countOverlap run in parallel
library(parallel)
library(GenomicAlignments)



# the function to generate a parral version of readGAlignments()
#this accepts a TxDB that has been processed by eByg<-exonsBy(txdb, by="gene")
mkWorker <- function(eByg) {
  # make sure values we need passed 
  # are available in this environment
  library(GenomicAlignments)
  force(eByg)
  # define any and every function our worker function 
  # needs in this environment
  #this subfunction takes the eByg and a bam record 
  #to generate a counts record for bam overlaps against exons in the genome summed by gene
  countOneBam <- function(eByg,aBam) {
    aligns <- GenomicAlignments::readGAlignments(aBam, index=aBam, use.names=T) # Substitute next two lines with this
    #print(aBam)
    GenomicRanges::countOverlaps(eByg, aligns)
  }
  # Finally: define and return our worker function.
  # 
  worker <- function(aBam) {
    countOneBam(eByg,aBam)
  }
  return(worker)
}
#define th4e dataframe (this is the total numbner of genes that will be aligned against availible form the eByg file)
df<-as.data.frame(matrix(nrow=56623))
#the first record
a<-1
#the last record
#mybam is a named vector of the file locations for each BAM file
b<-length(mybam)
#a while statment will feed chunks of BAM records to the counter
while(a<=b){
#open the parrallel cluster  
  parallelCluster <- parallel::makeCluster(parallel::detectCores())
  #print(parallelCluster) 
  print(a) #prints the record chunk start position, optional
#first bam record to read in chunk
  s<-a
#Last bam record to read in chunk  
  e<- a+6 # the number can be adjusted to determine the total number of records perchunk. 
  #I have 8 cores so I send 7 BAMS and one core to co-ordinate, 
  #this is not as efficent as sending everyone all at once, but saves RAM
  
  #this if statment makes sure that the records list is not exceeded
  if(e>b){
    e<-b
  }
  #set a for the next start elelment
  a<-e+1
  #send a chunk of BAM records  
  countedBams <- parallel::parLapply(parallelCluster,mybam[s:e], mkWorker(eByg))
  #change the list to a data.frame and concatenate
  df<-data.frame(df,(do.call(cbind, lapply(countedBams, data.frame, stringsAsFactors=FALSE))))
#close the cluster to release the memory allocation as
  # earlier version sent all records but the memmory seems to become slowly consumed, 
  #leading to a crash after about 1000 BAM records were sent
    stopCluster(parallelCluster)
}

# the df does not have the samples names but they are in the order of the bam records, so if can be renamed.


