# Recode the PLINK bimbam12 data into mean genotype format
# Code is not ran on local machines
 
 
#code to start R on cluster
qsub -I -l ncpus=1,mem=10gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
module load intel intelmpi R 
R
 
 
 
####################################################recode all chromosomes at once
# load bimfile to retreive alleles
bimfile <-read.table("/storage/nipm/kerimbae/ADNI/ADNI.bim", header = FALSE)

#read the bimbam12 genotype file
path = paste("/storage/nipm/kerimbae/ADNI/adni_allchr.recode.geno.txt" ,sep = "")
conn <- file(path,open="r")
lines <- readLines(conn)
 

# Assemble each line to in mgt format: rsID, minor, major, copies of minor allele for person1 etc.

for (i in 4:length(lines)){
   #get allele names from bimfile
   minor <- lapply(bimfile[i-3,5], as.character)
   major <- lapply(bimfile[i-3,6], as.character)
      
   tmp1 <-  paste(minor, major, sep = ",")
   tmp <- paste0(",",tmp1, ",")
      
   # insert allele names
      
   j0 <- sub("," , tmp, lines[i])
      
      
   #recode to 11 12 21 22 into 0,1,2 
   j <- gsub(",11",",0",j0)
   j1 <- gsub(",12",",1",j)
   j2 <- gsub(",22",",2",j1)
      
   #fill in missing genotypes with their mean (optional)
   m <-   as.factor(unlist(strsplit(j2, ",")))
   t <- table(m)
      
   sm <- round((t[3]+2*t[4])/sum(t[2:4]),2)
      
   j3 <- gsub("??",sm, j2, fixed = TRUE)
   
   #ADNI data had missing mean genotypes, drop them
   if(is.na(j3)) next 
   write(j3,file=paste("/storage/nipm/kerimbae/ADNI/adni_mgt_allchr.txt", sep=""),append=TRUE)
}
   
   
close(conn)
   
   










# recode each chromosome separately

for ( q in 1:22){
   
   allele <- subset(bimfile, bimfile$V1 == q, stringsAsFactors=FALSE)
   #see bimfile for current chr
   head(allele)
   
   # read genotype for single chromosome in bimbam12 format
   path = paste("/storage/nipm/kerimbae/ADNI/input/chr",q,".recode.geno.txt" ,sep = "")
   conn <- file(path,open="r")
   lines <- readLines(conn)
   
   for (i in 4:length(lines)){
     #get allele names
     minor <- lapply(allele[i-3,5], as.character)
     major <- lapply(allele[i-3,6], as.character)
     
     tmp1 <-  paste(minor, major, sep = ",")
     tmp <- paste0(",",tmp1, ",")
     
     # insert allele names
     j0 <- sub("," , tmp, lines[i])
     
     #recode to 0-2 
     j <- gsub(",11",",0",j0)
     j1 <- gsub(",12",",1",j)
     j2 <- gsub(",22",",2",j1)
     
     #fill in missing with their mean
     m <-   as.factor(unlist(strsplit(j2, ",")))
     t <- table(m)
     
     sm <- round((t[3]+2*t[4])/sum(t[2:4]),2)
     
     j3 <- gsub("??",sm, j2, fixed = TRUE)
     
     if(is.na(j3)) next 
     write(j3,file=paste("/storage/nipm/kerimbae/ADNI/input/chr",q,".txt", sep=""),append=TRUE)
   }
   
close(conn)

}
 