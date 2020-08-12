 adnifam <- read.table("dissertation/data/ADNI/ADNI.fam",sep="")
 dim(adnifam)= 757 by 6
 ricfile <- read.csv("dissertation/data/ADNI/phenotypes.csv", header=T)
 
 
 # Benazir Rowe
 # recoding all
 
 
 #code to start R on cluster
 qsub -I -l ncpus=1,mem=5gb,cput=5:0:0 -l walltime=5:0:0 /bin/bash
 module load intel intelmpi R 
 R
 
 
 
 #################################################### nee bimfile to get allele namesq
 
 #read bim file since it containes allele names
 
bimfile <-read.table("/storage/nipm/kerimbae/ADNI/ADNI.bim", header = FALSE)
 #subset the part for a chromosome in question
 
for ( q in 1:21){
   
   allele <- subset(bimfile, bimfile$V1 == q, stringsAsFactors=FALSE)
   #see bimfile for current chr
   head(allele)
   
   ###############################################
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
     
     j3
     if(is.na(j3)) next 
     write(j3,file=paste("/storage/nipm/kerimbae/ADNI/input/chr",q,".txt", sep=""),append=TRUE)
   }
   
   
   close(conn)
   
   
}
 
 
 
 
 ######### recode phenotype 2-affected, 1- unaffected fam  file format
 
 pheno = read.table( "swdrecode/chr1.recode.pheno.txt")
 
 for (i in 1: dim(pheno)[1])
 {
   if (pheno[i,] == "1")
   {
     pheno[i,] = 0
   }
   
   if (pheno[i,] == "2")
   {
     pheno[i,] = 1
   }
   
 }
 
 write.table(pheno,"swdrecode/phenoswd.txt", row.names = F, col.names = F, quote = F)
 
 