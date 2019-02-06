#!/usr/bin/env Rscript

# ----------------------------------------------
# Download and read in prokaryotic genome report
# ----------------------------------------------

cat("=> Fetching prokaryotic genome report from NCBI ...\n")
if (file.exists("prokaryotes.txt")){
  file.remove("prokaryotes.txt")
}
system("wget -q ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt")
cat("Done!\n")

cat("=> Reading prokaryotes.txt ...\n")
genome_report <- read.table(file="prokaryotes.txt", header=FALSE, fill = TRUE, stringsAsFactors = FALSE,
                            sep = "\t", row.names = NULL, quote = "")
col_names <- c("Organism_Name", "TaxID", "BioProject_Accession", "BioProject_ID", "Group",	
               "SubGroup", "Size_in_Mb", "GC_Content", "Replicons", "WGS",	"Scaffolds", 
               "Genes", "Proteins", "Release_Date", "Modify Date", "Status",	"Center", 
               "BioSample_Accession", "Assembly_Accession", "Reference", "FTP Path", 
               "Pubmed_ID", "Strain")
colnames(genome_report) <- col_names
cat("Done!\n")

# ----------------------------------------------
# Download and read in prokaryotic genome report
# ----------------------------------------------

reanno<-function(x){
	temp<-unlist(strsplit(as.character(x["Replicons"]),split=";"))
	temp<-gsub(".*:","",temp)
	accession<-c()
	refseq<-rep(NA,length(temp))
	
	for(i in 1:length(temp)){
		temp2<-unlist(strsplit(temp[i],split="/"))
		if(length(temp2)==2){
			accession<-c(accession,temp2[2])
			refseq[i]<-temp2[1]
			
		}
		if(length(temp2)==1){
			accession<-c(accession,temp2[1])
			
		}
	}
	out<-cbind(accession,refseq)
	out
}

no<-which(genome_report$Replicons=="-")
no1<-which(genome_report$Replicons=="")

no<-c(no,no1)
genome_report<-genome_report[-no,]
ref<-apply(genome_report,1,reanno)
ref<-do.call(rbind,ref)

save(ref, file = "accession_to_refseq")
cat("Done!\n")
