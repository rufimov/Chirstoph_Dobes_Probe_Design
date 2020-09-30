# control.R
# first/last update: 25.11.2019/26.11.2019
# reads in separated alignments and aligns them with the mRNAs and the original subject sequences
# first sequence in the alignment is the Malus mRNA, second the first subject sequence

# library(DECIPHER)

filtered_query_sequences <- "Pyrus_hits_loci_filtered_6_loci_with_chro00_02122019_queries_fully_compatible_both_genomes.fasta"  # name of file with filtered sequences for probe design; note that the sequence name and the sequence must be tab delimited and named "seq_no" and "sequence"; name without ">"

file_subject_sequences   <- "mRNAs subject sequences extracted from Malus and Pyrus.txt"

##########################################################################################################
########################################### read in data #################################################
##########################################################################################################

query_sequences   <- read.csv(filtered_query_sequences, sep="\t", stringsAsFactors = F)
subject_sequences <- read.csv(file_subject_sequences, sep="\t", stringsAsFactors = F)
subject_sequences <- subject_sequences[,c("query_ID","sequence")]
queries           <- sort(unique(subject_sequences[,"query_ID"]))

FileNames <- choose.files(".fa") # select sequences

len   <- nchar(FileNames)
query <- substring(FileNames,(nchar(getwd())+10),(nchar(getwd())+21))

for (i in 1:length(queries)) {
    query_name <- queries[i]
    pos1       <- rownames(subset(subject_sequences,query_ID==query_name))
    pos2       <- which(query==query_name)
    if(length(pos2)>0) {
        file1           <- as.data.frame(c(as.character(paste(">",query_name,sep="")),as.character(query_sequences[query_sequences$seq_no==query_name,"sequence"])))
        colnames(file1) <- "V1"
        file2           <- as.data.frame(c(as.character(paste(">",query_name,sep="")),as.character(subject_sequences[pos1[1],"sequence"])))
        colnames(file2) <- "V1"
        file_aligned    <- rbind(file1,file2)
        for(j in pos2) {
            file_aligned_j <- as.data.frame(read.csv(FileNames[j], sep="\t", stringsAsFactors = F, header=F))
            file_aligned   <- rbind(file_aligned,file_aligned_j)
        }
        FileName <- paste("control/aligned/",query_name,sep="")
        FileName <- paste(FileName,".fa",sep="")
        write.table(file_aligned,FileName,row.names=F,col.names=FALSE,quote=F)

        dna                  <- readDNAStringSet(FileName)
        aligned_dna          <- AlignSeqs(dna)
        aligned_adjusted_dna <- AdjustAlignment(aligned_dna)

        write.table(aligned_adjusted_dna,FileName,col.names=FALSE,quote=F)
        file     <- read.csv(FileName, sep=" ",  header = F, stringsAsFactors = F)
        file[,1] <- paste(">",file[,1],sep="")
        file     <- c(as.character(file[,1]),as.character(file[,2]))
        m        <- (length(file)/2)
        sort_key <- NULL
        for(j in 1:m) {
            sort_key_j <- c(j,j+m)
            sort_key   <- c(sort_key,sort_key_j)
        }
        file <- file[sort_key]
        write.table(file,FileName,row.names=FALSE,col.names=FALSE,quote=F)
    }
}

