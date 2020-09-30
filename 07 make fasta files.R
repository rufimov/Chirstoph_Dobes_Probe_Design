# make fasta files.R
# first/last update 8.11.2019/8.11.2019
# reads in the subject sequences of Malus and Pyrus (made by the loopXY.sh scripts) and names the sequences
# creates unaligned and aligned fasta files for groups of sequences
# this is an subscript of the grouping of query sequences.R script; note that the script requires the "subject sequences"; and the Loci_for_subject_extraction.R needs to be run first

# library(DECIPHER)

##########################################################################################################
################################## make files with unaligned sequence ####################################
##########################################################################################################

Malus_sequs     <- read.csv("malus_subject_sequences.txt", sep="\t", stringsAsFactors = F, header=F)
Pyrus_sequs     <- read.csv("pyrus_subject_sequences.txt", sep="\t", stringsAsFactors = F, header=F)

# make DNA complement for inversly orientated subjects
Malus_nts               <- c("A","C","G","T")
Malus_nts_complementary <- c("T","G","C","A")
for(i in seqs_reversed_Malus) {
    seq        <- rev(strsplit(Malus_sequs[i,],NULL)[[1]])
    position_h <- list()
    for(h in 1:length(Malus_nts)) {
        position_h_i    <- which(seq==Malus_nts[h])
        position_h[[h]] <- position_h_i
    }
    for(h in 1:length(Malus_nts)) {
        seq[position_h[[h]]] <- Malus_nts_complementary[h]
    }
    Malus_sequs[i,] <- paste(seq,collapse='')
}

Pyrus_nts               <- c("A","C","G","T","a","c","g","t")
Pyrus_nts_complementary <- c("T","G","C","A","T","G","C","A")
for(i in seqs_reversed_Pyrus) {
    seq        <- rev(strsplit(Pyrus_sequs[i,],NULL)[[1]])
    position_h <- list()
    for(h in 1:length(Pyrus_nts)) {
        position_h_i    <- which(seq==Pyrus_nts[h])
        position_h[[h]] <- position_h_i
    }
    for(h in 1:length(Pyrus_nts)) {
        seq[position_h[[h]]] <- Pyrus_nts_complementary[h]
    }
    Pyrus_sequs[i,] <- paste(seq,collapse='')
}

sequs           <- rbind(Malus_sequs,Pyrus_sequs)
m               <- rbind(data_Malus[,c("query_ID","locus_ID","chromosome")],data_Pyrus[,c("query_ID","locus_ID","chromosome")])
m[,"genome"]                       <- "Pyrus"
m[1:nrow(data_Malus),"genome"]     <- "Malus"
m[1:nrow(data_Malus),"chromosome"] <- paste("chr",m[1:nrow(data_Malus),"chromosome"],sep="")
sequ_names      <- paste(m[,"query_ID"],"_",m[,"genome"],"_",m[,"chromosome"],"_locus",m[,"locus_ID"],sep="")
sequs           <- cbind(m[,"query_ID"],sequ_names,sequs)
colnames(sequs) <- c("query_ID","sequence_name","sequence")

write.table(sequs,"mRNAs subject sequences extracted from Malus and Pyrus.txt",sep="\t",row.names=F,quote=F)

queries                   <- query_sequences[,c("seq_no","sequence")]
colnames(queries)         <- c("query_ID","sequence")
queries[,"sequence_name"] <- paste(query_sequences[,"seq_no"],"Malus_mRNA",sep="_")
queries                   <- queries[,c("query_ID","sequence_name","sequence")]

sequs                     <- rbind(queries, sequs)

queries_with_paralogs    <- as.character(query_accs_and_their_paralogs_unique)
queries_with_no_paralogs <- subset(N_loci_accessions,!query_ID %in% queries_with_paralogs)[,"query_ID"]

for(i in 1:nrow(query_accs_and_their_paralogs_unique)) {
    file <- subset(sequs, query_ID %in% query_accs_and_their_paralogs_unique[i,])[,c("sequence_name","sequence")]
    write.table(file, paste("unaligned sequences/",query_accs_and_their_paralogs_unique[i,1],"_",query_accs_and_their_paralogs_unique[i,2],"_unaligned.seq",sep=""),sep="\t",row.names=F,quote=F)
}

for(j in queries_with_no_paralogs) {
    file <- subset(sequs, query_ID == j)[,c("sequence_name","sequence")]
    write.table(file, paste("unaligned sequences/",j,"_unaligned.seq",sep=""),sep="\t",row.names=F,quote=F)
}

##########################################################################################################
################################ make fasta files and align sequences ####################################
##########################################################################################################

FileNames <- choose.files(".seq") # select unaligned sequences
for (i in 1:length(FileNames)) {
    file                   <- read.csv(FileNames[i], sep="\t", stringsAsFactors = F)
    file[,"sequence_name"] <- paste(">",file[,"sequence_name"],sep="")
    file                   <- c(as.character(file[,"sequence_name"]),as.character(file[,"sequence"]))
    m                      <- (length(file)/2)
    sort_key               <- NULL
    for(j in 1:m) {
        sort_key_j <- c(j,j+m)
        sort_key   <- c(sort_key,sort_key_j)
    }
    file                   <- file[sort_key]

    len      <- nchar(FileNames[i])  
    FileName <- paste("aligned sequences/",substring(FileNames[i],(nchar(getwd())+22),len-14),"_aligned.fa",sep="")
    write.table(file,FileName,row.names=F,col.names=FALSE,quote=F)

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
