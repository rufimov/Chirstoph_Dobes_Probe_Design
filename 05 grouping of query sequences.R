# 05 grouping of query sequences.R
# first/last update: 6.11.2019/22.11.2019
# the script reads in the result of a nucleotide blast search performed at NCBI
# mRNAs selected for bait design were either blasted against themselves or against the full set of Malus mRNAs
# the script identifies paralogous relationships among mRNAs and checks whether mRNAs which are member of a paralogous pair had also two hits in the blast search against the Malus genome
# based on the gained information mRNAs and extracted Malus and Pyrus sequences will be grouped and aligned in script ##

# names of input files
tsv                      <- "mRNAs blasted against each other.tsv"       # name of tsv-file
filtered_query_sequences <- "Pyrus_hits_loci_filtered_6_loci_with_chro00_02122019_queries_fully_compatible_both_genomes.fasta"  # name of file with filtered sequences for probe design; note that the sequence name and the sequence must be tab delimited and named "seq_no" and "sequence"; name without ">"
file_summary_statistics  <- "Pyrus_Malus_detailed_statistics_fully_compatible_queries_02122019.txt"                             # name of detailed summary statistics for the mRNAs fully compatible in both genomes

# name of output files
name_outputfile_data     <- "relationships among mRNAs 02122019"


##########################################################################################################
########################################### read in data #################################################
##########################################################################################################

query_sequences                  <- read.csv(filtered_query_sequences, sep="\t", stringsAsFactors = F)
query_sequences[,"query_length"] <- apply(query_sequences,1,nchar)[2,]
N_query_sequences                <- nrow(query_sequences)

results_blast_search             <- read.csv(tsv, sep="\t", dec=".", stringsAsFactors = F)
results_blast_search             <- subset(results_blast_search, query_acc %in% query_sequences[,1]) # restrict lines to query sequences
results_blast_search             <- results_blast_search[with(results_blast_search, order(query_acc, subject_acc, subject_start)), ]

summary_statistics               <- read.csv(file_summary_statistics, sep="\t", dec=",", stringsAsFactors = F)


##########################################################################################################
######################################### search for paralogs ############################################
##########################################################################################################

query_acc_with_paralogs           <- aggregate(results_blast_search[1:2], list(results_blast_search$query_acc), length)[,c(1,3)]
colnames(query_acc_with_paralogs) <- c("query_acc","N_hits")
query_acc_with_paralogs           <- subset(query_acc_with_paralogs, N_hits > 1)[,"query_acc"]

query_accs_and_their_paralogs     <- subset(results_blast_search, query_acc %in% query_acc_with_paralogs)[,c("query_acc","subject_acc")]

rows <- NULL
for(i in 1:nrow(query_accs_and_their_paralogs)) { # keep only rows for which query_acc != subject_acc
    row_i <- rownames(subset(query_accs_and_their_paralogs[i,], query_acc != query_accs_and_their_paralogs[i,2]))
    rows  <- c(rows, row_i)
}

query_accs_and_their_paralogs <- query_accs_and_their_paralogs[rows,]
query_accs_and_their_paralogs <- query_accs_and_their_paralogs[!duplicated(query_accs_and_their_paralogs), ]

if(sum(aggregate(query_accs_and_their_paralogs, list(query_accs_and_their_paralogs$query_acc), length)[,2]) != nrow(query_accs_and_their_paralogs)) {print("there are mRNAs with more than one paralogous mRNAs in the data")}

line <- NULL
for(i in 1:nrow(query_accs_and_their_paralogs)) { # sort accession numbers of paralogs columnwise (lower value in column 1)
    line_i <- sort(as.character(query_accs_and_their_paralogs[i,]))
    line  <- rbind(line, line_i)
}
rownames(line) <- 1:nrow(query_accs_and_their_paralogs)
colnames(line) <- c("seq1","seq2")

x <- query_accs_and_their_paralogs_unique <- line[!duplicated(line), ]
y <- line[duplicated(line), ]

if(nrow(x) != nrow(y)) {print("hits of paralogous mRNAs were NOT fully reciprocal")} else {print("hits of paralogous mRNAs were fully reciprocal")}


##########################################################################################################
########################### make statistics for the comparison of queries ################################
##########################################################################################################

x <- as.data.frame(cbind(sort(unique(results_blast_search[,"query_acc"])),1))
y <- subset(summary_statistics, genome == "Malus")[,c("query_ID","N_loci")]
y <- y[!duplicated(y), ]
y[,"N_mRNAs_in_data"]                                                                         <- 1 # number of loci expected from presence of paralogous mRNAs
y[,"N_mRNAs_in_data"][y[,"query_ID"] %in% as.character(query_accs_and_their_paralogs_unique)] <- 2
# subset(y, N_loci == 1) # queries which hat only locus in the Malus genome but for which a patalogous mRNA was found
N_loci_accessions              <- merge(y,query_accs_and_their_paralogs,by.x="query_ID",by.y="query_acc",all.x=T)
colnames(N_loci_accessions)[4] <- "ID_paralogous_query"

z                    <- summary_statistics[,c("genome","query_ID","locus_ID","N_loci","chromosome","query_start","query_end","subject_start","subject_end")]
z[,"length_subject_approx"] <- z[,"subject_end"]-z[,"subject_start"]
u                    <- aggregate(z[,"length_subject_approx"], list(z$genome,z$query_ID,z$N_loci,z$locus_ID,z$chromosome), sum)
colnames(u)          <- c("genome","query_ID","N_loci","locus_ID","chromosome","length_subject_approx") # that's only the approximate length
v                    <- aggregate(z, list(z$genome,z$query_ID,z$N_loci,z$locus_ID,z$chromosome), min)[,c("genome","query_ID","locus_ID","N_loci","chromosome","query_start")]
colnames(v)[6]       <- "query_start_minimum" # the lowest value observed
w                    <- aggregate(z, list(z$genome,z$query_ID,z$N_loci,z$locus_ID,z$chromosome), max)[,c("genome","query_ID","locus_ID","N_loci","chromosome","query_end")]
colnames(w)[6]       <- "query_end_maximum" # the highest value observed
x                    <- aggregate(z, list(z$genome,z$query_ID,z$N_loci,z$locus_ID,z$chromosome), min)[,c("genome","query_ID","locus_ID","N_loci","chromosome","subject_start")]
colnames(x)[6]       <- "subject_start_minimum" # the highest value observed
y                    <- aggregate(z, list(z$genome,z$query_ID,z$N_loci,z$locus_ID,z$chromosome), max)[,c("genome","query_ID","locus_ID","N_loci","chromosome","subject_end")]
colnames(y)[6]       <- "subject_end_maximum" # the highest value observed
x1                   <- aggregate(z, list(z$genome,z$query_ID,z$N_loci,z$locus_ID,z$chromosome), max)[,c("genome","query_ID","locus_ID","N_loci","chromosome","subject_start")]
colnames(x1)[6]      <- "subject_start_maximum" # the highest value observed
y1                   <- aggregate(z, list(z$genome,z$query_ID,z$N_loci,z$locus_ID,z$chromosome), min)[,c("genome","query_ID","locus_ID","N_loci","chromosome","subject_end")]
colnames(y1)[6]      <- "subject_end_minimum" # the highest value observed
start_end_Pyrus_Malus_loci <- merge(v,w)
start_end_Pyrus_Malus_loci <- merge(start_end_Pyrus_Malus_loci,x)
start_end_Pyrus_Malus_loci <- merge(start_end_Pyrus_Malus_loci,y)
start_end_Pyrus_Malus_loci <- merge(start_end_Pyrus_Malus_loci,x1)
start_end_Pyrus_Malus_loci <- merge(start_end_Pyrus_Malus_loci,y1)
start_end_Pyrus_Malus_loci <- merge(start_end_Pyrus_Malus_loci,query_sequences[,c("seq_no","query_length")],by.x="query_ID",by.y="seq_no")
start_end_Pyrus_Malus_loci <- merge(start_end_Pyrus_Malus_loci,u)
start_end_Pyrus_Malus_loci[,c("subject_start","subject_end","cut_subject_start","cut_subject_end")] <- 0
start_end_Pyrus_Malus_loci[,"subject_start"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0] <- start_end_Pyrus_Malus_loci[,"subject_start_minimum"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0]
start_end_Pyrus_Malus_loci[,"subject_end"]  [start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0] <- start_end_Pyrus_Malus_loci[,"subject_end_maximum"]  [start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0]
start_end_Pyrus_Malus_loci[,"subject_start"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0] <- start_end_Pyrus_Malus_loci[,"subject_start_maximum"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0]
start_end_Pyrus_Malus_loci[,"subject_end"]  [start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0] <- start_end_Pyrus_Malus_loci[,"subject_end_minimum"]  [start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0]
start_end_Pyrus_Malus_loci[,"cut_subject_start"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0] <- start_end_Pyrus_Malus_loci[,"subject_start"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0] - (start_end_Pyrus_Malus_loci[,"query_start_minimum"] [start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0] - 1)
start_end_Pyrus_Malus_loci[,"cut_subject_end"]  [start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0] <- start_end_Pyrus_Malus_loci[,"subject_end"]  [start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0] + (start_end_Pyrus_Malus_loci[,"query_length"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0] - start_end_Pyrus_Malus_loci[,"query_end_maximum"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] > 0])
start_end_Pyrus_Malus_loci[,"cut_subject_start"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0] <- start_end_Pyrus_Malus_loci[,"subject_end"]  [start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0] - (start_end_Pyrus_Malus_loci[,"query_length"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0] - start_end_Pyrus_Malus_loci[,"query_end_maximum"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0])
start_end_Pyrus_Malus_loci[,"cut_subject_end"]  [start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0] <- start_end_Pyrus_Malus_loci[,"subject_start"][start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0] + (start_end_Pyrus_Malus_loci[,"query_start_minimum"] [start_end_Pyrus_Malus_loci[,"length_subject_approx"] < 0] - 1)
start_end_Pyrus_Malus_loci[,"length_subject_real"] <- start_end_Pyrus_Malus_loci[,"cut_subject_end"] - start_end_Pyrus_Malus_loci[,"cut_subject_start"] + 1
start_end_Pyrus_Malus_loci <- start_end_Pyrus_Malus_loci[,c("genome","query_ID","N_loci","locus_ID","chromosome","query_start_minimum","query_end_maximum","query_length","subject_start","subject_start_minimum","cut_subject_start","subject_end","subject_end_maximum","cut_subject_end","length_subject_approx","length_subject_real")]

##########################################################################################################
################################################## save data #############################################
##########################################################################################################

write.table(query_accs_and_their_paralogs_unique,paste(name_outputfile_data,"accession numbers paralagous queries unique.txt",sep=" "),sep="\t",dec=",",row.names=F,quote=F)
write.table(N_loci_accessions                   ,paste(name_outputfile_data,"N loci queries observed and expected.txt"                ,sep=" "),sep="\t",dec=",",row.names=F,quote=F)
write.table(start_end_Pyrus_Malus_loci          ,paste(name_outputfile_data,"sequence coordinates.txt"                                ,sep=" "),sep="\t",dec=",",row.names=F,quote=F)
write.table(results_blast_search                ,paste(name_outputfile_data,"results blast search.txt"                                ,sep=" "),sep="\t",dec=",",row.names=F,quote=F)
source("06 loci for subject extraction.R") # makes input files for the extraction of subject sequences required by the Linux bash files loopMalus.sh and loopPyrus.sh
# source("07 make fasta files.R")
# source("08 sequence comparison.R")
