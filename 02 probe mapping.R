# probe_mapping.R
# first/last update: 9.4.2019/23.8.2019

##########################################################################################################
########################################### script desprition ############################################
##########################################################################################################

# The script first reads in a tsv-files which holds the output of a blast search (i.e. "hits" of "query sequences" against a "subject") created using the "Genome Database for Rosaceae": for instance a search of transiptome sequences from Malus (= "query sequences") against the assembled Malus genome (= "subject").
# The query sequences are also read in in fasta format.
# Queries suitable for probe design are selected in two steps (A and B):
# step A:
# (A1) Hits with a sequence identity relative to the subject lower than a preset value are excluded. The chosen value must be lower than the minimum sequences idenity allowing hybridisation between "probes" (an oligonuleotide used as a bait in the hybridization) and target sequences (around 85%).
# Next "hits" are assigned/combined to "loci" (a stretch of DNA with a particular position in the reference genome) based on the criterion that the length of "introns" separating the hits does not exceed a chosen limit. The limit should be several times the maximum read length (e.g. 3000 bp).
# (A2) Loci with hits ALL shorter than a preset value are excluded (e.g. 70bp), since hits that short are assumed to do not qualify as probes.
# (A3) Next only those query sequences are kept in the analyses which had hits in a preset maximum number of loci. This is to uncover whether loci are duplicated and to estimate the number of copies.
# step B:
# Loci (B1) with lengths of single hits lower than a minimum threshold (e.g. 120), (B2) a total length of hits (e.g. 240) NOT falling below a minimum threshold, and (B3) with intron lengths not exceeding a maximum value (e.g. 500 bp = two times the maximum read length) are marked. In addition, (B4) the number of loci observed for a query sequence, and (B5) the proportion of the query sequences in terms of its length matching the subject is noted.
# For queries with duplicated loci the maximum and minimum sequence divergence among loci is also calculated (actually this is a surrogate calculated from the identities between query and subject for the duplicated loci).
# Finally, the chosen settings are saved in a logfile, the basic results in a summary file (here only queries with with a more stringent maximum number of loci [B6; ideally 2] and loci with a sequence divergence above a preset threshold [B7] are kept), and the full information in another output file. The summary only will include queries and loci fulfilling all criteria for probe selection (A1-3 and B1-7), while the full output will report queries fulfilling criteria A1-3.


##########################################################################################################
############################################## settings ##################################################
##########################################################################################################
# unit is base pairs unless stated otherwise

# criteria for query sequence selection
min_sequence_identity              <- 80   # 80 # minimum identity between query and subject [unit is percent]
min_hit_length                     <- 70   # 70 # hits below this length are ignored
length_loci_significant            <- 10 # 10 exclude loci wihich match only a short proportion of the query sequence AND have low sequence identity
identity_loci_significant          <- 90 # 90
N_loci_first_allowed               <- 6    # 6 number of loci allowed per query sequence in a first selective step [4 to 6 is recommended] [N]
N_loci_finally_allowed             <- 2    # 2 maximum number of loci finally allowed per query sequence [2 is strongly recommended] [N]
intron_length_for_locus_assignment <- 3000 # 3000 # maximum distance between hits to be joined in a locus
required_length_single_hits        <- 120  # 120  # minimum length of hits required for probe design, only loci fulfilling the criterion will be kept
required_total_length_hits         <- 240  # 240 minimum joined length of hits (i.e. the sum of all hits for a locus) required for probe design, only loci fulfilling the criterion will be kept
required_prop_length_hits          <- 0.90 # 0.90 # minimum proportion of the query matching the subject; this setting compets with the previous one (required_total_length_hits)
intron_length_for_locus_selection  <- 500  # 500 only loci with maximum intron lengths lower than this are finally accepted
divergence_among_duplicated_loci   <- 6    # 6 # minimum sequence divergence among duplicated loci required [unit is percent]

# names of input files
tsv                      <- "2019May07_080415.blast.tsv"       # name of tsv-file
filtered_query_sequences <- "genes_mRNA_MalusGDDH13v1_1_for_GDR_Blast_08052019.fasta"  # name of file with filtered sequences for probe design; note that the sequence name and the sequence must be tab delimited and named "seq_no" and "sequence"; name without ">"
#scaffolds                <- "seq_contig_short.md"              # name of file with scaffold info; "seq_contig_short.md" = scaffolds from the MalDomGD1.0 assembly

# name of output files
name_outputfile_data     <- paste("Malus_hits_loci_filtered_",N_loci_first_allowed,"_loci_with_chro00_18072019",sep="")


##########################################################################################################
########################################### read in data #################################################
##########################################################################################################

query_sequences                  <- read.csv(filtered_query_sequences, sep="\t", stringsAsFactors = F)
query_sequences[,"query_length"] <- apply(query_sequences,1,nchar)[2,]
N_query_sequences                <- nrow(query_sequences)

results_blast_search             <- read.csv(tsv, sep="\t", dec=".", stringsAsFactors = F)
results_blast_search             <- subset(results_blast_search,query_acc_ver %in% query_sequences[,1]) # restrict lines to query sequences
results_blast_search             <- subset(results_blast_search, percent_identity > min_sequence_identity) # exclude hits below the preset required idenity of query and subject
results_blast_search             <- results_blast_search[with(results_blast_search, order(query_acc_ver,subject_acc_ver,subject_start)), ]

#scaffolds_on_chromosome          <- read.csv(scaffolds, sep="\t", dec=".", stringsAsFactors = F)


##########################################################################################################
############################################## statistics ################################################
##########################################################################################################

## on hits
min   (results_blast_search[,"percent_identity"])
mean  (results_blast_search[,"percent_identity"])
median(results_blast_search[,"percent_identity"])

## note length of chromosomes
#chromosome_lengths           <- aggregate(scaffolds_on_chromosome,list(scaffolds_on_chromosome$chromosome),max)[,c("Group.1","chr_stop")]
#colnames(chromosome_lengths) <- c("chromosome","length")

print("reading in data completed")
print(Sys.time())
##########################################################################################################
################### group / join hits based on their physical proximity = loci ###########################
##########################################################################################################

query_designation                          <- unique(results_blast_search[,"query_acc_ver"])
results_blast_search[,"dist_next_subject"] <- 99999
for(i in query_designation) {
    chromosome_i <- unique(results_blast_search[which(results_blast_search[,"query_acc_ver"]==i),"subject_acc_ver"])
    for(j in chromosome_i) {
        rows_i   <- which(results_blast_search[,"query_acc_ver"]==i & results_blast_search[,"subject_acc_ver"]==j)
        if(length(rows_i)>1) {
            results_blast_search[rows_i[-length(rows_i)],"dist_next_subject"] <- results_blast_search[rows_i[-1],"subject_start"] - results_blast_search[rows_i[-length(rows_i)],"subject_end"]
            }
        }
    }

print("joining of hits completed")
print(Sys.time())

counter        <- 1
locus          <- NULL
for(i in 1:nrow(results_blast_search)) {
   locus <- c(locus, counter)
   if(results_blast_search[i,"dist_next_subject"] > intron_length_for_locus_assignment) {counter <- counter + 1}
}
results_blast_search[,"locus_ID"] <- locus

print("counting of loci completed")
print(Sys.time())


##########################################################################################################
######################################## statistics loci properties ######################################
##########################################################################################################

## general statistics on hits / loci
alignment_length_total         <- aggregate(results_blast_search[-c(1:2)],list(results_blast_search$query_acc_ver,results_blast_search$locus_ID),sum,na.rm=T)[,c("Group.1","Group.2","alignment_length")] # total sequence length of hits per locus
N_hits                         <- aggregate(results_blast_search         ,list(results_blast_search$query_acc_ver,results_blast_search$locus_ID),length)[,"alignment_length"]
hit_length_min                 <- aggregate(results_blast_search         ,list(results_blast_search$query_acc_ver,results_blast_search$locus_ID),min,na.rm=T)[,"alignment_length"]
hit_length_max                 <- aggregate(results_blast_search         ,list(results_blast_search$query_acc_ver,results_blast_search$locus_ID),max,na.rm=T)[,"alignment_length"]
statistics_loci                <- cbind(alignment_length_total, N_hits, hit_length_min, hit_length_max)
colnames(statistics_loci)[1:3] <- c("query_ID","locus_ID","alignment_length_total")

## filter hits / loci based on the selection criteria
statistics_filtered_loci          <- subset(statistics_loci, hit_length_max > min_hit_length)        # exclude loci with short hits which do not qualify as probes, i.e. all hits shorter than the min_hit_length
ID_loci_filtered                  <- statistics_filtered_loci[,"locus_ID"]                           # note the ID of the kept loci
hits_loci_filtered                <- results_blast_search
hits_loci_filtered                <- merge(hits_loci_filtered,query_sequences[,c("seq_no","query_length")],by.x="query_acc_ver",by.y="seq_no",all.x=T)                     # add length of query sequences
hits_loci_filtered                <- merge(hits_loci_filtered,statistics_filtered_loci[,c("locus_ID","alignment_length_total", "N_hits", "hit_length_min")],by="locus_ID") # ,all.x=T add statistics on loci
hits_loci_filtered[,"prop_match"] <- hits_loci_filtered[,"alignment_length_total"]/hits_loci_filtered[,"query_length"]
m                                 <- aggregate(hits_loci_filtered[,c("locus_ID","percent_identity","prop_match")],list(hits_loci_filtered$locus_ID),mean,na.rm=T)[,-1]
ID_loci_filtered_significant      <- unique(subset(m,percent_identity>=identity_loci_significant | prop_match>=length_loci_significant/100)[,"locus_ID"])                      # filter out insignifcant loci
statistics_filtered_loci          <- statistics_filtered_loci[which(statistics_filtered_loci$locus_ID %in% ID_loci_filtered_significant),]
queries_with_few_loci             <- aggregate(statistics_filtered_loci,list(statistics_filtered_loci$query_ID),length)[,c("Group.1","locus_ID")]                          # estimate the copy number for the loci
queries_with_few_loci             <- subset(queries_with_few_loci, locus_ID <= N_loci_first_allowed) # keep only query sequences with copy numbers below the the threshold "N_loci_first_allowed"
colnames(queries_with_few_loci)   <- c("query_ID","N_loci")
ID_query_filtered                 <- queries_with_few_loci[,"query_ID"]                              # note the ID of the kept query sequences
hits_loci_filtered                <- hits_loci_filtered[which(hits_loci_filtered$query_acc_ver %in% ID_query_filtered & hits_loci_filtered$locus_ID %in% ID_loci_filtered_significant),] # limit hits to those from loci fulfilling the selection criteria min_hit_length and N_loci_first_allowed
hits_loci_filtered                <- merge(hits_loci_filtered,queries_with_few_loci,by.x="query_acc_ver",by.y="query_ID") # add number of loci
hits_loci_filtered[,"dist_next_subject"][hits_loci_filtered$dist_next_subject > intron_length_for_locus_assignment] <- NA
hits_loci_filtered[,"dist_next_subject"][hits_loci_filtered$dist_next_subject == 99999]                             <- NA
intron_length_max                 <- aggregate(hits_loci_filtered[,c("locus_ID","dist_next_subject")],list(hits_loci_filtered$locus_ID),max,na.rm=T)[,c("Group.1","dist_next_subject")]
colnames(intron_length_max)       <- c("locus_ID","length_intron_max")
hits_loci_filtered                <- merge(hits_loci_filtered,intron_length_max,by="locus_ID")  # add max intron length
hits_loci_filtered                <- hits_loci_filtered[,c("query_acc_ver", "locus_ID", "N_loci", "subject_acc_ver", "percent_identity", "mismatches", "gap_opens", "alignment_length", "query_start", "query_end", "subject_start", "subject_end", "dist_next_subject", "length_intron_max", "evalue", "query_length", "N_hits", "alignment_length_total", "prop_match", "hit_length_min")]
colnames(hits_loci_filtered)      <-                     c("query_ID",      "locus_ID", "N_loci", "chromosome",      "identity",         "mismatches", "gap_opens", "alignment_length", "query_start", "query_end", "subject_start", "subject_end", "length_intron"    , "length_intron_max", "evalue", "query_length", "N_hits", "hit_length_total",       "prop_match", "hit_length_min")
# name_Malus_full_statistics_file <- "Malus_hits_loci_filtered_6_loci_with_chro00_18072019.txt" # name full statistics file Malus, this spares running the whole script for some kinds of analyses under changed parameter settings
# hits_loci_filtered              <- read.csv(name_Malus_full_statistics_file, sep="\t", dec=",", stringsAsFactors = F, header = T)[,-which(colnames(hits_loci_filtered) %in% c("N_violations_locus", "N_violations_query"))]
hits_loci_filtered[,c("intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")]      <- 0
hits_loci_filtered[,"hit_lg_tot_critical"][hits_loci_filtered$alignment_length_total < required_total_length_hits]        <- 1
hits_loci_filtered[,"hit_prop_critical"][hits_loci_filtered$"prop_match"             < required_prop_length_hits]         <- 1
hits_loci_filtered[,"min_hit_length_critical"][hits_loci_filtered$hit_length_min     < required_length_single_hits]       <- 1
hits_loci_filtered[,"intron_length_critical"][hits_loci_filtered$length_intron_max   > intron_length_for_locus_selection] <- 1
rownames(hits_loci_filtered) <- seq(1,nrow(hits_loci_filtered),1)


## count violation of criteria for probe selection by loci
evaluation_queries                        <- aggregate(hits_loci_filtered[,c("locus_ID","intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")],list(hits_loci_filtered$locus_ID),max,na.rm=T)[,c("locus_ID","intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")]
evaluation_queries[,"N_violations_locus"] <- rowSums(evaluation_queries[,c("intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")])
hits_loci_filtered                        <- merge(hits_loci_filtered,evaluation_queries[,c( "locus_ID","N_violations_locus")],by="locus_ID",all.x=T)

## count violation of criteria for probe selection by queries
evaluation_queries                        <- aggregate(hits_loci_filtered[,c("query_ID","intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")],list(hits_loci_filtered$query_ID),max,na.rm=T)[,c("query_ID","intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")]
evaluation_queries[,"N_violations_query"] <- rowSums(evaluation_queries[,c("intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")])
hits_loci_filtered                        <- merge(hits_loci_filtered,evaluation_queries[,c( "query_ID","N_violations_query")],by="query_ID",all.x=T)


##########################################################################################################
######################################## statistics loci properties ######################################
##########################################################################################################

N_queries <- length(ID_query_filtered)
# N_queries <- length(unique((hits_loci_filtered[,"query_ID"]))) # use if fullstatistics file has been read in
N_loci    <- length(unique((hits_loci_filtered[,"locus_ID"])))

## extract queries with no violations of criteria for probe design and the the finally allowed number of loci in order to make a short statistics on their proberties
queries_loci_fully_compatible <- unique(subset(hits_loci_filtered, N_violations_query==0 & N_loci <= N_loci_finally_allowed)[,c("query_ID","locus_ID","N_loci","chromosome")])
m                             <-        subset(hits_loci_filtered, N_violations_query==0 & N_loci <= N_loci_finally_allowed)[,c("locus_ID","subject_start","N_hits","hit_length_total")]
n                             <-        subset(hits_loci_filtered, N_violations_query==0 & N_loci <= N_loci_finally_allowed)[,c("query_ID","locus_ID","identity")]
## calculate the minimum und maximum distances among duplicated loci
if(nrow(n)>0) {
     r                          <- aggregate(n,list(n$query_ID,n$locus_ID),mean,na.rm=T)[,c("Group.1","Group.2","identity")]
     colnames(r)                <- c("query_ID","locus_ID","identity")
     r[,c("div_min","div_max")] <- NA
     query_IDs                  <- unique(r[,"query_ID"])
     for(i in query_IDs) {
         s <- subset(r, query_ID == i)
         if(nrow(s) > 1) {
             t <- as.matrix(dist(s$identity))
             diag(t) <- NA
             r[,"div_min"][r$query_ID == i] <- apply(t,2,min,na.rm=T)
             r[,"div_max"][r$query_ID == i] <- apply(t,2,max,na.rm=T)
         }
     }
     r[,"symbol_plot"]                <- sample(c(0,1,3,4,8),nrow(r),replace=T) # randomly assign the duplicated loci pch symbols by which we try to distinguish in the plot below loci placed at the same locus (we expect that the transcript from paralogs map two both loci)
     queries_loci_fully_compatible    <- merge(queries_loci_fully_compatible,r[,c("locus_ID","identity","div_min","div_max","symbol_plot")],by="locus_ID",all.x=T)
}
p                             <- aggregate(m,list(m$locus_ID),mean,na.rm=T)[,c("locus_ID","subject_start","N_hits","hit_length_total")]
queries_loci_fully_compatible <- merge(queries_loci_fully_compatible,p,by="locus_ID",all.x=T)
queries_loci_fully_compatible <- subset(queries_loci_fully_compatible, div_min >= divergence_among_duplicated_loci | is.na(div_min)==T) # exclude insufficiently diverged loci (in case of N_loci_finally_allowed = 2 this means also exclusion of queries!)
queries_loci_fully_compatible[,c("identity","div_min","div_max")] <- round(queries_loci_fully_compatible[,c("identity","div_min","div_max")],2)
colnames(queries_loci_fully_compatible)[which(colnames(queries_loci_fully_compatible)=="identity")] <- "identity_mean"
queries_loci_fully_compatible <- queries_loci_fully_compatible[with(queries_loci_fully_compatible, order(query_ID,chromosome,subject_start)), ] 
queries_loci_fully_compatible[,"subject_start"] <- round(queries_loci_fully_compatible[,"subject_start"])
N_queries_fully_compatible    <- length(unique(queries_loci_fully_compatible[,"query_ID"]))
N_loci_fully_compatible       <- length(unique(queries_loci_fully_compatible[,"locus_ID"]))

query_sequences_for_probes    <- query_sequences[which(query_sequences$seq_no %in% unique(queries_loci_fully_compatible[,"query_ID"])),][,c("seq_no","sequence")]

# length(unique(subset(hits_loci_filtered, hit_prop_critical==0)[,"query_ID"])) # for outgroup evaluation

## make log file with settings and save N queries and N loci kept
settings           <- matrix(c(min_sequence_identity,  min_hit_length,   length_loci_significant,   identity_loci_significant, N_loci_first_allowed,  N_loci_finally_allowed,  intron_length_for_locus_assignment,  required_length_single_hits,  required_total_length_hits,   required_prop_length_hits,   intron_length_for_locus_selection,  divergence_among_duplicated_loci,  tsv,       filtered_query_sequences, N_query_sequences,   N_queries,   N_loci,   N_queries_fully_compatible,   N_loci_fully_compatible), ncol=1)
settings           <- cbind(c("min_sequence_identity","min_hit_length", "length_loci_significant", "identity_loci_significant", "N_loci_first_allowed","N_loci_finally_allowed","intron_length_for_locus_assignment","required_length_single_hits","required_total_length_hits", "required_prop_length_hits", "intron_length_for_locus_selection","divergence_among_duplicated_loci","tsv_file","file_query_sequences",   "N_query_sequences", "N_queries", "N_loci", "N_queries_fully_compatible", "N_loci_fully_compatible"), settings)
colnames(settings) <- c("parameter","value")


## save output
write.table(settings                     ,paste(name_outputfile_data,"_logfile.txt",sep=""),sep="\t",dec=",",row.names=F,quote=F)
write.table(hits_loci_filtered           ,paste(name_outputfile_data,".txt",sep=""),sep="\t",dec=",",row.names=F,quote=F)
write.table(queries_loci_fully_compatible,paste(name_outputfile_data,"_fully_compatible_summary.txt",sep=""),sep="\t",dec=",",row.names=F,quote=F)
write.table(query_sequences_for_probes   ,paste(name_outputfile_data,"_query_sequences_fully_compatible_.fasta",sep=""),sep="\t",row.names=F,quote=F)


##########################################################################################################
########################################### plot position of loci ########################################
##########################################################################################################

queries_loci_fully_compatible[,"symbol_plot"][queries_loci_fully_compatible$N_loci >2] <- 6  # sequences assigned to more than two loci are symbolized by a downright triangle
queries_loci_fully_compatible[,"symbol_plot"][queries_loci_fully_compatible$N_loci==1] <- 19 # sequences assigned to only one loci are symbolized by a black dot
plot(as.numeric(substr(queries_loci_fully_compatible[,"chromosome"],4,5)),queries_loci_fully_compatible[,"subject_start"],xlab=paste("chromosomes 0 to",as.numeric(max(substr(queries_loci_fully_compatible[,"chromosome"],4,5)))),ylab="position [bp]",pch=queries_loci_fully_compatible[,"symbol_plot"])