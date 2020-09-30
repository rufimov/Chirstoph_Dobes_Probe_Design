# probe_mapping Pyrus.R
# first/last update: 9.4.2019/23.8.2019

##########################################################################################################
########################################### script desprition ############################################
##########################################################################################################

# The script is an extention of the probe mapping.R script. It first reads in a tsv-files which holds the output of a blast search (i.e. "hits" of "query sequences" against a "subject") created using the "Genome Database for Rosaceae": for instance a search of transiptome sequences filtered from Malus (= "query sequences") against the assembled Pyrus bretschneideri genome (= "subject").
# The query sequences are also read in in fasta format.
# Queries suitable for probe design are selected in two steps (A and B):
# step A:
# (A1) Hits with a sequence identity relative to the subject lower than a preset value are excluded. 
# Next "hits" are assigned/combined to "loci" (a stretch of DNA with a particular position in the reference genome) based on the criterion that the length of "introns" separating the hits does not exceed a chosen limit. The limit should be considerably higher than that chosen for Malus (e.g. 10000 bp).
# (A2) Loci with hits ALL shorter than a preset value are excluded (e.g. 70bp), since hits that short are assumed to do not qualify as probes.
# (A3) Next only those query sequences are kept in the analyses which had hits in a preset maximum number of loci. This is to uncover whether loci are duplicated and to estimate the number of copies.
# step B:
# Loci (B1) with lengths of single hits lower than a minimum threshold (e.g. 120), (B2) a total length of hits (e.g. 240) NOT falling below a minimum threshold, and (B3) with intron lengths not exceeding a maximum value (e.g. 500 bp = two times the maximum read length) are marked (penalized). In addition, (B4) the number of loci observed for a query sequence, and (B5) the proportion of the query sequences in terms of its length matching the subject is noted. (B6) For queries with duplicated loci the sequence maximum and minimum divergence among loci is also calculated (actually this is a surrogate calculated from the identities between query and subject for the duplicated loci).
# Finally, the chosen settings are saved in a logfile, the basic results in summary files and the full information in another output file. Both only will include queries and loci fulfilling all criteria for probe selection (A1-3).
# Note that queries/loci are only penalized and saved in the two summary files: summary (Pyrus only) and summaries_both_genomes (this is a difference to the summary file created for the Blast against the Malus genome to which only fully compatible loci are saved). This allows to decide individually which query sequences are finally kept.
# However, a fasta file of those query sequences which are fully compatible for both genomes is saved.


##########################################################################################################
############################################## settings ##################################################
##########################################################################################################
# unit is base pairs unless stated otherwise

# criteria for query sequence selection
min_sequence_identity              <- 80   # minimum identity between query and subject [unit is percent]
min_hit_length                     <- 70   # hits below this length are ignored
length_loci_significant            <- 10 # exclude loci wihich match only a short proportion of the query sequence AND have low sequence identity [unit is percent]
identity_loci_significant          <- 90 # 90 [unit is percent]
N_loci_first_allowed               <- 6    # number of loci allowed per query sequence in a first selective step [4 to 6 is recommended] [N]
N_loci_finally_allowed             <- 2    # maximum number of loci finally allowed per query sequence [2 is strongly recommended] [N]
intron_length_for_locus_assignment <- 10000 # maximum distance between hits to be joined in a locus
required_length_single_hits        <- 120  # minimum length of hits required for probe design, only loci fulfilling the criterion will be kept
required_total_length_hits         <- 240  # minimum joined length of hits (i.e. the sum of all hits for a locus) required for probe desqign, only loci fulfilling the criterion will be kept
preferred_prop_length_hits         <- 0.90 # minimum proportion of the query matching the subject; this setting compets with the previous one (required_total_length_hits) [proportion]
required_prop_length_hits          <- 0.20 # minimum proportion of the query matching the subject in order be kept [proportion]
intron_length_for_locus_selection  <- 500  # only loci with maximum intron lengths lower than this are finally accepted
locus_sequence_divergence          <- T    # whether sequence divergence among loci should be considered as selection criterion for fully compatible loci: Yes (T), No (F)
    divergence_among_duplicated_loci   <- 6 # 6 minimum sequence divergence among duplicated loci required [unit is percent]
hits_in_both_genomes               <- F    # whether full compatibilty of queries should include that queries have hits in both genomes: Yes (T), No (F)

exclude_scaffolds                  <- F    # whether scaffolds, i.e. DNA blocks not assigned to chromosomes, should be excluded from the tsv file (T) or kept (F).

# names of input files
tsv                       <- "2019Aug21_034923.blast.Pyrus.tsv"                                            # name of tsv-file
filtered_query_sequences  <- "Malus_hits_loci_filtered_6_loci_with_chro00_18072019_query_sequences_for_Pyrus_RSearch.fasta"  # name of file with filtered sequences potentially suitable for probe design
summary_file_first_genome <- "Malus_hits_loci_filtered_6_loci_with_chro00_18072019_summary.txt"            # name of the file summarizing the results of the first genome
                              
# name of output files
name_outputfile_data     <- paste("Pyrus_hits_loci_filtered_",N_loci_first_allowed,"_loci_with_chro00_18072019",sep="")


##########################################################################################################
########################################### read in data #################################################
##########################################################################################################

query_sequences                  <- read.csv(filtered_query_sequences, sep="\t", stringsAsFactors = F)
query_sequences[,"query_length"] <- apply(query_sequences,1,nchar)[2,]
N_queries_blasted                <- nrow(query_sequences)

results_blast_search          <- read.csv(tsv, sep="\t", dec=".", stringsAsFactors = F)
results_blast_search          <- subset(results_blast_search,query_acc_ver %in% query_sequences[,1]) # restrict lines to query sequences
results_blast_search          <- subset(results_blast_search, percent_identity > min_sequence_identity) # exclude hits below the preset required idenity of query and subject
results_blast_search          <- results_blast_search[with(results_blast_search, order(query_acc_ver,subject_acc_ver,subject_start)), ]
subject_acc_ver_Pyrus         <- sort(unique(results_blast_search$subject_acc_ver))
subject_acc_ver_chromos_Pyrus <- subject_acc_ver_Pyrus[grep("chr", subject_acc_ver_Pyrus)]
if(exclude_scaffolds==T) {results_blast_search <- subset(results_blast_search, subject_acc_ver %in% subject_acc_ver_chromos_Pyrus)}
    # rename chromosomes
for(i in subject_acc_ver_chromos_Pyrus) {
   if(substr(i,5,5)=="-") {m <- substr(i,1,4)} else {m <- substr(i,1,5)}
   results_blast_search[,"subject_acc_ver"][results_blast_search$subject_acc_ver==i] <- m
   }


##########################################################################################################
############################################## statistics ################################################
##########################################################################################################

## on hits
min   (results_blast_search[,"percent_identity"])
mean  (results_blast_search[,"percent_identity"])
median(results_blast_search[,"percent_identity"])


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

counter        <- 1
locus          <- NULL
for(i in 1:nrow(results_blast_search)) {
   locus <- c(locus, counter)
   if(results_blast_search[i,"dist_next_subject"] > intron_length_for_locus_assignment) {counter <- counter + 1}
}
results_blast_search[,"locus_ID"] <- locus


##########################################################################################################
######################################## statistics loci properties ######################################
##########################################################################################################

## general statistics on hits / loci
alignment_length_total         <- aggregate(results_blast_search[-c(1:2)],list(results_blast_search$query_acc_ver,results_blast_search$locus_ID),sum,na.rm=T)[,c("Group.1","Group.2","alignment_length")] # total sequence length of hits per locus
N_hits                         <- aggregate(results_blast_search,         list(results_blast_search$query_acc_ver,results_blast_search$locus_ID),length)[,"alignment_length"]
hit_length_min                 <- aggregate(results_blast_search,         list(results_blast_search$query_acc_ver,results_blast_search$locus_ID),min,na.rm=T)[,"alignment_length"]
hit_length_max                 <- aggregate(results_blast_search,         list(results_blast_search$query_acc_ver,results_blast_search$locus_ID),max,na.rm=T)[,"alignment_length"]
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
hits_loci_filtered                <- merge(hits_loci_filtered,queries_with_few_loci,by.x="query_acc_ver",by.y="query_ID")
hits_loci_filtered[,"dist_next_subject"][hits_loci_filtered$dist_next_subject > intron_length_for_locus_assignment] <- NA
hits_loci_filtered[,"dist_next_subject"][hits_loci_filtered$dist_next_subject == 99999]                             <- NA
intron_length_max                 <- aggregate(hits_loci_filtered[,c("locus_ID","dist_next_subject")],list(hits_loci_filtered$locus_ID),max,na.rm=T)[,c("Group.1","dist_next_subject")]
colnames(intron_length_max)       <- c("locus_ID","length_intron_max")
hits_loci_filtered                <- merge(hits_loci_filtered,intron_length_max,by="locus_ID")
hits_loci_filtered                <- hits_loci_filtered[,c("query_acc_ver", "locus_ID", "N_loci", "subject_acc_ver", "percent_identity", "mismatches", "gap_opens", "alignment_length", "query_start", "query_end", "subject_start", "subject_end", "dist_next_subject", "length_intron_max", "evalue", "query_length", "N_hits", "alignment_length_total", "prop_match", "hit_length_min")]
colnames(hits_loci_filtered)      <-                     c("query_ID",      "locus_ID", "N_loci", "chromosome",      "identity",         "mismatches", "gap_opens", "alignment_length", "query_start", "query_end", "subject_start", "subject_end", "length_intron"    , "length_intron_max", "evalue", "query_length", "N_hits", "hit_length_total",       "prop_match", "hit_length_min")
hits_loci_filtered[,c("intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")]      <- 0
hits_loci_filtered[,"hit_lg_tot_critical"][hits_loci_filtered$alignment_length_total < required_total_length_hits]        <- 1
hits_loci_filtered[,"hit_prop_critical"][hits_loci_filtered$prop_match               < preferred_prop_length_hits]        <- 1
hits_loci_filtered[,"min_hit_length_critical"][hits_loci_filtered$hit_length_min     < required_length_single_hits]       <- 1
hits_loci_filtered[,"intron_length_critical"][hits_loci_filtered$length_intron_max   > intron_length_for_locus_selection] <- 1
rownames(hits_loci_filtered) <- seq(1,nrow(hits_loci_filtered),1)

## count violation of criteria for probe selection by queries
evaluation_queries                  <- aggregate(hits_loci_filtered[,c("query_ID","intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")],list(hits_loci_filtered$query_ID),max,na.rm=T)[,c("query_ID","intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")]
evaluation_queries[,"N_violations"] <- rowSums(evaluation_queries[,c("intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical")])
hits_loci_filtered                  <- merge(hits_loci_filtered,evaluation_queries[,c( "query_ID","N_violations")],by="query_ID",all.x=T)


##########################################################################################################
######################################## statistics loci properties ######################################
##########################################################################################################

## here shorter loci are excluded (applying the "required_prop_length_hits"-limit)
queries_summary_statistics <- unique(subset(hits_loci_filtered, prop_match >= required_prop_length_hits)[,c("query_ID","locus_ID","N_loci","chromosome")])
m                          <-        subset(hits_loci_filtered, prop_match >= required_prop_length_hits)[,c("locus_ID","subject_start","N_hits","hit_length_total")]
n                          <-        subset(hits_loci_filtered, prop_match >= required_prop_length_hits)[,c("query_ID","locus_ID","identity")]
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
     queries_summary_statistics       <- merge(queries_summary_statistics,r[,c("locus_ID","identity","div_min","div_max","symbol_plot")],by="locus_ID",all.x=T)
}
p                          <- aggregate(m,list(m$locus_ID),mean,na.rm=T)[,c("locus_ID","subject_start","N_hits","hit_length_total")]
queries_summary_statistics <- merge(queries_summary_statistics,p,by="locus_ID",all.x=T)
queries_summary_statistics <- queries_summary_statistics[with(queries_summary_statistics, order(query_ID,chromosome,subject_start)), ] 
queries_summary_statistics[,"subject_start"] <- round(queries_summary_statistics[,"subject_start"])
N_queries_for_penalization             <- length(unique(queries_summary_statistics[,"query_ID"]))
N_loci_for_penalization                <- length(unique(queries_summary_statistics[,"locus_ID"]))

## merge summaries of first and second genome
summary_first_genome                   <- read.csv(summary_file_first_genome, sep="\t", dec=",", stringsAsFactors = F)
summary_first_genome[,"genome"]        <- "Malus"
summary_first_genome[,c("intron_length_critical","hit_lg_tot_critical","hit_prop_critical","min_hit_length_critical","N_violations")]  <- 0
summary_second_genome                  <- queries_summary_statistics[,c("locus_ID","query_ID","N_loci","chromosome","identity","div_min","div_max","symbol_plot","subject_start","N_hits","hit_length_total")]
summary_second_genome[,"genome"]       <- "Pyrus"
summary_second_genome                  <- merge(summary_second_genome,evaluation_queries,by="query_ID",all.x=T)
colnames(summary_second_genome)[which(colnames(summary_second_genome)=="identity")] <- "identity_mean"
summary_both_genomes                   <- rbind(summary_first_genome, summary_second_genome)
summary_both_genomes[,c("identity_mean","div_min","div_max")] <- round(summary_both_genomes[,c("identity_mean","div_min","div_max")],2)
summary_both_genomes[,"chromosome_no"] <- as.numeric(substr(summary_both_genomes[,"chromosome"],4,5))
summary_both_genomes                   <- summary_both_genomes[with(summary_both_genomes, order(query_ID,genome,chromosome_no,subject_start)), ]
summary_both_genomes                   <- summary_both_genomes[,c("query_ID","genome","locus_ID","N_loci","chromosome","chromosome_no","subject_start","hit_length_total","identity_mean","div_min","div_max","N_hits","N_violations")]
m                                      <- aggregate(summary_both_genomes[,c("query_ID","N_loci")],list(summary_both_genomes$query_ID),max,na.rm=T)[,c("query_ID","N_loci")]
summary_both_genomes                   <- merge(summary_both_genomes,m,by="query_ID")
colnames(summary_both_genomes)[which(colnames(summary_both_genomes) %in% c("N_loci.x","N_loci.y"))] <- c("N_loci","N_loci_max_ges")
m                                      <- aggregate(summary_both_genomes[,c("query_ID","N_violations")],list(summary_both_genomes$query_ID),max,na.rm=T)[,c("query_ID","N_violations")]
summary_both_genomes                   <- merge(summary_both_genomes,m,by="query_ID")
colnames(summary_both_genomes)[which(colnames(summary_both_genomes) %in% c("N_violations.x","N_violations.y"))] <- c("N_violations","N_violations_both")

## make log file with settings and save N queries and N loci fully compatible
m                                        <- subset(summary_both_genomes, N_loci_max_ges <= N_loci_finally_allowed & N_violations_both == 0)
if(locus_sequence_divergence == T) {ID_m <- unique(subset(m, div_min <= divergence_among_duplicated_loci)[,"query_ID"])
m                                        <- subset(m, !query_ID %in% ID_m)}
if(hits_in_both_genomes == T) {ID_m <- unique(subset(m, genome == "Pyrus")[,"query_ID"])
m                                        <- subset(m, query_ID %in% ID_m)}
ID_m                                    <- unique(m[,"query_ID"])
N_queries_fully_compatible_both_genomes <- length(ID_m)
N_loci_fully_compatible_both_genomes    <- length(unique(m[,"locus_ID"]))
prop_queries_compatible_both_genomes <- round(N_queries_fully_compatible_both_genomes/length(unique(summary_first_genome$query_ID)),2)
settings                             <- matrix(c(min_sequence_identity,  min_hit_length,  length_loci_significant,  identity_loci_significant,  N_loci_first_allowed,  N_loci_finally_allowed,  intron_length_for_locus_assignment,  required_length_single_hits,  required_total_length_hits,   preferred_prop_length_hits,   intron_length_for_locus_selection,  locus_sequence_divergence,  divergence_among_duplicated_loci,  tsv,       filtered_query_sequences, N_queries_blasted,   N_queries_for_penalization,   N_loci_for_penalization,   hits_in_both_genomes,   N_queries_fully_compatible_both_genomes,   prop_queries_compatible_both_genomes,   N_loci_fully_compatible_both_genomes ), ncol=1)
settings                             <- cbind(c("min_sequence_identity","min_hit_length","length_loci_significant","identity_loci_significant","N_loci_first_allowed","N_loci_finally_allowed","intron_length_for_locus_assignment","required_length_single_hits","required_total_length_hits", "preferred_prop_length_hits", "intron_length_for_locus_selection","locus_sequence_divergence","divergence_among_duplicated_loci","tsv_file","file_query_sequences",   "N_queries_blasted", "N_queries_for_penalization", "N_loci_for_penalization", "hits_in_both_genomes", "N_queries_fully_compatible_both_genomes", "prop_queries_compatible_both_genomes", "N_loci_fully_compatible_both_genomes"), settings)
colnames(settings)                   <- c("parameter","value")

# extract query sequences fully compatible for both genomes
query_sequences_fully_compatible_both_genomes <- query_sequences[which(query_sequences$seq_no %in% ID_m),1:2]

## save output
write.table(settings                                     ,paste(name_outputfile_data,"_logfile.txt",sep=""),sep="\t",dec=",",row.names=F,quote=F)
write.table(hits_loci_filtered                           ,paste(name_outputfile_data,"_Pyrus.txt",sep=""),sep="\t",dec=",",row.names=F,quote=F)
write.table(queries_summary_statistics                   ,paste(name_outputfile_data,"Pyrus_all_queries_summary.txt",sep=""),sep="\t",dec=",",row.names=F,quote=F)
write.table(summary_both_genomes                         ,paste(name_outputfile_data,"_summaries_both_genomes_all_blasted_queries.txt",sep=""),sep="\t",dec=",",row.names=F,quote=F)
write.table(query_sequences_fully_compatible_both_genomes,paste(name_outputfile_data,"_queries_fully_compatible_both_genomes.fasta",sep=""),sep="\t",dec=",",row.names=F,quote=F)

##########################################################################################################
########################################### plot position of loci ########################################
##########################################################################################################

# visualize position of loci in the Pyrus genome (penalized sequences)
queries_summary_statistics[,"symbol_plot"][queries_summary_statistics$N_loci >2] <- 6  # sequences assigned to more than two loci are symbolized by a downright triangle
queries_summary_statistics[,"symbol_plot"][queries_summary_statistics$N_loci==1] <- 19 # sequences assigned to only one loci are symbolized by a black dot
plot(as.numeric(substr(queries_summary_statistics[,"chromosome"],4,5)),queries_summary_statistics[,"subject_start"],xlab=paste("chromosomes 0 to",max(as.numeric(substr(queries_summary_statistics[,"chromosome"],4,5)))),ylab="position [bp]",pch=queries_summary_statistics[,"symbol_plot"])
