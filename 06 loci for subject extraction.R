# Loci_for_subject_extraction.R
# first/last update: 8.11.2019/22.11.2019

data_Malus                <- subset(start_end_Pyrus_Malus_loci, genome=="Malus")[,c("query_ID","locus_ID","chromosome","cut_subject_start","cut_subject_end","length_subject_approx","length_subject_real")]
data_Malus[,"chromosome"] <- as.numeric(substring(data_Malus[,"chromosome"],4,5))
data_Malus[,"chromosome"][data_Malus[,"chromosome"]==0] <- 18 # chromosome 0 (the non anchored contigs) are last in the Malus chromosome file and therefore is designated the 18th

data_Pyrus                          <- subset(start_end_Pyrus_Malus_loci, genome=="Pyrus")[,c("query_ID","locus_ID","chromosome","cut_subject_start","cut_subject_end","length_subject_approx","length_subject_real")]
chro_order                          <- read.csv("order of chromosomes in Pyrus communis genome file.txt", sep="\t", stringsAsFactors = F)
chros                               <- sort(unique(data_Pyrus[,"chromosome"]))
data_Pyrus[,"order_in_genome_file"] <- 0
for(i in chro_order[,"chromosome"]) {
    data_Pyrus[,"order_in_genome_file"][data_Pyrus[,"chromosome"]==i] <- chro_order[chro_order$chromosome==i,"order_genome_file"]
}
data_Pyrus                          <- data_Pyrus[,c("query_ID","locus_ID","chromosome","order_in_genome_file","cut_subject_start","cut_subject_end","length_subject_approx","length_subject_real")]

rownames(data_Malus) <- seq(1,nrow(data_Malus),1)
rownames(data_Pyrus) <- seq(1,nrow(data_Pyrus),1)
seqs_reversed_Malus  <- rownames(subset(data_Malus, length_subject_approx <0))
seqs_reversed_Pyrus  <- rownames(subset(data_Pyrus, length_subject_approx <0))

write.table(data_Malus, "Malus_loci_for_subject_extraction.txt", sep="\t",row.names=F,quote=F)
write.table(data_Pyrus, "Pyrus_loci_for_subject_extraction.txt", sep="\t",row.names=F,quote=F)