# mRNA selection for Blast search.R
# first/last update: 7.5.2019/15.7.2019

###############################################
############### MalusGDDH1 ####################
###############################################

statistics                    <- read.csv("genes_mRNA_MalusGDDH13v1_1_from_GDR_table.txt", sep=",", stringsAsFactors = F)
statistics[,"length_on_chro"] <- statistics[,"Stop"] - statistics[,"Start"]
statistics                    <- subset(statistics, Type == "mRNA" & length_on_chro >= 900 & length_on_chro <= 5000)
seqIDs_anchored <- statistics$Name

sequences                <- read.csv("genes_mRNA_MalusGDDH13v1_1_from_GDR_fasta.txt", sep=";", stringsAsFactors = F)
sequences                <- sequences[,-2]
sequences[,"length_seq"] <- apply(sequences,1,nchar)[2,]
sequences[,"seqID"]      <- substr(sequences[,"seqID"],2,13)
sequences                <- subset(sequences,length_seq >= 800 & length_seq <= 3000 & seqID %in% seqIDs_anchored)

statistics <- merge(sequences[,c("seqID","length_seq")],statistics[,c("Name","Source","Start","Stop","length_on_chro")],by.x="seqID",by.y="Name",all.x=T)
statistics <- statistics[,c("seqID","Source","length_seq","Start","Stop","length_on_chro")]

write.table(sequences[,-3],"genes_mRNA_MalusGDDH13v1_1_for_GDR_Blast_09052019.fasta",sep=";",row.names=F)
write.table(statistics    ,"genes_mRNA_MalusGDDH13v1_1_for_GDR_Blast_statistics_09052019.txt",sep=";",row.names=F)


###############################################
################## HFTH1 ######################
###############################################

sequences                <- read.csv("genes_mRNA_MalusHFTH1_1_for_GDR_Blast_16072019_all_mRNAs.original.fasta", sep=";", header=T, stringsAsFactors=F)
sequences[,"length_seq"] <- apply(sequences,1,nchar)[2,]
sequences                <- subset(sequences,length_seq >= 800 & length_seq <= 3000)

write.table(sequences[,-3],"genes_mRNA_MalusHFTH1_1_for_GDR_Blast_16072019_length_filtered_mRNAs.original.fasta",sep=";",row.names=F,quote=F)
