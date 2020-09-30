# join_full_statistics_fully_compatible_genomes.R
# first/last update: 24.7.2019/24.7.2019
# first run the probe_mapping Pyrus.R script (the ID_m is required!)

malus <- read.csv("Malus_hits_loci_filtered_6_loci_with_chro00_18072019.txt", sep="\t", stringsAsFactors = F, header=T)[,-25]
pyrus <- read.csv("Pyrus_hits_loci_filtered_6_loci_with_chro00_18072019.txt", sep="\t", stringsAsFactors = F, header=T)
colnames(malus)[25] <- "N_violations"
column.names <- colnames(malus)
malus[,"genome"] <- "Malus"
pyrus[,"genome"] <- "Pyrus"
both  <- rbind(malus,pyrus)
both_fully_compatible <- subset(both, query_ID %in% ID_m & N_violations == 0)[,c("genome",column.names)] # shorter loci (applying the "required_prop_length_hits"-limit) were not excluded yet from the full tables
both_fully_compatible <- both_fully_compatible[with(both_fully_compatible, order(query_ID,genome)),]

#ID_m <- unique(m$query_ID)
#dim(both_fully_compatible)
#length(ID_m)
#length(unique(both_fully_compatible$query_ID))

#colnames(malus)
#colnames(pyrus)

write.table(both_fully_compatible,"Pyrus_Malus_detailed_statistics_fully_compatible_queries_18072019.txt",sep="\t",dec=",",row.names=F,quote=F)
