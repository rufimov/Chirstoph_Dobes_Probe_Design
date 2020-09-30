# 10 final select.R
# first/last update: 27.11.2019/05.12.2019
# final selection of loci

make directories "finally selected sequences", "finally selected sequences/exons", "finally selected sequences/introns"

name_statistics_file <- "statistics Malus and Pyrus run 2019-12-05.txt"

min_exon_length      <- 80 # minimum length of exons required

max_div_exons        <- 0.15 # maxial allowed divergence of exons: applies to single-copy loci only (in this case divergence refers to that among orthologs)

########################################################################################

statistics <- read.csv(name_statistics_file, sep="\t", dec=".", stringsAsFactors = F)

rows1 <- as.numeric(rownames(subset(statistics,  statistics[,"length.exon.min"] < min_exon_length)))
rows2 <- as.numeric(rownames(subset(statistics,  statistics[,"N.loci.Malus"]==1 & statistics[,"max.div.exons"] > max_div_exons)))
rows3 <- as.numeric(rownames(subset(statistics, (statistics[,"N.loci.Malus"]==1 & statistics[,"N.loci.Pyrus"]!=1) | (statistics[,"N.loci.Malus"]==2 & statistics[,"N.loci.Pyrus"]!=2))))
rows  <- sort(unique(c(rows1,rows2,rows3)))
cases_to_be_included <- substring(statistics[-rows,"queries"],1,12)
cases_excluded       <- statistics[rows,"queries"]

write.table(statistics[-rows,],paste(substring(name_statistics_file,1,(nchar(name_statistics_file)-4)),"final selection.txt",sep=" "),sep="\t",row.names=F,quote=F)

wd <- getwd()

########################################################################################
####################################### Exons ##########################################
########################################################################################

setwd(paste(wd,"/selected sequences/exons",sep=""))

current_folder_exons <- getwd()
new_exons_folder     <- paste(substring(getwd(),1,106),"/finally selected sequences/exons",sep="")

for(i in cases_to_be_included) {
    list.of.files <- list.files(current_folder_exons, i)
    file.copy(list.of.files, new_exons_folder)
}


########################################################################################
###################################### Introns #########################################
########################################################################################

setwd(paste(wd,"/selected sequences/introns",sep=""))

current_folder_introns <- getwd()
new_intron_folder      <- paste(substring(getwd(),1,106),"/finally selected sequences/introns",sep="")

for(i in cases_to_be_included) {
    list.of.files <- list.files(current_folder_introns, i)
    file.copy(list.of.files, new_intron_folder)
}

setwd(wd)
