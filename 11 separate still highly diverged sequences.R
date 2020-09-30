# 11 separate still highly diverged seuqnces.R
# first/last update: 28.11.2019/05.12.2019
# stores still highly diverged sequences in separate files
# and extracts the longest sequence of each alignment. The sequence is attached to a single file holding all sequences, which are named by the name of the alignment.
# alignments are read in from "finally selected sequences" and saved in "sequences to be send"

# set wd to "...\files for bait design\"
# make subdirectories sequences to be send, sequences to be send/exons, and sequences to be send/introns

# library(stringdist); library(plyr); library(DECIPHER)


##############################################################################################
######################################### settings ###########################################
##############################################################################################

max_div_exons   <- 0.08       # value above which sequences are separated
max_div_introns <- 0.08
M_cluster       <- "ward.D"   # clustering method of hclust (average, ward.D,...)

file_statistics_sequence_comparison <- "statistics Malus and Pyrus run 2019-12-05 final selection.txt" # name of the statistics file output by script sequence comparison.R


##########################################################################################################
######################################### execute program ################################################
##########################################################################################################

FileNamesExons   <- choose.files(".fa") # select finally selected exons
FileNamesIntrons <- choose.files(".fa") # select finally selected introns

wd       <- getwd()

data_out      <- NULL
sequence_file <- NULL
for(g in 1:2) {
    if(g==1) {FileNames <- FileNamesExons;   feature <- "exon"}
    if(g==2) {FileNames <- FileNamesIntrons; feature <- "intron"}
    for(h in 1:length(FileNames)) {
        save_as_is <- F
        variant    <- 1
        file       <- read.csv(FileNames[h], stringsAsFactors = F, header=F)
        if(nrow(file) > 2) {


##############################################################################################
################## calculate sequence divergence for all exons and introns ###################
##############################################################################################

            div     <- stringdistmatrix(toupper(file[seq(2,nrow(file),2),]))/nchar(file[2,])
            max_div <- max(div)
            if((max_div > max_div_exons & g==1) | (max_div > max_div_introns & g==2)) {
                fit        <- hclust(div, method = M_cluster)
                fit$labels <- file[seq(1,nrow(file),2),]
                branches   <- cutree(fit, 2)
                file_partition           <- as.data.frame(cbind(names(branches),as.numeric(branches)),stringsAsFactors=F)
                colnames(file_partition) <- c("seq_name","group")
                sequ_order               <- as.data.frame(cbind(file[seq(1,nrow(file),2),],seq(1,nrow(file),2)),stringsAsFactors=F)
                colnames(sequ_order)     <- c("seq_name","order")
                m                        <- merge(file_partition,sequ_order)
                rows                             <- as.numeric(m[,"order"][m$group==1])
                rows                             <- sort(c(rows,(rows+1)))
                seq_file1                        <- as.matrix(file[rows,],ncol=1)
                seq_file1[seq(1,length(rows),2)] <- file[rows[seq(1,length(rows),2)],]
                rows                             <- as.numeric(m[,"order"][m$group==2])
                rows                             <- sort(c(rows,(rows+1)))
                seq_file2                        <- as.matrix(file[rows,],ncol=1)
                seq_file2[seq(1,length(rows),2)] <- file[rows[seq(1,length(rows),2)],]
            } else { # sequences are not diverged beyond the set threshold
                save_as_is <- T
                seq_file1  <- file
                seq_file2  <- NULL
            }
        } else { # alignments holds one sequence
            save_as_is <- T
            seq_file1  <- file[1:2,]
            seq_file2  <- NULL
        }
        len   <- nchar(FileNames[h])
        if(g==1) query <- substring(FileNames[h],(nchar(getwd())+35),len-3)
        if(g==2) query <- substring(FileNames[h],(nchar(getwd())+37),len-3)
        if(save_as_is==F) {
            if(g==1) {file_name_1 <- paste(wd,"/sequences to be send/exons/",query,"_subset_1.fa",sep="");   file_name_2 <- paste(wd,"/sequences to be send/exons/",query,"_subset_2.fa",sep="");   seq_name_1 <- substring(file_name_1,(nchar(getwd())+29),len); seq_name_2 <- substring(file_name_2,(nchar(getwd())+29),len)}
            if(g==2) {file_name_1 <- paste(wd,"/sequences to be send/introns/",query,"_subset_1.fa",sep=""); file_name_2 <- paste(wd,"/sequences to be send/introns/",query,"_subset_2.fa",sep=""); seq_name_1 <- substring(file_name_1,(nchar(getwd())+31),len); seq_name_2 <- substring(file_name_2,(nchar(getwd())+31),len)}
            for(i in 1:2) # loop for file 1 and 2
            {
                if(i==1) {file <- seq_file1; file_name <- file_name_1;               seq_name <- paste(seq_name_1,"variant",variant,sep="_")}
                if(i==2) {file <- seq_file2; file_name <- file_name_2; variant <- 2; seq_name <- paste(seq_name_1,"variant",variant,sep="_")}
                for(j in seq(2,nrow(file),2)) {
                    m                <- strsplit(file[j,],NULL)[[1]]
                    coordinates_gaps <- which(m=="-")
                    if(length(coordinates_gaps)>0) {file[j,] <- paste(m[-coordinates_gaps],collapse="")} else {file[j,] <- paste(m,collapse="")}
                }
                n        <- file[seq(2,nrow(file),2),]
                o        <- nchar(n)
                sequence <- n[which(o == max(o))][1]
                if(nrow(file)>2) {
                    write.table(file,file_name,row.names=FALSE,col.names=FALSE,quote=F)
                    dna                  <- readDNAStringSet(file_name)
                    aligned_dna          <- AlignSeqs(dna)
                    aligned_adjusted_dna <- AdjustAlignment(aligned_dna)
                    write.table(aligned_adjusted_dna,file_name,col.names=FALSE,quote=F)
                    file                 <- read.csv(file_name,sep=" ",header=F,stringsAsFactors = F)
                    file[,1]             <- paste(">",file[,1],sep="")
                    file                 <- c(as.character(file[,1]),as.character(file[,2]))
                    m                    <- (length(file)/2)
                    sort_key             <- NULL
                    for(j in 1:m) {
                        sort_key_j <- c(j,j+m)
                        sort_key   <- c(sort_key,sort_key_j)
                    }
                    file            <- matrix(file[sort_key],ncol=1)
                    div_realigned   <- stringdistmatrix(toupper(file[seq(2,nrow(file),2),]))/nchar(file[2,1])
                    max_div_realign <- max(div_realigned)
                }
                write.table(file,file_name,row.names=FALSE,col.names=FALSE,quote=F)
                if(nrow(file)>2) {div <- round(as.numeric(max(stringdistmatrix(toupper(file[seq(2,nrow(file),2),]))/nchar(file[2,]))),2)} else {div <- NA}
                data_out      <- rbind(data_out,c(query,variant,feature,max(o),nrow(file)/2,div))
                sequence_file <- rbind(sequence_file,c(seq_name,sequence))
            }
        } else {
            if(g==1) {file_name <- paste(wd,"/sequences to be send/exons/",query,".fa",sep="");   seq_name <- substring(file_name,(nchar(getwd())+29),len-9)}
            if(g==2) {file_name <- paste(wd,"/sequences to be send/introns/",query,".fa",sep=""); seq_name <- substring(file_name,(nchar(getwd())+31),len-9)}
            write.table(seq_file1,file_name,row.names=FALSE,col.names=FALSE,quote=F)
            if(nrow(file)>2) {div <- round(as.numeric(max(stringdistmatrix(toupper(file[seq(2,nrow(file),2),]))/nchar(file[2,]))),2)} else {div <- NA}
            for(j in seq(2,nrow(file),2)) {
                m                <- strsplit(file[j,],NULL)[[1]]
                coordinates_gaps <- which(m=="-")
                if(length(coordinates_gaps)>0) {file[j,] <- paste(m[-coordinates_gaps],collapse="")} else {file[j,] <- paste(m,collapse="")}
            }
            n        <- file[seq(2,nrow(file),2),]
            o        <- nchar(n)
            sequence <- n[which(o == max(o))][1]
            data_out      <- rbind(data_out,c(query,variant,feature,max(o),nrow(file)/2,div))
            sequence_file <- rbind(sequence_file,c(seq_name,sequence))
        }
    } # end h-loop
} # end g-loop

data_out           <- as.data.frame(data_out)
colnames(data_out) <- c("alignment","variant","feature","length bp","N sequences","maximum_div")
write.table(data_out,     paste(substring(file_statistics_sequence_comparison,1,(nchar(file_statistics_sequence_comparison)-4)),"length alignments to be send.txt",sep=" "),sep="\t",dec=",",row.names=FALSE,quote=F)
write.table(sequence_file,paste(substring(file_statistics_sequence_comparison,1,(nchar(file_statistics_sequence_comparison)-4)),"sequences to be send.txt",sep=" ")        ,sep="\t",dec=",",row.names=FALSE,col.names=FALSE,quote=F)


