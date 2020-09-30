# 08 sequence comparison.R
# first/last update: 11.11.2019/22.11.2019
# identifies exon/inron limits, calculates sequence divergence among all exonic and intronic regions, partitions sequences according to a set divergence threshold, and saves exons and introns in separate files
# in addition, the start and end of alignments is trimmed
# introns below preset length are excluded as are introns which exceed as preset level of sequence divergence

# library(stringdist); library(plyr); library(DECIPHER)


##############################################################################################
######################################### settings ###########################################
##############################################################################################

# generell
min_length <- 80 # minimum length of introns for selection for bait design

# for alignment limits
size_window         <- 20   # size of the sliding window
threshold_div_exons <- 0.20 # the alignment starts and ends when this level of sequence divergence is untercut

# for sequence divergence
max_div_exons   <- 0.05 # value above which sequences are separated
max_div_introns <- 0.05

max_allowed_div_introns <- 0.20 # introns diverged within alignments (also within alignments which were separated) are written to warnings, they should be excluded from bait design

M_cluster       <- "ward.D"   # clustering method of hclust (average, ward.D,...)


##############################################################################################
################################### identify alignment limits ################################
##############################################################################################

FileNames <- choose.files(".fa") # select aligned sequences

data_out <- NULL
for(h in 1:length(FileNames)) {

file         <- read.csv(FileNames[h], stringsAsFactors = F, header=F)

len          <- nchar(FileNames[h])
queries      <- substring(FileNames[h],(nchar(getwd())+20),len-11)
N_loci_Malus <- unique(summary_statistics[,"N_loci"][summary_statistics[,"query_ID"]==substring(queries,1,12) & summary_statistics[,"genome"]=="Malus"])
N_loci_Pyrus <- unique(summary_statistics[,"N_loci"][summary_statistics[,"query_ID"]==substring(queries,1,12) & summary_statistics[,"genome"]=="Pyrus"])

seq_names       <- matrix(file[seq(1,nrow(file),2),],ncol=1)
seq_names_mRNA  <- substring(seq_names,nchar(seq_names)-3,nchar(seq_names))
N_mRNAs         <- length(seq_names_mRNA[,1][seq_names_mRNA[,1]=="mRNA"]) # number of mRNAs in the alignment

seq_names_Pyrus <- substring(seq_names,15,19)
N_Pyrus_seqs    <- length(seq_names_Pyrus[,1][seq_names_Pyrus[,1]=="Pyrus"]) # number of Pyrus subject sequences in the alignment

if(nrow(file) > 4) { # more than 2 sunject sequences are required for the trmming of the alignment, otherwise it is used unchanged (only mRNA and its subject are present)
    mean_d <- NULL
    for(i in 1:(nchar(file[2,])-size_window+1)) {
        mean_d_i <- mean(stringdistmatrix(toupper(substring(file[seq(2+N_mRNAs*2,nrow(file),2),],i,(i+size_window-1)))))/size_window
        mean_d   <- c(mean_d,mean_d_i)
    }
    a               <- which(mean_d < threshold_div_exons)
    start_alignment <- a[1]
    end_alignment   <- a[length(a)] + size_window-1
} else {
    start_alignment <- 1
    end_alignment   <- nchar(file[2,])
}

file[seq(2,nrow(file),2),] <- substring(file[seq(2,nrow(file),2),],start_alignment,end_alignment)


##############################################################################################
################################ identify intron_limits and exons limits #############################
##############################################################################################

# search for sequence which is most similar to the Malus mRNA, similarity based on exons only
mRNA <- file[2,]
mRNA_positions_without_gaps <- which(strsplit(mRNA,NULL)[[1]]!="-")

file_for_d <- NULL
for(i in seq(2,nrow(file),2)) {
    file_for_d_i <- matrix(strsplit(file[i,],NULL)[[1]],nrow=1)
    file_for_d_i <- paste(file_for_d_i[,mRNA_positions_without_gaps],collapse="")
    file_for_d   <- rbind(file_for_d, file_for_d_i)
}
d           <- as.matrix(stringdistmatrix(toupper(file_for_d)))
diag(d)     <- NA
seq_names   <- as.matrix(substr(file[seq(1,nrow(file),2),],2,nchar(file[seq(1,nrow(file),2),])),ncol=1)
d           <- as.data.frame(cbind(seq_names, d), stringsAsFactors = F)
colnames(d) <- "sequence"

min_dist            <- min(d[,2],na.rm=T) # we may better read in the distance from the summary file
name_Malus_ortholog <- subset(d,d[,2]==min_dist)[1,"sequence"]
row_Malus_ortholog  <- as.numeric(rownames(subset(d,d[,2]==min_dist)[1,]))*2
Malus_ortholog      <- file[row_Malus_ortholog,]
gaps_mRNA           <- as.matrix(which(strsplit(mRNA,NULL)[[1]]=="-"),ncol=1)
# gaps_Malus_ortholog <- as.matrix(which(strsplit(Malus_ortholog,NULL)[[1]]=="-"),ncol=1)
ident_positions     <- which(stringdist(toupper(strsplit(mRNA,NULL)[[1]]),toupper(strsplit(Malus_ortholog,NULL)[[1]]))==0)

if(nrow(gaps_mRNA) > 0) {
    if(nrow(gaps_mRNA) == 1) gap_limits <- c(gaps_mRNA,gaps_mRNA)
    if(nrow(gaps_mRNA)  > 1) {
        x          <- gaps_mRNA
        gap_limits <- x[1]
        for (i in 1:(length(x)-1)) {
           m          <- x[i]-x[i+1]
           if(m!=-1) {gap_limits_i <- c(x[i],x[i+1])
               gap_limits <- c(gap_limits,gap_limits_i)
           }
        }
    } # end if(nrow(gaps_mRNA)  > 1)
    gap_limits           <- c(gap_limits, x[length(x)])
    gap_limits           <- cbind(gap_limits[seq(1,length(gap_limits),2)],gap_limits[seq(2,length(gap_limits),2)])
    colnames(gap_limits) <- c("start_intron_limits","end_intron_limits")

    x            <- ident_positions
    ident_limits <- x[1]
    for (i in 1:(length(x)-1)) {
       m          <- x[i]-x[i+1]
       if(m!=-1) {ident_limits_i <- c(x[i],x[i+1])
           ident_limits <- c(ident_limits,ident_limits_i)
       }
    }
    ident_limits           <- c(ident_limits, x[length(x)])
    ident_limits           <- cbind(ident_limits[seq(1,length(ident_limits),2)],ident_limits[seq(2,length(ident_limits),2)])
    colnames(ident_limits) <- c("start_ident","end_ident")

    no_intron_limits <- NULL
    for(i in 1:nrow(gap_limits)) {
        for(j in 1:nrow(ident_limits)) {
            if(gap_limits[i,1] >= ident_limits[j,1] & gap_limits[i,2] <= ident_limits[j,2]) {no_intron_limits_i <- i; no_intron_limits <- c(no_intron_limits,no_intron_limits_i)}
        }
    }
    if(length(no_intron_limits)>0) {intron_limits <- gap_limits[-no_intron_limits,]} else {intron_limits <- gap_limits}
    if(length(intron_limits) ==2)  intron_limits  <- matrix(intron_limits,nrow=1)

    if(nrow(intron_limits) > 0) {
        exon_limits <- NULL
        if(intron_limits[1,1] == 1) { # alignment starts with intron
            if(intron_limits[nrow(intron_limits),2] == nchar(mRNA)) { # alignment starts with intron & alignment ends with intron
                for(j in 1:(nrow(intron_limits)-1)) {
                    exon_limits_1_j <- intron_limits[j,2]   + 1
                    exon_limits_2_j <- intron_limits[j+1,1] - 1
                    exon_limits     <- rbind(exon_limits,c(exon_limits_1_j,exon_limits_2_j))
                }
            } else {                                                  # alignment starts with intron & alignment ends with exon
                for(k in 1:nrow(intron_limits)) {
                    exon_limits_1_k <- intron_limits[k,2]   + 1
                    if(k!=nrow(intron_limits)) {exon_limits_2_k <- intron_limits[k+1,1] - 1} else {exon_limits_2_k <- nchar(mRNA)}
                    exon_limits     <- rbind(exon_limits,c(exon_limits_1_k,exon_limits_2_k))
                }
            }
        } else { # alignment starts with exon
            if(intron_limits[nrow(intron_limits),2] == nchar(mRNA)) { # alignment starts with exon & alignment ends with intron
                for(l in 1:nrow(intron_limits)) {
                    if(l==1) {exon_limits_1_l <- 1} else {exon_limits_1_l <- intron_limits[(l-1),2] + 1}
                    exon_limits_2_l <- intron_limits[l,1] - 1
                    exon_limits     <- rbind(exon_limits,c(exon_limits_1_l,exon_limits_2_l))
                }
            } else {                                                  # alignment starts with exon & alignment ends with exon
                for(m in 1:(nrow(intron_limits)+1)) {
                    if(m==1)                       {exon_limits_1_m <- 1} else {exon_limits_1_m <- intron_limits[(m-1),2] + 1}
                    if(m!=(nrow(intron_limits)+1)) {exon_limits_2_m <- intron_limits[m,1] - 1} else {exon_limits_2_m <- nchar(mRNA)}
                    exon_limits <- rbind(exon_limits,c(exon_limits_1_m,exon_limits_2_m))
                }
            }
        }
    } else {exon_limits <- matrix(c(1,nchar(mRNA)),nrow=1)} # no introns, i.e. only one exon
} else { # end if(nrow(gaps_mRNA) > 0)
exon_limits <- matrix(c(1,nchar(mRNA)),nrow=1); intron_limits  <- matrix(c(1,1),nrow=1)}

if(length(exon_limits)==2) exon_limits <- matrix(exon_limits,nrow=1)

intron_limits   <- subset(intron_limits, (intron_limits[,2] + 1 - intron_limits[,1]) >= min_length)
length_exon_min <- min(exon_limits[,2] + 1 - exon_limits[,1]); length_exon_max <- max(exon_limits[,2] + 1 - exon_limits[,1])


##############################################################################################
################## calculate sequence divergence for all exons and introns ###################
##############################################################################################

### Exons
max_exon_div         <- NULL
max_exon_realign_div <- NULL
mean_exon_div        <- NULL
length_exons_total   <- 0
seq_file2            <- NULL
warnings             <- NULL
for(i in 1:nrow(exon_limits)) {
    div_exons       <- stringdistmatrix(toupper(substr(file[seq(2,nrow(file),2),],exon_limits[i,1],exon_limits[i,2])))/(as.numeric(exon_limits[i,2]) + 1 - as.numeric(exon_limits[i,1]))
    max_div_exon_i  <- max(div_exons)
    mean_div_exon_i <- mean(div_exons) * (exon_limits[i,2]+1-exon_limits[i,1])/sum(exon_limits[,2]+1-exon_limits[,1])
    max_exon_div    <- c(max_exon_div, max_div_exon_i)
    mean_exon_div   <- c(mean_exon_div, mean_div_exon_i)
if(length(div_exons)>=3) {
    if(max_div_exon_i > max_div_exons) {
    fit        <- hclust(div_exons, method = M_cluster)
    fit$labels <- file[seq(1,nrow(file),2),]
    branches   <- cutree(fit, 2)
    file_partition           <- as.data.frame(cbind(names(branches),as.numeric(branches)),stringsAsFactors=F)
    colnames(file_partition) <- c("seq_name","group")
    sequ_order               <- as.data.frame(cbind(file[seq(1,nrow(file),2),],seq(1,nrow(file),2)),stringsAsFactors=F)
    colnames(sequ_order)     <- c("seq_name","order")
    m                        <- merge(file_partition,sequ_order)
    rows                             <- as.numeric(m[,"order"][m$group==1])
    rows                             <- sort(c(rows,(rows+1)))
    seq_file1                        <- as.matrix(substr(file[rows,],exon_limits[i,1],exon_limits[i,2]),ncol=1)
    seq_file1[seq(1,length(rows),2)] <- file[rows[seq(1,length(rows),2)],]
    rows                             <- as.numeric(m[,"order"][m$group==2])
    rows                             <- sort(c(rows,(rows+1)))
    seq_file2                        <- as.matrix(substr(file[rows,],exon_limits[i,1],exon_limits[i,2]),ncol=1)
    seq_file2[seq(1,length(rows),2)] <- file[rows[seq(1,length(rows),2)],]
    } else {seq_file1 <- as.matrix(substr(file[,1],exon_limits[i,1],exon_limits[i,2]),ncol=1)
            seq_file1[seq(1,nrow(file),2)] <- file[seq(1,nrow(file),2),]
            seq_file2 <- NULL}
} else { # the alignment only holds the mRNA and its corresponding intron
     seq_file1    <- as.matrix(substr(file[,],exon_limits[i,1],exon_limits[i,2]),ncol=1)
     seq_file1[seq(1,nrow(file),2)] <- file[seq(1,nrow(file),2),]
     seq_file2    <- NULL
}

file_name <- paste("aligned sequences/exons/",substring(FileNames[h],(nchar(getwd())+20),len-11),"_exon_",i,".fa",sep="")
write.table(seq_file1,file_name,row.names=FALSE,col.names=FALSE,quote=F)
file_name_selected <- paste("selected sequences/exons/",substr(file_name,25,nchar(file_name)),sep="")
if(is.null(seq_file1)==F & is.null(seq_file2)==T) {
    write.table(seq_file1,file_name_selected,row.names=FALSE,col.names=FALSE,quote=F)
}
## reaglin or remove gaps from alignments of diverged sequence groups (sequences have been partitioned to two aligments)
if(is.null(seq_file1)==F & is.null(seq_file2)==F) {
    for(j in seq(2,nrow(seq_file1),2)) {
        m                <- strsplit(seq_file1[j,],NULL)[[1]]
        coordinates_gaps <- which(m=="-")
        if(length(coordinates_gaps)>0) {seq_file1[j,] <- paste(m[-coordinates_gaps],collapse="")} else {seq_file1[j,] <- paste(m,collapse="")}
    }
    file_name_realigned <- paste("re",file_name,sep="")
    if(nrow(seq_file1)>2) {
    write.table(seq_file1,file_name_realigned,row.names=FALSE,col.names=FALSE,quote=F)

    dna                  <- readDNAStringSet(file_name_realigned)
    aligned_dna          <- AlignSeqs(dna)
    aligned_adjusted_dna <- AdjustAlignment(aligned_dna)

    write.table(aligned_adjusted_dna,file_name_realigned,col.names=FALSE,quote=F)
    seq_file1     <- read.csv(file_name_realigned, sep=" ",  header = F, stringsAsFactors = F)
    seq_file1[,1] <- paste(">",seq_file1[,1],sep="")
    seq_file1     <- c(as.character(seq_file1[,1]),as.character(seq_file1[,2]))
    m             <- (length(seq_file1)/2)
    sort_key <- NULL
    for(j in 1:m) {
        sort_key_j <- c(j,j+m)
        sort_key   <- c(sort_key,sort_key_j)
    }
    seq_file1 <- matrix(seq_file1[sort_key],ncol=1)
    div_exons_realigned    <- stringdistmatrix(toupper(seq_file1[seq(2,nrow(seq_file1),2),]))/nchar(seq_file1[2,1])
    max_div_exon_realign_i <- max(div_exons_realigned)
    max_exon_realign_div   <- c(max_exon_realign_div, max_div_exon_realign_i)
    } # if(nrow(seq_file1)>2)
write.table(seq_file1,file_name_realigned,row.names=FALSE,col.names=FALSE,quote=F)
write.table(seq_file1,file_name_selected, row.names=FALSE,col.names=FALSE,quote=F)
} # end if(is.null(seq_file1)==F & is.null(seq_file2)==F)
length_exons_total <- c(length_exons_total,nchar(seq_file1[2]))
if(is.na(nchar(seq_file1[2]))==F) {if(nchar(seq_file1[2])<38) {warning_h <- as.character(paste("Warning: sequences in file",file_name,"are too short!!")); warnings <- rbind(warnings,warning_h)}}

if(is.null(seq_file2)==F) {
    file_name <- paste(substring(file_name,1,nchar(file_name)-3),"_diverged_variant",".fa",sep="")
    write.table(seq_file2,file_name,row.names=FALSE,col.names=FALSE,quote=F)
    for(j in seq(2,nrow(seq_file2),2)) {
        m                <- strsplit(seq_file2[j,],NULL)[[1]]
        coordinates_gaps <- which(m=="-")
        if(length(coordinates_gaps)>0) {seq_file2[j,] <- paste(m[-coordinates_gaps],collapse="")} else {seq_file2[j,] <- paste(m,collapse="")}
   }
    file_name_realigned <- paste("re",file_name,sep="")
    if(nrow(seq_file2)>2) {
    write.table(seq_file2,file_name_realigned,row.names=FALSE,col.names=FALSE,quote=F)
    dna                  <- readDNAStringSet(file_name_realigned)
    aligned_dna          <- AlignSeqs(dna)
    aligned_adjusted_dna <- AdjustAlignment(aligned_dna)

    write.table(aligned_adjusted_dna,file_name_realigned,col.names=FALSE,quote=F)
    seq_file2     <- read.csv(file_name_realigned, sep=" ",  header = F, stringsAsFactors = F)
    seq_file2[,1] <- paste(">",seq_file2[,1],sep="")
    seq_file2     <- c(as.character(seq_file2[,1]),as.character(seq_file2[,2]))
    m             <- (length(seq_file2)/2)
    sort_key <- NULL
    for(j in 1:m) {
        sort_key_j <- c(j,j+m)
        sort_key   <- c(sort_key,sort_key_j)
    }
    seq_file2 <- matrix(seq_file2[sort_key],ncol=1)
    div_exons_realigned    <- stringdistmatrix(toupper(seq_file2[seq(2,nrow(seq_file2),2),]))/nchar(seq_file2[2,1])
    max_div_exon_realign_i <- max(div_exons_realigned)
    max_exon_realign_div   <- c(max_exon_realign_div, max_div_exon_realign_i)
    } # if(nrow(seq_file2)>2)
write.table(seq_file2,file_name_realigned,row.names=FALSE,col.names=FALSE,quote=F)
file_name_selected <- paste("selected sequences/exons/",substr(file_name,25,nchar(file_name)),sep="")
write.table(seq_file2,file_name_selected,row.names=FALSE,col.names=FALSE,quote=F)
length_exons_total <- c(length_exons_total,nchar(seq_file2[2]))
if(is.na(nchar(seq_file2[2]))==F) {if(nchar(seq_file2[2])<38) {warning_h <- as.character(paste("Warning: sequences in file",file_name,"are too short!!")); warnings <- rbind(warnings,warning_h)}}
} # end if(is.null(seq_file2)==F)
} # end exon-loop


### Introns
max_intron_div         <- NULL
max_intron_realign_div <- NULL
mean_intron_div        <- NULL
length_intron          <- NULL
length_introns_total   <- 0
seq_file1 <- seq_file2 <- NULL
counter                <- 0
N_subjects      <- (nrow(file) - N_mRNAs*2)/2 # introns are only present in the subject sequences, therefore the mRNAs are excluded
if(N_subjects>=2 & N_Pyrus_seqs>0 & nrow(intron_limits)>0) { # if there is only on subject sequence introns or no Pyrus subject sequences will no be extracted because divergence among introns cannot be (reliably) calculated
if(N_loci_Malus==1 | (N_loci_Malus==2 & N_subjects>2)) { # only groups with 2 sequences in case of 1 locus and 3 or more sequences in case of two lci remain (only for these homologous introns can be compared)
for(i in 1:nrow(intron_limits)) {
    div_introns       <- stringdistmatrix(toupper(substr(file[seq(2+N_mRNAs*2,nrow(file),2),],intron_limits[i,1],intron_limits[i,2])))/(as.numeric(intron_limits[i,2]) + 1 - as.numeric(intron_limits[i,1]))
if((max(div_introns) < max_allowed_div_introns & N_loci_Malus == 1) | N_loci_Malus == 2) { # for sequence groups with only one locus, homologous introns which are to much diverged are immediately excluded; for groups with two loci, the divergence is separately checked for the paralogous loci later
    max_div_intron_i  <- max(div_introns)
    counter           <- counter + 1
    mean_div_intron_i <- mean(div_introns) * (intron_limits[i,2]+1-intron_limits[i,1])/sum(intron_limits[,2]+1-intron_limits[,1])
    max_intron_div    <- c(max_intron_div, max_div_intron_i)
    mean_intron_div   <- c(mean_intron_div, mean_div_intron_i)
    length_intron_i   <- intron_limits[i,2] + 1 - intron_limits[i,1]
    length_intron     <- c(length_intron, length_intron_i)
    if(max_div_intron_i > max_div_introns) {
        fit        <- hclust(div_introns, method = M_cluster)
        fit$labels <- file[seq(1+N_mRNAs*2,nrow(file),2),]
        branches   <- cutree(fit, 2)
        file_partition           <- as.data.frame(cbind(names(branches),as.numeric(branches)),stringsAsFactors=F)
        colnames(file_partition) <- c("seq_name","group")
        sequ_order               <- as.data.frame(cbind(file[seq(1+N_mRNAs*2,nrow(file),2),],seq(1+N_mRNAs*2,nrow(file),2)),stringsAsFactors=F)
        colnames(sequ_order)     <- c("seq_name","order")
        m                        <- merge(file_partition,sequ_order)
        rows                     <- as.numeric(m[,"order"][m$group==1])
        rows                     <- sort(c(rows,(rows+1)))
        seq_file1                <- as.matrix(substr(file[rows,],intron_limits[i,1],intron_limits[i,2]),ncol=1)
        seq_file1[seq(1,length(rows),2)] <- file[rows[seq(1,length(rows),2)],]
        rows                     <- as.numeric(m[,"order"][m$group==2])
        rows                     <- sort(c(rows,(rows+1)))
        seq_file2                <- as.matrix(substr(file[rows,],intron_limits[i,1],intron_limits[i,2]),ncol=1)
        seq_file2[seq(1,length(rows),2)] <- file[rows[seq(1,length(rows),2)],]
    } else {seq_file1                             <- as.matrix(substr(file[seq(1+N_mRNAs*2,nrow(file),1),1],intron_limits[i,1],intron_limits[i,2]),ncol=1)
        seq_file1[seq(1,nrow(seq_file1),2),1] <- file[seq(1+N_mRNAs*2,nrow(file),2),]
        seq_file2 <- NULL
    }

## here alignments for which no homologous sequence comparisons are possible (only one unique subject sequence per locus i the alignment) are deleted; note that this is only necessary for sequence groups with two loci (for the groups with one locus this filter is already set before the introm loop starts); additionally note that for the 2 loci case subject sequences are always extracted twice (hits of both mRNAs!); see file comparability of subject sequences theoretical considerations.xlsx for expectations under differing numbers of Pyrus and Malus loci
                           if(N_loci_Malus==2 & N_mRNAs==1 & nrow(seq_file1) < 4)  seq_file1 <- NULL
if(is.null(seq_file2)==F) {if(N_loci_Malus==2 & N_mRNAs==1 & nrow(seq_file2) < 4) {seq_file2 <- NULL}}
if(is.null(seq_file1)==F) {if(N_loci_Malus==2 & N_mRNAs==2 & nrow(seq_file1) < 8) {seq_file1 <- NULL}}
if(is.null(seq_file2)==F) {if(N_loci_Malus==2 & N_mRNAs==2 & nrow(seq_file2) < 8) {seq_file2 <- NULL}}

file_name <- paste("aligned sequences/introns/",substring(FileNames[h],(nchar(getwd())+20),len-11),"_intron_",i,".fa",sep="")
if(is.null(seq_file1)==F) write.table(seq_file1,file_name,row.names=FALSE,col.names=FALSE,quote=F)
file_name_selected <- paste("selected sequences/introns/",substr(file_name,27,nchar(file_name)),sep="")
## reaglin or remove gaps from alignments of diverged sequence groups (sequences have been partitioned to two aligments)
if(is.null(seq_file1)==F & is.null(seq_file2)==T) {
    write.table(seq_file1,file_name_selected,row.names=FALSE,col.names=FALSE,quote=F)
}
if(is.null(seq_file1)==F & is.null(seq_file2)==F) {
    for(j in seq(2,nrow(seq_file1),2)) {
        m                <- strsplit(seq_file1[j,],NULL)[[1]]
        coordinates_gaps <- which(m=="-")
        if(length(coordinates_gaps)>0) {seq_file1[j,] <- paste(m[-coordinates_gaps],collapse="")} else {seq_file1[j,] <- paste(m,collapse="")}
    }
    file_name_realigned <- paste("re",file_name,sep="")
    if(nrow(seq_file1)>2) {
    write.table(seq_file1,file_name_realigned,row.names=FALSE,col.names=FALSE,quote=F)

    dna                  <- readDNAStringSet(file_name_realigned)
    aligned_dna          <- AlignSeqs(dna)
    aligned_adjusted_dna <- AdjustAlignment(aligned_dna)

    write.table(aligned_adjusted_dna,file_name_realigned,col.names=FALSE,quote=F)
    seq_file1     <- read.csv(file_name_realigned, sep=" ",  header = F, stringsAsFactors = F)
    seq_file1[,1] <- paste(">",seq_file1[,1],sep="")
    seq_file1     <- c(as.character(seq_file1[,1]),as.character(seq_file1[,2]))
    m             <- (length(seq_file1)/2)
    sort_key <- NULL
    for(j in 1:m) {
        sort_key_j <- c(j,j+m)
        sort_key   <- c(sort_key,sort_key_j)
    }
    seq_file1 <- matrix(seq_file1[sort_key],ncol=1)
    div_introns_realigned    <- stringdistmatrix(toupper(seq_file1[seq(2,nrow(seq_file1),2),]))/nchar(seq_file1[2,1])
    max_div_intron_realign_i <- max(div_introns_realigned)
    max_intron_realign_div   <- c(max_intron_realign_div, max_div_intron_realign_i)
    } # if(nrow(seq_file1)>2)
write.table(seq_file1,file_name_realigned,row.names=FALSE,col.names=FALSE,quote=F)
write.table(seq_file1,file_name_selected,row.names=FALSE,col.names=FALSE,quote=F)
} # end if(is.null(seq_file1)==F & is.null(seq_file2)==F)
length_introns_total <- c(length_introns_total,nchar(seq_file1[2]))
if(is.null(seq_file1)==F) {if(nchar(seq_file1[2])<38) {warning_h <- as.character(paste("Warning: sequences in file",file_name,"are too short!!")); warnings <- rbind(warnings,warning_h)}}

if(is.null(seq_file2)==F) {
    file_name <- paste(substring(file_name,1,nchar(file_name)-3),"_diverged_variant",".fa",sep="")
    write.table(seq_file2,file_name,row.names=FALSE,col.names=FALSE,quote=F)
    for(j in seq(2,nrow(seq_file2),2)) {
        m                <- strsplit(seq_file2[j,],NULL)[[1]]
        coordinates_gaps <- which(m=="-")
        if(length(coordinates_gaps)>0) {seq_file2[j,] <- paste(m[-coordinates_gaps],collapse="")} else {seq_file2[j,] <- paste(m,collapse="")}
   }
    file_name_realigned <- paste("re",file_name,sep="")
    if(nrow(seq_file2)>2) {
    write.table(seq_file2,file_name_realigned,row.names=FALSE,col.names=FALSE,quote=F)
    dna                  <- readDNAStringSet(file_name_realigned)
    aligned_dna          <- AlignSeqs(dna)
    aligned_adjusted_dna <- AdjustAlignment(aligned_dna)

    write.table(aligned_adjusted_dna,file_name_realigned,col.names=FALSE,quote=F)
    seq_file2     <- read.csv(file_name_realigned, sep=" ",  header = F, stringsAsFactors = F)
    seq_file2[,1] <- paste(">",seq_file2[,1],sep="")
    seq_file2     <- c(as.character(seq_file2[,1]),as.character(seq_file2[,2]))
    m             <- (length(seq_file2)/2)
    sort_key <- NULL
    for(j in 1:m) {
        sort_key_j <- c(j,j+m)
        sort_key   <- c(sort_key,sort_key_j)
    }
    seq_file2 <- matrix(seq_file2[sort_key],ncol=1)
    div_introns_realigned    <- stringdistmatrix(toupper(seq_file2[seq(2,nrow(seq_file2),2),]))/nchar(seq_file2[2,1])
    max_div_intron_realign_i <- max(div_introns_realigned)
    max_intron_realign_div   <- c(max_intron_realign_div, max_div_intron_realign_i)
    } # if(nrow(seq_file2)>2)
write.table(seq_file2,file_name_realigned,row.names=FALSE,col.names=FALSE,quote=F)
file_name_selected <- paste("selected sequences/introns/",substr(file_name,27,nchar(file_name)),sep="")
write.table(seq_file2,file_name_selected,row.names=FALSE,col.names=FALSE,quote=F)
length_introns_total <- c(length_introns_total,nchar(seq_file2[2]))
if(is.na(nchar(seq_file2[2]))==F) {if(nchar(seq_file2[2])<38) {warning_h <-  as.character(paste("Warning: sequences in file",file_name,"are too short!!")); warnings <- rbind(warnings,warning_h)}}
} # end if(is.null(seq_file2)==F)
} # end if((max(div_introns) < max_allowed_div_introns & N_loci_Malus == 1) | N_loci_Malus == 2)
} # end intron-loop
} # end if(N_loci_Malus == 1 | (N_loci_Malus == 2 & N_subjects>2))
} # end if(N_subjects>1 & N_Pyrus_seqs>0 & is.null(intron_limits)!=T)

max_exon_div       <- round(max(max_exon_div,na.rm=T),2)
mean_exon_div      <- round(sum(mean_exon_div,na.rm=T),2)
if(is.null(max_intron_div)==T)         {max_intron_div <- mean_intron_div <- NA} else {max_intron_div         <- round(max(max_intron_div,na.rm=T),2); mean_intron_div <- round(sum(mean_intron_div,na.rm=T),2)}
if(is.null(max_exon_realign_div)==T)   {max_exon_realign_div <- NA}              else {max_exon_realign_div   <- round(max(max_exon_realign_div,  na.rm=T),2)}
if(is.null(max_intron_realign_div)==T) {max_intron_realign_div <- NA}            else {max_intron_realign_div <- round(max(max_intron_realign_div,na.rm=T),2)}
if(is.null(length_intron) ==T) {length_intron_min <- length_intron_max <- NA} else {length_intron_min  <- min(length_intron, na.rm=T); length_intron_max  <- max(length_intron, na.rm=T)}
nrow_introns       <- counter
if(N_Pyrus_seqs==0) {with_Pyrus_subjects <- "No"} else {with_Pyrus_subjects <- "Yes"}
if(length(N_loci_Pyrus)==0) N_loci_Pyrus <- 0

data_out_h <- c(queries,N_loci_Malus,N_loci_Pyrus,N_mRNAs,N_subjects,with_Pyrus_subjects,nrow(exon_limits),nrow_introns,max_exon_div,mean_exon_div,max_intron_div,mean_intron_div,max_exon_realign_div,max_intron_realign_div,length_exon_min,length_exon_max,length_intron_min,length_intron_max,sum(length_exons_total,na.rm=T),sum(length_introns_total,na.rm=T))
data_out   <- rbind(data_out,data_out_h)
} # end File or h-loop

settings_header <- c("N seq groups","min length introns","size sliding window","threshold div alignment","max div exons","max div introns","max allowed divergence introns","cluster method") 
settings_values <- c(h,min_length,size_window,threshold_div_exons,max_div_exons,max_div_introns,max_allowed_div_introns,M_cluster)
settings        <- as.data.frame(cbind(settings_header,settings_values))

write.table(settings,   paste("settings run ",Sys.Date(),".txt",sep=""),sep="\t",row.names=FALSE,quote=F)
                        colnames(data_out) <- c("queries","N loci Malus","N loci Pyrus","NmRNAs","N subject sequences","Pyrus subjects present","N exons","N introns","max div exons","mean div exons","max div introns","mean div introns","max exon realign div","max intron realign div","length exon min","length exon max","length intron min","length intron max","total length exons","total length introns")
write.table(data_out,   paste("statistics Malus and Pyrus run ",Sys.Date(),".txt",sep=""),sep="\t",dec=",",row.names=FALSE,quote=F)

x <- paste(as.character(subset(data_out, data_out[,"N introns"]       > 0 & data_out[,"N loci Malus"]==1 & data_out[,"max div introns"]        > max_allowed_div_introns)[,"queries"]),collapse=",") # actually these cases are already excluded at the start of the intron loop
y <- paste(as.character(subset(data_out, data_out[,"N introns"]       > 0 & data_out[,"N loci Malus"]==2 & data_out[,"max intron realign div"] > max_allowed_div_introns)[,"queries"]),collapse=",")
z <- paste(as.character(subset(data_out, as.numeric(data_out[,"length exon min"]) < min_length)[,"queries"]),collapse=",")

if(nchar(x)>0) warnings <- rbind(warnings,paste("groups with single loci exceeding the allowed divergenve of introns",x))
if(nchar(y)>0) warnings <- rbind(warnings,paste("groups with two loci exceeding the allowed divergenve of introns",y))
if(nchar(z)>0) warnings <- rbind(warnings,paste("groups with exons shorter than the required length of baits (min_length)",z))
write.table(warnings,paste("warnings run ",Sys.Date(),".txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=F)

m                <- aggregate(summary_statistics[,"subject_start"], by=list(summary_statistics$genome, summary_statistics$query_ID, summary_statistics$chromosome),min)
colnames(m)      <- c("genome","queries","chromosome","subject start")
n                <- as.data.frame(data_out[,c("queries","N loci Malus","N loci Pyrus","max div exons","mean div exons","max div introns","mean div introns","max exon realign div","max intron realign div")],stringsAsFactors=F)
o                <- subset(n, nchar(n[,"queries"])==25)
o[,"queries"]    <- substring(o[,"queries"],14,25)
n[,"queries"][nchar(n[,"queries"])==25] <- substring(n[,"queries"][nchar(n[,"queries"])==25],1,12)
n                <- rbind(n,o)
rownames(n)      <- 1:nrow(n)
n                <- merge(n,m)
genomic_position <- n[,c("queries","genome","chromosome","subject start","N loci Malus","N loci Pyrus","max div exons","mean div exons","max div introns","mean div introns","max exon realign div","max intron realign div")]
genomic_position <- genomic_position[with(genomic_position, order(genomic_position[,"subject start"])), ]

chro_Malus <- sort(unique(genomic_position[,"chromosome"][genomic_position[,"genome"]=="Malus"]))[-1]
m          <- NULL
for(i in chro_Malus) {
    queries_chro <- unique(genomic_position[,"queries"][genomic_position[,"chromosome"]==i & genomic_position[,"N loci Malus"]==1])
    m_i          <- subset(genomic_position, queries %in% queries_chro & genomic_position[,"N loci Pyrus"]==1)
    m_i          <- m_i[with(m_i, order(m_i[,"queries"],m_i[,"genome"])),]
    m            <- rbind(m, m_i)
}
genomic_position_single_copy_loci <- m

write.table(genomic_position_single_copy_loci, paste("statistics Malus and Pyrus run",Sys.Date(),"sorted genomic position of single copy loci.txt",sep=" "),sep="\t",dec=",",row.names=FALSE,quote=F)
