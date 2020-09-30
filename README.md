Locus and sequence selection for bait design

The pipeline consists of the following R scripts which have to be consecutively excecuted. In some instances external files have to be provided (tsvs, genome sequences, order of chromosomes in the Pyrus genome).

It has been created within the Malinae project. Baits can be developed from mRNAs and modified Kew probes, respectively.
More details are provided within the scripts. See also strategy probe design v2.docx and the log.txt-files.
Genomes used were from the GDR: Malus x domestica GDDH13 V1.1 genome and Pyrus bretschneideri ‘DangshanSuli’ Genome Assembly v1.1 and finally, Pyrus communis Bartlett DH Genome v2.0 genome

Format of input files according to examples.

1. Transcripts are to be blasted against the Malus genome using the nBlast search of the GDR site. Output is a tsv file the columns of which need to be labeled. The tsv is read in by the script. The script 01 mRNA selection for Blast search.R is optional. It allows for selection of mRNAs based on length criteria of mRNAs and the area hit in the Malus genome.

2. 02 probe mapping.R. The script filters the transcripts (= queries) based on a host of settings, saves tables describing the properties of the queries as well as the selected queries.

3. Blast the selected queries against the Pyrus genome and feed them along with the obtained tsv-file to the script 03 probe mapping Pyrus.R. The script is similar to the former one and estimates the properties or suitability of the queries in the Pyrus genome.

4. The script 04 join statistics fully compatible genomes.R synthesizes the results of the searches/results in Malus and Pyrus.

5. Blast the finally selected queries against each other using the nBlast search of the NCBI site and save the hit table (choose the tab-delimited format). The script 05 grouping of query sequences.R. identifies paralogous relationships among mRNAs and checks whether mRNAs which are member of a paralogous pair had also two hits in the blast search against the Malus genome. Particularly important is the file with the ending "sequence coordinates.txt" which provides the coordinates for the extraction of Malus and Pyrus subject sequences. Inversion of subject sequences is noted and considered.

6. The script 06 loci for subject extraction.R actually makes the files holding the required information for the extraction of subject sequences from the genomes, saved in the output files _loci_for_subject_extraction.txt. The script requires the file order of chromosomes in Pyrus bretschneideri genome file.txt.

6. The _loci_for_subject_extraction.txt-output files need to be copied together with the genomes of Malus and Pyrus on the Linux GDDH server (directory /home/uni08/cdobes/pro/HybSeq/analysis ). Run the bash files loopMalus.sh and loopPyrus.sh (or loopMalusKewExons.sh). The Malus genome is saved in GDDH13_1-1_formatted.tab, the Pyrus genome in pbr.v1.1.chr.tab. Output of the bash files is the file XY_subject_sequences.txt. Copy this output file back to the working directory on the Windows system.

7. The script 07 make fasta files.R groups and aligns the subject as well as the query sequences. It requires the package Decipher. Before running the script some directories have to be created within the working directory as specified in the head of the script.

8. Finally, the 08 script sequence comparison.R separates the alignments into exonic and intronic regions, calculates pairwise distances among the sequences of each exon and intron and saves sequences in separate alignments if diverged beyond a set threshold. It also trims the ends of alignments and realigns the separated alignments. It provides a detailed statistics. Some specifics: The parameters "max div exons","mean div exons","max div introns","mean div introns" refer to sequence divergence among both exons and introns, i.e. in situations with only 1 locus they are the divergence among orthologs, in case of 2 loci among paralogs and orthologs. The parameters "max exon realign div","max intron realign div" have been separately calculated for the two alignments holding diverged sequences (max_div_exons/max_div_introns exceeded). Since in such cases usually paralogous sequences are separated from each other (orthologs are less diverged than paralogs) the realigned divergence refers to divergence among homologous sequences (in case of diverged sequences from groups with 1 locus the separated alignments each hold only one intron sequence, i.e. sequence divergence cannot be calculated). total length exons and total length introns are the sum of characters of the alignments (separated or not); i.e. are not the some of exons/introns in a single sequence. The alignments for bait design are saved in the subdirectories of directory \files for bait design\selected sequences. The script requires the packages stringdist and plyr.

9. The script 09 control.R aligns the first Malus mRNA and first subject sequence with the exons and introns. The exon and intron alignments (from the selected directory) have to be copied into a directory "control". Remove all gap symbols (-) from the aligments.

10. The script 10 final select.R allows to exclude loci based on some criteria like min exon lengths required or maximal divergence among sequences. It saves the selected exons/introns in the directory finally selected sequences.

11. The 11 script separate still highly diverged sequences.R reads in alignments from the directory finally selected sequences, calculated divergence among sequences and separates sequences if diverged above a set threshold (which can be different for exons and introns).
