# command to calculate the coverage of the deleted region for QC
# uses bedtools as this allows to deal with split reads!
# the bed file with coordinated of the deleted region is hardcoded in the perl script
#this command puts a bed file in each analysis folder with the counts of reads in the given region

find * -name "*.Aligned.sortedByCoord.out.bam" | perl del_count/count_reads.pl


