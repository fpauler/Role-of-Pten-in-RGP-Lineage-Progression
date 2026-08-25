# analyses read counts in genomic intervals using bedtools

#input bam files comes from pipe - usually with find * -name "*.Aligned.sortedByCoord.out.bam"

while (<STDIN>) {
	chomp($_);

	@parts = split(/\//, $_);
	$folder = "$parts[0]";
	for ($i=1; $i<$#parts; $i++) {
		$folder = $folder."/".$parts[$i];
	}

	$cmd = "bedtools intersect -c -b $_ -a del_count/Pten_ex5.bed -split -bed > $folder/Pten_del.counts.out";
	print "$cmd\n";
	system($cmd);
}



