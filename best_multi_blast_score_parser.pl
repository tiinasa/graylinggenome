#!usr/bin/perl

#This script takes a number of tabular blast outputs as input files (param 1...n) and creates a table with  the best score, number of best hits, input file name, and the blast line(s) of each query sequence as output for each query sequence. 


#run: perl best_multi_blast_score_parser.pl repet_graylingrep.word7-blastn.txt swissprot_graylingrep.word7-blastx.txt > best_multi_blast_score_grayling.txt

#blast columns (1-based):
#
#1. 	 qseqid 	 query (e.g., gene) sequence id
#2. 	 sseqid 	 subject (e.g., reference genome) sequence id
#3. 	 pident 	 percentage of identical matches
#4. 	 length 	 alignment length
#5. 	 mismatch 	 number of mismatches
#6. 	 gapopen 	 number of gap openings
#7. 	 qstart 	 start of alignment in query
#8. 	 qend 	 end of alignment in query
#9. 	 sstart 	 start of alignment in subject
#10. 	 send 	 end of alignment in subject
#11. 	 evalue 	 expect value
#12. 	 bitscore 	 bit score

my @inputfiles=@ARGV;
my %bestscores=(); #key: query name, value: score
my %bestscoringlines=(); #key: query name, value: array of lines

foreach my $filename (@inputfiles) {
	print STDERR "Reading file $filename... ";
	
	open(BLAST,$filename) || die("file $filename not found"); 
	while(my $line=<BLAST>) {
		$line=~s/\s*$//;
		my @cols=split("\t",$line); 
		#if best score not found so far or this score is better than the previous ones, save the score and make a new line array where the line is stored
		my $out=$filename."\t".$line;

		if(!exists($bestscores{$cols[0]}) || $bestscores{$cols[0]}<$cols[11]) {
			my @linearray=($out);
			$bestscores{$cols[0]}=$cols[11];
			$bestscoringlines{$cols[0]}=\@linearray;
						print STDERR "Made new array with len ".@{$bestscoringlines{$cols[0]}}."\n";

			my $newarraylen=@{$bestscoringlines{$cols[0]}}+0;
			print STDERR "New best score for query ".$cols[0]." with score ".$cols[11]."\n, array len for ".$cols[0]." now $newarraylen\nline: $out\n";
			
		}
		
		elsif($bestscores{$cols[0]}==$cols[11]) { 
		
			my $newarraylen=@{$bestscoringlines{$cols[0]}}+0;
			push(@{$bestscoringlines{$cols[0]}},$out); print STDERR "Added 1 array element with score ".$cols[11].", array len for ".$cols[0]." now $newarraylen\nline: $out\n"; 
			
		}
		else {next; }
	}
	print STDERR "Done! Wrote scores to ".keys(%bestscores)." hits\n";
}
close BLAST; 

#print best hits
foreach my $queryname (keys(%bestscores)) {
	my @lines=@{$bestscoringlines{$queryname}};
	foreach my $line(@lines) {
		print $queryname."\t".$bestscores{$queryname}."\t".@lines."\t".$line."\n";
	}
}

