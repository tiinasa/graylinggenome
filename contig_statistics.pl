#/usr/bin/perl -w
use strict;

#input: contigs in fasta or fastq format
#output: basic statistics of the sequence
#to run: perl contig_statistics.pl input.fa 
#Extra long chromosomes can cause error message "Substitution loop" while removing whitespaces


my ($len,$total)=(0,0);
my @x;
my ($min,$max)=(1000000000,0);
print("\nStatistics for $ARGV[0]:\n\n");
print("\nmaxrows: ");
my $isQuality=1;
my $row=0;

#print("lengths: ");
while(<>){
	$row++;
	if(/^[\+]/){			#if quality header line then flag=1, do nothing
		$isQuality=1;
	}	
	if(/^(\>|@)/){			#if header line then
		#print("head $_");
		$isQuality=0;
		if($len>0){
			$total+=$len;	#add length to total
			push @x,$len;	#save $len in array @x
			if($len<$min) { $min=$len; }
			if($len>$max) { $max=$len; print"$row-"; print "$max "; }
			#print "$len ";
		}
		$len=0;				#nollaa
	}
	else{
		if($isQuality==0) {
			s/\s*$//;
			if(m/\s/) {
				s/\s+//g;			
			}
			$len+=length($_);		#if sequence line then add length to $len
			
		}
	}
}
if ($len>0){				#save also the last data
	$total+=$len;
	push @x,$len;
	#print("$len\n");
	if($len<$min) { $min=$len; }
	if($len>$max) { $max=$len;  }

}
@x=sort{$b<=>$a} @x; 		#sort lengths 
my ($count,$half)=(0,0);
for (my $j=0;$j<@x;$j++){
	$count+=$x[$j];
	if (($count>=$total/2)&&($half==0)){
		print "\nN50: $x[$j]\n";
		print "L50: ".($j+1)."\n";
		$half=$x[$j]
	}elsif ($count>=$total*0.9){
		print "N90: $x[$j]\n";
		print "L90: ".($j+1)."\n";
		last;
	}
}
my $seqcount=@x;
my $ave=$total/$seqcount;

print("total length: $total\n");
print("count: $seqcount\n");
print("Average: $ave\n");
my $sumsq;
my $item;
my $dev;
my $sd;
foreach $item (@x) {

	$dev=($item-$ave);
	#print("dev $item - $ave = $dev squared: ");
	$dev=$dev**2;
	$sumsq=$sumsq+$dev;
	#print("$dev\n");
}
$sd=sqrt($sumsq/$seqcount);
print("sd: $sd\n");
print("min: $min\n");
print("max: $max\n");
