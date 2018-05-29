#!/usr/bin/perl
use strict;
use Getopt::Long;
my $VERSION=5.0;
#Last edited 19.12.2012

#example: perl barcodesplitter5.1.pl -i single -b barcodes.txt -s shorttest1.fq -d barcodeoutdir_shorttest -m 2 -t 2 -r TGCA
my $helptext=("

barcodesplitter$VERSION.pl
---------------------

use: perl barcodesplitter$VERSION.pl -b barcodefilename -i inputType -s inputfile1 (-p inputfile2) -d newOutDirectory (-m numberOfAllowedambiguousbases) (-t numberOfActionType) (-r restrictionsite)
help: perl barcodesplitter$VERSION.pl -h

input parameters:
-b string	barcodes file (tab delimited, id:s on 1st row and barcode 
		sequences on 2nd row. Reguired.
-i string	Input read type, 'single' or 'paired'. Reguired.
-s string	Single reads to be split in a fq file (upper case, N is a 
		special character describing any base).
-p string	In case of paired reads, also give pair file. IMPORTANT: Paired 
		reads are not assumed to contain barcodes in this version of Barcode-
		Splitter. 
are
		Input read file -s (and -p in case of paired) reguired.
-d string	output directory (must not exist yet). Reguired.
-m integer	allowed ambiguousbases (N bases) in barcodes (integer>=0, default=1).
            If restriction site exists, also include N:s in that sequence.
-e integer  allowed sequencing errors in sequence barcodes (integer>=0, default=0)
            If restriction site exists, also include N:s in that sequence.
-t integer	action type (integer 0-3, default=0):
		0: Barcode sizes are equalized with N:s at the end of barcode 
		sequences shorter than maximum barcode length. Barcode-matched 
		reads are outputted to one file, pair sequances to another file 
		and unmatched reads and pairs are outputted to their own files.		
		1: Barcodes are removed from -s reads and and barcode-matched 
		reads outputted to separate pairs of files. Unmatched pairs are 
		outputted to two unmatched files. To leave barcodes in place, add flag
                -l.
		2: Unrecognized barcodes are masked to N:s and all reads outputted 
		to single file. 
		IMPORTANT: This option currently doesn't have paired option! 
		3: same as (2) but also adds sampleid:s and barcodes to the beginning
		of headers and removes barcodes from sequences. (changed: v. 5.1)
		IMPORTANT: This option currently doesn't have paired option!
-l              If using -t action type 1, leave barcodes in reads. (added: v. 5.3)
-h 		Print help text and exit.
-r str	If given, then barcodesplitter also checks that restriction  
		site (str) is present after barcode, but doesn't cut that site 
		away from the output sequence (added v. 5.0)

");

my $barcodefile;
my $actionType=0; 			#default: action type 0
our $allowedambiguousbases=1; 	#default: 2 ambiguousbases allowed in barcodes
our $allowedseqerrors=0;  #allowed sequencing errors in barcodes (default: 0)
my $inputType;
my $readfile;
my $pairfile;
my $outputdir;
my $help;
my $restrictionsite;
my $leavebc=0;
GetOptions(	"b:s"=> \$barcodefile,
			"i:s"=> \$inputType,
			"s:s"=> \$readfile,
			"p:s"=> \$pairfile,
			"d:s"=> \$outputdir,
			"m:i"=> \$allowedambiguousbases,
			"e:i"=> \$allowedseqerrors,
			"t:i"=> \$actionType,
			"l"=> \$leavebc,
			"h"=> \$help,
			"r:s"=>\$restrictionsite
); 

my $readfileparsed=$readfile;  #remove paths from readfiles to use as output prefix
my $pairfileparsed=$pairfile;
$readfileparsed=~s/.*\///g;
$pairfileparsed=~s/.*\///g;
if($help) { print $helptext; exit; }

#check the validity of inputType and existence of read file(s)
if($inputType ne "single" & $inputType ne "paired") { die("Give valid input file type ('single' or 'paired' reads) as param -i!");}
if($readfile eq "") { die("Give read file as param -s");}
if($inputType eq "paired") { if($pairfile eq "") {die("Give pair read file as param -p");}}

mkdir($outputdir, 0777) or die("can't create outputdir, give output directory as param -d and check that the output directory doesn't exist before running the script!");

#barcode file is tab delimited file, id:s in first column, barcodes in second column
#reads barcode file into 
open(BARCODES,$barcodefile) or die("Can't open barcodefile. Give barcodefile as param -b");
my @row;
print("Reading barcodes...\nid\tbarcode\tfilename\tlength\n--\t-------\n");

my @ids;
my @barcodes;
my $filenames;
my @pairnames; #in case of paired reads use also these filenames
my $tempfile;
my @barcoderegexps;
my @barcodelengths;
my @barcodesubtractedlengths; #maxlength-barcode length
my @barcodesEqualized; #equalized length barcodes
my @filenames;
my @FILEHANDLES;
my @PAIRHANDLES;
my $maxBarcodeLength=0;				#after BARCODES file read this variable describes minimum barcode length
my $minBarcodeLength=1000000;		#after BARCODES file read this variable describes maximum barcode length
my $tempregexp;
my $tempLength;
my $tempequalbarcode;
my $bctotal=0;  #barcode counter (after BARCODES file read, this variable is set to barcode count -1 
my $restrictionsitelength=length($restrictionsite);
my @matchedreadsperbarcode;

#Read barcode file and 
while(my $line=<BARCODES>) {
	$line=~s/\s+$//;  #remove 
	if($line eq "") {next; }
	@row=split("\t",$line);
	
	push(@barcodes,@row[1]);
	push(@matchedreadsperbarcode,0);
	$tempregexp=@row[1];
	$tempregexp=~s/[A-Z]{1}/($&|N)/g;
	push(@barcoderegexps,$tempregexp); 
	
	push(@ids,@row[0]);
	
	push(@filenames,"$outputdir/@row[0]$readfileparsed"); #in case of single or paired reads 
	if($inputType=="paired") { 
		push(@pairnames,"$outputdir/@row[0]$pairfileparsed");  #in case of paired reads only
	} #in case of paired reads
	$tempLength=length(@row[1]);
	if($tempLength<$minBarcodeLength) { $minBarcodeLength=$tempLength; }
	if($tempLength>$maxBarcodeLength) { $maxBarcodeLength=$tempLength; }
	push(@barcodelengths,$tempLength);
	print("@ids[$bctotal]\t@barcodes[$bctotal]\t@filenames[$bctotal]\t@pairnames[$bctotal]\t@barcodelengths[$bctotal]\n");
	$bctotal++;

}#endwhile
close(BARCODES);

#length-equalized barcodes (actiontype 0, 2 and 3)
my $ii;
my $tempstring;
if($actionType==0||$actionType==2 ||$actionType==3) {
	print("Barcode lengths are equalized (add N:s to the end of shorter barcodes)\nif\tequalized barcode\n");
	for($ii=0;$ii<$bctotal;$ii++) { 				#if barcode lengths equalised then add this many N:s to the end of barcode
		@barcodesubtractedlengths[$ii]=$maxBarcodeLength-@barcodelengths[$ii];
		$tempstring="N" x @barcodesubtractedlengths[$ii];
		@barcodesEqualized[$ii]="@barcodes[$ii]$tempstring";
		print "@ids[$bctotal]\t@barcodesEqualized[$ii]\n";
	}#endfor 
}#endif
print("$bctotal barcodes read, min length $minBarcodeLength and max length $maxBarcodeLength.");

#If separate output files for each barcode then create filehandles for each barcode file
if($actionType==1) {
	my $iii;
	my $HANDLE;
	for($iii=0;$iii<$bctotal;$iii++) {
		local *HANDLE; #LOCALIZE FILEHANDLE SO IT BELONGS ONLY TO THIS LOOP!
		open(HANDLE,">>@filenames[$iii]");
		@FILEHANDLES[$iii]=*HANDLE;
	}#endfor
}#endif
#If separate output files for each barcode then create filehandles for each barcode file for paired reads
if($actionType==1 & $inputType eq "paired") {
	my $iii;
	my $HANDLE;
	for($iii=0;$iii<$bctotal;$iii++) {
		local *HANDLE; #LOCALIZE FILEHANDLE SO IT BELONGS ONLY TO THIS LOOP!
		open(HANDLE,">>@pairnames[$iii]");
		@PAIRHANDLES[$iii]=*HANDLE;
	}#endfor
}#endif


###############################
#sort reads to correct files
print("\rSorting reads...\n"); 


open(READS,$readfile) or die ("can't open read file. Give input read file as param -s");
my ($header,$seq,$qualheader,$qual);
my $barcodematched=0;
my $j;
my $ambiguousbases;
my $readcounter=0;
my $matchcounter=0; #counts barcode-identified reads or read pairs
my $mismatchcounter=0;  #counts number of reads or read pairs thrown to unmatched reads file because too many N:s or seqerrors 
						#found in barcode region

my($pairheader,$pairseq,$pairqualheader,$pairqual); 
if($inputType eq "paired") { open(PAIRS,$pairfile) or die ("can't open pair file. Give input read file as param -p"); }

####################
#if actionType==1, 
####################
#OUPUT READS TO SEPARATE FILES ACCORDING TO DIFFERENT BARCODES
if($actionType==1 & $inputType eq 'single') {
	my $FILEHANDLE;
	open(OUTFILE,">>$outputdir/unmatched_reads.fq");
	my $OUTFILE_UN=*OUTFILE;
	print("Use separate output files.\n");
	while(<READS>) {
		if($_ eq "\n") { next; }  #discard empty rows at the end
		$readcounter++;
		$header=$_;
		$seq=<READS>;
		$qualheader=<READS>;
		$qual=<READS>;
		
		$barcodematched=0;
		$FILEHANDLE=$OUTFILE_UN;
		#go through every barcode in @barcodes
		for($j=0;$j<$bctotal;$j++) { 
		#if barcode matches sequence, then  print to that file
		#print("test bc @barcodes[$j]...\n");
			if(matchbarcode($seq,@barcodes[$j])==1) { 
				$matchcounter++;
				$barcodematched=1;
				@matchedreadsperbarcode[$j]++;
				if($leavebc==0) {
					$seq=~s/^.{@barcodelengths[$j]}//;		#remove barcode
					$qual=~s/^\S{@barcodelengths[$j]}//;
				}
				$FILEHANDLE=@FILEHANDLES[$j];
				last; #match already found; stop searching matches for this read
			}
		}#endfor (barcodes list search)
		if($barcodematched==0) { $mismatchcounter++; }
		print $FILEHANDLE "$header$seq$qualheader$qual";
	}#endwhile
	close(READS);	

	my $iii;
	for($iii=0;$iii<$bctotal;$iii++) {
		close @FILEHANDLES[$iii];
	}#endfor
	close $OUTFILE_UN;
}#end if $actionType=1 & $inputType='single'


if($actionType==1 & $inputType eq 'paired') {
	my $FILEHANDLE1;
	my $FILEHANDLE2;
	open(OUTFILE1,">>$outputdir/unmatched_reads1.fq");
	open(OUTFILE2,">>$outputdir/unmatched_reads2.fq");
	my $OUTFILE_UN1=*OUTFILE1;
	my $OUTFILE_UN2=*OUTFILE2;
	print("Use separate output files for both pair files.\n");
	while(<READS>) {
		$readcounter++;
		$header=$_;
		$seq=<READS>;
		$qualheader=<READS>;
		$qual=<READS>;
		$pairheader=<PAIRS>;
		$pairseq=<PAIRS>;
		$pairqualheader=<PAIRS>;
		$pairqual=<PAIRS>;
		
		$barcodematched=0;
		$FILEHANDLE1=$OUTFILE_UN1;
		$FILEHANDLE2=$OUTFILE_UN2;
		#go through every barcode
		for($j=0;$j<$bctotal;$j++) { 
			#if barcode matches sequence, then  print to that file
			#print("test bc @barcodes[$j]...\n");
			if(matchbarcode($seq,@barcodes[$j])==1) { 
				$matchcounter++;
				$barcodematched=1;
				@matchedreadsperbarcode[$j]++;
				if($leavebc==0) {
					$seq=~s/^\S{@barcodelengths[$j]}//;		#remove barcode
					$qual=~s/^\S{@barcodelengths[$j]}//;
				}
				$FILEHANDLE1=@FILEHANDLES[$j];
				$FILEHANDLE2=@PAIRHANDLES[$j]; 
				last; #match already found; stop searching matches for this read
			}
		}#endfor (barcodes list search)
		if($barcodematched==0) { $mismatchcounter++; }
		print $FILEHANDLE1 "$header$seq$qualheader$qual";
		print $FILEHANDLE2 "$pairheader$pairseq$pairqualheader$pairqual";		
	}#endwhile
	close(READS);	
	close(PAIRS);
	
	my $iii;
	for($iii=0;$iii<$bctotal;$iii++) {
		close @FILEHANDLES[$iii];
		close @PAIRHANDLES[$iii];
	}#endfor
	close $OUTFILE_UN1;
	close $OUTFILE_UN2;

}#end if $actionType=1 & $inputType='paired'


############################
#IF SEPARATEFILES=0
############################
#IN CASE OF ONE OUTPUT FILE & SIZE EQUALIZATION
if($actionType==0 & $inputType eq "single") {
	my $FILEHANDLE;
	open(OUTFILE_UN,">>$outputdir/unmatched_reads.fq") || die("can't open outfile unmatched_reads.fq");
	my $OUTFILE_UN=*OUTFILE_UN;
	open(OUTFILE_MATCH,">>$outputdir/barcodetrimmed_$readfileparsed") or die("can't open $outputdir/barcodeEqualized_$readfileparsed to write to"); 
	my $OUTFILE_MATCH=*OUTFILE_MATCH;
	print("use one output file: $outputdir/barcodetrimmed_$readfileparsed\n");
	while(<READS>) {
		$header=$_;
		$seq=<READS>;
		$qualheader=<READS>;
		$qual=<READS>;
		$barcodematched=0;


##########
		$FILEHANDLE=$OUTFILE_UN;

		#go through every barcode
		for($j=0;$j<$bctotal;$j++) { 
			#if barcode matches sequence, then  print to that file

			#print "check bc:@barcodes[$j]...\n";
			if(matchbarcode($seq,@barcodes[$j])==1) { 
				$matchcounter++;
				$barcodematched=1;
				@matchedreadsperbarcode[$j]++;
				$seq=~s/^/'N' x @barcodesubtractedlengths[$j]/e; #add #:s to beginning to match with equalized barcode if barcode isn't max length				
				#$seq=~s/^@barcoderegexps[$j]/@barcodesEqualized[$j]/; #add N:s if this barcode isn't max length
				$qual=~s/^/'#' x @barcodesubtractedlengths[$j]/e; #add #:s to beginning to match with equalized barcode if barcode isn't max length
				$FILEHANDLE=$OUTFILE_MATCH;

				last; #match already found; stop searching matches for this read
			}
		}#endfor (barcodes list search)
		if($barcodematched==0) { $mismatchcounter++; }
		print $FILEHANDLE "$header$seq$qualheader$qual";
	}#endwhile		

###########



	close(READS);
	close(OUTFILE_UN);
	close(OUTFILE_MATCH);
}#END IF SEPARATEREADS==0 && inputType=="single"


#IF SEPARATEREADS==0 && inputType=="paired"
#reads _1 reads and corresponding _2 reads, then searches the barcodes from first file. If barcode is found from first file, then
#trims that barcode from both of the reads and outputs to correct output files.
if($actionType==0 & $inputType eq "paired") {
	my $FILEHANDLE1; #_1 file output
	my $FILEHANDLE2; #_2 file output
	open(OUTFILE_UN1,">>$outputdir/unmatched_reads1.fq") || die("can't open outfile unmatched_reads1.fq");
	open(OUTFILE_UN2,">>$outputdir/unmatched_reads2.fq") || die("can't open outfile unmatched_reads2.fq");
	my $OUTFILE_UN1=*OUTFILE_UN1;
	my $OUTFILE_UN2=*OUTFILE_UN2;
	open(OUTFILE_MATCH1,">>$outputdir/barcodetrimmed_$readfileparsed") or die("can't open $outputdir/barcodetrimmed_$readfileparsed to write to"); 
	open(OUTFILE_MATCH2,">>$outputdir/barcodetrimmed_$pairfileparsed") or die("can't open $outputdir/barcodetrimmed_$pairfileparsed to write to"); 

	my $OUTFILE_MATCH1=*OUTFILE_MATCH1;
	my $OUTFILE_MATCH2=*OUTFILE_MATCH2;

	print("use one output file for each of the read pair file: $outputdir/barcodetrimmed_$readfileparsed\n");
	while(<READS>) {
		$header=$_;
		$seq=<READS>;
		$qualheader=<READS>;
		$qual=<READS>;
		$pairheader=<PAIRS>;
		$pairseq=<PAIRS>;
		$pairqualheader=<PAIRS>;
		$pairqual=<PAIRS>;		
		
		$barcodematched=0;
		$FILEHANDLE1=$OUTFILE_UN1;	
		$FILEHANDLE2=$OUTFILE_UN2;	

		#go through every barcode
		for($j=0;$j<$bctotal;$j++) { 
			#if barcode matches sequence, then  print to that file

			#print "check bc:@barcodes[$j]...\n";
			if(matchbarcode($seq,@barcodes[$j])==1) { 
				$matchcounter++;
				$barcodematched=1;
				@matchedreadsperbarcode[$j]++;
				$seq=~s/^/'N' x @barcodesubtractedlengths[$j]/e; #add #:s to beginning to match with equalized barcode if barcode isn't max length				
				#$seq=~s/^@barcoderegexps[$j]/@barcodesEqualized[$j]/; #add N:s if this barcode isn't max length
				$qual=~s/^/'#' x @barcodesubtractedlengths[$j]/e; #add #:s to beginning to match with equalized barcode if barcode isn't max length
				$FILEHANDLE1=$OUTFILE_MATCH1; 
				$FILEHANDLE2=$OUTFILE_MATCH2; 

				last; #match already found; stop searching matches for this read
			}
		}#endfor (barcodes list search)
		if($barcodematched==0) { $mismatchcounter++; }
		print $FILEHANDLE1 "$header$seq$qualheader$qual";
		print $FILEHANDLE2 "$pairheader$pairseq$pairqualheader$pairqual";
	}#endwhile	

	close(READS);
	close(PAIRS);
	close(OUTFILE_UN1);
	close(OUTFILE_UN2);
	close(OUTFILE_MATCH1);
	close(OUTFILE_MATCH2);
}#END IF SEPARATEREADS==0 && inputType=="single"




############################
#IF SEPARATEREADS=2|3
############################
#IN CASE OF ONE OUTPUT FILE & SIZE EQUALIZATION

if($actionType==2 || $actionType==3) {
my $emptyRead='N'x$maxBarcodeLength; 
	my $FILEHANDLE;
	open(OUTFILE_MATCH,">>$outputdir/barcodetrimmed_$readfileparsed") or die("can't open $outputdir/barcodeEqualized_$readfileparsed to write to"); 
	my $OUTFILE_MATCH=*OUTFILE_MATCH;
	$FILEHANDLE=$OUTFILE_MATCH;	
	print("use one output file and mask unmatched barcodes with N:s; $outputdir/barcodetrimmed_$readfileparsed\n");
	while(<READS>) {
		$header=$_;
		$seq=<READS>;
		$qualheader=<READS>;
		$qual=<READS>;
		#go through every barcode
chomp $seq;
chomp $qual;
chomp $header;
		$barcodematched=0;



		#go through every barcode
		for($j=0;$j<$bctotal;$j++) { 
			#if barcode matches sequence, then  print to that file

			#print "check bc:@barcodes[$j]...\n";
			if(matchbarcode($seq,@barcodes[$j])==1) { 
				$matchcounter++;
				$barcodematched=1;
				@matchedreadsperbarcode[$j]++;
				if($actionType==2) {
					$seq=~s/^/'N' x @barcodesubtractedlengths[$j]/e; #add #:s to beginning to match with equalized barcode if barcode isn't max length				
					$qual=~s/^/'#' x @barcodesubtractedlengths[$j]/e; #add #:s to beginning to match with equalized barcode if barcode isn't max length
				}					
				if($actionType==3) {
					$header=substr($header,1);

					my $h=@barcodes[$j];
					my $tmpid=@ids[$j];
					$header="\@$tmpid"."_$h"."_$header";
					$seq=substr($seq,length($h));
					$qual=substr($qual,length($h));
					
				}
				last; #match already found; stop searching matches for this read
			}
		}#endfor (barcodes list search)
		if($barcodematched==0) { 
			$mismatchcounter++; 
			if($actionType==2) { $seq=~s/^[A-Z]{$maxBarcodeLength}/'N' x $maxBarcodeLength/e; }
			if($actionType==3) {
				$header=substr($header,1);
				my $h='N' x $maxBarcodeLength;
				$header="\@unknown_$h"."_$header";
			}
		}
		print $FILEHANDLE "$header\n$seq\n$qualheader$qual\n";

	}#endwhile	


	close(READS);
	close(OUTFILE_MATCH);
}#END IF SEPARATEREADS==2
############################
#IF SEPARATEREADS=3


my $allReads=$matchcounter+$mismatchcounter;
print("Matched $allReads reads/pairs with barcodes, of which $matchcounter barcodes and $mismatchcounter mismatches in reads/pairs.\n");
print("Summary:\nBarcode\tMatched reads/pairs\n"); 
my $bccount=@barcodes;
for(my $i=0;$i<$bccount;$i++) {
	print "@barcodes[$i]\t@matchedreadsperbarcode[$i]\n";

}

################################################################################################################################


#This sub gets a fastq file sequence as param1 and barcode as param2.
#Input sequence can contain hharacters A-Z. Char N is a special character
#which describes an ambiguous base (any base). N in input sequence barcode increases number of ambiguousbases count
#($ambiguousbases variable) in that input sequence, but doesn't affect to seqerrors count ($seqerrors variable)
#in the sequence. In addition to barcode, this script also checks if restriction site is found in the sequence,
#if -r parameter is used to define restriction site sequence (variable $restrictionsite).
#If number of ambiguous or seqerror bases is greater than $allowedAmbiguousbases or $allowedSeqerrors, then
#returns 0. Otherwise returns 1.
sub matchbarcode {
	my $inseq=shift;
	my $bc=shift; 
	my $bpindex_in_barcode=0;
	my $ambiguousbases=0; 	#N:s
	my $seqerrors=0; 		#wrong bases
	my $bclen=length($bc);
	for($bpindex_in_barcode=0;$bpindex_in_barcode<$bclen;$bpindex_in_barcode++) {
		my $seqbase=substr($inseq,$bpindex_in_barcode,1);
		my $bcbase=substr($bc,$bpindex_in_barcode,1);
		#print "ind $bpindex_in_barcode: seqbase $seqbase == bcbase $bcbase ???\n";
		if($seqbase eq 'N') { $ambiguousbases++; }
		else {
			if($seqbase ne $bcbase) { $seqerrors++; }
		}
		if($ambiguousbases>$allowedambiguousbases) { return 0; }
		if($seqerrors>$allowedseqerrors) { return 0; }		
	}#endfor
	
	#if restriction site check is in use (option -r / $restrictionsite defined
	#then compare these bases also
	if($restrictionsite ne "") {
		for(my $i=0;$i<$restrictionsitelength;$i++) {
			my $seqbase=substr($inseq,$bclen+$i,1);
			my $resbase=substr($restrictionsite,$i,1);
			if($seqbase eq 'N') { $ambiguousbases++; }
			else {
				if($seqbase ne $resbase) { $seqerrors++; }
			}
		}#endfor
			if($ambiguousbases>$allowedambiguousbases) { return 0; }
			if($seqerrors>$allowedseqerrors) { return 0; }
		
	}#endif
	return 1;
}#end sub
