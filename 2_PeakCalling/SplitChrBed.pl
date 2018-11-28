#SplitChrBed.pl
#splits a bed file into individual chromosome files, only for chromosomes selected below in @chrlist
#written by Nicolas Altemose
#2015

use warnings;
use strict;

my $usage = "USAGE: perl SplitChrBed.pl <input file> <output file prefix>";


my $infile = "";
if(defined $ARGV[0]){
	$infile = $ARGV[0];
	chomp($infile);
}
else{
	die "$usage\n";
}
my $outfilebase = "";
if(defined $ARGV[1]){
	$outfilebase = $ARGV[1];
	chomp($outfilebase);
}
else{
	die "$usage\n";
}

my @chrlist = qw(
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
);

my %chrs;
foreach my $chr(@chrlist){
	$chrs{$chr}=0;
}

my $prevchr='';
my $ct=0;
open(IN,$infile) or die "could not open $infile\n";
while(my $line = <IN>){
	chomp($line);
	$line=~/^(chr\S+)\s+/;
	my $chr=$1;
	next unless exists $chrs{$chr};
	if($chr ne $prevchr){
		if($ct>0){
			close OUT;
		}
		open(OUT,'>'.$outfilebase.".$chr.bed");
	}
	print OUT "$line\n";
	$ct++;
	$prevchr=$chr;
}
close IN;
close OUT;
