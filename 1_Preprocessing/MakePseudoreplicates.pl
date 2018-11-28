#MakePseudoreplicates.pl
#splits a fragment position bedfile into two pseudoreplicates
#by Nicolas Altemose
#2015

use warnings;
use strict;

my $usage = "USAGE: perl MakePseudoreplicates.pl <input file> <Proportion in Replicate 1 [default 0.5]> <chr size file [default hg19.chromsizes.tbl]> <bedtools path [default ~/bin]>";

my $infile;
my $bedtoolspath = '~/bin';
my $proprep1=0.5;
my $chrsizefile = 'hg19.chromsizes.tbl';


if(defined $ARGV[0]){
	$infile = $ARGV[0];
	chomp($infile);
}
else{
	die "$usage\n";
}
if(defined $ARGV[1]){
	$proprep1 = $ARGV[1];
	chomp($proprep1);
}
if(defined $ARGV[2]){
	$chrsizefile = $ARGV[2];
	chomp($chrsizefile);
}
if(defined $ARGV[3]){
	$bedtoolspath = $ARGV[3];
	chomp($bedtoolspath);
}
unless(-e $chrsizefile){
	die "ERROR: Could not find $chrsizefile in this directory\n"; 
}


my $outfile1 = "$infile.PR1.sorted.bed";
my $outfile2 = "$infile.PR2.sorted.bed";

my $outfile3 = "FragDepth.$infile.PR1.bed";
my $outfile4 = "FragDepth.$infile.PR2.bed";



my $ct1=0;
my $ct2=0;
open(IN,$infile);
open(OUT1,'>'.$outfile1);
open(OUT2,'>'.$outfile2);
while(my $line = <IN>){
	chomp($line);
	my $rando = rand;
	if($rando<$proprep1){
		print OUT1 "$line\n";
		$ct1++;
	}
	else{
		print OUT2 "$line\n";
		$ct2++;
	}
}
close IN;
close OUT1;
close OUT2;

my $tot = $ct2+$ct1;
my $prop1 = $ct1/$tot;
my $prop2 = $ct2/$tot;
print "\nprinted $ct1 ($prop1) to $outfile1\n";
print "\nprinted $ct2 ($prop2) to $outfile2\n";

print "getting coverage for FragDepth files\n";

system("$bedtoolspath/bedtools genomecov -bg -i $outfile1 -g $chrsizefile >$outfile3");
system("$bedtoolspath/bedtools genomecov -bg -i $outfile2 -g $chrsizefile >$outfile4");


