#GetFragmentDepth.pl
#Given a bam file with paired reads (piped in via samtools), outputs a bed file with fragment positions (after removing read pairs with both reads below a certain mapping quality threshold)
#by Nicolas Altemose
#2015

use warnings;
use strict;
my $tic = time;
print "\n\n";

my $usage = "USAGE: samtools view -F12 -q1 <SortedBAMfile.bam> | perl GetFragmentPositions.pl <Fragment Position Output File> <read length [default 51]> <max frag length [default 10000]>";

my $outfile1;
my $outfile2;
my $readlen = 51;
my $maxlen = 10000; #maximum allowable fragment length


if(defined $ARGV[0]){
	$outfile2 = $ARGV[0];
	chomp($outfile2);
}
else{
	die "$usage\n";
}
if(defined $ARGV[1]){
	$readlen = $ARGV[1];
	chomp($readlen);
}
if(defined $ARGV[2]){
	$maxlen = $ARGV[2];
	chomp($maxlen);
}


open(OUT2,'>'.$outfile2.'.bedtemp');
my @memorylist;
my %memorynames;

my $prevchr = '';
my $prevpos=0;
my $chrnum=0;
my $sort=0;

print "chr\n";
while(my $line = <STDIN>){
	chomp($line);
	if($line=~/^\S+\t\S+\t(\S+)\t(\S+)\t\S+\t\S+\t(\S+)\t(\S+)\t/){
		my $rchr = $1;
		next if($rchr eq '*');
		my $rpos = $2;
		my $pchr = $3;
		my $ppos = $4;

		if($rchr ne $prevchr){
			$chrnum++;
			@memorylist=();
			%memorynames=();
			$prevpos=0;
			print "reading $rchr\n";
		}

		if($rchr eq $pchr || $pchr eq '='){
			if(((abs($ppos-$rpos)+$readlen+1)<=$maxlen)){
				my $start=$rpos;
				if($ppos>$rpos){
					my $stop = $ppos+$readlen-1;
					my $start0 = $start-1;
					print OUT2 "$rchr\t$start0\t$stop\t1\t0\t+\n";
					push(@memorylist,$rpos);
					$memorynames{$rpos}=0;
				}
				else{
					unless(exists $memorynames{$ppos}){
						$start = $ppos;
						my $stop = $rpos+$readlen-1;
						my $start0 = $start-1;
						print OUT2 "$rchr\t$start0\t$stop\t1\t0\t+\n";
					}
				}
				$prevpos = $start;
			}
		}
		foreach my $pos(@memorylist){
			if($pos<($rpos-2*$maxlen)){
				shift(@memorylist);
				delete $memorynames{$pos};
			}
			else{
				last;
			}
		}
		$prevchr = $rchr;
	}
	elsif($line!~/^\@/){
		print "ERROR parsing line: $line\n";
	}
}
close OUT2;

print "\nsorting...\n";

my $sortedoutfile = "$outfile2.sorted.bed";
system("sort -k1,1V -k2,2n $outfile2.bedtemp >$sortedoutfile");
system("rm -f $outfile2.bedtemp");


#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
