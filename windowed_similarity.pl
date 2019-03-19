use strict;
use warnings;
use Getopt::Long;

my($vcf,$qName,$windowSize,$keepList,$pass,$pass1,$outName,$regList);

my $usage = "USAGE:\nperl $0 --in <vcf> --q <query name> --out <out prefix>\n";

GetOptions(
	"in=s" => \$vcf,
	"q=s" => \$qName,
	"out=s" => \$outName,
	"keep=s" => \$keepList,
	"bed=s" => \$regList,
	"window=s" => \$windowSize,
	"pass!" => \$pass,
	"pass1!" => \$pass1,
) or die $usage;

die $usage unless(defined $vcf and defined $qName and defined $outName);

unless(defined $windowSize){
	$windowSize = 10 * 1000;
}

# input keep samples
my %hash_sample;
if(defined $keepList){
	print "# Reading sample list...\n";
	open(IN,"<$keepList") or die $!;
	while(<IN>){
		chomp;
		$hash_sample{$_}{rank} = "";
	}
	close IN;
}

# input keep regions
my %hash_bed;
my %hash_chr;
if(defined $regList){
	open(IN,"<$regList") or die $!;
	while(<IN>){
		chomp;
		my($chr,$start,$end,$other) = split/\t/;
		next if(exists $hash_bed{$chr}{$start} and $hash_bed{$chr}{$start} > $end);
		$hash_bed{$chr}{$start} = $end;
	}
	close IN;
	# merge regions
	foreach my $chr(keys %hash_bed){
		$hash_chr{$chr} = "";
		my @starts = sort {$a <=> $b} keys %{$hash_bed{$chr}};
		for(my $i = 0; $i < @starts; $i++){
			my $starti = $starts[$i];
			my $endi = $hash_bed{$chr}{$starti};
			for(my $j = $i+1; $j < @starts; $j++){
				my $startj = $starts[$j];
				my $endj = $hash_bed{$chr}{$startj};
				last if($startj > $endi);
				if($endj > $endi){
					$endi = $endj;
					$hash_bed{$chr}{$starti} = $endj;
				}
				delete($hash_bed{$chr}{$startj});
				splice(@starts,$j,1);
				$j--;
			}
		}
	}
}

# read vcf file
open(IN,"<$vcf") or die $!;
open(OUT,">$outName.values");
# general settings
my @allSamples;
my @tRanks;
my @tSamples;
my @regions;
my $qRank;
# chr settings
my @remained_chrs;
# window settings
my $chr_tmp = "";
my $window_tmp = -1;
my $qCoverNum = 0;
my @tCoverNums;
my @identiNums;

while(<IN>){
	chomp;
	if($_ =~ /^#CHROM/){
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@allSamples) = split/\t/;
		for(my $i = 0; $i < @allSamples; $i++){
			next unless($allSamples[$i] eq $qName);
			$qRank = $i;
		}
		if(defined $keepList){
			for(my $i = 0; $i < @allSamples; $i++){
				next unless(exists $hash_sample{$allSamples[$i]});
				$hash_sample{$allSamples[$i]}{rank} = $i;
				push @tRanks, $i;
				push @tSamples, $allSamples[$i];
			}
		}else{
			for(my $i = 0; $i < @allSamples; $i++){
				$hash_sample{$allSamples[$i]}{rank} = $i;
				push @tRanks, $i;
				push @tSamples, $allSamples[$i];
			}
		}
		unless(defined $qRank){
			print "#ERROR: the query sample:$qName does not exist in the vcf file.\n";
			die;
		}
		unless(@tRanks > 0){
			print "#ERROR: no target sample found in the vcf file.\n";
			die;
		}
		print "\t# sample number=".@tRanks."\n";
		print OUT "#CHROM\tSTART\tEND\tMKNUM\t".join("\t",@tSamples)."\n";
		next;
	}
	next if($_ =~ /^#/);
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$datas_join) = split/\t/,$_,10;
	if(defined $pass1){
		next unless($filter eq "PASS");
	}
	if(defined $pass){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}
	# Check if the variant is a bi-allelic SNP
	next unless(length($ref) == 1 and length($alt) == 1);
	
	# Check if the position locates in the candidate regions
	goto NOBED unless(defined $regList);
	
	next unless(exists $hash_bed{$chr});
	if($chr ne $chr_tmp){
		@regions = ();
		if($qCoverNum > 0 and $window_tmp >= 0){
			my @outDatas;
			for(my $i = 0; $i < @tCoverNums; $i++){
				push @outDatas, "$identiNums[$i]|$tCoverNums[$i]";
			}
			my $windowStart = $windowSize * $window_tmp + 1;
			my $windowEnd = $windowSize * ($window_tmp + 1);
			print OUT "$chr_tmp\t$windowStart\t$windowEnd\t$qCoverNum\t".join("\t",@outDatas)."\n";
			$chr_tmp = $chr;
			$window_tmp = -1;
			$qCoverNum = 0;
			@tCoverNums = ();
			@identiNums = ();
		}
		foreach my $start(sort {$a <=> $b} keys %{$hash_bed{$chr}}){
			my $end = $hash_bed{$chr}{$start};	
			push @regions, "$start,$end";
		}
		$chr_tmp = $chr;
		if(exists $hash_chr{$chr}){
			delete($hash_chr{$chr});
		}
		@remained_chrs = keys %hash_chr;
		print "\t# Reading chr:$chr\n";
		if(@regions == 0){
			print "\t# Skipping chr:$chr\n";
		}
	}
	
	if(@remained_chrs == 0 and @regions == 0){
		last;
	}
	if(@regions == 0){
		next;
	}
	
	my $pickIt = 0;
	for(my $i = 0; $i < @regions; $i++){
		my($start,$end) = split/,/,$regions[$i];
		if($pos < $start){
			last;
		}elsif($pos > $end){
			splice(@regions, $i, 1);
			$i--;
			next;
		}else{
			$pickIt = 1;
			last;
		}
	}

	next unless($pickIt == 1);
	NOBED:
	
	# dissect
	my $windowRank = int($pos/$windowSize);
	if($chr_tmp ne $chr or $window_tmp ne $windowRank){
		if($window_tmp >= 0){
			my @outDatas;
			for(my $i = 0; $i < @tRanks; $i++){
				unless(defined $identiNums[$i]){
					$identiNums[$i] = 0;
				}
				unless(defined $tCoverNums[$i]){
					$tCoverNums[$i] = 0;
				}
				push @outDatas, "$identiNums[$i]|$tCoverNums[$i]";
			}
			my $windowStart = $windowSize * $window_tmp + 1;
			my $windowEnd = $windowSize * ($window_tmp + 1);
			print OUT "$chr_tmp\t$windowStart\t$windowEnd\t$qCoverNum\t".join("\t",@outDatas)."\n";
		}
		$chr_tmp = $chr;
		$window_tmp = $windowRank;
		$qCoverNum = 0;
		@tCoverNums = ();
		@identiNums = ();
	}
	
	my @datas = split/\t/, $datas_join;
	my $qGeno = "NA";
	if($datas[$qRank] =~ /^(\d+)\/(\d+)/ and $1 eq $2){
		$qGeno = $1;
		$qCoverNum++;
		for(my $i = 0; $i < @tRanks; $i++){
			my $tg = $datas[$tRanks[$i]];
			my $tGeno = "NA";
			if($tg =~ /^(\d+)\/(\d+)/ and $1 eq $2){
				$tGeno = $1;
				$tCoverNums[$i]++;
			}
			if($tGeno eq $qGeno){
				$identiNums[$i]++;
			}
		}
	}
}
close IN;
close OUT;
