use strict;
use warnings;

my($in,$out) = @ARGV;

open(IN,"<$in") or die $!;
open(OUT,">$out");
open(OUT1,">$out.stat");

my @samples;
my %hash_stat;

while(<IN>){
	chomp;
	my($chr,$start,$end,$mknum,@datas) = split/\t/;
	if($_ =~ /^#/){
		@samples = @datas;
		next;
	}
	next if($mknum < 1);
	my %hash_sample;
	for(my $i = 0; $i < @samples; $i++){
		my $sample = $samples[$i];
		my($idn,$cov) = split/\|/, $datas[$i];
		next unless($cov/$mknum > 1/3);
		$hash_sample{$sample}{idn} = $idn/$cov;
		$hash_sample{$sample}{cov} = $cov/$mknum;
	}
	my @samples_sorted = sort {$hash_sample{$b}{idn} <=> $hash_sample{$a}{idn} or $hash_sample{$b}{idn} <=> $hash_sample{$a}{idn}} keys %hash_sample;
	my @topSamples = ();
	my @topCovers = ();
	my $top;
	for(my $i = 0; $i < @samples_sorted; $i++){
		my $sample = $samples_sorted[$i];
		if($i == 0){
			push @topSamples, $sample;
			push @topCovers, $hash_sample{$sample}{cov};
			$top = $hash_sample{$sample}{idn};
			next;
		}
		last if($hash_sample{$sample}{idn} < $top);
		push @topSamples, $sample;
		push @topCovers, $hash_sample{$sample}{cov};
	}
	next if(@topSamples == 0);
	my $score = 1/@topSamples;
	foreach my $sample(@topSamples){
		$hash_stat{$sample}{abscore}++;
		$hash_stat{$sample}{rescore} += $score;
	}
	print OUT "$chr\t$start\t$end\t$mknum\t$top\t".join(",",@topSamples)."\t".join(",",@topCovers)."\n";
}
close IN;
close OUT;

foreach my $sample(@samples){
	my $abscore = 0;
	my $rescore = 0;
	if(exists $hash_stat{$sample}){
		$abscore = $hash_stat{$sample}{abscore};
		$rescore = $hash_stat{$sample}{rescore};
	}
	$rescore = sprintf("%.3f",$rescore);
	print OUT1 "$sample\t$abscore\t$rescore\n";
}
close OUT1;
