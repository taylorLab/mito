#!/usr/bin/perl -w
#
# produce frequency distributions
#
use strict;
my (@sp,$lim,%index,@range,$total);
my $numeric = 1;
my ($column,$pc) = @ARGV;
if (!$pc){
	$pc = 0;
}
$total = 0;
while (<STDIN>){
	chomp;
	@sp = split /\s+/;
	if (!defined $index{$sp[$column]}){
		$index{$sp[$column]} = 1;
		if ($sp[$column] =~ /\D/){
			$numeric = 0;
		}
		if  ($sp[$column] =~ /^00/){
		    $numeric = 0;
		}
	}
	else {
		$index{$sp[$column]}++;
	}
	$total++;
}
if ($pc == 1){
	foreach $lim (keys %index){
		$index{$lim} = ($index{$lim} / $total) * 100;
	}
}
elsif ($pc eq "s"){
	$numeric = 0;
}
if ($numeric == 1){
	@range = sort { $a <=> $b } keys %index;
	for $lim ($range[0] .. $range[$#range]){
		print "$lim\t";
		if (defined $index{$lim}){
			print "$index{$lim}\n";
		}
		else {
			print "0\n";
		}
	}
}
else {
	@range = sort { $a cmp $b } keys %index;
	for $lim (0 .. $#range){
		print "$range[$lim]\t";
		if (defined $index{$range[$lim]}){
			print "$index{$range[$lim]}\n";
		}
		else {
			print "0\n";
		}
	}
}
	
