# Martin Taylor June 2015
# Given two bed files, strip all lines from file B where the name matches a line from file A but print the line from A. (Remove a subset of redundently mapped reads).

use strict;
my ($fileA,$fileB) = @ARGV;
my %filter;
open(FA, "<$fileA") or die "Failed to read $fileA\n";
while(<FA>){
    chomp;
    my @sp = split /\s+/;
    $filter{$sp[3]} = 1;
    print;
    print "\n";
}
close(FA);
open(FB, "<$fileB") or die "Failed to read $fileB\n";
while(<FB>){
    chomp;
    my @sp = split /\s+/;
    next if($filter{$sp[3]});
    print;
    print "\n";
}

