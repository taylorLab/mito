#
# Martin Taylor, September 2015
# martin.taylor@igmm.ed.ac.uk
#
# Circularisation of DNA and bedtools doesn't play nice.
# Process bed files to clean for simple circularisatoin issues.
# Work on piped and assume that further processing on fasta file
# has a tandemly repeated circular chromosome sequence so
# grabbing over the end works correctly. If we need to grab
# before the start we can reset the coordinate to the concatamer
# junction.
# Need to tidy up the coordinate transforms at the end.

use strict;
# Set max size of circular chromosomes
my %circChrom = ('chrM' => $ENV{LENCIRC});
my ($slopSize) = @ARGV;

while (<STDIN>){
    chomp;
    my @sp = split /\s+/;
    $sp[1]-=$slopSize;
    $sp[2]+=$slopSize;
    if($sp[1] < 0){
	if($circChrom{$sp[0]}){
	    # the slop has pushed us round the corner of the circle!
	    $sp[1] = $circChrom{$sp[0]} + $sp[1];
	    $sp[2] = $circChrom{$sp[0]} + $sp[2];
	}
	else {
	    print STDERR "WARN: ".join"\t",@sp;
	    print STDERR "\n";
	    next;
	}
    }
    print join "\t", @sp;
    print "\n";
}
