my $offset = -1;
if (defined $ENV{OFFSET}){
    $offset = $ENV{OFFSET};
}
while (<STDIN>){
    chomp;
    my @F = split /\s+/;
    if($F[5]=~/\+/){
	$F[1] += $offset;
    }
    else{
	# The minus-one because we are converted an bed end to a bed start.
	$F[1]=($F[2] - $offset) - 1;
    }
    if($F[1] < 0){
	# because bedtools can't deal with negative numbers, but we can
	# tell it that the chromosome is longer than it really is and then
	# reconvert coordinates later.
	$F[1] = $ENV{LENCIRC} + $F[1];
    }
    # Apply the strand flip after offset as offests are relative to the
    # read 5' end, not the ribo-strand.
    if($ENV{FLIPSTRAND}){
	# if We are to flip the strand..
	$F[5]=~tr/\+-/-\+/;
    }
    $F[2] = $F[1] + 1;
    print join qq{\t}, @F;
    print "\n";
}
