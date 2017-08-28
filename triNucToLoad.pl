use strict;
while(<>){
    chomp;
    my @sp = split /\s+/;
    my @qq = split /:/, $sp[0];
    $sp[1] = uc $sp[1];
    if($qq[0] eq "chrM"){
	if($qq[1] > $ENV{LENCIRC}){
	    $qq[1] = $qq[1] - $ENV{LENCIRC}; 
	}
    }
    print "$qq[0]\t$qq[1]\t$qq[2]\t$qq[3]\t$sp[1]\n";
}
