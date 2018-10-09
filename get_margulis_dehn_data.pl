#!/usr/bin/env perl
use strict;
use warnings;
use Math::Trig;
use Math::Complex;
use List::Util qw[min max pairs];

my $type = '5';
my @mld_idx = (3, 4, 6, 7, 9, 10, 11, 15);

sub gcd {
    my ($a, $b) = (abs($_[0]), abs($_[1]));
    ($a,$b) = ($b,$a) if $a > $b; 
    if ($b == 0) { return $a; }
    while ($a) {
        ($a, $b) = ($b % $a, $a);
    } 
    return $b; 
}

my $s_len = 20;
my @slopes = (); 
for my $i (-$s_len .. $s_len) {
    for my $j (0 .. $s_len) {
        if (($j > 0 or $i > 0) and $i**2 + $j**2 <= $s_len**2 and gcd($i,$j) == 1) {
            my @p = ($i, $j);
            push @slopes, \@p;
        }   
    }   
}

#"\"name\",\"volume\",\"margulis\",\"left\",\"right\",\"ortho\"\n";

MAIN: 
for my $n (@mld_idx) {
    SLOPE:
    for my $s (@slopes) {
        my $name = "$type $n fill (@$s)";
        print STDERR "$name\n";
        `./margulis.py $type $n 0.8 @$s >> data_margulis_dehn_20.csv` 
    }
}
