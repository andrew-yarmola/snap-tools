#!/usr/bin/perl
use strict;
use warnings;
# use Data::Dumper;

sub get_geods {
    my $n = shift;
    my $l = shift;
    return  "echo 'read closed $n
                   print geodesics $l' | snap |";
}

sub get_ortho {
    my $n = shift;
    my $l = shift;
    my $d = shift;
    return "echo 'read closed $n
                  print geodesics $l
                  print ortholines $d @_' | snap |";
}

MAIN: 
for my $n (1 .. 11031) {
    print STDERR "$n\n";
    if ($n == 6109) { next MAIN; }
    my $min_len = -1;
    my %geods;
    my $l = 0.8;
    while (keys %geods < 8) {
        $l += 0.2;
        open(PS, get_geods($n,$l)) or die "Failed on : $!";
        while (<PS>) {
            if (/\[(\d+)\]  ([0-9\.]+)([0-9\.\+\-]+)*/) {
                my $m = $1;
                my $real = $2;
                my $imag = $3;
                if ($real == $min_len || keys %geods < 8) {
                    if ($min_len < 0) { $min_len = $real; }
                    $geods{$m} = "$real$imag*1j";
                } 
            }
            if (/.*Dirichlet.*/) {
                print STDERR "$n Dirichlet domain failed\n";
                next MAIN;
            }
        }
        close PS;
    }
#    print Dumper(\%geods);

    my $num_geods = keys %geods;
    my $num_pairs = int($num_geods*($num_geods - 1)/2);
    my %orthos;
    my $d = 1.8;
    while (keys %orthos < $num_pairs) {
        $d += 0.2;
        open(PS, get_ortho($n, $l, $d, keys %geods)) or die "Failed on : $!";
        while (<PS>) {
            if (/([0-9\.]+)([0-9\.\+\-]+)\*i  (\d+):.+  (\d+):.*/) {
                my $real = $1;
                my $imag = $2;
                my $left = $3;
                my $right = $4;
                if (! exists $orthos{"$left:$right"} && $real != 0.0) {
                    $orthos{"$left:$right"} = "$real$imag*1j"; 
                } 
            }
        }
        close PS;
    }
#    print Dumper(\%orthos);

    while (my ($k, $v) = each %orthos) {
        my ($left, $right) = split(/:/, $k);
        print "$n, $geods{$left}, $geods{$right}, $v\n";
    }
}
