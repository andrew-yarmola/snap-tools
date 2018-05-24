#!/usr/bin/perl
use strict;
use warnings;
# use Data::Dumper;

my $type = 'closed';
my $mfld_count  = 11031;
# my $mfld_count  = 1;
my %cusped_counts = ( 5 => 414, 6 => 961, 7 => 3551 );
# my %cusped_counts = ( 5 => 1, 6 => 1, 7 => 1 );
if (@ARGV > 0) {
    if ($ARGV[0] eq 'cusped') {
        my $tet_num = 7;
        if (@ARGV > 1) {
            $tet_num = int($ARGV[1]);
            if ($tet_num < 5 || 7 < $tet_num) { $tet_num = 7; }
        }
        $type = "census $tet_num";
        $mfld_count = $cusped_counts{$tet_num};
    }
}

sub get_geods {
    my $t = shift;
    my $n = shift;
    my $l = shift;
    return  "echo 'read $t $n
                   print geodesics $l' | snap |";
}

sub get_ortho {
    my $t = shift;
    my $n = shift;
    my $l = shift;
    my $d = shift;
    return "echo 'read $t $n
                  print geodesics $l
                  print ortholines $d @_' | snap |";
}

MAIN: 
for my $n (1 .. $mfld_count) {
    print STDERR "$type $n\n";
    if ($type eq 'closed' && $n == 6109) { next MAIN; }
    my $max_len = 0;
    my %geods;
    my $l = 0.8;
    while ($max_len < 0.8) {
        $l += 0.2;
        open(PS, get_geods($type, $n,$l)) or die "Failed on : $!";
        while (<PS>) {
            if (/\[(\d+)\]  ([0-9\.]+)([0-9\.\+\-]+)*/) {
                my $m = $1;
                my $real = $2;
                my $imag = $3;
                if ($max_len < $real) { $max_len = $real; }
                $geods{$m} = "$real$imag*1j";
            } elsif (/.*Dirichlet.*/) {
                print STDERR "$n Dirichlet domain failed\n";
                next MAIN;
            }
        }
        close PS;
    }
#    print Dumper(\%geods);

    my $num_geods = keys %geods;
    my $num_pairs = int($num_geods*($num_geods - 1)/2) + $num_geods;
    my %orthos;
    my $d = 1.8;
    while (keys %orthos < $num_pairs) {
        $d += 0.2;
        open(PS, get_ortho($type, $n, $l, $d, keys %geods)) or die "Failed on : $!";
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
        print "$type $n, $geods{$left}, $geods{$right}, $v\n";
    }
}
