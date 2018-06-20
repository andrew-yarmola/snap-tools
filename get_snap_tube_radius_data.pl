#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use Math::Complex;
use List::Util qw[min max];

my $length_bound = 2.0;
my $ortho_bound = 8.0;
my $type = 'closed';
my $mfld_count  = 11031;
my %cusped_counts = ( 5 => 414, 6 => 961, 7 => 3551 );
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
    my %geods;
    # Simple, but hang suseptable way of proecessing system command
    open(PS, get_geods($type, $n, $length_bound)) or die "Failed on : $!";
    while (<PS>) {
        if (/\[(\d+)\]  ([0-9\.]+)([0-9\.\+\-]+)*/) {
            my $m = $1;
            my $real = $2;
            my $imag = $3;
            $geods{$m} = $real + $imag * i;
        } elsif (/.*Dirichlet.*/) {
            print STDERR "Dirichlet domain failed for $type $n\n";
            next MAIN;
        }
    }
    close PS;

    if (keys %geods < 1) { next MAIN; }

    # Different/safer way of processing system command
    my $pid;
    my %orthos;
    INNER:
    foreach my $g (keys %geods) {
        my $ortho_len = max(Re($geods{$g}),2); 
        while ((! exists $orthos{$g}) and ($ortho_len < $ortho_bound)) { 
            eval {
                local $SIG{ALRM} = sub { die "alarm\n" }; # NB: \n required
                alarm 180;
                $pid = open(PS, get_ortho($type, $n, $length_bound, $ortho_len, $g)) or die "Failed on : $!";
                while (<PS>) {
                    if (/([0-9\.]+)([0-9\.\+\-]+)\*i  (\d+):.+  (\d+):.*/) {
                        my $real = $1;
                        my $imag = $2;
                        my $left = $3;
                        my $right = $4; 
                        if (($left == $right) and (! exists $orthos{$left})) {
                            if ($real != 0.0) {
                                $orthos{$left} = $real + $imag * i; 
                            } else {
                                print STDERR "Found 0 length ortholine for $type $n\n";
                            }
                        } 
                    }
                }
                close PS;
                alarm 0;
            };
            if ($@) {
                print STDERR "Timeout for $type $n\n";
                kill 9, $pid+2; # HACK!!! echo and snap in the command get their own pids.
                next INNER;
            }
            $ortho_len += 0.5;
        }
    }

    while (my ($g, $r) = each %orthos) {
        my $l = $geods{$g};
        my $d = 100.0;
        my $k = 1;
        while (Re($k*$l) < $d) {
            $d = min($d, acosh(cosh(Re($k*$l))*cosh(Re($r)/2)**2 - cos(Im($k*$l))*sinh(Re($r)/2)**2));
            $k += 1;
        }
        my $out = "$type $n, geod $g, length $l, ortho $r, move $d\n";
        $out =~ s/i/*1j/g;
        print "$out";
    }
}
