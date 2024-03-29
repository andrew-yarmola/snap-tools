#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use Math::Complex;
use List::Util qw[min max];

my $length_bound = 1.5;
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
    my ($t, $n, $l) = @_;
    return  "echo 'read $t $n
                   print name
                   print geodesics $l' | snap |";
}

sub get_ortho {
    my ($t, $n, $l, $d, @gs)  = @_;
    return "echo 'read $t $n
                  print geodesics $l
                  ortholines $d @gs' | snap |";
}

print "\"name\",\"geod\",\"length\",\"ortho\",\"move\"\n";

MAIN: 
for my $n (1 .. $mfld_count) {
    my $name = "$type $n";
    print STDERR "$name\n";
    my %geods;
    # Simple, but hang suseptable way of proecessing system command
    open(PS, get_geods($type, $n, $length_bound)) or die "Failed on : $!";
    while (<PS>) {
        if (/\[(\d+)\]  ([0-9\.]+)([0-9\.\+\-]+)*/) {
            my $m = $1;
            my $real = $2;
            my $imag = $3;
            $geods{$m} = $real + $imag * i;
        } elsif (/.*name : (.*)/) {
            $name = $1;
        } elsif (/.*Dirichlet.*/) {
            print STDERR "Dirichlet domain failed for $type $n\n";
            next MAIN;
        }
    }
    close PS;

    if (keys %geods < 1) { next MAIN; }

    # Different/safer way of processing system command
    my %orthos;
    INNER:
    foreach my $g (keys %geods) {
        my $ortho_len = max(Re($geods{$g}),2); 
        while ((! exists $orthos{$g}) and ($ortho_len < $ortho_bound)) { 
            my $pid;
            eval {
                local $SIG{ALRM} = sub { die "alarm\n" }; # NB: \n required
                alarm 240;
                $pid = open(PS, get_ortho($type, $n, $length_bound, $ortho_len, $g)) or die "Failed on : $!";
                while (<PS>) {
                    if (/([0-9\.]+)([0-9\.\+\-]+)\*i  (\d+):.+  (\d+):.*/) {
                        my $real = $1;
                        my $imag = $2;
                        my $left = $3;
                        my $right = $4; 
                        if (($left == $right) and ((! exists $orthos{$left}) or ($real < Re($orthos{$left})))) {
                            if ($real != 0.0) {
                                $orthos{$left} = $real + $imag * i; 
                            } else {
                                print STDERR "Found 0 length ortholine for $type $n\n";
                            }
                        } 
                    }
                }
                alarm 0;
                close PS;
            };
            alarm 0; # race contition 
            if ($@) {
                if ($@ eq "alarm\n") {
                    print STDERR "Timeout for $type $n\n";
                } else {
                    print STDERR "Failed for $type $n\n";
                }            
                if (defined $pid) { kill 9, $pid+2; } # HACK!!! echo and snap in the command get their own pids.
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
        my $out = "\"$name\",$g,$l,$r,$d\n";
        $out =~ s/i/j/g;
        print "$out";
    }
}
