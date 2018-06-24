#!/usr/bin/env perl
use strict;
use warnings;
use Math::Trig;
use Math::Complex;
use List::Util qw[min max pairs];

my $length_bound = 2.0;
my $ortho_bound = 8.0;
my $type = 'census 5';
my @mld_idx = (3, 4, 6, 7, 9, 10, 11, 15);
my $s_len = 10;
if (@ARGV > 0) {
   $s_len = $ARGV[0]; 
}

sub get_geods {
    my ($t, $n, $s, $l) = @_;
    return  "echo 'read $t $n
                   surgery @$s
                   print name
                   print geodesics $l' | snap |";
}

sub get_ortho {
    my ($t, $n, $s, $l, $d, @gs)  = @_;
    return "echo 'read $t $n
                  surgery @$s
                  print geodesics $l
                  ortholines $d @gs' | snap |";
}

sub gcd {
    my ($a, $b) = (abs($_[0]), abs($_[1]));
    ($a,$b) = ($b,$a) if $a > $b; 
    if ($b == 0) { return $a; }
    while ($a) {
        ($a, $b) = ($b % $a, $a);
    } 
    return $b; 
}

my @slopes = (); 
for my $i (-$s_len .. $s_len) {
    for my $j (0 .. $s_len) {
        if (($j > 0 or $i > 0) and $i**2 + $j**2 <= $s_len**2 and gcd($i,$j) == 1) {
            my @p = ($i, $j);
            push @slopes, \@p;
        }   
    }   
}

print "\"name\",\"geod\",\"length\",\"ortho\",\"move\"\n";

MAIN: 
for my $n (@mld_idx) {
    SLOPE:
    for my $s (@slopes) {
        my $name = "$type $n fill (@$s)";
        print STDERR "$name\n";
#        my @b = @$s;
#        if ($b[0] != 0 or $b[1] != 1) {
#            print STDERR "skip\n";
#            next SLOPE;
#        }
        my %geods;
        # Timeout way of processing system command
        my $g_pid;
        eval {
            local $SIG{ALRM} = sub { die "alarm\n" }; # NB: \n required
            alarm 180;
            $g_pid = open(PS, get_geods($type, $n, $s, $length_bound)) or die "Failed on : $!";
            while (<PS>) {
                if (/\[(\d+)\]  ([0-9\.]+)([0-9\.\+\-]+)*/) {
                    my $m = $1;
                    my $real = $2;
                    my $imag = $3;
                    $geods{$m} = $real + $imag * i;
                } elsif (/.*name : (.*)/) {
                    $name = $1;
                } elsif (/.*Dirichlet.*/) {
                    print STDERR "Dirichlet domain failed for $type $n fill (@$s)\n";
                    die "Dirichlet\n";
                } elsif (/.*(other|flat|degenerate).*/) {
                    print STDERR "Non-hyperbolic $1 solution for $type $n fill (@$s)\n";
                    die "nonhyp\n";
                }
            }
            alarm 0;
            close PS;
        };
        alarm 0; # race condition
        if ($@) { 
            if ($@ eq "alarm\n") {
                print STDERR "Timeout for $type $n fill (@$s) geod\n";
            } else {
                print STDERR "Failed for $type $n fill (@$s) geod\n";
            }            
            if (defined $g_pid) { kill 9, $g_pid+2; } # HACK!!! echo and snap in the command get their own pids.
            next SLOPE;
        }

        if (keys %geods < 1) { next SLOPE; }

        # Timeout way of processing system command
        my %orthos;
        INNER:
        foreach my $g (keys %geods) {
            my $ortho_len = max(Re($geods{$g}),2);
            while ((! exists $orthos{$g}) and ($ortho_len < $ortho_bound)) { 
                my $o_pid;
                eval {
                    local $SIG{ALRM} = sub { die "alarm\n" }; # NB: \n required
                    alarm 240;
                    $o_pid = open(PS, get_ortho($type, $n, $s, $length_bound, $ortho_len, $g)) or die "Failed on : $!";
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
                    alarm 0;
                    close PS;
                };
                alarm 0; # race contition 
                if ($@) { 
                    if ($@ eq "alarm\n") {
                        print STDERR "Timeout for $type $n (@$s) ortho for geod $g\n";
                    } else {
                        print STDERR "Failed for $type $n (@$s) ortho for geod $g\n";
                    }
                    if (defined $o_pid) { kill 9, $o_pid+2; }# HACK!!! echo and snap in the command get their own pids.
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
            $out =~ s/i/*1j/g;
            print "$out";
        }
    }
}
