#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use Math::Complex;

my $margulis_bound = 0.8;
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
                   print geodesics $l' | snap |";
}

sub get_ortho {
    my ($t, $n, $l, $d, @gs)  = @_;
    return "echo 'read $t $n
                  print geodesics $l
                  print ortholines $d @gs' | snap |";
}

MAIN: 
for my $n (1 .. $mfld_count) {
    print STDERR "$type $n\n";
    my %geods;
    # Simple, but hang suseptable way of proecessing system command
    open(PS, get_geods($type, $n, $margulis_bound)) or die "Failed on : $!";
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

    # Largest othroline length so that at midpoint the shortest geod can move by 0.8
    my $ortho_bound = 0.5;
    for my $w (values %geods) {
        my $z = $w; # copy to inciment
        while (Re($z) < $margulis_bound) {
            my $d = 2.*acosh(sqrt((cosh($margulis_bound)-cos(Im($z)))/(cosh(Re($z))-cos(Im($z)))));
            print STDERR "$z --> $d\n";
            if ($ortho_bound < $d) {
                $ortho_bound = $d;
            }
            $z += $w;
        }
    }
    $ortho_bound += 0.05; # just for sanity
    print STDERR "$ortho_bound\n";

    # Different/safer way of processing system command
    my $pid;
    my %orthos;
    eval {
        local $SIG{ALRM} = sub { die "alarm\n" }; # NB: \n required
        alarm 180;
        $pid = open(PS, get_ortho($type, $n, $margulis_bound, $ortho_bound, keys %geods)) or die "Failed on : $!";
        while (<PS>) {
            if (/([0-9\.]+)([0-9\.\+\-]+)\*i  (\d+):.+  (\d+):.*/) {
                my $real = $1;
                my $imag = $2;
                my $left = $3;
                my $right = $4;
                if (! exists $orthos{"$left:$right"}) {
                    if ($real != 0.0) {
                        $orthos{"$left:$right"} = $real + $imag * i; 
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
        next MAIN;
    }

    while (my ($k, $v) = each %orthos) {
        my ($left, $right) = split(/:/, $k);
#        print "    from [$left] to [$right] --> $v\n"; 
        my $out = "$type $n, $geods{$left}, $geods{$right}, $v\n";
        $out =~ s/i/*1j/g;
        print "$out";
    }
}
