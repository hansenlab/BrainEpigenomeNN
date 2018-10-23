#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Cwd;

my $dir = getcwd();
my @data_R1 = ();
my @data_R2 = ();
my $hiseq = $ARGV[0];

my @mbias = glob "$dir/*/*.M-bias.txt";

for(my $i=0;$i<scalar(@mbias);$i++){
    
    my @split = split("/",$mbias[$i]);
    my @filename = split(/\./,$split[scalar(@split)-1]);

    open (FILE,$mbias[$i]) or die "Could not open: $!";
    open (IN,">> ".$hiseq."-mbias.txt") or die "Could not write to: $!";
    print IN $filename[0]."\n";

    while (<FILE>) {
        push @data_R1, $_ if (/^CpG context \Q(R1)/ .. /^125/);
	push @data_R2, $_ if (/^CpG context \Q(R2)/ .. /^125/); 
    }
    for(my $j=0;$j<scalar(@data_R1);$j++){
	print IN $data_R1[$j];
    }
    print IN "\n";
    for(my $j=0;$j<scalar(@data_R2);$j++){
        print IN $data_R2[$j];
    }
    print IN "\n";
    @data_R1 = ();
    @data_R2 = ();
}

close(IN);
