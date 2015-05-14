#! /usr/bin/env perl
# Get tagged RO2 and radicals from Tim's tagged MOZART file
# Version 0: Jane Coates 23/4/2015

use strict;
use diagnostics;

my $RO2_file = "RO2_species.txt";
my $radicals_file = "radicals.txt";
my $spc_file = "mozart.spc";
my @loop = ($RO2_file, $radicals_file);

open my $spc_in, '<:encoding(utf-8)', $spc_file or die $!;
my @lines = <$spc_in>;
close $spc_in;
my %tagged_species;
foreach my $line (@lines) {
    next unless ($line =~ /IGNORE/ and $line =~ /_/);
    chomp $line;
    $line =~ s/ = IGNORE(.*?)$//;
    $tagged_species{$line} += 1;
}
my @tagged_species = keys %tagged_species;

foreach my $file (@loop) {
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    my @lines = <$in>;
    close $in;

    my @species;
    foreach my $line (@lines) {
        chomp $line;
        push @species, split /\s/, $line;
    }

    my @tagged_list;
    foreach my $species (@species) {
        foreach my $tagged (@tagged_species) {
            push @tagged_list, $tagged if ($tagged =~ /^${species}_/);
        }
    }
    push @species, @tagged_list;
    
    my $out_file = $file . ".tagged";
    open my $out, '>:encoding(utf-8)', $out_file or die $!;
    print $out "@species\n";
    close $out;
}
