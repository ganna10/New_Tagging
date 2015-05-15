#! /usr/bin/env perl
# Check which species have already been tagged and what tags they have, argument is .spc file
# Version 0: Jane Coates 15/5/2015

use strict;
use diagnostics;

my $file = $ARGV[0];
die "Specify .spc file" unless (defined $file);
open my $in, '<:encoding(utf-8)', $file or die $!;
my @lines = <$in>;
close $in;

my (%species, %tagged_species);
foreach my $line (@lines) {
    next unless ($line =~ /IGNORE/);
    my ($spc, $rest) = split / = /, $line;
    $species{$spc} = 1;
}

%tagged_species = (
    "INI"   => [],
    "CH4"   => [],
    "XTR"   => [],
    "notag" => [],
);

foreach my $spc (sort keys %species) {
    my ($species, $tag) = split /_/, $spc;
    if (defined $tag) {
        push @{$tagged_species{$tag}}, $species;
    } else {
        push @{$tagged_species{"notag"}}, $species;
    }
}

my($day, $month, $year) = (localtime)[3,4,5];
$month = sprintf '%02d', $month+1;
$day   = sprintf '%02d', $day;
$year = 1900 + $year;
my $out_file = $file . "tagged_species" . $year . $month . $day;
open my $out, '>:encoding(utf-8)', $out_file or die $!;
foreach my $tag (sort keys %tagged_species) {
    print $out "$tag : @{$tagged_species{$tag}}\n";
}
close $out;
