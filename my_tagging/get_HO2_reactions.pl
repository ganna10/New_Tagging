#! /usr/bin/env perl
# Get HO2 reactions from my mozart chemistry
# Version 0: Jane Coates 27/5/2015
#
### next version also output reaction number
### rate of HO2_X + HO2 = HO2 has to be multiplied by 2

use strict;
use diagnostics;
use KPP;

my $species = "HO2";
my $eqn = "my_mozart.eqn";
my $kpp = KPP->new($eqn);
my (%producers, %consumers);
my $producers = $kpp->producing($species);
my $consumers = $kpp->consuming($species);

foreach my $reaction (@$producers) {
    my $reaction_string = $kpp->reaction_string($reaction);
    next unless ($reaction_string =~ /HO2NO2/);
    my ($reactants, $products) = split / = /, $reaction_string;
    $reactants = "NO2HO2_X";
    $products = "HO2_X";
    my $rate_string = $kpp->rate_string($reaction);
    my $string = $reactants . " = " . $products . " : " . $rate_string . " ;";
    $producers{$string} += 1;
}

foreach my $reaction (@$consumers) {
    my $reaction_string = $kpp->reaction_string($reaction);
    my $rate_string = $kpp->rate_string($reaction);
    my ($reactants, $products) = split / = /, $reaction_string;
    $reactants =~ s/HO2/HO2_X/;
    if ($reactants eq "HO2_X \+ NO") {
        $products = "NO2_X + NO";
    } elsif ($reactants eq "HO2_X + NO2") {
        $products = "NO2HO2_X + NO2";
    } else {
        ($products = $reactants) =~ s/HO2_X//;
        $products =~ s/ \+ //;
    }
    my $string = $reactants . " = " . $products . " : " . $rate_string . " ;";
    $consumers{$string} += 1;
}

my $out_file = "ho2_tags.txt";
open my $out, '>:encoding(utf-8)', $out_file or die $!;
print $out "#HO2 producers\n";
print $out "$_\n" foreach (sort keys %producers);
print $out "#HO2 consumers\n";
print $out "$_\n" foreach (sort keys %consumers);
print $out "#NO2HO2 consumers\n";
print $out "NO2HO2_X + OH = OH : 3.2D-13*EXP(690/TEMP)*1.0 ;\n";
print $out "NO2HO2_X = UNITY : MOZART_VD(KPP_HO2NO2)/(zmbl*100.) ;\n";
close $out;
