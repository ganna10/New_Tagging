#! /usr/bin/env perl
# get Ox reactions for tagging
# Version 0: Jane Coates 29/5/2015
#
### next version also output reaction number

use strict;
use diagnostics;
use KPP;

my $eqn = "my_mozart.eqn";
my $kpp = KPP->new($eqn);
my (%weights, %reactions);
my @Ox = qw( O3 O O1D NO2 NO3 N2O5 HNO3 PAN MPAN ONIT ONITR ISOPNO3 HO2NO2 ) ;
my %families = ( 
    "Ox" => [ @Ox ],
);
$kpp->family({
        name    => "Ox",
        members => $families{"Ox"},
        weights => $weights{"Ox"},
});

#reactions with each Ox species
foreach my $ox (@Ox) {
    my $reactions = $kpp->consuming($ox);
    foreach my $reaction (@$reactions) {
        my $reaction_string = $kpp->reaction_string($reaction);
        my $rate_string = $kpp->rate_string($reaction);
        my ($reactants, $products) = split / = /, $reaction_string;
        $reactants =~ s/\b$ox\b/${ox}_X/;
        my ($non_tagged_reactant, $string);
        if ($reactants =~ /\+/ and $reactants !~ /hv/) { #include non-tagged reactant as product
            ($non_tagged_reactant = $reactants) =~ s/${ox}_X//;
            $non_tagged_reactant =~ s/ \+ //;
        }
        ##tag all Ox products
        my @products = split / \+ /, $products;
        my @final_products;
        foreach my $item (@products) {
            my ($yield, $spc) = split / /, $item;
            unless (defined $spc) {
                $spc = $yield;
                $yield = 1;
            }
            $yield = 1 unless (defined $yield);
            if ($spc ~~ @Ox) {
                $spc .= "_X";
            } elsif ($spc eq "UNITY" or $spc eq "NA") { ## do nothing
            } else {
                next;
            }
            if ($yield == "1") {
                push @final_products, $spc;
            } else {
                push @final_products, "$yield $spc";
            }
        }
        my $final_products = join ' + ', @final_products;
        $final_products = "UNITY" if ($final_products eq "");
        
        if (defined $non_tagged_reactant and $reactants eq "NO2 + O_X" and $final_products eq "NO3_X") {
            $string = "$reactants = $non_tagged_reactant : $rate_string ;";
        } elsif (defined $non_tagged_reactant and $reactants eq "NO2 + O3_X" and $final_products eq "NO3_X") {
            $string = "$reactants = $non_tagged_reactant : $rate_string ;";
        } elsif (defined $non_tagged_reactant and $reactants eq "NO2_X + NO3" and $final_products eq "NO2_X") {
            $string = "$reactants = $non_tagged_reactant : $rate_string ;";
        } elsif (defined $non_tagged_reactant) {
            $string = "$reactants = $non_tagged_reactant + $final_products : $rate_string ;";
        } else {
            $string = "$reactants = $final_products : $rate_string ;";
        }
        $reactions{$string} += 1;
    }
}

my $out_file = "ox_tags.txt";
open my $out, '>:encoding(utf-8)', $out_file or die $!;
print $out "# all reactions involving Ox\n";
foreach my $reaction (sort keys %reactions) {
    if ($reaction =~ /N2O5/) {
        print $out "Change to NO2NO3_X and NO3NO2_X : $reaction\n";
    } elsif ($reaction =~ /ISOPNO3_X \+ NO/ or $reaction =~ /ISOPNO3 \+ NO3_X/) {
        print $out "change NO2 yield: $reaction\n";
    } elsif ($reaction =~ /ONITR \+ NO3_X/ or $reaction =~ /ONITR_X \+ NO3/) {
        print $out "change NO2 yield: $reaction\n";
    } else {
        print $out "$reaction\n";
    }
}
close $out;
