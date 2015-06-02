#! /usr/bin/perl
# Verion 0 Jane Coates 29/07/2013: Extract a subset of MOZART, only enough reactions to fully oxidise the primary species given on the command line. The 48 "inorganic" reactions are always included
# Version 1 : Jane Coates 14/11/2013 including CO as a tagged species
# Version 2 : Jane Coates 5/6/2014 removing PANs from script, adding radicals and outputs to all files immediately, updating read & write from files

use strict;
use diagnostics;
use KPP;

my @non_chain = qw( OH HO2 NO NO2 NO3 H2O O2 O3 O1D CO2 H2 H2O2 HNO3 SO3 UNITY );
my %non_chain;
$non_chain{$_} = 1 for (@non_chain);

my $RO2_file = "RO2_species.txt";
open my $RO2_list, '<:encoding(utf-8)', $RO2_file or die $!;
my %ro2;
while (<$RO2_list>) {
    $ro2{$_} = 1 for (split /\s+/, $_);
}
close $RO2_list;

my $radicals_file = "radicals.txt";
open my $radicals_list, '<:encoding(utf-8)', $radicals_file or die $!;
my %radicals;
while (<$radicals_list>) {
    $radicals{$_} = 1 for (split /\s+/, $_);
}
close $radicals_list;

my $eqn_file = "gas.eqn";
die "$0: $eqn_file: Not found\n" unless (-f $eqn_file);
my $kpp = KPP->new($eqn_file); 
my $spc_file = "gas.spc";
die "$0: $spc_file: Not found\n" unless (-f $spc_file);

#include the first 48 (inorganic) reactions
my %chain_reactions;
for (1..48) {
    my $string = sprintf "R%03d", $_;
    $chain_reactions{$string} = $kpp->full_string($string);
}

# include emissions and all subsequent reactions of species specified in ARGV
my (%new_species, %new_ro2, %new_radicals);
foreach my $species (@ARGV) {
    my ($emissions) = $kpp->producing_from($species, "UNITY"); #get emission reactions
    foreach my $reaction (@$emissions) {
        $chain_reactions{$reaction} = $kpp->full_string($reaction);
    }
    my ($consumers) = $kpp->consuming($species);
    warn "No reactions consuming $species\n" unless (@$consumers > 0);
    &follow_chain($species);
}

#include NO and NO2 emissions
foreach my $species (@ARGV, 'NO', 'NO2') {
    my ($emissions) = $kpp->producing_from($species, "UNITY");
    foreach my $reaction (@$emissions) {
        $chain_reactions{$reaction} = $kpp->full_string($reaction);
    }
}

#include depostion reactions of some inorganic species, as well as the actual species in the tagged/replaced species
my @replaced_species = map {$_ =~ s/_.*?$//; $_} keys %new_species;
foreach my $species( qw( HNO3 SO2 NO2 O3 H2O2), @replaced_species) {
    my ($deposition) = $kpp->producing_from("UNITY", $species);
    foreach my $reaction (@$deposition) {
        $chain_reactions{$reaction} = $kpp->full_string($reaction) unless (exists $chain_reactions{"${reaction}_$species"});
    }
}

#read header of the untagged gas.eqn file, for copying into the new file
open my $eqn_input, '<:encoding(utf-8)', $eqn_file or die $!;
my @lines = ();
while (<$eqn_input>) {
    push @lines, $_;
    last if ($_ =~ "#EQUATIONS");
}
close $eqn_input;

#add comment line to header
splice (@lines, 0, 1, $lines[0], "// Tagged for @ARGV\n");

#add selected reactions
foreach my $reaction (sort keys %chain_reactions) {
    push @lines, $chain_reactions{$reaction};
}
my $eqn_text = join '', @lines;

#add the new tagged/replaced species to the definition of RO2, if appropriate
my $ro2_text = "";
foreach my $ro2 (keys %new_ro2) {
    $ro2_text .= "      IF (KPP_$ro2 /= 0) RO2 = RO2 + C(KPP_$ro2)\n"
}
$eqn_text =~ s/(RO2 = 0.\n)/$1$ro2_text/;

#write out new eqn file
my $suffix = ".tagged";
my $new_eqn_file = $eqn_file . $suffix;
open my $eqn_output, '>:encoding(utf-8)', $new_eqn_file or die $!;
print $eqn_output $eqn_text;
close $eqn_output;

#write out new spc file, with tagged/replaced species definitions, just appended to the old spc file text
open my $spc_input, '<:encoding(utf-8)', $spc_file or die $!;
@lines = <$spc_input>;
close $spc_input;
for (sort keys %new_species) {
    push @lines, "$_ = IGNORE ; {\@IGNORE} {}\n";
}
my $new_spc_file = $spc_file . $suffix;
open my $spc_output, '>:encoding(utf-8)', $new_spc_file or die $!;
print $spc_output @lines;
close $spc_output;

#write out a file containing a list of the new radicals species resulting from the tagging
my $new_radicals_file = $radicals_file . $suffix;
open my $radicals_input, '<:encoding(utf-8)', $radicals_file or die $!;
@lines = <$radicals_input>;
close $radicals_input;
for (sort keys %new_radicals) {
    push @lines, "$_\n";
}
open my $radicals_output, '>:encoding(utf-8)', $new_radicals_file or die $!;
print $radicals_output @lines;
close $radicals_output;

#write out a file containing a list of the new RO2 species resulting from the tagging
my $new_RO2_file = $RO2_file . $suffix;
open my $RO2_input, '<:encoding(utf-8)', $RO2_file or die $!;
@lines = <$RO2_input>;
close $RO2_input;
for (sort keys %new_ro2) {
    push @lines, "$_\n";
}
open my $RO2_output, '>:encoding(utf-8)', $new_RO2_file or die $!;
print $RO2_output @lines;
close $RO2_output;

#subroutine to recursively follow the oxidation of a parent species, generatin replacement reactions to tag all intermediates
sub follow_chain {
    my @history = @_;
    my $root_species = $history[0];
    my $species = $history[-1];
    my %history;
    $history{$_} = 1 for (@history);

    #loop over the reactions consuming the current species
    my ($consumers) = $kpp->consuming($species);
    foreach my $consumer (@$consumers) {
        my $label = "${consumer}_$root_species";
        next if (defined $chain_reactions{$label});
        my ($reactants) = $kpp->reactants($consumer);
        #tag the reactants with the root species
        my $new_reactants = join ' + ', map {
            if ($_ eq $species and $_ ne $root_species) {
                my $new_species = $_ . "_$root_species";
                $new_species{$new_species} = 1;
                $new_species;
            } else {
                $_;
            }
        } sort @$reactants;

        my ($products) = $kpp->products($consumer);
        #tag the products with the root species
        my @new_products = ();
        foreach my $product (sort @$products) {
            my $new_product;
            if (defined $non_chain{$product}) {
                $new_product = $product;
            } else {
                $new_product = $product . "_$root_species";
                $new_species{$new_product} = 1;
                $new_ro2{$new_product} = 1 if (defined $ro2{$product});
                $new_radicals{$new_product} = 1 if (defined $radicals{$product});
            }
            my $yield = $kpp->yield_of($product, [$consumer]);
            $yield = $yield->[0];
            if ($yield == 1) {
                push @new_products, $new_product;
            } else {
                push @new_products, "$yield $new_product";
            }
        }
        my $new_products = join ' + ', @new_products;

        #generate new reaction and save it, keyed by its reaction label
        my $rate_string = $kpp->rate_string($consumer);
        my $new_reaction_string = "{#$label} $new_reactants = $new_products : $rate_string;\n";
        $chain_reactions{$label} = $new_reaction_string;

        #perform the same action for all chain_products of this reaction
        foreach my $product (@$products) {
            unless ($non_chain{$product} or defined $history{$product}) {
                &follow_chain(@history, $product);
            }
        }
    }
}
