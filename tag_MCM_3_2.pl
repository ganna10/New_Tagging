#! /usr/bin/perl
# Jane Coates: Extract a subset of the MCM, only enough reactions to fully oxidise the primary species given on the command line. The 48 "inorganic" reactions are always included
# Version 1 : Jane Coates 17/11/2013 removed CO from non-tagged list and adapted tagging subroutine accordingly
# Version 2: Jane Coates 28/01/2014 included reactions of tagged CO with criegee biradicals
# Version 3: Jane Coates 07/05/2014 removing output of TOTPAN as it's not needed, updating output to files and including output of tagged radicals too
# Vresion 4: Jane Coates 4/7/2014 more diradicals included and sorted out the repeated criegeee reactions with CO

use strict;
use diagnostics;
use KPP;

my @non_chain = qw( OH HO2 NO NO2 NO3 H2O O2 O3 O3P O1D CO2 HCL HBR HI CLO BRO IO CL BR I CL2 BR2 I2 BRCL HOBR H2 H2O2 HNO3 SO3 UNITY );
my %non_chain;
$non_chain{$_} = 1 for (@non_chain);

my $RO2_file = "RO2_species.txt";
my %ro2 = read_slurp($RO2_file);

my $radicals_file = "radicals.txt";
my %radicals = read_slurp($radicals_file);

my $criegee_file = "diradicals.txt";
my %criegee = read_slurp($criegee_file);

my $eqn_file = "gas.eqn";
die "$0: $eqn_file: Not found\n" unless (-f $eqn_file);
my $kpp = KPP->new($eqn_file);

my $spc_file = "gas.spc";
die "$0: $spc_file: Not found\n" unless (-f $spc_file);

#include the first 48 (inorganic) reactions
my %chain_reactions;
for (1..48) {
    my $string = sprintf "R%05d", $_;
    $chain_reactions{$string} = $kpp->full_string($string);
}

#included tagged R19 - CO + OH = HO2
foreach my $parent (@ARGV) {
    my $number = "R00019";
    my $reaction = $kpp->full_string($number);
    $reaction =~ s/(CO)/$1_$parent/;
    $reaction =~ s/($number)/$1_$parent/;
    $chain_reactions{"${number}_$parent"} = $reaction;
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

#include reactions of tagged CO with criegee biradicals
my @tagged_CO;
push @tagged_CO, 'CO_'.$_ foreach (@ARGV);
foreach my $radical (sort keys %criegee) {
    foreach my $parent (@ARGV) {
        my $tagged_criegee = $radical . "_" . $parent;
        next unless (exists $new_species{$tagged_criegee});
        my $CO_reaction_nr = $kpp->reacting_with('CO', $radical);
        next unless (exists $CO_reaction_nr->[0]);
        my $CO_reaction = $kpp->full_string($CO_reaction_nr->[0]);
        my $counter = 1;
        foreach my $CO (@tagged_CO) { 
            my ($r_nr, $reactants, $products, $rate_constant) = $CO_reaction =~ /^{#(.*?)}(.*?)=(.*?):(.*?)$/;
            if ($reactants =~ /(CO \+)/) {
                $reactants =~ s/(CO \+ $radical)/$CO \+ $tagged_criegee/;
            } elsif ($reactants =~ /(\+ CO)/) {
                $reactants =~ s/($radical \+ CO)/$tagged_criegee \+ $CO/;
            }
            $products =~ s/\s$//;
            $chain_reactions{"${r_nr}_${parent}_$counter"} = "{#${r_nr}_${parent}_$counter}$reactants=${products}_${parent} :$rate_constant\n";
            $counter++;
        }
    }
}

#include NO and NO2 emissions
foreach my $species (@ARGV, 'NO', 'NO2') {
    my ($emissions) = $kpp->producing_from($species, "UNITY");
    $chain_reactions{$_} = $kpp->full_string($_) foreach (@$emissions);
}

#include depostion reactions of some inorganic species, as well as the actual species in the tagged/replaced species
my @replaced_species = map {$_ =~ s/_.*?$//; $_} keys %new_species;
foreach my $species( qw( HNO3 SO2 NO2 O3 H2O2), @replaced_species) {
    my ($deposition) = $kpp->producing_from("UNITY", $species);
    foreach my $reaction (@$deposition) {
        $chain_reactions{$reaction} = $kpp->full_string($reaction) unless (exists $chain_reactions{"${reaction}_$species"});
    }
}

#read header of the untagged gas.eqn file, for cpying into the new file
open my $eqn_in, '<:encoding(utf-8)', $eqn_file or die $!;
my @lines = ();
for (<$eqn_in>) {
    push @lines, $_;
    last if ($_ =~ "#EQUATIONS");
}
close $eqn_in;

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
open my $eqn_out, '>:encoding(utf-8)', $new_eqn_file or die $!;
print $eqn_out $eqn_text;
close $eqn_out;

#write out new spc file, with tagged/replaced species definitions, just appended to the old spc file text
open my $spc_in, '<:encoding(utf-8)', $spc_file or die $!;
@lines = <$spc_in>;
close $spc_in;
for (sort keys %new_species) {
    push @lines, "$_ = IGNORE ; {\@IGNORE} {}\n";
}
my $new_spc_file = $spc_file . $suffix;
open my $spc_out, '>:encoding(utf-8)', $new_spc_file or die $!;
print $spc_out @lines;
close $spc_out;

#write out a file containing a list of the new radicals species resulting from the tagging
open my $radicals_in, '<:encoding(utf-8)', $radicals_file or die $!;
@lines = <$radicals_in>;
close $radicals_in;
push @lines, "$_\n" for (sort keys %new_radicals);
my $new_radicals_file = $radicals_file . $suffix;
open my $radicals_out, '>:encoding(utf-8)', $new_radicals_file or die $!;
print $radicals_out @lines;
close $radicals_out;

#write out a file containing a list of the new RO2 species resulting from the tagging
open my $RO2_in, '<:encoding(utf-8)', $RO2_file or die $!;
@lines = <$RO2_in>;
close $RO2_in;
push @lines, "$_\n" for (sort keys %new_ro2);
my $new_ro2_file = $RO2_file . $suffix;
open my $RO2_out, '>:encoding(utf-8)', $new_ro2_file or die $!;
print $RO2_out @lines;
close $RO2_out;

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
            unless ($non_chain{$product} or defined $history{$product} or $product eq "CO") {
                &follow_chain(@history, $product);
            }
        }
    }
}

sub read_slurp {
    my ($file) = @_;

    my %data;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    while (<$in>) {
        $data{$_} = 1 for (split /\s+/, $_) ;
    }
    return %data;
}
