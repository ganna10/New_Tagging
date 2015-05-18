#! /usr/bin/env perl
# Tagging approach for KPP files using MCM v3.2 inorganic chemistry as per global model technique
# Version 0: Jane Coates 14/5/2015

use strict;
use diagnostics;
use Mechanism;

my $Ox_file = "MOZART_Ox_species.txt";
my @Ox_species; ###use later to tag these species as _X_TAG
open my $Ox_in, '<:encoding(utf-8)', $Ox_file or die $!;
my @Ox_lines = <$Ox_in>;
close $Ox_in;
foreach my $line (@Ox_lines) {
     chomp $line;
     $line =~ s/^\s+|\s+$//g;
     push @Ox_species, $line;
}

my @do_not_tag = qw( H2O2 N2 O2 H2O OH NO HONO UNITY NA hv SO2 SO3 HSO3 );
my %non_chain;
$non_chain{$_} += 1 foreach (@do_not_tag);

my %tag_species = ( ## add tags and the species that will be tagged
    CH4     => [ qw( CH4 ) ],
    #INI   => [ qw( TOLUENE ) ],
    #XTR   => [ qw( HO2 NO2 ) ],
);

my (@tagged_reactions, @non_tagged_reactions);
my $eqn_file = "test.eqn";
my $mech = Mechanism->new();
$mech->read_KPP($eqn_file);

my (%ro2, %radical);
open my $ro2_in, '<:encoding(utf-8)', "RO2_species.txt" or die $!;
my @ro2 = split /\s/, <$ro2_in>;
close $ro2_in;
$ro2{$_} = 1 foreach (@ro2);

open my $radical_in, '<:encoding(utf-8)', "radicals.txt" or die $!;
my @radical = split /\s/, <$radical_in>;
close $radical_in;
$radical{$_} = 1 foreach (@radical);

my (%chain_reactions, %new_species, %new_ro2, %new_radical);
foreach my $tag (keys %tag_species) {
    foreach my $species (@{$tag_species{$tag}}) {
        my ($consumers) = $mech->consuming($species);
        warn "No reactions consuming $species\n" unless @$consumers > 0;
        &follow_chain($tag, $species);
    }
}
foreach my $reaction (sort keys %chain_reactions) {
    print "$chain_reactions{$reaction}\n";
}

#include NO and NO2 emissions 
####include emissions of tagged species
#foreach my $species ('NO', 'NO2') {
#    my ($emissions) = $kpp->producing_from($species, "UNITY");
#    $chain_reactions{$_} = $kpp->full_string($_) foreach (@$emissions);
#}
#
##include depostion reactions of some inorganic species, as well as the actual species in the tagged/replaced species
#####include tagged species
#my @replaced_species = map {$_ =~ s/_.*?$//; $_} keys %new_species;
#foreach my $species( qw( HNO3 SO2 NO2 O3 H2O2), @replaced_species) {
#    my ($deposition) = $kpp->producing_from("UNITY", $species);
#    foreach my $reaction (@$deposition) {
#        $chain_reactions{$reaction} = $kpp->full_string($reaction) unless (exists $chain_reactions{"${reaction}_$species"});
#    }
#}

#read header of the untagged gas.eqn file, for copying into the new file
open my $eqn_in, '<:encoding(utf-8)', $eqn_file or die $!;
my @lines = ();
for (<$eqn_in>) {
    push @lines, $_;
    last if ($_ =~ "#EQUATIONS");
}
close $eqn_in;

#add comment line to header
splice (@lines, 0, 1, $lines[0], "// Tagged for testing\n");

#add selected reactions -> need to include all original inorganic and tagged species chemistry
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

##write out new eqn file
my $suffix = ".tagged";
#my $new_eqn_file = $eqn_file . $suffix;
#open my $eqn_out, '>:encoding(utf-8)', $new_eqn_file or die $!;
#print $eqn_out $eqn_text;
#close $eqn_out;
#
##write out new spc file, with tagged/replaced species definitions, just appended to the old spc file text
#open my $spc_in, '<:encoding(utf-8)', $spc_file or die $!;
#@lines = <$spc_in>;
#close $spc_in;
#for (sort keys %new_species) {
#    push @lines, "$_ = IGNORE ; {\@IGNORE} {}\n";
#}
#my $new_spc_file = $spc_file . $suffix;
#open my $spc_out, '>:encoding(utf-8)', $new_spc_file or die $!;
#print $spc_out @lines;
#close $spc_out;

#write out a file containing a list of the new radicals species resulting from the tagging
$radical{$_} = 1 foreach (sort keys %new_radical);
my $new_radicals_file = "radicals.txt" . $suffix;
open my $radicals_out, '>:encoding(utf-8)', $new_radicals_file or die $!;
print $radicals_out "$_\n" foreach (sort keys %radical);
close $radicals_out;

#write out a file containing a list of the new RO2 species resulting from the tagging
$ro2{$_} = 1 foreach (sort keys %new_ro2);
my $new_ro2_file = "RO2_species.txt" . $suffix;
open my $RO2_out, '>:encoding(utf-8)', $new_ro2_file or die $!;
print $RO2_out "$_\n" foreach (sort keys %ro2);
close $RO2_out; 

# Recursively follow the oxidation of a parent species, generating additional reactions to tag all intermediates
sub follow_chain {
    my @history = @_;
    my $tag = $history[0];
    my $species = $history[-1];
    my %history;
    $history{$_} = 1 for @history;
	#print "@history\n";
	# Loop over the reactions comsuming the current species
    my ($consumers) = $mech->consuming($species);
    foreach my $consumer (@$consumers) {
        my $label = "${consumer}_$tag";
		next if defined $chain_reactions{$label};
        my ($reactants) = $mech->reactants($consumer);
		my %reactants = map { $_ => 1 } @$reactants;
		# Tag the reactants
		my @new_reactants = ();
		my @new_products = ();
		foreach my $reactant (sort @$reactants) {
			if ($reactant eq $species) {
				my $new_species = $reactant . "_$tag";
				$new_species{$new_species} = 1;
				#$chain_species{$reactant} = 1;
				push @new_reactants, $new_species
			} else {
				#die "Reaction of two chain species" if defined $chain_species{$reactant};
				push @new_reactants, $reactant;
				push @new_products, $reactant unless $reactant eq 'hv';
			}
		}
		#if (@$reactants > 1 and $reactants->[0] eq $reactants->[1] and defined $chain_species{$reactants->[0]}) {
		#die "Self-reaction of a chain species";
		#}
		my $new_reactants = join ' + ', @new_reactants;
        my ($products) = $mech->products($consumer);
		my %products = map { $_ => 1 } @$products;
		# Tag the products with the root species
		foreach my $product (sort @$products) {
			my $new_product;
			if (defined $non_chain{$product}) {
				#$new_product = $product;
				next;
			} else {
				$new_product = $product . "_$tag";
				$new_species{$new_product} = 1;
				$new_ro2{$new_product} = 1 if defined $ro2{$product};
				$new_radical{$new_product} = 1 if defined $radical{$product};
			}
			my $yield = $mech->yield_of($product, [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_product;
			} else {
				push @new_products, "$yield $new_product";
			}
		}
        # Tag the NO to NO2 conversions
		if (defined $reactants{NO} and defined $products{NO2}) { #### jc: repeated NO2_tag products??
			my $new_no2 = "NO2_$tag";
			$new_species{$new_no2} = 1;
			my $yield = $mech->yield_of('NO2', [$consumer]);
			$yield = $yield->[0];
			if ($yield >= 1) { # Can't convert more than one NO to NO2(?)
				push @new_products, $new_no2;
			} else {
				push @new_products, "$yield $new_no2";
			}
		}
		# Tag the Ox reservoir species ONIT and ONITR
		if (defined $reactants{NO} and defined $products{ONIT}) {
			my $new_species = "ONIT_$tag";
			$new_species{$new_species} = 1;
			my $yield = $mech->yield_of('ONIT', [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_species;
			} else {
				push @new_products, "$yield $new_species";
			}
		}
		if (defined $reactants{NO} and not defined $reactants{ISOPNO3} and defined $products{ONITR}) {
			my $new_species = "ONITR_$tag";
			$new_species{$new_species} = 1;
			my $yield = $mech->yield_of('ONITR', [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_species;
			} else {
				push @new_products, "$yield $new_species";
			}
		}
		# Tag the produced HO2 -> repeating the produced yield * HO2_TAG
		if (defined $products{HO2}) {
			my $new_ho2 = "HO2_$tag";
			$new_species{$new_ho2} = 1;
			my $yield = $mech->yield_of('HO2', [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_ho2;
			} else {
				push @new_products, "$yield $new_ho2";
			}
		}
		# Tag the directly produced O3
		if (defined $products{O3} and not defined $reactants{O3}) { # re-produced O3 is handled elsewhere
			my $new_species = "O3_$tag";
			$new_species{$new_species} = 1;
			my $yield = $mech->yield_of('O3', [$consumer]);
			$yield = $yield->[0];
            if ($yield == 1) {
				push @new_products, $new_species;
			} else {
				push @new_products, "$yield $new_species";
			}
		}
		my $new_products = join ' + ', @new_products;
        # Generate the new reaction and save it, keyed by its new reaction label
        my $rate_string = $mech->rate_string($consumer);
        my $new_reaction_string = "{#$label} $new_reactants = $new_products : $rate_string;\n";
        print $new_reaction_string;
		$chain_reactions{$label} = $new_reaction_string;
		# Perform the same action for all chain-products of this reaction
        foreach my $product (@$products) {
            unless ($non_chain{$product} or defined $history{$product}) {
                &follow_chain(@history, $product);
            }
        }
    }
}
