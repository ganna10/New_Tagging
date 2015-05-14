#! /usr/bin/env perl
# Tagging approach for KPP files using MCM v3.2 inorganic chemistry as per global model technique
# Version 0: Jane Coates xx/xx/xxxx

use strict;
use diagnostics;
use Mechanism;

my @Ox_reservoirs = qw( O3 NO2 NO3 O O1D HO2 N2O5 NO2NO3 NO3NO2 HO2NO2 NO2HO2 HNO3 );
my @do_not_tag = qw( N2 O2 H2O OH NO HONO UNITY NA hv SO2 SO3 HSO3 );
my %non_chain;
$non_chain{$_} += 1 foreach (@do_not_tag);

my %tag_species = ( ## add tags and the species that will be tagged
    INI   => [ qw( O3 CO ) ],
    XTR   => [ qw( HO2 NO2 ) ],
);

my (@tagged_reactions, @non_tagged_reactions);
my $eqn_file = "test.eqn";
my $mech = Mechanism->new();
$mech->read_KPP($eqn_file);

my %chain_reactions;
#my %chain_species;
my %new_species;
foreach my $tag (keys %tag_species) {
    foreach my $species (@{$tag_species{$tag}}) {
        my ($consumers) = $mech->consuming($species);
        warn "No reactions consuming $species\n" unless @$consumers > 0;
        #&follow_chain($tag, $species);
    }
}

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
				#$new_ro2{$new_product} = 1 if defined $ro2{$product};
				#$new_pan{$new_product} = 1 if defined $pan{$product};
			}
			my $yield = $mech->yield_of($product, [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_product;
			} else {
				push @new_products, "$yield * $new_product";
			}
		}
		# Tag the NO to NO2 conversions
		if (defined $reactants{NO} and defined $products{NO2}) {
			my $new_no2 = "NO2_X_$tag";
			$new_species{$new_no2} = 1;
			my $yield = $mech->yield_of('NO2', [$consumer]);
			$yield = $yield->[0];
			if ($yield >= 1) { # Can't convert more than one NO to NO2(?)
				push @new_products, $new_no2;
			} else {
				push @new_products, "$yield * $new_no2";
			}
		}
		# Tag the Ox reservoir species ONIT and ONITR
		if (defined $reactants{NO} and defined $products{ONIT}) {
			my $new_species = "ONIT_X_$tag";
			$new_species{$new_species} = 1;
			my $yield = $mech->yield_of('ONIT', [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_species;
			} else {
				push @new_products, "$yield * $new_species";
			}
		}
		if (defined $reactants{NO} and not defined $reactants{ISOPNO3} and defined $products{ONITR}) {
			my $new_species = "ONITR_X_$tag";
			$new_species{$new_species} = 1;
			my $yield = $mech->yield_of('ONITR', [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_species;
			} else {
				push @new_products, "$yield * $new_species";
			}
		}
		# Tag the produced HO2
		if (defined $products{HO2}) {
			my $new_ho2 = "HO2_X_$tag";
			$new_species{$new_ho2} = 1;
			my $yield = $mech->yield_of('HO2', [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_ho2;
			} else {
				push @new_products, "$yield * $new_ho2";
			}
		}
		# Tag the directly produced O3
		if (defined $products{O3} and not defined $reactants{O3}) { # re-produced O3 is handled elsewhere
			my $new_species = "O3_X_$tag";
			$new_species{$new_species} = 1;
			my $yield = $mech->yield_of('O3', [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_species;
			} else {
				push @new_products, "$yield * $new_species";
			}
		}
		my $new_products = join ' + ', @new_products;
		# Generate the new reaction and save it, keyed by its new reaction label
        my $rate_string = $mech->rate_string($consumer);
		my ($rate1, $rate2) = split /\s*;\s*/, $rate_string;
		if (defined $rate1 and $rate1 ne '') { # generate a new tagged rate label
			$rate1 = tag_rate_label($rate1, $tag);
		}
		my $new_reaction_string = " $rate1 $new_reactants -> $new_products";
		$new_reaction_string .= " ; $rate2" if $rate2 =~ /\S/;
        $chain_reactions{$label} = $new_reaction_string;
		# Perform the same action for all chain-products of this reaction
        foreach my $product (@$products) {
            unless ($non_chain{$product} or defined $history{$product}) {
                &follow_chain(@history, $product);
            }
        }
    }
}
