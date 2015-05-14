#!/usr/bin/perl -w
#
use warnings;
use strict;

use Mechanism;

# Tags which are only applied to HO2 species, not VOCs (eg. inorganic production pathways)
my @ho2_tags = qw(XTR STR);
# Tags which are only applied to Ox species, not VOCs (eg. to account for initial ozone)
my @ox_tags = (@ho2_tags, qw(INI STR));
# Assign tags to VOCs
#my @common_voc_tags = qw(TAG);
my %species_tags = (
	BIGALK   => [ qw(ANT BMB) ],
	BIGENE   => [ qw(ANT BMB) ],
	C10H16   => [ qw(BIO) ],
	C2H2     => [ qw(ANT BMB) ],
	C2H4     => [ qw(ANT BMB BIO) ],
	C2H5OH   => [ qw(ANT BMB) ],
	C2H6     => [ qw(ANT BMB BIO INI) ],
	C3H6     => [ qw(ANT BMB BIO) ],
	C3H8     => [ qw(ANT BMB BIO) ],
	CH2O     => [ qw(ANT BMB) ],
	CH3CHO   => [ qw(ANT BMB) ],
	CH3CN    => [ qw(ANT BMB) ],
	CH3COCH3 => [ qw(ANT BMB BIO) ],
	CH3COOH  => [ qw(ANT BMB) ],
	CH3OH    => [ qw(ANT BMB BIO) ],
	CO       => [ qw(ANT BMB BIO INI) ],
	DMS      => [ qw(BIO) ],
	HCN      => [ qw(ANT BMB) ],
	HCOOH    => [ qw(ANT BMB) ],
	ISOP     => [ qw(BIO) ],
	MEK      => [ qw(ANT BMB) ],
	TOLUENE  => [ qw(ANT BMB) ],
	CH4      => [ qw(CH4) ],
	H2       => [ qw(INI) ],
);

# Get a list of all unique VOC tag names
my %voc_tags;
foreach my $species (keys %species_tags) {
	foreach my $tag (@{ $species_tags{$species} }) {
		$voc_tags{$tag} = 1;
	}
}
my @voc_tags = keys %voc_tags;
# ...and finally of all unique tag names
my %all_ox_tags = %voc_tags;
$all_ox_tags{$_} = 1 for @ox_tags;
my @all_ox_tags = keys %all_ox_tags;
my %all_ho2_tags = %voc_tags;
$all_ho2_tags{$_} = 1 for @ho2_tags;
my @all_ho2_tags = keys %all_ho2_tags;

# Species which do not inherit tags from VOC oxidation
my @non_chain = qw(M OH HO2 NO NO2 NO3 H2O O2 O3 O O1D CO2 H2O2 HNO3 SO2);
my %non_chain; $non_chain{$_} = 1 for @non_chain;

# Tagged species stubs to track the Ox and HO2 produced by VOC intermediates
my %tagged_ox = (
	NO2_X     => 'NO2',
	NO3_X     => 'NO3',
	HNO3_X    => 'HNO3',
	HO2NO2_X  => 'HO2NO2',
	NO3NO2_X  => 'N2O5',
	NO2NO3_X  => 'N2O5',
	PAN_X     => 'PAN',
	ONIT_X    => 'ONIT',
	MPAN_X    => 'MPAN',
	ISOPNO3_X => 'ISOPNO3',
	ONITR_X   => 'ONITR',
	NH4NO3_X  => 'NH4NO3',
	O_X       => 'O',
	O1D_X     => 'O1D',
	O3_X      => 'O3',
	HO2_X     => 'HO2',
	NO2HO2_X  => 'HO2NO2',
);

# Import the mechanism file into the Mechanism.pm module
my $file_base = "chem_mech";
my $mech_file = "${file_base}.in";
my $mech = Mechanism->new;
$mech->read_MOZART($mech_file);
# Also read in the mechanism file into a normal string
open FILE, $mech_file or die $!;
my $mech_text = join '', <FILE>;
close FILE;

# Read in the species mapping from the "Solution" section of the mechanism file
my %solution_mapping;
my ($solution_text) = $mech_text =~ /Solution\s*(.*?)\s*End Solution/s;
my @solution_species = split /,|\n/, $solution_text;
foreach my $species (@solution_species) {
	$species =~ s/\s+//g;
	my ($alias, $original) = split '->', $species;
	die "Problem reading solution string" unless defined $alias;
	$original = $alias unless defined $original;
	$solution_mapping{$alias} = $original;
}

# Follow the chain of VOC oxidation, adding tags to the intermediates
my %chain_reactions;
#my %chain_species;
my %new_species;
foreach my $species (keys %species_tags) {
	my ($consumers) = $mech->consuming($species);
	warn "No reactions consuming $species\n" unless @$consumers > 0;
	foreach my $tag (@{ $species_tags{$species} }) {
		&follow_chain($tag, $species);
	}
}

# Get a list of the new VOC tagged reactions
my @new_reactions;
foreach my $reaction (sort keys %chain_reactions) {
	push @new_reactions, $chain_reactions{$reaction};
}

# Add tagged reactions for the Ox which are specific for VOC tagging
my $ox_file = "${file_base}_ox_tags_voc.in";
push @new_reactions,  tag_non_voc($ox_file, \@all_ox_tags, \%tagged_ox, \%new_species);
# Add tagged reactions for the Ox produced by the VOC intermediates
$ox_file = "${file_base}_ox_tags.in";
push @new_reactions,  tag_non_voc($ox_file, \@all_ox_tags, \%tagged_ox, \%new_species);
# Add tagged reactions of the HO2 produced during VOC oxidation
my $ho2_file = "${file_base}_ho2_tags.in";
push @new_reactions,  tag_non_voc($ho2_file, \@all_ho2_tags, \%tagged_ox, \%new_species);

# A list of pre-tagged reactions which don't fit neatly into the above categories
my $extra_file = "${file_base}_voc_extra_tags.in";
push @new_reactions,  tag_extra($extra_file, \@all_ox_tags);

# Deal with long lines
my $maxchars = 120;
foreach my $l (0..$#new_reactions) {
	my $line = $new_reactions[$l];
	if (length $line > $maxchars) {
		my @newlines;
		my ($reaction, $rate) = split ';', $line;
		my ($reactants, $products) = split ' -> ', $reaction;
		die "Unable to split reaction into reactants and products" unless defined $reactants and defined $products;
		my @products = split '\+', $products;
		$newlines[0] = "$reactants -> $products[0]";
		$newlines[0] .= ";$rate" if defined $rate;
		die "Can't fit reaction into $maxchars characters:\n$line" if length($newlines[0]) > $maxchars;
		shift @products;
		while (@products) {
			my $newline = ' ' x 6;
			while (@products and length($newline) + length($products[0]) + 1 < $maxchars) {
				$newline .= '+' . shift @products;
			}
			push @newlines, $newline;
		}
		$new_reactions[$l] = join "\n", @newlines;
	}
}

# Split the new reactions into photolysis and thermal reactions
my (@new_photolysis, @new_thermal);
foreach my $reaction (@new_reactions) {
	if ($reaction =~ /hv/) {
		push @new_photolysis, $reaction;
	} else {
		push @new_thermal, $reaction;
	}
}

# Get a list of all newly-added species, and a mapping to the "original" species in the mechanism
my %tag_mapping;
foreach my $species (keys %new_species) {
	my $original = $species;
	$original =~ s/^(.*)_.*?$/$1/; # Remove the last underscore and everything following
	$original = $tagged_ox{$original} if defined $tagged_ox{$original};
	die "Don't know how to map species $species ($original)" unless defined $solution_mapping{$original};
	$tag_mapping{$species} = $solution_mapping{$original};
}
my @new_species = sort keys %tag_mapping;

# Generate the new strings to add to the mechanism file
my $solution_string = '';
$solution_string   .= "$_ -> $tag_mapping{$_}\n" for @new_species;
#my $implicit_string = join "\n", @new_species;
my $photo_string    = join "\n", @new_photolysis;
my $thermal_string  = join "\n", @new_thermal;

# Add the newly-generated strings to the mechanism file and write it out
$mech_text =~ s/^(\s*End Solution)/$solution_string$1/m;
#$mech_text =~ s/^(\s*End Implicit)/$implicit_string\n$1/m;
$mech_text =~ s/^(\s*End Photolysis)/$photo_string\n$1/m;
$mech_text =~ s/^(\s*End Reactions)/$thermal_string\n$1/m;
my $output_file = "${file_base}_voc_tagged.in";
open FILE, ">$output_file" or die $!;
print FILE $mech_text;
close FILE;

# Recursively follow the oxidation of a parent species, generating
# additional reactions to tag all intermediates
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

# Tag a rate inside square brackets
sub tag_rate_label {
	my ($rate, $tag) = @_;
	$rate =~ /\[(.*?)->,(.*?)\]/;
	if (defined $1 and defined $2) { # an aliased photo rate
		$rate = "[${1}_${tag}->,$2]";
	} elsif ($rate =~ /\[(j.*?)\]/) { # an unaliased photo rate
		$rate = "[${1}_$tag->,$1]";
	} else { # some other kind of rate
		$rate =~ s/\[(.*?)\]/[${1}_$tag]/;
	}
	return $rate;
}

# Add tagged reactions for the HO2/Ox produced by the VOC intermediates
sub tag_non_voc {
	my ($file, $tags, $tagged_species, $new_species) = @_;
	my @new_reactions;
	open FILE, $file or die $!;
	foreach my $line (<FILE>) {
		next if $line =~ /^\s*\*/;
		next unless $line =~ /\S/;
		die "$file: contains invalid reactions" unless $line =~ /->/;
		# Make sure potentially problematic characters are surrounded by space
		$line =~ s/\*/ * /g;
		$line =~ s/;/ ; /g;
		$line =~ s/\+/ + /g;
		my $reaction;
		foreach my $tag (@$tags) { # Make a copy of the ox tagged reactions for each tag
			$reaction = $line;
			foreach my $tagged_species (keys %$tagged_species) {
				$reaction =~ s/\s($tagged_species)\s/ $1_$tag /g;
				$new_species->{"$1_$tag"} = 1 if defined $1;
			}
			if ($reaction =~ s/^(\s*\[.*?\]\s*)//) {
				my $rate = tag_rate_label($1, "X_$tag");
				$reaction = " $rate $reaction";
			}
			# trim the spaces introduced above
			$reaction =~ s/\s+/ /g;
			push @new_reactions, $reaction;
		}
	}
	close FILE;
	return @new_reactions;
}

# Directly add reactions from the specified file, checking that the tagged species actually exist
sub tag_extra {
	my ($file, $ox_tags) = @_;
	my %ox_tags = map { $_ => 1 } @$ox_tags;
	my @new_reactions;
	open FILE, $file or die $!;
	foreach my $line (<FILE>) {
		next if $line =~ /^\s*\*/;
		next unless $line =~ /\S/;
		die "$file: contains invalid reactions" unless $line =~ /->/;
		my @tags = $line =~ /X_([A-Z]+)/g;
		foreach my $tag (@tags) {
			die "Tag \"$tag\" in file $file not previously defined" unless defined $ox_tags{$tag};
		}
		push @new_reactions, $line;
	}
	close FILE;
	return @new_reactions;
}
