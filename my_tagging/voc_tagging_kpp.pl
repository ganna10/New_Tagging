#! /usr/bin/env perl
# create tagged mechanism from KPP input and also output to KPP format
# Version 0: Jane Coates 1/6/2015

use strict;
use diagnostics; 
use Mechanism;

# Tags which are only applied to HO2 species, not VOCs (eg. inorganic production pathways)
my @ho2_tags = qw(XTR);
# Tags which are only applied to Ox species, not VOCs (eg. to account for initial ozone)
my @ox_tags = (@ho2_tags, qw(INI));
# Assign tags to VOCs
#my @common_voc_tags = qw(TAG);
my %species_tags = (
	CH4      => [ qw(CH4) ],
);
#	BIGALK   => [ qw(BIGALK) ],
#	BIGENE   => [ qw(BIGENE) ],
#	C2H4     => [ qw(C2H4) ],
#	C2H6     => [ qw(C2H6) ],
#	C3H6     => [ qw(C3H6) ],
#	C3H8     => [ qw(C3H8) ],
#	ISOP     => [ qw(ISOP) ],
#	TOLUENE  => [ qw(TOLUENE) ],

# Get a list of all unique VOC tag names
my (%voc_tags, @emitted_species);
foreach my $species (keys %species_tags) {
    push @emitted_species, $species;
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
my @non_chain = qw(M OH HO2 NO NO2 NO3 H2O O2 O3 O O1D CO2 H2O2 HNO3 SO2 HONO N2 UNITY HSO3 NA SA SO3);
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
my $file_base = "my_mozart";
my $mech_file = "${file_base}.eqn";
my $mech = Mechanism->new;
$mech->read_KPP($mech_file);

my (%RO2, %radicals);
#read in gas.spc, radicals.txt and RO2_species.txt
my $spc_file = "${file_base}.spc";

my $radicals_file = "radicals.txt";
open my $radicals_in, '<:encoding(utf-8)', $radicals_file or die $!;
my @radical_lines = <$radicals_in>;
close $radicals_in;
foreach my $line (@radical_lines) {
    chomp $line;
    $radicals{$line} += 1;
}

my $RO2_file = "RO2_species.txt";
open my $RO2_in, '<:encoding(utf-8)', $RO2_file or die $!;
my @RO2_lines = <$RO2_in>;
close $RO2_in;
foreach my $line (@RO2_lines) {
    chomp $line;
    $RO2{$line} += 1;
}

# Follow the chain of VOC oxidation, adding tags to the intermediates
my %chain_reactions;
#my %chain_species;
my (%new_species, %new_radicals, %new_RO2);
foreach my $species (keys %species_tags) {
	my ($consumers) = $mech->consuming($species);
	warn "No reactions consuming $species\n" unless (@$consumers > 0);
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
my $ox_file = "${file_base}_ox_tags.txt";
push @new_reactions,  tag_non_voc($ox_file, \@all_ox_tags, \%tagged_ox, \%new_species);
# Add tagged reactions of the HO2 produced during VOC oxidation
my $ho2_file = "${file_base}_ho2_tags.txt";
push @new_reactions,  tag_non_voc($ho2_file, \@all_ho2_tags, \%tagged_ox, \%new_species);

# A list of pre-tagged reactions which don't fit neatly into the above categories
my $extra_file = "${file_base}_extra_tags.txt";
push @new_reactions,  tag_extra($extra_file, \@all_ox_tags);

my %all_reactions;
my $rxns = $mech->all_reactions();
foreach my $reaction (@$rxns) {
    my $full_string =  $mech->full_string($reaction);
    $all_reactions{$full_string} = 1;
}
$all_reactions{$_} = 1 foreach (@new_reactions);

#emissions of tagged species
foreach my $spc (@emitted_species) {
    my ($emission_rxns) = $mech->producing_from($spc, "UNITY");
    foreach my $reaction (@$emission_rxns) {
        my $string = $mech->reaction_string($reaction);
        my $rate_string = $mech->rate_string($reaction);
        my $label = $reaction;
        foreach my $tag (@voc_tags) {
            $label .= "_$tag";
            $string =~ s/\b($spc)\b/$1_$tag/g;
            my $new_reaction = "{#$label} $string : $rate_string ;\n";
            $all_reactions{$new_reaction} = 1;
        } 
    }
}

#read header of the untagged gas.eqn file, for copying into the new file
open my $eqn_input, '<:encoding(utf-8)', $mech_file or die $!;
my @lines = ();
while (<$eqn_input>) {
    push @lines, $_;
    last if ($_ =~ "#EQUATIONS");
}
close $eqn_input;

foreach my $reaction (sort keys %all_reactions) {
    push @lines, $reaction;
}

#output to files
my $suffix = ".tagged";
my $eqn_output = $mech_file . $suffix;
open my $eqn_out, '>:encoding(utf-8)', $eqn_output or die $!;
print $eqn_out $_ foreach (@lines);
close $eqn_out;

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
foreach (sort keys %new_radicals) {
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
foreach (sort keys %new_RO2) {
    push @lines, "$_\n";
}
open my $RO2_output, '>:encoding(utf-8)', $new_RO2_file or die $!;
print $RO2_output @lines;
close $RO2_output;

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
		my $new_reactants = join ' + ', @new_reactants;
        my ($products) = $mech->products($consumer);
		my %products = map { $_ => 1 } @$products;
		# Tag the products with the root species
		foreach my $product (sort @$products) {
			my $new_product;
            if (defined $non_chain{$product} and $product eq "UNITY") { #include deposition reactions of tagged species
				$new_product = $product;
            } elsif (defined $non_chain{$product}) {
				next;
			} else {
				$new_product = $product . "_$tag";
				$new_species{$new_product} = 1;
				$new_RO2{$new_product} = 1 if (defined $RO2{$product});
				$new_radicals{$new_product} = 1 if (defined $radicals{$product});
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
		if (defined $reactants{NO} and defined $products{NO2}) {
			my $new_no2 = "NO2_X_$tag";
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
			my $new_species = "ONIT_X_$tag";
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
			my $new_species = "ONITR_X_$tag";
			$new_species{$new_species} = 1;
			my $yield = $mech->yield_of('ONITR', [$consumer]);
			$yield = $yield->[0];
			if ($yield == 1) {
				push @new_products, $new_species;
			} else {
				push @new_products, "$yield $new_species";
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
				push @new_products, "$yield $new_ho2";
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
				push @new_products, "$yield $new_species";
			}
		}
		my $new_products = join ' + ', @new_products;
		# Generate the new reaction and save it, keyed by its new reaction label
        my $rate_string = $mech->rate_string($consumer);
		my $new_reaction_string = "{#$label} $new_reactants = $new_products : $rate_string;\n";
		$chain_reactions{$label} = $new_reaction_string;
		# Perform the same action for all chain-products of this reaction
        foreach my $product (@$products) {
            unless ($non_chain{$product} or defined $history{$product}) {
                &follow_chain(@history, $product);
            }
        }
    }
}

# Add tagged reactions for the HO2/Ox produced by the VOC intermediates
sub tag_non_voc {
	my ($file, $tags, $tagged_species, $new_species) = @_;
	my @new_reactions;
    my $kpp = Mechanism->new;
    $kpp = $kpp->read_KPP($file);
    my $rxns = $kpp->all_reactions();
    foreach my $reaction (@$rxns) {
        my $reaction_string = $kpp->reaction_string($reaction);
        my $rate_string = $kpp->rate_string($reaction);
        my $label;
		foreach my $tag (@$tags) { # Make a copy of the ox tagged reactions for each tag
            $label = "${reaction}_$tag";
            foreach my $tagged_species (keys %$tagged_species) {
                $reaction_string =~ s/\b($tagged_species)\b/$1_$tag/g;
				$new_species->{"$1_$tag"} = 1 if defined $1;
			}
            my $new_reaction = "{#$label} $reaction_string : $rate_string ;\n";
            push @new_reactions, $new_reaction;
		}
    }
	return @new_reactions;
}

# Directly add reactions from the specified file, checking that the tagged species actually exist
sub tag_extra {
	my ($file, $ox_tags) = @_;
	my %ox_tags = map { $_ => 1 } @$ox_tags;
	my @new_reactions;

    my $xtr = Mechanism->new;
    $xtr = $xtr->read_KPP($file);
    my $rxns = $xtr->all_reactions();
    foreach my $reaction (@$rxns) {
        my $reaction_string = $xtr->reaction_string($reaction);
        my $rate_string = $xtr->rate_string($reaction);
		my @tags = $reaction_string =~ /X_([A-Z]+)/g;
        my $label;
		foreach my $tag (@tags) { 
            die "Tag \"$tag\" in file $file not previously defined" unless (defined $ox_tags{$tag});
            $label = "${reaction}_$tag";
		}
        my $new_reaction = "{#$label} $reaction_string : $rate_string ;\n";
        push @new_reactions, $new_reaction;
	}
	return @new_reactions;
}
