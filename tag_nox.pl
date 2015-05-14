#!/usr/bin/perl -w
#
use warnings;
use strict;

use Mechanism;

# Tags which are only applied to Ox species, not VOCs (eg. to account for initial ozone)
my @ox_tags = qw(INI XTR STR );
# Tags which are applied to emitted NOx
my @nox_tags = qw(STR LGT AIR BIO BMB ANT);

# Get a list of all tag names
my %nox_tags = map { $_ => 1 } @nox_tags;
my %all_ox_tags = %nox_tags;
$all_ox_tags{$_} = 1 for @ox_tags;
my @all_ox_tags = keys %all_ox_tags;

# Tagged species stubs to track the NOx
my %tagged_nox = (
	NO_ => 'NO',
	NO2_ => 'NO2',
	NO2_X => 'NO2',
	NO3_ => 'NO3',
	NO3_X => 'NO2',
	NO2NO3_ => 'N2O5',
	NO3NO2_ => 'N2O5',
	PAN_ => 'PAN',
	MPAN_ => 'MPAN',
	ISOPNO3_ => 'ISOPNO3',
	ONIT_ => 'ONIT',
	ONIT_X => 'ONIT',
	ONITR_=> 'ONITR',
	ONITR_X => 'ONITR',
	HO2NO2_ => 'HO2NO2',
	HNO3_ => 'HNO3',
);

# Tagged species stubs to track the Ox produced by reactions of NOx
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

my @new_reactions;
my %new_species;
# Tag the reactions of emitted NOx
my $nox_file = "${file_base}_nox_tags.in";
push @new_reactions,  tag_non_voc($nox_file, \@nox_tags, \%tagged_nox, \%new_species);
# Add tagged reactions for the Ox (reactions only relevant when tagging NOx)
my $ox_file = "${file_base}_ox_tags_nox.in";
push @new_reactions,  tag_non_voc($ox_file, \@all_ox_tags, \%tagged_ox, \%new_species);
# Add tagged reactions for the Ox (common to both NOx and VOC tagging)
$ox_file = "${file_base}_ox_tags.in";
push @new_reactions,  tag_non_voc($ox_file, \@all_ox_tags, \%tagged_ox, \%new_species);

# A list of pre-tagged reactions which don't fit neatly into the above categories
my $extra_file = "${file_base}_nox_extra_tags.in";
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
my %tagged_all = (%tagged_nox, %tagged_ox);
my %tag_mapping;
foreach my $species (keys %new_species) {
	my $original = $species;
	$original =~ s/^(.*)_.*?$/$1/; # Remove the last underscore and everything following
	$original .= '_' unless $original =~ /X$/; # Add the trailing underscore back for non-Ox NOx tagged species
	$original = $tagged_all{$original} if defined $tagged_all{$original};
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
my $output_file = "${file_base}_nox_tagged.in";
open FILE, ">$output_file" or die $!;
print FILE $mech_text;
close FILE;

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
				my $replacement;
				if ($tagged_species =~ /_$/) { # Annoying stuff with trailing underscores for non-Ox NOx species
					$replacement = "$tagged_species$tag";
				} else {
					$replacement = "${tagged_species}_$tag";
				}
				if ($reaction =~ s/\s($tagged_species)\s/ $replacement /g) {
					$new_species->{$replacement} = 1;
				}
			}
			if ($reaction =~ s/^(\s*\[.*?\]\s*)//) {
				my $rate = $1;
				if ($reaction =~ /_X_.*->/) {
					$rate = tag_rate_label($rate, "X_$tag");
				} else {
					$rate = tag_rate_label($rate, "$tag");
				}
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
