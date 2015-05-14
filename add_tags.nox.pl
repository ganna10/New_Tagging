#!/usr/bin/perl -w

# Add the specified tags to the CESM source code, based on template source
# code files
# Tim Butler <tim.butler@iass-potsdam.de> May 2014, adapted from the nox tagging script
# Tim Butler: August 2014, adapted to apply to *either* NOx or VOC tagging (using new fortran comment structure)

use warnings;
use strict;

# Choose NOX or VOC tagging
my $what_to_tag = 'NOX';
die "Can only tag NOX or VOC" unless $what_to_tag eq 'NOX' or $what_to_tag eq 'VOC';

my $model_version = "cesm1_1_1_ccmi23";

# These are the tags to add
my @uctags = qw(TAG LGT STR AIR BIO BMB ANT);
my @uc_x_tags = (@uctags, qw(INI XTR));

# Make a list of these in lowercase, too
my @lctags = @uctags;
my @lc_x_tags = @uc_x_tags;
@lctags = map { $_ =~ tr/A-Z/a-z/; $_ } @lctags;
@lc_x_tags = map { $_ =~ tr/A-Z/a-z/; $_ } @lc_x_tags;

# Get a list of all of the files to which we must add tags
my $template_dir = "./templates/$model_version";
opendir DIR, $template_dir or die $!;
my @template_sources = grep /src/, readdir DIR;
closedir DIR;
die "No template source directories found" unless @template_sources;
my @template_files = ();
foreach my $dir (@template_sources) {
	opendir DIR, "$template_dir/$dir" or die $!;
	push @template_files, map { "$dir/$_" } grep /F90$/, readdir DIR;
	closedir DIR;
}
die "No template files found" unless @template_files;

# Make sure the SourceMods structure is in place
my $sourcemods_dir = "./SourceMods";
if (! -d $sourcemods_dir) {
	mkdir $sourcemods_dir or die "Unable to create $sourcemods_dir: $!\n";
}
foreach my $dir (@template_sources) {
	if (! -d "$sourcemods_dir/$dir") {
		mkdir "$sourcemods_dir/$dir" or die "Unable to create $sourcemods_dir/$dir: $!\n";
	}
}

# Copy over each of the template files, adding the requested tags.
foreach my $file (@template_files) {
	open FILE, "$template_dir/$file" or die $!;
	my @lines = <FILE>;
	close FILE;
	open FILE, ">$sourcemods_dir/$file" or die $!;
	for (my $l = 0; $l < @lines; $l++) {
		my $line = $lines[$l];
		if ($line =~ /START $what_to_tag TAGGING CODE/) {
			my $template_text;
			for ($l++; $lines[$l] !~ /END $what_to_tag TAGGING CODE/; $l++) {
				$template_text .= $lines[$l];
			}
			$l--;
			$line .= replacement_text($template_text, \@uctags, \@lctags);
		} elsif ($line =~ /START TAGGING CODE/) {
			my $template_text;
			for ($l++; $lines[$l] !~ /END TAGGING CODE/; $l++) {
				$template_text .= $lines[$l];
			}
			$l--;
			$line .= replacement_text($template_text, \@uc_x_tags, \@lc_x_tags);
		}
		print FILE $line;
	}
	close FILE;
}

sub replacement_text {
	my ($template_text, $uctags, $lctags) = @_;
	my $replacement_text = "";
	foreach my $t (0..$#$uctags) {
		my $tmp = $template_text;
		$tmp =~ s/_TAG/_$uctags->[$t]/sg;
		$tmp =~ s/_tag/_$lctags->[$t]/sg;
		$replacement_text .= $tmp;
	}
	return $replacement_text;
}
