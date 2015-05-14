#!/usr/bin/env perl
#write out .eqn, .spc and mecca.f90 module from tagged MOZART mechanism using only inorganic and ch4 chemistry
#Version 0: Jane Coates 12/5/2015

use strict;
use diagnostics;

my @RO2_species = qw( HO2 HO2_INI HO2_XTR HO2_CH4 CH3O2 CH3O2_CH4 HOCH2OO HOCH2OO_CH4 );

my %species;

my @reactions;
my $photo_max = 0;
my (@photo_rate, @photo_labels, %photolysis);

#get tagged MCM inorganic reactions
open my $inorganic_in, "<:encoding(utf-8)", "MCM-3.2_inorganic_reactions_tagged.txt" or die $!;
my @mcm_inorganic = (<$inorganic_in>);
close $inorganic_in;
$_ =~ s/^\s+|\s+$//g foreach (@mcm_inorganic);
#print ":$_:\n" foreach @mcm_inorganic;
foreach my $reaction (@mcm_inorganic) {
    push @reactions, $reaction;
    my ($reactants, $products, $rate) = $reaction =~ /(.*?) = (.*?) : (.*?)$/;
    my @reactants = split / \+ /, $reactants;
    my @products = split / \+ /, $products;
    $species{$_} += 1 foreach (@reactants);
    $species{$_} += 1 foreach (@products);
    if ($reaction =~ /hv/) { ## process photolysis labels
        my ($label) = $rate =~ /JX\((.*?)\)/;
        $photolysis{$label} += 1;
    }
}

# get tagged and non-tagged CH4 chemistry
open my $ch4_in, '<:encoding(utf-8)', "ch4_chemistry_tagged_non-tagged.txt" or die $!;
my @ch4_in = <$ch4_in>;
close $ch4_in;
foreach my $reaction (@ch4_in) {
    chomp $reaction;
    push @reactions, $reaction;
    my ($reactants, $products, $rate) = $reaction =~ /(.*?) = (.*?) : (.*?)$/;
    my @reactants = split / \+ /, $reactants;
    my @products = split / \+ /, $products;
    $_ =~ s/^\d?\.?\d+\s// foreach (@products); #remove stoichiometeries
    $species{$_} += 1 foreach (@reactants);
    $species{$_} += 1 foreach (@products);
    if ($reaction =~ /hv/) { ## process photolysis labels
        my ($label) = $rate =~ /JX\((.*?)\)/;
        $photolysis{$label} += 1;
    }
}

delete $species{"hv"};
my @species = sort keys %species;

foreach my $label (sort keys %photolysis) {
    $photo_labels[$photo_max] = $label; #add JMCM<> value to photo_rate
    if ($label eq "ip_O3_O1D") {
        $photo_rate[$photo_max] = "J_MCM(1)";
    } elsif ($label eq "ip_O3_O") {
        $photo_rate[$photo_max] = "J_MCM(2)";
    } elsif ($label eq "ip_H2O2_OH_OH") {
        $photo_rate[$photo_max] = "J_MCM(3)";
    } elsif ($label eq "ip_NO2_NO_O") {
        $photo_rate[$photo_max] = "J_MCM(4)";
    } elsif ($label eq "ip_NO3_NO") {
        $photo_rate[$photo_max] = "J_MCM(5)";
    } elsif ($label eq "ip_NO3_NO2_O") {
        $photo_rate[$photo_max] = "J_MCM(6)";
    } elsif ($label eq "ip_HONO_NO_OH") {
        $photo_rate[$photo_max] = "J_MCM(7)";
    } elsif ($label eq "ip_HNO3_NO2_OH") {
        $photo_rate[$photo_max] = "J_MCM(8)";
    } elsif ($label eq "ip_CH2O_CO_H2") {
        $photo_rate[$photo_max] = "J_MCM(12)";
    } elsif ($label eq "ip_CH2O_CO_HO2") {
        $photo_rate[$photo_max] = "J_MCM(11)";
    } elsif ($label eq "ip_CH3OOH_CH2O_HO2") {
        $photo_rate[$photo_max] = "J_MCM(41)";
    } elsif ($label eq "ip_CH3OOH_CH2O_HO2_OH") {
        $photo_rate[$photo_max] = "J_MCM(41)";
    }
    $photo_max++;
}

# Extract photolysis expressions
my $photo_file = "MCM-3.1_photolysis.fac";
open my $photo_in, '<:encoding(utf-8)', $photo_file or die $!;
my @photo_lines = <$photo_in>;
close $photo_in;
my (@j_mcm, @jval_a, @jval_b, @jval_c);
foreach my $line (@photo_lines) {
    my ($j, $a, $b, $c) = $line =~ /([\dD.-]+)/g;
    next unless ( defined $j and defined $a and defined $b and defined $c ); 
    push @j_mcm, $j; 
    $jval_a[$j] = $a; 
    $jval_b[$j] = $b; 
    $jval_c[$j] = $c; 
}
my $n_j_mozart = $#jval_a + 1;

#add deposition
open my $dep_in, "<:encoding(utf-8)", "MOZART_deposition_velocities.txt.ch4" or die $!;
my @deposition_lines = grep /^%/, <$dep_in>;
close $dep_in;
my %deposition_velocities;

foreach my $reaction (@deposition_lines) {
    my ($rate, $reactants) = $reaction =~ /^%(.*?):(.*?)=/;
    $reactants =~ s/^\s+|\s+$//g;
    ($deposition_velocities{$reactants}) = $rate =~ /([\d.]+)/;
    my $product = 'UNITY';
    $rate = "MOZART_VD(KPP_$reactants)/(zmbl*100.)";
    push @reactions, "$reactants = $product : $rate";
} 
#print "$_ : $deposition_velocities{$_}\n" foreach keys %deposition_velocities;

#add emissions
#open my $primary_in, '<:encoding(utf-8)', "primary_species.txt" or die $!;
#my @primary_species = <$primary_in>;
#close $primary_in;
#chomp $_ foreach (@primary_species);
my @primary_species = qw( NO NO2 CO CH4 CH4_CH4 );
push @reactions, "UNITY = $_ : MOZART_EMIS(KPP_$_)/(zmbl*100.)" foreach (@primary_species) ;
#print "$_\n" foreach @reactions;

#get MCM rate expressions and to be included at beginning of file
open my $rate_in, '<:encoding(utf-8)', "mcm_rate_variables.txt" or die $!;
my @mcm_exps = <$rate_in>;
close $rate_in;
chomp $_ foreach (@mcm_exps);

my @mcm_vars;
foreach (@mcm_exps) {
    my ($var, $exp) = split / = /, $_;
    push @mcm_vars, $var;
}

my $no_of_reactions = @reactions;
my $label_length = length $no_of_reactions;

#output files
my $eqn_file = "mozart.eqn.ch4";
my $mecca_file = "mecca_mozart.f90.ch4";
my $spc_file = "mozart.spc.ch4";

#write spc file
open my $out_spc, ">:encoding(utf-8)", $spc_file or die $!;
print $out_spc "{ Created automatically from tagged inorganic and ch4 chemistry }\n";
print $out_spc "#DEFVAR\n";
print $out_spc "$_ = IGNORE ; {\@IGNORE} {}\n" for @species;
close $out_spc;

#Write eqn file
open my $eqn_out, '>:encoding(utf-8)', $eqn_file or die $!;
print $eqn_out "// inorganic and methane chemistry \n";
print $eqn_out "#INLINE F95_DECL\n";
print $eqn_out "      REAL(dp) :: temp\n";
print $eqn_out "      REAL(dp) :: press\n";
print $eqn_out "      REAL(dp) :: cair\n";
print $eqn_out "      REAL(dp) :: RO2\n";
print $eqn_out "      REAL(dp) :: $_\n" foreach (@mcm_vars);
print $eqn_out "\n";
print $eqn_out "      INTEGER, PARAMETER :: NJ_MOZART = $n_j_mozart\n";
print $eqn_out "      INTEGER, PARAMETER :: IP_MAX = $photo_max\n";
print $eqn_out "      INTEGER, PARAMETER :: $photo_labels[$_-1] = $_\n" for (1..$photo_max);
print $eqn_out "      REAL(dp), DIMENSION(IP_MAX) :: JX = 0.\n";
print $eqn_out "#ENDINLINE {above lines go into MODULE messy_mecca1_kpp_g_mem}\n";
print $eqn_out "\n";
print $eqn_out "#INLINE F95_RCONST\n";
print $eqn_out "      $_\n" foreach (@mcm_exps);
print $eqn_out "\n";
print $eqn_out "      RO2 = 0.\n";
print $eqn_out "      IF (KPP_$_ /= 0) RO2 = RO2 + C(KPP_$_)\n" for (@RO2_species);
print $eqn_out "#ENDINLINE {above lines go into the SUBROUTINE UPDATE_RCONST and UPDATE_PHOTO}\n";
print $eqn_out "\n";
print $eqn_out "#EQUATIONS\n";
for (1..$no_of_reactions) {
    my $label = '#R' . '0' x ($label_length - length $_) . $_; 
    print $eqn_out "{$label} $reactions[$_-1] ;\n";
}
close $eqn_out;

#write MECCA_MOZART module
open my $mecca_out, '>:encoding(utf-8)', $mecca_file or die $!;
print $mecca_out "!inorganic and methane chemistry\n";
print $mecca_out "MODULE MECCA_MOZART\n";
print $mecca_out "USE messy_mecca1_kpp_g_mem\n";
print $mecca_out "\n";
print $mecca_out "IMPLICIT NONE\n";
print $mecca_out "SAVE\n";
print $mecca_out "\n";
print $mecca_out "!Make sure these parameters are in  messy_mecca1_kpp_g_mem\n";
print $mecca_out "!REAL(dp), PRIVATE, DIMENSION(NJ_MOZART) :: J_MOZART = 0.\n";
print $mecca_out "!REAL(dp), PUBLIC, DIMENSION(IP_MAX) :: MOZART_PHOTO = 0.\n";
print $mecca_out "!REAL(dp), PUBLIC, DIMENSION(NSPEC) :: MOZART_PHOTO = 0.\n";
print $mecca_out "!REAL(dp), PUBLIC, DIMENSION(NSPEC) :: MOZART_VD = 0.\n";
print $mecca_out "!REAL(dp), PUBLIC, DIMENSION(NSPEC) :: MOZART_INI = 0.\n";
print $mecca_out "!REAL(dp), PUBLIC, DIMENSION(NSPEC) :: MOZART_EMIS = 0.\n";
print $mecca_out "!REAL(dp), PUBLIC, DIMENSION(NSPEC) :: MOZART_FIXED = 0.\n";
print $mecca_out "REAL(dp), PRIVATE :: INI_$_ = 0.\n" foreach (@species);
print $mecca_out "REAL(dp), PRIVATE :: EMIS_$_ = 0.\n" foreach (@species);
print $mecca_out "REAL(dp), PRIVATE :: FIXED_$_ = 0.\n" foreach (@species);
print $mecca_out "\n";
print $mecca_out "CONTAINS\n";
print $mecca_out "\n";
print $mecca_out "SUBROUTINE SET_MOZART_VD\n"; 
print $mecca_out "MOZART_VD(KPP_$_) = $deposition_velocities{$_}\n" foreach (keys %deposition_velocities);
print $mecca_out "END SUBROUTINE SET_MOZART_VD\n";
print $mecca_out "\n";
foreach my $thing ( qw(INI EMIS FIXED) ) {
    print $mecca_out "SUBROUTINE SET_MOZART_${thing}(iou)\n";
    print $mecca_out "INTEGER, INTENT(in) :: iou\n";
    print $mecca_out "LOGICAL :: lex = .false.\n";
    print $mecca_out "INTEGER :: fstat = 1\n";
    print $mecca_out "NAMELIST /${thing}/&\n";
    print $mecca_out "&${thing}_$_,&\n" for @species[0..$#species-1];
    print $mecca_out "&${thing}_$species[-1]\n";
    print $mecca_out "INQUIRE(file='MOZART_${thing}.nml', exist=lex)\n";
    print $mecca_out "IF (.NOT. lex) THEN\n";
    print $mecca_out "    WRITE (*,*) '*** WARNING: File MOZART_${thing}.nml NOT FOUND !'\n";
    print $mecca_out "    RETURN \n";
    print $mecca_out "ENDIF\n";
    print $mecca_out "OPEN(iou, file='MOZART_${thing}.nml')\n";
    print $mecca_out "WRITE (*,*) 'Reading namelist from MOZART_${thing}.nml'\n";
    print $mecca_out "READ(iou, NML=${thing}, IOSTAT=fstat)\n";
    print $mecca_out "IF (fstat /= 0) WRITE (*,*) '*** ERROR: reading MOZART_${thing}.nml!'\n";
    print $mecca_out "CLOSE(iou)\n";
    print $mecca_out "MOZART_${thing}(kpp_$_) = ${thing}_$_\n" foreach (@species);
    print $mecca_out "END SUBROUTINE SET_MOZART_${thing}\n";
    print $mecca_out "\n";
}

print $mecca_out "SUBROUTINE UPDATE_MOZART_PHOTO(COSX)\n";
print $mecca_out "REAL(dp) :: COSX\n";
print $mecca_out "J_MCM($_) = $jval_a[$_]*(COSX**$jval_b[$_])*EXP($jval_c[$_]/COSX)\n" foreach (@j_mcm);
print $mecca_out "MOZART_PHOTO($_) = $photo_rate[$_-1]\n" for (1..$photo_max);
print $mecca_out "END SUBROUTINE UPDATE_MOZART_PHOTO\n";
print $mecca_out "\n"; 
print $mecca_out "SUBROUTINE SET_MODEL_PARAMETERS(iou)\n";
print $mecca_out "INTEGER, INTENT(in) :: iou\n";
print $mecca_out "LOGICAL :: lex = .false.\n";
print $mecca_out "INTEGER :: fstat = 1\n";
print $mecca_out "NAMELIST /MODEL/&\n";
print $mecca_out "&end_model,&\n";
print $mecca_out "&scaling_factor,&\n";
print $mecca_out "&online_calc,&\n";
print $mecca_out "&temp,&\n";
print $mecca_out "&press,&\n";
print $mecca_out "&relhum,&\n";
print $mecca_out "&zmbl\n";
print $mecca_out "INQUIRE(file='MODEL_PARAMS.nml', exist=lex)\n";
print $mecca_out "IF (.NOT. lex) THEN\n";
print $mecca_out "    WRITE (*,*) '*** WARNING: File MODEL_PARAMS.nml NOT FOUND !'\n";
print $mecca_out "    RETURN \n";
print $mecca_out "ENDIF\n";
print $mecca_out "OPEN(iou, file='MODEL_PARAMS.nml')\n";
print $mecca_out "WRITE (*,*) 'Reading namelist from MODEL_PARAMS.nml'\n";
print $mecca_out "READ(iou, NML=MODEL, IOSTAT=fstat)\n";
print $mecca_out "IF (fstat /= 0) WRITE (*,*) '*** ERROR: reading MODEL_PARAMS.nml!'\n";
print $mecca_out "CLOSE(iou)\n";
print $mecca_out "WRITE (*,*) 'end_model = ', end_model, ' days'\n";
print $mecca_out "WRITE (*,*) 'scaling_factor = ', scaling_factor\n";
print $mecca_out "WRITE (*,*) 'temp = ', temp, ' K'\n";
print $mecca_out "WRITE (*,*) 'press = ', press, ' Pa'\n";
print $mecca_out "WRITE (*,*) 'relhum = ', relhum\n";
print $mecca_out "WRITE (*,*) 'zmbl = ', zmbl, ' m'\n";
print $mecca_out "WRITE (*,*) 'online_calc = ', online_calc\n";
print $mecca_out "END SUBROUTINE SET_MODEL_PARAMETERS\n";
print $mecca_out "\n";
print $mecca_out "END MODULE MECCA_MOZART\n";
close $mecca_out;
