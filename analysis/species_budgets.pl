#! /usr/bin/env perl
# Compare budgets in old and new tagging of species in @ARGV
# Version 0: Jane Coates 26/5/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use Statistics::R;

my $species = $ARGV[0];
die "Specify species\n" unless (defined $species);

my $base = "/local/home/coates";
my @runs = qw( MECCA New_tagging );
my %data;

my $mecca = MECCA->new("$base/$runs[0]/MOZART-4_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;
my $n_days = int $ntime / $n_per_day;

foreach my $run (@runs) {
    my $dir = "$base/$run";
    if ($run =~ /tagging/) {
        $dir .= "/MOZART-4_VOC_tagged";
    } else {
        $dir .= "/MOZART-4_tagged";
    }
    my $mecca = MECCA->new("$dir/boxmodel");
    my $kpp = KPP->new("$dir/gas.eqn");
    my $producers = $kpp->producing($species);
    my $producer_yields = $kpp->effect_on($species, $producers);
    my $consumers = $kpp->consuming($species);
    my $consumer_yields = $kpp->effect_on($species, $consumers);
    print "No consumers for $species\n" if (@$consumers == 0);
    print "No producers for $species\n" if (@$producers == 0);
    
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        $data{$run}{$reactants} += $rate;
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my $reaction_string = $kpp->reaction_string($reaction);
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        $data{$run}{$reactants} += $rate;
    }
}

foreach my $run (keys %data) {
    foreach my $reaction (keys %{$data{$run}}) {
        my $reshape = $data{$run}{$reaction}->reshape($n_per_day, $n_days);
        $data{$run}{$reaction} = $reshape->sumover;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `, 
);
$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` d = data.frame() `);

foreach my $run (sort keys %data) {
    $R->set('run', $run);
    $R->run(q` pre = data.frame(Time, Run = rep(run, length(Time))) `);
    foreach my $reaction (sort keys %{$data{$run}}) {
        $R->set('reaction', $reaction);
        $R->set('rate', [ map { $_ } $data{$run}{$reaction}->dog ]);
        $R->run(q` pre[reaction] = rate `);
    }
    $R->run(q` pre = gather(pre, Reaction, Rate, -Time, -Run) `,
            q` d = rbind(d, pre) `,
    );
}

$R->set('filename', "${species}_budget_by_reactions.pdf");
$R->run(q` p = ggplot(d, aes(x = Time, y = Rate, fill = Reaction, order = Reaction)) `,
        q` p = p + geom_bar(stat = "identity", data = subset(d, Rate < 0)) `,
        q` p = p + geom_bar(stat = "identity", data = subset(d, Rate > 0)) `,
        q` p = p + facet_wrap( ~ Run) `,
        q` p = p + theme(legend.title = element_blank()) `,
);

$R->run(q` CairoPDF(file = filename, width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();
