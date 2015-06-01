#! /usr/bin/env perl

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use Statistics::R;

my $base = "/local/home/coates";
my @runs = qw( MECCA New_tagging );
my %data;

my $mecca = MECCA->new("$base/$runs[0]/MOZART-4_tagged/boxmodel");
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;

foreach my $run (@runs) {
    my $dir = "$base/$run";
    if ($run =~ /tagging/) {
        $dir .= "/MOZART-4_VOC_tagged";
    } else {
        $dir .= "/MOZART-4_tagged";
    }
    my $mecca = MECCA->new("$dir/boxmodel");
    my $kpp = KPP->new("$dir/gas.eqn");

    my $reactions = $kpp->reacting_with("O3", "hv");
    foreach my $reaction (@$reactions) {
        my $reaction_string = $kpp->reaction_string($reaction);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        $data{$run}{$reaction_string} += $rate;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `, 
);
$R->set('Time', [ map { $_ } $time->dog ]);
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

$R->run(q` p = ggplot(d, aes(x = Time, y = Rate, colour = Reaction, order = Reaction)) `,
        q` p = p + geom_line() `,
        q` p = p + facet_wrap( ~ Run) `,
        q` p = p + theme(legend.title = element_blank()) `,
);

$R->run(q` CairoPDF(file = "O3_photolysis_rates_comparison.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();
