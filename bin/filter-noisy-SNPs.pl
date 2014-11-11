#!/usr/bin/env perl
# Mike Covington
# created: 2014-11-10
#
# Description: When only one parent is available during SNP identification
#              based on a population (e.g., Recombinant Inbred Lines), it isn't
#              possible to do noise-reduction based on parental genotyped.
#              Therefore, I chose a different approach: identify and ignore
#              positions that have an over-representation of heterozygosity
#              in individual lines across the entire population.
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use List::Util 'max';

my $cov_min          = 3;
my $homo_ratio_min   = 0.9;

my $options = GetOptions(
    "cov_min=i"          => \$cov_min,
    "homo_ratio_min=i"   => \$homo_ratio_min,
);

my $geno_file_list = \@ARGV;
my $homo_scores = score_snps( $geno_file_list, $cov_min, $homo_ratio_min );

exit;

sub score_snps {
    my ( $geno_file_list, $cov_min, $homo_ratio_min ) = @_;
    my $homo_scores;

    for my $geno_file (@$geno_file_list) {
        open my $geno_fh, "<", $geno_file;
        for (<$geno_fh>) {
            chomp;
            my ( $chr, $pos, $par1, $par2, $tot ) = split;
            my $score = homo_or_het( $par1, $par2, $tot, $cov_min,
                $homo_ratio_min );
            $$homo_scores{$chr}{$pos}{$score}++;
        }
        close $geno_fh;
    }

    return $homo_scores;
}

sub homo_or_het {
    my ( $par1, $par2, $tot, $cov_min, $homo_ratio_min ) = @_;

    my $score;

    if ( $tot == 0 || $tot < $cov_min ) {
        $score = 'NA';
    }
    else {
        my $homo_ratio = max( $par1, $par2 ) / $tot;
        $score = $homo_ratio >= $homo_ratio_min ? 'homo' : 'het';
    }

    return $score;
}

