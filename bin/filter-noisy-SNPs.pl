#!/usr/bin/env perl
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use v5.010_000;
use Getopt::Long;
use List::Util 'max';

=head1 NAME

filter-noisy-SNPs.pl - Noise-reduction for SNPs derived from a population

=head1 VERSION

Version 0.1.0

=cut

our $VERSION = '0.1.0';

=head1 SYNOPSIS

    perl filter-noisy-SNPs.pl genotyped/*.genotyped
        --cov_min           Minimum coverage to be considered [3]
        --homo_ratio_min    Minimum major allele ratio to be considered homozygous [0.9]
        --sample_ratio_min  Minimum ratio of homozygous samples to pass filter [0.9]
        --snp_dir           Directory containing [snp_master]
        --force             Overwrite previously generated noise-reduction output files

=head1 DESCRIPTION

When only one parent is available during SNP identification
based on a population (e.g., Recombinant Inbred Lines), it isn't
possible to do noise-reduction based on parental genotyped.
Therefore, I chose a different approach: identify and ignore
positions that have an over-representation of heterozygosity
in individual lines across the entire population.

=cut

my $cov_min          = 3;
my $homo_ratio_min   = 0.9;
my $sample_ratio_min = 0.9;
my $snp_dir          = "snp_master";
my $force;

my $options = GetOptions(
    "cov_min=i"          => \$cov_min,
    "homo_ratio_min=i"   => \$homo_ratio_min,
    "sample_ratio_min=i" => \$sample_ratio_min,
    "snp_dir=s"          => \$snp_dir,
    "force"              => \$force,
);

my $geno_file_list = \@ARGV;
my $homo_scores = score_snps( $geno_file_list, $cov_min, $homo_ratio_min );
filter_snps( $homo_scores, $sample_ratio_min );
write_snps( $homo_scores, $snp_dir, $force );

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

sub filter_snps {
    my ( $homo_scores, $sample_ratio_min ) = @_;

    for my $chr ( keys %{$homo_scores} ) {
        for my $pos ( keys %{ $$homo_scores{$chr} } ) {
            my $homo_count = $$homo_scores{$chr}{$pos}{'homo'} // 0;
            my $het_count  = $$homo_scores{$chr}{$pos}{'het'}  // 0;
            my $tot_count  = $homo_count + $het_count;
            my $sample_ratio = $tot_count > 0 ? $homo_count / $tot_count : 0;
            $$homo_scores{$chr}{$pos}{'keep'}
                = $sample_ratio >= $sample_ratio_min ? 1 : 0;
        }
    }
}

sub write_snps {
    my ( $homo_scores, $snp_dir, $force ) = @_;

    my @chr_list = keys %{$homo_scores};

    for my $chr (@chr_list) {
        my $snp_file = "$snp_dir/polyDB.$chr";
        my $snp_out_file = "$snp_file.nr";
        die "\nFile exists: $snp_out_file\nUse '--force' to overwrite.\n"
            if -e $snp_out_file && !$force;

        open my $snp_in_fh,  "<", $snp_file;
        open my $snp_out_fh, ">", $snp_out_file;

        my $header = <$snp_in_fh>;
        print $snp_out_fh $header;

        for my $snp_line (<$snp_in_fh>) {
            my ($pos) = ( split /\t/, $snp_line )[1];
            ($pos) = split /\./, $pos;    # for future support of inserts
            print $snp_out_fh $snp_line if $$homo_scores{$chr}{$pos}{'keep'};
        }

        close $snp_in_fh;
        close $snp_out_fh;
    }
}

=head1 AUTHOR

Michael F. Covington, <mfcovington@gmail.com>

=head1 SEE ALSO

L<https://github.com/mfcovington/snps-from-rils>

=head1 SOURCE AVAILABILITY

The source code is on Github:
L<https://github.com/mfcovington/noise-reduction-for-snps-from-pop>

=head1 BUGS

Please report any bugs or feature requests at
L<https://github.com/mfcovington/noise-reduction-for-snps-from-pop/issues>.

=head1 LICENSE

This is released under the Artistic
License. See L<perlartistic>.

=cut
