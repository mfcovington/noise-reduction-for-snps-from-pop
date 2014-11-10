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

