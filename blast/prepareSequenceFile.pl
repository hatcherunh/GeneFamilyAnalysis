#!/usr/bin/perl

# $Id: prepareSequenceFile.pl 35 2015-02-03 04:08:03Z pjh $
#
# Phil Hatcher, May 2007
#
# Written to process Nematode protein files to prepare to be processed
# by my findHomology pipeline. Only the FASTA headers are altered. The
# genome name is prepended to the gene name, and any annotation in the
# header is discarded.
#
# This script takes three arguments:
#   1. The genome name.
#   2. The input proteins file.
#   3. The output proteins file.
#
# pjh Jan 2015: Modified to generalize to support both nucleotide and
#               protein data. Actually just needed to change the
#               script name.

use strict;
use warnings;

my $separatorCharacter = "\$";

if (@ARGV != 3)
{
  die "Usage: prepareSequenceFile.pl genome inputFile outputFile\n";
}

my $genomeName = $ARGV[0];

# open input file
my $inputFile = $ARGV[1];
open(IN, "<", $inputFile) or
  die "cannot open input ($inputFile)\n";

# open output file
my $outputFile = $ARGV[2];
open(OUT, ">", $outputFile) or
  die "cannot open output ($outputFile)\n";

# read the input file
while (my $line = <IN>)
{
  chomp($line);

  # look for FASTA header
  if ($line =~ /^>/)
  {
    if ($line =~ /^>\s*([^\s]*)/)
    {
      my $geneName = $1;

      print OUT ">$genomeName$separatorCharacter$geneName\n";
    }
    else
    {
      die "can't parse FASTA header to get gene name!\n" .
          "  $line\n";
    }
  }
  else
  {
    print OUT "$line\n";
  }
}

