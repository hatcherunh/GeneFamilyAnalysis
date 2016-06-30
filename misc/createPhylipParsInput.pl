#!/usr/bin/perl

# $Id: createPhylipParsInput.pl 6 2013-11-16 01:42:31Z pjh $
#
# Phil Hatcher, May 2007
# modified May 2008 to take a fourth file containing the genome abbreviations
# modified May 2011 to make it clear what the input family file should contain
#
# Converts a ortholog family file produced by analyzeFamilies.pl to an
# input file for Phylip's PARS program. An ortholog family is a family
# with at most one gene from each genome.
#
# This script takes three command-line arguments:
#   1. name of ortholog family file to be processed
#   2. file containing the names of the genomes of interest
#   3. file containing genome name abbreviations to be used in the output
#   4. name of the output file to be created
#
# The genome list should have the genome names simply listed one per line.
# These genome names must match the genome names in the pan-ortholog family
# file, of course. The genome abbreviations should also be listed one per
# line and should be in the same order as the file of genome names. The
# genome abbreviations should be distinct within the first ten characters,
# given the constraints of the Phylip PARS input format. See below.
#
# The structure of the output file is:
#
# #genomes #families
# <abbrev1> 1 0 1 0 1 0...
# <abbrev2> 0 0 1 0 0 1...
# .
# .
# .
#
# where the <abbrev> field is 10 characters wide and there is 1 or 0 for
# whether the genome has a member in the family represented by that column.

use strict;
use warnings;

# genome name is separated from gene name by this character
my $separatorCharacter = "[\$]";
#my $separatorCharacter = "[_]";

if (@ARGV != 4)
{
  die "Usage: createPhylipParsInput.pl orthologFamiliesFile " .
        "genomeListFile genomeAbbreviationFile outputFile\n";
}

my $orthologFamiliesFile = $ARGV[0];
my $genomeListFile = $ARGV[1];
my $abbreviationFile = $ARGV[2];
my $outputFile = $ARGV[3];

# first read genome list and genome abbreviation files

open(LIST, "<", $genomeListFile) or
  die "cannot open input ($genomeListFile)\n";

open(ABBREV, "<", $abbreviationFile) or
  die "cannot open input ($abbreviationFile)\n";

my @genomes = <LIST>;
my @abbreviations = <ABBREV>;

if (@genomes != @abbreviations)
{
  die
    "genome list file and genome abbreviation file are not the same length\n";
}

for (my $i = 0; $i < @genomes; $i += 1) { chomp($genomes[$i]); }
for (my $i = 0; $i < @abbreviations; $i += 1) { chomp($abbreviations[$i]); }

my $genomeCount = @genomes;

close(LIST);
close(ABBREV);

# begin to format the output lines, one for each genome

my %outputLines = ();

for (my $i = 0; $i < @genomes; $i += 1)
{
  my $nameField = substr($abbreviations[$i], 0, 10);

  while (length($nameField) < 10)
  {
    $nameField .= " ";
  }
  $outputLines{$genomes[$i]} = $nameField;
}

# read the family file and build the output file in memory

my $familyCount = 0;

open(FAMILIES, "<", $orthologFamiliesFile) or
  die "cannot open input ($orthologFamiliesFile)\n";

while (my $line = <FAMILIES>)
{
  chomp($line);

  $familyCount += 1;

  print "processing family $familyCount...\n";

  my @pieces = split /\s/, $line;

  # discard family number
  shift(@pieces);

  # now process each member of the family
  foreach my $piece (@pieces)
  {
    my ($genome, $gene) = split /$separatorCharacter/, $piece;
    # print "  $genome\n";

    my $oldValue = $outputLines{$genome};
    if (!defined($oldValue))
    {
      #die "genome ($genome) does not have an output line?!\n";
      # skip genes from genomes not in the list
      # print "   no output line for $genome?\n";
      next;
    }

    # There may already have been a member of this genome in the family.
    # If so, this output line will already be of the proper length. In
    # that case just skip to the next member.
    if (length($oldValue) == (10 + (2 * $familyCount)))
    {
      # print "   skipping paralog?\n";
      next;
    }

    # If necessary, put 0's on output line for earlier families for which
    # this genome did not have a member. Current family should go in the
    # (10 + familyCount)th position (since the name field is of length 10).
    # If this genome's output line is less than this, then extend it with
    # 0's.
    while (length($oldValue) < (10 + 2*($familyCount - 1)))
    {
      $oldValue .= " 0";
    }

    # print "    marked!\n";
    # now tack on 1 for current family
    $oldValue .= " 1";

    # and put it back in the hash
    $outputLines{$genome} = $oldValue;
  }
}
close(FAMILIES);

# now print the output file

open(OUT, ">", $outputFile) or
  die "cannot open output ($outputFile)\n";

print OUT "$genomeCount $familyCount\n";

for (my $i = 0; $i < @genomes; $i += 1)
{
  my $value = $outputLines{$genomes[$i]};

  # Some lines may need 0's for the last families for which they did not have
  # members.
  while (length($value) < (10 + (2 * $familyCount)))
  {
    $value .= " 0";
  }
  print OUT "$value\n";
}
close(OUT);

