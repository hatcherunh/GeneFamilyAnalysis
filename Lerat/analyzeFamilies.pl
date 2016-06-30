#!/usr/bin/perl

# $Id: analyzeFamilies.pl 6 2013-11-16 01:42:31Z pjh $
#
# Phil Hatcher, April 2007
# modified in July 2008 to change terminology to use ortholog to mean
#   family has at most one gene from each genome.
#
# Analyzes the families computed by findHomologFamilies.pl.
#
# This scripts takes two command-line arguments:
#   1. name of family file to be processed
#   2. prefix to be used to build filenames for output files
#
# For each family it computes the number of genes from each genome represented.
# This and some other statistical data (e.g. size of largest family, etc.) is
# printed to <outputPrefix>.stats.
#
# In addition, it writes out families that contain at most one gene per genome.
# This output is written out to <outputPrefix>.orthologs. (Remember that
# a family contains at least two genes. Unique genes, which only appear in one
# genome, are reported by findUniques.pl.)
#
# Families with at least one represented genome with more than one gene will
# be written out to <outputPrefix>.paralogs. In addition, for each of these
# families a summary line is written containing the count of genes per genome.
# These summary lines are <outputPrefix>.paralogs-summary.
#

use strict;
use warnings;

my $startTime = time();

# genome name is separated from gene name by this character
my $separatorCharacter = "[\$]";
#my $separatorCharacter = "[_]";

if (@ARGV != 2)
{
  die "Usage: analyzeFamilies.pl familyFile outputPrefix\n";
}

my $inputFile = $ARGV[0];
my $outputPrefix = $ARGV[1];

# open file for printing statistics
# (June 2010: it should already exist so append to it)

open(STATS, ">>", $outputPrefix . ".stats") or
  die "cannot open output ($outputPrefix.stats)\n";

# read file once just to determine number of genomes present, largest
# family, average family size.

open(IN, "<", $inputFile) or
  die "cannot open input ($inputFile)\n";

my %genomeHash = ();

my $familyCount = 0;
my $orthologCount = 0;
my $panorthologCount = 0;

my $geneCount = 0;

my $largestFamily;
my $largestFamilySize = -1;

while (my $line = <IN>)
{
  chomp($line);

  $familyCount += 1;

  my @pieces = split / /, $line;

  # get family number
  my $family = shift(@pieces);
  $family =~ s/:$//;

  $geneCount += @pieces;

  if (@pieces > $largestFamilySize)
  {
    $largestFamily = $family;
    $largestFamilySize = @pieces;
  }

  foreach my $gene (@pieces)
  {
    (my $genome, my $discard) = split /$separatorCharacter/, $gene;

    $genomeHash{$genome} += 1;
  }
}
close(IN);

print STATS "families: $familyCount\n";
print STATS "genes in families: $geneCount\n";
print STATS "largest family size: $largestFamilySize\n";
print STATS "largest family number: $largestFamily\n";

my $avg = $geneCount / $familyCount;
print STATS "average family size: $avg\n";

print STATS "count of genes in families for each genome:\n";
my $genomeCount = 0;
while ((my $key, my $value) = each %genomeHash)
{
  $genomeCount += 1;
  print STATS "  $key $value\n";
}

# read families again and output for each family the number of genes present
# for each genome. if a family has at most one gene per genome, then output
# the family to the ortholog output file, otherwise output it to the
# paralog output file. And output a summary line for each paralog family
# to the paralog summary file. Also compute a basic histogram of family
# size. The histogram bins are: <= 2 * number of genomes
#                               <= 5 * number of genomes
#                               <= 10 * number of genomes
#                               <= 50 * number of genomes
#                                > 50 * number of genomes
# This histogram is written to the statistics file.

open(IN, "<", $inputFile) or
  die "cannot open input ($inputFile)\n";

open(ORTH, ">", $outputPrefix . ".orthologs") or
  die "cannot open output ($outputPrefix.orthologs)\n";

open(PAR, ">", $outputPrefix . ".paralog") or
  die "cannot open output ($outputPrefix.paralog)\n";

open(PARSUM, ">", $outputPrefix . ".paralog-summary") or
  die "cannot open output ($outputPrefix.paralog-summary)\n";

# used for histogram of family sizes
my $bin1 = 0;
my $bin2 = 0;
my $bin3 = 0;
my $bin4 = 0;
my $bin5 = 0;

while (my $line = <IN>)
{
  chomp($line);

  my @pieces = split / /, $line;

  # get family number
  my $family = shift(@pieces);
  $family =~ s/:$//;

  # compute basic histogram of family size
  if (@pieces <= 2 * $genomeCount)
  {
    $bin1 += 1;
  }
  elsif (@pieces <= 5 * $genomeCount)
  {
    $bin2 += 1;
  }
  elsif (@pieces <= 10 * $genomeCount)
  {
    $bin3 += 1;
  }
  elsif (@pieces <= 20 * $genomeCount)
  {
    $bin4 += 1;
  }
  else
  {
    $bin5 += 1;
  }

  %genomeHash = ();

  # count genes per genome
  foreach my $gene (@pieces)
  {
    (my $genome, my $discard) = split /$separatorCharacter/, $gene;

    $genomeHash{$genome} += 1;
  }

  my $cnt = 0;
  my $allOne = 1;
  my @lst = ();
  while ((my $key, my $value) = each %genomeHash)
  {
    if ($value != 1)
    {
      #print "DBG: family $family ($key, $value)\n";
      $allOne = 0;
    }
    push @lst, "$key $value";
    $cnt += 1;
  }

  if ($allOne)
  {
    # change from space-separated to tab-separated
    $line =~ s/ /\t/g;

    print ORTH "$line\n";

    $orthologCount += 1;
    if ($cnt == $genomeCount)
    {
      $panorthologCount += 1;
    }
  }
  else
  {
    print PAR "$line\n";

    print PARSUM "$family:";
    foreach my $x (sort @lst)
    {
      print PARSUM " $x";
    }
    print PARSUM "\n";
  }
}
close(IN);
close(ORTH);
close(PAR);
close(PARSUM);

print STATS "ortholog families: $orthologCount\n";
print STATS "panortholog families: $panorthologCount\n";

print STATS "family size histogram:\n";
my $lowBound = 2;
my $highBound = 2 * $genomeCount;
print STATS "  $lowBound to $highBound: $bin1\n";
$lowBound = $highBound + 1;
$highBound = 5 * $genomeCount;
print STATS "  $lowBound to $highBound: $bin2\n";
$lowBound = $highBound + 1;
$highBound = 10 * $genomeCount;
print STATS "  $lowBound to $highBound: $bin3\n";
$lowBound = $highBound + 1;
$highBound = 20 * $genomeCount;
print STATS "  $lowBound to $highBound: $bin4\n";
print STATS "  >$highBound: $bin5\n";
close(STATS);

my $stopTime = time();
my $executionTime = $stopTime - $startTime;
print "execution complete after $executionTime seconds.\n";

