#!/usr/bin/perl

# $Id: findMaximalPanorthologFamilies.pl 6 2013-11-16 01:42:31Z pjh $
#
# Phil Hatcher, June 2007
# updated Oct 2011 to allow paralog families to be in the input. they
#   are just skipped.
#
# Reads a family file and finds the families that have exactly one member
# from each genome in a given list of genomes
#
# This script takes two initial command-line arguments:
#   1. name of family file to be processed
#   2. output file
#
# These arguments are followed by a list of genome names that define the
# set of genomes to be used. This list must contain at least one genome.
#
# Instead of a list of genomes, -all can be specified. In this case the DONE
# file in the blast directory is consulted to get the list of genomes. The
# BLAST directory must be given as an extra argument after the -all in this
# case.


use strict;
use warnings;

my $startTime = time();

# genome name is separated from gene name by this character
my $separatorCharacter = "[\$]";
#my $separatorCharacter = "[_]";

if (@ARGV < 2)
{
  die "Usage: findMaximalPanorthologFamilies.pl panorthologFile outputFile " .
        "<list-of-genomes>\n";
}

my $inputFile = shift @ARGV;
my $outputFile = shift @ARGV;

my @genomes;
if ($ARGV[0] eq "-all")
{
  if (@ARGV != 2)
  {
    die "-all must be followed by (only) the blast directory!\n";
  }
  open(IN, "<", "$ARGV[1]/DONE") or
    die "cannot open input ($ARGV[1]/DONE)\n";
  while (my $line = <IN>)
  {
    chomp($line);
    push @genomes, $line;
  }
  close (IN);
}
else
{
  @genomes = @ARGV;
}
my $numberOfGenomes = @genomes;
print "number of genomes is $numberOfGenomes\n";

# put genomes in a hash
my %genomeHash = ();
foreach my $genome (@genomes)
{
  if (defined($genomeHash{$genome}))
  {
    die "duplicate genome ($genome) in the list of genomes!\n";
  }
  else
  {
    $genomeHash{$genome} = 0;
  }
}

# open output file
open(OUT, ">", "$outputFile") or
  die "cannot open output ($outputFile)\n";

# read families and output families that have exactly one gene from
# each genome in the list of genomes.

open(IN, "<", $inputFile) or
  die "cannot open input ($inputFile)\n";

my $maxCount = 0;

while (my $line = <IN>)
{
  chomp($line);

  my @pieces = split /\s/, $line;

  # get family number
  my $family = shift(@pieces);
  $family =~ s/:$//;

  # count how many genomes in the list are represented
  my $count = 0;
  foreach my $gene (@pieces)
  {
    (my $genome, my $discard) = split /$separatorCharacter/, $gene;

    #print "DBG: ($genome, $discard)\n";
    if (defined($genomeHash{$genome}))
    {
      # paralog found?
      if ($genomeHash{$genome} == 1)
      {
        $count = -1;
        next;
      }
      else
      {
        $genomeHash{$genome} = 1;
        $count += 1;
      }
    }
    else
    {
      die "$genome not found in genome hash?\n"
    }
  }
  #print "family $family: $count of $numberOfGenomes genomes represented\n";

  # if all genomes are represented, then output the family
  if ($count == $numberOfGenomes)
  {
    print OUT "$line\n";
    $maxCount += 1;
  }
  # re-intialize hash for next iteration of look
  foreach my $genome (@genomes)
  {
    $genomeHash{$genome} = 0;
  }
}
close(IN);
close(OUT);

print "$maxCount maximal pan-ortholog families found!\n";

my $stopTime = time();
my $executionTime = $stopTime - $startTime;
print "execution complete after $executionTime seconds.\n";

