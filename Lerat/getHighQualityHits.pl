#!/usr/bin/perl
# $Id: getHighQualityHits.pl 6 2013-11-16 01:42:31Z pjh $
#

# Phil Hatcher, May 2007
#
# Screen blastp data and keep the high-quality hits.
#
# There are two screening techniques:
#   1. E-value.
#   2. Lerat et al. technique that utilizes the ratio of the bitscore
#      to the maximal bit score.
#
# Input to this script are pair-wise BLAST results files. The filenames
# are <QueryGenome>-<TargetGenome>.blast. These files contain a line
# for each non-error gene in the query genome. The line starts with a
# query gene name followed by a list of hits. Each hit has the gene that
# was hit followed by bit score, e-value and alignment length, separated by
# exclamation points.
#
# All output is written to a single output file. This file contains a line
# for each gene, with the gene as the first thing on the line, followed by
# a space-separated list of genes that it hits above the threshold.
#
# Takes four initial command-line arguments:
#   1. directory that contains the BLAST results
#   2. -evalue or -lerat
#   3. either evalue threshold or the lerat ratio threshold
#   4. output file name
#
# These arguments are followed by a list of genome names that define the
# set of genomes being analyzed. This list must contain at least one genome.
#
# Instead of a list of genomes, -all can be specified. In this case the DONE
# file in the blast directory is consulted to get the list of genomes.
#

use strict;
use warnings;

my $startTime = time();

if (@ARGV < 5)
{
  die "Usage: getHighQualityHits.pl blastDirectory [-evalue | -lerat] " .
      "threshold outputFile <list of genomes>\n";
}

my $blastDirectory = shift @ARGV;
my $technique = shift @ARGV;
my $threshold = shift @ARGV;
my $outputFile = shift @ARGV;

print "Technique: $technique\n";
print "Threshold: $threshold\n";

if (($technique ne "-lerat") && ($technique ne "-evalue"))
{
  die "second argument must be either -lerat or -evalue";
}

my @genomes;
if ($ARGV[0] eq "-all")
{
  open(IN, "<", "$blastDirectory/DONE") or
    die "cannot open input ($blastDirectory/DONE)\n";
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

# open output file
open(OUT, ">", "$outputFile") or
  die "cannot open output ($outputFile)\n";

# process all pairs of genomes
foreach my $query (@genomes)
{
  print "Processing $query...\n";

  # for -lerat read self-hits into a hash
  my %selfHits = ();
  if ($technique eq "-lerat")
  {
    open(SELF, "<", "$blastDirectory/$query.self") or
      die "cannot open input ($blastDirectory/$query.self)\n";
    while (my $line = <SELF>)
    {
      chomp($line);
      my ($gene, $bitScore) = split / /, $line;
      if ($bitScore eq "")
      {
        die "null bit score for self-hit for $gene?\n";
      }
      if (defined($selfHits{$gene}))
      {
        die "$gene has more than one self-hit?\n";
      }
      $selfHits{$gene} = $bitScore;
    }
    close(SELF);
  }

  # keep a hash to collect together hits for the same query gene
  my %geneHash = ();

  foreach my $target (@genomes)
  {
    print "  Reading $query-$target.blast...\n";

    open(IN, "<", "$blastDirectory/$query-$target.blast") or
      die "cannot open input ($blastDirectory/$query-$target.blast)\n";

    my $line;
    while ($line = <IN>)
    {
      chomp($line);

      # line contains the query gene first followed by the genes it hit
      my @piece = split / /, $line;

      # grab the query gene
      my $key = shift @piece;
      my $keepers = "";

      # now screen for homologs (and remove self-hits too)
      for (my $i = 0; $i < @piece; $i += 1)
      {
        my @datum = split /!/, $piece[$i];
        if ($datum[0] eq $key)
        {
          next;
        }
        if ($technique eq "-lerat")
        {
          my $maxBitScore = $selfHits{$key};
          if (!defined($maxBitScore))
          {
            die "no self-hit for $key!\n";
          }
          if ($datum[1] >= ($threshold * $maxBitScore))
          {
            $keepers .= " $datum[0]";
          }
        }
        else
        {
          if ($datum[2] <= $threshold)
          {
            $keepers .= " $datum[0]";
          }
        }
      }
      if (defined($geneHash{$key}))
      {
        $geneHash{$key} = $geneHash{$key} . $keepers;
      }
      else
      {
        $geneHash{$key} = $keepers;
      }
    }
    close(IN);

    print "    Done.\n";
  }

  # dump the hash to the output file
  while (my ($key, $data) = each %geneHash)
  {
    print OUT "$key$data\n";
  }

  print "  Done.\n";
}
close(OUT);

my $stopTime = time();
my $executionTime = $stopTime - $startTime;
print "execution complete after $executionTime seconds.\n";
