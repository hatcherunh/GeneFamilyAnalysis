#!/usr/bin/perl

# $Id: findUniques.pl 6 2013-11-16 01:42:31Z pjh $
#
# Phil Hatcher, May 2007
#
# Using processed BLAST data for a set of genomes, find unique genes.
#
# A command-line switch allows either one-way or reciprocal hits to
# be used to determine uniqueness. That is, if the -oneway option is
# selected, then a one-way hit either to A or from A will eliminate
# A as a unique. If the -reciprocal option is selected, then a gene
# must have a reciprocal hit with another gene in order to be eliminated
# as a unique gene.
#
# Input to this script is a high-quality hits file and a high-quality
# reverse hits file.
#
# In the forward hits file, each line starts with a query gene name followed
# by a list of subject genes that were hit.
#
# In the reverse hits file, each line starts with a subject gene name followed
# by a list of query genes that hit the subject gene.
#
# Writes one output file that contains the unique genes, one per line.
#
# Takes four command-line arguments:
#   1. prefix for input and output files
#   2. -reciprocal or -oneway
#
# May 2011 - edits made to ensure each genome will be listed in the
#            list of genomes with unique gene counts. Genomes with
#            no unique genes will print a count of zero.

use strict;
use warnings;

# genome name is separated from gene name by this character
my $separatorCharacter = "[\$]";
#my $separatorCharacter = "[_]";

my $startTime = time();

if (@ARGV != 2)
{
  die "Usage: findUniques.pl filePrefix [-oneway | -reciprocal]\n";
}

my $prefix = $ARGV[0];
my $mode = $ARGV[1];
my $reverseHitsFile = "$prefix.reverse";
my $hitsFile = "$prefix.hits";
my $outputFile = "$prefix.unique";

# takes two space separated "lists" of genes and computes the intersection
# of the two lists. input could be either undef or empty string indicating
# the empty set.
sub intersect
{
  my $arg1 = $_[0];
  my $arg2 = $_[1];
  my @ret;

  if (!defined($arg1) || !defined($arg2))
  {
    @ret = ();
  }
  elsif ($arg1 eq "" || $arg2 eq "")
  {
    @ret = ();
  }
  else
  {
    my @list1 = split / /, $arg1;
    my @list2 = split / /, $arg2;
  
    my %tmp = ();
    foreach my $gene (@list1)
    {
      $tmp{$gene} = 0;
    }
    foreach my $gene (@list2)
    {
      if (defined($tmp{$gene}))
      {
        $tmp{$gene} = 1;
      }
    }
  
    @ret = ();
    while ((my $key, my $value) = each %tmp)
    {
      if ($value == 1)
      {
        push(@ret, $key);
      }
    }
  }
  return @ret;
}

open(HITS, "<", $hitsFile) or
  die "cannot open hits file for input ($hitsFile)\n";

open(REVERSE_HITS, "<", $reverseHitsFile) or
  die "cannot open reverse hits file for input ($reverseHitsFile)\n";

open(OUT, ">", $outputFile) or
  die "cannot open output file ($outputFile)\n";

my %genomeHash = ();

my $lineCount = 0;

while (my $line = <HITS>)
{
  chomp($line);

  $lineCount += 1;

  my @pieces = split / /, $line;

  my $gene = shift @pieces; 
  (my $genome, my $discard) = split /$separatorCharacter/, $gene;
  if (!exists $genomeHash{$genome})
  {
    $genomeHash{$genome} = 0;
  }

  my $hitsStr = join " ", @pieces;

  my $line2 = <REVERSE_HITS>;

  chomp($line2);

  if (!defined($line2))
  {
    die "unexpected EOF on reverse hits file at line $lineCount!\n";
  }

  my @pieces2 = split / /, $line2;
 
  my $gene2 = shift @pieces2;

  if ($gene ne $gene2)
  {
    die
     "gene mismatch ($gene, $gene2) in two input files at line $lineCount!\n";
  } 

  my $hitsStr2 = join " ", @pieces2;

  if ($mode eq "-oneway")
  {
    if ((@pieces == 0) && (@pieces2 == 0))
    {
      $genomeHash{$genome} += 1;

      print OUT "$gene\n";
    }
  }
  elsif ($mode eq "-reciprocal")
  {
    my @reciprocals = intersect($hitsStr, $hitsStr2);
    if (@reciprocals == 0)
    {
      $genomeHash{$genome} += 1;

      print OUT "$gene\n";
    }
  }
  else
  {
    die "bogus mode ($mode)!\n";
  }
}
close(HITS);
close(REVERSE_HITS);
close(OUT);

# June 2010: append unique counts to the *.stats file
open(STATS, ">>", $prefix . ".stats") or
  die "cannot open output ($prefix.stats)\n";
print STATS "count of unique genes for each genome:\n";
while ((my $key, my $value) = each %genomeHash)
{
  print STATS "  $key $value\n";
}
close(STATS);


my $stopTime = time();
my $executionTime = $stopTime - $startTime;
print "execution complete after $executionTime seconds.\n";
