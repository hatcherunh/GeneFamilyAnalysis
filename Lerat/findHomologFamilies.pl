#!/usr/bin/perl
# $Id: findHomologFamilies.pl 6 2013-11-16 01:42:31Z pjh $
#

#
# Phil Hatcher, July 2013
# modified to add an additional output file
#
# The additional output file is a comma-separated file that contains
# a line for each family, with a count for each genome of the number
# of genes from that genome in the family.
#
# The name of the output file is simply the user-given output file
# name (third command-line argument) with ".csv" appended.
#
# Phil Hatcher, May 2007
#
# Using processed BLAST data for a set of genomes, compute homolog
# families for the genes in those genomes.
#
# A command-line switch allows either one-way or reciprocal hits to
# be used to determine family membership.
#
# Primary input to this script is a high-quality hits file. And if
# the reciprocal option is being used, then a high-quality reverse hits
# file must also be input.
#
# In the forward hits file, each line starts with a query gene name followed
# by the list of subject genes that were hit.
#
# In the reverse hits file, each line starts with a subject gene name followed
# by the list of query genes that hit the subject gene.
#
# Writes one output file that contains the computed homolog families, one per
# line, where a pair of genes is added to a family if one is a high-quality hit
# of the other. (If reciprocal option, then they must be reciprocal
# high-quality hits.)
#
# Takes three or four command-line arguments:
#   1. input file for high-quality hits
#   2. -reciprocal or -oneway
#   3. output file
#   4. if -reciprocal then this is the input file for high-quality reverse hits
#

use strict;
use warnings;

my $startTime = time();

if (@ARGV < 3 || @ARGV > 4 ||
    ($ARGV[1] ne "-reciprocal" && $ARGV[1] ne "-oneway") ||
    (@ARGV == 4 && $ARGV[1] ne "-reciprocal"))
{
  die "Usage: findHomologFamilies.pl hitsInput -oneway output\n" .
      "       findHomologFamilies.pl hitsInput -reciprocal output " .
        "reverseHitsInput\n";
}

my $hitsFile = $ARGV[0];
my $mode = $ARGV[1];
my $outputFile = $ARGV[2];
my $outputFile2 = $ARGV[2] . ".csv";;
my $reverseHitsFile = "";

if ($mode eq "-reciprocal")
{
  $reverseHitsFile = $ARGV[3];
}

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

# list of families
my $familyCnt = 0;
my @family = ();

# hash to map genes to families that contain them
my %familyHash = ();

# hash to track genomes present in families
my %genomeHash = ();

# makes a new family
#   only is called with two genes
sub newFamily
{
  # takes space separated "list" of genes
  my $gene1 = $_[0];
  my $gene2 = $_[1];

  my $newFamilyCnt = $familyCnt;
  $familyCnt += 1;
  $family[$newFamilyCnt] = "$gene1 $gene2";

  $familyHash{$gene1} = $newFamilyCnt;
  $familyHash{$gene2} = $newFamilyCnt;
}

# add two genes to a family
#   genes might already be in separate families
#     if so families must be disjoint
#   genes might already be in the same family too!
sub updateFamilies
{
  my $gene1 = $_[0];
  my $gene2 = $_[1];

  my $family1 = $familyHash{$gene1};
  my $family2 = $familyHash{$gene2};

  if ((!defined($family1)) && (!defined($family2)))
  {
    newFamily($gene1, $gene2);
  }
  elsif (defined($family1) && (!defined($family2)))
  {
    my $new = $family[$family1] . " " . $gene2;
    $family[$family1] = $new;
    $familyHash{$gene2} = $family1;
  }
  elsif (defined($family2) && (!defined($family1)))
  {
    my $new = $family[$family2] . " " . $gene1;
    $family[$family2] = $new;
    $familyHash{$gene1} = $family2;
  }
  else
  {
    # if already in the same family then nothing to do
    if ($family1 != $family2)
    {
      # if in different families then they are disjoint
      my $new = $family[$family1] . " " . $family[$family2];
      $family[$family1] = $new;
      foreach my $g (split / /, $family[$family2])
      {
        $familyHash{$g} = $family1;
      }
      $family[$family2] = "";
    }
  }
}

open(HITS, "<", $hitsFile) or
  die "cannot open hits file for input ($hitsFile)\n";

if ($mode eq "-reciprocal")
{
  open(REVERSE_HITS, "<", $reverseHitsFile) or
    die "cannot open reverse hits file for input ($reverseHitsFile)\n";
}

my $lineCount = 0;
while (my $line = <HITS>)
{
  chomp($line);

  $lineCount += 1;

  my @pieces = split / /, $line;

  my $gene = shift @pieces; 

  if (($lineCount % 1000) == 0)
  {
    print "processing $gene ($lineCount):\n";
  }

  # remember the genome
  my ($genomeName, $geneName) = split /\$/, $gene;
  $genomeHash{$genomeName} = $genomeName;

  my $hitsStr = join " ", @pieces;

  if ($mode eq "-reciprocal")
  {
    my $line2 = <REVERSE_HITS>;

    if (!defined($line2))
    {
      die "unexpected EOF on reverse hits file at line $lineCount!\n";
    }

    chomp($line2);

    my @pieces2 = split / /, $line2;
 
    my $gene2 = shift @pieces2;

    if ($gene ne $gene2)
    {
      die
      "gene mismatch ($gene, $gene2) in two input files at line $lineCount!\n";
    } 

    my $hitsStr2 = join " ", @pieces2;

    my @reciprocals = intersect($hitsStr, $hitsStr2);

    # now use each reciprocal hit to update the families
    foreach my $reciprocalHit (@reciprocals)
    {
      updateFamilies($gene, $reciprocalHit);
    }
  }
  elsif ($mode eq "-oneway")
  {
    # just use the one-way hits to update the families
    foreach my $hit (split / /, $hitsStr)
    {
      updateFamilies($gene, $hit);
    }
  }
  else
  {
    die "bogus mode ($mode)!\n";
  }
}
close(HITS);
if ($mode eq "-reciprocal")
{
  close(REVERSE_HITS);
}

print "families complete.\n";

# now dump out the families
open(FAMILY, ">", $outputFile) or
  die "cannot open output ($outputFile)\n";
open(CSV, ">", $outputFile2) or
  die "cannot open output ($outputFile2)\n";
print CSV "family";
foreach my $genomeName ( keys(%genomeHash) ) {
  print CSV ",$genomeName";
}
print CSV "\n";
for (my $i = 0; $i < $familyCnt; $i += 1)
{
  my %countHash = ();
  if ($family[$i] ne "")
  {
    print FAMILY "$i:";
    print CSV "$i";
    my @genes = split / /, $family[$i];
    foreach my $gene (@genes)
    {
      my ($genomeName,$geneName) = split /\$/, $gene;
      my $n = $countHash{$genomeName};
      if (defined($n))
      {
        $countHash{$genomeName} = $n + 1;
      }
      else
      {
        $countHash{$genomeName} = 1;
      }
      print FAMILY " $gene";
    }
    print FAMILY "\n";
    foreach my $genomeName ( keys(%genomeHash) ) {
      my $cnt = $countHash{$genomeName};
      if (!defined($cnt))
      {
        $cnt = 0;
      }
      print CSV ",$cnt";
    }
    print CSV "\n";
  }
}
close(FAMILY);
close(CSV);

print "families dumped to file.\n";

my $stopTime = time();
my $executionTime = $stopTime - $startTime;
print "execution complete after $executionTime seconds.\n";
