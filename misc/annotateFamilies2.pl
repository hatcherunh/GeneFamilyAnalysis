#!/usr/bin/perl

# $Id: annotateFamilies2.pl 36 2015-02-09 17:34:16Z pjh $
#
# Phil Hatcher, February 2015
#
# Quick hack of annotateFamilies.pl.
# Reads an ortholog family file and adds annotation for each family. The
# annotation is obtained from the FASTA header of the proteins file of
# a genome in the family. A CSV file is generated.
#
# This script takes four command-line arguments:
#   1. name of family file to be annotated
#   2. file containing the names of the genomes of interest
#   3. directory containing the proteins files.
#   4. name of the output file to be created
#
# The genome list should have the genome names simply separated by
# whitespace. These genome names must match the genome names in the
# family file, of course. The annotation will be taken from the first
# genome in this list that is actually represented in the family.
#
#

use strict;
use warnings;

# genome name is separated from gene name by this character
my $separatorCharacterPattern = "[\$]";
my $separatorCharacter = "\$";

# what is the separator in the family file?
my $familyFileSeparatorPattern = "[\t]";

if (@ARGV != 4)
{
  die "Usage: annotateFamilies2.pl familiesFile genomeListFile " .
    "proteinsDirectory outputFile\n";
}

my $familiesFile = $ARGV[0];
my $genomeListFile = $ARGV[1];
my $proteinsDirectory = $ARGV[2];
my $outputFile = $ARGV[3];

# first read genome list file

open(LIST, "<", $genomeListFile) or
  die "cannot open input ($genomeListFile)\n";

my @inputLines = <LIST>;

for (my $i = 0; $i < @inputLines; $i += 1) { chomp($inputLines[$i]); }

my $genomesStr = join " ", @inputLines;
my @genomes = split /\s/, $genomesStr;

my $genomeCount = @genomes;

close(LIST);

# now read the proteins files of all the genomes and build a hash to
# map gene name to annotation string
my %geneHash = ();
foreach my $genome (@genomes)
{
  open(PROTEINS, "<", "$proteinsDirectory/$genome.proteins") or
    die "cannot open input ($proteinsDirectory/$genome.proteins)\n";

  while (my $line = <PROTEINS>)
  {
    chomp($line);
    if ($line =~ /^>([^\s]*)\s(.*)$/)
    {
      my $key = "$genome$separatorCharacter$1";
      my $annot = $2;

      # remove any commas in the annotation
      $annot =~ s/,//g;

      my $old = $geneHash{$key};
      if (defined($old))
      {
        print STDERR "$key is duplicate in geneHash?\n";
      }
      else
      {
        $geneHash{$key} = $annot;
      }
    }
  }
  close(PROTEINS);
}

# now read the family file and output annotation for each family

open(FAMILIES, "<", $familiesFile) or
  die "cannot open input ($familiesFile)\n";

open(OUT, ">", $outputFile) or
  die "cannot open output ($outputFile)\n";

# print header line

print OUT "family,size,annotation";
foreach my $genome (@genomes)
{
  print OUT ",$genome";
}
print OUT "\n";

while (my $line = <FAMILIES>)
{
  chomp($line);

  my @pieces = split /$familyFileSeparatorPattern/, $line;

  # grab family number
  my $familyNumber = shift(@pieces);
  if ($familyNumber =~ /^(\d*):$/)
  {
    $familyNumber = $1;
  }
  else
  {
    die "can't parse family number ($familyNumber)\n";
  }
    
  my $familySize = @pieces;

  print OUT "$familyNumber,$familySize";

  # now process each member of the family
  my %famHash = ();
  foreach my $piece (@pieces)
  {
    my ($genome, $gene) = split /$separatorCharacterPattern/, $piece;
    my $old = $famHash{$genome};
    if (defined($old))
    {
      print STDERR "$genome is duplicate in famHash?\n";
    }
    else
    {
      $famHash{$genome} = $gene;
    }
  }

  # print the annotation for the first genome represented in the family
  print OUT ",";
  foreach my $genome (@genomes)
  {
    my $gene = $famHash{$genome};
    if (defined($gene))
    {
      my $annot = $geneHash{"$genome$separatorCharacter$gene"};
      if (defined($annot))
      {
        print OUT "$annot";
        last;
      }
    }
  }

  # now print the family members in the specified order
  foreach my $genome (@genomes)
  {
    my $gene = $famHash{$genome};
    print OUT ",";
    # allow ortholog families here (i.e. not all genomes represented)
    if (defined($gene))
    {
      print OUT $gene;
    }
  }
  print OUT "\n";
}

  
close(FAMILIES);
close(OUT);

