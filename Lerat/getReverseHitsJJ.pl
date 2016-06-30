#!/usr/bin/perl

# $Id: getReverseHitsJJ.pl 31 2014-09-04 02:17:08Z pjh $
#
# Phil Hatcher, August 2014
#
# Takes a high-quality hits file and produces the corresponding reverse
# hits file.
#
# This script takes two arguments:
#   1. The input high-quality hits file.
#   2. The output high-quality reverse hits file.
#
# This version uses the notion of buckets from Jamie Jackson's MS thesis.

use strict;
use warnings;

# for safely creating temp files for the buckets
use File::Temp qw/ tempfile tempdir /;

# in units of numbers of genes
# i.e. divide gene number by this to get bucket number
my $BUCKET_SIZE = 500;

# maximum number of buckets
# this is driven by the limit on the number of files that can be open
#   at the same time.
my $MAX_NUMBER_OF_BUCKETS = 1000;

# temp file name prefix
my $TEMP_FILE_PREFIX = "bucket";

my $startTime = time();

if (@ARGV != 2)
{
  die "Usage: getReverseHitsJJ.pl inputFile outputFile\n";
}

# open input file
my $inputFile = $ARGV[0];
open(IN, "<", $inputFile) or
  die "cannot open input ($inputFile)\n";

# read the input file a first time to number the genes in order
my %geneToNumber = ();
my %numberToGene = ();
my $counter = 0;
while (my $line = <IN>)
{
  chomp($line);

  my @pieces = split / /, $line;

  my $gene = shift @pieces;

  $geneToNumber{$gene} = $counter;
  $numberToGene{$counter} = $gene;

  $counter += 1;
}
close(IN);
my $totalGeneCount = $counter;

# create the buckets
# i.e. create and open a temp file for each bucket

my $numberOfBuckets;
{ # limit scope of "use integer"
  use integer;
  $numberOfBuckets = $totalGeneCount / $BUCKET_SIZE;
  if ($totalGeneCount % $BUCKET_SIZE != 0)
  {
    $numberOfBuckets += 1;
  }
}
if ($numberOfBuckets > $MAX_NUMBER_OF_BUCKETS)
{
  $numberOfBuckets = $MAX_NUMBER_OF_BUCKETS;
  { # limit scope of "use integer"
    use integer;
    $BUCKET_SIZE = $totalGeneCount / $numberOfBuckets;
    if ($totalGeneCount % $numberOfBuckets != 0)
    {
      $BUCKET_SIZE += 1;
    }
  }
}
print "Total gene count = $totalGeneCount\n";
print "Number of buckets: $numberOfBuckets\n";
print "Bucket size: $BUCKET_SIZE\n";

my $i = 0;
my @buckets = ();
my @filenames = ();
while ($i < $numberOfBuckets)
{
  ($buckets[$i], $filenames[$i]) = tempfile($TEMP_FILE_PREFIX."XXXXXXXX");
  $i += 1;
}
print "Buckets created.\n";

# re-open input file
open(IN, "<", $inputFile) or
  die "cannot re-open input ($inputFile)\n";

# read the input file again
while (my $line = <IN>)
{
  chomp($line);

  my @pieces = split / /, $line;

  my $hitter = shift @pieces;

  foreach my $hitStr (@pieces)
  {
    my @data = split /!/, $hitStr;
    my $hittee = $data[0];
    my $number = $geneToNumber{$hittee};
    my $bucket;
    { # limit scope of "use integer"
      use integer;
      $bucket = $number / $BUCKET_SIZE;
    }
    my $file = $buckets[$bucket];
    print $file "$hittee $hitter\n";
  }
}
close(IN);

# close all the buckets
$i = 0;
while ($i < $numberOfBuckets)
{
  my $file = $buckets[$i];
  close($file);
  $i += 1;
}
print "Buckets filled.\n";

# open output file
my $outputFile = $ARGV[1];
open(OUT, ">", $outputFile) or
  die "cannot open output ($outputFile)\n";

# open and process the buckets in order
$i = 0;
while ($i < $numberOfBuckets)
{
  my $filename = $filenames[$i];
  open(BUCKET, "<", $filename) or
    die "cannot open bucket ($filename)\n";
  print "Processing bucket $i ($filename)\n";

  my %reverseHits = ();

  while (my $line = <BUCKET>)
  {
    chomp($line);

    my ($hittee, $hitter) = split / /, $line;

    my $old = $reverseHits{$hittee};
    if (defined($old))
    {
      $reverseHits{$hittee} = "$old $hitter";
    }
    else
    {
      $reverseHits{$hittee} = $hitter;
    }
  }
  close(BUCKET);

  # now print out the reverse hits in the same order as the input file
  my $geneNumber = $i * $BUCKET_SIZE;
  my $j = 0;
  while ($j < $BUCKET_SIZE && $geneNumber < $totalGeneCount)
  {
    my $geneName = $numberToGene{$geneNumber};
    
    my $data = $reverseHits{$geneName};
    if (defined($data))
    {
      print OUT "$geneName $data\n";
    }
    else
    {
      print OUT "$geneName\n";
    }
    delete $reverseHits{$geneName};

    $j += 1;
    $geneNumber += 1;
  }

  # sanity check: reverseHits hash should now be empty
  while (my ($key, $value) = each %reverseHits)
  {
    die "Bucket $i has extra genes: $key\n";
  }
  
  $i += 1;
}

close(OUT);

# delete all the buckets
$i = 0;
while ($i < $numberOfBuckets)
{
  my $filename = $filenames[$i];
  unlink($filename);
  $i += 1;
}

my $stopTime = time();
my $executionTime = $stopTime - $startTime;
print "execution complete after $executionTime seconds.\n";
