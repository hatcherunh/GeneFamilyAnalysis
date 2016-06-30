#!/usr/bin/perl -w
# $Id: makeFastaByFamily2.pl 42 2015-05-22 19:16:20Z pjh $
#
# THIS PROGRAM WAS NOT FINISHED!
#
# This script takes a family file input and generates a fasta file 
# for the family containing the DNA sequence of each member.
#
# This has been created from makeFastaByFamily.pl in order to allow
# it to be generalized and work with general families and not just
# panortholog families.
#
# The only challenge is to properly order the genes in the input
# to clustalw. (This ordering is only important when analyzing two
# different sets of families produced by two different gene family
# analysis, where we want the results to be consistent for two
# identical families.) I sort first by genome, using the gene
# order file, then I sort by gene name in each genome.
#
# The command-line arguments:
#   1. family file
#   2. genome order file (one genome per line), which will control the
#      order of genes in the output files. Clustal is apparently sensitive
#      to what order the genes are in the file passed to clustal.
#   3. file extension for input gene sequence files ("nuc" or "proteins")
#   4. directory where FASTA file for gene sequences are found.
#   5. prefix to use for output directory name

use warnings;
use strict;

sub usage()
{
  print STDERR
  "Usage: makeFastaByFamily2.pl family-file genome-order-file " .
  " <proteins|nuc> gene-dir output-prefix\n";
  exit;
}

if ( @ARGV != 5 || $ARGV[0] eq "-h" )
{
  &usage();
}

my $familyFile = $ARGV[0];
my $genomeOrderFile = $ARGV[1];
my $type = $ARGV[2];
my $geneDir = $ARGV[3];
my $outputPrefix = $ARGV[4];

my %taxaHash = ();

# Finds a gene from a fasta file. Requires filehandle OUT to work.
sub findGeneInFile () 
{
  # first argument is a string formatted like:
  #   taxa$genenum
  my ($taxa,$gene) = split (/\$/,$_[0]);
  #print STDERR "$taxa $gene\n";

  # second argument is the family number, which is only used for error
  # reporting
  my $familyNum = $_[1];

  my $ref = $taxaHash{$taxa};
  if (!$ref)
  {
    print "creating hash for $taxa...\n";

    open GENE,"$geneDir/$taxa\.$type" or
      die "cannot find gene file $taxa\.$type\n";

    my %geneHash = ();
    $taxaHash{$taxa} = \%geneHash;

    # go through each gene
    my $line;
    my $insert = "";
    my $geneName;
    while ($line = <GENE>)
    {
      # check if this line is a header
      if ($line =~ /^>\s*([^\s\t]*)/)
      {
        if ($insert ne "")
        {
          $geneHash{$geneName} = $insert;
        }
        $geneName = $1;
        $insert = "";
      }
      $insert .= $line;
    }
    if ($insert ne "")
    {
      $geneHash{$geneName} = $insert;
    }
    close(GENE);
    $ref = \%geneHash;
    print "  Done.\n";
  }
  my $out = ${$ref}{$gene};
  if ($out)
  {
    print OUT "$out";
  }
  else
  {
    print STDERR
      "family $familyNum: gene $gene not found in hash for taxa $taxa\n";
  }
}

sub padNum ()
{
  my $n = $_[0];
  my $res = "$n";
  $res = "0".$res if ($n < 10) ;
  $res = "0".$res if ($n < 100) ;
  $res = "0".$res if ($n < 1000) ;
  $res = "0".$res if ($n < 10000) ;
  $res = "0".$res if ($n < 100000) ;
  $res;
}

sub main ()
{
  # open and read the genome-order file
  my @order = ();
  my $numberOfTaxa = 0;
  open (ORDER,"$genomeOrderFile") or
    die "can't open $genomeOrderFile for input\n";
  while (my $genomeName = <ORDER>)
  {
    chomp($genomeName);
    $order[$numberOfTaxa] = $genomeName;
    $numberOfTaxa += 1;
  }
  close(ORDER);
  
  # open a family file
  open (FAMILIES,"$familyFile");
  
  # create a directory for results.
  my $familyDir = "$outputPrefix";
  if ($type eq "nuc") {
    $familyDir .= "-NA";
  }
  elsif ($type eq "proteins") {
    $familyDir .= "-AA";
  }
  else
  {
    die "unexpected type argument: $type\n";
  }

  mkdir "$familyDir", 0755 or die "Could not create directory";

  my $line = <FAMILIES>;

  while ($line) {
    # read a line
    my ($familyNum, @genes) = split (/:?\s+/,$line);

    # create an output file   
    my $file = &padNum($familyNum);
    open (OUT, ">$familyDir/$file\.$type");

    # put the genes into a hash by taxon name
    # could be multiple genes per taxon
    my %hash = ();
    foreach my $gene (@genes) {
      my ($taxa,$gene2) = split (/\$/,$gene);
      my $old = $hash{$taxa};
      if (defined($old))
      {
        $hash{$taxa} = $old . "\$" . $gene2;
      }
      else
      {
        $hash{$taxa} = $gene2;
      }
    }

    # now process genes in the specified order
    # first, in genome order
    # second, sort gene names for genes in the same genome
    for (my $i = 0; $i < $numberOfTaxa; $i += 1)
    {
      my $genes = $hash{$order[$i]};
      if (defined($genes))
      {
        my @geneList = split (/\$/,$genes);
        my @sortedList = sort @geneList;
        foreach my $gene (@sortedList) {
          &findGeneInFile($order[$i] . "\$" . $gene, $familyNum);
        }
      }
    }

    close OUT;
    $line = <FAMILIES>;
  }

  # ta da!
}

&main();

