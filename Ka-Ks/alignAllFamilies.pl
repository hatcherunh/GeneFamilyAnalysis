#!/usr/bin/perl -w
#
# $Id: alignAllFamilies.pl 6 2013-11-16 01:42:31Z pjh $
#
# This script opens the directory 'family' (the results of the 
# makeFastaByFamily script) and runs clustalw2 multiple sequence alignment
# on all gene files.
#
# Sam Vohr (svohr@unh.edu)
# March 11, 2008
#
# Phil Hatcher, August 2008
# edits to better integrate the script into my pipeline
#
# Phil Hatcher, October 2008
# add command-line args to choose format of output files
# either "phylip" or "fasta"
#
# Note: this script should be passed the protein family directory
#

use warnings;
use strict;

sub usage()
{
  print STDERR
    "Usage: alignAllFamilies [phylip | fasta | clustal] protein-family-dir\n" .
    "note: clustalw2 must be in path\n";
  exit;
}

if ( @ARGV != 2 || $ARGV[0] eq "-h" )
{
  &usage();
}


sub main ()
{
  my $format;
  my $fileExtension;

  if ($ARGV[0] eq "phylip")
  {
    $format = "phylip";
    $fileExtension = "phy";
  }
  elsif ($ARGV[0] eq "fasta")
  {
    $format = "fasta";
    $fileExtension = "fasta";
  }
  elsif ($ARGV[0] eq "clustal")
  {
    $format = "clustal";
    $fileExtension = "aln";
  }
  else
  {
    usage();
  }

  # open a family file
  my $familyDir = $ARGV[1];
  opendir (FAMILIES,"$familyDir");

  if ( substr($familyDir, -1) eq "/" ) {
    chop $familyDir;
  }
  # create a directory for results.
  my $resultsDir = $familyDir . "-aligned";
  mkdir "$resultsDir", 0755 or die "Could not create directory";

  # get all the family names
  my @files = grep (/^[^\.]/,readdir(FAMILIES));
  #print STDERR "$_\n" foreach (@files);

  foreach my $file (@files) {
    # parse the file name
    my $fileStart;
    if ($file =~ /([0-9]+)[.]/)
    {
      $fileStart = $1;
    }
    else
    {
      die "cannot parse the input file name ($file)\n";
    }

    # do the alignment
    my @clustal_cmd = ("clustalw2 -OUTPUT=$format -OUTORDER=INPUT" .
                       " -INFILE=$familyDir/$file");
    my $err = system(@clustal_cmd);
    if ($err != 0)
    {
      die "Failed: clustalw2 -OUTPUT=$format -OUTORDER=INPUT" .
          " -INFILE=$familyDir/$file";
    }
    
    # move the resulting file to the results directory
    my $outFile = "$fileStart.$fileExtension";
    my @move = ("mv $familyDir/$outFile $resultsDir/");
    $err = system(@move);
    if ($err != 0)
    {
      die "Failed: mv $familyDir/$outFile $resultsDir/";
    }
    
    # remove the *.dnd file produced by clustal
    my $dndFile = "$fileStart.dnd";
    my @rm = ("rm $familyDir/$dndFile");
    $err = system(@rm);
    if ($err != 0)
    {
      die "Failed: rm $familyDir/$dndFile";
    }
  }

  closedir FAMILIES;

  # ta da!
}

&main();

