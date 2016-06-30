#!/usr/bin/perl -w
#
# $Id: alignAllNAbyAAFamilies.pl 6 2013-11-16 01:42:31Z pjh $
#
# This script opens the directory 'family' (the results of the 
# makeFastaByFamily script) and runs alignNAbyAA using the corresponding
# aligned protein file.
# on all gene files.
#
# Sam Vohr (svohr@unh.edu)
# March 11, 2008
#
# Phil Hatcher, August 2008
# made some edits to work better with my pipeline
# Note: you pass this script the nucleotide family directory
# and the aligned protein directory.
#

use warnings;
use strict;

sub usage()
{
  print STDERR
    "Usage: alignAllNAbyAAFamilies nuc-family-dir aligned-protein-dir\n" ;
  exit;
}

if ( @ARGV < 2 || $ARGV[0] eq "-h" )
{
  &usage();
}


sub main ()
{
  # open a family file
  my $familyDir = $ARGV[0];
  my $proteinDir = $ARGV[1];
  opendir (FAMILIES,"$familyDir");

  if ( substr($familyDir, -1) eq "/" ) {
    chop $familyDir;
  }

  # create a directory for results.
  my $resultsDir = $familyDir . "-aligned";
  mkdir "$resultsDir", 0755 or die "Could not create directory";

  # get all the family names
  my @files = grep (/^[^\.]*.nuc/,readdir(FAMILIES));
 #print STDERR "$_\n" foreach (@files);

  foreach my $file (@files) {
    # do the alignment
    my $AAfile = substr($file, 0, -4) . ".fasta"; 
    my @script = ("alignNAbyAA.pl", "$proteinDir/$AAfile", "$familyDir/$file");
    my $err = system(@script);
    if ($err != 0)
    {
      die "alignNAbyAA.p failed for $proteinDir/$AAfile $familyDir/$file\n";
    }
    
    # move the resulting file to the results directory
    my $aligned = substr($file, 0, -4) . "-aligned.fna";
    my @move = ("mv", 
                "$familyDir/$aligned", 
                "$resultsDir/$file");
    $err = system(@move);
    if ($err != 0)
    {
      die
        "move of results ($familyDir/$aligned to $resultsDir/$file) failed\n";
    }
  }

  closedir FAMILIES;

  # ta da!
}

&main();

