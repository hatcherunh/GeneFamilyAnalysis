#!/usr/bin/perl -w
#
# $Id: consensusAllFamilies.pl 41 2015-05-20 19:59:27Z pjh $
#
# This script opens the directory 'family' and runs cons from
# EMBOSS on a set of aligned genes.
#
# Sam Vohr (svohr@unh.edu)
# March 11, 2008
# edited by Phil Hatcher, October 2010 (including the above comment
# which was bogus before).
#

use warnings;
use strict;

sub usage()
{
  print STDERR "Usage: consensusAllFamilies family-dir\n" .
               "note: cons from EMBOSS is required in path\n";
  exit;
}

if ( @ARGV < 1 || $ARGV[0] eq "-h" )
{
  &usage();
}


sub main ()
{
  # open a family file
  my $familyDir = $ARGV[0];
  opendir (FAMILIES,"$familyDir");

  # create a directory for results.
  chop $familyDir if (substr($familyDir, -1) eq "/");
  my $resultsDir = $familyDir . "-consensus";
  mkdir "$resultsDir", 0755 or die "Could not create directory";

  # get all the family names
  #my @files = grep (/^[^\.]*\.faa/,readdir(FAMILIES));
  # pjh Feb 2015: generalize and accept either .faa or .fasta
  my @files = grep (/^[^\.]*\.fa/,readdir(FAMILIES));
  #print STDERR "$_\n" foreach (@files);

  foreach my $file (@files) {
    # get the consensus sequence
    system("cons $familyDir/$file $resultsDir/$file");
  }

  closedir FAMILIES;

  # ta da!
}

&main();

