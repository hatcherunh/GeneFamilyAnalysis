#!/usr/bin/perl -w
#
# $Id$
#
# This script reads a directory of FASTA files and edits the gene names to
# be generic genome names (ie genome1, genome2, etc). This is to work around
# the limitation that phylip-formatted files only have 10 characters for
# the gene name. Also it works around codeml's dislike for numeric
# gene names, which result in this kind of error message:
#   species number 638203616 outside range
#

use warnings;
use strict;

sub usage()
{
  print STDERR "Usage: replaceGeneNames2.pl file-extension family-dir\n";
  exit;
}

if ( @ARGV < 2 || $ARGV[0] eq "-h" )
{
  &usage();
}

sub main ()
{
  # open a family file
  my $fileExtension = $ARGV[0];
  my $familyDir = $ARGV[1];

  opendir (FAMILIES,"$familyDir");

  if ( substr($familyDir, -1) eq "/" )
  {
    chop $familyDir;
  }

  # create results directories
  my $resultsDir = $familyDir . "-renamed";
  mkdir "$resultsDir", 0755 or die "Could not create directory $resultsDir";

  # get all the family names
  my @files = grep (/^[^\.]*\.$fileExtension/,readdir(FAMILIES));

  foreach my $file (@files)
  {
    # i.e. the file name is expected to be based on the family number
    # with $fileExtension for the extension. The "plus 1" is for the
    # dot before the extension.
    my $familyNum = substr($file,0,-(length($fileExtension)+1));
    
    # read old FASTA file, write new FASTA file, and
    # build hash to map gene names to generic genome names
    open (SEQ, ">$resultsDir/$file");
    open (IN, "$familyDir/$file");
    my %geneHash = ();
    my $i = 0;
    while (my $line = <IN>)
    {
      chomp($line);

      # FASTA header?
      if ($line =~ /^>\s*([^\s]*)/)
      {
        $i += 1;
        $geneHash{$1} = "genome$i";
        print SEQ ">genome$i\n";
      }
      else
      {
        print SEQ "$line\n";
      }
    }
    close SEQ;
    close IN;
  }

  closedir FAMILIES;

  # ta da!
}

&main();

