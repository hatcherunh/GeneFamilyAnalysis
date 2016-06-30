#!/usr/bin/perl -w
# $Id$
#
#
# This script opens a "family" directory of nucleotide sequences and runs
# clustalw2 to create an alignment in Phylip format. Then dnaml from Phylip
# is run in order to generate a Phylogenetic tree.
#
# This version assumes the gene names have already been shortened.

use warnings;
use strict;

sub usage()
{
  print STDERR "Usage: phylipFormatAllFamilies.pl family-dir\n";
  exit;
}

if ( @ARGV == 0 || $ARGV[0] eq "-h" )
{
  &usage();
}

# pjh Dec 2013
# what is this function doing? and why?
# it is just counting off inner parentheses in the tree file
# apparently codeml wants this
# but with some datasets, it does nothing because there are not
# any nested parentheses.
sub addBranchLabels ()
{
  open TREE,"outtree" || die "addBranchLabels: can't open outtree for input\n";
  
  my @treelines = <TREE>;
  
  close TREE;

  my $tree = "";

  foreach my $line (@treelines) {
    $tree .= $line;
  }

  my $result;
  my $stack; # pjh Dec 2013: only set and never read

  my $branchLabel = 1;
  for ( my $i = 0; $i < length $tree; $i++ ) {
    my $tmp = substr ($tree, $i, 1);
    if ($tmp eq "(") { 
        $stack .= "(";
        $result .= "("; 
    }
    elsif ($tmp eq ")") { 
      if ( substr ($tree, $i+1, 1) ne ";" ) {
        $result .= ")\#$branchLabel";
        $stack .= ")";
        $branchLabel += 1;
      }
      else {
        $stack .= ")";
        $result .= ")";
      }
    }
    else {
      $result .= substr ($tree, "$i", 1);
    }
  }
  
  open TREE,"+>outtree" ||
    die "addBranchLabels: can't open outtree for output\n";
  
  print TREE "$result";

  close TREE;
}    


sub main ()
{
  # open a directory of family files
  my $familyDir = $ARGV[0];
  opendir (FAMILIES,"$familyDir") || die "can't open $familyDir\n";
  if ( substr($familyDir, -1) eq "/" ) {
    chop $familyDir;
  }

  # get all the family names
  my @files = grep (/^[^\.].+\.nuc/,readdir(FAMILIES));
 #print STDERR "$_\n" foreach (@files);
  
  my $resultsDir = $familyDir . "-trees";
    
  mkdir "$resultsDir", 0755 or die "Could not create directory";

  # get all the family names

  foreach my $file (@files) {
    
    # do the alignment
    my @clustal_phylip = (  
      "clustalw2 -OUTPUT=phylip -INFILE=$familyDir/$file -OUTORDER=INPUT");
    system(@clustal_phylip);
    
    $file = substr($file, 0, -4) . ".phy";
 
    my $copy_from = "cp $familyDir/$file infile";
    system($copy_from);

    # generate a tree
    my $dnaml = "yes | dnaml"; 
    system($dnaml);

    &addBranchLabels;
    
    # move the resulting file to the results directory
    my $newtree = substr($file, 0, -4);
    my $copy_back = "mv outtree $resultsDir/$newtree.tree";
    system($copy_back);
    $copy_back = "mv outfile $resultsDir/$newtree.file";
    system($copy_back);
    
  }

  closedir FAMILIES;

  # ta da!
}

&main();

