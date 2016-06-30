#!/usr/bin/perl -w
#
#$Id: codemlAllFamilies.pl 17 2013-11-26 18:18:13Z pjh $
#
# This script runs codeml on all families contained the the specified directory.
#
# Sam Vohr (svohr@unh.edu)
# July 15, 2008
# Edited by Phil Hatcher, Oct 2010.

use warnings;
use strict;


# This is the configuration file template for codeml. the seqfile and resfile
# placeholders will be replaced with actual filenames.

my $configTemplate =
q[      seqfile = <seqfile>
     treefile = <treefile>
      outfile = <resfile>           * main result file name

        noisy = 3  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        model = 2
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 4.54006   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates

        getSE = 0   * 0: do n0t want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)

  fix_blength = 1  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0   * 0: simultaneous; 1: one branch at a time

* Specifications for duplicating results for the small data set in table 1
* of Yang (1998 MBE 15:568-573).
* see the tree file lysozyme.trees for specification of node (branch) labels ];

 



sub usage()
{
  print STDERR "Usage: codemlAllFamilies family-dir tree-dir\n";
  exit;
}

if ( @ARGV < 2 || $ARGV[0] eq "-h" )
{
  &usage();
}

sub generateInput () 
{
  my $file = $_[0];
  my $input = "";

  # keep track of the length of each sequence, and number of sequences
  my $length = 0;
  my $count = 0;

  # read in the contents of the file
  open (IN,$file);

  my $in = <IN>;

  # read in the file
  while ($in) {
    $input .= $in;
    $count += 1;
    $length = 0;
    
    # read the next sequence
    $in = <IN>;
    while ( $in && substr ($in, 0, 1) ne ">" ) {
      $length += length($in) - 1;
      
      if (length $in == 4) {
        # For some reason, the stop codon cannot have it's own 
        chomp $input;
      }
      $input .= $in;
      $in = <IN>;
    }
  }

  # don't count the stop codon.
  $length -= 3;  
  $input = "$count $length\n" . $input;

  return $input;
}

sub main ()
{
  # open a family file
  my $familyDir = $ARGV[0];
  my $treeDir = $ARGV[1];

  opendir (FAMILIES,"$familyDir");

  if ( substr($familyDir, -1) eq "/" ) {
    chop $familyDir;
  }
  if ( substr($treeDir, -1) eq "/" ) {
    chop $treeDir;
  }

  # create results directory
  my $resultsDir = $familyDir . "-codeml";
  mkdir "$resultsDir", 0755 or die "Could not create directory";

  # get all the family names
  my @files = grep (/^[^\.]*\.nuc/,readdir(FAMILIES));
 #print STDERR "$_\n" foreach (@files);

  foreach my $file (@files) {
    my $familyNum = substr($file,0,-4);
    mkdir "$resultsDir/$familyNum", 0755;
    
    # build the input file.
    open (SEQ, ">$resultsDir/$familyNum/$file");
    my $input = &generateInput("$familyDir/$file");
    print SEQ "$input";
    close SEQ;

    # copy the tree file
    my $treeFile = "$familyNum.tree";
    system "cp $treeDir/$treeFile $resultsDir/$familyNum/";
    
    # create the configuration file
    open (CFG, ">$resultsDir/$familyNum/codeml.ctl");
    
    my $configText = $configTemplate;
    $configText =~ s/<seqfile>/$file/g;
    $configText =~ s/<treefile>/$treeFile/g;
    $configText =~ s/<resfile>/$familyNum-result.txt/g;
    
    print CFG "$configText";
    close CFG;

    # move to the results directory and execute PAML codeml.
    chdir "$resultsDir/$familyNum";
    system ("codeml");
    
    if ( $? == -1 ) {
      print STDERR "codeml failed to execute: $!\n";
    }
    elsif ( ($? >> 8) != 0) {
      print STDERR "ERROR: codeml did not successfully complete on family " .
            "$familyNum.\n";
    }

    chdir "../..";
  }

  closedir FAMILIES;

  # ta da!
}

&main();

