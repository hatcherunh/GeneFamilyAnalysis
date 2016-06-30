#! /usr/bin/perl
# $Id: maxDiffFromConsensus.pl 38 2015-02-10 13:32:07Z pjh $
#
#
# Compares every gene in each family against its consensus sequence. It
# records the maximum consecutive difference found for each family. This
# can be used later to screen out families from the analysis.
#
# maxDiffFromConsensus.pl <AA-aligned> <AA-consensus> <outputFile>
#
# The format of the output file is simply a line for each family with
#   familyNum maxDiff
#
# Created by Phil Hatcher, October 2010, from Sam Vohr's binByResidue.pl.
#

use warnings;
use strict;

sub usage()
{
  print STDERR
    "Usage: maxDiffFromConsensus.pl <AA-aligned> <AA-consensus> <outFile>\n";
  exit;
}

if ( @ARGV != 3 || $ARGV[0] eq "-h" )
{
  &usage();
}

# reads the consensus sequence for a family file
sub getConsensus()
{
  my $file = $_[0];
  my $dir = $_[1];

  # open the file and read the header
  open CONS,"$dir/$file" or die "missing consensus file for $file\n" .
                                "$dir/$file\n";
  my $header = <CONS>;

  # read the sequence
  my $seq = "";
  my $in = <CONS>;
  while ($in && substr($in, 0, 1) ne '>') { 
    chomp $in;
    $seq .= $in;
    $in = <CONS>;
  }

  close CONS;

  return (uc $seq);
}

sub main ()
{
  # for curiousity keep the "max of the max" and report to stdout
  my $maxOfMax = 0;

  # process directory names
  my $familyDir = $ARGV[0];
  my $consensusDir = $ARGV[1];
  my $outputFile =  $ARGV[2];  

  # open output file
  open (OUT, ">$outputFile") or
    die "cannot open output file: $outputFile\n";

  # clean up the directory names if necessary.
  if (substr($familyDir, -1) eq "/")
  {
    $familyDir = substr($familyDir, 0, -1);
  }
  if (substr($consensusDir, -1) eq "/")
  {
    $consensusDir = substr($consensusDir, 0, -1);
  }

  # open the aligned directory and read each family
  opendir (FAMILIES,"$familyDir") or
    die "can't open $familyDir\n";
  #my @families = grep (/^[^\.].*\.faa/, readdir(FAMILIES));
  # pjh Feb 2015: generalize to allow either .fasta or .faa
  my @families = grep (/^[^\.].*\.fa/, readdir(FAMILIES));

  foreach my $file (@families)
  {
    # track the largest mismatch in the family
    my $maxMismatch = 0;

    # get the consensus sequence
    my $consensusSeq = &getConsensus($file,$consensusDir);

    # compare all sequences
    open (IN,"$familyDir/$file") or
      die "can't open $familyDir/$file";

    print STDERR "reading first line from $familyDir/$file\n";

    my $in = <IN>;

    if (substr ($in,0,1) ne ">")
    {
      die "expected fasta formatting\n";
    }

    while ($in)
    {
      # build the sequence
      $in = <IN>;
      my $seq = "";
      while ( $in && substr ($in, 0, 1) ne ">" )
      {
        chomp($in);
        $seq .= $in;
        $in = <IN>;
      }
      
      # got the sequence. compare it to the consensus
      my $mismatch = 0;
      my $matchstr = "";
      my $conSeq = $consensusSeq;
      while (length ($seq) > 0)
      {
        my $seqResidue = substr($seq, 0, 1);
        my $conResidue = substr($conSeq, 0, 1);

        # check for a match. update count
        if ( $seqResidue eq $conResidue )
        {
          $matchstr .= "=";

          # once I have a match, reset my counter
          # because we are interested in the maximum consecutive mismatch
          $mismatch = 0;
        }
        else
        {
          $matchstr .= $seqResidue;
          $mismatch += 1;

          # is this the biggest consecutive mismatch that I have seen so far?
          if ($mismatch > $maxMismatch)
          {
            $maxMismatch = $mismatch;
          }
        }

        $seq = substr ($seq, 1);
        $conSeq = substr ($conSeq, 1);
      }

      print "$file: $maxMismatch\n";
      print "$matchstr\n";
    }

    # get family number: discard the ".faa" from the file name
    # pjh Feb 2015
    #my $familyNum = substr($file, 0, -4);
    my $familyNum;
    if ($file =~ /^(\d*)/)
    {
      $familyNum = $1;
    }
    else
    {
      die "can't parse filename ($file)\n";
    }

    # output the max difference from the consensus
    print OUT "$familyNum $maxMismatch\n";
    print "** $familyNum: $maxMismatch\n\n";

    # is this the biggest max I have seen so far
    if ($maxMismatch > $maxOfMax)
    {
      $maxOfMax = $maxMismatch;
    }

    close IN;
  }

  close OUT;

  print "maximum difference seen across all families is $maxOfMax\n";

  # and we're done!
}


&main();

