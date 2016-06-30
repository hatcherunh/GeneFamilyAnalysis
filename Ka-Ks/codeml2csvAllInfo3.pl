#!/usr/bin/perl
#
#$Id: codeml2csvAllInfo3.pl 42 2015-05-22 19:16:20Z pjh $
#
# May 2012, Phil Hatcher
# This is a hack of codeml2csvAllInfo2.pl to be used when there are
# no details files.
#
# Creates a .csv file for the codeml results.
#
# Sam Vohr (svohr@unh.edu)
# July 21, 2009
# Edited by Phil Hatcher in October 2010.
# Edited by Phil Hatcher in March 2012 to add a commandline argument to
#   specify where the details files are.

use strict;
use warnings;

#---------------------------------------------------------------------
# Constants
#my @taxa = ("Methanobrevibacter-ruminantium-M1","Methanobrevibacter-smithii-ATCC-35061","Methanothermobacter-thermautotrophicus-Delta-H","Methanosphaera-stadtmanae-DSM-3091");
#my $familyFile = "group2.panorthologs";

#---------------------------------------------------------------------

 
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

  if ( @ARGV != 4 || $ARGV[0] eq "-h" )
  {
    die
  "Usage: codeml2csvAllInfo3.pl codemlDir familyFile orderFile maxDiffFile\n";
  }

  my $pamlDir = $ARGV[0];
  my $familyFile = $ARGV[1];
  my $orderFile = $ARGV[2];
  my $maxDiffFile = $ARGV[3];


  open (ORDER,"$orderFile") or die "Cannot open $orderFile.";
  my @taxa = <ORDER>;
  close (ORDER);

  my %Genes = ();

  open (FAMILIES,"$familyFile") or die "Cannot open $familyFile.";
  
  foreach my $line (<FAMILIES>) {
    my ($num, @genes) = split (/:?\s+/,$line);
    $num = &padNum($num);
  
    foreach my $gene (@genes) {
      my ($taxon,$genenum) = split (/\$/,$gene);
      my $key = "$taxon\_$num";

      $Genes{ $key } = $genenum;
    }
  }
  close FAMILIES;
 
  my %AAdiff = ();
  open (DIFF,"$maxDiffFile") or die "Cannot open maxDiff file.";

  foreach my $line (<DIFF>)
  {
    chomp($line);

    my ($num, $max) = split (/\s/, $line);

    $AAdiff{$num} = $max;
  }

  opendir (PAML,$pamlDir);

  my @families = grep (/^[^\.].*/, readdir(PAML));

  print STDOUT "Family, MaxDiff, Ka, Ks";
  print STDOUT "\n";

  foreach my $family (sort @families) {
    my $dn = 0;
    my $ds = 0;
    my $wasMatch = 0;

    my $maxDiff = $AAdiff{$family};
    if (!defined $maxDiff)
    {
      print STDERR "family $family has no maxDiff?\n";
      $maxDiff = "";
    }

    open (IN,"$pamlDir/$family/$family-result.txt");
    foreach my $line (<IN>) {
      if ($line =~ m/tree\ length\ for\ dN\:\s+(\d+\.\d+)/) {
        $dn = $1;
        $wasMatch += 1;
      }
      if ($line =~ m/tree\ length\ for\ dS\:\s+(\d+\.\d+)/) {
        $ds = $1;
        $wasMatch += 1;
      }
    }
    close IN;

    if ( $wasMatch < 2 ) {
      print STDERR "$family did not have results ($wasMatch)\n";
    }
    else {
      print STDOUT "$family, $maxDiff, $dn, $ds";
      print STDOUT "\n";
    }
  }
}

&main;

