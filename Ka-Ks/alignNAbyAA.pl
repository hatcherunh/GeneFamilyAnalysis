#! /usr/bin/perl
#
# $Id: alignNAbyAA.pl 6 2013-11-16 01:42:31Z pjh $
#
# alignNAbyAA aligned_protein_file unaligned_nucleotide_file 
#
# This script takes a protein alignment in FASTA format, and using a FASTA
# file containing the nucleotide sequences of the aligned genes, produces
# a FASTA file containing the aligned nucleotide sequences.
#
# Sam Vohr (svohr@unh.edu)
# March 28, 2008
#

use warnings;
use strict;

sub usage()
{
    print STDERR "Usage: alignNAbyAA "
                 . "aligned_protein_file unaligned_nucleotide_file\n"; 
    exit;
}


# Hash table to store codons associated with amino acids. This is used to 
# check the correctness of the output data. Adapted from "Perl for 
# Exploring DNA", page 185.
#
# U -> A
# A -> T
# C -> G
# G -> C
#
my %amino_acids = 
  ( "TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
    "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
    "TAT" => "Y", "TAC" => "Y", "TAA" => "(STOP)", "TAG" => "(STOP)",
    "TGT" => "C", "TGC" => "C", "TGA" => "(STOP)", "TGG" => "W",
    "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
    "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
    "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
    "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
    "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
    "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
    "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
    "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
    "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
    "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
    "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
    "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G" );

my %na_sequences;

# builds a table with a nucleotide sequence for each gene number.
sub build_NA_table()
{
    my $header = <NA_FILE>;
    
    while ($header)
    {
    
        # pjh: always a problem: what is the structure of the fasta header? 
        if (!($header =~ m/^>([^\s\t]+)/))
        {
          die "can't parse FASTA header for nucleotide sequence file!\n";
        }
        my $genename = $1; 

        # pjh: mega-hack to address fact that clustalw2 re-writes colons in
        # gene names to underscores
        $genename =~ tr/:/_/;
        
        # read all lines of the sequence
        my $seq = "";
        my $in = <NA_FILE>;
        while ( $in && substr( $in, 0, 1 ) ne '>' )
        {
            chomp($in);
            $seq .= $in;
            $in = <NA_FILE>;
        }
 
        # pjh March 2009
        # promote nucleotides to uppercase
        $seq =~ tr/acgtnysvwhkdbrm/ACGTNYSVWHKDBRM/;

        $header = $in;
        $na_sequences{$genename} = $seq;
        #print STDERR "hashing nucleotide sequence for $genename\n";
    }
    #print STDERR "=================================================\n";
}

sub main()
{
    # check arguments
    if ( @ARGV < 2 || $ARGV[0] eq "-h" )
    {
        &usage();
    }

    # open files
    my $aa_filename = $ARGV[0];
    my $na_filename = $ARGV[1];
    my $outfile = substr ($na_filename, 0, -4) . "-aligned.fna";

    open ( AA_FILE,$aa_filename ) or die "Unable to open: ".$aa_filename ;
    open ( NA_FILE,$na_filename ) or die "Unable to open: ".$na_filename ;
    open (OUT, ">$outfile");

    # build hash table for nucleic acid sequences.
    &build_NA_table();

    # for each protein sequence build an nucleic acid alignment.
    my $header = <AA_FILE>;
    
    while ($header)
    {
        # pjh: always a problem: what is the structure of the fasta header? 
        if (!($header =~ m/^>([^\s\t]+)/))
        {
          die "can't parse FASTA header for aligned protein file!\n";
        }
        my $genename = $1; 
        my $seq = "";

        my $in = <AA_FILE>;

        # read the entire amino acid sequence.
        while ($in  && substr ($in, 0, 1) ne ">")
        {
            chomp($in);
            $seq .= $in;
            $in = <AA_FILE>;
        }
        
        $header = $in;

        # initialize the aligned nucleic acid sequence and fetch the
        # unaligned sequence.
        my $aligned_na = "";
        my $unaligned_na = $na_sequences{$genename};
        if (!($unaligned_na))
        {
          die "can't get nucleotides for $genename!\n";
        }
        my $pos = 1;
        
        # for each amino acid in the aligned sequence...
        while (length ($seq) > 0) 
        {
            my $aacid = substr($seq, 0, 1);
           
            # if there is a gap in the AA sequence...
            if ($aacid eq "-")
            {
                # ...place 3 gaps in the NA sequence.
                $aligned_na .= "---";
            } 
            else
            {
                my $codon = substr ($unaligned_na, 0, 3);
                # pjh Oct 2008: pad out last codon, if necessary
                while (length($codon) < 3)
                {
                  $codon .= "N";
                  $unaligned_na .= "N";
                }
                my $translatedCodon = $amino_acids{$codon};
                if (!($translatedCodon))
                {
                  if ($codon eq "NNN")
                  {
                    $translatedCodon = "*";
                  }
                  else
                  {
                    # pjh Oct 2008: deal with wildcards ala checkInputs.pl?
                    my $modifiedCodon = $codon;
                    $modifiedCodon =~ s/N/A/g;
                    $modifiedCodon =~ s/Y/C/g;
                    $modifiedCodon =~ s/S/C/g;
                    $modifiedCodon =~ s/V/C/g;
                    $modifiedCodon =~ s/W/A/g;
                    $modifiedCodon =~ s/H/A/g;
                    $modifiedCodon =~ s/K/G/g;
                    $modifiedCodon =~ s/D/G/g;
                    $modifiedCodon =~ s/B/G/g;
                    $modifiedCodon =~ s/R/A/g;
                    $modifiedCodon =~ s/M/A/g;

                    $translatedCodon = $amino_acids{$modifiedCodon};
                    if (!($translatedCodon))
                    {
                      die "can't get amino acid for $codon in $genename!\n";
                    }
                  }
                }
                if ($translatedCodon eq $aacid)
                {
                    $aligned_na .= $codon;
                    $unaligned_na = substr ($unaligned_na, 3);
                }
                else
                {
                    $aligned_na .= $codon;
                    $unaligned_na = substr ($unaligned_na, 3);

                    # pjh: check for wildcards
                    # Oct 2008: not sure this is still necessary
                    if ($aacid ne "X" && $codon ne "NNN")
                    {
                      if ($aacid eq "M")
                      {
                        if ($pos == 1)
                        {
                          print STDERR
                          "$na_filename: start codon problem ($aacid $codon)\n";
                        }
                        else
                        {
                          print STDERR
                          "$na_filename: codon $codon does not match M\n";
                        }
                      }
                      else
                      {
                        print STDERR
                          "$na_filename $genename $pos: $codon $aacid\n";
                      }
                    }
                }
            }

            $seq = substr ($seq, 1);
            $pos += 1;
        }

        if ($unaligned_na ne "")
        {
            if ($amino_acids{$unaligned_na} eq "(STOP)") {
                $aligned_na .= $unaligned_na;
            }
            else {
                print STDERR "weird: NA sequence longer than AA sequence\n";
                print STDERR "$na_filename $genename: $unaligned_na\n";
            }
        }

        # print the header.
        print OUT ">$genename\n";
        
        # print the sequence.
        foreach my $l ( unpack( "(A60)*", $aligned_na ) )
        {
            print OUT "$l\n";
        }
    }

    close OUT;

    # done.
}

# call the main function.
&main();

