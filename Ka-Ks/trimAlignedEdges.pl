#! /usr/bin/perl
#
# $Id: trimAlignedEdges.pl 6 2013-11-16 01:42:31Z pjh $
#
# trimAlignedEdges <aligned_NA_dir> <aligned_AA_dir>
#
# This script takes an alignment in FASTA format and removes the leading 
# characters of all sequences equal corresponding to the longest sequence of
# '-' characters appearing at the start of the sequences. Does the same thing
# for the end of the sequence. This is kinda confusing...
# 
# The aligned sequences:
#   AGGA-CTTA --GAGCTTA -GGACCAT- 
# will be trimmed to:
#   GA-CTT GAGCTT GACCAT
#
# if that makes any sense.
#
# Sam Vohr (svohr@unh.edu)
# March 30, 2008
# Edited by Phil Hatcher, Oct 2010.
# Edited by Phil Hatcher, Jul 2013 in order to maintain the order of
# the genes in the output. This is important downstream, for instance
# to properly identify the outgroup by position, as Phylip requires.
# The code only needed to be modified for the amino acid sequences, as
# the nucleotide sequences were not being scrambled.

use warnings;
use strict;

sub usage()
{
    print STDERR "trimAlignedEdges <aligned_AA_dir> <aligned_NA_dir>\n";
    exit;
}

if ( @ARGV < 2 || $ARGV[0] eq "-h" )
{
    &usage();
} 

# global path variables so we don't need to pass 4 args into cutEdges().
my $AAdir = $ARGV[0];
my $NAdir = $ARGV[1];

if ( substr($AAdir,-1) eq "/" ) {
    chop $AAdir;
}
if ( substr($NAdir,-1) eq "/" ) {
    chop $NAdir;
}

my $AAresdir = $AAdir . "-edged";
my $NAresdir = $NAdir . "-edged";

# sorry :(


sub cutEdges ()
{
    # open files
    my $AAfilename = $_[0];
    my $AAout = substr ($AAfilename, 0, -6) . ".faa";
    
    #my $NAfilename = substr ($AAfilename, 0, -6) . ".fna";
    my $NAfilename = substr ($AAfilename, 0, -6) . ".nuc";
    my $NAout = $NAfilename;

    open ( AASEQ, "$AAdir/$AAfilename" ) or die "Unable to open: ".$AAfilename; 
    open ( AAOUT, ">$AAresdir/$AAout" );

    # hash to store the gene sequences
    my %table;

    # store the max sequence of '-' chars
    my $cut_lead = 0;
    my $cut_tail = 0;

    # need to work to make the order of the genes in the output the
    # same as in the input
    my $header = <AASEQ>;
    my @order = ();
    my $cnt = 0;
    
    while ($header)
    {
        $order[$cnt] = $header;
        $cnt += 1;

        # get the gene's number
        my $seq = "";
        my $in = <AASEQ>;

        # read all lines of the sequence
        while ( $in && substr( $in, 0, 1 ) ne '>' )
        {
            chomp($in);
            $seq .= $in;
            $in = <AASEQ>;
        }
       
        # if there a leading '-' see if there are more than
        # already have been seen.
        if ($seq =~ m/(^-+)/)
        {
            my $len = length $1;
            if ($len > $cut_lead)
            {
                $cut_lead = $len;
            }
        }
        # do the same for the trailing '-'
        if ($seq =~ m/(-+$)/)
        {
            my $len = length $1;
            if ($len > $cut_tail)
            {
                $cut_tail = $len;
            }
        }
       
        $table{$header} = $seq;
        $header = $in;
    }

    print STDERR "$cut_lead $cut_tail\n";
    for (my $i = 0; $i < $cnt; $i += 1)
    {
        my $gene = $order[$i];
        my $seq = $table{$gene};
        # print the header.
        print AAOUT "$gene"; 
        
        # print the sequence.
        $seq = substr ($seq, $cut_lead);
        $seq = substr ($seq, 0, (length($seq) - $cut_tail));
        foreach my $l ( unpack( "(A60)*", $seq ) )
        {
            print AAOUT "$l\n";
        }
    }
    
    close AAOUT;
    close AASEQ;

    # trim the edges of the nucleotide sequences now
    open ( NASEQ, "$NAdir/$NAfilename" ) or die "Unable to open: ".$NAfilename; 
    open ( NAOUT, ">$NAresdir/$NAout" );
    
    $header = <NASEQ>;
   
    # for nucleotides triple the number of characters to trim
    # (1 amino acid = 1 codon = 3 nucleotides)
    $cut_lead *= 3;
    $cut_tail *= 3;
 
    while ($header)
    {
        # get the gene's number
        my $seq = "";
        my $in = <NASEQ>;

        # read all lines of the sequence
        while ( $in && substr( $in, 0, 1 ) ne '>' )
        {
            chomp($in);
            $seq .= $in;
            $in = <NASEQ>;
        }

        # print the header
        print NAOUT "$header";
        
        # print the sequence.
        $seq = substr ($seq, $cut_lead);
        $seq = substr ($seq, 0, (length($seq) - $cut_tail));
        foreach my $l ( unpack( "(A60)*", $seq ) )
        {
            print NAOUT "$l\n";
        }
        
        $header = $in;
    }
 
    close NAOUT;
    close NASEQ;

    # done.
}

sub main()
{
    # create new result directories
    mkdir "$AAresdir", 0755;
    mkdir "$NAresdir", 0755;

    # open AA dir
    opendir (AADIR, "$AAdir");

    my @files = grep (/^[^\.].*\.fasta/,readdir(AADIR));

    foreach my $file (@files) {
        &cutEdges($file);
    }
    
    # done?
    close AADIR;
}


# call the main function.
&main();

