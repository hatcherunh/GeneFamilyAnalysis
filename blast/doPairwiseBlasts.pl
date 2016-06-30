#!/usr/bin/perl

# $Id: doPairwiseBlasts.pl 57 2015-07-22 21:41:52Z pjh $
#
#
# Phil Hatcher, June 2007
#
# Do pair-wise genome BLASTs involving a set of new genomes and a set
# of genomes that have already been processed.
#
# It takes six initial arguments:
#   1. optional -n (indicating that the FASTA sequences are nucleotides)
#   2. directory containing FASTA sequence files
#   3. directory containing BLAST results
#   4. evalue threshold to be passed via -e argument to blastall
#   5. number of MPI processes to use
#   6. MPI machine file
#
# The additional arguments are the new genomes to process. There must be at
# least one genome specified.
#
# The FASTA sequence files are assumed to be named either <GENOME>.proteins
# or <GENOME>.nuc.
#
# The blast directory should contain a file called "DONE", which gives
# list of genomes that have already been processed. This file will be
# updated by this script. This file must exist, even if it is empty.
#
# The BLAST results for each pairwise comparison are stored in a file
# called <GENOME1>-<GENOME2>.blast.
#
# And for each new genome these is a file created called <GENOME>.self that
# contains a line for each non-error gene, with the line containing that
# gene's self-hit. Error genes are ones which do not have self-hits.
# BLAST results are also filtered to remove these genes. A file <GENOME>.errors
# is also created that contains these error gene names.
#
# pjh Nov. 2007: Updated to use NFS-based cluster
#
# pjh Sep. 2009: Modified to make the use of MPI optionable. These changes
#                were suggested by Greg Schwendimann.
#
# pjh Sep. 2011: Modified to use MPI machinefile
#
# pjh Jan. 2014: Fixed a bug in that the script assumed that the
#                full path is specified for the MPI machinefile.
#
# pjh Aug. 2014: Now using Jamie's mpiBlast to do the BLAST-ing. This script
#                also now always uses MPI.
#
# pjh Dec. 2014: Now using blast 2.2.30, which uses makeblastdb rather than
#                formatdb.
#
# pjh Jan. 2015: Added support for optionally using blastn.
#
# pjh Jun. 2015: Adapt to changes to mpiBlast. In particular now must
#                explicitly provide the output format (6).
#
# pjh Jul. 2015: Adapt to changes to mpiBlast. In particular now must
#                explicitly provide the number of hits to keep (500).

use strict;
use warnings;

use POSIX;

if (@ARGV < 5)
{
  die "Usage: doPairwiseBlasts.pl [-n] sequenceDirectory blastDirectory " .
    "evalueThreshold numberOfProcessors machinefile " .
    "<list of genome names>\n";
}

my $useBlastn = 0;
my $sequenceExt = "proteins";

if ($ARGV[0] eq "-n")
{
  $sequenceExt = "nuc";
  $useBlastn = 1;
  shift @ARGV;
}

my $sequenceDirectory = shift @ARGV;
my $blastDirectory = shift @ARGV;
my $evalueThreshold = shift @ARGV;
my $numberOfProcessors = shift @ARGV;
my $machinefile = shift @ARGV;;
my $tempMachinefile;

# check to see if machinefile can be opened
open (MF, "<", $machinefile) or
  die "cannot open machinefile provided for MPI\n";
close (MF);

my @newGenomes = @ARGV;

# Read an errors file into a hash
sub readErrorsFile
{
  my $filename = $_[0];
  my %errorHash = ();
  open(XYZ, "<", $filename) or
    die "cannot open input ($filename)\n";
  while (my $line = <XYZ>)
  {
    chomp($line);

    $errorHash{$line} = $line;
  }
  close(XYZ);

  return %errorHash;
}

# trim trailing and leading spaces from a string
sub trimSpaces
{
  my $str = $_[0];

  if ($str =~ /^\s*([^\s]*)\s*$/)
  {
    $str = $1;
  }
  else
  {
    die("embedded space in input: $str");
  }

  return $str;
}

# Given the tabular output of an NCBI blast call, parse it into a list of
# (gene, bit score, e-value, alignment length) tuples. Return the list.
sub getHits
{
  my $arrayRef = shift;

  my @lines = @$arrayRef;

  # build a hash of hits to filter out multiple hits with same subject gene
  my %ret = ();

  # process each line of output
  foreach my $line (@lines)
  {
    chomp($line);

    my @data = split /\t/, $line;
    my $queryId = trimSpaces($data[0]);
    my $subjectId = trimSpaces($data[1]);
    my $alignLength = trimSpaces($data[3]);
    my $evalue = trimSpaces($data[10]);
    my $bitScore = trimSpaces($data[11]);

    # could be multiple hits for same subject; keep the best bit score
    my $oldStr = $ret{$subjectId};
    if (defined($oldStr))
    {
      my @oldData = split /!/, $oldStr;
      my $oldBitScore = $oldData[1];
      if ($bitScore > $oldBitScore)
      {
        $ret{$subjectId} = "$subjectId!$bitScore!$evalue!$alignLength";
      }
    }
    else
    {
      $ret{$subjectId} = "$subjectId!$bitScore!$evalue!$alignLength";
    }
  }

  return values %ret;
}


# This performs one BLAST operation between two genomes. The
# first genome is the set of query genes and the two second genome acts
# as the database to be searched. It produces two output files:
#   1. <genome>-<db>.blast: the raw blast results.
#   2. <genome>.errors: list of genes that do not have any hits.
#
sub doOneBlast
{
  my $genome = $_[0];
  my $db = $_[1];

  print "  Executing the BLASTs using MPI...\n";

  # format the db
  #system "formatdb -p -i $db.prepared";
  my $blastType;
  if ($useBlastn)
  {
    $blastType = "blastn";
    system "makeblastdb -dbtype nucl -in $db.prepared";
  }
  else
  {
    $blastType = "blastp";
    system "makeblastdb -dbtype prot -in $db.prepared";
  }

  if($? != 0)
  {
      die("Could not format database $db.prepared");
  }

  # run mpiBlast (which takes two extra processes (scheduler and writer)
  # keep up to 500 blast hits
  # use output format 6
  my $actualProcessCount = $numberOfProcessors + 2;
  system "mpiexec -n $actualProcessCount -f $tempMachinefile mpiBlast $blastType -query $genome.prepared -db $db.prepared -evalue $evalueThreshold -max_target_seqs 500 -outfmt 6 -out $genome-$db.temp 2>&1 </dev/null";

  if($? != 0)
  {
      die("mpiexec of mpiBlast failed for $genome.prepared $db.prepared");
  }

  # categorize the results as either "good" or "error"

  open(OUTPUT, ">", "$genome-$db.blast") or
    die("Could not open blast output file $genome-$db.blast");

  open(OUTPUT_ERRORS, ">", "$genome-$db.errors") or
    die("Could not open blast output file $genome-$db.errors");

  open (INPUT, "<", "$genome-$db.temp") or
    die("Could not open temp file $genome-$db.temp");

  my @outputLines = <INPUT>;
  close INPUT;

  my %sequenceHash = ();

  my @queryHits = ();

  foreach my $line (@outputLines)
  {
    chomp($line);

    my @data = split /\t/, $line;
    my $queryId = $data[0];

    if (defined($sequenceHash{$queryId}))
    {
      push(@queryHits, $line);
    }
    else
    {
      if (@queryHits > 0)
      {
        my @tempArray = split /\t/, $queryHits[0];

        my @hits = getHits(\@queryHits);

        if (@hits > 0)
        {
          print OUTPUT "$tempArray[0]";

          foreach my $hit (@hits)
          {
            print OUTPUT " $hit";
          }
          print OUTPUT "\n";

          @queryHits = ();
        }
        else
        {
        }
      }

      $sequenceHash{$queryId} = $queryId;
      push(@queryHits, $line);
    }
  }

  if(@queryHits > 0)
  {
    my @tempArray = split /\t/, $queryHits[0];

    my @hits = getHits(\@queryHits);

    if (@hits > 0)
    {
      print OUTPUT "$tempArray[0]";

      foreach my $hit (@hits)
      {
        print OUTPUT " $hit";
      }
      print OUTPUT "\n";

      @queryHits = ();
    }
  }

  close OUTPUT;
                                                              
  open (INPUT, "<", "$genome.prepared") or
    die("Could not open $genome.prepared");

  while(<INPUT>)
  {
    chomp;

    if (/^>/)
    {
      my $seq = substr $_, 1;

      if (!defined($sequenceHash{$seq}))
      {
        print OUTPUT_ERRORS "$seq\n";
      }
    }
  }

  close INPUT;

  # cleanup temp files
  system "rm *.prepared.* *.temp";

}

# prepare the pretein file for each new genome
print "Preparing sequences for new genomes...\n";
my %newGenomeHash = ();
foreach my $genome (@newGenomes)
{
  if (defined($newGenomeHash{$genome}))
  {
    die "duplicate genome ($genome) in list of new genomes!\n";
  }
  else
  {
    $newGenomeHash{$genome} = $genome;
  }
  my $err = system "prepareSequenceFile.pl $genome " .
    "$sequenceDirectory/$genome.$sequenceExt $blastDirectory/$genome.prepared";
  if ($err != 0)
  {
    die "prepareSequenceFile.pl failed for $genome\n";
  }
  else
  {
    print "  $genome.\n";
  }
}
print "  Done.\n";

# copy the MPI machinefile into the blast directory so that
# it will be accessible in case we don't have its full path
$tempMachinefile = "tmp-" . POSIX::getpid() . "-mf";
my $cpErr = system "cp $machinefile $blastDirectory/$tempMachinefile ";
if ($cpErr != 0)
{
  die "FAILED: copy of machinefile!\n";
}

# make the blast directory the working directory
chdir ("$blastDirectory") or
  die "FAILED: cd $blastDirectory\n";

# identify and validate previously processed genomes
print "Identifying and validating previously processed genomes...\n";
my @oldGenomes = ();
open(IN, "<", "DONE") or
  die "Cannot open file containing list of previously processed genomes.\n" .
      "  The file should be named \"DONE\".\n";
while (my $genome = <IN>)
{
  chomp($genome);
  push @oldGenomes, $genome;
  if (defined($newGenomeHash{$genome}))
  {
    die "genome ($genome) has already been processed?\n";
  }
}
close(IN);
foreach my $old (@oldGenomes)
{
  open(TMP, "<", "$old.self") or
    die "self-hit file for $old was not found!\n";
  close(TMP);
  open(TMP, "<", "$old.prepared") or
    die "prepared FASTA file for $old was not found!\n";
  close(TMP);
  open(TMP, "<", "$old.errors") or
    die "errors file for $old was not found!\n";
  close(TMP);

  foreach my $otherOld (@oldGenomes)
  {
    open(TMP, "<", "$old-$otherOld.blast") or
      die "blast results file $old-$otherOld.blast was not found!\n";
    close(TMP);
  }
  print "  $old.\n";
}
print "  Done.\n";

# Now loop through each new genome and BLAST it against itself and create the
# self-hit file and errors file.
foreach my $new (@newGenomes)
{
  # count the genes in this genome
  my $geneCount = 0;
  open(IN, "<", "$new.prepared") or
    die "cannot open input ($new.prepared)\n";
  while (my $line = <IN>)
  {
    chomp($line);

    if ($line =~ /^>/)
    {
      $geneCount += 1;
    }
  }
  close(IN);

  # BLAST new x new
  print "BLAST $new against itself...\n";
  doOneBlast($new, $new);
  print "  Done.\n";

  # Read the error file produced by this BLAST into a hash
  my %errorHash = readErrorsFile("$new-$new.errors");
  my $bogusGenes = keys %errorHash;

  # Find the self-hit and filter the BLAST results to remove hits to error genes
  print "Finding the self-hits and the error genes...\n";
  my $missingSelfHitCount = 0;
  open(IN, "<", "$new-$new.blast") or
    die "cannot open input ($new.blast)\n";
  open(OUT, ">", "$new.self") or
    die "cannot open output ($new.self)\n";
  while (my $line = <IN>)
  {
    chomp($line);

    # line contains the query gene first followed by the genes it hit
    my @piece = split / /, $line;

    # grab the query gene
    my $key = shift @piece;

    # sanity check: query gene should never be an error gene
    if (defined($errorHash{$key}))
    {
      die "query gene ($key) is in error list?\n";
    }

    # iterate through the hits to find self-hit
    my $foundIt = 0;
    for (my $i = 0; $i < @piece; $i += 1)
    {
      my @datum = split /!/, $piece[$i];
      if ($datum[0] eq $key)
      {
        print OUT "$key $datum[1]\n";
        $foundIt = 1;
      }
    }
    # if no self-hit, then add it to the error hash
    if ($foundIt == 0)
    {
      $missingSelfHitCount += 1;
      $errorHash{$key} = $key;
    }
  }
  close(IN);
  close(OUT);

  # Filter the BLAST results to remove error genes
  open(IN, "<", "$new-$new.blast") or
    die "cannot open input ($new-$new.blast)\n";
  open(OUT, ">", "$new-$new.tmp") or
    die "cannot open output ($new-$new.tmp)\n";
  while (my $line = <IN>)
  {
    chomp($line);

    # line contains the query gene first followed by the genes it hit
    my @piece = split / /, $line;

    # grab the query gene
    my $key = shift @piece;
    if (defined($errorHash{$key}))
    {
      next;  # skip if the query gene is on the error list (missing self-hit)
    }
    else
    {
      print OUT "$key";
    }

    # iterate through the hits to filter out errors
    for (my $i = 0; $i < @piece; $i += 1)
    {
      my @datum = split /!/, $piece[$i];
      if (defined($errorHash{$datum[0]}))
      {
        next;  # skip if the subject gene is on the error list
      }
      else
      {
        print OUT " $piece[$i]";
      }
    }
    print OUT "\n";
  }
  close(IN);
  close(OUT);

  # rename filtered tmp file to be the BLAST result file
  rename("$new-$new.tmp", "$new-$new.blast") or
    die "rename of $new-$new.tmp to $new-$new.blast failed!";

  # dump error hash to a file
  my $errorCount = 0;
  open(OUT, ">", "$new.errors") or
    die "cannot open output ($new.errors)\n";
  foreach my $gene (keys %errorHash)
  {
    print OUT "$gene\n";
    $errorCount += 1;
  }
  close(OUT);
  print "  $errorCount ($bogusGenes + $missingSelfHitCount) errors in " .
           "$geneCount genes ";
  my $percent = ($errorCount / $geneCount) * 100;
  printf "(%.1f%%)\n", $percent;

  # can remove the errors file produced by the BLAST
  unlink ("$new-$new.errors") or
    die "delete of $new-$new.errors failed!";

  print "  Done.\n";
}

# Now loop through each new genome and BLAST it against the other new genomes.
foreach my $new (@newGenomes)
{
  # read the new genome's error file into a hash
  my %errorHash = readErrorsFile("$new.errors");

  foreach my $otherNew (@newGenomes)
  {
    if ($new eq $otherNew)
    {
      next;
    }

    print "BLAST $new against $otherNew...\n";
    doOneBlast($new, $otherNew);
    print "  Done.\n";

    # the error file produced by this step can be discarded
    unlink "$new-$otherNew.errors" or
      die "cannot delete $new-$otherNew.errors";

    print "Filtering error genes from the BLAST results...\n";

    # read in error genes for other new genome
    my %otherErrorHash = readErrorsFile("$otherNew.errors");

    open(IN, "<", "$new-$otherNew.blast") or
      die "cannot open input ($new-$otherNew.blast)\n";
    open(OUT, ">", "$new-$otherNew.tmp") or
      die "cannot open output ($new-$otherNew.tmp)\n";
    while (my $line = <IN>)
    {
      chomp($line);

      # line contains the query gene first followed by the genes it hit
      my @piece = split / /, $line;

      # grab the query gene
      my $key = shift @piece;
      if (defined($errorHash{$key}))
      {
        next;  # skip if the query gene is on the error list
      }
      else
      {
        print OUT "$key";
      }

      # iterate through the hits to filter out errors
      for (my $i = 0; $i < @piece; $i += 1)
      {
        my @datum = split /!/, $piece[$i];
        if (defined($otherErrorHash{$datum[0]}))
        {
          next;  # skip if the subject gene is on the error list
        }
        else
        {
          print OUT " $piece[$i]";
        }
      }
      print OUT "\n";
    }
    close(IN);
    close(OUT);

    # rename filtered tmp file to be the BLAST result file
    rename ("$new-$otherNew.tmp", "$new-$otherNew.blast") or
      die "rename of $new-$otherNew.tmp to $new-$otherNew.blast failed!";

    print "  Done.\n";
  }
}

# Now loop through each new genome and BLAST it against each old genome
foreach my $new (@newGenomes)
{
  # read in error genes for new genome
  my %errorHash = readErrorsFile("$new.errors");

  foreach my $old (@oldGenomes)
  {
    print "BLAST $new against $old...\n";
    doOneBlast($new, $old);
    print "  Done.\n";

    # the error file produced by this step can be discarded
    unlink "$new-$old.errors" or
      die "cannot delete $new-$old.errors";

    print "Filtering error genes from the BLAST results...\n";

    # read in error genes for old genome
    my %otherErrorHash = readErrorsFile("$old.errors");

    open(IN, "<", "$new-$old.blast") or
      die "cannot open input ($new-$old.blast)\n";
    open(OUT, ">", "$new-$old.tmp") or
      die "cannot open output ($new-$old.tmp)\n";
    while (my $line = <IN>)
    {
      chomp($line);

      # line contains the query gene first followed by the genes it hit
      my @piece = split / /, $line;

      # grab the query gene
      my $key = shift @piece;
      if (defined($errorHash{$key}))
      {
        next;  # skip if the query gene is on the error list
      }
      else
      {
        print OUT "$key";
      }

      # iterate through the hits to filter out errors
      for (my $i = 0; $i < @piece; $i += 1)
      {
        my @datum = split /!/, $piece[$i];
        if (defined($otherErrorHash{$datum[0]}))
        {
          next;  # skip if the subject gene is on the error list
        }
        else
        {
          print OUT " $piece[$i]";
        }
      }
      print OUT "\n";
    }
    close(IN);
    close(OUT);

    # rename filtered tmp file to be the BLAST result file
    rename ("$new-$old.tmp", "$new-$old.blast") or
      die "rename of $new-$old.tmp to $new-$old.blast failed!";

    print "  Done.\n";
  }
}

# Now loop through each new genome and BLAST each old genome against it
foreach my $new (@newGenomes)
{
  # read in error genes for new genome
  my %errorHash = readErrorsFile("$new.errors");

  foreach my $old (@oldGenomes)
  {
    print "BLAST $old against $new...\n";
    doOneBlast($old, $new);
    print "  Done.\n";

    # the error file produced by this step can be discarded
    unlink "$old-$new.errors" or
      die "cannot delete $old-$new.errors";

    print "Filtering error genes from the BLAST results...\n";

    # read in error genes for other new genome
    my %otherErrorHash = readErrorsFile("$old.errors");

    open(IN, "<", "$old-$new.blast") or
      die "cannot open input ($old-$new.blast)\n";
    open(OUT, ">", "$old-$new.tmp") or
      die "cannot open output ($old-$new.tmp)\n";
    while (my $line = <IN>)
    {
      chomp($line);

      # line contains the query gene first followed by the genes it hit
      my @piece = split / /, $line;

      # grab the query gene
      my $key = shift @piece;
      if (defined($otherErrorHash{$key}))
      {
        next;  # skip if the query gene is on the error list
      }
      else
      {
        print OUT "$key";
      }

      # iterate through the hits to filter out errors
      for (my $i = 0; $i < @piece; $i += 1)
      {
        my @datum = split /!/, $piece[$i];
        if (defined($errorHash{$datum[0]}))
        {
          next;  # skip if the subject gene is on the error list
        }
        else
        {
          print OUT " $piece[$i]";
        }
      }
      print OUT "\n";
    }
    close(IN);
    close(OUT);

    # rename filtered tmp file to be the BLAST result file
    rename ("$old-$new.tmp", "$old-$new.blast") or
      die "rename of $old-$new.tmp to $old-$new.blast failed!";

    print "  Done.\n";
  }
}

# Finally, update the DONE file
open(OUT, ">>", "DONE") or
  die "Cannot open DONE for appending!\n";
foreach my $new (@newGenomes)
{
  print OUT "$new\n";
}
close(OUT);

# clean up copy of the machinefile
$tempMachinefile = "tmp-" . getpid() . "-mf";
my $rmErr = system "rm $tempMachinefile";
if ($rmErr != 0)
{
  die "FAILED: remove of temporary machinefile!\n";
}

print "\nProcess Complete.\n";

