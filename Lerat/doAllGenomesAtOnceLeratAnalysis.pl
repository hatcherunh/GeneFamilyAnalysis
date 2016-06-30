#!/usr/bin/perl

# $Id: doAllGenomesAtOnceLeratAnalysis.pl 57 2015-07-22 21:41:52Z pjh $
#
#
# Phil Hatcher, July 2008
#
# Do Lerat genome homolog analysis for a set of genomes
# on which pairwise BLASTs have already been done.
#
# The analysis performed for all genomes includes:
#   1. identify high-quality hits --> <prefix>.hits (according to some Lerat
#        or e-value threshold)
#   2. reverse the high-quality hits --> <prefix>.reverse
#   3. identify homolog families --> <prefix>.family 
#   4. identify orthologs --> <prefix>.orthologs
#   5. identify panorthologs --> <prefix>.panorthologs
#   6. identify paralog families (ie families that have at least one genome
#        with more than one gene present in the family) --> <prefix>.paralogs
#   7. identify unique genes --> <prefix>.unique
#   8. basic statistics (e.g. largest family size) --> <prefix>.stats
#   9. summary (genes per genome) for each paralog family -->
#        <prefix>.paralogs-summary
#
# It takes five initial arguments:
#   1. -lerat or -evalue
#   2. either evalue threshold or the lerat ratio threshold
#   3. -oneway or -reciprocal
#   4. name of directory that contains the pairwise BLAST results
#   5. prefix to use for output files and as directory name that will
#        be created for storing these files
#
# These arguments are followed by a list of genome names that define the
# set of genomes being analyzed. This list must contain at least two genomes.
#
# Instead of a list of genomes, -all can be specified. In this case the DONE
# file in the blast directory is consulted to get the list of genomes.
#
# Modified by pjh in June 2010 by request of Nancy Garnhart to make the output
# files more amenable for further processing.
#

use strict;
use warnings;
use Cwd;
use Sys::Hostname;

if (@ARGV < 6)
{
  die "Usage: doAllGenomesAtOnceLeratAnalysis.pl [-lerat | -evalue] " .
    "threshold [-oneway | -reciprocal] blastDirectory prefix " .
    "<list of genome names>\n";
}

my $technique = shift @ARGV;
my $threshold = shift @ARGV;
my $membership = shift @ARGV;
my $blastDirectory = shift @ARGV;
my $prefix = shift @ARGV;

my @genomes = @ARGV;

if ($genomes[0] eq "-all")
{
  if (@genomes > 1)
  {
    die "-all should be the last argument!";
  }

  # discard -all
  shift @genomes;

  open(IN, "<", "$blastDirectory/DONE") or
    die "cannot open input ($blastDirectory/DONE)\n";

  while (my $line = <IN>)
  {
    chomp($line);
    push @genomes, $line;
  }
  close (IN);
}
else
{
  if (@genomes < 2)
  {
    die "must provide at least two genome names!\n";
  }
}

if ($technique ne "-lerat" && $technique ne "-evalue")
{
  die "first argument must be either -lerat or -evalue\n";
}

if ($membership ne "-oneway" && $membership ne "-reciprocal")
{
  die "third argument must be either -oneway or -reciprocal\n";
}

my $exit;

# create directory to contain result files
print "Create directory (./$prefix) to contain results\n";
$exit = system "mkdir $prefix";
if ($exit == 0)
{ 
  print "  Done.\n";
}
else
{
  print "  Failed with exit code $exit. Aborting....\n";
  die "";
}

# now update the prefix variable to include the directory
$prefix = $prefix . "/" . $prefix;

# create and open the stats file
open(STATS, ">", $prefix . ".stats") or
  die "cannot open output ($prefix.stats)\n";

# get and format the local time
my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek,
  $dayOfYear, $daylightSavings) = localtime();
my $year = 1900 + $yearOffset;
my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] " .
              "$dayOfMonth, $year";

# get the hostname
my $host = hostname;

# get the current working directory
my $dir = getcwd;

# get the login running this script
my $login = getlogin;
if (!defined($login))
{
  $login = "PBS";
}

# write output basic info about the run to the stats file
print STATS "date: $theTime\n";
print STATS "host: $host\n";
print STATS "login: $login\n";
print STATS "technique: $technique all-at-once\n";
print STATS "threshold: $threshold\n";
print STATS "membership: $membership\n";
print STATS "cwd: $dir\n";
print STATS "blast results: $blastDirectory\n";

my $genomeString = join " ", @genomes;
print "Genomes being processed: $genomeString\n";
print STATS "genomes: $genomeString\n";

# close the stats file as it will be appended to by other scripts
close(STATS);

# validate that all the genomes have had their BLAST processing done
print "Validating that all genomes have had their BLAST processing done...\n";
foreach my $genome (@genomes)
{
  open(TMP, "<", "$blastDirectory/$genome.self") or
    die "self-hit file for $genome was not found!\n";
  close(TMP);

  foreach my $otherGenome (@genomes)
  {
    open(TMP, "<", "$blastDirectory/$genome-$otherGenome.blast") or
      die "blast results file $genome-$otherGenome.blast was not found!\n";
    close(TMP);
  }
  print "    $genome.\n";
}
print "  Done.\n";

# get the high-quality hits
print "Get the high-quality hits ($technique $threshold)...\n";
$exit = system "getHighQualityHits.pl $blastDirectory $technique " .
                 "$threshold $prefix.hits $genomeString";
if ($exit == 0)
{ 
  print "  Done.\n";
}
else
{
  print "  Failed with exit code $exit. Aborting....\n";
  die "";
}

# reverse the high-quality hits
print "Reverse the high-quality hits...\n";
$exit = system "getReverseHitsJJ.pl $prefix.hits $prefix.reverse";
if ($exit == 0)
{ 
  print "  Done.\n";
}
else
{
  print "  Failed with exit code $exit. Aborting....\n";
  die "";
}

# find the homolog families
print "Find the homolog families...\n";
if ($membership eq "-oneway")
{
  $exit = system "findHomologFamilies.pl $prefix.hits $membership " .
                   "$prefix.family";
}
else
{
  $exit = system "findHomologFamilies.pl $prefix.hits $membership " .
                   "$prefix.family $prefix.reverse";
}
if ($exit == 0)
{ 
  print "  Done.\n";
}
else
{
  print "  Failed with exit code $exit. Aborting....\n";
  die "";
}

# find the unique genes
print "Find the unique genes...\n";
$exit = system "findUniques.pl $prefix $membership";
if ($exit == 0)
{ 
  print "  Done.\n";
}
else
{
  print "  Failed with exit code $exit. Aborting....\n";
  die "";
}

# analyze the families
print "Analyze the families...\n";
$exit = system "analyzeFamilies.pl $prefix.family $prefix";
if ($exit == 0)
{ 
  print "  Done.\n";
}
else
{
  print "  Failed with exit code $exit. Aborting....\n";
  die "";
}

# find the panorthologs
print "Find the panorthologs...\n";
$exit = system "findMaximalPanorthologFamilies.pl $prefix.orthologs " .
  "$prefix.panorthologs $genomeString";
if ($exit == 0)
{ 
  print "  Done.\n";
}
else
{
  print "  Failed with exit code $exit. Aborting....\n";
  die "";
}

print "\nProcess Complete.\n";

