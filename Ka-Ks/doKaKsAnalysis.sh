#!/bin/sh -v

# $Id: doDeadMicePipe2.sh 46 2015-06-02 00:09:58Z pjh $
#
# use this version if you do not have "details" files.
#
# Phil Hatcher, December 2013
# reordered steps to replace gene names before running phylip
#
# Phil Hatcher, November 2013
# edits to add a second argument to give the file extension
# for the family file and a third argument to specify the genome
# order file.
#
# so this script takes three arguments:
# 1. prefix of the family file name
# 2. extension of the family file name
# 3. the taxon order file
#
# This script still expects the "nuc" and "proteins" directory to
# be "up one level" at ../nuc and ../proteins.

makeFastaByFamily2.pl $1.$2 $3 nuc ../nuc $1
makeFastaByFamily2.pl $1.$2 $3 proteins ../proteins $1
alignAllFamilies.pl fasta $1-AA
alignAllNAbyAAFamilies.pl $1-NA $1-AA-aligned
trimAlignedEdges.pl $1-AA-aligned $1-NA-aligned
consensusAllFamilies.pl $1-AA-aligned-edged
maxDiffFromConsensus.pl $1-AA-aligned-edged $1-AA-aligned-edged-consensus $1-maxdiff
replaceGeneNames2.pl nuc $1-NA-aligned-edged
phylipFormatAllFamilies3.pl $1-NA-aligned-edged-renamed
codemlAllFamilies.pl $1-NA-aligned-edged-renamed $1-NA-aligned-edged-renamed-trees
codeml2csvAllInfo3.pl $1-NA-aligned-edged-renamed-codeml $1.$2 $3 $1-maxdiff >$1.csv
