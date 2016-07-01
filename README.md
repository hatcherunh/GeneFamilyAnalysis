# Gene Family Analysis

This is a collection of code to compute gene families.

The foundation of the analysis is an efficient mechanism for executing
all-by-all BLASTS between the genes of a collection of genomes.

The family analysis follows the approach described by Lerat et al.,
"From Gene Trees to Organismal Phylogeny in Prokaryotes: The
Case of the γ-Proteobacteria", *PLoS Biology*, 2004,
which identifies homologs as those gene pairs that had BLAST hits
in both directions within a given scaled bit score threshold.

There are also scripts to estimate the evolutionary rates of
the families.

For an example of using these scripts, see Cooper et al.,
"Why Genes Evolve Faster on Secondary Chromosomes in Bacteria",
*PLoS Computational Biology*, 2010.

INSTALLATION
--

**Notes**

- The code has been designed for, and only tested on, Linux systems.

- The code is organized into four directories:
 - *blast*: perform all-by-all BLASTs.
 - *Lerat*: compute gene families.
 - *Ka-Ks*: estimate evolutionary rates.
 - *misc*: other miscellaneous scripts.

**Dependencies**

- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (tested with 2.2.30)

- [MPICH](http://www.mpich.org/downloads/) (tested with 3.1.1)

- [PHYLIP](http://evolution.genetics.washington.edu/phylip/getme.html) (tested with 3.66)

- [ClustalW2](http://www.clustal.org/clustal2/) (tested with 2.1)

- [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html#download) (tested with 4.7a)

**To Do**

1. Install the above dependencies.

2. Include in your PATH the directories containing the executables for the
above dependencies.

3. Download all the Perl and bash scripts and place them into a directory
that is in your PATH.
Add execute permission to all the scripts.

4. Download the single C file (blast/mpiBlast.c), compile it using mpicc
from MPICH (```mpicc -o mpiBlast mpiBlast.c```), and place the executable
in a directory that is in your PATH.
Add execute permission to this file, if necessary.

USER GUIDE
--

**Input Files**

The ortholog detection algorithms use the amino acid sequences (proteins)
for your genomes.
If you plan to compute Ka and Ks, then you will also need the nucleotide
sequences for your genes.
And you should verify that these two sequences match for each gene.

The sequences should be in FASTA format, with the sequences for each
genome being in its own file.
The file name for each amino acid sequence file should be *[genome].proteins*,
where *[genome]* is the name of the genome.
The file name for each nucleotide sequence file should be *[genome].nuc*,
where *[genome]* is the name of the genome.
The filenames for the two sequence files for a genome should use the
same name for the genome.

Each FASTA header should start with a unique identifier for the gene.
This identifier is assumed to be what is in between the *>* and the
first space.
Text after the first space is ignored.
The identifier for a gene should be the same in the amino acid file
and in the nucleotide file.

Place all the amino acid sequence files in a directory called
*proteins".
Place all the nucleotide sequence files in a directory called
*nuc".

**Doing the BLASTs**

MPI (as implemented by MPICH) is the mechanism used to run the
BLAST analysis in parallel.
You must create an MPI *machine file* to indicate how many
processes MPI should utilize when executing the BLAST analysis.
The typical user will want to use a single machine with multiple cores,
which would require a one-line machine file, which lists the
name of the machine and the number of processes (e.g. *compserver:60*).

Create a directory for the BLAST results.
Put an empty file called *DONE* into this directory.
(This file will track which genomes have been processed.)

Run *doPairwiseBlasts.pl*.
It takes five initial arguments:

1. directory containing protein FASTA files.

2. directory to contain BLAST results.

3. evalue threshold for which hits to keep.

4. number of processes MPI should use.

5. machine file MPI should use.

For the remaining arguments, list the genomes you want to process.
The names you list as arguments must match the names you used for
the amino acid files.

This step will likely take some time if you are running with a large
set of genomes, since it performs all-by-all blasts for all pairs
of genomes.

The BLASTs can be done incrementally.
You can run *doPairwiseBlasts.pl* for some of the genomes, and
then run it again later to add more genomes to the mix.
For the second run, you only need to list the new genomes.

**Compute Gene Families**

Run *doAllGenomesAtOnceLeratAnalysis.pl*.
It takes five initial arguments:

1. "-lerat" (the alternative is "-evalue" which will select homologs based
upon BLAST hits below a certain evalue threshold, rather than using
scaled bitscores).

2. homolog selection threshold (for scaled bitscores, a value in the range
0-1; for "-evalue", an evalue threshold).

3. "-reciprocal" (requires the gene hits to satisfy the threshold requirement
in both directions; the alternative is "-oneway", which requires the
threshold to be satisfied in only one direction).

4. the directory containing the BLAST results from *doPairwiseBlasts.pl*.

5. a prefix to be used for the output files (see below for a description
of the output files), and which is also used as the directory name for a
directory that is created and where the output files will be placed.

After the fifth argument, either list the genomes you want to analyze
or say "-all" to indicate you want all the genomes analyzed for which
there is BLAST data.

The bit scores are scaled by the bit score of the self hit of the query gene.
Gene families are formed by including two genes in a family if they had
been identified as homologs.

The script will produce a number of output files:
- *[prefix]*.panorthologs: These are the gene families that have exactly
one member from each genome.
- *[prefix]*.orthologs: These are the gene families that have at most
one member from each genome.
- *[prefix]*.paralog: These are the gene families with at least two members
from one genome.
- *[prefix]*.family: All gene families.
- *[prefix]*.unique: The genes that did not end up in any family.
- *[prefix]*.paralog-summary: Counts of genes per genome for paralog
families.
- *[prefix]*.stats: Interesting statistics about the families.
- *[prefix]*.hits and [prefix].reverse: BLAST hit data used to generate the
families. These are not usually useful.

The files containing gene families have one family per line, with each
line beginning with a unique numeric family identifier.
The genes in the family are represented as *[genome]$[geneID]*.

**Estimate Evolutionary Rates**

Run *doKaKsAnalysis.sh*.
This script expects three arguments:

1. prefix of the family file to be processed.

2. extension of the family file to be processed.

3. the genome order file.

Therefore, the family file to be processed is *[prefix].[extension]*.

The genome order file is a file that contains the names of the genomes,
one genome per line.
 
The script assumes that the amino acid sequence files are available
at "../proteins" and that the nucleotide sequence files are available
at "../nuc".

The script performs the following steps:

1. The amino acid sequences of each family are aligned using ClustalW2.

2. The codon boundaries are used to align the nucleotide sequences.

3. The leading and trailing edges of each amino acid sequence in every family
is trimmed to generate consensus edges, and then the nucleotide sequences are
trimmed to match.

4. From this trimmed file, a consensus sequence for the family is found,
using the *cons* utility from the EMBOSS suite.

5. Each sequence in the family is compared against the consensus sequence
and the maximum number of amino acid differences for a gene in the family is
computed.

6. Phylogenetic trees are constructed for each family using DNAML
(maximum likelihood) in PHYLIP using default settings and the Newick
formatted trees are saved.

7. dN and dS are calculated from the trimmed nucleic acid alignment and
the DNAML tree as a guide using codeml in the PAML package].
Codeml model 0, which allows for a single dN and dS value throughout
the phylogeny, is used.

8. A CSV (comma-separated value) file is written containing the following
fields for each family: family identifier, maximum number of amino acid
differences from the consensus sequence, dN and dS.
The file is named *[prefix].csv*, where *[prefix]* is the first
argument supplied to the script.
(Numerous intermediate directories and files are also created, which
can be ignored.)

**Note:** The script *doKaKsAnalysis.sh* is designed to be run with
families of size four or bigger.

**Other Scripts**

1. *annotateFamilies2.pl*: Reads a family file and adds annotation for
each family.
The annotation is obtained from the FASTA header of the amino acid
sequence file.
A CSV file is output.
See the comments at the top of the script for directions on how to
run it.

2. *createPhylipParsInput.pl*: Converts an ortholog family file to
an input file for PHYLIP's PARS program.
See the comments at the top of the script for directions on how to
run it.

SAMPLE RUN
--

The *sample-run* directory contains input data and output data from
a sample run.
The *proteins* and *nuc* subdirectories contain the amino acid and
nucleotide sequences, respectively, for four bacteria genomes.
The *genome-order* file contains the names of the four genomes.
The *genome-abbrev* file contains the abbreviated (10 or fewer
characters) names of the four genomes, required by PARS.
These are the necesary input files.

The following sequence of commands was executed:
```
mkdir blast
cd blast
touch DONE
cd ..
doPairwiseBlasts.pl proteins blast 1.0 4 ~/mf.c4 Aliivibrio-salmonicida-LFI1238 Photobacterium-profundum-SS9 Vibrio-cholerae-2740-80 Vibrio-cholerae-M66-2
doAllGenomesAtOnceLeratAnalysis.pl -lerat .7 -reciprocal blast point7 -all
mkdir KaKs-point7-panorthologs
cd KaKs-point7-panorthologs
cp ../point7/point7.panorthologs .
doKaKsAnalysis.sh point7 panorthologs ../genome-order
cd ..
mkdir annotate
cd annotate
annotateFamilies2.pl ../point7/point7.panorthologs ../genome-order ../proteins annotated-point7-panorthologs
cd ..
mkdir PARS-input
cd PARS-input
createPhylipParsInput.pl ../point7/point7.orthologs ../genome-order ../genome-abbrev point7-orthologs-PARS-input
cd ..
```

This generates the results of the family analysis on the four genomes using
a .7 scaled bitscale threshold into the *point7* subdirectory.

The evolutionary rate estimates for the .7 panortholog families
are generated into the *KaKs-point7-panorthologs* subdirectory.

The annotated .7 panortholog families are generated into the
*annotate* subdirectory.

A PARS input for the .7 ortholog families is generated into the
*PARS-input* subdirectory.

**Note**: Numerous generated intermediate files and subdirectories have
been deleted.

CONTRIBUTORS
--

The code was written by Phil Hatcher, Jamie Jackson and Sam Vohr.
Key design ideas were provided by Kelley Thomas and Vaughn Cooper.
Phil Hatcher (hatcher@unh.edu) maintains the code.

