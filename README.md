# BashRegion
This function takes as input an Interval-format file and output all sgRNA
sites in the interval in both CSV and FASTA format.

Author: Neville Sanjana (nsanjana@nygenome.org)

CHANGELOG:
--------------------
* December 2015
	- Support for 5' PAMs (e.g. Cpf1)
	
* June 2015
	- Consolidated modules
	- Prepared for public release


* November 2014
	- Initial version
	
USAGE REQUIREMENTS
---------------------
Python 2.7 with the following libraries installed:
* Default Python libraries (sys, os, csv, urllib2, math, ElementTree)
* BioConductor (for SeqIO, SeqUtils, Seq, IUPAC)
* numpy
* pandas

Please make sure all of these libraries are installed and working before proceeding.

To run, place the BashRegion.py in a directory that also includes the input file 
(specified in intervalFile_INPUT). Then run the following command
> ``python BashRegion.py``





INPUT/PARAMETERS TO SET:
-------------------------
1) ``intervalFile_INPUT`` : A tab-delimited file (such as the included "SampleInput.interval" file) with a list of 
coordinates. See comment in the PARAMETERS section for an example of the format.

2) ``referenceGenomeForDAS`` : The reference genome to use for translating the genomic coordinates in ``intervalFile_INPUT``
to genomic DNA sequence. For the included file ("SampleInput.interval"), make sure to 
set ``referenceGenomeForDAS = "hg19"`` since the coordinates are from the hg19 human
reference genome.

3) ``outDirectory`` : The name of the directory for output files. It will be created if not already present.

4) ``outShortName`` : A short string to prefix the output files with.

5) ``cutoffSpacing`` : The minimum spacing between selected sgRNAs on the same strand. For example, if 
``cutoffSpacing = 2``, then the next sgRNA must be > 2 bp from the previous sgRNA chosen. To
find all sgRNAs in a region, set ``cutoffSpacing = 0``.

6) ``PAM`` : A string that contains the PAM sequence (e.g. ``"NGG"`` for SpCas9)


FORMAT OF **intervalFile_INPUT**:
-------------------------
The interval file should have 5 columns and be tab-delimited

``<chromosome> <strand (optional but needs column to be there)> <base start> <base end> <geneID>``

Example of **intervalFile_INPUT**:
>``chr12	.	51057788	51173928	ATF1``

>``chr1		+	8821058		8921418		ENO1``

Note that the strand parameter is not used by BashRegion.py but placed here for compatibility with standard interval file formats.



OUTPUT FILES:
-------------------------
Output is all sgRNAs on both strands in CSV and FASTA format. 

The CSV format output includes the sgRNA sequences and the following parameters:
* Spacer number
* Strand
* Start Pos
* End Pos
* Chromosome
* ChrStart
* ChrEnd
* Cut Site
* Distance from Previous Spacer
* Spacer (i.e. sgRNA)
* PAM
* Spacer GC
* Seq ID


