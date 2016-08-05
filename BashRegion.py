# BashRegion.py

'''
This function takes as input an Interval-format file and output all sgRNA
sites in the interval in both CSV and FASTA format
Author: Neville Sanjana (Zhang Lab)

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
	Default Python libraries (sys, os, csv, urllib2, math, ElementTree)
	BioConductor (for SeqIO, SeqUtils, Seq, IUPAC)
	numpy
	pandas

Please make sure all of these libraries are installed and working before proceeding.

To run, place the BashRegion.py in a directory that also includes the input file 
(specified in intervalFile_INPUT). Then run the following command
>python BashRegion.py





INPUT/PARAMETERS TO SET:
-------------------------
1) intervalFile_INPUT 
A tab-delimited file (such as the included "SampleInput.interval" file) with a list of 
coordinates. See comment in the PARAMETERS section for an example of the format.

2) referenceGenomeForDAS
The reference genome to use for translating the genomic coordinates in intervalFile_INPUT 
to genomic DNA sequence. For the included file ("SampleInput.interval"), make sure to 
set referenceGenomeForDAS = "hg19" since the coordinates are from the hg19 human
reference genome.

3) outDirectory
The name of the directory for output files. It will be created if not already present.

4) outShortName
A short string to prefix the output files with.

5) cutoffSpacing
The minimum spacing between selected sgRNAs on the same strand. For example, if 
cutoffSpacing = 2, then the next sgRNA must be > 2 bp from the previous sgRNA chosen. To
find all sgRNAs in a region, set cutoffSpacing = 0.

6) PAM
A string that contains the PAM sequence (e.g. "NGG" for SpCas9)





OUTPUT FILES:
-------------------------
Output is all sgRNAs on both strands in CSV and FASTA format. 

The CSV format includes the sgRNA sequences and the following parameters:
	Spacer number
	Strand
	Start Pos
	End Pos
	Chromosome
	ChrStart
	ChrEnd
	Cut Site
	Distance from Previous Spacer
	Spacer (i.e. sgRNA)
	PAM
	Spacer GC
	Seq ID

'''

#####################
#####################
#### START PARAMS ###
#####################
#####################

# intervalFile_INPUT specifies the input file name
intervalFile_INPUT = "NONessentialGene_interval.txt"
'''
 The interval file should have 5 columns and be tab-delimited
<chromosome> <strand (optional but needs column to be there)> <base start> <base end> <geneID>

Example:
chr12	.	51057788	51173928	ATF1
chr1	+	8821058	8921418	ENO1

Note that the strand parameter is not used but placed here for compatibility with standard interval file formats.

'''

referenceGenomeForDAS = "hg19"
#referenceGenomeForDAS = "mm9"
#referenceGenomeForDAS = "hg38"
#referenceGenomeForDAS = "mm10"

# outDirectory is the name of the directory for output files
# outDirectory will be created if it doesn't already exist
# Make sure that outDirectory doesn't contain spaces
outDirectory = "BashRegion_output";  # no trailing slash here

# outShortName is the prefix for the output files (CSV and FASTA)
# Make sure that outShortName doesn't contain spaces
outShortName = "cpf1_NONessential_sgRNA"

# cutoffSpacing is the minimum distance allowed between 2 sgRNAs
# set cutoffSpacing to 0 if you want to output all sgRNAs in the region
cutoffSpacing=0

# PAM is the sequence of the CRISPR PAM to look for
# All IUPAC ambiguous bases are allowed in addition to ATCG
# For example, for S. pyogenes Cas9 (SpCas9), use PAM="NGG"
# For example, for S. aureus Cas9 (SaCas9), use PAM="NNGRRT"
# For example, for AsCpf1, use PAM="TTTN"
PAM="TTTN"

# PAMside indicates if PAM is on the 5' or 3' side
# For Cas9, use PAMside = 3
# For Cpf1, use PAMside = 5
PAMside = 5;

# spacerLength is the length of the PAM-adjacent region to select for each target
# For example, for S. pyogenes Cas9 (SpCas9), use spacerLength=20
# For example, for S. aureus Cas9 (SaCas9), use spacerLength=21 or 22
# For example, for AsCpf1, use spacerLength=20
spacerLength=20

# distanceToCutSiteFromPAM_bp indicates how many bases from the PAM 
# the Cas9 cuts. For SpCas9 and SaCas9, this value is 3. For Cpf1, this value is 19.
# (Technically, for Cpf1, it is 19bp on targeted strand and 23 bp on non-targeted strand)
# Direction does not matter here since the program already knows the PAMside.
distanceToCutSiteFromPAM_bp = 19

#####################
#####################
###### END PARAMS ###
#####################
#####################

# Libraries to import
import csv
import math
import os.path
import sys
import urllib2
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from xml.etree import ElementTree


'''

FUNCTION: GalaxyIntervalParser 

'''
def GalaxyIntervalParser(inputFile_galaxy):
	chrom=[]
	endRef=[]
	startRef=[]
	assocRefSeq=[]
		
	if not os.path.isfile(inputFile_galaxy):
		sys.exit("BashRegion: Input file not found! Please make sure the input file exists and try again.")
	
	with open(inputFile_galaxy) as f:
		reader=csv.reader(f,delimiter='\t')
		for row in reader:
			chrom.append(row[0])
			startRef.append(row[2])
			endRef.append(row[3])
			assocRefSeq.append(row[4])
	
	return {'chrom':chrom, 'startRefPos':startRef,'endRefPos':endRef,'assocRefSeq':assocRefSeq }		



'''

FUNCTION: spacerCSV2fa 

'''
def spacerCSV2fa(csv_file, outFile):

	if not os.path.isfile(csv_file):
		sys.exit(sys._getframe().f_code.co_name+": Input file not found!!")

	df = pd.read_csv(csv_file)
	# saved_column = df.column_name #you can also use df['column_name']

	# Create unique ID on the fly from columns of CSV file
	# fasta_ids = [">%s_%s_%s_%s" % t for t in zip(df['Seq ID'], df['Chromosome'],df['Strand'],df['ChrStart'])];
	
	# Use single column from CSV file as unique ID
	fasta_ids = [">%s" % t for t in zip(df['Seq ID'])];
	#print fasta_ids[1]
	fasta_seqs = df['Spacer']
	rowsFASTA = ["%s\n%s" % t for t in zip(fasta_ids, fasta_seqs)]
	
	with open(outFile, 'w') as thefile:
		for item in rowsFASTA:
  			thefile.write("%s\n" % item)



'''

FUNCTION: WriteColumns2CSV 

'''
def WriteColumns2CSV(listOfRows,outputFile_csv, newCSV, header):
		
	if newCSV:
		f=open(outputFile_csv,'w')
		csvOut=csv.writer(f)
		csvOut.writerow(header)    	
	else:
		f=open(outputFile_csv,'a')
		csvOut=csv.writer(f)

	for idx, item in enumerate (listOfRows):
		csvOut.writerow(item)
	f.close()



'''

FUNCTION: coords2fa 

'''	
def coords2fa(chromosome, start, end, referenceGenomeForDAS):
	base = 'http://genome.ucsc.edu/cgi-bin/das/' + referenceGenomeForDAS +'/dna?segment='
	url = base + chromosome + ':' + str(start) + ',' + str(end)
	# Check on proper URL formation
	#print "as={ur}".format(ur=url)
	request = urllib2.Request(url, headers={"Accept" : "application/xml"})
	u = urllib2.urlopen(request)
	tree = ElementTree.parse(u)
	doc = tree.getroot()
	
	dict={}
	for i in doc.iter():
		if i.text!=None:
			dict[i.tag]=i.text.strip()
		else:
			dict[i.tag]=""
 	#print dict['DNA']
 	# Sequence from UCSC has newline characters to remove
 	sequence = dict['DNA'].replace('\n','')
 	 
 	# print out key-value pairs
	#for key in sorted(dict.keys(), key=lambda v: v.upper()):
	#	print key+":"+dict[key]


	if sequence == '':
		sequence = 'THE SEQUENCE DOES NOT EXIST FOR GIVEN COORDINATES'
	#print "here you go: {myseq}".format(myseq=sequence)
	return sequence


	
'''

FUNCTION: FindSpacersInterval 

'''	
def FindSpacersInterval(chromPos, chromStartRG, chromEndRG, seqStr, cutoff_spacing, referenceGenomeForDAS,spacerLength, distanceToCutSiteFromPAM_bp):
	from Bio import SeqFeature	

	if PAMside == 3:
		distanceToCutSiteFrom5pEnd = spacerLength - distanceToCutSiteFromPAM_bp
		# For SpCas9 (as an example): spacerLength = 20bp ; distanceToCutSiteFromPAM_bp = 3bp; distanceToCutSiteFrom5pEnd = 17bp
	else:
		distanceToCutSiteFrom5pEnd = distanceToCutSiteFromPAM_bp-1
		# For AsCpf1 (as an example): spacerLength = 20bp ; distanceToCutSiteFromPAM_bp = 19bp; distanceToCutSiteFrom5pEnd = 18bp		
		
	s = coords2fa(chromPos, chromStartRG, chromEndRG, referenceGenomeForDAS);	
	
	s=s.upper();
	PAM = Seq(seqStr, IUPAC.ambiguous_dna)
	PAM_length = len(seqStr);
	if seqStr == str(PAM.reverse_complement()):
		DoRevComp=0
		forwardNameString = "{name}_{num:0{width}}"
	else:
		DoRevComp=1
		forwardNameString = "{name}_F{num:0{width}}"
	listSpacer=[]
	listDistBetweenSpacers=[]
	
	spacerNum=0
	prevStartLocInRefSeq=-9999
	if PAMside == 3:
		gbStringForSearch = s[spacerLength:];	# Cas9
	else:
		gbStringForSearch = s[:-spacerLength];   # Cpf1, get all but last ~20 bases of sequence
				
	spacerInds = SeqUtils.nt_search(gbStringForSearch, str(PAM))
	if len(spacerInds) > 1:	# matches found 
		del spacerInds[0] # first result from nt_search is regexp expansion
		#print "len line below {fname}".format(fname=len(spacerInds))
		formatDigitsN = int(math.ceil(math.log(len(spacerInds),10)));
		print "Plus strand sgRNAs found: {num}".format(num=len(spacerInds)) 

		for idx, item in enumerate(spacerInds):
			startPos = SeqFeature.ExactPosition(item)	# start and end pos of PAM
			endPos = SeqFeature.ExactPosition(item+PAM_length)  	

			if PAMside == 3:		# Cas9-like
				startLocInRefSeq = startPos+1
				endLocInRefSeq = startLocInRefSeq+spacerLength-1
			else:					# Cpf1-like
				startLocInRefSeq = endPos  #Starts immediately after PAM
				endLocInRefSeq = startLocInRefSeq+spacerLength  

			startLocInRefGenome = chromStartRG+startLocInRefSeq
			endLocInRefGenome = chromStartRG+endLocInRefSeq-1
			cutSiteInRefGenome = startLocInRefGenome+distanceToCutSiteFrom5pEnd

			# Only add the spacer if it is a certain distance from the previous spacer
			if (startLocInRefSeq-prevStartLocInRefSeq) > cutoff_spacing: 
				spacerNum += 1
				strand="+"
				if spacerNum > 1:
					distFromPrevSpacer = startLocInRefSeq-prevStartLocInRefSeq
				else:
					distFromPrevSpacer = 0
				if PAMside == 3:
					spacerAsStr = str(s[startLocInRefSeq-1:endLocInRefSeq])
					exactPAM = s[endLocInRefSeq:endLocInRefSeq+PAM_length];
				else:
					spacerAsStr = str(s[startLocInRefSeq:endLocInRefSeq])
					exactPAM = s[startLocInRefSeq-PAM_length:startLocInRefSeq];  # Python slices: second index is first char you *DON'T* want

				GCcontent = SeqUtils.GC(spacerAsStr);
				listItem = [spacerNum, strand, startLocInRefSeq, endLocInRefSeq, chromPos, startLocInRefGenome, endLocInRefGenome, 
							cutSiteInRefGenome, distFromPrevSpacer, spacerAsStr, exactPAM, GCcontent]				
				listSpacer.append(listItem)
				listDistBetweenSpacers.append(float(distFromPrevSpacer))
				prevStartLocInRefSeq=startLocInRefSeq
	
	
	print "Plus strand sgRNAs included after minimum spacing (> {limit}bp between sgRNAs): {num}".format(limit=cutoff_spacing,num=spacerNum) 
	spacerNumTotal=spacerNum
			
	# Search rev complement of PAM
	# print PAM
	# print PAM.reverse_complement()
	prevStartLocInRefSeq=-9999
	spacerNum=0
	if DoRevComp:
		if PAMside == 3:
			gbStringForSearch = s[:-spacerLength];   # get all but last ~20 bases of sequence
		else:
			gbStringForSearch = s[spacerLength:];
			
		spacerInds = SeqUtils.nt_search(gbStringForSearch,str(PAM.reverse_complement()))
		if len(spacerInds) > 1:	# matches found 
			del spacerInds[0] # first result from nt_search is regexp expansion
			#print "len line below {fname}".format(fname=len(spacerInds))
			formatDigitsN = int(math.ceil(math.log(len(spacerInds),10)));                                                                                                                                                                          
			print "Minus strand sgRNAs found: {num}".format(num=len(spacerInds))

			for idx, item in enumerate(spacerInds): 
				startPos = SeqFeature.ExactPosition(item) 
				endPos = SeqFeature.ExactPosition(item+PAM_length)   
				#print "Start pos: {num}  End pos: {num2}".format(num=startPos,num2=endPos)
				 			
				# Start and end locations are flipped here due to reverse strand
				if PAMside == 3:
					endLocInRefSeq = endPos+1  #flipped for reverse strand
					startLocInRefSeq = endLocInRefSeq+spacerLength-1  #flipped for reverse strand
				else:
					# startLocInRefSeq is 5' end of spacer on PAM-containing strand
					# endLocInRefSeq is 3' end of spacer on PAM-containing strand
					# Hence endLocInRefSeq <  startLocInRefSeq since this is reverse strand
					startLocInRefSeq = startPos + spacerLength # b/c spacer length is the offset between gbStringForSearch to RefSeq
					endLocInRefSeq = startLocInRefSeq - spacerLength +1

				startLocInRefGenome = chromStartRG+startLocInRefSeq-1
				endLocInRefGenome = chromStartRG+endLocInRefSeq-1
				cutSiteInRefGenome = startLocInRefGenome-distanceToCutSiteFrom5pEnd
												
				# Only add the spacer if it is a certain distance from the previous spacer
				if (startLocInRefSeq-prevStartLocInRefSeq) > cutoff_spacing: 
					spacerNum += 1
					strand="-"
					if spacerNum > 1:
						distFromPrevSpacer = startLocInRefSeq-prevStartLocInRefSeq
					else:
						distFromPrevSpacer = 0
					if PAMside == 3:# Cas9-like
						spacerRC = Seq(str(s[endLocInRefSeq-1:startLocInRefSeq]), IUPAC.ambiguous_dna)
						spacerAsStr = str(spacerRC.reverse_complement())
						exactPAM = str(Seq(str(s[endLocInRefSeq-(PAM_length+1):endLocInRefSeq-1]), IUPAC.ambiguous_dna).reverse_complement())
					else:	# Cpf1-like
						spacerRC = Seq(str(s[endLocInRefSeq-1:startLocInRefSeq]), IUPAC.ambiguous_dna)
						spacerAsStr = str(spacerRC.reverse_complement())
						exactPAM = str(Seq(str(s[startLocInRefSeq:startLocInRefSeq+PAM_length]), IUPAC.ambiguous_dna).reverse_complement())
						

					GCcontent = SeqUtils.GC(spacerAsStr);
					listItem = [spacerNum, strand, startLocInRefSeq, endLocInRefSeq, chromPos, startLocInRefGenome, endLocInRefGenome, 
								cutSiteInRefGenome, distFromPrevSpacer, spacerAsStr, exactPAM, GCcontent]				

					listSpacer.append(listItem)
					listDistBetweenSpacers.append(float(distFromPrevSpacer))
					prevStartLocInRefSeq=startLocInRefSeq		

		print "Minus strand sgRNAs included after minimum spacing (> {limit}bp between sgRNAs): {num}".format(limit=cutoff_spacing,num=spacerNum) 
		spacerNumTotal=spacerNumTotal+spacerNum;
	
	arrDistBetweenSpacers = np.asarray(listDistBetweenSpacers)
	meanDist = np.mean(arrDistBetweenSpacers)
	return (listSpacer, spacerNumTotal, meanDist)
	
	
'''

MAIN PROGRAM 

'''	
# Read Galaxy intervals
ginterval = GalaxyIntervalParser(intervalFile_INPUT)
# print ginterval['chrom']

# For each interval, 
for idx, currentGenomicInterval in enumerate(ginterval['chrom']):
	print idx
	print currentGenomicInterval
	print int(ginterval['startRefPos'][idx])
	# get list of NGGs within the interval 
	listSpacer, spacerNumTotal, meanDist=FindSpacersInterval(ginterval['chrom'][idx],int(ginterval['startRefPos'][idx]),int(ginterval['endRefPos'][idx]),PAM,cutoffSpacing, referenceGenomeForDAS, spacerLength, distanceToCutSiteFromPAM_bp);
	
	# Make unique IDs for each	
	for item in listSpacer:
		item.append(ginterval['assocRefSeq'][idx]+"_"+item[4]+"_"+str(item[5])+"_"+item[1])
	#print listSpacer[0]
	
	# Write separate CSV file containing NGGs for each gene/region
	
	if not os.path.exists(outDirectory):
		os.makedirs(outDirectory)
	
	#outShortName = "sg_" + ginterval['chrom'][idx] + '_' + str(ginterval['startRefPos'][idx]) + '_' + str(ginterval['startRefPos'][idx])
	csvFilename = outDirectory +"/" + outShortName+"_out.csv"
	header=("Spacer number","Strand","Start Pos","End Pos","Chromosome","ChrStart","ChrEnd","Cut Site","Distance from Previous Spacer","Spacer","PAM","Spacer GC","Seq ID");	

	if idx == 0:
		WriteColumns2CSV(listSpacer,csvFilename, True, header)
	else:
		WriteColumns2CSV(listSpacer,csvFilename, False, header)

print "BashRegion: Finished finding spacers!"
# After CSV file is complete, output FASTA version with UIDs and seqs	
spacerCSV2fa(csvFilename,outDirectory +"/" + outShortName+"_out.fa");
