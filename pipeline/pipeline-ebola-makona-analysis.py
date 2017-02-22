#################################################################
#################################################################
###############  ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, os, gzip, glob
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import pandas.rpy.common as com

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineEbolaMakonaAnalysis as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
### 1. Scripts
kallistoPath = '/Users/denis/Documents/Projects/scripts/tools/kallisto/kallisto'
fastqcPath = '/Applications/FastQC.app/Contents/MacOS/fastqc'

### 2. Files
macacaTranscriptomeFastq = 'rawdata.dir/transcriptome/Macaca_mulatta.Mmul_8.0.1.cdna.all.fa.gz'
sequencingFastq = 'rawdata.dir/fastq/*.fastq.gz'
sampleAnnotationFile = 'rawdata.dir/annotations/Animal_IDs_PBMC_Makona_.txt'
sampleClinicalDataFile = 'rawdata.dir/annotations/ZEBOV_makona_clinicalparameters.xls'

##### 2. R Connection #####
# Local source
rSource = 'pipeline/scripts/pipeline-ebola-makona-analysis.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Transcriptome processing
#######################################################
#######################################################

#############################################
########## 1. Build transcriptome index
#############################################

@follows(mkdir('f1-transcriptome.dir'))

@transform(macacaTranscriptomeFastq,
		   regex(r'.*/(.*).fa.gz'),
		   r'f1-transcriptome.dir/\1.idx')

def buildTranscriptomeIndex(infile, outfile):

	# Prepare statement
	statement = kallistoPath + ' index -i %(outfile)s %(infile)s' % locals()

	# Run
	os.system(statement)

#############################################
########## 2. Create annotations
#############################################

@transform(macacaTranscriptomeFastq,
		   regex(r'.*/(.*).fa.gz'),
		   r'f1-transcriptome.dir/\1_annotation.txt')

def makeTranscriptomeAnnotationTable(infile, outfile):

	# Read infile
	with gzip.open(infile,'r') as openfile:        
	    transcriptData = [x.strip() for x in openfile.readlines() if x[0] == '>']

	# Process data
	transcriptDataProcessed = [x.replace(' cdna', '').replace('>', '').split(' ', 6) for x in transcriptData]

	# Create dict
	transcriptDataDict = {x[0]:{y.split(':', 1)[0]:y.split(':', 1)[1] for y in x[1:]} for x in transcriptDataProcessed}

	# Select rows
	selectedRows = ['gene', 'gene_biotype', 'transcript_biotype', 'chromosome', 'gene_symbol' ,'description']

	# Convert to dataframe
	transcriptAnnotationDataframe = pd.DataFrame(transcriptDataDict).loc[selectedRows].T

	# Write file
	transcriptAnnotationDataframe.to_csv(outfile, sep='\t', index_label='ensembl_transcript_id')

#######################################################
#######################################################
########## S2. Read Quality Control
#######################################################
#######################################################

#############################################
########## 1. FastQC
#############################################

@follows(mkdir('f2-qc.dir'))

@transform(sequencingFastq,
		   regex(r'.*/(.*).fastq.gz'),
		   r'f2-qc.dir/\1_fastqc.html')

def runFastQC(infile, outfile):

	# Get FastQC output files
	outputHtmlFile = infile.replace('.fastq.gz', '_fastqc.html')
	outputZipFile = infile.replace('.fastq.gz', '_fastqc.zip')

	# Prepare statement
	statement = fastqcPath + ' %(infile)s; mv %(outputHtmlFile)s %(outfile)s; rm %(outputZipFile)s;' % locals()

	# Run statement
	os.system(statement)

#######################################################
#######################################################
########## S3. Read quantification
#######################################################
#######################################################

#############################################
########## 1. Quantification
#############################################

@follows(mkdir('f3-expression_quantification.dir'))

@transform(sequencingFastq,
		   regex(r'.*/(.*).fastq.gz'),
		   add_inputs(buildTranscriptomeIndex),
		   r'f3-expression_quantification.dir/\1')

def quantifyExpressionData(infiles, outfile):

	# Split infiles
	fastqFile, transcriptomeFile = infiles

	# Prepare statement
	statement = kallistoPath + ' quant -i %(transcriptomeFile)s -o %(outfile)s -b 100 --single -l 100 -s 20 %(fastqFile)s' % locals()

	# Run
	os.system(statement)

#######################################################
#######################################################
########## S4. Annotation Data
#######################################################
#######################################################

#############################################
########## 1. Sample annotations
#############################################

@follows(mkdir('f4-annotation_data.dir'))

@files([sampleAnnotationFile,
		sampleClinicalDataFile],
	   'f4-annotation_data.dir/ebola_pbmc-sample_annotations.txt')

def makeSampleAnnotationTable(infiles, outfile):

	# Split infiles
	sampleInfile, clinicalInfile = infiles

	# Read sample table
	sampleDataframe = pd.read_table(sampleInfile)

	# Read clinical table
	clinicalDataframe = pd.read_excel(clinicalInfile).dropna()

	# Fix Animal ID
	clinicalDataframe['Animal ID'] = [int(x) for x in clinicalDataframe['Animal ID']]

	# Add sample column
	sampleDataframe['FASTQ_name'] = [x.split('.')[0] for x in sampleDataframe['FASTQ_file']]

	# Add sample name
	sampleDataframe['sample_name'] = ['%(Tissue)s_%(animal_ID)s_day%(day)s' % locals() for Tissue, animal_ID, day in sampleDataframe[['Tissue', 'animal_ID', 'day']].as_matrix()]

	# Merge dataframes
	mergedDataframe = sampleDataframe.merge(clinicalDataframe, left_on=['animal_ID', 'day'], right_on=['Animal ID', 'DPI'])

	# Drop columns
	mergedDataframe.drop(['FASTQ_file', 'Animal ID', 'DPI'], axis=1, inplace=True)

	# Save file
	mergedDataframe.to_csv(outfile, sep='\t', index=False)	

#######################################################
#######################################################
########## S5. Annotation Data
#######################################################
#######################################################

#############################################
########## 1. Transcript expression
#############################################

@follows(mkdir('f5-expression_data.dir'))

@merge([makeSampleAnnotationTable,
		quantifyExpressionData],
	   'f5-expression_data.dir/ebola_pbmc-transcript_rawcounts.txt')

def getTranscriptRawcounts(infiles, outfile):

	# Split infiles
	annotationFile = infiles.pop(0)
	abundanceFiles = [os.path.join(x, 'abundance.tsv') for x in infiles]

	# Read annotation dataframe
	annotationDataframe = pd.read_table(annotationFile)

	# Get FASTQ-name conversion
	nameConversionDict = annotationDataframe[['FASTQ_name', 'sample_name']].set_index('FASTQ_name').to_dict()['sample_name']

	# Initialize empty dataframe
	mergedDataframe = pd.DataFrame()

	# Loop through files
	for abundanceFile in abundanceFiles:
	    
	    # Read abundance dataframe
	    abundanceDataframe = pd.read_table(abundanceFile, usecols=['target_id', 'est_counts'])
	    
	    # Get sample name
	    fastqFile = abundanceFile.split('/')[1]
	    
	    # Add sample name
	    abundanceDataframe['FASTQ_name'] = nameConversionDict[fastqFile]

	    # Concatenate
	    mergedDataframe = pd.concat([mergedDataframe, abundanceDataframe])

	# Cast dataframe
	mergedDataframeCast = mergedDataframe.pivot(index='target_id', columns='FASTQ_name', values='est_counts').astype(int)

	# Write file
	mergedDataframeCast.to_csv(outfile, sep='\t', index_label='ensembl_transcript_id')

#############################################
########## 2. Gene expression
#############################################

@transform('f5-expression_data.dir/ebola_pbmc-transcript_rawcounts.txt',#getTranscriptRawcounts,
		   suffix('transcript_rawcounts.txt'),
		   add_inputs(makeTranscriptomeAnnotationTable),
		   'gene_rawcounts.txt')

def getGeneRawcounts(infiles, outfile):

	# Split infiles
	transcriptInfile, annotationInfile = infiles

	# Get rawcount dataframe
	transcriptRawcountDataframe = pd.read_table(transcriptInfile)

	# Get annotation dataframe
	annotationDataframe = pd.read_table(annotationInfile, usecols=['ensembl_transcript_id', 'gene_symbol'])

	# Merge
	mergedDataframe = transcriptRawcountDataframe.merge(annotationDataframe, on='ensembl_transcript_id', how='inner')

	# Melt
	mergedDataframeMelt = pd.melt(mergedDataframe.drop('ensembl_transcript_id', axis=1), id_vars='gene_symbol')

	# Sum
	mergedDataframeCombined = mergedDataframeMelt.groupby(['gene_symbol', 'variable']).sum().reset_index()

	# Cast
	resultDataframe = mergedDataframeCombined.pivot(index='gene_symbol', columns='variable', values='value')

	# Save
	resultDataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S6. Normalize Expression Data
#######################################################
#######################################################

#############################################
########## 1. Run VST
#############################################

@follows(mkdir('f6-normalized_expression_data.dir'))

@transform(getGeneRawcounts,
		   regex(r'.*/(.*)rawcounts.txt'),
		   r'f6-normalized_expression_data.dir/\1vst.txt')

def runVst(infile, outfile):

	# Read expression dataframe
	expressionDataframe = pd.read_table(infile, index_col='gene_symbol')

	# Run function
	vstMatrix = r.runVST(com.convert_to_r_dataframe(expressionDataframe))

	# Convert to dataframe
	vstDataframe = com.convert_robj(vstMatrix)

	# Write file
	vstDataframe.to_csv(outfile, sep='\t', index_label='gene_symbol')

#######################################################
#######################################################
########## S7. Differential Expression Analysis
#######################################################
#######################################################

#############################################
########## 1. Characteristic Direction
#############################################

@follows(mkdir('f7-differential_expression.dir'))

@transform(runVst,
		   regex(r'.*/(.*)vst.txt'),
		   add_inputs(makeSampleAnnotationTable),
		   r'f7-differential_expression.dir/\1differential_expression.txt')

def runCharacteristicDirection(infiles, outfile):

	# Split infiles
	vstFile, annotationFile = infiles

	# Read expression data
	vstDataframe = pd.read_table(vstFile, index_col='gene_symbol')

	# Read annotation data
	annotationDataframe = pd.read_table(annotationFile, index_col='sample_name')

	# Get gene variance
	geneVariance = vstDataframe.apply(np.var, 1).sort_values(ascending=False)

	# Get variable genes
	variableGenes = geneVariance.index[geneVariance > 3].tolist()

	# Filter expression dataframe
	vstDataframeFiltered = vstDataframe.loc[variableGenes]

	# Get columns
	experimentColumns = vstDataframeFiltered.columns[:10].tolist()
	controlColumns = vstDataframeFiltered.columns[10:21].tolist()

	# Run characteristic direction
	cdResults = r.runCharacteristicDirection(com.convert_to_r_dataframe(vstDataframeFiltered), experimentColumns, controlColumns)

	print cdResults

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=4, verbose=1)
print('Done!')
