#def bootstrap_lib(lib_path):
#  import os, sys
#  from glob import glob#
#
#  exist_eggs = [os.path.basename(i) for i in sys.path];
#
# sys.path.insert(1,lib_path);#

#  for d in glob(os.path.join(lib_path,"*egg")):
#    if os.path.basename(d) in exist_eggs:
#      continue;
#    else:
#      sys.path.insert(1,d);

#bootstrap_lib("/net/snowwhite/home/welchr/lib/python2.7/site-packages");

import sys
import os
import re
import getopt
import os.path
import subprocess
import gzip
import numpy
import pandas
import time
import matplotlib
import pylab
import random

#FUNCTIONS
def calculatemaf(inputvcf,outprefix,samplestokeep,sepchr=False):
	#get the header from the vcf and figure out sample columns to keep
	sampleline = subprocess.Popen('zcat '+ inputvcf +' | grep CHROM', shell=True,stdout=subprocess.PIPE)
	sampleline = sampleline.communicate()[0].rstrip().split()
	colstokeep = "1,2"
	samplestokeeplist = list(samplestokeep)
	for columnnum in range(10,len(sampleline)):
		if sampleline[columnnum] in samplestokeeplist:
			colstokeep = colstokeep + ',' + str(columnnum)
	command = 'zcat ' + inputvcf + ' | grep -v "#" | cut -f ' + colstokeep + "| sed 's/:" + '[^\\t]*//g' + "' | sed '" + 's/0\/0/0/g' + "' | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' | sed 's/\.\/\./-1/g' | sed 's/0\/./1/g' | sed 's/.\/./2/g' > " + inputvcf.replace('.gz','') + '_GENOTYPES'
	print(command)
	os.system(command)
	genotypefile = pandas.read_table(inputvcf.replace('.gz','') + '_GENOTYPES',header=None)
	#remove first two columns chr and position
	genotypefile = genotypefile.set_index(genotypefile[genotypefile.columns[0]].map(str) + "\t" + genotypefile[genotypefile.columns[1]].map(str))
	del genotypefile[genotypefile.columns[0]]
	del genotypefile[genotypefile.columns[0]]
	print("getting genotype frequencies")
	frq1 = []
	frq2 = []
	nalleles = []
	n_chr = []
	counts1 = []
	counts2 = []
	for index,row in genotypefile.iterrows():
		counts = row.value_counts()
		n0 = 0
		n1 = 0
		n2 = 0
		if 0 in counts.index:
			n0 = counts.loc[0]
		if 1 in counts.index:
			n1 = counts.loc[1]
		if 2 in counts.index:
			n2 = counts.loc[2]
		if (n0 + n1 + n2) > 0:
			frq2.append((0.5*n1+n2)/(n0+n1+n2))
			frq1.append((0.5*n1+n0)/(n0+n1+n2))
			counts1.append(2*n0 + n1)
			counts2.append(2*n2 + n1)
		else:
			frq2.append(None)
			frq1.append(None)
			counts1.append(None)
			counts2.append(None)
		nalleles.append(2)
		n_chr.append(n0+n1+n2)
	genotypefile['N_ALLELES'] = nalleles
	genotypefile['N_CHR'] = n_chr
	genotypefile['FREQ1'] = frq1
	genotypefile['FREQ2'] = frq2
	colnames = genotypefile.columns
	for i in range(0,(len(colnames) - 4)):
		del genotypefile[colnames[i]]
	genotypefile.to_csv(outprefix + '.frq',sep="\t",index=True,quotechar='')
	os.system("cat " + outprefix + '.frq | sed \'s/"//g\' > ' + outprefix + 'TEMP')
	os.system('mv ' + outprefix + 'TEMP ' + outprefix + '.frq')
	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Variant Frequency file created" + "\n")

#this function just performs a filter on a column (or multiple columns) of a file by reading it in as a numpy table.
def evalexpression(exprs, filetofilter):
	exprs = exprs.strip()
	print exprs
	ex = ''
	open = False
	i = 0
	if '(' not in exprs:
		#base case
		if '==' in exprs:
			splitex = exprs.split('==')
			splitex[1] = splitex[1].strip()
			colname = splitex[0].strip()
			condition = splitex[1]
			if isinstance(filetofilter[colname][0],int):
				condition = int(condition)	
			if isinstance(filetofilter[colname][0],float):
				condition = float(condition)
			rowstokeep = numpy.where(filetofilter[colname] == condition)
		elif '!=' in exprs:
			splitex = exprs.split('!=')
			splitex[1] = splitex[1].strip()
			colname = splitex[0].strip()
			condition = splitex[1]
			if isinstance(filetofilter[colname][0],int):
				condition = int(condition)	
			if isinstance(filetofilter[colname][0],float):
				condition = float(condition)
			rowstokeep = numpy.where(filetofilter[colname] != condition)
		elif '<=' in exprs:
			splitex = exprs.split('<=')
			splitex[1] = float(splitex[1])
			colname = splitex[0].strip()
			condition = splitex[1]
			rowstokeep = numpy.where(filetofilter[colname] <= condition)			
		elif '<' in exprs:
			splitex = exprs.split('<')
			splitex[1] = float(splitex[1])
			colname = splitex[0].strip()
			condition = splitex[1]
			rowstokeep = numpy.where(filetofilter[colname] < condition)			
		elif '>=' in exprs:
			splitex = exprs.split('>=')
			splitex[1] = float(splitex[1])
			colname = splitex[0].strip()
			condition = splitex[1]
			rowstokeep = numpy.where(filetofilter[colname] >= condition)			
		elif '>' in exprs:
			splitex = exprs.split('>')
			splitex[1] = float(splitex[1])
			colname = splitex[0].strip()
			condition = splitex[1]
			rowstokeep = numpy.where(filetofilter[colname] > condition)
		elif exprs == '':
			rowstokeep = NULL
		else:
			rowstokeep = exprs
		return rowstokeep
	else:
		#find closed loop
		firstopen = exprs.index('(')
		numopen = 1
		numclose = 0
		index = firstopen
		while numopen != numclose:
			index = index + 1
			if exprs[index] == '(':
				numopen = numopen + 1
			if exprs[index] == ')':
				numclose = numclose + 1
		firstclose = index
		if ((index+1) < len(exprs)) and (exprs[index +1] == ' '):
			exprslist = list(exprs)
			del(exprslist[index+1])
			exprs = "".join(exprslist)
		if (index + 1) < len(exprs):
			nextsign = exprs[index+1].strip()
			if nextsign == '&':
				parta = evalexpression(exprs[firstopen+1:index],pedfile)
				partb = evalexpression(exprs[firstclose:].strip(), pedfile)
				return set(list(parta[0])).intersection(set(list(partb[0])))
			elif nextsign == '|':
				parta = evalexpression(exprs[firstopen+1:index],pedfile)
				partb = evalexpression(exprs[firstclose:].strip(), pedfile)			
				return set(list(parta[0])).union(set(list(partb[0])))
			else:
				sys.exit("Invalid expression in maf filter!")
		else:
			return evalexpression(exprs[firstopen+1:index],pedfile)

#this function creates the required group file for epacts
def create_group_file(outfilename, vcffilename, annofile, maffilterfile, maffilter, DELIM_MAF, DELIM_ANNO, samplestokeep):
	vcfname = vcffilename
	trimmedvcf = options['OUTPREFIX'] + vcfname.replace('.vcf.gz', '').replace(options['INPUTDIR'] + '/',"")+'_trimmed.txt'
	tofilter = False
	
	if maffilterfile != 'NA':
		tofilter = True
		maffile = pandas.read_table(maffilterfile,header=0,sep="\t")
		snprowstokeep = evalexpression(maffilter, maffile)
		passsnps = maffile.iloc[list(snpsrowstokeep),:]
		passsnps = passsnps['SNP']
	elif (maffilterfile == 'NA') and (maffilter != "NA"):	#in this case, we have to create an MAF file
		tofilter = True
		print("filtering by MAF")
		LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t filtering by MAF\n")
		if 'SEPCHR' not in options:
#			command = 'vcftools --gzvcf ' + options['INPUTDIR'] + "/" + options['VCFFILE'] + ' --keep '+ options['OUTPREFIX'] + '.samplestokeep.txt ' + '--freq2 --out ' + options['OUTPREFIX'] + '> /dev/null'
#			os.system(command)
#			LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + command + "\n")
			command = 'zcat ' + options['INPUTDIR'] + "/" + options['VCFFILE'] + ' | grep -v "^##" | cut -f 1,2,3 > ' + options['OUTPREFIX'] + '.variantids.txt'
			os.system(command)
			calculatemaf(vcffilename, options['OUTPREFIX'],samplestokeep)	# this creates an MAF file
			
		else:
			vcffilenametomatch = options['VCFFILE'].split('chr1')
			for f in os.listdir(options['INPUTDIR']):
				match = True
				for f2 in vcffilenametomatch:
					if f2 not in f:
						match = False
						break
				if match:
					calculatemaf(options['INPUTDIR'] + '/' + f, options['OUTPREFIX'],samplestokeep)	# this creates an MAF file
					command = 'zcat ' + options['INPUTDIR'] + "/" + f + ' | grep -v "^##" | cut -f 1,2,3 >> ' + options['OUTPREFIX'] + '.variantids.txt'
					os.system(command)

		maffile = pandas.read_table(options['OUTPREFIX'] + ".frq",skiprows=1,header=None)
		allidsfile = pandas.read_table(options['OUTPREFIX'] + '.variantids.txt',header=0)
		r1 = set(maffile.index[maffile[maffile.columns[5]] < float(maffilter)])
		r2 = set(maffile.index[maffile[maffile.columns[5]] > (1 - float(maffilter))])
		rowstokeep = r1.union(r2)
		r3 = set(maffile.index[maffile[maffile.columns[5]] > 0])
		rowstokeep = rowstokeep.intersection(r3)
		r4 = set(maffile.index[maffile[maffile.columns[5]] < 1])
		rowstokeep = rowstokeep.intersection(r4)
		maffile = maffile.iloc[list(rowstokeep)]
		passsnps = pandas.merge(maffile, allidsfile,left_on = [maffile.columns[0], maffile.columns[1]], right_on = [allidsfile.columns[0], allidsfile.columns[1]],how='inner')
		passsnps = list(passsnps['ID'])
		
	#make a bed file containing only the positions you need
	genes = []
	snps = []
	numline = 0
	bedfile = open(options['OUTPREFIX'] + vcfname.replace('.vcf.gz', '').replace(options['INPUTDIR'] + '/', '') + '_trimmed.bed',"w")

	#this is the annotation file with the group names (gene names)
	if 'gz' in annofile:
		filename = gzip.open(annofile,'rb')
	else:
		filename = open(annofile)
	#read in annotation file and store the gene and snp in 2 arrays
	line = filename.readline()
	for line in filename:
		line = line.rstrip()
		temp = line.split()
		if (not tofilter) or (temp[defaults['ANNOTVARCOL']] in passsnps):
			genes.append(temp[defaults['ANNOTGENECOL']])
			snps.append(temp[defaults['ANNOTVARCOL']])
			if '-' in temp[defaults['ANNOTPOSCOL']]:
				bedfile.write(temp[defaults['ANNOTPOSCOL']].split(':')[0] + "\t" + temp[defaults['ANNOTPOSCOL']].split(':')[1].split('-')[0] + "\t" + temp[defaults['ANNOTPOSCOL']].split(':')[1].split('-')[1] + "\n")
			else:
				bedfile.write(temp[defaults['ANNOTPOSCOL']].split(':')[0] + "\t" + temp[defaults['ANNOTPOSCOL']].split(':')[1] + "\t" + str(int(temp[defaults['ANNOTPOSCOL']].split(':')[1])+1) + "\n")
	filename.close()
	bedfile.close()
	tabixcommand = 'tabix ' + vcfname + ' -B ' + options['OUTPREFIX'] + vcfname.replace('.vcf.gz', '').replace(options['INPUTDIR'] + '/', '') + '_trimmed.bed | grep -v "#" > ' + trimmedvcf 
	print(tabixcommand)
	os.system(tabixcommand)
	
	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "finished creating trimmed vcf\n")
	genes = numpy.asarray(genes)
	snps = numpy.asarray(snps)
	geneinfo = numpy.column_stack((genes, snps))

	#opening in the trimmed vcf file and reading in snp information and storing it in the dictionary snpinfo
	vcffile = open(trimmedvcf)
	snpinfo = dict()
	for line in vcffile:
		temp = line.split()
		snpinfo[temp[2]] = temp[0]+":"+temp[1]+"_"+temp[3]+'/'+temp[4]
	vcffile.close()

	print("SNP Info obtained\n")
	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "SNP info obtained\n")
	#now sort geneinfo by gene names
	if genes.size > 0:
		geneinfo = geneinfo[geneinfo[:,0].argsort()]

	variants = dict() #storing the gene annotation as it should be written to the output file
	i = 0
	#go through each row of geneinfo, and if successive rows have the same gene, concatenate the snps together in the variants dictionary using the gene name as key
	while i < genes.size:
		if snpinfo.has_key(geneinfo[i,1]):
			variants[geneinfo[i,0]] = snpinfo[geneinfo[i,1]]
		else:
			variants[geneinfo[i,0]] = ""
		while ((i+1) < genes.size)and(geneinfo[i+1,0] == geneinfo[i,0]):
			i = i+1
			if snpinfo.has_key(geneinfo[i,1]):
				variants[geneinfo[i,0]] = variants[geneinfo[i,0]]+"\t"+snpinfo[geneinfo[i,1]]
		i = i + 1

	outname = outfilename
	print outname
	out = open(outname, 'w') #file to which to write annotations
	for eachkey in variants.keys():
		variants[eachkey] = variants[eachkey].rstrip("\t")
		temp = variants[eachkey].split()
		temp.sort()
		variants[eachkey] = ""
		for i in temp:
			variants[eachkey] = variants[eachkey]+i+"\t"
		variants[eachkey] = variants[eachkey].rstrip("\t")
		if variants[eachkey] != "":
				out.write(eachkey+"\t"+variants[eachkey]+"\n")
	out.close()

def check_options(user_option):
#Check right combination of user options is provided
#	required['group'] = ['VCF', 'PED', 'OUT','GROUPF','PHENO']
#	required['single'] = ['VCF', 'PED', 'OUT','GROUPF']
	return 0
	
def check_file_exists(filename):
	if os.path.exists(filename):
		return 0
	else:
		print("Error: File " + filename + " does not exist\n")
		sys.exit()
##END FUNCTIONS

#Starting main

#if no options are specified, then print out message
if not len(sys.argv) > 1:
	print("runepacts.py\n" + "Option: --config required\n" + "Usage: python runepacts.py --config config.txt\n")
	sys.exit()

#help message
if sys.argv[1] == '--help':
	print "runepacts options: --config"
	print "Config File:"
	print "Required Options"
	print "INPUTDIR\tDirectory containing input files\nVCF\tFile name of gzipped vcf\nPEDFILE\tSample information file with column headers"
	print "MODEL\tPhenotype and covariates to use\nTEST\tEPACTS test to use in the format group=testname or single=testname"
	print "OUTPREFIX\tOutput folder/prefix name prefixed to output"
	print "Use python runepacts.py --example to view example config file"
	sys.exit()

try:
	options,args = getopt.getopt(sys.argv[1:],'c', ["config"])
	if not os.path.exists(args[0]):
		print("Config file doesn't exist!")
		sys.exit()

except getopt.GetoptError,err:
	print str(err)
	print("RunEPACTS options: --config\n")
	sys.exit()

#Then, read in config file
options = dict()
print("Reading config file...")

#This section counts the number of lines in the config file, so we know how many lines to read in
configfile = open(args[0])
lines_in_file = subprocess.Popen('cat '+ args[0] +' | wc -l', shell=True,stdout=subprocess.PIPE)
lines_in_file = lines_in_file.communicate()[0]
line_num = 0
line = configfile.readline()

#Set defaults
defaults = dict()
defaults['VCFFILE'] = '999'
defaults['PEDFILE'] = '999'
defaults['FILTERPED'] = '999'
defaults['ANNOTFILE'] = '999'
defaults['FILTERMAF'] = 'NA'
defaults['MAFFILE'] = '999'
defaults['GROUPFILE'] = '999'
defaults['KINSHIPFILE'] = '999'
defaults['MODEL'] = '999'
defaults['PVALUETHRESHOLD'] = 0.05
defaults['OUTPREFIX'] = ''
defaults['INPUTDIR'] = ''
defaults['SINGLEMARKERTEST'] = 'q.linear'
defaults['MINMAF'] = 0
defaults['GENEMINMAC'] = 0
defaults['MAXMAF'] = 1
defaults['EPACTSGROUPFILE'] = 'NA'
defaults['VERBOSE'] = 'OFF'
defaults['GENELIST'] = 'NA'
defaults['MINMAC'] = 0
defaults['ANNOTGENECOL'] = 4
defaults['ANNOTVARCOL'] = 0
defaults['ANNOTPOSCOL'] = 1

##open log file to write to
logfilename = "runepacts_" +  time.asctime(time.localtime(time.time())).replace(" ","_") + ".log"
LOGFILE = open(logfilename, "w")
numtests = 0

#This while loop reads in each line, till the end, and if the line has a 'PROCESS', it does the required processing, and run the epacts test
while line_num < int(lines_in_file):
	#ignore empty lines
	while len(line) < 3:
		line = configfile.readline()
		line_num = line_num + 1
	TESTS = []
	#keep reading in lines, and storing the parameters till you hit a 'PROCESS'
	while line[0:7] != "PROCESS":
		line = line.rstrip()
		splitline = line.split("\t")
		if len(splitline) < 2:
			print("Argument " + splitline[0] + 'has no value\n')
			LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Argument " + splitline[0] + "has no value\n")
			sys.exit()
		if splitline[0] == 'TEST':
			TESTS.append(splitline[1])
		else:
			options[splitline[0]] = splitline[1]
		line = configfile.readline()
		line_num = line_num + 1
	
	#Now, line has 'PROCESS'
	print("Running test "+ str(numtests + 1))
	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Running test " + str(numtests) + "\n")
	numtests = numtests + 1
	
	#TO COMPLETE LATER: CHECK IF ALL THE REQUIRED OPTIONS FOR THIS TEST ARE DEFINED, and the ones which are defined are valid
#	check_file_exists(options['INPUTDIR']+'/'+options['VCFFILE'])
#	check_file_exists(options['INPUTDIR']+'/'+options['PEDFILE'])
#	check_file_exists(options['INPUTDIR']+'/'+options['ANNOTFILE'])

	#if folder to write files in does not exist, make it
	if not os.path.isdir(re.sub('/.*', '', options['OUTPREFIX'])):
		os.system('mkdir ' + re.sub('/.*', '', options['OUTPREFIX']))

	#set the defaults as the input options
	if 'VCFFILE' in options:
		defaults['VCFFILE'] = options['VCFFILE']
	if 'PEDFILE' in options:
		defaults['PEDFILE'] = options['PEDFILE']
	if 'MODEL' in options:
		defaults['MODEL'] = options['MODEL']
	if 'GROUPFILE' in options:
		defaults['GROUPFILE'] = options['GROUPFILE']
		check_file_exists(options['INPUTDIR'] + '/' + options['GROUPFILE'])
	if 'ANNOTFILE' in options.keys():
		defaults['ANNOTFILE'] = options['ANNOTFILE']
#		check_file_exists(options['INPUTDIR'] + '/' + options['ANNOTFILE'])
	if 'MAFFILE' in options.keys():
		defaults['MAFFILE'] = options['MAFFILE']
	if 'KINSHIPFILE' in options.keys():
		defaults['KINSHIPFILE'] = options['KINSHIPFILE']
		check_file_exists(options['INPUTDIR'] + '/' + options['KINSHIPFILE'])
	if 'INPUTDIR' in options:
		defaults['INPUTDIR'] = options['INPUTDIR']
	if 'SINGLEMARKERTEST' in options.keys():
		defaults['SINGLEMARKERTEST'] = options['SINGLEMARKERTEST']
	if 'FILTERMAF' in options.keys():
		defaults['FILTERMAF'] = options['FILTERMAF']
	if 'PVALUETHRESHOLD' in options.keys():
		PVALUETHRESHOLD = options['PVALUETHRESHOLD']
	else:
		PVALUETHRESHOLD = 0.05
	if 'MINVARS' in options.keys():
		MINVARS = options['MINVARS']
	else:
		MINVARS = 2
	if 'GENEMINMAC' in options.keys():
		GENEMINMAC = options['GENEMINMAC']
	else:
		GENEMINMAC = 5
	if 'ANNOTGENECOL' in options.keys():
		defaults['ANNOTGENECOL'] = int(options['ANNOTGENECOL']) - 1
	if 'ANNOTVARCOL' in options.keys():
		defaults['ANNOTVARCOL'] = int(options['ANNOTVARCOL']) - 1
	if 'ANNOTPOSCOL' in options.keys():
		defaults['ANNOTPOSCOL'] = int(options['ANNOTPOSCOL']) - 1
	#now you have all the options for the test
	#this block sets the default options (which change if specified by the user)
	kinshipcommand = ''

	if 'EPACTSDIR' not in options.keys():
		epacts = 'epacts'
	else:
		epacts = options['EPACTSDIR'] + '/epacts'	
		defaults['EPACTSDIR'] = options['EPACTSDIR'] + '/epacts'

	OUTDIR = re.sub(".*/",'',options['OUTPREFIX'])

	##This section has to be worked on. The user can give a list of delimiters, one for each file. This section parses that option.
	DELIM_PED = '\t'
	DELIM_ANNO = '\t'
	DELIM_MAF = '\t'
	if 'DELIMITER' in options.keys():
		DELIM = options['DELIMITER']
		delimsplit = DELIM.split(",")
		for i in delimsplit:
			a = i.strip()
			splut = a.split("=")
			splut[0] = splut[0].strip()
			splut[1] = splut[1].strip()
			if splut[0] == 'PEDFILE':
				DELIM_PED = splut[1].strip().replace('COMMA', ',').replace('TAB', '\t').replace('SPACE', '\s').replace('WHITESPACE','\s+')
			elif splut[0] == 'ANNOTFILE':
				DELIM_ANNO = splut[1].strip().replace('COMMA', ',').replace('TAB', '\t').replace('SPACE', '\s').replace('WHITESPACE','\s+')
			elif splut[0] == 'MAFFILE':
				DELIM_MAF = splut[1].strip().replace('COMMA', ',').replace('TAB', '\t').replace('SPACE', '\s').replace('WHITESPACE','\s+')
			else:
				print(splut[0] + " not a valid delimiter!")

	#the user can ask for only some columns to be used from the input 'ped' file. this section parses that input
	if 'PEDCOLUMNS' in options.keys():
		pedcolumns = options['PEDCOLUMNS'].split(",")
	else:
		pedcolumns = []

	#The '~' in the model option indicates that there are covariates to be used, since the MODEL option is MODEL = PHENO ~ COV1 + COV2
	if '~' in options['MODEL']:
		covariates = options['MODEL'].split('~')[1].split('+')
	else:
		covariates = []
	
	#the user can allow for extra epacts commands to be included in the run
	if 'EPACTSCMD' in options.keys():
		usercommand = ' ' + options['EPACTSCMD'].strip("'") + ' '
	else:
		usercommand = ''

	if 'SEPCHR' in options.keys():
		if options['SEPCHR'] == 'ON':
			sepchr = ' -sepchr '
		else:
			sepchr = ''
	else:
		sepchr = ''

	if 'GENOTYPE' in options.keys():
		field = ' -field ' + options['GENOTYPE']
	else:
		field = ' '
	
	if 'MINMAF' in options.keys():
		minmaf = ' -min-maf ' + str(options['MINMAF'])
	else:
		minmaf = ' -min-maf ' + str(defaults['MINMAF'])
	
	if 'MINMAC' in options.keys():
		minmac = ' -min-mac ' + str(options['MINMAC'])
	else:
		minmac = ' -min-mac ' + str(defaults['MINMAC'])
	
	if 'MAXMAF' in options.keys():
		maxmaf = ' -max-maf ' + str(options['MAXMAF'])
	else:
		maxmaf = ' -max-maf ' + str(defaults['MAXMAF'])
	############################################################################
	
	#First step: Index vcffile if indexed file does not already exist
	if not os.path.exists(defaults['INPUTDIR'] + '/' + options['VCFFILE'] + '.tbi'):
		print("Indexing vcf. Warning: This may take longer for large VCFs\n")
		LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Indexing vcf\n")
		tabix_command = 'tabix -pvcf -f ' + defaults['INPUTDIR'] + '/' + options['VCFFILE'] + " > /dev/null"
		os.system(tabix_command)
		LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + tabix_command + "\n")

	#Then, make a ped file after using the filters provided
	print "Creating ped file...\n"
	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Creating ped file...\n")

	pedfile = pandas.read_table(defaults['INPUTDIR'] + '/' + defaults['PEDFILE'] , sep="\t",header=0)
	headersplit = pedfile.columns
	#the variable filter contains the awk filter to use for the ped file. we are replacing the filter-column name with the column
	#number to get the awk command
	print "Filtering ped file...\n"
	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Filtering ped file...\n")
	filter = ''
	if 'FILTERPED' in options.keys():
		filter = options['FILTERPED']
		filteredrows = evalexpression(filter, pedfile)
		pedfile = pedfile.iloc[list(filteredrows),:]
	pedfile[pedfile.columns[1]].to_csv(options['OUTPREFIX']+ ".samplestokeep.txt", sep="\t",index=False, na_rep='NA')
	samplestokeep = pedfile[pedfile.columns[1]]
	if len(pedcolumns) > 0:
		if len(covariates) > 0:
			for i in range(0,len(covariates)):
				if covariates[i].strip() not in pedcolumns:
					pedcolumns.append(covariates[i].strip())
		pedfile = pedfile.loc[:,pedcolumns]

	col0 = pedfile.columns[0]
	#default columns
	pedfile.rename(columns={col0:'#FAM_ID'}, inplace=True)
	col1 = pedfile.columns[1]
	pedfile.rename(columns={col1:'IND_ID'}, inplace=True)
	col2 = pedfile.columns[2]
	pedfile.rename(columns={col2:'FAT_ID'}, inplace=True)
	col3 = pedfile.columns[3]
	pedfile.rename(columns={col3:'MOT_ID'}, inplace=True)
	col4 = pedfile.columns[4]
	pedfile.rename(columns={col4:'SEX'}, inplace=True)

	#Checking if there are categorical covariates. If yes, convert to dummy variables
	print "Adding covariates...\n"
	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Adding covariates...\n")
	covariatecommand = ''

	for i in range(0,len(covariates)):
		#if covariate is categorical, then create dummy variables for it
		if 'as.factor' in covariates[i]:
			covariate = covariates[i].replace('as.factor(', '')
			covariate.strip(')')
			values = list(set(list(pedfile[covariate])))			#get the categories
			numvalues = len(set(list(pedfile[covariate])))			#get the number of categories
			covcolumns = numpy.empty([len(pedfile),numvalues-1])
			for row in range(0,len(pedfile)):
				for col in range(0,numvalues-1):
					if pedfile[covariate][row] == values[col]:
						covcolumns[row][col] = 1
					else:
						covcolumns[row][col] = 0
			num = 0
			for numvalue in range(0,len(values)-1):
				pedfile[covariate + '_' + str(num)] = covcolumns[:,num]
				num = num + 1
				covariatecommand = covariatecommand + ' -cov ' + covariate + '_' + str(num)

		else:		#no dummy variable has to be added to the file
			covariatecommand = covariatecommand + ' -cov ' + covariates[i] + ' '

	pedfile.to_csv(options['OUTPREFIX'] + '.pheno.ped', sep="\t",index=False, na_rep='NA')
	#Covariates have now been added to the ped file and it has been written out

	###create groupfile if one is not already given by the user
	if 'GROUPFILE' not in options.keys():
		print "Creating groupfile...\n"
		LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + "Creating groupfile...\n")

		if 'ANNOTFILE' not in options:
			if 'SEPCHR' not in options:
				annotatecommand = epacts + ' anno -in ' + defaults['INPUTDIR'] + '/' + defaults['VCFFILE'] + ' -out ' + options['OUTPREFIX'] + '.anno.vcf.gz'
				os.system(annotatecommand)
				annotatecommand = epacts + ' make-group --vcf ' + options['OUTPREFIX'] + '.anno.vcf.gz --out ' + options['OUTPREFIX'] + '.annotated.txt --format epacts --nonsyn'
				print(annotatecommand)
				os.system(annotatecommand)
				finalgroupfilename = options['OUTPREFIX'] + '.annotated.txt'
			else:
				vcffilenametomatch = options['VCFFILE'].split('chr1')
				for f in os.listdir(options['INPUTDIR']):
					match = True
					for f2 in vcffilenametomatch:
						if f2 not in f:
							match = False
							break
					if match:
						annotatecommand = epacts + ' anno -in ' + defaults['INPUTDIR'] + '/' + f + ' -out ' + options['OUTPREFIX'] + '.' + f
						os.system(annotatecommand)
						annotatecommand = epacts + ' make-group --vcf ' + options['OUTPREFIX'] + '.' + f + ' --out ' + options['OUTPREFIX'] + '.annotated.' + f  + ' --format epacts'
						os.system(annotatecommand)
						annotatecommand = 'cat ' + options['OUTPREFIX'] + '.' + f + ' >> ' + options['OUTPREFIX'] + '.annotated.txt'
						os.system(annotatecommand)
				finalgroupfilename = options['OUTPREFIX'] + '.annotated.txt' 
		
		else:
			if '.gz' in options['ANNOTFILE']:
				cat = 'zcat '
			else:
				cat = 'cat '

			if 'FILTERANNOT' in options.keys():
				print "Filtering annotations...\n"
				LOGFILE.write("Filtering annotations...\n")
				filter = options['FILTERANNOT']
				annofile = pandas.read_table(defaults['INPUTDIR'] + '/' + options['ANNOTFILE'],	header=0,sep="\t",compression='gzip')
				annoheadersplit = annofile.columns
				
				annorowstokeep = evalexpression(filter, annofile)
				annofile = annofile.iloc[list(annorowstokeep),:]
				annofile.to_csv(defaults['INPUTDIR'] + '/' + options['ANNOTFILE'].replace('.gz', '') + '.filtered.txt',sep="\t", index=False)
				print("still filtering...\n")
				if 'MAFFILE' not in options.keys():
					create_group_file(options['OUTPREFIX'] + '_groups.txt', defaults['INPUTDIR'] + '/' + defaults['VCFFILE'], defaults['INPUTDIR'] + '/' + options['ANNOTFILE'].replace('.gz', '') + '.filtered.txt', 'NA', defaults['FILTERMAF'], DELIM_MAF, DELIM_ANNO,samplestokeep)
				else:
					create_group_file(options['OUTPREFIX'] + '_groups.txt',defaults['INPUTDIR'] + '/' + options['VCFFILE'], defaults['INPUTDIR'] + '/' + options['ANNOTFILE'].replace('.gz', '') + '.filtered.txt', \
								defaults['INPUTDIR'] + '/' + options['MAFFILE'] , defaults['FILTERMAF'], DELIM_MAF, DELIM_ANNO,samplestokeep)

			else:
				if 'MAFFILE' not in options.keys():
					create_group_file(options['OUTPREFIX'] + '_groups.txt', defaults['INPUTDIR'] + '/' + defaults['VCFFILE'], \
								defaults['INPUTDIR'] + '/' + defaults['ANNOTFILE'], 'NA', 'NA', DELIM_MAF, DELIM_ANNO,samplestokeep)
				else:
					create_group_file(options['OUTPREFIX'] + '_groups.txt', defaults['INPUTDIR'] + '/' + defaults['VCFFILE'], \
								defaults['INPUTDIR'] + '/' + defaults['MAFFILE'], defaults['FILTERMAF'], DELIM_MAF, DELIM_ANNO,samplestokeep)

			finalgroupfilename = options['OUTPREFIX'] + '_groups.txt'
	else:
		finalgroupfilename = defaults['INPUTDIR'] + '/' + defaults['GROUPFILE']

	if 'GENELIST' in options:
		outputname = random.randrange(1,1000000)
		getonlygenelistcommand = 'cat ' + finalgroupfilename + ' | grep -wf ' + options['INPUTDIR'] + '/' + options['GENELIST'] + ' > ' + outputname
		os.system(getonlygenelistcommand)
		os.system('mv ' + outputname + ' ' + finalgroupfilename)
		
	phenotype = options['MODEL'].split('~')[0]

	for TEST in TESTS:
		#Check if test includes an emmax test. Then check if kinship file has been provided. If not, then make one.
		if ('emmax' in TEST.split('=')[1]) or ('mmskat' in TEST.split('=')[1]) or ('q.emmax' in options['SINGLEMARKERTEST']):
			if 'KINSHIPFILE' in options.keys():
				if not os.path.exists(defaults['INPUTDIR'] + '/' + defaults['KINSHIPFILE']):
					print("Kinship file specified does not exist")
					LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Kinship file specified does not exist\n")
					sys.exit()
				kinshipcommand = ' -kin ' + defaults['INPUTDIR'] + '/' + defaults['KINSHIPFILE'] + ' '

			else:
				#make kinship file from input using emmax
				makekinshipcommand = epacts + ' make-kin '  + ' -vcf ' + defaults['INPUTDIR'] + '/' + defaults['VCFFILE'] \
										+ ' -ped ' + options['OUTPREFIX'] + '.pheno.ped ' + ' -out ' + options['OUTPREFIX'] \
										+ '.kinf' \
										+ ' -min-maf 0.01 -min-callrate 0.95' + ' -run 1'
				#print makekinshipcommand
				os.system(makekinshipcommand)
				LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + makekinshipcommand + "\n")
				kinshipcommand = ' -kin ' + options['OUTPREFIX']  + '.kinf '

		#now create the epacts command
		epacts_command = epacts + ' ' + TEST.split('=')[0] + ' -test ' + TEST.split('=')[1] + ' -vcf ' \
										+ defaults['INPUTDIR'] + '/' + defaults['VCFFILE'] + ' -pheno ' + phenotype + ' -ped ' \
										+ options['OUTPREFIX'] + '.pheno.ped ' + ' -groupf ' \
										+ finalgroupfilename + ' -out ' + options['OUTPREFIX'] + '.' + TEST.split('=')[1] \
										+ minmaf + ' ' + minmac + ' ' + covariatecommand + usercommand + sepchr + kinshipcommand + ' -run 1'

		if ('RUNEPACTS' in options.keys()) and (options['RUNEPACTS'] == 'FALSE'):
			print epacts_command
			LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + epacts_command + "\n")
		else:
			print epacts_command
			os.system(epacts_command)
			LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + epacts_command + "\n")
			print("Epacts run! Formatting Output...\n")
	############
	
	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Epacts run! Formatting output..." + "\n")

	for TEST in TESTS:
		outputtest = pandas.read_table(options['OUTPREFIX'] + '.' + TEST.split('=')[1] + '.epacts',sep="\t",header=0)	#reading in output file for test
		if 'BEGIN' not in outputtest.columns:
			outputtest.rename(columns={'BEG':'BEGIN'}, inplace=True)
		
		outputtest = outputtest.sort(['#CHROM','BEGIN'])
		if TEST == TESTS[0]:
			output = outputtest
		else:
			output = output.append(outputtest)
	
	#the output dataframe contains concatenated results from all EPACTS tests run
	output = output.sort(columns=[output.columns[0], 'BEGIN'],inplace=False)
	allgenes = output['MARKER_ID']
	allgenesoutput = output
	output = output[output.PVALUE <= float(PVALUETHRESHOLD)]
	numvars = output[['MARKER_ID','PASS_MARKERS']]
	numvars = numvars.drop_duplicates()
	
	#if output is empty, then ??
	if len(output) == 0:
		print("Epacts groupwise test returned no significant results")
		continue

	singlemarkertest = 'false'
	#get minvars and minmac
	genespassingfilters = dict()
	print("Creating input for single marker test...")
	if ('SINGLEMARKERTEST' not in options) or (options['SINGLEMARKERTEST'] != 'FALSE'):
		singlemarkertest = 'true'
		#get the variants in each gene if test if not VT
		genes = output['MARKER_ID']
		genes = genes.drop_duplicates()
		print(finalgroupfilename)
		groupfile = open(finalgroupfilename)
		markerlistforgenes = dict()
		print("Reading groupfile\n")
		#read in groupfile, and get the list of variants belonging to each gene. Store this in markerlistforgenes
		for line in groupfile:
			line = line.rstrip()
			temp = line.split()
			genename = temp[0]
			temp = temp[1:]

			for tempindex in range(0,len(temp)):
				temp[tempindex] = re.sub(r"_.*","",temp[tempindex])
			markerlistforgenes[genename] = temp
		groupfile.close()

		#use VCFtools to get the mac for each variant from the vcffile

		if 'GENEMINMAC' in options:
#			maccommand = 'vcftools --gzvcf ' + options['INPUTDIR'] + '/' + defaults['VCFFILE'] + ' --keep ' + options['OUTPREFIX'] + ".samplestokeep.txt"  + ' --counts --out ' + options['OUTPREFIX'] + '.allelecounts.temp'
#			os.system(maccommand)
#			LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + maccommand + "\n" )
#			maccommand = 'cat ' + options['OUTPREFIX'] + ".allelecounts.temp.frq.count | sed 's/[A|T|G|C]://g' > " + options['OUTPREFIX'] + '.variantcounts.txt'
#			os.system(maccommand)
#			maccommand = 'rm ' + options['OUTPREFIX'] + '.allelecounts.temp.frq.count > /dev/null'
#			os.system(maccommand)
#			maccommand = 'cat ' + options['OUTPREFIX'] + ".variantcounts.txt | sed 's/[A|T|G|C]://g'| awk '{min=$5;if($6 < $5){min = $6};print $1\"\t\"$2\"\t\"min}' > " + options['OUTPREFIX'] + '.maccounts.txt'
#			os.system(maccommand)
			if os.path.isfile(options['OUTPREFIX'] + '.frq'):
				command = 'cat ' + options['OUTPREFIX'] + '.frq | ' + "awk '{print $1\"\\t\"$2\"\\t\"$4*$5\"\\t\"$4*$6}' | awk '{min=$3;if($4 < min){min=$4;};print $1\"\\t\"$2\"\\t\"min}' > " + options['OUTPREFIX'] + '.maccounts.txt'
				os.system(command)
			else:
				calculatemaf(options['INPUTDIR'] + '/' + options['VCFFILE'], options['OUTPREFIX'], samplestokeep)
				command = 'cat ' + options['OUTPREFIX'] + '.frq | ' + "awk '{print $1\"\\t\"$2\"\\t\"$4*$5\"\\t\"$4*$6}' | awk '{min=$3;if($4 < min){min=$4;};print $1\"\\t\"$2\"\\t\"min}' > " + options['OUTPREFIX'] + '.maccounts.txt'
				os.system(command)
			macs = pandas.read_table(options['OUTPREFIX'] + '.maccounts.txt',skiprows=1,header=None)
			macs['MARKER'] = macs[macs.columns[0]].map(str) + ':'  + macs[macs.columns[1]].map(str)
			del macs[macs.columns[0]]
			del macs[macs.columns[1]]

			#For all genes (not only the significant ones), get the minor allele count for that gene. Filter by GENEMINMAC and MINVARS, output to genespassingfilters
			LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Filtering genes by MINMAC/MINVARS" + "\n")
			for gene in allgenes:
				mac_gene = 0
				genename = re.sub(".*_", "", gene)
				for ms in markerlistforgenes[genename]:
					macstoadd = macs[macs['MARKER'] == ms]
					macstoadd = macstoadd[macs.columns[0]]
					mac_gene = mac_gene + int(macstoadd)
				if mac_gene > int(GENEMINMAC):
					varnums = numvars[numvars['MARKER_ID'] == gene]
					if (len(varnums) > 0) and (varnums['PASS_MARKERS'].loc[varnums.index[0]] >= int(MINVARS)):
						genespassingfilters[gene] = re.sub(".*_", "", gene)

			genespassingfilters = pandas.DataFrame(genespassingfilters.keys())
			genespassingfilters['GENENAME'] = genespassingfilters[genespassingfilters.columns[0]].replace(r".*_","")
			genespassingfilters.to_csv(options['OUTPREFIX'] + ".genespassingfilters.txt",index=False,index_label=False,sep="\t")
		elif 'MINVARS' in options:
			for gene in allgenes:
				genename = re.sub(".*_", "", gene)
				if len(markerlistforgenes[genename]) > int(MINVARS):
					genespassingfilters[gene] = re.sub(".*_", "", gene)
			genespassingfilters = pandas.DataFrame(genespassingfilters.keys())
			genespassingfilters['GENENAME'] = genespassingfilters[genespassingfilters.columns[0]].replace(r".*_","")
			genespassingfilters.to_csv(options['OUTPREFIX'] + ".genespassingfilters.txt",index=False,index_label=False,sep="\t")			
		else:
			genespassingfilters = pandas.DataFrame(allgenes)
			genespassingfilters['GENENAME'] = genespassingfilters[genespassingfilters.columns[0]].replace(r".*_","")
			genespassingfilters.to_csv(options['OUTPREFIX'] + ".genespassingfilters.txt",index=False,index_label=False,sep="\t")

		os.system('zcat ' + defaults['INPUTDIR'] + '/' + defaults['VCFFILE'] + ' | grep -m 1 "#CHROM" > ' + options['OUTPREFIX'] + '.singlemarkers.vcf')
		LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + 'zcat ' + defaults['INPUTDIR'] + '/' + defaults['VCFFILE'] + ' | grep -m 1 "#CHROM" > ' + options['OUTPREFIX'] + '.singlemarkers.vcf\n')
		notthere = 0
		
		for eachgeneindex in genes.index:
			gene = genes[eachgeneindex]
			genename = gene.split('_')[1]
			if genename not in markerlistforgenes:
				notthere = notthere + 1
				continue
			markers = markerlistforgenes[genename]
			
			for marker in markers:
				markerchr = marker.split('_')[0].split(':')[0]
				markerposition = marker.split('_')[0].split(':')[1]
				tabixcommand = 'tabix ' + options['INPUTDIR'] + '/' + options['VCFFILE'] + ' ' + markerchr + ':' + markerposition + '-' + markerposition + ' >> ' + options['OUTPREFIX'] + '.singlemarkers.vcf'
				os.system(tabixcommand)
				LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + tabixcommand + "\n")
		

		#tabix and zip
		print("Indexing file...")
		
		os.system('cat ' + options['OUTPREFIX'] + '.singlemarkers.vcf | sort -g -k1 -k2 > TEMPXXX')
		os.system('mv TEMPXXX ' + options['OUTPREFIX'] + '.singlemarkers.vcf')
		os.system('bgzip -f ' + options['OUTPREFIX'] + '.singlemarkers.vcf')
		LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + 'bgzip ' + options['OUTPREFIX'] + '.singlemarkers.vcf'+ "\n")
		os.system('tabix -pvcf ' + options['OUTPREFIX'] + '.singlemarkers.vcf.gz')
		
		#run epacts for single markers to get single marker p values
		print("Running Single Marker EPACTS test...")
		epactscommand = epacts + ' single --vcf ' + options['OUTPREFIX'] + '.singlemarkers.vcf.gz -ped ' + options['OUTPREFIX'] + '.pheno.ped ' + '-pheno ' + phenotype + ' ' + covariatecommand + ' -test ' + defaults['SINGLEMARKERTEST']+ ' ' + kinshipcommand + ' -out ' + options['OUTPREFIX'] + '.singlemarker --run 1 > /dev/null'
		print(epactscommand)
		LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + epactscommand + "\n")
		os.system(epactscommand)
		
		#read annot file if it exists
		annotations = dict()
		if ('ANNOTCOLUMNS' in options) and ('ANNOTFILE' in options):
			annotation = pandas.read_table(defaults['INPUTDIR'] + '/' + defaults['ANNOTFILE'],compression="gzip")
			annotation = annotation.rename(columns={annotation.columns[4]:"GROUP"})
		
		#read single marker epacts results
		print("Processing results...")
		LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\tProcessing results...\n")
		singlemarkervcfs = pandas.read_table(options['OUTPREFIX'] + '.singlemarkers.vcf.gz', compression="gzip")
		singlemarkervcfs = singlemarkervcfs.rename(columns={singlemarkervcfs.columns[0]:'CHROM'})
		singlemarkervcfs['INDEX'] = singlemarkervcfs['CHROM'].map(str) + ":" + singlemarkervcfs['POS'].map(str)
		singlemarkervcfs = singlemarkervcfs.set_index('INDEX')
		del singlemarkervcfs['CHROM']
		del singlemarkervcfs['POS']
		del singlemarkervcfs['REF']
		del singlemarkervcfs['ALT']
		del singlemarkervcfs['FILTER']
		del singlemarkervcfs['QUAL']
		del singlemarkervcfs['ID']
		del singlemarkervcfs['FORMAT']
		del singlemarkervcfs['INFO']

		singlemarkeroutput = pandas.read_table(options['OUTPREFIX'] + '.singlemarker.epacts.gz', compression="gzip")
		singlemarkeroutput['MARKER_ID'] = singlemarkeroutput['MARKER_ID'].str.replace("_.\/._","|")
		singlemarkerresults = pandas.DataFrame(singlemarkeroutput.MARKER_ID.str.split('|').tolist(),columns="MARKER_ID MARKER".split())
		
		markeranno = singlemarkerresults[['MARKER', 'MARKER_ID']]
		singlemarkerresults['MAF'] = singlemarkeroutput['MAF']
		singlemarkerresults['PVALUE'] = singlemarkeroutput['PVALUE']
		singlemarkerresults['BETA'] = singlemarkeroutput['BETA']
		singlemarkerresults['MAC'] = singlemarkeroutput['NS']*2*singlemarkeroutput['MAF']
		markerstoadd = []
		pvaluestoadd = []
		betastoadd = []
		genestoadd = []
		individualstoadd = []
		annottoadd = []
		genotypestoadd = []
		mafstoadd = []
		macstoadd = []
		num = 0
		macsofgenes = dict()
		for eachgeneindex in genes.index:
			gene = genes[eachgeneindex]
			genename = gene.split('_')[1]
			mac = 0
			if genename not in markerlistforgenes:
				continue
			for marker in markerlistforgenes[genename]:
				markerresult = singlemarkerresults[singlemarkerresults.MARKER_ID == marker]
				mac = mac + markerresult.MAC[markerresult.MAC.index[0]]
				
				if (len(markerresult) > 0) and not (pandas.isnull(markerresult.PVALUE[markerresult.index[0]])):
					genotypes = singlemarkervcfs.loc[marker]
					genotypes = genotypes.str.replace(':.*','')
					genotypes = genotypes[genotypes != 0]
					if len(genotypes) > 0:
						for sample in range(0,len(genotypes)):
							individualstoadd.append(genotypes.index[sample])
							genotypestoadd.append(genotypes[sample])
							genestoadd.append(gene)
							markerstoadd.append(markerresult['MARKER_ID'][markerresult.index[0]])
							pvaluestoadd.append(markerresult['PVALUE'][markerresult.index[0]])
							betastoadd.append(markerresult['BETA'][markerresult.index[0]])
							mafstoadd.append(markerresult['MAF'][markerresult.index[0]])
							macstoadd.append(markerresult['MAC'][markerresult.index[0]])
					num = num + 1
				else:
					donothing = 1
			macsofgenes[gene] = mac

		variants = pandas.DataFrame(genestoadd)
		variants['MARKER'] = markerstoadd
		variants['PVALUE'] = pvaluestoadd
		variants['GENOTYPE'] = genotypestoadd
		variants['IND'] = individualstoadd
		variants['BETA'] = betastoadd
		variants['MAF'] = mafstoadd
		variants['MAC'] = macstoadd
		variants = variants.rename(columns={0:"GROUP"})
		allvariants = variants
	#	allvariants.to_csv(options['OUTPREFIX'] + ".allvariants.txt",sep="\t",index=False, index_label=False)
		#filter by mac of gene
		geneswithvalidmac = variants[['GROUP']]
		for macind in geneswithvalidmac.index:
			a = geneswithvalidmac.GROUP.iloc[macind]
			if macsofgenes[a] >= int(GENEMINMAC):
				geneswithvalidmac.iloc[macind] = 1
			else:
				geneswithvalidmac.iloc[macind] = 0
		
		variants = variants[geneswithvalidmac.GROUP == 1]
		#now you have genes filtered by minMAC
		
		variants = variants.rename(columns={"MARKER":"MARKER_ID"})
		#filter by numver of variants per gene
		
		variants.to_csv(options['OUTPREFIX'] + ".variants.txt",sep="\t", index=False, index_label=False)
		
		if 'ANNOTCOLUMNS' in options:
			columns = options['ANNOTCOLUMNS'].split(",")
			annonum = 0
			colnum = 0
			for col in annotation.columns:
				if ((col not in columns)) and (colnum > 0):
					del annotation[col]
				colnum = colnum + 1
			merged1 = pandas.merge(markeranno, annotation, left_on='MARKER', right_on = annotation.columns[0], how="left")
			del merged1[annotation.columns[0]]
			merged2 = pandas.merge(allvariants, merged1, left_on = 'MARKER', right_on='MARKER_ID',how="left")
			del merged2['MARKER_x']
			merged2 = merged2.rename(columns={"MARKER_y":"MARKERNAME"})
			allvariants = merged2
		allvariants = allvariants.rename(columns={"MARKER":"MARKER_ID"})
		allvariants.to_csv(options['OUTPREFIX'] + ".allvariants.txt",sep="\t",index=False, index_label=False)		
	
	LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + "Creating output plots for each significant gene..." + "\n")
	SUMMARYFILE = open(options['OUTPREFIX'] + ".summary.txt","w")
	SUMMARYFILE.write("DATE\t" + str(time.asctime(time.localtime(time.time()))) + "\n")
	SUMMARYFILE.write("MODEL\t" + options['MODEL'] + "\n")
	SUMMARYFILE.write("VCFFILE\t" + options['VCFFILE'] + "\n")
	SUMMARYFILE.write("PEDFILE\t" + options['PEDFILE'] + "\n")
	teststowrite = ''
	for test in TESTS:
		SUMMARYFILE.write("TEST\t" + test + "\n")
		teststowrite = teststowrite + test.split('=')[1] + ','
	if 'FILTERMAF' in options:
		SUMMARYFILE.write("FILTERMAF\t" + options['FILTERMAF'] + "\n")
	if 'ANNOTFILE' in options:
		SUMMARYFILE.write("ANNOTFILE\t" + options['ANNOTFILE'] + "\n")
	if 'FILTERANNOT' in options:
		SUMMARYFILE.write("FILTERANNOT\t" + options['FILTERANNOT'] + "\n")
	if 'PVALUETHRESHOLD' in options:
		SUMMARYFILE.write("PVALUETHRESHOLD\t" + options['PVALUETHRESHOLD'] + "\n")
	if 'ANNOTCOLUMNS' in options:
		SUMMARYFILE.write("ANNOTCOLUMNS\t" + options['ANNOTCOLUMNS'] + "\n")
	if 'PVALUETHRESHOLD' in options:
		SUMMARYFILE.write("PVALUETHRESHOLD\t" + options['PVALUETHRESHOLD'] + "\n")
	if 'MINVARS' in options:
		SUMMARYFILE.write("PLOT.MINVARS\t" + options['MINVARS'] + "\n")
	if 'GENEMINMAC' in options:
		SUMMARYFILE.write("PLOT.GENEMINMAC\t" + options['GENEMINMAC'] + "\n")
	
	SUMMARYFILE.close()
	teststowrite = teststowrite.strip(",")
	
	if len(covariates) > 0:
		os.system("R --vanilla --slave --args "+ options['OUTPREFIX'] + '.variants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(PVALUETHRESHOLD) + " " + options['MODEL'].split('~')[1] + ' < lattice.multiplegenes2.R')
		LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + "R --vanilla --slave --args "+ options['OUTPREFIX'] + '.variants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(PVALUETHRESHOLD) + ' "' + options['MODEL'].split('~')[1] + '" < lattice.multiplegenes2.R')
		print("R --vanilla --slave --args "+ options['OUTPREFIX'] + '.variants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(PVALUETHRESHOLD) + ' "' + options['MODEL'].split('~')[1] + '" < lattice.multiplegenes2.R')
	else:
		os.system("R --vanilla --slave --args "+ options['OUTPREFIX'] + '.variants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(PVALUETHRESHOLD) + " " + 'NA' + ' < lattice.multiplegenes2.R')
		LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + "R --vanilla --slave --args "+ options['OUTPREFIX'] + '.variants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(PVALUETHRESHOLD) + ' ' + 'NA '+ ' < lattice.multiplegenes2.R' + "\n")
		print("R --vanilla --slave --args "+ options['OUTPREFIX'] + '.variants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(PVALUETHRESHOLD) + ' ' + 'NA '+ ' < lattice.multiplegenes2.R')
	############
	line = configfile.readline()
	line_num = line_num + 1

LOGFILE.close()

command = 'mv ' + logfilename + ' ' + options['OUTPREFIX'] + '.log'
os.system(command)
print("Output written to " + options['OUTPREFIX'] + '.topgenes.pdf')
