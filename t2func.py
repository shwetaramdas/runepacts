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

######
def set_default(options = 'NA'):
	#Set defaults
	if options == 'NA':
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
		defaults['VCFDIR'] = None
		defaults['MARKERMINMAF'] = 0
		defaults['MARKERMAXMAF'] = 1
		defaults['MARKERMINMAC'] = 1
		defaults['MINVARS'] = 1
		return defaults
	else:
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
			check_file_exists(os.path.join(options['INPUTDIR'],options['KINSHIPFILE']))
		if 'INPUTDIR' in options:
			defaults['INPUTDIR'] = options['INPUTDIR']
		if 'SINGLEMARKERTEST' in options.keys():
			defaults['SINGLEMARKERTEST'] = options['SINGLEMARKERTEST']
		if 'FILTERMAF' in options.keys():
			defaults['FILTERMAF'] = options['FILTERMAF']
		if 'PVALUETHRESHOLD' in options.keys():
			defaults['PVALUETHRESHOLD ']= options['PVALUETHRESHOLD']
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
		if 'MARKERMINMAF' in options.keys():
			defaults['MARKERMINMAF'] = float(options['MARKERMINMAF'])
		if 'MARKERMINMAC' in options.keys():
			defaults['MARKERMINMAC'] = int(options['MARKERMINMAC'])
		
######
def get_vcf_header(vcf):
	if vcf.endswith(".gz"):
		f = gzip.open(vcf);
	else:
		f = open(vcf);
		
	header = None;
	with f:
		for line in f:
			if line.startswith("#CHROM"):
				header = line.rstrip().split("\t");
				break;
				
	return header;

#####
def calculatemaf(inputvcf,outprefix,samplestokeep,sepchr=False):
#	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Starting to calculate MAF\n")
	header = get_vcf_header(inputvcf);
		#get the header from the vcf and figure out sample columns to keep

	colstokeep = []
	samplestokeeplist = list(samplestokeep)
	for columnnum in range(10,len(header)):
		if header[columnnum] in samplestokeeplist:
			colstokeep.append(columnnum)	
	if inputvcf.endswith(".gz"):
		f = gzip.open(inputvcf);
	else:
		f = open(inputvcf);

	rows = []
	n_alleles = 2
	for line in f:
		if line.startswith("#"):
			continue;

		# Split the line on tabs (VCF must be tab-delimited)
		ls = line.split("\t");
		ls[-1] = ls[-1].rstrip();

		# Find which element of the genotype field is the genotype itself		
		gt_ind = 0;
		
		# Grab genotypes for only those samples we care about.
		# gt_ind is the index within the genotype field for the genotype itself
		# Remember sometimes that field can contain dosage, genotype likelihoods, etc.
		# Example: GT:EC:DS 0/0:1.3314141:8
		genos = [ls[i].split(":")[gt_ind] for i in colstokeep];

		# Count alleles per genotype.
		count = Counter();
		map(count.update,genos);
		n_chr = count['0'] + count['1'];
		freq_0 = count['0'] / float(n_chr);
		freq_1 = count['1'] / float(n_chr);
		mac = min(count['0'],count['1']);
		maf = min(freq_0,freq_1);

		# Store information on this variant
		chr, pos, id, ref, alt = ls[0:5];

		rows.append([chr,pos,id,ref,alt,n_alleles,n_chr,freq_0,freq_1,mac,maf]);
#	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Calculated maf!\n")	
	df = pandas.DataFrame(
		rows, columns=["CHROM","POS","ID","REF","ALT","N_ALLELES","N_CHR","FREQ0","FREQ1","MAC","MAF"]
	);
	df.to_csv(outprefix + ".maffile.txt",sep="\t",index=False,index_label=False)
#	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Variant Frequency file created" + "\n")
	return df
	
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
def create_group_file(outprefix, vcffilename, annofile, maffilterfile, maffilter, DELIM_MAF, DELIM_ANNO, samplestokeep,annotvarcol=0, annotgenecol=4, annotposcol= 1, sepchr=''):
	vcfname = vcffilename
	vcfdir = vcffilename.split('/')
	vcfdir = '/'.join(vcffilename.split('/')[0:(len(vcfdir)-1)])
	outfilename = outprefix + '_groups.txt'
	trimmedvcf = outprefix + vcfname.replace('.vcf.gz', '').replace(vcfdir + '/',"")+'_trimmed.txt'
	tofilter = False

	#make a bed file containing only the positions you need
	genes = []
	snps = []
	numline = 0
	bedfile = open(outprefix + vcfname.replace('.vcf.gz', '').replace(vcfdir + '/', '') + '_trimmed.bed',"w")

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
		if (not tofilter) or (temp[annotvarcol] in passsnps):
			genes.append(temp[annotgenecol])
			snps.append(temp[annotvarcol])
			if '-' in temp[annotposcol]:
				bedfile.write(temp[annotposcol].split(':')[0] + "\t" + temp[annotposcol].split(':')[1].split('-')[0] + "\t" + temp[annotposcol].split(':')[1].split('-')[1] + "\n")
			else:
				bedfile.write(temp[annotposcol].split(':')[0] + "\t" + temp[annotposcol].split(':')[1] + "\t" + str(int(temp[annotposcol].split(':')[1])+1) + "\n")
	filename.close()
	bedfile.close()
	tabixcommand = 'tabix ' + vcfname + ' -B ' + outprefix + vcfname.replace('.vcf.gz', '').replace(vcfdir + '/', '') + '_trimmed.bed | grep -v "#" > ' + trimmedvcf 
	os.system(tabixcommand)
	
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
#	LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "SNP info obtained\n")
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

######
def check_options(user_option):
#Check right combination of user options is provided
#	required['group'] = ['VCF', 'PED', 'OUT','GROUPF','PHENO']
#	required['single'] = ['VCF', 'PED', 'OUT','GROUPF']
	return 0
	
######
def check_file_exists(filename):
	if os.path.exists(filename):
		return 0
	else:
		print("Error: File " + filename + " does not exist\n")
		sys.exit()
##END FUNCTIONS
