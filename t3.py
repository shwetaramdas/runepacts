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
from t2func import *
from collections import deque, Counter

def main():
	DELIM_MAF = "\t"
	DELIM_ANNO = "\t"
	
	#if no options are specified, then print out message
	if not len(sys.argv) > 1:
		print("runepacts.py\n" + "Option: --config required\n" + "Usage: python runepacts.py --config config.txt\n")
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
	print(lines_in_file)
	line_num = 0
	line = configfile.readline()

	#set default parameters to be used
	defaults = set_default()
	
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
			print(line)
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
		print("Running process "+ str(numtests + 1))
		LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Running process " + str(numtests) + "\n")
		numtests = numtests + 1
		if len(TESTS) == 0:
			print("Error: Must specify test.")
			sys.exit()
		
		#if folder to write files in does not exist, make it
		if not os.path.isdir(re.sub('/.*', '', options['OUTPREFIX'])):
			os.system('mkdir ' + re.sub('/.*', '', options['OUTPREFIX']))		

		#set the default parameters as the input options
		defaults['VCFFILE'] = options['VCFFILE']

		if '/' in options['VCFFILE']:		#vcffile is in a different folder than the inputdir
			vcfdir = options['VCFFILE'].split('/')
			defaults['VCFDIR'] = "/".join(vcfdir[0:(len(vcfdir)-1)])
			options['VCFDIR'] = defaults['VCFDIR']
			vcffilename = options['VCFFILE']
		else:
			vcffilename = os.path.join(options['INPUTDIR'], options['VCFFILE'])
			
########## Setting defaults to input parameters

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
			if options['FILTERMAF'] == 'NA':
				defaults['FILTERMAF'] = 0.05
			else:
				defaults['FILTERMAF'] = options['FILTERMAF']
		if 'PVALUETHRESHOLD' in options.keys():
			defaults['PVALUETHRESHOLD ']= options['PVALUETHRESHOLD']
		else:
			PVALUETHRESHOLD = 0.05
		if 'MINMAF' in options.keys():
			defaults['MINMAF'] = options['MINMAF']
		if 'MINVARS' in options.keys():
			defaults['MINVARS'] = options['MINVARS']
		else:
			defaults['MINVARS'] = 2
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

		kinshipcommand = ''

		if 'EPACTSDIR' not in options.keys():
			epacts = 'epacts'
		else:
			epacts = options['EPACTSDIR'] + '/epacts'
			defaults['EPACTSDIR'] = options['EPACTSDIR'] + '/epacts'

		OUTDIR = re.sub(".*/",'',options['OUTPREFIX'])

#############################################
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
		
		if 'FILTERMAF' in options.keys():
			maxmaf = ' -max-maf ' + str(options['FILTERMAF'])
		else:
			maxmaf = ' -max-maf ' + str(defaults['FILTERMAF'])

		#for single marker test options
		if 'MARKERMINMAF' in options.keys():
			markerminmaf = ' -min-maf ' + str(options['MARKERMINMAF'])
		else:
			markerminmaf = ' -min-maf ' + str(defaults['MARKERMINMAF'])

		if 'MARKERMINMAC' in options.keys():
			markerminmac = ' -min-mac ' + str(options['MARKERMINMAC'])
		else:
			markerminmac = ' -min-mac ' + str(defaults['MARKERMINMAC'])

		if 'MARKERMAXMAF' in options.keys():
			markermaxmaf = ' -max-maf ' + str(options['MARKERMAXMAF'])
		else:
			markermaxmaf = ' -max-maf ' + str(defaults['MARKERMAXMAF'])

		phenotype = options['MODEL'].split('~')[0]
		#################################################

		#First step: Index vcffile if indexed file does not already exist
		if not os.path.exists(vcffilename + '.tbi'):
			print("Indexing vcf. Warning: This may take longer for large VCFs\n")
			LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Indexing vcf\n")
			tabix_command = 'tabix -pvcf -f ' + vcffilename + " > /dev/null"
			os.system(tabix_command)
			LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + tabix_command + "\n")
		
		#if sepchr is on, then tabix all other chromosomes too
		if ('SEPCHR' in options) and (options['SEPCHR'] == 'ON'):
			vcffilenametomatch = vcffilename.split('/')[len(vcffilename.split('/'))-1].split('chr1')
			for f in os.listdir(options['VCFDIR']):
				if '.gz.tbi' in f:
					continue
				match = True
				for f2 in vcffilenametomatch:
					if f2 not in f:
						match = False
						break
				if match:
					tabix_command = 'tabix -pvcf -f ' + os.path.join(options['VCFDIR'],f) + " > /dev/null"
					os.system(tabix_command)
	
		#Then, make a ped file after using the filters provided
		print "Creating ped file...\n"
		LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Creating ped file...\n")

		pedfile = pandas.read_table(os.path.join(defaults['INPUTDIR'],defaults['PEDFILE']) , sep="\t",header=0)
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
			
			#no annotation file provided by the user, so we have to use epacts to create the groupfile
			if 'ANNOTFILE' not in options:
				if 'SEPCHR' not in options:
					annotatecommand = epacts + ' anno -in ' + vcffilename + ' -out ' + options['OUTPREFIX'] + '.anno.vcf.gz'
					os.system(annotatecommand)
					annotatecommand = epacts + ' make-group --vcf ' + options['OUTPREFIX'] + '.anno.vcf.gz --out ' + options['OUTPREFIX'] + '.annotated.txt --format epacts --nonsyn'
					print(annotatecommand)
					os.system(annotatecommand)
					finalgroupfilename = options['OUTPREFIX'] + '.annotated.txt'
				
				#SEPCHR, so we annotate each chromosome separately
				else:	
					vcffilenametomatch = options['VCFFILE'].split('/')[len(options['VCFFILE'].split('/')) - 1].split('chr1')
					for f in os.listdir(options['VCFDIR']):
						match = True
						for f2 in vcffilenametomatch:
							if f2 not in f:
								match = False
								break
						if match:
							annotatecommand = epacts + ' anno -in ' + os.path.join(defaults['VCFDIR'],f) + ' -out ' + options['OUTPREFIX'] + '.' + f
							os.system(annotatecommand)
							annotatecommand = epacts + ' make-group --vcf ' + options['OUTPREFIX'] + '.' + f + ' --out ' + options['OUTPREFIX'] + '.annotated.' + f  + ' --format epacts'
							os.system(annotatecommand)
							annotatecommand = 'cat ' + options['OUTPREFIX'] + '.' + f + ' >> ' + options['OUTPREFIX'] + '.annotated.txt'
							os.system(annotatecommand)
					finalgroupfilename = options['OUTPREFIX'] + '.annotated.txt' 
			
			#annotation file is in options 
			else:
				if '.gz' in options['ANNOTFILE']:
					cat = 'zcat '
				else:
					cat = 'cat '
				
				if 'FILTERANNOT' in options.keys():
					print "Filtering annotations...\n"
					LOGFILE.write("Filtering annotations...\n")
					filter = options['FILTERANNOT']
					annofile = pandas.read_table(os.path.join(defaults['INPUTDIR'],options['ANNOTFILE']), header=0,sep="\t",compression='gzip')
					annoheadersplit = annofile.columns
					
					annorowstokeep = evalexpression(filter, annofile)
					annofile = annofile.iloc[list(annorowstokeep),:]
					annofile.to_csv(os.path.join(defaults['INPUTDIR'], options['ANNOTFILE'].replace('.gz', '') + '.filtered.txt'),sep="\t", index=False)
					create_group_file(options['OUTPREFIX'], vcffilename, defaults['INPUTDIR'] + '/' + options['ANNOTFILE'].replace('.gz', '') + '.filtered.txt', 'NA', defaults['FILTERMAF'], DELIM_MAF, DELIM_ANNO,samplestokeep,defaults['ANNOTVARCOL'],defaults['ANNOTGENECOL'],defaults['ANNOTPOSCOL'],sepchr)
					
				else:
					create_group_file(options['OUTPREFIX'], vcffilename, os.path.join(defaults['INPUTDIR'], defaults['ANNOTFILE']), 'NA', defaults['FILTERMAF'], DELIM_MAF, DELIM_ANNO,samplestokeep,defaults['ANNOTVARCOL'],defaults['ANNOTGENECOL'],defaults['ANNOTPOSCOL'],sepchr)
				finalgroupfilename = options['OUTPREFIX'] + '_groups.txt'
		else:
			finalgroupfilename = os.path.join(defaults['INPUTDIR'], defaults['GROUPFILE'])

		if 'GENELIST' in options:
			outputname = random.randrange(1,1000000)
			getonlygenelistcommand = 'cat ' + finalgroupfilename + ' | grep -wf ' + os.path.join(options['INPUTDIR'], options['GENELIST']) + ' > ' + outputname
			os.system(getonlygenelistcommand)
			os.system('mv ' + outputname + ' ' + finalgroupfilename)		
		
		############
	
		LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Epacts run! Formatting output..." + "\n")		

		#Create kinship matrix if needed
		for TEST in TESTS:
			#Check if test includes an emmax test. Then check if kinship file has been provided. If not, then make one.
			if ('emmax' in TEST.split('=')[1]) or ('emmaxVT' in TEST.split('=')[1]) or ('emmaxCMC' in TEST.split('=')[1])  or ('emmaxSKAT' in TEST.split('=')[1]) or ('mmskat' in TEST.split('=')[1]) or ('q.emmax' in options['SINGLEMARKERTEST']):
				if 'KINSHIPFILE' in options.keys():
					if not os.path.exists(os.path.join(defaults['INPUTDIR'],defaults['KINSHIPFILE'])):
						print("Kinship file specified does not exist")
						LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + "Kinship file specified does not exist\n")
						sys.exit()
					kinshipcommand = ' -kin ' + os.path.join(defaults['INPUTDIR'],defaults['KINSHIPFILE']) + ' '
				else:
					#make kinship file from input using emmax
					makekinshipcommand = epacts + ' make-kin '  + ' -vcf ' + vcffilename \
											+ ' -ped ' + options['OUTPREFIX'] + '.pheno.ped ' + ' -out ' + options['OUTPREFIX'] \
											+ '.kinf' \
											+ ' -min-maf 0.01 -min-callrate 0.95' + ' -run 1 > /dev/null'
					#print makekinshipcommand
					os.system(makekinshipcommand)
					LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + makekinshipcommand + "\n")
					kinshipcommand = ' -kin ' + options['OUTPREFIX']  + '.kinf '
		

		#Run single marker epacts test
		print("Running Single Marker EPACTS test...")
		epactscommand = epacts + ' single --vcf ' + vcffilename + sepchr + ' -ped ' + options['OUTPREFIX'] + '.pheno.ped ' + '-pheno ' + phenotype + ' ' + covariatecommand + ' -test ' + defaults['SINGLEMARKERTEST']+ ' ' + kinshipcommand + markerminmaf + ' ' + markermaxmaf + ' ' + markerminmac +  ' -out ' + options['OUTPREFIX'] + '.singlemarker --run 1 > /dev/null'
		print(epactscommand)
		LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + epactscommand + "\n")
		os.system(epactscommand)		

		#read in single marker results
		singlemarkeroutput = pandas.read_table(options['OUTPREFIX'] + '.singlemarker.epacts.gz', compression="gzip")
		#singlemarkeroutput['MARKER_ID'] = singlemarkeroutput['MARKER_ID'].str.replace("_.\/._","|")
		singlemarkerresults = pandas.DataFrame(singlemarkeroutput.MARKER_ID.str.split('_').tolist())
		if len(singlemarkerresults.columns) > 3:
			ncol = len(singlemarkerresults.columns)
			for i in range((ncol-1),2,-1):
				singlemarkerresults[singlemarkerresults.columns[2]] = singlemarkerresults[singlemarkerresults.columns[2]].map(str) + "_" + singlemarkerresults[singlemarkerresults.columns[i]].map(str)
				singlemarkerresults[singlemarkerresults.columns[2]] = singlemarkerresults[singlemarkerresults.columns[2]].str.replace("_None","")
				del singlemarkerresults[singlemarkerresults.columns[i]]
			
		singlemarkerresults[singlemarkerresults.columns[0]] = singlemarkerresults[singlemarkerresults.columns[0]].map(str) + "_" + singlemarkerresults[singlemarkerresults.columns[1]].map(str)
		del singlemarkerresults[singlemarkerresults.columns[1]]
		
		
		singlemarkerresults.columns = "MARKER_ID MARKER".split()
		
		markeranno = singlemarkerresults[['MARKER', 'MARKER_ID']]
		singlemarkerresults['MAF'] = singlemarkeroutput['MAF']
		singlemarkerresults['PVALUE'] = singlemarkeroutput['PVALUE']
		singlemarkerresults['BETA'] = singlemarkeroutput['BETA']
		singlemarkerresults['MAC'] = singlemarkeroutput['NS']*2*singlemarkeroutput['MAF']
		
		#get variants that pass the MAF/MAC filter
		passsnps = singlemarkerresults[singlemarkerresults['MAF'] < float(defaults['FILTERMAF'])]
		passsnps = passsnps[passsnps['MAF'] >= defaults['MINMAF']]
		passsnps = passsnps[passsnps['MAC'] >= defaults['MARKERMINMAC']]
		macs = passsnps[['MARKER_ID', 'MAC']].drop_duplicates()
		passsnps = list(passsnps['MARKER_ID'])
		pandas.DataFrame(passsnps).to_csv("PASSSNPS.txt",sep="\t",index=False,index_label=False)
		
		#get all markers belonging to a gene
		groupfile = open(finalgroupfilename)
		newgroupfile = open(finalgroupfilename + '.temp', 'w')
		markerlistforgenes = dict()
		macsofgenes = dict()
		genespassingfilters = dict()
		markernames = []
		
		print("Reading groupfile\n")
		#read in groupfile, and get the list of variants belonging to each gene. Store this in markerlistforgenes.
		#Also get the genemac and the numvars for each gene and get the genespassingfilters
		for line in groupfile:
			line = line.rstrip()
			temp = line.split()
			originaltemp = line.split()
			genename = temp[0]
			temp = temp[1:]
			originaltemp = originaltemp[1:]
			numvars = 0
			mac = 0 
			
			towrite = genename
			#iterating over all markers in gene, find if gene passes filters. also get 
			for tempindex in range(0,len(temp)):
				original = temp[tempindex]
				#temp[tempindex] = re.sub(r"_.*","",temp[tempindex])
				if temp[tempindex] in passsnps:
					towrite = towrite + "\t" + original
					numvars = numvars + 1
					mac = mac + float(macs[macs['MARKER_ID'] == temp[tempindex]]['MAC'])
			
			if (mac >= int(defaults['GENEMINMAC'])) and (numvars >= int(defaults['MINVARS'])):
				genespassingfilters[genename] = 1
				for tempindex in range(0,len(temp)):
					if temp[tempindex] in passsnps:
						markernames.append(originaltemp[tempindex])
			
			newgroupfile.write(towrite + "\n")
			markerlistforgenes[genename] = temp
			macsofgenes[genename] = mac
			
		groupfile.close()
		newgroupfile.close()
		os.system('mv ' + finalgroupfilename + '.temp ' + finalgroupfilename)
		genespassingfilters = pandas.DataFrame(genespassingfilters.keys())
		genespassingfilters.to_csv(options['OUTPREFIX'] + '.genespassingfilters.txt',index=False,index_label=False)
		
		for TEST in TESTS:
			#now create the epacts command
			epacts_command = epacts + ' ' + TEST.split('=')[0] + ' -test ' + TEST.split('=')[1] + ' -vcf ' \
											+ vcffilename + ' -pheno ' + phenotype + ' -ped ' \
											+ options['OUTPREFIX'] + '.pheno.ped ' + ' -groupf ' \
											+ finalgroupfilename + ' -out ' + options['OUTPREFIX'] + '.' + TEST.split('=')[1] \
											+ minmaf + ' ' + markerminmac + ' ' + covariatecommand + usercommand + sepchr + kinshipcommand + ' -run 1  >/dev/null'

			if ('RUNEPACTS' in options.keys()) and (options['RUNEPACTS'] == 'FALSE'):
				print epacts_command
				LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + epacts_command + "\n")
			else:
				print epacts_command
				os.system(epacts_command)
				LOGFILE.write(str(time.asctime( time.localtime(time.time()))) + "\t" + epacts_command + "\n")
				print("Epacts run! \n")
		
		#reading in group test results
		for TEST in TESTS:
			outputtest = pandas.read_table(options['OUTPREFIX'] + '.' + TEST.split('=')[1] + '.epacts',sep="\t",header=0)	#reading in output file for test
			if 'NUM_PASS_VARS' in outputtest.columns:
				outputtest.rename(columns={'NUM_PASS_VARS':'PASS_MARKERS'}, inplace=True)
			if 'BEGIN' not in outputtest.columns:
				outputtest.rename(columns={'BEG':'BEGIN'}, inplace=True)
			
			outputtest = outputtest.sort(['#CHROM','BEGIN'])
			if TEST == TESTS[0]:
				output = outputtest
			else:
				output = output.append(outputtest)
		
		#the output dataframe contains concatenated results from all EPACTS tests run
		output = output.sort(columns=[output.columns[0], 'BEGIN'],inplace=False)
		allgenesoutput = output
		allgenes = allgenesoutput['MARKER_ID'].drop_duplicates()
		numvars = output[['MARKER_ID','PASS_MARKERS']]
		numvars = numvars.drop_duplicates()

		output = output[output.PVALUE <= float(defaults['PVALUETHRESHOLD'])]
		genes = output['MARKER_ID'].drop_duplicates()
		
		#if output is empty, then ??
		if len(output) == 0:
			print("Epacts groupwise test returned no significant results")
			continue		
		
		markernames = pandas.DataFrame(markernames)
		markernames[markernames.columns[0]] = markernames[markernames.columns[0]].str.replace("_",":")
		markernames = pandas.DataFrame(markernames[markernames.columns[0]].str.split(":").tolist())
		markernames[markernames.columns[1]] = markernames[markernames.columns[1]].astype(int)
		
		markernameswithalleles = markernames
		del markernames[markernames.columns[2]]
		
		markernames['END'] = markernames[markernames.columns[1]].astype(int) + numpy.repeat(1,len(markernames))
		markernames[markernames.columns[1]] = markernames[markernames.columns[1]].astype(int) - numpy.repeat(1,len(markernames))
		markernames = markernames.sort(columns=[markernames.columns[0], markernames.columns[1]])
		markernames.to_csv(options['OUTPREFIX'] + '.markersfromsignificantgenes.txt',sep="\t",index=False,index_label=False)
		
		markernames = markernameswithalleles
		
		#now extract these markers from the full vcf
		singlemarkervcfs = pandas.DataFrame()
		if 'SEPCHR' not in options:
#			tabixcommand = "vcftools --gzvcf " + vcffilename + ' --bed ' + options['OUTPREFIX'] + '.markersfromsignificantgenes.txt --recode --out ' + options['OUTPREFIX'] + '.markersfromsignificantgenes.temp'
#			os.system(tabixcommand)
#			tabixcommand = 'cat ' + options['OUTPREFIX'] + '.markersfromsignificantgenes.temp.recode.vcf | grep -v "##" > ' + options['OUTPREFIX'] + '.markersfromsignificantgenes.recode.vcf'
#			os.system(tabixcommand)
			tabixcommand = 'tabix ' + vcffilename + ' -B ' + options['OUTPREFIX'] + '.markersfromsignificantgenes.txt | grep -v "##" > ' + options['OUTPREFIX'] + '.markersfromsignificantgenes.recode.vcf'
			os.system(tabixcommand)
			singlemarkervcfs = pandas.read_table(options['OUTPREFIX'] + '.markersfromsignificantgenes.recode.vcf',header=None)
		else:
			vcffilenametomatch = vcffilename.split('/')[len(vcffilename.split('/'))-1].split('chr1')
			for f in os.listdir(options['VCFDIR']):
				if '.gz.tbi' in f:
					continue
				match = True
				for f2 in vcffilenametomatch:
					if f2 not in f:
						match = False
						break
				if match:
					tabixcommand = 'tabix ' + os.path.join(options['VCFDIR'], f) + ' -B ' + options['OUTPREFIX'] + '.markersfromsignificantgenes.txt | grep -v "##" > ' + options['OUTPREFIX'] + '.markersfromsignificantgenes.recode.vcf'
					os.system(tabixcommand)
					singlemarkervcfs = singlemarkervcfs.append(pandas.read_table(options['OUTPREFIX'] + '.markersfromsignificantgenes.recode.vcf',header=None))
						
		header = subprocess.Popen('zcat '+ vcffilename +' | grep -m1 CHROM', shell=True,stdout=subprocess.PIPE)
		header = header.communicate()[0]
		singlemarkervcfs.columns = header.rstrip().split()
		
		singlemarkervcfs = singlemarkervcfs.rename(columns={singlemarkervcfs.columns[0]:'CHROM'})
		singlemarkervcfs = singlemarkervcfs.rename(columns={singlemarkervcfs.columns[1]:'POS'})
		singlemarkervcfs = singlemarkervcfs.rename(columns={singlemarkervcfs.columns[3]:'REF'})
		singlemarkervcfs = singlemarkervcfs.rename(columns={singlemarkervcfs.columns[4]:'ALT'})
		singlemarkervcfs['INDEX'] = singlemarkervcfs['CHROM'].map(str) + ":" + singlemarkervcfs['POS'].map(str) + "_" + singlemarkervcfs['REF'].map(str) + '/' + singlemarkervcfs['ALT'].map(str)
		singlemarkervcfs = singlemarkervcfs.set_index('INDEX')		

		for sample in range(len(singlemarkervcfs.columns)-1,0,-1):
			if singlemarkervcfs.columns[sample] not in list(samplestokeep):
				del singlemarkervcfs[singlemarkervcfs.columns[sample]]
		
		#now for each significant gene, get the genotypes for each sample, and 
		individualstoadd = []
		genotypestoadd = []
		genestoadd = []
		markerstoadd = []
		mafstoadd = []
		macstoadd = []
		betastoadd = []
		pvaluestoadd = []

		for geneindex in genes.index:
			gene = genes[geneindex]
			genename = gene.split('_')[1]
			
			if genename not in list(genespassingfilters[0]):
				continue

			for marker in markerlistforgenes[genename]:
				if marker not in singlemarkervcfs.index:
					continue

				markerresult = singlemarkerresults[singlemarkerresults.MARKER_ID == marker]
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
		
		
		#read annot file if it exists
		annotations = dict()
		if ('ANNOTCOLUMNS' in options) and ('ANNOTFILE' in options):
			annotation = pandas.read_table(defaults['INPUTDIR'] + '/' + defaults['ANNOTFILE'],compression="gzip")
			annotation = annotation.rename(columns={annotation.columns[4]:"GROUP"})
			
		if 'ANNOTCOLUMNS' in options:
			columns = options['ANNOTCOLUMNS'].split(",")
			annonum = 0
			colnum = 0
			for col in annotation.columns:
				if ((col not in columns)) and (colnum > 0):
					del annotation[col]
				colnum = colnum + 1
			merged1 = pandas.merge(markeranno, annotation, left_on='MARKER', right_on = annotation.columns[defaults['ANNOTVARCOL']], how="left")
			print(merged1)
			print(merged1.iloc[0])
			del merged1[annotation.columns[defaults['ANNOTVARCOL']]]
			merged2 = pandas.merge(allvariants, merged1, left_on = 'MARKER', right_on='MARKER_ID',how="left")
			print(merged2)
			print(merged1.iloc[0])
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
		if 'MINVARS' in options:
			SUMMARYFILE.write("PLOT.MINVARS\t" + options['MINVARS'] + "\n")
		if 'GENEMINMAC' in options:
			SUMMARYFILE.write("PLOT.GENEMINMAC\t" + options['GENEMINMAC'] + "\n")
		
		SUMMARYFILE.close()
		teststowrite = teststowrite.strip(",")
		
		if len(covariates) > 0:
			os.system("R --vanilla --slave --args "+ options['OUTPREFIX'] + '.allvariants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(defaults['PVALUETHRESHOLD']) + " " + options['MODEL'].split('~')[1] + ' < l3.R')
			LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + "R --vanilla --slave --args "+ options['OUTPREFIX'] + '.variants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(defaults['PVALUETHRESHOLD']) + ' "' + options['MODEL'].split('~')[1] + '" < l3.R')
			print("R --vanilla --slave --args "+ options['OUTPREFIX'] + '.allvariants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(defaults['PVALUETHRESHOLD']) + ' "' + options['MODEL'].split('~')[1] + '" < l3.R')
		else:
			os.system("R --vanilla --slave --args "+ options['OUTPREFIX'] + '.allvariants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(defaults['PVALUETHRESHOLD']) + " " + 'NA' + ' < l3.R')
			LOGFILE.write(str(time.asctime(time.localtime(time.time()))) + "\t" + "R --vanilla --slave --args "+ options['OUTPREFIX'] + '.allvariants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(defaults['PVALUETHRESHOLD']) + ' ' + 'NA '+ ' < l3.R' + "\n")
			print("R --vanilla --slave --args "+ options['OUTPREFIX'] + '.allvariants.txt ' + phenotype + ' ' + options['OUTPREFIX'] + ' ' + teststowrite + ' ' + str(defaults['PVALUETHRESHOLD']) + ' ' + 'NA '+ ' < l3.R')
		############


		line = configfile.readline()
		line_num = line_num + 1
		
	LOGFILE.close()
	command = 'mv ' + logfilename + ' ' + options['OUTPREFIX'] + '.log'
	os.system(command)
	print("Output written to " + options['OUTPREFIX'] + '.topgenes.pdf')
	
if __name__ == "__main__":
	main()
