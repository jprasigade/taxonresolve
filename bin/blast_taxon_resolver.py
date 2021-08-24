#!/usr/bin/env python3

"""
blast_overlap_resolver	 V1.0	 anais.barray@gmail.com

The aim of this program is to generate a matrix of read counts to be used for diversity & ecological analyses, based on BLAST results.
If a BLAST search returns several best hits for a read and cannot resolve the most likely taxon to assign, a group of species is formed. 
The tool resolves the eventual taxonomical overlaps between groups, by generating minimal groups based on their intersections.
It then computes the weight of minimal groups within all samples, to distribute more sparingly the read counts when they are shared between several minimal groups.


Input file is a manifest file listing all samples and their file path, with the following format :
   e.g. manifest.txt content : 
	 sample1:/path/to/file1.blast
	 sample2:/path/to/file2.blast
	 
The BLAST files must be generated respecting the following formatting options : 
   1) -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' 
   2) -max_target_seqs >= 1 (e.g. -max_target_seqs 20)

This script allows the possibility to dereplicate reads prior to BLAST in order to accelerate computing times. e.g. with vsearch --derep_fulllength --sizeout command.
 
#### ON TAXONOMY REFERENCE ######################################################################
#  This script makes use of the new_taxdump files from the NCBI, download beforehand with :
#     wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip & unzip new_taxdump.zip
#################################################################################################
 
#### ON 16S NORMALIZATION ######################################################################################################################## 
#  If dealing with 16S metabarcoding data, we suggest to normalize the input files' reads content based on the 16S copy numbers of the taxonomies
#  This script only requires the genus 16S copy number.
#  The 16S copy number NCBI database can be acquired and generated as follow :
#	  wget https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.5_pantaxa_stats_NCBI.tsv.zip & unzip rrnDB-5.5_pantaxa_stats_NCBI.tsv.zip
#	  head -n 1 rrnDB-5.5_pantaxa_stats_NCBI.tsv > rrnDB-5.5_pantaxa_stats_NCBI_genus.tsv
#	  grep "genus" rrnDB-5.5_pantaxa_stats_NCBI.tsv | grep -v "species" >> rrnDB-5.5_pantaxa_stats_NCBI_genus.tsv
#  Manual of the database content : https://rrndb.umms.med.umich.edu/help/#help-download-pan-taxa
#  Some database entries were manually updated with their new taxonomy or Blast database nomenclature (e.g. Actinomyces was duplicated as Schaalia, 
#  or manually added for missing bacteria, mocks and spikes (e.g. Dolosigranulum, Imtechella & Allobacillus)
##################################################################################################################################################

Calling example : 
blast_taxon_resolver.py -i manifest.txt -dr /path/to/rankedlineage.dmp -n yes -dn /path/to/16S_copynb_db
"""

import sys, argparse, os, re, string, math
import numpy as np
from operator import itemgetter
from joblib import Parallel, delayed

#### Parsing options
parser = argparse.ArgumentParser(description='Creates a Species contingency table from BLAST taxonomies obtained in metagenomics projects. Output displays the read counts of each Species (rows) for each sample (column). Bacterial read numbers can be normalized to represent Species counts instead, provided a 16S copy number database ; in such case both contingency tables will be produced. This algorithm aims to improve Species discrimination by resolving ambiguous BLAST taxonomic assignations.', add_help=True)

parser.add_argument('-i', metavar='Input sample list', type=str, dest="infile", required=True,
					help="Manifest file with each line = sample_name:/path/to/file.blast" )
parser.add_argument('-o', metavar='Output read count matrix', type=str, dest="countfile",
					help="Outputs a matrix of taxonomy counts found in each sample. Default : readcount_table.tsv", default="readcount_table.tsv")
parser.add_argument('-a', metavar='Min abundance threshold in percent', type=str, dest='abundthreshold', 
					help='Abundance threshold (INT percent) above which a taxonomy must be represented in at least one sample. Default : 0', required=False, default=0)
parser.add_argument('-id', metavar='Identity threshold', dest='id', type=str, help='BLAST Top hits identity treshold. Default : 90', required=False, default='90')
parser.add_argument('-cov', metavar='Coverage threshold', dest='cov', type=str,	help='BLAST Top hits coverage treshold. Default : 80', required=False, default='80')
parser.add_argument('-sep', metavar='Taxonomy separator', dest='sep', type=str,	help='Taxonomy format separator between quotes, e.g. : "\t", ";"', required=False, default=";")
parser.add_argument('-dr', metavar='Ranked lineages database', dest='rankedlineage', type=str, help='Reference json taxonomy file', required=True)
parser.add_argument('-n', metavar='Normalize 16S copy numbers', type=str, dest='normalize', 
					help='Normalize 16S copy numbers from read counts with the NCBI database. Default : no', required=False, choices=['no', 'yes'], nargs='?', default='no', const='yes')
parser.add_argument('-dn', metavar='Normalization 16S copy number database', type=str, dest="normdb", required=False, help="Path to the 16S copy number genus NCBI database")	
parser.add_argument('-t', metavar='Number of threads', type=int, dest='threads', 
					help='Number of CPU threads for parallelisation. Default : -1 (all available)', required=False, default="-1")
args = parser.parse_args()

#### Variables & checking input files/databases
input_file = args.infile
count_file = args.countfile
refdb = args.rankedlineage
abund_threshold = float(args.abundthreshold)
nb_threads = args.threads
norm_bool = True if args.normalize == "yes" else False
if args.sep == "\\t" : 
	sep = "\t"
else :
	sep = args.sep


if not os.path.isfile(input_file):
	print("ERROR : input file is not recognized, option -i\n")
	parser.print_help()
	sys.exit()
	
if not os.path.isfile(refdb):
	print("ERROR : database file is not recognized, option -dr\n")
	parser.print_help()
	sys.exit()
	
if norm_bool :
	db16S = args.normdb
	norm_countfile = "normalized_"+str(os.path.basename(count_file))
	if os.path.dirname(count_file) : norm_countfile = str(os.path.dirname(count_file)+"/"+norm_countfile) 
	if not db16S or not os.path.isfile(db16S):
		print("ERROR : No recognized input 16S copy number database, option -dn\n")
		parser.print_help()
		sys.exit()

def getDictFiles(F):
	# From the input file, returns a dictionnary of samples, where the keys are sample names, and values are path/to/file.

	dict_files = {}
	fofs = open(F,'r')
	
	while(True):
		line = fofs.readline()
		if not line: break	#EOF
		if not line.strip(): continue # avoid empty lines
		sample_name = line.split(":")[0]
		path = line.split(":")[1].strip()
		dict_files[sample_name] = path
	fofs.close()
	
	return dict_files


def linecount(sample_path):
	# Count number of lines in input
	
    i = 0
    with open(sample_path, "r") as s:
        for i, l in enumerate(s):
            pass
    return i


def get_bestHitTaxid(queryLines):
	# Get list of taxid of all best bitscores BLAST hits
	
	list_taxid = []
	best_bitscore = 0

	for i in range(0, len(queryLines)) :
		if float(queryLines[i][4]) >= float(args.id) and float(queryLines[i][5]) >= float(args.cov):
			try : 
				bitscore = int(queryLines[i][7])
				if bitscore >= best_bitscore :
					best_bitscore = bitscore
					taxid = queryLines[i][3]
					if taxid not in list_taxid : list_taxid.append(taxid) 
			except : # blast bitscore as float for some entries, but are always lower than best_bitscore
				pass

	return list_taxid
	
	
def getListTaxid(sample_path):
	# Returns a list of taxa found in all samples
	
	list_taxid = []
	try : 
		
		fofs = open(sample_path,'r')
		while(True):
			line = fofs.readline()
			if not line: break	#EOF
			if not line.strip(): continue # avoid empty lines
			taxid = line.split("\t")[3]
			if taxid not in list_taxid : list_taxid.append(taxid)
		fofs.close()

	except FileNotFoundError: 
		print("err:" , sample_path, "\n>FILE NOT FOUND< ... Provide correct file path and try again ... ")
		parser.print_help()
		sys.exit()
		
	return list_taxid
	
	
def getSampleNbReads(dictTaxidToCount, sample) :
	# Returns the total number of reads assigned to a taxonomy in a sample.

	sample_nbread = sum(dictTaxidToCount.values())
	
	return (sample_nbread, sample)
	

def getDictGeneraCopynb(db16S):
	""" 
	  Make a dictionary with key : genus, value : 16S copy number
	  Also capture the mean of 16S copy number to normalize unclassified or unlisted entries.
	"""
	
	dictGeneraCopynb = {}
	
	db = open(db16S,'r')
	db.readline() # skip header
	while(True):
		line = db.readline()
		if not line: break	#EOF
		if not line.strip(): continue # avoid empty lines
		list_line = line.split()
		genus = list_line[2]

		# Very rare cases where Genus is composed of 2 names (e.g. Candidatus Bipolaricaulis). Only the first name will be kept. 
		# If same first name is encountered again (e.g. Candidatus Nasuia), the old value is replaced by the new, but they are usually identical.
		check = list_line[3] 
		try:
			c = int(check)
		except ValueError:
			del list_line[3] 
			
		copynb = float(list_line[6])
		dictGeneraCopynb[genus] = copynb
	db.close()
	
	mean_norm = float(round(np.array(list(dictGeneraCopynb.values())).mean(),1))
	
	return dictGeneraCopynb, mean_norm


def getDictReferenceTaxonomy(list_unique_taxids):
	# Generate a dictionary of taxids (key) to their taxonomy (value in a dictionary format) based on NCBI's nomenclatures for the list of taxids found in all samples.
	
    dict_taxonomy = {}
    with open(refdb) as rankedlineage:
        for tax in rankedlineage:
            tax = tax.split("|")
            taxonid = tax[0].strip()
            if taxonid in list_unique_taxids :
	            species = tax[1].strip() if tax[1].strip() else "unknown species"
	            genus = tax[3].strip() if tax[3].strip() else "unknown genus"
	            family = tax[4].strip() if tax[4].strip() else "unknown family"
	            order = tax[5].strip() if tax[5].strip() else "unknown order"
	            classe = tax[6].strip() if tax[6].strip() else "unknown class"
	            phylum = tax[7].strip() if tax[7].strip() else "unknown phylum"
	            kingdom = tax[8].strip() if tax[8].strip() else "unknown kingdom"
	            superkingdom = tax[9].strip() if tax[9].strip() else "unknown superkingdom"
	            dict_taxonomy[taxonid] = {"species":species, "genus":genus, "family":family, "order":order, "class":classe, "phylum":phylum, "kingdom":kingdom,"superkingdom":superkingdom}
	    
    return dict_taxonomy


def getDictTaxo(dict_taxonomy) :
	# Make a dictionary of taxid (key) and its taxonomy separated by ';' (value)
	
	dict_taxidToTaxo = {}
	
	for taxid in dict_taxonomy:
		kingdom = dict_taxonomy[taxid]["kingdom"]
		superkingdom = dict_taxonomy[taxid]["superkingdom"]
		
		# Keep only genus and species names (no spp) and remove accessory formating.
		sp_list = dict_taxonomy[taxid]["species"].split()
		genus_species = "_".join(sp_list[0:2])
		genus_species = re.sub('\[|\]', '', genus_species)
		
		if kingdom and kingdom != "unknown kingdom" :
			taxonomy = sep.join([dict_taxonomy[taxid]["kingdom"],dict_taxonomy[taxid]["phylum"],dict_taxonomy[taxid]["class"],dict_taxonomy[taxid]["order"],dict_taxonomy[taxid]["family"],dict_taxonomy[taxid]["genus"],genus_species])
		elif superkingdom and superkingdom != "unknown superkingdom" :
			taxonomy = sep.join([dict_taxonomy[taxid]["superkingdom"],dict_taxonomy[taxid]["phylum"],dict_taxonomy[taxid]["class"],dict_taxonomy[taxid]["order"],dict_taxonomy[taxid]["family"],dict_taxonomy[taxid]["genus"],genus_species])
		else : # Taxonomy not found
			continue
		
		set_taxid = set()
		set_taxid.add(taxid)
		for idtax, taxo in dict_taxidToTaxo.items():
			if taxonomy == taxo :	# if several taxids share the same taxonomy, the key of this taxonomy will be a frozenset of all taxids.
				set_taxid = set_taxid.union(idtax)
				del dict_taxidToTaxo[idtax]
				break
		
		dict_taxidToTaxo[frozenset(set_taxid)] = taxonomy

	return dict_taxidToTaxo

	
							
def getDictTaxidToCount(sample_path, dict_taxidToTaxo, sample) :
	""" 
	  Make a dictionary of taxid or group of taxid (key) and its number of reads (value) in a sample.
	  If multiple taxa show a best Blast bitscores for a read, a set is generated from their taxids
	  e.g. dictTaxidToCount[frozenset({'1234680'}), frozenset({'889205', '257758'})] = 28.0
	"""
	
	dictTaxidToCount = {}
	lastLineCount = linecount(sample_path)
	derep = True
	
	with open(sample_path, "r") as s:
		
		for num, line in enumerate(s): 
			query = line.split("\t")[0]

			if num == 0 : # Start first query
				queryList = []
				queryLines = []
				queryList.append(query)
				queryLines.append(line.rstrip().split("\t"))
				try :
					nread = float(query.split("=")[1]) # get original read count before dereplication, which follows the query name with a "=" if using vsearch.
				except IndexError : # read not dereplicated
					derep = False
					nread = 1
				continue
							
			if query in queryList :
				queryLines.append(line.rstrip().split("\t"))
			
			else : # new query, store previous query results in dict.
	
				list_taxid = get_bestHitTaxid(queryLines)
				if list_taxid :
					set_taxid = set()
					for taxid in dict_taxidToTaxo.keys():#?????
						if taxid.intersection(set(list_taxid)): 
							set_taxid.add(taxid)
						
					try : # sum with nbread value
						dictTaxidToCount[frozenset(set_taxid)] += nread
					except KeyError : # new dict entry
						dictTaxidToCount[frozenset(set_taxid)] = nread

				queryList = [] # start new query
				queryList.append(query)
				queryLines = []
				queryLines.append(line.rstrip().split("\t"))
				if derep : nread = float(query.split("=")[1])
								
			if num == lastLineCount : # store last query
				list_taxid = get_bestHitTaxid(queryLines)
				if list_taxid :
					set_taxid = set()
					for taxid in dict_taxidToTaxo.keys():#?????
						if taxid.intersection(set(list_taxid)): 
							set_taxid.add(taxid)
					try : # sum with nbread value
						dictTaxidToCount[frozenset(set_taxid)] += nread
					except KeyError : # new dict entry
						dictTaxidToCount[frozenset(set_taxid)] = nread
				
	return (dictTaxidToCount, sample)


def getDictTaxidToCountTotal(list_dictTaxidToCount) :
	# Make a dictionary of taxid (key) and their total number of reads (value) found in all samples
	
	dict_taxidToCountTotal = {}
	
	for (dictTaxidToCount, sample) in list_dictTaxidToCount:
		for (taxid, count) in dictTaxidToCount.items() :
			try :
				dict_taxidToCountTotal[taxid] += count
			except KeyError : 
				dict_taxidToCountTotal[taxid] = count
				
	return(dict_taxidToCountTotal)
		
				
def getSetMinimalGroups(unique_set_taxid):
	"""
	  Step 1) Search all overlaps between taxids groups, i.e. when groups share at least 1 spp. 
	  Step 2) Create minimal groups (mg), i.e. groups defined by the smallest intersection(s) of an overlap with all other overlaps. If no intersection is found, the overlap is its own minimal group.
	"""
	
	# Step 1.
	set_overlap = set() 
	
	for set_taxid1 in unique_set_taxid :
		for set_taxid2 in unique_set_taxid : 
			overlap = set_taxid1.intersection(set_taxid2)
			if (overlap) : 
				set_overlap.add(overlap)

	# Step 2.
	set_all_mg = set() 
	
	for overlap1 in set_overlap :
		set_mg = set()
		if (len(overlap1) == 1) :
			set_all_mg.add(overlap1)
			continue
	
		for overlap2 in set_overlap :
			current_mg = overlap1.intersection(overlap2)
			if (current_mg) :
				if (set_mg) : # if new mg has intersection with previously found mg, replace by newest intersection to form a smaller mg.
					for mg in set_mg.copy() :
						if (current_mg.intersection(mg)):
							current_mg = current_mg.intersection(mg)
							set_mg.remove(mg)
							break
					set_mg.add(current_mg)
				else :
					set_mg.add(current_mg)
		
		for mg in set_mg : 
			set_all_mg.add(mg)
	
	return(set_all_mg)


def getDictMgToTaxo(set_all_mg, dict_taxidToTaxo):
	"""
	  Rebuild the taxonomy of the minimal groups in a dictionnary of set_mg (key) and their taxonomies (values)
	  e.g. dict_taxidToTaxo[frozenset({frozenset({'1121326, 1121327'}), frozenset({'1650662'}), frozenset({'382954'}), frozenset({'536227'}), frozenset({'339861'}), frozenset({'332101'})})] = 'Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium_aciditolerans-carboxidivorans-drakei-magnum-nitrophenolicum__Desnuesiella_massiliensis'
	"""
	
	dict_MgToTaxo = {}
	
	for mg in set_all_mg :
		list_mg_taxo = []
		
		for set_taxid in mg :
			list_mg_taxo.append(dict_taxidToTaxo[set_taxid])
			
		list_mg_taxo.sort()
		taxo = sep.join(list_mg_taxo[0].split(sep)[:-1])
		
		set_genera = set([taxo.split(sep)[-1].split("_")[0] for taxo in list_mg_taxo])
		
		list_genera_spp = []
		for genera in set_genera :
			list_spp = [mg_taxo.split(sep)[-1].split("_")[1] for mg_taxo in list_mg_taxo if genera in mg_taxo.split(";")[-1]]
			list_spp.sort()
			list_genera_spp.append(genera + "_" + "-".join(list_spp))
		
		mg_taxo = taxo + sep + "__".join(sorted(list_genera_spp))
		dict_MgToTaxo[mg] = mg_taxo

	return(dict_MgToTaxo)




def getDictMgToWeight(dict_taxidToCountTotal, set_all_mg): 
	""" 
	Make a dictionary of minimal groups (key) and their weight (value) accross all samples
	The weight is defined by the number of reads belonging to a mg divided by the total number of reads
	"""
	
	dict_MgToCountTotal = {}
	dict_MgToWeight = {}


	for set_mg in set_all_mg :
		for (taxid, count) in dict_taxidToCountTotal.items() :
			if set_mg.intersection(taxid):
				try :
					dict_MgToCountTotal[set_mg] += count
				except KeyError : 
					dict_MgToCountTotal[set_mg] = count
	
	total_count = sum(dict_MgToCountTotal.values())

	for (set_mg, total_mg_count) in dict_MgToCountTotal.items() :
		dict_MgToWeight[set_mg] = total_mg_count/total_count
		
	return(dict_MgToWeight)
			
			

def getDictMgToCount(dictTaxidToCount, sample, dict_MgToTaxo, dict_MgToWeight) :
	#  Make a dictionary of minimal groups' taxonomy (key) and their number of reads (value) found in a sample.
	
	dictMgToCount = {}
	
	for set_taxid, count in dictTaxidToCount.items() : 
		
		list_mg = []
		sum_weight = 0
		
		for mg, mg_weight in dict_MgToWeight.items() :
			if (set_taxid.intersection(mg)) : 
				list_mg.append(mg)
				sum_weight += mg_weight
		
		for mg in list_mg :
			mg_taxo = dict_MgToTaxo[mg]
			mg_weight = dict_MgToWeight[mg]
			
			try :
				dictMgToCount[mg_taxo] += float((count*mg_weight)/sum_weight)
			except KeyError : 
				dictMgToCount[mg_taxo] = float((count*mg_weight)/sum_weight)


	return(dictMgToCount, sample)
	
		
def makeOutputTables(list_dictMgToCount, dict_MgToTaxo, list_sampleNbReads) :
	"""
	  Generate output tables of read counts, where each column is a sample and each row is a minimal group taxonomy. In last column the full taxonomy in written based on taxid and separated by ';'
	  16S copy number normalization is performed if option & database are informed in input.
	"""

	list_samples = [sample for (dict_MgToCount, sample) in list_dictMgToCount]
	list_ordered_nbreads = [nbreads for sampleA in list_samples for (nbreads, sampleB) in list_sampleNbReads if sampleA == sampleB] # reorder samples in case multithreading changed it
	
	if sep == "\t" :
		first_line = "\t".join(list_samples) + "\t" + "\t".join(["Kingdom","Phylum","Class","Order","Family","Genus","Species"]) + "\n"
	else : 
		first_line = "\t".join(list_samples) + "\t" + "taxonomy" + "\n"
	
	with open(count_file, "w") as c:
		c.write(first_line)
		for set_mg, mg_taxo in dict_MgToTaxo.items() : 
			
			if abund_threshold > 0 :
				write_bool = False # output writing flag if abundance meets input threshold. e.g. write if MG represents at least 1% of reads in at least 1 sample.
			else :  
				write_bool = True 
				
			output_count_list = []
			
			for sample in range(0, len(list_samples)):
				sample_count_mg = 0
				sample_abund_mg = 0
				
				try : 
					sample_count_mg = list_dictMgToCount[sample][0][mg_taxo]
					sample_abund_mg = round((sample_count_mg/list_ordered_nbreads[sample]*100),2)
				except KeyError : # MG not in sample
					pass

				output_count_list.append(int(sample_count_mg))
				
				if list_ordered_nbreads[sample] != 0 and not write_bool: # check if sample is not empty or has no reads passing alignment thresholds
					if (sample_abund_mg >= abund_threshold) : # if abundance threshold requirement is met by at least one sample, taxon OK to write in output.
						write_bool = True
				
			if write_bool :
				c.write("\t".join(map(str, output_count_list)) + "\t" + mg_taxo + "\n")
	
	
	if norm_bool : 
		dictGeneraCopynb, mean_norm = getDictGeneraCopynb(db16S) 
		
		# Creating new lists with updated samples' & taxa's nbreads after count normalization.
		list_normnbread_allsamples = []
		list_normcount_allsamples = []

		with open(norm_countfile, "w") as cn :
			cn.write(first_line)
			list_taxa = []
			
			for set_mg, mg_taxo in dict_MgToTaxo.items() :
				list_taxa.append(mg_taxo)
				output_normcount_list = []
				genus = mg_taxo.split(sep)[-2]

				try :
					copynb = dictGeneraCopynb[genus]
				except KeyError : # unclassified or entry not in 16Sdb
					copynb = mean_norm
				
				for sample in range(0, len(list_samples)):
					sample_normcount_mg = 0
					try : 
						sample_normcount_mg = (list_dictMgToCount[sample][0][mg_taxo]/copynb)
					except KeyError : # MG not in sample
						pass
					output_normcount_list.append(int(sample_normcount_mg))
					
					# update total number of reads in each sample after accounting for normalization
					if (len(list_normnbread_allsamples) > sample) :
						list_normnbread_allsamples[sample] += sample_normcount_mg
					else : 
						list_normnbread_allsamples.append(sample_normcount_mg)
				list_normcount_allsamples.append(output_normcount_list)
				
				
			for taxon in range(0, len(list_taxa)):
				if ( sum(list_normcount_allsamples[taxon]) >= 1 ) : # Some taxa now have 0 count after 16S copy number normalization, keep only if > 0.
					write_bool = False
					for sample in range(0, len(list_samples)) :
						if (list_ordered_nbreads[sample] != 0): # check if sample is not empty or have no reads passing alignment thresholds
							sample_normabund_taxon = round((list_normcount_allsamples[taxon][sample]/list_normnbread_allsamples[sample]*100),2) 
							if (sample_normabund_taxon >= abund_threshold) : # if abundance threshold requirement is met by at least one sample, taxon OK to write in output.
								write_bool = True
					if write_bool :
						cn.write("\t".join(map(str, list_normcount_allsamples[taxon])) + "\t" + list_taxa[taxon] + "\n")
				


def makeMatrices():
	#### Generating dictionaries of samples, read counts, taxonomies and output tables
	
	dict_samples = getDictFiles(input_file)
	list_taxidAllSamples = Parallel(n_jobs=nb_threads)(delayed(getListTaxid)(dict_samples[sample]) for sample in dict_samples.keys())
	list_unique_taxid = list(set([taxid for list_taxidsample in list_taxidAllSamples for taxid in list_taxidsample]))
	dict_taxonomy = getDictReferenceTaxonomy(list_unique_taxid)
	dict_taxidToTaxo = getDictTaxo(dict_taxonomy)
	list_dictTaxidToCount = Parallel(n_jobs=nb_threads)(delayed(getDictTaxidToCount)(dict_samples[sample], dict_taxidToTaxo, sample) for sample in sorted(dict_samples.keys()))
	dict_taxidToCountTotal = getDictTaxidToCountTotal(list_dictTaxidToCount)
	unique_set_taxid = set([set_taxid for (dictTaxidToCount, sample) in list_dictTaxidToCount for set_taxid in dictTaxidToCount])
	list_sampleNbReads = Parallel(n_jobs=nb_threads)(delayed(getSampleNbReads)(dictTaxidToCount, sample) for (dictTaxidToCount, sample) in list_dictTaxidToCount)
	totalNbReads = sum(nbReads for (nbReads, sample) in list_sampleNbReads)
	set_all_mg = getSetMinimalGroups(unique_set_taxid)
	dict_MgToTaxo = getDictMgToTaxo(set_all_mg, dict_taxidToTaxo)
	dict_MgToWeight = getDictMgToWeight(dict_taxidToCountTotal, set_all_mg)
	list_dictMgToCount = Parallel(n_jobs=nb_threads)(delayed(getDictMgToCount)(dictTaxidToCount, sample, dict_MgToTaxo, dict_MgToWeight) for (dictTaxidToCount, sample) in list_dictTaxidToCount)
	makeOutputTables(list_dictMgToCount, dict_MgToTaxo, list_sampleNbReads) 
	
	
def main():
	makeMatrices()
	
	
if __name__ == "__main__":
	main()
