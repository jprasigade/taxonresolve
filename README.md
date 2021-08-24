# Blast taxon resolver
## Introduction
### 1. Purpose
Classifiers aim to perform taxonomic analyses in order to uncover microbial diversity within samples, and are often a necessary step in metagenomics and metabarcoding projects. Many tools have been developed to accomodate increasingly larger projects, but species resolution trade-offs are often needed in order to perform faster and large-scale analyses. To this day, the gold standard remains BLAST (Altschul et al. 1990) as the most accurate alignment algorithm in most situations. For a read query, BLAST outputs several alignment scores against a database of reference sequences, offering an insight into which microbe the read may originate from, at the cost of long analyses run time.

Nonetheless, several technical limitations may hinders those taxonomic assignations, inducing BLAST to wrongly or ambiguously output a taxon as the best candidate(s). These issues can spawn from PCR and sequencing errors, limited information from genomic content (e.g. short or scarce reads), genetic proximity between several species (especially with metabarcoding samples), non-exhaustive reference database, *etc*.

This project stems from the need to improve BLAST's resolution at the species level by taking into consideration the overall taxa distribution within samples. Its intent is to correct, if at least minimise, the limitations previously mentionned.

### 2. How
The aim of this program is to generate a table of read counts to be used for diversity & ecological analyses, based on corrected BLAST results. 

First, if a BLAST search returns several best hits with identical bitscores for a read and cannot resolve the most likely taxon to assign, a group of species is formed. 
The code then resolves the eventual taxonomical overlaps between groups, by generating minimal groups based on their intersections.
It then computes the weight of minimal groups within all samples, to distribute more sparingly the read counts when they are shared between several minimal groups.

The output displays the read counts of each species (rows) for each sample (column). Bacterial read numbers can be normalized to represent species counts instead, provided a 16S copy number database ; in such case 2 contingency tables will be produced, with and without 16S normalization. 

### 3. Prerequisite
The BLAST files must be generated respecting the following formatting options : 
   1) `-outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore'`
   2) `-max_target_seqs >= 1 (e.g. -max_target_seqs 20)`

This script allows the possibility to dereplicate reads prior to BLAST in order to accelerate computing times ; \
e.g. with the `vsearch --derep_fulllength --sizeout` command.

## Download
### All necessary scripts, databases and test datasets are provided in this projet.
A **run example** is provided as well in the test directory, along with the expected outputs.

To note, the code makes use of 2 databases, which can be occasionally updated as follow :
- A taxonomy reference, compiled in the new_taxdump files from the NCBI
>wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip & unzip new_taxdump.zip
- A rRNA 16S copy number database at the genus level, provided by the NCBI 
>wget https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.5_pantaxa_stats_NCBI.tsv.zip & unzip rrnDB-5.5_pantaxa_stats_NCBI.tsv.zip
>
>head -n 1 rrnDB-5.5_pantaxa_stats_NCBI.tsv > rrnDB-5.5_pantaxa_stats_NCBI_genus.tsv
>
>grep "genus" rrnDB-5.5_pantaxa_stats_NCBI.tsv | grep -v "species" >> rrnDB-5.5_pantaxa_stats_NCBI_genus.tsv

Some genera may be missing, a good practice is to check if the genus of your species of interest are in this table. Otherwise an average of the sample's 16S copy numbers will be used. You can manually add the missing entries by checking a species 16S copy numbers here : https://rrndb.umms.med.umich.edu/

If the scripts provided do not run, you may need to give execution rights, for example `chmod +x bin/blast_taxon_resolver.py`.

## Dependencies
This code runs in Python3 and requires NumPy's library.

## Run examples

### 1. Running Blast taxon resolver code

**Step 1**

The code requires a manifest file in order to locate and analyze multiple samples at once.\
&nbsp;&nbsp;&nbsp;&nbsp;e.g. manifest.txt content :
>sample1:/path/to/file1.blast\
>sample2:/path/to/file2.blast

Which can be manually done or automatically generated like in the following example :

>for i in \`ls ../data/\*.blastn\`; do j=$(basename -- "$i") ; k="${j%.*}"; echo $k:$i >> manifest.txt; done

**Step 2** 

A basic Blast taxon resolver analysis runs as follow :

`python3 ../bin/blast_taxon_resolver.py -i manifest.txt -dr ../database/rankedlineage.dmp`

As an advanced analysis example, let's perform rRNA 16S copy normalization, be more stringent on Blast's identity scores with a threshold of 97%, and keep only the taxonomies which are abundant at minimum 1% in at least one sample, using 40 threads :

`python3 ../bin/blast_taxon_resolver.py -i manifest.txt -dr ../database/rankedlineage.dmp -n yes -dn ../database/rrnDB-5.5_pantaxa_stats_NCBI_genus.tsv -a 1 -id 97 -t 40 -o macotra_table.tsv`

### 2. Testing the package

You can first try your hand on a test dataset using this command :

`test/example_run.sh`
 

## Test dataset
The dataset was obtained from the MACOTRA project work package 4, consisting of 16 patients. Each patient's nasal microbiota was sampled at 5 different time points to evaluate the effect of antibiotic therapy on bacterial diversity and success of *Staphylococcus aureus* decolonization.

Publication on the dataset and its analysis is underway.
## Contributors
* Anaïs Barray anais.barray@gmail.com (Initiator, coding, testing, documentation, evaluation)
* Jean-Philippe Rasigade jean-philippe.rasigade@univ-lyon1.fr (Feature suggestions, evaluation)
## Licence
The Blast taxon resolver code is licensed  under the GNU Affero General Public License v3.0. Please see LICENSE.txt for details.
## Bugs and Getting help
All bug reports are highly appreciated. You may send an email to anais.barray@gmail.com.

Similarly, any help can be inquired by email if the `blast_taxon_resolver.py --help` command doesn't provide the necessary information.
