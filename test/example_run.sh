# Step 1 : data files extraction ; optional if already performed
cd "${0%/*}"
cd ../data
7z e macotra.7z -y
# Step 2 : create manifest file ; contains sample names and their corresponding Blast file path
cd "${0%/*}"
rm manifest.txt
for i in `ls ../data/*.blastn`; do j=$(basename -- "$i") ; k="${j%.*}"; echo $k:$i >> manifest.txt; done
# Step 3 : run taxon resolver algorithm
python3 ../bin/blast_taxon_resolver.py -i manifest.txt -dr ../database/rankedlineage.dmp -n yes -dn ../database/rrnDB-5.5_pantaxa_stats_NCBI_genus.tsv -o macotra_table.tsv 
