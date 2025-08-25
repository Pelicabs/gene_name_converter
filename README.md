# gene_name_converter
Converts gene names of multiple formats into one universal format

# Version 2:
Primarily addresses issues of genes mapping to multiple records. This was combatted via the following:
1. For genes that could be identified as ENSEMBL/NCBI/HGNC ID, only their respective column was searched for increased speed and accuracy
2. For genes that still returned multiple records, all instances are logged with corresponding match types.