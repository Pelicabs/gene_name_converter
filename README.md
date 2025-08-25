# gene_name_converter
Converts gene names of multiple formats into one universal format. In this case, the tool accepts: Ensembl IDs, NCBI gene IDs, HGNC IDs, Approved symbols, Previous symbols, and Alias symbols. In this iteration, the tool outputs the approved symbol, approved name, alias symbols, and previous symbols.

An example use case is shown in gene_lookup_example.ipynb

# Version 2:
Primarily addresses issues of genes mapping to multiple records. This was combatted via the following:
1. For genes that could be identified as ENSEMBL/NCBI/HGNC ID, only their respective column was searched for increased speed and accuracy
2. For genes that still returned multiple records, all instances are logged with corresponding match types.
