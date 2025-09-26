# gene_name_converter
Converts gene names of multiple formats into one universal format. In this case, the tool accepts: Ensembl IDs, NCBI gene IDs, HGNC IDs, Approved symbols, Previous symbols, and Alias symbols. In this iteration, the tool outputs the approved symbol, approved name, alias symbols, and previous symbols.

An example use case is shown in gene_lookup_example.ipynb

# Version 2:
Primarily addresses issues of genes mapping to multiple records. This was combatted via the following:
1. For genes that could be identified as ENSEMBL/NCBI/HGNC ID, only their respective column was searched for increased speed and accuracy
2. For genes that still returned multiple records, all instances are logged with corresponding match types.

# Version 3:
1. Added handling of different delimiters in search_single_gene, allowing ',' and ', '
2. Added a column for match status and original user input
3. Modified logging to account for unmatched genes
### Changes to functions
- Modified addColumns to be faster by only searching until the categories are found
- Modified createDownloadURL to be faster by avoiding repeated string concat
- Rewrote transform_string to be cleaner and slightly faster by using partition for splitting
- Modified search_single_gene to be faster by only opening the file once and uses sets for faster lookup
- getData and convert_gene_names rewritten for less indentation

# Version 4:
1. Modified handling of NCBI/Entrez to just use .isdigit() instead of ReGex
2. Got rid of code block in findAPI function that was unsafe