# RNAseq
Expression analysis following induced expression in *Bacillus subtilis* in exponential growth

The raw sequencing data was deposited at NCBI's Gene Expression Omnibus under accession number [GSE187004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE187004).

Sequencing data was processed by the Center for Genomics and Bioinformatics at Indiana University, as described in the manuscript. The results of that analysis are differential expression (DE) values which are the strating point for the analysis here. DE data can be found in [RNAseq/data/DEseq](data/DEseq).

This analysis had two, independent parts


* **volcano_plot.R**
Per-treatment plots of differential expression, with enrichment test for sporulation genes. 

* **corr_plot.R**
Correlation between gene expression when phage genes are induced and host genes are induced.
