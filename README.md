# gtex_enrichment

This is a project for testing expression enrichment in tissues from the whole body (tissue samples from GTEx), starting with a query of a gene set of interest.

## Data
Raw data was downloaded from the [GTEx portal](https://gtexportal.org/home/datasets).

- Sample annotations metadata -- GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
- RNASeq data -- GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

Expression data was extracted and processed with [cmapPy](https://github.com/cmap/cmapPy) in python2.7

Downstream processing was performed using python3 and pandas, numpy and scipy.

