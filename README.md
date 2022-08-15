# BreastCancerR

A script based on Seurat Package (https://satijalab.org/seurat/) to interrogate publicly available scRNAseq data on breast cancer. The objective is to interrogate the expression of SPARC gene in the cellular subsets found in the breast tumor (as well as test the significance of differential expression in the individual cell populations). The script can be simply adapted to interrogate the expression of any other gene present in the dataset.

Following datasets should be loaded to re-create the figures of the article (Alcaraz L. et al. pending):

Figure 3: article: PMID: 32790115; https://www.embopress.org/doi/full/10.15252/embj.2019104063 [data: https://ega-archive.org/datasets/EGAD00001006981]


Figure S6: article: PMID: 30181541; https://www.nature.com/articles/s41467-018-06052-0 [data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118390]


Figure S7: article: PMID: 32434947; https://aacrjournals.org/cancerdiscovery/article/10/9/1330/2752/Single-Cell-Analysis-Reveals-Fibroblast-Clusters [data:https://ega-archive.org/search-results.php?query=EGAS00001004030]


N.B.: Installation of specific R Packages (incl. Seurat) is required before using the script; See the "Libraries" section of the script. Questions/issues regarding the availability of the abovementioned data, please address to the corresponding author of the respective study.


