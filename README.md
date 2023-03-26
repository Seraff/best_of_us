# Best Of Us
The script is designed to pick out sugnificant genes for differntial expression analysis.
It uses DEseq2 and edgeR to determine the most signifficant genes.

## Input

### Counts Table

The standard tab-delimided table generated by featureCounts tool.

First line is a comment.

The structure of the table should be like this:

| Geneid  | Chr | Start | End | Length | /path/to/mapping_1.bam | /path/to/mapping_2.bam
| ------- | --- | ----- | --- | ------ | ---------------------- | ----------------------
| ------- | --- | ----- | --- | ------ | ---------------------- | ----------------------

Please, prepare your table to this format.

The most important columns which should be there are `Geneid` and "bam" columns with coverage.
If you don't have `Chr`, `Start`, etc additional columns, add empty ones instead.

### Metadata

The table should be generated manmually.

The structure of the table should be in the follwing format:

|  | sample_1  | sample_2   | sample_3 | sample_4 | ...
| - | -  | - | - | - | -
| state | anaerobic | anaerobic  | aerobic  | aerobic | ...

The name of row `state` is important. If you metadata table has a different name, rename it.

## Output

In the outpu folder there will be the following fies:

* `significant_genes.txt` - list of significant genes
* `vienn.png` - Vienn diagram of genes determined as significant for different analysis
* `heatmap.png` - plot of differential expression of significant genes
* `pca.png` - PCA plot of significant genes
* `ma.png` - MA plot of significant genes
* `volcano.png` - Volcano plot of significant genes
