# Applying ES-PGS

## Example pipeline
We provide an example implementation and output of ES-PGS in the `code` directory.

This example has 6 steps:
1. Gather evolutionary annotations and identify "background" regions with `BEDTools`
2. Generate matched control regions for your evolutionary annotations using `R`
3. Compute ES-PGS with `PRSet`
4. Adjust ES-PGS for genetic PCs and reformat data in `R`
5. ES-PGS analysis in `R`
6. Visualize results with `R`

The primary pipeline script (`code/ES-PGS_pipeline.sh`) and accompanying R scripts (`code/Rscripts/`) should hopefully be easily adapted to your own data.

You will not be able to run this example without access to the EpiSLI / 1000 Genomes merged sample described in our manuscript. To run the example you would need to download the EpiSLI data from dbGaP and 1000 Genomes Phase 3 data and merge them. You should be able to easily swap out the file paths with your own data in these scripts to run an ES-PGS analysis.

Additionally, you will likely want to use imputed SNP genotypes or high-coverage whole genome sequencing data. Many evolutionary annotations make up a very small proportion of the genome, relying only on directly genotyped SNPs from SNP array data may not provide enough coverage for ES-PGS analysis.

## Citations
If you use this method, please cite our manuscript describing [ES-PGS](https://www.biorxiv.org/content/10.1101/2025.03.07.641231v2).

If you use the recommended software to compute and run the ES-PGS analysis, cite the following papers:
- [BEDTools](https://pubmed.ncbi.nlm.nih.gov/20110278/)
- [PRSet](https://pubmed.ncbi.nlm.nih.gov/36749789/)
- `R` (in R run `citation()` and `citation("tidyverse")` to get the appropriate info for your version)