# Understanding how human language evolved with ES-PGS

## Overview
Evolutionary Stratified Polygenic Score (ES-PGS) analysis is a novel method to partition polygenic scores based on evolutionary annotations, enabling the investigation of how genetic variants from different evolutionary periods contribute to complex traits. We used ES-PGS to understand which evolutionary events contributed to human language, code for these analyses can be found in the `manuscript` directory. 

To get the most out of ES-PGS you will need a decent understanding of polygenic scores, [this review](https://pmc.ncbi.nlm.nih.gov/articles/PMC7612115/) should cover a lot of the necessary background information. Papers describing "pathway based PGS" may also be helpful, [like this one](https://pubmed.ncbi.nlm.nih.gov/36749789/).

Below is a graphical abstract showcasing how ES-PGS can be used to better understand evolutionary events influencing complex traits (in our case, human language).
![graphical abstract](manuscript/figures/graphical_abstract.png)

## Requirements
In this repository we adapt established tools for our ES-PGS calculation, which should minimize needing to download and compile new software. Basic requirements can be seen below:

- `BEDTools` (for partitioning the genome), [download here](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- `PRSet` (for ES-PGS calculation), [download here](https://choishingwan.github.io/PRSice/)
- `R` (for ES-PGS modeling) and the `tidyverse` package (for data manipulation and visualization), [download here](https://cran.r-project.org/)

## Example usage
We provide an example implementation and output of ES-PGS in the `ES-PGS_example` directory, we hope this code can be easily adapted to a wide range of research areas.

A general ES-PGS analysis has 6 steps, which are covered in our example:
1. Gather evolutionary annotations and identify "background" regions with `BEDTools`
2. Generate matched control regions for your evolutionary annotations using R
3. Compute ES-PGS with `PRSet`
4. Adjust ES-PGS for genetic PCs and reformat data in `R`
5. ES-PGS analysis with your phenotype(s) in `R`
6. Visualize results with `R`

## Citations
If you use this method, please cite our manuscript describing [ES-PGS](https://www.biorxiv.org/content/10.1101/2025.03.07.641231v1).

If you use the recommended software to compute and run the ES-PGS analysis, cite the following papers:
- [BEDTools](https://pubmed.ncbi.nlm.nih.gov/20110278/)
- [PRSet](https://pubmed.ncbi.nlm.nih.gov/36749789/)
- `R` (in R run `citation()` and `citation("tidyverse")` to get the appropriate info for your version)