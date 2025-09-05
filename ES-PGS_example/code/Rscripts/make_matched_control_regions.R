## -----------------------------------------
## load all dependencies (there are a lot)
## -----------------------------------------
## required
library(tidyverse)
library(regioneR)  # This has randomizeRegions
library(nullranges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(GenomicFeatures)  # Add this for cds(), exons(), genes(), promoters()
## recommended: for parallelization of random sampling 
library(future)
library(furrr)
# Set up parallel processing
plan(multisession, workers = 6)  # Adjust based on your CPU cores

# chromosomes of interest (tweak to match your data - sex chromosomes will probably be dropped in PGS calculations)
target_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
hg19_granges <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg19))
hg19_filtered <- hg19_granges[seqnames(hg19_granges) %in% target_chroms]
hg19_filtered <- keepSeqlevels(hg19_filtered, target_chroms, pruning.mode = "coarse")

## get annotation file list to loop over
files <- c("ES-PGS_example/example_data/HAQERs.bed", "ES-PGS_example/example_data/HARs.bed")

## for each annotation, find matched random regions
for(f in files) {
    cat('\n\n\n\n\n\n')
    message(rep('<>', times = 15))
    message('Working on: ', f)
    message(rep('<>', times = 15))
    message("Reading / prepping data...")
    # Load your regions and make sure we're only working with chromosomes of interest
    original_regions <- import(f)
    original_regions <- original_regions[seqnames(original_regions) %in% target_chroms]

    # Create buffer zones around regions of interest
    buffer_distance <- 100000  # 100kb buffer around ROIs
    original_with_buffer <- resize(original_regions, 
                                width = width(original_regions) + 2 * buffer_distance, 
                                fix = "center")
    excluded_fraction <- sum(width(reduce(original_with_buffer))) / sum(seqlengths(hg19_filtered)[target_chroms])
    message("Excluded fraction of genome (after 100Kb mask around ROIs): ", round(excluded_fraction * 100, 4), "%\n")
    # Generate many random regions for the sequence matching pool of the same size on the same chromosomes of our ROIs
    message("Generating random chromosome/size matched sequences so we can make a matched control set, might take a few minutes...")

    # Parallelized random sequence generation
    ## this will sample 1,000 random size matched sequences for each input region (later we will identify the closest matching one)
    random_sets <- future_map(1:1000, function(i) {
                                            set.seed(i + 42)
                                            randomizeRegions(
                                                            A = original_regions,
                                                            genome = hg19_filtered,
                                                            allow.overlaps = FALSE,
                                                            per.chromosome = TRUE,
                                                            mask = ## mask regions of interest + buffer zone
                                                            )  
                                                            }, .options = furrr_options(seed = TRUE))

    # Combine and deduplicate
    random_pool <- do.call(c, random_sets)
    rm(random_sets)
    gc()
    pool_regions <- random_pool[!duplicated(random_pool)]
    # Make sure we exclude any that overlap with your original regions (though, this should be taken care of in randomizeRegions())
    pool_regions <- random_pool[!overlapsAny(random_pool, original_regions)]

    # Add required annotations
    message("Annotating genomic regions so we can accurately match them...")
    original_regions <- original_regions %>%
      mutate(
          # GC content
          gc = getSeq(BSgenome.Hsapiens.UCSC.hg19, .) %>% 
              letterFrequency(letters = "GC") %>% 
              rowSums() / width(.),
          
          # CDS density (per kb)
          cds_density = countOverlaps(., cds(txdb)) / (width(.) / 1000),
          
          # Exon density (per kb)
          exon_density = countOverlaps(., exons(txdb)) / (width(.) / 1000),
          
          # Promoter density (TSS +/- 2kb, per kb)
          promoter_density = countOverlaps(., promoters(txdb, upstream = 2000, downstream = 2000)) / (width(.) / 1000),
          
          # Repeat content (N content as proxy)
          repeat_content = getSeq(BSgenome.Hsapiens.UCSC.hg19, .) %>%
                          letterFrequency(letters = "N") %>%
                          rowSums() / width(.),
          
          # Distance to nearest gene
          nearest_gene_dist = distanceToNearest(., genes(txdb)) %>% mcols() %>% .$distance,
          
          # Number of genes within 1 megabase
          gene_density_1mb = countOverlaps(resize(., width=1e6, fix="center"), genes(txdb)),

          # Make binary indicators and log gene distance
          cds_overlap = ifelse(cds_density > 0, 1, 0),
          exon_overlap = ifelse(exon_density > 0, 1, 0),
          promoter_overlap = ifelse(promoter_density > 0, 1, 0),
          nearest_gene_dist_log10 = log10(nearest_gene_dist + 1)
        )

    pool_regions <- pool_regions %>%
      mutate(
          # GC content
          gc = getSeq(BSgenome.Hsapiens.UCSC.hg19, .) %>% 
              letterFrequency(letters = "GC") %>% 
              rowSums() / width(.),
          
          # CDS density (per kb)
          cds_density = countOverlaps(., cds(txdb)) / (width(.) / 1000),
          
          # Exon density (per kb)
          exon_density = countOverlaps(., exons(txdb)) / (width(.) / 1000),
          
          # Promoter density (TSS +/- 2kb, per kb)
          promoter_density = countOverlaps(., promoters(txdb, upstream = 2000, downstream = 2000)) / (width(.) / 1000),
          
          # Repeat content (N content as proxy)
          repeat_content = getSeq(BSgenome.Hsapiens.UCSC.hg19, .) %>%
                          letterFrequency(letters = "N") %>%
                          rowSums() / width(.),
          
          # Distance to nearest gene
          nearest_gene_dist = distanceToNearest(., genes(txdb)) %>% mcols() %>% .$distance,

          # Number of genes within 1 megabase
          gene_density_1mb = countOverlaps(resize(., width=1e6, fix="center"), genes(txdb)),

          # Make binary indicators and log gene distance
          cds_overlap = ifelse(cds_density > 0, 1, 0),
          exon_overlap = ifelse(exon_density > 0, 1, 0),
          promoter_overlap = ifelse(promoter_density > 0, 1, 0),
          nearest_gene_dist_log10 = log10(nearest_gene_dist + 1)
        )
    # Add 'width' as a metadata column
    mcols(original_regions)$range_width <- width(original_regions)
    mcols(pool_regions)$range_width <- width(pool_regions)

    # Clean data 
    # Function to completely clean a GRanges object
    clean_granges_completely <- function(gr) {
    # Extract all components
    chr <- as.character(seqnames(gr))
    start_pos <- start(gr)
    end_pos <- end(gr)
    strand_info <- as.character(strand(gr))
    
    # Extract metadata as clean data.frame
    if (ncol(mcols(gr)) > 0) {
        meta_data <- as.data.frame(mcols(gr))
        rownames(meta_data) <- NULL
        names(meta_data) <- make.names(names(meta_data))  # Ensure valid column names
    } else {
        meta_data <- data.frame()
    }
    
    # Create new GRanges from scratch
    new_gr <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = start_pos, end = end_pos),
        strand = strand_info
    )
    
    # Add metadata if it exists
    if (nrow(meta_data) > 0) {
        mcols(new_gr) <- meta_data
    }
    
    return(new_gr)
    }
    pool_regions <- clean_granges_completely(pool_regions)

    # Remove exact duplicate ranges (same chr, start, end) (REDUNDANT CODE)
    pool_regions <- pool_regions[!duplicated(pool_regions)]

    # Match using the covariates
    message("Matching genomic regions...")
    set.seed(1010)
    matched_regions <- matchRanges(
      focal = original_regions,
      pool = pool_regions,
      covar = ~ range_width + gc + repeat_content + cds_overlap + promoter_overlap + nearest_gene_dist_log10 + gene_density_1mb,
      method = "stratified", ## package default method ("rejection") fails, this one seems more robust
      replace = FALSE
    )

    ## plot matched sequences vs actual sequences on covariate info
    # plotCovariate(matched_regions, covar = "gc")
    # plotCovariate(matched_regions, covar = "cds_density")
    # plotCovariate(matched_regions, covar = "promoter_density")
    # plotCovariate(matched_regions, covar = "repeat_content")
    # plotCovariate(matched_regions, covar = "nearest_gene_dist")

    # Save matched regions
    message("Done with data processing for this annotation...")
    tmp <- as.data.frame(matched_regions) %>% 
        as_tibble() %>% 
        mutate(chr = as.numeric(str_remove_all(seqnames, pattern = 'chr'))) %>% 
        arrange(chr, start, end) 
    tmp$chr <- NULL

    cat("\n\n\n")
    message("# rows in original: ", nrow(as.data.frame(original_regions)))
    message("# rows in matched control set: ", nrow(tmp))
    message("# unique rows in matched control set: ", nrow(distinct(tmp[,1:3])))

    ##
    message("Saving files now...")
    outf <- str_c("ES-PGS_example/example_data/", str_replace_all(basename(f), '.bed', '.matched_control_sequences.bed'))
    tmp[,1:3] %>% 
        write_tsv(outf, col_names = FALSE)
    # Save genomic info about regions of interest
    outf2 <- str_c("ES-PGS_example/example_data/", str_replace_all(basename(f), '.bed', '_genomic_annotation_info.csv'))
    original_df <- as.data.frame(original_regions) %>% 
        as_tibble() %>% 
        mutate(chr = as.numeric(str_remove_all(seqnames, pattern = 'chr'))) %>% 
        arrange(chr, start, end) %>% 
        dplyr::select(-chr)
    original_df %>% 
        write_csv(outf2)

    # Save genomic info about matched control regions
    outf3 <- str_c("ES-PGS_example/example_data/matched_control.", str_replace_all(basename(f), '.bed', '_genomic_annotation_info.csv'))
    matched_df <- as.data.frame(matched_regions) %>% 
        as_tibble() %>% 
        mutate(chr = as.numeric(str_remove_all(seqnames, pattern = 'chr'))) %>% 
        arrange(chr, start, end) %>% 
        dplyr::select(-chr)
    matched_df %>% 
        write_csv(outf3)

    # Clean up workspace
    rm(original_regions, matched_df, original_df, random_pool, pool_regions, tmp, outf, outf2, outf3)
    gc()
}
## done