## load packages
import sys
from Bio import AlignIO, Align
from Bio.Align import MultipleSeqAlignment
import gzip
from collections import defaultdict
import csv

## get chromosome to work on
i = sys.argv[1]
print("Chromosome number passed to script: " + str(i))

## function to read in and parse a BED file
def read_bed_file(bed_file):
    """
    Reads a BED file and extracts regions as a list of tuples.
    Args:
        bed_file (str): Path to the BED file.
    Returns:
        list: List of regions as (chrom, start, end) tuples.
    """
    regions = []
    with open(bed_file, 'r') as file:
        for line in file:
            # Skip empty lines or lines starting with a comment
            if line.strip() and not line.startswith("#"):
                fields = line.strip().split("\t")
                ## if no "chr" prefix uncomment the lines below
                # chr_pre = fields[0]
                # chrom = chr_pre.replace('chr', '')
                chrom = fields[0]
                start = int(fields[1])  # Start position (0-based)
                end = int(fields[2])    # End position (exclusive)
                regions.append((chrom, start, end))
    return regions


#########################################
## function to trim alignments just to sequences of interest
def subset_maf_by_species(maf_file, regions, output_file, species_id):
    """
    Subset a MAF file to include alignments overlapping specific regions based on one species (e.g., hg38).
    Alignments are trimmed based on the reference genome indices, keeping aligned bases for all species. 
    Args:
        maf_file (str): Path to the input MAF file (gzip compressed).
        regions (list): List of regions as tuples (chrom, start, end).
        output_file (str): Path to the output MAF file.
        species_id (str): Species identifier in the MAF file (e.g., "hg38").
    """
    with gzip.open(maf_file, 'rt') as infile, open(output_file, 'w') as outfile:
        for alignment in AlignIO.parse(infile, "maf"):
            # Find the reference species record
            ref_record = None
            for record in alignment:
                if record.id.startswith(species_id):
                    ref_record = record
                    break
            if not ref_record:
                continue  # Skip if no reference species in this alignment
            ref_start = ref_record.annotations["start"]
            ref_size = ref_record.annotations["size"]
            ref_end = ref_start + ref_size
            ref_chrom = ref_record.id.split('.')[1]
            for region in regions:
                chrom, region_start, region_end = region
                # Ensure chromosome matches and check for overlap
                if ref_chrom == chrom and max(region_start, ref_start) < min(region_end, ref_end):
                    # Calculate the overlap range relative to the reference species
                    overlap_start = max(region_start, ref_start)
                    overlap_end = min(region_end, ref_end)
                    # Find alignment column indices corresponding to this range
                    ref_alignment_indices = []
                    ref_pos = ref_start
                    for col_idx, base in enumerate(ref_record.seq):
                        if base != '-':  # Count only non-gap positions for the reference
                            if ref_pos >= overlap_start and ref_pos < overlap_end:
                                ref_alignment_indices.append(col_idx)
                            ref_pos += 1
                    if not ref_alignment_indices:
                        continue  # No valid overlap in this alignment
                    # Create a trimmed alignment for all species based on these indices
                    trimmed_alignment = Align.MultipleSeqAlignment([])
                    for record in alignment:
                        # Extract the corresponding aligned bases
                        trimmed_seq = "".join(record.seq[idx] for idx in ref_alignment_indices)  
                        # Create a new trimmed record
                        trimmed_record = record[:]
                        trimmed_record.seq = trimmed_seq
                        trimmed_record.annotations["start"] = overlap_start
                        trimmed_record.annotations["size"] = len(trimmed_seq)
                        trimmed_alignment.append(trimmed_record)
                    # Write the trimmed alignment if valid
                    if len(trimmed_alignment) > 0:
                        AlignIO.write(trimmed_alignment, outfile, "maf")


######################################
## function to compute sequence similarity with hg38 as the reference
def calculate_similarity_per_species(alignment):
    # Initialize dictionaries to store results for each species
    species_results = defaultdict(lambda: {"conserved_bases": 0, "sequence_length": 0, "similarity": 0})
    # Extract species names from sequence IDs
    species_list = [rec.id.split(".")[0] for rec in alignment]
    # Iterate through columns of the alignment to calculate conserved bases and similarity
    for i in range(len(alignment[0])):
        column = alignment[:, i]
        # Find the reference base (first species in the column)
        ref_base = column[0]
        for j, base in enumerate(column):
            species = species_list[j]
            if base.lower() == ref_base.lower():
                species_results[species]["conserved_bases"] += 1
            species_results[species]["sequence_length"] += 1
    # Calculate similarity for each species
    for species in species_results:
        conserved = species_results[species]["conserved_bases"]
        length = species_results[species]["sequence_length"]
        species_results[species]["similarity"] = conserved / length if length > 0 else 0
    # Include region metadata
    region_metadata = [
        f"{rec.id}:{rec.annotations['start']}-{rec.annotations['start'] + rec.annotations['size']}"
        for rec in alignment
    ]
    return {
        "species_results": species_results,
        "regions": region_metadata,
    }

## function to reformat and save results to a CSV file
def save_results_to_csv(maf_file, output_file):
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # Write the header
        writer.writerow(["region", "species", "conserved_bases", "seq_length", "similarity_pct"])
        # Parse the MAF file and calculate similarity for each alignment
        with open(maf_file, 'rt') as infile: ## gzip.
            for alignment in AlignIO.parse(infile, "maf"):
                results = calculate_similarity_per_species(alignment)
                species_results = results["species_results"]
                regions = results["regions"]
                # Write data for each species and region
                for region in regions:
                    for species, stats in species_results.items():
                        if region.split('.')[0] == 'hg38':
                            writer.writerow([
                                region, 
                                species, 
                                stats["conserved_bases"], 
                                stats["sequence_length"], 
                                f"{stats['similarity'] * 100:.2f}"
                            ])

# Read in HAQERs
regions = read_bed_file('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/git/language_evo/manuscript/supplemental_materials/HAQER.hg38.sorted_autosomes.bed')
print(f"There are {len(regions)} regions of interest in this annotation")

# Loop over HAQERs
print("\n\n\n\n")
print(f"Working on chromosome: {i}")
maf_file = '/Dedicated/jmichaelson-genome/Zoonomia/447-way_mammalian_alignment_2023/maf/chr' + str(i) + '.maf.gz' ## input multiple sequence alignment
# maf_file = '/Dedicated/jmichaelson-genome/Zoonomia/447-way_mammalian_alignment_2023/tst.maf.gz' ## input multiple sequence alignment
maf_file_subset = '/Dedicated/jmichaelson-genome/Zoonomia/447-way_mammalian_alignment_2023/tmp/chr' + str(i) + '.HAQERs.maf' ## maf subset to regions of interest
output_file = "/Dedicated/jmichaelson-genome/Zoonomia/447-way_mammalian_alignment_2023/tmp/HAQER_hg38_sequence_similarity_results.chr" + str(i) + ".csv"  # Output CSV file
## subset to regions of interest on this chromosome
chromosome = "chr" + str(i)
roi = [region for region in regions if region[0] == chromosome]
print(f"There are {len(roi)} regions of interest on this chromosome")
subset_maf_by_species(maf_file, roi, maf_file_subset, 'hg38')
## loop over regions of interest and compute sequence similarity
save_results_to_csv(maf_file_subset, output_file)
