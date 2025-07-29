from Bio import AlignIO
import numpy as np
import pandas as pd

## read in HAQER bed file
# data = pd.read_csv ("/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.bed", sep = '\t', names = ["chr", "start", "end", "id"])
data = pd.read_csv ("/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/complete_primate_genomes-Yoo-Nature2025/HAQERs_v2.hg38.bed", sep = '\t', names = ["chr", "start", "end"])

## loop over each chromosome to get sequence of interest
for chr in range(1, 23):
    chr = str(chr)
    ## filepath to alignment
    mafpath = "/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/5wayPrimate/chr" + chr + ".fa"
    print(f'\n\n\n\n')
    print("========================================================")
    print(f"Working on chromosome: {chr}")
    ## read in alignment
    aln = AlignIO.read(mafpath, 'fasta')
    print("Number of rows (species): %i" % len(aln))
    # dir(aln)
    # print(aln[:,1:100])
    ## extract sequences 
    hg38=aln[0,:]
    panpan = aln[1,:]
    pantro = aln[2,:]
    gor = aln[3,:]
    ponabe = aln[4,:]
    ##
    # dir(hg38)
    # dir(hg38.seq)
    ## count up n insertions into hg38 and see if removing those leaves sequence length of known hg38 chr length
    print(f'Number of gaps in hg38: {hg38.seq.count('-')} (these will be dropped)')
    # hg38_chr1_len = 248956422
    # len(hg38.seq) - hg38.seq.count('-') ## alignment length without gaps
    # hg38_chr1_len ## actual length
    ## find all insertions into hg38 so the coordinates match
    def find(s, ch):
        return [i for i, ltr in enumerate(s) if ltr == ch]
    drop_idx = find(hg38.seq, "-")
    # Convert string to a list to make it mutable
    string_list = list(hg38.seq)
    string_list_panpan = list(panpan.seq)
    string_list_pantro = list(pantro.seq)
    string_list_gor = list(gor.seq)
    string_list_ponabe = list(ponabe.seq)
    # Use list comprehension to exclude the elements at the removal indices
    # filtered_string_list = [char for i, char in enumerate(string_list) if i not in drop_idx]
    # filtered_string_list_panpan = [char for i, char in enumerate(string_list_panpan) if i not in drop_idx]
    # filtered_string_list_pantro = [char for i, char in enumerate(string_list_pantro) if i not in drop_idx]
    # filtered_string_list_gor = [char for i, char in enumerate(string_list_gor) if i not in drop_idx]
    # filtered_string_list_ponabe = [char for i, char in enumerate(string_list_ponabe) if i not in drop_idx]
    ## extract sequence of interest from the MSA for each species and remove hg38 gap positions
    arr = np.array(string_list)
    new_arr = np.delete(arr, drop_idx)
    arr = np.array(string_list_panpan)
    new_arr_panpan = np.delete(arr, drop_idx)
    arr = np.array(string_list_pantro)
    new_arr_pantro = np.delete(arr, drop_idx)
    arr = np.array(string_list_gor)
    new_arr_gor = np.delete(arr, drop_idx)
    arr = np.array(string_list_ponabe)
    new_arr_ponabe = np.delete(arr, drop_idx)
    print(f'Length of hg38 sequence: {len(new_arr)}')
    print(f'Length of orangutaun sequence: {len(new_arr_ponabe)}')
    ## read in HAQERs are reformat
    chrs = "chr" + chr
    df = data[data['chr'] == chrs] 
    df_sort = df.sort_values("start")
    df_sort["start"] = df_sort["start"] - 10001
    df_sort["end"] = df_sort["end"] + 10000
    ## subset full chromosome to alleles of interest
    subset_array = np.concatenate([new_arr[start:end] for start, end in zip(df_sort['start'], df_sort['end'])])
    subset_array_panpan = np.concatenate([new_arr_panpan[start:end] for start, end in zip(df_sort['start'], df_sort['end'])])
    subset_array_pantro = np.concatenate([new_arr_pantro[start:end] for start, end in zip(df_sort['start'], df_sort['end'])])
    subset_array_gor = np.concatenate([new_arr_gor[start:end] for start, end in zip(df_sort['start'], df_sort['end'])])
    subset_array_ponabe = np.concatenate([new_arr_ponabe[start:end] for start, end in zip(df_sort['start'], df_sort['end'])])
    ## check number of extracted alleles matches HAQER length
    print(f'Length of expected HAQER sequence on this chromosome: {sum(df_sort["end"] - df_sort["start"])}')
    print(f'Length of actual hg38 HAQER sequence on this chromosome: {len(subset_array)}')
    ## get genome positions
    index_values = [i for start, end in zip(df_sort['start'], df_sort['end']) for i in range(start, end)]
    index_values_mapped = new_list = [x + 1 for x in index_values]
    # Create a new DataFrame with all extracted info and save to file
    new_df = pd.DataFrame({'chr': chrs, 'pos_hg38': index_values_mapped, "human_allele": subset_array, "bonobo_allele": subset_array_panpan, "chimp_allele": subset_array_pantro, "gorilla_allele": subset_array_gor, "orangatuan_allele": subset_array_ponabe})
    outf="/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/5wayPrimate/HAQER_alleles/" + chrs + "_v2.csv"
    new_df.to_csv(outf)

# import subprocess
# import sys
# import sys
# import subprocess
# def install(package):
#     subprocess.check_call([sys.executable, "-m", "pip", "install", package])
# install('pandas')
