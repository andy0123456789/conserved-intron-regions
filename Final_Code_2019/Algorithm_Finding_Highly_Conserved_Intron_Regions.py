'''dataframe'''
import pandas as pd
'''array processing package with Python'''
import numpy as np
'''2-dimensional Python plotting library'''
from matplotlib import pyplot as plt
'''Python data visualization library based on matplotlib'''
import seaborn as sns
'''Python module for unix style pathnames'''
from glob import glob
'''dataframes'''
from pandas import DataFrame
'''module provides the Python interfaces to stream handling'''
import io
'''The library parses JSON into a Python dictionary or list'''
import json
'''Puthonic interface to the HDF5 binary data format'''
import h5py
'''Find peaks in a 1-D array with wavelet transformation.'''
from scipy.signal import medfilt,find_peaks_cwt, fftconvolve, general_gaussian

sns.set_style("white")
#_____________________________________________________________________________________________________________________________________________________
'''defines glob object which has all of the files to be under data folder and have filepath of "data/intronic_conservation_bedgraphs/huandy_*"
'''
files = glob("data/intronic_conservation_bedgraphs/huandy_*")
'''Initialize the empty h5 file'''

with h5py.File('data/aggregate_regions_bp_resolution2.h5', 'w') as h5file:

    '''For each bed file, read and append to the h5 collection'''
    for idx,f in enumerate(files):

        '''Log the current file # to keep track of progress'''
        if idx % 1000 == 0:
            print(idx, end= " ")

        region = pd.read_table(f,names=["chromosome","start","end","conservation"])
        
        '''Convert to correct datatypes'''
        region["start"] = region["start"].astype(int)
        region["end"] = region["end"].astype(int)
        region["conservation"] = region["conservation"].astype(float)

        '''Identify the region by the file name.
         -1 means taking the last element of each file path.'''
        region_id =f.split("/")[-1]
        
        '''Create separate datasets for start, end, and conservation'''
        
        starts = region["start"]
        ends = region["end"]
        conservations = region["conservation"]
        
        '''Expand conserved mini-regions into base-pair resolution array'''
        exploded_conservations = np.repeat(conservations,ends-starts)

        '''Input the base-pair resolution array into a new h5'''
        h5file.create_dataset(region_id,  data=exploded_conservations)
        def v(ranges):
    """
    Merge overlapping and adjacent ranges and yield the merged ranges
    in order. The argument must be an iterable of pairs (start, stop).

    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []
    """
    ranges = iter(sorted(ranges))
    current_start, current_stop = next(ranges)
    for start, stop in ranges:
        if start > current_stop:
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop
   UCSC_exons = pd.read_table("data/UCSC_exons.bed",names=["chrom","start","end","id","0_1","strand"])
exons_info = UCSC_exons[["chrom","start","end"]]

non_overlapping_exons = pd.DataFrame()

for chrom in list(set(exons_info["chrom"])):
    chrom_specific = exons_info[exons_info["chrom"]==chrom]
    
    exon_regions = list(zip(chrom_specific["start"],chrom_specific["end"]))
    
    exon_regions = np.array(list(merge_ranges(exon_regions)))
    
    specific_regions = pd.DataFrame(data=exon_regions,columns=["start","stop"])
    
    specific_regions["chrom"] = chrom
    
    non_overlapping_exons = non_overlapping_exons.append(specific_regions)
    
    


non_overlapping_exons = non_overlapping_exons.sort_values(by=["chrom","start"])
def has_exon(row):
    
    chrom = row["chrom"]
    start = row["start"]+1
    end = row["end"]-1
    
    chrom_specific = non_overlapping_exons[non_overlapping_exons["chrom"]==str(chrom)]
    
    exon_starts = np.array(chrom_specific["start"])
    exon_ends = np.array(chrom_specific["stop"])
    
    left_exon_end = np.searchsorted(exon_ends, start)
    right_exon_start = np.searchsorted(exon_starts, end)
    
    if left_exon_end == right_exon_start-1:
        return True
    else:
        return False
 intron_regions_df = pd.DataFrame()

chroms = np.asarray([x.split("_")[-3] for x in regions])
starts = np.asarray([int(x.split("_")[-2]) for x in regions])
ends = np.asarray([int(x.split("_")[-1]) for x in regions])

intron_regions_df["chrom"] = chroms
intron_regions_df["start"] = starts
intron_regions_df["end"] = ends

intron_regions_df["start"] = intron_regions_df["start"].astype(int)
intron_regions_df["end"] = intron_regions_df["end"].astype(int)

intron_regions_df = intron_regions_df.set_index(regions)
contains_exons = intron_regions_df.apply(has_exon, axis=1)

contains_exons.to_csv("contains_exons.csv")
with h5py.File('data/aggregate_regions_bp_resolution.h5', 'r') as h5file:
    regions = list(h5file.keys())
    '''Define count'''
    count = 0
    
    peak_df = pd.DataFrame()
    peaks_files = []
    
    chroms = np.asarray([x.split("_")[-3] for x in regions])
    starts = np.asarray([int(x.split("_")[-2]) for x in regions])
    ends = np.asarray([int(x.split("_")[-1]) for x in regions])
    
#     select_regions = np.where((chroms == "chr3") & (starts>170582665) & (starts<170588045))

    select_regions = np.arange(len(regions))

   # select_regions = np.random.choice(len(regions),100)

    
    print(len(select_regions))
    
    region_lengths = []
    
    for idx,region in enumerate(np.asarray(regions)[select_regions]):
        
        if idx%100==0:
            print(idx,end=" ")
        
        conservation_scores = h5file[region][:]
        
        region_length = len(conservation_scores)
        if region_length < 10**4:
        

            window_width = 9

            if len(conservation_scores) > 32:

                
                convolve_length = 16
                
                window = general_gaussian(64, p=0.5, sig=24)
                filtered = fftconvolve(window, conservation_scores)
                filtered = (np.average(conservation_scores) / np.average(filtered)) * filtered
                filtered = np.roll(filtered, -25)
                
                running_mean = np.convolve(conservation_scores, np.ones((convolve_length,))/convolve_length, mode='valid')
                peaks = find_peaks_cwt(running_mean,widths = np.arange(8,64),min_snr=2)

                if len(peaks) > 0:

                    peak_conservations = conservation_scores[peaks]
                    peaks = peaks[peak_conservations > 6]

                    if len(peaks) > 2 and len(peaks) < 8:
                        count+=1

                        for bp in peaks:
                            plt.axvline(bp)

                       # plt.step(np.arange(len(conservation_scores)),conservation_scores,c="grey")
                       # plt.step(np.arange(len(filtered))+16,filtered,c="orange")
                       # print(region)
                       # plt.show()

                        peaks_files.append(region)
                        
    print(len(peaks_files))

    peak_df["file_region"] = peaks_files
    peak_df.to_csv("peak_df.csv")
