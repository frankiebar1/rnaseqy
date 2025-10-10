#!/usr/bin/env python

## Imports needed
import argparse
import os
import glob
import pandas as pd


## Functions
# merge counts in one table
def merge_featurecounts(input_dir, output_file):
    # Collect files
    files = glob.glob(os.path.join(input_dir, "*.txt"))
    if not files:
        raise ValueError("No files found in", input_dir)
    
    merged = None
    for f in files:
        df = pd.read_csv(f, sep="\t")
        #print(df.head)
        df = df.iloc[:, [0, -1]]  # keep geneid + counts column
        sample = os.path.basename(f).replace(".featureCounts.txt", "")
        df.columns = ["Geneid", sample]
        merged = df if merged is None else merged.merge(df, on="Geneid", how="outer")

    #merged.fillna(0, inplace=True)
    merged.to_csv(output_file, sep="\t", index=False)
    print("merged")
    return output_file

# plot PCA


## Main
parser = argparse.ArgumentParser()

parser.add_argument("--indir", required=True, help="Directory with *.featureCounts.tsv files")
parser.add_argument("--out", default="featureCounts_matrix.tsv", help="Output merged matrix file")
parser.add_argument("--plot", default="featureCounts_PCA.png", help="Output PCA plot file")
args = parser.parse_args()

print("Script is used")
merge_featurecounts(args.indir, args.out)
