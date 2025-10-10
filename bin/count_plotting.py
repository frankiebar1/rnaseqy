#!/usr/bin/env python

## Imports needed
import argparse
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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

def plot_pca(counts_file, plot_file):

    # Load counts data
    df = pd.read_csv(counts_file, sep="\t")
    df.set_index("Geneid", inplace=True)
    
    # Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df.T)
    scaled_data = pd.DataFrame(scaled_data, index=df.columns, columns=df.index)
    scaled_data = scaled_data.dropna(axis=1)  # Drop genes with NaN values after scaling

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)

    # Create a DataFrame for PCA results
    pca_df = pd.DataFrame(data=pca_result, columns=["PC1", "PC2"])
    pca_df.index = df.columns

    # Plot PCA
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_df["PC1"], pca_df["PC2"])

    for i, sample in enumerate(pca_df.index):
        plt.text(pca_df["PC1"][i], pca_df["PC2"][i], sample)

    plt.title("PCA of FeatureCounts Data")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
    plt.grid()
    plt.savefig(plot_file)
    print("PCA plot saved to", plot_file)


## Main
parser = argparse.ArgumentParser()

parser.add_argument("--indir", required=True, help="Directory with *.featureCounts.tsv files")
parser.add_argument("--outdir", default="featureCounts_matrix.tsv", help="Output merged matrix file")
parser.add_argument("--plot", default="featureCounts_PCA.png", help="Output PCA plot file")
args = parser.parse_args()

print("Script is used")
output_file = os.path.join(args.outdir, "merged_counts.tsv")
output_plot = os.path.join(args.outdir, "pca.png")

out_file = merge_featurecounts(args.indir, output_file)
plot_pca(out_file, output_plot)