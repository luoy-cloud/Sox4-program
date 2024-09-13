#%% Read arguments
# [sc counts] [sc types] [spatial counts] [spatial spots] [spatial genes] [spatial locations] [output file]
import sys, argparse
msg = "Command-line wrapper for SpaDecon"
# Initialize parser
parser = argparse.ArgumentParser(description = msg)
parser.add_argument('sc_counts', help = "Gene expression data (cells x genes) from scRNA-seq (csv)")
parser.add_argument('sc_types', help = "Cell type of each cell in scRNA-seq data (csv)")
parser.add_argument('spatial_counts', help = "Gene expression data (spots x genes) from spatial dataset (mtx)")
parser.add_argument('spot_names', help = 'List of spot names in spatial dataset (csv)')#
parser.add_argument('gene_names', help = 'List of gene names in spatial dataset (csv)')#
parser.add_argument('spatial_locations', help = "Locations of each spot in spatial dataset (csv)")#
parser.add_argument('output_file', help = "Where to export deconvolution results")

args = parser.parse_args()
print(args.sc_counts)
print(args.spatial_counts)

#%% Load packages
import SpaDecon as spd
import scanpy as sc
import pandas as pd
import numpy as np
from skimage import io
import os
import scipy as sp

#%% Load single-cell data
#Read scRNA-seq GE data (rows = cells, columns = genes)
sc_ge = pd.read_csv(args.sc_counts, index_col = 0)

#Read scRNA-seq cell-type labels
sc_types = pd.read_csv(args.sc_types, index_col = 0)

#Convert scRNA-seq GE data to AnnData object
adata_sc = sc.AnnData(sc_ge)
print('sc data loaded')
#Insert cell-type labels into "celltype" column adata_sc.obs
adata_sc.obs['celltype'] = sc_types['celltype']
print('celltypes assigned')

#%% Load spatial data
#Read SRT GE data (rows = cells, columns = genes)
srt_ge = sp.io.mmread(args.spatial_counts)
srt_ge = sp.sparse.csc_matrix(srt_ge)
srt_spots = pd.read_csv(args.spot_names, header = None).set_index(0)
srt_spots.index.name = None
srt_genes = pd.read_csv(args.gene_names, header = None).set_index(0)
srt_genes.index.name = None

#Convert SRT GE data to AnnData object
adata_srt = sc.AnnData(srt_ge, obs=srt_spots, var=srt_genes)

#Read file with SRT spatial locations
locations = pd.read_csv(args.spatial_locations,header=None,index_col=0) 
locations = locations.loc[adata_srt.obs.index]

#%% Run SpaDecon
clf = spd.SpaDecon()
clf.deconvolution(source_data=adata_sc, target_data=adata_srt, spatial_locations=locations, histology=False)
spadecon_proportions = clf.props
spadecon_proportions.to_csv(args.output_file) 