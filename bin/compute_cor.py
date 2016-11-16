from numpy import *
from scipy.stats import *
from matplotlib import pyplot
import os, argparse, sys, gzip, collections

#########################
# Load data sources
#########################

labels = ["miRanda",'miRDB','targetscan', 'starmirdb',"PITA" ]
ifiles = [ "master_files/%s.tsv"%x for x in labels ]
no_files = len(labels)

resources = collections.defaultdict(list)

#for idx, ifile in enumerate(args.resource_masterfiles):
for idx, ifile in enumerate(ifiles):
    handle = gzip.open(ifile) if ifile.endswith(".gz") else open(ifile)
    data = [ x.strip("\n").split("\t") for x in handle ]
    for record in data:
        key = tuple(record[:3])
        val = float(record[5])
        if not key in resources:
            resources[key] = [ 0.0 for x in range(no_files) ]
        resources[key][idx] = val

resources = array([ x for x in resources.values()],dtype=float)

#########################
# Determine intersection, compute correlations, and plot
#########################

spearman_cor = array([ [None for x in range(no_files)]
                       for y in range(no_files) ], dtype=float)
spearman_pval = array([ [None for x in range(no_files)]
                       for y in range(no_files) ], dtype=float)

pearson_cor = array([ [None for x in range(no_files)]
                       for y in range(no_files) ], dtype=float)
pearson_pval = array([ [None for x in range(no_files)]
                       for y in range(no_files) ], dtype=float)

# Plot individual distributions?

SKIP_PROCESSED = True
nbins = 1000
for i in range(no_files-1):
    for j in range(i+1,no_files):
        idx2use =  logical_and( resources[:,i] > 0,  resources[:,j] > 0 )
        pearson_cor[i,j], pearson_pval[i,j] = pearsonr( resources[idx2use,i],resources[idx2use,j] )
        spearman_cor[i,j], spearman_pval[i,j] = spearmanr( resources[idx2use,i],resources[idx2use,j] )

savetxt( "spearman_cor.txt" , spearman_cor, delimiter="\t"  )
savetxt( "spearman_pval.txt" , spearman_pval, delimiter="\t"  )
savetxt( "pearson_cor.txt" , pearson_cor, delimiter="\t"  )
savetxt( "pearson_pval.txt", pearson_pval, delimiter="\t"  )
