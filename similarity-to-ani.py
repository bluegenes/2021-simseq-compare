###############################################
# Convert sourmash distances into ANI estimates
# Author: N. Tessa Pierce, ntpierce@ucdavis.edu
###############################################

import os
import sys
import argparse

import numpy as np
import pandas as pd


# jaccard, containment --> mash distance ANI, AAI
def similarity_to_mashD_ANI(similarity, ksize, b1=1.0, b2=1.0, return_ANI=False):
    # orig used with jaccard; trying with jaccard, max containment, anchor containment
    if similarity == 0:
        return np.nan
   # first, get proportion of observed differences
    p = 1 - np.power(2*similarity/(similarity + 1),(1/float(ksize)))
    # apply poisson correction
    d = -(b1*np.log((1-p)/b2))
    # if want ANI, return 1-d
    if return_ANI:
        return 1-d
    return d


## use mutation rate intervals equations instead
#python ../mutation-rate-intervals/r1-from-minhash-jaccard.py L=5000000 k=21 C=0.95 J=0.209971 S=5000
def containment_to_ANI_mri(containment, ksize, num_unique_kmers, conf=0.95, return_ANI=False):
    # to be implemented
    d = np.nan
    if return_ANI:
        return 1-d
    return d




def estimate_ANI_AAI(row):
    # use mash distance equations
    row["jaccard_ANI-AAI"] = similarity_to_mashD_ANI(row["jaccard"], row["ksize"], return_ANI=True)
    row["mc_ANI-AAI"] = similarity_to_mashD_ANI(row["max_containment"], row["ksize"], return_ANI=True)
    row["ac_ANI-AAI"] = similarity_to_mashD_ANI(row["anchor_containment"], row["ksize"], return_ANI=True)
    # now estimate ANI/AAI using mutation-rate-intervals

    # get average unique k-mers
    #sigA_unique = row["anchor_hashes"] * row['scaled']
    #sigB_unique = row["query_hashes"] * row['scaled']
    #avg_unique_kmers= (sigA_unique + sigB_unique) /2
    #row["jaccard_ANI-AAI_mri"] = similarity_to_mashD_ANI(row["jaccard"], row["ksize"], avg_unique_kmers, return_ANI=True)
    #row["mc_ANI-AAI_mri"] = similarity_to_mashD_ANI(row["max_containment"], row["ksize"], avg_unique_kmers, return_ANI=True)
    #row["ac_ANI-AAI_mri"] = similarity_to_mashD_ANI(row["anchor_containment"], row["ksize"], avg_unique_kmers, return_ANI=True)
    return row



def main(args):
    # read in distance csv
    sim_cols = ["jaccard", "anchor_containment", "max_containment"]
    simDF = pd.read_csv(args.similarity_csv, sep=",")
    simDF[sim_cols] = simDF[sim_cols].replace({0:np.nan})
    simDF = simDF.apply(estimate_ANI_AAI, axis=1)
    simDF.to_csv(args.output_csv, index=False)

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--similarity-csv")
    p.add_argument("--output-csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
