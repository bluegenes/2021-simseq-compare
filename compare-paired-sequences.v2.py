import os
import sys
import argparse
import glob
import pprint

import pandas as pd

import screed
import sourmash
from sourmash import load_one_signature
from sourmash.sourmash_args import load_file_as_signatures

from collections import defaultdict, namedtuple

CompareResult = namedtuple('CompareResult',
                           'comparison_name, sig1_name, sig2_name, alphabet, ksize, scaled, jaccard, max_containment, sig1_containment, sig2_containment, sig1_hashes, sig2_hashes, num_common')

def compare_sigs(sigA, sigB, comparison_name, alpha, ksize, scaled):
    # handle scaled
    if scaled != sigA.minhash.scaled:
        sigA.minhash = sigA.minhash.downsample(scaled=scaled)
    if scaled != sigB.minhash.scaled:
        sigB.minhash = sigB.minhash.downsample(scaled=scaled)
    # compare
    sigA_numhashes = len(sigA.minhash.hashes)
    sigB_numhashes = len(sigB.minhash.hashes)
    intersect_numhashes = sigA.minhash.count_common(sigB.minhash)
    jaccard = sigA.jaccard(sigB)
    containA = sigA.contained_by(sigB)
    containB = sigA.contained_by(sigA)
    max_contain = sigA.max_containment(sigB)
    return CompareResult(comparison_name, str(sigA).split(" ")[0], str(sigB).split(" ")[0], alpha, ksize, scaled, jaccard, max_contain, containA, containB, sigA_numhashes, sigB_numhashes, intersect_numhashes)

def main(args):
    ksize=args.ksize
    scaled_vals=args.scaled
    alphabet=args.alphabet
    if alphabet == "nucleotide":
        moltype = "DNA"
    else:
        moltype = alphabet

    # names, info for simulated fasta pairs to compare
    #compare_info = pd.read_csv(args.simreads_info_csv, dtype=str, sep=",", header=0)
    # read list of sig paths
    # build dictionary signame::sigfile so we don't need to load all sigs into memory at once (only load when needed)
    siglist = [x.rstrip() for x in open(args.siglist)]
    sigD={}
    compare_names = set()
    for sigF in siglist:
        name = os.path.basename(sigF).rsplit(args.sig_suffix)[0]
        # specific to these paired seqs
        compare_name = name.rsplit("-seq")[0]
        compare_names.add(compare_name)
        if not os.path.exists(sigF):
            full_sigF = os.path.join(args.sigdir, sigF)
            if not os.path.exists(full_sigF):
                print(f"sig {name} cannot be found at {sigF} or within sigdir {args.sigdir}")
                continue
            else:
                sigF=full_sigF
        sigD[name] = sigF

    # get comparison list
    #names = list(compare_info["name"])
    seq_comparisons = []
    for n, comparison_name in enumerate(compare_names):
        if n !=0 and n % 50 == 0:
            print(f"... assessing {n}th pair, {name}\n")
        # for paired simulated sequences, sigfiles should just be f"{name}-seq1" and f"{name}-seq2"
        seq1_name = f"{comparison_name}-seq1"
        seq2_name = f"{comparison_name}-seq2"
        # select and load sigs
        #selector = load_file_as_signatures(sigD[seq1_name], ksize=ksize, select_moltype=moltype)
        #sig1 = next(selector)
        #selector = load_file_as_signatures(sigD[seq2_name], ksize=ksize, select_moltype=moltype)
        #sig2 = next(selector)
        sig1 = load_one_signature(sigD[seq1_name], ksize=ksize, select_moltype=moltype)
        sig2 = load_one_signature(sigD[seq2_name], ksize=ksize, select_moltype=moltype)

        # can only compare at largest scaled
        orig_scaled = max(sig1.minhash.scaled, sig2.minhash.scaled) # can only compare at largest scaled

        # if scaled is not specified, use original scaled
        if not scaled_vals:
            scaled_vals=[orig_scaled]

        # compare at each scaled val
        for sc in scaled_vals:
            if sc < orig_scaled:
                print(f"Can't downsample: desired scaled {sc} is smaller than original scaled, {orig_scaled}. Ignoring scaled {sc}...")
                continue
                # store signature info
            comparison = compare_sigs(sig1, sig2, comparison_name, alphabet, ksize, sc)
            seq_comparisons.append(comparison)


    # convert path comparison info to pandas dataframe
    comparisonDF = pd.DataFrame.from_records(seq_comparisons, columns = CompareResult._fields)

    # print to csv
    comparisonDF.to_csv(args.output_csv, index=False)
    print(f"done! simread comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    #p.add_argument("--simreads-info-csv", default="simreads-info.csv.gz")
    p.add_argument("--siglist", default="simreads.signatures.txt")
    p.add_argument("--sigdir", default="signatures")
    p.add_argument("--sig-suffix", default=".sig")
    p.add_argument("--alphabet", default="protein")
    p.add_argument("--ksize", default=10, type=int)
    p.add_argument("-s", "--scaled", action="append", type=int, help= "provide additional scaled values for downsampling")
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
