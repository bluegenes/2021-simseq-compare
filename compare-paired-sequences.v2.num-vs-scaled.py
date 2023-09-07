import os
import sys
import argparse
import glob
import pprint
import csv

import pandas as pd

import screed
import sourmash
from sourmash import load_one_signature
from sourmash.sourmash_args import load_file_as_signatures

from collections import defaultdict, namedtuple

NumScaledCompareResult = namedtuple('NumScaledCompareResult',
                           'comparison_name, sig1_name, sig2_name, alphabet, ksize, scaled, scaled_jaccard, scaled_intersect, num, num_jaccard, num_intersect')

def compare_mh(scaled_mhA, scaled_mhB, num_mhA, num_mhB, A_name, B_name, comparison_name, num, scaled, alpha, ksize):
    # compare scaled mh
    scaled_intersect = scaled_mhA.count_common(scaled_mhB)
    scaled_jaccard = scaled_mhA.jaccard(scaled_mhB)
    num_intersect    = num_mhA.count_common(num_mhB)
    num_jaccard = num_mhA.jaccard(num_mhB)
    return NumCompareResult(comparison_name, A_name, B_name, alpha, ksize, scaled, scaled_jaccard, scaled_intersect, num, num_jaccard, num_intersect)

def main(args):
    ksize=args.ksize
    num=args.num
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
    num_vals = num
    num_comparisons = len(compare_names)
    output_fieldnames = NumCompareResult._fields # may need to be list?
    with open(args.output_csv, 'w') as outF:
        w = csv.DictWriter(outF, fieldnames=output_fieldnames)

        for n, comparison_name in enumerate(compare_names):
            if n % 10 == 0:
                print(f"... assessing {n}th pair, {name} of {num_comparisons} total\n")
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

            # get all hashvals for these sigs
            sig1_hashvals = sig1.minhash.hashes.keys()
            sig2_hashvals = sig2.minhash.hashes.keys()
            sig1_name = str(sig1).split(" ")[0]
            sig2_name = str(sig2).split(" ")[0]

            minimum_hashes= min(len(sig1_hashvals), len(sig2_hashvals)) # can only compare at minimum num hashes between the two

            orig_scaled = max(sig1.minhash.scaled, sig2.minhash.scaled) # can only compare at largest scaled
            ## if num is not specified, use original num
            #if not num_vals:
            #    num_vals=[orig_num]

            # compare at each num val
            seq_comparisons = []
            for sc in scaled_vals:
                if sc < orig_scaled:
                    print(f"Can't downsample: desired scaled {sc} is smaller than original scaled, {orig_scaled}. Ignoring scaled {sc}...")
                    continue
                    # store signature info

                # build new scaled minhashes (downsample)
                if scaled != sigA.minhash.scaled:
                    sc_mh1 = sigA.minhash.downsample(scaled=scaled)
                if scaled != sigB.minhash.scaled:
                    sc_mh2 = sigB.minhash.downsample(scaled=scaled)

                numhashes_1 = len(sc_mh1.hashes.keys())
                numhashes_2 = len(sc_mh2.hashes.keys())
                import pdb;pdb.set_trace()
                average_num = round(np.mean(numhashes_1,numhashes_2))

                # build new num minhashes
                mh1_num = sourmash.MinHash(n=average_num, ksize=ksize)
                mh2_num = sourmash.MinHash(n=average_num, ksize=ksize)
                mh1_num.add_many(sig1_hashvals)
                mh2_num.add_many(sig2_hashvals)

                comparison = compare_mh(mh1_num, mh2_num, sig1_name, sig2_name, comparison_name, nu, alphabet, ksize)
                seq_comparisons.append(comparison)

           # print(f"writing comparison {comparison_name}\n")
            for c in seq_comparisons:
                w.writerow(c._asdict())


    # print to csv
    #comparisonDF.to_csv(args.output_csv, index=False)
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
    p.add_argument("-s", "--num", action="append", type=int, help= "provide additional num values for downsampling")
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
