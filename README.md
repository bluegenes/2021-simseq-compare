# Evaluate sourmash distance estimation using simulated sequences

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bluegenes/2021-simseq-compare/main)

**click on the notebooks folder, then choose the notebook you want to run.**


## Simulated Sequences 

This repo uses a series of simulated sequences from [Criscuolo 2020](https://doi.org/10.12688/f1000research.26930.1) to evaluate distance estimation using sourmash.
The simulated sequences are archived on [Zenodo](https://zenodo.org/record/4034462) and are nucleotide sequences containing "--" stretches to represent indels.
Each fasta is provided in a column of a tsv.xz file with the simulation information included. As stated in the [Zenodo](https://zenodo.org/record/4034462) reference:

"""
This repository contains 24,000 pairs of nucleotide sequences (and associated parameters) that have been simulated for testing alignment-free genome distance estimates. Given an evolutionary distance d varying from 0.05 to 1.00 nucleotide substitutions per character (step = 0.05), the program INDELible was used to simulate the evolution of 200 nucleotide sequence pairs with d substitution events per character under the models GTR and GTR+Γ. Each model was adjusted with three different equilibrium frequencies:

f1: equal frequencies, i.e. freq(A) = freq(C) = freq(G) = freq(T) = 0.25,
f2: GC-rich, i.e. freq(A) = 0.1, freq(C) = 0.3, freq(G) = 0.4, freq(T) = 0.2,
f3: AT-rich, i.e. freq(A) = freq(T) = 0.4, freq(C) = freq(G) = 0.1.
For each simulated sequence pair, model parameters (i.e. GTR: six relative rates of nucleotide substitution; GTR+Γ: six rates and one Γ shape parameter) were randomly drawn from 142 sets of parameters derived from real-case data (see file GTR.params.trees.tsv at https://zenodo.org/record/4034261). **Initial sequence length was 5 Mbs**, and an indel rate of 0.01 was set with indel length drawn from [1, 50000] according to a Zipf distribution with parameter 1.5 (see INDELible manual).
"""

I used the following names to generate a unique name for each simulated pair of reads: `data-<DISTANCE>-<NUCL_SUBSTITUION>-<EVOLUTIONARY-MODEL>`
  
  - Nucleotide Substitution:
    - f1
    - f2
    - f3
  - Evolutionary Models:
    - `nogam`: GTR (six relative rates of nucleotide substitution)
    - `gamma`: GTR+Γ (six rates and one Γ shape parameter)
  - Evolutionary distances:
    - d0.05 --> d1.0 (stepping by 0.05)
  - Seed value provided to INDELIBLE:
    - seedXX

So, `data-d0.05-f1-nogam-seed36` represents a sequence simulation at 0.05 evolutonary distance, using equal frequence nucleotide substitution, without using a Γ shape parameter.

I generated a single csv with all comparison names and the corresponding simulation info: [simreads-info.csv.gz](https://github.com/bluegenes/2021-simseq-compare/raw/main/simreads-info.csv.gz)

## Analysis Workflow

### Analysis workflow:
- Indels (`--`) were removed (adjacent sections were concatenated), and fasta files were written for each simulated sequence.
- **For DNA comparisons:**
    - each fasta was sketched with `sourmash v4.0`, dna ksizes `21`,`31`,`51`, scaled=`1` (retains all unique k-mers)
    - sketches from sequence pairs were compared using `compare-paired-sequences.v2.py` for each set of simulation parameters, for scaled values: `1`,`100`,`1000`,`2000`
    - **results** from all simulation parameters were aggregated into a single results file, [simreads-compare.dnainput.csv.gz](https://osf.io/xn7vt/download)
- **For protein comparisons:**
    - PRODIGAL translation:
      - each fasta was translated using `PRODIGAL v2.6.3`
      - prodigal-translated reads were sketched with `sourmash sketch protein` in `sourmash v4.0`, with protein k-sizes `7-12`, scaled=`1` (retains all unique k-mers)
      - sketches from sequence pairs were compared using `compare-paired-sequences.v2.py` for each set of simulation parameters, for scaled values: `1`,`100`,`1000`,`2000`
      - **results** from all simulation parameters were aggregated into a single results file, `simreads-compare.prodigal.csv.gz` (to be posted)
    - 6-frame translation:
      - DNA fastas were sketched with `sourmash sketch translate` (6-frame translation) in `sourmash v4.0`, with protein k-sizes `7-12`, scaled=`1` (retains all unique k-mers)
      - sketches from sequence pairs were compared using `compare-paired-sequences.v2.py` for each set of simulation parameters, for scaled values: `1`,`100`,`1000`,`2000`
      - **results** from all simulation parameters were aggregated into a single results file, `simreads-compare.translate.csv.gz` (to be posted)

### Results format:

Each `simreads-compare*.csv.gz` file consists of the following columns:

- `comparison_name` - name containing all simulation information, e.g. `data-d0.05-f1-nogam-seed36`
- `sig1_name` - comparison name with `-seq1`
- `sig2_name` - comparison name with `-seq2`
- `alphabet` - alphabet used for sketching (`nucleotide`/`dna`, `protein`, `dayhoff`, or `hp`)
- `ksize` - kmer size used for sketching. Corresponds to the alphabet (k=21 for nucleotide is 21 nucleotides; k=7 for protein is 7 amino acids)
- `scaled` - scaled value. Initial scaled=1 sketches were downsampled to this value for this comparison.
- `jaccard` - `sourmash` estimate of the Jaccard Index
- `max_containment` - `sourmash` estimate of Maximum Containment (maximum Containment Index between both Containment directions)
- `sig1_containment` - `sourmash` estimate of the Containment Index for sig1
- `sig2_containment` - `sourmash` estimate of the Containment Index for sig2
- `sig1_hashes` - number of unique k-mers in sig1 sketch 
- `sig2_hashes` - number of unique k-mers in sig2 sketch
- `num_common` - number of intersected k-mers between sig1,sig2
- `alpha-ksize` - combined column for alphabet and ksize

  
> Note that many of these columns can be regenerated from other columns, I just found them useful for plotting, comparisons, etc.
> I may also upload csvs containing a pared down set of columns.

_Note to self: this repo was split from 2020-dist-compare; version 1 analysis remains there._
