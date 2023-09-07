"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Usage: Use pyANI to Compare distance estimation on pairs of unaligned nucleotide sequences
Installation: mamba env install -f envs/simseq-env.yml; conda activate simseq
Run: snakemake -s pyani-compare.snakefile -n
"""
#####
####

import pandas as pd

configfile: conf/pyani-compare.yml
out_dir = config["output_dir"]
logs_dir = out_dir + "/logs"
basename = config.get("basename", "simreads-compare")


# just do nucl input
FASTA_TYPES = ["dnainput"]

rule all: 
    input:
        expand(os.path.join(out_dir, "compare", "{siminfo}.{fasta_type}.fastalist.txt"), siminfo=simulation_info, fasta_type=FASTA_TYPES),
        os.path.join(out_dir, f"{basename}.dnainput.csv.gz")
        #expand(os.path.join(out_dir, "compare", "{siminfo}.{fasta_type}.compare.csv.gz"), siminfo=simulation_info, fasta_type = FASTA_TYPES)

rule download_tsv:
    output: os.path.join(out_dir,"data/{siminfo}.tsv.xz")
    params:
        url=lambda w: f"https://zenodo.org/record/4034462/files/{w.siminfo}.tsv.xz?download=1"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    shell:
        """
        wget -O {output} {params.url}
        """

# extract & aggregate *just simulation info* from multiple tsv.xz files --> single csv file
rule write_siminfo_filelist:
    input: ancient(expand(os.path.join(out_dir,"data/{siminfo}.tsv.xz"), siminfo=simulation_info))
    output: os.path.join(out_dir, "data/simulation-info.filelist.txt")
    run: 
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


rule aggregate_simulation_info:
    input: os.path.join(out_dir, "data/simulation-info.filelist.txt")
    output: os.path.join(out_dir, "data/simulation-info.csv.gz")
    conda: "envs/pdist-env.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=600,
    shell:
        """
        python aggregate-simulation-info.py {input} {output}
        """

# extract fastas from the simulation tsv.xz files
rule simreads_to_fasta_f1:
    input: ancient(os.path.join(out_dir,"data/{siminfo}.tsv.xz"))
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f1_seednums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """

rule simreads_to_fasta_f2:
    input: ancient(os.path.join(out_dir,"data/{siminfo}.tsv.xz"))
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f2_seednums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """

rule simreads_to_fasta_f3:
    input: ancient(os.path.join(out_dir,"data/{siminfo}.tsv.xz"))
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f3_seednums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """

localrules: write_fasta_filelist
rule write_fasta_filelist:
    input: 
        lambda w: ancient(expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = seedinfo[w.siminfo], seq=[1,2])),
    output: os.path.join(out_dir, "compare", "{siminfo}.dnainput.fastalist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    scaled = min(scaled) #take minimum value of scaled list
    return f"{ks},scaled={scaled},abund"

rule sourmash_sketch:
    input:
        ancient(os.path.join(out_dir, "data/simreads", "{siminfo}-seed{seed}-seq{seq}.fasta")),
    output:
        os.path.join(out_dir, "signatures", "{siminfo}-seed{seed}-seq{seq}.sig"),
    params:
        sketch_params=make_param_str(config["ksize"], config["scaled"]),
        signame = lambda w: f"{w.siminfo}-seed{w.seed}-seq{w.seq}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch", "{siminfo}-seed{seed}-seq{seq}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch", "{siminfo}-seed{seed}-seq{seq}.sketch.benchmark")
    conda: "envs/smash-compare.yml"
    shell:
        """
        sourmash sketch dna -p {params.sketch_params} -o {output} --name {params.signame:q} {input} 2> {log}
        """

localrules: signames_to_file
rule signames_to_file:
    input:
        lambda w: expand(os.path.join(out_dir, "signatures", "{{siminfo}}-seed{seed}-seq{seq}.sig"), seed = seedinfo[w.siminfo], seq=[1,2]),
    output: os.path.join(out_dir, "compare", "{siminfo}.dnainput.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule sourmash_compare:
    input:
        #simulation_info = ancient(os.path.join(out_dir,"data/{siminfo}.tsv.xz")),
        #simulation_info = os.path.join(out_dir, "data/simulation-info.csv.gz"),
        siglist = os.path.join(out_dir, "compare", "{siminfo}.dnainput.siglist.txt")
    output:
        os.path.join(out_dir, "compare", "{siminfo}.{alphabet}-k{ksize}.dnainput.compare.csv.gz"),
    params:
        sig_dir = os.path.join(out_dir, "signatures"),
        scaled_cmd = " --scaled " + " --scaled ".join(map(str, config["scaled"]))
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{siminfo}.{alphabet}-k{ksize}.dnainput.log")
    benchmark: os.path.join(logs_dir, "compare", "{siminfo}.{alphabet}-k{ksize}.dnainput.benchmark")
    conda: "envs/smash-compare.yml"
    shell:
        """
        python compare-paired-sequences.v2.py \
               --siglist {input.siglist} --sigdir {params.sig_dir} \
               --alphabet {wildcards.alphabet} --ksize {wildcards.ksize} \
               {params.scaled_cmd} --output-csv {output} 2> {log}
        """
 #--simreads-info-csv {input.simulation_info} \
localrules: aggregate_sourmash_compare
rule aggregate_sourmash_compare:
    input:
        expand(os.path.join(out_dir, "compare", "{siminfo}.{alpha}-k{k}.dnainput.compare.csv.gz"), siminfo=simulation_info, alpha=config["alphabet"], k=config["ksize"]),
    output:
        os.path.join(out_dir, "{basename}.dnainput.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] = aggDF["alphabet"] + "-" + aggDF["ksize"].astype(str)
        aggDF.to_csv(str(output), index=False)

