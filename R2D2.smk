# pileline for pre-process hires data on nanopore platform
#@author skelviper 20220721
# 

SAMPLES = [i.split(sep='.fastq.gz')[0] for i in os.listdir("./Rawdata")]
barcodes = "R2D2/barcodes.txt"
import pandas as pd 
BARCODES = pd.read_csv(barcodes).values.T[0].tolist()

configfile: "R2D2/config.yaml"

#localrules:
#    demultiplex,

rule all:
    input:
        expand("processed/demultiplex/{sample}.R1.{barcode}.fastq.gz", sample=SAMPLES, barcode=BARCODES),
        expand("processed/demultiplex/{sample}.R2.{barcode}.fastq.gz", sample=SAMPLES, barcode=BARCODES),
        expand("result/RNA_Res/counts.{type}.{genome}.format.tsv", type = ["gene","exon"],genome=["genome1","genome2","total"]),
    threads: 30
    shell:"""
        set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate py3; set -u;

        mkdir -p stat
        #bash ./hires3rd/scripts/generateStat.sh
        echo "All done!"

        set +u; conda deactivate; set -u;
    """

rule trim_ligation:
    input:
        R1 = "Rawdata/{sample}/{sample}_R1.fq.gz",
        R2 = "Rawdata/{sample}/{sample}_R2.fq.gz",
    output:
        trimmed_R1 = "processed/trimmed/{sample}.trim.R1_fastq.gz",
        trimmed_R2 = "processed/trimmed/{sample}.trim.R2_fastq.gz",
    threads:
        10
    conda:"envs/main_env.yaml"
    shell:"""

        cutadapt -j {threads} --discard-untrimmed -m 20 -g "NNNNNNNNGTAGGTGTGAGTGATGGTTGAGGTAGT;o=35;e=0.2;" --rename='{{r1.match_sequence}}_{{header}}' {input.R1} {input.R2} -o {output.trimmed_R1} -p {output.trimmed_R2}
    """    

rule demultiplex:
    input:
        trimmed_R1 = rules.trim_ligation.output.trimmed_R1,
        trimmed_R2 = rules.trim_ligation.output.trimmed_R2,
        barcodes = barcodes,
    output:
        demultiplex_R1 = expand("processed/demultiplex/{{sample}}.R1.{barcode}.fastq.gz", barcode=BARCODES),
        demultiplex_R2 = expand("processed/demultiplex/{{sample}}.R2.{barcode}.fastq.gz", barcode=BARCODES),
    params:
        prefixR1 = "processed/demultiplex/{sample}.R1.",
        prefixR2 = "processed/demultiplex/{sample}.R2.",
    threads: 10
    conda:"envs/demultiplex.yaml"
    shell:"""

        python3 ./R2D2/demultiplex.py -b {input.barcodes} -i {input.trimmed_R1} -o {params.prefixR1} -t {threads}
        python3 ./R2D2/demultiplex.py -b {input.barcodes} -i {input.trimmed_R2} -o {params.prefixR2} -t {threads}

    """

# demultiplex output order is not the same as input order, so we need to sort it
rule sortfastq:
    input:
        RNA_R1="processed/demultiplex/{sample}.R1.{barcode}.fastq.gz",
        RNA_R2="processed/demultiplex/{sample}.R2.{barcode}.fastq.gz"
    output:
        RNAsort_R1 = temp("processed/demultiplex/{sample}.R1.{barcode}.sort.fastq.gz"),
        RNAsort_R2 = temp("processed/demultiplex/{sample}.R2.{barcode}.sort.fastq.gz"),
    conda:  "envs/main_env.yaml"
    shell:"""

        zcat {input.RNA_R1} | seqkit sort -n - | sed  "s/@[ATCGN]\+\_/@/g" | gzip > {output.RNAsort_R1}
        zcat {input.RNA_R2} | seqkit sort -n - | sed  "s/@[ATCGN]\+\_/@/g" | gzip > {output.RNAsort_R2}

    """

rule extract_umi:
    input:
        RNA_R1 = rules.sortfastq.output.RNAsort_R1,
        RNA_R2 = rules.sortfastq.output.RNAsort_R2,
    output:
        umi1="processed/umi/umi.{sample}.{barcode}.rna.R1.fq.gz",
        umi2="processed/umi/umi.{sample}.{barcode}.rna.R2.fq.gz",
        unzip_umi1="processed/RNA_all/umibycell.{sample}.{barcode}.rna.R1.fq"
    resources:
        nodes = 1
    params:
        pattern=r"NNNNNNNN",
    conda:"envs/main_env.yaml"
    shell: """
        umi_tools extract -p {params.pattern} -I {input.RNA_R2} -S {output.umi2} --read2-in={input.RNA_R1} --read2-out={output.umi1}

        gunzip --force -c {output.umi1} > processed/umi/umi.{wildcards.sample}.{wildcards.barcode}.rna.R1.fq

        sed 's/_/_{wildcards.sample}.{wildcards.barcode}_/' processed/umi/umi.{wildcards.sample}.{wildcards.barcode}.rna.R1.fq > {output.unzip_umi1}

        """

#output UMIs matrix.

rule RNAmerge:
    input:
        expand("processed/RNA_all/umibycell.{sample}.{barcode}.rna.R1.fq",sample=SAMPLES,barcode=BARCODES),
    output:
        rnaAll = temp("processed/RNA_all/rnaAll.fq"),
    conda:"envs/main_env.yaml",
    shell:"""
        cat {input} > {output.rnaAll}
    """

rule RNAclean:
    input:
        rnaAll = rules.RNAmerge.output.rnaAll,
    output:
        rnaAllclean="processed/RNA_all/rnaAll.clean.fq"
    threads: config["resources"]["cutadapt_cpu_threads"]
    resources:
        nodes = config["resources"]["cutadapt_cpu_threads"]
    params:
        adapter=r"AAAAAAAAAAAAAANNNNNNNNCATTGCGCAATACTACCTCAACCCTGTCTCTTATA",
    conda:"envs/main_env.yaml",
    shell:"""
        # cut nextera transposes and clean empty line.
        cutadapt -a {params.adapter} {input.rnaAll} -j {threads} --minimum-length 1 > {output.rnaAllclean}
    """

rule star_mapping:
    input:
        fastqIn = rules.RNAclean.output.rnaAllclean,
        starIndex = config["refs"][config["ref_genome"]]["star_index"]
    output:
        bamOut = "processed/RNA_all/starOut/star.Aligned.sortedByCoord.out.bam"
    threads: config["resources"]["star_cpu_threads"]
    resources:
        nodes = config["resources"]["star_cpu_threads"]
    conda:"envs/main_env.yaml"
    shell:"""
        STAR --runThreadN {threads} \
        --genomeDir {input.starIndex} \
        --readFilesIn {input.fastqIn} \
        --outFileNamePrefix processed/RNA_all/starOut/star. \
        --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes NH HI NM MD
        """

rule makeFolder:
    input:
        bam = rules.star_mapping.output.bamOut,
    output:
        totalBamOut = "processed/RNA_all/bamToCount/star.Aligned.sortedByCoord.out.total.bam",
    conda:"envs/main_env.yaml",
    shell:"""
        mkdir -p processed/RNA_all/bamToCount
        cp {input.bam} {output.totalBamOut}
    """

rule RNA_snp_split:
    input:
        RNAsnp = config["refs"][config["ref_genome"]]["RNAsnp"],
        bam = rules.makeFolder.output.totalBamOut,
        bamToProcess = rules.star_mapping.output.bamOut,
    output:
        dirOut = directory("processed/RNA_all/snpSplitRes"),
        genome1BamOut = "processed/RNA_all/snpSplitRes/star.Aligned.sortedByCoord.out.genome1.bam",
        genome2BamOut = "processed/RNA_all/snpSplitRes/star.Aligned.sortedByCoord.out.genome2.bam",
        genome1Copied = "processed/RNA_all/bamToCount/star.Aligned.sortedByCoord.out.genome1.bam",
        genome2Copied = "processed/RNA_all/bamToCount/star.Aligned.sortedByCoord.out.genome2.bam",
    conda:"envs/main_env.yaml"
    shell:"""
        
        mkdir -p processed/RNA_all/snpSplitRes
        SNPsplit --snp_file {input.RNAsnp} {input.bamToProcess} -o {output.dirOut}

        cp {output.genome1BamOut} {output.genome1Copied}
        cp {output.genome2BamOut} {output.genome2Copied}

        """

rule count:
    input:
        bamIn = "processed/RNA_all/bamToCount/star.Aligned.sortedByCoord.out.{genome}.bam",
        geneAnnotations = config["refs"][config["ref_genome"]]["annotations"],
    output:
        CountMatrix = "result/RNA_Res/counts.{type}.{genome}.tsv",
    params:
        countParams = r"--per-gene --per-cell --gene-tag=XT --wide-format-cell-counts --assigned-status-tag=XS",
    threads: config["resources"]["count_cpu_threads"]
    resources:
        nodes = config["resources"]["count_cpu_threads"]
    conda:"envs/main_env.yaml"
    shell:"""

        mkdir -p ./processed/RNA_all/feature_{wildcards.type}_{wildcards.genome}/

        featureCounts -a {input.geneAnnotations} -o ./processed/RNA_all/feature_{wildcards.type}_{wildcards.genome}/gene_assigned -R BAM {input.bamIn} -T {threads} -Q 30 -t {wildcards.type} -g gene_name -O -s 1

        #OUTPUT countMatix by gene
        samtools sort ./processed/RNA_all/feature_{wildcards.type}_{wildcards.genome}/star.Aligned.sortedByCoord.out.{wildcards.genome}.bam.featureCounts.bam -o ./processed/RNA_all/feature_{wildcards.type}_{wildcards.genome}/samsort.bam
        samtools index ./processed/RNA_all/feature_{wildcards.type}_{wildcards.genome}/samsort.bam
        umi_tools count {params.countParams} -I ./processed/RNA_all/feature_{wildcards.type}_{wildcards.genome}/samsort.bam -S ./result/RNA_Res/counts.{wildcards.type}.{wildcards.genome}.tsv

        """

rule convertCountFormat:
    input:
        countMatrix = rules.count.output.CountMatrix,
    output:
        convertedCountMatrix = "result/RNA_Res/counts.{type}.{genome}.format.tsv",
    conda:"envs/main_env.yaml"
    shell:"""
        
        Rscript R2D2/fraction.R {input.countMatrix} {output.convertedCountMatrix}
        
    """