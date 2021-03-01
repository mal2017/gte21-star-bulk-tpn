rule concat_fqs:
    input:
        lambda wc: [determine_resource(x) for x in config["samples"][wc.s]["fastq"][wc.end]]
    output:
        temp("results/fastq/{s}_{end}.fq.gz")
    shell:
        "cat {input} > {output}"

def get_fqs_for_trim(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return ["results/fastq/{s}_r1.fq.gz"]
    else:
        return ["results/fastq/{s}_r1.fq.gz", "results/fastq/{s}_r2.fq.gz"]

def get_fqs_for_aln(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return ["results/fastq/{s}_r1.trimmed.fq.gz"]
    else:
        return ["results/fastq/{s}_r1.trimmed.fq.gz", "results/fastq/{s}_r2.trimmed.fq.gz"]


def get_proper_ended_fastp_call(x):
    fqs = get_fqs_for_trim(x)
    if len(fqs) == 1:
        return "--in1 {r1}".format(r1=fqs[0].format(s=x))
    else:
        return "--in1 {r1} --in2 {r2}".format(r1=fqs[0].format(s=x), r2=fqs[1].format(s=x))

def get_proper_ended_fastp_out(x):
    fqs = get_fqs_for_aln(x)
    if len(fqs) == 1:
        return "--out1 {r1}".format(r1=fqs[0].format(s=x))
    else:
        return "--out1 {r1} --out2 {r2}".format(r1=fqs[0].format(s=x), r2=fqs[1].format(s=x))

rule trim_se:
    input:
        fq = lambda wc: get_fqs_for_trim(wc.s)
    output:
        r1 = "results/fastq/{s}_r1.trimmed.fq.gz",
        html = "results/fastq/{s}_fastp.html",
        json = "results/fastq/{s}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.s),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.s)
    conda:
        "../envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.s}_fastp"

rule trim_pe:
    input:
        fq = lambda wc: get_fqs_for_trim(wc.s)
    output:
        r1 = "results/fastq/{s}_r1.trimmed.fq.gz",
        r2 = "results/fastq/{s}_r2.trimmed.fq.gz",
        html = "results/fastq/{s}_fastp.html",
        json = "results/fastq/{s}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.s),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.s)
    conda:
        "../envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.s}_fastp"


rule generate_star_idx:
    """
    https://kb.10xgenomics.com/hc/en-us/articles/115004415786-What-parameters-are-used-for-STAR-alignment-
    https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/align_reads/__init__.py
    https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/reference.py
    """
    input:
        genome = custom_genome('results/custom-genome/combined.fasta'),
        gtf = custom_genome('results/custom-genome/combined.fixed.gtf')
    output:
        "results/index/chrLength.txt",
        "results/index/chrNameLength.txt",
        "results/index/chrName.txt",
        "results/index/chrStart.txt",
        "results/index/exonGeTrInfo.tab",
        "results/index/exonInfo.tab",
        "results/index/Genome",
        "results/index/geneInfo.tab",
        "results/index/genomeParameters.txt",
        "results/index/SA",
        "results/index/SAindex",
        "results/index/sjdbList.fromGTF.out.tab",
        "results/index/sjdbInfo.txt",
        "results/index/sjdbList.out.tab",
        "results/index/transcriptInfo.tab",
    threads:
        12
    conda:
        "../envs/star.yaml"
    shadow:
        "minimal"
    resources:
        time=60,
        mem=20000,
        cpus=12
    #singularity:
    #    "docker://quay.io/biocontainers/star:2.7.3a--0"
    shell:
        "STAR --runMode genomeGenerate "
        "--genomeDir results/index/ --runThreadN {threads} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} "
        "--genomeSAindexNbases 11"

rule star_align:
    input:
        rules.generate_star_idx.output,
        reads = lambda wc: get_fqs_for_aln(wc.s),
        gtf = custom_genome('results/custom-genome/combined.fixed.gtf')
    output:
        directory("results/star/{s}/")
    threads:
        12
    conda:
        "../envs/star.yaml"
    resources:
        time=240,
        mem=20000,
        cpus=12
    #shadow:
    #    "minimal"
    #singularity:
    #    "docker://quay.io/biocontainers/star:2.7.3a--0"
    shell:
        "STAR "
        "--genomeDir results/index --runThreadN {threads} "
        "--readFilesIn {input.reads} "
        "--outSAMmultNmax 1 "
        "--outSAMtype BAM Unsorted "
        "--quantMode GeneCounts "
        "--sjdbGTFfile {input.gtf} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix {output}/"
