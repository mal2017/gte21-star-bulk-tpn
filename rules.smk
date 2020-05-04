rule concat_fqs:
    input:
        lambda wc: [determine_resource(x) for x in config["samples"][wc.s]["fastq"][wc.end]]
    output:
        temp("fastq/{s}_{end}.fq.gz")
    shell:
        "cat {input} > {output}"

def get_fqs_for_trim(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return ["fastq/{s}_r1.fq.gz"]
    else:
        return ["fastq/{s}_r1.fq.gz", "fastq/{s}_r2.fq.gz"]

def get_fqs_for_aln(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return ["fastq/{s}_r1.trimmed.fq.gz"]
    else:
        return ["fastq/{s}_r1.trimmed.fq.gz", "fastq/{s}_r2.trimmed.fq.gz"]


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
        r1 = "fastq/{s}_r1.trimmed.fq.gz",
        html = "fastq/{s}_fastp.html",
        json = "fastq/{s}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.s),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.s)
    conda:
        "envs/fastp.yaml"
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
        r1 = "fastq/{s}_r1.trimmed.fq.gz",
        r2 = "fastq/{s}_r2.trimmed.fq.gz",
        html = "fastq/{s}_fastp.html",
        json = "fastq/{s}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.s),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.s)
    conda:
        "envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.s}_fastp"

rule get_gtf:
    output:
        "ref/tx.gtf"
    params:
        fl=config.get("gtf",None)
    shell:
        "curl -J -L {params.fl} > {output[0]}"

rule generate_star_idx:
    """
    https://kb.10xgenomics.com/hc/en-us/articles/115004415786-What-parameters-are-used-for-STAR-alignment-
    https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/align_reads/__init__.py
    https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/reference.py
    """
    input:
        genome = determine_resource(config.get("fasta",None)),
        gtf = determine_resource(config.get('gtf', None))
    output:
        "index/chrLength.txt",
        "index/chrNameLength.txt",
        "index/chrName.txt",
        "index/chrStart.txt",
        "index/exonGeTrInfo.tab",
        "index/exonInfo.tab",
        "index/Genome",
        "index/geneInfo.tab",
        "index/genomeParameters.txt",
        "index/SA",
        "index/SAindex",
        "index/sjdbList.fromGTF.out.tab",
        "index/sjdbInfo.txt",
        "index/sjdbList.out.tab",
        "index/transcriptInfo.tab",
    threads:
        config.get("MAX_THREADS", 1)
    conda:
        "envs/star.yaml"
    shadow:
        "minimal"
    #singularity:
    #    "docker://quay.io/biocontainers/star:2.7.3a--0"
    shell:
        "STAR --runMode genomeGenerate "
        "--genomeDir index/ --runThreadN {threads} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} "
        "--genomeSAindexNbases 11"
