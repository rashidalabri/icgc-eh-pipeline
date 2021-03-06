import pandas as pd
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

AWS_SAMPLES = pd.read_table("data/manifest.aws.tsv").set_index("object_id", drop=False)
PDC_SAMPLES = pd.read_table("data/manifest.pdc.tsv").set_index("guid", drop=False)

rule aws:
    input:
        S3.remote(expand("icgc-eh-bucket/results/aws/{object_id}.json", object_id=AWS_SAMPLES.index)),
        S3.remote(expand("icgc-eh-bucket/results/aws/{object_id}.vcf", object_id=AWS_SAMPLES.index)),
        S3.remote(expand("icgc-eh-bucket/results/aws/{object_id}_realigned.bam", object_id=AWS_SAMPLES.index)),

rule pdc:
    input:
        S3.remote(expand("icgc-eh-bucket/results/pdc/{guid}.json", guid=PDC_SAMPLES.index)),
        S3.remote(expand("icgc-eh-bucket/results/pdc/{guid}.vcf", guid=PDC_SAMPLES.index)),
        S3.remote(expand("icgc-eh-bucket/results/pdc/{guid}_realigned.bam", guid=PDC_SAMPLES.index))

rule test_pdc:
    input:
        S3.remote(f"icgc-eh-bucket/results/pdc/39dfdb0b-4a7d-4b35-ba0b-abd85b39f846.json"),
        S3.remote(f"icgc-eh-bucket/results/pdc/39dfdb0b-4a7d-4b35-ba0b-abd85b39f846.vcf"),
        S3.remote(f"icgc-eh-bucket/results/pdc/39dfdb0b-4a7d-4b35-ba0b-abd85b39f846_realigned.bam")

rule test_aws:
    input:
        S3.remote(f"icgc-eh-bucket/results/aws/ccf30f0f-0ec1-5aa2-8eb3-14cb1f4f58cd.json"),
        S3.remote(f"icgc-eh-bucket/results/aws/ccf30f0f-0ec1-5aa2-8eb3-14cb1f4f58cd.vcf"),
        S3.remote(f"icgc-eh-bucket/results/aws/ccf30f0f-0ec1-5aa2-8eb3-14cb1f4f58cd_realigned.bam")


rule genotype_aws:
    input:
        fa="data/ref/hs37d5.fa.gz",
        fai="data/ref/hs37d5.fa.gz.fai",
        var="data/variant_catalog_hg19.json"
    params:
        sex=lambda wildcards: AWS_SAMPLES.loc[wildcards["object_id"], "sex"],
        file_name=lambda wildcards: AWS_SAMPLES.loc[wildcards["object_id"], "file_name"],
        prefix=lambda wildcards, output: output.json[:-5]
    output:
        json=S3.remote("icgc-eh-bucket/results/aws/{object_id}.json"),
        vcf=S3.remote("icgc-eh-bucket/results/aws/{object_id}.vcf"),
        realigned_bam=S3.remote("icgc-eh-bucket/results/aws/{object_id}_realigned.bam")
    envmodules:
        "samtools",
        "oracle-java"
    log:
        score=S3.remote("icgc-eh-bucket/logs/aws/{object_id}-eh.log"),
        eh=S3.remote("icgc-eh-bucket/logs/aws/{object_id}-score-client.log")
    resources:
        mem=25,
        time=24
    threads: 16
    shell:
        "/home/rashid/score-client-5.8.1/bin/score-client "
        "--profile collab download --output-dir {resources.tmpdir} "
        "--object-id {wildcards.object_id} > {log.score} && "
        "ExpansionHunter --reads {resources.tmpdir}/{params.file_name} "
        "--reference {input.fa} "
        "--variant-catalog {input.var} "
        "--output-prefix {params.prefix} "
        "--sex {params.sex} "
        "--analysis-mode streaming "
        "--threads {threads} > {log.eh}"

rule genotype_pdc:
    input:
        fa="data/ref/hs37d5.fa.gz",
        fai="data/ref/hs37d5.fa.gz.fai",
        var="data/variant_catalog_hg19.json"
    params:
        sex=lambda wildcards: PDC_SAMPLES.loc[wildcards["guid"], "sex"],
        file_name=lambda wildcards: PDC_SAMPLES.loc[wildcards["guid"], "file_name"],
        prefix=lambda wildcards, output: output.json[:-5],
        out_dir=lambda wildcards, output: "/".join(output.json.split("/")[:-1])
    output:
        json=S3.remote("icgc-eh-bucket/results/pdc/{guid}.json"),
        vcf=S3.remote("icgc-eh-bucket/results/pdc/{guid}.vcf"),
        realigned_bam=S3.remote("icgc-eh-bucket/results/pdc/{guid}_realigned.bam")
    envmodules:
        "samtools"
    log:
        S3.remote("icgc-eh-bucket/logs/pdc/{guid}.log")
    resources:
        mem=25,
        time=24
    threads: 16
    shell:
        "gen3-client download-single --profile=icgc --guid {wildcards.guid} "
        "--download-path {tmpdir} --no-prompt --skip-completed && "
        "ExpansionHunter --reads {tmpdir}/{params.file_name} "
        "--reference {input.fa} "
        "--variant-catalog {input.var} "
        "--output-prefix {params.prefix} "
        "--sex {params.sex} "
        "--analysis-mode streaming "
        "--threads {threads} "
        "> {log} && "
        "rm -f {tmpdir}/{params.file_name}"