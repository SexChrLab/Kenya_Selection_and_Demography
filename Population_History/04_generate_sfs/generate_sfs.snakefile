import os

configfile: "generate_sfs_config.json"

rule all:
    input:
        expand("sfs/chr{chrm_n}_sfs.txt", chrm_n=config["autosomes"]),
        expand("sfs/chr{chrm_n}_raw_sfs.txt", chrm_n=config["autosomes"])

rule generate_sfs_autosomes:
    input:
        os.path.join(config["vcf_directory"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.biallelic.pass.filtered_DP_neutral.vcf.gz")
    output:
        "sfs/chr{chrm_n}_sfs.txt"
    params:
        script = config["generate_sfs_script"]
    shell:
        """
        python {params.script} --vcf_file {input} --sfs_all --sfs_all_out {output} --ploidy diploid
        """

# Troubleshooting
rule generate_sfs_autosomes_raw:
    input:
        os.path.join(config["vcf_directory"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.biallelic.pass.vcf.gz")
    output:
        "sfs/chr{chrm_n}_raw_sfs.txt"
    params:
        script = config["generate_sfs_script"]
    shell:
        """
        python {params.script} --vcf_file {input} --sfs_all --sfs_all_out {output} --ploidy diploid
        """
