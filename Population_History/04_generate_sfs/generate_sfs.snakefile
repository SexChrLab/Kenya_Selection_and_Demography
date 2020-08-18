import os

configfile: "generate_sfs_config.json"

rule all:
    input:
        expand("sfs/chr{chrm_n}_sfs.txt", chrm_n=config["autosomes"])

rule generate_sfs_autosomes:
    input:
        os.path.join(config["vcf_directory"], "chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.neutral.vcf.gz")
    output:
        "sfs/chr{chrm_n}_sfs.txt"
    params:
        script = config["generate_sfs_script"]
    shell:
        """
        python {params.script} --vcf_file {input} --sfs_all --sfs_all_out {output} --ploidy diploid
        """
