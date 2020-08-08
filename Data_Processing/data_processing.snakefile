import os

configfile: "data_processing_config.json"

rule all:
    input:
        expand("{chrm_n}_missingness.smiss", chrm_n=config["autosomes"])

rule compute_individual_missingness_plink:
    input:
    output:
        "{chrm_n}_missingness.smiss"
    params:
        input_basename = "/data/CEM/wilsonlab/from_collaborators/princeton_kenya/plink/high_low_coverage_combined/impute1panel.{chrm_n}.merged_LDfilt_plink",
        output_basename = "{chrm_n}_missingness"
    shell:
        """
        plink2 --bfile {params.input_basename} --missing --out {params.output_basename}
        """
