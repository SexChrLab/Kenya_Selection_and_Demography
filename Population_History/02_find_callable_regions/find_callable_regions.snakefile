# Environment: turkana_demography
import os

configfile: "find_callable_regions_config.json"

rule all:
    input:
        expand("DP/chr{chrm_n}_DP.txt", chrm_n=config["autosomes"]),
        expand("DP/chr{chrm_n}_DP_sum.txt", chrm_n=config["autosomes"]),
        "DP/autosomes_DP_sum.txt"

rule extract_DP:
    input:
        "/data/CEM/wilsonlab/from_collaborators/princeton_kenya/vcf/high_coverage_all_sites/chr{chrm_n}.all_high_cov.emit_all.vcf.gz"
    output:
        "DP/chr{chrm_n}_DP.txt"
    params:
        script = config["extract_stats_from_vcf_script"],
        stat = "DP"
    shell:
        """
        python {params.script} {params.stat} --vcf {input} --outfile {output}
        """

rule compute_sum_DP:
    input:
        "DP/chr{chrm_n}_DP.txt"
    output:
        "DP/chr{chrm_n}_DP_sum.txt"
    params:
        script = config["compute_sum_DP_script"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule concat_sum_DP:
    input:
        expand("DP/chr{chrm_n}_DP_sum.txt", chrm_n=config["autosomes"])
    output:
        "DP/autosomes_DP_sum.txt"
    shell:
        """
        cat {input} >> {output}
        """
