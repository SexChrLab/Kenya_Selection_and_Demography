# Environment: turkana_demography
import os

configfile: "find_callable_regions_config.json"

rule all:
    input:
        expand("DP/chr{chrm_n}_DP.txt", chrm_n=config["autosomes"]),
        expand("DP/chr{chrm_n}_DP_sum.txt", chrm_n=config["autosomes"]),
        "DP/autosomes_DP_sum.txt",
        expand(os.path.join(config["all_sites_vcf_dir"], "chr{chrm_n}.all_high_cov.emit_all.vcf.gz.tbi"), chrm_n=config["autosomes"]),
        expand("callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.vcf.gz", chrm_n=config["autosomes"]),
        expand("callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf.gz.tbi", chrm_n=config["autosomes"]), #filter by AN, bgzip and tabix
        expand("count_num_sites/all_sites/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #count number of all_sites
        expand("count_num_sites/post_DP_filter/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #count number of sites post DP filter
        expand("count_num_sites/post_AN_filter/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #count number of sites post AN filter
        expand("count_num_sites/putatively_neutral/chr{chrm_n}_num_sites.txt", chrm_n=config["autosomes"]), #filter for putatively neutral and count number of sites after neutral regions filtering
        expand("count_num_sites/nre_neutral/chr{chrm_n}_num_sites.txt", chrm_n=config["autosomes"]) #filter for nre neutral and count the number of sites after neutral regions filtering

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

rule index_vcfs:
    input:
        os.path.join(config["all_sites_vcf_dir"], "chr{chrm_n}.all_high_cov.emit_all.vcf.gz")
    output:
        os.path.join(config["all_sites_vcf_dir"], "chr{chrm_n}.all_high_cov.emit_all.vcf.gz.tbi")
    shell:
        """
        tabix -p vcf {input}
        """

rule DP_filter:
    input:
        ref = config["ref"],
        vcf = os.path.join(config["all_sites_vcf_dir"], "chr{chrm_n}.all_high_cov.emit_all.vcf.gz")
    output:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.vcf.gz"
    params:
        gatk = config["gatk4_path"],
        DP = 130.0
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """-select "DP >= {params.DP}" """

rule AN_filter:
    input:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.vcf.gz"
    output:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf"
    params:
        basename = "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN"
    shell:
        """
        vcftools --gzvcf {input} --max-missing-count 55 --out {params.basename} --recode
        """

rule bgzip_vcfs_post_filter:
    input:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf"
    output:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf.gz"
    shell:
        """
        bgzip -c {input} > {output}
        """

rule index_vcfs_post_filter:
    input:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf.gz"
    output:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input}
        """

rule count_number_of_all_sites:
    input:
        os.path.join(config["all_sites_vcf_dir"], "chr{chrm_n}.all_high_cov.emit_all.vcf.gz")
    output:
        "count_num_sites/all_sites/chr{chrm_n}_num_variants.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """

rule count_number_of_sites_post_DP_filter:
    input:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.vcf.gz"
    output:
        "count_num_sites/post_DP_filter/chr{chrm_n}_num_variants.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """

rule count_number_of_sites_post_AN_filter:
    input:
        "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf.gz"
    output:
        "count_num_sites/post_AN_filter/chr{chrm_n}_num_variants.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """

# ----------------------------------
# Tanya's putatively neutral regions
# ----------------------------------
rule neutral_regions_filter:
    input:
        vcf = "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf.gz",
        interval = os.path.join(config["neutral_regions_dir"], "chr{chrm_n}/chr{chrm_n}_putatively_neutral_regions_112719.bed")
    output:
        "neutral_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.neutral.vcf"
    shell:
        """
        bcftools view --regions-file {input.interval} {input.vcf} -output-type z --output-file {output}
        """


rule count_number_of_sites_post_neutral_regions_filter:
    input:
        "neutral_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.neutral.vcf"
    output:
        "count_num_sites/putatively_neutral/chr{chrm_n}_num_sites.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """

# -------------------
# NRE neutral regions
# -------------------
rule nre_neutral_regions_filter:
    input:
        vcf = "callable_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.vcf.gz",
        interval = os.path.join(config["nre_neutral_regions_dir"], "neutral_regions_nre_autosomes_liftoverhg38_chr{chrm_n}.bed")
    output:
        "neutral_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.nre.neutral.vcf"
    shell:
        """
        bcftools view --regions-file {input.interval} {input.vcf} -output-type z --output-file {output}
        """


rule count_number_of_sites_post_nre_neutral_regions_filter:
    input:
        "neutral_regions/chr{chrm_n}.all_high_cov.emit_all.filtered_DP.AN.recode.nre.neutral.vcf"
    output:
        "count_num_sites/nre_neutral/chr{chrm_n}_num_sites.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """
