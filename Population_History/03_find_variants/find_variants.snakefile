import os

configfile: "find_variants_config.json"

rule all:
    input: #Amanda
        "pass_variants/high_cov_chr22.SNP1.vqsr.recode_select.pass.biallelic.snps.vcf.gz"
    input: # from all sites
        "pass_variants/chr22.all_high_cov.emit_all_select.biallelic.variants_hard.filter_select.pass_snp.vcf.gz"
    input:
        expand(os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.vcf.gz"), chrm_n=config["autosomes"]),
        expand(os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.vcf.gz.tbi"), chrm_n=config["autosomes"])


rule split_by_chromosomes:
    input:
        config["variant_vcf_from_amanda"]
    output:
        os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.vcf")
    params:
        chrm_n = "chr{chrm_n}",
        out_basename = os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr")
    shell:
        """
        vcftools --gzvcf {input} --chr {params.chrm_n} --recode --recode-INFO-all --out {params.out_basename}
        """

rule bgzip_vcfs:
    input:
        os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.vcf")
    output:
        os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.vcf.gz")
    shell:
        """
        bgzip -c {input} > {output}
        """

rule index_vcfs:
    input:
        os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.vcf.gz")
    output:
        os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.vcf.gz.tbi")
    shell:
        """
        tabix -p vcf {input}
        """

rule select_pass_variant:
    input:
        os.path.join(config["variant_vcf_dir"], "high_cov_chr22.SNP1.vqsr.recode.vcf.gz")
    output:
        "pass_variants/high_cov_chr22.SNP1.vqsr.recode_select.pass.vcf"
    shell:
        """
        python /home/tphung3/softwares/tanya_repos/vcfhelper/obtain_pass_variants.py --input_vcf {input} --output_vcf {output}
        """

rule select_biallelic_snps:
    input:
        ref = config["ref"],
        vcf = "pass_variants/high_cov_chr22.SNP1.vqsr.recode_select.pass.vcf.gz"
    output:
        "pass_variants/high_cov_chr22.SNP1.vqsr.recode_select.pass.biallelic.snps.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """--select-type-to-include SNP """
        """--restrict-alleles-to BIALLELIC """
        """-O {output} """

# rule select_pass_variants:
#     input:
#         ref = config["ref"],
#         vcf = os.path.join(config["variant_vcf_dir"], "high_cov_chr{chrm_n}.SNP1.vqsr.recode.vcf.gz")
#     output:
#         "pass_variants/high_cov_chr{chrm_n}.SNP1.vqsr.recode.biallelic.pass.vcf.gz"
#     params:
#         gatk = config["gatk4_path"]
#     shell:
#         """{params.gatk} SelectVariants """
#         """-R {input.ref} """
#         """-V {input.vcf} """
#         """-O {output} """
#         """--exclude-filtered """
#         """--restrict-alleles-to BIALLELIC """
#
# rule filter_by_DP:
#     input:
#         ref = config["ref"],
#         vcf = "pass_variants/high_cov_chr{chrm_n}.SNP1.vqsr.recode.biallelic.pass.vcf.gz"
#     output:
#         "pass_variants/high_cov_chr{chrm_n}.SNP1.vqsr.recode.biallelic.pass.filtered_DP.vcf.gz"
#     params:
#         gatk = config["gatk4_path"],
#         DP = 130.0
#     shell:
#         """{params.gatk} SelectVariants """
#         """-R {input.ref} """
#         """-V {input.vcf} """
#         """-O {output} """
#         """-select "DP >= {params.DP}" """
#
# rule select_variants_in_neutral_regions:
#     input:
#         ref = config["ref"],
#         vcf = "pass_variants/high_cov_chr{chrm_n}.SNP1.vqsr.recode.biallelic.pass.filtered_DP.vcf.gz",
#         neutral = os.path.join(config["neutral_regions_dir"], "chr{chrm_n}/chr{chrm_n}_putatively_neutral_regions_112719.bed")
#     output:
#         "pass_variants/high_cov_chr{chrm_n}.SNP1.vqsr.recode.biallelic.pass.filtered_DP_neutral.vcf.gz"
#     params:
#         gatk = config["gatk4_path"]
#     shell:
#         """{params.gatk} SelectVariants """
#         """-R {input.ref} """
#         """-V {input.vcf} """
#         """-O {output} """
#         """--intervals {input.neutral} """

# ---------------
# Troubleshooting
# ---------------

# Hard filter from all sites output
rule select_variants:
    input:
        ref = config["ref"],
        vcf = "/data/CEM/wilsonlab/from_collaborators/princeton_kenya/vcf/high_coverage_all_sites/chr22.all_high_cov.emit_all.vcf.gz"
    output:
        "pass_variants/chr22.all_high_cov.emit_all_select.biallelic.variants.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--restrict-alleles-to BIALLELIC """

rule hard_filter_chr22:
    input:
        ref = config["ref"],
        vcf = "pass_variants/chr22.all_high_cov.emit_all_select.biallelic.variants.vcf.gz"
    output:
        "pass_variants/chr22.all_high_cov.emit_all_select.biallelic.variants_hard.filter.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} VariantFiltration """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-filter "QD < 2.0" --filter-name "QD2" """
        """-filter "QUAL < 30.0" --filter-name "QUAL30" """
        """-filter "SOR > 3.0" --filter-name "SOR3" """
        """-filter "FS > 60.0" --filter-name "FS60" """
        """-filter "MQ < 40.0" --filter-name "MQ40" """
        """-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" """
        """-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" """
        """-O {output} """

rule select_pass_variants_post_hard_filter_chr22:
    input:
        ref = config["ref"],
        vcf = "pass_variants/chr22.all_high_cov.emit_all_select.biallelic.variants_hard.filter.vcf.gz"
    output:
        "pass_variants/chr22.all_high_cov.emit_all_select.biallelic.variants_hard.filter_select.pass.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--exclude-filtered """

rule select_pass_variants_post_hard_filter_snp_chr22:
    input:
        ref = config["ref"],
        vcf = "pass_variants/chr22.all_high_cov.emit_all_select.biallelic.variants_hard.filter_select.pass.vcf.gz"
    output:
        "pass_variants/chr22.all_high_cov.emit_all_select.biallelic.variants_hard.filter_select.pass_snp.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """--select-type-to-include SNP """
        """-O {output} """
