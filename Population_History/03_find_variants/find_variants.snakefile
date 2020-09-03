import os

# Overall description: from the all sites VCF, filter to obtain variants that are high quality and in putatively neutral regions

configfile: "find_variants_config.json"

rule all:
    input:
        expand("biallelic_snps/chr{chrm_n}.all_high_cov.emit_all.biallelic_snps.vcf.gz", chrm_n=config["autosomes"]), #select_biallelic_snps
        expand("count_num_variants/biallelic_snps/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #count_number_of_biallelic_snps
        expand("vqsr/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.vcf.gz", chrm_n=config["autosomes"]), #vqsr
        expand("count_num_variants/post_vqsr/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #count_number_of_snps_post_vqsr
        expand("filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.vcf.gz.tbi", chrm_n=config["autosomes"]), #hwe, bgzip, tabix
        expand("count_num_variants/post_hwe_filter/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #count_number_of_snps_post_hwe_filter
        expand("count_num_variants/post_pilot_masks/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #pilot masks filter and count the number of variants
        expand("count_num_variants/putatively_neutral/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #neutral regions filter and count the number of variants
        expand("count_num_variants/nre_neutral/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"]), #nre neutral regions filter and count the number of variants
        "filtered_variants/autosomes.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.vcf.gz",
        expand("filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11.recode.vcf.gz", chrm_n=config["autosomes"]), #remove A1 and A11
        expand("filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11.recode.vcf.gz.tbi", chrm_n=config["autosomes"]), #remove A1 and A11
        expand("count_num_variants/post_remove_individuals/chr{chrm_n}_num_variants.txt", chrm_n=config["autosomes"])

rule select_biallelic_snps:
    input:
        ref = config["ref"],
        vcf = os.path.join(config["all_sites_vcf_dir"], "chr{chrm_n}.all_high_cov.emit_all.vcf.gz")
    output:
        "biallelic_snps/chr{chrm_n}.all_high_cov.emit_all.biallelic_snps.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--restrict-alleles-to BIALLELIC """
        """--select-type-to-include SNP """

rule count_number_of_biallelic_snps:
    input:
        "biallelic_snps/chr{chrm_n}.all_high_cov.emit_all.biallelic_snps.vcf.gz"
    output:
        "count_num_variants/biallelic_snps/chr{chrm_n}_num_variants.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """

rule gatk_variantrecalibrator:
    input:
        ref = config["ref"],
        vcf = "biallelic_snps/chr{chrm_n}.all_high_cov.emit_all.biallelic_snps.vcf.gz",
        hapmap = "/data/CEM/shared/public_data/validated_variant_resources/hapmap_3.3.hg38.vcf.gz",
        omni = "/data/CEM/shared/public_data/validated_variant_resources/1000G_omni2.5.hg38.vcf.gz",
        thousandG = "/data/CEM/shared/public_data/validated_variant_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    output:
        recal = "vqsr/chr{chrm_n}/chr{chrm_n}_output.recal",
        tranches = "vqsr/chr{chrm_n}/chr{chrm_n}_output.tranches"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} --java-options "-Xmx16g" VariantRecalibrator """
        """-R {input.ref} -V {input.vcf}  """
        """--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} """
        """--resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} """
        """--resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} """
        """-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff """
        """-mode SNP """
        """-O {output.recal} """
        """--tranches-file {output.tranches} """
        """--max-gaussians 4 """

rule gatk_applyvqsr:
    input:
        ref = config["ref"],
        vcf = "biallelic_snps/chr{chrm_n}.all_high_cov.emit_all.biallelic_snps.vcf.gz",
        tranches = "vqsr/chr{chrm_n}/chr{chrm_n}_output.tranches",
        recal = "vqsr/chr{chrm_n}/chr{chrm_n}_output.recal"
    output:
        "vqsr/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} --java-options "-Xmx16g" ApplyVQSR """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--truth-sensitivity-filter-level 90.0 """
        """--tranches-file {input.tranches} """
        """--recal-file {input.recal} """
        """-mode SNP """

rule gatk_selectvariants:
    input:
        ref = config["ref"],
        vcf = "vqsr/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.vcf.gz"
    output:
        "vqsr/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} --java-options "-Xmx16g" SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """--exclude-filtered """
        """-O {output} """

rule count_number_of_snps_post_vqsr:
    input:
        "vqsr/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.vcf.gz"
    output:
        "count_num_variants/post_vqsr/chr{chrm_n}_num_variants.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """

rule hwe_filter:
    input:
        "vqsr/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.vcf.gz"
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.vcf"
    params:
        basename = "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe"
    shell:
        """
        vcftools --gzvcf {input} --hwe 10e-6 --out {params.basename} --recode
        """

rule bgzip_vcfs_post_hwe_filter:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.vcf"
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.vcf.gz"
    shell:
        """
        bgzip -c {input} > {output}
        """

rule index_vcfs_post_hwe_filter:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.vcf.gz"
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input}
        """

rule count_number_of_snps_post_hwe_filter:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.vcf.gz"
    output:
        "count_num_variants/post_hwe_filter/chr{chrm_n}_num_variants.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """

rule pilot_masks_filter:
    input:
        ref = config["ref"],
        vcf = "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.vcf.gz",
        interval = "/data/CEM/wilsonlab/from_collaborators/princeton_kenya/reference/1000_genomes_pilot_masks/20160622.allChr.pilot_mask.chr{chrm_n}.bed"
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} VariantFiltration """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-L {input.interval} """
        """-O {output} """

rule count_number_of_snps_post_masks_filter:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.vcf.gz"
    output:
        "count_num_variants/post_pilot_masks/chr{chrm_n}_num_variants.txt"
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
        ref = config["ref"],
        vcf = "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.vcf.gz",
        interval = os.path.join(config["neutral_regions_dir"], "chr{chrm_n}/chr{chrm_n}_putatively_neutral_regions_112719.bed")
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.neutral.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} VariantFiltration """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-L {input.interval} """
        """-O {output} """


rule count_number_of_snps_post_neutral_regions_filter:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.neutral.vcf.gz"
    output:
        "count_num_variants/putatively_neutral/chr{chrm_n}_num_variants.txt"
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
        ref = config["ref"],
        vcf = "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.vcf.gz",
        interval = os.path.join(config["nre_neutral_regions_dir"], "neutral_regions_nre_autosomes_liftoverhg38_chr{chrm_n}.bed")
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.vcf.gz"
    params:
        gatk = config["gatk4_path"]
    shell:
        """{params.gatk} VariantFiltration """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-L {input.interval} """
        """-O {output} """


rule count_number_of_snps_post_nre_neutral_regions_filter:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.vcf.gz"
    output:
        "count_num_variants/nre_neutral/chr{chrm_n}_num_variants.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """

# ---------------
# Concat variants
# ---------------
rule concat_vcfs_autosomes:
    input:
        vcfs = expand(
			"filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.vcf.gz", chrm_n=config["autosomes"])
    output:
        "filtered_variants/autosomes.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.vcf.gz"
    run:
        variant_files = []
        for i in input.vcfs:
        	variant_files.append("-I " + i)
        variant_files = " ".join(variant_files)
        shell(
        	"""picard MergeVcfs {variant_files} -O {output} """)

# -----------------
# Remove A7 and A11
# -----------------
rule remove_individuals:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.vcf.gz"
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11.recode.vcf"
    params:
        basename = "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11"
    shell:
        """
        vcftools --gzvcf {input} --remove-indv A7 --remove-indv A11 --out {params.basename} --recode
        """

rule bgzip_vcfs_post_remove_individuals:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11.recode.vcf"
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11.recode.vcf.gz"
    shell:
        """
        bgzip -c {input} > {output}
        """

rule index_vcfs_post_remove_individuals:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11.recode.vcf.gz"
    output:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11.recode.vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input}
        """

rule count_number_of_snps_post_remove_individuals:
    input:
        "filtered_variants/chr{chrm_n}/chr{chrm_n}.gatk.called.vqsr.select.variants.hwe.recode.pilot.masks.nre.neutral.rmA7A11.recode.vcf.gz"
    output:
        "count_num_variants/post_remove_individuals/chr{chrm_n}_num_variants.txt"
    params:
        script = config["calc_num_sites_in_vcf_script"],
        id = "{chrm_n}"
    shell:
        """
        python {params.script} --vcf {input} --id {params.id} > {output}
        """
