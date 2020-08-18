## 01_find_neutral_regions
- The (putatively) neutral regions are defined as:
    1. Outside of genic regions, conserved regions, CpG islands, and repeats
    2. Bscores are greater or equal to 0.9
    3. #TODO: genetic distance?
- (Putatively) neutral regions were previously identified for a separate project. The bed files representing the (putatively) neutral regions are in: `/scratch/tphung3/Turkana_Demography/04_find_neutral_regions/results/`
- Briefly, below is how the (putatively) neutral regions were generated:

1. Find regions from UCSC tracks
    1. UCSC tracks:
        1. https://genome.ucsc.edu/cgi-bin/hgTables
    2. **Genic regions**: clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: Genes and Gene Predictions; track: GENCODE v32; table: knownGene; region: genome; output format: BED - browser extensible data; output file: GRCh38_gencodev32_genome_knownGene; file type returned: gzip compressed; Create one BED record per: Whole Gene
    3. **Conserved regions**: clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: Comparative Genomics; track: Conservation; table: 100 Vert. El (phastConsElements100way); region: genome; output format: BED - browser extensible data; output file: GRCh38_phastConsElements100way; file type returned: gzip compressed
    4. **Repeats**: clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: Repeats; track: RepeatMasker; table: msk; region: genome; output format: BED - browser extensible data; output file: GRCh38_repeats; file type returned: gzip compressed
    5. **CpG islands**: clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: regulation; track: CpG islands; table: cpgIslandExt; region: genome; output format: BED - browser extensible data; output file: GRCh38_cpgislands; file type returned: gzip compressed

2. B values
    1. Obtaining B values:
        - The file `bscores.tsv.gz` was obtained from Adriana Arneson (UCLA) who downloaded it from the CADD annotation.
        - The format of this file is chr, start, end, b values
    2. Process the bscore file:
        - We want to convert the naming of chromosome from 1 to chr1.
        - Also, we just want to keep the sites where bscore is greater or equal to 0.9
        - Lift over the coordinates for bscores from hg19 to GRCh38

3. Create bed files representing putatively neutral regions using bedtools

## 02_find_callable_regions
- The callable regions are defined as:
    1. DP is between 50% and 150% of the mean
    2. AN is equal or greater than 110
        - Rationale: there are 110 individuals
        - For the first pass, I am using 110 as a threshold so this mean that there are at least half of the total alleles being genotyped
    2. Overlapping with 1000 genome masks (pilot)

- Snakefile: `find_callable_regions.snakefile`
    - config file: `find_callable_regions_config.json`

- Steps:
    1. Extract DP
        1. rule `extract_DP` (line 10)
    2. Compute sum DP and number of sites
        1. rule `compute_sum_DP` (line 23) to compute sum per chromosome
        2. rule `concat_sum_DP` (line 36) to concat across all chromosomes genome-wide
    3. Compute genome-wide mean DP and 50% of 150% of the mean
        1. `Rscript compute_DP_threshold.R`
        2. Genome-wide: mean DP is 259.4879.
        3. Threshold: filter out sites where DP is less than 50% of the mean. In other words, remove any sites where DP is less than 129.7439 (rounded to 130).
        4. **Notes/To be considered: a more strigent threshold for DP (having an upper threshold for DP)?**
    4. Filter sites based on DP threshold:
        1. rule `find_callable_regions` (line 59)


## 03_find_variants
- From the vcf file for variants that Amanda sent, do the following:
    1. Select "PASS" variants
    2. Select variants based on DP and AN thresholds
    3. Select variants in putatively neutral regions


## 04_generate_sfs
- Aggregate across all chromosomes:
```
python tabulate_foldedSFS.py --num_bin 110 --directory sfs/ --out_filename sfs/autosomes_sfs.txt
```

## 05_demographic_inference
1. Convert SFS to dadi file format
  ```
  grep -v af_bin ../05_generate_sfs/results/chrX/chrX.females.sfs.txt > ../05_generate_sfs/results/chrX/chrX.females.sfs.clean.txt

  python convert_sfs_to_dadi_format.py --num_bin 93 --folded_or_unfolded folded --population_name Turkana --sfs_filename ../05_generate_sfs/results/chrX/chrX.females.sfs.clean.txt --num_individuals 46 --out_filename chrX_females_sfs_02082020_dadi_format.sfs
  ```

  ```
  grep -v af_bin ../05_generate_sfs/results/chrX/chrX.males.sfs.txt > ../05_generate_sfs/results/chrX/chrX.males.sfs.clean.txt

  python convert_sfs_to_dadi_format.py --num_bin 73 --folded_or_unfolded folded --population_name Turkana --sfs_filename ../05_generate_sfs/results/chrX/chrX.males.sfs.clean.txt --num_individuals 36 --out_filename chrX_males_sfs_03052020_dadi_format.sfs
  ```

2. Run dadi
  1. chrX females
  ```
  for i in {1..50}; do python 1D.2Epoch.dadi.py --runNum ${i} --pop Turkana --mu 1.5e-8 --L 10526249 --sfs chrX_females_sfs_02082020_dadi_format.sfs --outdir 2Epoch_chrX_females_run_${i}; done;
  ```

  2. chrX males
  ```

  ```

  ```
  python merge_out.py --directory . --run_base_directory 2Epoch_chrX_males --file_basename Turkana.dadi.inference.1D.2Epoch.runNum --out_filename Turkana.dadi.inference.1D.2Epoch.chrX.males.50.runs.csv --out_filename_sorted Turkana.dadi.inference.1D.2Epoch.chrX.males.50.runs.sorted.csv
  ```

  ```
  python parse_dadi_expsfs.py --dadi_expsfs 2Epoch_run_3/Turkana.dadi.inference.1D.2Epoch.runNum.3.expSFS --num_individuals 18 --theta 9028.43395222795 --out_filename 2Epoch_chrX_males_run_3_sfs_normalized_by_theta.txt
  ```



## Notes:
- In the vcf files with all sites output, there are 110 individuals
- In the vcf files with variants only after vqsr, there are 108 individuals
