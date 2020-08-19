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
### Without including monomorphic sites
1. Convert SFS to dadi file format
```
python ~/softwares/tanya_repos/dadi_tutorial/convert_sfs_to_dadi_format.py --num_bin 221 --folded_or_unfolded folded --population_name Turkana --sfs_filename ../04_generate_sfs/sfs/autosomes_sfs.txt --num_individuals 110 --out_filename autosomes_sfs_dadi_format.sfs --num_monomorphic 0 --include_monomorphic no
```

2. Run dadi
```
for i in {1..50}; do python ~/softwares/tanya_repos/dadi_tutorial/1D.2Epoch.dadi.py --runNum ${i} --pop Turkana --mu 1.5e-8 --L 215622954 --sfs autosomes_sfs_dadi_format.sfs --outdir 1D_2Epoch/1D_2Epoch_autosomes_run_${i}; done;
```

```
python ~/softwares/tanya_repos/dadi_tutorial/merge_out.py --directory . --run_base_directory 1D_2Epoch_autosomes --file_basename Turkana.dadi.inference.1D.2Epoch.runNum --out_filename Turkana.dadi.inference.1D.2Epoch.autosomes.50.runs.csv --out_filename_sorted Turkana.dadi.inference.1D.2Epoch.autosomes.50.runs.sorted.csv
```

```
python ~/softwares/tanya_repos/dadi_tutorial/parse_dadi_expsfs.py --dadi_expsfs 1D_2Epoch_autosomes_run_19/Turkana.dadi.inference.1D.2Epoch.runNum.19.expSFS --num_individuals 101 --theta 207120.395651646 --out_filename 1D_2Epoch_autosomes_run_19/Turkana.dadi.inference.1D.2Epoch.runNum.19.expSFS.normalized.by.theta
```

### With including monomorphic sites
1. Convert SFS to dadi file format
```
python ~/softwares/tanya_repos/dadi_tutorial/convert_sfs_to_dadi_format.py --num_bin 221 --folded_or_unfolded folded --population_name Turkana --sfs_filename ../04_generate_sfs/sfs/autosomes_sfs.txt --num_individuals 110 --out_filename autosomes_sfs_dadi_format_include_monomorphic.sfs --num_monomorphic 211484690  --include_monomorphic yes
```

2. Run dadi
```
for i in {1..50}; do python ~/softwares/tanya_repos/dadi_tutorial/1D.2Epoch.dadi.py --runNum ${i} --pop Turkana --mu 1.5e-8 --L 213127407 --sfs autosomes_sfs_dadi_format_include_monomorphic.sfs --outdir 1D_2Epoch_include_monomorphic/1D_2Epoch_autosomes_run_${i}; done;
```

```
python ~/softwares/tanya_repos/dadi_tutorial/merge_out.py --directory . --run_base_directory 1D_2Epoch_autosomes --file_basename Turkana.dadi.inference.1D.2Epoch.runNum --out_filename Turkana.dadi.inference.1D.2Epoch.autosomes.50.runs.csv --out_filename_sorted Turkana.dadi.inference.1D.2Epoch.autosomes.50.runs.sorted.csv
```

```
python ~/softwares/tanya_repos/dadi_tutorial/parse_dadi_expsfs.py --dadi_expsfs 1D_2Epoch_autosomes_run_3/Turkana.dadi.inference.1D.2Epoch.runNum.3.expSFS --num_individuals 110 --theta 207137.229608562 --out_filename 1D_2Epoch_autosomes_run_3/Turkana.dadi.inference.1D.2Epoch.runNum.3.expSFS.normalized.by.theta
```
