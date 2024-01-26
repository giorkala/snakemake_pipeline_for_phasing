## Snakemake pipeline for phasing rare variants with SHAPEIT5

In brief, to phase rare WES variants (rule-3) we first need to phase array genotypes (rule-2), combined with common WES variants, and use those as scaffold for the rare ones. Apart from that, we can select a few trios to phase by Mendelian transmission (rule-4) and thus assess the statistical phasing (rule-5). Rule-4 is not necessary but good to have. Computationally, the tough ones are 2 and 3.
#### Pipeline overview (rules)
1. prepare BCF files (common & rare)
2. phase common for the whole cohort (the scaffold)
3. phase rare for the whole cohort, in chunks
4. prepare and phase trios (just common variants; phase using the script for step-2)
5. assess phasing and get files with the PP distribution (needed for plotting).

### Required input ###
* {TAG}.100trios.pedigree: a list of {child,father,mother} triplets for trio phasing. See [here](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/#usage2-phasing-related-samples).
* Several lists of sample IDs, such as mapping between array and WES (optional).
* [Genetic maps](https://github.com/odelaneau/shapeit5/tree/main/resources/maps/b38) and lists of [chunks](https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38/4cM), per chromosome (see the corresponding scripts; what's out for build mismatches)
* A few parameters for Shapeit5 which are already set.
* Data assumed at hand: BED for array; VCF for WES. Fix the corresponding paths in `Snakefile`.

### Important notes ####
* Software dependencies: bcftools (+HTSlib), SHAPEIT5, snakemake (+Python).
* The working directory should be the one from where the scripts will be executed. 
* A few lists of sample IDs -- such as a map between array and WES indices, list of parents, etc -- are required in several parts of the pipeline. These can be generated once, e.g. with the `smk_00_prep_sample_lists.sh`, and are assumed to be in the sample_lists folder. But this file is GNH-specific, so please make any such lists appropriately. In fact, this step could be skipped completely!
* Also, if no trios are sequenced, "phase_trios" and "assess_phasing" are not possible, thus the corresponding rules can be left out.
* Each script requires a few params, such as paths to genetic maps, or chunk lists. These are not pre-defined and need to be passed as arguments.
* For trio phasing, as min-AC ~ 1/300, no rare variants exist, thus phasing is just one step. This is performed with `smk_02_phase_common.sh`, but with the appropriate input (a BCF for trios-only and the pedigree).
* [TODO] The SER analysis is performed using the `switch` tool provided by SHAPEIT5 - this might need to change , as we can similarly work with `bcftools +trio-switch-rate ...`, which requires simpler input.
* Now using ligate instead of concat within phase-rare.
* We need to use Shapeit v5.0.0 ("[preprint version](https://github.com/odelaneau/shapeit5/releases/tag/v1.0.0)" so as to get PP values for singletons. See the discussion [here](https://github.com/odelaneau/shapeit5/issues/56).
* In case you face the "Invalid CONTIG id -1", see [here](https://github.com/odelaneau/shapeit5/issues/34), though the above solution should work anyway. These are the reasons why a combination of binaries and modules for Shapeit5; please choose what's best.

#### How to - LSF ####
* One way to make the pipeline work is to comment the rule "run_all" appropriately in `Snakefile`.
* the following will start the pipeline and submit all jobs with all the jobs required:
`snakemake --cluster "bsub -M 16G -R 'select[mem>16G] rusage[mem=16G]' -n10 -G team281 -q normal -o run_all.stdout -e run_all.stderr" --jobs 10 --latency-wait 10 -T 2 run_all`
* first run the following the initialise the required folders:
`mkdir -p phased_genotypes_common phased_genotypes_rare phasing_assessment sandbox`

GK - Sep 26th, 2023
