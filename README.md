## Snakefile pipeline for the phasing common and rare variants with SHAPEIT5

The working directory will be the one from where this script will be executed
The pipeline consists of 4 small steps (small conceptually, might be tough computationally)
    1) prepare BCF files (common & rare)
    2) phase common for the whole cohort (the scaffold)
    3) phase rare for the whole cohort, in chunks
    4) prepare and phase trios (just common variants; phase using the script for step-2)
    5) assess phasing and get files with the PP distribution

#### Required input ####
1. {TAG}.100trios.pedigree: a list of {child,father,mother} triplets for trio phasing
2. Several lists of sample IDs, such as mapping between array and WES (optional)
3. Genetic maps and lists of chunks, per chromosome (see the corresponding scripts)
4. A few parameters for Shapeit5 which are already set

### Important notes ####
* [TODO] A few lists of sample IDs -- such as a map between array and WES indices, list of parents, etc -- are required in several parts of the pipeline. These can be generated once, e.g. with the `smk_00_prep_sample_lists.sh`, and are assumed to be in the sample_lists folder. But this file is GNH-specific, so please make any such lists appropriately.
* Each script requires a few params, such as paths to genetic maps, or chunk lists. These are not pre-defined and need to be passed as arguments.
* Software dependencies: bcftools (+HTSlib), SHAPEIT5, snakemake (+Python).
* For trio phasing, as min-AC ~ 1/300, no rare variants exist, thus phasing is just one step. This is performed with `smk_02_phase_common.sh`, but with the appropriate input (a BCF for trios-only and the pedigree).
* I'm currently using a combination of binaries and modules for Shapeit5 (to deal with issues 34,56), please choose what's best.
* [TODO] The SER analysis is performed using the `switch` tool provided by SHAPEIT5 - this might need to change , as we can similarly work with `bcftools +trio-switch-rate ...`, which requires simpler input.
* Now using ligate instead of concat within phase-rare.

#### How to - LSF ####
* the following will start the pipeline and submit all jobs with all the jobs required:
`snakemake --cluster "bsub -M 16G -R 'select[mem>16G] rusage[mem=16G]' -n10 -G team281 -q normal -o run_all.stdout -e run_all.stderr" --jobs 10 --latency-wait 10 -T 2 run_all`
* first run the following the initialise the required folders:
`mkdir -p phased_genotypes_common phased_genotypes_rare phasing_assessment sandbox`

GK - Sep 26th, 2023