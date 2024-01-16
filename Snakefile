"""
Snakefile pipeline for the phasing common and rare variants with SHAPEIT5

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

Check the `README.md` for more details.

# for LSF:
`snakemake --cluster "bsub -M 16G -R 'select[mem>16G] rusage[mem=16G]' -n10 -G team281 -q normal -o run_all.stdout -e run_all.stderr" --jobs 22 --latency-wait 10 -T 2 run_all`

GK - Sep 26th, 2023
"""

## General input and parameters ##
input_wes_prefix="/FIXTHIS/filtered_vcfs"
input_bed_prefix="/FIXTHIS/GSA_QC/GNH_GSAv3EAMD_DupExcl_autosome_atcgSNP"
genet_map_prefix="/FIXTHIS/genetic_maps"
chunk_list_prefix="/FIXTHIS/chunks_b38_4cM/chunks_"

WORK_DIR="/FIXTHIS/WORK_DIR"
TAG="GNH_39k"
PEDIGREE=WORK_DIR+"/GNH_39k.100trios.pedigree"

rule prepare_bcfs:
    resources:
        mem_mb=16000,
        threads=5
    input:
        raw_wes = input_wes_prefix+"/chr{chrom}_hard_filters.vcf.gz",
    output:
        bcf_common = "{wd}/sandbox/{tag}.notrios.common_merged.chr{chrom}.bcf",
        bcf_rare   = "{wd}/sandbox/{tag}.notrios.rare_prepared.chr{chrom}.bcf"
    run:
        shell("bash smk_01_prep_all_bcf.sh {wildcards.chrom} {TAG} {WORK_DIR} {input.raw_wes} "+input_bed_prefix)

rule phase_common:
    resources:
        mem_mb=16000,
        threads=10
    input:
        input_bcf = rules.prepare_bcfs.output.bcf_common,
        genet_map = genet_map_prefix + "/chr{chrom}.b38.gmap.gz"
    output:
        phased_bcf = "{wd}/phased_genotypes_common/{tag}.notrios.phased.chr{chrom}.bcf"
    shell:
        """
        bash smk_02_phase_common.sh {wildcards.chrom} {input.input_bcf} {TAG}.notrios {input.genet_map}
        """

rule phase_rare:
    resources:
        mem_mb=16000,
        threads=5
    input:
        input_bcf = rules.prepare_bcfs.output.bcf_rare,
        scaffold = rules.phase_common.output.phased_bcf,
        genet_map = genet_map_prefix + "/chr{chrom}.b38.gmap.gz",
        chunk_list = chunk_list_prefix + "chr{chrom}.txt"
    output:
        phased_bcf = "{wd}/phased_genotypes_rare/{tag}.notrios.phased_all.chr{chrom}.bcf"
    shell:
        """
        bash smk_03_phase_rare.sh {wildcards.chrom} {input.scaffold} {input.input_bcf} {TAG}.notrios {input.genet_map} {input.chunk_list}
        """

rule phase_trios:
    resources:
        mem_mb=16000,
        threads=10
    input:
        raw_wes = input_wes_prefix + "/chr{chrom}_hard_filters.vcf.gz",
        genet_map = genet_map_prefix + "/chr{chrom}.b38.gmap.gz"
    output:
        trios_prep = "{wd}/sandbox/{tag}.trios.prepared.chr{chrom}.bcf",
        phased_bcf = "{wd}/phased_genotypes_common/{tag}.trios.phased.chr{chrom}.bcf"
    shell:
        """
        bash smk_04_phase_trios.sh {wildcards.chrom} {TAG} {WORK_DIR} {input.raw_wes} {input_bed_prefix}
        bash smk_02_phase_common.sh {wildcards.chrom} {output.trios_prep} {TAG}.trios {input.genet_map} {PEDIGREE}
        """

rule assess_phasing:
    resources:
        threads=5
    input:
        phased_common= rules.phase_common.output.phased_bcf,
        phased_rare  = rules.phase_rare.output.phased_bcf,
        phased_trios = rules.phase_trios.output.phased_bcf
    output:
        ser_comn = "{wd}/phasing_assessment/{tag}.assess_common.chr{chrom}.variant.switch.txt.gz",
        ser_rare = "{wd}/phasing_assessment/{tag}.assess_rare.pp0.50.chr{chrom}.variant.switch.txt.gz"
    run:
        prefix = WORK_DIR + "/phasing_assessment/" + TAG
        shell("bash smk_05_assess_phasing.sh common {wildcards.chrom} {input.phased_trios} {input.phased_common} {PEDIGREE} {prefix}")
        shell("bash smk_05_assess_phasing.sh rare {wildcards.chrom} {input.phased_trios} {input.phased_rare} {PEDIGREE} {prefix}")

rule get_pp_distribution:
    resources:
        threads=4
    input:
        phased_rare  = rules.phase_rare.output.phased_bcf
    output:
        pp_mac0_pp90 = "{wd}/phasing_assessment/{tag}.pp.chr{chrom}.maf00015.pp0.90.gz",
        pp_mac0_pp50 = "{wd}/phasing_assessment/{tag}.pp.chr{chrom}.maf00015.pp0.50.gz"
    run:
        prefix = WORK_DIR + "/phasing_assessment/" + TAG
        shell("bash smk_05_assess_phasing.sh getpp {wildcards.chrom} {input.phased_rare} {prefix}.pp")


rule run_all:
    input:
        expand(WORK_DIR+"/phased_genotypes_common/{tag}.trios.phased.chr{chrom}.bcf", chrom=range(1,23), tag={TAG}),
        expand(WORK_DIR+"/phased_genotypes_rare/{tag}.notrios.phased_all.chr{chrom}.bcf", chrom=range(1,23), tag={TAG}),   
        expand(WORK_DIR+"/phasing_assessment/{tag}.assess_rare.pp0.50.chr{chrom}.variant.switch.txt.gz", chrom=range(1,23), wd={WORK_DIR}, tag={TAG}),
        # expand(WORK_DIR+"/phasing_assessment/{tag}.assess_common.chr{chrom}.variant.switch.txt.gz", chrom=range(21,22), tag={TAG}),
        expand(WORK_DIR+"/phasing_assessment/{tag}.pp.chr{chrom}.maf00015.pp0.90.gz", chrom=range(1,23), tag={TAG}),
