#!/usr/bin/env bash
#
# phase common variants using Shapeit5
# minimal script to be used within the Snakemake pipeline
#
# GK - Jul 14th, 2023, 26/9/23

chr=$1 # number of chrom to process
to_phase=$2 # the *preprocessed* BCF to phase
tag=$3 # a file identifier
gmap=$4 # the genetic map
pedigree=$5 # if given, the pedigree for phasing trios

out_prefix="./phased_genotypes_common/$tag.phased.chr${chr}"
phased="${out_prefix}.bcf"

module load common-apps/bcftools/1.16
# module load HGI/common/shapeit/contig
SHAPEIT_phase_common='/software/team281/bin/shapeit5/phase_common_static'
threads=10

### hyper-parameters ###
pbwt_modulo=0.1 # deault is 0.1, in shapeit4 with sequencing it is 0.0002
pbwt_depth=4 # deault is 4
pbwt_mac=5 # deafult is 5
pbwt_mdr=0.1 # default is 0.1
min_maf=0.001 # no default value

if [ ! -f ${to_phase}.csi ]; then
    echo "Generating the index for $to_phase..."
    bcftools index $to_phase
fi

SECONDS=0

if [[ -z $5 ]]; then
    echo -e "\nPhasing without trio information.\n"
    ${SHAPEIT_phase_common} \
      --input ${to_phase} \
      --map ${gmap} \
      --region "chr${chr}" \
      --thread ${threads} \
      --output "${phased}" \
      --log "${out_prefix}.log" \
      --pbwt-modulo ${pbwt_modulo} \
      --pbwt-depth ${pbwt_depth} \
      --pbwt-mac ${pbwt_mac} \
      --pbwt-mdr ${pbwt_mdr} \
      --filter-maf ${min_maf} \
      && echo "Finished phasing variants for chr${chr}, out: ${phased}, duration: ${SECONDS}" \
      || echo "ERROR: Phasing variants failed for chr${chr}, duration: ${SECONDS}"

else
    echo -e "\nPhasing using the pedigree in $pedigree.\n"
    ${SHAPEIT_phase_common} \
      --input ${to_phase} \
      --map ${gmap} \
      --region "chr${chr}" \
      --thread ${threads} \
      --output "${phased}" \
      --log "${out_prefix}.log" \
      --pbwt-modulo ${pbwt_modulo} \
      --pbwt-depth ${pbwt_depth} \
      --pbwt-mac ${pbwt_mac} \
      --pbwt-mdr ${pbwt_mdr} \
      --filter-maf ${min_maf} \
      --pedigree $pedigree \
      && echo "Finished phasing variants for chr${chr}, out: ${phased}, duration: ${SECONDS}" \
      || echo "ERROR: Phasing variants failed for chr${chr}, duration: ${SECONDS}"
fi
