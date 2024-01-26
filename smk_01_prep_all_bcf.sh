#!/bin/bash
#
# prepare all the BCF files required for phasing - this is treaky as several files are needed
# which should be included in ./sample_lists/
# In brief, to merge WES and array variants we proceed as follows:
# 1. extract "common" variants (AF>0.001)
# 2. we filter-out anything tagged as "ExcessHet" (you might have to remove this filter)
# 3. we filter-out indels longer than 25bp (if any; see discussion with high numbers of pLoF genotypes)
# we save those to a VCF. Then we prepare array variants and make another VCF
# For any overlapping variants, we chose those in the WES dataset. Then merge the two VCFs, and save as BCF.
# To prepare the rare variants we follow a similar procedure but simpler.
# Note that if no trios are available, `samples.$tag.100trios` should be an empty file to get things working.
# Also, should exclude variants in LCRs, but that can be dealt with in the annotation.

### resources ###
module load common-apps/bcftools/1.16
plink=/software/team281/bin/plink
mem=16000
threads=5

### input ###
chr=$1
tag=$2
work_dir=$3
input_wes=$4
input_snp=$5

# lists of samples required:
samples_trios="$work_dir/sample_lists/samples.$tag.100trios"
samples_wes_array="$work_dir/sample_lists/samples.WES_in_GSA.fam"
samples_final="$work_dir/sample_lists/samples.$tag.final" # EDIT THIS ACCORDINGLY
# samples_final1="$work_dir/samples.$tag.final" # this file is WITHOUT underscores, might not exist at first
# update_ids_wes="$work_dir/sample_lists/samples.update_ids_wes.txt" # not needed in general
# update_ids_array="$work_dir/sample_lists/samples.update_ids_gsa.txt" # not needed in general
echo "Samples w. both WES and array genotypes:" $(wc -l $samples_wes_array | awk '{print $1}')
echo "Samples among trios to run SER analysis:" $(wc -l $samples_trios | awk '{print $1}')
echo "Samples to keep for downstream analyses:" $(wc -l $samples_final | awk '{print $1}')
# note that `samples_final` should contain all samples except the parents of the trios, if any!

# intermediate files
wes_new_prefix="$work_dir/sandbox/$tag.notrios.chr$chr.exome"
snp_new_prefix="$work_dir/sandbox/$tag.notrios.chr$chr.array"

### output ###
out_bcf_full_common="${work_dir}/sandbox/${tag}.notrios.common_merged.chr${chr}.bcf"
out_bcf_full_rare="${work_dir}/sandbox/${tag}.notrios.rare_prepared.chr${chr}.bcf"

SECONDS=0

### A. Extract common genotypes from raw WES ###
# skip rare variants (MAF<0.001), those that are ExcessHet, or those being long indels
if [ ! -f $wes_new_prefix.vcf.gz ]; then
    echo -e "\nCalling BCFtools to prepare the WES input.\n"
    bcftools \
        view $input_wes -S $samples_final -Ou | bcftools \
        filter --include 'MAF>0.001 & FILTER!="ExcessHet" & MAX(STRLEN(ALT))<25 & MAX(STRLEN(REF))<25' -Oz -o $wes_new_prefix.vcf.gz
    bcftools index $wes_new_prefix.vcf.gz
else
    echo -e "\nWES input already exists - moving on.\n"
fi

# if [ ! -f $samples_final1 ]; then
#     echo "Exctract the final list of available samples:"
#     bcftools query $wes_new_prefix.vcf.gz --list-samples > $samples_final1
# fi

### B. Extract array genotypes from raw GSA ###

echo -e "\nCalling PLINK to prepare the GSA input."
# subselect individuals, replace "22" with "chr22", and make the vcf in one go
$plink --bfile $input_snp --chr $chr \
    --maf 0.001 \
    --keep $samples_wes_array \
    --real-ref-alleles \
    --recode vcf-fid bgz \
    --output-chr chr26 \
    --out $snp_new_prefix \
    --threads $threads --memory $mem
    # --update-ids $update_ids_array \
bcftools index $snp_new_prefix.vcf.gz

### C. Merge common WES with array ###

# first, make a list of overlapping variants to exclude from GSA
echo -e "\nExtracting overlapping variants between WES and array..."
bcftools view -HG $wes_new_prefix.vcf.gz | awk '{ print $1":"$2}' > $wes_new_prefix.snpid
bcftools view -HG $snp_new_prefix.vcf.gz | awk '{ print $1":"$2}' > $snp_new_prefix.snpid
comm -23 <( sort $snp_new_prefix.snpid) <( sort $wes_new_prefix.snpid) > $snp_new_prefix.snpid_unique
sed -i 's/:/\t/g' $snp_new_prefix.snpid_unique
a=$(wc -l $snp_new_prefix.snpid | awk '{print $1}')
b=$(wc -l $snp_new_prefix.snpid_unique | awk '{print $1}')
echo "Overlapping variants to remove:" $((a-b))

echo -e "\nProceeding with BCFtools and final merging."
bcftools view $snp_new_prefix.vcf.gz -S $samples_final -R $snp_new_prefix.snpid_unique -Oz -o $snp_new_prefix.sorted.vcf.gz
bcftools index $snp_new_prefix.sorted.vcf.gz
echo -e "\nAlmost there, perform merging..."
bcftools concat -a -Oz -o $out_bcf_full_common.nosnpID.vcf.gz $snp_new_prefix.sorted.vcf.gz $wes_new_prefix.vcf.gz --threads 5

# finally, update the SNP IDs to a uniform index
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' $out_bcf_full_common.nosnpID.vcf.gz -Ob -o $out_bcf_full_common
bcftools index -f $out_bcf_full_common

rm $out_bcf_full_common.nosnpID* $snp_new_prefix.sorted*
echo "Done with preparation of common, duration: ${SECONDS}."

### D. Extract rare genotypes from raw WES ###

echo -e "\n################\n"
echo "Preparing a new BCF for rare variants..."
SECONDS=0
# here we have select rare variants and QC them, then select samples
bcftools \
    view -i 'MAF<0.001 & MAX(STRLEN(ALT))<25 & MAX(STRLEN(REF))<25 & FILTER!="ExcessHet"' -S $samples_final | bcftools \
    annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' -Ob -o $out_bcf_full_rare
bcftools index $out_bcf_full_rare
echo "Done with preprocessing rare, duration: ${SECONDS}."

echo "All done!"
