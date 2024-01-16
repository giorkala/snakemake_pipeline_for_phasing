#!/bin/bash
#
# prepape a file with genotypes for trios to phase for SER analysis
# no need to discriminate between rare/common here, as the lowest maf will be >1/600
# input: raw_vcf, array genotypes, list of trios
# + samples.update_ids_{gsa,wes}.txt (a file to change sample IDs, optional)
# + samples.WES_in_GSA.fam (an index of WES samples who are also genotyped, optional)
# Note: for phasing, use smk_02_phase_common.sh with the appropriate input

module load common-apps/bcftools/1.16
plink=/software/team281/bin/plink
mem=16000

### input ###
chr=$1
tag=$2
work_dir=$3
input_wes=$4
input_snp=$5

# lists of samples required:
samples_wes_array="$work_dir/sample_lists/samples.WES_in_GSA.fam"
samples_final="$work_dir/sample_lists/samples.$tag.final_" # this files includes IDs with underscores (GNH-related)
samples_trios="$work_dir/sample_lists/samples.$tag.100trios"
ids_wes_to_array="$work_dir/sample_lists/samples.update_ids_wes.txt"
ids_array_to_wes="$work_dir/sample_lists/samples.update_ids_gsa.txt"
echo "Samples w. both WES and array genotypes:" $(wc -l $samples_wes_array | awk '{print $1}')
echo "Samples among trios to run SER analysis:" $(wc -l $samples_trios | awk '{print $1}')
echo "Samples to keep for downstream analyses:" $(wc -l $samples_final | awk '{print $1}')

wes_new_prefix="$work_dir/sandbox/$tag.trios.chr$chr.exome"
gsa_new_prefix="$work_dir/sandbox/$tag.trios.chr$chr.array"

### output ###
final_out="$work_dir/sandbox/$tag.trios.prepared.chr$chr.bcf"

SECONDS=0
# make a new VCF focusing on the new set of samples and MAC>0 ##
if [ ! -f $wes_new_prefix.vcf.gz ]; then
    echo -e "\nCalling BCFtools to prepare the WES input.\n"
    bcftools view $input_wes -S $samples_trios -Ou | bcftools \
    reheader --samples $work_dir/sample_lists/samples.update_ids_wes.txt | bcftools \
    filter --exclude 'MAF==0.0 | MAX(STRLEN(ALT))>25 | MAX(STRLEN(REF))>25 | FILTER=="ExcessHet"' -Oz -o $wes_new_prefix.vcf.gz
else
    echo -e "\nWES input already exists - moving on.\n"
fi

# prepare the array genotypes
# focusing on the new set of samples (on the original index) - but extracting all samples (not just the trios)
# now subselect individuals, change from the original index to the VCF one, replace "22" with "chr22", and make the vcf in one go
$plink --bfile $input_snp --chr $chr \
    --maf 0.001 \
    --update-ids $work_dir/sample_lists/samples.update_ids_gsa.txt \
    --keep $work_dir/sample_lists/samples.WES_in_GSA.fam \
    --real-ref-alleles \
    --recode vcf-fid bgz \
    --output-chr chr26 \
    --out $gsa_new_prefix \
    --threads 5 --memory $mem

# make a list of overlapping variants to exclude from GSA
echo -e "\nDealing with overlapping markers to exclude."
bcftools view -HG $wes_new_prefix.vcf.gz | awk '{ print $1":"$2}' > $wes_new_prefix.snpid
bcftools view -HG $gsa_new_prefix.vcf.gz | awk '{ print $1":"$2}' > $gsa_new_prefix.snpid
comm -23 <( sort $gsa_new_prefix.snpid) <( sort $wes_new_prefix.snpid) > $gsa_new_prefix.snpid_unique
sed -i 's/:/\t/g' $gsa_new_prefix.snpid_unique
a=$(wc -l $gsa_new_prefix.snpid | awk '{print $1}')
b=$(wc -l $gsa_new_prefix.snpid_unique | awk '{print $1}')
echo "Overlapping variants to remove:" $((a-b))

# merge the two VCFs 
echo -e "\nProceeding with BCFtools and final merging."
bcftools index $gsa_new_prefix.vcf.gz
bcftools index $wes_new_prefix.vcf.gz
# sort individuals in GSA's VCF according to the former
bcftools view $gsa_new_prefix.vcf.gz -S $samples_trios -R $gsa_new_prefix.snpid_unique -Oz -o $gsa_new_prefix.sorted.vcf.gz
bcftools index $gsa_new_prefix.sorted.vcf.gz
echo -e "\nAlmost there, last step in progress..."
bcftools concat -a -Oz -o $final_out.nosnpID $gsa_new_prefix.sorted.vcf.gz $wes_new_prefix.vcf.gz --threads 5

# finally, update the SNP IDs to a uniform index
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' $final_out.nosnpID -Ob -o $final_out
bcftools index -f $final_out

echo "Done with preparation, duration: ${SECONDS}. Cleaning files..."
rm $gsa_new_prefix* $final_out.nosnpID
echo "All done!"
