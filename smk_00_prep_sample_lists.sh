#!/bin/bash
#

module load common-apps/bcftools/1.16
tag='GNH_39k'
work_dir="./sample_lists"

input_wes="/FIXTHIS/filtered_vcfs/chr21_hard_filters.vcf.gz"
input_gsa="/FIXTHIS/GSA_QC/GNH_GSAv3EAMD_DupExcl_autosome_atcgSNP"
pedigree="/FIXTHIS/NEW_RUN/GNH_39k.100trios.pedigree"
link_OrageneID="/FIXTHIS/link_OrageneID_all-WES_GSA.txt"

# main output (+more files)
samples_trios="$work_dir/samples.$tag.100trios"
samples_final="$work_dir/samples.$tag.final_"

mkdir -p sample_lists
# create an index of samples in the WES cohort
bcftools query $input_wes --list samples > $work_dir/samples.WES_all.txt
# find the intersection of those with the GSA cohort - the link file contains more
join -j 1 $work_dir/samples.WES_all.txt <(awk 'NR>1{if($2>0){print $0} }' $link_OrageneID | sort) > $work_dir/samples.common_WES_GSA.link 
# for t in 2 3 4; do sed -i -E "s/\_$t//g" $work_dir/samples.common_WES_GSA.link; done
awk '{print $1}' $work_dir/samples.common_WES_GSA.link > $work_dir/samples.WES_in_GSA.txt
# fix names with underscores for WES
awk '{print $1}' $work_dir/samples.common_WES_GSA.link | grep "_" > $work_dir/samples.WES_with_underscore.txt
paste $work_dir/samples.WES_with_underscore.txt <(sed -e 's/\_2//g' -e 's/\_3//g' -e 's/\_4//g' $work_dir/samples.WES_with_underscore.txt) > $work_dir/samples.update_ids_wes.txt
# also create a new linkage file according to PLINK's required format for updating IDs
join -1 2 -2 1 <(cat $work_dir/samples.common_WES_GSA.link | sort -t' ' -k2 ) <(cat $input_gsa.fam | sort) | awk '{print $1,$3,$2,$2}' > $work_dir/samples.update_ids_gsa.txt
sed -i -e 's/\_2//g' -e 's/\_3//g' -e 's/\_4//g' $work_dir/samples.update_ids_gsa.txt
awk '{print $1,$1}' $work_dir/samples.common_WES_GSA.link | sed -e 's/\_2//g' -e 's/\_3//g' -e 's/\_4//g' > $work_dir/samples.WES_in_GSA.fam 

echo "Generate a list of $tag and the corresponding parents"
# NOTE: This is now performed with python according to mendel-errors
cut -f2,3 $pedigree | tr "\t" "\n" |  sort | uniq > $work_dir/samples.$tag.parents 
cat $pedigree | tr "\t" "\n" | sort  | uniq > $samples_trios

comm -23 $work_dir/samples.WES_in_GSA.txt $work_dir/samples.$tag.parents > $work_dir/samples.$tag.final_

echo "Samples w. both WES and array genotypes:" $(wc -l $work_dir/samples.WES_in_GSA.fam  | awk '{print $1}')
echo "Samples among trios to run SER analysis:" $(wc -l $samples_trios | awk '{print $1}')
echo "Samples to keep for downstream analyses:" $(wc -l $samples_final | awk '{print $1}')
