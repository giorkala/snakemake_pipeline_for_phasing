#!/bin/bash

# GK - Jul 28th, 2023, Sep 26/9/23
# based on https://odelaneau.github.io/shapeit5/docs/documentation/switch/
# note, all files should be BCF for higher efficiency
# see README.md for more details

mode=$1
chr=$2
phased_true=$3
phased_estd=$4
pedigree=$5
out_prefix=$6

# module load HGI/common/shapeit/contig
SHAPEIT_switch='/software/team281/bin/shapeit5/switch_static'
module load common-apps/bcftools/1.16
threads=5

if [ $mode = "common" ]; then

    echo -e "\nWork for common variants at $out_prefix"
    # out_prefix="$out_dir/$tag.assess_common.chr$chr"
    $SHAPEIT_switch --validation $phased_true --estimation $phased_estd \
    --region chr$chr --pedigree $pedigree \
    --output $out_prefix.assess_common.chr$chr --thread $threads
    zcat $out_prefix.assess_common.chr$chr.sample.switch.txt.gz | awk 'BEGIN { e=0; t=0; } { e+=$2; t+=$3; } END { print "SER =", e*100/t; }'

elif [ $mode = "rare" ]; then
    
    echo -e "\nWork for rare variants at $out_prefix"

    # get list of the rare variants phased for each cohort
    bcftools query -f '%CHROM:%POS:%REF:%ALT\t%POS\t%ID\t%AC\t%AN\n' $phased_true > $out_prefix.trios.chr$chr.snpinfo
    bcftools query -f '%CHROM:%POS:%REF:%ALT\t%POS\t%ID\t%AC\t%AN\n' $phased_estd > $out_prefix.notrios.chr$chr.snpinfo

    for PP in 0.50 0.80 0.90; do 
        new_prefix="$out_prefix.assess_rare.pp$PP.chr$chr"
        $SHAPEIT_switch --validation $phased_true --estimation $phased_estd --min-pp $PP \
        --singleton --pedigree $pedigree \
        --region chr$chr --output $new_prefix --thread $threads --log $new_prefix.log
    done

    zcat $new_prefix.sample.switch.txt.gz | awk 'BEGIN { e=0; t=0; } { e+=$2; t+=$3; } END { print "SER =", e*100/t; }'

elif [ $mode = "getpp" ]; then

    BCF=$3
    out_prefix=$4
    echo "Extracting PP summaries for $BCF (multiple tasks)"

    # MAF<0.001 + epsilon
    bcftools view $BCF --max-af 0.0015 -Ou | bcftools query -i'GT="het"' -f'[%CHROM:%POS:%REF:%ALT %GT %PP \n]' | awk '($3>=0.50)' | gzip > $out_prefix.chr$chr.maf00015.pp0.50.gz &
    bcftools view $BCF --max-af 0.0015 -Ou | bcftools query -i'GT="het"' -f'[%CHROM:%POS:%REF:%ALT %GT %PP \n]' | awk '($3>=0.90)' | gzip > $out_prefix.chr$chr.maf00015.pp0.90.gz &
    # MAC==1
    bcftools view $BCF -C 1 -Ou | bcftools query -i'GT="het"' -f'[%CHROM:%POS:%REF:%ALT %GT %PP \n]' | gzip > $out_prefix.chr$chr.mac1.gz &
    # MAC==2
    bcftools view $BCF -c 2 -C 2 -Ou | bcftools query -i'GT="het"' -f'[%CHROM:%POS:%REF:%ALT %GT %PP \n]' | gzip > $out_prefix.chr$chr.mac2.gz &
    # MAC==3:
    bcftools view $BCF -c 3 -C 3 -Ou | bcftools query -i'GT="het"' -f'[%CHROM:%POS:%REF:%ALT %GT %PP \n]' | gzip > $out_prefix.chr$chr.mac3.gz

else
    echo "Error: unrecognised argument: mode needs to be one of `common`,`rare`,`getpp`."
fi
