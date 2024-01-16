#!/usr/bin/env bash
#
# phase rare variants using SHAPEIT5
# minimal script to be used within the Snakemake pipeline
#
# GK - Jul 21st, 2023, 26/9/23
# based on https://github.com/frhl/wes_ko_ukbb/blob/main/scripts/phasing/phasing/02_phase_chunks
#   and on https://odelaneau.github.io/shapeit5/docs/documentation/phase_rare/

# genetic maps downloaded from https://github.com/odelaneau/shapeit5/tree/main/resources/maps
#  chunk lists downloaded from https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38/4cM
# see README.md for more details

chr=$1 # number of chrom to process
phased_scaffold=$2 # the BCF of phased common variants
to_phase=$3 # the *preprocessed* BCF to phase
tag=$4 # just an identifier
gmap=$5 # the genetic map
chunk_list=$6 # a file with a list of chunks

### hyper-parameters ###
threads=5
# we need to both deal with "invalid CONTIG" issues and obtain PPs for singletons. this is enabled by 'shapeit/contig'
module load HGI/common/shapeit/contig
module load common-apps/bcftools/1.16 

### Output prefices ###
work_dir="."
final_out="${work_dir}/phased_genotypes_rare/${tag}.phased_all.chr${chr}.bcf"
out_prefix="${work_dir}/phased_genotypes_rare/chunks/${tag}.phased_all.chr${chr}"

if [ ! -f ${to_phase}.csi ]; then
    echo "Generating the index for $to_phase..."
    bcftools index $to_phase
fi

mkdir -p ${work_dir}/phased_genotypes_rare/chunks/

## paramters for phasing with shapeit
phased_set_error="0.0001" # this is not needed 
pbwt_min_mac=2 # for shapeit5r
pbwt_depth=4 # 5
pbwt_modulo=0.1 # default is 0.1 but 0.0004 ( 0.2 / 50 ) is default value when using --sequencing arugment
pbwt_mdr=0.1
pop_effective_size=15000

SECONDS=0
rm -f $work_dir/sandbox/files.$chr
cat $chunk_list | while read LINE; do
	CHK=$(echo $LINE | awk '{ print $1; }')
	SRG=$(echo $LINE | awk '{ print $3; }')
	IRG=$(echo $LINE | awk '{ print $4; }')

    phased_chunk=$out_prefix.chunk$CHK.bcf
    echo "Proceeding with chunk $CHK."
    echo $phased_chunk >> $work_dir/sandbox/files.$chr

    if [ ! -f ${phased_chunk} ]; then
        phase_rare \
            --input $to_phase \
            --input-region $IRG \
            --scaffold $phased_scaffold \
            --scaffold-region $SRG \
            --map $gmap \
            --pbwt-mac $pbwt_min_mac \
            --pbwt-depth-rare $pbwt_depth \
            --pbwt-modulo $pbwt_modulo \
            --pbwt-mdr $pbwt_mdr \
            --effective-size $pop_effective_size \
            --output $phased_chunk \
            --log $out_prefix.log \
            --thread $threads
            # --phased-set-error $phased_set_error 
    else
        echo "Nothing to phase, $phased_chunk already exists."
    fi
done

echo "Done with phasing, duration: ${SECONDS}."

# second, concatenate the phased chunks
if [ -f ${phased_chunk} ]; then
    # we can either use bcftools concat-ligate (which has an issue), or ligate by SHAPEIT5
    # ls -1v $out_prefix.chunk*.bcf > $work_dir/sandbox/files.chr$chr
    # bcftools concat -n -Ob -o $final_out -f $work_dir/sandbox/files.chr$chr
    ligate --input $work_dir/sandbox/files.$chr --output $final_out --thread 2
    bcftools index $final_out
    rm $work_dir/sandbox/files.$chr
fi
