# put all the scaffolds in different files

for i in $( cut -f1 GCA_014839835.1_bDryPub1.pri_genomic.fna.fai ); do 
echo $i > test.txt; 
seqtk subseq GCA_014839835.1_bDryPub1.pri_genomic.fna test.txt > scaffolds/$i.fasta
done

grep "CM" GCA_014839835.1_bDryPub1.pri_genomic.fna.fai > rm.fai


# run the following script in a submission job

#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=RM
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=18
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-47

input=$( head -n${SLURM_ARRAY_TASK_ID} rm.fai | tail -n1 | cut -f1 )

# use a custom repeatmasker database to annotate the genome for TEs 
# includes the:
# RepBase vertebrate database v24.03 sequences
# certhia americana custom repeatmodeler sequences (doi: 10.1093/gbe/evab120)
# colaptes auratus custom repeatmodeler sequences (doi:10.1093/g3journal/jkaa026)

# run repeat masker v1.332
cd /lustre/scratch/jmanthey/09_dryobates_repeatmasker/scaffolds
RepeatMasker -pa 18 -s -lib ~/RepeatMasker/Libraries/custom_library_certhia_colaptes.fa ${input}.fasta

