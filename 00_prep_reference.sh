cd 

cd references

samtools faidx GCA_014839835.1_bDryPub1.pri_genomic.fna

bwa index GCA_014839835.1_bDryPub1.pri_genomic.fna

cd ../

java -jar picard.jar CreateSequenceDictionary R=/home/jmanthey/references/GCA_014839835.1_bDryPub1.pri_genomic.fna \
O=/home/jmanthey/references/GCA_014839835.1_bDryPub1.pri_genomic.dict
