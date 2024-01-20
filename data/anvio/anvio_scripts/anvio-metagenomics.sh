# kjkibler 2022-Nov-22
# Go through anvio metagenomics tutorial to generate contig database and summary profile statistic

### --- Anvio installation --- ###
# https://anvio.org/install/

### Code ###
conda --version
 #
conda update conda

conda create -y --name anvio-7.1 python=3.6

conda activate anvio-7.1

conda install -y -c bioconda "sqlite>=3.31.1"
conda install -y -c bioconda prodigal
conda install -y -c bioconda mcl
conda install -y -c bioconda muscle=3.8.1551
conda install -y -c bioconda hmmer
conda install -y -c bioconda diamond
conda install -y -c bioconda blast
conda install -y -c bioconda megahit
conda install -y -c bioconda spades
conda install -y -c bioconda bowtie2 tbb=2019.8
conda install -y -c bioconda bwa
conda install -y -c bioconda samtools=1.9
conda install -y -c bioconda centrifuge
conda install -y -c bioconda trimal
conda install -y -c bioconda iqtree
conda install -y -c bioconda trnascan-se
conda install -y -c bioconda r-base
conda install -y -c bioconda r-stringi
conda install -y -c bioconda r-tidyverse
conda install -y -c bioconda r-magrittr
conda install -y -c bioconda r-optparse
conda install -y -c bioconda bioconductor-qvalue
conda install -y -c bioconda fasttree
conda install -y -c bioconda vmatch
conda install -y -c bioconda fastani

curl -L https://github.com/merenlab/anvio/releases/download/v7.1/anvio-7.1.tar.gz \
        --output anvio-7.1.tar.gz

pip install anvio-7.1.tar.gz

anvi-self-test --suite mini --no-interactive

anvi-setup-scg-taxonomy
anvi-setup-ncbi-cogs --reset
anvi-setup-kegg-kofams --reset

### --- Anvio Metagenomics --- ###
# https://merenlab.org/2016/06/22/anvio-tutorial-v2/

### Files ###
# Concatenated fasta file of genomes
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags.fna

# Bam files
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2/*.bam

### Code ###
# Rename raw bam files to name-RAW.bam #
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2
for f in *.bam
do
  mv $f ${f/.bam/-RAW.bam}
done

# Make/Move to working directory #
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data
mkdir anvio
mkdir anvio/anvi-metagenomics
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics

# Generate Contig Database #
#$ anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'An example contigs database'
anvi-gen-contigs-database --help

anvi-gen-contigs-database -f /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags.fna \
  -o 16cyanoMags.db -n 'LM 16 Cyano Mags' --split-length 0

# Run Hmms
#$ anvi-run-hmms -c contigs.db

anvi-run-hmms -c 16cyanoMags.db --num-threads 4

# NCBI cogs/functions
#$ anvi-run-ncbi-cogs -c contigs.db

anvi-run-ncbi-cogs -c 16cyanoMags.db --num-threads 4
# offending protien ids that were unable to be annotated
#'WP_071782595.1, WP_026788925.1, WP_015109478.1, ABW27134.1, WP_127054283.1, BAG00033.1, WP_042151191.1, AVZ30803.1, AFY39127.1, AFY33857.1, WP_071782561.1'

# Profiling BAM files #

# initialize bam biles
#$ anvi-init-bam SAMPLE-01-RAW.bam -o SAMPLE-01.bam
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2
for f in *.bam
do
anvi-init-bam $f -o ${f/-RAW.bam/.bam}
done

# profiling bam files
#$ anvi-profile -i SAMPLE-01.bam -c contigs.db
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics
mkdir profile
for f in /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2/*OR.bam
do
  file=`basename $f`
  file=${file/.OR.bam/}
  anvi-profile -i $f -c 16cyanoMags.db \
    --output-dir profile/$file \
    --sample-name metag_$file \
    --profile-SCVs \
    --num-threads 8
done

# merge profiles
#$ anvi-merge */PROFILE.db -o SAMPLES-MERGED -c contigs.db

cd profile
anvi-merge */PROFILE.db -o 16cyanoMags-MERGED -c ../16cyanoMags.db

# import bin info

# fix genomes.stb to have non-number at beginning of bin names
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags
parse_stb.py --reverse -f *-fixed.fna  -o genomes.stb
sed -i 's/-fixed.fna//g' genomes.stb
sed -i 's/\.//g' genomes.stb
sed -i 's/-/_/g' genomes.stb
awk '$2="mag_"$2 {print}' genomes.stb > genomes1.stb
mv genomes1.stb genomes.stb
cat genomes.stb | tr [:blank:] \\t > genomes1.stb
mv genomes1.stb genomes.stb


#$ anvi-import-collection binning_results.txt -p SAMPLES-MERGED/PROFILE.db -c contigs.db --source "SOURCE_NAME"
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile

anvi-import-collection /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/genomes.stb \
-p 16cyanoMags-MERGED/PROFILE.db -c ../16cyanoMags.db -C LM_16cyanoMags --contigs-mode


### Summary and refinement ###
#$ anvi-summarize -p SAMPLES-MERGED/PROFILE.db -c contigs.db -o SAMPLES-SUMMARY -C CONCOCT
anvi-summarize -p 16cyanoMags-MERGED/PROFILE.db -c ../16cyanoMags.db -o ../Original-SUMMARY -C LM_16cyanoMags

#refine bins
# 4 bins that had above 5% redundancy, refining to below 5% while maintaining completion
anvi-refine -p 16cyanoMags-MERGED/PROFILE.db -c ../16cyanoMags.db -C LM_16cyanoMags \
    -b mag_2582580577

anvi-refine -p 16cyanoMags-MERGED/PROFILE.db -c ../16cyanoMags.db -C LM_16cyanoMags \
    -b mag_3300020490_bin6

anvi-refine -p 16cyanoMags-MERGED/PROFILE.db -c ../16cyanoMags.db -C LM_16cyanoMags \
    -b mag_3300020542_bin8

anvi-refine -p 16cyanoMags-MERGED/PROFILE.db -c ../16cyanoMags.db -C LM_16cyanoMags \
    -b mag_GEODES117_bin31


# new summary with refined bins
anvi-summarize -p 16cyanoMags-MERGED/PROFILE.db -c ../16cyanoMags.db -o ../Refined-SUMMARY -C LM_16cyanoMags

# move refined bins to mag data folder
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/Refined-SUMMARY/bin_by_bin
for f in *
do
cd $f
filename=`basename *-contigs.fa`
cp $filename /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/Refined-SUMMARY/bin_by_bin
done
