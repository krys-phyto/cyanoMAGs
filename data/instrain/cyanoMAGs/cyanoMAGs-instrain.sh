# kjkibler 2022-Dec-10
# Installing instrain and then running profile on all bins


### Instrain ###


# Install
conda install -c conda-forge -c bioconda -c defaults instrain


### Files ###
# Bam files
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2/*.bam

# Refined bins
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/*-contigs.fa

### Code ###
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/
mkdir instrain
mkdir instrain/profile

# making a genomes.stb file since it hates my scaffolds_to_bin file
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/
parse_stb.py --reverse -f *-contigs.fa  -o genomes-refined.stb

# prodigal
prodigal -i 16cyanoMags-refined.fna -d 16cyanoMags-refined.fna.genes.fna


# inStrain profile
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2/

for file in *.bam
do
  echo $file
    inStrain profile $file /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna  \
    -o /home/glbrc.org/kjkibler/paper-1-cyanomags/data/instrain/profile/${file/.bam/.IS} -p 6 \
    -g /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna.genes.fna \
    -s /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/genomes-refined.stb \
    --min_genome_coverage 0.01
done

for file in *.bam
do
  echo $file
  inStrain profile_genes  \
    -i /home/glbrc.org/kjkibler/paper-1-cyanomags/data/instrain/profile/${file/.bam/.IS} -p 6 \
    -g /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna.genes.fna
done

for file in *.bam
do
  echo $file
  inStrain genome_wide  \
    -i /home/glbrc.org/kjkibler/paper-1-cyanomags/data/instrain/profile/${file/.bam/.IS} -p 6 \
    -s /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/genomes-refined.stb
done

for file in *.bam
do
  echo $file
  inStrain plot -i /home/glbrc.org/kjkibler/paper-1-cyanomags/data/instrain/profile/${file/.bam/.IS} -p 6 \
  -p1 -p2 -p3 -p4 -p5 -p6 -p7 -p8 -p9
done


# rerun with genome coverage min 5
for file in *.OR.bam
do
  echo $file
    inStrain profile $file /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna  \
    -o /home/glbrc.org/kjkibler/paper-1-cyanomags/data/instrain/profile-5/${file/.bam/.IS} -p 6 \
    -g /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna.genes.fna \
    -s /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/genomes-refined.stb \
    --min_genome_coverage 5
done

# rerun with genome coverage min 3
for file in 2012*.OR.bam
do
  echo $file
    inStrain profile $file /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna  \
    -o /home/glbrc.org/kjkibler/paper-1-cyanomags/data/instrain/profile-3/${file/.bam/.IS} -p 6 \
    -g /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna.genes.fna \
    -s /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/genomes-refined.stb \
    --min_genome_coverage 3
done

# rerun with genome coverage min 1
for file in *.OR.bam
do
  echo $file
    inStrain profile $file /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna  \
    -o /home/glbrc.org/kjkibler/paper-1-cyanomags/data/instrain/profile-1/${file/.bam/.IS} -p 6 \
    -g /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags-refined.fna.genes.fna \
    -s /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/genomes-refined.stb \
    --min_genome_coverage 1
done





# COMPARE
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/instrain/profile-5
inStrain compare -i 20080813.OR.IS 20081017.OR.IS 20090609.OR.IS 20100613.OR.IS 20100716.OR.IS 20100830.OR.IS 20110628.OR.IS 20120713.OR.IS 20120713.OR.IS 20120720.OR.IS 20121105.OR.IS 20121116.OR.IS \
  -o compare-select  \
  -s /home/glbrc.org/kjkibler/paper-1-cyanomags/data/all-gene/instrain/genome-list.txt
