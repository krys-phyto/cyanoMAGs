# Krys Kibler
# 2023-07-24

# InStrain

cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/aph-instrain/


prodigal -i mag_3300020490_bin6-contigs.fa -d mag_3300020490_bin6-contigs.fa.genes.fna

# mapping
mkdir bt2

bowtie2-build mag_3300020490_bin6-contigs.fa \
  bt2/mag_3300020490_bin6-contigs.fa

for f in /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/fastq/*.fastq
do
  echo $f
  NAME="$(basename $f .fastq)"
  bowtie2 -p 6 -x bt2/mag_3300020490_bin6-contigs.fa --very-sensitive \
    -q $f > bt2/$NAME.sam
done

# profile

mkdir instrain
mkdir instrain/profile

for file in bt2/*.sam
do
  echo $file
  NAME="$(basename $file .OR.sorted.bam)"
  inStrain profile $file mag_3300020490_bin6-contigs.fa  \
      -o instrain/profile/${NAME}.IS -p 6 \
      -g mag_3300020490_bin6-contigs.fa.genes.fna
done


# compare

inStrain compare -i profile/* -s ../mag_3300020490_bin6-contigs.fa

# plot compare

inStrain plot -i instrainComparer/ -p10 
