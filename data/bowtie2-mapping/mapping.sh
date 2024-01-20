# kjkibler 2022-Nov-21
# Generating new mapping files with bowtie2 for both anvio and instrain analyses

### Files ###
# concatenated mags
/home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags.fna
# mapping fastq files
/home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/fastq/*.OR.fastq


### Bowtie2 ###
#install current version of bowtie2
conda install bowtie2
# update bt2
conda update bowtie2

# version
#/home/glbrc.org/kjkibler/miniconda3/bin/bowtie2-align-s version 2.5.0
#64-bit
#Built on fv-az123-980
#Tue Nov  1 03:50:31 UTC 2022
#Compiler: gcc version 10.4.0 (conda-forge gcc 10.4.0-19)
#Options: -O3 -msse2 -funroll-loops -g3 -fvisibility-inlines-hidden -std=c++17 -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/glbrc.org/kjkibler/miniconda3/include -fdebug-prefix-map=/opt/conda/conda-bld/bowtie2_1667273633358/work=/usr/local/src/conda/bowtie2-2.5.0 -fdebug-prefix-map=/home/glbrc.org/kjkibler/miniconda3=/usr/local/src/conda-prefix -std=c++11 -DPOPCNT_CAPABILITY -DNO_SPINLOCK -DWITH_QUEUELOCK=1 -DWITH_ZSTD
#Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}

# mv to working directory
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping

# generate index
mkdir bt2

bowtie2-build /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/16cyanoMags.fna \
 bt2/16cyanoMags.fa

# change working dir and map with bt2
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/fastq/
for f in *.fastq
do
  echo $f
  bowtie2 -p 6 -x /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2/16cyanoMags.fa --very-sensitive \
  -q $f > /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2/${f/.fastq/.sam}
done

# convert sam to bam
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/bt2
for f in *.sam
do
    echo $f
    samtools view -S -b $f > ${f/.sam/.bam}
done
