# Krys Kibler
# 2023-12-07

# InStrain

cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/maer-instrain/
mkdir maer-instrain/
cd maer-instrain/

prodigal -i mag_3300020573_bin18-contigs.fa -d mag_3300020573_bin18-contigs.fa.genes.fna

# mapping
mkdir bt2

bowtie2-build mag_3300020573_bin18-contigs.fa \
  bt2/mag_3300020573_bin18-contigs.fa

for f in /home/glbrc.org/kjkibler/paper-1-cyanomags/data/bowtie2-mapping/fastq/*.fastq
do
  echo $f
  NAME="$(basename $f .fastq)"
  bowtie2 -p 6 -x bt2/mag_3300020573_bin18-contigs.fa --very-sensitive \
    -q $f > bt2/$NAME.sam
done

# profile

mkdir instrain
mkdir instrain/profile

for file in bt2/*.sam
do
  echo $file
  NAME="$(basename $file .OR.sorted.bam)"
  inStrain profile $file mag_3300020573_bin18-contigs.fa  \
      -o instrain/profile/${NAME}.IS -p 6 \
      -g mag_3300020573_bin18-contigs.fa.genes.fna
done


# compare

inStrain compare -i /profile/* -s ../mag_3300020573_bin18-contigs.fa

inStrain compare -i profile/* #these two profile broke the code #sadge#
#profile/20100615.OR.sam.IS
#profile/20100621.OR.sam.IS

inStrain compare -i profile/20080626.OR.sam.IS \
profile/20100505.OR.sam.IS \
profile/20111101.OR.sam.IS \
profile/20080703.OR.sam.IS \
profile/20100518.OR.sam.IS \
profile/20111130.OR.sam.IS \
profile/20080709.OR.sam.IS \
profile/20100520.OR.sam.IS \
profile/20120305.OR.sam.IS \
profile/20080718.OR.sam.IS \
profile/20100602.OR.sam.IS \
profile/20120402.OR.sam.IS \
profile/20080719.OR.sam.IS \
profile/20100613.OR.sam.IS \
profile/20120505.OR.sam.IS \
profile/20080721.OR.sam.IS \
profile/20120517.OR.sam.IS \
profile/20080723.OR.sam.IS \
profile/20120602.OR.sam.IS \
profile/20080805.OR.sam.IS \
profile/20100706.OR.sam.IS \
profile/20120608.OR.sam.IS \
profile/20080813.OR.sam.IS \
profile/20100715.OR.sam.IS \
profile/20120615.OR.sam.IS \
profile/20080820.OR.sam.IS \
profile/20100716.OR.sam.IS \
profile/20120622.OR.sam.IS \
profile/20080827.OR.sam.IS \
profile/20100727.OR.sam.IS \
profile/20120629.OR.sam.IS \
profile/20080912.OR.sam.IS \
profile/20100805.OR.sam.IS \
profile/20120706.OR.sam.IS \
profile/20080913.OR.sam.IS \
profile/20100817.OR.sam.IS \
profile/20120713.OR.sam.IS \
profile/20080925.OR.sam.IS \
profile/20100830.OR.sam.IS \
profile/20120717.OR.sam.IS \
profile/20081008.OR.sam.IS \
profile/20100831.OR.sam.IS \
profile/20120720.OR.sam.IS \
profile/20081017.OR.sam.IS \
profile/20100914.OR.sam.IS \
profile/20120803.OR.sam.IS \
profile/20090422.OR.sam.IS \
profile/20100926.OR.sam.IS \
profile/20120817.OR.sam.IS \
profile/20090429.OR.sam.IS \
profile/20101013.OR.sam.IS \
profile/20120824.OR.sam.IS \
profile/20090609.OR.sam.IS \
profile/20101029.OR.sam.IS \
profile/20120831.OR.sam.IS \
profile/20090618.OR.sam.IS \
profile/20101119.OR.sam.IS \
profile/20120907.OR.sam.IS \
profile/20090626.OR.sam.IS \
profile/20110503.OR.sam.IS \
profile/20120913.OR.sam.IS \
profile/20090707.OR.sam.IS \
profile/20110518.OR.sam.IS \
profile/20120921.OR.sam.IS \
profile/20090730.OR.sam.IS \
profile/20110601.OR.sam.IS \
profile/20120927.OR.sam.IS \
profile/20090810.OR.sam.IS \
profile/20110613.OR.sam.IS \
profile/20121008.OR.sam.IS \
profile/20090826.OR.sam.IS \
profile/20110628.OR.sam.IS \
profile/20121012.OR.sam.IS \
profile/20090913.OR.sam.IS \
profile/20110712.OR.sam.IS \
profile/20121022.OR.sam.IS \
profile/20090914.OR.sam.IS \
profile/20110725.OR.sam.IS \
profile/20121026.OR.sam.IS \
profile/20090927.OR.sam.IS \
profile/20110809.OR.sam.IS \
profile/20121105.OR.sam.IS \
profile/20091007.OR.sam.IS \
profile/20110822.OR.sam.IS \
profile/20121109.OR.sam.IS \
profile/20091026.OR.sam.IS \
profile/20110904.OR.sam.IS \
profile/20121116.OR.sam.IS \
profile/20091114.OR.sam.IS \
profile/20110921.OR.sam.IS \
profile/20100420.OR.sam.IS \
profile/20111003.OR.sam.IS -s ../mag_3300020573_bin18-contigs.fa

# plot compare

inStrain plot -i instrainComparer/ -p10
