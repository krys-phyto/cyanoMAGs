# kjkibler 2022-Dec-02
# Determine taxonomy using anvio's various programs and functions


### --- SCG Taxonomy --- ###

# Files
# contig.db
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db

# profile.db
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db

### SCG Taxonomy and bonus coverages ###

# code
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/
mkdir anvio-taxonomy
anvi-run-scg-taxonomy -c ../anvi-metagenomics/16cyanoMags.db --min-percent-identity 95

anvi-estimate-scg-taxonomy -c ../anvi-metagenomics/16cyanoMags.db \
                           -p ../anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                           -C LM_16cyanoMags \
                           -o SCG-TAXONOMY.txt

# quick coverage
mkdir ../anvio-coverages
anvi-estimate-scg-taxonomy -c ../anvi-metagenomics/16cyanoMags.db \
                           -p ../anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                           -C LM_16cyanoMags \
                           --compute-scg-coverages -o SCG-coverages.txt
mv SCG-coverages.txt ../anvio-coverages/SCG-coverages.txt


### Phylogenomics ###

# Files
#refined bins
# /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/*-contigs.fa


# code
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-taxonomy
mkdir refined-contigdb/
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags

# generate a contigdb for every refined 16cyanomag genome
for file in *-contigs.fa
do
  anvi-gen-contigs-database -f $file \
      -o /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-taxonomy/refined-contigdb/${file/-contigs.fa/.db} \
      -T 4
  anvi-run-hmms -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-taxonomy/refined-contigdb/${file/-contigs.fa/.db}
done

# generate a contigdb for every reference genome
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/mags/reference-genomes
for file in *.fna
do
  anvi-script-reformat-fasta $file -o ${file/.fna/-fixed.fna} -l 0 --simplify-names
  mv ${file/.fna/-fixed.fna} $file
  anvi-gen-contigs-database -f $file \
      -o /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-taxonomy/refined-contigdb/${file/_genomic.fna/.db} \
      -T 4
  anvi-run-hmms -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-taxonomy/refined-contigdb/${file/_genomic.fna/.db}
done


# create external-genomes.txt
cd /home/glbrc.org//kjkibler/paper-1-cyanomags/data/anvio/anvio-taxonomy/refined-contigdb
for file in GC*.db
do
        filename=`basename $file | sed 's/\.1/_1/g'`
        mv $file $filename
done

echo "name contigs_db_path" > external-genomes.txt
for file in *.db
do
        echo "${file%%.db} $file" >> external-genomes.txt
done
cat external-genomes.txt | tr [:blank:] \\t > external-genomes1.txt
mv external-genomes1.txt external-genomes.txt

# get shared hmms
anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
      --hmm-source Bacteria_71 \
      --list-available-gene-names

# get list of shared genes from terminal and do this to create a txt file
genenames="Ribosomal_L1, Ribosomal_L13, Ribosomal_L14,Ribosomal_L16, Ribosomal_L17, Ribosomal_L18p, Ribosomal_L19, Ribosomal_L2,Ribosomal_L20, Ribosomal_L21p, Ribosomal_L22, Ribosomal_L23, Ribosomal_L27,Ribosomal_L27A, Ribosomal_L28, Ribosomal_L29, Ribosomal_L3, Ribosomal_L32p,Ribosomal_L35p, Ribosomal_L4, Ribosomal_L5, Ribosomal_L6, Ribosomal_L9_C,Ribosomal_S10, Ribosomal_S11, Ribosomal_S13, Ribosomal_S15, Ribosomal_S16,Ribosomal_S17, Ribosomal_S19, Ribosomal_S2, Ribosomal_S20p, Ribosomal_S3_C, Ribosomal_S6, Ribosomal_S7, Ribosomal_S8, Ribosomal_S9"
echo $genenames | awk -F, -v OFS="\n" '{$1=$1; print}' | sed 's/ //g' > gene-names.txt

# get all the sequences of the shared hmms
anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
      -o ../concatenated-proteins.fa \
      --hmm-source Bacteria_71 \
      --gene-names  gene-names.txt \
      --return-best-hit \
      --get-aa-sequences \
      --concatenate

# rename generate concatenated protein file to something else informative

# generate tree
cd ../
anvi-gen-phylogenomic-tree -f concatenated-proteins_ref-16cyanoMags.fna \
                           -o concatenated-proteins_ref-16cyanoMags.fna-tree.txt


anvi-interactive -p phylogenomic-profile.db \
                 -t phylogenomic-tree.txt \
                 --title "Phylogenomics Tutorial Example #1" \
                 --manual



### Bonus Code ###
# split taxonomy
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-taxonomy

anvi-export-splits-taxonomy -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                            -o splits-taxonomy-txt
 # nevermind it basically requires running kaiju, importing that information in, then you can export it again lol


### Phycocyanin Maybe ###
cd /home/glbrc.org//kjkibler/paper-1-cyanomags/data/anvio/anvio-taxonomy
anvi-export-locus -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                   --num-genes 2 \
                   -o Cpc \
                   -O Phycocyanin \
                   --search-term “Phycocyanin/phycoerythrin beta chain, CpcB/CpeB”
