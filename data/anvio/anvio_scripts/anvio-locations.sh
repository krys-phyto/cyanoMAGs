# Krystyn Kibler 2023-01-16
# Purpose find where the genes are on contigs to then compare to where snvs are located on contigs
# then can calculate snv density in regions

### Gene locations ###
# https://anvio.org/help/main/programs/anvi-export-gene-calls/
# where genes are on contigs
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-locations
anvi-export-gene-calls -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                       --gene-caller prodigal \
                       --skip-sequence-reporting \
                       -o gene-calls.txt

# what scg gene calls are
# https://merenlab.org/2019/10/08/anvio-scg-taxonomy/
anvi-estimate-scg-taxonomy -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                      --debug --just-do-it --output-file scg-gene-ids.txt

# Get NRPS regions from antismash



# confirm that gene ids are the same between gene-calls.txt, scg-gene-ids.txt, and 16cyanoMags-refined_FUNCTIONS.txt
#snvs / length of gene = snv density per that 'gene type'

#hopefully, nrps regions will have a greater snv density than say scg genes

#observation to find and basically the whole meat of my paper that ill write and be done
#with this project would be

#nrps snv density > regular genes for shit > scg genes

#maybe therell be significant differences between different groups of cyanos
