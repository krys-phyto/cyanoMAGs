# kjkibler 2022-Dec-06
# Export coverages

### Export coverages ###

# Files
# contig.db
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db

# profile.db
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db


### code ###
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-coverages

anvi-export-splits-and-coverages -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                                 -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                                 -O 16cyanoMags-refined


anvi-export-splits-and-coverages -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                                 -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                                 -O 16cyanoMags-refined \
                                 --use-Q2Q3-coverages

# mv things back to anvio-coverages/ vs in the profile directory *eye roll*


anvi-export-gene-coverage-and-detection -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                                        -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                                        -O 16cyanoMags-refined


anvi-export-functions -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                      -o 16cyanoMags-refined_FUNCTIONS.txt
