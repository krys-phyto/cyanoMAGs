# kjkibler 2022-12-07
# Analyzing SNVs and SCAAVs with anvio variability programs

# Files #
# contig.db
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db

# profile.db
#/home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db

# Collection name
#LM_16cyanoMags

### Variability Profiles for SNVs and SCAAVs ###
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-snvs


anvi-script-add-default-collection -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                             -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db


anvi-gen-variability-profile -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                             -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --quince-mode


anvi-gen-variability-profile -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                             -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --engine AA \
                             --quince-mode


anvi-gen-variability-profile -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                             -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --engine CDN \
                             --quince-mode



### Variability profiles for selected genes
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-snvs

    # genes-of-interest file: scg-genes-of-interest.txt
# text file of scg genes
anvi-gen-variability-profile  -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                              -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                              --genes-of-interest scg-genes-of-interest.txt \
                              --quince-mode \
                              --include-contig-names \
                              --compute-gene-coverage-stats \
                              --quince-mode \
                              --output-file scg-snv-variability-profile-NUC.txt

  # strings-of-interest file: nrps_contigs-of-interest.txt
sed -i 's/$/_split_00001/' nrps_contigs-of-interest.txt
sed -i '1d' nrps_contigs-of-interest.txt
sed -i 's/mag_//' nrps_contigs-of-interest.txt



anvi-gen-variability-profile  -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                              -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                              --splits-of-interest nrps_contigs-of-interest.txt \
                              --quince-mode \
                              --include-contig-names \
                              --compute-gene-coverage-stats \
                              --output-file nrps-snv-variability-profile-NUC.txt


# Selected genes -- gene keys from paper-1-cyanoMags/data/anvio/anvio-snvs/gene-regions_keys.txt
anvi-gen-variability-profile  -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                              -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                              --genes-of-interest gene-regions_keys.txt \
                              --quince-mode \
                              --include-contig-names \
                              --compute-gene-coverage-stats \
                              --output-file gene-regions_snv-variability-profile_NUC.txt

Total number of variable positions in samples .................: metag_20080626: 19034;
                                                                 metag_20080703: 132070;
                                                                 metag_20080709: 27090;
                                                                 metag_20080718: 93536;
                                                                 metag_20080719: 8896;
                                                                 metag_20080721: 45190;
                                                                 metag_20080723: 22976;
                                                                 metag_20080805: 9079;
                                                                 metag_20080813: 102651;
                                                                 metag_20080820: 3830;
                                                                 metag_20080827: 17437;
                                                                 metag_20080912: 89568;
                                                                 metag_20080913: 14287;
                                                                 metag_20080925: 7290;
                                                                 metag_20081008: 48491;
                                                                 metag_20081017: 31235;
                                                                 metag_20090422: 98499;
                                                                 metag_20090429: 9468;
                                                                 metag_20090609: 14370;
                                                                 metag_20090618: 22336;
                                                                 metag_20090626: 23814;
                                                                 metag_20090707: 51889;
                                                                 metag_20090730: 2115;
                                                                 metag_20090810: 109202;
                                                                 metag_20090826: 29934;
                                                                 metag_20090913: 10237;
                                                                 metag_20090914: 0;
                                                                 metag_20090927: 188;
                                                                 metag_20091007: 19546;
                                                                 metag_20091026: 1713;
                                                                 metag_20091114: 188;
                                                                 metag_20100420: 6559;
                                                                 metag_20100505: 13880;
                                                                 metag_20100518: 303;
                                                                 metag_20100520: 179721;
                                                                 metag_20100602: 3756;
                                                                 metag_20100613: 60662;
                                                                 metag_20100615: 12;
                                                                 metag_20100621: 12;
                                                                 metag_20100706: 19;
                                                                 metag_20100715: 20638;
                                                                 metag_20100716: 151640;
                                                                 metag_20100727: 3827;
                                                                 metag_20100805: 0;
                                                                 metag_20100817: 25366;
                                                                 metag_20100830: 45615;
                                                                 metag_20100831: 5494;
                                                                 metag_20100914: 19797;
                                                                 metag_20100926: 98527;
                                                                 metag_20101013: 1;
                                                                 metag_20101029: 269;
                                                                 metag_20101119: 21876;
                                                                 metag_20110503: 154;
                                                                 metag_20110518: 1;
                                                                 metag_20110601: 39304;
                                                                 metag_20110613: 7385;
                                                                 metag_20110628: 8956;
                                                                 metag_20110712: 2323;
                                                                 metag_20110725: 290;
                                                                 metag_20110809: 23276;
                                                                 metag_20110822: 9957;
                                                                 metag_20110904: 50998;
                                                                 metag_20110921: 75990;
                                                                 metag_20111003: 15758;
                                                                 metag_20111101: 5212;
                                                                 metag_20111130: 26095;
                                                                 metag_20120305: 4132;
                                                                 metag_20120402: 49144;
                                                                 metag_20120505: 2444;
                                                                 metag_20120517: 5617;
                                                                 metag_20120602: 45336;
                                                                 metag_20120608: 717;
                                                                 metag_20120615: 7162;
                                                                 metag_20120622: 1299;
                                                                 metag_20120629: 73632;
                                                                 metag_20120706: 23849;
                                                                 metag_20120713: 34949;
                                                                 metag_20120717: 3848;
                                                                 metag_20120720: 170490;
                                                                 metag_20120803: 6433;
                                                                 metag_20120817: 2422;
                                                                 metag_20120824: 14894;
                                                                 metag_20120831: 10732;
                                                                 metag_20120907: 26542;
                                                                 metag_20120913: 3431;
                                                                 metag_20120921: 5898;
                                                                 metag_20120927: 13092;
                                                                 metag_20121008: 9103;
                                                                 metag_20121012: 45452;
                                                                 metag_20121022: 5364;
                                                                 metag_20121026: 113538;
                                                                 metag_20121105: 38960;
                                                                 metag_20121109: 4348;
                                                                 metag_20121116: 36394


# 2023-12-5
cd /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvio-snvs
anvi-gen-variability-profile -p /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/profile/16cyanoMags-MERGED/PROFILE.db \
                             -c /home/glbrc.org/kjkibler/paper-1-cyanomags/data/anvio/anvi-metagenomics/16cyanoMags.db \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --engine AA \
                             --quince-mode
