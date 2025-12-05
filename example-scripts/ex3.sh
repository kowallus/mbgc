#Listeria monocytogenes bacterial genomes were downloaded from US National Center for Biotechnology Information:
#https://www.ncbi.nlm.nih.gov/pathogens

echo "1. compress small collection of genomes"
mbgc c seqlist.txt archive3.mbgc

echo
echo "2. repack archive in max mode with files containing 182 substring"
mbgc r -m3 -e 182 archive3.mbgc archive3_max.mbgc

echo
echo "3. decompress to 'out' folder with original EOLs (after every 80 chars of DNA)"
mbgc d archive3_max.mbgc out

echo
echo "4. validation (differences expected due to lack of one file)"
cmp GCA_000585755.1_Lm1823_genomic.fna out/GCA_000585755.1_Lm1823_genomic.fna
cmp GCA_000585775.1_Lm1824_genomic.fna out/GCA_000585775.1_Lm1824_genomic.fna
cmp GCA_000585795.1_Lm1840_genomic.fna out/GCA_000585795.1_Lm1840_genomic.fna

echo
echo "5. append archive (files with duplicate names are skipped)"
mbgc a seqlist.txt archive3_max.mbgc

echo
echo "6. decompress to 'out' folder with original EOLs (after every 80 chars of DNA)"
echo "  (overwrites files extracted in step 3.)"
mbgc d -f archive3_max.mbgc out

echo
echo "7. validation"
cmp GCA_000585755.1_Lm1823_genomic.fna out/GCA_000585755.1_Lm1823_genomic.fna
cmp GCA_000585775.1_Lm1824_genomic.fna out/GCA_000585775.1_Lm1824_genomic.fna
cmp GCA_000585795.1_Lm1840_genomic.fna out/GCA_000585795.1_Lm1840_genomic.fna

echo
echo "8. cleanup"
rm -R out
rm archive3.mbgc
rm archive3_max.mbgc
