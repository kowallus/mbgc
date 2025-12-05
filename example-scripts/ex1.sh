#Listeria monocytogenes bacterial genomes were downloaded from US National Center for Biotechnology Information:
#https://www.ncbi.nlm.nih.gov/pathogens

echo "1. compress small collection of genomes"
mbgc c seqlist.txt archive1.mbgc

echo
echo "2. decompress to 'out' folder with original EOLs (after every 80 chars of DNA)"
mbgc d archive1.mbgc out

echo
echo "3. validation"
cmp GCA_000585755.1_Lm1823_genomic.fna out/GCA_000585755.1_Lm1823_genomic.fna
cmp GCA_000585775.1_Lm1824_genomic.fna out/GCA_000585775.1_Lm1824_genomic.fna
cmp GCA_000585795.1_Lm1840_genomic.fna out/GCA_000585795.1_Lm1840_genomic.fna

echo
echo "4. decompress to 'out' folder a file whose name contains 1824 substring,"
echo "  without EOLs in DNA sequences (overwrites files extracted in step 2.)"
mbgc d -f -l 0 -e 1824 archive1.mbgc out

echo
echo "5. validation (differences expected due to lack of EOLs)"
cmp GCA_000585775.1_Lm1824_genomic.fna out/GCA_000585775.1_Lm1824_genomic.fna

echo
echo "6. cleanup"
rm -R out
rm archive1.mbgc
