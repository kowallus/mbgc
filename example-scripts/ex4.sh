#Listeria monocytogenes bacterial genomes were downloaded from US National Center for Biotechnology Information:
#https://www.ncbi.nlm.nih.gov/pathogens

echo "1. compress a single genome"
mbgc c -i GCA_000585755.1_Lm1823_genomic.fna archive4.mbgc

echo
echo "2. decompress to 'out' folder with original EOLs (after every 80 chars of DNA)"
mbgc d archive4.mbgc out

echo
echo "3. validation"
cmp GCA_000585755.1_Lm1823_genomic.fna out/GCA_000585755.1_Lm1823_genomic.fna

echo
echo "4. append archive with contents of remaining two genomes"
mbgc a -i GCA_000585775.1_Lm1824_genomic.fna archive4.mbgc
mbgc a -i GCA_000585795.1_Lm1840_genomic.fna archive4.mbgc

echo
echo "5. decompress to 'out' folder with original EOLs (after every 80 chars of DNA)"
mbgc d -f archive4.mbgc out

echo
echo "6. validation (against concatenated multi-FASTA)"
cmp GCA_lm_concat.fna out/GCA_000585755.1_Lm1823_genomic.fna

echo
echo "7. cleanup"
rm -R out
rm archive4.mbgc
