#Listeria monocytogenes bacterial genomes were downloaded from US National Center for Biotechnology Information:
#https://www.ncbi.nlm.nih.gov/pathogens

echo "1. compress small collection of genomes from a single file"
mbgc c -i GCA_lm_concat.fna archive2.mbgc

echo
echo "2. decompress to 'out' folder with original EOLs (after every 80 chars of DNA)"
mbgc d archive2.mbgc out

echo
echo "3. validation"
cmp GCA_lm_concat.fna out/GCA_lm_concat.fna

echo
echo "4. decompress to 'out' folder without EOLs in DNA sequences"
echo "  (overwrites files extracted in step 2.)"
mbgc d -f -l 0 archive2.mbgc out

echo
echo "5. validation (differences expected due to lack of EOLs)"
cmp GCA_lm_concat.fna out/GCA_lm_concat.fna

echo
echo "6. cleanup"
rm -R out
rm archive2.mbgc
