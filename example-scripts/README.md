# MBGC usage examples

The data used in the examples are Listeria monocytogenes bacterial genomes in FASTA format:
* GCA_000585755.1_Lm1823_genomic.fna
* GCA_000585775.1_Lm1824_genomic.fna
* GCA_000585795.1_Lm1840_genomic.fna

The genomes were downloaded from [National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/pathogens).

The folder contains two scripts performing two basic scenarios of MBGC usage:
* [ex1.sh](ex1.sh) - (de)compression of 3 FASTA files representing a collection of genomes
* [ex2.sh](ex2.sh) - (de)compression of a small collection of genomes stored in a single file

The list of FASTA files used in `ex2.sh` is in:
* [seqlist.txt](seqlist.txt)

The simple concatenation of the Listeria monocytogenes bacterial genomes used in `ex2.sh` is in:
* GCA_lm_concat.fna

### Executing the examples on Linux

For the examples to work, the `mbgc` binary needs to be located in a programs' folder 
included in the system path, e.g. `/usr/local/bin`.
```
cp mbgc /usr/local/bin 
```

Then you can start the tests by running:
```
bash ex1.sh
bash ex2.sh
```
