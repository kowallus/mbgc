# MBGC: Multiple Bacteria Genome Compressor
 
Multiple Bacteria Genome Compressor (MBGC) is a tool for compressing 
genomes in FASTA (or gzipped FASTA) input format. 
It is tailored for fast and efficient compression of bacteria species collections.  

The implementation:
* supports gz archives as input, 
* preserves folder structure of compressed files, and 
* decompresses DNA streams 80 bases per line,
* handles FASTA with non-bacteria species as well.
### Installation on Linux
mbgc requires [libdeflate](https://github.com/ebiggers/libdeflate) library to work (an example 
[install howto](https://pkgs.org/search/?q=libdeflate-dev) for Debian and Ubuntu).

The following steps create an *mbgc* executable. 
On Linux *mbgc* build requires installed cmake version >= 3.4 (check using ```cmake --version```):
```bash
git clone https://github.com/kowallus/mbgc.git
cd mbgc
mkdir build
cd build
cmake ..
make mbgc
```

### Basic usage

```
Usage for compression: mbgc [-c level] [-t noOfThreads] [-m] <sequencesListFile> <archiveFile>

Usage for decompression: mbgc -d [-t noOfThreads] <archiveFile> [<outputPath>]

-c compression level (1 - fast; 2 - default; 3 - max)
-t number of threads used (8 - default)
-m compression of mixed species collection
```

compression of fasta files listed in *samplelist.txt* file (one fasta or gz per line):
```
./mbgc samplelist.txt -o comp.mbgc
```
decompression of fasta files to *out* folder (which is created if it does not exist):
```
./mbgc -d comp.mbgc out
```
Please note, that decompression overwrites existing files.

<!--
## Publications
[Szymon Grabowski, Tomasz M. Kowalski: MBGC: Multiple Bacteria Genome Compressor (2021).]()

[supplementary data]()
-->
