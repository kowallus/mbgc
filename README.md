# MBGC: Multiple Bacteria Genome Compressor
 
Multiple Bacteria Genome Compressor (MBGC) is a tool for compressing 
genomes in FASTA (or gzipped FASTA) input format. 
It is tailored for fast and efficient compression of bacteria species collections.  

The implementation:
* supports gz archives as input, 
* preserves folder structure of compressed files,
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
Usage for compression: mbgc [-t noOfThreads] <sequencesListFile> <archiveFile>

Usage for decompression: mbgc -d [-t noOfThreads] [-f pattern] <archiveFile> [<outputPath>]

-t number of threads used (8 - default)
-d decompression mode
-f decompress files with names containing the given pattern
```

compression of fasta files listed in *samplelist.txt* file (one fasta or gz per line):
```
./mbgc samplelist.txt comp.mbgc
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
