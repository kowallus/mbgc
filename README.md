# MBGC: Multiple Bacteria Genome Compressor

Multiple Bacteria Genome Compressor (MBGC) is a tool for compressing
genomes in FASTA (or gzipped FASTA) input format.
It is tailored for fast and efficient compression of bacteria species collections.

The implementation:
* supports gzipped (.gz) archives as input,
* preserves folder structure of compressed files,
* decompresses DNA streams without EOLs or with a fixed number of bases per line,
* handles FASTA with non-bacteria species as well,
* supports standard input and output during (de)compression.

### Installation on Linux
mbgc requires [libdeflate](https://github.com/ebiggers/libdeflate) library to work (an example
[install howto](https://pkgs.org/search/?q=libdeflate-dev) for Debian and Ubuntu).

The following steps create an *mbgc* executable.
On Linux *mbgc* build requires cmake version >= 3.4 installed (check using ```cmake --version```):
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
Usage for compression using files list as input: mbgc [-t noOfThreads] <sequencesListFile> <archiveFile>
Usage for single file compression: mbgc [-t noOfThreads] -i <inputFastaFile> <archiveFile>

Usage for decompression: mbgc -d [-t noOfThreads] [-f pattern] [-l dnaLineLength] <archiveFile> [<outputPath>]

-t number of threads used (default: 8)
-d decompression mode
-l format decompressed DNA (i.e., sets the number of bases per row)
-f decompress files with names containing the given pattern

<sequencesListFile> name of text file containing a list of FASTA files (raw or in gz archives)
        (given in separate lines) for compression
<inputFastaFile> name of a FASTA file (raw or in gz archive) for compression
<archiveFile> mbgc archive filename
<outputPath> extraction target path root (if skipped the root path is the current directory)
```

compression of FASTA files (raw or gzipped) listed in *seqlist.txt* file (one FASTA file per line):
```
./mbgc seqlist.txt comp.mbgc
```
compression of a single FASTA file (in FASTA or gzipped FASTA format):
```
./mbgc -i input.fasta comp.mbgc
```
decompression to *out* folder (which is created if it does not exist) using 80 bases per row DNA formatting:
```
./mbgc -l 80 -d comp.mbgc out
```
Please note that decompression overwrites existing files!

Exemplary data and scripts demonstrating usages of MBGC in basic compression scenarios are located in
[example-scripts](example-scripts) folder.

### Basic usage with standard I/O

Following POSIX convention, a single hyphen character can be used to specify input from
or output to the standard input and output streams.

```
for standard input set <inputFastaFile> to -
for standard input (resp. output) in compression (resp. decompression) set <archiveFile> to -
for standard output set <outputPath> to - (all files are concatenated)
```

compression of FASTA in standard input data stream (in raw or gzipped FASTA format):
```
./mbgc -i - comp.mbgc
```
decompression to standard output (without EOLs symbols within DNA sequences):
```
./mbgc -d comp.mbgc -
```

<!--
## Publications
[Szymon Grabowski, Tomasz M. Kowalski: MBGC: Multiple Bacteria Genome Compressor (2021).]()

[supplementary data]()
-->
