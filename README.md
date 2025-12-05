# MBGC: Multiple Bacteria Genome Compressor

[![GitHub downloads](https://img.shields.io/github/downloads/kowallus/mbgc/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/kowallus/mbgc/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/mbgc.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/mbgc)

Multiple Bacteria Genome Compressor (MBGC) is a tool for compressing
genomes in FASTA (or gzipped FASTA) input format.
It performs efficiently in terms of compression ratio and speed for various collections 
but is tailored for bacteria species. The implementation obtained >6.6 GB/s (both compression- and decompression-wise)
squeezing ~4.9GB collection of 1k *E. coli* genomes to ~4.5MB, tested in RAM disk on: 
Intel Core i9-10940X (14 cores) 3.3 GHz CPU, 128 GB of DDR4-RAM (2666 MHz, CL 16).

Major features:
* supports gzipped (.gz) archives as input and output,
* preserves folder structure of compressed files,
* preserves "well-formed" sequence formatting (defined in [NAF](https://kirill-kryukov.com/study/naf/) specification),
and offers [lossy compression](#lossless-and-lossy-compression-modes) when sequences are malformed,
* supports standard input and output during (de)compression,
* enables appending FASTA to existing archives,
* supports decompression (or repacking to new archive) of selected files,
* handles FASTA with non-bacteria species and amino acids as well.

*mbgc* is available through [bioconda repository](https://anaconda.org/bioconda/mbgc). More information can be found in the [bioconda page](https://bioconda.github.io/recipes/mbgc/README.html). 


### Manual build

The following steps create *mbgc* executable.
*mbgc* build requires cmake version >= 3.5 installed 
(check using ```cmake --version```).
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
Usage for multiple file compression (list of files given as input):
        mbgc c [-m <compressionMode>] [-L] <fastaFileList> <archiveFile>
Usage for single file compression:
        mbgc c [-m <compressionMode>] [-L] -i <inputFastaFile> <archiveFile>
Usage for decompression:
        mbgc d [-z <gzLevel>] <archiveFile> [<outputPath>]
Usage for partial decompression (list of patterns given as input):
        mbgc d [-E <patternsFile>] <archiveFile> [<outputPath>]

<fastaFileList> name of text file with a list of FASTA files (raw or gz)
        (given in separate lines) for compression
<inputFastaFile> name of a FASTA file (raw or gz) for compression
<archiveFile> mbgc archive filename
<patternsFile> name of text file with list of patterns (in separate lines)
        excludes files not matching any pattern (does not invalidate -e option)
<outputPath> extraction target path root (current directory by default)

Basic options (for compression, decompression and commons):
        [-m <compressionMode>] (speed: 0; default: 1; repo: 2; max: 3)
        [-L] allow lossy compression
        [-z <gzLevel>] extract FASTA files to gz archives
                (compression level: 1 <= z <= 12, recommended: 3)
        [-l <basesPerRow>] custom format of decoded sequences (0 - unlimited)
        [-f] overwrite an existing output file
        [-e <pattern>] exclude files with names not containing pattern
        [-E <patternsFile>] exclude files not matching any pattern 
        [-t <noOfThreads>] set limit of used threads
        [-I] ignore FASTA files paths (use only filenames)
        [-h] print full command help and exit
        [-v] print version number and exit

Compression modes description:
        (0) speed - for speed (fastest compression and decompression)
        (1) default - regular mode (good ratio, fast)
        (2) repo - for public repositories (better ratio, good speed)
        (3) max - for long-term storage (best ratio, memory-frugal)
```
compression of FASTA files (raw or gzipped) listed in *fastaList.txt* file 
(one FASTA file per line):
```
./mbgc c fastaList.txt comp.mbgc
```
compression of a single FASTA file (in FASTA or gzipped FASTA format):
```
./mbgc c -i input.fasta comp.mbgc
```
decompression to *out* folder (which is created if it does not exist):
```
./mbgc d comp.mbgc out
```
decompression to gz archives of files containing at least one pattern 
specified in *patterns.txt* file (one pattern per line) 
to *out* folder (which is created if it does not exist):
```
./mbgc d -z2 -E patterns.txt comp.mbgc out
```
Please note that since version 2.1 MBGC does not overwrite existing files!

To force overwriting use ```-f``` option, for example:
```
./mbgc d -f -comp.mbgc out
```

Exemplary data and scripts demonstrating usages of 
MBGC in basic compression scenarios are located in
[example-scripts](example-scripts) folder.
More test datasets and scripts are available [here](https://coach.iis.p.lodz.pl/mbgc2-datasets/).


### Lossless and lossy compression modes

Since version 2.1 MBGC ensures lossless compression by default.
Lossless compression is supported for multi-FASTA files:
- with headers starting with the '>' symbol,
- with a fixed line length for sequence lines within each file (except that the last line of a sequence may be shorter),
- without empty lines (yet, empty sequences are supported),
- without comments starting with the ';' symbol.

If a FASTA file cannot be compressed losslessly, lossy compression mode can be enabled using ```-L``` option.
```
./mbgc c -L fastaList.txt comp.mbgc
```
When sequences are not "well-formed" (see [NAF](https://kirill-kryukov.com/study/naf/) specification)
and line-length cannot be preserved losslessly,
MBGC will store the length of the longest sequence line within each file.
For lossy single-file compression this may differ, because the file is compressed in blocks 
and each block stores the length of the longest sequence line within that block. 


### Basic usage with standard I/O

Following POSIX convention, 
a single hyphen character can be used to specify input from
or output to the standard input and output streams.
```
for standard input set <inputFastaFile> to -
for standard input (resp. output) 
        in compression (resp. decompression) set <archiveFile> to -
for standard output set <outputPath> to - (all files are concatenated)
```
compression of FASTA in standard input data stream 
(in raw or gzipped FASTA format):
```
./mbgc c -i - comp.mbgc
```
decompression to standard output without EOLs symbols within sequences:
```
./mbgc d -l0 comp.mbgc -
```
### Usage with other commands

MBGC offers following commands:
```
Available commands (i - default):
        c       compress FASTA file(/s) into archive
        d       decompress FASTA file(/s) from archive
        i       info about contents (FASTA file names & headers) of archive
        a       append FASTA file(/s) to the given archive
        r       repack selected FASTA files from existing archive to new archive
```
listing filenames in given archive:
```
./mbgc i comp.mbgc
```
or using default command syntax: 
```
./mbgc comp.mbgc
```
listing headers (using convention: ">sequencename>filename") 
in filenames containing ASM17 pattern:
```
./mbgc i -e ASM17 comp.mbgc
```
appending FASTA files (raw or gzipped) listed in fastaList.txt file 
(one FASTA file per line) to archive:
```
./mbgc a fastalist.txt comp.mbgc
```
Please note that appending ignores FASTA with file names 
already existing in the archive.

repacking archive (in current or older mbgc version) to max compression mode:
```
./mbgc r -m3 comp.mbgc max.mbgc
```
repacking archive (using default compression mode) 
skipping files not matching any pattern specified in *patterns.txt* file 
(one pattern per line):
```
./mbgc r -E patterns.txt comp.mbgc part.mbgc
```
## Additional remarks

* FASTA files are compressed and stored by *mbgc* (and later extracted)
  in the order defined by *fastaFileList*.
* Reducing number of threads below 4 may result in ratio improvement (in a single-threaded runs
  from 1% up to 24% for a cross-pathogen dataset in the *speed* mode) at the expense of serious
  performance penalty (up to a factor of ~6 for human genome collections).
  The impact of threads on the ratio does not apply to *max* mode and is weakest in *repo* mode
  (up to ~6% for human genomes collections)
* Multithreaded compression ratio (mainly in *speed* and *default* modes) may vary up to ~1% from run to run.

## Publications

[Szymon Grabowski, Tomasz M. Kowalski: MBGC: Multiple Bacteria Genome Compressor (2022). *GigaScience*, Volume 11, 2022, giab099](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giab099/6515740) (concerns first version of MBGC)
