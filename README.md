# kbbq

K-mer based base quality score recalibration.

## downloading, building, and installing

```
git clone https://github.com/adamjorr/kbbq.git
cd kbbq/
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cmake --install .
```
Use the `--prefix` option to change where CMake installs the executable if you don't have permission to use the default:
```
cmake --install . --prefix ~/.local/
```
will put the executable in `~/.local/bin/`, a common place for such things. Alternatively, move the executable to somewhere on your path or call it from the build directory. 

The only dependency is [HTSLib](https://github.com/samtools/htslib/). Currently `kbbq` requires POSIX-compliance for posix_memalign and getopt_long. This means you may struggle to compile on Windows systems; compiling on Cygwin for Windows should be OK.

It also requires a processor that supports AVX1 instructions. If you run into a problem compiling because your processor doesn't support AVX, open an issue and I will write a workaround.

## quickstart

`kbbq --threads 6 --genomelen NUM_BASEPAIRS --coverage COVERAGE INPUT.fq > RECALIBRATED.fq.gz`

For example, use `--genomelen 1000000 --coverage 20` for a 1 megabase region sequenced to a depth of 20X. If the coverage is not provided, `kbbq` will take an extra pass through the data to estimate it. If the genome length is not provided, it can be estimated using the headers in a BAM file, but it will fail with FASTQ input.

`kbbq` accepts reads in SAM/BAM or FASTQ format, gzipped or not. It will automatically detect the format of the input. If the input is SAM/BAM, none of the alignment information is used; it is merely supported for convenience. The output will be gzipped.

## options

Parameter | Short Option | Default Value | Summary
--- | --- | --- | ---
`--ksize` | `-k` | 32 | Size of k-mer to use for correction
`--use-oq` | `-u` | Off | Use BAM OQ tag values as quality scores
`--set-oq` | `-s` | Off | Set BAM OQ tag values before recalibration
`--genomelen` | `-g` | Estimated for BAM, required for FASTQ | The approximate size of the sequenced region in base-pairs.
`--coverage` | `-c` | Estimated from data | Approximate sequencing coverage
`--fixed` | `-f` | Off | Treat changes to reads in the given file as errors and recalibrate. 
`--alpha` | `-a` | 7 / coverage | Rate to sample k-mers
`--threads` | `-t` | 1 | Number of CPU threads to use

## details

`kbbq` uses an error correction method similar to [Lighter](www.github.com/mourisl/Lighter) to find errors in whole-genome sequencing reads, then applies a hierarchical model similar to GATK's BaseRecalibrator to recalibrate quality scores.
