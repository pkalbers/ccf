# CCF
Dynamic programming method to infer the cumulative coalescent function (CCF) from haplotype data and knowledge about allele age.

This document describes the usage of the dynamic programming technique used to infer the CCF between two genomic sequences, driven by estimates of allele age.

> **Dating genomic variants and shared ancestry in population-scale sequencing data**  
> *Patrick K. Albers and Gil McVean*  
> doi: https://doi.org/10.1101/416610
###### See manuscript on bioRxiv: https://www.biorxiv.org/content/early/2018/09/13/416610

The CCF algorithm is described in detail in the Supplementary Text.

### Compilation ...
... is straightforward, because the implementation is very simple. To compile the program, type
```
g++ -std=c++11 -O2 -o ccf ccf.cpp
```
on the command line, or execute the provided `./build.sh` script.
This generates a single executable called `ccf`.


## Input
The program is executed on the command line:
```
./ccf /path/to/INPUT.txt
```
where `/path/to/INPUT.txt` refers to an input file, `INPUT.txt`, which is created separately by the user.  
The description below demonstrates how such an input file is generated.

Consider two haplotype sequences (variant alleles; encoded as 0s and 1s); one from a given **target** individual, and one from another, **comparator** individual.
As usual, variation is initially sorted by variant position along the chromosomal sequence.
We only take the variants carried by the target genome (1s), and compare those to the same positions in the comparator, such that each variant can be labelled as either being *shared* (**S**) or *unshared* (**U**). Only the shared/unshared states are retained to form an observation sequence.
This sequence is then sorted by the age of the variants (indicated by *a, b, c, ..., k* in the example below); from youngest to oldest. The sorted sequence of states is the input of the CCF algorithm (for simplicity encoded as 0s and 1s, but not to be confused with 0/1-encoded alleles; see example below).

###### Example:
```
Position:     0 -------------------------------- 0.5 -------------------------------- 1

Age:*                f     b  h           d  c     g     j  a  i  k           e      
                     |     |  |           |  |     |     |  |  |  |           |
  Target:      0  0  1  0  1  1  0  0  0  1  1  0  1  0  1  1  1  1  0  0  0  1  0  0
  Comparator:  1  0  1  0  1  0  0  1  0  1  1  1  0  0  1  0  0  1  0  0  0  1  1  0

State:               S     S  U           S  S     U     S  U  U  S           S

* Age of variants:     a < b < c < d < e < f < g < h < i < j < k
Observation sequence:  U   S   S   S   S   S   U   U   U   S   S
   (Encoded sequence:  0   1   1   1   1   1   0   0   0   1   1 )
```

In the above, the variants that are shared or unshared between the two haplotypes are used to create an encoded sequence of states (bottom). This sequence is written into an input file (`INPUT.txt`), where each state is given on a new line, which is then supplied to run the `ccf` executable.  
The input file generated from the example above contains the following:
```
0
1
1
1
1
1
0
0
0
1
1
```
That's it.

## Output
The inferred path sequence is printed to the standard output. This can be captured, for example, by typing:
```
./ccf /path/to/INPUT.txt > OUTPUT.txt
```

The output (e.g. captured in `OUTPUT.txt`) is the inferred fraction of the genome shared at each position along the time-sorted input sequence.
To arrive at the fraction of the genome shared at different points in **time**, simply map the ages of the age-sorted variants to the output sequence of states.


## Comment
Fractions are inferred at a discrete number of "hidden" states. By default, this is set to 200, which is hardcoded.
To change this value, simply modify
```C++
#define NSTATES 200
```
in the `ccf.cpp` file before compilation.
