# biu

Utilities for work with DNA or RNA sequenses and FASTQ.

## Installation and Requirments

Just copy the folder and work with "main.py". It contains imports and two functions "run_dna_rna_tools" and "filer_fastq". There are no dependencies - only the Python standard library.

Python 3.x (tested on Python 3.8 and above).

## Usage

### 1. Work with DNA or RNA

"run_dna_rna_tools" takes as input a DNA or RNA sequence (str), as well as the name of the procedure to be performed ('is_nucleic_acid', 'transcribe', 'reverse', 'complement', 'reverse_complement')

### Validation RNA or DNA. Return bool result.

run_dna_rna_tools('TTUU', 'is_nucleic_acid') \# → False run_dna_rna_tools('TTAA', 'is_nucleic_acid') \# → TRUE

### Make a T to U substitution. Returns the transcribed sequence (DNA → RNA)

run_dna_rna_tools('ATGC', 'transcribe') \# → 'AUGC'

### Return the reversed (mirror) sequence.

run_dna_rna_tools('ATG', 'reverse') \# → 'GTA'

### Return the complementary sequence.

run_dna_rna_tools('AtG', 'complement') \# → 'TaC'

### Return the reverse complementary sequence.

run_dna_rna_tools('ATg', 'reverse_complement') \# → 'cAT'

###When you have two or more sequenses.

run_dna_rna_tools('ATG', 'aT', 'reverse') \# → ['GTA', 'Ta']

### 2. Work with FASTQ

"filer_fastq" takes 4 arguments: seqs, gc_bounds, length_bounds, quality_threshold)

seqs = dict.

gc_bounds=(0, 100), by default all reads are preserved. If a single number is passed as an argument, it is considered the upper limit.

length_bounds=(0, 2\*\*32), by default all reads are preserved. If a single number is passed as an argument, it is considered the upper limit.

quality_threshold=0.0, the average read quality threshold for filtering is 0 by default (phred33 scale). Reads with an average quality for all nucleotides below the threshold are discarded.

## Team

Loginova Olga and web :)

## Citation

If you use theese tools in your work, please cite [<https://github.com/jallopukki/biu.git>].

## Feedback and bug reports

If you have any troubles running tools, please attach params.txt
