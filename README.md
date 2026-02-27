# **BIU**

A set of utilities for processing biological data: FASTA transformation, BLAST parsing, FASTQ filtering, and nucleotide sequence tools.




### **Installation and requirments**

Copy the folder and work with `biu`. Python 3.x (tested on Python 3.8 and above).





###  **`main.py`** 


1. **`run_dna_rna_tools`** takes as input a DNA or RNA sequence (str), as well as the name of the procedure to be performed ('is_nucleic_acid', 'transcribe', 'reverse', 'complement', 'reverse_complement')

    - `is_nucleic_acid` Validation RNA or DNA. Return bool result.
    - `transcribe` Make a T to U substitution. Returns the transcribed sequence.
    - `reverse` Return the reversed (mirror) sequence. 
    - `complement` Return the complementary sequence. 
    - `reverse_complement` Return the reverse complementary sequence. 

    *example of usage: run_dna_rna_tools('TTUU', 'is_nucleic_acid')*


2. **`filter_fastq`** uses *.fastq and *.fasta files and save them into filtered data ./filtered.


### **`bio_files_processor.py`**

1) FASTA: Convert multiline format to single line per sequence.
2) BLAST: Extract unique annotations from text report.
3) FASTQ: Filter reads by GC content, length, and quality.


### **Autor**

Loginova Olga




### **Feedback and bug reports**

If you have any troubles running tools, please attach params.txt
