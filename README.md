---
editor_options: 
  markdown: 
    wrap: 72
---

# **BIU**

Tool for work with biological sequences: FASTA transformation, BLAST
parsing, FASTQ filtering, and nucleotide sequence tools.\
**biu**/\
├── main.py\
├── requirements.txt\
├── data/\
└── filtered/\

**Installation**

-   Clone repository `bash`\
    `git clone <repository_url>`\
    `cd biu`\
-   Create or/and activate virtual environment\
    `python3 -m venv venv`\
    `source venv/bin/activate`\
-   Install dependencies\
    `pip install -r requirements.txt`\
-   Run\
    `python main.py`\

**Usage**

==== For DNA ====\
`from main import DNASequence`\
dna = DNASequence("ATGCATGC")\
dna.complement() \# TACGTACG\
dna.reverse() \# CGTACGTA\
dna.reverse_complement() \# GCATGCAT\
dna.transcribe() \# RNASequence("AUGCAUGC")\
\
==== For RNA ====\
`from main import RNASequence`\
rna = RNASequence("AUGCAUGC")\
rna.complement() \# UACGUACG\
\
==== For Proteins ====\
`from main import AminoAcidSequence`\
protein = AminoAcidSequence("ACDEFGHIK")\
protein.count_amino_acids() \# =\> AA number\
\
==== Filtering reads by GC content, length, and quality: ====\
`from main import filter_fastq`

result = filter_fastq(input_fastq="data/example_fastq.fastq",
output_fastq="filtered/output_fastq.fastq", gc_bounds=(40, 60),
length_bounds=(100, 10000), quality_threshold=20.0, )\

==== FASTA and BLAS processing: ====\
`from main import convert_multiline_fasta_to_oneline, parse_blast_output`

FASTA conversion: multiline → one line per sequence
convert_multiline_fasta_to_oneline("data/input.fasta",
"filtered/oneline.fasta")

BLAST parsing: extracting unique protein descriptions
parse_blast_output("data/blast_results.txt",
"filtered/blast_parsed.txt")\

**Requirements**

-   Python 3.8+

-   Biopython 1.81+

**Author**

Loginova Olga

**Feedback and bug reports**

If you have any troubles running tools, please attach params.txt
