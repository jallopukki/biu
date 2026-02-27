def mean_quality(quality: str) -> float:
    "Calculate mean quality of read (Phred+33)."
    if not quality:
        return 0.0
    return sum(ord(char) - 33 for char in quality) / len(quality)


def normalize_bounds(bounds):
    "Normalize bounds to a tuple of floats (min, max)."
    if isinstance(bounds, (int, float)):
        return (0.0, float(bounds))
    return (float(bounds[0]), float(bounds[1]))


def calculate_gc(seq: str) -> float:
    "Calculate GC content percentage of a DNA sequence."
    if not seq:
        return 0.0
    gc_count = seq.upper().count("G") + seq.upper().count("C")
    return (gc_count / len(seq)) * 100


def read_fastq(input_fastq: str) -> dict:
    "Read a FASTQ file and return a dict (read_name: sequence, quality)."
    reads = {}
    with open(input_fastq, "r") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            read_name = header[1:]  # remove '@'
            reads[read_name] = (sequence, quality)
    return reads


def write_fastq(output_fastq: str, filtered_reads: dict):
    "Write filtered reads to a FASTQ file."
    with open(output_fastq, "w") as f:
        for name, (seq, qual) in filtered_reads.items():
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def filter_fastq(
    input_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0.0,
) -> dict:
    "Reads FASTQ and returns filtered dictionary."
    gc_min, gc_max = normalize_bounds(gc_bounds)
    len_min, len_max = normalize_bounds(length_bounds)

    reads = read_fastq(input_fastq)
    filtered = {}

    for name, (seq, quality) in reads.items():
        if not seq:
            continue

        gc_percentage = calculate_gc(seq)
        length = len(seq)
        avg_qual = mean_quality(quality)

        if (
            gc_min <= gc_percentage <= gc_max
            and len_min <= length <= len_max
            and avg_qual >= quality_threshold
        ):
            filtered[name] = (seq, quality)

    return filtered
