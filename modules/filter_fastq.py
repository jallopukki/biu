def mean_quality(quality: str) -> float:
    """Calculate mean quality of read (Phred+33)."""
    if not quality:
        return 0.0
    return sum(ord(char) - 33 for char in quality) / len(quality)
