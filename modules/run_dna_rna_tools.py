def transcribe(seq: str) -> str:
    """Makes a T to U substitution. Returns the transcribed sequence."""
    return seq.replace("T", "U").replace("t", "u")


def reverse(seq: str) -> str:
    """Returns the reversed (mirror) sequence."""
    return seq[::-1]


def complement(seq: str) -> str:
    """Returns the complementary sequence."""
    upper_seq = seq.upper()
    if "U" in upper_seq:
        table = str.maketrans("AUGCaugc", "UACGuacg")
    else:
        table = str.maketrans("ATGCatgc", "TACGtacg")
    return seq.translate(table)


def reverse_complement(seq: str) -> str:
    """Returns the reverse complementary sequence."""
    return complement(reverse(seq))
