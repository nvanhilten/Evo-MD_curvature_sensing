hydropathy = { # Kyte J, Doolittle RF (May 1982). "A simple method for displaying the hydropathic character of a protein". Journal of Molecular Biology. 157 (1): 105–32. CiteSeerX 10.1.1.458.454. doi:10.1016/0022-2836(82)90515-0. PMID 7108955.
    "A" : 1.8, "G" :-0.4, "M" : 1.9, "S" :-0.8,
    "C" : 2.5, "H" :-3.2, "N" :-3.5, "T" :-0.7,
    "D" :-3.5, "I" : 4.5, "P" :-1.6, "V" : 4.2,
    "E" :-3.5, "K" :-3.9, "Q" :-3.5, "W" :-0.9,
    "F" : 2.8, "L" : 3.8, "R" :-4.5, "Y" :-1.3
}

CGsize = { # Number of CG Martini beads per residue.
    "A" : 1, "G" : 1, "M" : 2, "S" : 2,
    "C" : 2, "H" : 4, "N" : 2, "T" : 2,
    "D" : 2, "I" : 2, "P" : 2, "V" : 2,
    "E" : 2, "K" : 3, "Q" : 2, "W" : 5,
    "F" : 4, "L" : 2, "R" : 3, "Y" : 4
}

index = {
    "A" : 0, "G" : 5, "M" : 10, "S" : 15,
    "C" : 1, "H" : 6, "N" : 11, "T" : 16,
    "D" : 2, "I" : 7, "P" : 12, "V" : 17,
    "E" : 3, "K" : 8, "Q" : 13, "W" : 18,
    "F" : 4, "L" : 9, "R" : 14, "Y" : 19
}

def compute_hydropathy_sequence(sequence):
    hydro_seq = []
    for element in sequence:
        hydro_seq.append(hydropathy[element])

    return hydro_seq

def compute_hydropathy_seq_norm(sequence):
    return [i/max(map(abs, hydropathy.values())) for i in compute_hydropathy_sequence(sequence)]

def compute_CGsize_sequence(sequence):
    CGsize_seq = []
    for element in sequence:
        CGsize_seq.append(CGsize[element])

    return CGsize_seq

def compute_CGsize_seq_norm(sequence):
    return [i/max(map(abs, CGsize.values())) for i in compute_CGsize_sequence(sequence)]

