"""
Reference sequences for building pMHC single-chain trimers (SCTs).

Pulled out of analyses/05_SCT_generation.ipynb so the same constants back both
the notebook and analyses/06_HLA_test.py. Everything here is reference data
only; the codon optimisation itself lives in generate_pMHC_SCT.py.

Nucleotide sequences are stored lowercase, amino acid sequences uppercase,
matching how they were recorded in the notebook.
"""

# Flexible glycine-serine linkers joining the SCT components.
# L1 joins peptide to B2M, L2 joins B2M to the HLA heavy chain.
linker_nuc = {
    "L1": "ggtggaggaggttctggaggtggtggtagtggtggtggtggttcc",
    "L2": "ggtggtggtggtagtggtggtggtggttcaggtggtggtggttccggtggtggtggttcc",
}

# Beta-2 microglobulin, mature chain (leader removed).
B2M_nuc = "atacaaagaactccaaagatccaagtttacagtagacatcctgctgaaaacggtaaatctaatttcttgaactgttacgtctccggtttccacccaagtgatatagaagttgacttgttgaaaaatggtgaaagaatcgaaaaggttgaacattcagatttgtctttttctaaggactggtccttctatttgttgtactacacagaattcactccaactgaaaaggatgaatacgcttgcagagttaatcatgtaaccttgtctcaacctaaaatcgttaagtgggatagagacatg"

# Homology arms for Gibson assembly into each destination vector.
destination_vector_nuc = {
    "pHIV-EGFP": {
        "5_addon": "tcgtgagcggccgctgagtt",
        "3_addon": "aactattctagagtacccgggctaggatcc",
    },
    "pHIV-dTomato": {
        "5_addon": "gaactgaattcatcgacgtt",
        "3_addon": "aactattctagagtacccgggctaggatcc",
    },
}

# Leader (signal) peptides, as amino acids: these are codon optimised at run
# time rather than stored as nucleotides.
leader_peptide_aa = {
    "B2M": "MSRSVALAVLALLSLSGLEA",
    "HGH": "MATGSRTSLLLAFGLLCLPWLQEGSA",
}

# HLA heavy chains (alpha1-alpha3, no leader, no transmembrane domain) nucleotide sequence.
HLA_nuc = {
    "HLA-A2": "catagtatgagatatttctttacttctgtttcaagaccaggtagaggtgaacctagattcatcgcagtcggttacgttgatgacacacaatttgtaagattcgattccgacgctgcaagtcaaagaatggaaccaagagcaccttggattgaacaagaaggtccagaatattgggatggtgaaactagaaaagttaaggcccattctcaaactcacagagtagatttgggtacattaagaggtgcttataatcaatctgaagcaggttcacatacagtacaaagaatgtacggttgtgatgtcggttcagactggagatttttgagaggttatcaccaatatgcttacgatggtaaagactacattgcattgaaggaagatttgagatcctggaccgccgctgacatggcagcccaaactacaaaacataagtgggaagctgcacacgtagcagaacaattgagagcctatttggaaggtacatgtgtcgaatggttgagaagatacttagaaaacggtaaagaaacattgcaaagaaccgatgctccaaagactcatatgacacatcacgccgttagtgatcacgaagctactttgagatgctgggcattatctttttaccctgccgaaatcacattgacctggcaaagagatggtgaagaccaaacccaagatactgaattagttgaaaccagaccagcaggtgacggtactttccaaaaatgggccgctgttgtagtcccttcaggtcaagaacaaagatacacatgccatgtccaacacgaaggtttaccaaagccattgacattgagatgggaaccatcc",
}

# HLA heavy chains (alpha1-alpha3, no leader, with transmembrane domain) amino acid sequence
HLA_aa = {
    "HLA-A201_BM": "HSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGAYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSS",
    "HLA-A201_RS": "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGAYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSS",
}

# Order the SCT components are joined in:
# signal peptide - peptide - L1 - B2M - L2 - HLA heavy chain
SCT_CHAIN_ORDER = ("5_addon", "SP", "peptide", "L1", "B2M", "L2", "HLA", "3_addon")
