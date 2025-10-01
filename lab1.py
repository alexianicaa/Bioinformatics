
#Make an application that is able to find the alphabet of a sequence of text. THis seq may be an ARN seq or ADN seq or protein seq.
def detect_seq(S: str) -> str:

    dna = set("ACGT")
    rna = set("ACGU")
    protein = set("ACDEFGHIKLMNPQRSTVWY")  

    letters = set(S)

    if letters.issubset(dna):
        return "DNA"
    elif letters.issubset(rna):
        return "RNA"
    elif letters.issubset(protein):
        return "Protein"
    else:
        return "Unknown sequence type"


print(detect_seq("ATGCGTACG")) #DNA
print(detect_seq("AUGCGAU"))     #RNA
print(detect_seq("MTEYKLVVVG")) #protein
print(detect_seq("XYZ")) #unknown       
