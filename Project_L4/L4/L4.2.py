#2 Download from NCBI the FASTA files containing the covid-19 genom and the influenza genom. Use the AI to compare the codons frequencies between the two. Make a chart that shows the top 10 most frequent codons for covid-19.
#b. Make a chart that shows the top 10 most frequent codons for influenza.
#c. Compare the two results and show the most frequent codons between the two
#d. Show in the output of the console top three amino acids for each genom


from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt

# --- Step 1: Read FASTA files ---
covid_file = "covid-19.fasta"       # SARS-CoV-2 complete genome
influenza_file = "influenza.fasta" # Influenza A genome

def read_sequence(fasta_file):
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.seq).upper().replace("T", "U")  # Convert DNA to RNA

# --- Step 2: Count codons ---
def get_codons(seq):
    return [seq[i:i+3] for i in range(0, len(seq)-2, 3) if len(seq[i:i+3]) == 3]

def codon_frequencies(seq):
    return Counter(get_codons(seq))

# --- Step 3: Genetic code dictionary ---
genetic_code = {
    'UUU':'Phe', 'UUC':'Phe', 'UUA':'Leu', 'UUG':'Leu',
    'CUU':'Leu', 'CUC':'Leu', 'CUA':'Leu', 'CUG':'Leu',
    'AUU':'Ile', 'AUC':'Ile', 'AUA':'Ile', 'AUG':'Met',
    'GUU':'Val', 'GUC':'Val', 'GUA':'Val', 'GUG':'Val',
    'UCU':'Ser', 'UCC':'Ser', 'UCA':'Ser', 'UCG':'Ser',
    'CCU':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro',
    'ACU':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr',
    'GCU':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala',
    'UAU':'Tyr', 'UAC':'Tyr', 'UAA':'STOP', 'UAG':'STOP',
    'CAU':'His', 'CAC':'His', 'CAA':'Gln', 'CAG':'Gln',
    'AAU':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys',
    'GAU':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu',
    'UGU':'Cys', 'UGC':'Cys', 'UGA':'STOP', 'UGG':'Trp',
    'CGU':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg',
    'AGU':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGG':'Arg',
    'GGU':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly'
}

def translate_codons_to_amino_acids(codons):
    amino_acids = [genetic_code.get(c, '') for c in codons if c in genetic_code]
    return Counter(amino_acids)

# --- Step 4: Foods low in each amino acid ---
low_aa_foods = {
    'Ala': ['apple', 'lettuce', 'rice'],
    'Arg': ['milk', 'cheese', 'apple', 'white rice'],
    'Asn': ['pear', 'lettuce', 'white rice'],
    'Asp': ['apple', 'cucumber', 'rice'],
    'Cys': ['orange', 'grapes', 'lettuce'],
    'Gln': ['banana', 'lettuce', 'rice'],
    'Glu': ['apple', 'rice', 'tomato'],
    'Gly': ['cucumber', 'lettuce', 'apple'],
    'His': ['rice', 'bread', 'pear'],
    'Ile': ['apple', 'white rice', 'lettuce'],
    'Leu': ['apple', 'rice', 'lettuce'],
    'Lys': ['apple', 'rice', 'cucumber'],
    'Met': ['apple', 'lettuce', 'white rice'],
    'Phe': ['apple', 'rice', 'cucumber'],
    'Pro': ['apple', 'lettuce', 'rice'],
    'Ser': ['berries', 'cucumber', 'white rice', 'bread'],
    'Thr': ['apple', 'rice', 'lettuce'],
    'Trp': ['apple', 'lettuce', 'rice'],
    'Tyr': ['apple', 'lettuce', 'rice'],
    'Val': ['apple', 'lettuce', 'rice']
}

def recommend_foods(amino_acids):
    recommendations = set()
    for aa in amino_acids:
        if aa in low_aa_foods:
            recommendations.update(low_aa_foods[aa])
    return recommendations

# --- Step 5: Run analysis ---
covid_seq = read_sequence(covid_file)
influenza_seq = read_sequence(influenza_file)

covid_codons = codon_frequencies(covid_seq)
influenza_codons = codon_frequencies(influenza_seq)

# --- Step 6: Plot Top 10 Codons in Two Separate Charts ---
def plot_top_codons(counter, title):
    top = counter.most_common(10)
    codons, freqs = zip(*top)
    plt.figure(figsize=(8, 4))
    plt.bar(codons, freqs, color='skyblue', edgecolor='black')
    plt.title(title)
    plt.xlabel("Codons")
    plt.ylabel("Frequency")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

plot_top_codons(covid_codons, "Top 10 Codons - COVID-19")
plot_top_codons(influenza_codons, "Top 10 Codons - Influenza")

# --- Step 7: Compare most frequent codons ---
covid_top = {c for c, _ in covid_codons.most_common(10)}
influenza_top = {c for c, _ in influenza_codons.most_common(10)}
common_codons = covid_top.intersection(influenza_top)

print("\nCommon Top Codons between COVID-19 and Influenza:")
print(common_codons)

# --- Step 8: Amino acid analysis ---
covid_amino = translate_codons_to_amino_acids(get_codons(covid_seq))
influenza_amino = translate_codons_to_amino_acids(get_codons(influenza_seq))

print("\nTop 3 Amino Acids for COVID-19:")
for aa, freq in covid_amino.most_common(3):
    print(f"{aa}: {freq}")

print("\nTop 3 Amino Acids for Influenza:")
for aa, freq in influenza_amino.most_common(3):
    print(f"{aa}: {freq}")

# --- Step 9: Recommend foods based on amino acids ---
covid_foods = recommend_foods([aa for aa, _ in covid_amino.most_common(3)])
flu_foods = recommend_foods([aa for aa, _ in influenza_amino.most_common(3)])

print("\nRecommended foods low in COVID-19 top amino acids:")
print(sorted(covid_foods))

print("\nRecommended foods low in Influenza top amino acids:")
print(sorted(flu_foods))
