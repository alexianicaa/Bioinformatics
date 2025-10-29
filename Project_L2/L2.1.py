from itertools import product

S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

letters = [('A','C','G','T'), ('A','C','G','T')]
dinucleotide = [''.join(x) for x in list(product(*letters))]
#print(dinucleotide)

letters = [('A','C','G','T'), ('A','C','G','T'), ('A','C','G','T')]
trinucleotide = [''.join(x) for x in list(product(*letters))]
#print(trinucleotide)

def count_subs(seq, sub):
    count = 0
    start = 0
    while True:
        start = seq.find(sub, start)
        if start == -1:
            break
        count += 1
        start += 1
    return count

total_dinuc = len(S) - 1
total_trinuc = len(S) - 2 

dinuc_percentages = {}
for d in dinucleotide:
    count = count_subs(S, d)
    percent = (count / total_dinuc) * 100
    dinuc_percentages[d] = percent
    
trinuc_percentages = {}
for t in trinucleotide:
    count = count_subs(S, t)
    percent = (count / total_trinuc) * 100
    trinuc_percentages[t] = percent

print("=== Dinucleotide percentages ===")
for d, pct in dinuc_percentages.items():
    print(f"{d}: {pct:.2f}%")

print("\n=== Trinucleotide percentages ===")
for t, pct in trinuc_percentages.items():
    print(f"{t}: {pct:.2f}%")

