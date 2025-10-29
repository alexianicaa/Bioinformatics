
# Find in sequence S only the dinucleotides and trinucleotides that exists, without the use of the brute force engine. In order to acieve the results one must check these combinations starting from the begining of the seq. until the end of the seq. Ex. S = "ABAA" 
S="ABAA"

dinucleotides = set()
for i in range(len(S) - 1):
    dinucleotides.add(S[i:i+2])

trinucleotides = set()
for i in range(len(S) - 2):
    trinucleotides.add(S[i:i+3])

print("Dinucleotides found:", sorted(dinucleotides))
print("Trinucleotides found:", sorted(trinucleotides))