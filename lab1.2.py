
#A DNA seq is given: S = "ACGGGCATATGCGC". Make an app which is able to show the percentage of the components from the alphabet of the sequence S. In other words, the input of the seq S and output is the alphabet of the seq and the percentage of each letter in the alphabet found in seq S


def percentage (S: str):
    length = len(S)
    X=sorted(set(S))
    P = {char: (S.count(char) / length) * 100 for char in X}
    return P
    
S = "ACGGGCATATGCGC"

P = percentage(S)

print("Percentages: ")

for char, pct in P.items():
    print(f"{char}: {pct:.2f}%")