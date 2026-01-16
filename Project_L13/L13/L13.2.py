
#2.Use a DNA sequence of 50 letetrs in oreder to calucalte the transition matrix that represent the sequence and store it in a JSON file.

import json
import numpy as np

def dna_transition_matrix(dna_sequence):
    # MEMORY TYPE 1: SYMBOL MEMORY
    symbols = ['A', 'C', 'G', 'T']
    symbol_to_id = {s: i for i, s in enumerate(symbols)}

    # MEMORY TYPE 2: TRANSITION MEMORY
    matrix = np.zeros((4, 4))

    for i in range(len(dna_sequence) - 1):
        current = symbol_to_id[dna_sequence[i]]
        next_state = symbol_to_id[dna_sequence[i + 1]]
        matrix[current][next_state] += 1

    # Normalize rows
    row_sums = matrix.sum(axis=1, keepdims=True)
    matrix = np.divide(matrix, row_sums, where=row_sums != 0)

    return {
        "symbol_memory": symbol_to_id,
        "transition_matrix": matrix.tolist()
    }


dna_sequence = "ACGTGCTAGCTAGGCTAACGTAGCTAGCTAGGCTAACGTAGCTAGC"  # 50 letters

dna_data = dna_transition_matrix(dna_sequence)

with open("dna_transition_matrix.json", "w") as f:
    json.dump(dna_data, f, indent=4)

print("DNA transition matrix saved.")
