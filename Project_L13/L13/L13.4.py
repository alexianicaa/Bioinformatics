# Use the JSON file to create an engine that is able to use the transition porbabilities in order to produce an output similar to the original text from which the training was made. one has to create 2 separate apps, one producing DNA sequences and ine english text.
import json
import random

def generate_dna(json_file, length=50):
    with open(json_file, "r") as f:
        data = json.load(f)

    matrix = data["transition_matrix"]
    symbols = ["A", "C", "G", "T"]

    current = random.randint(0, 3)
    sequence = [symbols[current]]

    for _ in range(length - 1):
        probs = matrix[current]
        current = random.choices(range(4), probs)[0]
        sequence.append(symbols[current])

    return "".join(sequence)


if __name__ == "__main__":
    dna = generate_dna("dna_transition_matrix.json", 50)
    print(dna)
