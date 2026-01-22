import math

S1 = "ATCGATTCGATATCATACACGTAT"
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"
S = "CAGGTTGGAAACGTAA"

alphabet = ['A', 'C', 'G', 'T']


def build_model(sequence, alphabet):
    counts = {
        a: {b: 1 for b in alphabet}
        for a in alphabet
    }

    for i in range(len(sequence) - 1):
        counts[sequence[i]][sequence[i + 1]] += 1

    probs = {}
    for a in alphabet:
        total = sum(counts[a].values())
        probs[a] = {
            b: counts[a][b] / total for b in alphabet
        }

    return probs


def show_matrix(matrix, title):
    print(f"\n--- {title} ---")
    print("    " + " ".join(f"{b:>8}" for b in alphabet))
    for a in alphabet:
        row = f"{a:>3} "
        for b in alphabet:
            row += f"{matrix[a][b]:8.3f}"
        print(row)


def log_odds_matrix(model_pos, model_neg):
    beta = {}
    for a in alphabet:
        beta[a] = {}
        for b in alphabet:
            beta[a][b] = math.log2(
                model_pos[a][b] / model_neg[a][b]
            )
    return beta


def score(seq, beta):
    total = 0.0
    for i in range(len(seq) - 1):
        total += beta[seq[i]][seq[i + 1]]
    return total


model_plus = build_model(S1, alphabet)
model_minus = build_model(S2, alphabet)

show_matrix(model_plus, "CpG Island (+) Transition Probabilities")
show_matrix(model_minus, "Non-CpG (−) Transition Probabilities")

beta = log_odds_matrix(model_plus, model_minus)
show_matrix(beta, "Log-Likelihood Ratio Matrix (β)")

final_score = score(S, beta)

print("\nSequence:", S)
print("Total log-likelihood ratio:", round(final_score, 4))

if final_score > 0:
    print("Result: ISLAND (+)")
else:
    print("Result: NON ISLAND (−)")
