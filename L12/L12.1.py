import numpy as np
import pandas as pd

motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT"
]

# Count Matrix
def count_matrix(motifs):
    length = len(motifs[0])
    counts = pd.DataFrame(0, index=['A','C','G','T'], columns=range(length))
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            if nucleotide in counts.index:
                counts.loc[nucleotide, i] += 1
    return counts

counts = count_matrix(motifs)
print("Count Matrix:\n", counts)

# Weight Matrix
def weight_matrix(counts, pseudocount=1):
    total = counts.sum(axis=0) + 4*pseudocount
    pwm = (counts + pseudocount) / total
    return pwm

pwm = weight_matrix(counts)
print("\nWeight Matrix (PWM):\n", pwm)

# Relative Frequencies Matrix
rel_freq = pwm
print("\nRelative Frequencies:\n", rel_freq)

# Log-likelihoods Matrix
def log_likelihood_matrix(pwm, background=0.25):
    return np.log(pwm / background)

log_pwm = log_likelihood_matrix(pwm)
print("\nLog-likelihood Matrix:\n", log_pwm)

# Analyze sequence S
S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

def score_sequence(seq, log_pwm):
    k = log_pwm.shape[1]  # motif length
    scores = []
    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        score = 0
        for j, nucleotide in enumerate(window):
            if nucleotide in log_pwm.index:
                score += log_pwm.loc[nucleotide, j]
            else:
                score += 0
        scores.append((i, window, score))
    return scores

scores = score_sequence(S, log_pwm)

scores_df = pd.DataFrame(scores, columns=['Start','Window','Score'])
print("\nTop scoring windows:")
print(scores_df.sort_values(by='Score', ascending=False).head(10))
# The sequence S shows signals of exon-intron boundaries. 
# The highest scoring windows are at positions 0–8 (CAGGTTGGA) and 9–17 (AACGTAATC), suggesting these regions are likely exon-intron borders.