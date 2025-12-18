# Download 10 influenza genomes. Adapt your application from the previous assignment in order to scan each genome for possible motifs. For each genome make a chart that shows the signal with the most likely locations of real functional motifs.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.signal import find_peaks

def get_log_likelihood(motifs):
    nucleotides = ['A', 'C', 'G', 'T']
    k = len(motifs[0])
    counts = np.zeros((4, k))
    ntIndex = {nt: i for i, nt in enumerate(nucleotides)}
    for motif in motifs:
        for i, nt in enumerate(motif):
            if nt in ntIndex:
                counts[ntIndex[nt], i] += 1
    
    pseudo = 1
    pwm = (counts + pseudo) / (len(motifs) + 4 * pseudo)
    return pd.DataFrame(np.log(pwm / 0.25), index=nucleotides)

def score_sequence(seq, log_pwm):
    k = log_pwm.shape[1]
    score_dict = log_pwm.to_dict()
    seq_chars = list(seq)
    scores = np.zeros(len(seq) - k + 1)
    for i in range(len(seq) - k + 1):
        window_score = 0
        for j in range(k):
            nt = seq_chars[i + j]
            if nt in score_dict[j]:
                window_score += score_dict[j][nt]
        scores[i] = window_score
    return scores

motifs = ["GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
          "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"]

log_pwm = get_log_likelihood(motifs)
folder = "influenza"

if not os.path.exists(folder):
    print(f"Folder '{folder}' not found.")
else:
    genome_files = sorted([f for f in os.listdir(folder) if f.endswith(".fasta")])
    num_genomes = len(genome_files)
    
    fig, axes = plt.subplots(num_genomes, 1, figsize=(15, 3 * num_genomes), sharex=False)
    
    if num_genomes == 1: axes = [axes]

    for idx, genome_file in enumerate(genome_files):
        path = os.path.join(folder, genome_file)
        record = SeqIO.read(path, "fasta")
        sequence = str(record.seq).upper()
        
        raw_scores = score_sequence(sequence, log_pwm)
        
        scores_norm = (raw_scores - np.min(raw_scores)) / (np.max(raw_scores) - np.min(raw_scores))
        peaks, _ = find_peaks(scores_norm, height=0.85, distance=30)
        
        ax = axes[idx]
        ax.plot(scores_norm, label='Score Signal', color='steelblue', linewidth=0.7)
        ax.plot(peaks, scores_norm[peaks], "ro", markersize=4, label='Predicted Border')
        
        ax.set_title(f"Influenza {idx + 1}", fontweight='bold')
        ax.set_ylabel("Norm Score")
        if idx == num_genomes - 1:
            ax.set_xlabel("Genome Position (bp)")
        
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize='small')

    plt.tight_layout()
    plt.show()