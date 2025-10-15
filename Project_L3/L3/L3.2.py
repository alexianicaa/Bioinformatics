#Design an app that uses a sliding window method in order to read the tm over the seq S use a sliding winow of 8 positions and chose a FASTA file as input.
import math
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog

def tm_1(seq, Na = 0.0001):
        A = seq.count("A")
        T = seq.count("T")
        G = seq.count("G")
        C = seq.count("C")
        
        tm = 4 * (G + C) + 2 * (A + T)
        return tm

def tm_2(seq, Na = 0.0001):
        G = seq.count("G")
        C = seq.count("C")
        length = len(seq)
        
        tm = -(81.5 + 16.6 * math.log10(Na) + 0.41 * ((G + C)/length*100) - 600/length)
        return tm

def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper()

def sliding_window_tm(sequence, window_size=8, method=1, Na=0.001):
    results = []
    for i in range(len(sequence) - window_size + 1):
        window_seq = sequence[i:i + window_size]
        if method == 1:
            tm_value = tm_1(window_seq, Na)
        elif method == 2:
            tm_value = tm_2(window_seq, Na)
        else:
            raise ValueError("Method must be 1 or 2.")
        results.append({
            "start": i + 1,
            "end": i + window_size,
            "window": window_seq,
            "tm": tm_value
        })
    return results

def plot_tm(results, method):
    positions = [r["start"] for r in results]
    tm_values = [r["tm"] for r in results]

    plt.figure(figsize=(10, 5))
    plt.plot(positions, tm_values, marker='.', linestyle='-', linewidth=1.5)
    plt.title(f"Sliding Window Tm (window size = 8, method = {method})")
    plt.xlabel("Window start position along sequence")
    plt.ylabel("Melting Temperature (°C)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def main():
    Tk().withdraw()
    fasta_file = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
    )

    if not fasta_file:
        print("❌ No file selected. Exiting.")
        return

    print(f"✅ Selected file: {fasta_file}")

    try:
        seq = read_fasta(fasta_file)
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return

    print(f"Loaded sequence length: {len(seq)} bases")

    try:
        method = int(input("Choose Tm calculation method (1 or 2): ").strip() or "1")
    except ValueError:
        method = 1

    results = sliding_window_tm(seq, window_size=8, method=method, Na=0.001)

    print("\nSliding Window Tm Results:")
    for r in results:
        print(f"{r['start']:>4}-{r['end']:>4} | {r['window']} | Tm = {r['tm']:.2f} °C")

    plot_tm(results, method)

if __name__ == "__main__":
    main()