
#Design an application using the AI, wich contains the GUI that allows the user to choose a FASTA file. The comntent of the file should be analysed by using a sliding window of 30 positions. The content of each sliding window should be used in oreder to extract the procenteges, the relative freq o the symbols found in the alphabet of the seq. Thus your input will be the input DNA seq from the FASTA file and the output should be the values of the rel freq. of the each symbol in the alphabet of the seq. Translate the lines on a chart.
import re
from collections import Counter
import tkinter as tk
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

# -------------------------
# Helper functions
# -------------------------
def read_fasta(file_path):
    """Read FASTA file and return sequence containing only A/C/G/T"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    seq = ''.join([line.strip() for line in lines if not line.startswith('>')])
    seq = seq.upper()
    cleaned_seq = re.sub('[^ACGT]', '', seq)
    return cleaned_seq

def calculate_relative_freq(window_seq):
    counts = Counter(window_seq)
    length = len(window_seq)
    if length == 0:
        return {'A':0, 'C':0, 'G':0, 'T':0}
    return {b: counts.get(b,0)/length for b in 'ACGT'}

def sliding_window_analysis(seq, window_size=30, step=1):
    results = []
    L = len(seq)
    if L == 0:
        return results
    if L <= window_size:
        results.append(calculate_relative_freq(seq))
        return results
    for i in range(0, L - window_size + 1, step):
        window_seq = seq[i:i+window_size]
        results.append(calculate_relative_freq(window_seq))
    return results

def smooth_series(data, smooth_window=10):
    """Apply moving average smoothing"""
    if len(data) < smooth_window:
        return data
    smoothed = []
    half = smooth_window // 2
    for i in range(len(data)):
        start = max(0, i - half)
        end = min(len(data), i + half + 1)
        smoothed.append(sum(data[start:end]) / (end-start))
    return smoothed

# -------------------------
# Plotting function (embedded)
# -------------------------
def plot_results_in_tkinter(frame, results):
    for widget in frame.winfo_children():
        widget.destroy()

    if not results:
        messagebox.showinfo("Info", "No data to plot.")
        return

    x = list(range(1, len(results)+1))
    A = smooth_series([r['A'] for r in results], smooth_window=10)
    C = smooth_series([r['C'] for r in results], smooth_window=10)
    G = smooth_series([r['G'] for r in results], smooth_window=10)
    T = smooth_series([r['T'] for r in results], smooth_window=10)

    fig, ax = plt.subplots(figsize=(12,5))
    ax.plot(x, A, label='A', color='green', linewidth=2)
    ax.plot(x, C, label='C', color='blue', linewidth=2)
    ax.plot(x, G, label='G', color='orange', linewidth=2)
    ax.plot(x, T, label='T', color='red', linewidth=2)
    ax.set_xlabel("Sliding Window Index")
    ax.set_ylabel("Relative Frequency")
    ax.set_title("Relative Frequencies of A, C, G, T (Sliding Window = 30)")
    ax.set_ylim(0,1)
    ax.legend()
    ax.grid(True)

    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# -------------------------
# GUI actions
# -------------------------
def select_fasta_file():
    file_path = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files","*.fasta *.fa *.txt"), ("All files","*.*")]
    )
    if not file_path:
        return
    try:
        seq = read_fasta(file_path)
        if len(seq) == 0:
            messagebox.showwarning("Warning", "Sequence contains no valid bases (A,C,G,T).")
            return
        results = sliding_window_analysis(seq, window_size=30, step=1)
        plot_results_in_tkinter(plot_frame, results)
    except Exception as e:
        messagebox.showerror("Error", f"Failed to process file: {e}")

# -------------------------
# Main window
# -------------------------
root = tk.Tk()
root.title("DNA Sliding Window Analyzer")
root.geometry("1000x500")

top_frame = tk.Frame(root)
top_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)

label = tk.Label(top_frame, text="Analyze DNA sequence from a FASTA file (Sliding window = 30, Step = 1)", wraplength=900)
label.pack(side=tk.LEFT, padx=5)

button = tk.Button(top_frame, text="Choose FASTA File", command=select_fasta_file)
button.pack(side=tk.RIGHT, padx=5)

plot_frame = tk.Frame(root)
plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

root.mainloop()
