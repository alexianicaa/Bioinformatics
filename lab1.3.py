
#Use the AI to adapt your current algorithm in order to make an app that takes a FASTA file and read the sew content from it and display the rel percenteges for the symbolic present in the alphabet of seq. Note: FASTA represents a file format that contains DNA, ARN or proteins seq. Thus, it contains the information for your input

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

def read_fasta(file_path: str) -> str:
    """Reads a FASTA file and returns the sequence as a single string."""
    sequence = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):  # ignore headers and empty lines
                continue
            sequence.append(line.upper())
    return "".join(sequence)


def percentage(S: str):
    """Calculates the relative percentage of each symbol in the sequence."""
    length = len(S)
    if length == 0:
        return {}
    X = sorted(set(S))
    P = {char: (S.count(char) / length) * 100 for char in X}
    return P


def detect_seq(S: str) -> str:
    """Detects whether the sequence is DNA, RNA, Protein, or Unknown."""
    dna = set("ACGT")
    rna = set("ACGU")
    protein = set("ACDEFGHIKLMNPQRSTVWY")

    letters = set(S)

    if letters.issubset(dna):
        return "DNA"
    elif letters.issubset(rna):
        return "RNA"
    elif letters.issubset(protein):
        return "Protein"
    else:
        return "Unknown sequence type"

def browse_file():
    file_path = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
    )
    if not file_path:
        return
    
    try:
        S = read_fasta(file_path)
        seq_type = detect_seq(S)
        P = percentage(S)

        # Clear previous results
        for row in tree.get_children():
            tree.delete(row)

        # Insert new results
        for char, pct in P.items():
            tree.insert("", "end", values=(char, f"{pct:.2f}%"))

        label_seq_type.config(text=f"Sequence Type: {seq_type}")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to read file:\n{e}")

# --- GUI Setup ---
root = tk.Tk()
root.title("FASTA Sequence Analyzer")
root.geometry("400x400")

btn_browse = tk.Button(root, text="Browse FASTA File", command=browse_file)
btn_browse.pack(pady=10)

label_seq_type = tk.Label(root, text="Sequence Type: -", font=("Arial", 12, "bold"))
label_seq_type.pack(pady=5)

# Table for percentages
tree = ttk.Treeview(root, columns=("Symbol", "Percentage"), show="headings", height=10)
tree.heading("Symbol", text="Symbol")
tree.heading("Percentage", text="Percentage")
tree.pack(expand=True, fill="both", padx=10, pady=10)

root.mainloop()
