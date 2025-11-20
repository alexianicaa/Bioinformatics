## download from NCBS a total of 3 bacteria genom use these genoms as the input for your application, also modify your application in order to be able to handle the ammount of information of this genome. your application must detect transposable elements and the results from the output must show their position and their length. in this case the inverted repeated sequences are unknown unlike the previous assignment in which these are known. NOTE: the following cases must be taken in consideration 1. transposome envolving 2. transposome overlapping 3. the minumum size of the inveted repeats must be of 4 bases and the maximum size of the inverted repeats must be of 6
import random
import os

# --- HELPER FUNCTIONS ---

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]

def read_fasta(path):
    """Read a FASTA file and return sequence and header."""
    header = None
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:]
            else:
                seq_parts.append(line)
    return "".join(seq_parts).upper(), (header if header else os.path.basename(path))

def find_inverted_repeats(sequence, min_len=4, max_len=6, max_spacer=100):
    """Detect candidate Transposable Elements by searching for inverted repeats."""
    repeats = []
    n = len(sequence)
    for repeat_len in range(min_len, max_len + 1):
        for i in range(n - repeat_len):
            left = sequence[i:i + repeat_len]
            if set(left) - set('ACGT'):
                continue
            rc = reverse_complement(left)
            search_start = i + repeat_len
            search_end = min(i + repeat_len + max_spacer, n - repeat_len + 1)
            for j in range(search_start, search_end):
                right = sequence[j:j + repeat_len]
                if right == rc:
                    spacer = j - (i + repeat_len)
                    repeats.append({
                        'left_seq': left,
                        'right_seq': right,
                        'left_pos': i,
                        'right_pos': j,
                        'length': repeat_len,
                        'spacer': spacer
                    })
    return repeats

def filter_repeats(repeats, max_results=50):
    """Filter overlapping repeats, keeping the longest and most spaced TEs first."""
    sorted_repeats = sorted(repeats, key=lambda x: (-x['length'], -x['spacer']))
    filtered = []
    used_positions = set()
    for repeat in sorted_repeats:
        if len(filtered) >= max_results:
            break
        all_pos = set(range(repeat['left_pos'], repeat['right_pos'] + repeat['length']))
        if not all_pos.intersection(used_positions):
            filtered.append(repeat)
            used_positions.update(all_pos)
    return sorted(filtered, key=lambda x: x['left_pos'])

def analyze_genome(filename):
    """Analyze a genome file and detect candidate TEs."""
    print(f"\nAnalyzing: {filename}")
    sequence, genome_id = read_fasta(filename)
    print(f"  Header: {genome_id}")
    print(f"  Genome length: {len(sequence):,} bp")

    print(f"  Searching for inverted repeats (length 4-6, spacer ≤100)...")
    repeats = find_inverted_repeats(sequence)
    print(f"  Found {len(repeats):,} raw inverted repeat pairs.")

    filtered = filter_repeats(repeats)
    print(f"  Selected {len(filtered)} top non-overlapping candidates.")

    counts = {4: 0, 5: 0, 6: 0}
    for r in filtered:
        counts[r['length']] += 1

    return {
        'filename': filename,
        'id': genome_id,
        'length': len(sequence),
        'total_repeats': len(repeats),
        'filtered_repeats': filtered,
        'counts': counts
    }

def generate_report(results):
    """Create a text report of detected TEs."""
    print("\nGenerating text report...")
    filepath = 'transposon_report.txt'
    with open(filepath, 'w') as f:
        f.write("="*80 + "\n")
        f.write("TRANSPOSON DETECTION REPORT - De Novo IR Analysis\n")
        f.write("="*80 + "\n\n")
        f.write("METHODOLOGY:\n")
        f.write("-"*80 + "\n")
        f.write("Candidate Transposable Elements were identified by searching for Terminal Inverted Repeats (TIRs) in bacterial genomes.\n\n")
        f.write("Parameters:\n")
        f.write("  - Inverted repeat length: 4-6 bp\n")
        f.write("  - Max spacer distance: 100 bp\n")
        f.write("  - Overlapping repeats filtered to keep top candidates\n\n")

        for i, result in enumerate(results, 1):
            f.write(f"\n{'='*80}\n")
            f.write(f"GENOME {i}: {result['id']} ({result['filename']})\n")
            f.write(f"{'='*80}\n")
            f.write(f"Genome length: {result['length']:,} bp\n")
            f.write(f"Total raw IR pairs found: {result['total_repeats']:,}\n")
            f.write(f"Top non-redundant TE candidates: {len(result['filtered_repeats'])}\n\n")
            f.write("Distribution of detected IR lengths (top candidates):\n")
            for length in [4,5,6]:
                f.write(f"  {length} bp repeats: {result['counts'][length]}\n")

            f.write("\n--- SELECTED TE CANDIDATES (1-Based Positions) ---\n\n")
            f.write(f"{'#':<4} {'LIR_POS':<9} {'RIR_POS':<9} {'TE_START':<9} {'TE_END':<9} {'IR_LEN':<6} {'SPACER':<8} {'TE_LEN':<8} {'LIR_SEQ':<8}\n")
            f.write("-"*80 + "\n")

            for j, repeat in enumerate(result['filtered_repeats'], 1):
                te_start = repeat['left_pos'] + 1
                te_end = repeat['right_pos'] + repeat['length']
                te_len = te_end - te_start + 1
                f.write(f"{j:<4} {repeat['left_pos']+1:<9} {repeat['right_pos']+1:<9} {te_start:<9} {te_end:<9} "
                        f"{repeat['length']:<6} {repeat['spacer']:<8} {te_len:<8} {repeat['left_seq']:<8}\n")

        f.write("\n" + "="*80 + "\n")
        f.write("CONCLUSIONS\n")
        f.write("="*80 + "\n")
        f.write("Detected inverted repeats may represent Transposable Element candidates.\n")

    print(f"  Text report saved: {filepath}")

def main():
    print("="*80)
    print("TRANSPOSON DETECTION - INVERTED REPEAT ANALYSIS")
    print("="*80)

    fasta_files = ["Escherichia.fasta", "Pseudomonas.fasta", "Bacillus.fasta"]
    print(f"\nAnalyzing {len(fasta_files)} genomes: {', '.join(fasta_files)}")

    results = []
    for fasta_file in fasta_files:
        try:
            result = analyze_genome(fasta_file)
            results.append(result)
        except FileNotFoundError as e:
            print(f"  ERROR: {e}")
        except Exception as e:
            print(f"  ERROR processing {fasta_file}: {e}")

    if results:
        generate_report(results)
    else:
        print("No genomes successfully analyzed.")

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)

if __name__ == "__main__":
    main()
