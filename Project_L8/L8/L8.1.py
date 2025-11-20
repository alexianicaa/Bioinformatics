import random

MIN_LEN = 300
MAX_LEN = 400
NUM_TES = 4
MUTATION_RATE = 0.05
SIMILARITY_THRESHOLD = 0.9
TE_SEQ = "TAAGGATCCGTTAAGCGATGCGATCCGTTA"

def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))

def mutate(seq, rate):
    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < rate:
            seq[i] = random.choice([b for b in 'ACGT' if b != seq[i]])
    return "".join(seq)

def detect_tes(dna, te, min_score):
    hits = []
    for i in range(len(dna)-len(te)+1):
        window = dna[i:i+len(te)]
        score = sum(a==b for a,b in zip(window, te))
        if score >= min_score:
            hits.append({'start':i+1, 'end':i+len(te), 'score':score, 'seq':window})
    return hits

seq_len = random.randint(MIN_LEN, MAX_LEN)
sequence = random_dna(max(50, seq_len - NUM_TES*len(TE_SEQ)))

te_locations = []
offset = 0
positions = [random.randint(0, len(sequence)) for _ in range(NUM_TES)]

for i, pos in enumerate(positions):
    te = mutate(TE_SEQ, MUTATION_RATE)
    real_pos = pos + offset
    sequence = sequence[:real_pos] + te + sequence[real_pos:]
    te_locations.append({'name':f'TE{i+1}', 'start':real_pos+1, 'end':real_pos+len(te), 'seq':te})
    offset += len(te)

print(f"\nSequence length: {len(sequence)} bp")
for te in te_locations:
    print(f"{te['name']} | {te['start']}-{te['end']} | {te['seq']}")

print("\n>synthetic_sequence")
for i in range(0, len(sequence), 60):
    print(sequence[i:i+60])

min_score = int(len(TE_SEQ) * SIMILARITY_THRESHOLD)
detected = detect_tes(sequence, TE_SEQ, min_score)

for d in detected:
    print(f"{d['start']}-{d['end']} | Score: {d['score']}/{len(TE_SEQ)} | {d['seq']}")

for i in range(len(detected)-1):
    if detected[i]['end'] >= detected[i+1]['start']:
        overlap = detected[i]['end'] - detected[i+1]['start'] + 1
        print(f"Overlap {detected[i]['start']}-{detected[i]['end']} and {detected[i+1]['start']}-{detected[i+1]['end']} ({overlap} bp)")

gt_starts = {te['start'] for te in te_locations}
TP = sum(1 for d in detected if d['start'] in gt_starts)
FP = len(detected) - TP
FN = NUM_TES - TP

print(f"\nTP: {TP} | FP: {FP} | FN: {FN}")
