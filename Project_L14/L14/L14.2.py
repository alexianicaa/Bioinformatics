import math
import string
import matplotlib.pyplot as plt
import numpy as np


text_eminescu_train = """
Somnoroase păsărele
Pe la cuiburi se adună,
Se ascund în ramurele —
Noapte bună!

Doar izvoarele suspină,
Pe când codrul negru tace;
Dorm și florile-n grădină —
Dormi în pace!

Trece lebăda pe ape
Între trestii să se culce —
Fie-ți îngerii aproape,
Somnul dulce!

Peste-a nopții feerie
Se ridică mândra lună,
Totu-i vis și armonie —
Noapte bună!

Când răsare zorile
Și codrul se întinde lin,
Florile-și deschid petalele,
Iar lumea se trezește blând.
"""

text_stanescu_train = """
Cuvintele se ating unele pe altele
ca niște ființe fără trup,
iar gândul meu se îndoaie
în forma aerului dintre ele.

Timpul merge desculț
prin propoziții neterminate,
iar eu respir sensuri
care nu au fost încă inventate.

Umbrele se înalță în mijlocul zilei,
iar sunetul tăcerii vorbește mai tare
decât vocea oricărui cuvânt.
Fiecare idee își caută forma,
iar eu o prind în palme invizibile
și o las să zboare înapoi în lumină.

Și noaptea cade, dar fără întuneric,
iar lumina se împletește cu gândurile,
creând un univers de idei
fără margini, fără început și fără sfârșit.
"""

text_mihai_test = """
Noaptea cade peste cuvinte
iar păsările gândului tac.
Sensurile adorm pe ramuri,
dar timpul încă respiră.

Izvoarele vorbesc despre idei,
despre liniște și despre aer,
iar propozițiile se rup
în somn și în lumină.

Lebăda trece pe ape
și codrul se întinde ușor.
Florile-și deschid petalele,
iar gândurile zboară libere.

Cuvintele se ating, se împletesc,
iar ideile plutesc printre umbre.
Noaptea și ziua se amestecă
într-un dans de lumină și tăcere.

"""


def preprocess(text):
    text = text.lower().replace('\n', ' ')
    table = str.maketrans('', '', string.punctuation)
    text = text.translate(table)
    return text.split()

words_e = preprocess(text_eminescu_train)
words_s = preprocess(text_stanescu_train)
words_test = preprocess(text_mihai_test)

vocabulary = sorted(list(set(words_e) | set(words_s) | set(words_test)))
word_to_id = {word: i for i, word in enumerate(vocabulary)}
id_to_word = {i: word for word, i in word_to_id.items()}

ids_e = [word_to_id[w] for w in words_e]
ids_s = [word_to_id[w] for w in words_s]
ids_test = [word_to_id[w] for w in words_test]

def train_markov_model(ids, vocab_size):
    counts = np.ones((vocab_size, vocab_size))
    
    for k in range(len(ids) - 1):
        counts[ids[k]][ids[k+1]] += 1
            
    probs = counts / counts.sum(axis=1, keepdims=True)
    return probs

vocab_size = len(vocabulary)
prob_e = train_markov_model(ids_e, vocab_size)
prob_s = train_markov_model(ids_s, vocab_size)

window_size = 5
scores = []
segments = []

for i in range(len(ids_test) - window_size + 1):
    window = ids_test[i : i + window_size]
    window_score = 0
    
    for j in range(len(window) - 1):
        p_e = prob_e[window[j]][window[j+1]]
        p_s = prob_s[window[j]][window[j+1]]
        window_score += math.log2(p_e / p_s)
    
    scores.append(window_score)
    segments.append(" ".join([id_to_word[idx] for idx in window]))


plt.figure(figsize=(12, 6))
plt.plot(scores, color='gray', linestyle='--', alpha=0.5)
plt.axhline(0, color='black', linewidth=1)

x = np.arange(len(scores))
plt.fill_between(x, scores, 0, where=(np.array(scores) > 0), color='green', alpha=0.7, label='Stil Eminescu')
plt.fill_between(x, scores, 0, where=(np.array(scores) < 0), color='red', alpha=0.7, label='Stil Stănescu')

plt.title("Stylometric Analysis (Mihai)")
plt.xlabel("Word Position in Test Text")
plt.ylabel("Log-Likelihood Score (Eminescu vs Stănescu)")
plt.legend()
plt.grid(alpha=0.2)
plt.show()

print(f"{'TEXT SEGMENT':<50} | {'SCORE':<8} | {'PREDICTION'}")
print("-" * 75)
for i in range(0, len(scores), 3): 
    v = scores[i]
    if v > 1: res = "EMINESCU"
    elif v < -1: res = "STANESCU"
    else: res = "NEITHER"
    
    print(f"{segments[i][:48]:<50} | {v:>8.2f} | {res}")