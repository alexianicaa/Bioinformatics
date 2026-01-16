import json
import random

def generate_text(json_file, word_count=100):
    with open(json_file, "r") as f:
        data = json.load(f)

    id_to_word = {str(v): k for k, v in data["symbol_memory"].items()}
    transitions = data["transition_probabilities"]

    current = random.choice(list(transitions.keys()))
    words = [id_to_word[current]]

    for _ in range(word_count - 1):
        if current not in transitions:
            break

        next_ids = list(transitions[current].keys())
        probabilities = list(transitions[current].values())

        current = random.choices(next_ids, probabilities)[0]
        words.append(id_to_word[current])

    return " ".join(words)


if __name__ == "__main__":
    text = generate_text("word_transition_probabilities.json", word_count=120)
    print("Generated text:")
    print(text)
