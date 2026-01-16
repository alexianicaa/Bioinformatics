#3. Use an english text of a length od 300 letters and calculate the transition probabilities between the words inside the text and store that as a JSON file.
#Represent each new word with a unique symbol
import json
import re
from collections import defaultdict

def text_transition_probabilities(text):
    words = re.findall(r'\b[a-z]+\b', text.lower())

    # MEMORY TYPE 1: SYMBOL MEMORY
    unique_words = sorted(set(words))
    word_to_id = {word: i for i, word in enumerate(unique_words)}

    # MEMORY TYPE 2: TRANSITION MEMORY
    transitions = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)

    for i in range(len(words) - 1):
        current = word_to_id[words[i]]
        next_word = word_to_id[words[i + 1]]
        transitions[current][next_word] += 1
        totals[current] += 1

    # Convert to probabilities
    transition_probabilities = {}
    for current_id, next_words in sorted(transitions.items()):
        transition_probabilities[current_id] = {
            str(next_id): count / totals[current_id]
            for next_id, count in next_words.items()
        }

    return {
        "symbol_memory": word_to_id,
        "transition_probabilities": transition_probabilities
    }


text = """
Kindergarten - I was five years old when I started kindergarten, but with a September birthday, I turned six soon after. Mrs. Baker was round and wore dresses every day. The first day of school, a little girl in my class cried and cried. She cried every day for weeks. I was curious about her. I watched her come in with her mother. Outside the door, she didn’t shed a tear. But once her mother guided her inside our classroom door, the water works began. She cried all morning. When she wasn’t shriek-crying, she was sobbing, her shoulders raising and falling dramatically. I don’t remember tears, but I do remember that she had a hyper-productive snot gland. She stood there helplessly while Mrs. Baker tried to calm her. Every morning it was the same. Every morning, Mrs. Baker assured her mother that Miss Cries-a-lot just needed a little more time. Months passed. Suddenly one day the classroom was quiet. She was gone and I knew exactly where she was before Mrs. Baker even made the announcement. She was next door, in another teacher’s room. You see, that little crybaby didn’t cry because she missed her mommy. She cried because she wanted to be in the other kinder class – the one with more toys. That’s what she howled about every single day. But that’s not the story Mrs. Baker told. During circle time that morning, she announced that the little girl was now in the other classroom because there were too many kids in her class. I suspect that every single one of my classmates had the same thought…”move me…move me.”

After circle time, all the kids were playing when I approached Mrs. Baker. “She really got to move because she cried, right?” I asked. She didn’t answer. “Go play,” she said.
She would later tell my mother how smart I was. “She figured out exactly what happened.” She shared.

It didn’t take a rocket scientist.

Now it's your turn. Open up a page on your word procesor and begin with Kindergarten. Write just a paragraph or two of the first memory that comes to mind. Tie it up quickly and stop.
Then challenge yourself to add to it regularly. As a writer, I find it is a good "warm-up" to the writing I must get done on a particular day. Or you may want to use it if you are stuck for what to write. And in the mean time you are doing something even more amazing, you are taking something that is yours and yours alone- a memory- and recording it for generations to come. Good for you! Write Now, Karen
"""

text_data = text_transition_probabilities(text)

with open("word_transition_probabilities.json", "w") as f:
    json.dump(text_data, f, indent=4)

print("Word transition probabilities saved.")
