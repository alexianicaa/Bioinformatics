import numpy as np

def predict_steps(matrix, initial_vector, total_steps=5):

    A = np.array(matrix)
    v = np.array(initial_vector)
    
    if A.shape[0] != A.shape[1]:
        raise ValueError("The transition matrix must be square.")
    if A.shape[1] != v.shape[0]:
        raise ValueError("Matrix columns must match vector rows.")

    history = [v]
    current_v = v

    print(f"Initial State (Step 0): {current_v}")
    
    for step in range(1, total_steps + 1):
        current_v = np.dot(current_v, A)
        history.append(current_v)
        print(f"Step {step}: {current_v}")
        
    return history

A_matrix = [[0.3, 0.35, 0.35],
            [0, 0, 1],
            [0.9, 0, 0.1]]

v_start = [1,0,0]

results = predict_steps(A_matrix, v_start, total_steps=5)