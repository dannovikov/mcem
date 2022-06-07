import math
import random

def perturb_sequences(seqs, p):
    """
    Input: 
    Seqs: List of sequences
    p: Probability of mutation
    """

    N = len(seqs)
    L = len(seqs[0])

    letter_counts = [None] * L

    for k in range(L):
        letter_counts[k] = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}

    for seq in seqs:
        for k in range(L):
            letter = seq[k]
            letter_counts[k][letter] += 1

    R = [0] * L
    R_total = 0
    consensus = ['N'] * L

    chars = ['A','C','G','T']
    for k in range(L):
        counts = letter_counts[k]
        if counts['N'] == N:
            R[k] = 0
        else:
            max_val = 0
            random.shuffle(chars)
            for char in chars:
                if counts[char] > max_val:
                    max_val = counts[char]
                    consensus[k] = char
            R[k] = (N - counts['N'] - max_val) / (N - counts['N'])
            R_total += R[k]

    # Should we check if R_total == 0? I.e., is it realistic for 'N" to be the consensus everywhere?

    c = p * N * L / R_total
    #cR = [math.floor(c * r) for r in R]
    cR = [round(c * r) for r in R]

    for k in range(L):
        if cR[k] > 0:
            rows = [i for i in range(N) if (seqs[i][k] != consensus[k]) and (seqs[i][k] != 'N')]
            n_mut = min(len(rows),cR[k])
            rows_mut = random.sample(rows, n_mut)
            for i in rows_mut:
                chars_mut = [char for char in ['A','C','G','T'] if char != seqs[i][k]]
                seqs[i] = seqs[i][:k] + random.choice(chars_mut) + seqs[i][(k+1):]

    return seqs

if __name__ == '__main__':
    seqs = ["AAAACAAA","CAAACAAA","ACCACCAC","ANCATNAN","TTTTTTTT","GGGGGGGG"]
    p = 0.25 # Intended to be 0.01
    seqs = perturb_sequences(seqs, p)
    print(seqs)