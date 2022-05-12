def create_seqs_dict(fasta):
    """
    From a fasta file, the function builds a dictionary {label:sequence}.
    In fasta format, the even # lines are labels that start with ">",
    the odd # lines are sequences.
    """
    seqs = {}
    with open(fasta, "r") as f:
        last_label = ""
        for line, text in enumerate(f):
            if line % 2 == 0:  # is label
                # Create new entry in dictionary for upcoming sequence
                seqs[text.strip()[1:]] = ""
                last_label = text.strip()[1:]
            else:
                # Add sequence to newly created entry
                seqs[last_label] = text.strip()
    return seqs


def create_seqs_matrix(seqs, SEQ_LEN):
    """
    ONE-HOT ENCODING
    seqs is a dict mapping id:sequence
    return N x M x 5 tensor where each row is a M x 5 one-hot encoding of an M-length sequence
    """
    seqs_m = np.zeros(shape=(len(seqs), SEQ_LEN, 5))
    seqs_index = {}
    # for i, seq_id in enumerate(tqdm(seqs)):
    for i, seq_id in enumerate(seqs):
        seqs_index[seq_id] = i
        for j, nucl in enumerate(seqs[seq_id]):
            if nucl == "A":
                seqs_m[i][j][0] = 1
            elif nucl == "C":
                seqs_m[i][j][1] = 1
            elif nucl == "T":
                seqs_m[i][j][2] = 1
            elif nucl == "G":
                seqs_m[i][j][3] = 1
            else:
                assert nucl == "-" or nucl == "N", f"nucl: {nucl}"
                seqs_m[i][j][4] = 1
    return seqs_m, seqs_index


def create_seqs_matrix_numerical(seqs, SEQ_LEN):
    """
    seqs is a dict mapping id:sequence
    return index dictionary and an N x M matrix where each row encodes a sequence as {A,C, T, G, -} -> {0, 1, 2, 3, 4}
    """
    seqs_m = np.zeros(shape=(len(seqs), SEQ_LEN))
    seqs_index = {}
    # for i, seq_id in enumerate(tqdm(seqs)):
    for i, seq_id in enumerate(seqs):
        seqs_index[seq_id] = i
        for j, nucl in enumerate(seqs[seq_id]):
            if nucl == "A":
                seqs_m[i][j] = 0
            elif nucl == "C":
                seqs_m[i][j] = 1
            elif nucl == "T":
                seqs_m[i][j] = 2
            elif nucl == "G":
                seqs_m[i][j] = 3
            else:
                assert nucl == "-" or nucl == "N"
                seqs_m[i][j] = 4
    return seqs_m, seqs_index