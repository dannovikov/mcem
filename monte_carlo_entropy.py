"""
    nwk, fasta, out_nwk, num_iter = sys.argv[1:5]
    try: existing_tree = sys.argv[5]
    except: existing_tree = False
"""

import dendropy
import copy
from tqdm import tqdm
import numpy as np
import itertools
import sys
from matplotlib import pyplot as plt

from utils.seq_utils import create_seqs_dict

sys.setrecursionlimit(10000)

LABEL_COUNT = 0
SEQ_LEN = 0

def compute_initial_tree_entropy(tree, counts):
    """
    This function computes entropy of every cluster (node) in the tree.
    It is designed to run be used once to obtain initial entropies each cluster.
    """
    tree_entropy = 0
    for cluster in tree.nodes():
        clust_entropy = compute_cluster_entropy(cluster, counts)
        tree_entropy += clust_entropy * counts[cluster.label]["size"]
        counts[cluster.label]["entropy"] = clust_entropy
    return -1 * tree_entropy


def compute_tree_entropy(tree, counts):
    tree_entropy = 0
    for node in tree.nodes():
        entropy = counts[node.label]["entropy"]
        size = counts[node.label]["size"]
        tree_entropy += entropy * size
    return -1 * tree_entropy


def compute_tree_entropy_dividebyparent(tree, counts):
    tree_entropy = 0
    for cluster in tree.nodes():
        if cluster == tree.seed_node:
            continue
        entropy = counts[cluster.label]["entropy"]
        size = counts[cluster.label]["size"]
        parent = cluster.parent_node
        parent_size = counts[parent.label]["size"]
        tree_entropy += (entropy * size) / parent_size
        tree_entropy += entropy
    return -1 * tree_entropy


def compute_cluster_entropy(cluster, counts):
    col_entropy = 0
    size = counts[cluster.label]["size"]
    for i in range(4):  # {A,C,T,G}
        for pos in range(SEQ_LEN):
            count = counts[cluster.label]["counts"][pos, i]  # of A's @ pos i
            p_i = count / size
            # assert 0 <= p_i <= 1, f"ERROR: mc encountered {p_i=}, {count=}, {pos=}, {SEQ_LEN=}"
            entropy = p_i * np.log2(p_i) if p_i != 0 else 0
            col_entropy += entropy
    return col_entropy / SEQ_LEN


def create_counts_matrices(tree, seqs_m, seqs_index, SEQ_LEN):
    counts = {}
    for i in tree.nodes():
        leaves = [k.taxon.label for k in i.leaf_iter()]
        count = np.zeros(shape=(SEQ_LEN, 5))
        for leaf in leaves:
            idx = seqs_index[leaf.replace(" ", "_")]
            seq = seqs_m[idx]
            count += seq
        counts[i.label] = {"counts": count, "size": len(leaves), "entropy": 0}
    return counts


def update_counts(counts, i, i_parent, j_parent, lca):
    """
    After a move, this function is called to update
        1. the nucleotide counts
        2. the entropy
    of all affected ancestor clusters; those until the LCA of i and j

    Params:
        i: node i being moved
        i_parent: grandparent of i
        j_parent: destination parent
    """
    p = i_parent
    while p and p != lca:
        counts[p.label]["counts"] -= counts[i.label]["counts"]
        assert all(j >= 0 for i in counts[p.label]['counts'] for j in i)
        counts[p.label]["size"] -= counts[i.label]["size"]
        counts[p.label]["entropy"] = compute_cluster_entropy(p, counts)
        p = p.parent_node

    p = j_parent
    while p and p != lca:
        counts[p.label]["counts"] += counts[i.label]["counts"]
        assert all(j >= 0 for i in counts[p.label]['counts'] for j in i), "Error in update counts adding"
        counts[p.label]["size"] += counts[i.label]["size"]
        counts[p.label]["entropy"] = compute_cluster_entropy(p, counts)
        p = p.parent_node

    return counts


def update_counts_after_reroot(counts, i, j_parent, new_root):
    p = j_parent
    while p:
        counts[p.label]["counts"] += counts[i.label]["counts"]
        counts[p.label]["size"] += counts[i.label]["size"]
        counts[p.label]["entropy"] = compute_cluster_entropy(p, counts)
        if p == new_root:
            break
        p = p.parent_node
    return counts


def move(tree, node_1, node_2, counts):
    n1_parent= node_1.parent_node
    n1_sibling = node_1.sibling_nodes()[0]
    n2_parent = node_2.parent_node
    tree.prune_subtree(node_1, suppress_unifurcations=False)
    tree.prune_subtree(node_2, suppress_unifurcations=False)
    new_internal = create_node()
    insert_internal_node(new_internal, parent=n2_parent, children=[node_1, node_2])
    if reroot:=(n1_parent == tree.seed_node):
        tree.reroot_at_node(n1_sibling, suppress_unifurcations=False)
        tree.prune_subtree(n1_parent, suppress_unifurcations=False)
    else:
        n1_grandparent = n1_parent.parent_node
        tree.prune_subtree(n1_sibling, suppress_unifurcations=False)
        tree.prune_subtree(n1_parent, suppress_unifurcations=False)
        n1_grandparent.add_child(n1_sibling)
    counts = add_counts_entry(new_internal, node_1, node_2, counts)
    return counts, reroot


def find_valid_move(tree):
    non_root_nodes = tree.nodes(lambda x: x != tree.seed_node)
    i = np.random.choice(non_root_nodes)
    j = np.random.choice(non_root_nodes)
    keep_searching = len(tree.nodes()) * 10
    while (
        j == i.parent_node
        or j.parent_node == i.parent_node 
        or j in i.preorder_iter()
    ) and keep_searching:
        i = np.random.choice(non_root_nodes)
        j = np.random.choice(non_root_nodes)
        keep_searching -= 1
    if not keep_searching:
        found = False
    else:
        found = True
    return i, j, found


def add_counts_entry(new_internal, node_1, node_2, counts):
    counts[new_internal.label] = {
        "counts": counts[node_1.label]["counts"] + counts[node_2.label]["counts"],
        "size": counts[node_1.label]["size"] + counts[node_2.label]["size"],
    }
    counts[new_internal.label]["entropy"] = compute_cluster_entropy(new_internal, counts)
    return counts


def create_node():
    global LABEL_COUNT
    nd = dendropy.Node()
    nd._set_label(str(LABEL_COUNT))
    LABEL_COUNT += 1
    return nd


def insert_internal_node(node, parent, children):
    node.set_child_nodes(children)
    parent.add_child(node)


def create_seqs_matrix(seqs, SEQ_LEN):
    """
    ONE-HOT ENCODING
    seqs is a dict mapping id:sequence
    return N x M x 5 tensor where each row is a M x 5 one-hot encoding of an M-length sequence
    """
    seqs_m = np.zeros(shape=(len(seqs), SEQ_LEN, 5))
    seqs_index = {}
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


def get_lca(tree, i, j):
    """
    returns lowest common ancestor of nodes i and j
    """
    i_leaves = [k.taxon.label for k in i.leaf_iter()]
    j_leaves = [k.taxon.label for k in j.leaf_iter()]
    leaf_taxa = set([taxa for taxa in itertools.chain(i_leaves, j_leaves)])
    lca = tree.mrca(taxon_labels=(leaf_taxa))
    return lca


def assign_internal_node_labels(tree):
    label = 0
    for i in tree.nodes():
        i._set_label(str(label))
        label += 1
    return label


def finish_node(node):
    if node.taxon is None:
        return node
    curr_label = node.taxon.label
    try:
        taxon, label = curr_label.split(" ")
    except ValueError:
        taxon, label = curr_label.split("_")

    node.taxon.label = taxon
    node.label = label
    return node


def read_args():
    try:
        return sys.argv[5], sys.argv[2:5]
    except IndexError:
        return sys.argv[1:5]


def tree_copy(tree):
    tree_str = tree.as_string(
        schema="newick",
        suppress_leaf_node_labels=False,
        suppress_internal_node_labels=False,
        suppress_edge_lengths=True,
        suppress_rooting=True
    )
    tree_new = dendropy.Tree.get(
        data=tree_str,
        schema="newick",
        suppress_edge_lengths=True,
        preserve_underscores=True,
        finish_node_fn=finish_node,
        rooting='force-rooted'
    )
    return tree_new


def read_tree(path):
    tree = dendropy.Tree.get(
        path=path,
        schema="newick",
        suppress_edge_lengths=True,
        preserve_underscores=True,
        rooting='force-rooted')
    return tree


def main():
    global LABEL_COUNT
    global SEQ_LEN

    nwk, fasta, out_nwk, num_iter = read_args()

    tree_0 = read_tree(nwk)
    LABEL_COUNT = assign_internal_node_labels(tree_0)

    seqs_d = create_seqs_dict(fasta)  
    SEQ_LEN = len(list(seqs_d.values())[0])
    seqs_m, seqs_index = create_seqs_matrix(seqs_d, SEQ_LEN)
    counts = create_counts_matrices(tree_0, seqs_m, seqs_index, SEQ_LEN)
    starting_entropy = compute_initial_tree_entropy(tree_0, counts)

    current_entropy = starting_entropy
    total_moves = 0

    for k in tqdm(range(int(num_iter))):
        # TODO: Retain copies between rejected moves
        tree = tree_copy(tree_0)
        new_counts = copy.deepcopy(counts)

        i, j, found = find_valid_move(tree_0)
        if not found:
            print("no suitable moves found.")
            break

        node_i = tree.find_node_with_label(i.label)
        node_j = tree.find_node_with_label(j.label)
        j_parent = node_j.parent_node

        new_counts, rerooted = move(tree, node_i, node_j, new_counts)

        if rerooted:
            new_counts = update_counts_after_reroot(new_counts, node_i, j_parent, tree.seed_node)
        else:
            lca = get_lca(tree, i, j)
            i_grandparent = i.parent_node.parent_node
            new_counts = update_counts(new_counts, node_i, i_grandparent, j_parent, lca)

        entropy = compute_tree_entropy(tree, new_counts)
        if entropy < current_entropy:
            tree_0 = tree
            counts = new_counts
            current_entropy = entropy
            total_moves += 1
            tree_0.write(path=out_nwk,schema="newick",
                suppress_internal_node_labels=True,suppress_edge_lengths=True, suppress_rooting=True)
            

    print(f" Entropy:   {starting_entropy:.2f} --> {current_entropy:.2f} after {total_moves}/{num_iter} accepted moves.")
    #Save final tree
    tree_0.write(
        path=out_nwk,
        schema="newick",
        suppress_internal_node_labels=True,
        suppress_edge_lengths=True,
        suppress_rooting=True,
    )


if __name__ == "__main__":
    main()
