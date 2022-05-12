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


def assign_internal_node_labels(tree):
    label = 0
    for i in tree.nodes():
        i._set_label(str(label))
        label += 1
    return label


def create_counts_matrices(tree, seqs_m, seqs_index, SEQ_LEN):
    counts = {}
    # for i in tqdm(tree.nodes()):
    for i in tree.nodes():
        leaves = [k.taxon.label for k in i.leaf_iter()]
        count = np.zeros(shape=(SEQ_LEN, 5))
        for leaf in leaves:
            idx = seqs_index[leaf.replace(" ", "_")]
            seq = seqs_m[idx]
            count += seq
        counts[i.label] = {"counts": count, "size": len(leaves), "entropy": 0}
    return counts


def update_counts(counts, tree, i, i_parent, j_parent, lca):
    """
    After a move, this function is called to update
        1. the nucleotide counts
        2. the entropy
    of all affected ancestor clusters; those until the LCA of i and j

    Params:
        tree: dendropy Tree
        i: node i being moved
        i_parent: parent above where i was removed
        j_parent: destination parent
    """
    p = i_parent
    while p and p != lca:
        counts[p.label]["counts"] -= counts[i.label]["counts"]
        assert all(j >= 0 for i in counts[p.label]['counts'] for j in
                   i), f"Error in update counts subtracting {[x for x in counts[p.label]['counts'] if -1 in x]}, {p.label=}"
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


def get_lca(tree, i, j):
    """
    params:
        tree - dendropy tree
        i - node being moved
        j - destination node

    returns lowest common ancestor of i and j
    """
    i_leaves = [k.taxon.label for k in i.leaf_iter()]
    j_leaves = [k.taxon.label for k in j.leaf_iter()]
    leaf_taxa = set([taxa for taxa in itertools.chain(i_leaves, j_leaves)])
    try:
        lca = tree.mrca(taxon_labels=(leaf_taxa))
    except Exception as e:
        print("error in LCA")
        print(e)
        print(
            leaf_taxa,
            len(leaf_taxa),
            tree.taxon_namespace.get_taxa(labels=leaf_taxa),
            len(tree.taxon_namespace.get_taxa(labels=leaf_taxa)),
            [
                i
                for i in leaf_taxa
                if i not in tree.taxon_namespace.get_taxa(labels=leaf_taxa)
            ],
            len(tree.taxon_namespace),
        )
        exit()
    return lca


def compute_cluster_entropy(cluster, counts):
    col_entropy = 0
    size = counts[cluster.label]["size"]
    for i in range(4):  # {A,C,T,G}
        for pos in range(SEQ_LEN):
            count = counts[cluster.label]["counts"][pos, i]  # of A's @ pos i
            p_i = count / size
            assert 0 <= p_i <= 1, f"ERROR: mc encountered {p_i=}, {count=}, {pos=}, {SEQ_LEN=}"
            entropy = p_i * np.log2(p_i) if p_i != 0 else 0
            col_entropy += entropy
    return col_entropy / SEQ_LEN


def compute_tree_entropy(tree, counts):
    tree_entropy = 0
    # for node in tqdm(tree.nodes(), leave=False):
    for node in tree.nodes():
        entropy = counts[node.label]["entropy"]
        size = counts[node.label]["size"]
        tree_entropy += entropy * size
    return -1 * tree_entropy


def compute_tree_entropy_dividebyparent(tree, counts):
    tree_entropy = 0
    # for cluster in tqdm(tree.nodes(), leave=False):
    for cluster in tree.nodes():
        if cluster == tree.seed_node:
            continue
        entropy = counts[cluster.label]["entropy"]
        size = counts[cluster.label]["size"]
        parent = cluster.parent_node
        parent_size = counts[parent.label]["size"]
        # tree_entropy += (entropy * size) / parent_size
        # tree_entropy += entropy
    return -1 * tree_entropy


def compute_initial_tree_entropy(tree, counts):
    """
    This function computes entropy of every cluster (node) in the tree.
    It is designed to run be used once to obtain initial entropies each cluster.
    """
    tree_entropy = 0
    # for cluster in tqdm(tree.nodes(), leave=False):
    for cluster in tree.nodes():
        clust_entropy = compute_cluster_entropy(cluster, counts)
        tree_entropy += clust_entropy * counts[cluster.label]["size"]
        counts[cluster.label]["entropy"] = clust_entropy
    return -1 * tree_entropy


def move(tree, node_1, node_2, counts):
    # print(' moving', node_1.label, 'to', node_2.label)
    global label_counter
    target_parent = node_2.parent_node
    original_parent = node_1.parent_node
    try:
        sibling = node_1.sibling_nodes()[0]
    except:
        print("error, can't find sibling, node being moved is likely the root.")
        print(f"{node_1 == tree.seed_node}")
        tree.print_plot()
        print(tree, f"\n{node_1=}{len(node_1)=}\n{tree.seed_node=}")
        exit()

    tree.prune_subtree(node_1, suppress_unifurcations=False)
    tree.prune_subtree(node_2, suppress_unifurcations=False)

    new_internal = dendropy.Node()
    new_internal.set_child_nodes([node_1, node_2])
    target_parent.add_child(new_internal)

    # assign new label, add new entry to counts matrix
    new_internal._set_label(str(label_counter))
    counts[new_internal.label] = {
        "counts": counts[node_1.label]["counts"] + counts[node_2.label]["counts"],
        "size": counts[node_1.label]["size"] + counts[node_2.label]["size"],
    }
    counts[new_internal.label]["entropy"] = compute_cluster_entropy(
        new_internal, counts
    )
    if original_parent == tree.seed_node:
        tree.reroot_at_node(sibling, suppress_unifurcations=True)
        try:
            tree.prune_subtree(original_parent, suppress_unifurcations=False)
        except:
            pass
    else:
        ogp_sibling = original_parent.sibling_nodes()[0]
        original_parent.parent_node.set_child_nodes([sibling, ogp_sibling])

    # counts.pop(original_parent.label)
    label_counter += 1
    return counts, original_parent


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


label_counter = 0
SEQ_LEN = 0


def main():
    global label_counter
    global SEQ_LEN

    nwk, fasta, out_nwk, num_iter = sys.argv[1:5]
    try:
        existing_tree = sys.argv[5]
    except:
        existing_tree = False
    if existing_tree:
        print("Using existing tree at", existing_tree)
        nwk = existing_tree

    history = []
    seqs_d = create_seqs_dict(fasta)  # {seq_id: seq}
    SEQ_LEN = len(list(seqs_d.values())[0])

    tree_0 = dendropy.Tree.get(
        path=nwk, schema="newick", suppress_edge_lengths=True, preserve_underscores=True
    )
    label_counter = assign_internal_node_labels(tree_0)

    # print("Creating matrix representaion of sequences... ")
    seqs_m, seqs_index = create_seqs_matrix(seqs_d, SEQ_LEN)

    # print("Creating counts matrices for all nodes... ")
    counts = create_counts_matrices(tree_0, seqs_m, seqs_index, SEQ_LEN)

    # print("Computing starting entropy...")
    current_entropy = compute_initial_tree_entropy(tree_0, counts)
    starting_entropy = current_entropy
    # print(f"Starting entropy = {starting_entropy}")

    # print("Running moves experiment: ")

    total_moves = 0
    # for k in tqdm(range(int(num_iter))):
    for k in range(int(num_iter)):
        non_root_nodes = tree_0.nodes(lambda x: x != tree_0.seed_node)

        # Select qualifying source node and dest node for move
        i = np.random.choice(non_root_nodes)
        j = np.random.choice(non_root_nodes)
        stop = len(tree_0.nodes()) * 10
        while stop != 0 and not (
                j not in i.preorder_iter()
                and j.parent_node != i.parent_node
                and j != i.parent_node
                # and i.parent_node != tree_0.seed_node
        ):
            i = np.random.choice(non_root_nodes)
            j = np.random.choice(non_root_nodes)
            stop -= 1

        if stop == 0:
            print("no suitable moves found.")
            break

        # Make the move and compute entropy
        # This is where we used to make deep copies of the tree
        tree_0_str = tree_0.as_string(
            schema="newick",
            suppress_leaf_node_labels=False,
            suppress_internal_node_labels=False,
            suppress_edge_lengths=True,
        )
        tree = dendropy.Tree.get(
            data=tree_0_str,
            schema="newick",
            suppress_edge_lengths=True,
            preserve_underscores=True,
            finish_node_fn=finish_node,
        )
        tree_mod = dendropy.Tree.get(
            data=tree_0_str,
            schema="newick",
            suppress_edge_lengths=True,
            preserve_underscores=True,
            finish_node_fn=finish_node,
        )

        node_i = tree.find_node_with_label(i._get_label())
        node_j = tree.find_node_with_label(j._get_label())

        # get grandparent of source node b.c. parent will become a singleton & collapse w/ grandparent
        i_parent = (node_i.parent_node.parent_node)
        j_parent = node_j.parent_node

        new_counts = copy.deepcopy(counts)
        lca = get_lca(tree_mod, i, j)
        new_counts, original_parent = move(tree, node_i, node_j, new_counts)
        new_counts = update_counts(new_counts, tree, node_i, i_parent, j_parent, lca)
        entropy = compute_tree_entropy(tree, new_counts)

        if entropy < current_entropy:
            # print(f"### Overwriting tree: new tree has entropy {entropy} compared to prev {current_entropy}")
            counts = new_counts
            tree_0 = tree
            counts.pop(original_parent.label)
            history.append(current_entropy)
            current_entropy = entropy
            tree_0.write(
                path=out_nwk,
                schema="newick",
                suppress_internal_node_labels=True,
                suppress_edge_lengths=True,
            )
            total_moves += 1

    print(
        f" Entropy:   {starting_entropy:.2f} --> {current_entropy:.2f} after {total_moves}/{num_iter} accepted moves."
    )
    tree_0.write(
        path=out_nwk,
        schema="newick",
        suppress_internal_node_labels=True,
        suppress_edge_lengths=True,
    )
    # print(f"Tree saved to: {out_nwk}")
    history.append(current_entropy)
    plt.plot(history)
    plt.ylabel("entropy")
    plt.savefig("history_g1a_sphere_new.png")


if __name__ == "__main__":
    main()
