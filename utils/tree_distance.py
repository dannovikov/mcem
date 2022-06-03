"""
This program computes distances along a phylogenetic tree between two sequences.
It uses Fitch's algorithm to assign sequences to internal nodes, then for each
of the two sequences, it computes the sum of edge lengths to their common
ancestor. Edge lengths are the hamming distance between the nodes. 

Usage:

    Input: tree_path, fasta_path 
    
    tree, taxa = read_tree(tree_path)
    tree = infer_internal_node_sequences(tree, taxa, fasta_path)
    d = average_distance_along_tree(tree, fasta_[ath])

    s1 = "EPI_ISL_402124"
    s2 = "EPI_ISL_650723"

    d = distance_along_tree(tree, s1, s2)

"""

from distutils.command import check
from sys import intern
import dendropy
from dendropy.model.parsimony import fitch_down_pass, fitch_up_pass
import itertools
from numpy import hamming
from tqdm import tqdm


def compute_avg_dist_along_tree(tree_path, fasta_path):
    tree, taxa = read_tree(tree_path)
    tree = infer_internal_node_sequences(tree, taxa, fasta_path)
    d,h,dv = average_distance_along_tree(tree, fasta_path)
    return d,h,dv


def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2), f"{len(seq1)} != {len(seq2)}, {seq1} != {seq2}"
    d = 0
    for i in range(len(seq1)):
        if 'N' in (seq1[i], seq2[i]): continue
        if seq1[i] != seq2[i]:
            d += 1
    return d

def hamming_distance_by_label(tree, l1, l2):
    seq1 = tree.find_node_with_taxon_label(l1).seq
    seq2 = tree.find_node_with_taxon_label(l2).seq
    return hamming_distance(seq1, seq2)

def get_lca(tree, i, j):
    """
    Returns lowest common ancestor of nodes i and j
    params:
        tree - dendropy tree
        i, j - nodes in tree

    """
    i_leaves = [k.taxon.label for k in i.leaf_iter()]
    j_leaves = [k.taxon.label for k in j.leaf_iter()]
    leaf_taxa = set([taxa for taxa in itertools.chain(i_leaves, j_leaves)])
    lca = tree.mrca(taxon_labels=(leaf_taxa))
    return lca


def assign_seqs_to_nodes(tree, seq_len=None):
    '''
    After fitch up and down passes, this function assigns
    node.seq to the sequence inferred by Fitch algorithm.

    It does this by converting the numerical states assigned by Fitch
    into strings in the normal {A,C,T,G} alphabet.
    '''
    alphabet = dendropy.DNA_STATE_ALPHABET
    for nd in tree:
        ssl = nd.state_sets         
        seq = ''
        for ss in ssl:
            s = list(ss)[0]             # For multistate sets, choose first state
            nucl = str(alphabet[s])
            seq += nucl
        nd.seq = seq
        if seq_len is not None:
            assert len(seq) == seq_len

    return tree
    

def extract_sequence_labels(fasta):
    seqs = []
    with open(fasta, 'r') as f:
        for i, line in enumerate(f):
            if i % 2 == 0:
                seqs.append(line.strip()[1:])
    return seqs


def read_tree(tree_path):
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get_from_path(tree_path,"newick",taxon_namespace=taxa,preserve_underscores=True)
    return tree, taxa


def ensure_internal_seqs_exist(tree):
        assert hasattr(tree.nodes()[0], 'seq'), \
        "all nodes, including internal, must have a .seq property with a DNA sequence. "\
        "Use 'tree = infer_internal_node_sequences(tree, taxon_namespace, fasta)'"


def infer_internal_node_sequences(tree, taxa, fasta):
    
    data = dendropy.DnaCharacterMatrix.get(path=fasta, schema="fasta", taxon_namespace=taxa)
    taxon_state_sets_map = data.taxon_state_sets_map(gaps_as_missing=True)

    try:
        fitch_down_pass(tree.postorder_node_iter(), taxon_state_sets_map=taxon_state_sets_map)
        fitch_up_pass(tree.preorder_node_iter())
    except ValueError as e:
        print("May have encountered something \
            other than a bifurcation at an internal node.")
        raise e

    tree = assign_seqs_to_nodes(tree)
    return tree
        

def distance_to_ancestor(node, ancestor):
    assert node in ancestor.preorder_iter()
    nd = node
    p  = node.parent_node
    d = hamming_distance(nd.seq, p.seq)
    while p!= ancestor:
        nd = p
        p = p.parent_node
        d += hamming_distance(nd.seq, p.seq)
    return d


def distance_along_tree(tree, seq1, seq2):
    """
        tree = Dendropy Tree on which infer_internal_node_sequences has been called
            meaning ALL nodes have a .seq property containing their DNA sequence.

        node1, node2 = nodes in tree
        --
        returns distance along the tree between nodes
    """
    ensure_internal_seqs_exist(tree)
    node1 = tree.find_node_with_taxon_label(seq1)
    node2 = tree.find_node_with_taxon_label(seq2)
    lca = get_lca(tree, node1, node2)
    # d1 = distance_to_ancestor(node1, lca)
    # d2 = distance_to_ancestor(node2, lca)
    d1 = hamming_distance(node1.seq, lca.seq)
    d2 = hamming_distance(node2.seq, lca.seq)
    return d1 + d2


def average_distance_along_tree(tree, fasta):
    ensure_internal_seqs_exist(tree)
    seqs = extract_sequence_labels(fasta)
    dists = []
    hdists = []
    deviations = []
    for i, s1 in enumerate(tqdm(seqs, leave=False)):
        for j, s2 in enumerate(seqs):
            if i != j:
                d = distance_along_tree(tree, s1, s2)
                dists.append(d)
                hdists.append(hamming_distance_by_label(tree, s1, s2))
                deviations.append((d - hdists[-1])**2)
    return sum(dists)/len(dists), hdists, sum(deviations)/len(deviations)
