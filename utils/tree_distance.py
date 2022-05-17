import dendropy
from dendropy.model.parsimony import fitch_down_pass, fitch_up_pass
import itertools

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
            s = list(ss)[0]             # if multistate set, choose first state
            nucl = str(alphabet[s])
            seq += nucl
        nd.seq = seq
        if seq_len is not None:
            assert len(seq) == seq_len

    return tree


def infer_internal_node_sequences(tree_path, fasta):

    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get_from_path(tree_path,"newick",taxon_namespace=taxa,preserve_underscores=True)
    data = dendropy.DnaCharacterMatrix.get(path=fasta, schema="fasta", taxon_namespace=taxa)

    taxon_state_sets_map = data.taxon_state_sets_map(gaps_as_missing=True)
    try:
        fitch_down_pass(tree.postorder_node_iter(), taxon_state_sets_map=taxon_state_sets_map)
    except ValueError as e:
        print("Please ensure tree bifurcates at every internal node.")
        raise e
    fitch_up_pass(tree.preorder_node_iter())

    tree = assign_seqs_to_nodes(tree)
    
    return tree
        

def hamming_distance(seq1, seq2):
    d = 0
    for i in range(len(seq1)):
        if 'N' in (seq1[i], seq2[i]): continue
        if seq1[i] != seq2[i]:
            d += 1
    return d


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


def distance_along_tree(tree_path, seq1, seq2, fasta):
    """
        tree = path to dendropy tree 
        node1, node2 = nodes in tree
        returns distance along tree between nodes
    """
    tree = infer_internal_node_sequences(tree_path, fasta)
    node1 = tree.find_node_with_taxon_label(seq1)
    node2 = tree.find_node_with_taxon_label(seq2)
    lca = get_lca(tree, node1, node2)
    d1 = distance_to_ancestor(node1, lca)
    d2 = distance_to_ancestor(node2, lca)
    print(f"Hamming distance: {hamming_distance(node1.seq, node2.seq)}")
    return d1 + d2
