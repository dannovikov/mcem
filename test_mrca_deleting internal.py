import dendropy
import copy 
from tqdm import tqdm
from numba import jit
import numpy as np
import itertools
import sys

def assign_internal_node_labels(tree):
    label = 0
    for i in tree.nodes():
        print(f'adding label {label} to node with descendants {[nd.taxon for nd in i.preorder_iter()]}')
        i._set_label(str(label))
        label += 1
        
    return label


nwk, fasta, num_iter = sys.argv[1:]

tree_0 = dendropy.Tree.get(path=nwk, schema="newick", suppress_edge_lengths=True, preserve_underscores=True)


label_counter = 0
label_counter = assign_internal_node_labels(tree_0)

tree = copy.deepcopy(tree_0)

print(tree)

lca = dendropy.Tree.mrca(tree, taxon_labels=['A', 'B'])
print(tree)
print(lca)

print(tree_0)
tree_0.print_plot()
