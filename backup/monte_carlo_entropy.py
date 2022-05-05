import dendropy
import copy 
from tqdm import tqdm
from numba import jit
import numpy as np
import itertools
import sys
from matplotlib import pyplot as plt
import os

sys.setrecursionlimit(10000)


def create_seqs_dict(fasta):
    seqs = {}
    with open(fasta, "r") as f:
        last_label = ""
        for line, text in enumerate(f):
            if line % 2 == 0:  # is label
                seqs[text.strip()[1:]] = ""
                last_label = text.strip()[1:]
            else:
                seqs[last_label] = text.strip()
    return seqs


def create_seqs_matrix(seqs, SEQ_LEN):
    """
    seqs is a dict mapping id:sequence
    return N x M x 5 tensor where each row is a M x 5 one-hot encoding of an M-length sequence
    """
    seqs_m = np.zeros(shape=(len(seqs), SEQ_LEN, 5))
    seqs_index = {}
    for i, seq_id in enumerate(tqdm(seqs)):
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
                assert nucl == "-" or nucl == 'N'
                seqs_m[i][j][4] = 1
    return seqs_m, seqs_index


def assign_internal_node_labels(tree):
    label = 0
    for i in tree.nodes():
        #print(f'adding label {label} to node with descendants {[nd.taxon for nd in i.preorder_iter()]}')
        i._set_label(str(label))
        label += 1
        
    return label


def create_counts_matrices(tree, seqs_m, seqs_index, SEQ_LEN):
    counts = {}
    for i in tqdm(tree.nodes()):
        leaves = [k.taxon.label for k in i.leaf_iter()]
        count = np.zeros(shape=(SEQ_LEN, 5))
        for leaf in leaves:
            idx = seqs_index[leaf.replace(" ", "_")]
            seq = seqs_m[idx]
            count += seq
        counts[i.label] = {'counts': count, 'size': len(leaves), 'entropy': 0}
    return counts


def update_counts(counts, tree, i, i_parent, j_parent, lca):
    """
    After a move, this function updates 
        1. the nucleotide counts 
        2. the entropy
    of all affected ancestor clusters; those until the LCA of i and j

    Inputs:
    tree: dendropy Tree
    i: node i being moved
    i_parent: parent above where i was removed
    j_parent:  destination parent
    """

    p = i_parent
    while p and p != lca:
        counts[p.label]['counts'] -= counts[i.label]['counts']
        counts[p.label]['size'] -= counts[i.label]['size']
        counts[p.label]['entropy'] = compute_cluster_entropy(p, counts)
        p = p.parent_node

    p = j_parent
    while p and p != lca:
        counts[p.label]['counts'] += counts[i.label]['counts']
        counts[p.label]['size'] += counts[i.label]['size']
        counts[p.label]['entropy'] = compute_cluster_entropy(p, counts)
        p = p.parent_node
    
    return counts


def get_lca(tree, i, j):
    # node_i = tree.find_node_with_label(i._get_label())
    # node_j = tree.find_node_with_label(j._get_label())
    #print(f'{i=}, {j=}, {tree.seed_node}')
    #tree.print_plot()
    i_leaves = [k.taxon.label for k in i.leaf_iter()]
    j_leaves = [k.taxon.label for k in j.leaf_iter()]
    leaf_taxa = set([taxa for taxa in itertools.chain(i_leaves, j_leaves)])
    try:
        lca = tree.mrca(taxon_labels=(leaf_taxa))
    except Exception as e:
        print('error in LCA')
        print(e)
        print(leaf_taxa, len(leaf_taxa), tree.taxon_namespace.get_taxa(labels=leaf_taxa), len(tree.taxon_namespace.get_taxa(labels=leaf_taxa)), 
            [i for i in leaf_taxa if i not in tree.taxon_namespace.get_taxa(labels=leaf_taxa)], len(tree.taxon_namespace))  
        exit()
    return lca


def compute_cluster_entropy(cluster, counts):
    # entropy: for i in a,c,t,g: sum += counts[i]/len(cluster)
    col_entropy = 0
    size = counts[cluster.label]['size']
    for i in range(4):  # not 5, don't count dashes
        count = counts[cluster.label]['counts'][:,i]
        p_i = sum(count) / size
        entropy = p_i * np.log(p_i) if p_i != 0 else 0
        col_entropy += entropy
    return col_entropy / SEQ_LEN


# def compute_tree_entropy(tree, counts):
#     tree_entropy = 0
#     for cluster in tqdm(tree.nodes(), leave = False):
#         clust_entropy = compute_cluster_entropy(cluster, counts)
#         tree_entropy += clust_entropy * counts[cluster.label]['size']
#     return tree_entropy

def compute_tree_entropy(tree, counts):
    tree_entropy=0
    for node in tqdm(tree.nodes(), leave=False):
        entropy = counts[node.label]['entropy']
        size = counts[node.label]['size']
        tree_entropy += entropy * size
    return tree_entropy


def compute_initial_tree_entropy(tree, counts):
    '''
    This function computes entropy of every cluster (node) in the tree.
    It is designed to run once to obtain initial entropy of each cluster.

    '''
    tree_entropy = 0
    for cluster in tqdm(tree.nodes(), leave = False):
        clust_entropy = compute_cluster_entropy(cluster, counts)
        tree_entropy += clust_entropy * counts[cluster.label]['size']
        counts[cluster.label]['entropy'] = clust_entropy
    return tree_entropy




def move(tree, node_1, node_2, counts):
    #print('moving', node_1.label, 'to', node_2.label)
    global label_counter
    target_parent = node_2.parent_node
    original_parent = node_1.parent_node
    if original_parent == tree.seed_node:
        print('Skipping reroot move')
        return counts, original_parent
    counts = copy.deepcopy(counts)
    try:
        sibling = node_1.sibling_nodes()[0]
    except:
        print('error finding sibling, means node_1 is root')
        tree.print_plot()
        print(tree,f'\n{node_1=}{len(node_1)=}\n{tree.seed_node=}')
        exit()


    tree.prune_subtree(node_1, suppress_unifurcations=False) 
    tree.prune_subtree(node_2, suppress_unifurcations=False)

    new_internal = dendropy.Node()
    new_internal.set_child_nodes([node_1, node_2])
    target_parent.add_child(new_internal)

    #assign new label, and add its entry to counts matrix
    new_internal._set_label(str(label_counter))
    counts[new_internal.label] = {'counts': counts[node_1.label]['counts'] + counts[node_2.label]['counts'], 'size': counts[node_1.label]['size'] + counts[node_2.label]['size']}
    counts[new_internal.label]['entropy'] = compute_cluster_entropy(new_internal, counts)
    # if original_parent == tree.seed_node:
    #     print('REROOT: original parent', original_parent, 'rerooting at', sibling)
    #     tree.reroot_at_node(sibling, suppress_unifurcations=True)
    #     tree.prune_subtree(original_parent,suppress_unifurcations=False)
    # else:
    ogp_sibling = original_parent.sibling_nodes()[0]
    original_parent.parent_node.set_child_nodes([sibling, ogp_sibling])

    #counts.pop(original_parent.label)
    label_counter += 1
    return counts, original_parent


def dict_print(d):
    for k,v in d.items():
        print(k,v)

def move_and_compute_entropy(tree_0, i, j, counts):
    """
    i, j are nodes in t0
    """

    #print('tree_0 before cloning=', tree_0)
    #tree = tree_0.clone(depth=2)
    tree = copy.deepcopy(tree_0)
    tree_mod = copy.deepcopy(tree_0)
    #lca = tree.find_node_with_label(lca_label)
    # for new tree
    node_i = tree.find_node_with_label(i._get_label())
    node_j = tree.find_node_with_label(j._get_label())
    #print(f'\nMoving {i=}, \n{node_i=}, \n to {j=}, \n{node_j=}\n\n')
    i_parent = node_i.parent_node.parent_node #getting grandparent because parent will become a singleton and collapse w/ grandparent
    j_parent = node_j.parent_node

    
    #print('tree_0 before lca', tree_0)
    lca = get_lca(tree_mod, i, j)
    #print('after lca', tree_0, tree)

    #print('tree_0 before move=', tree_0)
    counts, original_parent = move(tree, node_i, node_j, counts)
    #print('tree_0 before move=', tree_0)
    #print('counts after move:')
    #dict_print(counts)
    counts = update_counts(counts, tree, node_i, i_parent, j_parent, lca)
    #print('counts after update counts:')
    #dict_print(counts)
    #print('tree_0 before cte=', tree_0)
    # entropy = newickEntropy.compute_entropy(tree=tree, seqs=seqs)
    entropy = compute_tree_entropy(tree, counts)
    #print('tree_0 after cte=', tree_0)
    return tree, entropy, counts 


label_counter = 0
SEQ_LEN = 0

def main():
    global label_counter
    global SEQ_LEN

    nwk, fasta, out_nwk, num_iter = sys.argv[1:5]
    try: existing_tree = sys.argv[5]
    except: existing_tree = False

    if existing_tree:
        print('Using existing tree at', existing_tree)
        nwk = existing_tree

    history = []
    seqs_d = create_seqs_dict(fasta)  # {seq_id: seq}
    SEQ_LEN = len(list(seqs_d.values())[0])

    tree_0 = dendropy.Tree.get(path=nwk, schema="newick", suppress_edge_lengths=True, preserve_underscores=True)
    label_counter = assign_internal_node_labels(tree_0)

    print("Creating matrix representaion of sequences... ")
    seqs_m, seqs_index = create_seqs_matrix(seqs_d, SEQ_LEN)


    print("Creating counts matrices for all nodes... ")
    counts = create_counts_matrices(tree_0, seqs_m, seqs_index, SEQ_LEN)

    print("Computing starting entropy...")
    current_entropy = compute_initial_tree_entropy(tree_0, counts)
    print(f"Starting entropy = {current_entropy}")

    print("Running moves experiment: ")
    for k in tqdm(range(int(num_iter))):
        #print('Iter: ', k)
        non_root_nodes = tree_0.nodes(lambda x: x != tree_0.seed_node)
        
        #Choose source and target for move
        i = np.random.choice(non_root_nodes)
        j = np.random.choice(non_root_nodes)

        stop = len(tree_0.nodes())* 10
        while stop!=0 and not (
            j not in i.preorder_iter()
            and j.parent_node != i.parent_node
            and j != i.parent_node
        ):
            i = np.random.choice(non_root_nodes)
            j = np.random.choice(non_root_nodes)
            stop -=1 

        if not stop: 
            print('no suitable moves found.')
            break
        
        #TODO: Memory leak. We add the new internal node's counts to the counts table even if we don't keep the move
        #Solution: seperate the move from the counts update. but then how do we compuote entropy? keep a serpeate counts table?
        # print('tree_0=', tree_0, '')
        #tree, entropy, counts = move_and_compute_entropy(tree_0, i, j, counts)

        #Make the move and compute entropy
        tree = copy.deepcopy(tree_0)
        tree_mod = copy.deepcopy(tree_0)
        node_i = tree.find_node_with_label(i._get_label())
        node_j = tree.find_node_with_label(j._get_label())
        i_parent = node_i.parent_node.parent_node #getting grandparent because parent will become a singleton and collapse w/ grandparent
        j_parent = node_j.parent_node
        lca = get_lca(tree_mod, i, j)
        new_counts, original_parent = move(tree, node_i, node_j, counts)
        new_counts = update_counts(new_counts, tree, node_i, i_parent, j_parent, lca)
        entropy = compute_tree_entropy(tree, new_counts)

        # print('tree_0 after mace method=', tree_0, '')
        # print('\n\ntree returned by mace=', tree, '')
        # print('counts after mace='), counts, '\n\n'
        if entropy < current_entropy:
            #print(f'##############overwriting tree. new tree has entropy {entropy} compared to prev {current_entropy}')
            counts = new_counts
            tree_0 = tree
            counts.pop(original_parent.label)
            #tree_0.print_plot()
            #print(tree_0)
            history.append(current_entropy)
            current_entropy = entropy
            #tqdm.write(f'Modified entropy: {entropy}')
            tree_0.write(
                path = out_nwk,
                schema= 'newick',
                suppress_internal_node_labels=True,
                suppress_edge_lengths=True
            ) 
        else:
            #print(f'##############SKIPPING tree with entropy {entropy} compared to prev {current_entropy}')
            pass
    print(current_entropy)
    history.append(current_entropy)
    plt.plot(history)
    plt.ylabel('entropy')
    plt.savefig('history.png')



if __name__ == "__main__":
    main()








# def make_clusters(tree, seqs):
#     """
#     Reads dendropy tree object
#     returns clusters (descendant set for each internal node) of seq IDs
#     """
#     clusters = {}
#     for node in tree:
#         if node.is_internal():#and node != tree.seed_node:
#             cluster_ids = set()
#             for child in node.leaf_iter(lambda x: x!=node):
#                 cluster_ids.add(str(child.taxon).replace(' ', '_'))
#             if 'None' in cluster_ids:
#                 cluster_ids.remove('None')
#             cluster = []
#             for seq_id in cluster_ids:
#                 cluster.append(seqs[seq_id[1:-1]])
#             clusters[node] = cluster
#     return clusters

# Here, let's loop and random choose a source and target cluster
# make a move and recompute entropy, decide keep, revert
# also let's create a table mapping cluster labels to counts matrix and size


# Create all tasks to run in multiproc pool
# tasks = []
# print("Creating tasks...")
# seqs = newickEntropy.make_seqs_dict(fasta)
# for i in non_root_nodes:
#     for j in non_root_nodes:
#         if j not in i.preorder_iter() and j.parent_node != i.parent_node and j != i.parent_node:
#             tasks.append((tree_0, i, j, seqs))

# #Create mp pool and run move_entropy tasks
# print("Starting MP pool...")
# pool = mp.Pool(parallel_num_procs)
# entropies = pool.starmap(move_entropy, tqdm(tasks, total=len(tasks)), chunksize=1)
# print(entropies)
# pool.close()
# # pool.join()


# num_procs += 1
# proc = mp.Process(move_entropy(
#         tree_0,
#         i,
#         j,
#         entropies,
#         f'raw_tree_{i._get_label()}_{j._get_label()}.nwk',
#         f'tree_{i._get_label()}_{j._get_label()}.nwk' ))
# proc.start()
# procs.append(proc)
# print('starting process...')
# if num_procs ==4:# mp.cpu_count():
#     for p in procs:
#         if p.is_alive():
#             p.join()
#     procs.clear()
#     num_procs = 0
#     print('queue cleared')


# output tree to newick file for entropy program.
# instead of writing to file, let's do everything internally.
# Let's
#  import newickEntropy
# tree.write(
#             path = raw_tree_file,
#             schema= 'newick',
#             suppress_internal_node_labels=True,
#             suppress_edge_lengths=True
#         )
# with open(raw_tree_file, 'r') as f:
#     with open(tree_file, 'w') as o:
#         text = f.readlines()[0]
#         text=text.replace(',,', ',')
#         o.write(text[5:])
# run_entropy = f"python newickEntropy.py {tree_file}".split(' ')
# result=subprocess.check_output(run_entropy)
# entropy = float(str(result.strip())[2:-2])
# print(entropy)

# print('after move original:', tree_0)
# print('after move cloned:', tree, '\n')
