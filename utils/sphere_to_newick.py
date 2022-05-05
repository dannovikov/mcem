import networkx as nx
import sys


def inorder_newick(G, root):
    """
    Inorder recursive traversal and newick construction.
    Builds the children newick first,
    then builds and adds parent newicks as siblings.
    Collapsed nodes are handled as siblings.
    """
    if len(list(G.neighbors(root))) != 0:
        # process children
        output = ""
        output += "("
        for child in G.neighbors(root):
            output += inorder_newick(G, child)
            output += ","
        output = output[:-1]  # remove last comma
        output += ")"
        # process parent
        output += ","
        for seq in nodes[root]:
            output += seq
            output += ","
        output = output[:-1]  # remove last comma
        # return '(' + output + '))'
        return "(" + output + ")"
    else:
        # base case for children
        return ",".join(nodes[root])


def build_node_seqs_dict(nodes_file):
    nodes = {}
    with open(nodes_file, "r") as f:
        for index, line in enumerate(f.readlines()):
            if index != 0:
                seq, node = line.strip().split(",")
                node = int(node)
                if node in nodes:
                    nodes[node].append(seq)
                else:
                    nodes[node] = [seq]
    return nodes


def build_edge_list(edges_file):
    edges = []
    with open(edges_file, "r") as f:
        for index, line in enumerate(f.readlines()):
            if index != 0:
                src, trg, nmux = line.strip().split(",")
                edges.append((int(src), int(trg), int(nmux)))
    return edges


try:
    nodes_file = sys.argv[1]
    edges_file = sys.argv[2]
    seqs_file = sys.argv[3]
    out_file = sys.argv[4]
    nodes = build_node_seqs_dict(nodes_file)
    edges = build_edge_list(edges_file)

except:
    print(
        """
    nodes_file = sys.argv[1]
    edges_file = sys.argv[2]
    are missing. Running default test case.
        """
    )

    nodes = {
        0: ["a", "b", "c"],
        1: ["e"],
        2: ["d"],
        3: ["g", "h"],
        4: ["x"],
    }

    edges = [
        (0, 1, 1),
        (0, 2, 1),
        (2, 4, 1),
        (4, 3, 1),
    ]


G = nx.DiGraph()

for n in nodes.items():
    G.add_node(n[0])

for e in edges:
    G.add_edge(e[0], e[1])

nwk = inorder_newick(G, 0) + ";"

with open(out_file, "w") as f:
    f.write(nwk)
