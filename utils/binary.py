"""
Input: Newick with polytomies/multifurcations. (Like from Sphere)
Output: Binary newick; polytomies resolved with arbitrary dichotomic structure

"""

from ete3 import Tree
import sys

# t = Tree("gisaid_1a_sphere.nwk")
t = Tree(sys.argv[1])
t.resolve_polytomy(recursive=True)
t.standardize(delete_orphan=False)  # deletes internal nodes with 1 child
t.write(outfile=sys.argv[2], format=9)  # 9: Leaf labels only
