import dendropy
import sys

nwk = sys.argv[1]
taxa = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(path=nwk, schema='newick', taxon_namespace=taxa, preserve_underscores=True, suppress_edge_lengths=True)
tree.print_plot()
