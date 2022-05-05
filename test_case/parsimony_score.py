import dendropy
from dendropy.model.parsimony import parsimony_score
import sys
import os

nwk,fasta= sys.argv[1:]
taxa = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(path=nwk, schema='newick', taxon_namespace=taxa, preserve_underscores=True)
chars = dendropy.DnaCharacterMatrix.get(path=fasta, schema='fasta', taxon_namespace=taxa)

print(parsimony_score(tree, chars))
