import dendropy
from dendropy.model.parsimony import parsimony_score
import sys
import os

# nwk = r'.\sphere_binary\standardized\standardized_binary_sphere_gisaid_1a.nwk'
# fasta = r'.\sequences_2021-03-04_08-34_trimmed_50_l_r_no_amb_aligned_before_march_5_3757_total_sorted.fasta'


def compute_parsimony_score(nwk, fasta):
    """
    nwk - PATH to nwk file
    fasta - PATH to fasta file
    """

    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(
        path=nwk, schema="newick", taxon_namespace=taxa, preserve_underscores=True, suppress_edge_lengths=True
    )
    chars = dendropy.DnaCharacterMatrix.get(
        path=fasta, schema="fasta", taxon_namespace=taxa
    )

    return parsimony_score(tree, chars)
