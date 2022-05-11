# mcem
Monte Carlo Entropy Minimization 
----
`monte_carlo_entropy.py` takes a fasta file of aligned sequences and a phylogenetic tree of the sequences given in the newick format, and randomly moves subtrees around, accepting changes that reduce the total tree entropy. It returns a locally-minimum entropy phylogeny tree.

The definition of tree entropy is not settled, but currently we treat each internal node as a cluster of the seqeuences in its descendant leaf nodes, and define tree entropy as `sum(entropy(cluster) * size(cluster))` where size is the number of descendant leaves. Entropy of a cluster of aligned sequences is computed column-wise on the nucleotide frequencies.
----
`tree_stability.py` is an experiment to see how far away trees get from themselves when you randomly perturb a small percentage of nucleotides.
