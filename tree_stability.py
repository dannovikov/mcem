"""This program tests the stability of a phylogeny tree by randomly mutating a small percentage of nucleotides,
recomputing the tree, and comparing the differences between before and after.

First, a fasta file of aligned sequences will be perturbed, changing nucleotides randomly. 

Once the perturbed file is created, 

For each method being tested:
    construct two trees, one on original data and one on perturbed data
    (Note: for random methods, we need to find the tree mode in 100 runs.)
    Compute distance from tree_original to tree_new

"""
import sys
import subprocess
import random
import json
import os
import pandas as pd
from tqdm import tqdm

import dendropy
from dendropy.calculate import treecompare
from utils.parsimony_score import compute_parsimony_score

MC_DIR = "/home/dnovikov1/dan/entropy"  # Monte carlo file
SPHERE_DIR = "/home/dnovikov1/dan/sphere"  # Sphere
RAXML_DIR = "/home/dnovikov1/dan/tools/standard-RAxML"
UTILS_DIR = "/home/dnovikov1/dan/entropy/utils"
WORK_DIR = "/home/dnovikov1/dan/entropy/utils/test_outputs"  # Working directory
OUT_DIR = "/home/dnovikov1/dan/entropy/utils/results"  # Output directory


def create_seqs_dict(fasta):
    """
    From a fasta file, the function builds a dictionary {label:sequence}.
    In fasta format, the even # lines are labels that start with ">",
    the odd # lines are sequences.
    """
    seqs = {}
    with open(fasta, "r") as f:
        last_label = ""
        for line, text in enumerate(f):
            if line % 2 == 0:  # is label
                # Create new entry in dictionary for upcoming sequence
                seqs[text.strip()[1:]] = ""
                last_label = text.strip()[1:]
            else:
                # Add sequence to newly created entry
                seqs[last_label] = text.strip()
    return seqs


def perturb_sequences(seqs, p=0.05, keep_same_nucl_allowed=False):
    """
    This function randomly perturbs nucleotides -> changing to any other nucleotide
    with equal probability, even to the same nucl if keep_same_nucl_allowed=True.

    parameters:
        seqs: dictionary mapping {labels:sequences}
        p:  probability that nucleotide changes
        keep_same_nucl_allowed:  allow changing nucleotide to keep same state

    returns:
        petrubred_seqs: dictionary mapping {label: perturbed sequences}
    """
    perturbed_seqs = {}
    for s_id in seqs:
        seq = seqs[s_id]
        for i in range(len(seq)):
            nucl = seq[i]
            if random.random() < p:
                options = ["A", "C", "G", "T"]
                if not keep_same_nucl_allowed:
                    options.remove(nucl)
                seq = seq[:i] + random.choice(options) + seq[i+1:]
        perturbed_seqs[s_id] = seq
    return perturbed_seqs


def get_pairwise_distances(computed_trees):
    # Compute pairwise distances
    results = {}
    taxa = dendropy.TaxonNamespace()
    for i, tree in tqdm(enumerate(computed_trees)):
        for j, other_tree in enumerate(computed_trees):
            if i != j:
                t0 = dendropy.Tree.get(
                    path=tree, schema="newick", taxon_namespace=taxa
                )
                t1 = dendropy.Tree.get(
                    path=other_tree, schema="newick", taxon_namespace=taxa
                )
                d = treecompare.symmetric_difference(t0, t1)
                if i in results:
                    results[i][j] = d
                else:
                    results[i] = {j: d}
    return results


def write_seqs(seqs, out_file):
    with open(out_file, "w") as f:
        for s_id in seqs:
            f.write(f">{s_id}\n{seqs[s_id]}\n")


def cat(fasta, ref):
    # Combine fasta and reference sequence
    with open(fasta, "r") as f:
        with open(ref, "r") as r:
            with open(f"{WORK_DIR}/combined.fasta", "w") as c:
                for line, text in enumerate(r):
                    if text != "\n":
                        c.write(text)
                c.write("\n")
                for line, text in enumerate(f):
                    if text != "\n":
                        c.write(text)
    return f"{WORK_DIR}/combined.fasta"


def matrix_dict_to_df(results):
    # results is a dict of {tree : {other_tree: distance}}
    # returns a pandas dataframe whos values are distances between the trees on the axes
    df = pd.DataFrame()
    for tree in results:
        df.loc[tree, tree] = 0
        for other_tree in results[tree]:
            df.loc[tree, other_tree] = results[tree][other_tree]
    return df


def clear_working_dir():
    os.system(f"rm {WORK_DIR}/*")


def run_sphere(fasta, custom_ref=None, perturbed=False, tree_index=None):
    if tree_index is None:
        raise ValueError("Please provide an index to name the sphere output file")

    # Construct commands to run sphere and produce binary newick
    if perturbed:
        nwk_path = f"{WORK_DIR}/sphere_binary_p_{tree_index}.nwk"
    else:
        nwk_path = f"{WORK_DIR}/sphere_binary_{tree_index}.nwk"

    if custom_ref:
        #print("Using reference: ", custom_ref)
        sphere_cmd = (
            f"java -jar {SPHERE_DIR}/sphere/sphere.jar "
            f"-i {fasta} "
            f"-r {custom_ref} "
            f"-e {WORK_DIR}/sphere_edges.txt "
            f"-v {WORK_DIR}/sphere_nodes.txt "
            f"-s {WORK_DIR}/sphere_seqs.txt"
        )
    else:
        sphere_cmd = (
            f"java -jar {SPHERE_DIR}/sphere/sphere.jar "
            f"-i {fasta} "
            f"-r {SPHERE_DIR}/sample_inputs/ref.fas "
            f"-e {WORK_DIR}/sphere_edges.txt "
            f"-v {WORK_DIR}/sphere_nodes.txt "
            f"-s {WORK_DIR}/sphere_seqs.txt"
        )

    sphere_to_nwk_cmd = (
        f"python {UTILS_DIR}/sphere_to_newick.py "
        f"{WORK_DIR}/sphere_nodes.txt "
        f"{WORK_DIR}/sphere_edges.txt "
        f"{WORK_DIR}/sphere_seqs.txt "
        f"{WORK_DIR}/sphere_output.nwk"
    )

    sphere_to_binary_nwk_cmd = (
        f"python {UTILS_DIR}/binary.py "
        f"{WORK_DIR}/sphere_output.nwk "
        f"{nwk_path}"
    )

    # Run commands
    subprocess.run(sphere_cmd.split(" "), stdout=subprocess.DEVNULL)
    subprocess.run(sphere_to_nwk_cmd.split(" "), stdout=subprocess.DEVNULL)
    subprocess.run(sphere_to_binary_nwk_cmd.split(" "), stdout=subprocess.DEVNULL)

    # Read output file to return newick as string
    nwk_file = open(nwk_path, "r")
    nwk = nwk_file.read()
    nwk_file.close()
    return nwk, nwk_path


def run_monte_carlo(nwk, fasta, index, num_iter=10):
    nwk_out_path = f"{WORK_DIR}/mc_output_{index}.nwk"
    mc_cmd = (
        f"python {MC_DIR}/monte_carlo_entropy.py "
        f"{nwk} {fasta} {nwk_out_path} {num_iter}"
    )
    subprocess.run(mc_cmd.split(" "), stdout=subprocess.DEVNULL)
    with open(nwk_out_path, "r") as nwk_file:
        nwk = nwk_file.read()
    return nwk, nwk_out_path


def run_raxml(fasta, index):
    rax_outpath = f"{WORK_DIR}/RAxML_bestTree.{index}"
    nwk_outpath = f"{WORK_DIR}/RAxML_bestTree_standardized.{index}"
    raxml_cmd = f"{RAXML_DIR}/raxmlHPC -m GTRCAT -V -s {fasta} -n {index} -w {WORK_DIR} -p 12345 "
    subprocess.run(raxml_cmd.split(" "), stdout=subprocess.DEVNULL)
    standardize_raxml_cmd = (
        f"python {UTILS_DIR}/binary.py "
        f"{rax_outpath} "
        f"{nwk_outpath}"
    )
    subprocess.run(standardize_raxml_cmd.split(" "), stdout=subprocess.DEVNULL)
    with open(nwk_outpath, "r") as nwk_file:
        nwk = nwk_file.read()
    return nwk, nwk_outpath


def run_experiment(fasta, ref, method, tree_index, perturbed, mc_iter=10):
    print(f"Running experiment {fasta=}, {ref=}, {perturbed=}, {method=}, {tree_index=}")
    fasta_with_ref = cat(fasta, ref)

    if perturbed:
        mc_out = f'{method}_p_{tree_index}'
        raxml_out = f'tree_p_{tree_index}'
    else:
        mc_out = f'{method}_{tree_index}'
        raxml_out = f'tree_{tree_index}'

    if method.lower() == 'sphere':
        nwk, nwk_path = run_sphere(fasta, ref, perturbed, tree_index)
    elif method.lower() == 'raxml':
        nwk, nwk_path = run_raxml(fasta_with_ref, raxml_out)
    else:
        raise ValueError(f'{method} is an invalid choice for method. Please use one of ["sphere", "raxml"]')

    orig_pscore = compute_parsimony_score(nwk_path, fasta_with_ref)
    mc_nwk, mc_nwk_path = run_monte_carlo(nwk_path, fasta_with_ref, mc_out, num_iter=mc_iter)
    mc_pscore = compute_parsimony_score(mc_nwk_path, fasta_with_ref)
    return nwk_path, mc_nwk_path, orig_pscore, mc_pscore


def main():
    """
    Required arguments:
        fasta:  path to fasta file WITHOUT REFERENCE SEQUENCE
    Optional arguments:
        ref:  path to reference fasta file [default: EPI_ISL_405124]
    """
    num_trees = 10
    mc_iter = 25
    p = 0.01

    fasta_path = sys.argv[1]
    try:
        ref_path = sys.argv[2]
    except IndexError:
        ref_path = f"{SPHERE_DIR}/sample_inputs/ref.fas"

    tree_ids = {
        0: "sphere",
        1: "sphere_mc",
        2: "raxml",
        3: "raxml_mc",
        4: "sphere_pert",
        5: "sphere_pert_mc",
        6: "raxml_pert",
        7: "raxml_pert_mc"
    }

    seqs = create_seqs_dict(fasta_path)

    tree_index = 0
    running_results = {"dists": [], "pscores": []}
    while tree_index < num_trees:
        print(f"\ntree_index: {tree_index}")
        clear_working_dir()

        # Create perturbed sequences
        pert_seqs = perturb_sequences(seqs, p=p)
        pert_fasta_path = f"{WORK_DIR}/pert_seqs.fasta"
        write_seqs(pert_seqs, pert_fasta_path)

        # Run experiments
        computed_trees = []
        pscores = {}
        tree_id = 0
        for fasta in [fasta_path, pert_fasta_path]:
            perturbed = False if fasta == fasta_path else True
            for method in ['sphere', 'raxml']:
                nwk_path, mc_nwk_path, orig_pscore, mc_pscore = run_experiment(
                    fasta,
                    ref_path,
                    method,
                    tree_index,
                    perturbed,
                    mc_iter
                )
                print(f"{orig_pscore=}, {mc_pscore=}")
                computed_trees.append(nwk_path)
                computed_trees.append(mc_nwk_path)

                pscores[tree_ids[tree_id]] = orig_pscore
                tree_id += 1
                pscores[tree_ids[tree_id]] = mc_pscore
                tree_id += 1

        distances = get_pairwise_distances(computed_trees)
        running_results["dists"].append(matrix_dict_to_df(distances))
        running_results["pscores"].append(pscores)
        tree_index += 1

    final_results = sum(running_results["dists"]) / len(running_results["dists"])
    final_results.rename(index=tree_ids, columns=tree_ids, inplace=True)
    final_results.to_csv(f"{OUT_DIR}/distances.csv")

    pscore_df = pd.DataFrame(running_results["pscores"]).T
    pscore_df.to_csv(f"{OUT_DIR}/pscores.csv")

    print(final_results)
    print(pscore_df)


if __name__ == "__main__":
    main()