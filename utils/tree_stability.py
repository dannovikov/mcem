"""
This program tests the stability of a phylogeny tree by randomly mutating a small percentage of nucleotides, recomputing the tree, and comparing the differences between before and after.

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
from parsimony_score import compute_parsimony_score

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


def perturb_sequences(seqs, p=0.01, keep_same_nucl_allowed=True):
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
                seq = seq[:i] + random.choice(options) + seq[i + 1 :]
        perturbed_seqs[s_id] = seq
    return perturbed_seqs


def run_sphere(fasta, custom_ref=None):
    # Construct commands to run sphere and produce binary newick
    if custom_ref:
        print("Using referece: ", custom_ref)
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
        f"{WORK_DIR}/sphere_output_binary.nwk"
    )

    # Run commands
    subprocess.run(sphere_cmd.split(" "), stdout=subprocess.DEVNULL)
    subprocess.run(sphere_to_nwk_cmd.split(" "), stdout=subprocess.DEVNULL)
    subprocess.run(sphere_to_binary_nwk_cmd.split(" "), stdout=subprocess.DEVNULL)

    # Read output file to return newick as string
    nwk_path = f"{WORK_DIR}/sphere_output_binary.nwk"
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
    outpath = f"{WORK_DIR}/RAxML_bestTree.{index}"
    delete_old_files = False
    for f in os.listdir(WORK_DIR):
        if "RAxML_" in f:
            delete_old_files = True
            break
    if delete_old_files:
        os.system(f"rm {WORK_DIR}/RAxML_*")

    raxml_cmd = f"{RAXML_DIR}/raxmlHPC -m GTRCAT -V -s {fasta} -n {index} -w {WORK_DIR} -p 12345 "

    subprocess.run(raxml_cmd.split(" "), stdout=subprocess.DEVNULL)
    with open(outpath, "r") as nwk_file:
        nwk = nwk_file.read()
    return nwk, outpath


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


def main():
    """
    Required arguments:
        fasta:  path to fasta file WITHOUT REFERENCE SEQUENCE
    Optional arguments:
        ref:  path to reference fasta file [default: EPI_ISL_405124]


    Future:
    [X] p:  probability that nucleotide changes
    [X] keep_same_nucl_allowed:  allow changing nucleotide to keep same state

    """

    fasta = sys.argv[1]
    try:
        ref = sys.argv[2]
    except:
        ref = f"{SPHERE_DIR}/sample_inputs/ref.fas"

    num_iter = 20

    seqs = create_seqs_dict(fasta)

    running_results = {"dists": [], "pscores": []}
    tree_index = 0

    while tree_index < num_iter:
        print(f"\ntree_index: {tree_index}")

        clear_working_dir()

        pscores = {}

        # Create perturbed sequences
        pert_seqs = perturb_sequences(seqs)
        pert_fasta_path = "pert_seqs.fasta"
        write_seqs(pert_seqs, pert_fasta_path)

        # Concatenate fasta and reference
        fasta_with_ref = cat(fasta, ref)
        pert_fasta_with_ref = cat(pert_fasta_path, ref)

        # Run Sphere and MC on original sequences
        print("Running Sphere on original sequences")
        sphere_nwk, sphere_nwk_path = run_sphere(fasta, ref)
        pscores["sphere"] = compute_parsimony_score(sphere_nwk_path, fasta_with_ref)

        print("Running Monte Carlo on sphere output")
        sphere_mc_nwk, sphere_mc_nwk_path = run_monte_carlo(
            sphere_nwk_path, fasta_with_ref, f"orig_{tree_index}", num_iter=25
        )
        pscores["sphere_mc"] = compute_parsimony_score(
            sphere_mc_nwk_path, fasta_with_ref
        )

        print("[SPHERE] Parsimony of original:\t", pscores["sphere"])
        print("[SPHERE] Parsimony of MC min_ent:\t ", pscores["sphere_mc"])

        # Run Sphere and MC on perturbed sequences
        print("Running Sphere on perturbed sequences")
        sphere_pert_nwk, sphere_pert_nwk_path = run_sphere(pert_fasta_path, ref)
        pscores["sphere_pert"] = compute_parsimony_score(
            sphere_pert_nwk_path, pert_fasta_with_ref
        )

        print("Running Monte Carlo on sphere output")
        sphere_mc_pert_nwk, sphere_mc_pert_nwk_path = run_monte_carlo(
            sphere_pert_nwk_path, pert_fasta_with_ref, f"pert_{tree_index}", num_iter=25
        )
        pscores["sphere_mc_pert"] = compute_parsimony_score(
            sphere_mc_pert_nwk_path, pert_fasta_with_ref
        )

        print("[SPHERE] Parsimony of original (perturbed):\t", pscores["sphere_pert"])
        print(
            "[SPHERE] Parsimony of MC min_ent (perturbed):\t ",
            pscores["sphere_mc_pert"],
        )

        # Run RAxML and MC on original sequences
        print("Running RAxML on original sequences")
        raxml_nwk, raxml_nwk_path = run_raxml(fasta_with_ref, f"orig_{tree_index}")
        pscores["raxml"] = compute_parsimony_score(raxml_nwk_path, fasta_with_ref)

        print("Running Monte Carlo on RAxML output")
        raxml_mc_nwk, raxml_mc_nwk_path = run_monte_carlo(
            raxml_nwk_path, fasta_with_ref, f"orig_{tree_index}", num_iter=25
        )
        pscores["raxml_mc"] = compute_parsimony_score(raxml_mc_nwk_path, fasta_with_ref)

        print("[RAxML] Parsimony of original:\t", pscores["raxml"])
        print("[RAxML] Parsimony of MC min_ent:\t ", pscores["raxml_mc"])

        # Run RAxML and MC on perturbed sequences
        raxml_pert_nwk, raxml_pert_nwk_path = run_raxml(
            pert_fasta_with_ref, f"pert_{tree_index}"
        )
        pscores["raxml_pert"] = compute_parsimony_score(
            raxml_pert_nwk_path, pert_fasta_with_ref
        )
        raxml_mc_pert_nwk, raxml_mc_pert_nwk_path = run_monte_carlo(
            raxml_pert_nwk_path, pert_fasta_with_ref, f"pert_{tree_index}", num_iter=25
        )
        pscores["raxml_mc_pert"] = compute_parsimony_score(
            raxml_mc_pert_nwk_path, pert_fasta_with_ref
        )

        print("[RAxML] Parsimony of original (perturbed):\t", pscores["raxml_pert"])
        print(
            "[RAxML] Parsimony of MC min_ent (perturbed):\t ",
            pscores["raxml_mc_pert"],
        )

        computed_trees = []

        for tree in [
            sphere_nwk,
            sphere_mc_nwk,
            sphere_pert_nwk,
            sphere_mc_pert_nwk,
            raxml_nwk,
            raxml_mc_nwk,
            raxml_pert_nwk,
            raxml_mc_pert_nwk,
        ]:
            computed_trees.append(tree)

        # Compute pairwise distances
        results = {}
        taxa = dendropy.TaxonNamespace()
        for i, tree in tqdm(enumerate(computed_trees)):
            for j, other_tree in enumerate(computed_trees):
                if i != j:
                    t0 = dendropy.Tree.get_from_string(
                        tree, "newick", taxon_namespace=taxa
                    )
                    t1 = dendropy.Tree.get_from_string(
                        other_tree, "newick", taxon_namespace=taxa
                    )
                    d = treecompare.symmetric_difference(t0, t1)
                    if i in results:
                        results[i][j] = d
                    else:
                        results[i] = {j: d}

        # save results
        with open(f"{WORK_DIR}/results_{tree_index}.json", "w") as f:
            json.dump(results, f)

        running_results["dists"].append(matrix_dict_to_df(results))
        running_results["pscores"].append(pscores)

        tree_index += 1

    final_results = sum(running_results["dists"]) / len(running_results["dists"])
    tree_ids = {
        0: "sphere",
        1: "sphere_mc",
        2: "sphere_pert",
        3: "sphere_mc_pert",
        4: "raxml",
        5: "raxml_mc",
        6: "raxml_pert",
        7: "raxml_mc_pert",
    }
    final_results.rename(index=tree_ids, columns=tree_ids, inplace=True)
    final_results.to_csv(f"{OUT_DIR}/distances.csv")
    print(final_results)

    pscore_df = pd.DataFrame(running_results["pscores"]).T
    pscore_df.to_csv(f"{OUT_DIR}/pscores.csv")
    print(pscore_df)


if __name__ == "__main__":
    main()
