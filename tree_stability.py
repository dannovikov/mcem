"""
Required arguments:
    1. fasta:  path to fasta file WITHOUT REFERENCE SEQUENCE
Optional arguments:
    2. ref:  path to reference fasta file [default: EPI_ISL_405124]
"""

import cProfile
import io
import pstats
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
from utils.seq_utils import create_seqs_dict
import utils.tree_distance as td

MC_DIR = "/home/dnovikov1/dan/entropy"  # Monte carlo file
SPHERE_DIR = "/home/dnovikov1/dan/sphere"  # Sphere
RAXML_DIR = "/home/dnovikov1/dan/tools/standard-RAxML"
UTILS_DIR = "/home/dnovikov1/dan/entropy/utils"
WORK_DIR = "/home/dnovikov1/dan/entropy/utils/test_outputs"  # Working directory
OUT_DIR = "/home/dnovikov1/dan/entropy/utils/results"  # Output directory


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
                options = ["A", "C", "G", "T", "N"]
                if not keep_same_nucl_allowed:
                    try:
                        options.remove(nucl)
                    except:
                        print(f"missing nucl = {nucl}")
                        raise ValueError
                seq = seq[:i] + random.choice(options) + seq[i + 1 :]
        perturbed_seqs[s_id] = seq
    return perturbed_seqs

from utils.perturb_sequences import perturb_sequences as perturb_seq_list
def perturb_sequences(seqs, p=0.05):
    s = list(seqs.values())
    s = perturb_seq_list(s, p)
    d = dict(zip(seqs.keys(), s))
    return d


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
                    c.write(f"{text.strip()}\n")
                for line, text in enumerate(f):
                    c.write(f"{text.strip()}\n")
    return f"{WORK_DIR}/combined.fasta"


def compute_pairwise_rf_distances(computed_trees):
    # Compute pairwise distances
    results = {}
    taxa = dendropy.TaxonNamespace()
    for i, tree in tqdm(enumerate(computed_trees)):
        for j, other_tree in enumerate(computed_trees):
            if i != j:
                t0 = dendropy.Tree.get(path=tree, schema="newick", taxon_namespace=taxa)
                t1 = dendropy.Tree.get(
                    path=other_tree, schema="newick", taxon_namespace=taxa
                )
                d = treecompare.symmetric_difference(t0, t1)
                if i in results:
                    results[i][j] = d
                else:
                    results[i] = {j: d}
    return results


def matrix_dict_to_df(results):
    # results is a dict of {tree : {other_tree: distance}}
    # returns a pandas dataframe with trees as the axes and distances as the values
    df = pd.DataFrame()
    for tree in results:
        df.loc[tree, tree] = 0
        for other_tree in results[tree]:
            df.loc[tree, other_tree] = results[tree][other_tree]
    return df


def clear_working_dir():
    os.system(f"rm {WORK_DIR}/*")


def run_sphere(fasta, custom_ref=None, perturbed=False, tree_index=None):
    # Construct commands to run sphere and produce binary newick
    if tree_index is None:
        raise ValueError(
            "Please provide an index to uniquely name the sphere output file"
        )

    if perturbed:
        nwk_path = f"{WORK_DIR}/sphere_binary_p_{tree_index}.nwk"
    else:
        nwk_path = f"{WORK_DIR}/sphere_binary_{tree_index}.nwk"

    if custom_ref:
        ref_path = custom_ref
    else:
        ref_path = f"{SPHERE_DIR}/sample_inputs/ref.fas"

    sphere_cmd = (
        f"java -jar {SPHERE_DIR}/sphere/sphere.jar "
        f"-i {fasta} "
        f"-r {ref_path} "
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
        f"python {UTILS_DIR}/binary.py " f"{WORK_DIR}/sphere_output.nwk " f"{nwk_path}"
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
    mc_cmd = f"python {MC_DIR}/monte_carlo_entropy.py {nwk} {fasta} {nwk_out_path} {num_iter}"
    subprocess.run(mc_cmd.split(" "))
    with open(nwk_out_path, "r") as nwk_file:
        nwk = nwk_file.read()
    return nwk, nwk_out_path


def run_raxml(fasta, index):
    rax_outpath = f"{WORK_DIR}/RAxML_bestTree.{index}"
    final_outpath = f"{WORK_DIR}/RAxML_bestTree_standardized.{index}"

    # Run RAxML
    raxml_cmd = f"{RAXML_DIR}/raxmlHPC -m GTRCAT -V -s {fasta} -n {index} -w {WORK_DIR} -p 12345 "
    subprocess.run(raxml_cmd.split(" "), stdout=subprocess.DEVNULL)

    # Standardize output tree to ensure it is binary
    standardize_raxml_cmd = (
        f"python {UTILS_DIR}/binary.py {rax_outpath} {final_outpath}"
    )
    subprocess.run(standardize_raxml_cmd.split(" "), stdout=subprocess.DEVNULL)

    with open(final_outpath, "r") as nwk_file:
        nwk = nwk_file.read()
    return nwk, final_outpath


def run_experiment(fasta, ref, method, tree_index, perturbed, mc_iter=10):
    print(f"[{method}] {perturbed=}\t{fasta}\t{ref}")

    # Combine fasta and reference
    fasta_with_ref = cat(fasta, ref)

    # Set output files
    if perturbed:
        mc_out = f"{method}_p_{tree_index}"
        raxml_out = f"tree_p_{tree_index}"
    else:
        mc_out = f"{method}_{tree_index}"
        raxml_out = f"tree_{tree_index}"

    # Run method
    if method.lower() == "sphere":
        nwk, nwk_path = run_sphere(fasta, ref, perturbed, tree_index)
    elif method.lower() == "raxml":
        nwk, nwk_path = run_raxml(fasta_with_ref, raxml_out)
    else:
        raise ValueError(
            f'{method} is an invalid choice for method. Please use one of ["sphere", "raxml"]'
        )

    # Minimize tree entropy with monte carlo
    mc_nwk, mc_nwk_path = run_monte_carlo(
        nwk_path, fasta_with_ref, mc_out, num_iter=mc_iter
    )

    print('Computing parsimony score...')
    orig_pscore = compute_parsimony_score(nwk_path, fasta_with_ref)
    mc_pscore = compute_parsimony_score(mc_nwk_path, fasta_with_ref)
    
    print('Computing distance along tree...')
    orig_tdist, orig_hdist, orig_dev = td.compute_avg_dist_along_tree(nwk_path, fasta_with_ref)
    mc_tdist, mc_hdist, mc_dev = td.compute_avg_dist_along_tree(mc_nwk_path, fasta_with_ref)

    return nwk_path, mc_nwk_path, orig_pscore, mc_pscore, orig_tdist, mc_tdist, orig_dev, mc_dev, orig_hdist, mc_hdist


def read_args():
    fasta_path = sys.argv[1]
    try:
        ref_path = sys.argv[2]
    except IndexError:
        ref_path = f"{SPHERE_DIR}/sample_inputs/ref.fas"
    return fasta_path, ref_path


def main():
    num_trees = 2
    mc_iter = 50
    p = 0.01

    fasta_path, ref_path = read_args()
    seqs = create_seqs_dict(fasta_path)

    tree_ids = {
        0: "sphere",
        1: "sphere_mc",
        2: "raxml",
        3: "raxml_mc",
        4: "sphere_pert",
        5: "sphere_pert_mc",
        6: "raxml_pert",
        7: "raxml_pert_mc",
    }

    running_results = {"dists": [], "pscores": [], "tdists": [], 'deviations': [], 'hdists': []}
    tree_index = 0

    while tree_index < num_trees:
        clear_working_dir()
        print(f"\nTree: {tree_index}")

        # Create perturbed sequences
        pert_seqs = perturb_sequences(seqs, p=p)
        pert_fasta_path = f"{WORK_DIR}/pert_seqs.fasta"
        write_seqs(pert_seqs, pert_fasta_path)

        # Run experiments
        computed_trees = []
        tdists = []
        deviations = []
        pscores = {}
        tree_id = 0
        for fasta in [fasta_path, pert_fasta_path]:
            perturbed = False if fasta == fasta_path else True
            for method in ["sphere", "raxml"]:
                nwk_path, mc_nwk_path, orig_pscore, mc_pscore, orig_tdist, mc_tdist, orig_dev, mc_dev, orig_hdist, mc_hdist = run_experiment(
                    fasta, ref_path, method, tree_index, perturbed, mc_iter
                )

                print(f" Parsimony: {orig_pscore}  --> {mc_pscore}")
                pscores[tree_ids[tree_id]] = orig_pscore
                tree_id += 1
                pscores[tree_ids[tree_id]] = mc_pscore
                tree_id += 1
                computed_trees.append(nwk_path)
                computed_trees.append(mc_nwk_path)
                tdists.append(orig_tdist)
                tdists.append(mc_tdist)
                deviations.append(orig_dev)
                deviations.append(mc_dev)

        distances = compute_pairwise_rf_distances(computed_trees)

        running_results["dists"].append(matrix_dict_to_df(distances))
        running_results["pscores"].append(pscores)
        running_results['tdists'].append(tdists)
        running_results['deviations'].append(deviations)

        tree_index += 1

    # Accumulate results
    avg_dists = sum(running_results["dists"]) / len(running_results["dists"])
    avg_dists.rename(index=tree_ids, columns=tree_ids, inplace=True)
    avg_dists.to_csv(f"{OUT_DIR}/distances.csv")
    print(avg_dists)

    pscore_df = pd.DataFrame(running_results["pscores"]).T
    pscore_df.to_csv(f"{OUT_DIR}/pscores.csv")
    print(pscore_df)

    tdist_df = pd.DataFrame(running_results['tdists']).T
    tdist_df.rename(index=tree_ids, inplace=True)
    tdist_df.to_csv(f"{OUT_DIR}/tdists.csv")
    print(tdist_df)

    dev_df = pd.DataFrame(running_results['deviations']).T
    dev_df.rename(index=tree_ids, inplace=True)
    dev_df.to_csv(f"{OUT_DIR}/deviations.csv")
    print(dev_df)

if __name__ == "__main__":
    main()
    # with cProfile.Profile() as pr:
    #     main()
    # stream = io.StringIO()
    # stats = pstats.Stats(pr, stream=stream)
    # stats.sort_stats("cumtime")
    # stats.print_stats()

    # with open('stats.txt', 'w') as f:
    #     f.write(stream.getvalue())
    
