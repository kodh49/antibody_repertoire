import networkx as nx
from tqdm import trange
from itertools import combinations
import pickle
import os
from loguru import logger

# Configure loguru to log to both stdout and a file
logger.add("logs/hamming.log", format="{time} {message}", rotation="10 MB", retention="10 days", backtrace=True, diagnose=True, level="INFO")

def hamming_distance(str1, str2):
    """
    Compute the Hamming distance between two strings of equal length.

    Parameters:
        str1 (str): The first string.
        str2 (str): The second string.

    Returns:
        int: The Hamming distance (number of differing characters).
    """
    return sum(symb1 != symb2 for symb1, symb2 in zip(str1, str2))

def init(igblast_result, use_cache: bool = False):
    """
    Generate a Hamming graph from igBlast result data.

    In this graph:
    - Each vertex represents a CDR3 sequence.
    - Two vertices are connected if:
      1. They have the same length.
      2. Their Hamming distance (nucleotide sequences) is below 10% of the CDR3 length.
      3. They represent VDJ sequences with identical V and J genes (ignoring alleles).

    Parameters:
        igblast_result (dict): A dictionary of dictionaries where each key is a sequence ID 
                               and the value is a dictionary containing:
                               - 'cdr3': CDR3 sequence (str)
                               - 'cdr_length': Length of the CDR3 sequence (int)
                               - 'v_call': Set of V genes (set)
                               - 'j_call': Set of J genes (set)
        use_cache (bool): If True, attempts to load the graph from a cached file. If False,
                          computes the graph from scratch and saves it for future use.

    Returns:
        nx.Graph: A NetworkX graph object representing the computed Hamming graph.
    """
    cache_path = "cache/HammingGraph.pkl"

    # Check if the graph exists in cache and load if requested
    if use_cache and os.path.isfile(cache_path):
        logger.info("Loading Hamming graph from cache...")
        try:
            with open(cache_path, 'rb') as f:
                HammingGraph = pickle.load(f)
            logger.info("Hamming graph successfully loaded from cache.")
            return HammingGraph
        except Exception as e:
            logger.error(f"Failed to load cache file: {e}")
            raise

    logger.info("Initializing Hamming graph computation...")

    # Initialize an empty graph
    HammingGraph = nx.Graph()
    igblast_result = igblast_result.to_dict('index')

    # Generate all unique pairs of VDJ sequences
    vdj_seq_pairs = list(combinations(igblast_result.keys(), 2))
    logger.info(f"Total sequence pairs to process: {len(vdj_seq_pairs)}")

    # Loop over each pair of sequences
    for i in trange(len(vdj_seq_pairs), desc="Processing sequence pairs"):
        vdj_seq_pair = vdj_seq_pairs[i]

        # Retrieve sequence data for the current pair
        seq_a = igblast_result[vdj_seq_pair[0]]
        seq_b = igblast_result[vdj_seq_pair[1]]

        # Add nodes for the CDR3 sequences to the graph
        HammingGraph.add_node(seq_a['cdr3'])
        HammingGraph.add_node(seq_b['cdr3'])

        # Check conditions for adding an edge between these nodes
        try:
            # (1) Ensure both sequences have the same CDR3 length
            assert seq_a['cdr_length'] == seq_b['cdr_length'], "CDR3 lengths differ"

            # (2) Ensure Hamming distance is below 10% of the CDR3 length
            hamming_dist = hamming_distance(seq_a['cdr3'], seq_b['cdr3'])
            assert hamming_dist <= 0.1 * seq_a['cdr_length'], f"Hamming distance {hamming_dist} exceeds threshold"

            # (3) Ensure sequences have identical V and J genes (ignoring alleles)
            assert len(seq_a['v_call'].intersection(seq_b['v_call'])) > 0, "No common V genes"
            assert len(seq_a['j_call'].intersection(seq_b['j_call'])) > 0, "No common J genes"

            # Add an edge between the nodes if all conditions are met
            HammingGraph.add_edge(seq_a['cdr3'], seq_b['cdr3'])

        except AssertionError as e:
            continue

    # Save the computed graph to cache
    os.makedirs("cache", exist_ok=True)
    try:
        with open(cache_path, 'wb') as f:
            pickle.dump(HammingGraph, f)
        logger.info("Hamming graph successfully saved to cache.")
    except Exception as e:
        logger.error(f"Failed to save Hamming graph to cache: {e}")
        raise

    logger.success("Hamming graph computation completed.")
    return HammingGraph