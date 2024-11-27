import networkx as nx
import pandas as pd
import pickle
import math
import os
from loguru import logger

# Configure logger
logger.add("logs/lineage_analysis.log", format="{time} {message}", level="DEBUG", rotation="10 MB")

def get_clonal_lineages(HammingGraph):
    """
    Identify and report clonal lineages from a Hamming graph.

    Each connected component in the Hamming graph corresponds to a clonal lineage.

    Parameters:
        HammingGraph (networkx.Graph): The Hamming graph representing relationships between CDR3 sequences.

    Returns:
        list: A list of connected components, sorted by size in descending order.
    """
    logger.info("Extracting clonal lineages from the Hamming graph.")
    try:
        lineages = sorted(nx.connected_components(HammingGraph), key=len, reverse=True)
        logger.info(f"Identified {len(lineages)} clonal lineages.")
        return lineages
    except Exception as e:
        logger.error(f"Failed to extract clonal lineages: {e}")
        raise


def get_usage_stats(igblast_result, HammingGraph, use_cache=False, cache_path="cache/UsageCounts.pkl"):
    """
    Compute usage statistics of V genes based on clonal lineages.

    The usage statistics are the number of clonal lineages that each V gene is part of. Results can be cached for future use.

    Parameters:
        igblast_result (pd.DataFrame): The processed IgBLAST result containing V gene calls and CDR3 sequences.
        HammingGraph (networkx.Graph): The Hamming graph used to identify clonal lineages.
        use_cache (bool, optional): If True, load results from the cache if available. Defaults to False.
        cache_path (str, optional): Path to the cache file. Defaults to "cache/UsageCounts.pkl".

    Returns:
        dict: A dictionary where keys are V genes and values are their usage counts.
    """
    logger.info("Computing usage statistics for V genes.")
    if use_cache and os.path.isfile(cache_path):
        logger.info("Loading usage statistics from cache.")
        try:
            with open(cache_path, 'rb') as f:
                usage_counts = pickle.load(f)
            logger.info("Loaded usage statistics from cache.")
            return usage_counts
        except Exception as e:
            logger.warning(f"Failed to load cache file: {e}")

    # Map nodes to their clonal lineage
    logger.info("Mapping nodes to their respective clonal lineages.")
    try:
        node_to_component = {}
        for i, component in enumerate(nx.connected_components(HammingGraph)):
            for node in component:
                node_to_component[node] = i
    except Exception as e:
        logger.error(f"Failed to map nodes to clonal lineages: {e}")
        raise

    # Compute usage counts for V genes
    v_genes = set(gene for genes in igblast_result['v_call'] for gene in genes)
    usage_counts = {v_gene: 0 for v_gene in v_genes}

    logger.info("Counting V gene usage.")
    try:
        for v_call, cdr3 in zip(igblast_result['v_call'], igblast_result['cdr3']):
            if cdr3 in node_to_component:
                for v_gene in v_call:
                    usage_counts[v_gene] += 1
    except Exception as e:
        logger.error(f"Error during V gene usage counting: {e}")
        raise

    # Cache the results
    logger.info("Caching usage statistics.")
    try:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, 'wb') as f:
            pickle.dump(usage_counts, f)
        logger.info("Cached usage statistics.")
    except Exception as e:
        logger.warning(f"Failed to cache usage statistics: {e}")

    return usage_counts


def get_aaseq_from_lcl(igblast_result, ClonalLineage, filename=None):
    """
    Extract amino acid sequences of CDR3s from the largest clonal lineage and optionally save them to a file.

    Parameters:
        igblast_result (pd.DataFrame): The processed IgBLAST result containing CDR3 amino acid sequences.
        ClonalLineage (list): A list of connected components from the Hamming graph.
        filename (str, optional): If provided, the sequences will be saved to this file. Defaults to None.

    Returns:
        str: A string containing amino acid sequences, one per line.
    """
    logger.info("Extracting amino acid sequences from the largest clonal lineage.")
    try:
        # Extract the largest clonal lineage
        largest_clonal_lineage = list(ClonalLineage[0])
        cdr3_df = pd.DataFrame(largest_clonal_lineage, columns=['cdr3'])

        # Merge with igblast_result to get amino acid sequences
        aaseq_df = pd.merge(igblast_result, cdr3_df, how='inner', on='cdr3')
        aa_seqs = aaseq_df['cdr3_aa'].str.replace(' ', '*').tolist()

        # Compute average length and pad sequences
        avg_length = math.ceil(sum(map(len, aa_seqs)) / len(aa_seqs))
        aa_seqs = [seq.ljust(avg_length, '*') for seq in aa_seqs]

        # Generate the query string
        query_string = '\n'.join(aa_seqs)

        # Save to file if a filename is provided
        if filename:
            output_dir = "data"
            os.makedirs(output_dir, exist_ok=True)
            with open(os.path.join(output_dir, f"{filename}.txt"), 'w') as f:
                f.write(query_string)
            logger.info(f"Amino acid sequences saved to {filename}.txt in the data directory.")

        return query_string
    except Exception as e:
        logger.error(f"Failed to extract amino acid sequences: {e}")
        raise
