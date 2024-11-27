import networkx as nx
import pandas as pd
import pickle
import math
import os
from loguru import logger

# Configure logger
logger.add("logs/lineage_analysis.log", format="{time} {level} {message}", level="DEBUG", rotation="10 MB")

def get_clonal_lineages(HammingGraph):
    """
    Identify and report clonal lineages from a Hamming graph.

    Each connected component in the Hamming graph corresponds to a clonal lineage.
    The lineages are sorted by size in descending order.

    Parameters:
        HammingGraph (networkx.Graph): The Hamming graph representing relationships between CDR3 sequences.

    Returns:
        list: A list of connected components (clonal lineages), sorted by size in descending order.

    Logs:
        - Number of clonal lineages identified.
        - Errors encountered during the process.
    """
    logger.info("Extracting clonal lineages from the Hamming graph.")
    try:
        lineages = sorted(nx.connected_components(HammingGraph), key=len, reverse=True)
        logger.info(f"Identified {len(lineages)} clonal lineages.")
        return lineages
    except Exception as e:
        logger.error(f"Failed to extract clonal lineages: {e}")
        raise


def get_lineage_from_cdr3(row, ClonalLineage):
    """
    Find and return the clonal lineage to which a given CDR3 sequence belongs.

    Parameters:
        row (pd.Series): A row from the DataFrame containing the 'cdr3' field.
        ClonalLineage (list): A list of clonal lineages (sets of CDR3 sequences).

    Returns:
        set or None: The clonal lineage (set of CDR3 sequences) the sequence belongs to, or None if not found.
    """
    for comp in ClonalLineage:
        if row['cdr3'] in comp:
            return comp
    return None  # Return None if the CDR3 does not belong to any lineage.


def get_lineages_stats(igblast_result, ClonalLineage):
    """
    Compute and report statistics related to clonal lineages.

    The function calculates several statistics, including the number of lineages,
    details about the largest lineage, and lineages with significant representation.

    Parameters:
        igblast_result (pd.DataFrame): Processed IgBLAST result containing CDR3 sequences.
        ClonalLineage (list): A list of clonal lineages (sets of CDR3 sequences).

    Returns:
        dict: A dictionary containing statistics about clonal lineages.

    Logs:
        - Number of lineages and statistics about the largest lineage.
        - Errors during clonal lineage assignment or grouping.
    """
    logger.info("Starting computation of clonal lineage statistics.")

    # Total number of clonal lineages
    num_of_lineages = len(ClonalLineage)

    # Largest clonal lineage
    lcl = ClonalLineage[0] if ClonalLineage else set()
    num_of_cdr3_lcl = len(lcl)  # Total number of unique CDR3 regions in the largest clonal lineage

    # Number of clonal lineages with at least 10 unique CDR3s
    num_of_lineages_cdr3_cutoff = len([x for x in ClonalLineage if len(x) >= 10])

    # Assign clonal lineages to each sequence
    logger.info("Assigning clonal lineages to sequences in the IgBLAST result.")
    try:
        igblast_result['lineages'] = igblast_result.apply(get_lineage_from_cdr3, axis=1, args=(ClonalLineage,))
        igblast_result['lineages'] = igblast_result['lineages'].apply(lambda x: tuple(x) if x is not None else None)
    except Exception as e:
        logger.error(f"Failed to assign clonal lineages: {e}")
        raise

    # Mark whether each sequence belongs to the largest clonal lineage
    igblast_result['in_lcl'] = igblast_result['lineages'].apply(lambda x: x == tuple(lcl))
    num_of_seqs_lcl = len(igblast_result.loc[igblast_result['in_lcl']].to_numpy())

    # Group by clonal lineage and filter those with at least 10 sequences
    logger.info("Grouping sequences by clonal lineage and filtering lineages with at least 10 sequences.")
    try:
        lineage_counts = igblast_result.groupby('lineages').size().reset_index(name='count')
        filtered_lineages = lineage_counts[lineage_counts['count'] >= 10].to_numpy()
        num_of_lineages_seq_cutoff = len(filtered_lineages)
    except Exception as e:
        logger.error(f"Failed to compute lineage groupings: {e}")
        raise

    # Prepare the result dictionary
    stats = {
        'number of clonal lineages': num_of_lineages,
        'number of unique CDR3s in the largest clonal lineage': num_of_cdr3_lcl,
        'number of sequences in the largest clonal lineage': num_of_seqs_lcl,
        'number of lineages represented by at least 10 unique CDR3s': num_of_lineages_cdr3_cutoff,
        'number of lineages represented by at least 10 sequences': num_of_lineages_seq_cutoff,
    }

    logger.info("Completed computation of clonal lineage statistics.")
    return stats


def get_usage_stats(igblast_result, HammingGraph, use_cache=False, cache_path="cache/UsageCounts.pkl"):
    """
    Compute usage statistics of V genes based on clonal lineages.

    Parameters:
        igblast_result (pd.DataFrame): Processed IgBLAST result containing V gene calls and CDR3 sequences.
        HammingGraph (networkx.Graph): The Hamming graph used to identify clonal lineages.
        use_cache (bool, optional): If True, load results from cache if available. Defaults to False.
        cache_path (str, optional): Path to the cache file. Defaults to "cache/UsageCounts.pkl".

    Returns:
        dict: A dictionary where keys are V genes and values are their usage counts.

    Logs:
        - Usage statistics for each V gene.
        - Errors during computation or cache operations.
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
        node_to_component = {node: i for i, component in enumerate(nx.connected_components(HammingGraph)) for node in component}
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

    Logs:
        - Details about the amino acid sequence extraction and file saving process.
        - Errors encountered during the process.
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
