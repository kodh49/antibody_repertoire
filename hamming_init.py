import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
import igblast_init as ig

# generate the Hamming graph from the data
"""
Compute a Hamming graph on the computed CDR3s where each CDR3 represents a vertex and two CDR3s are connected with an edge if (i) they have same length, (ii) the Hamming distance between their nucleotide sequences is below 10%, (iii) they represent VDJ sequences with identical V and J genes. Ignore gene alleles while checking the third condition: e.g., alleles IGHV1-69*01 and IGHV1-69*02 represent the same gene.
"""

igblast_result = ig.igblast_preprocess()
HammingGraph = nx.Graph() # initiate the hamming graph

def init(igblast_result):
    """
    Create a Hamming graph on the computed CDR3s from the sequencing data
    """
    vdj_seq_pairs = list(combinations(igblast_result, 2)) 
    for vdj_seq_pair in vdj_seq_pairs:
        # select two VDJ sequences
        seq_a = vdj_seq_pair[0]
        seq_b = vdj_seq_pair[1]
        # add CDR3s of selected sequences to the Hamming Graph
        HammingGraph.add_node(seq_a['cdr3'])
        HammingGraph.add_node(seq_b['cdr3'])
        # check the condition for external edges
        try:
            assert len(seq_a['cdr3']) == len(seq_b['cdr3']) # Two CDR3s have the same length
            assert hamming_distance(seq_a['cdr3'], seq_b['cdr3']) < 0.1 * len(seq_a['cdr3']) # hamming distance between nucleotide identity is below 10%
            assert seq_a['v'] == seq_b['v'] and seq_a['j'] == seq_b['j'] # represent identical V and J genes
            # add an edge between two CDR3s
            HammingGraph.add_edge(seq_a['cdr3'], seq_b['cdr3'])
        except:
            # skip the addition of an edge between two selected sequences
            continue


def hamming_distance(str1, str2):
    """
    Computes the Hamming Distance between two strings of equal length
    """
    return sum(symb1 != symb2 for symb1, symb2 in zip(str1, str2))