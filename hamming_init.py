import networkx as nx
from tqdm import trange
from itertools import combinations
import igblast_init as ig
import pickle

# generate the Hamming graph from the data
def init(igblast_result):
    """
    Compute a Hamming graph on the computed CDR3s where each CDR3 represents a vertex and two CDR3s are connected with an edge if (i) they have same length, (ii) the Hamming distance between their nucleotide sequences is below 10%, (iii) they represent VDJ sequences with identical V and J genes. Ignore gene alleles while checking the third condition: e.g., alleles IGHV1-69*01 and IGHV1-69*02 represent the same gene.
    """
    HammingGraph = nx.Graph() # initiate the hamming graph
    vdj_seq_pairs = list(combinations(igblast_result.keys(), 2)) 
    for i in trange(len(vdj_seq_pairs)):
        vdj_seq_pair = vdj_seq_pairs[i]
        # select two VDJ sequence object from igBlast result set
        seq_a = igblast_result[vdj_seq_pair[0]]
        seq_b = igblast_result[vdj_seq_pair[1]]
        # add CDR3s of selected sequences to the Hamming Graph
        HammingGraph.add_node(seq_a['cdr3'])
        HammingGraph.add_node(seq_b['cdr3'])
        # check the condition for external edges
        try:
            assert seq_a['cdr_length'] == seq_b['cdr_length'] # Two CDR3s have the same length
            assert hamming_distance(seq_a['cdr3'], seq_b['cdr3']) <= 0.1 * seq_a['cdr_length'] # hamming distance between nucleotide identity is below 10%
            assert len(seq_a['v_call'].intersection(seq_b['v_call'])) != 0 # represent identical V genes (non-empty intersection)
            assert len(seq_a['j_call'].intersection(seq_b['j_call'])) != 0 # represent identical J genes (non-empty intersection)
            # add an edge between two CDR3s
            HammingGraph.add_edge(seq_a['cdr3'], seq_b['cdr3'])
        except:
            continue # skip the addition of an edge between two selected sequences
    return HammingGraph


def hamming_distance(str1, str2):
    """
    Computes the Hamming Distance between two strings of equal length
    """
    return sum(symb1 != symb2 for symb1, symb2 in zip(str1, str2))

igblast_result = ig.igblast_preprocess()
HammingGraph = init(igblast_result=igblast_result)

# save graph object to file
pickle.dump(HammingGraph, open('HammingGraph.pkl', 'wb'))

# load graph object from file
HammingGraph = pickle.load(open('HammingGraph.pkl', 'rb'))

print(HammingGraph)