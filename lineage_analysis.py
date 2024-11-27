import networkx as nx
import pickle
import igblast_init as ig
import pandas as pd
import utils
import math

def get_clonal_lineage(HammingGraph):
    """
    Report connected components of the computed Hamming graph as clonal lineages.
    """
    ConnectedComponents = list(filter(lambda x: len(x)>1,
                                  [c for c in sorted(nx.connected_components(HammingGraph), key=len, reverse=True)]
                                  ))
    return ConnectedComponents

def get_usage_stats(igblast_result, HammingGraph, verbose):
    """
    Create a usage plot of the computed V genes (x axis = V genes, y axis = the number of clonal lineages formed by each of V genes).
    For each V gene in the igBlast analysis, compute the number of clonal lineages, i.e. number of connected components that this particular V gene is part of.
    """
    v_cdr3 = igblast_result.loc[:, ['v_call', 'cdr3']]
    v_genes = set()
    for v_gene in v_cdr3['v_call']:
        v_genes.update(v_gene)
    
    # Create a dictionary with keys from all V genes and initial value 0
    usage_counts = {v_gene: 0 for v_gene in v_genes}

    for _, sequence in v_cdr3.iterrows():
        cdr3 = sequence['cdr3']
        component = nx.node_connected_component(HammingGraph, cdr3)
        if len(component) != 0: # selected cdr3 is in the clonal lineage
            for v_gene in sequence['v_call']:
                usage_counts[v_gene] += 1

    # cache the result for future use
    pickle.dump(usage_counts, open('cache/usage_counts.pkl', 'wb'))

    return usage_counts


def get_aaseq_from_lcl(igblast_result, clonal_lineage, filename):
    """
    Create a list of amino acid sequences of CDR3s from the largest clonal lineage and save it as a text file with the given filename
    """
    lcl = list(clonal_lineage[0]) # largest clonal lineage
    cdr3_df = pd.DataFrame(lcl, columns=['cdr3']) # table of amino-acid sequences
    aaseq_df = pd.merge(igblast_result, cdr3_df, how='inner', on='cdr3')
    aa_seqs = list(map(lambda x: x.replace(' ', '*'), aaseq_df['cdr3_aa'].tolist()))
    amino_acid_length = math.ceil(sum(map(len, aa_seqs))/len(aa_seqs))
    for i in range(len(aa_seqs)):
        amino_acid_seq = aa_seqs[i]
        if len(amino_acid_seq) < amino_acid_length:
            aa_seqs[i] = amino_acid_seq + "*"
    query_string = '\n'.join(aa_seqs)

    if filename is not None:
        txt_file = open(f'data/{filename}.txt', 'w')
        txt_file.write(query_string)
        txt_file.close()

    return query_string

# load graph object from file
HammingGraph = pickle.load(open('HammingGraph.pkl', 'rb'))
ConnectedComponents = get_clonal_lineage(HammingGraph=HammingGraph)
igblast_result = ig.igblast_preprocess(to_dict=False)
usage_stats = get_usage_stats(igblast_result=igblast_result, HammingGraph=HammingGraph)
query_string = get_aaseq_from_lcl(igblast_result=igblast_result, clonal_lineage=ConnectedComponents)
print(query_string)

utils.plot_usage_stats(usage_data=usage_stats)