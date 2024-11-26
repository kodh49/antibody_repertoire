import networkx as nx
import pickle

# load graph object from file
HammingGraph = pickle.load(open('HammingGraph.pkl', 'rb'))


def get_clonal_lineage(HammingGraph):
    """
    Report connected components of the computed Hamming graph as clonal lineages.
    """
    ConnectedComponents = list(filter(lambda x: len(x)>1,
                                  [c for c in sorted(nx.connected_components(HammingGraph), key=len, reverse=True)]
                                  ))
    return ConnectedComponents

def get_usage_stats(igblast_result, ConnectedComponents):
    """
    Create a usage plot of the computed V genes (x axis = V genes, y axis = the number of clonal lineages formed by each of V genes).
    For each V gene in the igBlast analysis, compute the number of clonal lineages, i.e. number of connected components that this particular V gene is part of.
    """
    v_cdr3 = igblast_result.loc[['v_call', 'cdr3']]
    v_genes = set()
    for v_gene in v_cdr3['v_call']:
        v_genes.update(v_gene)
    
    # Create a dictionary with keys from all V genes and initial value 0
    usage_counts = {v_gene: 0 for v_gene in v_genes}

    for sequence in v_cdr3:
        cdr3 = sequence['cdr3']
        component = nx.node_connected_component(ConnectedComponents, cdr3)
        if len(component) != 0: # selected cdr3 is in the clonal lineage
            for v_gene in sequence['v_call']:
                usage_counts[v_gene] += 1
    
    return usage_counts


"""
Create a [web logo](https://weblogo.berkeley.edu/logo.cgi)[Links to an external site.](https://weblogo.berkeley.edu/logo.cgi) plot of amino acid sequences of CDR3s from the largest clonal lineage.
"""