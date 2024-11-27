import hamming_init
import igblast_init
import lineage_analysis
import utils

"""
This is the driver code of the immunogemonics analysis
"""

if __name__ == '__main__':

    # preprocess the igblast result
    igblast_result_dataframe = igblast_init.igblast_preprocess(to_dict=False)
    igblast_result_dict = igblast_init.igblast_preprocess(to_dict=True)

    # create Hamming graph from the igblast result
    HammingGraph = hamming_init.init(igblast_result=igblast_result_dict)

    # report clonal lineage
    ClonalLineage = lineage_analysis.get_clonal_lineage(HammingGraph=HammingGraph)

    # report usage stats
    UsageStats = lineage_analysis.get_usage_stats(igblast_result=igblast_result_dataframe, HammingGraph=HammingGraph)
    
    # Generate a list of Amino-acid sequences of CDR3s from the largest clonal lineage
    AminoAcidSeqs = lineage_analysis.get_aaseq_from_lcl(igblast_result=igblast_result_dataframe, clonal_lineage=ClonalLineage)
    
    # Plot the graph
    utils.plot_usage_stats(usage_data=UsageStats)