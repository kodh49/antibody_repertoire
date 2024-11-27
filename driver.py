import hamming
import igblast
import lineage_analysis
import utils
import argparse

"""
This is the driver code of the immunogemonics analysis
"""

def add_arguments(parser):
    parser.add_argument(
        "--cache",
        type=str,
        choices=["y", "n"],
        help="Select whether to use cached data for the computation.",
        required=True,
    )
    parser.add_argument(
        "--usage_plot",
        type=str,
        help="Filename of the V gene usage plot.",
        required=True,
    )
    parser.add_argument(
        "--weblogo_query",
        type=str,
        help="Filename of the query for generating weblogo plot.",
        required=True,
    )
    parser.add_argument(
        "--igblast",
        type=str,
        help="Path to the NCBI igBLAST analysis result tsv file.",
        required=True,
    )
    parser.add_argument(
        "--num_threads",
        type=int,
        help="Number of threads to be used for the computation.",
        required=False,
        default=8,
    )


def main(args):
    usage_plot = args.usage_plot # filename of the usage plot
    weblogo_query = args.weblogo_query # filename of the query for generating weblogo
    num_threads = args.num_threads # number of threads to use for the computation
    igblast_result_path = args.igblast  # path to NCBI igBLAST analysis result
    if args.cache == 'y':
        use_cache = True
    else:
        use_cache = False

    # preprocess the igblast result
    igblast_result_dataframe = igblast.igblast_preprocess(to_dict=False, filepath=igblast_result_path)
    igblast_result_dict = igblast.igblast_preprocess(to_dict=True, filepath=igblast_result_path)

    # create Hamming graph from the igblast result
    HammingGraph = hamming.init(igblast_result=igblast_result_dict, use_cache=use_cache)

    # report clonal lineage
    ClonalLineage = lineage_analysis.get_clonal_lineages(HammingGraph=HammingGraph)

    # report usage stats
    UsageStats = lineage_analysis.get_usage_stats(igblast_result=igblast_result_dataframe,
                                                  HammingGraph=HammingGraph,
                                                  use_cache=use_cache)
    
    # Plot the graph of usage statistics
    utils.plot_usage_stats(usage_data=UsageStats, plot=True, filename=usage_plot)
    
    # Generate a list of Amino-acid sequences of CDR3s from the largest clonal lineage
    AminoAcidSeqs = lineage_analysis.get_aaseq_from_lcl(igblast_result=igblast_result_dataframe,
                                                        ClonalLineage=ClonalLineage,
                                                        filename=weblogo_query)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script computes classical optimal n-couplings from multiple marginal distributions and a cost tensor.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)