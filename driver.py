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
        required=False,
        default="y",
    )
    parser.add_argument(
        "--igblast",
        type=str,
        help="Path to the NCBI igBLAST analysis result tsv file.",
        required=True,
    )
    parser.add_argument(
        "--usage_plot",
        type=str,
        help="Filename of the V gene usage plot.",
        required=False,
    )
    parser.add_argument(
        "--weblogo_query",
        type=str,
        help="Filename of the query for generating weblogo plot.",
        required=False,
    )


def main(args):
    usage_plot = args.usage_plot # filename of the usage plot
    weblogo_query = args.weblogo_query # filename of the query for generating weblogo
    igblast_result_path = args.igblast  # path to NCBI igBLAST analysis result
    if args.cache == 'y':
        use_cache = True
    else:
        use_cache = False

    # preprocess the igblast result
    igblast_result = igblast.igblast_preprocess(filepath=igblast_result_path)

    # create Hamming graph from the igblast result
    HammingGraph = hamming.init(igblast_result=igblast_result, use_cache=use_cache)

    # report clonal lineage
    ClonalLineage = lineage_analysis.get_clonal_lineages(HammingGraph=HammingGraph)

    LineageStats = lineage_analysis.get_lineages_stats(igblast_result=igblast_result, ClonalLineage=ClonalLineage)

    # Plot the graph of lineage statistics
    utils.plot_stats(LineageStats)

    # report usage stats
    UsageStats = lineage_analysis.get_usage_stats(igblast_result=igblast_result, HammingGraph=HammingGraph, use_cache=use_cache)
    
    # Plot the graph of usage statistics
    utils.plot_usage_stats(usage_data=UsageStats, plot=True, filename=usage_plot)
    
    # Generate a list of Amino-acid sequences of CDR3s from the largest clonal lineage
    AminoAcidSeqs = lineage_analysis.get_aaseq_from_lcl(igblast_result=igblast_result, ClonalLineage=ClonalLineage, filename=weblogo_query)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script performs clonal lineage analysis based on the NCBI-igBlast processed immunoglobin data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)