import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from loguru import logger
import warnings
import os

warnings.filterwarnings("ignore")

# Configure logger
logger.add("logs/plots.log", format="{time} {message}", level="DEBUG", rotation="5 MB")

def plot_usage_stats(usage_data, plot: bool = False, filename: str = None) -> None:
    """
    Plot usage statistics of V genes based on clonal lineages.

    Parameters:
        usage_data (dict): Dictionary where keys are V genes and values are their usage counts.
        plot (bool): Whether to display the plot interactively. Defaults to False.
        filename (str): If provided, saves the plot to a file with this name in the 'plots' directory.

    Returns:
        None

    Logs:
        - Logs the process of plotting and saving the file.
    """
    logger.info("Starting to plot usage statistics.")

    # Validate input data
    if not usage_data or not isinstance(usage_data, dict):
        logger.error("Invalid usage data provided. Expected a non-empty dictionary.")
        raise ValueError("usage_data must be a non-empty dictionary.")

    # Sort data by values for better visualization
    sorted_data = dict(sorted(usage_data.items(), key=lambda item: item[1], reverse=True))

    try:
        # Set up the plot
        plt.figure(figsize=(14, 8))
        sns.barplot(x=list(sorted_data.keys()), y=list(sorted_data.values()), palette="viridis")

        # Add titles and labels
        plt.title('V Gene Usage Statistics', fontsize=16, fontweight='bold', color='darkblue')
        plt.xlabel('V Genes', fontsize=14, fontweight='bold')
        plt.ylabel('Clonal Lineage Count', fontsize=14, fontweight='bold')

        # Improve X-axis readability
        plt.xticks(rotation=90, fontsize=10, color='darkgreen')
        plt.yticks(fontsize=10, color='darkgreen')

        # Add grid for better clarity
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        # Optimize layout
        plt.tight_layout()

        # Show the plot if requested
        if plot:
            logger.info("Displaying the plot interactively.")
            plt.show()

        # Save the plot if a filename is provided
        if filename:
            output_dir = "plots"
            os.makedirs(output_dir, exist_ok=True)
            file_path = os.path.join(output_dir, f"{filename}.png")
            plt.savefig(file_path)
            logger.info(f"Plot saved to file: {file_path}")

        # Clear the figure after saving or showing
        plt.close()

    except Exception as e:
        logger.error(f"An error occurred while plotting: {e}")
        raise

    logger.success("Plotting completed successfully.")



def plot_stats(stats):
    """
    Convert the stats dictionary into a DataFrame and create a bar plot.

    Parameters:
        stats (dict): A dictionary containing clonal lineage statistics.

    Returns:
        None: Displays a plot and saves it as an image.
    """
    # Convert the dictionary to a DataFrame
    df_stats = pd.DataFrame(list(stats.items()), columns=['Metric', 'Value'])

    # Set up the plot
    plt.figure(figsize=(12, 6))
    sns.barplot(data=df_stats, x='Value', y='Metric', palette='viridis', orient='h')

    # Add titles and labels
    plt.title('Clonal Lineage Statistics', fontsize=16, fontweight='bold', color='darkblue')
    plt.xlabel('Value', fontsize=14, fontweight='bold')
    plt.ylabel('Metric', fontsize=14, fontweight='bold')

    # Add values to the bars for better readability
    for index, value in enumerate(df_stats['Value']):
        plt.text(value + 0.1, index, str(value), va='center', fontsize=12, color='black')

    # Optimize layout and display
    plt.tight_layout()

    # Save and show the plot
    output_dir = "plots"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "clonal_lineage_stats.png"))
    plt.show()