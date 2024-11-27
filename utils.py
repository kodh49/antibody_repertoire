import seaborn as sns
import matplotlib.pyplot as plt

def plot_usage_stats(usage_data):    

    # Sort data by values for better visualization
    sorted_data = dict(sorted(usage_data.items(), key=lambda item: item[1], reverse=True))

    # Enhanced plotting with seaborn and matplotlib
    plt.figure(figsize=(14, 8))
    sns.barplot(x=list(sorted_data.keys()), y=list(sorted_data.values()), palette="viridis")

    # Add titles and labels with improved styling
    plt.title('V Gene Usage Statistics', fontsize=16, fontweight='bold', color='darkblue')
    plt.xlabel('V Genes', fontsize=14, fontweight='bold')
    plt.ylabel('Clonal Lineage Count', fontsize=14, fontweight='bold')

    # Improve X-axis readability
    plt.xticks(rotation=90, fontsize=10, color='darkgreen')
    plt.yticks(fontsize=10, color='darkgreen')

    # Add grid for better visual clarity
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Tight layout for better spacing
    plt.tight_layout()
    plt.show()
