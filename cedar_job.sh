#!/bin/bash
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=00:38:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

srun --ntasks $SLURM_NNODES --tasks-per-node=1 bash << EOF

module load python/3.12

virtualenv --no-download $SLURM_TMPDIR/test_env
source $SLURM_TMPDIR/env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt

EOF

# Define input and output paths
QUERY="data/PRJNA324093_Dnr4_10k.fasta"
OUTPUT="data/igblast_results.tsv"
DB_DIR="database"
USAGE_PLOT="usage_plot"
WEBLOGO_QUERY="weblogo_query"
AUXILIARY="optional_file/human_gl.aux"
CORES=48

cd ncbi-igblast-1.22.0

# run ncbi-igblash-1.22.0
bin/igblastn -germline_db_V $DB_DIR/my_seq_V -germline_db_D $DB_DIR/my_seq_D -germline_db_J $DB_DIR/my_seq_J -organism human -domain_system imgt -query ../$QUERY -outfmt 19 -out ../$OUTPUT -auxiliary_data $AUXILIARY

cd ..

# run analysis script
python3 driver.py --cache n --igblast $OUTPUT --usage_plot $USAGE_PLOT --weblogo_query $WEBLOGO_QUERY
