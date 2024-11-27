import pandas as pd
import dask.dataframe as dd
from loguru import logger
import os

# Configure logger
logger.add("logs/preprocess.log", format="{time} {message}", level="DEBUG", rotation="10 MB")

def igblast_preprocess(filepath='data/igblast_results.tsv'):
    """
    Preprocess the NCBI IgBLAST results to extract relevant columns, normalize gene calls,
    compute CDR3 lengths, and optionally convert to dictionary format.

    Parameters:
        to_dict (bool): Whether to return the processed DataFrame as a dictionary. Defaults to False.
        filepath (str): Path to the IgBLAST results file. Defaults to 'data/igblast_results.tsv'.

    Returns:
        pd.DataFrame or dict: Processed DataFrame or its dictionary representation.

    Logs:
        - Logs the preprocessing steps, including file loading, transformations, and errors.
    """
    logger.info("Starting IgBLAST preprocessing.")
    
    # Define columns to retain from the input file
    assorted_columns = ['v_call', 'd_call', 'j_call', 'cdr3', 'cdr3_aa', 'cdr3_start', 'cdr3_end']

    # Validate file existence
    if not os.path.isfile(filepath):
        logger.error(f"File not found: {filepath}")
        raise FileNotFoundError(f"Input file '{filepath}' does not exist.")

    try:
        # Load the file and select only the relevant columns
        logger.info(f"Loading file: {filepath}")
        df = dd.read_csv(filepath, sep='\t', usecols=assorted_columns)
    except Exception as e:
        logger.error(f"Failed to load file '{filepath}': {e}")
        raise

    try:
        # Normalize gene calls by removing alleles
        logger.info("Normalizing gene calls by removing alleles.")
        for col in ['v_call', 'd_call', 'j_call']:
            df[col] = df[col].str.replace(r"\*\d+", '', regex=True)

        # Compute CDR3 region lengths
        logger.info("Computing CDR3 region lengths.")
        df['cdr_length'] = df['cdr3_end'] - df['cdr3_start'] + 1

        # Drop unnecessary columns
        df = df.drop(columns=['cdr3_start', 'cdr3_end'])

        # Convert comma-separated genes into sets
        logger.info("Converting comma-separated gene calls into sets.")
        for col in ['v_call', 'd_call', 'j_call']:
            df[col] = df[col].apply(lambda x: set(str(x).split(",")) if pd.notnull(x) else set())

        logger.success("Preprocessing completed successfully.")
        return df

    except Exception as e:
        logger.error(f"Error during preprocessing: {e}")
        raise