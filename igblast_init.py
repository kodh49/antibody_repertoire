import pandas as pd

# preprocess the ncbi-igblast result

def igblast_preprocess(to_dict: False):
    df = pd.read_csv('data/igblast_results.tsv', sep='\t', header=0)
    assorted_columns = ['v_call', 'd_call', 'j_call', 'cdr3', 'cdr3_start', 'cdr3_end']
    df.drop(df.columns.difference(assorted_columns), axis=1, inplace=True)

    # normalize matching gene information by ignoring alleles 
    df['v_call'] = df['v_call'].replace( { r"\*\d+" : '' }, regex = True) # ignore alleles
    df['d_call'] = df['d_call'].replace( { r"\*\d+" : '' }, regex = True) # ignore alleles
    df['j_call'] = df['j_call'].replace( { r"\*\d+" : '' }, regex = True) # ignore alleles
    
    # compute CDR3 region length
    df['cdr_length'] = df['cdr3_end'] - df['cdr3_start'] + 1
    df.drop(['cdr3_start', 'cdr3_end'], axis=1, inplace=True)

    # convert comma separated matching genes into list
    df['v_call'] = df['v_call'].apply(lambda x: set(str(x).split(",")))
    df['d_call'] = df['d_call'].apply(lambda x: set(str(x).split(",")))
    df['j_call'] = df['j_call'].apply(lambda x: set(str(x).split(",")))

    if to_dict == True:
        # convert to dictionary format
        dictionary = df.to_dict('index')
        return dictionary
    
    # return original dataframe otherwise
    return df