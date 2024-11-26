import pandas as pd

# preprocess the ncbi-igblast result

def igblast_preprocess():
    df = pd.read_csv('data/igblast_results.tsv', sep='\t', header=0)
    assorted_columns = ['v_call', 'd_call', 'j_call', 'cdr3', 'cdr3_start', 'cdr3_end']
    df.drop(df.columns.difference(assorted_columns), axis=1, inplace=True)      
    df['v_call'] = df['v_call'].replace( { r"\*\d+" : '' }, regex = True) # ignore alleles
    df['d_call'] = df['d_call'].replace( { r"\*\d+" : '' }, regex = True) # ignore alleles
    df['j_call'] = df['j_call'].replace( { r"\*\d+" : '' }, regex = True) # ignore alleles
    df['cdr_length'] = df['cdr3_end'] - df['cdr3_start'] + 1

    # drop cdr3 start and end
    df.drop(['cdr3_start', 'cdr3_end'], axis=1, inplace=True)
    dictionary = df.to_dict('index')
    print(df)
    return dictionary

igblast_preprocess()