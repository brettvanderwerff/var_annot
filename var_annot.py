import numpy as np
import os
import pandas as pd
import requests
import sys
import json

'''
    var_annot
    ~~~~~~~~~

    A simple VCF annotation tool.
    
    Parses a VCF file and outputs an annotated VCF file with the following information for each variant:
    
    - 'TYPE' (The type of allele, either snp, mnp, ins, del, or complex)
    - 'DP' (Depth of sequence coverage at the site of variation)
    - 'AO' (Alternate allele observation count)
    - 'AF' (Estimated allele frequency in the range (0,1))
    - '#CHROM' (Chromosome that the variant is on)
    - 'POS' (Position of the variant on a chromosome)
    - 'REF' (Allele of the reference genome)
    - 'ALT' (Alternative allele detected by the variant caller)
    - 'ALLELE_FREQUENCY(ExAC)' (Allele frequency as reported by the ExAC Browser API (http://exac.hms.harvard.edu/))
    
'''

def convert_col_to_df(df, col, key_value_sep, elem_sep):
    '''
    Converts a column of a DataFrame that contains delimited key value pairs into a new
    DataFrame where the keys are column headers and the values are the rows.

    :param df: input DataFrame
    :param col: column to be converted into a DataFrame
    :param key_value_sep: key value pair delimiter
    :param elem_sep: column delimiter
    :return: DataFrame resulting from converting the column
    '''
    df = df[col].str.split(elem_sep, expand=True)
    key_values = df.iloc[0]
    header = [item.split(key_value_sep)[0] for item in key_values]
    df.columns = header
    df = df.applymap(lambda x: x.split(key_value_sep)[1])
    return df

def expand_df(df, columns, sep):
    '''

    Takes columns of a DataFrame that contain multiple delimited values and inserts these individual values as new rows.
    Mating values in columns that do not contain delimited values are copied directly to the new rows.

    :param df: input DataFrame
    :param columns: list of columns that contain delimited values
    :param sep: value delimiter
    :return: DataFrame with new rows for each delimited value
    '''
    print('Processing VCF data. This may take several minutes.')
    new_df = pd.DataFrame(columns=columns)
    for name, row in df.iterrows():
        columns_dict = {}
        for column in columns:
            columns_dict[column] = row[column].split(sep)
        zipped = zip(*columns_dict.values())
        for item in zipped:
            for counter, column in enumerate(columns_dict.keys()):
                row[column] = item[counter]
            new_df = new_df.append(row.apply(str)) #ToDo Pandas converts dtype int to float at this step, str is workaround
    return new_df.reset_index(drop=True)

def build_cprv_col(df, chromosome, position, reference, variant):
    '''
    Builds url endpoints to communicate with The Exome Aggregation Consortium (ExAC) Browser API
    (http://exac.hms.harvard.edu/). url endpoints need to be in the format: 'CHROMOSOME-POSITION-REFERENCE-VARIANT'
    (CPRV).

    :param df: input DataFrame
    :param chromosome: DataFrame column for the chromosome value
    :param position: DataFrame column for the position value
    :param reference: DataFrame column for the reference allele value
    :param variant: DataFrame column for the alternative allele value
    :return: a DataFrame with a new 'CPRV' column for the 'CHROMOSOME-POSITION-REFERENCE-VARIANT' data
    '''
    columns = [chromosome, position, reference, variant]
    df['CPRV'] = df[columns[0]]
    for column in columns[1:]:
        df['CPRV'] += '-' + df[column].apply(str)
    return df

def call_exac_server(df):
    '''
    Makes a call to the The Exome Aggregation Consortium (ExAC) Browser API (http://exac.hms.harvard.edu/). Selects the
     'CPRV' column from the input DataFrame. Converts the information in this column to JSON. Then makes a bulk
     query to the ExAC Browser API.

    :param df: input DataFrame
    :return: JSON data containing variant information from the ExAC Browser API
    '''
    print('Getting variant data from the ExAC server. This may take take several minutes')
    endpoint_column = df['CPRV'].to_list()
    endpoint_column_json = json.dumps(endpoint_column)
    response = requests.post('http://exac.hms.harvard.edu/rest/bulk/variant/variant', data=endpoint_column_json)
    return response.json()

def get_allele_freqs(json_response):
    '''
    Parses JSON data returned from a query of The Exome Aggregation Consortium (ExAC) Browser API
     (http://exac.hms.harvard.edu/). Isolates allele frequency data from this response.

    :param json_response: JSON data from the API call
    :return: a list where each element represents allele frequencies
    '''
    allele_freqs = []
    for item in json_response.values():
        try:
            allele_freq = item['allele_freq']
            allele_freqs.append(allele_freq)
        except KeyError as e:
            allele_freqs.append(np.nan)
    return allele_freqs

def get_cprv(json_response):
    '''
    Parses JSON data returned from a query of The Exome Aggregation Consortium (ExAC) Browser API
     (http://exac.hms.harvard.edu/). Isolates the 'CHROMOSOME-POSITION-REFERENCE-VARIANT' (CPRV) data that was used to
     make the call.

    :param json_response: JSON data from the API call
    :return: a list where each element represents 'CHROMOSOME-POSITION-REFERENCE-VARIANT' (CPRV) data
    '''
    return [item for item in json_response]

def get_vcf_filename(file_path):
    '''
    Gets the name of a VCF file from a path.

    :param file_path: Path to VCF file
    :return: VCF filename without extension
    '''
    base = os.path.basename(file_path)
    return os.path.splitext(base)[0]

def write_output(file_name, df):
    '''

    Writes an annotated VCF file to the output folder.

    :param file_name: file name of the VCF file
    :param df: input DataFrame
    '''
    print('Writing annotated VCF file to the "output" folder.')
    prefix = './output'
    filename_with_extension = file_name + '_annotated.vcf'
    write_path  = '/'.join([prefix, filename_with_extension])
    df = df[['#CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'DP', 'AO', 'AF', 'ALLELE_FREQUENCY(ExAC)']] # Re-order to reflect original VCF file
    df.to_csv(write_path, sep='\t', index=False)

def subset_vcf(df):
    '''

    Subsets a DataFrame containing VCF data to retain only the following information:

    - 'TYPE' (The type of allele, either snp, mnp, ins, del, or complex)
    - 'DP' (Depth of sequence coverage at the site of variation)
    - 'AO' (Alternate allele observation count)
    - 'AF' (Estimated allele frequency in the range (0,1))
    - '#CHROM' (Chromosome that the variant is on)
    - 'POS' (Position of the variant on a chromosome)
    - 'REF' (Allele of the reference genome)
    - 'ALT' (Alternative allele detected by the variant caller)

    :param df: input DataFrame
    :return: subset of input DataFrame
    '''
    info = convert_col_to_df(df, 'INFO', '=', ';')
    info = info[['DP', 'AO', 'AF', 'TYPE']]
    df = df[['#CHROM', 'POS', 'REF', 'ALT']]
    df = df.join(info)
    df = expand_df(df, ['TYPE', 'ALT', 'AF'], ',')
    return df

def get_exac_data(df):
    '''
    Handles making a call to The Exome Aggregation Consortium (ExAC) Browser API (http://exac.hms.harvard.edu/). This
    call is made to obtain allele frequency data.

    :param df: input DataFrame
    :return: DataFrame with a new ALLELE_FREQUENCY(ExAC) column containing allele frequency data.
    '''
    df = build_cprv_col(df, '#CHROM', 'POS', 'REF', 'ALT')
    json_response = call_exac_server(df)
    allele_freqs = get_allele_freqs(json_response)
    cprv = get_cprv(json_response)
    response_df = pd.DataFrame({'ALLELE_FREQUENCY(ExAC)': allele_freqs,
                           'CPRV': cprv})
    df = df.merge(response_df, on='CPRV')
    df = df.drop(['CPRV'], axis=1)  # CPRV column was only needed to ensure correct merging.
    return df

def run(file_path, header_row):
    '''
    Main function to run the the var_annot program.

    :param file_path: Path to the VCF file that is to be annotated.
    :param header_row: Header row for the VCF file, usually this is the first row following the meta-information rows.
    '''
    df = pd.read_csv(file_path, skiprows=header_row, sep='\t')
    df = subset_vcf(df)
    df = get_exac_data(df)
    file_name = get_vcf_filename(file_path)
    write_output(file_name, df)

if __name__ ==  '__main__':
    file_path = sys.argv[1] # Get first command line argument for the file_path
    header_row = int(sys.argv[2]) # Get second command line argument for the header_row
    run(file_path, header_row)


