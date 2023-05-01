"""
author: Elisa Verhofstadt
studentnumber: 2261793
"""
import networkx as nx
import pandas as pd
from typing import List


def read_csv(name: str):
    '''Reads in csv file and returns pandas DataFrame'''
    with open(name, mode='r') as dna_file:
        # every reading as string in a list
        df_as_list = dna_file.readlines()

        # turn strings into lists of values
        for index_r, reading in enumerate(df_as_list):
            # split values, remove whitespace
            df_as_list[index_r] = reading.strip().split(',')
            # record every values as integers
            for index_v, value in enumerate(df_as_list[index_r]):
                df_as_list[index_r][index_v] = eval(value)

        # return list of lists as pandas DataFrame
        return pd.DataFrame(df_as_list, columns=[
            'SegmentNr', 'Position', 'A', 'C', 'G', 'T'])


# Clean the dataframe

# Generate JSON sequences from the dataframe

# Construct de Bruijn graph

# Plot the de Bruijn graph

# Check whether the de Bruijn graph can be sequenced

# Construct DNA sequence

# Save DNA sequence or write the error message

if __name__ == "__main__":
    print(read_csv('DNA_1_5.csv'))
