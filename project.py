"""
author: Elisa Verhofstadt
studentnumber: 2261793
"""
import networkx as nx
import pandas as pd
from typing import List


def read_csv(name: str) -> pd.DataFrame:
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
def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    # check everything that involves the position of the segments: missing, duplicate, errors
    df, to_remove_segment = _clean_positions(df)
    # check for double occurring segments
    to_remove_segment = _clean_segments(df, to_remove_segment)

    # remove selected segments
    df = df.query('SegmentNr not in @to_remove_segment')
    return df


def _clean_positions(df: pd.DataFrame):
    # identify all the unique sequence numbers
    segments = df['SegmentNr'].unique()

    # for each segment check if all positions are included
    to_remove_index = []
    to_remove_segment = []

    index = 0
    for segment in segments:
        segment_index = 0
        position_count = 1

        # select df of all positions of this segment
        segment_df = df[df['SegmentNr'] == segment]

        # go through each position in a segment and do the necessary checks
        for _, row in segment_df.iterrows():
            # indicates missing position, e.g. 1, 2, 4, 5
            if (row['Position'] > position_count) and row['SegmentNr'] not in to_remove_segment:
                to_remove_segment += [row['SegmentNr']]
            # indicated duplicate position, e.g. 1, 2, 2, 3
            elif row['Position'] < position_count:
                # check whether duplicate positions have equal values
                if list(segment_df.iloc[segment_index - 1]) == list(row):
                    # same information, remove one of the positions
                    to_remove_index += [index]
                elif row['SegmentNr'] not in to_remove_segment:
                    # different information, discard segment
                    to_remove_segment += [row['SegmentNr']]
            else:  # correct position
                # check for errors in every row (all zeros or multiple ones)
                # select all zero and one values for one entry
                values = list(row[['A', 'C', 'G', 'T']])
                if (not any(values) or values.count(1) > 1) and row['SegmentNr'] not in to_remove_segment:
                    to_remove_segment += [row['SegmentNr']]
                position_count += 1
            index += 1
            segment_index += 1
    # remove selected positions within segments
    df = df.drop(to_remove_index)

    return df, to_remove_segment


def _clean_segments(df: pd.DataFrame, to_remove_segment: list) -> list:
    #
    segments = df['SegmentNr'].unique()
    segment_list = []
    for segment in segments:
        # for each segment make list of lists with every measurement, don't include segment number
        segment_list += [df.query('SegmentNr == @segment')
                         [['A', 'C', 'G', 'T']].values.tolist()]

    # make sure that for every combination of entries is only counted once
    for index, list_values in enumerate(segment_list):
        if segment_list.count(list_values) > 1:
            # record the segment number
            to_remove_segment += [segments[index]]
            # remove it from the list such that it will not be removed twice
            segment_list[index] = ''
    return to_remove_segment


# Generate JSON sequences from the dataframe

# Construct de Bruijn graph

# Plot the de Bruijn graph

# Check whether the de Bruijn graph can be sequenced

# Construct DNA sequence

# Save DNA sequence or write the error message
# testing
if __name__ == "__main__":
    #df_1 = read_csv('DNA_1_5.csv')
    df_clean = read_csv('DNA_1_5_to_clean.csv')
    print(clean_data(df_clean))
