"""
author: Elisa Verhofstadt
studentnumber: 2261793
"""
import networkx as nx
from networkx.drawing import planar_layout
import pandas as pd
import matplotlib.pyplot as plt
import json
from typing import List

# read the csv file given its name


def read_csv(name: str) -> pd.DataFrame:
    '''Reads in csv file and returns pandas DataFrame'''
    with open(name, mode='r') as dna_file:
        # every reading as string in a list
        df_as_list = dna_file.readlines()

        # turn strings into lists of values
        for index_r, reading in enumerate(df_as_list):
            # split values, remove whitespace
            df_as_list[index_r] = reading.strip().split(',')
            # record every value as integers
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
    # reset index as several segments or positions have been removed
    return df.reset_index(drop=True)


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
        # for each segment make list of lists with every measurement,
        # don't include segment number
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

def generate_sequences(df: pd.DataFrame) -> str:
    return df.to_json(indent=0, orient='records')


# Construct de Bruijn graph


def construct_graph(json_data: str, k: int) -> nx.MultiDiGraph:
    ''' TO FINISH '''
    # initiate graph
    de_Bruij_G = nx.MultiDiGraph()

    # turn json object into dictionary
    dna_dict = json.loads(json_data)

    # dna_dict is now list with dict for each entry
    # select all segment numbers
    segment_nrs = set([position['SegmentNr'] for position in dna_dict])
    # add nodes for each sequence
    for segment_nr in segment_nrs:
        # reconstruct DNA structure single segment
        single_segment = [
            position for position in dna_dict if position['SegmentNr'] == segment_nr]
        dna_str = _get_dna_string(single_segment)

        # get all k-mers
        k_mer_list = _generate_k_mers(dna_str, k)

        # for each k-mer get two (k-1) mers -> L and R
        for k_mer in k_mer_list:
            left = k_mer[:-1]  # L
            right = k_mer[1:]  # R

            # for each L and R add a node and edge from left to right
            de_Bruij_G.add_edge(left, right)
    return de_Bruij_G


def _generate_k_mers(dna_str: str, k: int) -> list:
    ''' TO FINISH'''
    k_mers = []
    # determine how many k-mers we will return
    n_k_mer = len(dna_str) - k + 1

    # slice dna_str n_k_mer times, each time changing one place
    n_begin = 0
    n_end = len(dna_str) - n_k_mer + 1
    for _ in range(n_k_mer):
        k_mers += [dna_str[n_begin:n_end]]
        # update slicers
        n_begin += 1
        n_end += 1
    return k_mers


def _get_dna_string(dna_data: List[dict]) -> str:
    ''' TO FINISH'''
    dna_str = []
    # loop through all positions
    for line in dna_data:
        # each line is a dictionary
        # select only the keys where value is 1 for dna letter
        dna_str += [key for (key, value) in line.items() if value == 1 and
                    (key != 'SegmentNr' and key != 'Position')]
    return ''.join(dna_str)


# Plot the de Bruijn graph
def plot_graph(graph: nx.MultiDiGraph, filename: str) -> None:
    pos = nx.planar_layout(graph)
    # use matplotlib make to plot
    plt.figure()
    nx.draw_networkx(graph, pos, with_labels=True)
    # Save the plot to the output file
    plt.savefig(filename)
    plt.close()

# Check whether the de Bruijn graph can be sequenced


def is_valid_graph(graph: nx.MultiDiGraph) -> bool:
    ''' TO FINISH'''
    connectivity_check = True
    # store whether or not a first or last node has been identified
    first_node = False
    last_node = False
    # keep track of how many nodes with different degrees have been encountered
    differ_degree = 0
    for node in list(graph.nodes()):
        # to pass for the connectivity test the result of the dfs should be equal
        # to all the edges in the graph
        if sorted(_dfs_recursive(graph, node, [])) != sorted(list(graph.nodes())):
            connectivity_check = False
        # check in and out degree of this node
        if graph.in_degree(node) != graph.out_degree(node):
            differ_degree += 1
            if graph.in_degree(node) - graph.out_degree(node) == 1:
                last_node = True
            elif graph.out_degree(node) - graph.in_degree(node) == 1:
                first_node = True
    # second condition check: if nodes with different degrees are present one should be first and other one last
    degree_check = (last_node and first_node and differ_degree ==
                    2) or differ_degree == 0
    # only returns true is both conditions are met
    return connectivity_check and degree_check


def _dfs_recursive(graph: nx.MultiDiGraph, starting_node: str, visited_list: list) -> list:
    '''TO FINISH'''
    if visited_list is None:
        visited_list = []
    visited_list.append(starting_node)
    # loop through all neigbours (both directions -> weakly connected)
    for neigbour_node in set(list(graph.neighbors(starting_node)) + list(graph.predecessors(starting_node))):
        if neigbour_node not in visited_list:
            _dfs_recursive(graph, neigbour_node, visited_list)
    # only reaches here when all possible nodes are traversed
    return visited_list

# Construct DNA sequence


def construct_dna_sequence(graph: nx.MultiDiGraph):
    # check if graph is valid
    if is_valid_graph(graph):
        # add extra edge if graph is not eulerian yet
        _make_eulerian_graph(graph)
    else:
        # idk yet what i should do here/if this should be here
        print('you enetred a graph where no path can be found')


def _make_eulerian_graph(graph: nx.MultiDiGraph):
    '''TO FINISH'''
    # we know that graph is valid, so differ_degree is either 0 or 2
    # if 0, we have to do nothing
    # if 2, we have to connect beginning to ending node
    begin = ''
    end = ''
    for node in graph.nodes():
        if node.in_degree() - node.out_degree() == 1:
            end = node
        elif node.in_degree() - node.out_degree() == -1:
            begin = node
    # if begin and end are defined, if not they are still '' and return False
    if begin and end:
        graph.add_edge(begin, end)

# Save DNA sequence or write the error message


# testing
if __name__ == "__main__":
    # df_1 = read_csv('DNA_1_5.csv')
    df_clean = pd.DataFrame(data=[
        [1, 1, 0, 0, 0, 1],
        [1, 2, 0, 0, 0, 1],
        [1, 3, 0, 1, 0, 0],
        [2, 1, 0, 0, 0, 1],
        [2, 2, 0, 1, 0, 0],
        [2, 3, 0, 1, 0, 0],
        [2, 4, 1, 0, 0, 0],
        [3, 1, 0, 0, 1, 0],
        [3, 2, 1, 0, 0, 0],
        [3, 3, 0, 1, 0, 0],
        [3, 4, 1, 0, 0, 0],
        [4, 1, 0, 0, 0, 1],
        [5, 1, 0, 0, 0, 1],
        [5, 2, 0, 1, 0, 0]],
        columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
    json_file = generate_sequences(clean_data(df_clean))
    G = construct_graph(json_file, 3)
    json_file2 = '[{"SegmentNr":1,"Position":1,"A":1,"C":0,"G":0,"T":0},' + '{"SegmentNr":1,"Position":2,"A":0,"C":0,"G":0,"T":1},' + '{"SegmentNr":1,"Position":3,"A":0,"C":0,"G":1,"T":0},' + '{"SegmentNr":1,"Position":4,"A":1,"C":0,"G":0,"T":0},' + \
        '{"SegmentNr":1,"Position":5,"A":0,"C":1,"G":0,"T":0},' + '{"SegmentNr":1,"Position":6,"A":0,"C":0,"G":0,"T":1},' + \
        '{"SegmentNr":1,"Position":7,"A":0,"C":0,"G":1,"T":0},' + \
        '{"SegmentNr":1,"Position":8,"A":1,"C":0,"G":0,"T":0},' + \
        '{"SegmentNr":1,"Position":9,"A":1,"C":0,"G":0,"T":0}]'
    G2 = construct_graph(json_file2, 3)
