"""
author: Elisa Verhofstadt
studentnumber: 2261793
"""
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from typing import List
import sys

# read the csv file given its name


def read_csv(name: str) -> pd.DataFrame:
    '''Read in csv file and returns pandas DataFrame of the same information'''
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
    '''
    Filter out specified mistakes in the input DataFrame related to segments and positions.

    Input: DataFrame to be cleaned.

    Returns cleaned DataFrame without the specified errors.

    The function applies a series of checks and cleaning operations to remove the identified errors.
    It returns a new DataFrame that excludes the prespecified mistakes.

    '''
    # check everything that involves the position of the segments: missing, duplicate, errors
    df, to_remove_segment = _clean_positions(df)
    # check for double occurring segments
    to_remove_segment = _clean_segments(df, to_remove_segment)

    # remove selected segments
    df = df.query('SegmentNr not in @to_remove_segment')
    # reset index as several segments or positions have been removed
    return df.reset_index(drop=True)


def _clean_positions(df: pd.DataFrame):
    '''
    Clean positions in a DataFrame by performing checks for identifying and removing erroneous data.

    Input: DataFrame containing position information.

    Returns a cleaned DataFrame (in function of positions) and a list of segments to be removed.

    '''
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
    '''
    Identify duplicate segments from the input DataFrame.

    Input: DataFrame containing segment information and list of segments to be removed.

    Returns an updated list of segments to be removed.
    '''
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
    '''
    Convert the input DataFrame to a JSON string representation.

    Input: DataFrame containing sequences.

    Returns SON string representation of the DataFrame.
    '''
    return df.to_json(indent=0, orient='records')


# Construct de Bruijn graph

def construct_graph(json_data: str, k: int) -> nx.MultiDiGraph:
    ''' 
    Construct a de Bruijn graph from the DNA sequences provided in JSON format.

    Input: JSON string representing the DNA sequences and
    length of k-mers to be used in constructing the graph.

    Returns the constructed de Bruijn graph.
    '''
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
    '''
    Generate k-mers from the given DNA string.

    Input: DNA string from which k-mers are generated
    and length of the k-mers to be generated.

    Returns a list of generated k-mers.
    '''
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
    '''
    Retrieve the DNA string from the given DNA data.

    Input: list of dictionaries representing DNA data.

    Returns DNA string extracted from the DNA data.
    '''
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
    '''
    Plot the de Bruijn graph using NetworkX and save it to a file.

    Input: de Bruijn graph to be plotted and output filename to save the plot.
    '''
    pos = nx.planar_layout(graph)
    # use matplotlib make to plot
    plt.figure()
    nx.draw_networkx(graph, pos, with_labels=True)
    # Save the plot to the output file
    plt.savefig(filename)
    plt.close()

# Check whether the de Bruijn graph can be sequenced


def is_valid_graph(graph: nx.MultiDiGraph) -> bool:
    '''
    Check if a given de Bruijn graph contains a valid Euler's path based on connectivity and node degree.

    Input: de Bruijn graph to be checked.

    Returns True if the graph is valid, False if not.
    '''
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
    '''
    Perform a depth-first search (DFS) traversal on a given graph starting from a specified node.

    Input: Graph to perform DFS on, starting node for DFS traversal 
    and List to store visited nodes during DFS traversal.

    Returns list of visited nodes in the order they were visited.
    '''
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
    '''
    Construct a DNA sequence from a given graph using Eulerian path/circuit.

    Input: the input graph representing DNA sequences.

    Returns the constructed DNA sequence.
    '''

    # add extra edge if graph is not eulerian yet
    _make_eulerian_graph(graph)
    # apply hierhozer algorithm to get sequence
    sequence_index_list = _hierholzer_algorithm(
        _from_edges_to_matrix(graph), 0, [])
    # translate result into corresponding names of the nodes
    sequence_list = [list(graph.nodes())[index]
                     for index in sequence_index_list]

    # combine all nodes into full sequence
    full_sequence = ''
    for index, node in enumerate(sequence_list):
        if index == 0:
            full_sequence += node
        else:
            full_sequence += node[-1]
    return full_sequence


def _make_eulerian_graph(graph: nx.MultiDiGraph):
    '''
    Make the given graph Eulerian by adding an edge to connect the beginning and ending nodes if necessary.

    Input: Graph to check and possibly make Eulerian.

    Returns Eulertian graph
    '''
    # we know that graph is valid, so differ_degree is either 0 or 2
    # if 0, we have to do nothing
    # if 2, we have to connect beginning to ending node
    begin = ''
    end = ''
    for node in graph.nodes():
        if graph.in_degree(node) - graph.out_degree(node) == 1:
            end = node
        elif graph.in_degree(node) - graph.out_degree(node) == -1:
            begin = node
    # if begin and end are defined, if not they are still '' and return False
    if begin and end:
        graph.add_edge(end, begin)


def _from_edges_to_matrix(graph: nx.MultiDiGraph):
    '''
    Convert the given graph's edge representation into a matrix.

    Input: graph to convert.

    Returns matrix representing the connections between nodes.
    '''
    # make matrix with overview of connections between nodes
    edges_matrix = np.zeros(
        (graph.number_of_nodes(), graph.number_of_nodes()), dtype=int)
    for start_node, end_node in graph.edges():
        start_index = list(graph.nodes()).index(start_node)
        end_index = list(graph.nodes()).index(end_node)
        # directed graph: only conncection from start to end
        # this way rows represent nodes where the edge starts and
        # columns represents nodes where the edge ends
        # += 1 because one edge can be traversed multiple times
        edges_matrix[start_index][end_index] += 1
    return edges_matrix


def _hierholzer_algorithm(matrix, starting_row_index: int, full_sequence: list) -> list:
    '''
    Apply the Hierholzer's algorithm to find Eulerian circuits or paths in a graph.

    Input: Matrix representing the graph's connections, index of the row to start the algorithm
    and list to store the sequence of nodes.

    Returns list representing the Eulerian path.
    '''
    # identify all values in this row that have a 1; meaning that they are an
    # ending node of this edge
    ending_nodes_index = np.where(matrix[starting_row_index] > 0)[0]
    # use ending node as beginning node by using recursive function
    for ending_node_index in ending_nodes_index:
        # remove edge from the edge_matrix
        matrix[starting_row_index][ending_node_index] -= 1
        # recursive function: look at edges or ending node
        _hierholzer_algorithm(matrix, ending_node_index, full_sequence)

        # if we have completed a cycle and need to backtrack we are at the end
        # so we can add this node to our full_sequence list
        full_sequence.insert(0, starting_row_index)

        # stopping condition: if we have a matrix with all zeroes, all edges have been added
        if (matrix == np.zeros(matrix.shape, dtype=int)).all():
            # full_sequence.append(ending_node_index)
            return full_sequence


# Save DNA sequence or write the error message

def save_output(s: str, filename: str):
    '''
    Save the given string `s` to a text file with the specified `filename`.
    '''
    with open(f'{filename}.txt', mode='w') as dna_file:
        dna_file.write(s)


# for running the file from the command line
if __name__ == "__main__":
    input_file = sys.argv[1]
    _, x, k = input_file.split('_')
    k, _ = k.split('.')

    dna_dataframe = read_csv(input_file)
    dna_clean_dataframe = clean_data(dna_dataframe)
    dna_json = generate_sequences(dna_clean_dataframe)
    db_graph = construct_graph(dna_json, int(k))
    plot_graph(db_graph, f'DNA_{x}.png')
    if is_valid_graph(db_graph):
        to_write_string = construct_dna_sequence(db_graph)
        print(to_write_string)
    else:
        to_write_string = 'DNA sequence can not be constructed.'
    save_output(to_write_string, f'DNA_{x}')
