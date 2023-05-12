from project import clean_data
from pytest import mark
import pandas as pd
import networkx as nx


@mark.parametrize(
    'dna_df, expected',
    [
        (
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 1],
                [1, 2, 0, 0, 0, 1],
                [2, 1, 1, 0, 0, 0],
                [2, 2, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
                [2, 1, 1, 0, 0, 0],
                [2, 2, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
        )
    ],

)
def test_clean_data(dna_df: pd.DataFrame, expected: pd.DataFrame) -> None:
    assert clean_data(dna_df).equals(expected)


# @mark.parametrize(
#     'dna_df, expected_json_str',
#     [
#         # not given, you have to choose your own JSON structure
#     ], )
# def test_generate_sequences(dna_df: pd.DataFrame, expected_json_str: str) -> None:
#     assert (generate_sequences(dna_df) == expected_json_str)


# @mark.parametrize(
#     'json_data, k,  expected_edge_list',
#     [
#         # create a JSON in your own format that contains one segment with sequence "ATTACTC"
#         # for k = 5, it should output a graph with the edges: [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('TACT', 'ACTC')]
#     ])
# def test_construct_graph(json_data: str, k: int, expected_edge_list: list) -> None:
#     pass

# @mark.parametrize(
#     'DNA_edge_list,  expected_validity, tname',
#     [
#         (
#                 [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('TACT', 'ACTC'), ('ACTC', 'ATTA')],
#                 True
#         )
#     ])
# def test_is_valid_graph(DNA_edge_list: list, expected_validity: bool) -> None:
#     debruijn_graph = nx.MultiDiGraph()
#     for edge in DNA_edge_list:
#         debruijn_graph.add_edge(edge[0], edge[1])

#     assert is_valid_graph(debruijn_graph) is expected_validity


# @mark.parametrize(
#     'DNA_edge_list,  possible_dna_sequence',
#     [
#         (
#                 [('AAA', 'AAC'), ('AAC', 'ACA'), ('ACA', 'CAC')],
#                 ["AAACAC"]
#         )
#     ])
# def test_construct_dna_sequence(DNA_edge_list: list, possible_dna_sequence) -> None:
#     debruijn_graph = nx.MultiDiGraph()
#     for edge in DNA_edge_list:
#         debruijn_graph.add_edge(edge[0], edge[1])

#     assert construct_dna_sequence(debruijn_graph) in possible_dna_sequence
