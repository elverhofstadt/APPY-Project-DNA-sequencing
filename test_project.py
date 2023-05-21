from project import clean_data, generate_sequences, construct_graph, is_valid_graph
from pytest import mark
import pandas as pd
import networkx as nx


@mark.parametrize(
    'dna_df, expected',
    [
        (
            # wrong position: more than one value 1
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
        ),
        (
            # wrong position: all values are 0's
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 0],
                [1, 2, 0, 0, 0, 1],
                [2, 1, 0, 0, 0, 0],
                [2, 2, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 0],
                [1, 2, 0, 0, 0, 1]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
        ),
        (
            # missing position: skipping one
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 0],
                [1, 2, 0, 0, 0, 1],
                [1, 4, 0, 1, 0, 0],
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])

        ),
        (
            # duplicate position: same values
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 0],
                [1, 2, 0, 0, 0, 1],
                [1, 3, 0, 1, 0, 0],
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 1, 0, 0],
                [2, 2, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 0],
                [1, 2, 0, 0, 0, 1],
                [1, 3, 0, 1, 0, 0],
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
        ),
        (
            # duplicate position: different values
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 0],
                [1, 2, 0, 0, 0, 1],
                [1, 3, 0, 1, 0, 0],
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 1, 0, 0],
                [2, 2, 1, 0, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 0],
                [1, 2, 0, 0, 0, 1],
                [1, 3, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
        ),
        (
            # duplicate segments (segment 1 and 5)
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 0],
                [1, 2, 0, 0, 0, 1],
                [1, 3, 0, 1, 0, 0],
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 1, 0, 0],
                [2, 3, 1, 0, 0, 0],
                [3, 1, 0, 0, 0, 1],
                [4, 1, 0, 1, 0, 0],
                [4, 2, 1, 0, 0, 0],
                [5, 1, 1, 0, 0, 0],
                [5, 2, 0, 0, 0, 1],
                [5, 3, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 1, 0, 0],
                [2, 3, 1, 0, 0, 0],
                [3, 1, 0, 0, 0, 1],
                [4, 1, 0, 1, 0, 0],
                [4, 2, 1, 0, 0, 0],
                [5, 1, 1, 0, 0, 0],
                [5, 2, 0, 0, 0, 1],
                [5, 3, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
        ),
        (

            # combination of different error (each segment has one error)
            pd.DataFrame(data=[
                [1, 1, 1, 0, 0, 1],
                [1, 2, 0, 0, 0, 1],
                [1, 3, 0, 1, 0, 0],
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 0, 0, 0],
                [2, 3, 1, 0, 0, 0],
                [3, 1, 0, 0, 0, 1],
                [3, 1, 0, 0, 0, 1],
                [4, 1, 0, 1, 0, 0],
                [4, 2, 1, 0, 0, 0],
                [4, 2, 0, 1, 0, 0],
                [4, 3, 0, 0, 0, 1],
                [5, 1, 1, 0, 0, 0],
                [5, 2, 0, 0, 0, 1],
                [5, 3, 0, 1, 0, 0],
                [6, 1, 1, 0, 0, 0],
                [6, 2, 0, 0, 0, 1],
                [6, 3, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
                [3, 1, 0, 0, 0, 1],
                [6, 1, 1, 0, 0, 0],
                [6, 2, 0, 0, 0, 1],
                [6, 3, 0, 1, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
        ),

        (
            # segments with multiple errors
            pd.DataFrame(data=[
                [1, 1, 1, 1, 0, 1],
                [1, 2, 0, 0, 0, 1],
                [1, 2, 0, 1, 0, 0],
                [2, 1, 0, 0, 0, 1],
                [2, 2, 0, 1, 0, 0],
                [2, 2, 0, 1, 0, 0],
                [2, 3, 1, 0, 0, 0],
                [3, 1, 0, 0, 0, 1],
                [3, 2, 0, 1, 0, 0],
                [3, 2, 0, 1, 0, 0],
                [3, 3, 1, 0, 0, 0],
                [4, 5, 0, 0, 0, 1],
                [4, 5, 0, 0, 0, 1]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
                [3, 1, 0, 0, 0, 1],
                [3, 2, 0, 1, 0, 0],
                [3, 3, 1, 0, 0, 0]],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
        ),

        (
            # empty dataframe
            pd.DataFrame(data=[],
                         columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[],
                         columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])
        ),
        (
            # no errors
            pd.DataFrame(data=[
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
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
            pd.DataFrame(data=[
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
        )
    ],

)
def test_clean_data(dna_df: pd.DataFrame, expected: pd.DataFrame) -> None:
    assert clean_data(dna_df).equals(expected)


@mark.parametrize(
    'dna_df, expected_json_str',
    [(
        pd.DataFrame(data=[
            [1, 1, 0, 0, 0, 1],
            [1, 2, 0, 0, 0, 1]],
            columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
        '[{"SegmentNr":1,"Position":1,"A":0,"C":0,"G":0,"T":1},' +
        '{"SegmentNr":1,"Position":2,"A":0,"C":0,"G":0,"T":1}]'
    ),
        (
        # empty dataframe
        pd.DataFrame(data=[],
                     columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
        '[]'
    ),
        (
        # multiple segments
        pd.DataFrame(data=[
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
            columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
        '[{"SegmentNr":1,"Position":1,"A":0,"C":0,"G":0,"T":1},' +
        '{"SegmentNr":1,"Position":2,"A":0,"C":0,"G":0,"T":1},' +
        '{"SegmentNr":1,"Position":3,"A":0,"C":1,"G":0,"T":0},' +
        '{"SegmentNr":2,"Position":1,"A":0,"C":0,"G":0,"T":1},' +
        '{"SegmentNr":2,"Position":2,"A":0,"C":1,"G":0,"T":0},' +
        '{"SegmentNr":2,"Position":3,"A":0,"C":1,"G":0,"T":0},' +
        '{"SegmentNr":2,"Position":4,"A":1,"C":0,"G":0,"T":0},' +
        '{"SegmentNr":3,"Position":1,"A":0,"C":0,"G":1,"T":0},' +
        '{"SegmentNr":3,"Position":2,"A":1,"C":0,"G":0,"T":0},' +
        '{"SegmentNr":3,"Position":3,"A":0,"C":1,"G":0,"T":0},' +
        '{"SegmentNr":3,"Position":4,"A":1,"C":0,"G":0,"T":0},' +
        '{"SegmentNr":4,"Position":1,"A":0,"C":0,"G":0,"T":1},' +
        '{"SegmentNr":5,"Position":1,"A":0,"C":0,"G":0,"T":1},' +
        '{"SegmentNr":5,"Position":2,"A":0,"C":1,"G":0,"T":0}]'

    )
    ],
)
def test_generate_sequences(dna_df: pd.DataFrame, expected_json_str: str) -> None:
    assert (generate_sequences(dna_df) == expected_json_str)


@mark.parametrize(
    'json_data, k,  expected_edge_list',
    [
        # create a JSON in your own format that contains one segment with sequence "ATTACTC"
        # for k = 5, it should output a graph with the edges: [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('TACT', 'ACTC')]
        (
            '[{"SegmentNr":1,"Position":1,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":2,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":3,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":4,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":5,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":6,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":7,"A":0,"C":1,"G":0,"T":0}]',
            5,
            [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('TACT', 'ACTC')]
        ),
        # multiple sequences without overlapping k-mers
        # ATGCTAA, CTCAATA
        (
            '[{"SegmentNr":1,"Position":1,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":2,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":3,"A":0,"C":0,"G":1,"T":0},' +
            '{"SegmentNr":1,"Position":4,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":5,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":6,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":7,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":1,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":2,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":2,"Position":3,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":4,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":5,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":6,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":2,"Position":7,"A":1,"C":0,"G":0,"T":0}]',
            4,
            [('ATG', 'TGC'), ('TGC', 'GCT'), ('GCT', 'CTA'), ('CTA', 'TAA'), ('CTC', 'TCA'),
             ('TCA', 'CAA'), ('CAA', 'AAT'), ('AAT', 'ATA')]
        ),
        # multiple sequences with overlapping sequences
        # ATGCTAA, CTCGCTA
        (
            '[{"SegmentNr":1,"Position":1,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":2,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":3,"A":0,"C":0,"G":1,"T":0},' +
            '{"SegmentNr":1,"Position":4,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":5,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":6,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":7,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":1,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":2,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":2,"Position":3,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":4,"A":0,"C":0,"G":1,"T":0},' +
            '{"SegmentNr":2,"Position":5,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":6,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":2,"Position":7,"A":1,"C":0,"G":0,"T":0}]',
            4,
            [('ATG', 'TGC'), ('TGC', 'GCT'), ('GCT', 'CTA'), ('CTA', 'TAA'), ('CTC', 'TCG'),
             ('TCG', 'CGC'), ('CGC', 'GCT'), ('GCT', 'CTA')]
        ),
        # k of 2 (nodes of one letter) + node with connection to itself
        # ATGCG, TAA
        (
            '[{"SegmentNr":1,"Position":1,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":2,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":3,"A":0,"C":0,"G":1,"T":0},' +
            '{"SegmentNr":1,"Position":4,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":5,"A":0,"C":0,"G":1,"T":0},' +
            '{"SegmentNr":2,"Position":1,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":2,"Position":2,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":2,"Position":3,"A":1,"C":0,"G":0,"T":0}]',
            2,
            [('A', 'T'), ('T', 'G'), ('G', 'C'),
             ('C', 'G'), ('T', 'A'), ('A', 'A')]
        ),
        # segment with repeating edges within itself
        # ATGACTGAA
        (
            '[{"SegmentNr":1,"Position":1,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":2,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":3,"A":0,"C":0,"G":1,"T":0},' +
            '{"SegmentNr":1,"Position":4,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":5,"A":0,"C":1,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":6,"A":0,"C":0,"G":0,"T":1},' +
            '{"SegmentNr":1,"Position":7,"A":0,"C":0,"G":1,"T":0},' +
            '{"SegmentNr":1,"Position":8,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":9,"A":1,"C":0,"G":0,"T":0}]',
            3,
            [('AT', 'TG'), ('TG', 'GA'), ('GA', 'AC'), ('AC', 'CT'),
             ('CT', 'TG'), ('TG', 'GA'), ('GA', 'AA')]
        ),
        # segment consisting of same letter
        # AAAAAAA
        (
            '[{"SegmentNr":1,"Position":1,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":2,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":3,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":4,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":5,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":6,"A":1,"C":0,"G":0,"T":0},' +
            '{"SegmentNr":1,"Position":7,"A":1,"C":0,"G":0,"T":0}]',
            4,
            [('AAA', 'AAA'), ('AAA', 'AAA'), ('AAA', 'AAA'), ('AAA', 'AAA')]
        )
    ])
def test_construct_graph(json_data: str, k: int, expected_edge_list: list) -> None:
    assert (sorted(list(construct_graph(json_data, k).edges()))
            ) == sorted(expected_edge_list)


@mark.parametrize(
    'DNA_edge_list,  expected_validity',
    [
        (
            [('ATTA', 'TTAC'), ('TTAC', 'TACT'),
             ('TACT', 'ACTC'), ('ACTC', 'ATTA')],
            True
        ),
        # weakly connected, but in and out degree is not correct (4 different in and out)
        (
            [('ATT', 'TTA'), ('TTA', 'TAG'), ('TAG', 'AGT'),
             ('TAG', 'AGA'), ('AGA', 'GAA'), ('AGA', 'GAT'),
             ('GAT', 'ATT')],
            False
        ),

        # not weakly connected, but in and out degree is correct
        (
            [('ATT', 'TTA'), ('TTA', 'TAG'), ('TAG', 'GAT'), ('GAT', 'ATT'),
             ('GTC', 'TCA'), (('TCA', 'CAC'))],
            False
        ),

        # not weakly connected, in and out degree is not correct
        (
            [('ATT', 'TTA'), ('TTA', 'TAG'), ('TTA', 'TAT'),
             ('TAG', 'GAT'), ('GAT', 'ATT'),
             ('GTC', 'TCA'), (('TCA', 'CAC'))],
            False
        ),

        # segment that repeats within itself
        (
            [('ATT', 'TTA'), ('TTA', 'TAG'), ('TAG', 'AGT'),
             ('AGT', 'GTA'), ('GTA', 'TTA')],
            True
        )

    ])
def test_is_valid_graph(DNA_edge_list: list, expected_validity: bool) -> None:
    debruijn_graph = nx.MultiDiGraph()
    for edge in DNA_edge_list:
        debruijn_graph.add_edge(edge[0], edge[1])

    assert is_valid_graph(debruijn_graph) is expected_validity


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
