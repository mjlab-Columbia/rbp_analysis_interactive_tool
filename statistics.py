from pandas import DataFrame
from data import InteractionTypeValue
from networkx import betweenness_centrality, Graph

from typing import List, Tuple


def get_number_of_nodes(df: DataFrame) -> int:
    baits = set(df["Baits"].unique())  # Set[str]
    preys = set(df["Preys"].unique())  # Set[str]

    # Total nodes is union of baits and preys
    num_nodes = len(baits | preys)  # int

    return num_nodes


def get_number_of_interactions(df: DataFrame) -> int:
    # Total number of interactions = number of bait-prey entries in df
    return df.shape[0]


def get_number_of_direct_interactions(df: DataFrame) -> int:
    direct = InteractionTypeValue.direct.value  # str
    direct_df = df[df["IP_interaction_type"] == direct]  # DataFrame
    direct_interactions = direct_df.shape[0]  # int
    return direct_interactions


def get_number_of_shared_interactions(df: DataFrame) -> int:
    return 0


def get_top_betweenness_proteins(df: DataFrame,
                                 graph: Graph,
                                 num_proteins: int) -> List[Tuple[str, float]]:
    bw_dict = betweenness_centrality(graph)  # Dict[str, float]
    sorted_bw = sorted(bw_dict.items(),
                       key=lambda x: x[1])  # List[Tuple[str, float]]
    sorted_bw_subset = sorted_bw[:num_proteins]  # List[Tuple[str, float]]
    return sorted_bw_subset


def get_top_betweenness_complexes(df: DataFrame,
                                  num_complexes: int = 5) -> int:
    return 0


def get_graph_statistics(df: DataFrame,
                         graph: Graph,
                         num_proteins: int = 5,
                         num_complexes: int = 5) -> str:
    return ""
