from pandas import DataFrame
from data import InteractionTypeValue
from networkx import betweenness_centrality, Graph

from typing import List, Tuple


def get_number_of_nodes(df: DataFrame) -> int:
    baits = set(df["Bait"].unique())  # Set[str]
    preys = set(df["Prey"].unique())  # Set[str]

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


def get_number_of_mediated_interactions(df: DataFrame) -> int:
    mediated = InteractionTypeValue.mediated.value  # str
    mediated_df = df[df["IP_interaction_type"] == mediated]  # DataFrame
    mediated_interactions = mediated_df.shape[0]  # int
    return mediated_interactions


def get_number_of_shielded_interactions(df: DataFrame) -> int:
    shielded = InteractionTypeValue.shielded.value  # str
    shielded_df = df[df["IP_interaction_type"] == shielded]  # DataFrame
    shielded_interactions = shielded_df.shape[0]  # int
    return shielded_interactions


def get_number_of_shared_interactions(df: DataFrame) -> int:
    col = "Interaction_support"  # str
    subset_df = df[df[col] == "both"]  # DataFrame
    return subset_df.shape[0]


def get_number_of_ip_interactions(df: DataFrame) -> int:
    col = "Interaction_support"  # str
    subset_df = df[df[col] == "IP"]  # DataFrame
    return subset_df.shape[0]


def get_top_betweenness_proteins(df: DataFrame,
                                 graph: Graph,
                                 num_proteins: int) -> List[Tuple[str, float]]:
    bw_dict = betweenness_centrality(graph)  # Dict[str, float]
    sorted_bw = sorted(bw_dict.items(),
                       key=lambda x: x[1])  # List[Tuple[str, float]]
    sorted_bw_subset = sorted_bw[:num_proteins]  # List[Tuple[str, float]]
    return sorted_bw_subset


# TODO: Implement if needed (ask Lena)
def get_top_betweenness_complexes(df: DataFrame,
                                  num_complexes: int = 5) -> int:
    return 0


def get_graph_statistics(df: DataFrame,
                         graph: Graph,
                         num_proteins: int = 5,
                         num_complexes: int = 5) -> str:
    """
    Calculates and formats summary statistics for the interactome graph
    """
    # Get all of the statistics we care about
    num_nodes = get_number_of_nodes(df)  # int
    num_interactions = get_number_of_interactions(df)  # int

    direct = get_number_of_direct_interactions(df)  # int
    rna_mediated = get_number_of_mediated_interactions(df)  # int
    rna_shielded = get_number_of_shielded_interactions(df)  # int

    shared_interactions = get_number_of_shared_interactions(df)  # int
    ip_interactions = get_number_of_ip_interactions(df)  # int

    bw_proteins = get_top_betweenness_proteins(df,
                                               graph,
                                               num_proteins=5)

    top_proteins = [protein for protein, _ in bw_proteins]  # List[str]

    # Format all the statistics into single string to print
    stats_string = f"""
    Number of Nodes: {num_nodes}
    Number of Interactions: {num_interactions}

    Direct: {direct}
    RNA Mediated: {rna_mediated}
    RNA Shielded: {rna_shielded}

    Shared Interactions: {shared_interactions}
    IP Interactions: {ip_interactions}

    Top Five Betweeness Proteins:
    {top_proteins[0]}
    {top_proteins[1]}
    {top_proteins[2]}
    {top_proteins[3]}
    {top_proteins[4]}
    """

    return stats_string
