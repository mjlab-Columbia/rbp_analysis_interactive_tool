from pandas import DataFrame, read_excel
from data import InteractionTypeValue, Dataset
from networkx import betweenness_centrality, Graph
from pdb import set_trace

from typing import List, Tuple, Union

from buttons import dataset_checkbox_button

Clusters = List[Tuple[str]]


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
    subset_df = df[df[col] == "both"].copy()  # DataFrame
    return subset_df.shape[0]


def get_number_of_sec_interactions(df: DataFrame) -> int:
    col = "Interaction_support"  # str
    subset_df = df[df[col] == "SEC"].copy()  # DataFrame
    return subset_df.shape[0]


def get_number_of_ip_interactions(df: DataFrame) -> int:
    col = "Interaction_support"  # str
    subset_df = df[df[col] == "IP"].copy()  # DataFrame
    return subset_df.shape[0]


def get_top_betweenness_proteins(df: DataFrame,
                                 graph: Graph,
                                 num_proteins: int) -> List[Tuple[str, float]]:
    bw_dict = betweenness_centrality(graph)  # Dict[str, float]
    sorted_bw_asc = sorted(bw_dict.items(),
                           key=lambda x: x[1])  # List[Tuple[str, float]]
    sorted_bw_desc = list(reversed(sorted_bw_asc))  # List[Tuple[str, float]]
    sorted_bw_subset = sorted_bw_desc[:num_proteins]
    return sorted_bw_subset


def get_num_corum_complexes(df: DataFrame) -> int:
    unique_complexes = df["CORUM_complex_2022"].dropna().unique()  # ndarray
    return len(unique_complexes)


def get_top_betweenness_complexes(df: DataFrame,
                                  num_complexes: int = 5) -> int:

    return 0


def get_graph_statistics(df: DataFrame,
                         graph: Graph,
                         num_proteins: int,
                         num_complexes: int,
                         clusters: Union[Clusters, None]) -> str:
    """
    Calculates and formats summary statistics for the interactome graph
    """
    active_datasets = dataset_checkbox_button.active

    # Get all of the statistics we care about
    num_nodes = get_number_of_nodes(df)  # int
    num_interactions = get_number_of_interactions(df)  # int

    direct = get_number_of_direct_interactions(df)  # int
    rna_mediated = get_number_of_mediated_interactions(df)  # int
    rna_shielded = get_number_of_shielded_interactions(df)  # int

    shared_interactions = get_number_of_shared_interactions(df)  # int
    sec_interactions = get_number_of_sec_interactions(df)  # int
    ip_interactions = get_number_of_ip_interactions(df)  # int

    bw_proteins = get_top_betweenness_proteins(df=df,
                                               graph=graph,
                                               num_proteins=5)

    top_proteins = [protein for protein, _ in bw_proteins]  # List[str]

    # Format all the statistics into single string to print
    num_nodes_stats = f"Number of Nodes: {num_nodes}"
    num_interactions_stats = f"Number of Interactions: {num_interactions}"
    node_stats = "\n".join([num_nodes_stats, num_interactions_stats])

    direct_stats = f"Direct: {direct}"  # str
    mediated_stats = f"RNA Mediated: {rna_mediated}"  # str
    shielded_stats = f"RNA Shieleded: {rna_shielded}"  # str
    interaction_type_stats = "\n".join([direct_stats,
                                        mediated_stats,
                                        shielded_stats])  # str

    # TODO: Ask Lena about shortening this by doing SEC Interactions: N/A
    # Create each string separately to avoid extra indentations
    if Dataset.SEC.value in active_datasets:
        shared = f"Shared Interactions: {shared_interactions}"
        ip = f"IP Interactions: {ip_interactions}"
        sec = f"SEC Interactions: {sec_interactions}"
        interaction_origin_stats = "\n".join([shared, ip, sec])
    else:
        shared = f"Shared Interactions: {shared_interactions}"
        ip = f"IP Interactions: {ip_interactions}"
        interaction_origin_stats = "\n".join([shared, ip])

    preamble = "Top Five Betweeness Proteins:"
    protein_list = "\n".join(top_proteins[:5])
    top_five_proteins = f"{preamble}\n{protein_list}"

    if clusters is None:
        stats_list = [node_stats, interaction_type_stats,
                      interaction_origin_stats, top_five_proteins]  # List[str]
        stats_string = "\n\n".join(stats_list)  # str
    else:
        num_complexes = get_num_corum_complexes(df)  # int

        norm_bw_filepath = "SupplementalTable_8.xlsx"  # str
        norm_bw_df = read_excel(norm_bw_filepath,
                                sheet_name=1,
                                index_col=0)  # DataFrame
        norm_bw_df.reset_index(inplace=True)  # None
        norm_bw_df.columns = ["protein",
                              "bw_centrality", "is_bait"]  # List[str]

        # Populate dict with most common CORUM complex per cluster and median
        # betweenness centrality of proteins in cluster
        complex_bw = {}  # Dict[str, float]
        for cluster in clusters:
            bait_mask = df["Bait"].isin(cluster)  # Series
            prey_mask = df["Prey"].isin(cluster)  # Series
            df_subset = df[bait_mask | prey_mask].copy()  # DataFrame

            # Find the most common CORUM complex in the cluster
            corum_complexes = df_subset["CORUM_complex_2022"]  # Series
            corum_complexes_no_nan = corum_complexes.dropna()
            complex_counts = corum_complexes_no_nan.value_counts()  # Series

            # Only add items to complex_bw if there's >= 1 complex
            if complex_counts.shape[0] != 0:
                most_common_complex = complex_counts.idxmax()  # str

                # Find the median betweenness centrality of proteins in cluster
                proteins = df_subset[["Bait", "Prey"]].values.ravel()
                protein_mask = norm_bw_df["protein"].isin(proteins)  # Series
                median_bw = norm_bw_df[protein_mask]["bw_centrality"].median()
                complex_bw[most_common_complex] = median_bw

        # Sort complexes by betweenness centrality
        sorted_complex_bw = sorted(complex_bw.items(),
                                   key=lambda x: x[1],
                                   reverse=True)  # List[Tuple[str, str]]
        top_complexes = [name for name, _ in sorted_complex_bw]  # List[str]
        num_complexes_stats = f"Number of Complexes: {num_complexes}"
        preamble = "Top Five Betweenness Complexes:"
        complex_list = "\n".join(top_complexes[:5])
        top_five_complexes = "\n".join([preamble, complex_list])

        stats_list = [node_stats, interaction_type_stats,
                      interaction_origin_stats, top_five_proteins,
                      num_complexes_stats, top_five_complexes]
        stats_string = "\n\n".join(stats_list)

    return stats_string
