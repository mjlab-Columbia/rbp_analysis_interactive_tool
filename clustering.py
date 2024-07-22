from networkx.algorithms.community import louvain_communities
import networkx as nx
from scipy.sparse import csr_matrix
import markov_clustering as mc
from os import system
from typing import List, Tuple, Dict, Union, Hashable, Any, Set
from itertools import combinations
from pdb import set_trace

# Mypy type aliases for brevity
NodeColoring = Dict[Hashable, Dict[str, str]]
Clusters = List[Tuple[str]]
ClusteredColoring = Tuple[Dict[Hashable, Dict[str, str]], Union[Clusters, None]]


def generate_clusters(graph: nx.Graph,
                      method: str,
                      res: float) -> Clusters:
    """
    Create node clusters based on a graph, clustering method, and resolution
    """
    print(f"Performing {method} clustering")

    if method == "louvain":
        res = 1 + (5 * res)
        raw_clusters = louvain_communities(graph,
                                           resolution=res)  # List[Set[str]]
        clusters = [tuple(cluster)
                    for cluster in raw_clusters]  # Clusters
    elif method == "markov":
        array = nx.to_scipy_sparse_array(graph)  # csr_array
        matrix = csr_matrix(array)  # csr_matrix

        # TODO: Figure out why these parameters were chosen
        inflation = 5.0 + (res * 1.5)
        result = mc.run_mcl(matrix, inflation=inflation)  # csc_matrix
        index_clusters = mc.get_clusters(result)  # List[Tuple[int]]

        # Output is by node index so we need to convert to node names
        node_map = dict(enumerate(graph.nodes()))  # Dict[int, str]
        clusters = [tuple(node_map[idx] for idx in c)
                    for c in index_clusters]  # Clusters
    elif method == "clusterONE":
        # Need to write edgelist in bytes for clusterONE
        input_file = "data.txt"  # str
        output_file = "clusters.txt"  # str
        with open(input_file, "wb") as f:
            nx.write_edgelist(graph, f, data=False)  # None

        base_command = "java -jar cluster_one-1.0.jar"  # str
        full_command = f"{base_command} {input_file} > {output_file}"  # str
        system(full_command)  # None

        # Since clusterONE is a CLI program, we need to parse the output
        with open("clusters.txt", "r") as f:
            lines = f.readlines()  # List[str]
            clusters = [tuple(line.split("\t")) for line in lines]  # Clusters
    else:
        message = "Supported clustering values: louvain, markov, clusterONE"
        raise Exception(message)

    return clusters


def get_edges_between_subgraphs(graph: nx.DiGraph,
                                group1: Tuple[str],
                                group2: Tuple[str]) -> List[Tuple[Any, Any]]:
    subgraph1 = graph.subgraph(group1)
    subgraph2 = graph.subgraph(group2)

    shared_edges = [(source, target) for source, target in graph.edges()
                    if source in subgraph1 and target in subgraph1
                    or source in subgraph2 and target in subgraph2]

    return shared_edges


def create_cluster_graph(graph: nx.DiGraph,
                         method: str,
                         res: float,
                         node_coloring: NodeColoring) -> nx.DiGraph:

    # Only cluster nodes for which we have color information
    colored_nodes: Set[str] = set(node_coloring.keys())
    graph_nodes: Set[str] = set(graph.nodes())
    nodes_to_cluster: Set[str] = colored_nodes.intersection(graph_nodes)
    colored_subgraph: nx.DiGraph = graph.subgraph(nodes_to_cluster)

    # Generate list of node clusters based on the method
    clusters: List[Tuple[str]] = generate_clusters(graph=colored_subgraph, res=res, method=method)
    cluster_indices: List[int] = [i for i in range(len(clusters))]
    cluster_pairs: List[Tuple[int, int]] = [combo for combo in combinations(cluster_indices, 2)]
    cluster_graph: nx.DiGraph = nx.DiGraph()

    # For every cluster pair in n Choose 2, determine the edge and its attributes
    for group1_idx, group2_idx in cluster_pairs:
        shared_edges = get_edges_between_subgraphs(graph,
                                                   clusters[group1_idx],
                                                   clusters[group2_idx])

        if len(shared_edges) != 0:
            cluster_graph.add_edge(group1_idx, group2_idx, weight=len(shared_edges))

            group1_coloring = [node_coloring[node]["node_color"] for node in clusters[group1_idx]]
            group2_coloring = [node_coloring[node]["node_color"] for node in clusters[group2_idx]]

            group1_majority_color = max(set(group1_coloring), key=group1_coloring.count)
            group2_majority_color = max(set(group2_coloring), key=group2_coloring.count)

            size_attr = {group1_idx: len(clusters[group1_idx]),
                         group2_idx: len(clusters[group2_idx])}

            color_attr = {group1_idx: group1_majority_color,
                          group2_idx: group2_majority_color}

            nx.set_node_attributes(cluster_graph, size_attr, "size")
            nx.set_node_attributes(cluster_graph, color_attr, "color")
        else:
            cluster_graph.add_node(group1_idx)
            cluster_graph.add_node(group2_idx)

    return cluster_graph
