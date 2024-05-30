from networkx.algorithms.community import louvain_communities
import networkx as nx
from scipy.sparse import csr_matrix
import markov_clustering as mc
from os import system
from typing import List, Tuple, Dict, Union, Hashable
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


def get_clustered_coloring(graph: nx.Graph,
                           method: str,
                           res: float,
                           node_coloring: NodeColoring) -> ClusteredColoring:
    """
    Clusters graph based on one of several methods and returns
    coloring dictionary with majority color within each cluster
    """

    # Only recolor nodes if we're clustering
    if method != "no clustering":
        clusters = generate_clusters(graph, method, res)  # List[Tuple[str]]
        subgraph_colorings = []  # List[Dict[str, Dict[str, str]]]

        for cluster in clusters:
            subgraph = graph.subgraph(cluster)  # nx.Graph
            nodes = set(subgraph.nodes())  # Set[str]
            node_colors = {k: v["node_color"]
                           for k, v in node_coloring.items()
                           if k in nodes}

            # Find the color that appears the most often in the cluster
            all_colors = list(node_colors.values())  # List[str]
            majority_color = max(set(all_colors),
                                 key=all_colors.count)  # str
            majority_coloring = {k: {"node_color": majority_color}
                                 for k in subgraph.nodes()}  # Dict[str, str]
            subgraph_colorings.append(majority_coloring)

        # TODO: Resolve the "Need type annotation error" for new_node_colors
        new_node_colors: Dict[Hashable, Dict[str, str]] = {}
        for coloring in subgraph_colorings:
            new_node_colors |= coloring

        return new_node_colors, clusters
    else:
        # TODO: Visualize the clusters even if "Color nodes by" is none
        # New node colors = old node colors if we're not clustering
        new_node_colors = {k: {"node_color": v["node_color"]}
                           for k, v in node_coloring.items()}
        return new_node_colors, None
