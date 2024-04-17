from os.path import dirname, join
from os import system

import networkx as nx
from networkx.algorithms.community import louvain_communities
from pandas import DataFrame, read_csv

from bokeh.io import curdoc
from bokeh.models import Div, PreText, ColumnDataSource
from bokeh.models import Legend, LegendItem, LabelSet
from bokeh.models import CustomJS
from bokeh.layouts import column, row
from bokeh.plotting import figure

import controls
from buttons import create_interactome_button, download_button
from colors import NodeColors, EdgeColors
from data import subset_by_protein, subset_by_edge_type, subset_by_node_type
import markov_clustering as mc

# GLOBAL
graph_viewer = figure(
    height=600,
    width=800,
    tools=[
        "pan",
        "wheel_zoom",
        "save",
        "reset"
    ],
    active_scroll="wheel_zoom"
)

# Create Column Data Source that will be used by the plot
empty_edges_dict = dict(xs=[], ys=[], style=[], color=[],
                        label=[])  # Dict[str, List]
source_edges = ColumnDataSource(empty_edges_dict)
empty_nodes_dict = dict(xs=[], ys=[], color=[],
                        name=[], node_size=[], label=[])  # Dict[str, List]
source_nodes = ColumnDataSource(empty_nodes_dict)

graph_viewer.multi_line(
    'xs',
    'ys',
    line_color='color',
    line_dash='style',
    source=source_edges,
    alpha=0.5,
    width=1
)

graph_viewer.scatter(
    'xs',
    'ys',
    fill_color='color',
    line_color='black',
    size='node_size',
    source=source_nodes
)

graph_viewer.add_layout(
    LabelSet(
        x='xs',
        y='ys',
        text='name',
        source=source_nodes,
        background_fill_color='white',
        text_font_size='8px',
        background_fill_alpha=.7
    )
)

# Create legend for edges
edge_type_map = {
    "PPI Support": {
        "line_color": EdgeColors.ppi_support.value,
        "line_dash": 'solid'
    },
    "Direct": {
        "line_color": EdgeColors.direct.value,
        "line_dash": 'solid'
    },
    "RNA Shielded": {
        "line_color": EdgeColors.rna_shielded.value,
        "line_dash": 'dashed'
    },
    "RNA Mediated": {
        "line_color": EdgeColors.rna_mediated.value,
        "line_dash": 'dotted'
    }
}

# Define node type colormap
node_type_colormap = {
    "IP": NodeColors.IP.value,
    "SEC": NodeColors.SEC.value,
    "Both": NodeColors.Both.value,
}

node_type_map = {
    "IP": {
        "line_color": NodeColors.IP.value,
        "line_dash": "solid"
    },
    "SEC": {
        "line_color": NodeColors.SEC.value,
        "line_dash": "solid"
    },
    "Both": {
        "line_color": NodeColors.Both.value,
        "line_dash": "solid"
    },
}

# TODO: Try to resolve (BAD_COLUMN_NAME) error without ruining functionality
legend_items = []
for label, attr in edge_type_map.items():
    renderer = graph_viewer.line(
        [0],
        [0],
        line_color=attr["line_color"],
        line_dash=attr["line_dash"]
    )
    legend_item = LegendItem(label=label, renderers=[renderer])
    legend_items.append(legend_item)

for label, attr in node_type_map.items():
    renderer = graph_viewer.line(
        [0],
        [0],
        line_color=attr["line_color"],
        line_dash=attr["line_dash"]
    )
    legend_item = LegendItem(label=label, renderers=[renderer])
    legend_items.append(legend_item)

legend = Legend(
    items=legend_items,
    location="right",
    orientation='vertical',
    label_text_font_size='6pt'
)

# Add static legends around the graph
graph_viewer.add_layout(legend, 'right')

# Other
statistics = PreText(text="Graph Statistics:", width=300)
global_df = read_csv("SupplementalTable_2.csv")
subset_df = global_df.copy()


def subset_dataframe() -> DataFrame:
    """
    Returns portion of Supplemental Table 2 based on controls tagged
    """
    df = global_df.copy()
    node_filtered_df = subset_by_node_type(df)
    edge_filtered_df = subset_by_edge_type(node_filtered_df)
    protein_filtered_df = subset_by_protein(edge_filtered_df)

    print("Original", df.head())
    print("Node filtered", node_filtered_df.head())
    print("Edge filtered", edge_filtered_df.head())
    print("Protein filtered", protein_filtered_df.head())
    return protein_filtered_df


def cluster_graph(graph: nx.Graph, method: str) -> nx.Graph:
    """
    Clusters nodes in an undirected graph
    """
    res = controls.clustering_resolution.value

    if method == "louvain":
        clusters = louvain_communities(graph, resolution=res)  # List[Set[str]]
    elif method == "markov":
        matrix = nx.to_scipy_sparse_matrix(graph)
        result = mc.run_mcl(matrix, inflation=res)
        clusters = mc.get_clusters(result)  # List[List[str]]
    elif method == "clusterONE":
        nx.write_edgelist(graph, "tmp.txt", data=False)  # None
        command = 'java -jar cluster_one-1.0.jar tmp.txt > clusters.txt'  # str
        system(command)

        with open("clusters.txt", "r") as f:
            clusters = f.readlines()
    else:
        message = "Supported values: louvain, markov, clusterONE"
        raise Exception(message)

    return clusters


def update() -> None:
    df = subset_dataframe()  # DataFrame
    subset_df = df.copy()  # DataFrame
    graph = nx.from_pandas_edgelist(df, source="Bait", target="Prey",
                                    edge_attr=["edge_color", "edge_style"])
    node_sizes = [8 for _ in range(len(graph.nodes()))]  # List[int]
    new_edges = list(graph.edges())  # List[Tuple[str, str]]

    # Layout is a mapping of nodes --> coordinates
    # Dict[str, np.ndarray]
    layout = nx.spring_layout(graph)

    new_names = list(layout.keys())  # List[str]
    # List[np.ndarray]
    coordinates = list(layout.values())

    # Tuple[Tuple[np.float64]]
    node_x, node_y = zip(*coordinates)

    # For each new edge, get the respective x coordinates for source -> target
    edge_x = [
        [layout[edge[0]][0], layout[edge[1]][0]]
        for edge in new_edges
    ]  # List[List[float, float]]

    # For each new edge, get the respective x coordinates for source -> target
    edge_y = [
        [layout[edge[0]][1], layout[edge[1]][1]]
        for edge in new_edges
    ]  # List[List[float, float]]

    # TODO: Implement color by logic here
    new_node_colors = ["black" for _ in range(len(graph.nodes()))]

    # TODO: Implement interaction type coloring/styling
    new_edge_colors = ["black" for _ in range(len(graph.nodes()))]
    new_edge_styles = ["solid" for _ in range(len(graph.nodes()))]

    source_nodes.data = dict(
        xs=node_x,
        ys=node_y,
        names=new_names,
        color=new_node_colors,
        node_size=node_sizes,
        label=new_names
    )

    source_edges.data = dict(
        xs=edge_x,
        ys=edge_y,
        style=new_edge_styles,
        color=new_edge_colors,
        label=new_names
    )


# Arranging final output
primary_div = Div(
    text=open(
        join(dirname(__file__), "Interactome.html")
    ).read(),
    sizing_mode="stretch_width"
)

# Triggers the creation of the interactome graph
create_interactome_button.on_click(update)

# Download button which triggers Supplemental Table 2 subset download
download_button.js_on_click(
    CustomJS(
        args=dict(
            source=subset_df.to_csv(index=False)
        ),
        code=open(join(dirname(__file__), "download.js")).read()
    )
)

# User-defined toggles, parameters, cutoffs, etc.
control_inputs = column(*controls.all_controls, width=300)

# Displaying graph statistics
statistics_output = column(*[statistics], width=300)

# Parent column where everything is laid out
root_column = column(
    primary_div,
    row(control_inputs, graph_viewer, statistics_output),
    sizing_mode="scale_both"
)

# curdoc == current document
curdoc().add_root(root_column)
curdoc().title = "Interactome"
