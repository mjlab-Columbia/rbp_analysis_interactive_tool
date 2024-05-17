from os.path import dirname, join
from tqdm import tqdm

# TODO: Only import functions that are used from networkx
import networkx as nx
from pandas import DataFrame, read_csv

from clustering import get_clustered_coloring

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import Div, PreText, ColumnDataSource
from bokeh.models import Legend, LegendItem, LabelSet
from bokeh.models import CustomJS
from bokeh.plotting import figure

from buttons import create_interactome_button, download_button
from buttons import dataset_checkbox_button
import controls
from colors import EdgeColors, GRAY, PPI_SUPPORT, LifecycleColorsDict
from data import subset_by_protein, subset_by_edge_type, subset_by_node_type
from data import remap_df, determine_node_coloring, Dataset

from numpy import ndarray

from statistics import get_graph_statistics

from typing import List, Tuple, Dict

from pdb import set_trace

# Type aliases for readability
Edges = Tuple[List[List[float]], List[List[float]]]
EdgeStyles = List[List[int]]

# TODO: Turn type comments into mypy type annotations

# GLOBAL
graph_viewer = figure(
    height=600,
    width=800,
    tools=["pan", "wheel_zoom", "save", "reset"],
    active_scroll="wheel_zoom"
)

# Create Column Data Source that will be used by the plot
empty_edges_dict: Dict[str, List[str]] = dict(xs=[], ys=[],
                                              line_dash=[], line_color=[])
source_edges = ColumnDataSource(empty_edges_dict)
ppi_edges = ColumnDataSource(dict(xs=[], ys=[], line_dash=[]))
empty_nodes_dict: Dict[str, List] = dict(xs=[], ys=[], color=[], name=[],
                                         node_size=[], label=[])
source_nodes = ColumnDataSource(empty_nodes_dict)

graph_viewer.multi_line(
    'xs',
    'ys',
    line_color='line_color',
    line_dash='line_dash',
    source=source_edges,
    alpha=0.8,
    width=1
)

graph_viewer.multi_line(
    'xs',
    'ys',
    line_color=PPI_SUPPORT,
    line_dash='line_dash',
    source=ppi_edges,
    alpha=0.8,
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
        text='label',
        source=source_nodes,
        background_fill_color='white',
        text_font_size='8px',
        background_fill_alpha=.7
    )
)

# Create legend for edges
edge_legend_map = {
    "PPI Support": {"line_color": PPI_SUPPORT, "line_dash": "solid"},
    "Direct": {"line_color": GRAY, "line_dash": "solid"},
    "RNA Shielded": {"line_color": GRAY, "line_dash": "dashed"},
    "RNA Mediated": {"line_color": GRAY, "line_dash": "dotted"},
    "IP": {"line_color": EdgeColors.IP.value, "line_dash": "solid"},
    "SEC": {"line_color": EdgeColors.SEC.value, "line_dash": "solid"},
    "Both": {"line_color": EdgeColors.Both.value, "line_dash": "solid"},
}

edge_legend_items = []
for label, attr in edge_legend_map.items():
    renderer = graph_viewer.line(
        [0],
        [0],
        line_color=attr["line_color"],
        line_dash=attr["line_dash"]
    )
    legend_item = LegendItem(label=label, renderers=[renderer])
    edge_legend_items.append(legend_item)

edge_legend = Legend(
    items=edge_legend_items,
    location="right",
    orientation='vertical',
    label_text_font_size='6pt'
)

# Add static legend for edges to the right of the graph display
graph_viewer.add_layout(edge_legend, 'right')

# Create legend for nodes
node_legend_map = {
    "3' end processing": {"node_color": LifecycleColorsDict["3' end processing"]},
    "decay": {"node_color": LifecycleColorsDict["decay"]},
    "export": {"node_color": LifecycleColorsDict["export"]},
    "localization": {"node_color": LifecycleColorsDict["localization"]},
    "modification": {"node_color": LifecycleColorsDict["modification"]},
    "splicing": {"node_color": LifecycleColorsDict["splicing"]},
    "transcription": {"node_color": LifecycleColorsDict["transcription"]},
    "translation": {"node_color": LifecycleColorsDict["translation"]},
    "other": {"node_color": LifecycleColorsDict["undetermined"]},
}  # Dict[Hashable, Dict[str, str]]

node_legend_items = []  # List[GlyphRenderer]
for label, attr in node_legend_map.items():
    renderer = graph_viewer.scatter(
        [0],
        [0],
        [0],
        marker=['circle'],
        color=attr["node_color"]
    )
    legend_item = LegendItem(label=label, renderers=[renderer])
    node_legend_items.append(legend_item)

node_legend = Legend(
    items=node_legend_items,
    location="bottom_center",
    orientation="horizontal",
    label_text_font_size="6pt"
)

# Add static legend for node colors below the graph display
graph_viewer.add_layout(node_legend, "below")

# Other
statistics_info = PreText(text="Graph Statistics:", width=300)
global_df = read_csv("SupplementalTable_2.csv")


def subset_dataframe() -> DataFrame:
    """
    Returns portion of Supplemental Table 2 based on controls tagged
    """
    df = global_df.copy()
    remapped_df = remap_df(df)
    node_filtered_df = subset_by_node_type(remapped_df)
    edge_filtered_df = subset_by_edge_type(node_filtered_df)
    protein_filtered_df = subset_by_protein(edge_filtered_df)

    return protein_filtered_df


def get_edges(df: DataFrame,
              graph: nx.Graph,
              layout: Dict[str, ndarray],
              alpha: float) -> Tuple[Edges, Edges, List[str], EdgeStyles]:
    edge_x = []  # List[List[float64]]
    edge_y = []  # List[List[float64]]
    ppi_x = []  # List[List[float64]]
    ppi_y = []  # List[List[float64]]
    edge_colors = []  # List[str]
    edge_styles = []  # List[List[int]]

    progress_bar = tqdm(list(graph.edges()), total=len(graph.edges()))

    # If above statement didn't return anything, we render CORUM edges
    for edge in progress_bar:
        # Undirected edge is: (s_x, s_y) <--> (t_x, t_y)
        source = edge[0]
        target = edge[1]
        source_coordinates = layout[source]
        target_coordinates = layout[target]
        s_x, s_y = source_coordinates
        t_x, t_y = target_coordinates

        # Subset to bait-prey pair in CORUM (in_CORUM_2022 is a boolean Series)
        bait_mask = df["Bait"] == source
        prey_mask = df["Prey"] == target
        corum_mask = df["in_CORUM_2022"]
        df_edge = df[bait_mask & prey_mask & corum_mask]

        if df_edge.shape[0] == 1:
            # Find x and y components of the slope
            slope_x = t_x - s_x
            slope_y = t_y - s_y

            # Translate original edge coordinates to make room for PPI edge
            new_s_x = s_x + (alpha * slope_y)
            new_s_y = s_y - (alpha * slope_x)
            new_t_x = t_x + (alpha * slope_y)
            new_t_y = t_y - (alpha * slope_x)

            # Add coordinates, colors, and styles to lists
            edge_x.append([new_s_x, new_t_x])
            edge_y.append([new_s_y, new_t_y])
            edge_color = df_edge["edge_color"].values[0]
            edge_style = df_edge["edge_style"].values[0]
            edge_colors.append(edge_color)
            edge_styles.append(edge_style)

            # Create source (s) and target (t) coordinates for PPI edge
            ppi_s_x = s_x - (alpha * slope_y)
            ppi_s_y = s_y + (alpha * slope_x)
            ppi_t_x = t_x - (alpha * slope_y)
            ppi_t_y = t_y + (alpha * slope_x)

            # Add coordinates to list of ppi edges
            ppi_x.append([ppi_s_x, ppi_t_x])
            ppi_y.append([ppi_s_y, ppi_t_y])
        elif df_edge.shape[0] == 0:
            # No PPI support --> don't translate coordinates
            edge_x.append([s_x, t_x])
            edge_y.append([s_y, t_y])
            edge_color = df[bait_mask & prey_mask]["edge_color"].values[0]
            edge_style = df[bait_mask & prey_mask]["edge_style"].values[0]
            edge_colors.append(edge_color)
            edge_styles.append(edge_style)
        else:
            message = "Each edge should correspond to at most 1 entry"
            raise Exception(message)

    return (edge_x, edge_y), (ppi_x, ppi_y), edge_colors, edge_styles


def update() -> None:
    df = subset_dataframe()  # DataFrame
    graph = nx.from_pandas_edgelist(df,
                                    source="Bait",
                                    target="Prey",
                                    edge_attr=["edge_color", "edge_style"],
                                    create_using=nx.DiGraph)

    # TODO: Relocate this to initial ColumnDataSource declaration
    node_sizes = [12 for _ in range(len(graph.nodes()))]  # List[int]

    # Layout is a mapping of nodes --> coordinates
    layout = nx.spring_layout(graph)  # Dict[str, ndarray]

    # Show protein names only if the LABELS checkbox is active
    if controls.apply_labels_checkbox.active:
        new_node_names = list(layout.keys())  # List[str]
    else:
        new_node_names = ["" for _ in layout.keys()]

    node_coordinates = zip(*layout.values())
    node_x, node_y = node_coordinates

    datasets_toggled = dataset_checkbox_button.active
    if Dataset.CORUM.value in datasets_toggled:
        # TODO: Find way to dynamically set alpha based on node_size
        alpha = 5e-5
        new_edges_out = get_edges(df, graph, layout, alpha)
        new_edges, new_ppi_edges, new_edge_colors, new_edge_styles = new_edges_out
        edge_x, edge_y = new_edges
        ppi_x, ppi_y = new_ppi_edges
    else:
        edge_x = [[layout[source][0], layout[target][0]]
                  for source, target in graph.edges()]
        edge_y = [[layout[source][1], layout[target][1]]
                  for source, target in graph.edges()]

        edge_data = list(graph.edges(data=True))
        new_edge_colors = [attr["edge_color"] for _, _, attr in edge_data]
        new_edge_styles = [attr["edge_style"] for _, _, attr in edge_data]

        ppi_x, ppi_y = [], []

    # Account for the "color by" selection menu
    clustering_method = str(controls.graph_clustering_selection.value)  # str
    resolution = float(controls.clustering_resolution.value)  # float
    init_node_coloring = determine_node_coloring(df)
    clustering_output = get_clustered_coloring(graph,
                                               method=clustering_method,
                                               res=resolution,
                                               node_coloring=init_node_coloring)
    new_node_coloring, clusters = clustering_output

    # Get new summary statistics based on new subsetting of the dataframe/graph
    new_stats = get_graph_statistics(df,
                                     graph,
                                     num_proteins=5,
                                     num_complexes=5,
                                     clusters=clusters)  # str
    statistics_info.text = new_stats

    new_node_colors = []  # List[str]
    for node in layout.keys():
        if node in new_node_coloring:
            color = new_node_coloring[node]["node_color"]
            new_node_colors.append(color)
        else:
            new_node_colors.append(GRAY)

    source_nodes.data = dict(
        xs=node_x,
        ys=node_y,
        names=new_node_names,
        color=new_node_colors,
        node_size=node_sizes,
        label=new_node_names,
    )

    source_edges.data = dict(
        xs=edge_x,
        ys=edge_y,
        line_dash=new_edge_styles,
        line_color=new_edge_colors,
    )

    ppi_edges.data = dict(
        xs=ppi_x,
        ys=ppi_y,

        # TODO: Replace with "solid" in ColumnDataSource declaration
        line_dash=[[] for _ in ppi_x]
    )

    # Set callback on each update to ensure only subset_df is downloaded
    download_button.js_on_click(
        CustomJS(
            args=dict(source=subset_dataframe().to_csv(index=False)),
            code=open(join(dirname(__file__), "download.js")).read()
        )
    )


# Arranging final output
primary_div = Div(
    text=open(join(dirname(__file__), "Interactome.html")).read(),
    sizing_mode="stretch_width"
)

create_interactome_button.on_click(update)
control_inputs = column(*controls.all_controls, width=300)
statistics_output = column(*[statistics_info], width=300)

root_column = column(
    primary_div,
    row(control_inputs, graph_viewer, statistics_output),
    sizing_mode="scale_both"
)

# curdoc == current document
curdoc().add_root(root_column)
curdoc().title = "Interactome"
