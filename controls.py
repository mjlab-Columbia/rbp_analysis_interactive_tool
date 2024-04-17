from bokeh.models.widgets.inputs import Checkbox, TextInput
from bokeh.models.widgets import Slider, Select
import buttons

# Toggle for whether to display the labels on nodes
apply_labels_checkbox = Checkbox(active=False, label='LABELS')  # Checkbox

# Field for entering a specific protein name by user
protein_text_input = TextInput(title="Protein")  # TextInput

# Clustering options
graph_clustering_selection = Select(
    title="Clustering method",
    options=[
        'no clustering',
        'louvain',
        'markov',
        'clusterONE'
    ],
    value="no clustering"
)

# Clustering resolution
clustering_resolution = Slider(
    title="Clustering resolution",
    value=0.0,
    start=0.0,
    end=1.0,
    step=0.1
)

# Overlay options
node_coloring_selection = Select(
    title="Color node by",
    options=[
        'none',
        'lifecycle stage',
        'location',
        'disease'
    ],
    value="none"
)

num_neighbors_selection = Select(
    title="Number of Neighbors",
    options=[
        "1",
        "2",
        "3"
    ],
    value="1"
)

all_controls = [
    buttons.dataset_checkbox_button,
    buttons.dataset_combined_checkbox_button,
    buttons.interaction_type_checkbox_button,
    protein_text_input,
    num_neighbors_selection,
    graph_clustering_selection,
    clustering_resolution,
    node_coloring_selection,
    apply_labels_checkbox,
    buttons.create_interactome_button,
    buttons.download_button,
]
