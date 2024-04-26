# from bokeh.core.properties import DashPattern
from enum import Enum
from pandas import DataFrame, concat
from matplotlib.pyplot import cm
from numpy import nan

from typing import Dict, Hashable, List

from buttons import dataset_combined_checkbox_button, dataset_checkbox_button
from buttons import interaction_type_checkbox_button
from colors import EdgeColors, LifecycleColorsDict, LifecycleColors
import controls


class Dataset(Enum):
    IP: int = 0
    SEC: int = 1
    CORUM: int = 2


class CombinedDataset(Enum):
    IP_and_SEC: int = 0


class InteractionTypeCheckbox(Enum):
    direct: int = 0
    mediated: int = 1
    shielded: int = 2
    undetermined: int = 3


class InteractionTypeValue(Enum):
    direct: str = "direct"
    mediated: str = "mediated"
    shielded: str = "shielded"
    undetermined: str = "undetermined"


def remap_df(df: DataFrame) -> DataFrame:
    new_bait_lifecycle_dict = {
        nan: nan,
        "decay": "decay",
        "export": "export",
        "splicing": "splicing",
        "new": "undetermined",
        "localization": "localization",
        "3' end processing": "3' end processing",
        "translation": "translation",
        "new": "undetermined",
        "transcription": "transcription",
        "negative control": "negative control",
        "modification": "modification"
    }

    new_prey_lifecycle_dict = {
        "undetermined": "undetermined",
        "splicing": "splicing",
        "degradation": "decay",
        "translation": "translation",
        "transcription": "transcription",
        "export": "export",
        "3' end processing": "3' end processing",
        "localization": "localization",
        "modification": "modification"
    }

    bait_col = "Bait_lifecycle_step"
    prey_col = "prey_lifecycle_stage_main_by_most_common_bait_stage"
    df_copy = df.copy()
    df_copy[bait_col] = df_copy[bait_col].map(new_bait_lifecycle_dict)
    df_copy[prey_col] = df_copy[prey_col].map(new_prey_lifecycle_dict)
    return df_copy


def rgba_to_hex(r, g, b, a) -> str:
    r_int = int(r * 255)
    g_int = int(g * 255)
    b_int = int(b * 255)
    a_int = int(a * 255)
    hex_color = "#{:02x}{:02x}{:02x}{:02x}".format(r_int, g_int, b_int, a_int)
    return hex_color


def subset_by_node_type(df: DataFrame) -> DataFrame:
    """
    Return the baits corresponding to the dataset the user selected via
    `dataset_checkbox_button` and `dataset_combined_checkbox_button`
    """

    # Each type of interaction has a color
    df["edge_color"] = df["IP_interaction_type"].map({
        "IP": EdgeColors.IP.value,
        "SEC": EdgeColors.SEC.value,
        "both": EdgeColors.Both.value
    })
    datasets = dataset_checkbox_button.active  # List[int]
    # List[int]
    ip_and_sec_toggle = dataset_combined_checkbox_button.active

    # Group combinations of toggles
    ip_bool = Dataset.IP.value in datasets  # bool
    sec_bool = Dataset.SEC.value in datasets  # bool

    ip_only = (ip_bool) and (not sec_bool)  # bool
    sec_only = (sec_bool) and (not ip_bool)  # bool
    ip_or_sec = (ip_bool) and (sec_bool)  # bool
    ip_and_sec = CombinedDataset.IP_and_SEC.value in ip_and_sec_toggle  # bool

    # A mask tells us which rows of the dataframe (df) to return
    ip_mask = df["Interaction_support"] == "IP"  # Series
    sec_mask = df["Interaction_support"] == "SEC"  # Series
    both_mask = df["Interaction_support"] == "both"  # Series

    # Depending on the toggle, we return a different union of masks
    if ip_only:
        return df[ip_mask | both_mask]
    elif sec_only:
        return df[sec_mask | both_mask]
    elif ip_or_sec:
        return df[ip_mask | sec_mask | both_mask]
    elif ip_and_sec:
        return df[both_mask]
    else:
        return df


def subset_by_edge_type(df: DataFrame) -> DataFrame:
    """
    Subset df to bait-prey pairs based on protein interaction type
    """
    df_copy = df.copy()

    # Note: DashPatterns take iterable as [dash, gap] ([] is a solid line)
    df_copy["edge_style"] = df_copy["IP_interaction_type"].map({
        InteractionTypeValue.direct.value: [],
        InteractionTypeValue.mediated.value: [8, 4],
        InteractionTypeValue.shielded.value: [5, 20],
        InteractionTypeValue.undetermined.value: [1, 1],
        nan: []
    })

    df_copy["edge_color"] = df_copy["Interaction_support"].map({
        "IP": EdgeColors.IP.value,
        "SEC": EdgeColors.SEC.value,
        "both": EdgeColors.Both.value
    })

    interaction_type = interaction_type_checkbox_button.active  # List[int]
    types_to_include = []

    if len(interaction_type) == 0:
        return df

    # TODO: This logic can be simpler, so figure out a simplification
    if InteractionTypeCheckbox.direct.value in interaction_type:
        types_to_include.append(InteractionTypeValue.direct.value)

    if InteractionTypeCheckbox.mediated.value in interaction_type:
        types_to_include.append(InteractionTypeValue.mediated.value)

    if InteractionTypeCheckbox.shielded.value in interaction_type:
        types_to_include.append(InteractionTypeValue.shielded.value)

    # Only consider edges with IP or IP & SEC for IP interaction type filtering
    ip_mask = df_copy["Interaction_support"].isin(["IP", "both"])  # Series
    ip_type_mask = df_copy[ip_mask]["IP_interaction_type"].isin(
        types_to_include)
    ip_rows = df_copy[ip_mask][ip_type_mask].copy()  # DataFrame

    # Get SEC rows separately
    sec_mask = df_copy["Interaction_support"] == "SEC"  # Series
    sec_rows = df_copy[sec_mask].copy()  # DataFrame

    filtered_df = concat([ip_rows, sec_rows], axis=0)  # DataFrame
    return filtered_df


def filter_by_dfs(df: DataFrame,
                  nodes: List[str],
                  n: int) -> DataFrame:
    """
    Performs depth first search on bait-prey dataframe
    """
    indices = []
    visited = set()

    def dfs(node, depth):
        visited.add(node)
        if depth <= n:
            neighbors = df[df['Bait'] == node].index
            for neighbor_idx in neighbors:
                indices.append(neighbor_idx)
                neighbor = df.loc[neighbor_idx, 'Prey']
                if neighbor not in visited:
                    dfs(neighbor, depth + 1)

    for node in nodes:
        dfs(node, 0)

    filtered_df = df.loc[indices]
    return filtered_df


def subset_by_protein(df: DataFrame) -> DataFrame:
    """
    Subset df to bait proteins list specified in protein_text_input
    """
    protein_str = str(controls.protein_text_input.value).strip()  # str
    protein_list = [p.strip() for p in protein_str.split(",")]  # List[str]

    # If there's no proteins entered, we return all edges
    if (len(protein_list) == 1) and (protein_list[0] == ""):
        return df
    else:
        num_neighbors_value = controls.num_neighbors_selection.value  # str
        num_neighbors = int(num_neighbors_value) - 1
        filtered_df = filter_by_dfs(df, protein_list, num_neighbors)
        return filtered_df


def lifecycle_node_attributes(df: DataFrame) -> Dict[Hashable, Dict[str, str]]:
    bait_col = "Bait_lifecycle_step"  # str
    prey_col = "prey_lifecycle_stage_main_by_most_common_bait_stage"  # str

    # Colors for bait proteins
    bait_lifecycle_df = df[["Bait", bait_col]].copy()  # DataFrame
    bait_lifecycle_df.dropna(inplace=True)  # None
    bait_lifecycle_df.drop_duplicates(keep="first", inplace=True)  # None
    bait_lifecycle_df.set_index("Bait", inplace=True)  # DataFrame

    # TODO: Fix autoformatter putting type hints on new line when >80 chars
    # Dict[str, Dict[str, str]]
    bait_lifecycle_dict = bait_lifecycle_df.to_dict()
    bait_lifecycle_dict = bait_lifecycle_dict[bait_col]  # Dict[str, str]

    # Colors for prey proteins
    prey_lifecycle_df = df[["Prey", prey_col]].copy()
    prey_lifecycle_df.dropna(inplace=True)
    prey_lifecycle_df.drop_duplicates(keep="first", inplace=True)
    prey_lifecycle_df.set_index("Prey", inplace=True)

    # Dict[str, Dict[str, str]]
    prey_lifecycle_dict = prey_lifecycle_df.to_dict()
    prey_lifecycle_dict = prey_lifecycle_dict[prey_col]  # Dict[str, str]

    # The | operation merges two dictionaries by keys
    protein_dict = bait_lifecycle_dict | prey_lifecycle_dict  # Dict[str, str]

    node_attributes = {}  # Dict[Hashable, Dict[str, str]]
    for protein, lifecycle_step in protein_dict.items():
        if lifecycle_step in LifecycleColorsDict:
            node_color = LifecycleColorsDict[lifecycle_step]
        else:
            node_color = LifecycleColors.undetermined.value

        attributes = {
            "legend_label": lifecycle_step,
            "node_color": node_color
        }  # Dict[str, str]
        node_attributes[protein] = attributes

    return node_attributes


def location_node_attributes(df: DataFrame) -> Dict[Hashable, Dict[str, str]]:
    bait_col = "Bait_main_location_HPA"
    prey_col = "Prey_main_location_HPA"

    # Get location annotations for baits
    bait_location_df = df[["Bait", bait_col]].copy()
    bait_location_df.dropna(inplace=True)
    bait_location_df.drop_duplicates(keep="first", inplace=True)
    bait_location_df.set_index("Bait", inplace=True)

    bait_location_dict = bait_location_df.to_dict()
    bait_location_dict = bait_location_dict[bait_col]
    unique_bait_locations = set(bait_location_dict.values())

    # Get location annotations for preys
    prey_location_df = df[["Prey", prey_col]].copy()  # DataFrame
    prey_location_df.dropna(inplace=True)
    prey_location_df.drop_duplicates(keep="first", inplace=True)
    prey_location_df.set_index("Prey", inplace=True)

    # Dict[str, Dict[str, str]]
    prey_location_dict = prey_location_df.to_dict()
    prey_location_dict = prey_location_dict[prey_col]  # Dict[str, str]
    unique_prey_locations = set(prey_location_dict.values())  # Set[str]

    # The | operator merges dictionaries based on keys
    # Dict[str, str]
    combined_location_dict = bait_location_dict | prey_location_dict

    unique_locations = list(unique_bait_locations |
                            unique_prey_locations)  # List[str]

    # integer --> color in RGBA
    color_map = cm.get_cmap("tab10", len(unique_locations))  # Colormap

    # location --> color
    color_dict = {location: color_map(i)
                  for i, location in enumerate(unique_locations)}

    # TODO: Make a legend for the disease coloring

    # node --> {"node_attribute": color}
    node_attributes = {
        key: {"node_color": rgba_to_hex(*color_dict[value])}
        for key, value in combined_location_dict.items()
    }  # Dict[Hashable, Dict[str, str]]

    return node_attributes


def disease_node_attribute(df: DataFrame) -> Dict[Hashable, Dict[str, str]]:
    bait_col = "Bait_Disgenet_disease_.1"
    prey_col = "Prey_Disgenet_disease_.1"

    # Get disease annotations for baits
    bait_disease_df = df[["Bait", bait_col]].copy()
    bait_disease_df.dropna(inplace=True)
    bait_disease_df.drop_duplicates(keep="first", inplace=True)
    bait_disease_df.set_index("Bait", inplace=True)

    bait_disease_dict = bait_disease_df.to_dict()  # Dict[str, Dict[str, str]]
    bait_disease_dict = bait_disease_dict[bait_col]  # Dict[str, str]
    unique_bait_diseases = set(bait_disease_dict.values())

    # Get disease annotations for preys
    prey_disease_df = df[["Prey", prey_col]].copy()
    prey_disease_df.dropna(inplace=True)
    prey_disease_df.drop_duplicates(keep="first", inplace=True)
    prey_disease_df.set_index("Prey", inplace=True)

    prey_disease_dict = prey_disease_df.to_dict()  # Dict[str, Dict[str, str]]
    prey_disease_dict = prey_disease_dict[prey_col]  # Dict[str, str]
    unique_prey_diseases = set(prey_disease_dict.values())

    unique_diseases = list(unique_bait_diseases |
                           unique_prey_diseases)  # List[str]

    # The | merges two dictionaries by keys
    # Dict[str, str]
    combined_disease_dict = bait_disease_dict | prey_disease_dict

    # integer --> color in RGBA
    color_map = cm.get_cmap("tab10", len(unique_diseases))  # Colormap

    # disease --> color
    color_dict = {disease: color_map(i)
                  for i, disease in enumerate(unique_diseases)}

    # TODO: Make a legend for the disease coloring

    # node --> {"node_color": disease color}
    # Note: RGBA tuple -> rgba_to_hex(*RGBA) -> rgba_to_hex(R, G, B, A) -> hex
    node_attributes = {
        key: {"node_color": rgba_to_hex(*color_dict[value])}
        for key, value in combined_disease_dict.items()
    }  # Dict[Hashable, Dict[str, str]]

    return node_attributes


def determine_node_coloring(df: DataFrame) -> Dict[Hashable, Dict[str, str]]:
    coloring_selection = controls.node_coloring_selection.value  # str

    if coloring_selection == "none":
        baits = set(df["Bait"].unique())
        preys = set(df["Prey"].unique())
        proteins = baits | preys

        no_attributes = {
            protein: {"node_color": LifecycleColors.undetermined.value,
                      "legend_label": "Other"}
            for protein in proteins
        }  # Dict[str, Dict[str, str]]
        return no_attributes
    elif coloring_selection == "lifecycle stage":
        lifecycle_attributes = lifecycle_node_attributes(df)
        return lifecycle_attributes
    elif coloring_selection == "location":
        location_attributes = location_node_attributes(df)
        return location_attributes
    elif coloring_selection == "disease":
        disease_attributes = disease_node_attribute(df)
        return disease_attributes
    else:
        message = "supported values: lifecycle stage, location, disease, none"
        raise Exception(message)
