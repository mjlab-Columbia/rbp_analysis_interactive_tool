from bokeh.core.properties import DashPattern

from buttons import dataset_combined_checkbox_button, dataset_checkbox_button
from buttons import interaction_type_checkbox_button
from enum import Enum

from pandas import DataFrame, Series

from matplotlib.pyplot import cm

from colors import NodeColors, LifecycleColors
from colors import LifecycleColorsDict

import controls

from typing import Dict, Tuple, Any

from pdb import set_trace as st

all_global_data_columns = [
    'Bait',
    'Prey',
    'IP',
    'Bait_lifecycle_step',
    'IP_experiment',
    'IP_reps',
    'IP_bait_gene',
    'IP_prey_gene',
    'IP_volc_pval_cutoff',
    'IP_volc_log2FC',
    'IP_volc_pvalue',
    'IP_bait_gene_uniprotID',
    'IP_prey_gene_uniprotID',
    'SEC_gene_A',
    'SEC_gene_B',
    'SEC_gene_A_uniprotID',
    'SEC_gene_B_uniprotID',
    'IP_interaction_type',
    'Interaction_support',
    'IP_prey_over_bait_intensity_stoichimoetry',
    'stoich_bait_uniprotID',
    'stoich_prey_uniprotID',
    'totalProteome_prey_over_bait_intensity_stoichimoetry',
    'in_CORUM_2020',
    'CORUM_complex_2020',
    'in_CORUM_2022',
    'CORUM_complex_2022',
    'in_Mentha',
    'Mentha_score',
    'in_HI',
    'in_HIPPIE',
    'HIPPIE_score',
    'in_recY2H',
    'recY2H_score',
    'in_BioID',
    'in_OpenCell',
    'OpenCell_enrichment',
    'in_BioPlex3',
    'BioPlex3_score',
    'PPI_support',
    'bait_eclipped',
    'prey_eclipped',
    'Bait_targets_cnt',
    'Prey_targets_cnt',
    'shared_targets',
    'Bait_main_location_HPA',
    'Prey_main_location_HPA',
    'Bait_main_location.nucl.cyto.',
    'Prey_main_location.nucl.cyto.',
    'Prey_is_RBP_HydRa',
    'Prey_HydRa_score',
    'Prey_PSP_PScore',
    'Bait_Disgenet_disease_.1',
    'Bait_Disgenet_disease_.2',
    'Prey_Disgenet_disease_.1',
    'Prey_Disgenet_disease_.2',
    'Prey_SAFE_localization_BioID',
    'Prey_NMF_localization_BioID',
    'Bait_SAFE_localization_BioID',
    'Bait_NMF_localization_BioID',
    'Prey_NMF_localization_2018_BioID',
    'prey_lifecycle_stage_main_by_most_common_bait_stage',
    'bait_prey_pctLocalizationOverlap_HPA'
]


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


def subset_by_node_type(df: DataFrame) -> DataFrame:
    """
    Return the baits corresponding to the dataset the user selected via
    `dataset_checkbox_button` and `dataset_combined_checkbox_button`
    """

    # Each type of interaction has a color
    df["edge_color"] = df["IP_interaction_type"].map({
        "IP": NodeColors.IP.value,
        "SEC": NodeColors.SEC.value,
        "both": NodeColors.Both.value
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
        raise Exception("Improper combination of dataset checkboxes")


def subset_by_edge_type(df: DataFrame) -> DataFrame:
    """
    Subset df to bait-prey pairs based on protein interaction type
    """

    # Note: DashPatterns take iterable as [dash, gap] ([] is a solid line)
    df["edge_style"] = df["IP_interaction_type"].map({
        InteractionTypeValue.direct.value: DashPattern([]),
        InteractionTypeValue.mediated.value: DashPattern([8, 4]),
        InteractionTypeValue.shielded.value: DashPattern([5, 20]),
        InteractionTypeValue.undetermined.value: DashPattern([1, 1])
    })
    interaction_type = interaction_type_checkbox_button.active  # List[int]
    types_to_include = []

    if len(interaction_type) == 0:
        message = "Supported values: Direct, RNA mediated, or RNA shielded"
        raise Exception(message)

    # TODO: This logic can be simpler, so figure out a simplification
    if InteractionTypeCheckbox.direct.value in interaction_type:
        types_to_include.append(InteractionTypeValue.direct.value)

    if InteractionTypeCheckbox.mediated.value in interaction_type:
        types_to_include.append(InteractionTypeValue.mediated.value)

    if InteractionTypeCheckbox.shielded.value in interaction_type:
        types_to_include.append(InteractionTypeValue.shielded.value)

    st()
    mask = df["IP_interaction_type"].isin(types_to_include)
    filtered_df = df[mask].copy()
    return filtered_df


def subset_by_protein(df: DataFrame) -> DataFrame:
    """
    Subset df to bait proteins list specified in protein_text_input
    """
    protein_str = str(controls.protein_text_input.value).strip()  # str
    protein_list = [p.strip() for p in protein_str.split(",")]  # List[str]
    if (len(protein_list) == 0) and (protein_list[0] == ""):
        return df
    else:
        return df[df["Bait"].isin(protein_list)]


def lifecycle_node_attributes(df: DataFrame) -> Dict[str, Dict[str, str]]:
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

    # Map each bait node to its lifecycle color
    bait_colors = {
        key: str(LifecycleColorsDict[value])
        for key, value in bait_lifecycle_dict.items()
    }  # Dict[str, str]

    # Colors for prey proteins
    prey_lifecycle_df = df[["Prey", prey_col]].copy()
    prey_lifecycle_df.dropna(inplace=True)
    prey_lifecycle_df.drop_duplicates(keep="first", inplace=True)
    prey_lifecycle_df.set_index("Prey", inplace=True)

    # Dict[str, Dict[str, str]]
    prey_lifecycle_dict = prey_lifecycle_df.to_dict()
    prey_lifecycle_dict = prey_lifecycle_dict[prey_col]  # Dict[str, str]

    # Map each prey node to its lifecycle color
    prey_colors = {
        key: str(LifecycleColorsDict[value])
        for key, value in prey_lifecycle_dict.items()
    }  # Dict[str, str]

    # The | operation merges two dictionaries by keys
    protein_colors = bait_colors | prey_colors

    node_attributes = {
        key: {"node_color": value}
        for key, value in protein_colors
    }  # Dict[str, Dict[str, LifecycleColors]]

    return node_attributes


def location_node_attributes(df: DataFrame) -> Dict[str, Dict[str, Tuple[float, float, float, float]]]:
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

    # integer --> color
    color_map = cm.get_cmap("tab10", len(unique_locations))

    # location --> color
    color_dict = {location: color_map(i)
                  for i, location in enumerate(unique_locations)}

    # node --> {"node_attribute": color}
    node_attributes = {
        key: {"node_attribute": color_dict[value]}
        for key, value in combined_location_dict.items()
    }  # Dict[str, Dict[str, Tuple[float, float, float, float]]]

    return node_attributes


def disease_node_attribute(df: DataFrame) -> Dict[str, Dict[str, Tuple[float, float, float, float]]]:
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

    # integer --> color
    color_map = cm.get_cmap("tab10", len(unique_diseases))

    # disease --> color
    color_dict = {disease: color_map(i)
                  for i, disease in enumerate(unique_diseases)}

    # node --> {"node_attribute": disease color}
    node_attributes = {
        key: {"node_attributes": color_dict[value]}
        for key, value in combined_disease_dict.items()
    }  # Dict[str, Dict[str, Tuple[float, float, float, float]]]

    return node_attributes


# TODO: Create type aliases to replace `Any` type annotation
def determine_node_coloring(df: DataFrame) -> Tuple[DataFrame, Any]:
    coloring_selection = controls.node_coloring_selection

    if coloring_selection.value == "none":
        baits = set(df["Baits"].unique())
        preys = set(df["Preys"].unique())
        proteins = baits | preys

        no_attributes = {
            protein: {"node_color": ""}
            for protein in proteins
        }  # Dict[str, Dict[str, str]]
        return df, no_attributes
    elif coloring_selection.value == "lifecycle stage":
        lifecycle_attributes = lifecycle_node_attributes(df)
        return df, lifecycle_attributes
    elif coloring_selection.value == "location":
        location_attributes = location_node_attributes(df)
        return df, location_attributes
    elif coloring_selection.value == "disease":
        disease_attributes = disease_node_attribute(df)
        return df, disease_attributes
    else:
        message = "supported values: lifecycle stage, location, disease, none"
        raise Exception(message)
