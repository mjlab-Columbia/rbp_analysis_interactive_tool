from bokeh.models import CheckboxButtonGroup, Button
from enum import Enum

# Generates the interactome based on the user's input
create_interactome_button = Button(
    label="Create Interactome",
    button_type="success"
)  # Button

# Download button
download_button = Button(
    label="Download",
    button_type="success"
)  # Button

# Checkbox buttons for selecting the dataset
dataset_checkbox_button = CheckboxButtonGroup(
    labels=['IP', 'SEC', 'CORUM'])  # CheckboxButtonGroup
dataset_combined_checkbox_button = CheckboxButtonGroup(
    labels=['IP+SEC'])  # CheckboxButtonGroup
interaction_type_checkbox_button = CheckboxButtonGroup(
    labels=["Direct", "RNA mediated", "RNA shielded"])  # CheckboxButtonGroup

all_buttons = [
    create_interactome_button,
    download_button,
    dataset_checkbox_button,
    dataset_combined_checkbox_button,
    interaction_type_checkbox_button
]  # List[CheckboxButtonGroup]
