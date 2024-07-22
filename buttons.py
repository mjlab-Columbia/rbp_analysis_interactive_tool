from bokeh.models import CheckboxButtonGroup, Button
from enum import Enum

# Generates the interactome based on the user's input
create_interactome_button = Button(
    label="Create Interactome",
    button_type="success",
    width_policy="max"
)  # Button

# Download button
download_button = Button(
    label="Download",
    button_type="success",
    width_policy="max"
)  # Button

# Checkbox buttons for selecting the dataset
dataset_checkbox_button = CheckboxButtonGroup(labels=['IP', 'SEC'],
                                              width_policy="max")
interaction_support_checkbox_button = CheckboxButtonGroup(labels=['IP+SEC', 'CORUM'],
                                                          width_policy="max")
interaction_type_checkbox_button = CheckboxButtonGroup(labels=["Direct", "RNA mediated", "RNA shielded"],
                                                       width_policy="max")

all_buttons = [
    create_interactome_button,
    download_button,
    dataset_checkbox_button,
    interaction_support_checkbox_button,
    interaction_type_checkbox_button
]  # List[CheckboxButtonGroup]
