# IP, SEC, CORUM, and IP+SEC filtering

- If dataset is IP, subset by: if Interaction_support is IP or both
- If dataset is SEC, subset by: if Interaction_support is SEC or both
- If dataset is IP+SEC, subset by: if Interaction_support is both
- If dataset is CORUM subset by: if in_CORUM_2022 is True
- For combinations of these take the intersection of baits?

## Extras

- IP+SEC might be useful to change to IP & SEC if time permits

# Direct, RNA mediated, RNA shielded filtering

- If it's Direct, subset by: if IP_interaction_type is direct
- If it's RNA mediated, subset by: if IP_interaction_type is mediated
- If it's RNA shielded, subset by: if IP_interaction_type is shielded
- If it's Undetermined, subset by: if IP_interaction_type is undetermined

## Extras

- Add a specific line type for undetermined (can use dotted)

# Protein text input field

- If a user types in a protein name, subset to baits with that protein name
  - If IP only, then filtering baits is fine
  - If SEC is involved, need to look in Bait and Prey columns

# Clustering method

- Self-explanatory, take the graph we have and cluster with one of the methods

# Color node by

- if lifecycle stage, each stage should have following colors: Bait_lifecycle_step for bait and prey_lifecycle_stage_main_by_most_common_bait_stage for prey
- if location, each stage should have following colors: Bait_main_location_HPA for bait and Prey_main_location_HPA for prey
- if disease, each stage should have following colors: use Bait_Disgenet_disease\_.1 for bait and .1 for prey as well (i.e. should be matching)

## Extras

- By default nodes are uncolored (only edges are styled and colored by default)

# Other Notes

- include a basic README.md with instructions of how to run the code
