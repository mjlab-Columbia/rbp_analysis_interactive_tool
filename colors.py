from enum import Enum

GRAY = "#808080"
PPI_SUPPORT = "#39B54A"


class EdgeColors(Enum):
    IP: str = "#2B3990"
    SEC: str = "#BCBEC0"
    Both: str = "#F26B21"


# class EdgeStyles(Enum):
#     direct: DashPattern = DashPattern([])
#     rna_mediated: DashPattern = DashPattern([8, 4])
#     rna_shielded: DashPattern = DashPattern([5, 20])
#     undetermined: DashPattern = DashPattern([1, 1])


class LifecycleColors(Enum):
    decay: str = "#BE1E2D"
    export: str = "#2CB56B"
    splicing: str = "#FBDC06"
    localization: str = "#F4B063"
    three_prime_end_processing: str = "#00AEEF"
    translation: str = "#354089"
    transcription: str = "#D43C96"
    negative_control: str = "#939598"
    modification: str = "#69477F"
    new: str = "#939598"
    undetermined: str = "#939598"


LifecycleColorsDict = {
    "decay": LifecycleColors.decay.value,
    "export": LifecycleColors.export.value,
    "splicing": LifecycleColors.splicing.value,
    "localization": LifecycleColors.localization.value,
    "3' end processing": LifecycleColors.three_prime_end_processing.value,
    "translation": LifecycleColors.translation.value,
    "new": LifecycleColors.new.value,
    "transcription": LifecycleColors.transcription.value,
    "negative control": LifecycleColors.negative_control.value,
    "modification": LifecycleColors.modification.value,
    "undetermined": LifecycleColors.undetermined.value
}  # Dict[str, str]
