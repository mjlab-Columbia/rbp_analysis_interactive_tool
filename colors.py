from enum import Enum


class NodeColors(Enum):
    IP: str = "#2B3990"
    SEC: str = "#BCBEC0"
    Both: str = "#F26B21"


class EdgeColors(Enum):
    direct: str = "#808080"
    rna_mediated: str = "#808080"
    rna_shielded: str = "#808080"
    undetermined: str = "#808080"
    ppi_support: str = "#39B54A"


class LifecycleColors(Enum):
    decay: str = "#BE1E2D"
    export: str = "#2CB56B"
    splicing: str = "#FBDC06"
    localization: str = "#F4B063"
    three_prime_end_processing: str = "#00AEEF"
    translation: str = "#354089"
    transcription: str = "#D43C96"
    negative_control: str = "#FF0000"
    modification: str = "#69477F"
    new: str = "#939598"


LifecycleColorsDict = {
    "decay": LifecycleColors.decay,
    "export": LifecycleColors.export,
    "splicing": LifecycleColors.splicing,
    "localization": LifecycleColors.localization,
    "3' end processing": LifecycleColors.three_prime_end_processing,
    "translation": LifecycleColors.translation,
    "new": LifecycleColors.new,
    "transcription": LifecycleColors.transcription,
    "negative control": LifecycleColors.negative_control,
    "modification": LifecycleColors.modification
}  # Dict[str, LifecycleColors]
