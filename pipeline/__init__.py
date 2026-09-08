"""
pipeline/ — Modules du pipeline NanoVar

    utils.py      — Configuration, logging, run_command
    concat.py     — Étape 1 : concaténation FASTQ
    align.py      — Étape 2 : Minimap2
    picard.py     — Étape 3 : Picard MarkDuplicates
    variants.py   — Étape 4 : Clair3 (Docker)
    annotation.py — Étape 5 : ANNOVAR
    export.py     — Étapes 6-7 : Excel + IGV
"""

from pipeline.utils      import log, run_command, check_file_exists
from pipeline.utils      import check_existing_outputs, safe_cleanup, validate_configuration
from pipeline.concat     import run_concat
from pipeline.align      import run_minimap2
from pipeline.picard     import run_mark_duplicates
from pipeline.variants   import run_clair3
from pipeline.annotation import run_annovar
from pipeline.export     import convert_to_excel, run_igv_snapshot

__all__ = [
    "log", "run_command", "check_file_exists",
    "check_existing_outputs", "safe_cleanup", "validate_configuration",
    "run_concat", "run_minimap2", "run_mark_duplicates",
    "run_clair3", "run_annovar", "convert_to_excel", "run_igv_snapshot",
]
