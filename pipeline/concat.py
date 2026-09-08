"""
pipeline/concat.py
ETAPE 1/7 — Concaténation des fichiers FASTQ Nanopore.

Les données Nanopore sortent du basecalling (Guppy/Dorado) sous forme
de nombreux petits fichiers .fastq.gz par barcode. Cette étape les
fusionne en un seul fichier avant l'alignement.

La recherche est récursive pour gérer toutes les structures possibles :
    barcode01/*.fastq.gz
    barcode01/fastq_pass/*.fastq.gz
    run_folder/barcode01/*.fastq.gz
"""

import os
import glob

from pipeline.utils import log, run_command, check_file_exists


def run_concat(fastq_dir: str, concat_fastq: str) -> str:
    """Concatène tous les fichiers FASTQ d'un dossier barcode en un seul fichier.

    La recherche est récursive — tous les sous-dossiers sont explorés.
    Accepte les formats .fastq.gz et .fastq (mixtes).
    Les fichiers sont triés alphabétiquement pour un ordre reproductible.

    Args:
        fastq_dir:    Chemin vers le dossier contenant les fichiers FASTQ Nanopore
        concat_fastq: Chemin du fichier de sortie concaténé (.fastq.gz)

    Returns:
        Chemin du fichier FASTQ concaténé

    Raises:
        FileNotFoundError si aucun fichier FASTQ n'est trouvé dans fastq_dir
    """
    log("[ETAPE 1/7] Concaténation FASTQ")
    log(f"   → Recherche dans: {fastq_dir}")

    # Recherche récursive dans tous les sous-dossiers
    # Gère toutes les structures de sortie Nanopore :
    #   barcode01/*.fastq.gz
    #   barcode01/fastq_pass/*.fastq.gz
    #   run_folder/barcode01/fastq_pass/*.fastq.gz
    fastq_files = sorted(
        glob.glob(os.path.join(fastq_dir, "**", "*.fastq.gz"), recursive=True) +
        glob.glob(os.path.join(fastq_dir, "**", "*.fastq"),    recursive=True) +
        glob.glob(os.path.join(fastq_dir, "*.fastq.gz")) +
        glob.glob(os.path.join(fastq_dir, "*.fastq"))
    )

    # Dédupliquer (un fichier à la racine peut être trouvé deux fois)
    fastq_files = sorted(set(fastq_files))

    if not fastq_files:
        raise FileNotFoundError(
            f"❌ Aucun fichier FASTQ trouvé dans: {fastq_dir}\n"
            f"   Vérifiez que le dossier contient des fichiers .fastq ou .fastq.gz"
        )

    log(f"   → {len(fastq_files)} fichiers FASTQ détectés")
    for f in fastq_files[:5]:
        log(f"      · {os.path.relpath(f, fastq_dir)}")
    if len(fastq_files) > 5:
        log(f"      · ... et {len(fastq_files) - 5} autres")

    # Concaténation avec cat
    files_str = ' '.join(f'"{f}"' for f in fastq_files)
    run_command(
        f"cat {files_str} > \"{concat_fastq}\"",
        "Erreur concaténation FASTQ",
        timeout=3600
    )

    check_file_exists(concat_fastq, "FASTQ concaténé")
    log(f"✅ Concaténation terminée → {os.path.basename(concat_fastq)}")
    return concat_fastq
