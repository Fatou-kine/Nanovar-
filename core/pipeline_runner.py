import os
import sys
import time
import subprocess

from config.settings import BASE_DIR, OUTPUT_DIR
from core.genes import generate_bed_from_gene
from core.jobs import (
    update_job,
    get_log_file_path,
    add_log,
)
from core.results import get_sample_results


def parse_progression(line):
    """Extrait la progression depuis une ligne de log du pipeline."""

    etapes = {
        "ETAPE 1/7": (1, "Concaténation FASTQ"),
        "ETAPE 2/7": (2, "Alignement Minimap2"),
        "ETAPE 3/7": (3, "Picard MarkDuplicates"),
        "ETAPE 4/7": (4, "Clair3 - Appel de variants"),
        "ETAPE 5/7": (5, "Annotation ANNOVAR"),
        "ETAPE 6/7": (6, "Conversion Excel"),
        "ETAPE 7/7": (7, "IGV Snapshots"),
    }

    for key, (num, label) in etapes.items():
        if key in line:
            progression = int((num / 7) * 100)
            return progression, label

    return None, None


def run_pipeline_job(
    sample_name,
    fastq_dir,
    mode,
    gene_name=None,
    bed_file=None
):
    """
    Fonction exécutée dans un thread séparé pour chaque échantillon.

    Lance NANOPORE_PIPELINE.py via subprocess et lit les logs
    en temps réel.
    """

    update_job(
        sample_name,
        statut="en_cours",
        debut=time.strftime("%H:%M:%S"),
        progression=0,
        etape="Initialisation..."
    )

    try:

        # Préparer le dossier de sortie de l'échantillon
        sample_output_dir = os.path.join(
            OUTPUT_DIR,
            sample_name
        )

        os.makedirs(
            sample_output_dir,
            exist_ok=True
        )

        # Réinitialiser le fichier de log
        log_path = get_log_file_path(sample_name)

        with open(
            log_path,
            "w",
            encoding="utf-8"
        ) as f:
            f.write(
                f"=== NanoVar — Démarrage {sample_name} — "
                f"{time.strftime('%Y-%m-%d %H:%M:%S')} ===\n"
            )

        # Construire la commande
        cmd = [
            sys.executable,
            os.path.join(
                BASE_DIR,
                "NANOPORE_PIPELINE.py"
            ),
            "--sample",
            sample_name,
            "--fastq",
            fastq_dir,
        ]

        if mode == "gene" and gene_name:

            bed_path, ctg_name = generate_bed_from_gene(
                gene_name,
                sample_output_dir
            )

            if bed_path:
                cmd += [
                    "--bed",
                    bed_path,
                    "--ctg_name",
                    ctg_name
                ]

                add_log(
                    sample_name,
                    f"[INFO] Gène: {gene_name} → {ctg_name}"
                )

            else:
                add_log(
                    sample_name,
                    f"[WARNING] Gène {gene_name} non trouvé "
                    f"dans la base — analyse sur génome complet"
                )

        elif mode == "panel" and bed_file:

            cmd += [
                "--bed",
                bed_file
            ]

            add_log(
                sample_name,
                f"[INFO] Panel BED: {os.path.basename(bed_file)}"
            )

        elif mode == "chromosome" and gene_name:

            cmd += [
                "--ctg_name",
                gene_name
            ]

            add_log(
                sample_name,
                f"[INFO] Chromosome: {gene_name}"
            )

        add_log(
            sample_name,
            f"[INFO] Lancement: {' '.join(cmd)}"
        )

        # Lancer le pipeline
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
            universal_newlines=True,
            encoding="utf-8",
            errors="replace"
        )

        # Stocker le PID
        update_job(
            sample_name,
            pid=process.pid
        )

        # Lire les logs en temps réel
        for line in process.stdout:

            line = line.rstrip()

            if not line:
                continue

            add_log(
                sample_name,
                line
            )

            progression, etape = parse_progression(line)

            if progression is not None:

                update_job(
                    sample_name,
                    progression=progression,
                    etape=etape
                )

            # Détecter les erreurs critiques
            if "ERREUR CRITIQUE" in line or "💥" in line:

                update_job(
                    sample_name,
                    statut="erreur",
                    etape="Erreur critique"
                )

        process.wait()

        # Pipeline terminé
        if process.returncode == 0:

            update_job(
                sample_name,
                statut="termine",
                progression=100,
                etape="Terminé ✓",
                fin=time.strftime("%H:%M:%S"),
                resultats=get_sample_results(sample_name)
            )

            add_log(
                sample_name,
                f"✅ Pipeline terminé avec succès pour {sample_name}"
            )

        else:

            update_job(
                sample_name,
                statut="erreur",
                etape="Erreur — voir les logs",
                fin=time.strftime("%H:%M:%S")
            )

            add_log(
                sample_name,
                f"❌ Pipeline échoué pour {sample_name} "
                f"(code: {process.returncode})"
            )

    except Exception as e:

        update_job(
            sample_name,
            statut="erreur",
            etape=f"Erreur: {str(e)}",
            fin=time.strftime("%H:%M:%S")
        )

        add_log(
            sample_name,
            f"❌ Erreur: {str(e)}"
        )
