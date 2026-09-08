"""
NANOPORE_PIPELINE.py — Point d'entrée principal NanoVar
Appelé par app.py (Flask) ou directement en ligne de commande.

Usage:
    python3 NANOPORE_PIPELINE.py --sample NOM --fastq /chemin/fastq
    python3 NANOPORE_PIPELINE.py --sample NOM --fastq /chemin --ctg_name chr17
    python3 NANOPORE_PIPELINE.py --sample NOM --fastq /chemin --bed /chemin/panel.bed
"""
import os
import sys
import time
import argparse

from pipeline import (
    log, validate_configuration, check_existing_outputs, safe_cleanup,
    run_concat, run_minimap2, run_mark_duplicates, run_clair3,
    run_annovar, convert_to_excel, run_igv_snapshot,
)
from pipeline.utils import OUTPUT_DIR


def process_sample(sample_name, fastq_dir, ctg_name=None, bed_file=None):
    """Lance le pipeline complet pour un échantillon."""
    start_time = time.time()

    log(f"{'='*60}")
    log(f"🚀 TRAITEMENT: {sample_name}")
    log(f"   Modèle Clair3 : r1041_e82_400bps_sup_v430 (P2 Solo)")
    log(f"   BQSR          : exclu (non adapté Nanopore)")
    if ctg_name: log(f"   Chromosome    : {ctg_name}")
    if bed_file:  log(f"   BED           : {os.path.basename(bed_file)}")
    log(f"{'='*60}")

    sample_dir = os.path.join(OUTPUT_DIR, sample_name)
    os.makedirs(sample_dir, exist_ok=True)

    log("🔍 Vérification des étapes déjà terminées...")
    done = check_existing_outputs(sample_dir, sample_name)

    concat_fastq   = os.path.join(sample_dir, f"{sample_name}.fastq.gz")
    bam_file       = os.path.join(sample_dir, f"{sample_name}.bam")
    markdup_bam    = os.path.join(sample_dir, f"{sample_name}_markdup.bam")
    dedup_bam      = os.path.join(sample_dir, f"{sample_name}_dedup.bam")
    metrics_file   = os.path.join(sample_dir, f"{sample_name}_markdup_metrics.txt")
    clair3_out_dir = os.path.join(sample_dir, "clair3_output")
    annovar_prefix = os.path.join(sample_dir, f"{sample_name}_annovar")
    excel_file     = os.path.join(sample_dir, f"{sample_name}_variants.xlsx")

    # ── ÉTAPE 1 ─────────────────────────────────────────────────
    if 'concat' not in done:
        concat_fastq = run_concat(fastq_dir, concat_fastq)
    else:
        concat_fastq = done['concat']
        log("1/7 Concaténation — Déjà terminé ✓")

    # ── ÉTAPE 2 ─────────────────────────────────────────────────
    if 'minimap2' not in done:
        bam_file = run_minimap2(concat_fastq, bam_file, sample_name)
    else:
        bam_file = done['minimap2']
        log("2/7 Minimap2 — Déjà terminé ✓")

    # ── ÉTAPE 3 ─────────────────────────────────────────────────
    if 'dedup' not in done:
        markdup_bam, dedup_bam, metrics_file = run_mark_duplicates(
            bam_file, markdup_bam, dedup_bam, metrics_file
        )
    else:
        dedup_bam = done['dedup']
        log("3/7 Picard MarkDuplicates — Déjà terminé ✓")

    # ── ÉTAPE 4 ─────────────────────────────────────────────────
    if 'clair3' not in done:
        clair3_vcf = run_clair3(
            dedup_bam, clair3_out_dir, sample_name,
            ctg_name=ctg_name, bed_file=bed_file
        )
    else:
        clair3_vcf = done['clair3']
        log("4/7 Clair3 — Déjà terminé ✓")

    # ── ÉTAPE 5 ─────────────────────────────────────────────────
    if 'annovar' not in done:
        annovar_txt = run_annovar(clair3_vcf, annovar_prefix)
    else:
        annovar_txt = done['annovar']
        log("5/7 ANNOVAR — Déjà terminé ✓")

    # ── ÉTAPE 6 ─────────────────────────────────────────────────
    excel_file = convert_to_excel(annovar_txt, excel_file)

    # ── ÉTAPE 7 ─────────────────────────────────────────────────
    run_igv_snapshot(dedup_bam, clair3_vcf, sample_dir, sample_name)

    safe_cleanup(sample_dir, sample_name)

    elapsed = time.time() - start_time
    log(f"✅ PIPELINE TERMINÉ: {sample_name} en {elapsed/60:.1f} minutes")

    return {
        "sample":  sample_name,
        "annovar": annovar_txt,
        "excel":   excel_file,
        "bam":     dedup_bam,
        "vcf":     clair3_vcf,
        "statut":  "termine"
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NanoVar — Pipeline Nanopore")
    parser.add_argument("--sample",   required=True)
    parser.add_argument("--fastq",    required=True)
    parser.add_argument("--ctg_name", required=False, default=None)
    parser.add_argument("--bed",      required=False, default=None)
    args = parser.parse_args()

    validate_configuration()

    if not os.path.exists(args.fastq):
        log(f"❌ Dossier FASTQ non trouvé: {args.fastq}")
        sys.exit(1)

    try:
        result = process_sample(
            args.sample, args.fastq,
            ctg_name=args.ctg_name, bed_file=args.bed
        )
        log(f"🎉 Succès: {result['sample']}")
        sys.exit(0)
    except Exception as e:
        log(f"💥 ERREUR CRITIQUE: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
