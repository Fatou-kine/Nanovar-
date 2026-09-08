"""
pipeline/picard.py
ETAPE 3/7 — Picard MarkDuplicates + suppression des duplicats
Note: BQSR intentionnellement exclu — non adapté aux données Nanopore
"""
import os
from pipeline.utils import log, run_command, check_file_exists, THREADS, PICARD_JAR


def run_mark_duplicates(bam_file, markdup_bam, dedup_bam, metrics_file):
    log("[ETAPE 3/7] Picard MarkDuplicates + Deduplication")
    log("   ℹ️  BQSR exclu — non adapté Nanopore")
    check_file_exists(bam_file, "Fichier BAM trié")

    run_command(
        f"java -jar \"{PICARD_JAR}\" MarkDuplicates "
        f"I=\"{bam_file}\" O=\"{markdup_bam}\" M=\"{metrics_file}\" "
        f"CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT",
        "Erreur Picard MarkDuplicates", timeout=14400
    )
    check_file_exists(markdup_bam, "BAM MarkDuplicates")
    _log_metrics(metrics_file)

    run_command(
        f"samtools view -@ {THREADS} -F 1024 -b \"{markdup_bam}\" "
        f"| samtools sort -@ {THREADS} -o \"{dedup_bam}\" "
        f"&& samtools index -@ {THREADS} \"{dedup_bam}\"",
        "Erreur suppression duplicats", timeout=7200
    )
    check_file_exists(dedup_bam, "BAM dédupliqué")
    log(f"✅ Picard terminé → {os.path.basename(dedup_bam)}")
    return markdup_bam, dedup_bam, metrics_file


def _log_metrics(metrics_file):
    try:
        with open(metrics_file) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if 'PERCENT_DUPLICATION' in line and i + 1 < len(lines):
                h = line.strip().split('\t')
                v = lines[i+1].strip().split('\t')
                if len(h) == len(v):
                    pct = float(dict(zip(h,v)).get('PERCENT_DUPLICATION', 0)) * 100
                    log(f"   📊 Taux de duplication: {pct:.2f}%")
                break
    except Exception:
        pass
