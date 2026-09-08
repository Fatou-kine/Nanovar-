"""
pipeline/align.py
ETAPE 2/7 — Alignement Minimap2 (map-ont) + SAM→BAM + tri + indexation
"""
import os
import subprocess
from pipeline.utils import log, run_command, check_file_exists, THREADS, REFERENCE_GENOME


def run_minimap2(concat_fastq, bam_file, sample_name):
    """Aligne les reads Nanopore sur hg38 avec Minimap2 en mode map-ont.
    Pipeline en une passe : minimap2 → samtools view → sort → index
    """
    log("[ETAPE 2/7] Minimap2 alignement + SAM→BAM + tri + index")
    check_file_exists(concat_fastq, "FASTQ concaténé")
    check_file_exists(REFERENCE_GENOME, "Génome de référence hg38")

    rg_tag = (
        f"@RG\\tID:{sample_name}\\tSM:{sample_name}"
        f"\\tPL:ONT\\tLB:lib1\\tPU:unit1"
    )
    cmd = (
        f"minimap2 -a -x map-ont -t {THREADS} -R '{rg_tag}' "
        f"\"{REFERENCE_GENOME}\" \"{concat_fastq}\" "
        f"| samtools view -@ {THREADS} -b -S "
        f"| samtools sort -@ {THREADS} -o \"{bam_file}\" "
        f"&& samtools index -@ {THREADS} \"{bam_file}\""
    )
    run_command(cmd, "Erreur Minimap2 / Samtools", timeout=28800)
    check_file_exists(bam_file, "Fichier BAM trié et indexé")
    _log_alignment_stats(bam_file)
    log(f"✅ Alignement terminé → {os.path.basename(bam_file)}")
    return bam_file


def _log_alignment_stats(bam_file):
    try:
        result = subprocess.run(
            f"samtools flagstat \"{bam_file}\"",
            shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=120
        )
        if result.returncode == 0:
            for line in result.stdout.decode('utf-8', errors='replace').split('\n')[:4]:
                if line.strip():
                    log(f"   📊 {line.strip()}")
    except Exception:
        pass
