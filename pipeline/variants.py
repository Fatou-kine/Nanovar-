"""
pipeline/variants.py
ETAPE 4/7 — Appel de variants avec Clair3 via Docker.

Image  : hkubal/clair3:v1.0.10 (validée et testée)
Modèle : r1041_e82_400bps_sup_v430 (P2 Solo, R10.4.1, SUP 400bps)

IMPORTANT — Solution au deadlock Docker sur Mac (OSError Errno 35) :
Docker Desktop sur Mac génère un Resource deadlock quand plusieurs volumes
sont montés simultanément et qu'un gros fichier indexé (.fai) est lu en
parallèle par un sous-processus.

Solution retenue :
  - Copier BAM, BAI et hg38.fa.fai dans un dossier de travail temporaire
  - Créer un hard link vers hg38.fa (instantané, zéro espace supplémentaire)
  - Monter UNIQUEMENT ce dossier de travail dans Docker (1 seul volume de données)
  - Récupérer les résultats avec shutil.copytree puis nettoyer
"""

import os
import shutil

from pipeline.utils import (
    log, run_command, check_file_exists,
    THREADS, REFERENCE_GENOME,
    CLAIR3_IMAGE, CLAIR3_MODEL_LOCAL
)


def run_clair3(
    dedup_bam: str,
    clair3_output_dir: str,
    sample_name: str,
    ctg_name: str = None,
    bed_file: str = None
) -> str:
    """Appelle les variants avec Clair3 via Docker.

    Args:
        dedup_bam:         BAM dédupliqué en entrée
        clair3_output_dir: Dossier de sortie final des résultats Clair3
        sample_name:       Nom de l'échantillon
        ctg_name:          Chromosome cible (ex: 'chr17') — optionnel
        bed_file:          Fichier BED pour panel de gènes — optionnel

    Returns:
        Chemin du fichier merge_output.vcf.gz

    Raises:
        FileNotFoundError si les fichiers requis sont manquants
        RuntimeError si Docker échoue
    """
    log("[ETAPE 4/7] Clair3 appel de variants via Docker")
    log(f"   Modèle : r1041_e82_400bps_sup_v430 (P2 Solo)")

    check_file_exists(dedup_bam,         "Fichier BAM dédupliqué")
    check_file_exists(REFERENCE_GENOME,  "Génome de référence hg38")
    check_file_exists(CLAIR3_MODEL_LOCAL,"Modèle Clair3 local")

    # ── Préparer le dossier de travail unique ──────────────────────
    work_dir   = os.path.join(os.path.dirname(dedup_bam), "clair3_workdir")
    output_dir = os.path.join(work_dir, "output")
    os.makedirs(work_dir,   exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.chmod(work_dir,   0o777)
    os.chmod(output_dir, 0o777)

    bam_filename = os.path.basename(dedup_bam)
    bai_file     = dedup_bam + ".bai"
    ref_filename = os.path.basename(REFERENCE_GENOME)
    fai_file     = REFERENCE_GENOME + ".fai"

    # ── Copier BAM + BAI ──────────────────────────────────────────
    log("   → Copie BAM vers dossier de travail Docker...")
    shutil.copy2(dedup_bam, os.path.join(work_dir, bam_filename))
    if os.path.exists(bai_file):
        shutil.copy2(bai_file, os.path.join(work_dir, bam_filename + ".bai"))
    else:
        run_command(
            f"samtools index \"{os.path.join(work_dir, bam_filename)}\"",
            "Erreur indexation BAM pour Docker",
            timeout=600
        )

    # ── Copier le .fai (19KB — instantané) ───────────────────────
    if os.path.exists(fai_file):
        shutil.copy2(fai_file, os.path.join(work_dir, ref_filename + ".fai"))
        log("   → hg38.fa.fai copié")
    else:
        raise FileNotFoundError(f"❌ Index hg38.fa.fai manquant: {fai_file}")

    # ── Hard link vers hg38.fa ────────────────────────────────────
    # Instantané, zéro espace supplémentaire, Docker le voit comme un vrai fichier
    # Les symlinks et montages fichier unique échouent sur VirtioFS Mac
    ref_link = os.path.join(work_dir, ref_filename)
    if not os.path.exists(ref_link):
        try:
            os.link(REFERENCE_GENOME, ref_link)
            log("   → hg38.fa lié (hard link)")
        except OSError:
            log("   → hg38.fa : hard link impossible, copie en cours (3GB)...")
            shutil.copy2(REFERENCE_GENOME, ref_link)
            log("   → hg38.fa copié")

    # ── Copier le fichier BED si fourni ──────────────────────────
    bed_in_docker = ""
    if bed_file and os.path.exists(bed_file):
        bed_dest = os.path.join(work_dir, os.path.basename(bed_file))
        shutil.copy2(bed_file, bed_dest)
        bed_in_docker = f"--bed_fn=/data/{os.path.basename(bed_file)}"
        log(f"   → Fichier BED: {os.path.basename(bed_file)}")

    # ── Option chromosome ─────────────────────────────────────────
    ctg_opt = f"--ctg_name={ctg_name}" if ctg_name and not bed_in_docker else ""
    if ctg_name:
        log(f"   → Chromosome ciblé: {ctg_name}")

    # ── Construire les options Clair3 ────────────────────────────
    clair3_opts = (
        f"--bam_fn=/data/{bam_filename} "
        f"--ref_fn=/data/{ref_filename} "
        f"--threads={THREADS} "
        f"--platform=ont "
        f"--model_path=/opt/models/r1041_e82_400bps_sup_v430 "
        f"--output=/data/output "
        f"--enable_phasing "
        f"--use_whatshap_for_final_output_phasing "
        f"{ctg_opt} "
        f"{bed_in_docker}"
    ).strip()

    # ── Commande Docker avec UN SEUL montage de données ──────────
    cmd = (
        f'docker run --rm '
        f'--user $(id -u):$(id -g) '
        f'-v "{work_dir}":/data '
        f'-v "{CLAIR3_MODEL_LOCAL}":/opt/models/r1041_e82_400bps_sup_v430 '
        f'{CLAIR3_IMAGE} '
        f'/bin/bash -c "source activate clair3 && /opt/bin/run_clair3.sh {clair3_opts}"'
    )

    run_command(cmd, "Erreur Clair3 Docker", timeout=57600)

    # ── Récupérer les résultats ───────────────────────────────────
    vcf_in_workdir = os.path.join(output_dir, "merge_output.vcf.gz")
    check_file_exists(vcf_in_workdir, "VCF Clair3 (merge_output.vcf.gz)")

    # Copier tout le dossier output vers la destination finale
    if os.path.exists(clair3_output_dir):
        shutil.rmtree(clair3_output_dir)
    shutil.copytree(output_dir, clair3_output_dir)

    # ── Nettoyer le dossier de travail temporaire ─────────────────
    try:
        shutil.rmtree(work_dir)
        log("   → Dossier de travail Docker nettoyé")
    except Exception as e:
        log(f"   ⚠️  Nettoyage workdir: {e}")

    vcf_output = os.path.join(clair3_output_dir, "merge_output.vcf.gz")
    check_file_exists(vcf_output, "VCF Clair3 final")

    log(f"✅ Clair3 terminé → {vcf_output}")
    return vcf_output
