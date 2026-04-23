import subprocess
import os
import sys
import time
import glob
import argparse
import yaml

# -----------------------------
# CHARGEMENT DE LA CONFIGURATION
# -----------------------------

def load_config():
    """Charge config.yaml depuis le même dossier que ce script"""
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.yaml")
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"❌ config.yaml non trouvé: {config_path}")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

CONFIG             = load_config()
THREADS            = CONFIG['threads']
REFERENCE_GENOME   = CONFIG['reference_genome']
OUTPUT_DIR         = CONFIG['output_dir']
PICARD_JAR         = CONFIG['picard_jar']
ANNOVAR_DIR        = CONFIG['annovar_dir']
HUMANDB            = CONFIG['humandb']
CLAIR3_IMAGE       = CONFIG['clair3_docker_image']
CLAIR3_MODEL       = CONFIG['clair3_model']
CLAIR3_MODEL_LOCAL = CONFIG['clair3_model_local']
IGV_PATH           = CONFIG.get('igv_path', '')

os.makedirs(OUTPUT_DIR, exist_ok=True)

# -----------------------------
# FONCTIONS UTILITAIRES
# -----------------------------

def log(message):
    """Affiche un message horodaté — lu par Flask pour la progression"""
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)

def check_file_exists(file_path, description):
    """Vérifie qu'un fichier existe et n'est pas vide"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"❌ {description} non trouvé: {file_path}")
    if os.path.getsize(file_path) == 0:
        raise ValueError(f"❌ {description} est vide: {file_path}")
    log(f"✅ {description} vérifié")
    return True

def run_command(cmd, error_message, timeout=14400):
    """Exécute une commande shell avec gestion d'erreur et timeout"""
    if isinstance(cmd, list):
        full_command = ' '.join(cmd)
    else:
        full_command = cmd

    log(f"🔧 {full_command}")

    try:
        start_time = time.time()
        result = subprocess.run(
            full_command,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            timeout=timeout,
            executable='/bin/bash'
        )

        stderr_output = result.stderr.decode('utf-8') if result.stderr else ""
        stdout_output = result.stdout.decode('utf-8') if result.stdout else ""
        elapsed = time.time() - start_time

        log(f"✅ Terminé en {elapsed:.2f}s")

        if result.returncode != 0:
            log(f"🔍 STDERR: {stderr_output}")
            log(f"🔍 STDOUT: {stdout_output}")
            raise RuntimeError(f"{error_message}: {stderr_output}")

        return result

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Timeout après {timeout}s: {error_message}")

def check_existing_outputs(sample_dir, sample_name):
    """Vérifie les étapes déjà terminées pour reprise en cas d'interruption"""
    expected_files = {
        'concat':   os.path.join(sample_dir, f"{sample_name}.fastq.gz"),
        'minimap2': os.path.join(sample_dir, f"{sample_name}.bam"),
        'dedup':    os.path.join(sample_dir, f"{sample_name}_dedup.bam"),
        'clair3':   os.path.join(sample_dir, "clair3_output", "merge_output.vcf.gz"),
        'annovar':  os.path.join(sample_dir, f"{sample_name}_annovar.hg38_multianno.txt")
    }

    completed_steps = {}
    for step, file_path in expected_files.items():
        if os.path.exists(file_path) and os.path.getsize(file_path) > 1000:
            completed_steps[step] = file_path
            log(f"   ✅ Étape '{step}' déjà terminée")
        else:
            log(f"   ⏳ Étape '{step}' à exécuter")

    return completed_steps

def safe_cleanup(sample_dir, sample_name):
    """Supprime les fichiers intermédiaires après succès"""
    log("🧹 Nettoyage des fichiers intermédiaires...")

    files_to_keep = [
        os.path.join(sample_dir, f"{sample_name}.fastq.gz"),
        os.path.join(sample_dir, f"{sample_name}_dedup.bam"),
        os.path.join(sample_dir, f"{sample_name}_dedup.bam.bai"),
        os.path.join(sample_dir, f"{sample_name}_markdup_metrics.txt"),
        os.path.join(sample_dir, f"{sample_name}_annovar.hg38_multianno.txt"),
        os.path.join(sample_dir, f"{sample_name}_annovar.hg38_multianno.vcf"),
        os.path.join(sample_dir, f"{sample_name}_annovar.avinput"),
        os.path.join(sample_dir, f"{sample_name}_variants.xlsx"),
    ]

    clair3_dir = os.path.join(sample_dir, "clair3_output")
    deleted = 0

    for file_path in glob.glob(os.path.join(sample_dir, "*")):
        if file_path == clair3_dir:
            continue
        if file_path not in files_to_keep and os.path.isfile(file_path):
            try:
                os.remove(file_path)
                deleted += 1
            except Exception as e:
                log(f"⚠️  Impossible de supprimer {os.path.basename(file_path)}: {e}")

    log(f"✅ Nettoyage terminé — {deleted} fichiers supprimés")

def validate_configuration():
    """Valide tous les chemins au démarrage"""
    log("🔍 Validation de la configuration...")

    checks = [
        (REFERENCE_GENOME,   "Génome de référence hg38"),
        (PICARD_JAR,         "Picard JAR"),
        (ANNOVAR_DIR,        "Répertoire ANNOVAR"),
        (HUMANDB,            "Base de données ANNOVAR humandb"),
        (CLAIR3_MODEL_LOCAL, "Modèle Clair3 local"),
    ]

    for path, description in checks:
        if not os.path.exists(path):
            raise FileNotFoundError(f"❌ {description} non trouvé: {path}")
        log(f"   ✅ {description}")

    # Vérifier que Docker est disponible et démarré
    result = subprocess.run(
        "docker info",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    if result.returncode != 0:
        raise RuntimeError("❌ Docker n'est pas disponible ou n'est pas démarré")
    log("   ✅ Docker disponible")

    log("✅ Configuration validée")

# -----------------------------
# FONCTIONS DU PIPELINE
# -----------------------------

def run_concat(fastq_dir, concat_fastq):
    """ETAPE 1/7 — Concaténation des fichiers FASTQ Nanopore"""
    log("[ETAPE 1/7] Concaténation FASTQ")

    fastq_files = sorted(
        glob.glob(os.path.join(fastq_dir, "*.fastq.gz")) +
        glob.glob(os.path.join(fastq_dir, "*.fastq"))
    )

    if not fastq_files:
        raise FileNotFoundError(f"❌ Aucun fichier FASTQ dans: {fastq_dir}")

    log(f"   → {len(fastq_files)} fichiers FASTQ détectés")

    files_str = ' '.join(fastq_files)
    run_command(
        f"cat {files_str} > {concat_fastq}",
        "Erreur concaténation FASTQ",
        timeout=3600
    )

    check_file_exists(concat_fastq, "FASTQ concaténé")
    return concat_fastq

def run_minimap2(concat_fastq, bam_file, sample_name):
    """ETAPE 2/7 — Alignement Minimap2 (map-ont) + SAM→BAM + tri + index"""
    log("[ETAPE 2/7] Minimap2 alignement + SAM→BAM + tri + index")

    check_file_exists(concat_fastq, "FASTQ concaténé")
    check_file_exists(REFERENCE_GENOME, "Génome de référence")

    rg_tag = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ONT\\tLB:lib1\\tPU:unit1"

    cmd = (
        f"minimap2 -a -x map-ont -t {THREADS} -R '{rg_tag}' "
        f"{REFERENCE_GENOME} {concat_fastq} "
        f"| samtools view -@ {THREADS} -b -S "
        f"| samtools sort -@ {THREADS} -o {bam_file} "
        f"&& samtools index -@ {THREADS} {bam_file}"
    )

    run_command(cmd, "Erreur Minimap2 / Samtools", timeout=28800)
    check_file_exists(bam_file, "Fichier BAM trié et indexé")
    return bam_file

def run_mark_duplicates(bam_file, markdup_bam, dedup_bam, metrics_file):
    """ETAPE 3/7 — Picard MarkDuplicates + suppression des duplicats
    Note: BQSR intentionnellement exclu — non adapté aux données Nanopore
    """
    log("[ETAPE 3/7] Picard MarkDuplicates + Deduplication")

    check_file_exists(bam_file, "Fichier BAM trié")

    # Marquage des duplicats
    markdup_cmd = (
        f"java -jar {PICARD_JAR} MarkDuplicates "
        f"I={bam_file} O={markdup_bam} M={metrics_file} "
        f"CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT"
    )
    run_command(markdup_cmd, "Erreur Picard MarkDuplicates", timeout=14400)
    check_file_exists(markdup_bam, "Fichier BAM MarkDuplicates")

    # Suppression des duplicats — flag -F 1024 exclut les reads PCR_DUPLICATE
    dedup_cmd = (
        f"samtools view -@ {THREADS} -F 1024 -b {markdup_bam} "
        f"| samtools sort -@ {THREADS} -o {dedup_bam} "
        f"&& samtools index -@ {THREADS} {dedup_bam}"
    )
    run_command(dedup_cmd, "Erreur suppression duplicats", timeout=7200)
    check_file_exists(dedup_bam, "Fichier BAM dédupliqué")

    return markdup_bam, dedup_bam, metrics_file

def run_clair3(dedup_bam, clair3_output_dir, sample_name, ctg_name=None, bed_file=None):
    """ETAPE 4/7 — Appel de variants avec Clair3 via Docker
    Image  : hkubal/clair3:v1.0.10 (validée et testée)
    Modèle : r1041_e82_400bps_sup_v430 (P2 Solo, R10.4.1, SUP 400bps)
    Monté depuis le disque local vers /opt/models/ dans le conteneur
    """
    log("[ETAPE 4/7] Clair3 appel de variants via Docker")

    check_file_exists(dedup_bam, "Fichier BAM dédupliqué")
    check_file_exists(REFERENCE_GENOME, "Génome de référence")
    check_file_exists(CLAIR3_MODEL_LOCAL, "Modèle Clair3 local")

    os.makedirs(clair3_output_dir, exist_ok=True)
    os.chmod(clair3_output_dir, 0o777)

    sample_dir   = os.path.dirname(dedup_bam)
    ref_dir      = os.path.dirname(REFERENCE_GENOME)
    ref_filename = os.path.basename(REFERENCE_GENOME)
    bam_filename = os.path.basename(dedup_bam)

    # Construction des options Clair3
    clair3_opts = (
        f"--bam_fn=/data/sample/{bam_filename} "
        f"--ref_fn=/data/ref/{ref_filename} "
        f"--threads={THREADS} "
        f"--platform=ont "
        f"--model_path=/opt/models/r1041_e82_400bps_sup_v430 "
        f"--output=/data/output "
        f"--no_phasing_for_fa"
    )

    # Limiter à un chromosome (gène unique)
    if ctg_name:
        clair3_opts += f" --ctg_name={ctg_name}"
        log(f"   → Chromosome ciblé: {ctg_name}")

    # Utiliser un fichier BED (panel de gènes)
    if bed_file and os.path.exists(bed_file):
        bed_filename = os.path.basename(bed_file)
        clair3_opts += f" --bed_fn=/data/sample/{bed_filename}"
        log(f"   → Fichier BED: {bed_filename}")

    # Commande Docker complète — validée et testée
    cmd = (
        f'docker run --rm '
        f'--user $(id -u):$(id -g) '
        f'-v {sample_dir}:/data/sample '
        f'-v {ref_dir}:/data/ref '
        f'-v {clair3_output_dir}:/data/output '
        f'-v {CLAIR3_MODEL_LOCAL}:/opt/models/r1041_e82_400bps_sup_v430 '
        f'{CLAIR3_IMAGE} '
        f'/bin/bash -c "source activate clair3 && /opt/bin/run_clair3.sh {clair3_opts}"'
    )

    run_command(cmd, "Erreur Clair3 Docker", timeout=57600)

    vcf_output = os.path.join(clair3_output_dir, "merge_output.vcf.gz")
    check_file_exists(vcf_output, "VCF Clair3 (merge_output.vcf.gz)")

    log(f"✅ Clair3 terminé: {vcf_output}")
    return vcf_output

def run_annovar(vcf_file, annovar_prefix):
    """ETAPE 5/7 — Annotation des variants avec ANNOVAR
    Bases: refGene, ClinVar 2024, gnomAD genome, dbNSFP 4.7a, avsnp151
    """
    log("[ETAPE 5/7] ANNOVAR annotation")

    check_file_exists(vcf_file, "Fichier VCF Clair3")

    cmd = (
        f"perl {os.path.join(ANNOVAR_DIR, 'table_annovar.pl')} "
        f"{vcf_file} {HUMANDB} "
        f"-buildver hg38 "
        f"-out {annovar_prefix} "
        f"-protocol refGene,clinvar_20240611,gnomad_genome,dbnsfp47a,avsnp151 "
        f"-operation g,f,f,f,f "
        f"-nastring . "
        f"-vcfinput"
    )

    run_command(cmd, "Erreur ANNOVAR", timeout=14400)

    annovar_txt = f"{annovar_prefix}.hg38_multianno.txt"
    annovar_vcf = f"{annovar_prefix}.hg38_multianno.vcf"

    check_file_exists(annovar_txt, "Fichier TXT ANNOVAR")
    check_file_exists(annovar_vcf, "Fichier VCF ANNOVAR")

    log("✅ ANNOVAR terminé")
    return annovar_txt

def convert_to_excel(annovar_txt, excel_file):
    """ETAPE 6/7 — Conversion ANNOVAR → Excel (.xlsx)"""
    log("[ETAPE 6/7] Conversion Excel")

    check_file_exists(annovar_txt, "Fichier TXT ANNOVAR")

    try:
        import pandas as pd

        df = pd.read_csv(annovar_txt, sep='\t', low_memory=False)

        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df.to_excel(writer, index=False, sheet_name='Variants')
            worksheet = writer.sheets['Variants']
            for col in worksheet.columns:
                max_length = max(
                    len(str(cell.value)) if cell.value else 0 for cell in col
                )
                worksheet.column_dimensions[col[0].column_letter].width = min(
                    max_length + 2, 50
                )

        log(f"✅ Excel généré: {excel_file}")
        return excel_file

    except ImportError:
        log("⚠️  pandas/openpyxl non installé — pip install pandas openpyxl")
        return None

def run_igv_snapshot(dedup_bam, vcf_file, sample_dir, sample_name):
    """ETAPE 7/7 — Snapshots IGV (optionnel, non bloquant)"""
    log("[ETAPE 7/7] IGV snapshots")

    if not IGV_PATH or not os.path.exists(IGV_PATH):
        log("⚠️  IGV non configuré — étape ignorée")
        return None

    snapshot_dir = os.path.join(sample_dir, "igv_snapshots")
    os.makedirs(snapshot_dir, exist_ok=True)

    igv_batch = os.path.join(sample_dir, f"{sample_name}_igv_batch.txt")
    with open(igv_batch, 'w') as f:
        f.write("new\n")
        f.write("genome hg38\n")
        f.write(f"load {dedup_bam}\n")
        f.write(f"load {vcf_file}\n")
        f.write(f"snapshotDirectory {snapshot_dir}\n")
        f.write(f"snapshot {sample_name}_overview.png\n")
        f.write("exit\n")

    try:
        run_command(f"bash {IGV_PATH} --batch {igv_batch}", "Erreur IGV", timeout=3600)
        log(f"✅ IGV snapshots: {snapshot_dir}")
    except Exception as e:
        log(f"⚠️  IGV échoué (non bloquant): {e}")

    return snapshot_dir

# -----------------------------
# PIPELINE PRINCIPAL
# -----------------------------

def process_sample(sample_name, fastq_dir, ctg_name=None, bed_file=None):
    """Lance le pipeline complet pour un échantillon"""
    start_time = time.time()

    log(f"{'='*60}")
    log(f"🚀 TRAITEMENT: {sample_name}")
    log(f"   Modèle Clair3 : r1041_e82_400bps_sup_v430 (P2 Solo)")
    log(f"   BQSR          : exclu (non adapté Nanopore)")
    if ctg_name:
        log(f"   Chromosome    : {ctg_name}")
    if bed_file:
        log(f"   BED           : {os.path.basename(bed_file)}")
    log(f"{'='*60}")

    sample_dir = os.path.join(OUTPUT_DIR, sample_name)
    os.makedirs(sample_dir, exist_ok=True)

    log("🔍 Vérification des étapes déjà terminées...")
    completed_steps = check_existing_outputs(sample_dir, sample_name)

    # Chemins des fichiers
    concat_fastq   = os.path.join(sample_dir, f"{sample_name}.fastq.gz")
    bam_file       = os.path.join(sample_dir, f"{sample_name}.bam")
    markdup_bam    = os.path.join(sample_dir, f"{sample_name}_markdup.bam")
    dedup_bam      = os.path.join(sample_dir, f"{sample_name}_dedup.bam")
    metrics_file   = os.path.join(sample_dir, f"{sample_name}_markdup_metrics.txt")
    clair3_out_dir = os.path.join(sample_dir, "clair3_output")
    annovar_prefix = os.path.join(sample_dir, f"{sample_name}_annovar")
    excel_file     = os.path.join(sample_dir, f"{sample_name}_variants.xlsx")

    # ── ÉTAPE 1 : Concaténation ─────────────────────────────────
    if 'concat' not in completed_steps:
        concat_fastq = run_concat(fastq_dir, concat_fastq)
    else:
        concat_fastq = completed_steps['concat']
        log("1/7 Concaténation — Déjà terminé ✓")

    # ── ÉTAPE 2 : Minimap2 ──────────────────────────────────────
    if 'minimap2' not in completed_steps:
        bam_file = run_minimap2(concat_fastq, bam_file, sample_name)
    else:
        bam_file = completed_steps['minimap2']
        log("2/7 Minimap2 — Déjà terminé ✓")

    # ── ÉTAPE 3 : Picard ────────────────────────────────────────
    if 'dedup' not in completed_steps:
        markdup_bam, dedup_bam, metrics_file = run_mark_duplicates(
            bam_file, markdup_bam, dedup_bam, metrics_file
        )
    else:
        dedup_bam = completed_steps['dedup']
        log("3/7 Picard MarkDuplicates — Déjà terminé ✓")

    # ── ÉTAPE 4 : Clair3 ────────────────────────────────────────
    if 'clair3' not in completed_steps:
        clair3_vcf = run_clair3(
            dedup_bam, clair3_out_dir, sample_name,
            ctg_name=ctg_name, bed_file=bed_file
        )
    else:
        clair3_vcf = completed_steps['clair3']
        log("4/7 Clair3 — Déjà terminé ✓")

    # ── ÉTAPE 5 : ANNOVAR ───────────────────────────────────────
    if 'annovar' not in completed_steps:
        annovar_txt = run_annovar(clair3_vcf, annovar_prefix)
    else:
        annovar_txt = completed_steps['annovar']
        log("5/7 ANNOVAR — Déjà terminé ✓")

    # ── ÉTAPE 6 : Excel ─────────────────────────────────────────
    excel_file = convert_to_excel(annovar_txt, excel_file)

    # ── ÉTAPE 7 : IGV ───────────────────────────────────────────
    run_igv_snapshot(dedup_bam, clair3_vcf, sample_dir, sample_name)

    # ── NETTOYAGE ───────────────────────────────────────────────
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

# -----------------------------
# POINT D'ENTRÉE
# -----------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline Nanopore — NanoVar")

    parser.add_argument(
        "--sample", required=True,
        help="Nom de l'échantillon (ex: SAMPLE_01)"
    )
    parser.add_argument(
        "--fastq", required=True,
        help="Chemin vers le dossier FASTQ Nanopore"
    )
    parser.add_argument(
        "--ctg_name", required=False, default=None,
        help="Chromosome cible (ex: chr7) — gène unique"
    )
    parser.add_argument(
        "--bed", required=False, default=None,
        help="Chemin vers le fichier BED — panel de gènes"
    )

    args = parser.parse_args()

    validate_configuration()

    if not os.path.exists(args.fastq):
        log(f"❌ Dossier FASTQ non trouvé: {args.fastq}")
        sys.exit(1)

    try:
        result = process_sample(
            args.sample,
            args.fastq,
            ctg_name=args.ctg_name,
            bed_file=args.bed
        )
        log(f"🎉 Succès: {result['sample']}")
        sys.exit(0)
    except Exception as e:
        log(f"💥 ERREUR CRITIQUE: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
