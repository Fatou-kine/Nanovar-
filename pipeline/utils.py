"""
pipeline/utils.py
Fonctions utilitaires partagées par tous les modules du pipeline NanoVar.
Chargement de la configuration, logging, exécution de commandes.
"""

import subprocess
import os
import sys
import time
import glob
import yaml


# -----------------------------
# CHARGEMENT DE LA CONFIGURATION
# -----------------------------

def load_config():
    """Charge config.yaml depuis la racine du projet (dossier parent de pipeline/)"""
    # Remonte d'un niveau depuis pipeline/ vers la racine du projet
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_path = os.path.join(root_dir, "config.yaml")
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"❌ config.yaml non trouvé: {config_path}")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


# Chargement unique au démarrage — importé par tous les modules
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
# LOGGING
# -----------------------------

def log(message):
    """Affiche un message horodaté — lu par Flask pour la progression"""
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


# -----------------------------
# EXÉCUTION DE COMMANDES
# -----------------------------

def run_command(cmd, error_message, timeout=14400):
    """Exécute une commande shell avec gestion d'erreur et timeout.

    Args:
        cmd: str ou list — la commande à exécuter
        error_message: str — message d'erreur en cas d'échec
        timeout: int — timeout en secondes (défaut 4h)

    Returns:
        subprocess.CompletedProcess

    Raises:
        RuntimeError si la commande échoue ou timeout
    """
    full_command = ' '.join(cmd) if isinstance(cmd, list) else cmd
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

        stderr_output = result.stderr.decode('utf-8', errors='replace') if result.stderr else ""
        stdout_output = result.stdout.decode('utf-8', errors='replace') if result.stdout else ""
        elapsed = time.time() - start_time

        log(f"✅ Terminé en {elapsed:.2f}s")

        if result.returncode != 0:
            log(f"🔍 STDERR: {stderr_output}")
            log(f"🔍 STDOUT: {stdout_output}")
            raise RuntimeError(f"{error_message}: {stderr_output}")

        return result

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Timeout après {timeout}s: {error_message}")


# -----------------------------
# VÉRIFICATION DE FICHIERS
# -----------------------------

def check_file_exists(file_path, description):
    """Vérifie qu'un fichier existe et n'est pas vide.

    Raises:
        FileNotFoundError si le fichier n'existe pas
        ValueError si le fichier est vide
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"❌ {description} non trouvé: {file_path}")
    if os.path.getsize(file_path) == 0:
        raise ValueError(f"❌ {description} est vide: {file_path}")
    log(f"✅ {description} vérifié")
    return True


# -----------------------------
# GESTION DES ÉTAPES
# -----------------------------

def check_existing_outputs(sample_dir, sample_name):
    """Vérifie quelles étapes sont déjà terminées pour reprendre en cas d'interruption.

    Returns:
        dict — clés = nom de l'étape, valeurs = chemin du fichier de sortie
    """
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
    """Supprime les fichiers intermédiaires après succès du pipeline."""
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


# -----------------------------
# VALIDATION DE LA CONFIGURATION
# -----------------------------

def validate_configuration():
    """Valide tous les chemins et ressources au démarrage.

    Raises:
        FileNotFoundError si un chemin obligatoire est manquant
        RuntimeError si Docker n'est pas disponible
    """
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

    # Vérifier que Docker est disponible
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
