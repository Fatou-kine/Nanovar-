import os
import sys
import time
import json
import yaml
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from flask import Flask, request, jsonify, render_template, send_file
from flask_cors import CORS

# -----------------------------
# INITIALISATION FLASK
# -----------------------------

app = Flask(__name__)
CORS(app)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# -----------------------------
# CHARGEMENT CONFIGURATION
# -----------------------------

def load_config():
    config_path = os.path.join(BASE_DIR, "config.yaml")
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"config.yaml non trouvé: {config_path}")
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

CONFIG = load_config()
MAX_PARALLEL = int(CONFIG.get('max_parallel', 2))
OUTPUT_DIR   = CONFIG['output_dir']

# -----------------------------
# CHARGEMENT BASE DE DONNÉES GÈNES
# -----------------------------

def load_genes():
    genes_path = os.path.join(BASE_DIR, "genes_hg38.json")
    if not os.path.exists(genes_path):
        return []
    with open(genes_path, 'r') as f:
        data = json.load(f)
    return data.get('genes', [])

GENES_DB = load_genes()

# -----------------------------
# JOB STORE — MÉMOIRE PARTAGÉE
# -----------------------------
# Dictionnaire global qui stocke l'état de chaque job en temps réel
# Clé = sample_name, Valeur = dict avec statut, étape, progression, logs

jobs = {}
jobs_lock = threading.Lock()

# Pool de threads limité à MAX_PARALLEL échantillons simultanés
executor = ThreadPoolExecutor(max_workers=MAX_PARALLEL)

# -----------------------------
# FONCTIONS DE GESTION DES JOBS
# -----------------------------

def update_job(sample_name, **kwargs):
    """Met à jour l'état d'un job dans le job store"""
    with jobs_lock:
        if sample_name in jobs:
            jobs[sample_name].update(kwargs)

def add_log(sample_name, message):
    """Ajoute une ligne de log à un job"""
    with jobs_lock:
        if sample_name in jobs:
            jobs[sample_name]['logs'].append(message)
            # Garder uniquement les 200 dernières lignes
            if len(jobs[sample_name]['logs']) > 200:
                jobs[sample_name]['logs'] = jobs[sample_name]['logs'][-200:]

def parse_progression(line):
    """Extrait la progression depuis une ligne de log du pipeline
    Le pipeline émet des messages [ETAPE X/7] que Flask intercepte
    """
    etapes = {
        "ETAPE 1/7": (1,  "Concaténation FASTQ"),
        "ETAPE 2/7": (2,  "Alignement Minimap2"),
        "ETAPE 3/7": (3,  "Picard MarkDuplicates"),
        "ETAPE 4/7": (4,  "Clair3 - Appel de variants"),
        "ETAPE 5/7": (5,  "Annotation ANNOVAR"),
        "ETAPE 6/7": (6,  "Conversion Excel"),
        "ETAPE 7/7": (7,  "IGV Snapshots"),
    }
    for key, (num, label) in etapes.items():
        if key in line:
            progression = int((num / 7) * 100)
            return progression, label
    return None, None

def generate_bed_from_gene(gene_name, sample_dir, margin=500):
    """Génère un fichier BED à partir du nom d'un gène
    Ajoute une marge de 500bp de chaque côté pour ne pas rater
    les variants aux extrémités de la région
    """
    gene = next((g for g in GENES_DB if g['name'].upper() == gene_name.upper()), None)
    if not gene:
        return None, None

    # Appliquer la marge
    start = max(0, gene['start'] - margin)
    end   = gene['end'] + margin

    # Créer le fichier BED
    bed_content = f"{gene['chr']}\t{start}\t{end}\t{gene['name']}\n"
    bed_path = os.path.join(sample_dir, f"{gene_name}.bed")

    os.makedirs(sample_dir, exist_ok=True)
    with open(bed_path, 'w') as f:
        f.write(bed_content)

    return bed_path, gene['chr']

def run_pipeline_job(sample_name, fastq_dir, mode, gene_name=None, bed_file=None):
    """Fonction exécutée dans un thread séparé pour chaque échantillon
    Lance le pipeline Python via subprocess et lit les logs en temps réel
    """
    update_job(sample_name,
        statut="en_cours",
        debut=time.strftime("%H:%M:%S"),
        progression=0,
        etape="Initialisation..."
    )

    try:
        # Préparer le dossier de sortie de l'échantillon
        sample_output_dir = os.path.join(OUTPUT_DIR, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)

        # Construire la commande selon le mode
        cmd = [
            sys.executable,
            os.path.join(BASE_DIR, "NANOPORE_PIPELINE.py"),
            "--sample", sample_name,
            "--fastq",  fastq_dir,
        ]

        ctg_name = None

        if mode == "gene" and gene_name:
            # Générer le BED depuis le nom du gène
            bed_path, ctg_name = generate_bed_from_gene(gene_name, sample_output_dir)
            if bed_path:
                cmd += ["--bed", bed_path]
                cmd += ["--ctg_name", ctg_name]
                add_log(sample_name, f"[INFO] Gène: {gene_name} → {ctg_name}")
            else:
                add_log(sample_name, f"[WARNING] Gène {gene_name} non trouvé dans la base — analyse sur génome complet")

        elif mode == "panel" and bed_file:
            # BED fourni directement par l'utilisateur
            cmd += ["--bed", bed_file]
            add_log(sample_name, f"[INFO] Panel BED: {os.path.basename(bed_file)}")

        elif mode == "chromosome" and gene_name:
            # gene_name contient le chromosome (ex: chr7)
            cmd += ["--ctg_name", gene_name]
            add_log(sample_name, f"[INFO] Chromosome: {gene_name}")

        add_log(sample_name, f"[INFO] Lancement: {' '.join(cmd)}")

        # Lancer le pipeline et lire les logs en temps réel
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
            universal_newlines=True
        )

        # Stocker le PID pour permettre l'annulation
        update_job(sample_name, pid=process.pid)

        # Lire chaque ligne de sortie du pipeline
        for line in process.stdout:
            line = line.rstrip()
            if line:
                add_log(sample_name, line)

                # Mettre à jour la progression selon l'étape détectée
                progression, etape = parse_progression(line)
                if progression is not None:
                    update_job(sample_name,
                        progression=progression,
                        etape=etape
                    )

                # Détecter les erreurs critiques
                if "ERREUR CRITIQUE" in line or "💥" in line:
                    update_job(sample_name, statut="erreur", etape="Erreur critique")

        process.wait()

        if process.returncode == 0:
            # Pipeline terminé avec succès
            update_job(sample_name,
                statut="termine",
                progression=100,
                etape="Terminé ✓",
                fin=time.strftime("%H:%M:%S"),
                resultats=get_sample_results(sample_name)
            )
            add_log(sample_name, f"✅ Pipeline terminé avec succès pour {sample_name}")
        else:
            update_job(sample_name,
                statut="erreur",
                etape="Erreur — voir les logs",
                fin=time.strftime("%H:%M:%S")
            )
            add_log(sample_name, f"❌ Pipeline échoué pour {sample_name} (code: {process.returncode})")

    except Exception as e:
        update_job(sample_name,
            statut="erreur",
            etape=f"Erreur: {str(e)}",
            fin=time.strftime("%H:%M:%S")
        )
        add_log(sample_name, f"❌ Erreur: {str(e)}")

def get_sample_results(sample_name):
    """Retourne les chemins des fichiers de résultats d'un échantillon"""
    sample_dir = os.path.join(OUTPUT_DIR, sample_name)
    results = {}

    files = {
        "excel":   f"{sample_name}_variants.xlsx",
        "vcf_txt": f"{sample_name}_annovar.hg38_multianno.txt",
        "vcf":     f"{sample_name}_annovar.hg38_multianno.vcf",
        "bam":     f"{sample_name}_dedup.bam",
        "metrics": f"{sample_name}_markdup_metrics.txt",
        "vcf_gz":  os.path.join("clair3_output", "merge_output.vcf.gz"),
    }

    for key, filename in files.items():
        filepath = os.path.join(sample_dir, filename)
        if os.path.exists(filepath):
            results[key] = filepath

    return results

# -----------------------------
# ROUTES FLASK
# -----------------------------

@app.route('/')
def index():
    """Page principale — sert l'interface utilisateur"""
    return render_template('index.html')

@app.route('/api/genes', methods=['GET'])
def get_genes():
    """Retourne la liste des gènes pour l'autocomplétion
    Supporte la recherche par query ?q=BRCA
    """
    query = request.args.get('q', '').upper().strip()

    if not query or len(query) < 2:
        return jsonify([])

    # Filtrer les gènes qui correspondent à la recherche
    matching = [
        {
            "name":        g['name'],
            "chr":         g['chr'],
            "description": g['description']
        }
        for g in GENES_DB
        if query in g['name'].upper()
    ]

    # Limiter à 10 résultats
    return jsonify(matching[:10])

@app.route('/api/chromosomes', methods=['GET'])
def get_chromosomes():
    """Retourne la liste des chromosomes disponibles"""
    chromosomes = (
        [f"chr{i}" for i in range(1, 23)] +
        ["chrX", "chrY", "chrM"]
    )
    return jsonify(chromosomes)

@app.route('/api/lancer', methods=['POST'])
def lancer_pipeline():
    """Lance le pipeline pour un ou plusieurs échantillons
    Body JSON attendu:
    {
        "echantillons": [
            {
                "sample_name": "SAMPLE_01",
                "fastq_dir": "/chemin/vers/barcode01",
                "mode": "gene",          // "gene", "panel", "chromosome"
                "gene_name": "BRCA1",    // si mode=gene ou mode=chromosome
                "bed_file": "/chemin/vers/panel.bed"  // si mode=panel
            }
        ]
    }
    """
    data = request.get_json()

    if not data or 'echantillons' not in data:
        return jsonify({"erreur": "Données manquantes"}), 400

    echantillons = data['echantillons']

    if len(echantillons) > 5:
        return jsonify({"erreur": "Maximum 5 échantillons simultanés"}), 400

    launched = []
    errors   = []

    for ech in echantillons:
        sample_name = ech.get('sample_name', '').strip()
        fastq_dir   = ech.get('fastq_dir', '').strip()
        mode        = ech.get('mode', 'gene')
        gene_name   = ech.get('gene_name', '').strip()
        bed_file    = ech.get('bed_file', '').strip()

        # Validations
        if not sample_name:
            errors.append("Nom d'échantillon manquant")
            continue

        if not fastq_dir or not os.path.exists(fastq_dir):
            errors.append(f"{sample_name}: Dossier FASTQ non trouvé: {fastq_dir}")
            continue

        # Vérifier que l'échantillon n'est pas déjà en cours
        with jobs_lock:
            if sample_name in jobs and jobs[sample_name].get('statut') == 'en_cours':
                errors.append(f"{sample_name}: Déjà en cours d'analyse")
                continue

            # Initialiser le job dans le store
            jobs[sample_name] = {
                "sample_name": sample_name,
                "fastq_dir":   fastq_dir,
                "mode":        mode,
                "gene_name":   gene_name,
                "statut":      "en_attente",
                "etape":       "En attente...",
                "progression": 0,
                "debut":       None,
                "fin":         None,
                "logs":        [],
                "resultats":   {},
                "pid":         None
            }

        # Soumettre le job au pool de threads
        executor.submit(
            run_pipeline_job,
            sample_name,
            fastq_dir,
            mode,
            gene_name if gene_name else None,
            bed_file  if bed_file  else None
        )

        launched.append(sample_name)

    return jsonify({
        "lances":  launched,
        "erreurs": errors,
        "message": f"{len(launched)} échantillon(s) soumis"
    })

@app.route('/api/statut', methods=['GET'])
def get_statut():
    """Retourne l'état de tous les jobs en cours
    Appelé toutes les 3 secondes par le frontend (polling)
    """
    with jobs_lock:
        statuts = {}
        for sample_name, job in jobs.items():
            statuts[sample_name] = {
                "sample_name": job['sample_name'],
                "statut":      job['statut'],
                "etape":       job['etape'],
                "progression": job['progression'],
                "debut":       job['debut'],
                "fin":         job['fin'],
                "gene_name":   job.get('gene_name', ''),
                "mode":        job.get('mode', ''),
                "resultats":   job.get('resultats', {}),
            }
    return jsonify(statuts)

@app.route('/api/logs/<sample_name>', methods=['GET'])
def get_logs(sample_name):
    """Retourne les logs d'un échantillon spécifique"""
    with jobs_lock:
        if sample_name not in jobs:
            return jsonify({"erreur": "Échantillon non trouvé"}), 404
        logs_list = jobs[sample_name]['logs'].copy()

    return jsonify({"logs": logs_list})

@app.route('/api/annuler/<sample_name>', methods=['POST'])
def annuler_job(sample_name):
    """Annule un job en cours"""
    with jobs_lock:
        if sample_name not in jobs:
            return jsonify({"erreur": "Échantillon non trouvé"}), 404

        job = jobs[sample_name]

        if job['statut'] != 'en_cours':
            return jsonify({"erreur": "Ce job n'est pas en cours"}), 400

        pid = job.get('pid')

    # Tuer le processus si on a son PID
    if pid:
        try:
            import signal
            os.kill(pid, signal.SIGTERM)
            update_job(sample_name,
                statut="annule",
                etape="Annulé par l'utilisateur",
                fin=time.strftime("%H:%M:%S")
            )
            return jsonify({"message": f"{sample_name} annulé"})
        except ProcessLookupError:
            update_job(sample_name, statut="annule")
            return jsonify({"message": f"{sample_name} annulé"})

    return jsonify({"erreur": "Impossible d'annuler"}), 500

@app.route('/api/telecharger/<sample_name>/<file_type>', methods=['GET'])
def telecharger_fichier(sample_name, file_type):
    """Télécharge un fichier de résultat
    file_type: excel, vcf, vcf_txt, bam, metrics, vcf_gz
    """
    with jobs_lock:
        if sample_name not in jobs:
            return jsonify({"erreur": "Échantillon non trouvé"}), 404
        resultats = jobs[sample_name].get('resultats', {})

    if file_type not in resultats:
        # Chercher directement dans le dossier de sortie
        resultats = get_sample_results(sample_name)

    if file_type not in resultats:
        return jsonify({"erreur": f"Fichier {file_type} non disponible"}), 404

    file_path = resultats[file_type]

    if not os.path.exists(file_path):
        return jsonify({"erreur": "Fichier introuvable sur le disque"}), 404

    return send_file(
        file_path,
        as_attachment=True,
        download_name=os.path.basename(file_path)
    )

@app.route('/api/supprimer/<sample_name>', methods=['DELETE'])
def supprimer_job(sample_name):
    """Supprime un job du job store (ne supprime pas les fichiers)"""
    with jobs_lock:
        if sample_name not in jobs:
            return jsonify({"erreur": "Échantillon non trouvé"}), 404
        if jobs[sample_name]['statut'] == 'en_cours':
            return jsonify({"erreur": "Impossible de supprimer un job en cours"}), 400
        del jobs[sample_name]

    return jsonify({"message": f"{sample_name} supprimé du tableau de bord"})

@app.route('/api/config', methods=['GET'])
def get_config_info():
    """Retourne les informations de configuration (sans données sensibles)"""
    return jsonify({
        "max_parallel":   MAX_PARALLEL,
        "output_dir":     OUTPUT_DIR,
        "genes_count":    len(GENES_DB),
        "clair3_image":   CONFIG.get('clair3_docker_image', ''),
        "clair3_model":   os.path.basename(CONFIG.get('clair3_model_local', '')),
    })

@app.route('/api/verifier_docker', methods=['GET'])
def verifier_docker():
    """Vérifie que Docker est disponible et l'image Clair3 est présente"""
    try:
        result = subprocess.run(
            "docker info",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10
        )
        if result.returncode != 0:
            return jsonify({"docker": False, "message": "Docker non disponible"}), 500

        # Vérifier que l'image Clair3 est présente
        image = CONFIG.get('clair3_docker_image', '')
        result2 = subprocess.run(
            f"docker image inspect {image}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10
        )
        image_ok = result2.returncode == 0

        return jsonify({
            "docker":   True,
            "image":    image,
            "image_ok": image_ok,
            "message":  "Docker disponible" + (f" — image {image} présente" if image_ok else f" — image {image} non trouvée")
        })

    except Exception as e:
        return jsonify({"docker": False, "message": str(e)}), 500

# -----------------------------
# DÉMARRAGE DU SERVEUR
# -----------------------------

if __name__ == '__main__':
    print("=" * 55)
    print("🧬 NanoVar — Plateforme d'analyse Nanopore")
    print("=" * 55)
    print(f"📁 Output         : {OUTPUT_DIR}")
    print(f"⚡ Max parallèle  : {MAX_PARALLEL} échantillons")
    print(f"🧬 Gènes chargés  : {len(GENES_DB)}")
    print(f"🐳 Clair3 image   : {CONFIG.get('clair3_docker_image', '')}")
    print("=" * 55)
    print("🌐 Ouvrir dans le navigateur : http://localhost:5000")
    print("=" * 55)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    app.run(
        host='0.0.0.0',
        port=5000,
        debug=False,
        threaded=True
    )
