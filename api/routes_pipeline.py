import os
import time
import signal

from flask import Blueprint, request, jsonify

from config.settings import OUTPUT_DIR
from core.jobs import (
    jobs,
    jobs_lock,
    executor,
    get_log_file_path,
)
from core.pipeline_runner import run_pipeline_job


pipeline_bp = Blueprint(
    "pipeline",
    __name__
)


@pipeline_bp.route("/api/lancer", methods=["POST"])
def lancer_pipeline():
    """
    Lance le pipeline pour un ou plusieurs échantillons.

    Body JSON attendu :

    {
        "echantillons": [
            {
                "sample_name": "SAMPLE_01",
                "fastq_dir": "/chemin/vers/fastq",
                "mode": "gene",
                "gene_name": "BRCA1",
                "bed_file": ""
            }
        ]
    }
    """

    data = request.get_json()

    if not data or "echantillons" not in data:
        return jsonify({
            "erreur": "Données manquantes"
        }), 400

    echantillons = data["echantillons"]

    if not isinstance(echantillons, list):
        return jsonify({
            "erreur": "echantillons doit être une liste"
        }), 400

    if len(echantillons) > 5:
        return jsonify({
            "erreur": "Maximum 5 échantillons simultanés"
        }), 400

    launched = []
    errors = []

    for ech in echantillons:

        sample_name = ech.get(
            "sample_name",
            ""
        ).strip()

        fastq_dir = ech.get(
            "fastq_dir",
            ""
        ).strip()

        mode = ech.get(
            "mode",
            "gene"
        )

        gene_name = ech.get(
            "gene_name",
            ""
        ).strip()

        bed_file = ech.get(
            "bed_file",
            ""
        ).strip()

        # -----------------------------
        # VALIDATION
        # -----------------------------

        if not sample_name:
            errors.append(
                "Nom d'échantillon manquant"
            )
            continue

        if not fastq_dir:
            errors.append(
                f"{sample_name}: dossier FASTQ manquant"
            )
            continue

        if not os.path.exists(fastq_dir):
            errors.append(
                f"{sample_name}: Dossier FASTQ non trouvé: {fastq_dir}"
            )
            continue

        # -----------------------------
        # INITIALISATION DU JOB
        # -----------------------------

        with jobs_lock:

            if (
                sample_name in jobs
                and jobs[sample_name].get("statut") == "en_cours"
            ):
                errors.append(
                    f"{sample_name}: Déjà en cours d'analyse"
                )
                continue

            jobs[sample_name] = {
                "sample_name": sample_name,
                "fastq_dir": fastq_dir,
                "mode": mode,
                "gene_name": gene_name,
                "statut": "en_attente",
                "etape": "En attente...",
                "progression": 0,
                "debut": None,
                "fin": None,
                "logs": [],
                "resultats": {},
                "pid": None,
            }

        # -----------------------------
        # SOUMISSION AU POOL
        # -----------------------------

        executor.submit(
            run_pipeline_job,
            sample_name,
            fastq_dir,
            mode,
            gene_name if gene_name else None,
            bed_file if bed_file else None,
        )

        launched.append(sample_name)

    return jsonify({
        "lances": launched,
        "erreurs": errors,
        "message": f"{len(launched)} échantillon(s) soumis",
    })


@pipeline_bp.route("/api/statut", methods=["GET"])
def get_statut():
    """
    Retourne l'état de tous les jobs.
    """

    with jobs_lock:

        statuts = {}

        for sample_name, job in jobs.items():

            statuts[sample_name] = {
                "sample_name": job["sample_name"],
                "statut": job["statut"],
                "etape": job["etape"],
                "progression": job["progression"],
                "debut": job["debut"],
                "fin": job["fin"],
                "gene_name": job.get(
                    "gene_name",
                    ""
                ),
                "mode": job.get(
                    "mode",
                    ""
                ),
                "resultats": job.get(
                    "resultats",
                    {}
                ),
            }

    return jsonify(statuts)


@pipeline_bp.route(
    "/api/logs/<sample_name>",
    methods=["GET"]
)
def get_logs(sample_name):
    """
    Retourne les logs d'un échantillon.

    Priorité :
    1. mémoire
    2. fichier pipeline.log
    """

    with jobs_lock:

        if sample_name in jobs:

            logs_list = jobs[
                sample_name
            ].get(
                "logs",
                []
            ).copy()

            return jsonify({
                "logs": logs_list,
                "source": "memoire",
            })

    # -----------------------------
    # REPLI SUR LE DISQUE
    # -----------------------------

    log_path = get_log_file_path(
        sample_name
    )

    if os.path.exists(log_path):

        try:

            with open(
                log_path,
                "r",
                encoding="utf-8"
            ) as f:

                lines = [
                    line.rstrip("\n")
                    for line in f.readlines()
                ]

            return jsonify({
                "logs": lines[-300:],
                "source": "disque",
            })

        except Exception as e:

            return jsonify({
                "erreur": f"Erreur lecture logs: {str(e)}"
            }), 500

    return jsonify({
        "erreur": "Échantillon non trouvé"
    }), 404


@pipeline_bp.route(
    "/api/annuler/<sample_name>",
    methods=["POST"]
)
def annuler_job(sample_name):
    """
    Annule un job en cours.
    """

    with jobs_lock:

        if sample_name not in jobs:

            return jsonify({
                "erreur": "Échantillon non trouvé"
            }), 404

        job = jobs[sample_name]

        if job["statut"] != "en_cours":

            return jsonify({
                "erreur": "Ce job n'est pas en cours"
            }), 400

        pid = job.get("pid")

    # -----------------------------
    # TERMINER LE PROCESSUS
    # -----------------------------

    if not pid:

        return jsonify({
            "erreur": "Impossible d'annuler : PID introuvable"
        }), 500

    try:

        os.kill(
            pid,
            signal.SIGTERM
        )

        from core.jobs import update_job

        update_job(
            sample_name,
            statut="annule",
            etape="Annulé par l'utilisateur",
            fin=time.strftime("%H:%M:%S"),
        )

        return jsonify({
            "message": f"{sample_name} annulé"
        })

    except ProcessLookupError:

        from core.jobs import update_job

        update_job(
            sample_name,
            statut="annule",
            etape="Annulé",
            fin=time.strftime("%H:%M:%S"),
        )

        return jsonify({
            "message": f"{sample_name} annulé"
        })

    except Exception as e:

        return jsonify({
            "erreur": f"Impossible d'annuler le job: {str(e)}"
        }), 500
