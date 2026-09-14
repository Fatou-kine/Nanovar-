import os

from flask import Blueprint, jsonify, send_file

from core.jobs import jobs, jobs_lock
from core.results import get_sample_results


files_bp = Blueprint(
    "files",
    __name__
)


@files_bp.route(
    "/api/telecharger/<sample_name>/<file_type>",
    methods=["GET"]
)
def telecharger_fichier(sample_name, file_type):
    """
    Télécharge un fichier de résultat.

    Types disponibles :
    excel
    vcf
    vcf_txt
    bam
    metrics
    vcf_gz
    """

    with jobs_lock:

        if sample_name not in jobs:

            # Le job peut ne plus être en mémoire
            # mais les résultats peuvent toujours exister
            resultats = get_sample_results(sample_name)

        else:

            resultats = jobs[
                sample_name
            ].get(
                "resultats",
                {}
            )

    # Si le résultat n'est pas présent en mémoire,
    # chercher directement sur disque
    if file_type not in resultats:

        resultats = get_sample_results(
            sample_name
        )

    if file_type not in resultats:

        return jsonify({
            "erreur": (
                f"Fichier {file_type} "
                "non disponible"
            )
        }), 404

    file_path = resultats[file_type]

    if not os.path.exists(file_path):

        return jsonify({
            "erreur": "Fichier introuvable sur le disque"
        }), 404

    return send_file(
        file_path,
        as_attachment=True,
        download_name=os.path.basename(file_path)
    )


@files_bp.route(
    "/api/supprimer/<sample_name>",
    methods=["DELETE"]
)
def supprimer_job(sample_name):
    """
    Supprime un job du job store.

    Les fichiers de résultats ne sont PAS supprimés.
    """

    with jobs_lock:

        if sample_name not in jobs:

            return jsonify({
                "erreur": "Échantillon non trouvé"
            }), 404

        if jobs[
            sample_name
        ].get("statut") == "en_cours":

            return jsonify({
                "erreur": (
                    "Impossible de supprimer "
                    "un job en cours"
                )
            }), 400

        del jobs[sample_name]

    return jsonify({
        "message": (
            f"{sample_name} supprimé "
            "du tableau de bord"
        )
    })
