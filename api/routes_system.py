import os
import subprocess

from flask import Blueprint, jsonify

from config.settings import (
    CONFIG,
    MAX_PARALLEL,
    OUTPUT_DIR,
)

from core.genes import GENES_DB


system_bp = Blueprint(
    "system",
    __name__
)


@system_bp.route(
    "/api/config",
    methods=["GET"]
)
def get_config_info():
    """
    Retourne les informations de configuration
    sans exposer de données sensibles.
    """

    return jsonify({
        "max_parallel": MAX_PARALLEL,
        "output_dir": OUTPUT_DIR,
        "genes_count": len(GENES_DB),
        "clair3_image": CONFIG.get(
            "clair3_docker_image",
            ""
        ),
        "clair3_model": os.path.basename(
            CONFIG.get(
                "clair3_model_local",
                ""
            )
        ),
    })


@system_bp.route(
    "/api/verifier_docker",
    methods=["GET"]
)
def verifier_docker():
    """
    Vérifie que Docker est disponible
    et que l'image Clair3 est présente.
    """

    image = CONFIG.get(
        "clair3_docker_image",
        ""
    )

    try:

        docker = subprocess.run(
            ["docker", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )

        if docker.returncode != 0:

            return jsonify({
                "docker": False,
                "image": image,
                "message": (
                    "Docker n'est pas disponible"
                )
            })

        result = subprocess.run(
            [
                "docker",
                "image",
                "inspect",
                image
            ],
            capture_output=True,
            text=True,
            timeout=10
        )

        image_present = (
            result.returncode == 0
        )

        return jsonify({
            "docker": True,
            "image": image,
            "image_present": image_present,
            "docker_version": docker.stdout.strip(),
            "message": (
                "Docker et l'image Clair3 "
                "sont disponibles"
                if image_present
                else
                "Docker est disponible mais "
                "l'image Clair3 est absente"
            )
        })

    except FileNotFoundError:

        return jsonify({
            "docker": False,
            "image": image,
            "message": (
                "Commande Docker introuvable"
            )
        })

    except Exception as e:

        return jsonify({
            "docker": False,
            "image": image,
            "message": str(e)
        })
