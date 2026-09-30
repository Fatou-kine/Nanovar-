import os
import threading
from concurrent.futures import ThreadPoolExecutor

from config.settings import OUTPUT_DIR, MAX_PARALLEL


# Stockage des jobs en mémoire
jobs = {}

# Protection du dictionnaire partagé
jobs_lock = threading.Lock()

# Pool de threads
executor = ThreadPoolExecutor(
    max_workers=MAX_PARALLEL
)


def update_job(sample_name, **kwargs):
    """Met à jour l'état d'un job."""

    with jobs_lock:
        if sample_name in jobs:
            jobs[sample_name].update(kwargs)


def get_log_file_path(sample_name):
    """Retourne le chemin du fichier de log."""

    return os.path.join(
        OUTPUT_DIR,
        sample_name,
        "pipeline.log"
    )


def add_log(sample_name, message):
    """
    Ajoute un log en mémoire et sur disque.
    """

    with jobs_lock:

        if sample_name in jobs:

            jobs[sample_name]["logs"].append(
                message
            )

            # Garder les 300 dernières lignes
            if len(jobs[sample_name]["logs"]) > 300:
                jobs[sample_name]["logs"] = (
                    jobs[sample_name]["logs"][-300:]
                )

    # Persistance sur disque
    try:

        log_path = get_log_file_path(
            sample_name
        )

        os.makedirs(
            os.path.dirname(log_path),
            exist_ok=True
        )

        with open(
            log_path,
            "a",
            encoding="utf-8"
        ) as f:
            f.write(message + "\n")

    except Exception:
        pass
