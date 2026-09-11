import os
import yaml


# Racine du projet NanoVar
BASE_DIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)


def load_config():
    """Charge le fichier config.yaml."""

    config_path = os.path.join(
        BASE_DIR,
        "config.yaml"
    )

    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"config.yaml non trouvé: {config_path}"
        )

    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


CONFIG = load_config()

MAX_PARALLEL = int(
    CONFIG.get("max_parallel", 2)
)

OUTPUT_DIR = CONFIG["output_dir"]
