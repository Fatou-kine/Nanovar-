import os
import json

from config.settings import BASE_DIR


GENES_PATH = os.path.join(
    BASE_DIR,
    "genes_hg38.json"
)


def load_genes():
    """Charge la base de données des gènes."""

    if not os.path.exists(GENES_PATH):
        return []

    with open(
        GENES_PATH,
        "r",
        encoding="utf-8"
    ) as f:
        data = json.load(f)

    return data.get("genes", [])


GENES_DB = load_genes()


def generate_bed_from_gene(
    gene_name,
    sample_dir,
    margin=500
):
    """
    Génère un fichier BED à partir du nom d'un gène.

    Une marge de 500 bp est ajoutée de chaque côté.
    """

    gene = next(
        (
            g for g in GENES_DB
            if g["name"].upper() == gene_name.upper()
        ),
        None
    )

    if not gene:
        return None, None

    start = max(
        0,
        gene["start"] - margin
    )

    end = gene["end"] + margin

    bed_content = (
        f"{gene['chr']}\t"
        f"{start}\t"
        f"{end}\t"
        f"{gene['name']}\n"
    )

    os.makedirs(
        sample_dir,
        exist_ok=True
    )

    bed_path = os.path.join(
        sample_dir,
        f"{gene_name}.bed"
    )

    with open(
        bed_path,
        "w",
        encoding="utf-8"
    ) as f:
        f.write(bed_content)

    return bed_path, gene["chr"]
