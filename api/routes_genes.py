from flask import Blueprint, request, jsonify

from core.genes import GENES_DB


genes_bp = Blueprint(
    "genes",
    __name__
)


@genes_bp.route("/api/genes", methods=["GET"])
def get_genes():
    """
    Retourne la liste des gènes pour l'autocomplétion.

    Exemple :
        /api/genes?q=BRCA
    """

    query = request.args.get(
        "q",
        ""
    ).upper().strip()

    if not query or len(query) < 2:
        return jsonify([])

    matching = [
        {
            "name": g["name"],
            "chr": g["chr"],
            "description": g["description"]
        }
        for g in GENES_DB
        if query in g["name"].upper()
    ]

    return jsonify(matching[:10])


@genes_bp.route("/api/chromosomes", methods=["GET"])
def get_chromosomes():
    """Retourne la liste des chromosomes disponibles."""

    chromosomes = (
        [f"chr{i}" for i in range(1, 23)]
        + ["chrX", "chrY", "chrM"]
    )

    return jsonify(chromosomes)
