import os

from config.settings import OUTPUT_DIR


def get_sample_results(sample_name):
    """Retourne les fichiers de résultats disponibles."""

    sample_dir = os.path.join(
        OUTPUT_DIR,
        sample_name
    )

    results = {}

    files = {
        "excel": f"{sample_name}_variants.xlsx",
        "vcf_txt": (
            f"{sample_name}_annovar.hg38_multianno.txt"
        ),
        "vcf": (
            f"{sample_name}_annovar.hg38_multianno.vcf"
        ),
        "bam": f"{sample_name}_dedup.bam",
        "metrics": f"{sample_name}_markdup_metrics.txt",
        "vcf_gz": os.path.join(
            "clair3_output",
            "merge_output.vcf.gz"
        ),
    }

    for key, filename in files.items():

        filepath = os.path.join(
            sample_dir,
            filename
        )

        if os.path.exists(filepath):
            results[key] = filepath

    return results
