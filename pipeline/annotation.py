"""
pipeline/annotation.py
ETAPE 5/7 — Annotation des variants avec ANNOVAR

Bases: refGene, ClinVar 2024, gnomAD genome, dbNSFP 4.7a, avsnp151

On passe par convert2annovar.pl séparément (avec --minqual 0) avant
table_annovar.pl pour éviter l'erreur 'last column should start with
Otherinfo' avec les VCF produits par Clair3.
"""
import os
from pipeline.utils import log, run_command, check_file_exists, ANNOVAR_DIR, HUMANDB


def run_annovar(vcf_file, annovar_prefix):
    """Annote les variants avec ANNOVAR en deux étapes.

    Étape 5a: convert2annovar.pl → VCF vers .avinput (--minqual 0)
    Étape 5b: table_annovar.pl   → annotation avec les bases hg38
    """
    log("[ETAPE 5/7] ANNOVAR annotation")
    check_file_exists(vcf_file, "Fichier VCF Clair3")

    avinput = f"{annovar_prefix}.avinput"

    # ── Étape 5a : Conversion VCF → avinput ──────────────────────
    log("   → Conversion VCF → avinput (--minqual 0)...")
    run_command(
        f"perl \"{os.path.join(ANNOVAR_DIR, 'convert2annovar.pl')}\" "
        f"-format vcf4 \"{vcf_file}\" "
        f"-includeinfo -withfreq -allsample -minqual 0 "
        f"> \"{avinput}\"",
        "Erreur convert2annovar", timeout=3600
    )

    if not os.path.exists(avinput) or os.path.getsize(avinput) == 0:
        log("⚠️  Aucun variant dans l'avinput — VCF Clair3 vide")
        _create_empty_outputs(annovar_prefix)
        return f"{annovar_prefix}.hg38_multianno.txt"

    n = sum(1 for _ in open(avinput) if _.strip())
    log(f"   → {n} variant(s) à annoter")

    # ── Étape 5b : Annotation ──────────────────────────────────────
    log("   → Annotation avec refGene, ClinVar, gnomAD, dbNSFP, avsnp151...")
    run_command(
        f"perl \"{os.path.join(ANNOVAR_DIR, 'table_annovar.pl')}\" "
        f"\"{avinput}\" \"{HUMANDB}\" "
        f"-buildver hg38 -out \"{annovar_prefix}\" "
        f"-protocol refGene,clinvar_20240611,gnomad_genome,dbnsfp47a,avsnp151 "
        f"-operation g,f,f,f,f -nastring . -otherinfo",
        "Erreur table_annovar", timeout=14400
    )

    annovar_txt = f"{annovar_prefix}.hg38_multianno.txt"
    check_file_exists(annovar_txt, "Fichier TXT ANNOVAR")
    log(f"✅ ANNOVAR terminé → {os.path.basename(annovar_txt)}")
    return annovar_txt


def _create_empty_outputs(annovar_prefix):
    for suffix in ['.hg38_multianno.txt', '.hg38_multianno.vcf']:
        path = f"{annovar_prefix}{suffix}"
        if not os.path.exists(path):
            open(path, 'w').close()
    log("   → Fichiers vides créés (pipeline non bloqué)")
