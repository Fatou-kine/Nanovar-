"""
pipeline/export.py
ETAPE 6/7 — Conversion du tableau ANNOVAR en Excel (.xlsx)
ETAPE 7/7 — Snapshots IGV (optionnel, non bloquant)
"""
import os
from pipeline.utils import log, run_command, IGV_PATH


def convert_to_excel(annovar_txt, excel_file):
    """Convertit le fichier hg38_multianno.txt en Excel avec colonnes ajustées."""
    log("[ETAPE 6/7] Conversion Excel")
    if not os.path.exists(annovar_txt) or os.path.getsize(annovar_txt) == 0:
        log("⚠️  Fichier ANNOVAR vide — export Excel ignoré")
        return None
    try:
        import pandas as pd
        df = pd.read_csv(annovar_txt, sep='\t', low_memory=False)
        if df.empty:
            log("⚠️  Aucun variant à exporter")
            return None
        log(f"   → {len(df)} variant(s) à exporter")
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df.to_excel(writer, index=False, sheet_name='Variants')
            ws = writer.sheets['Variants']
            for col in ws.columns:
                max_len = max(len(str(c.value)) if c.value else 0 for c in col)
                ws.column_dimensions[col[0].column_letter].width = min(max_len + 2, 50)
            ws.freeze_panes = "A2"
        log(f"✅ Excel généré → {os.path.basename(excel_file)}")
        return excel_file
    except ImportError:
        log("⚠️  pandas/openpyxl manquant — pip install pandas openpyxl")
        return None
    except Exception as e:
        log(f"⚠️  Erreur Excel: {e}")
        return None


def run_igv_snapshot(dedup_bam, vcf_file, sample_dir, sample_name):
    """Génère des snapshots IGV automatiques (optionnel, non bloquant)."""
    log("[ETAPE 7/7] IGV snapshots")
    if not IGV_PATH or not os.path.exists(IGV_PATH):
        log("⚠️  IGV non configuré — étape ignorée")
        return None

    snapshot_dir = os.path.join(sample_dir, "igv_snapshots")
    os.makedirs(snapshot_dir, exist_ok=True)

    igv_batch = os.path.join(sample_dir, f"{sample_name}_igv_batch.txt")
    with open(igv_batch, 'w') as f:
        f.write("new\n")
        f.write("genome hg38\n")
        f.write(f"load {dedup_bam}\n")
        f.write(f"load {vcf_file}\n")
        f.write(f"snapshotDirectory {snapshot_dir}\n")
        f.write(f"snapshot {sample_name}_overview.png\n")
        f.write("exit\n")

    try:
        run_command(f"bash \"{IGV_PATH}\" --batch \"{igv_batch}\"",
                    "Erreur IGV", timeout=3600)
        log(f"✅ IGV snapshots → {snapshot_dir}")
    except Exception as e:
        log(f"⚠️  IGV échoué (non bloquant): {e}")
    return snapshot_dir
