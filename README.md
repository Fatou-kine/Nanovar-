# NanoVar 🧬

**Plateforme d'analyse de données de séquençage Nanopore**

NanoVar est une plateforme web locale permettant l'analyse complète de données de séquençage Oxford Nanopore (ONT) pour des applications cliniques et de recherche sur gène unique ou panel de gènes.

---

## Fonctionnalités

- 📁 **Import simplifié** — l'utilisateur fournit uniquement le dossier de reads Nanopore
- 🔄 **Pipeline automatisé** — concaténation, alignement, déduplication, appel de variants, annotation
- 🧬 **Sélection de gène par nom** — base de données de 120+ gènes cliniquement pertinents (hg38)
- ⚡ **Traitement parallèle** — jusqu'à 5 échantillons simultanément
- 📊 **Export multiple** — VCF annoté, Excel, BAM, rapport PDF
- 🖥️ **Interface web locale** — accessible via navigateur, aucune installation complexe

---

## Pipeline

```
Dossier FASTQ (Nanopore)
        ↓
1. Concaténation des reads         cat
        ↓
2. Alignement                      Minimap2 (map-ont)
        ↓
3. Conversion SAM→BAM + tri        Samtools
        ↓
4. Marquage des duplicats          Picard MarkDuplicates
        ↓
5. Appel de variants               Clair3 (r1041_e82_400bps_sup_v430)
        ↓
6. Annotation                      ANNOVAR (hg38)
        ↓
7. Export                          VCF + Excel + BAM
```

---

## Prérequis

### Outils bioinformatiques
- [Minimap2](https://github.com/lh3/minimap2)
- [Samtools](http://www.htslib.org/)
- [Picard](https://broadinstitute.github.io/picard/)
- [ANNOVAR](https://annovar.openbioinformatics.org/)
- [Docker](https://www.docker.com/) — pour Clair3

### Image Docker Clair3
```bash
docker pull hkubal/clair3:v1.0.10
```

### Modèle Clair3
```bash
# Télécharger le modèle r1041_e82_400bps_sup_v430 (P2 Solo)
git clone https://github.com/nanoporetech/rerio.git
cd rerio
python3 download_model.py --clair3 clair3_models/r1041_e82_400bps_sup_v430_model
```

### Python
```bash
pip install -r requirements.txt
```

---

## Installation

```bash
# Cloner le dépôt
git clone https://github.com/Fatou-kine/NanoVar.git
cd NanoVar

# Installer les dépendances Python
pip install -r requirements.txt

# Configurer les chemins
cp config.example.yaml config.yaml
# Éditer config.yaml avec vos chemins locaux
```

---

## Structure du projet

```
NanoVar/
├── config.example.yaml       # Template de configuration
├── config.yaml               # Configuration locale (non versionné)
├── NANOPORE_PIPELINE.py      # Pipeline principal
├── app.py                    # Backend Flask
├── genes_hg38.json           # Base de données des gènes (120+ gènes)
├── requirements.txt          # Dépendances Python
├── templates/
│   └── index.html            # Interface utilisateur
├── static/
│   ├── app.js                # JavaScript frontend
│   └── style.css             # Styles
├── tools/
│   └── picard.jar            # Picard (non versionné)
└── annovar/                  # ANNOVAR (non versionné)
    └── humandb/
```

---

## Utilisation

### Lancer la plateforme
```bash
cd /chemin/vers/NanoVar
python app.py
```
Ouvrir le navigateur à l'adresse : **http://localhost:5000**

### Utilisation en ligne de commande (sans interface)
```bash
python NANOPORE_PIPELINE.py \
  --sample NOM_ECHANTILLON \
  --fastq /chemin/vers/dossier/reads
```

---

## Données d'entrée

NanoVar accepte les données issues du basecalling **Guppy** ou **Dorado** en mode SUP (Super Accuracy). Le dossier d'entrée doit contenir les fichiers `.fastq.gz` ou `.fastq` du barcode correspondant à l'échantillon.

```
run_nanopore/
└── barcode01/
    ├── FAT12345_pass_barcode01_0.fastq.gz
    ├── FAT12345_pass_barcode01_1.fastq.gz
    └── ...
```

---

## Données de sortie

Pour chaque échantillon, NanoVar génère dans le dossier `OUTPUT/` :

| Fichier | Description |
|---|---|
| `{sample}_dedup.bam` | BAM final aligné et dédupliqué |
| `clair3_output/merge_output.vcf.gz` | VCF brut Clair3 |
| `{sample}_annovar.hg38_multianno.vcf` | VCF annoté ANNOVAR |
| `{sample}_annovar.hg38_multianno.txt` | Tableau annoté (TSV) |
| `{sample}_variants.xlsx` | Variants en format Excel |
| `{sample}_markdup_metrics.txt` | Métriques Picard |

---

## Gènes supportés

NanoVar inclut une base de données de **120+ gènes** couvrant :
- Cancers héréditaires (BRCA1/2, MLH1, MSH2, TP53, APC...)
- Cardiomyopathies (MYH7, MYBPC3, SCN5A, KCNQ1...)
- Maladies neurologiques (HTT, LRRK2, PSEN1, FMR1...)
- Maladies hématologiques (HBB, JAK2, FLT3, NPM1...)
- Maladies métaboliques (PAH, CFTR, GBA, HEXA...)
- Et bien d'autres...

---

## Notes importantes

- ⚠️ Le **BQSR** (BaseRecalibrator) est intentionnellement exclu — non adapté aux données Nanopore
- ⚠️ Le **BCFtools filter** est exclu — Clair3 intègre son propre système de scoring
- ✅ Modèle Clair3 validé : `r1041_e82_400bps_sup_v430` (P2 Solo, R10.4.1, SUP 400bps)

---

## Auteur

**Fatou-kine** — [@Fatou-kine](https://github.com/Fatou-kine)

---

## Licence

Ce projet est à usage académique et de recherche.
