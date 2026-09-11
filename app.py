import os
import threading
import time
import webbrowser

from flask import Flask, render_template
from flask_cors import CORS

from config.settings import OUTPUT_DIR

from api.routes_genes import genes_bp
from api.routes_pipeline import pipeline_bp
from api.routes_files import files_bp
from api.routes_system import system_bp


# ============================================================
# INITIALISATION FLASK
# ============================================================

app = Flask(__name__)

CORS(app)


# ============================================================
# BLUEPRINTS
# ============================================================

app.register_blueprint(genes_bp)
app.register_blueprint(pipeline_bp)
app.register_blueprint(files_bp)
app.register_blueprint(system_bp)


# ============================================================
# PAGE PRINCIPALE
# ============================================================

@app.route("/")
def index():
    """Interface principale NanoVar."""
    return render_template("index.html")


# ============================================================
# DÉMARRAGE
# ============================================================

if __name__ == "__main__":

    print("=" * 60)
    print("🧬 NanoVar — Plateforme d'analyse Nanopore")
    print("=" * 60)
    print(f"📁 Output         : {OUTPUT_DIR}")
    print("🌐 Serveur Flask  : http://127.0.0.1:5000")
    print("=" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    def open_browser():
        time.sleep(1.5)
        webbrowser.open_new(
            "http://127.0.0.1:5000/"
        )

    threading.Thread(
        target=open_browser,
        daemon=True
    ).start()

    app.run(
        host="127.0.0.1",
        port=5000,
        debug=False,
        threaded=True,
        use_reloader=False
    )