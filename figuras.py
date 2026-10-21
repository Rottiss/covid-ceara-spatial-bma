"""Generate seven spatial figures from the included result tables."""
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
language = os.environ.get("FIGURE_LANGUAGE", "en").lower()
if language not in {"en", "pt"}:
    raise SystemExit("FIGURE_LANGUAGE must be 'en' or 'pt'.")
env = os.environ.copy()
env["FIGURE_LANGUAGE"] = language
env["MPLBACKEND"] = "Agg"
env["FIGURAS_SAIDA"] = str(ROOT / "execucoes_locais" / "figuras" / language)
env["DADOS_FIGURAS"] = str(ROOT / "execucoes_locais" / "dados_figuras")
if not (Path(env["DADOS_FIGURAS"]) / "pips_queen.csv").is_file():
    raise SystemExit("Run first: Rscript reproduzir.R dados-figuras")
for script in ("51_fig3_pip_especificacoes.py", "52_figuras_v3.py", "53_fig7_priori_sem.py"):
    subprocess.run([sys.executable, str(ROOT / "espacial" / script)], env=env, check=True)
print("Figures saved to:", env["FIGURAS_SAIDA"])
