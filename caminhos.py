"""Project directories, relative to this file's location."""
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
DADOS = RAIZ / "dados"
ESTIMACOES = RAIZ / "estimacoes"
ESPACIAL = RAIZ / "espacial"
MANUSCRITO = RAIZ / "manuscrito"
FIGURAS = RAIZ / "figuras"
DOCS = RAIZ / "docs"

for pasta in (ESTIMACOES, ESPACIAL):
    if not pasta.is_dir():
        raise FileNotFoundError(pasta)
