"""Check Python syntax, encoding, and personal paths in text files.

The --normalizar-utf8 option converts Windows-1252 files to UTF-8.
Dependencies computed at runtime are not resolved.
"""
import argparse
import ast
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--normalizar-utf8", action="store_true")
args = parser.parse_args()
issues = []
converted = []
inventory = []
for folder in (ROOT, ROOT / "espacial", ROOT / "estimacoes", ROOT / "scripts"):
    for path in sorted(folder.iterdir()):
        if path.suffix not in (".R", ".py", ".md") or not path.is_file():
            continue
        rel = path.relative_to(ROOT).as_posix()
        raw = path.read_bytes()
        try:
            text = raw.decode("utf-8-sig")
        except UnicodeDecodeError:
            if not args.normalizar_utf8:
                issues.append(f"{rel}: encoding is not UTF-8")
                continue
            text = raw.decode("cp1252")
            path.write_bytes(text.encode("utf-8"))
            converted.append(rel)
        if path.resolve() != Path(__file__).resolve() and re.search(r"[A-Za-z]:[/\\]Users[/\\]|/home/[^/]+/|/Users/[^/]+/", text):
            issues.append(f"{rel}: absolute personal path")
        if path.suffix == ".py":
            try:
                ast.parse(text, filename=rel)
            except SyntaxError as exc:
                issues.append(f"{rel}: {exc}")
        if path.suffix == ".R":
            reads = [line.strip() for line in text.splitlines()
                     if re.search(r"readRDS\(|read[.]csv\(|read_excel\(|source\(|load\(", line)]
            inventory.append({"script": rel, "leituras_estaticas": reads})
report = {"problemas": issues, "convertidos_utf8": converted,
          "leituras_r": inventory,
          "limite": "Static inspection; does not guarantee that all dependencies exist."}
out = ROOT / "execucoes_locais" / "auditoria_codigo.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"R scripts inventoried: {len(inventory)}; files converted: {len(converted)}")
for item in converted + issues:
    print(item)
print("Report:", out)
raise SystemExit(bool(issues))
