# -*- coding: utf-8 -*-
import os, re, sys

root = r".\data"
# Heurísticas rápidas:
# - muy pequeño => probablemente HTML/404
SMALL = 100_000  # 100 KB (ajústalo: 50–200 KB)
# - patrón de float en las primeras líneas no vacías
FLOAT = re.compile(r'^[\+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][\+\-]?\d+)?\s*$')

bad, good = [], []
for name in os.listdir(root):
    if not name.startswith("zeros_") or not name.endswith(".dat"):
        continue
    p = os.path.join(root, name)
    try:
        sz = os.path.getsize(p)
        if sz < SMALL:
            bad.append(name); continue

        # lee ~200 líneas como máximo
        k, ok = 0, 0
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                s = line.strip()
                if not s or s[0] in "#;<":        # ignora vacías/comentarios/HTML
                    continue
                if FLOAT.match(s.split()[0]):     # primer token parece número
                    ok += 1
                k += 1
                if k >= 200: break

        if ok < 20:               # menos de 20 líneas numéricas entre las primeras ~200
            bad.append(name)
        else:
            good.append(name)
    except Exception:
        bad.append(name)

print("GOOD", len(good))
print("BAD", len(bad))
open("bad_scan.txt","w",encoding="utf-8").write("\n".join(bad))
