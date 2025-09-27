#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, hashlib, os, re, sys, gzip
from typing import Dict, List, Optional, Tuple
from decimal import Decimal, InvalidOperation

NUM_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")
HTML_HINTS = ("<!doctype html","<!DOCTYPE html","<html","cloudflare","captcha",
              "please enable javascript","access denied","verifying you are human",
              "<title>lmfdb","beta.lmfdb.org")

def detect_html_gate(s:str)->bool:
    l=s.lower()
    return any(h in l for h in HTML_HINTS)

def md5_of(path:str, buf:int=1<<20)->str:
    h=hashlib.md5()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(buf), b""): h.update(b)
    return h.hexdigest()

def load_md5(md5path:str)->Dict[str,str]:
    m={}
    with open(md5path,"rt",encoding="utf-8",errors="ignore") as f:
        for r in f:
            line=r.strip()
            if not line or line[0] in "#;": continue
            parts=line.split()
            if len(parts)<2: continue
            md5=parts[0]
            star=line.find('*')
            name=line[star+1:].strip() if star>=0 else parts[-1]
            m[name]=md5
    return m

def basename_height(fn:str)->Optional[int]:
    m=re.search(r"zeros_(\d+)\.dat$", fn)
    return int(m.group(1)) if m else None

def is_int_decimal(x:Decimal)->bool:
    try: return x==x.to_integral_value()
    except: return False

def parse_nums(line:str)->List[Decimal]:
    out=[]
    for tok in NUM_RE.findall(line):
        try: out.append(Decimal(tok))
        except (InvalidOperation, ValueError): continue
    return out

def candidate_list(nums:List[Decimal], start_thr:Decimal) -> List[Decimal]:
    thr = max(Decimal(10), start_thr)
    return [v for v in nums if v >= thr]

def open_text_auto(path:str):
    # Si el fichero es GZIP (1F 8B), abrimos en modo texto via gzip; si no, como texto normal
    with open(path,"rb") as fb:
        head = fb.read(2)
    if head == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    else:
        return open(path, "rt", encoding="utf-8", errors="ignore")

def fast_seek_first_gamma(path:str, start_thr:Decimal, max_scan:int=200000) -> Tuple[Optional[Decimal], int]:
    f = open_text_auto(path)
    try:
        for ln, raw in enumerate(f,1):
            line = raw.strip()
            if not line: continue
            if detect_html_gate(line): continue
            nums = parse_nums(line)
            if not nums: continue
            cands = candidate_list(nums, start_thr)
            if cands:
                return (min(cands), ln)
            if ln >= max_scan:
                break
    finally:
        f.close()
    return (None, -1)

def read_dat_file(path:str, max_skips:int, start_frac:float, max_step:float, badlog:List[str])->List[str]:
    fn=os.path.basename(path)
    base=basename_height(fn)
    if base and base>=1000:
        start_thr = Decimal(base) * Decimal(str(start_frac))
    else:
        start_thr = Decimal(10)

    gammas:List[str]=[]
    last:Optional[Decimal]=None
    max_step_dec = Decimal(str(max_step))

    g0, ln0 = fast_seek_first_gamma(path, start_thr, max_scan=max_skips)
    if g0 is None:
        raise ValueError(f"No se encontró ningún valor ≥ start_thr en {path} tras {max_skips} líneas (start_thr={start_thr})")

    f = open_text_auto(path)
    try:
        for _ in range(ln0-1):
            f.readline()

        gammas.append(str(g0))
        last = g0
        skips = 0

        for ln_off, raw in enumerate(f, ln0):
            line=raw.strip()
            if not line or line.startswith(('#',';')): continue
            if detect_html_gate(line): continue

            nums = parse_nums(line)
            if not nums:
                skips += 1
                if skips>max_skips:
                    snippet = (line[:120] if line else "")
                    raise ValueError(f"Demasiadas líneas sin candidato válido en {path} (>{max_skips}); última línea: '{snippet}' con last={last}")
                continue

            cands:List[Decimal] = []
            if len(nums)>=2 and is_int_decimal(nums[0]):
                g_direct = nums[1]
                if g_direct >= max(Decimal(10), start_thr):
                    cands.append(g_direct)
            cands.extend([v for v in nums if v >= max(Decimal(10), start_thr)])

            if not cands:
                skips += 1
                if skips>max_skips:
                    snippet = (line[:120] if line else "")
                    raise ValueError(f"Demasiadas líneas sin candidato válido en {path} (>{max_skips}); última línea: '{snippet}' con last={last}")
                continue

            deltas = [(g, g - last) for g in cands if g > last]
            deltas = [(g, d) for (g, d) in deltas if d <= max_step_dec]
            if not deltas:
                skips += 1
                if skips>max_skips:
                    snippet = (line[:120] if line else "")
                    raise ValueError(f"Demasiadas no-monotonías/ausencias en {path} (>{max_skips}); última: '{snippet}' con last={last}")
                continue

            g, d = min(deltas, key=lambda t: t[1])
            if g <= last:
                skips += 1
                badlog.append(f"{fn}:{ln_off}: non-mono {g} after {last} | '{line}'")
                if skips>max_skips:
                    raise ValueError(f"Demasiadas no-monotonías en {path} (>{max_skips}); última: {g} tras {last}")
                continue

            gammas.append(str(g))
            last=g
    finally:
        f.close()

    if not gammas:
        raise ValueError(f"Sin datos numéricos en {path}")
    return gammas

def main():
    ap=argparse.ArgumentParser(description="Une zeros_*.dat (LMFDB) en lista de γ creciente (soporta GZIP).")
    ap.add_argument("--in", dest="indir", required=True, help="directorio con zeros_*.dat")
    ap.add_argument("--out", dest="out", required=True, help="fichero de salida (lista de gammas)")
    ap.add_argument("--md5", dest="md5file", help="md5.txt para verificar")
    ap.add_argument("--pattern", default="zeros_", help="prefijo de archivos (por defecto zeros_)")
    ap.add_argument("--max-nonmono", type=int, default=200000, help="límite de líneas problemáticas por fichero")
    ap.add_argument("--start-frac", type=float, default=0.35, help="fracción de la altura base para arrancar (p.ej. 0.35)")
    ap.add_argument("--max-step", type=float, default=5.0, help="salto máximo permitido entre γ consecutivos")
    ap.add_argument("--only", default="", help="si se da, procesa solo ese nombre de fichero")
    args=ap.parse_args()

    files=[fn for fn in os.listdir(args.indir) if fn.startswith(args.pattern) and fn.endswith(".dat")]
    files.sort(key=lambda s:(len(s),s))
    if args.only:
        files=[fn for fn in files if fn==args.only]
        if not files:
            print(f"[ERROR] --only {args.only} no encontrado en {args.indir}", file=sys.stderr); sys.exit(1)
    if not files:
        print(f"[WARN] No hay '{args.pattern}*.dat' en {args.indir}", file=sys.stderr); sys.exit(1)

    if args.md5file:
        md5map=load_md5(args.md5file)
        ok=True; chk=0
        for fn in files:
            if fn not in md5map: continue
            p=os.path.join(args.indir, fn)
            got=md5_of(p).lower(); exp=md5map[fn].lower()
            print(f"[MD5] {fn:>20}  {'OK' if got==exp else 'MISMATCH'}" + ("" if got==exp else f" (got {got}, expected {exp})"))
            ok = ok and (got==exp); chk+=1
        print(f"[MD5] Checked {chk} files; overall: {'OK' if ok else 'MISMATCH'}")

    out_all:List[str]=[]
    badlog:List[str]=[]

    for fn in files:
        p=os.path.join(args.indir, fn)
        try:
            g=read_dat_file(p,
                            max_skips=args.max_nonmono,
                            start_frac=args.start_frac,
                            max_step=args.max_step,
                            badlog=badlog)
        except Exception as e:
            print(f"[ERROR] {fn}: {e}", file=sys.stderr); sys.exit(2)

        for s in g:
            v=Decimal(s)
            if out_all and v <= Decimal(out_all[-1]):
                badlog.append(f"{fn}: global non-mono {v} after {out_all[-1]} (descartado)")
                continue
            out_all.append(str(v))

    with open(args.out,"wt",encoding="utf-8") as w:
        for s in out_all: w.write(s+"\n")
    print(f"[DONE] wrote {len(out_all)} γ to {args.out}")
    if badlog:
        bl=args.out+".badlines.log"
        with open(bl,"wt",encoding="utf-8") as w: w.write("\n".join(badlog))
        print(f"[INFO] líneas ignoradas: {len(badlog)} (ver {bl}); start_frac={args.start_frac}, max_step={args.max_step}")

if __name__=="__main__":
    main()
