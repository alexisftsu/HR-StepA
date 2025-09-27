
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
li_verify_full.py

Loads gamma ordinates from a file (one γ per line),
builds zeros ρ = 1/2 + iγ, plugs them into Criterio.LiCriterion,
and verifies λ_n > 0 for n ≤ Nmax with interval arithmetic used in Criterio.

Assumptions:
- The provided γ are certified to lie on the critical line (e.g., LMFDB verified range).
- Tail bounds beyond the last γ use the analytic bound in Criterio.

Usage:
  python li_verify_full.py --zeros zeros.txt --nmax 500 --t0_from_data

Options:
  --zeros PATH      Path to file with one gamma per line.
  --nmax N         Verify λ_n for 1..N (default 100)
  --prec P         Decimal precision (mpmath) (default 200)
  --t0_from_data   Use last |γ| as T0; otherwise T0=100.0

Outputs JSON with summary and writes a human-readable log.
"""

import argparse, json, time, sys
from pathlib import Path

# Import Criterio from same directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
import Criterio as C

def load_gammas(path: Path):
    gammas = []
    with path.open('r', encoding='utf-8', errors='ignore') as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"): 
                continue
            try:
                gammas.append(float(ln.split()[0]))
            except:
                pass
    if not gammas:
        raise RuntimeError("No gamma values parsed.")
    return gammas

class FileBackedZeros(C.ZetaZeros):
    """Provide LMFDB-backed zeros to Criterio by overriding constructor."""
    def __init__(self, gammas, T0=None):
        self.zeros = [complex(0.5, g) for g in gammas]
        # Ensure sorted by |Im|
        self.zeros.sort(key=lambda z: abs(z.imag))
        self.T0 = float(T0 if T0 is not None else abs(self.zeros[-1].imag))
    def get_zeros_up_to(self, T):
        # Return zeros with |Im| ≤ T (both positive and negative imaginary parts mirrored)
        # Criterio only uses positive imag part list; but provide both signs for completeness
        res = []
        for z in self.zeros:
            if z.imag <= T:
                res.append(z)
        return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zeros", required=True)
    ap.add_argument("--nmax", type=int, default=100)
    ap.add_argument("--prec", type=int, default=200)
    ap.add_argument("--t0_from_data", action="store_true")
    args = ap.parse_args()

    path = Path(args.zeros)
    gammas = load_gammas(path)
    T0 = abs(gammas[-1]) if args.t0_from_data else 100.0

    # Monkey-patch Criterio to use FileBackedZeros
    fbz = FileBackedZeros(gammas, T0=T0)
    li = C.LiCriterion(c=0.1, T0=T0)
    li.zeta_zeros = fbz  # replace provider

    C.mp.mp.dps = args.prec

    start = time.time()
    ok_all = True
    results = []
    for n in range(1, args.nmax+1):
        ok, inter = li.verify_lambda_positive(n)
        ok_all &= ok
        results.append({"n": n, "ok": bool(ok), "interval": [float(inter.a), float(inter.b)] if hasattr(inter, 'a') else None})
    elapsed = time.time() - start

    summary = {
        "zeros_file": str(path),
        "count_zeros": len(gammas),
        "T0": T0,
        "nmax": args.nmax,
        "prec": args.prec,
        "ok_all": ok_all,
        "elapsed_s": elapsed,
    }
    out_json = Path("li_verify_full_summary.json")
    out_json.write_text(json.dumps({"summary": summary, "results": results}, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
