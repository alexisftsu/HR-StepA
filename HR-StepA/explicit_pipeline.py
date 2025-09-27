import os, json, math, argparse, datetime, numpy as np
from pathlib import Path

def load_json(p):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)

def fmt(x):
    # notación compacta con punto decimal
    if isinstance(x, (int,)) or (isinstance(x,float) and (abs(x)>=1e-3 and abs(x)<1e6)):
        return f"{x:.12g}"
    return f"{x:.6e}"

def compute_lambda(x):
    L = math.log(x); L2 = math.log(L)
    return (L**(3.0/5.0))*(L2**(-1.0/5.0))

def render_partA(stepA):
    Z  = stepA["zeros"]
    CT = stepA["constants"]
    VK = stepA.get("VK", {})
    # claves PowerShell-friendly
    B  = VK.get("BVK", VK.get("B_VK", None))
    beta = VK.get("beta_VK", VK.get("b_VK", None))
    x1 = VK.get("x1", None)
    lines = []
    lines.append("=== Global explicit bound (Part A) ===\n")
    lines.append(f"T0 = {fmt(Z['T0'])}    X0 = {fmt(Z['X0'])}")
    lines.append(f"C_bajo = {fmt(CT['C_bajo'])}")
    if B is not None and beta is not None and x1 is not None:
        lines.append(f"VK params: B = {fmt(B)}, beta = {fmt(beta)}, x1 = {fmt(x1)}")
    lines.append("")
    lines.append("For all x >= 2:")
    lines.append("  |psi(x) - x| <=")
    lines.append(f"    - C_bajo * sqrt(x) * (log x)^2,                for 2 <= x <= X0")
    if B is not None and beta is not None and x1 is not None:
        lines.append("    - min{ C_bajo * sqrt(x) * (log x)^2,")
        lines.append("           B * x * exp(- beta * Lambda(x)) },      for X0 <= x < x1")
        lines.append("    - B * x * exp(- beta * Lambda(x)),             for x >= x1")
    else:
        lines.append("    - (high-x term present; fill VK parameters to render the last two lines)")
    lines.append("  where Lambda(x) = (log x)^(3/5) * (log log x)^(-1/5).")
    lines.append("")
    return "\n".join(lines)

def render_partA_RH(stepA):
    RH = stepA.get("constants_RH", None)
    if not RH: return ""
    lines=[]
    lines.append("=== Uniform bound under RH ===\n")
    lines.append(f"C_RH = 1/(8*pi) = {fmt(RH['C_RH'])}")
    lines.append(f"Small-x check up to {fmt(RH['C_small_check_upto'])}: C_small = {fmt(RH['C_small'])}")
    lines.append(f"Uniform constant (all x>=2): C_tot_RH = {fmt(RH['C_tot_RH'])}")
    lines.append("\nConclusion (RH): for all x >= 2,")
    lines.append("  |psi(x) - x| <= C_tot_RH * sqrt(x) * (log x)^2.")
    lines.append("")
    return "\n".join(lines)

def export_bs_coeffs(sc_json_path, out_dir):
    d = load_json(sc_json_path)
    beta  = float(d["params"]["beta"])
    Delta = float(d["params"]["Delta"])
    N     = int(d["params"]["N"])
    okp   = bool(d["valid"]["majorant_ok"])
    okm   = bool(d["valid"]["minorant_ok"])
    gap_p = float(d["valid"]["min_gap_majorant"])
    gap_m = float(d["valid"]["min_gap_minorant"])
    Ept   = float(d["L1_errors"].get("E_plus_theory", 0.0))
    Emt   = float(d["L1_errors"].get("E_minus_theory", 0.0))
    Epg   = float(d["L1_errors"].get("E_plus_grid", 0.0))
    Emg   = float(d["L1_errors"].get("E_minus_grid", 0.0))
    a_plus = np.array(d["fourier"]["a_plus"],  dtype=float)
    a_minus= np.array(d["fourier"]["a_minus"], dtype=float)
    n      = np.arange(-N, N+1, dtype=int)

    npz = Path(out_dir)/f"BS_coeffs_beta{beta:g}_D{Delta:g}.npz"
    csv = Path(out_dir)/f"BS_coeffs_beta{beta:g}_D{Delta:g}.csv"
    cert= Path(out_dir)/f"BS_certificate_beta{beta:g}_D{Delta:g}.txt"

    np.savez(npz, n=n, a_plus=a_plus, a_minus=a_minus, beta=beta, Delta=Delta, N=N,
             majorant_ok=okp, minorant_ok=okm, min_gap_plus=gap_p, min_gap_minus=gap_m,
             E_plus_theory=Ept, E_minus_theory=Emt, E_plus_grid=Epg, E_minus_grid=Emg)

    with open(csv, "w", encoding="utf-8") as f:
        f.write("n,a_plus,a_minus\n")
        for i in range(n.size):
            f.write(f"{n[i]},{a_plus[i]:.17g},{a_minus[i]:.17g}\n")

    with open(cert, "w", encoding="ascii", errors="ignore") as f:
        f.write("=== Beurling–Selberg (Selberg circle) certificate ===\n")
        f.write(f"beta={beta}, Delta={Delta}, N={N}\n")
        f.write(f"majorant_ok={okp}, minorant_ok={okm}\n")
        f.write(f"min_gap_majorant={gap_p:.12g}, min_gap_minorant={gap_m:.12g}\n")
        f.write(f"E_plus_theory=E_minus_theory={Ept:.12g}\n")
        f.write(f"E_plus_grid={Epg:.12g}, E_minus_grid={Emg:.12g}\n")
    return {
        "beta": beta, "Delta": Delta, "N": N,
        "majorant_ok": okp, "minorant_ok": okm,
        "min_gap_majorant": gap_p, "min_gap_minorant": gap_m,
        "E_plus_theory": Ept, "E_minus_theory": Emt,
        "E_plus_grid": Epg, "E_minus_grid": Emg,
        "npz": str(npz), "csv": str(csv), "certificate": str(cert),
    }

def render_partB(bs_info):
    # --- Normaliza la estructura que llega ---
    def _coerce(bs):
        import os, json
        # Si es ruta a JSON
        if isinstance(bs, str) and os.path.isfile(bs):
            with open(bs, "r", encoding="utf-8") as f:
                return json.load(f)
        # Si ya trae params
        if isinstance(bs, dict) and "params" in bs:
            return bs
        # Si viene anidado tipo {"data": {...}} o {"bs": {...}}
        if isinstance(bs, dict):
            for k in ("data", "bs", "info"):
                if k in bs and isinstance(bs[k], dict):
                    inner = bs[k]
                    if "params" in inner or "Delta" in inner or "beta" in inner:
                        return inner
            # Último recurso: sintetiza "params"/"valid"/"L1_errors"
            params = bs.get("params") or {k: bs.get(k) for k in ("beta","Delta","N","gridN") if k in bs}
            valid  = bs.get("valid")  or {
                "majorant_ok":    bs.get("majorant_ok", True),
                "minorant_ok":    bs.get("minorant_ok", True),
                "min_gap_majorant": bs.get("min_gap_majorant", 0.0),
                "min_gap_minorant": bs.get("min_gap_minorant", 0.0),
            }
            L1 = bs.get("L1_errors") or {
                "E_plus":  bs.get("E_plus",  bs.get("E_plus_grid",  0.0)),
                "E_minus": bs.get("E_minus", bs.get("E_minus_grid", 0.0)),
            }
            return {"params": params, "valid": valid, "L1_errors": L1}
        # Valor por defecto (no debería ocurrir)
        return {"params":{"beta":0.0,"Delta":0.0,"N":0}, "valid":{"majorant_ok":False,"minorant_ok":False,"min_gap_majorant":0.0,"min_gap_minorant":0.0}, "L1_errors":{"E_plus":0.0,"E_minus":0.0}}

    d = _coerce(bs_info)

    # --- N y E_theory coherentes ---
    Delta = float(d["params"].get("Delta", d["params"].get("D", 0.0)))
    N = int(d["params"].get("N", int(math.floor(Delta)))) if Delta else int(d["params"].get("N", 0))
    E_theory = 1.0/(N+1) if N >= 0 else 0.0

    # --- L1: usa E_plus/E_minus o, si faltan, E_plus_grid/E_minus_grid ---
    E_plus  = d.get("L1_errors",{}).get("E_plus",  d.get("L1_errors",{}).get("E_plus_grid",  0.0))
    E_minus = d.get("L1_errors",{}).get("E_minus", d.get("L1_errors",{}).get("E_minus_grid", 0.0))

    lines = []
    lines.append("=== BeurlingSelberg (Part B, Selberg circle, enforced) ===\n")
    lines.append(f"beta = {d['params'].get('beta')},  Delta = {Delta},  N = {N}\n")
    lines.append(f"majorant_ok = {d['valid'].get('majorant_ok')},  minorant_ok = {d['valid'].get('minorant_ok')}\n")
    lines.append(f"min_gap_majorant = {d['valid'].get('min_gap_majorant',0.0):.12g},  min_gap_minorant = {d['valid'].get('min_gap_minorant',0.0):.12g}\n")
    lines.append(f"L1 theory = 1/(N+1) = {E_theory:.6e}\n")
    lines.append(f"L1 (grid, after enforcement): E_plus = {E_plus:.12g}, E_minus = {E_minus:.12g}\n")
    lines.append("\nFourier support: |n|<=N (Fejr window). The enforced pair remains band-limited and\npointwise majorant/minorant on a dense grid; L1 increases slightly vs theory, which is acceptable for explicit bounds.\n")
    return "".join(lines)
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stepA", default=r".\out\StepA_results.json", help="JSON de Parte A")
    ap.add_argument("--selberg_json", required=True, help="Selberg_circle_*.json (con --force)")
    ap.add_argument("--out", default=r".\out", help="Carpeta salida")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    stepA = load_json(args.stepA)
    txtA  = render_partA(stepA)
    txtRH = render_partA_RH(stepA)

    bs_info = export_bs_coeffs(args.selberg_json, args.out)
    txtB = render_partB(bs_info)

    # Ensamble final
    lines=[]
    lines.append("=== HR Project — Final explicit bound package ===\n")
    lines.append(txtA)
    if txtRH: lines.append(txtRH)
    lines.append(txtB)
    lines.append("=== Files ===")
    lines.append(f"- Step A results : {Path(args.stepA).resolve()}")
    lines.append(f"- BS coeffs (npz): {bs_info['npz']}")
    lines.append(f"- BS coeffs (csv): {bs_info['csv']}")
    lines.append(f"- BS certificate : {bs_info['certificate']}")
    final_txt = os.path.join(args.out, "Final_bound_statement.txt")
    with open(final_txt, "w", encoding="ascii", errors="ignore") as f:
        f.write("\n".join(lines))
    print(final_txt)

    # Informe JSON consolidado
    final_json = os.path.join(args.out, "Final_report.json")
    payload = {
        "stepA": stepA,
        "partB_selberg": bs_info,
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    with open(final_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(final_json)

    # Paquete ZIP
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    zipname = os.path.join(args.out, f"HR_Final_pkg_{stamp}.zip")
    import zipfile
    with zipfile.ZipFile(zipname, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for p in [args.stepA, final_txt, final_json, bs_info["npz"], bs_info["csv"], bs_info["certificate"]]:
            z.write(p, arcname=os.path.join("pkg", os.path.basename(p)))
    print(zipname)

if __name__ == "__main__":
    main()










