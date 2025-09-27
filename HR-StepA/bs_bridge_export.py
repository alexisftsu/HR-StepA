import os, json, argparse, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--sc-json", required=True, help="Selberg_circle_*.json")
ap.add_argument("--out", required=True, help="Carpeta de salida")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

d = json.load(open(args.sc_json, "r", encoding="utf-8"))
beta  = float(d["params"]["beta"])
Delta = float(d["params"]["Delta"])
N     = int(d["params"]["N"])
okp   = bool(d["valid"]["majorant_ok"])
okm   = bool(d["valid"]["minorant_ok"])
gap_p = float(d["valid"]["min_gap_majorant"])
gap_m = float(d["valid"]["min_gap_minorant"])
Ept   = float(d["L1_errors"]["E_plus_theory"])
Emt   = float(d["L1_errors"]["E_minus_theory"])
Epg   = float(d["L1_errors"]["E_plus_grid"])
Emg   = float(d["L1_errors"]["E_minus_grid"])

n      = np.arange(-N, N+1, dtype=int)
a_plus = np.array(d["fourier"]["a_plus"],  dtype=float)
a_minus= np.array(d["fourier"]["a_minus"], dtype=float)

# Guardar en formatos prácticos
np.savez(os.path.join(args.out, f"BS_coeffs_beta{beta:g}_D{Delta:g}.npz"),
         n=n, a_plus=a_plus, a_minus=a_minus, beta=beta, Delta=Delta, N=N,
         majorant_ok=okp, minorant_ok=okm, min_gap_plus=gap_p, min_gap_minus=gap_m,
         E_plus_theory=Ept, E_minus_theory=Emt, E_plus_grid=Epg, E_minus_grid=Emg)

with open(os.path.join(args.out, f"BS_coeffs_beta{beta:g}_D{Delta:g}.csv"), "w", encoding="utf-8") as f:
    f.write("n,a_plus,a_minus\n")
    for i in range(n.size):
        f.write(f"{n[i]},{a_plus[i]:.17g},{a_minus[i]:.17g}\n")

with open(os.path.join(args.out, f"BS_certificate_beta{beta:g}_D{Delta:g}.txt"), "w", encoding="ascii", errors="ignore") as f:
    f.write("=== Beurling–Selberg (Selberg circle version) certificate ===\n")
    f.write(f"beta={beta}, Delta={Delta}, N=floor(Delta)={N}\n")
    f.write(f"majorant_ok={okp}, minorant_ok={okm}\n")
    f.write(f"min_gap_majorant={gap_p:.12g}, min_gap_minorant={gap_m:.12g}\n")
    f.write(f"E_plus_theory=E_minus_theory=1/(N+1)={Ept:.12g}\n")
    f.write(f"E_plus_grid={Epg:.12g}, E_minus_grid={Emg:.12g}\n")
print("OK")
