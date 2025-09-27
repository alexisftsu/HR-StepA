import os, json, math, argparse
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt

# ---------- Utils ----------
def K_delta(x, Delta):
    # KΔ(x) = (1/Δ) * (sin(πΔ x)/(π x))^2, KΔ >= 0, ∫KΔ = 1,  FT(KΔ)=Λ(|ξ|/Δ)
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    # manejar x=0
    eps = 1e-30
    z = np.where(np.abs(x)<eps, eps, x)
    s = np.sin(np.pi*Delta*z)
    out = (1.0/Delta)*(s/(np.pi*z))**2
    out[np.abs(x)<eps] = (1.0/Delta)*(Delta)**2  # límite: (1/Δ)*(Δ/π)^2 * π^2? simplifica a Δ (pero ya lo controla abajo)
    # Corrige el valor en 0 analíticamente:  KΔ(0) = Δ/π^2 * π^2?  En realidad: lim x->0 (sin(πΔ x)/(π x))^2 / Δ = Δ
    out[np.abs(x)<eps] = Delta
    return out

def tri_hat(ksi, Delta):
    # Λ(|ξ|/Δ) = max(0, 1 - |ξ|/Δ)
    r = np.maximum(0.0, 1.0 - np.abs(ksi)/Delta)
    return r

def hat_chi(ksi, beta):
    # \hat{chi}_{[-β,β]}(ξ) = ∫_{-β}^β e^{-2π i ξ y} dy = 2β * sinc(2β ξ)  (sinc normalizado sin(πx)/(πx))
    return 2.0*beta*np.sinc(2.0*beta*ksi)

def conv_chi_K(x, beta, Delta, grid_step=None, U_extra=None):
    # Convolución precisa: opción "tight" vía grid_step y U_extra desde variables globales
    global __GRID_STEP_FACTOR, __U_EXTRA_FACTOR
    if grid_step is None:
        grid_step = 1.0/( (__GRID_STEP_FACTOR or 20.0) * Delta )
    if U_extra is None:
        U_extra = (__U_EXTRA_FACTOR or 12.0)/Delta
    U = 2.0*beta + U_extra
    uu = np.arange(-U, U+grid_step, grid_step, dtype=float)
    Ku = K_delta(uu, Delta)
    G = np.concatenate(([0.0], np.cumsum(0.5*(Ku[1:]+Ku[:-1])*grid_step)))
    def antideriv(t):
        tt = np.clip(t, uu[0], uu[-1])
        idx = (tt - uu[0]) / grid_step
        i0 = np.floor(idx).astype(int)
        i1 = np.minimum(i0+1, len(uu)-1)
        frac = idx - i0
        return G[i0] + frac * (G[i1]-G[i0])
    a = x - beta; b = x + beta
    return antideriv(b) - antideriv(a)

def build_selberg(beta, Delta, c=0.5, x_min=None, x_max=None, Nx=4001):
    # rango para chequear: un poco más allá del intervalo
    if x_min is None: x_min = -beta - 6.0/Delta
    if x_max is None: x_max =  beta + 6.0/Delta
    xs = np.linspace(x_min, x_max, Nx)
    conv = conv_chi_K(xs, beta, Delta)
    bump = (c/Delta)*(K_delta(xs-beta, Delta) + K_delta(xs+beta, Delta))
    S_plus  = conv + bump
    S_minus = conv - bump
    chi = ((xs >= -beta) & (xs <= beta)).astype(float)
    # márgenes
    gap_plus  = S_plus - chi
    gap_minus = chi - S_minus
    min_gap_plus  = float(np.min(gap_plus))
    min_gap_minus = float(np.min(gap_minus))
    ok_plus  = (min_gap_plus  >= -1e-12)
    ok_minus = (min_gap_minus >= -1e-12)
    return xs, chi, S_plus, S_minus, ok_plus, ok_minus, min_gap_plus, min_gap_minus

def hat_S(beta, Delta, c, ksis):
    # \hat S^±(ξ) = Λ(|ξ|/Δ) * [ \hat χ(ξ) ± (2c/Δ) cos(2πβξ) ]
    lam = tri_hat(ksis, Delta)
    hchi = hat_chi(ksis, beta)
    bump = (2.0*c/Delta)*np.cos(2.0*np.pi*beta*ksis)
    Sphat = lam*(hchi + bump)   # mayorante
    Smhat = lam*(hchi - bump)   # minorante
    return Sphat, Smhat, lam, hchi, bump

# ---------- CLI ----------
ap = argparse.ArgumentParser()
ap.add_argument("--beta", type=float, required=True, help="Semiancho del intervalo [-beta, beta]")
ap.add_argument("--Delta", type=float, required=True, help="Bandlimit (soporte de Fourier en |ξ|<=Delta)")
ap.add_argument("--out", type=str, required=True, help="Carpeta de salida (se crea si no existe)")
ap.add_argument("--gridN", type=int, default=4001, help="Puntos de malla en x para validar mayoración/menoración")
ap.add_argument("--auto-c", action="store_true", help="Autoajusta c>=0.5 hasta que S^-<=chi<=S^+ pase en malla")
ap.add_argument("--plot", action="store_true", help="Guarda PNG de S^± y χ")
ap.add_argument("--c0", type=float, default=0.5, help="Valor inicial de c (por defecto 0.5)")
ap.add_argument("--tight", action="store_true", help="Convolución más precisa (malla fina y cola más larga)")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

beta  = args.beta
Delta = args.Delta
# control de precisión global para conv
__GRID_STEP_FACTOR = 40.0 if args.tight else 20.0
__U_EXTRA_FACTOR   = 20.0 if args.tight else 12.0
c = args.c0
xs, chi, Splus, Sminus, okp, okm, gpp, gpm = build_selberg(beta, Delta, c=c, Nx=args.gridN)

if args.auto_c:
    k=0
    while (not okp or not okm) and c < 5.0:
        c *= 1.05
        xs, chi, Splus, Sminus, okp, okm, gpp, gpm = build_selberg(beta, Delta, c=c, Nx=args.gridN)
        k+=1

# L1-errors aproximados (trapecios)
def trapz(y,x): return float(np.trapezoid(y, x))
E_plus  = trapz(Splus-chi, xs)  # ≈ 1/Δ si c≈0.5
E_minus = trapz(chi-Sminus, xs)

# Fourier: malla en [-Delta, Delta]
K = 2001
ksis = np.linspace(-Delta, Delta, K)
Sphat, Smhat, lam, hchi, bump = hat_S(beta, Delta, c, ksis)

# Salida JSON
payload = {
  "params": {"beta": beta, "Delta": Delta, "c_bump": c, "gridN": args.gridN, "auto_c_used": args.auto_c},
  "valid": {"majorant_ok": bool(okp), "minorant_ok": bool(okm),
            "min_gap_majorant": gpp, "min_gap_minorant": gpm},
  "L1_errors": {"E_plus": E_plus, "E_minus": E_minus, "target_1_over_Delta": 1.0/Delta},
  "fourier_grid": {"ksi_min": -Delta, "ksi_max": Delta, "K": K},
  "fourier": {
      "hat_S_plus":  Sphat.tolist(),
      "hat_S_minus": Smhat.tolist(),
      "Lambda":      lam.tolist(),
      "hat_chi":     hchi.tolist(),
      "bump_term":   bump.tolist()
  }
}
out_json = os.path.join(args.out, f"BS_beta{beta:g}_D{Delta:g}.json")
with open(out_json,"w",encoding="utf-8") as f: json.dump(payload,f,indent=2)
print(out_json)

# Plot opcional
if args.plot:
    fig1,ax=plt.subplots(figsize=(7,4))
    ax.plot(xs, chi,   lw=2, label="chi_{[-β,β]}")
    ax.plot(xs, Splus, lw=1.5, label="S_Δ^+(x)")
    ax.plot(xs, Sminus,lw=1.5, label="S_Δ^-(x)")
    ax.set_title(f"BS majorant/minorant  β={beta}, Δ={Delta}, c={c:.3f}")
    ax.legend(); ax.grid(True, alpha=0.25)
    fig1.tight_layout()
    png = os.path.join(args.out, f"BS_beta{beta:g}_D{Delta:g}.png")
    fig1.savefig(png, dpi=140)






