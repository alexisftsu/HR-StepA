import os, json, math, argparse
import numpy as np
import matplotlib.pyplot as plt

# --- Núcleos y FT auxiliares ---
def K_delta(x, Delta):
    x = np.asarray(x, float)
    out = np.empty_like(x)
    eps=1e-18
    z = np.where(np.abs(x)<eps, eps, x)
    s = np.sin(np.pi*Delta*z)
    out = (1.0/Delta)*(s/(np.pi*z))**2
    out[np.abs(x)<eps] = Delta
    return out

def tri_hat(ksi, Delta):
    return np.maximum(0.0, 1.0 - np.abs(ksi)/Delta)

def hat_chi(ksi, beta):
    return 2.0*beta*np.sinc(2.0*beta*ksi)

# --- Construcción base: (chi * KΔ) y bumps de borde ---
def build_initial(beta, Delta, xs):
    grid_step = 1.0/(60.0*Delta)
    U_extra   = 32.0/Delta
    U = beta + U_extra
    uu = np.arange(-U, U+grid_step, grid_step, float)
    Ku = K_delta(uu, Delta)
    # primitiva de KΔ (trapecios)
    G = np.concatenate(([0.0], np.cumsum(0.5*(Ku[1:]+Ku[:-1])*grid_step)))
    def antideriv(t):
        tt = np.clip(t, uu[0], uu[-1])
        idx = (tt - uu[0]) / grid_step
        i0  = np.floor(idx).astype(int)
        i1  = np.minimum(i0+1, len(uu)-1)
        frac= idx - i0
        return G[i0] + frac*(G[i1]-G[i0])
    def conv(x):
        a = x - beta; b = x + beta
        return antideriv(b) - antideriv(a)
    base = conv(xs)
    bump_edge = (K_delta(xs-beta, Delta) + K_delta(xs+beta, Delta))
    return base, bump_edge

# --- Refuerzo estable (mayorante/menorante) ---
def enforce_majorant(xs, S, chi, bump, itmax=800):
    S = S.copy()
    for _ in range(itmax):
        gap = S - chi
        mg  = float(np.min(gap))
        if mg >= -1e-11: break
        i   = int(np.argmin(gap))
        den = max(float(bump[i]), 1e-3)
        alpha = -1.05*mg/den
        alpha = min(max(alpha, 0.0), 0.2)  # tope por iteración
        S += alpha*bump
    return S

def enforce_minorant(xs, S, chi, bump, itmax=800):
    S = S.copy()
    for _ in range(itmax):
        gap = chi - S
        mg  = float(np.min(gap))
        if mg >= -1e-11: break
        i   = int(np.argmin(gap))
        den = max(float(bump[i]), 1e-3)
        alpha = -1.05*mg/den
        alpha = min(max(alpha, 0.0), 0.2)
        S -= alpha*bump
    return S

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta",  type=float, required=True)
    ap.add_argument("--Delta", type=float, required=True)
    ap.add_argument("--out",   type=str,   required=True)
    ap.add_argument("--gridN", type=int,   default=30001)
    ap.add_argument("--plot",  action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    beta = args.beta; Delta = args.Delta
    # ventana amplia en x para certificar
    xs   = np.linspace(-beta-24.0/Delta, beta+24.0/Delta, args.gridN)
    chi  = ((xs>=-beta)&(xs<=beta)).astype(float)

    # base y bump compuesto (bordes ±β + centro) normalizado por Δ
    base, bump_edge = build_initial(beta, Delta, xs)
    bump_center = K_delta(xs, Delta)
    B = (bump_edge + 0.5*bump_center) / Delta

    # arranque c=0.5 y refuerzo estable
    Splus  = base + 0.5*B
    Sminus = base - 0.5*B
    Splus  = enforce_majorant(xs, Splus,  chi, B)
    Sminus = enforce_minorant(xs, Sminus, chi, B)

    # DC nudge si quedara residuo negativo
    gapP = Splus - chi
    gapM = chi   - Sminus
    okp = bool(np.min(gapP) >= -1e-12)
    okm = bool(np.min(gapM) >= -1e-12)
    if not okp or not okm:
        need_up   = max(0.0, -float(np.min(gapP)))
        need_down = max(0.0, -float(np.min(gapM)))
        eps = 1.02*max(need_up, need_down)
        if eps > 0.0:
            Splus  = Splus  + eps
            Sminus = Sminus - eps
            gapP = Splus - chi
            gapM = chi   - Sminus
            okp = bool(np.min(gapP) >= -1e-12)
            okm = bool(np.min(gapM) >= -1e-12)

    # L1 en malla (objetivo ~1/Δ)
    Epg = float(np.trapezoid(np.maximum(gapP,0.0), xs))
    Emg = float(np.trapezoid(np.maximum(gapM,0.0), xs))

    # FT (solo diagnóstico)
    K=3001
    ksis = np.linspace(-Delta, Delta, K)
    lam  = tri_hat(ksis, Delta)
    hchi = hat_chi(ksis, beta)

    payload = {
        "model":"PW-real",
        "params":{"beta":beta,"Delta":Delta,"gridN":args.gridN},
        "valid":{"majorant_ok":okp,"minorant_ok":okm,
                 "min_gap_majorant":float(np.min(gapP)),
                 "min_gap_minorant":float(np.min(gapM))},
        "L1_errors":{"E_plus_grid":Epg,"E_minus_grid":Emg,"target_approx":1.0/Delta},
        "fourier_grid":{"ksi_min":-Delta,"ksi_max":Delta,"K":K},
        "fourier":{"Lambda":lam.tolist(),"hat_chi":hchi.tolist()}
    }
    outj = os.path.join(args.out, f"BS_PW_beta{beta:g}_D{Delta:g}.json")
    with open(outj,"w",encoding="utf-8") as f: json.dump(payload,f,indent=2)
    print(outj)

    if args.plot:
        import matplotlib.ticker as mt
        fig,ax = plt.subplots(figsize=(7.5,3.6))
        ax.plot(xs, chi,   lw=2, label="χ_{[-β,β]}")
        ax.plot(xs, Splus, lw=1.4, label="S^+ (PW, enforced)")
        ax.plot(xs, Sminus,lw=1.4, label="S^- (PW, enforced)")
        ax.set_title(f"PW real-line: β={beta}, Δ={Delta}")
        ax.grid(alpha=0.25); ax.legend()
        ax.xaxis.set_major_locator(mt.MaxNLocator(11))
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, f"BS_PW_beta{beta:g}_D{Delta:g}.png"), dpi=140)
