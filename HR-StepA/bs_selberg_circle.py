import os, json, math, argparse
import numpy as np
import matplotlib.pyplot as plt

# =========================
#  Selberg (circle version)
#  - Fourier support: |n| <= N = floor(Delta)
#  - Baseline coefficients (Fejer-weighted)
#  - Optional "force" step using Fejer bumps at ±beta to enforce pointwise
#    majorant/minorant on a dense grid (adds a small numerical slack)
# =========================

def selberg_coeffs(beta: float, Delta: float):
    N = int(math.floor(Delta))
    if N < 1:
        raise ValueError("Delta must be >= 1")
    a0_plus  = 2.0*beta + 1.0/(N+1)
    a0_minus = 2.0*beta - 1.0/(N+1)
    n = np.arange(-N, N+1, dtype=int)
    coeff = np.zeros((2, 2*N+1), dtype=float)  # [0]=plus, [1]=minus

    # n = 0
    coeff[0, N] = a0_plus
    coeff[1, N] = a0_minus

    # n != 0
    for k in range(1, N+1):
        w = 1.0 - k/(N+1.0)         # Fejer weight
        s = math.sin(2.0*math.pi*beta*k)/(math.pi*k)
        c = w*s
        coeff[0, N+k] =  c
        coeff[0, N-k] =  c
        coeff[1, N+k] =  c
        coeff[1, N-k] =  c
    return n, coeff  # integer freqs, coeffs for S^+, S^-

def eval_trig(n, a, xs):
    # S(x) = sum_{|n|<=N} a_n e^{2π i n x}, with real-even coeffs -> cosine sum
    N = (len(n)-1)//2
    a0 = a[N]
    val = np.full_like(xs, a0, dtype=float)
    for k in range(1, N+1):
        ak = a[N+k]  # = a[N-k]
        val += 2.0*ak * np.cos(2.0*np.pi*k*xs)
    return val

# ---------- Fejer bump helpers (for "force") ----------
def _Fejer(N, x):
    s = math.sin(math.pi*(N+1)*x)
    d = math.sin(math.pi*x)
    if abs(d) < 1e-14:
        return N+1.0
    return (1.0/(N+1.0)) * (s/d)**2

def _bump_fejer(N, beta, xs):
    # symmetric bumps centered at ±beta
    xs = np.asarray(xs)
    return np.array([_Fejer(N, x-beta) + _Fejer(N, x+beta) for x in xs], dtype=float)

def _enforce_majorant(xs, S, chi, N, beta, itmax=30):
    B = _bump_fejer(N, beta, xs)
    S = S.copy()
    for _ in range(itmax):
        gap = S - chi
        mg = float(np.min(gap))
        if mg >= -1e-10:
            break
        i = int(np.argmin(gap))
        alpha = -1.01*mg/max(B[i], 1e-16)  # lift at worst point +1%
        S += alpha*B
    return S

def _enforce_minorant(xs, S, chi, N, beta, itmax=30):
    B = _bump_fejer(N, beta, xs)
    S = S.copy()
    for _ in range(itmax):
        gap = chi - S
        mg = float(np.min(gap))
        if mg >= -1e-10:
            break
        i = int(np.argmin(gap))
        alpha = -1.01*mg/max(B[i], 1e-16)  # push down at worst point +1%
        S -= alpha*B
    return S

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--Delta", type=float, required=True, help="degree N = floor(Delta)")
    ap.add_argument("--out", type=str, required=True)
    ap.add_argument("--gridN", type=int, default=20001)
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--force", action="store_true", help="Enforce S^± to be pointwise majorant/minorant on the grid")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    beta = args.beta
    Delta = args.Delta
    n, coeffs = selberg_coeffs(beta, Delta)
    a_plus, a_minus = coeffs[0], coeffs[1]
    N = int(math.floor(Delta))

    # validate on one period [-1/2, 1/2]
    xs = np.linspace(-0.5, 0.5, args.gridN)
    Splus  = eval_trig(n, a_plus, xs)
    Sminus = eval_trig(n, a_minus, xs)
    chi = ((xs >= -beta) & (xs <= beta)).astype(float)

    if args.force:
        Splus  = _enforce_majorant(xs, Splus,  chi, N, beta)
        Sminus = _enforce_minorant(xs, Sminus, chi, N, beta)

    # Gaps (after optional enforcement)
    gap_plus_now  = Splus - chi
    gap_minus_now = chi   - Sminus
    okp = bool(np.min(gap_plus_now)  >= -1e-12)
    okm = bool(np.min(gap_minus_now) >= -1e-12)

    # L1 errors:
    #   - theory on circle: 1/(N+1) (for the unforced Selberg polynomials)
    #   - grid numerical after enforcement (may be >= theory)
    E_theory = 1.0/(N+1)
    dx = (xs[-1]-xs[0])/(len(xs)-1)
    E_plus_grid  = float(np.trapezoid(np.maximum(gap_plus_now,  0.0), xs))
    E_minus_grid = float(np.trapezoid(np.maximum(gap_minus_now, 0.0), xs))

    payload = {
      "model": "Selberg-circle",
      "params": {"beta": beta, "Delta": Delta, "N": N, "gridN": args.gridN, "force": bool(args.force)},
      "valid": {
        "majorant_ok": okp, "minorant_ok": okm,
        "min_gap_majorant": float(np.min(gap_plus_now)),
        "min_gap_minorant": float(np.min(gap_minus_now))
      },
      "L1_errors": {
        "E_plus_theory":  E_theory,
        "E_minus_theory": E_theory,
        "E_plus_grid":    E_plus_grid,
        "E_minus_grid":   E_minus_grid
      },
      "fourier": {
        "n_min": int(n[0]), "n_max": int(n[-1]),
        "a_plus":  a_plus.tolist(),
        "a_minus": a_minus.tolist()
      }
    }
    out_json = os.path.join(args.out, f"Selberg_circle_beta{beta:g}_D{Delta:g}.json")
    with open(out_json, "w", encoding="utf-8") as f: json.dump(payload, f, indent=2)
    print(out_json)

    if args.plot:
        import matplotlib.ticker as mt
        fig,ax = plt.subplots(figsize=(7.5,3.6))
        ax.plot(xs, chi,   lw=2, label="chi_{[-beta,beta]} (periodic)")
        ax.plot(xs, Splus, lw=1.5, label="S^+ (Selberg)")
        ax.plot(xs, Sminus,lw=1.5, label="S^- (Selberg)")
        ax.set_title(f"Selberg (circle): beta={beta}, N=floor(Delta)={N}, force={bool(args.force)}")
        ax.xaxis.set_major_locator(mt.MaxNLocator(11))
        ax.grid(alpha=0.25); ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, f"Selberg_circle_beta{beta:g}_D{Delta:g}.png"), dpi=140)
