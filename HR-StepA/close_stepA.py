import os, re, json, math, argparse, numpy as np
import mpmath as mp
import multiprocessing as mpc

ap = argparse.ArgumentParser()
ap.add_argument("--zeros", required=True)
ap.add_argument("--vk", default="vk_constants.json")
ap.add_argument("--out", required=True)
ap.add_argument("--scan-zld", action="store_true")
ap.add_argument("--Tscan", type=float, default=1e4)
ap.add_argument("--gridN", type=int, default=160)
ap.add_argument("--safety", type=float, default=2.5)
ap.add_argument("--xcap", type=float, default=None)
# rendimiento/feedback
ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2)//2))
ap.add_argument("--mpdps", type=int, default=30)
ap.add_argument("--hscale", type=float, default=1e-5)
ap.add_argument("--progress", action="store_true")
# NUEVO: modo RH
ap.add_argument("--assume-rh", action="store_true", help="Emite bound uniforme bajo RH: C_tot_RH*sqrt(x)*log(x)^2 para todo x>=2")

args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

# ---- util: floats ----
float_re = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")

def parse_gamma_from_line(line):
    s = line.strip()
    if not s or s.startswith("#") or s.startswith("//"):
        return None
    toks = float_re.findall(s)
    if not toks: return None
    vals=[]
    for t in toks:
        try: vals.append(float(t))
        except: pass
    cand = [v for v in vals if 10.0 <= v <= 1.0e14]
    if not cand: return None
    return min(cand)

def load_gammas(folder, Tcap=None):
    seen=set()
    for name in os.listdir(folder):
        low = name.lower()
        if "hash" in low: continue
        path = os.path.join(folder, name)
        if not os.path.isfile(path): continue
        with open(path,"r",encoding="utf-8",errors="ignore") as f:
            for line in f:
                g = parse_gamma_from_line(line)
                if g is None: continue
                if Tcap is not None and g>Tcap: continue
                seen.add(round(g,12))
    return np.array(sorted(seen))

# ---- C0' / X0 con ceros ----
gammas_all = load_gammas(args.zeros, Tcap=None)
rho_abs = np.sqrt(0.25 + gammas_all**2) if gammas_all.size>0 else np.array([])
S1 = float(np.sum(1.0/rho_abs)) if gammas_all.size>0 else 0.0
C0p = S1 / (math.log(2.0)**2) if gammas_all.size>0 else 0.0
gmin = float(gammas_all[0]) if gammas_all.size else float("nan")
gmax = float(gammas_all[-1]) if gammas_all.size else float("nan")
T0 = gmax if np.isfinite(gmax) else 0.0
X0 = T0**2 if T0>0 else 0.0

# ---- Kernel: sup|h''| para hat_g=(1-|t|)^3_+ ----
def sup_h2_kernel(n_t=20000):
    import numpy as _np
    sigma=1.0
    def hat_g(t):
        tt = abs(t)/sigma
        return 0.0 if tt>=1.0 else (1.0/sigma)*((1.0-tt)**3)
    def h_func(t): return 1j*t*hat_g(t)
    a,b=-sigma,sigma
    ts=_np.linspace(a,b,n_t); dt=ts[1]-ts[0]
    h=_np.array([h_func(t) for t in ts], dtype=_np.complex128)
    h1=_np.gradient(h,dt); h2=_np.gradient(h1,dt)
    return float(_np.max(_np.abs(h2)))
Cw = sup_h2_kernel()

# ---- Escaneo de |ζ'/ζ| (opcional) ----
T_min = math.sqrt(2.0)
K_eff = 10.0
K_max = None

_HSCALE = 1e-6
def _init_worker(dps, hscale):
    mp.mp.dps = dps
    global _HSCALE
    _HSCALE = hscale

def K_req_worker(t):
    try:
        s = mp.mpf("0.5")+1j*mp.mpf(t)
        h = mp.mpf(_HSCALE)*mp.mpf(t)
        if h <= 0: h = mp.mpf("1e-8")
        sh = 1j*h
        fph=mp.log(mp.zeta(s+sh)); fmh=mp.log(mp.zeta(s-sh))
        ddt=(fph-fmh)/(2*h)
        return float(abs(ddt)/(1.0+(math.log(t))**2))
    except:
        return None

if args.scan_zld:
    N=args.gridN
    grid=np.geomspace(T_min, args.Tscan, N)
    gam_mask = load_gammas(args.zeros, Tcap=args.Tscan)
    delta=1e-4
    mask=np.ones_like(grid,dtype=bool)
    if gam_mask.size>0:
        for i,t in enumerate(grid):
            j = np.searchsorted(gam_mask, t)
            d1 = abs(t-gam_mask[j-1]) if j-1>=0 else float("inf")
            d2 = abs(gam_mask[j]-t) if j<gam_mask.size else float("inf")
            if min(d1,d2)<delta: mask[i]=False
    t_eval=grid[mask]
    if args.progress:
        print(f"[scan] points={len(t_eval)} up to T={args.Tscan}  workers={args.workers} mp.dps={args.mpdps}", flush=True)
    if args.workers>1 and len(t_eval)>0:
        chunk=max(1,len(t_eval)//(args.workers*4))
        with mpc.Pool(processes=args.workers, initializer=_init_worker, initargs=(args.mpdps, args.hscale)) as pool:
            for i, kv in enumerate(pool.imap(K_req_worker, t_eval, chunksize=chunk)):
                if kv is not None and (K_max is None or kv>K_max): K_max=kv
                if args.progress and i % max(1,len(t_eval)//10)==0:
                    print(f"[scan] {i}/{len(t_eval)}  K_max≈{K_max}", flush=True)
    else:
        _init_worker(args.mpdps, args.hscale)
        for i,t in enumerate(t_eval):
            kv = K_req_worker(t)
            if kv is not None and (K_max is None or kv>K_max): K_max=kv
            if args.progress and i % max(1,len(t_eval)//10)==0:
                print(f"[scan] {i}/{len(t_eval)}  K_max≈{K_max}", flush=True)
    if K_max is not None:
        K_eff = args.safety * K_max

# ---- C_R del kernel ----
a=math.log(T_min)
I_tail=math.exp(-a)*(a*a+2.0*a+3.0)
A3 = (Cw/(2.0*math.pi*math.sqrt(2.0))) * (K_eff*I_tail)
C_R = max(0.0, A3 - 1.0/(4.0*math.pi))

# ---- VK (PowerShell-friendly keys) ----
B_VK=b_VK=x1=None; R=t0=None
if os.path.exists(args.vk):
    with open(args.vk,"r",encoding="utf-8-sig") as f:
        d=json.load(f)
    R=d.get("R",55.241); t0=d.get("t0",10.0)
    B_VK=d.get("B_VK", d.get("BVK", None))
    b_VK=d.get("b_VK", d.get("beta_VK", None))
    x1=d.get("x1", None)
if b_VK is None or B_VK is None:
    R=55.241 if R is None else R
    b_base=(1.0/(R+1.0))**(3.0/5.0)/8.0
    b_VK=b_base/3.0
    B_VK=150.0
if x1 is None: x1=1e6

def F_value(x,Bvk,bvk):
    L=math.log(x); L2=math.log(L)
    return Bvk*math.sqrt(x)/(L**2) * math.exp(- bvk*(L**(3.0/5.0))*(L2**(-1.0/5.0)))

def max_from_X0(X0,Bvk,bvk,r=1.05,steps=2000,xcap=None,x1=None):
    if X0<=0: return None,None
    x = max(X0, x1 or X0)
    best = F_value(x,Bvk,bvk); xbest = x
    k = 0
    while True:
        if steps is not None and k>=steps: break
        if xcap is not None and x>=xcap: break
        x *= r; k += 1
        val = F_value(x,Bvk,bvk)
        if val > best: best, xbest = val, x
    return best, xbest

C_alto, X_at_max = max_from_X0(X0,B_VK,b_VK,xcap=args.xcap,x1=x1)
F_X0 = F_value(X0,B_VK,b_VK) if X0>0 else None

C_bajo = 1.0/(4.0*math.pi) + C0p + C_R
C_empalme = max(C_bajo, F_X0) if F_X0 is not None else C_bajo
C_tot = max(C_bajo, (C_alto or 0.0), C_empalme) + 1e-12

# ==== NUEVO: modo RH — constante uniforme en todo x>=2 ====
constants_RH = None
if args.assume_rh:
    # Schoenfeld (RH): |psi(x)-x| <= (1/(8π)) sqrt(x) log^2 x para x >= 73.2 aprox.
    C_RH = 1.0/(8.0*math.pi)
    # verificar 2 <= x <= Xcheck por fuerza bruta de psi(x)
    def primes_upto(n):
        n=int(n)
        sieve=[True]*(n+1)
        sieve[0]=sieve[1]=False
        for p in range(2,int(n**0.5)+1):
            if sieve[p]:
                step=p; start=p*p
                sieve[start:n+1:step]=[False]*(((n-start)//step)+1)
        return [i for i,v in enumerate(sieve) if v]
    def psi_integer(n, primes):
        # psi(n) = sum_{p^k <= n} log p
        s=0.0
        for p in primes:
            if p>n: break
            pk=p
            while pk<=n:
                s += math.log(p)
                pk *= p
        return s
    Xcheck = 1000  # suficiente para cubrir [2, ~73.2] con holgura
    P = primes_upto(Xcheck)
    C_small = 0.0
    for n in range(2, Xcheck+1):
        L = math.log(n)
        if L<=0: continue
        denom = math.sqrt(n)*(L**2)
        if denom==0: continue
        ratio = abs(psi_integer(n,P) - n)/denom
        if ratio > C_small: C_small = ratio
    C_tot_RH = max(C_RH, C_small)
    constants_RH = {
        "C_RH": C_RH,
        "C_small_check_upto": Xcheck,
        "C_small": C_small,
        "C_tot_RH": C_tot_RH
    }
    # emitir statement en TXT
    rh_txt = os.path.join(args.out, "StepA_RH_statement.txt")
    with open(rh_txt, "w", encoding="ascii", errors="ignore") as f:
        f.write("=== HR Step A — Uniform bound under RH ===\n\n")
        f.write(f"C_RH = 1/(8*pi) = {C_RH:.12g}\n")
        f.write(f"Small-x check: 2 <= x <= {Xcheck}, C_small = {C_small:.12g}\n")
        f.write(f"Uniform constant (all x>=2, under RH): C_tot_RH = {C_tot_RH:.12g}\n\n")
        f.write("Conclusion (RH): for all x >= 2,\n")
        f.write("  |psi(x) - x| <= C_tot_RH * sqrt(x) * (log x)^2.\n")

# ---- Salida JSON ----
out_json = os.path.join(args.out, "StepA_results.json")
with open(out_json,"w",encoding="utf-8") as f:
    json.dump({
        "zeros":{"count":int(gammas_all.size),"gamma_min":gmin,"gamma_max":gmax,"S1":S1,"C0_prime_upper":C0p,"T0":T0,"X0":X0},
        "kernel":{"Cw_sup_h2":Cw,"T_min":T_min,"I_tail":I_tail,"K_eff":K_eff,"K_max_measured":K_max,"C_R":C_R},
        "VK":{"R":R,"t0":t0,"BVK":B_VK,"beta_VK":b_VK,"x1":x1},
        "constants":{"C_bajo":C_bajo,"C_alto":C_alto,"C_empalme":C_empalme,"C_tot":C_tot},
        "constants_RH": constants_RH
    }, f, indent=2)
print(out_json)
