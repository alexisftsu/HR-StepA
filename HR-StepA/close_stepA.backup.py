import os, re, json, math, argparse, numpy as np
import mpmath as mp

ap = argparse.ArgumentParser()
ap.add_argument("--zeros", required=True)
ap.add_argument("--vk", default="vk_constants.json")
ap.add_argument("--out", required=True)
ap.add_argument("--scan-zld", action="store_true")
ap.add_argument("--Tscan", type=float, default=1e4)
ap.add_argument("--gridN", type=int, default=160)
ap.add_argument("--safety", type=float, default=2.5)
ap.add_argument("--xcap", type=float, default=None)
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

float_re = re.compile(r"[-+]?\d+(?:[.,]\d+)?")
def parse_gamma_from_line(line):
    s = line.strip()
    if not s or s.startswith("#") or s.startswith("//"): return None
    toks = float_re.findall(s)
    if not toks: return None
    cand=[]
    for t in toks:
        tt = t.replace(",", ".")
        try: cand.append(float(tt))
        except: pass
    if not cand: return None
    g = max([x for x in cand if x>0.0], default=None)
    if g is None or g<10.0: return None
    return g

def load_gammas(folder, Tcap=None):
    seen=set()
    for name in os.listdir(folder):
        path = os.path.join(folder, name)
        if not os.path.isfile(path): continue
        with open(path,"r",encoding="utf-8",errors="ignore") as f:
            for line in f:
                g = parse_gamma_from_line(line)
                if g is None: continue
                if Tcap is not None and g>Tcap: continue
                seen.add(round(g,12))
    return np.array(sorted(seen))

gammas_all = load_gammas(args.zeros, Tcap=None)
rho_abs = np.sqrt(0.25 + gammas_all**2) if gammas_all.size>0 else np.array([])
S1 = float(np.sum(1.0/rho_abs)) if gammas_all.size>0 else 0.0
C0p = S1 / (math.log(2.0)**2) if gammas_all.size>0 else 0.0
gmin = float(gammas_all[0]) if gammas_all.size else float("nan")
gmax = float(gammas_all[-1]) if gammas_all.size else float("nan")
T0 = gmax if np.isfinite(gmax) else 0.0
X0 = T0**2 if T0>0 else 0.0

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

T_min = math.sqrt(2.0)
K_eff = 10.0
K_max = None
if args.scan_zld:
    mp.mp.dps = 60
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
    def K_req(t):
        s = mp.mpf("0.5")+1j*mp.mpf(t)
        h = max(1e-6, float(t)*1e-6)
        sh = 1j*mp.mpf(h)
        try:
            fph=mp.log(mp.zeta(s+sh)); fmh=mp.log(mp.zeta(s-sh))
            ddt=(fph-fmh)/(2*mp.mpf(h))
            return float(abs(ddt)/(1.0+(math.log(t))**2))
        except: return None
    Ks=[]
    for t in t_eval:
        kv=K_req(t)
        if kv is not None and math.isfinite(kv): Ks.append(kv)
    if Ks:
        K_max=max(Ks); K_eff=args.safety*K_max

a=math.log(T_min)
I_tail=math.exp(-a)*(a*a+2.0*a+3.0)
A3 = (Cw/(2.0*math.pi*math.sqrt(2.0))) * (K_eff*I_tail)
C_R = max(0.0, A3 - 1.0/(4.0*math.pi))

B_VK=b_VK=x1=None; R=t0=None
if os.path.exists(args.vk):
    with open(args.vk,"r") as f: d=json.load(f)
    R=d.get("R",55.241); t0=d.get("t0",10.0)
    B_VK=d.get("B_VK",None); b_VK=d.get("b_VK",None); x1=d.get("x1",None)
if b_VK is None or B_VK is None:
    R=55.241 if R is None else R
    b_base=(1.0/(R+1.0))**(3.0/5.0)/8.0
    b_VK=b_base/3.0
    B_VK=150.0
if x1 is None: x1=1e6

def F_value(x,Bvk,bvk):
    L=math.log(x); L2=math.log(L)
    return Bvk*math.sqrt(x)/(L**2) * math.exp(- bvk*(L**(3.0/5.0))*(L2**(-1.0/5.0)))
def max_from_X0(X0,Bvk,bvk,r=1.05,steps=2000):
    if X0<=0: return None,None
    x=X0; best=F_value(x,Bvk,bvk); xbest=x
    for _ in range(steps):
        x*=r; val=F_value(x,Bvk,bvk)
        if val>best: best=val; xbest=x
    return best,xbest
C_alto, X_at_max = max_from_X0(X0,B_VK,b_VK)
F_X0 = F_value(X0,B_VK,b_VK) if X0>0 else None

C_bajo = 1.0/(4.0*math.pi) + C0p + C_R
C_empalme = max(C_bajo, F_X0) if F_X0 is not None else C_bajo
C_tot = max(C_bajo, (C_alto or 0.0), C_empalme) + 1e-12

out_json=os.path.join(args.out,"StepA_results.json")
with open(out_json,"w") as f:
    json.dump({
        "zeros":{"count":int(gammas_all.size),"gamma_min":gmin,"gamma_max":gmax,"S1":S1,"C0_prime_upper":C0p,"T0":T0,"X0":X0},
        "kernel":{"Cw_sup_h2":Cw,"T_min":T_min,"I_tail":I_tail,"K_eff":K_eff,"K_max_measured":K_max,"C_R":C_R},
        "VK":{"R":R,"t0":t0,"BVK":B_VK,"beta_VK":b_VK,"x1":x1},
        "constants":{"C_bajo":C_bajo,"C_alto":C_alto,"C_empalme":C_empalme,"C_tot":C_tot}
    }, f, indent=2)
print(out_json)



