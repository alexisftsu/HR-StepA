# -*- coding: utf-8 -*-
# explicit_compare_policy.py
# Igual que explicit_compare_npyfree.py, con modo T_mode="policy"

import argparse, math

def read_gammas(path):
    gs=[]
    with open(path,"r",encoding="utf-8-sig") as f:
        for s in f:
            s=s.strip()
            if s: gs.append(float(s))
    gs.sort()
    return gs

def psi_exact_upto(xmax):
    spf = list(range(xmax+1))
    for i in range(2,int(xmax**0.5)+1):
        if spf[i]==i:
            for j in range(i*i, xmax+1, i):
                if spf[j]==j: spf[j]=i
    logp=[0.0]*(xmax+1)
    for i in range(2,xmax+1):
        if spf[i]==i: logp[i]=math.log(i)
    lam=[0.0]*(xmax+1)
    for n in range(2,xmax+1):
        p=spf[n]; m=n
        while m%p==0: m//=p
        if m==1: lam[n]=logp[p]
    psi=[0.0]*(xmax+1); run=0.0
    for n in range(1,xmax+1):
        run+=lam[n]; psi[n]=run
    return psi

def psi_explicit_truncated(x, gammas, T):
    L = math.log(x)
    sq = math.sqrt(x)
    S = 0.0
    for g in gammas:
        if g>T: break
        den = 0.25 + g*g
        cg = math.cos(g*L)
        sg = math.sin(g*L)
        S += 2.0*sq*((0.5*cg + g*sg)/den)
    tail = math.log(2.0*math.pi) + 0.5*math.log(1.0 - x**-2)
    return x - S - tail

def make_x_points(xmax, k, xmin=2):
    xmin = max(2, int(xmin))
    if k<=1: return [max(xmin, xmax)]
    xs=[]
    L0, L1 = math.log(float(xmin)), math.log(float(xmax))
    for i in range(k):
        t=i/(k-1)
        x=int(round(math.exp(L0 + t*(L1-L0))))
        x=max(xmin, min(x, xmax))
        if not xs or x!=xs[-1]:
            xs.append(x)
    return xs

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--xmax", type=int, required=True)
    ap.add_argument("--points", type=int, default=40)
    ap.add_argument("--zeros", required=True)
    ap.add_argument("--T_mode", choices=["constant","sqrtx","policy"], default="constant")
    ap.add_argument("--T_const", type=float, default=1000.0)
    ap.add_argument("--T_min", type=float, default=100.0, help="m?nimo T si usas sqrtx")
    ap.add_argument("--xmin", type=int, default=2, help="x m?nimo a evaluar")
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    gammas=read_gammas(args.zeros)
    xmax=args.xmax
    xs=make_x_points(xmax, args.points, args.xmin)
    psi=psi_exact_upto(xmax)

    with open(args.out,"w",encoding="utf-8") as f:
        f.write("x,psi_exact,psi_explicit,remainder psi_exact - explicit,ratio = |remainder|/(sqrt(x)*log(x)^2),T_used\n")
        for x in xs:
            # Selecci?n de T
            if args.T_mode=="constant":
                T = args.T_const
            elif args.T_mode=="sqrtx":
                T = max(args.T_min, math.sqrt(x))
            else:  # policy
                if x < 500:
                    T = max(700.0, math.sqrt(x))
                elif x < 3000:
                    T = 2000.0
                else:
                    T = 5000.0
            if gammas and T>gammas[-1]:
                T=gammas[-1]

            pe = psi[x]
            px = psi_explicit_truncated(x, gammas, T)
            rem = pe - px
            denom = math.sqrt(x)*(math.log(x)**2)
            ratio = abs(rem)/denom if denom>0 else 0.0
            f.write(f"{x},{pe:.12f},{px:.12f},{rem:.12f},{ratio:.12e},{T:.0f}\n")

if __name__=="__main__":
    main()
