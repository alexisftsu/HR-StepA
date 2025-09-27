# explicit_adaptive_target.py
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
    spf=list(range(xmax+1))
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

def make_x_points(xmax, k, xmin=2):
    xmin=max(2,int(xmin))
    if k<=1: return [max(xmin,xmax)]
    xs=[]; L0, L1 = math.log(float(xmin)), math.log(float(xmax))
    for i in range(k):
        t=i/(k-1)
        x=int(round(math.exp(L0 + t*(L1-L0))))
        x=max(xmin, min(x,xmax))
        if not xs or x!=xs[-1]: xs.append(x)
    return xs

def explicit_terms_for_x(x, gammas):
    L=math.log(x); sq=math.sqrt(x)
    terms=[]  # T_i = 2*sq*((0.5*cos(gL)+g*sin(gL))/(0.25+g^2))
    for g in gammas:
        den=0.25+g*g
        cg=math.cos(g*L); sg=math.sin(g*L)
        terms.append( 2.0*sq*((0.5*cg + g*sg)/den) )
    # prefijos S_k
    pref=[0.0]*len(terms); acc=0.0
    for i,t in enumerate(terms):
        acc+=t; pref[i]=acc
    return pref

def psi_explicit_from_prefix(x, gammas, pref, k):
    # usa todas las gammas con ?ndice <= k
    tail = math.log(2.0*math.pi) + 0.5*math.log(1.0 - x**-2)
    S = pref[k] if k>=0 else 0.0
    return x - S - tail

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--xmax", type=int, required=True)
    ap.add_argument("--points", type=int, default=200)
    ap.add_argument("--xmin", type=int, default=100)
    ap.add_argument("--zeros", required=True)
    ap.add_argument("--target", type=float, default=1e-3, help="umbral para el ratio")
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    gammas=read_gammas(args.zeros)
    psi=psi_exact_upto(args.xmax)
    xs=make_x_points(args.xmax, args.points, args.xmin)

    with open(args.out,"w",encoding="utf-8") as f:
        f.write("x,psi_exact,psi_explicit,remainder psi_exact - explicit,ratio = |remainder|/(sqrt(x)*log(x)^2),T_chosen,k\n")
        for x in xs:
            pe=psi[x]
            pref=explicit_terms_for_x(x, gammas)
            # b?squeda secuencial/binary sobre k (0..len(gammas)-1)
            lo, hi = -1, len(gammas)-1
            best_k = hi
            while lo < hi:
                mid = (lo+hi)//2
                px = psi_explicit_from_prefix(x, gammas, pref, mid)
                rem = pe - px
                denom = math.sqrt(x)*(math.log(x)**2)
                ratio = abs(rem)/denom if denom>0 else 0.0
                if ratio <= args.target:
                    best_k = mid; hi = mid
                else:
                    lo = mid+1
            # resultado con best_k
            px = psi_explicit_from_prefix(x, gammas, pref, best_k)
            rem = pe - px
            denom = math.sqrt(x)*(math.log(x)**2)
            ratio = abs(rem)/denom if denom>0 else 0.0
            T = gammas[best_k] if best_k>=0 else 0.0
            f.write(f"{x},{pe:.12f},{px:.12f},{rem:.12f},{ratio:.12e},{T:.6f},{best_k+1}\n")

if __name__=="__main__":
    main()
