gs=[]
with open(r".\canonical_idxgamma.txt","r",encoding="utf-8-sig") as f:
    for s in f:
        s=s.strip()
        if s: gs.append(float(s))
print("count=",len(gs),"strict_inc=",all(gs[i]<gs[i+1] for i in range(len(gs)-1)))
