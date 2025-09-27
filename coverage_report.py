import math
with open(r".\canonical_idxgamma.txt","r",encoding="utf-8-sig") as f:
    gs=[float(s) for s in f if s.strip()]
def N(T): q=T/(2*math.pi); return q*math.log(q)-q+0.875
for T in (1000,3162):
    have=sum(1 for g in gs if g<=T); exp=N(T)
    print(f"T={T:5.0f}  have={have:6d}  RvM≈{exp:9.1f}  gap={have-exp:8.1f}")
