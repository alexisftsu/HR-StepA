p=r".\canonical_idxgamma.txt"
with open(p,"r",encoding="utf-8-sig") as f:
    data=f.read()
with open(p,"w",encoding="utf-8") as f:
    f.write(data)
print("rewritten without BOM")
