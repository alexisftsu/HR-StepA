\# Proyecto global HR — Estado y reproducibilidad (Policy E)



\## TL;DR

\- Fórmula explícita truncada para ψ(x) con política \*\*E\*\* de recorte:

&nbsp; - \*\*low\*\*: `x < 300` → `T = max(2000, √x)` (en la práctica 2000)

&nbsp; - \*\*mid\*\*: `300 ≤ x < 3000` → `T = 3000`

&nbsp; - \*\*high\*\*: `x ≥ 3000` → `T = 5000`

&nbsp; - (capado por disponibilidad de ceros: `T ≤ 10000` de `canonical\_idxgamma\_T10000.txt`)

\- Corrida de referencia: \*\*500 puntos\*\* log-espaciados entre `x=100` y `x=10^6`.

\- Métricas (archivo `policy\_E\_p500.csv`):

&nbsp; - \*\*Points\*\*: 500  

&nbsp; - \*\*MaxRatio\*\*: 9.108e-03 @ x=127  

&nbsp; - \*\*Median\*\*: 1.170e-04  

&nbsp; - \*\*P90\*\*: 4.081e-04  

&nbsp; - \*\*P95\*\*: 9.478e-04  

&nbsp; - \*\*Mean\*\*: 3.291e-04  

&nbsp; - \*\*<5e-3\*\*: 492 / 500  

&nbsp; - \*\*<1e-3\*\*: 475 / 500  

\- `T\_used` \*\*verificado\*\* contra la política: ✓ 500/500 filas.



---



\## Cómo reproducir



\### 0) Prerrequisitos

\- `canonical\_idxgamma\_T10000.txt` (ceros positivos hasta T≤10000, UTF-8 sin BOM).

\- Scripts parcheados para escribir `T\_used` en el CSV:

&nbsp; - `explicit\_compare\_policy\_param.py`



\### 1) Ejecutar Policy E con 500 puntos

```powershell

python .\\explicit\_compare\_policy\_param.py --xmax 1000000 --points 500 --xmin 100 --zeros .\\canonical\_idxgamma\_T10000.txt --b1 300 --b2 3000 --Tmin\_low 2000 --T\_mid 3000 --T\_high 5000 --out .\\policy\_E\_p500.csv



