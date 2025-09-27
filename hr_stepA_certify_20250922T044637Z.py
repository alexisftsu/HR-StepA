# HR Step A certification script (generated 20250922T044637Z UTC)
import math, json, os, numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime

T0 = 3000000000000.0
X0 = T0**2
C0_prime = 0.0
C_R = 0.0
B_VK = 10.0
b_VK = 0.01
x1 = 10.0
EPSILON = 1e-12

def C_bajo_value(C0p: float, CR: float) -> float:
    return 1.0/(4.0*math.pi) + C0p + CR

def F_value(x: float, Bvk: float, bvk: float) -> float:
    L = math.log(x); L2 = math.log(L)
    return Bvk * math.sqrt(x) / (L**2) * math.exp(- bvk * (L**(3.0/5.0)) * (L2**(-1.0/5.0)))

def maximize_F_from_X0(X0: float, Bvk: float, bvk: float, r: float = 1.05, steps: int = 2000):
    x_vals = [X0]; F_vals = [F_value(X0, Bvk, bvk)]
    decay_win = 8; ctr = 0
    for k in range(1, steps+1):
        x_vals.append(x_vals[-1] * r)
        F_vals.append(F_value(x_vals[-1], Bvk, bvk))
        if len(F_vals) >= decay_win+1 and all(F_vals[-j-1] <= 0.9*F_vals[-j-2] for j in range(decay_win)):
            ctr = k; break
    C_alto = max(F_vals); X_at_max = x_vals[int(np.argmax(F_vals))]
    return C_alto, x_vals, F_vals, ctr, X_at_max

def main(out_dir="/mnt/data"):
    C_bajo = C_bajo_value(C0_prime, C_R)
    F_X0   = F_value(X0, B_VK, b_VK)
    C_alto, x_vals, F_vals, steps_used, X_at_max = maximize_F_from_X0(X0, B_VK, b_VK)
    C_empalme = max(C_bajo, F_X0)
    C_tot = max(C_bajo, C_alto, C_empalme) + EPSILON
    timestamp = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    json_path = os.path.join(out_dir, f"HR_StepA_constants_{timestamp}.json")
    pdf_path  = os.path.join(out_dir, f"HR_StepA_report_{timestamp}.pdf")
    payload = {
        "generated_at_utc": timestamp,
        "T0": T0, "X0": X0,
        "C_bajo": C_bajo,
        "C0_prime": C0_prime, "C_R": C_R,
        "VK_region_constant": 55.241,
        "B_VK": B_VK, "b_VK": b_VK, "x1": x1,
        "F_X0": F_X0, "C_alto": C_alto, "C_empalme": C_empalme,
        "C_tot": C_tot, "epsilon": EPSILON,
        "grid_info": {"r": 1.05, "steps_used": steps_used, "X_at_max": X_at_max, "sampled_points": len(x_vals)},
        "placeholders_warning": "B_VK, b_VK, C0_prime, C_R are placeholders. Replace with certified values."
    }
    with open(json_path, "w") as f: json.dump(payload, f, indent=2)
    with PdfPages(pdf_path) as pdf:
        fig1 = plt.figure(figsize=(8.27, 11.69)); fig1.clf(); plt.axis("off")
        text = f"HR Project — Step A\n\nC_tot = max{C_bajo, C_alto, C_empalme} + ε\n\n" \               f"T0={T0:.3e}, X0={X0:.3e}\n" \               f"C_bajo={C_bajo:.12e}, F(X0)={F_X0:.12e}, C_alto={C_alto:.12e}, C_empalme={C_empalme:.12e}\n" \               f"ε={EPSILON:.1e}, C_tot={C_tot:.12e}\n"
        plt.text(0.03, 0.98, text, va="top", ha="left", family="monospace")
        pdf.savefig(fig1); plt.close(fig1)
        fig2 = plt.figure(figsize=(8.27, 5.0))
        plt.plot(x_vals, F_vals); plt.xscale("log"); plt.yscale("log")
        plt.xlabel("x (log scale)"); plt.ylabel("F(x)  [log scale]")
        plt.title("VK-based upper bound for F(x)"); plt.tight_layout(); pdf.savefig(fig2); plt.close(fig2)
    return json_path, pdf_path

if __name__ == "__main__":
    paths = main()
    print("Wrote:", paths)
