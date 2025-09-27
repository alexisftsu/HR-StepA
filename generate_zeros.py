# generate_zeros.py
# Requiere: pip install mpmath

from mpmath import *
import time

# Configurar precisión
mp.dps = 50  # 50 decimales de precisión

print("Generando zeros de la función zeta de Riemann...")
print(f"Precisión: {mp.dps} decimales")

zeros = []
batch_size = 1000
max_zeros = 10000

start_time = time.time()

for n in range(1, max_zeros + 1):
    try:
        zero = zetazero(n)
        zeros.append(float(zero.imag))
        
        # Progreso cada 1000 zeros
        if n % batch_size == 0:
            elapsed = time.time() - start_time
            rate = n / elapsed
            eta = (max_zeros - n) / rate
            print(f"Calculados: {n}/{max_zeros} ({n/max_zeros*100:.1f}%) - ETA: {eta/60:.1f} min")
    
    except Exception as e:
        print(f"Error en zero {n}: {e}")
        continue

print(f"\nCalculados {len(zeros)} zeros en {(time.time() - start_time)/60:.1f} minutos")

# Guardar en archivo
output_file = "zeros_python_generated.dat"
with open(output_file, 'w') as f:
    for i, z in enumerate(zeros, 1):
        f.write(f"{z:.15f}\n")

print(f"Zeros guardados en: {output_file}")
print(f"Primer zero: {zeros[0]:.15f}")
print(f"Último zero: {zeros[-1]:.15f}")

# Verificación básica
print("\nVerificación de primeros 5 zeros:")
known_zeros = [14.134725141734693, 21.022039638771554, 25.010857580145688, 
               30.424876125859513, 32.935061587739189]

for i in range(min(5, len(zeros))):
    calculated = zeros[i]
    known = known_zeros[i] if i < len(known_zeros) else "N/A"
    if known != "N/A":
        diff = abs(calculated - known)
        print(f"Zero {i+1}: {calculated:.15f} (diff: {diff:.2e})")
    else:
        print(f"Zero {i+1}: {calculated:.15f}")