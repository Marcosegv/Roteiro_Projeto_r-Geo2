#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os
import re
from lmfit import Model

# -------------------- CONFIGURAÇÕES --------------------
SIM_DIR = "Simulations_fit"
OUT_PNG = "eos_fit_vinet_lmfit.png"
POSCAR = "arquivos/POSCAR"
# ------------------------------------------------------

# Modelo Vinet
def vinet(V, E0, V0, B0, Bp):
    eta = (V / V0) ** (1/3)
    return E0 + 9 * V0 * B0 / 16 * (((1 - eta) ** 3) * Bp + (1 - eta) ** 2 * (6 - 4 * eta))

# Função para ler volume do POSCAR
def ler_volume_poscar(poscar_path):
    with open(poscar_path) as f:
        lines = f.readlines()
    scale = float(lines[1])
    a = np.array([float(x) for x in lines[2].split()])
    b = np.array([float(x) for x in lines[3].split()])
    c = np.array([float(x) for x in lines[4].split()])
    volume = np.abs(np.dot(a, np.cross(b, c))) * scale**3
    return volume

# Função para ler volumes e energias
def ler_dados(sim_dir, V0):
    volumes = []
    energies = []
    for folder in sorted(os.listdir(sim_dir)):
        m = re.match(r"scale_([0-9.]+)", folder)
        if not m:
            continue
        factor = float(m.group(1))
        V = V0 * factor**3
        outcar = os.path.join(sim_dir, folder, "OUTCAR")
        E = None
        if os.path.exists(outcar):
            with open(outcar) as f:
                for line in reversed(f.readlines()):
                    if "free  energy   TOTEN" in line:
                        E = float(line.split()[-2])
                        break
        if E is None:
            print(f"⚠️  Energia não encontrada em {sim_dir}/{folder}")
            continue
        volumes.append(V)
        energies.append(E)
        print(f"✓ {folder}: V = {V:.3f} Å³ | E = {E:.6f} eV")
    return np.array(volumes), np.array(energies)

# Ajuste com lmfit
def ajustar_vinet(volumes, energies):
    model = Model(vinet)
    params = model.make_params(
        E0=min(energies),
        V0=volumes[np.argmin(energies)],
        B0=0.5,    # GPa
        Bp=4.0
    )
    result = model.fit(energies, params, V=volumes)
    return result

# ---------- MAIN ----------
if not os.path.exists(POSCAR):
    print(f"❌ Arquivo POSCAR não encontrado na pasta atual.")
    exit(1)

V0 = ler_volume_poscar(POSCAR)
print(f"🔹 V0 (volume do POSCAR) = {V0:.6f} Å³")

if not os.path.exists(SIM_DIR):
    print(f"❌ Pasta {SIM_DIR} não encontrada.")
    exit(1)

volumes, energies = ler_dados(SIM_DIR, V0)
if len(volumes) == 0:
    print(f"⚠️  Nenhum dado encontrado para ajuste.")
    exit(1)

result = ajustar_vinet(volumes, energies)

E0 = result.params['E0'].value
V0_fit = result.params['V0'].value
B0 = result.params['B0'].value
Bp = result.params['Bp'].value

# Conversão para GPa
B0_GPa = B0 * 160.21766208  # eV/Å³ → GPa

print("\n=== Parâmetros ajustados ===")
print(f"E0 = {E0:.6f} eV")
print(f"V0 = {V0_fit:.3f} Å³")
print(f"B0 = {B0_GPa:.3f} GPa")
print(f"B' = {Bp:.3f}")

# Plot
Vfit = np.linspace(min(volumes)*0.98, max(volumes)*1.02, 300)
Efit = vinet(Vfit, E0, V0_fit, B0, Bp)

plt.figure(figsize=(8,6))
plt.plot(Vfit, Efit, 'r-', label="Ajuste Vinet (lmfit)")
plt.scatter(volumes, energies, color='b', label="Dados VASP")
plt.xlabel("Volume (Å³)")
plt.ylabel("Energia (eV)")
plt.title("EOS Fit: Vinet (com B')")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
print(f"✅ Gráfico salvo como: {OUT_PNG}")
plt.show()
