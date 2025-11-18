#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun
from scipy.stats import linregress
import re

# ----------------------------
# CONFIGURAÇÕES
# ----------------------------
SIM_DIR = "Simulations"
OUT_CSV = "deformation_potential.csv"
OUT_PNG_V = "deformation_potential_vs_volume.png"
OUT_PNG_S = "deformation_potential_vs_strain.png"

# ----------------------------
# Leitura de gaps
# ----------------------------
volumes = []
gap_direct = []
gap_indirect = []

print(f"Lendo pastas dentro de: {SIM_DIR}")

for folder in sorted(os.listdir(SIM_DIR)):
    m = re.match(r"scale_([0-9.]+)", folder)
    if not m:
        continue

    bands_dir = os.path.join(SIM_DIR, folder, "bands")
    vasprun_path = os.path.join(bands_dir, "vasprun.xml")
    outcar_path = os.path.join(bands_dir, "OUTCAR")

    if not os.path.isfile(vasprun_path) or not os.path.isfile(outcar_path):
        print(f"⚠️  Arquivos faltando em {folder}, pulando.")
        continue

    # Lê volume do OUTCAR
    V = None
    with open(outcar_path) as f:
        for line in f:
            if 'volume of cell' in line:
                V = float(line.split()[-1])
                break

    if V is None:
        print(f"⚠️  Volume não encontrado em {folder}, pulando.")
        continue

    # Lê vasprun.xml
    try:
        vasprun = Vasprun(vasprun_path, parse_projected_eigen=False)
        band_gap = vasprun.get_band_structure().get_band_gap()

        gap_ind = band_gap["energy"]
        is_direct = band_gap["direct"]

        gap_dir = gap_ind if is_direct else vasprun.get_band_structure().get_direct_band_gap()

        print(f"✓ {folder}: V = {V:.3f} Å³ | Gap indirect = {gap_ind:.6f} eV | Gap direct = {gap_dir:.6f} eV")

        volumes.append(V)
        gap_indirect.append(gap_ind)
        gap_direct.append(gap_dir)

    except Exception as e:
        print(f"⚠️  Erro lendo {vasprun_path}: {e}")
        continue

volumes = np.array(volumes)
gap_indirect = np.array(gap_indirect)
gap_direct = np.array(gap_direct)

# ----------------------------
# Salva CSV
# ----------------------------
df = pd.DataFrame({
    "Volume (Å³)": volumes,
    "Gap_indirect (eV)": gap_indirect,
    "Gap_direct (eV)": gap_direct
})
df.to_csv(OUT_CSV, index=False)
print(f"✅ CSV salvo: {OUT_CSV}")

# ----------------------------
# Leitura de V0 do ajuste Vinet (melhor prática)
# ----------------------------
fit_vinet_path = os.path.join("fit-EOS-Vinet.out")
if not os.path.exists(fit_vinet_path):
    raise FileNotFoundError(f"Arquivo {fit_vinet_path} não encontrado!")

with open(fit_vinet_path) as f:
    texto = f.read()

match = re.search(r'V0 = ([\d.]+) Å³', texto)
if match:
    V0 = float(match.group(1))
    print(f"🔹 V₀ (referência) = {V0:.3f} Å³")
else:
    raise ValueError("❌ V0 não encontrado no arquivo fit-EOS-Vinet.out")

strain = (volumes - V0) / V0

# ----------------------------
# Ajuste linear em Volume
# ----------------------------
slope_V_ind, intercept_V_ind, r_V_ind, _, _ = linregress(volumes, gap_indirect)
slope_V_dir, intercept_V_dir, r_V_dir, _, _ = linregress(volumes, gap_direct)

# ----------------------------
# Ajuste linear em Strain
# ----------------------------
slope_S_ind, intercept_S_ind, r_S_ind, _, _ = linregress(strain, gap_indirect)
slope_S_dir, intercept_S_dir, r_S_dir, _, _ = linregress(strain, gap_direct)

# ----------------------------
# Prints
# ----------------------------
print("\n=== Deformation potential (V) ===")
print(f"Indirect gap: slope = {slope_V_ind:.6f} eV/Å³, R² = {r_V_ind**2:.4f}")
print(f"Direct gap:   slope = {slope_V_dir:.6f} eV/Å³, R² = {r_V_dir**2:.4f}")

print("\n=== Deformation potential (strain) ===")
print(f"Indirect gap: slope = {slope_S_ind:.6f} eV, R² = {r_S_ind**2:.4f}")
print(f"Direct gap:   slope = {slope_S_dir:.6f} eV, R² = {r_S_dir**2:.4f}")

# ----------------------------
# Plot: E_gap vs Volume
# ----------------------------
plt.figure(figsize=(8,6))
plt.scatter(volumes, gap_indirect, label="Indirect gap", marker="o")
plt.scatter(volumes, gap_direct, label="Direct gap", marker="s")
Vfit = np.linspace(min(volumes), max(volumes), 100)
plt.plot(Vfit, slope_V_ind * Vfit + intercept_V_ind, 'b--', label=f'Indirect fit: slope={slope_V_ind:.4f}')
plt.plot(Vfit, slope_V_dir * Vfit + intercept_V_dir, 'g--', label=f'Direct fit: slope={slope_V_dir:.4f}')
plt.xlabel("Volume (Å³)")
plt.ylabel("Band gap (eV)")
plt.title("Deformation Potential: E_gap vs Volume")
plt.grid(True, ls='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(OUT_PNG_V, dpi=300)
print(f"✅ Gráfico salvo: {OUT_PNG_V}")
plt.show()

# ----------------------------
# Plot: E_gap vs Strain
# ----------------------------
plt.figure(figsize=(8,6))
plt.scatter(strain, gap_indirect, label="Indirect gap", marker="o")
plt.scatter(strain, gap_direct, label="Direct gap", marker="s")
S_fit = np.linspace(min(strain), max(strain), 100)
plt.plot(S_fit, slope_S_ind * S_fit + intercept_S_ind, 'b--', label=f'Indirect fit: slope={slope_S_ind:.4f}')
plt.plot(S_fit, slope_S_dir * S_fit + intercept_S_dir, 'g--', label=f'Direct fit: slope={slope_S_dir:.4f}')
plt.xlabel("Strain ε = (V - V₀)/V₀")
plt.ylabel("Band gap (eV)")
plt.title("Deformation Potential: E_gap vs Strain")
plt.grid(True, ls='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(OUT_PNG_S, dpi=300)
print(f"✅ Gráfico salvo: {OUT_PNG_S}")
plt.show()

# ----------------------------
# Info
# ----------------------------
print(f"\n🔹 V₀ (referência) = {V0:.3f} Å³")
print(f"Strain range: {strain.min():.4f} → {strain.max():.4f}")
