# Roteiro Projeto r-GeO2 - Parte I
*** 
## 0. Organização geral do projeto
Estrutura sugerida do diretório do projeto:
```bash
projeto/
 ├── 0.INPUTS_BASE/
 │    ├── POSCAR
 │    ├── POTCAR
 │    ├── KPOINTS
 │    ├── KPOINTS-bands
 │    ├── INCAR-opt_rede
 │    ├── INCAR-opt_completa
 │    ├── INCAR-bands_scf
 │    ├── INCAR-bands_nscf
 │    ├── INCAR-optics
 │    └──minushalf_potfiles/
 │        ├── POTCAR.ge
 │        └── POTCAR.o
 ├── LDA/
 ├── PBE/
 └── PBEsol/
```
> Cada aluno escolhe uma pasta funcional (LDA, PBE ou PBEsol) e fazer tudo lá dentro, sempre usando os arquivos base de `0.INPUTS_BASE`, alterando apenas o necessário para seu sistema.

## 1. Otimização
Entrar na pasta do projeto e criar o diretório do funcional:
```bash
cd /home/ff281/FF281-2025/Geo2-Projeto_Final

# Ex.: se você ficou com PBEsol
mkdir PBEsol
cd PBEsol
```

Criar uma cópia local da pasta de inputs base (opcional, mas deixa tudo organizado):
```bash
cp -r ../0.INPUTS_BASE ./inputs
```
> A partir daqui, sempre que for preciso um arquivo base, usar `./inputs/...` ou `../0.INPUTS_BASE/...` (como você preferir).

### 1.1. Otimização de rede (volume)
```bash
mkdir 1.optimize_rede
cd 1.optimize_rede

cp ../inputs/POSCAR .
cp ../inputs/POTCAR .
cp ../inputs/KPOINTS .
cp ../inputs/INCAR-opt_rede INCAR
```

Rode o vasp via slurm. Para isso, crie o arquivo `run.srm`:
```bash
#!/bin/bash
#SBATCH --job-name=vasp      # Escolha um nome para o Job!
#SBATCH --partition=local    # Nome da fila
#SBATCH --nodes=1            # nó einstein
#SBATCH --ntasks=8           # 8 processos MPI; se quiser rodar com menos cores altere aqui!
#SBATCH --time=7-00:00:00    # Tempo limite do cálculo rodando, opcional!

# se for puro MPI:
export OMP_NUM_THREADS=1

# aqui vai o mpirun que você já está acostumado
mpirun -np $SLURM_NTASKS vasp.6.5.1_std
```
Em seguida rode:
```bash
# Veja se tem recursos disponíveis
sinfo

# Rode o vasp:
sbatch run.srm

# Acompanhe o cálculo na fila ou rodando
squeue
```
> **Atenção: Não se deve rodar contas fora do slurm, ou seja, nada de nohup ou executar o executável diretamente! O slurm é responsável por fazer o balanceamento de carga da máquina, assim como gerenciar a fila de jobs! Se não tiver recursos disponíveis o job ficaria na fila por ordem de prioridade!

**Objetivo desta etapa:** encontrar o volume/lado de rede otimizado para o funcional do aluno.

**Quando terminar:**
Usaremos o CONTCAR desta etapa como ponto de partida para a próxima.

### 1.2. Otimização completa da estrutura
```bash
cd ..
mkdir 2.optimize_completa
cd 2.optimize_completa

cp ../1.optimize_rede/CONTCAR POSCAR
cp ../inputs/POTCAR .
cp ../inputs/KPOINTS .
cp ../inputs/INCAR-opt_completa INCAR
cp ../1.optimize_rede/run.srm .
```

Rode o vasp via slurm:
```bash
sbatch run.srm
```

**Objetivo:** relaxar todos os graus de liberdade (átomos + eventualmente forma da célula, se for o caso) com o volume inicial vindo da etapa anterior.

## 2. DFT - 1/2

Crie uma pasta para rodar o dft-1/2. Em seguida, é necessário preparar os arquivos abaixo:

Arquivos necessários:
- `INCAR`
- `KPOINTS`
- `POSCAR`
- `POTCAR`
- `minushalf.yaml`
- `minushalf_potfiles/POTCAR.ge`
- `minushalf_potfiles/POTCAR.o`

Arquivo `INCAR`:
```bash
SYSTEM = r-GeO2 ENCUT convergence test
ISTART = 0
ICHARG = 2
PREC   = Accurate
EDIFF  = 1E-8
ENCUT = 520.0
NSW    = 0
IBRION = -1
ISMEAR = 0
SIGMA  = 0.05
LREAL=.FALSE.
ALGO   = Normal
LWAVE  = .FALSE.
LCHARG = .FALSE.
LASPH  = .TRUE.
NPAR   = 4
LORBIT = 11
```

Arquivo `minushalf.yaml`:
```bash
software: VASP
vasp:
    command: ['mpirun','-np','4','vasp.6.5.1_std']

correction:
    correction_code: v
```

Basta rodar o script abaixo (`sbatch run-dft-1_2`):
```bash
#!/bin/bash
#SBATCH --job-name=minushalf
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=1-00:00:00

export OMP_NUM_THREADS=1

# Ativando o ambiente do minushalf
conda activate minushalf_env

# Executa o minushalf (que internamente chama o VASP várias vezes)
minushalf execute > minushalf.out 2>&1
```

## 3. Variação de volume (RASCUNHO)

```bash
cd ../ ; mkdir volume ; mkdir volume/dft 

cd volume/dft
mkdir arquivos
cp ../2.optimize_completa/INCAR arquivos/INCAR_rlx
cp ../2.optimize_completa/POSCAR arquivos/.
cp ../2.optimize_completa/POTCAR arquivos/.
cp ../2.optimize_completa/KPOINTS arquivos/.

cp ~/Documentos/Mestrado/vasp/scripts/propriedades/2.job-EOS_inputs.srm .
cp ~/Documentos/Mestrado/vasp/scripts/propriedades/4.run_relax_dos_bands.srm .

~/Documentos/Mestrado/1.vasp/scripts/2.propriedades/1.generate_EOS_inputs.sh > saida-generate_EOS_inputs.out

cp ~/Documentos/Mestrado/1.vasp/scripts/2.propriedades/2.job-EOS_inputs.srm .
python ~/Documentos/Mestrado/vasp/scripts/propriedades/3.fit_eos_vinet.py > fit-EOS-Vinet.out
python ~/Documentos/Mestrado/vasp/scripts/propriedades/5.extract_gap_plot_with_slope_and_strain_v3.py > saida-slope_and_strain.out
```

`1.generate_EOS_inputs.sh`:
```bash
#!/bin/bash
set -euo pipefail

# -------------------------
# Configurações
# -------------------------
POSCAR="POSCAR"
ARQUIVOS="./arquivos"
#ESCALAS=($(seq 0.90 0.01 1.10)) #(0.94 0.96 0.98 1.00 1.02 1.04 1.06)

# ----------------------------------------------
# Escolha das escalas: para variação de ±1% no volume
# ----------------------------------------------
# O volume V varia com o cubo do fator de escala linear s: V = V0 * s^3
# Para garantir que a variação de volume fique dentro de ±1%:
#   → s = (1 ± 0.01)^(1/3) ≈ 0.9967 a 1.0033
# Assim, a escala linear varia apenas ±0.33% para obter uma variação de ±1% no volume
ESCALAS=(0.9967 0.9983 1.0000 1.0017 1.0033)


echo "Gerando POSCARs escalados com fatores: ${ESCALAS[*]}"
echo "Usando arquivos de entrada de: $ARQUIVOS"

# -------------------------
# Calcula V0 a partir dos vetores da matriz
# -------------------------

# Extrai os vetores da célula
read -r a1 a2 a3 <<< $(awk 'NR==3 {print $1, $2, $3}' $POSCAR)
read -r b1 b2 b3 <<< $(awk 'NR==4 {print $1, $2, $3}' $POSCAR)
read -r c1 c2 c3 <<< $(awk 'NR==5 {print $1, $2, $3}' $POSCAR)

# Calcula volume V0 = |a . (b x c)|
V0=$(echo "scale=8; \
  ($a1*($b2*$c3 - $b3*$c2) + \
   $a2*($b3*$c1 - $b1*$c3) + \
   $a3*($b1*$c2 - $b2*$c1))" | bc -l | awk '{print ($1<0)?-$1:$1}')

echo "Volume original V0 = $V0 Å³"

# -------------------------
# Loop sobre escalas
# -------------------------
for s in "${ESCALAS[@]}"; do
    # Calcula novo volume: V = V0 * s^3
    V=$(echo "scale=8; $V0 * ($s^3)" | bc -l)

    # Cria pasta
    DIR="Simulations/scale_${s}"
    mkdir -p "$DIR"

    # Gera POSCAR com fator de escala negativo
    cp "$POSCAR" "$DIR/POSCAR"
    sed -i "2s/.*/-$V/" "$DIR/POSCAR"

    # Copia arquivos de entrada
    cp "$ARQUIVOS/INCAR"    "$DIR/"
    cp "$ARQUIVOS/KPOINTS" "$DIR/"
    cp "$ARQUIVOS/POTCAR"  "$DIR/"

    echo "Gerado: $DIR/POSCAR (Volume: -$V Å³)"
done

echo "Todos os inputs foram gerados com sucesso."
```
`2.job-EOS_inputs.srm`:
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH -p par32
#SBATCH -J EOS_CsPbSnI3

echo "Nós alocados: $SLURM_JOB_NODELIST"
echo "CPUs: $SLURM_NTASKS"

cd $SLURM_SUBMIT_DIR
echo "Diretório de trabalho: $SLURM_SUBMIT_DIR"

module load vasp/6.3.0

VASP_CMD="mpirun -np $SLURM_NTASKS vasp_std"

SIM_DIR="Simulations"

scale_dirs=( $(ls -1d "$SIM_DIR"/scale_* | sort -V) )

for dir in "${scale_dirs[@]}"; do
  echo ">>> Rodando VASP em: $dir"
  (
    cd "$dir" || exit 1
    time $VASP_CMD > saida-EOS_inputs.out
    if [ $? -ne 0 ]; then
      echo "Erro no VASP em $dir"
    fi
  )
done

echo "Tudo concluído. Verifique os OUTCARs em cada pasta."
```

`3.fit_eos_vinet.py`:
```bash
#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os
import re
from lmfit import Model

# -------------------- CONFIGURAÇÕES --------------------
SIM_DIR = "Simulations"
OUT_PNG = "eos_fit_vinet_lmfit.png"
POSCAR = "POSCAR"
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
```

`4.run_relax_dos_bands.srm`
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH -p par32
#SBATCH -J deformation_potential

echo "Nós alocados: $SLURM_JOB_NODELIST"
cd $SLURM_SUBMIT_DIR
echo "Diretório: $SLURM_SUBMIT_DIR"

module load vasp/6.3.0

VASP_CMD="mpirun -np $SLURM_NTASKS vasp_std"
SIM_DIR="Simulations"
ARQ_DIR="arquivos"

for scale in ${SIM_DIR}/scale_*; do
  echo ">>> Processando $scale"

  # -----------------------------
  # 1) RELAXAÇÃO
  # -----------------------------
  mkdir -p $scale/rlx
  cp $scale/POSCAR $scale/rlx/
  cp $scale/POTCAR $scale/rlx/
  cp $ARQ_DIR/INCAR_rlx $scale/rlx/INCAR
  cp $ARQ_DIR/KPOINTS $scale/rlx/

  echo " • Relaxação..."
  cd $scale/rlx
    time $VASP_CMD > saida-rlx.out
  cd - >/dev/null

  # -----------------------------
  # 2) DOS
  # -----------------------------
  mkdir -p $scale/dos
  cp $scale/rlx/CONTCAR $scale/dos/POSCAR
  cp $scale/POTCAR $scale/dos/
  cp $ARQ_DIR/dos_INCAR $scale/dos/INCAR
  cp $ARQ_DIR/dos_KPOINTS $scale/dos/KPOINTS

  echo " • DOS..."
  cd $scale/dos
    time $VASP_CMD > saida-dos.out
  cd - >/dev/null

  # -----------------------------
  # 3) BANDAS
  # -----------------------------
  mkdir -p $scale/bands
  cp $scale/rlx/CONTCAR $scale/bands/POSCAR
  cp $scale/POTCAR $scale/bands/
  cp $ARQ_DIR/bands_INCAR $scale/bands/INCAR
  cp $ARQ_DIR/KPOINTS_band $scale/bands/KPOINTS
  cp $scale/dos/CHGCAR $scale/bands/

  echo " • Bandas..."
  cd $scale/bands
    time $VASP_CMD > saida-bands.out
  cd - >/dev/null

done

echo ">>> Todo o fluxo relaxação + DOS + bandas concluído!"
```



`5.extract_gap_plot_with_slope_and_strain_v3.py`:
```bash
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
```


Refazer para dft-1/2
