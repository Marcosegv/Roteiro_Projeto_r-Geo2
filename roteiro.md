# Roteiro Projeto r-GeO2
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

## 3. Variação de volume e potencial de deformação (strain)

## 3.1 Objetivo físico

Buscamos duas coisas:

1. **Ajustar a Equação de Estado (EoS)** para obter o **volume de equilíbrio** $V_0$, o **módulo volumétrico** $B_0$ e sua derivada $B'$ a partir da curva de energia total $E(V)$.
2. **Medir o potencial de deformação do gap**, isto é, como a **largura do gap** $E_g$ varia com a **deformação volumétrica** $\varepsilon$ nas vizinhanças do equilíbrio.

Definimos a deformação volumétrica em torno do volume de equilíbrio como

$$
\varepsilon \equiv \frac{V - V_0}{V_0},
$$

e o **potencial de deformação do gap** como a derivada de $E_g$ em relação à deformação no ponto de equilíbrio:

$$
\Xi \equiv \left.\frac{dE_g}{d\varepsilon}\right|_{\varepsilon=0}.
$$

Na prática, procedemos em duas etapas complementares:

- **(i) Fit amplo de $E(V)$** usando uma EoS (neste projeto, **Vinet**) com cálculos **estáticos a volume fixo** ($\text{NSW}=0$, $\text{ISIF}=2$). Dessa etapa obtemos $V_0$, $B_0$ e $B'$.
- **(ii) Varredura fina ao redor de $V_0$** (tipicamente $\pm 1\%$ em volume) para extrair $E_g(\varepsilon)$. Fazemos **relax iônico a volume fixo** ($\text{ISIF}=2$, $\text{NSW}>0$) e depois **SCF + NSCF de bandas** para cada volume. Ajustamos uma reta a $E_g(\varepsilon)$ e o coeficiente angular fornece $\Xi$.

> Observação: usar o **$V_0$ do fit** como referência de $\varepsilon$ garante consistência entre os pontos (mesmos parâmetros eletrônicos) e evita misturar efeitos de relaxações geométricas na definição de deformação.

### 3.1 Preparando os arquivos necessários

Voltem para a pasta base do seu pseudopotencial e crie uma pasta chamada `Volume`:

```bash
mkdir -p volume ; mkdir -p volume/dft

cd volume/dft

mkdir -p arquivos
```

Copie os arquivos necessários para a seguinte pasta:
```bash
cp ../../2.optimize_completa/INCAR arquivos/INCAR_rlx
cp ../../2.optimize_completa/POSCAR arquivos/.
cp ../../2.optimize_completa/POTCAR arquivos/.
cp ../../2.optimize_completa/KPOINTS arquivos/.
cp ~/FF281-2025/Geo2-Projeto_Final/0.INPUTS_BASE/INCAR-bands* arquivos/.
cp ~/FF281-2025/Geo2-Projeto_Final/0.INPUTS_BASE/KPOINTS_bands arquivos/.
```

** Atenção: No arquivo `INCAR_rlx`, mude de `ISIF = 3` para `ISIF = 2`. 

Crie o arquivo `INCAR` no diretório `arquivos` (ou seja, `arquivos/INCAR`):
```bash
SYSTEM = r-GeO2

Electronic relaxation:
ENCUT   = 520.000   ! Cutoff energy
ALGO    = Normal
NELMIN  = 6
NELM    = 80
NELMDL  = -4
EDIFF   = 1.0E-5

Calculation mode:
PREC    = Normal
ADDGRID = .TRUE.
LASPH   = .TRUE.

# Static calculation (no relaxation):
NSW     = 0         ! No ionic steps
EDIFFG  = -0.010    ! (não será usado pois NSW=0)
IBRION  = -1        ! No ionic relaxation
ISIF    = 2         ! Stresses calculated, but no cell relaxation
POTIM   = 0.00

# Integration over the Brillouin zone:
ISMEAR  = 0
SIGMA   = 0.010

OUTCAR size:
NWRITE  = 1
LWAVE   = .FALSE.
LCHARG  = .FALSE.

Key for parallel mode calculation:
NCORE = 4
LPLANE = .TRUE.
```

**Atenção: É importante copiar todos os scripts que estão em anexo (pasta: `scripts-Secao_3`) para o diretório atual. Pode-se também deixar em outra pasta e digitar todo o caminho, no entanto, caso não tenha experiencia com script bash, recomendo não fazer.**

### 3.2. Fase de FIT (varredura ampla)
Gerar inputs para o ajuste (em `Simulations_fit`)
```bash
POSCAR=arquivos/POSCAR ARQUIVOS=arquivos \
OUTDIR=Simulations_fit SCALES="$(seq 0.90 0.01 1.10)" \
bash 1.generate_EOS_inputs.sh
```

Rodar EOS estático nos volumes do fit:
```bash
sbatch --export=ALL,SIM_DIR=Simulations_fit 2.job-EOS_inputs.srm
```

Ajustar Vinet e plotar:
```bash
python 3.fit_eos_vinet.py > fit-EOS-Vinet.out
```


### 3.3. Fase fina (±1% ao redor de V0) → gaps/strain
Gerar inputs para ±1% (em Simulations)
```bash
POSCAR=arquivos/POSCAR ARQUIVOS=arquivos \
OUTDIR=Simulations SCALES="0.9967 0.9983 1.0000 1.0017 1.0033" \
bash 1.generate_EOS_inputs.sh
```

Relax interno + SCF + NSCF (bandas) nos volumes finos!
```bash
sbatch --export=ALL,SIM_DIR=Simulations 4.run_relax_scf_nscf_bands.srm
```

Extrair gaps vs V e vs strain:
```bash
python 5.extract_gap_plot_with_slope_and_strain_v3.py > saida-slope_and_strain.out
```

Agora, volte para a pasta volume. Crie um diretório para `dft-1_2` e repita toda a seção 3.3. **Lembre-se de pegar o POTCAR do dft-1/2!!!**

## 4. Varredura anisotrópica \(a,c\) e propriedades eletrônicas

Nesta etapa, vamos usar o resultado “melhor” da Seção 3 (volume de referência) para:
1. Fazer uma **varredura em \(a\)** (parâmetro in-plane), e para cada \(a\),  
2. Fazer uma **varredura em \(c\)** (via `c_scan.sh`) para achar o \(c^\*\) que minimiza a energia,  
3. Calcular **bandas + DOS** no ponto \((a, c^\*)\).  

Tudo isso já está automatizado nos scripts `1.run.sh` e `c_scan.sh`.

### 4.1 Preparar a pasta e os arquivos de base

Na pasta do funcional (por exemplo, `LDA/`):

```bash
cd /home/ff281/FF281-2025/Geo2-Projeto_Final/LDA

mkdir -p secao_4
cd secao_4

# Copiar os scripts prontos
cp ~/FF281-2025/Fisica_Materia_Condensada-FF281/scripts-secao_4/1.run.sh .
cp ~/FF281-2025/Fisica_Materia_Condensada-FF281/scripts-secao_4/c_scan.sh .

chmod +x 1.run.sh c_scan.sh

# Pasta de arquivos de entrada de referência
mkdir -p BASE_FILES
```

Agora é fundamental montar corretamente a pasta `BASE_FILES/` usando os melhores resultados da seção anterior (Seção 3 – variação de volume / V₀, bandas e DOS):

```bash
BASE_FILES/
 ├── POSCAR
 ├── POTCAR
 ├── INCAR           
 ├── KPOINTS      
 ├── 2.INCAR-bands_SCF
 ├── 2.KPOINTS-bands_NSCF
 ├── 3.INCAR-bands_NSCF
 ├── 3.KPOINTS-dos
 └── 4.INCAR-dos
```

## 4.2 Rodar a varredura a,c e as propriedades
Rode o comando com

```bash
sbatch 1.run.sh
```
