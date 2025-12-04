#!/bin/bash
#SBATCH --job-name=vasp      # Escolha um nome para o Job!
#SBATCH --partition=local    # Nome da fila
#SBATCH --nodes=1            # nó einstein
#SBATCH --ntasks=8           # 8 processos MPI; se quiser rodar com menos cores altere aqui!
#SBATCH --time=7-00:00:00    # Tempo limite do cálculo rodando, opcional!

echo "Nós alocados: $SLURM_JOB_NODELIST"
echo "CPUs: $SLURM_NTASKS"

cd "$SLURM_SUBMIT_DIR"
echo "Diretório: $SLURM_SUBMIT_DIR"

set -euo pipefail
export LC_ALL=C

# ========== CONFIG ==========
A_START=4.1
A_END=5.0
A_STEP=0.1

VASP_CMD="mpirun -np $SLURM_NTASKS vasp_std"
BASE_DIR="BASE_FILES"
C_SCAN_SCRIPT="./c_scan.sh"         # seu scanner de c (ISIF=2)
RESULTS_CSV="results_all.csv"
# ============================

# 🔧 torne o caminho absoluto (funciona mesmo após cd)
if command -v realpath >/dev/null 2>&1; then
  C_SCAN_SCRIPT="$(realpath "$C_SCAN_SCRIPT")"
else
  C_SCAN_SCRIPT="$(cd "$(dirname "$C_SCAN_SCRIPT")" && pwd)/$(basename "$C_SCAN_SCRIPT")"
fi
[[ -f "$C_SCAN_SCRIPT" ]] || { echo "[erro] não encontrei $C_SCAN_SCRIPT"; exit 1; }

# --------- Sanidade ----------
for f in POSCAR POTCAR INCAR; do
  [[ -f "$BASE_DIR/$f" ]] || { echo "[erro] falta $BASE_DIR/$f em $BASE_DIR"; exit 1; }
done
[[ -x "$C_SCAN_SCRIPT" || -f "$C_SCAN_SCRIPT" ]] || { echo "[erro] não encontrei $C_SCAN_SCRIPT"; exit 1; }
command -v awk >/dev/null || { echo "[erro] awk não encontrado"; exit 1; }
# -----------------------------

# Vetores |a0|,|b0|,|c0| do POSCAR de referência (para strain)
read -r a1 a2 a3 < <(awk 'NR==3{print $1,$2,$3}' "$BASE_DIR/POSCAR")
read -r b1 b2 b3 < <(awk 'NR==4{print $1,$2,$3}' "$BASE_DIR/POSCAR")
read -r c1 c2 c3 < <(awk 'NR==5{print $1,$2,$3}' "$BASE_DIR/POSCAR")
a0=$(awk -v x=$a1 -v y=$a2 -v z=$a3 'BEGIN{printf "%.8f", sqrt(x*x+y*y+z*z)}')
b0=$(awk -v x=$b1 -v y=$b2 -v z=$b3 'BEGIN{printf "%.8f", sqrt(x*x+y*y+z*z)}')
c0=$(awk -v x=$c1 -v y=$c2 -v z=$c3 'BEGIN{printf "%.8f", sqrt(x*x+y*y+z*z)}')

echo "a_A,c_star_A,E_eV,eps_a,eps_c,sxx_GPa,syy_GPa,szz_GPa,workdir" > "$RESULTS_CSV"

# lista de 'a'
mapfile -t A_LIST < <(awk -v a0="$A_START" -v a1="$A_END" -v da="$A_STEP" 'BEGIN{for(x=a0;x<=a1+1e-9;x+=da) printf "%.3f\n", x}')
echo "[info] a-list: ${A_LIST[*]}"

# --- helpers ---
make_poscar_with_a() { # reescala |a|=|b|=a_target mantendo direções
  local poscar_in="$1" poscar_out="$2" a_target="$3"
  awk -v at="$a_target" '
    NR==1{t=$0;next} NR==2{s=$0;next}
    NR==3{a1=$1;a2=$2;a3=$3;next}
    NR==4{b1=$1;b2=$2;b3=$3;next}
    NR==5{c1=$1;c2=$2;c3=$3;next}
    {rest=rest $0 "\n";next}
    END{
      an=sqrt(a1*a1+a2*a2+a3*a3); af=at/an
      bn=sqrt(b1*b1+b2*b2+b3*b3); bf=at/bn
      printf "%s\n%s\n",t,s
      printf "  %.16f  %.16f  %.16f\n", a1*af,a2*af,a3*af
      printf "  %.16f  %.16f  %.16f\n", b1*bf,b2*bf,b3*bf
      printf "  %.16f  %.16f  %.16f\n", c1,c2,c3
      printf "%s",rest
    }' "$poscar_in" > "$poscar_out"
}
c_length_from_poscar(){ awk 'NR==5{printf "%.8f\n", sqrt($1*$1+$2*$2+$3*$3)}' "$1"; }
min_from_summary(){ # devolve "scale,c_A,E_eV,dir"
  local summary="$1" scandir="$2"
  awk -v sd="$scandir" -F, 'NR==1{next}{E=$3+0;if(NR==2||E<Emin){Emin=E;S=$1;C=$2}} END{printf "%.3f,%s,%.10f,%s/c_%.3f\n",S,C,Emin,sd,S}' "$summary"
}
stress_from_outcar(){
  awk '
    /in kB/ {
      getline; xx=$1;
      getline; yy=$2;
      getline; zz=$3;
    }
    END{
      if(xx!="") printf "%.6f,%.6f,%.6f\n", xx*0.1, yy*0.1, zz*0.1;
      else print "NaN,NaN,NaN"
    }
  ' "$1"
}
# ----------------

for a in "${A_LIST[@]}"; do
  wdir=$(printf "a_%.2f" "$a")
  echo "===================="
  echo "[a] alvo = $a Å  →  $wdir"
  mkdir -p "$wdir/BASE_FILES"

  # POSCAR de entrada com a,b reescalados
  make_poscar_with_a "$BASE_DIR/POSCAR" "$wdir/BASE_FILES/POSCAR" "$a"

  # insumos do scan-c + propriedades
  cp "$BASE_DIR/INCAR"                 "$wdir/BASE_FILES/"
  [[ -f "$BASE_DIR/KPOINTS" ]] && cp "$BASE_DIR/KPOINTS" "$wdir/BASE_FILES/" || true
  cp "$BASE_DIR/POTCAR"                "$wdir/BASE_FILES/"
  cp "$BASE_DIR/2.INCAR-bands_SCF"     "$wdir/BASE_FILES/"
  cp "$BASE_DIR/2.KPOINTS-bands_NSCF"  "$wdir/BASE_FILES/"
  cp "$BASE_DIR/3.INCAR-bands_NSCF"    "$wdir/BASE_FILES/"
  cp "$BASE_DIR/3.KPOINTS-dos"         "$wdir/BASE_FILES/"
  cp "$BASE_DIR/4.INCAR-dos"           "$wdir/BASE_FILES/"

  # roda scanner de c
  ( cd "$wdir" && bash "$C_SCAN_SCRIPT" )

  # mínimo de energia
  summary="$wdir/SCAN_c/summary.csv"
  [[ -s "$summary" ]] || { echo "[erro] sem summary em $wdir"; exit 1; }
  read -r sc cstar Eev bestdir < <(min_from_summary "$summary" "$wdir/SCAN_c" | awk -F, '{print $1,$2,$3,$4}')
  echo "[ok] min em $bestdir | scale=$sc  c*=$cstar Å  E=$Eev eV"

  # stress no ponto ótimo
  outcar="$bestdir/OUTCAR"
  read -r sxx syy szz < <(stress_from_outcar "$outcar" | awk -F, '{print $1,$2,$3}')

  # strains
  eps_a=$(awk -v a="$a" -v a0="$a0" 'BEGIN{printf "%.8f", (a-a0)/a0}')
  cstar_true=$(c_length_from_poscar "$bestdir/POSCAR")
  eps_c=$(awk -v c="$cstar_true" -v c0="$c0" 'BEGIN{printf "%.8f", (c-c0)/c0}')

  # registra
  echo "$a,$cstar_true,$Eev,$eps_a,$eps_c,$sxx,$syy,$szz,$wdir" >> "$RESULTS_CSV"

  # ===== Propriedades para este (a,c*) =====
  prop="$wdir/properties"; mkdir -p "$prop/bands" "$prop/dos"

  # usa CONTCAR relaxado do melhor c
  cp "$bestdir/CONTCAR"           "$prop/bands/POSCAR"
  cp "$wdir/BASE_FILES/POTCAR"    "$prop/bands/POTCAR"

  # ---------- SCF (usa KPOINTS de SCF) ----------
  cp "$wdir/BASE_FILES/2.INCAR-bands_SCF" "$prop/bands/INCAR"
  cp "$wdir/BASE_FILES/KPOINTS"           "$prop/bands/KPOINTS"
  echo "[run] SCF → $prop/bands"
  ( cd "$prop/bands" && $VASP_CMD > vasp_scf.out )

  # guarda inputs usados
  mv "$prop/bands/INCAR"   "$prop/bands/INCAR-SCF"
  mv "$prop/bands/KPOINTS" "$prop/bands/KPOINTS-SCF"

  # ---------- Bands (NSCF) ----------
  cp "$wdir/BASE_FILES/3.INCAR-bands_NSCF"   "$prop/bands/INCAR"
  cp "$wdir/BASE_FILES/2.KPOINTS-bands_NSCF" "$prop/bands/KPOINTS"
  echo "[run] Bands (NSCF) → $prop/bands"
  ( cd "$prop/bands" && $VASP_CMD > vasp_bands_nscf.out )

  # ---------- DOS ----------
  cp "$prop/bands/POSCAR"  "$prop/dos/POSCAR"
  cp "$prop/bands/POTCAR"  "$prop/dos/POTCAR"
  cp "$prop/bands/CHGCAR"  "$prop/dos/CHGCAR"
  cp "$prop/bands/WAVECAR" "$prop/dos/WAVECAR"
  cp "$wdir/BASE_FILES/4.INCAR-dos"  "$prop/dos/INCAR"
  cp "$wdir/BASE_FILES/3.KPOINTS-dos" "$prop/dos/KPOINTS"
  echo "[run] DOS → $prop/dos"
  ( cd "$prop/dos" && $VASP_CMD > vasp_dos.out )

done

echo "[fim] Tabela mestre: $RESULTS_CSV"
