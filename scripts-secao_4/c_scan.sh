#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# -------- Config --------
BASE_DIR="BASE_FILES"          # INCAR, POTCAR, POSCAR (referência) e opcional KPOINTS
VASP_CMD="${VASP_CMD:-mpirun -np ${SLURM_NTASKS:-4} vasp_std}"
OUT_ROOT="SCAN_c"
# Janela de ±1.0% em passos de 0.2%  → 0.990 .. 1.010 passo 0.002
SCAN_MIN=0.900
SCAN_MAX=1.100
SCAN_STEP=0.002
# ------------------------

# Sanidade
command -v awk >/dev/null || { echo "[erro] awk não encontrado"; exit 1; }
for f in INCAR POTCAR POSCAR; do
  [[ -f "$BASE_DIR/$f" ]] || { echo "[erro] falta $BASE_DIR/$f"; exit 1; }
done
mkdir -p "$OUT_ROOT"

# Aviso se não houver KPOINTS e INCAR não tiver KSPACING
if [[ ! -f "$BASE_DIR/KPOINTS" ]]; then
  if ! grep -qi '^[[:space:]]*KSPACING[[:space:]]*=' "$BASE_DIR/INCAR"; then
    echo "[aviso] Não encontrei $BASE_DIR/KPOINTS e o INCAR não tem KSPACING definido."
    echo "        Defina KSPACING no INCAR OU forneça um arquivo KPOINTS."
  fi
fi

# Lista de fatores (geração robusta com awk)
mapfile -t FACTORS < <(
  awk -v min="$SCAN_MIN" -v max="$SCAN_MAX" -v step="$SCAN_STEP" '
    BEGIN{
      for (x=min; x<=max+1e-9; x+=step) printf "%.3f\n", x;
    }'
)
echo "[info] Fatores: ${FACTORS[*]}"

# Funções utilitárias
c_length_from_poscar() {
  awk 'NR==5{printf "%.8f\n", sqrt($1*$1+$2*$2+$3*$3)}' "$1"
}
sigma_zz_from_outcar() {
  # último bloco de tensões "in kB" → pega s_zz e converte kB→GPa (×0.1)
  awk '
    /in kB/ {inblk=1; next}
    inblk && NF==6 {szz=$3; inblk=0}
    END{if(szz!="") printf "%.6f\n", szz*0.1; else print "NaN"}
  ' "$1"
}
energy_from_oszicar() {
  awk '/F=/{E=$5} END{if(E!="") print E; else print "NaN"}' "$1"
}

# Loop principal
for s in "${FACTORS[@]}"; do
  D="$OUT_ROOT/c_${s}"
  mkdir -p "$D"

  # Pula se já convergiu
  if [[ -f "$D/OUTCAR" ]] && grep -q "Voluntary" "$D/OUTCAR"; then
    echo "[skip] $D já finalizado."
    continue
  fi

  # Copia base
  cp "$BASE_DIR/INCAR"  "$D/"
  [[ -f "$BASE_DIR/KPOINTS" ]] && cp "$BASE_DIR/KPOINTS" "$D/" || true
  cp "$BASE_DIR/POTCAR" "$D/"

  # Gera POSCAR escalando APENAS o vetor c (linha 5)
  {
    read -r line1
    read -r line2
    read -r a1 a2 a3
    read -r b1 b2 b3
    read -r c1 c2 c3
    echo "$line1"
    echo "$line2"
    printf "  %.16f  %.16f  %.16f\n" "$a1" "$a2" "$a3"
    printf "  %.16f  %.16f  %.16f\n" "$b1" "$b2" "$b3"
    # escala c
    awk -v c1="$c1" -v c2="$c2" -v c3="$c3" -v s="$s" 'BEGIN{
      printf "  %.16f  %.16f  %.16f\n", c1*s, c2*s, c3*s
    }'
    cat
  } < "$BASE_DIR/POSCAR" > "$D/POSCAR"

  echo "[run] $D"
  ( cd "$D" && $VASP_CMD > vasp.out )
done

# Pós-processamento
OUT_CSV="$OUT_ROOT/summary.csv"
echo "scale,c_A,E_eV,sigmaZZ_GPa" > "$OUT_CSV"
for s in "${FACTORS[@]}"; do
  D="$OUT_ROOT/c_${s}"
  POS="$D/POSCAR"
  OSZ="$D/OSZICAR"
  OUT="$D/OUTCAR"
  cA="$(c_length_from_poscar "$POS" 2>/dev/null || echo NaN)"
  EeV="$(energy_from_oszicar "$OSZ" 2>/dev/null || echo NaN)"
  szz="$(sigma_zz_from_outcar "$OUT" 2>/dev/null || echo NaN)"
  echo "$s,$cA,$EeV,$szz" >> "$OUT_CSV"
done

echo "[ok] Resultados em: $OUT_CSV"
