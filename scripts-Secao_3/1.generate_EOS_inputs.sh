#!/bin/bash
set -euo pipefail

# Parâmetros por variável de ambiente (com padrões)
POSCAR="${POSCAR:-POSCAR}"
ARQUIVOS="${ARQUIVOS:-./arquivos}"
OUTDIR="${OUTDIR:-Simulations}"          # <- para escolher a pasta de saída
SCALES_STR="${SCALES:-}"                  # <- lista de escalas; se vazio usa seq 0.90..1.10

# Default de escalas (varredura ampla)
if [[ -z "$SCALES_STR" ]]; then
  SCALES_STR="$(seq 0.90 0.01 1.10)"
fi

# Calcula V0 do POSCAR
read -r a1 a2 a3 <<< $(awk 'NR==3 {print $1, $2, $3}' "$POSCAR")
read -r b1 b2 b3 <<< $(awk 'NR==4 {print $1, $2, $3}' "$POSCAR")
read -r c1 c2 c3 <<< $(awk 'NR==5 {print $1, $2, $3}' "$POSCAR")
V0=$(echo "scale=8; ($a1*($b2*$c3-$b3*$c2) + $a2*($b3*$c1-$b1*$c3) + $a3*($b1*$c2-$b2*$c1))" | bc -l | awk '{print ($1<0)?-$1:$1}')

echo "V0 = $V0 Å^3"
echo "OUTDIR = $OUTDIR"
echo "ESCALAS = $SCALES_STR"

mkdir -p "$OUTDIR"

# Gera entradas para cada escala
for s in $SCALES_STR; do
  V=$(echo "scale=8; $V0 * ($s^3)" | bc -l)
  DIR="$OUTDIR/scale_${s}"
  mkdir -p "$DIR"

  cp "$POSCAR" "$DIR/POSCAR"
  sed -i "2s/.*/-$V/" "$DIR/POSCAR"     # escala negativa -> volume absoluto

  cp "$ARQUIVOS/INCAR"    "$DIR/"
  cp "$ARQUIVOS/KPOINTS"  "$DIR/"
  cp "$ARQUIVOS/POTCAR"   "$DIR/"

  echo "• $DIR  (V = -$V Å³)"
done

echo "[OK] Inputs gerados em: $OUTDIR"

