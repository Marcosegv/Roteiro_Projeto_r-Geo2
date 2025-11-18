#!/bin/bash
set -euo pipefail

# ---------- Config / parâmetros por variável de ambiente ----------
export LC_NUMERIC=C                           # garante ponto decimal
POSCAR="${POSCAR:-POSCAR}"                    # POSCAR de referência
ARQUIVOS="${ARQUIVOS:-./arquivos}"            # pasta com INCAR/KPOINTS/POTCAR
OUTDIR="${OUTDIR:-Simulations}"               # pasta de saída (ex.: Simulations_fit / Simulations)
SCALES_STR="${SCALES:-}"                      # lista de escalas (ex.: "0.90 0.91 ... 1.10")
SCALES_STR="${SCALES_STR//,/.}"               # troca vírgula por ponto se vier de pt_BR

# Defaults de escalas (varredura ampla) se nada for informado
if [[ -z "$SCALES_STR" ]]; then
  SCALES_STR="$(LC_NUMERIC=C seq 0.90 0.01 1.10)"
fi

# ---------- Sanity checks ----------
[[ -f "$POSCAR" ]] || { echo "ERRO: POSCAR '$POSCAR' não encontrado."; exit 1; }
for f in INCAR KPOINTS POTCAR; do
  [[ -f "$ARQUIVOS/$f" ]] || { echo "ERRO: Arquivo '$ARQUIVOS/$f' não encontrado."; exit 1; }
done

# ---------- Calcula V0 do POSCAR (suporta escala negativa e notação científica) ----------
V0=$(awk '
  NR==2 { scale=$1; next }
  NR==3 { a1=$1; a2=$2; a3=$3; next }
  NR==4 { b1=$1; b2=$2; b3=$3; next }
  NR==5 { c1=$1; c2=$2; c3=$3; next }
  END {
    det = a1*(b2*c3 - b3*c2) + a2*(b3*c1 - b1*c3) + a3*(b1*c2 - b2*c1);
    if (scale < 0) { vol = -scale; }         # escala negativa => 2a linha = |V|
    else           { vol = det * scale*scale*scale; }
    if (vol < 0) vol = -vol;
    printf("%.8f", vol);
  }
' "$POSCAR")

echo "V0 = $V0 Å^3"
echo "OUTDIR = $OUTDIR"
echo "ESCALAS = $SCALES_STR"

mkdir -p "$OUTDIR"

# ---------- Gera entradas para cada escala ----------
for s in $SCALES_STR; do
  # V = V0 * s^3  (feito no awk para aceitar floats em qualquer locale)
  V=$(awk -v V0="$V0" -v s="$s" 'BEGIN { printf("%.8f", V0 * s * s * s) }')

  DIR="$OUTDIR/scale_${s}"
  mkdir -p "$DIR"

  # POSCAR com volume absoluto (linha 2 negativa = V fixo)
  cp "$POSCAR" "$DIR/POSCAR"
  sed -i "2s/.*/-$V/" "$DIR/POSCAR"

  # Inputs eletrônicos
  cp "$ARQUIVOS/INCAR"   "$DIR/"
  cp "$ARQUIVOS/KPOINTS" "$DIR/"
  cp "$ARQUIVOS/POTCAR"  "$DIR/"

  echo "• $DIR  (V = -$V Å³)"
done

echo "[OK] Inputs gerados em: $OUTDIR"
