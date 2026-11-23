#!/usr/bin/env bash
set -euo pipefail

: "${TEMPERAMENT_ANALYSIS_ROOT:?Set TEMPERAMENT_ANALYSIS_ROOT}"
cd "$TEMPERAMENT_ANALYSIS_ROOT"

OUT="results/ldsc/developmental_rg"
REF="data/ldsc_reference/eur_w_ld_chr/"

mkdir -p "$OUT"

temperament=(EMO ACT SHY)
outcomes=(CINT CADHD AINT ASUD ARISK)

for trait in "${temperament[@]}"; do
    for outcome in "${outcomes[@]}"; do

        if [[ "$outcome" == C* ]]; then
            group="childhood"
        else
            group="adult"
        fi

        prefix="$OUT/${trait}_vs_${outcome}"

        if [[ -f "${prefix}.log" ]] &&
           grep -q "Genetic Correlation:" "${prefix}.log"; then
            echo "SKIPPING $trait vs $outcome"
            continue
        fi

        echo "RUNNING $trait vs $outcome"

        if python tools/ldsc/ldsc.py \
            --rg "data/ldsc_munged/temperament/${trait}.sumstats.gz,data/ldsc_munged/${group}/${outcome}.sumstats.gz" \
            --ref-ld-chr "$REF" \
            --w-ld-chr "$REF" \
            --out "$prefix" \
            > "${prefix}.stdout" 2>&1; then

            grep "Genetic Correlation:" "${prefix}.log" ||
                echo "NO ESTIMATE: $trait vs $outcome"
        else
            echo "FAILED: $trait vs $outcome"
            tail -15 "${prefix}.stdout"
        fi
    done
done
