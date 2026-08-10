#!/bin/bash

start=$1
end=$2
outfile=$3

datafiles="/scratch/pawsey0272/zwan4817/J1752/addsearch/links/*.sf"

if [[ -e "$outfile" ]]; then
    echo "Error: $outfile already exists. Aborting." >&2
    exit 1
fi

for f in $datafiles; do
    fname="$(basename "$f")"
    [[ "$fname" =~ ^[0-9]{8}T[0-9]{2}:[0-9]{2}:[0-9]{2}_ ]] || continue
    ts="${fname%%_*}"
    if [[ "$ts" > "$start" && "$ts" < "$end" ]]; then
        echo "$f"
    fi
done | sort > "$outfile"

echo "Done. $(wc -l < "$outfile") files matched, written to $outfile"