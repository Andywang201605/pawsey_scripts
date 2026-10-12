#!/bin/bash -l
#SBATCH --account=pawsey0272
#SBATCH --output=./logs/%x_%A_%a.out
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --job-name=download
#SBATCH --array=1-10

FULL_LIST=full_list.txt
N_TASKS=10

wgetcmd="wget --tries=0 --retry-connrefused --waitretry=5 --timeout=30 --read-timeout=60 --dns-timeout=10  --retry-on-host-error  --retry-on-http-error=401,403,408,429,500,502,503,504  --connect-timeout=10 -c"

N_LINES=$(wc -l < "$FULL_LIST")
CHUNK_SIZE=$(( (N_LINES + N_TASKS - 1) / N_TASKS ))

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))

sed -n "${START},${END}p" "$FULL_LIST" | \
while IFS= read -r line; do
    echo "Processing: $line"
    $wgetcmd "$line"
done