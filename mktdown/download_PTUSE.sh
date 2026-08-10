#!/bin/bash -l
#SBATCH --account=pawsey0272
#SBATCH --partition=copy
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --mem=10G

wgetcmd="wget --tries=0 --retry-connrefused --waitretry=5 --timeout=30 --read-timeout=60 --dns-timeout=10  --retry-on-host-error  --retry-on-http-error=401,403,408,429,500,502,503,504  --connect-timeout=10 -c"

while IFS= read -r line; do
    echo "Processing: $line"
    $wgetcmd "$line"
done < full_list.txt

