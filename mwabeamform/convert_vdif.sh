
#!/bin/bash

# Usage:
# ./convert_vdif.sh <hdr_file_or_list> [-o output_dir] [-b nbit] [-n nchan] [-d dm]
#
# Examples:
# ./convert_vdif.sh observation.hdr
# ./convert_vdif.sh observation.hdr -o ./filterbank
# ./convert_vdif.sh hdr_list.txt -o ./filterbank -b 8 -n 64
# ./convert_vdif.sh hdr_list.txt -o ./filterbank -b 8 -n 64 -d 120

set -e

# Default parameters
OUTDIR="./"
NBIT=8
NCHAN=64
DM=""

# Check positional argument
if [[ $# -lt 1 ]]; then
    echo "Usage: ${SCRIPT_NAME} <hdr_list> [-b nbit] [-n nchan] [-d dm] [-o output_base]"
    exit 1
fi

# Get input list
INPUT="$1"
shift

# Parse optional arguments
while getopts ":b:n:d:o:" opt; do
    case "$opt" in
        b)
            NBIT="$OPTARG"
            ;;
        n)
            NCHAN="$OPTARG"
            ;;
        d)
            DM="$OPTARG"
            ;;
        o)
            OUTDIR="$OPTARG"
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument."
            exit 1
            ;;
        \?)
            echo "Error: Invalid option -$OPTARG"
            exit 1
            ;;
    esac
done

# Check input list
if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input list not found: $INPUT"
    exit 1
fi

SCRIPT_NAME=$(basename "$0")

# Create output directory
mkdir -p "$OUTDIR"

# Convert a single HDR file
convert_file() {
    local hdr="$1"

    [[ -z "$hdr" || "$hdr" == \#* ]] && return

    if [[ ! -f "$hdr" ]]; then
        echo "Warning: File not found: $hdr"
        return
    fi

    local filename
    filename=$(basename "$hdr" .hdr)
    local outfile="${OUTDIR}/${filename}.fil"

    echo "Converting: $hdr"
    echo "Output: $outfile"

    # if [[ -z "$DM" ]]; then
    #     echo "No DM provided."
    #     echo "Running: $psraly digifil -b $NBIT -F$NCHAN -o $outfile $hdr"
    # else
    #     echo "DM provided: $DM"
    #     echo "Running: $psraly digifil -b $NBIT -F${NCHAN}:D -D $DM -o $outfile $hdr"
    # fi


    if [[ -z "$DM" ]]; then
        $psraly digifil -b "$NBIT" \
                -F"$NCHAN" \
                -o "$outfile" \
                -d 1 \
                "$hdr"
    else
        $psraly digifil -b "$NBIT" \
                -F"${NCHAN}:D" \
                -D "$DM" \
                -o "$outfile" \
                -d 1 \
                "$hdr"
    fi

    echo "Finished: $outfile"
    echo "-----------------------------------"
}

# Process input file or list
if [[ "$INPUT" == *.hdr ]]; then
    convert_file "$INPUT"
elif [[ -f "$INPUT" ]]; then
    while IFS= read -r hdr || [[ -n "$hdr" ]]; do
        convert_file "$hdr"
    done < "$INPUT"
else
    echo "Error: Input file or list not found: $INPUT"
    exit 1
fi

echo "All conversions completed!"
