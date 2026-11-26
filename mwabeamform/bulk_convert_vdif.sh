
#!/bin/bash

# Usage:
# ./bulk_convert_vdif.sh <hdr_list> [-b nbit] [-n nchan] [-d dm] [-o output_base]
#
# Examples:
# ./bulk_convert_vdif.sh ./hdr_list.txt
# ./bulk_convert_vdif.sh ./hdr_list.txt -d 181
# ./bulk_convert_vdif.sh ./hdr_list.txt -b 8 -n 64 -d 181 -o ./filterbank
#
# Output example:
# ./1468789376/1468789376_ch145_beam02.fil

set -e

SCRIPT_NAME=$(basename "$0")
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# Path to the existing single-file conversion script
CONVERT_SCRIPT="${SCRIPT_DIR}/convert_vdif.sh"

# Default parameters
NBIT=8
NCHAN=64
DM=""
OUTBASE="./"

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
            OUTBASE="$OPTARG"
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

# Check conversion script
if [[ ! -f "$CONVERT_SCRIPT" ]]; then
    echo "Error: Conversion script not found: $CONVERT_SCRIPT"
    exit 1
fi

# Create output base directory
mkdir -p "$OUTBASE"

echo "==================================="
echo "Bulk VDIF to filterbank conversion"
echo "Input list : $INPUT"
echo "Output base: $OUTBASE"
echo "NBIT       : $NBIT"
echo "NCHAN      : $NCHAN"
echo "DM         : ${DM:-not specified}"
echo "==================================="

# Process each HDR file
while IFS= read -r hdr || [[ -n "$hdr" ]]; do

    # Remove carriage return (for Windows-formatted lists)
    hdr="${hdr%$'\r'}"

    # Skip empty lines and comments
    [[ -z "$hdr" || "$hdr" == \#* ]] && continue

    # Extract observation ID from parent directory
    OBS_ID=$(basename "$(dirname "$hdr")")

    # Create observation-specific output directory
    OUTDIR="${OUTBASE}/${OBS_ID}"
    mkdir -p "$OUTDIR"

    echo ""
    echo "-----------------------------------"
    echo "Observation: $OBS_ID"
    echo "Input      : $hdr"
    echo "Output dir : $OUTDIR"
    echo "-----------------------------------"

    # Build conversion command
    CMD=(bash "$CONVERT_SCRIPT" "$hdr" -o "$OUTDIR" -b "$NBIT" -n "$NCHAN")

    if [[ -n "$DM" ]]; then
        CMD+=(-d "$DM")
    fi

    echo "Executing command: ${CMD[@]}"

    # Execute conversion
    "${CMD[@]}"

done < "$INPUT"

echo ""
echo "==================================="
echo "All conversions completed!"
echo "==================================="
