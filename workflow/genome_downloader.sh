#!/bin/bash

# Ensure correct usage
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <TaxonID_or_ScientificName> <OutputDirectory>"
    exit 1
fi

QUERY="$1"       # Taxon ID or scientific name
OUTDIR="$2"      # Output directory

# Create the output directory if it doesn't exist
mkdir -p "$OUTDIR"

echo "Fetching genome for: $QUERY"
echo "Saving files in: $OUTDIR"

################################################################################
# Try multiple assembly levels in descending order of completeness
#   Complete Genome -> Chromosome -> Scaffold -> Contig
################################################################################

BEST_ASSEMBLY=""
BEST_LEVEL=""

for LEVEL in "Complete" "Chromosome" "Scaffold" "Contig"; do
  # Attempt to find the first assembly at this level
  CANDIDATE=$(datasets summary genome taxon "$QUERY" \
    | jq -r "[.reports[] | select(.assembly_info.assembly_level==\"$LEVEL\")][0].accession")

  if [ -n "$CANDIDATE" ] && [ "$CANDIDATE" != "null" ]; then
    BEST_ASSEMBLY="$CANDIDATE"
    BEST_LEVEL="$LEVEL"
    echo "Found a $LEVEL assembly: $BEST_ASSEMBLY"
    break
  fi
done

# If nothing found, exit with a message
if [ -z "$BEST_ASSEMBLY" ] || [ "$BEST_ASSEMBLY" == "null" ]; then
    echo "No valid genome assembly (Complete/Chromosome/Scaffold/Contig) found for: $QUERY"
    exit 1
fi

# Convert spaces in the level to underscores for a filename-safe string
LEVEL_SAFE="${BEST_LEVEL// /_}"

# Download the best available genome into the output directory
ZIPFILE="$OUTDIR/genome_${LEVEL_SAFE}_${QUERY}.zip"

echo "Downloading accession: $BEST_ASSEMBLY"
datasets download genome accession "$BEST_ASSEMBLY" \
--include genome \
--filename "$ZIPFILE"

# Extract the FASTA file and save it in the output directory
echo "Extracting FASTA file for $QUERY"
unzip -p -v "$ZIPFILE" ncbi_dataset/data/*/*.fna > "$OUTDIR/${QUERY}.fa"

# Remove the downloaded ZIP file
rm -v "$ZIPFILE"

echo "Download complete: $OUTDIR/${QUERY}.fa (Assembly: $BEST_ASSEMBLY, Level: $BEST_LEVEL)"
