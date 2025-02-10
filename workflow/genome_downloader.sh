#!/bin/bash

# Ensure correct usage
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <AccessionCode_or_TaxonID_or_ScientificName> <OutputDirectory>"
    exit 1
fi

QUERY="$1"       # Accession Code, Taxon ID, or Scientific Name
OUTDIR="$(realpath -m "$2")"  # Normalize and clean output directory

echo "Fetching genome for: $QUERY"

################################################################################
# Check if input is an Accession Code (GCF_ or GCA_)
################################################################################
if [[ "$QUERY" =~ ^G[CA]F_[0-9]+\.[0-9]+$ ]]; then
    echo "Detected Accession Code: $QUERY"
    BEST_ASSEMBLY="$QUERY"
    BEST_LEVEL="Direct_Accession"
    
    # Replace "." with "_" in filenames
    QUERY_SAFE="${QUERY//./_}"
else
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

    # If not an accession code, keep QUERY as is
    QUERY_SAFE="$QUERY"
fi

# Convert spaces in the level to underscores for a filename-safe string
LEVEL_SAFE="${BEST_LEVEL// /_}"

# Normalize all filenames using `realpath -m`
ZIPFILE="$(realpath -m "$OUTDIR/genome_${LEVEL_SAFE}_${QUERY_SAFE}.zip")"
FASTA_FILE="$(realpath -m "$OUTDIR/${QUERY_SAFE}.fa")"

echo "Downloading accession: $BEST_ASSEMBLY"
datasets download genome accession "$BEST_ASSEMBLY" \
--include genome \
--filename "$ZIPFILE"

# Get the actual uncompressed size of the FASTA file (sum sizes of extracted files)
UNCOMPRESSED_SIZE=$(unzip -Z -1 "$ZIPFILE" ncbi_dataset/data/*/*.fna | xargs -I{} unzip -p "$ZIPFILE" {} | wc -c)

# Extract the FASTA file with an accurate progress bar
echo "Extracting FASTA file for $QUERY..."
unzip -p "$ZIPFILE" ncbi_dataset/data/*/*.fna | pv -s "$UNCOMPRESSED_SIZE" > "$FASTA_FILE"

# Remove the ZIP file with a progress bar
echo "Removing ZIP file: $ZIPFILE"
echo "$ZIPFILE" | pv -l -s 1 | xargs -d '\n' rm

echo "Download complete: $FASTA_FILE (Assembly: $BEST_ASSEMBLY, Level: $BEST_LEVEL)"
