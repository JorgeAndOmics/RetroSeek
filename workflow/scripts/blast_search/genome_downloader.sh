#!/bin/bash
# =============================================================================
# genome_downloader.sh
# =============================================================================
# Download one genome from NCBI Datasets: an assembly accession directly, or the
# most complete assembly (Complete > Chromosome > Scaffold > Contig) for a taxon
# name or ID.
#
# Usage: genome_downloader.sh <accession|taxon> <output dir> <download log>
#
# Messages follow the pipeline's line contract (ADR-021) on stderr:
#   HH:MM:SS LEVEL genome_downloader <query> | message
# =============================================================================

QUERY="$1"                     # Accession, BioProject, taxon ID or name

# One contract line on stderr; LEVEL is INFO, OK, WARN or ERROR.
say() {
    printf '%s %s genome_downloader %s | %s\n' "$(date +%H:%M:%S)" "$1" "${QUERY:-all}" "$2" >&2
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    say ERROR "usage: $0 <accession|taxon> <output dir> <download log>"
    exit 1
fi

OUTDIR="$(realpath -m "$2")"
LOGFILE="$(realpath -m "$3")"  # one tab-separated line per download

if [ -z "$NCBI_API_KEY" ]; then
    say INFO "no NCBI_API_KEY set; downloads will be slower (export NCBI_API_KEY=...)"
    API_KEY_FLAG=""
else
    API_KEY_FLAG="--api-key $NCBI_API_KEY"
fi

if [[ "$QUERY" =~ ^GC[AF]_[0-9]+(\.[0-9]+)?$ ]]; then
    BEST_ASSEMBLY="$QUERY"
    BEST_LEVEL="Direct_Accession"
else
    BEST_ASSEMBLY=""
    BEST_LEVEL=""
    for LEVEL in "Complete" "Chromosome" "Scaffold" "Contig"; do
        # The client's "New version of client" notice would break the JSON.
        DATASETS_OUTPUT="$(datasets summary genome taxon "$QUERY" $API_KEY_FLAG 2>&1 \
                           | sed '/^New version of client (/d')"

        if echo "$DATASETS_OUTPUT" | grep -q "The taxonomy name"; then
            say ERROR "NCBI does not know the name '$QUERY' (ambiguous or invalid); use an assembly accession in species:"
            exit 1
        fi

        CANDIDATE="$(echo "$DATASETS_OUTPUT" \
                     | jq -r "[.reports[] | select(.assembly_info.assembly_level==\"$LEVEL\")][0].accession")"
        if [ -n "$CANDIDATE" ] && [ "$CANDIDATE" != "null" ]; then
            BEST_ASSEMBLY="$CANDIDATE"
            BEST_LEVEL="$LEVEL"
            say INFO "best assembly is $LEVEL level: $BEST_ASSEMBLY"
            break
        fi
    done

    if [ -z "$BEST_ASSEMBLY" ]; then
        say ERROR "no Complete, Chromosome, Scaffold or Contig assembly found for '$QUERY'"
        exit 1
    fi
fi

ZIPFILE="$(realpath -m "$OUTDIR/genome_${BEST_LEVEL// /_}_${QUERY}.zip")"
FASTA_FILE="$(realpath -m "$OUTDIR/${QUERY}.fa")"

say INFO "downloading $BEST_ASSEMBLY"
if ! datasets download genome accession "$BEST_ASSEMBLY" $API_KEY_FLAG \
        --include genome --assembly-version latest --exclude-atypical \
        --filename "$ZIPFILE"; then
    say ERROR "datasets could not download $BEST_ASSEMBLY; check the network and the accession"
    exit 1
fi

if ! unzip -p "$ZIPFILE" 'ncbi_dataset/data/*/*.fna' > "$FASTA_FILE"; then
    say ERROR "could not unpack $ZIPFILE; the download may be incomplete, delete it and retry"
    exit 1
fi
rm -f "$ZIPFILE"

printf '%s\t%s\t%s\n' "$QUERY" "$BEST_ASSEMBLY" "$BEST_LEVEL" >> "$LOGFILE"
say OK "$BEST_ASSEMBLY ($BEST_LEVEL) downloaded to $FASTA_FILE"
