# -------------------
# DEPENDENCIES
# -------------------

# Suppress messages from loading libraries, so the console remains uncluttered.
suppressMessages({
  library(yaml)           # For reading YAML configuration file
  library(arrow)          # Provides tools for reading and writing Parquet files
  library(tidyverse)      # R packages for data manipulation and visualization
  library(argparse)       # Command-line argument parsing
  library(GenomicRanges)  # Genomic interval operations
  library(plyranges)      # "Tidyverse"-style GRanges operations
  library(rtracklayer)    # Reading/writing genome annotation files (e.g., GFF, BED)
})

# -------------------------------
# 1b. SOURCE PURE PAIRING HELPERS (get_5prime + find_pairs; unit-tested)
# -------------------------------
.resolve_script_dir <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(dirname(ofile))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd_args[grepl("^--file=", cmd_args)]
  if (length(file_arg) > 0L) return(dirname(sub("^--file=", "", file_arg[1])))
  "."
}
.script_dir <- .resolve_script_dir()
source(file.path(.script_dir, "utils", "log.R"))  # line contract, run_main (ADR-021)
source(file.path(.script_dir, "pair_detector", "pairing.R"))

# -------------------------------
# 2. PARSE COMMAND-LINE ARGUMENTS
# -------------------------------
parser <- ArgumentParser(
  description = "Process tBLASTn and LTRdigest integration overlaps"
)

# Define expected command-line arguments
parser$add_argument("--config", required = TRUE, help = "Configuration YAML file")
parser$add_argument("--gff", required = TRUE,
                    help = "GFF3 input: Ranges file to analyse")
parser$add_argument("--parquet_dir", required = TRUE,
                    help = "Directory for pipeline-internal Parquet output")
parser$add_argument("--csv_dir", required = TRUE,
                    help = "Directory for user-facing CSV output")
parser$add_argument("--log", default = NULL,
                    help = "job log file; the Snakemake log: path")

args <- parser$parse_args()
log_job(args$log, "pair_detector")

# main(): sections 3 to 5, run through run_main() so warnings are logged and
# every ending is recorded the same way (ADR-021).
main <- function(args) {
  # Log which tBLASTn file is being processed
  file_basename <- tools::file_path_sans_ext(basename(args$gff))
  log_info("processing pairs for %s", file_basename)

  # -----------------------------
  # 3. CONFIGURATION PARAMETERS
  # -----------------------------
  # Read YAML configuration file
  config <- yaml::read_yaml(args$config)

  # Import thresholds and merging behavior from config
  probe_to_pair <- config$parameters$probe_to_pair

  # -----------------------------
  # 3. GFF3 FILE IMPORT
  # -----------------------------
  # Import the GFF3 file containing validated hits
  ranges <- rtracklayer::import(args$gff, format = "gff3")

  # Assert the existence of the required fields
  required_fields <- c("probe", "label", "virus")
  missing_fields <- setdiff(required_fields, names(mcols(ranges)))

  if (length(missing_fields) > 0) {
    abort_hint(
      sprintf("%s lacks the GFF3 attribute(s) %s", args$gff,
              paste(missing_fields, collapse = ", ")),
      "rebuild the element-hit tracks with ./RetroSeek --ranges-analysis"
    )
  }

  # Assert the existence of the required probe from config within the `probe` field
  if (!any(grepl(probe_to_pair, ranges$probe))) {
    abort_hint(
      sprintf("probe %s (parameters.probe_to_pair) has no hits in %s", probe_to_pair,
              args$gff),
      "set parameters.probe_to_pair to a probe name from the probe CSV"
    )
  }

  # -----------------------------
  # 4. DATA PREPARATION
  # -----------------------------
  # Remove duplicates from ranges
  ranges <- ranges[!duplicated(ranges)]

  # Maximum gap (bp) allowed between paired ranges, from config
  max_gap <- as.integer(config$parameters$pair_max_gap)

  # Pairing logic (get_5prime + find_pairs) lives in pair_detector/pairing.R
  pair_df <- find_pairs(ranges, probe_to_pair, max_gap)


  # ------------------------------
  # 5. PLOT GENERATION
  # ------------------------------
  # ------------------------------
  # EXPORT RESULTS
  # ------------------------------
  write.csv(pair_df, file.path(args$csv_dir, paste0(file_basename, ".csv")),
            row.names = FALSE)
  arrow::write_parquet(pair_df,
                       file.path(args$parquet_dir, paste0(file_basename, ".parquet")))

  log_ok("%s %s pairs", format(nrow(pair_df), big.mark = ","), probe_to_pair)
}

run_main(function() main(args))
