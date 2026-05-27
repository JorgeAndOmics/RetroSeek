.PHONY: help env env-update lint format format-check typecheck test test-py test-r test-snakemake check clean

ENV_FILE := data/config/environment.yml
ENV_NAME := retroseek
PY_SRC   := workflow/scripts
PY_TESTS := tests workflow/tests

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

# ── environment ──────────────────────────────────────────
env: ## Create the conda/mamba env from data/config/environment.yml
	mamba env create -f $(ENV_FILE) -n $(ENV_NAME) || conda env create -f $(ENV_FILE) -n $(ENV_NAME)

env-update: ## Update the conda/mamba env in place
	mamba env update -f $(ENV_FILE) -n $(ENV_NAME) --prune || conda env update -f $(ENV_FILE) -n $(ENV_NAME) --prune

# ── linting / formatting / typing ───────────────────────
lint: ## Run ruff lint on Python sources
	ruff check $(PY_SRC) $(PY_TESTS)

format: ## Auto-format Python and run ruff --fix
	ruff format $(PY_SRC) $(PY_TESTS)
	ruff check --fix $(PY_SRC) $(PY_TESTS)

format-check: ## Verify formatting without modifying
	ruff format --check $(PY_SRC) $(PY_TESTS)

typecheck: ## mypy strict on Python sources
	mypy $(PY_SRC)

# ── testing ─────────────────────────────────────────────
test: test-py test-r ## Run Python and R tests

test-py: ## pytest
	pytest

test-r: ## R testthat suite
	Rscript -e 'testthat::test_dir("workflow/tests/testthat")'

test-snakemake: ## Dry-run the genome-prep DAG (no probe CSV / network needed)
	snakemake --configfile data/config/config.yaml -n --cores 1 genome_downloader blast_db_generator ltr_index_generator

# ── combined gate ───────────────────────────────────────
check: lint format-check typecheck test ## Run all quality gates (use before commit)

# ── cleanup ─────────────────────────────────────────────
clean: ## Remove Python and R caches
	rm -rf .pytest_cache .mypy_cache .ruff_cache htmlcov .coverage coverage.xml
	find . -type d -name __pycache__ -not -path '*/results/*' -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.pyc' -delete 2>/dev/null || true
	rm -f workflow/scripts/.RData workflow/scripts/.Rhistory
