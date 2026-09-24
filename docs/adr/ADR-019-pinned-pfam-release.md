# ADR-019: Pin the Pfam release and record it

- **Status**: Accepted
- **Date**: 2026-09-23
- **Deciders**: Jorge González García
- **Builds on**: [ADR-015](ADR-015-domain-evidence-by-symmetric-scan.md) (the domain scan and its curated table)

## Context

`pfam_hmm_downloader` fetched `current_release`, whatever Pfam that was on the day,
and recorded nothing about it. Two problems followed:

- **A stale library stopped a run.** The curated class table
  (`data/config/pfam_domain_classes.tsv`) was checked against Pfam 38.2 and names
  families that are only a few releases old (PF29337, PF29688, PF29843, PF29845,
  PF30679, PF30843). A server holding an older download failed in
  `pfam_subset_builder` with a bare list of accessions, and the only way out was to
  delete the library by hand and download it again.
- **Nobody could say which release produced a set of domain calls.** Pfam-A.hmm
  carries no release number.

A third trap sits beside them: `Pfam-A.hmm` is an input of `ltr_digester_setup`.
A fresh download is newer than every LTRdigest output, so Snakemake reruns
LTRdigest on every genome (about a day each). Timestamps alone trigger this, so
`--rerun-triggers mtime` does not help.

## Decision

- **Pin the release.** `input.pfam_release` (default `'38.2'`) selects
  `releases/Pfam<release>/` on the EBI server instead of `current_release`. The
  downloader's outputs are unchanged, so existing installs keep their files.
- **Record it.** A new rule, `pfam_release_recorder`, writes `Pfam.version` beside
  the curated subset: the pinned release, EBI's own version file for it, and
  whether the local download's checksum matches that release. It is demanded only
  by the aggregate targets `domain_scanner` and `taxonomy_classify`, and it is never
  an input of LTRdigest.
- **Explain failures.** When the library lacks curated accessions, the error names
  them, the likely cause (an older release) and the command that fixes it.
  `subset_pfam.missing_accessions` lets the launcher run the same check before any
  stage starts.

## Consequences

- Positive:
  - Every install downloads the same library, so domain calls are reproducible.
  - Each run's domain evidence can be traced to a Pfam release.
  - A stale library is reported with its cause and its fix.
- Negative:
  - On an install whose downloader already ran under Snakemake, the new URLs count
    as changed parameters, so Snakemake wants to rerun the downloader once, and
    LTRdigest after it. Clear that record once, before the next run:

    ```bash
    ./RetroSeek --download-hmm --configfile <config> -skp \
        --cleanup-metadata <accessory>/hmm_profiles/Pfam-A.hmm <accessory>/hmm_profiles/md5sum.txt
    ```

    The files and their timestamps stay as they are; only the parameter record is
    forgotten.
  - Moving to a newer release is a deliberate act (see below).
- Neutral:
  - LTRdigest keeps reading the full `Pfam-A.hmm`, as ADR-015 and ADR-016 require.

## Moving to a newer release

1. Check the curated table against the new library (`subset_pfam.missing_accessions`)
   and update the table if families were retired or merged.
2. Set `input.pfam_release` to the new release.
3. Run `./RetroSeek --download-hmm --forcerun pfam_hmm_downloader`.
4. Keep LTRdigest from rerunning: backdate the new files to before the LTRdigest
   outputs (`touch -d 2000-01-01 Pfam-A.hmm md5sum.txt`), unless LTRdigest should
   deliberately be redone with the new library.

## Alternatives considered

- **Keep `current_release` and only improve the error**: fixes the message but not
  the drift between installs.
- **Add the version file as a downloader output**: would rerun the downloader on
  every existing install and, through it, LTRdigest.
- **Check the release inside the HMM file**: the file carries no release number.

## Revisit trigger

- Pfam retiring or renaming a family in the curated table.
- The domain scan moving off Pfam, or LTRdigest being dropped (then the timestamp
  trap disappears).

## References

- [ADR-015](ADR-015-domain-evidence-by-symmetric-scan.md), [ADR-016](ADR-016-one-domain-classification.md)
- Pfam release archive: https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/
