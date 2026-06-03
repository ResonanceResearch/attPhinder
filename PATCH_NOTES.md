# Patch notes: serine-integrase-aware attP/attB search

## Changed files

- `README.md`
- `index.html`
- `src/app.js`
- `src/styles.css`
- `src/worker.js`
- `tools/build_host_index_cli.py`
- `tools/build_host_index_interactive.py`

## Summary of changes

This patch broadens the tool from a mainly tyrosine-integrase/tRNA-style screen into a more useful mycobacteriophage screen that also handles large-serine-integrase patterns.

Key changes:

1. Added serine-aware host indexing.
   - The host-index builders still index tRNA/tmRNA/rRNA neighborhoods by default.
   - They now also index selected CDS/gene neighborhoods likely to be relevant for serine integrases, especially `groL`, `groEL`, `groEL1`, `cpn60`, chaperonin 60/Hsp60, and glutamyl-tRNA amidotransferase-like genes.
   - This can be disabled with `--no-serine-cds-hotspots`.

2. Added small-core/stem candidate detection.
   - The browser worker still performs the original longer exact-core search.
   - It now also detects short shared cores, default 8 bp, when both the phage-side and host-side contexts have inverted-complement arm support around the core.
   - This was added for patterns like `GAGACCGC` with complementary arms around the core.

3. Added result annotations.
   - Results now include `integraseFamily`, `matchModel`, `serineEvidence`, `hostFeatureClass`, and `candidateReason` in the table and downloads.

4. Updated scoring and filtering.
   - Serine-integrase CDS hotspots receive explicit scoring support.
   - Low-complexity exact matches are filtered to reduce false positives from long A/T or simple-repeat runs.

5. Updated demo.
   - The demo now shows a serine-integrase/groL-style short-core/stem example rather than only a long exact-core tRNA-style example.

## Recommended first-pass host index command

```bash
python tools/build_host_index_cli.py /path/to/host.gb --flank 5000 --no-global-chunks --out host_index.json
```

This keeps global chunks off for compactness but keeps serine CDS hotspots enabled.

If the relevant host target is not a tRNA/tmRNA/rRNA and is not one of the currently recognized CDS hotspots, rebuild with global chunks enabled:

```bash
python tools/build_host_index_cli.py /path/to/host.gb --flank 5000 --chunk-size 20000 --chunk-step 10000 --out host_index.json
```
