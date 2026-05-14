# attP / attB Finder

A static GitHub Pages webapp for heuristic prediction of candidate phage attP and bacterial host attB sites.

The app is designed around a practical two-step workflow:

1. Upload one or more annotated phage GenBank files in the browser.
2. For large bacterial genomes, upload a prebuilt `host_index.json` instead of the full bacterial GenBank file.

This avoids asking the browser to parse a large bacterial `.gb` file directly.

## What the app does

For each phage GenBank file, the browser searches annotated CDS features for likely integrases. It extracts the integrase gene plus a configurable flank, 750 bp by default. It then compares that phage integrase-proximal region with the host side.

For the host side, the recommended input is `host_index.json`, generated with the included Python script. The host index contains compact tRNA/tmRNA/rRNA neighborhoods and, optionally, broader genome chunks. The browser searches feature neighborhoods first, then optional broader chunks.

The output is a ranked list of candidate shared core sequences with approximate phage and host coordinates, strand, locus label, core length, seed support, and a heuristic prediction score.

This is a prioritization tool, not definitive experimental proof of attP/attB identity.

## Create `host_index.json` from a bacterial GenBank file

Use the interactive script if you want the script to ask for paths:

```bash
python tools/build_host_index_interactive.py
```

For a large bacterial GenBank file, the safest first answers are:

- flank: `5000`
- include broader whole-genome chunks: `no`
- optional GFF/GFF3 path: leave blank
- pretty-print JSON: `no`

That creates a compact first-pass index focused on tRNA/tmRNA/rRNA neighborhoods.

You can also run the command-line version directly:

```bash
python tools/build_host_index_cli.py /path/to/host.gb --flank 5000 --no-global-chunks --out host_index.json
```

To include broader indexed genome chunks:

```bash
python tools/build_host_index_cli.py /path/to/host.gb --flank 5000 --chunk-size 20000 --chunk-step 10000 --out host_index.json
```

If the JSON becomes too large, omit global chunks or generate feature neighborhoods only.

## Browser usage

Open the GitHub Pages site. Drop in:

- phage `.gb`, `.gbk`, or `.genbank` file containing an annotated integrase CDS
- host `host_index.json`

Then click **Run attP/attB prediction**.

Small bacterial GenBank or FASTA files can still be uploaded directly, but the recommended workflow for real bacterial genomes is `host_index.json`.

## Deploy on GitHub Pages

Create a GitHub repository and upload the contents of this folder. Then enable Pages from the repository settings. You can either publish from the root of the main branch or use the included `.github/workflows/pages.yml` workflow.

## Limitations

The tool uses exact k-mer/core matching and heuristic scoring. It does not run BLAST, minimap2, HMMER, tRNAscan-SE, or local alignment refinement. False negatives can occur when attP is far from the annotated integrase, when the integrase is not annotated, when the host genome is not the true host, or when the shared core is short or divergent. False positives can occur around conserved tRNAs, repeats, mobile elements, and low-complexity sequence.

For publication-quality claims, inspect candidate loci manually and validate experimentally.
