# attP / attB Finder

A static GitHub Pages webapp for heuristic prediction of candidate phage attP and bacterial host attB sites.

The app is designed around a practical two-step workflow:

1. Upload one or more annotated phage GenBank files in the browser.
2. For large bacterial genomes, upload a prebuilt `host_index.json` instead of the full bacterial GenBank file.

This avoids asking the browser to parse a large bacterial `.gb` file directly.

## What the app does

For each phage GenBank file, the browser searches annotated CDS features for likely integrases. It extracts the integrase gene plus a configurable flank, 750 bp by default. This intentionally includes attP sites that sit adjacent to the integrase and attP sites that overlap or fall inside the integrase CDS.

For the host side, the recommended input is `host_index.json`, generated with the included Python script. The host index contains compact tRNA/tmRNA/rRNA neighborhoods, optional serine-integrase CDS hotspot neighborhoods such as `groL`/`groEL`/chaperonin 60 and glutamyl-tRNA amidotransferase-like genes, and optionally broader genome chunks. The browser searches feature neighborhoods first, then optional broader chunks.

The output is a ranked list of candidate shared core sequences with approximate phage and host coordinates, strand, locus label, core length, seed support, match model, stem evidence when relevant, integrase-family model fit, and a heuristic prediction score.

This is a prioritization tool, not definitive experimental proof of attP/attB identity.

## Tyrosine versus serine integrase support

The original version was most useful for tyrosine-integrase-like searches because it emphasized longer exact shared cores near tRNA/tmRNA/rRNA loci. This patched version adds serine-aware and family-aware behaviour:

- Host indexing now optionally includes CDS hotspots for serine integrase systems, especially `groL`/`groEL`/chaperonin 60 and glutamyl-tRNA amidotransferase-like genes. This matters because not all mycobacteriophage attB sites are near tRNAs.
- The web worker can now run a small-core/stem search. This mode allows short shared cores, default 8 bp, but only upgrades them when both the phage-side and host-side sites have local inverted-complement arm support around that core. This is intended for cases with a short shared core and complementary arms on either side; no particular core sequence is assumed.
- The browser now asks for the integrase type. If you select **Known serine integrase**, the app suppresses tRNA/tmRNA/rRNA-biased candidates unless they have unusually strong exact-core evidence. If you select **Known tyrosine integrase**, the app suppresses serine-CDS/stem-biased candidates unless the exact-core evidence is unusually strong. **Auto** uses the GenBank annotation when it explicitly identifies serine or tyrosine.

Long exact-core matching is still retained and remains preferable for tyrosine-integrase candidates. Use **Unknown / run both models without prior** when the integrase family is genuinely unclear.

## Create `host_index.json` from a bacterial GenBank file

Use the interactive script if you want the script to ask for paths:

```bash
python tools/build_host_index_interactive.py
```

For a large bacterial GenBank file, the safest first answers are:

- flank: `5000`
- include serine-integrase CDS hotspots: `yes`
- include broader whole-genome chunks: `no`
- optional GFF/GFF3 path: leave blank
- pretty-print JSON: `no`

That creates a compact first-pass index focused on tRNA/tmRNA/rRNA neighborhoods and selected CDS hotspots such as `groL`/`groEL`.

You can also run the command-line version directly:

```bash
python tools/build_host_index_cli.py /path/to/host.gb --flank 5000 --no-global-chunks --out host_index.json
```

To disable CDS hotspot indexing:

```bash
python tools/build_host_index_cli.py /path/to/host.gb --flank 5000 --no-serine-cds-hotspots --no-global-chunks --out host_index.json
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

Before clicking **Run attP/attB prediction**, choose the integrase type if you know it. This is especially useful for serine integrases because otherwise short-core/stem logic can produce plausible but biologically off-model tRNA hits. Leave the setting at **Auto** if the phage GenBank annotation clearly says serine or tyrosine, or choose **Unknown / run both models without prior** if you do not want family-based filtering.

Small bacterial GenBank or FASTA files can still be uploaded directly, but the recommended workflow for real bacterial genomes is `host_index.json`.

## Deploy on GitHub Pages

Create a GitHub repository and upload the contents of this folder. Then enable Pages from the repository settings. You can either publish from the root of the main branch or use the included `.github/workflows/pages.yml` workflow.

## Limitations

The tool uses exact k-mer/core matching, a simple serine small-core/stem model, and heuristic scoring. It does not run BLAST, minimap2, HMMER, tRNAscan-SE, covariance models, structural RNA prediction, or local alignment refinement. False negatives can occur when attP is far from the annotated integrase, when the integrase is not annotated, when the host genome is not the true host, when the shared core is very short without clear flanking-arm support, or when the relevant host CDS hotspot is not indexed. False positives can occur around repeats, conserved genes, conserved RNAs, mobile elements, and low-complexity sequence.

For publication-quality claims, inspect candidate loci manually and validate experimentally.
