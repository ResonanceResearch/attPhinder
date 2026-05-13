# attP / attB Finder

A static GitHub Pages web app for heuristic identification of candidate phage attP and bacterial attB sites from uploaded phage and host genomes.

The app runs entirely in the browser. Genome files are not uploaded to a server.

## What changed in the large-GenBank version

This version is patched for large bacterial GenBank files, including files around 20 MB or larger.

Key changes:

- Heavy parsing and comparison now run in a Web Worker, so the page should remain responsive during large analyses.
- Host/bacterial GenBank files are parsed in host-optimized mode: the app extracts the genome sequence and only `tRNA`/`tmRNA` features, rather than loading every bacterial CDS annotation into memory.
- Phage GenBank files are still parsed for CDS annotations so annotated integrase genes can be detected.
- A Stop button was added for long runs.
- The app warns when a large host file is selected.
- Broader whole-genome search remains optional and should be kept modest for large genomes.

## Intended workflow

1. Upload one or more annotated phage `.gb`, `.gbk`, or `.genbank` files.
2. Upload one or more bacterial host genomes as `.gb/.gbk/.genbank` or FASTA.
3. For host GenBank files, the app first screens regions around annotated bacterial tRNAs/tmRNAs.
4. If enabled, the app then performs a broader k-mer scan across the host genome.
5. Candidate shared core sequences are ranked and can be exported as TSV, JSON, or FASTA.

## Input notes

### Phage files

Annotated GenBank is strongly preferred. The app searches CDS qualifiers such as `gene`, `product`, `note`, `function`, and `locus_tag` for integrase-like terms:

- integrase
- int
- recombinase
- tyrosine
- serine
- site-specific

The integrase CDS plus the configured flanking region is used as the likely attP-containing search window.

### Bacterial files

GenBank is preferred when tRNA-first screening is important. For bacterial GenBank files, only tRNA/tmRNA features are parsed to avoid browser memory problems with large files.

FASTA is also supported, but tRNA-first screening is unavailable unless tRNA annotations are present in a GenBank file. FASTA input can still be used for the broader global search.

## Recommended settings for a 20 MB bacterial GenBank file

Start with:

- Flank around integrase: 750 bp
- Flank around bacterial tRNAs: 750 bp
- k-mer size: 12
- Minimum exact core: 18 bp
- Maximum global candidates: 10–25
- Broader whole-genome search: off for the first run, then on if the tRNA-first search does not identify convincing candidates

If the browser is still slow, use a bacterial FASTA for global-only screening or run the page in Chrome/Edge rather than Safari.

## Limitations

This is a heuristic prioritization tool, not a definitive attP/attB caller. False negatives can occur when:

- the integrase is missing or misannotated;
- attP is not near the integrase;
- the shared core is very short or diverged;
- the supplied bacterium is not the true host;
- the bacterial genome is incomplete;
- the integration site is not near a tRNA/tmRNA.

False positives can occur in conserved tRNAs, repeats, mobile elements, or low-complexity regions.

Candidate sites should be manually inspected and validated experimentally.

## Deploying to GitHub Pages

Upload the contents of this repository to GitHub. Then either:

1. Enable Pages from the repository root, or
2. Use the included `.github/workflows/pages.yml` GitHub Actions workflow.

The app is fully static and does not require a backend.
