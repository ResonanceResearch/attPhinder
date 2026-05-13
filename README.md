# attP / attB Finder

A static GitHub Pages web app for heuristic prediction of phage `attP` and bacterial `attB` sites from uploaded genomes.

## What it does

The app runs completely in the browser:

1. Reads annotated phage GenBank files and bacterial GenBank/FASTA files.
2. Locates likely phage integrase CDS features from GenBank qualifiers such as `integrase`, `int`, `recombinase`, `tyrosine`, or `serine`.
3. Extracts the integrase gene plus a user-defined flank, default ±750 bp.
4. Searches bacterial tRNA regions first, using annotated `tRNA` features plus a user-defined flank, default ±750 bp.
5. Optionally performs a broader genome-wide k-mer search.
6. Reports candidate shared core sequences, coordinates, strand, exact core length, seed support, and a prioritization score.
7. Exports TSV, JSON, and FASTA files.

## Important limitations

This is not a definitive att-site caller. It is a prioritization tool.

It will perform best when:

- the phage GenBank file contains a correctly annotated integrase;
- the real `attP` is close to the integrase;
- the host genome is complete or contains the relevant locus;
- the bacterial GenBank file contains annotated tRNA features;
- the shared att core has enough exact identity to seed k-mer detection.

It can miss sites if the attP is far from the integrase, if the host is wrong, if the bacterial genome is fragmented, or if the shared core is short or divergent. It can produce false positives in repetitive/mobile elements or conserved tRNA regions.

## Deploying on GitHub Pages

1. Create a new GitHub repository.
2. Upload these files to the repository root.
3. In GitHub, go to **Settings → Pages**.
4. Select **Deploy from a branch**.
5. Choose the `main` branch and `/root` folder.
6. Save. GitHub Pages will publish the static app.

Alternatively, the included workflow can deploy the static site when GitHub Pages is configured for GitHub Actions.

## Local use

Open `index.html` directly in a browser, or run a small local server:

```bash
python3 -m http.server 8000
```

Then visit `http://localhost:8000`.

## Suggested validation workflow

For each high-scoring candidate:

1. Inspect whether the core overlaps a tRNA or lies at the tRNA 3′ end.
2. Manually align the phage and bacterial regions around the core.
3. Check whether the predicted core is plausible in orientation relative to the integrase.
4. Compare related phages to see if the same core is conserved.
5. Experimentally validate integration junctions if biological confirmation is required.

## Files

- `index.html` — main app page.
- `src/app.js` — GenBank/FASTA parsing and att-site prediction logic.
- `src/styles.css` — app styling.
- `.github/workflows/pages.yml` — optional GitHub Actions deployment workflow.
- `.nojekyll` — disables Jekyll processing on GitHub Pages.
