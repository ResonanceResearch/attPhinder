#!/usr/bin/env python3
"""
build_host_index.py

Prepare a compact host_index.json file from a bacterial GenBank or FASTA file
for browser-based attP/attB prediction tools.

Purpose
-------
Large bacterial GenBank files are often too heavy for a browser to parse.
This script extracts only the information needed by the webapp:

- genome/contig sequences
- tRNA/tmRNA/rRNA neighborhood sequences and optional serine-integrase CDS hotspot neighborhoods
- feature coordinates and annotations
- optional broader-genome search chunks
- metadata and basic sequence statistics

Typical use
-----------
python build_host_index.py host.gb --flank 5000 --chunk-size 20000 --chunk-step 10000 --out host_index.json

For FASTA-only input:
python build_host_index.py host.fasta --out host_index.json

FASTA-only input will not contain tRNA annotations unless you also provide a GFF file:

python build_host_index.py host.fasta --gff host.gff --out host_index.json

Inputs supported
----------------
- GenBank: .gb, .gbk, .genbank
- FASTA: .fa, .fasta, .fna
- Optional GFF3: .gff, .gff3

No Biopython dependency is required.

Output
------
A JSON file designed to be uploaded into a static browser webapp.

Author: Generated with ChatGPT
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Optional, Tuple


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_text_auto(path: str):
    """Open plain text or gzipped text files."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def clean_dna(seq: str) -> str:
    """Keep only IUPAC-ish DNA letters and uppercase them."""
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1].upper()


def sha256_short(seq: str, n: int = 16) -> str:
    return hashlib.sha256(seq.encode("ascii", errors="ignore")).hexdigest()[:n]


def guess_format(path: str) -> str:
    lower = path.lower()
    if lower.endswith((".gb", ".gbk", ".genbank", ".gb.gz", ".gbk.gz", ".genbank.gz")):
        return "genbank"
    if lower.endswith((".fa", ".fasta", ".fna", ".faa", ".fa.gz", ".fasta.gz", ".fna.gz")):
        return "fasta"
    # Fall back to first non-empty line
    with open_text_auto(path) as handle:
        for line in handle:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                return "fasta"
            if s.startswith("LOCUS"):
                return "genbank"
            break
    raise ValueError(f"Could not guess input format for {path!r}. Use --format.")


def parse_fasta(path: str) -> List[Dict]:
    records = []
    current_id = None
    current_desc = ""
    parts = []

    with open_text_auto(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seq = clean_dna("".join(parts))
                    records.append({
                        "id": current_id,
                        "description": current_desc,
                        "sequence": seq,
                        "length": len(seq),
                        "topology": "unknown",
                    })
                header = line[1:].strip()
                current_desc = header
                current_id = header.split()[0] if header else f"contig_{len(records)+1}"
                parts = []
            else:
                parts.append(line.strip())

    if current_id is not None:
        seq = clean_dna("".join(parts))
        records.append({
            "id": current_id,
            "description": current_desc,
            "sequence": seq,
            "length": len(seq),
            "topology": "unknown",
        })

    if not records:
        raise ValueError("No FASTA records found.")

    return records


def split_genbank_records(text: str) -> List[str]:
    chunks = []
    current = []
    for line in text.splitlines():
        current.append(line)
        if line.strip() == "//":
            chunks.append("\n".join(current))
            current = []
    if current and any(x.strip() for x in current):
        chunks.append("\n".join(current))
    return chunks


def parse_location(raw: str) -> Optional[Dict]:
    """
    Parse common GenBank location strings.

    Supports:
    - 123..456
    - complement(123..456)
    - join(123..200,300..456)
    - complement(join(...))
    - <123..>456

    Returns 1-based inclusive start/end plus strand.
    For complex joins, returns the span and stores parts.
    """
    loc = raw.strip()
    strand = -1 if loc.startswith("complement(") else 1

    # Strip wrappers while preserving coordinate text.
    cleaned = loc
    cleaned = cleaned.replace("complement(", "")
    cleaned = cleaned.replace("join(", "")
    cleaned = cleaned.replace("order(", "")
    cleaned = cleaned.replace(")", "")
    cleaned = cleaned.replace("<", "").replace(">", "")

    nums = [int(x) for x in re.findall(r"\d+", cleaned)]
    if not nums:
        return None

    pairs = []
    if ".." in cleaned or "," in cleaned:
        # approximate pairs from all numbers
        for i in range(0, len(nums) - 1, 2):
            a, b = nums[i], nums[i + 1]
            pairs.append([min(a, b), max(a, b)])
    else:
        # single-base feature
        pairs.append([nums[0], nums[0]])

    start = min(p[0] for p in pairs)
    end = max(p[1] for p in pairs)

    return {
        "raw": raw.strip(),
        "start": start,
        "end": end,
        "strand": strand,
        "parts": pairs,
        "is_complex": len(pairs) > 1,
    }


def parse_qualifiers(lines: List[str]) -> Dict[str, List[str]]:
    quals: Dict[str, List[str]] = {}
    current_key = None
    current_val_parts = []

    def flush():
        nonlocal current_key, current_val_parts
        if current_key is not None:
            val = " ".join(current_val_parts).strip()
            val = val.strip('"')
            quals.setdefault(current_key, []).append(val)
        current_key = None
        current_val_parts = []

    for line in lines:
        s = line.strip()
        if not s.startswith("/"):
            if current_key is not None:
                current_val_parts.append(s.strip('"'))
            continue

        flush()
        body = s[1:]
        if "=" in body:
            key, val = body.split("=", 1)
            current_key = key.strip()
            current_val_parts = [val.strip().strip('"')]
        else:
            current_key = body.strip()
            current_val_parts = ["true"]

    flush()
    return quals


def extract_genbank_sequence(record_text: str) -> str:
    if "\nORIGIN" not in record_text and "\nORIGIN " not in record_text:
        return ""
    origin_part = re.split(r"\nORIGIN\b", record_text, maxsplit=1)[1]
    origin_part = origin_part.split("//")[0]
    letters = []
    for line in origin_part.splitlines():
        letters.append(re.sub(r"[^A-Za-z]", "", line))
    return clean_dna("".join(letters))


def parse_genbank_features(record_text: str, wanted_types: Optional[set] = None) -> List[Dict]:
    """
    Lightweight GenBank FEATURES parser.

    Only returns features in wanted_types if provided.
    """
    lines = record_text.splitlines()
    in_features = False
    feature_blocks = []
    current = None

    feature_start_re = re.compile(r"^     (\S+)\s+(.+)$")

    for line in lines:
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if in_features and line.startswith("ORIGIN"):
            if current:
                feature_blocks.append(current)
            break
        if not in_features:
            continue

        m = feature_start_re.match(line)
        if m:
            if current:
                feature_blocks.append(current)
            ftype = m.group(1)
            loc = m.group(2).strip()
            current = {"type": ftype, "location": loc, "qualifier_lines": []}
        elif current is not None:
            # Continuation line. It may be part of location or qualifier.
            content = line[21:] if len(line) > 21 else line.strip()
            stripped = content.strip()
            if stripped.startswith("/"):
                current["qualifier_lines"].append(stripped)
            elif current["qualifier_lines"]:
                current["qualifier_lines"].append(stripped)
            else:
                current["location"] += stripped

    features = []
    for block in feature_blocks:
        ftype = block["type"]
        if wanted_types and ftype not in wanted_types:
            continue
        loc = parse_location(block["location"])
        if loc is None:
            continue
        quals = parse_qualifiers(block["qualifier_lines"])
        features.append({
            "type": ftype,
            "location": loc,
            "qualifiers": quals,
        })

    return features


def parse_genbank(path: str, wanted_feature_types: Optional[set] = None) -> Tuple[List[Dict], Dict]:
    with open_text_auto(path) as handle:
        text = handle.read()

    gb_records = split_genbank_records(text)
    records = []
    meta = {"organism": None, "source": None, "accessions": []}

    for i, rec_text in enumerate(gb_records, start=1):
        locus_match = re.search(r"^LOCUS\s+(\S+).*$", rec_text, flags=re.MULTILINE)
        contig_id = locus_match.group(1) if locus_match else f"record_{i}"

        topology = "circular" if re.search(r"^LOCUS\s+.*\bcircular\b", rec_text, flags=re.MULTILINE) else "linear"

        definition_match = re.search(r"^DEFINITION\s+(.+?)(?=\n[A-Z][A-Z ]+\s|\Z)", rec_text, flags=re.MULTILINE | re.DOTALL)
        definition = " ".join(definition_match.group(1).split()) if definition_match else ""

        accession_match = re.search(r"^ACCESSION\s+(.+)$", rec_text, flags=re.MULTILINE)
        if accession_match:
            meta["accessions"].extend(accession_match.group(1).split())

        organism_match = re.search(r"^\s+ORGANISM\s+(.+)$", rec_text, flags=re.MULTILINE)
        if organism_match and not meta["organism"]:
            meta["organism"] = organism_match.group(1).strip()

        source_match = re.search(r"^SOURCE\s+(.+)$", rec_text, flags=re.MULTILINE)
        if source_match and not meta["source"]:
            meta["source"] = source_match.group(1).strip()

        seq = extract_genbank_sequence(rec_text)
        features = parse_genbank_features(rec_text, wanted_feature_types)

        records.append({
            "id": contig_id,
            "description": definition,
            "sequence": seq,
            "length": len(seq),
            "topology": topology,
            "features": features,
        })

    if not records:
        raise ValueError("No GenBank records found.")

    return records, meta


def parse_gff_attributes(attr_text: str) -> Dict[str, str]:
    attrs = {}
    for part in attr_text.split(";"):
        if not part.strip():
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
        else:
            k, v = part, ""
        attrs[k.strip()] = v.strip().strip('"')
    return attrs


def parse_gff_features(path: str, wanted_types: set) -> List[Dict]:
    features = []
    with open_text_auto(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            if ftype not in wanted_types:
                continue
            q = parse_gff_attributes(attrs)
            features.append({
                "contig_id": seqid,
                "type": ftype,
                "location": {
                    "raw": f"{start}..{end}",
                    "start": int(start),
                    "end": int(end),
                    "strand": -1 if strand == "-" else 1,
                    "parts": [[int(start), int(end)]],
                    "is_complex": False,
                },
                "qualifiers": {
                    "product": [q.get("product") or q.get("Name") or q.get("gene") or q.get("ID") or ""],
                    "gene": [q.get("gene", "")],
                    "locus_tag": [q.get("locus_tag") or q.get("ID") or ""],
                    "note": [q.get("Note", "")],
                },
            })
    return features


def circular_slice(seq: str, start0: int, end0: int, circular: bool) -> str:
    """
    Return slice using 0-based half-open coordinates.
    If circular is true, wrap around contig ends.
    """
    n = len(seq)
    if n == 0:
        return ""
    if not circular:
        start0 = max(0, start0)
        end0 = min(n, end0)
        if end0 <= start0:
            return ""
        return seq[start0:end0]

    # circular
    length = end0 - start0
    if length <= 0:
        return ""
    if length >= n:
        return seq
    start0 %= n
    end0 = start0 + length
    if end0 <= n:
        return seq[start0:end0]
    return seq[start0:] + seq[:end0 - n]


def feature_label(feature: Dict) -> str:
    q = feature.get("qualifiers", {})
    for key in ["product", "gene", "locus_tag", "note"]:
        vals = q.get(key) or []
        for v in vals:
            if v:
                return v
    return feature.get("type", "feature")



def feature_text(feature: Dict) -> str:
    q = feature.get("qualifiers", {})
    parts = [feature.get("type", "")]
    for key in ["product", "gene", "locus_tag", "note", "function", "standard_name"]:
        value = q.get(key) or []
        if isinstance(value, str):
            parts.append(value)
        else:
            parts.extend(str(v) for v in value if v)
    return " ".join(parts).lower().replace("_", " ").replace("-", " ")


def serine_cds_hotspot_reason(feature: Dict) -> Optional[str]:
    """Return a reason if a CDS/gene feature is a plausible serine-integrase attB hotspot."""
    if feature.get("type") not in {"CDS", "gene"}:
        return None
    text = feature_text(feature)
    if re.search(r"\b(groel1?|grol1?|cpn60|hsp60|chaperonin\s*60|60\s*kda\s+chaperonin|gro\s*el)\b", text):
        return "serine-integrase CDS hotspot: groL/groEL/chaperonin 60"
    if re.search(r"\b(glutamyl\s*trna|glutamyl\s+trna|gln\s*trna|gat[abc]|amidotransferase)\b", text) and re.search(r"\b(amidotransferase|glutamyl|gln|gat[abc])\b", text):
        return "reported serine-integrase/vector attB hotspot: glutamyl-tRNA amidotransferase-related CDS"
    return None


def feature_is_indexed(feature: Dict, feature_types: set, include_serine_cds_hotspots: bool) -> bool:
    if feature.get("type") in feature_types:
        return True
    return include_serine_cds_hotspots and serine_cds_hotspot_reason(feature) is not None


def build_neighborhoods(records: List[Dict], feature_types: set, flank: int, max_neighborhood_bp: int, include_serine_cds_hotspots: bool = True) -> List[Dict]:
    neighborhoods = []
    feature_counter = 0

    for rec in records:
        seq = rec.get("sequence", "")
        n = len(seq)
        circular = rec.get("topology") == "circular"
        for feat in rec.get("features", []):
            if not feature_is_indexed(feat, feature_types, include_serine_cds_hotspots):
                continue
            candidate_reason = serine_cds_hotspot_reason(feat) or ""
            loc = feat["location"]
            start = int(loc["start"])
            end = int(loc["end"])
            strand = int(loc.get("strand", 1))

            # Convert 1-based inclusive to 0-based half-open and add flank.
            win_start0 = start - 1 - flank
            win_end0 = end + flank
            raw_seq = circular_slice(seq, win_start0, win_end0, circular=circular)

            # Safety cap in case someone gives very large flanks or odd annotations.
            if len(raw_seq) > max_neighborhood_bp:
                extra = len(raw_seq) - max_neighborhood_bp
                trim_left = extra // 2
                raw_seq = raw_seq[trim_left:trim_left + max_neighborhood_bp]

            feature_counter += 1
            nid = f"{rec['id']}:{feat['type']}:{start}-{end}:{feature_counter}"

            neighborhoods.append({
                "id": nid,
                "contig_id": rec["id"],
                "feature_type": feat["type"],
                "feature_label": feature_label(feat),
                "feature_start": start,
                "feature_end": end,
                "feature_strand": strand,
                "window_start_approx": max(1, start - flank),
                "window_end_approx": min(n, end + flank),
                "window_is_circular_wrapped": circular and (win_start0 < 0 or win_end0 > n),
                "sequence": raw_seq,
                "length": len(raw_seq),
                "sequence_sha256_16": sha256_short(raw_seq),
                "qualifiers": feat.get("qualifiers", {}),
                "candidate_reason": candidate_reason,
            })

    return neighborhoods


def build_global_chunks(records: List[Dict], chunk_size: int, chunk_step: int, include_sequence: bool) -> List[Dict]:
    chunks = []
    if chunk_size <= 0 or chunk_step <= 0:
        return chunks

    for rec in records:
        seq = rec.get("sequence", "")
        n = len(seq)
        if n == 0:
            continue

        # Ensure at least one chunk for short contigs.
        if n <= chunk_size:
            chunks.append({
                "id": f"{rec['id']}:1-{n}",
                "contig_id": rec["id"],
                "start": 1,
                "end": n,
                "length": n,
                "sequence": seq if include_sequence else None,
                "sequence_sha256_16": sha256_short(seq),
            })
            continue

        pos = 0
        while pos < n:
            end = min(n, pos + chunk_size)
            chunk_seq = seq[pos:end]
            chunks.append({
                "id": f"{rec['id']}:{pos+1}-{end}",
                "contig_id": rec["id"],
                "start": pos + 1,
                "end": end,
                "length": len(chunk_seq),
                "sequence": chunk_seq if include_sequence else None,
                "sequence_sha256_16": sha256_short(chunk_seq),
            })
            if end == n:
                break
            pos += chunk_step

    return chunks


def gc_fraction(seq: str) -> Optional[float]:
    clean = re.sub(r"[^ACGT]", "", seq.upper())
    if not clean:
        return None
    return round((clean.count("G") + clean.count("C")) / len(clean), 5)


def summarize_records(records: List[Dict]) -> List[Dict]:
    out = []
    for rec in records:
        seq = rec.get("sequence", "")
        out.append({
            "id": rec.get("id"),
            "description": rec.get("description", ""),
            "length": len(seq),
            "topology": rec.get("topology", "unknown"),
            "gc_fraction": gc_fraction(seq),
            "sequence_sha256_16": sha256_short(seq),
        })
    return out


def remove_null_sequence_fields(obj):
    """
    Recursively remove sequence:null keys from chunks to keep JSON tidy.
    """
    if isinstance(obj, dict):
        return {
            k: remove_null_sequence_fields(v)
            for k, v in obj.items()
            if not (k == "sequence" and v is None)
        }
    if isinstance(obj, list):
        return [remove_null_sequence_fields(x) for x in obj]
    return obj


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build compact host_index.json for attP/attB webapp searches."
    )
    parser.add_argument("input", help="Host bacterial GenBank or FASTA file. Gzipped files are accepted.")
    parser.add_argument("--format", choices=["genbank", "fasta"], default=None, help="Input format. Default: auto-detect.")
    parser.add_argument("--gff", default=None, help="Optional GFF/GFF3 file with host feature annotations, useful with FASTA input.")
    parser.add_argument("--out", default="host_index.json", help="Output JSON file. Default: host_index.json")
    parser.add_argument("--flank", type=int, default=5000, help="Flank around indexed host features. Default: 5000 bp.")
    parser.add_argument("--feature-types", default="tRNA,tmRNA,rRNA", help="Comma-separated host features to index. Default: tRNA,tmRNA,rRNA")
    parser.add_argument("--no-serine-cds-hotspots", action="store_true", help="Do not add CDS/gene neighborhoods for serine-integrase hotspots such as groL/groEL and glutamyl-tRNA amidotransferase.")
    parser.add_argument("--chunk-size", type=int, default=20000, help="Broader-search chunk size. Default: 20000 bp.")
    parser.add_argument("--chunk-step", type=int, default=10000, help="Broader-search chunk step. Default: 10000 bp.")
    parser.add_argument("--no-global-chunks", action="store_true", help="Do not include broader-search chunks.")
    parser.add_argument("--omit-global-sequences", action="store_true", help="Include global chunk metadata only, not sequence. Smaller JSON but disables broad browser search.")
    parser.add_argument("--max-neighborhood-bp", type=int, default=25000, help="Safety cap for each feature neighborhood sequence. Default: 25000 bp.")
    parser.add_argument("--organism", default=None, help="Override organism name in metadata.")
    parser.add_argument("--name", default=None, help="Friendly host name for metadata.")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON. Larger but human-readable.")

    args = parser.parse_args(argv)

    if not os.path.exists(args.input):
        raise FileNotFoundError(args.input)

    input_format = args.format or guess_format(args.input)
    wanted_feature_types = {x.strip() for x in args.feature_types.split(",") if x.strip()}
    parse_feature_types = set(wanted_feature_types)
    if not args.no_serine_cds_hotspots:
        parse_feature_types.update({"CDS", "gene"})

    metadata = {
        "name": args.name or os.path.basename(args.input),
        "organism": args.organism,
        "source_file": os.path.basename(args.input),
        "source_format": input_format,
    }

    if input_format == "genbank":
        records, gb_meta = parse_genbank(args.input, parse_feature_types)
        metadata["organism"] = args.organism or gb_meta.get("organism") or gb_meta.get("source")
        metadata["accessions"] = gb_meta.get("accessions", [])
    else:
        records = parse_fasta(args.input)
        metadata["accessions"] = []

    # If GFF is provided, add feature annotations to matching contigs.
    if args.gff:
        gff_features = parse_gff_features(args.gff, parse_feature_types)
        by_contig = {}
        for feat in gff_features:
            by_contig.setdefault(feat["contig_id"], []).append(feat)

        for rec in records:
            rec.setdefault("features", [])
            rec["features"].extend(by_contig.get(rec["id"], []))

        # If exact IDs don't match, try a forgiving match using first token.
        unmatched = sum(len(v) for k, v in by_contig.items() if k not in {r["id"] for r in records})
        if unmatched:
            id_map = {r["id"].split()[0]: r for r in records}
            for contig_id, feats in by_contig.items():
                if contig_id in {r["id"] for r in records}:
                    continue
                short = contig_id.split()[0]
                if short in id_map:
                    id_map[short].setdefault("features", []).extend(feats)

    neighborhoods = build_neighborhoods(
        records=records,
        feature_types=wanted_feature_types,
        flank=args.flank,
        max_neighborhood_bp=args.max_neighborhood_bp,
        include_serine_cds_hotspots=not args.no_serine_cds_hotspots,
    )

    include_global = not args.no_global_chunks
    global_chunks = []
    if include_global:
        global_chunks = build_global_chunks(
            records=records,
            chunk_size=args.chunk_size,
            chunk_step=args.chunk_step,
            include_sequence=not args.omit_global_sequences,
        )

    total_bp = sum(len(r.get("sequence", "")) for r in records)

    host_index = {
        "schema": "att-site-host-index-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "metadata": metadata,
        "parameters": {
            "flank": args.flank,
            "feature_types": sorted(wanted_feature_types),
            "serine_cds_hotspots": not args.no_serine_cds_hotspots,
            "chunk_size": args.chunk_size,
            "chunk_step": args.chunk_step,
            "global_chunks_included": include_global,
            "global_chunk_sequences_included": include_global and not args.omit_global_sequences,
            "max_neighborhood_bp": args.max_neighborhood_bp,
        },
        "summary": {
            "contig_count": len(records),
            "total_bp": total_bp,
            "gc_fraction": gc_fraction("".join(r.get("sequence", "") for r in records)),
            "feature_neighborhood_count": len(neighborhoods),
            "global_chunk_count": len(global_chunks),
        },
        "contigs": summarize_records(records),
        "feature_neighborhoods": neighborhoods,
        "global_chunks": global_chunks,
    }

    host_index = remove_null_sequence_fields(host_index)

    indent = 2 if args.pretty else None
    with open(args.out, "w", encoding="utf-8") as out:
        json.dump(host_index, out, indent=indent, separators=None if args.pretty else (",", ":"))

    size_mb = os.path.getsize(args.out) / (1024 * 1024)

    print(f"Created: {args.out}")
    print(f"Input format: {input_format}")
    print(f"Contigs: {len(records)}")
    print(f"Total bp: {total_bp:,}")
    print(f"Feature neighborhoods: {len(neighborhoods)}")
    print(f"Global chunks: {len(global_chunks)}")
    print(f"Output size: {size_mb:.2f} MB")

    if len(neighborhoods) == 0:
        print()
        print("WARNING: No tRNA/tmRNA/rRNA or serine-CDS-hotspot neighborhoods were found.")
        print("If you used FASTA input, provide a GFF file with --gff.")
        print("If you used GenBank input, check whether feature types and annotations include tRNA/tmRNA/rRNA, groL/groEL, or glutamyl-tRNA amidotransferase-like CDS entries.")

    if size_mb > 50:
        print()
        print("WARNING: host_index.json is large.")
        print("Consider using --no-global-chunks or --omit-global-sequences for a smaller file.")
        print("For tyrosine-integrase first-pass prediction, tRNA/tmRNA/rRNA neighborhoods are often sufficient; for serine integrases, keep serine CDS hotspots enabled and consider global chunks if the host target is unknown.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Cancelled.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(1)
