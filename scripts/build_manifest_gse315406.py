#!/usr/bin/env python3
"""
build_manifest_gse315406.py
Fetches all 228 sample records for GSE315406 (ENCORE) from NCBI SRA,
groups them into experiments, and writes config/manifest_gse315406.csv.

Usage (run from the eclip_pipeline root on Puma):
    python3 scripts/build_manifest_gse315406.py

Requirements: Python 3.6+, standard library only (urllib, xml).
Tested with the eclip conda environment.
"""
import re
import sys
import time
import json
import urllib.request
import urllib.error
from pathlib import Path
from xml.etree import ElementTree as ET
from collections import defaultdict

# ── Config ────────────────────────────────────────────────────────────────────
GEO_ACC       = "GSE315406"
EUTILS        = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ADAPTER_INPUT = "/home/u11/maarowosegbe/eclip_pipeline/config/adapters/eclip_input_adapter.fasta"
ADAPTER_IP    = "/home/u11/maarowosegbe/eclip_pipeline/config/adapters/eclip_ip_adapter.fasta"
FASTQ_DIR     = "/xdisk/haining/maarowosegbe/eclip_pipeline/data/fastq/gse315406"
OUT_MANIFEST  = Path(__file__).parent.parent / "config" / "manifest_gse315406.csv"
BATCH_SIZE    = 50   # IDs per efetch request
SLEEP_SEC     = 0.4  # polite rate-limiting


def fetch(url, retries=4):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8", errors="replace")
        except urllib.error.URLError as e:
            if attempt == retries - 1:
                raise
            wait = 2 ** attempt
            print(f"  [retry {attempt+1}] {e} — waiting {wait}s", file=sys.stderr)
            time.sleep(wait)


def esearch_ids(db, term, retmax=300):
    url = f"{EUTILS}/esearch.fcgi?db={db}&term={term}&retmax={retmax}&retmode=json"
    data = json.loads(fetch(url))
    return data["esearchresult"]["idlist"]


# ── Sample title format ───────────────────────────────────────────────────────
# "YTHDC1 (5096) eCLIP in K562 cells, IP replicate 2"
# "YTHDC1 (5096) eCLIP in K562 cells, INPUT replicate 2"
TITLE_RE = re.compile(
    r"^(\S+)\s+\((\d+)\)\s+eCLIP\s+in\s+(\S+)\s+cells,\s+(IP|INPUT)\s+replicate\s+(\d+)$",
    re.IGNORECASE,
)


def parse_sra_xml_batch(xml_text):
    """
    Parse a batch of SRA XML records. Returns a list of dicts with keys:
    srx, srr, rbp, enc_id, cell, stype ('IP'/'INPUT'), rep (int).
    """
    records = []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        print(f"  [warn] XML parse error: {e}", file=sys.stderr)
        return records

    for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
        exp = pkg.find("EXPERIMENT")
        run_set = pkg.find("RUN_SET")
        sample = pkg.find(".//SAMPLE")

        if exp is None or run_set is None:
            continue

        srx = exp.get("accession", "")
        # First RUN accession
        run = run_set.find("RUN")
        srr = run.get("accession", "") if run is not None else ""

        # Sample title lives in SAMPLE > TITLE
        title = ""
        if sample is not None:
            title_el = sample.find("TITLE")
            if title_el is not None and title_el.text:
                title = title_el.text.strip()

        # Fallback: LIBRARY_NAME
        if not title:
            ln = exp.find(".//LIBRARY_NAME")
            if ln is not None and ln.text:
                title = ln.text.strip()

        m = TITLE_RE.match(title)
        if not m:
            print(f"  [skip] unrecognised title: {title!r}", file=sys.stderr)
            continue

        rbp, enc_id, cell, stype, rep_str = m.groups()
        records.append({
            "srx":    srx,
            "srr":    srr,
            "rbp":    rbp,
            "enc_id": enc_id,
            "cell":   cell,
            "stype":  stype.upper(),
            "rep":    int(rep_str),
        })
    return records


def main():
    print(f"[1/3] Finding SRA IDs for {GEO_ACC} ...")
    sra_ids = esearch_ids("sra", GEO_ACC, retmax=300)
    print(f"      {len(sra_ids)} SRA experiments found")

    print("[2/3] Fetching SRA metadata in batches ...")
    all_samples = []
    for i in range(0, len(sra_ids), BATCH_SIZE):
        chunk = sra_ids[i : i + BATCH_SIZE]
        url = f"{EUTILS}/efetch.fcgi?db=sra&id={','.join(chunk)}&retmode=xml"
        xml_text = fetch(url)
        batch_records = parse_sra_xml_batch(xml_text)
        all_samples.extend(batch_records)
        print(f"      batch {i//BATCH_SIZE + 1}: {len(batch_records)} records parsed")
        time.sleep(SLEEP_SEC)

    print(f"      total samples parsed: {len(all_samples)}")
    if len(all_samples) != len(sra_ids):
        print(f"  [warn] expected {len(sra_ids)}, got {len(all_samples)}", file=sys.stderr)

    print("[3/3] Building manifest ...")

    # Group: (rbp, enc_id, cell) → {ip: {rep: record}, input: {rep: record}}
    experiments = defaultdict(lambda: {"IP": {}, "INPUT": {}})
    for s in all_samples:
        key = (s["rbp"], s["enc_id"], s["cell"])
        experiments[key][s["stype"]][s["rep"]] = s

    # Write CSV
    header = (
        "Experiment,Sample,Cells,"
        "Input_replicate,Input_adapter,Input_fastq,"
        "CLIP_replicate,CLIP_adapter,CLIP_fastq,Notes"
    )
    rows = []
    for (rbp, enc_id, cell), exp in sorted(experiments.items()):
        # Experiment name: RBP_CELL_ENCID  e.g. YTHDC1_K562_5096
        exp_name = f"{rbp}_{cell}_{enc_id}"
        all_reps = sorted(set(list(exp["IP"].keys()) + list(exp["INPUT"].keys())))
        for rep_n in all_reps:
            ip_s  = exp["IP"].get(rep_n)
            inp_s = exp["INPUT"].get(rep_n)
            if ip_s is None or inp_s is None:
                print(
                    f"  [warn] incomplete rep{rep_n} for {exp_name}: "
                    f"IP={'yes' if ip_s else 'NO'}, INPUT={'yes' if inp_s else 'NO'}",
                    file=sys.stderr,
                )
                continue
            ip_fastq  = f"{FASTQ_DIR}/{exp_name}_ip_rep{rep_n}.fastq.gz"
            inp_fastq = f"{FASTQ_DIR}/{exp_name}_input_rep{rep_n}.fastq.gz"
            notes     = f"{ip_s['srr']}/{inp_s['srr']}"
            rows.append(",".join([
                exp_name, exp_name, f"{cell}_cells",
                str(rep_n), ADAPTER_INPUT, inp_fastq,
                str(rep_n), ADAPTER_IP,    ip_fastq,
                notes,
            ]))

    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_MANIFEST, "w") as f:
        f.write(f"# eCLIP sample manifest — ENCORE, {GEO_ACC}\n")
        f.write(f"# 57 RBPs (K562 and HepG2 cells), ENCODE4 eCLIP, 10 nt UMI, single-end\n")
        f.write(f"# Auto-generated by scripts/build_manifest_gse315406.py\n")
        f.write(f"#\n")
        f.write(header + "\n")
        for row in rows:
            f.write(row + "\n")

    print(f"Wrote {len(rows)} rows ({len(experiments)} experiments) → {OUT_MANIFEST}")

    # Print the SRR accessions needed for the download script
    print("\n# SRR accessions for download script:")
    for s in sorted(all_samples, key=lambda x: (x["rbp"], x["enc_id"], x["cell"], x["stype"], x["rep"])):
        exp_name = f"{s['rbp']}_{s['cell']}_{s['enc_id']}"
        stype_lbl = "ip" if s["stype"] == "IP" else "input"
        label = f"{exp_name}_{stype_lbl}_rep{s['rep']}"
        print(f"download_se {label} {s['srr']}")


if __name__ == "__main__":
    main()
