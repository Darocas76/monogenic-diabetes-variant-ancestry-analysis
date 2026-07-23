#!/usr/bin/env python3
"""
Provenance script: build the ClinVar classification lookup and the global
ClinVar composition (manuscript Table 1) from the ClinVar GRCh38 VCF.

Input : ClinVar GRCh38 VCF (bgzipped), release 2026-06-28 (accessed 2026-06-29)
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
Output: clinvar_lookup.json  -> {"chr-pos-ref-alt": "VUS|PLP|BLB|Other", ...}
        table1_global.csv     -> global composition (VUS / B-LB / P-LP / Other; total)

CLNSIG mapping (identical rule used for the gene-level and global analyses):
  s = CLNSIG.lower().replace("_", " ")
  "conflicting" in s              -> Other
  "uncertain significance" in s   -> VUS
  else: p = "pathogenic" in s; b = "benign" in s
        p and not b -> PLP;  b and not p -> BLB;  else -> Other

Note: clinvar_lookup.json shipped with this repo is the intersection of this
mapping with the gnomAD v4 target-gene variant keys (merge step 02); running
this script over the full VCF also reproduces the Table 1 global counts.
"""
import gzip, json, csv, re, sys

VCF = sys.argv[1] if len(sys.argv) > 1 else "clinvar.vcf.gz"

def classify(clnsig: str) -> str:
    s = clnsig.lower().replace("_", " ")
    if "conflicting" in s:
        return "Other"
    if "uncertain significance" in s:
        return "VUS"
    p = "pathogenic" in s
    b = "benign" in s
    if p and not b:
        return "PLP"
    if b and not p:
        return "BLB"
    return "Other"

lookup = {}
counts = {"VUS": 0, "BLB": 0, "PLP": 0, "Other": 0}
total = 0
op = gzip.open(VCF, "rt")
for line in op:
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    chrom, pos, _id, ref, alt, _q, _fl, info = f[:8]
    m = re.search(r'CLNSIG=([^;]+)', info)
    if not m:
        continue
    cat = classify(m.group(1))
    total += 1
    counts[cat] += 1
    key = f"{chrom}-{pos}-{ref}-{alt}"   # GRCh38, no 'chr' prefix (matches gnomAD variant_id)
    lookup[key] = cat
op.close()

json.dump(lookup, open("clinvar_lookup.json", "w"))
with open("table1_global.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["classification", "N", "pct"])
    for k, lab in [("VUS", "Uncertain significance (VUS)"),
                   ("BLB", "Benign / Likely benign"),
                   ("PLP", "Pathogenic / Likely pathogenic"),
                   ("Other", "Other / Conflicting")]:
        w.writerow([lab, counts[k], f"{100*counts[k]/total:.1f}%"])
    w.writerow(["Total records", total, "-"])
print(f"records={total}  " + "  ".join(f"{k}={100*v/total:.1f}%" for k, v in counts.items()))
