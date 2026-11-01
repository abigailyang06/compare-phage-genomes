#!/usr/bin/env python3
"""
compare_phage_genbank.py

Compare the CDS / protein content of two annotated phage genomes (GenBank files).

Pipeline
--------
1. Parse CDS features from both GenBank files (nucleotide sequence + protein).
2. Pair up orthologous proteins between the genomes
   (default: reciprocal best hits by protein alignment; or by gene / locus_tag).
3. Globally align each protein pair (BLOSUM62) and walk the alignment codon by codon.
4. Classify every aligned position:
     - identical codon              
     - same amino acid, new codon (SYNONYMOUS difference)
     - different amino acid       (NONSYNONYMOUS difference)
   (Residues that fall in alignment gaps are counted separately as indel residues.)
5. Report per-gene results, genome-wide totals, and a per-site list of differences.

Usage
-----
    python compare_phage_genbank_dr3_(final).py phageA.gb phageB.gb
    python compare_phage_genbank_dr3_(final).py phageA.gb phageB.gb -o results --min-identity 0.4
    python compare_phage_genbank_dr3_(final).py phageA.gb phageB.gb --match-by locus_tag

Requires: biopython >= 1.80
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Align import PairwiseAligner, substitution_matrices


# --------------------------------------------------------------------------- #
# Data structures
# --------------------------------------------------------------------------- #
@dataclass
class CDS:
    uid: str          # unique id used in reports
    locus_tag: str
    gene: str
    product: str
    nt: str           # coding nucleotide sequence (in frame, stop codon removed)
    protein: str      # translated protein (no stop)
    start: int        # 1-based genomic start (for reference)
    end: int
    strand: int


@dataclass
class GeneComparison:
    a: CDS
    b: CDS
    aligned_residues: int = 0
    aa_identical: int = 0          # same amino acid (includes identical codons)
    identical_codons: int = 0      # same amino acid AND same codon
    synonymous: int = 0            # same amino acid, different codon
    nonsynonymous: int = 0         # different amino acid
    synonymous_nt_changes: int = 0
    nonsynonymous_nt_changes: int = 0
    indel_residues: int = 0        # residues in A or B not aligned to anything
    ambiguous_codons: int = 0      # skipped (X or N-containing codons)
    conservative: int = 0          # aa difference with BLOSUM62 score > 0
    nonconservative: int = 0       # aa difference with BLOSUM62 score <= 0

    @property
    def total_nt_changes(self):
        return self.synonymous_nt_changes + self.nonsynonymous_nt_changes

    @property
    def pct_nt_identity(self):
        """Nucleotide identity over aligned (non-gap, non-ambiguous) codons."""
        n = 3 * self.aligned_residues
        return 100.0 * (n - self.total_nt_changes) / n if n else 0.0

    @property
    def pct_aa_identity(self):
        denom = max(len(self.a.protein), len(self.b.protein))
        return 100.0 * self.aa_identical / denom if denom else 0.0


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def load_cds(path, label):
    """Read all CDS features from a GenBank file."""
    cds_list = []
    seen_ids = Counter()
    n_skipped = 0

    for record in SeqIO.parse(path, "genbank"):
        for n, feat in enumerate(record.features):
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            if "pseudo" in q or "pseudogene" in q:
                n_skipped += 1
                continue
            try:
                nt = feat.extract(record.seq)
            except Exception as exc:  # noqa: BLE001
                print(f"[{label}] WARNING: could not extract CDS #{n}: {exc}", file=sys.stderr)
                n_skipped += 1
                continue

            codon_start = int(q.get("codon_start", ["1"])[0])
            table = int(q.get("transl_table", ["11"])[0])  # 11 = bacterial/phage
            nt = nt[codon_start - 1:]
            if len(nt) % 3:
                nt = nt[: len(nt) - (len(nt) % 3)]

            protein = str(nt.translate(table=table))
            nt_str = str(nt).upper()
            # Alternative start codons (GTG/TTG...) are decoded as Met (initiator tRNA),
            # so translate the first codon as M when it is a valid start codon.
            if codon_start == 1 and protein and nt_str[:3] in CodonTable.unambiguous_dna_by_id[table].start_codons:
                protein = "M" + protein[1:]
            if protein.endswith("*"):
                protein = protein[:-1]
                nt_str = nt_str[:-3]
            if "*" in protein or not protein:
                print(f"[{label}] WARNING: skipping CDS with internal stop / empty "
                      f"protein ({q.get('locus_tag', q.get('gene', ['?']))[0]})", file=sys.stderr)
                n_skipped += 1
                continue

            locus_tag = q.get("locus_tag", [""])[0]
            gene = q.get("gene", [""])[0]
            product = q.get("product", [""])[0]
            base = locus_tag or q.get("protein_id", [""])[0] or gene or f"CDS_{feat.location.start + 1}"
            seen_ids[base] += 1
            uid = base if seen_ids[base] == 1 else f"{base}_{seen_ids[base]}"

            cds_list.append(CDS(
                uid=uid, locus_tag=locus_tag, gene=gene, product=product,
                nt=nt_str, protein=protein,
                start=int(feat.location.start) + 1, end=int(feat.location.end),
                strand=feat.location.strand or 1,
            ))

    print(f"[{label}] {path}: {len(cds_list)} CDS loaded"
          + (f" ({n_skipped} skipped)" if n_skipped else ""))
    return cds_list


# --------------------------------------------------------------------------- #
# Alignment helpers
# --------------------------------------------------------------------------- #
def make_aligner(blast_like=False):
    """blast_like=True: local alignment with BLASTP default gap costs
    (existence 11, extension 1 -> a gap of length k costs 11 + k)."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    if blast_like:
        aligner.mode = "local"
        aligner.open_gap_score = -12.0   # first gapped position = 11 + 1
        aligner.extend_gap_score = -1.0
    else:
        aligner.mode = "global"
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5
    return aligner


def sanitize(protein, alphabet):
    """Replace any letter BLOSUM62 doesn't know (e.g. U, O, J) with X."""
    return "".join(c if c in alphabet else "X" for c in protein)


# --------------------------------------------------------------------------- #
# Ortholog pairing
# --------------------------------------------------------------------------- #
def kmer_set(p, k=3):
    return {p[i:i + k] for i in range(len(p) - k + 1)}


def pair_by_alignment(cds_a, cds_b, aligner, top_n=25, max_full_pairs=250_000):
    """One-to-one pairing of proteins by alignment score (greedy best-first).

    1. Score protein pairs with a global BLOSUM62 alignment. For phage-sized
       genomes EVERY A-vs-B pair is scored (no prefilter); only for very large
       inputs is a 3-mer prefilter used to pick candidates.
    2. Sort all scored pairs from best to worst and accept a pair whenever
       neither protein is already taken. Unlike reciprocal-best-hit, a gene whose
       best partner was claimed by a closer relative still gets its next-best partner.

    Returns (pairs, best_a, best_b); best_* map protein id -> (best partner id, score)
    and are used to explain why a protein ended up unpaired.
    """
    alphabet = set(aligner.substitution_matrix.alphabet)
    prot_a = [sanitize(c.protein, alphabet) for c in cds_a]
    prot_b = [sanitize(c.protein, alphabet) for c in cds_b]

    if len(cds_a) * len(cds_b) <= max_full_pairs:
        candidates = {(i, j) for i in range(len(cds_a)) for j in range(len(cds_b))}
    else:
        km_a = [kmer_set(p) for p in prot_a]
        km_b = [kmer_set(p) for p in prot_b]

        def index(kms):
            idx = defaultdict(list)
            for i, ks in enumerate(kms):
                for k in ks:
                    idx[k].append(i)
            return idx

        idx_a, idx_b = index(km_a), index(km_b)
        candidates = set()
        for i, ks in enumerate(km_a):
            hits = Counter(j for k in ks for j in idx_b.get(k, ()))
            candidates.update((i, j) for j, _ in hits.most_common(top_n))
        for j, ks in enumerate(km_b):
            hits = Counter(i for k in ks for i in idx_a.get(k, ()))
            candidates.update((i, j) for i, _ in hits.most_common(top_n))

    self_a = [aligner.score(p, p) for p in prot_a]
    self_b = [aligner.score(p, p) for p in prot_b]

    scored = []
    for i, j in candidates:
        raw = aligner.score(prot_a[i], prot_b[j])
        scored.append((raw / (self_a[i] * self_b[j]) ** 0.5, i, j))
    scored.sort(reverse=True)

    best_a, best_b = {}, {}
    used_a, used_b, pairs = set(), set(), []
    for norm, i, j in scored:
        best_a.setdefault(cds_a[i].uid, (cds_b[j].uid, norm))
        best_b.setdefault(cds_b[j].uid, (cds_a[i].uid, norm))
        if i not in used_a and j not in used_b:
            used_a.add(i)
            used_b.add(j)
            pairs.append((cds_a[i], cds_b[j]))
    return pairs, best_a, best_b


def pair_by_qualifier(cds_a, cds_b, field):
    """Pair genes that share the same locus_tag or gene name."""
    lookup = {getattr(c, field): c for c in cds_b if getattr(c, field)}
    return [(c, lookup[getattr(c, field)]) for c in cds_a
            if getattr(c, field) and getattr(c, field) in lookup]


# --------------------------------------------------------------------------- #
# Codon-level comparison
# --------------------------------------------------------------------------- #
VALID_CODON = set("ACGT")


def compare_pair(a, b, aligner, min_conservative_score=1):
    """Align proteins, walk the alignment, classify each aligned codon.

    Returns (GeneComparison, list_of_site_dicts for sites that differ).
    """
    alphabet = set(aligner.substitution_matrix.alphabet)
    pa, pb = sanitize(a.protein, alphabet), sanitize(b.protein, alphabet)
    aln = aligner.align(pa, pb)[0]

    res = GeneComparison(a=a, b=b)
    sites = []

    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(a1 - a0):
            i, j = a0 + k, b0 + k
            aa_a, aa_b = a.protein[i], b.protein[j]
            cod_a, cod_b = a.nt[3 * i:3 * i + 3], b.nt[3 * j:3 * j + 3]

            if aa_a == "X" or aa_b == "X" or not (set(cod_a) <= VALID_CODON and set(cod_b) <= VALID_CODON):
                res.ambiguous_codons += 1
                continue

            res.aligned_residues += 1
            nt_diff = sum(x != y for x, y in zip(cod_a, cod_b))

            if aa_a == aa_b:
                res.aa_identical += 1
                if cod_a == cod_b:
                    res.identical_codons += 1
                else:
                    res.synonymous += 1
                    res.synonymous_nt_changes += nt_diff
                    kind = "synonymous"
                    sites.append(dict(gene_A=a.uid, gene_B=b.uid, aa_pos_A=i + 1, aa_pos_B=j + 1,
                                      aa_A=aa_a, aa_B=aa_b, codon_A=cod_a, codon_B=cod_b,
                                      nt_changes=nt_diff, type=kind,
                                      blosum62="", aa_change_class=""))
            else:
                res.nonsynonymous += 1
                res.nonsynonymous_nt_changes += nt_diff
                # Default (BLAST "positives" convention): score >= 1 is conservative.
                # Use --conservative-min-score 0 to also count score-0 pairs as conservative.
                score = int(aligner.substitution_matrix[aa_a, aa_b])
                if score >= min_conservative_score:
                    res.conservative += 1
                    cls = "conservative"
                else:
                    res.nonconservative += 1
                    cls = "non-conservative"
                sites.append(dict(gene_A=a.uid, gene_B=b.uid, aa_pos_A=i + 1, aa_pos_B=j + 1,
                                  aa_A=aa_a, aa_B=aa_b, codon_A=cod_a, codon_B=cod_b,
                                  nt_changes=nt_diff, type="nonsynonymous",
                                  blosum62=score, aa_change_class=cls))

    n_aligned_positions = sum(a1 - a0 for a0, a1 in aln.aligned[0])
    res.indel_residues = len(a.protein) + len(b.protein) - 2 * n_aligned_positions
    return res, sites


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
GENE_COLS = [
    "gene_A", "gene_B", "product_A", "product_B", "len_aa_A", "len_aa_B",
    "aligned_codons", "aa_same", "aa_different(nonsyn)", "synonymous_codon_changes",
    "identical_codons", "indel_residues", "ambiguous_codons_skipped",
    "synonymous_nt_changes", "nonsynonymous_nt_changes", "pct_aa_identity",
    "conservative_aa_changes", "nonconservative_aa_changes",
    "total_nt_changes", "pct_nt_identity",
]
SITE_COLS = ["gene_A", "gene_B", "aa_pos_A", "aa_pos_B", "aa_A", "aa_B",
             "codon_A", "codon_B", "nt_changes", "type", "blosum62", "aa_change_class"]


def gene_row(r):
    return [r.a.uid, r.b.uid, r.a.product, r.b.product, len(r.a.protein), len(r.b.protein),
            r.aligned_residues, r.aa_identical, r.nonsynonymous, r.synonymous,
            r.identical_codons, r.indel_residues, r.ambiguous_codons,
            r.synonymous_nt_changes, r.nonsynonymous_nt_changes, f"{r.pct_aa_identity:.2f}",
            r.conservative, r.nonconservative, r.total_nt_changes, f"{r.pct_nt_identity:.2f}"]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genbank_a", help="First GenBank file (.gb)")
    ap.add_argument("genbank_b", help="Second GenBank file (.gb)")
    ap.add_argument("-o", "--output-prefix", default="phage_comparison",
                    help="Prefix for output files (default: phage_comparison)")
    ap.add_argument("--match-by", choices=["alignment", "locus_tag", "gene"], default="alignment",
                    help="How to pair proteins between genomes (default: alignment = reciprocal best hit)")
    ap.add_argument("--min-identity", type=float, default=0.30,
                    help="Minimum fraction amino-acid identity to keep a pair (default: 0.30)")
    ap.add_argument("--conservative-min-score", type=int, default=1,
                    help="Minimum BLOSUM62 score for an amino acid change to count as "
                         "conservative (default: 1, i.e. score > 0; use 0 to include score-0 pairs)")
    ap.add_argument("--blast-like", action="store_true",
                    help="Compare each protein pair with a BLASTP-style local alignment "
                         "(gap open 11 / extend 1) instead of a global alignment. "
                         "Unaligned protein ends are then reported as indel residues.")
    args = ap.parse_args()

    cds_a = load_cds(args.genbank_a, "A")
    cds_b = load_cds(args.genbank_b, "B")
    if not cds_a or not cds_b:
        sys.exit("ERROR: no CDS features found in one of the files.")

    aligner = make_aligner()                              # used for ortholog pairing
    cmp_aligner = make_aligner(blast_like=args.blast_like)  # used for the comparison

    if args.match_by == "alignment":
        pairs, best_a, best_b = pair_by_alignment(cds_a, cds_b, aligner)
    else:
        pairs = pair_by_qualifier(cds_a, cds_b, args.match_by)
        best_a, best_b = {}, {}
    print(f"Paired {len(pairs)} proteins ({args.match_by})")

    results, all_sites = [], []
    rej_a, rej_b = {}, {}   # uid -> (partner uid, aa identity) for pairs failing the identity cutoff
    for a, b in pairs:
        res, sites = compare_pair(a, b, cmp_aligner, args.conservative_min_score)
        ident = res.pct_aa_identity / 100.0
        if ident < args.min_identity:
            rej_a[a.uid] = (b.uid, ident)
            rej_b[b.uid] = (a.uid, ident)
            continue
        results.append(res)
        all_sites.extend(sites)
    if rej_a:
        print(f"Dropped {len(rej_a)} pairs below {args.min_identity:.0%} amino-acid identity "
              f"(see unmatched file; lower --min-identity to keep more divergent pairs)")

    results.sort(key=lambda r: r.a.start)

    paired_a = {r.a.uid for r in results}
    paired_b = {r.b.uid for r in results}
    only_a = [c for c in cds_a if c.uid not in paired_a]
    only_b = [c for c in cds_b if c.uid not in paired_b]

    # ----- write files -----
    gene_path = f"{args.output_prefix}_per_gene.tsv"
    site_path = f"{args.output_prefix}_differences.tsv"
    unmatched_path = f"{args.output_prefix}_unmatched.tsv"
    write_tsv(gene_path, GENE_COLS, [gene_row(r) for r in results])
    write_tsv(site_path, SITE_COLS, [[s[c] for c in SITE_COLS] for s in all_sites])
    holder_a = {r.a.uid: r.b.uid for r in results}   # A id -> its partner in B
    holder_b = {r.b.uid: r.a.uid for r in results}   # B id -> its partner in A

    def why(uid, side):
        best = (best_a if side == "A" else best_b).get(uid)
        rej = (rej_a if side == "A" else rej_b).get(uid)
        taken = holder_b if side == "A" else holder_a
        parts = []
        if rej:
            parts.append(f"assigned {rej[0]} but only {rej[1]:.0%} aa identity")
        if best:
            note = f"best hit {best[0]} (score {best[1]:.2f})"
            if best[0] in taken:
                note += f" is paired with {taken[best[0]]}"
            parts.append(note)
        if not parts:
            parts.append("no partner left / no candidate" if best_a else "no shared " + args.match_by)
        return "; ".join(parts)

    unmatched_rows = [["A", c.uid, c.gene, c.product, len(c.protein), why(c.uid, "A")] for c in only_a] + \
                     [["B", c.uid, c.gene, c.product, len(c.protein), why(c.uid, "B")] for c in only_b]
    write_tsv(unmatched_path, ["genome", "id", "gene", "product", "length_aa", "reason_unpaired"],
              unmatched_rows)

    # ----- totals -----
    tot = lambda attr: sum(getattr(r, attr) for r in results)  # noqa: E731
    aligned = tot("aligned_residues")
    print("\n" + "=" * 62)
    print("GENOME-WIDE SUMMARY (over paired proteins)")
    print("=" * 62)
    print(f"Paired proteins compared             : {len(results)}")
    print(f"Unpaired proteins  A / B             : {len(only_a)} / {len(only_b)}")
    print(f"Aligned codon positions              : {aligned}")
    print(f"  Amino acids SAME                   : {tot('aa_identical')}")
    print(f"      identical codon                : {tot('identical_codons')}")
    print(f"      synonymous differences         : {tot('synonymous')}")
    print(f"  Amino acids DIFFERENT              : {tot('nonsynonymous')}"
          "   (= nonsynonymous differences)")
    print(f"      conservative (BLOSUM62 >= {args.conservative_min_score})   : {tot('conservative')}")
    print(f"      non-conservative (< {args.conservative_min_score})        : {tot('nonconservative')}")
    print(f"Indel residues (in gaps)             : {tot('indel_residues')}")
    print(f"Nucleotide changes, synonymous       : {tot('synonymous_nt_changes')}")
    print(f"Nucleotide changes, nonsynonymous    : {tot('nonsynonymous_nt_changes')}")
    print(f"Nucleotide changes, TOTAL in CDS     : {tot('total_nt_changes')}")
    if aligned:
        print(f"Overall aa identity (aligned sites)  : {100.0 * tot('aa_identical') / aligned:.2f}%")
        print(f"Overall nt identity (coding, aligned): "
              f"{100.0 * (3 * aligned - tot('total_nt_changes')) / (3 * aligned):.2f}%")

    if unmatched_rows:
        print("\nUNPAIRED PROTEINS")
        for g_, uid, _gene, prod, n, why_ in unmatched_rows:
            print(f"  [{g_}] {uid:<24} {n:>4} aa  {prod[:28]:<28} {why_}")
    print(f"\nWrote: {gene_path}\n       {site_path}\n       {unmatched_path}")


if __name__ == "__main__":
    main()
