# mcr-chromosomal-origins

Analysis code and intermediate results for a study of chromosomal phosphoethanolamine
transferases as repeated and independent sources of mobilised colistin resistance in
*Enterobacterales*.

The question is where each *mcr* family came from, whether the families share one origin,
and which chromosomal genes sit close enough to the mobilised ones to be worth watching.
Everything here reproduces the figures, tables and numbers in the manuscript.

## What the analysis does

1. Retrieves the reference allele of every *mcr* family with an entry in the NCBI
   Bacterial Antimicrobial Resistance Reference Gene Catalog, together with chromosomal
   phosphoethanolamine transferases from RefSeq across 44 genera and from UniProtKB for
   the genera implicated as donors.
2. Aligns them with FAMSA, trims with trimAl, and infers a maximum-likelihood tree with
   IQ-TREE under a model chosen by ModelFinder, with 1000 ultrafast bootstrap replicates
   and 1000 SH-aLRT replicates.
3. Tests monophyly of the mobilised families against the unconstrained tree with the
   approximately unbiased test.
4. Reconstructs ancestral states at the node each mobilised lineage shares with its
   chromosomal relatives, separating the positions that are invariant across the whole
   dataset from those that carry information.
5. Estimates synonymous and non-synonymous divergence with codeml, first between families,
   where synonymous sites turn out to be saturated, then within families, where they are not.
6. Examines genomic context around six mobilised deposits and four chromosomal donor loci.
7. Screens the chromosomal proteins against profile hidden Markov models built from the
   catalogued allele sets, to test whether reference-based surveillance would find them.

## Headline results

| Result | Value |
| --- | --- |
| Best-fit model, primary alignment | LG+F+I+R6 |
| Unconstrained ML tree | lnL = -44382.42 |
| Tree constrained to *mcr* monophyly | lnL = -45207.22, deltaL = 824.8 |
| Approximately unbiased test | p = 8.1 x 10<sup>-6</sup>, monophyly rejected |
| Independent mobilised lineages | 9 in the 161-sequence analysis, 10 in the 191-sequence rooted analysis |
| Synonymous saturation between families | dS > 1 in 102 of 105 pairs, median 63.95 |
| dN/dS within families | 0.080 to 0.640, not estimable for *mcr-5* |
| Reference-based detection | 3 of 178 chromosomal proteins reached a family threshold |

## Layout

```
alignments/   FASTA, untrimmed and trimmed, for both datasets
trees/        ML trees, consensus trees, the constraint, and the topology test output
data/         metadata, catalytic site table, clade assignments, synteny, base composition
results/      every numeric result quoted in the manuscript, plus codeml runs
scripts/      the analysis, in the order it was run
```

`results/results.json` holds every number that appears in the manuscript, keyed by section.
`results/Table_S1_candidate_homologues.csv` is Supplementary Table S1.

## Running it

Python 3.11. Install the dependencies with

```bash
pip install -r requirements.txt
```

Three external binaries are needed and are not on PyPI:

- IQ-TREE 2.3.6, from <http://www.iqtree.org>
- PAML 4.10.10, for `codeml`, from <http://abacus.gene.ucl.ac.uk/software/paml.html>
- the ISEScan transposase profile library, from <https://github.com/xiezhq/ISEScan>

Scripts are numbered in run order within each folder. The primary analysis runs first,
then the revision scripts, then the figures, then the document builders. Paths at the top
of each script point at the directory layout used during the work and will need adjusting
if you move things around.

Two steps need network access that a locked-down environment may not grant: retrieval from
NCBI and the European Nucleotide Archive in `scripts/01_primary`, and the UniProtKB queries
in `scripts/02_revision`. The retrieved sequences are included here, so the rest of the
pipeline runs without them.

## Figures

| Figure | Script |
| --- | --- |
| 1, phylogeny | `scripts/03_figures/50_fig1.py` |
| 2, identity | `scripts/03_figures/51_fig2.py` |
| 3, domains and active site | `scripts/03_figures/53_fig3.py` |
| 4, genomic context | `scripts/03_figures/54_fig4.py` |
| 5, selection and base composition | `scripts/03_figures/52_fig5.py` |

## Data sources

All sequences are public records held by NCBI, the European Nucleotide Archive and
UniProtKB, and are identified by accession in `results/Table_S1_candidate_homologues.csv`
and in the manuscript. Nothing in this repository is original sequence data.

## Licence

MIT, see `LICENSE`. The sequence records redistributed here remain under the terms of the
databases they came from.
