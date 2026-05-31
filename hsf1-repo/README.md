# HSF1 cross-species domain analysis

Analysis code and data for the cross-species comparison of Heat Shock Factor 1
(HSF1) domain architecture in *Homo sapiens* (UniProt Q00613), *Caenorhabditis
elegans* (G5EFT5) and *Saccharomyces cerevisiae* (P10961), supporting the review
article *[YOUR TITLE HERE]* ([YOUR NAME], [YEAR]).

This repository regenerates the supplementary figures and the per-domain
quantitative analyses from the raw inputs.

## What's here

```
scripts/
  make_figureS1_coiledcoil.py     Figure S1 — multi-algorithm coiled-coil prediction
  make_figureS2_superposition.py  Figure S2 — DBD structural superposition + RMSD
data/
  AF_human.pdb, AF_worm.pdb, AF_yeast.pdb        AlphaFold DB models
  HSF1_*_{marcoil,multicoil2,ncoils,paircoil2}.txt   Waggawagga coiled-coil output (window 21)
  hsf1_msa.fasta                                  MAFFT alignment of the three sequences
  coils_*.csv                                     per-residue coiled-coil probabilities
  hsf1_perdomain_identity.xlsx                    per-domain pairwise identity
```

## Reproduce the figures

```
pip install -r requirements.txt
python scripts/make_figureS1_coiledcoil.py      # -> figureS1_coiledcoil.png
python scripts/make_figureS2_superposition.py   # -> figureS2_superposition.png (prints RMSD)
```

Run on Google Colab, a local Python 3.8+ environment, or any Jupyter setup.

## Data provenance

All primary analyses were run on public web services; this repository contains
the downstream parsing, calculation and plotting only.

| Input | Source |
|---|---|
| Sequences | UniProtKB: Q00613, G5EFT5, P10961 |
| Alignment | MAFFT, EMBL-EBI Job Dispatcher (job mafft-I20260523-213607-0461-26747537-p1m) |
| Coiled-coil | Waggawagga server (Marcoil, Multicoil2, nCOILS, Paircoil2; window 21) |
| Structures | AlphaFold Protein Structure Database |

## Domain coordinates used

DBD: human 15–120, worm 89–196, yeast 170–259.
HR-A/B: human 130–203, worm 206–253, yeast 350–403.
(See the article's Table 1 / Table S1 for the full set and evidence.)

## Methods summary

- **Per-domain identity:** computed from the MAFFT alignment over matched (non-gap)
  positions within each domain, human coordinates as reference frame.
- **Coiled-coil:** Waggawagga heptad-register output parsed to per-residue
  coiled-coil segments; window 21 primary.
- **Superposition:** Cα superposition (Kabsch via Biopython) over alignment-matched
  positions, restricted to the DBD and HR-A/B (mean pLDDT > 82); disordered regions
  excluded. RMSD in Å.

## Citation

If you use this code or data, please cite the article above. [Add Zenodo DOI here
once archived.]

## License

MIT (see LICENSE).
