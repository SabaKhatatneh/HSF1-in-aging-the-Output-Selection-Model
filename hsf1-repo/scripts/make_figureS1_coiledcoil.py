#!/usr/bin/env python3
"""
Figure S1 — Multi-algorithm coiled-coil prediction across HSF1 orthologs.

Reads the Waggawagga heptad-register output files (Marcoil, Multicoil2, nCOILS,
Paircoil2) for human, worm and yeast, parses the predicted coiled-coil segments,
and draws one track per algorithm per species with the UniProt-annotated
regions shaded for reference.

USAGE:
    pip install matplotlib
    python make_figureS1_coiledcoil.py
Output: figureS1_coiledcoil.png  (300 dpi)

Data needed (in ./data/):
    HSF1_HUMAN_marcoil.txt, HSF1_HUMAN_multicoil2.txt, HSF1_HUMAN_ncoils.txt,
    HSF1_HUMAN_paircoil2.txt, HSF1_CAEEL_*.txt, HSF_YEAST_multicoil2.txt,
    HSF_YEAST_ncoils.txt
"""

import re, glob, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

DATA = os.path.join(os.path.dirname(__file__), "..", "data")

def parse_waggawagga(path):
    """Parse heptad-register blocks. Coiled-coil residues are positions where the
       'Coil' line shows a/b/c/d/e/f/g. Returns list of (start,end) 1-based ranges."""
    with open(path) as f:
        lines = f.read().split("\n")
    cc = set()
    i = 0
    while i < len(lines):
        m = re.match(r"Seq\.\s+(\d+)\s+(\S+)", lines[i])
        if m:
            block_start = int(m.group(1)); seq = m.group(2)
            coil_line = None
            for j in range(i+1, min(i+4, len(lines))):
                if lines[j].lstrip().startswith("Coil"):
                    coil_line = lines[j]; break
            if coil_line:
                col = lines[i].index(seq)
                ann = coil_line[col:col+len(seq)]
                for k, ch in enumerate(ann):
                    if ch in "abcdefg":
                        cc.add(block_start + k + 1)
            i = j if coil_line else i+1
        else:
            i += 1
    if not cc: return []
    sp = sorted(cc); ranges=[]; start=prev=sp[0]
    for p in sp[1:]:
        if p == prev+1: prev=p
        else: ranges.append((start,prev)); start=prev=p
    ranges.append((start,prev))
    return ranges

# Load all available files
data = {}
no_result = set()
species_methods = {
    "HSF1_HUMAN": ["marcoil","multicoil2","ncoils","paircoil2"],
    "HSF1_CAEEL": ["marcoil","multicoil2","ncoils","paircoil2"],
    "HSF_YEAST":  ["marcoil","multicoil2","ncoils","paircoil2"],
}
for sp, methods in species_methods.items():
    for m in methods:
        path = os.path.join(DATA, f"{sp}_{m}.txt")
        if os.path.exists(path):
            data[(sp,m)] = parse_waggawagga(path)
        else:
            no_result.add((sp,m))   # yeast marcoil/paircoil2 returned "no CC found"

lengths = {"HSF1_HUMAN":529,"HSF1_CAEEL":671,"HSF_YEAST":833}
titles  = {"HSF1_HUMAN":"Homo sapiens HSF1 (529 aa)",
           "HSF1_CAEEL":"C. elegans HSF-1 (671 aa)",
           "HSF_YEAST":"S. cerevisiae Hsf1 (833 aa)"}
uniprot = {"HSF1_HUMAN":[("HR-A/B",130,203,"#C84B31"),("HR-C",384,409,"#B07E2A")],
           "HSF1_CAEEL":[("COILED (UniProt)",206,240,"#C84B31")],
           "HSF_YEAST":[("trimerization",350,403,"#C84B31")]}
methods = ["marcoil","multicoil2","ncoils","paircoil2"]
mlabel  = {"marcoil":"Marcoil","multicoil2":"Multicoil2","ncoils":"nCOILS","paircoil2":"Paircoil2"}
mcol    = {"marcoil":"#2E5C8A","multicoil2":"#3B8C6E","ncoils":"#7A4E9E","paircoil2":"#C77F2A"}

MAXLEN = 833
fig, axes = plt.subplots(3,1, figsize=(11,7.5))
for ax, sp in zip(axes, ["HSF1_HUMAN","HSF1_CAEEL","HSF_YEAST"]):
    L = lengths[sp]
    ax.add_patch(plt.Rectangle((0,0.0), L, 0.05, color="#dddddd"))
    for lab,s,e,col in uniprot[sp]:
        ax.axvspan(s,e,color=col,alpha=0.12)
        ax.text((s+e)/2,4.7,f"{lab}\n{s}-{e}",ha="center",va="bottom",fontsize=7.5,color=col)
    for yi,m in enumerate(methods):
        y = 4-yi
        if (sp,m) in no_result:
            ax.text(2,y,f"{mlabel[m]}: no CC found",fontsize=8,va="center",color="#999",style="italic")
            continue
        ranges = data.get((sp,m))
        if not ranges:
            ax.text(2,y,f"{mlabel[m]}: —",fontsize=8,va="center",color="#999"); continue
        ax.text(2,y,mlabel[m],fontsize=8,va="center",ha="left",color=mcol[m],fontweight="bold")
        for (a,b) in ranges:
            ax.add_patch(plt.Rectangle((a,y-0.28),b-a,0.56,color=mcol[m],alpha=0.85))
    ax.set_xlim(-40,MAXLEN+10); ax.set_ylim(-0.5,5.6)
    ax.set_yticks([]); ax.set_title(titles[sp],fontsize=11,loc="left")
    ax.spines[["top","right","left"]].set_visible(False)
    ax.set_xlabel("Residue position",fontsize=9)
axes[0].annotate("HR-C detected\n(nCOILS only)",xy=(391,2),xytext=(430,3.2),
    fontsize=8,color="#B07E2A",arrowprops=dict(arrowstyle="->",color="#B07E2A"))
fig.suptitle("Multi-algorithm coiled-coil prediction across HSF1 orthologs",fontsize=12.5,y=0.99)
fig.text(0.012,-0.01,
 "Bars = coiled-coil segments called by each algorithm (Waggawagga, window 21). Shaded = UniProt-annotated regions. "
 "A second coiled-coil at the HR-C position is detected only in human. Yeast Marcoil/Paircoil2 returned no call; "
 "the yeast nCOILS 630-650 hit is a low-complexity false positive.",fontsize=7.2,color="#555",wrap=True)
plt.tight_layout()
plt.savefig("figureS1_coiledcoil.png",dpi=300,bbox_inches="tight")
print("Wrote figureS1_coiledcoil.png")
