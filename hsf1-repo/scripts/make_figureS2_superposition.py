#!/usr/bin/env python3
"""
Figure S2 — Structural superposition of the HSF1 DNA-binding domain.

Loads the three AlphaFold models, extracts the DBD of each, superposes worm and
yeast DBD onto the human DBD over alignment-matched Cα positions (Kabsch via
Biopython), reports pairwise RMSD, and draws the three Cα traces overlaid.

USAGE:
    pip install biopython matplotlib numpy
    python make_figureS2_superposition.py
Output: figureS2_superposition.png (300 dpi) + prints RMSD values

Data needed (in ./data/):
    AF_human.pdb, AF_worm.pdb, AF_yeast.pdb   (AlphaFold DB PDB files)
    hsf1_msa.fasta                            (MAFFT alignment, 3 sequences)

NOTE: This is a Cα-trace figure. For publication-quality cartoon ribbons,
open the same PDBs in PyMOL/ChimeraX (commands in README). The RMSD values
are identical either way — they come from this script.
"""

import warnings; warnings.simplefilter("ignore")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser, Superimposer
from Bio import SeqIO

import os
DATA = os.path.join(os.path.dirname(__file__), "..", "data")

# Domain boundaries (residue numbers in each native sequence)
DBD  = {"human":(15,120),  "worm":(89,196),  "yeast":(170,259)}
HRAB = {"human":(130,203), "worm":(206,253), "yeast":(350,403)}
files = {"human":f"{DATA}/AF_human.pdb","worm":f"{DATA}/AF_worm.pdb","yeast":f"{DATA}/AF_yeast.pdb"}
namemap = {"human":"HSF1_HUMAN","worm":"HSF1_CAEEL","yeast":"HSF_YEAST"}

p = PDBParser(QUIET=True)
structs = {n: p.get_structure(n,f) for n,f in files.items()}
aln = {r.id.split()[0]: str(r.seq) for r in SeqIO.parse(f"{DATA}/hsf1_msa.fasta","fasta")}

def resnum_to_col(a):
    m={}; rn=0
    for col,ch in enumerate(a):
        if ch!='-': rn+=1; m[rn]=col
    return m
def col_map(a):
    m={}; rn=0
    for col,ch in enumerate(a):
        if ch!='-': rn+=1; m[col]=rn
    return m

def ca_atoms(struct,s,e):
    d={}
    for r in struct.get_residues():
        if r.id[0]==' ' and s<=r.id[1]<=e:
            for at in r:
                if at.get_name()=='CA': d[r.id[1]]=at
    return d
def ca_coords(struct,s,e):
    d={}
    for r in struct.get_residues():
        if r.id[0]==' ' and s<=r.id[1]<=e:
            for at in r:
                if at.get_name()=='CA': d[r.id[1]]=at.get_coord()
    return d

def matched_pairs(spA,spB,domA,domB):
    aA=aln[namemap[spA]]; aB=aln[namemap[spB]]
    mA=col_map(aA); mB=col_map(aB); pairs=[]
    sA,eA=domA; sB,eB=domB
    for col in range(len(aA)):
        if aA[col]!='-' and aB[col]!='-':
            ra=mA[col]; rb=mB[col]
            if sA<=ra<=eA and sB<=rb<=eB: pairs.append((ra,rb))
    return pairs

def rmsd(spA,spB,domdict):
    caA=ca_atoms(structs[spA],*domdict[spA]); caB=ca_atoms(structs[spB],*domdict[spB])
    fixed=[];moving=[]
    for ra,rb in matched_pairs(spA,spB,domdict[spA],domdict[spB]):
        if ra in caA and rb in caB: fixed.append(caA[ra]); moving.append(caB[rb])
    si=Superimposer(); si.set_atoms(fixed,moving)
    return si.rms, len(fixed), si.rotran

print("=== DBD RMSD (Cα, alignment-matched) ===")
for a,b in [("human","worm"),("human","yeast"),("worm","yeast")]:
    r,n,_=rmsd(a,b,DBD); print(f"  {a}-{b}: {r:.2f} Å over {n} atoms")
print("=== HR-A/B RMSD ===")
for a,b in [("human","worm"),("human","yeast"),("worm","yeast")]:
    r,n,_=rmsd(a,b,HRAB); print(f"  {a}-{b}: {r:.2f} Å over {n} atoms")

# ---- 3D figure: DBD traces superposed on human ----
fig=plt.figure(figsize=(9,7)); ax=fig.add_subplot(111,projection='3d')
colors={"human":"#2E5C8A","worm":"#C84B31","yeast":"#7A9E50"}
hc=ca_coords(structs["human"],*DBD["human"]); hk=sorted(hc)
arr=np.array([hc[k] for k in hk])
ax.plot(arr[:,0],arr[:,1],arr[:,2],color=colors["human"],lw=2,label="Human DBD (ref)")
for sp in ["worm","yeast"]:
    r,n,(rot,tran)=rmsd("human",sp,DBD)
    bc=ca_coords(structs[sp],*DBD[sp]); bk=sorted(bc)
    barr=np.array([bc[k] for k in bk]); barr_t=np.dot(barr,rot)+tran
    ax.plot(barr_t[:,0],barr_t[:,1],barr_t[:,2],color=colors[sp],lw=2,
            label=f"{sp.capitalize()} DBD (RMSD {r:.2f} Å)")
ax.set_title("DBD structural superposition — HSF1 orthologs\n(AlphaFold models, Cα trace)",fontsize=12)
ax.legend(fontsize=9,loc="upper left")
ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)")
ax.view_init(elev=18,azim=45)
plt.tight_layout()
plt.savefig("figureS2_superposition.png",dpi=300,bbox_inches="tight")
print("Wrote figureS2_superposition.png")
