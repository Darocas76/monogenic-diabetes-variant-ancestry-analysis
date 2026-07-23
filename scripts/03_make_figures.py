import json, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
recs=[x for x in json.load(open("merged_v4.json")) if x["gene"]!="KLF11"]  # primary 16
EUR={"nfe","fin","asj"}; NONEUR={"afr","amr","eas","sas","mid"}
def pct(a,b): return 100*a/b if b else 0.0

# overall population-private classified
def grpstats(grp):
    cl=[x for x in recs if x["grp"]==grp and x["cat"] in ("VUS","PLP","BLB")]
    n=len(cl); c=collections.Counter(x["cat"] for x in cl)
    return n, pct(c["VUS"],n), pct(c["PLP"],n), pct(c["BLB"],n)
ne_n,ne_v,ne_p,ne_b = grpstats("nonEUR-private")
e_n,e_v,e_p,e_b = grpstats("EUR-private")

# Fig 1 — headline: classification of population-private variants by ancestry
fig,ax=plt.subplots(figsize=(7,4.5))
cats=["VUS","P/LP","B/LB"]; x=np.arange(3); w=0.38
eur=[e_v,e_p,e_b]; non=[ne_v,ne_p,ne_b]
b1=ax.bar(x-w/2,eur,w,label=f"European-private (n={e_n})",color="#d9772b")
b2=ax.bar(x+w/2,non,w,label=f"Non-European-private (n={ne_n})",color="#3a78b5")
ax.set_xticks(x); ax.set_xticklabels(cats); ax.set_ylabel("% of classified variants")
ax.set_title("Classification of population-private variants by ancestry\n(16 monogenic diabetes genes, gnomAD v4, ClinVar 2026-06-29)",fontsize=10)
for bars in (b1,b2):
    for r in bars: ax.text(r.get_x()+r.get_width()/2,r.get_height()+0.6,f"{r.get_height():.1f}%",ha="center",fontsize=8)
ax.legend(fontsize=8); ax.set_ylim(0,60); fig.tight_layout(); fig.savefig("Figure_2.png",dpi=300); plt.close(fig)

# Fig 2 — annotation gap near-universal (gap-rate by ancestry)
by=collections.defaultdict(lambda:{"tot":0,"ann":0})
for x in recs:
    if x["grp"]=="unassigned": continue
    by[x["grp"]]["tot"]+=1; by[x["grp"]]["ann"]+= x["cat"] is not None
tot_all=len(recs); ann_all=sum(1 for x in recs if x["cat"] is not None)
labels=["EUR-private","nonEUR-private","shared","Overall"]
vals=[pct(by["EUR-private"]["tot"]-by["EUR-private"]["ann"],by["EUR-private"]["tot"]),
      pct(by["nonEUR-private"]["tot"]-by["nonEUR-private"]["ann"],by["nonEUR-private"]["tot"]),
      pct(by["shared"]["tot"]-by["shared"]["ann"],by["shared"]["tot"]),
      pct(tot_all-ann_all,tot_all)]
fig,ax=plt.subplots(figsize=(6.5,4.2))
cols=["#3a78b5","#d9772b","#6aa84f","#777777"]
b=ax.bar(labels,vals,color=cols)
for r in b: ax.text(r.get_x()+r.get_width()/2,r.get_height()+0.8,f"{r.get_height():.1f}%",ha="center",fontsize=9)
ax.set_ylabel("% gnomAD variants without ClinVar classification"); ax.set_ylim(0,100)
ax.set_title("Annotation gap is near-universal across ancestry\n(not disproportionately non-European)",fontsize=10)
fig.tight_layout(); fig.savefig("Figure_1.png",dpi=300); plt.close(fig)

# Fig 3 — per-gene heatmap of P/LP% by ancestry
GENES=["WFS1","KCNJ11","NEUROD1","HNF1A","ABCC8","GCK","HNF1B","PDX1","INS","PPARG","HNF4A","CEL","BLK","PAX4","PTF1A","APPL1"]
M=[]; ann=[]
bg=collections.defaultdict(list)
for x in recs: bg[x["gene"]].append(x)
for g in GENES:
    row=[]; arow=[]
    for grp in ("EUR-private","nonEUR-private"):
        cl=[x for x in bg[g] if x["grp"]==grp and x["cat"] in ("VUS","PLP","BLB")]
        n=len(cl); p=pct(sum(1 for x in cl if x["cat"]=="PLP"),n)
        row.append(p); arow.append(f"{p:.0f}%\n(n={n})")
    M.append(row); ann.append(arow)
M=np.array(M)
fig,ax=plt.subplots(figsize=(5.2,7))
im=ax.imshow(M,cmap="YlOrRd",vmin=0,vmax=70,aspect="auto")
ax.set_xticks([0,1]); ax.set_xticklabels(["European-\nprivate","Non-European-\nprivate"],fontsize=9)
ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES,fontsize=8,style="italic")
for i in range(len(GENES)):
    for j in range(2): ax.text(j,i,ann[i][j],ha="center",va="center",fontsize=6.5)
ax.set_title("P/LP rate among classified\npopulation-private variants (%)",fontsize=10)
fig.colorbar(im,ax=ax,shrink=0.6,label="P/LP %"); fig.tight_layout(); fig.savefig("Figure_3.png",dpi=300); plt.close(fig)
print("figures written:", "Figure_2.png Figure_1.png Figure_3.png")
