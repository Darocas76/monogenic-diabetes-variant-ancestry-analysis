import json, collections
gnv=json.load(open("gnomad_v4_all.json"))
cl=json.load(open("clinvar_lookup.json"))
EUR={"nfe","fin","asj"}; NONEUR={"afr","amr","eas","sas","mid"}
def grp(pops):
    eur=any(pops.get(p,{}).get("ac",0)>0 for p in EUR)
    non=any(pops.get(p,{}).get("ac",0)>0 for p in NONEUR)
    if eur and non: return "shared"
    if eur: return "EUR-private"
    if non: return "nonEUR-private"
    return "unassigned"
# build merged records
recs=[]
for g,rows in gnv.items():
    for r in rows:
        key=r["variant_id"]
        recs.append({"gene":g,"key":key,"cons":r["consequence"],
                     "cat":cl.get(key),            # None if not in ClinVar
                     "grp":grp(r["pops"])})
def analyze(records, label):
    tot=len(records)
    ann=[x for x in records if x["cat"] is not None]
    gap=[x for x in records if x["cat"] is None]
    classified=[x for x in ann if x["cat"] in ("VUS","PLP","BLB")]
    print(f"\n===== {label}  (genes: {len(set(x['gene'] for x in records))}) =====")
    print(f"gnomAD variants: {tot} | ClinVar-annotated: {len(ann)} ({100*len(ann)/tot:.1f}%) | annotation gap: {len(gap)} ({100*len(gap)/tot:.1f}%)")
    print(f"classified (VUS/PLP/BLB): {len(classified)}")
    # (2) annotation gap stratified by ancestry
    gd=collections.Counter(x["grp"] for x in gap)
    print("ANNOTATION GAP by ancestry:")
    asg=sum(v for k,v in gd.items() if k!='unassigned')
    for k in ("EUR-private","nonEUR-private","shared","unassigned"):
        v=gd.get(k,0); pe=f"{100*v/asg:.1f}% of assigned" if k!='unassigned' and asg else ""
        print(f"   {k:<16} {v:>6}  {pe}")
    # (1) population-private classified comparison
    print("POPULATION-PRIVATE classified comparison (VUS / PLP / BLB):")
    for grpname in ("EUR-private","nonEUR-private","shared"):
        sub=[x for x in classified if x["grp"]==grpname]
        n=len(sub)
        if not n: continue
        c=collections.Counter(x["cat"] for x in sub)
        print(f"   {grpname:<16} n={n:<6} VUS={100*c['VUS']/n:4.1f}%  PLP={100*c['PLP']/n:4.1f}%  BLB={100*c['BLB']/n:4.1f}%")
    return tot,len(ann),len(gap)

analyze(recs, "ALL 17 GENES")
analyze([x for x in recs if x["gene"]!="KLF11"], "PRIMARY: 16 GENES (KLF11 excluded)")
# save merged for later
json.dump(recs, open("merged_v4.json","w"))
print("\nsaved merged_v4.json")
