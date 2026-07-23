import json, time, os, sys
import requests
API="https://gnomad.broadinstitute.org/api"
DS="gnomad_r4"
OUT="gnomad_v4"
GENES=["HNF1A","HNF4A","HNF1B","GCK","KCNJ11","ABCC8","INS","PDX1","NEUROD1","PTF1A","CEL","PPARG","APPL1","BLK","KLF11","PAX4","WFS1"]
POPS={"afr","amr","eas","sas","nfe","fin","asj","mid"}
Q='query($g:String!,$ds:DatasetId!){gene(gene_symbol:$g,reference_genome:GRCh38){variants(dataset:$ds){variant_id consequence exome{ac an populations{id ac an}} genome{ac an populations{id ac an}}}}}'
def fetch(gene):
    r=requests.post(API,json={"query":Q,"variables":{"g":gene,"ds":DS}},
                    headers={"Content-Type":"application/json"},timeout=40)
    if r.status_code==429:
        time.sleep(12); r=requests.post(API,json={"query":Q,"variables":{"g":gene,"ds":DS}},headers={"Content-Type":"application/json"},timeout=40)
    r.raise_for_status(); d=r.json()
    if "errors" in d: raise RuntimeError(str(d["errors"])[:300])
    return d["data"]["gene"]["variants"]
def combine(v):
    pops={}
    for src in ("exome","genome"):
        blk=v.get(src) or {}
        for p in (blk.get("populations") or []):
            pid=p["id"].lower()
            if pid not in POPS: continue
            a=pops.setdefault(pid,{"ac":0,"an":0})
            a["ac"]+=p.get("ac",0) or 0; a["an"]+=p.get("an",0) or 0
    gac=((v.get("exome") or {}).get("ac",0) or 0)+((v.get("genome") or {}).get("ac",0) or 0)
    gan=((v.get("exome") or {}).get("an",0) or 0)+((v.get("genome") or {}).get("an",0) or 0)
    return {"variant_id":v["variant_id"],"consequence":v.get("consequence",""),"pops":pops,"ac":gac,"an":gan}
start=time.time()
done=[]
for g in GENES:
    fp=f"{OUT}/{g}.json"
    if os.path.exists(fp): continue
    if time.time()-start>38:
        print("time guard"); break
    try:
        vs=fetch(g); rows=[combine(v) for v in vs]
        json.dump(rows,open(fp,"w"))
        print(f"{g}: {len(rows)} variants saved"); done.append(g)
    except Exception as e:
        print(f"{g} FAILED: {e}")
    time.sleep(4)
have=[g for g in GENES if os.path.exists(f"{OUT}/{g}.json")]
print(f"have {len(have)}/17: {have}")
