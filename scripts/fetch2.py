import sys, json, subprocess, os, pandas as pd
S = sys.argv[1]
raw, red = f"{S}/data", f"{S}/reduced2"
api = "https://repo-prod.prod.sagebase.org/repo/v1"
kids = json.loads(subprocess.check_output(["curl","-sS","-X","POST",f"{api}/entity/children",
  "-H","Content-Type: application/json","-d",'{"parentId":"syn20613593","includeTypes":["file"]}']))["page"]
keep = ["treatment","cell_line","time","fileID","p.ERK","p.MEK","p.p90RSK","p.Akt.Ser473.","p.AKT.Thr308.",
        "p.S6","p.S6K","p.4EBP1","cleavedCas","Ki.67","p.RB","GAPDH",
        "p.p38","p.MAPKAPK2","p.MKK3.MKK6","p.GSK3b","b.CATENIN","p.JNK","p.CREB"]
for k in kids:
    name = k["name"].replace(".csv","")
    out = f"{red}/{name}.parquet"
    if os.path.exists(out): continue
    f = f"{raw}/{k['name']}"
    for attempt in range(5):
        if os.path.exists(f) and os.path.getsize(f) > 1_000_000: break
        u = subprocess.check_output(["curl","-sS","-o","/dev/null","-w","%{redirect_url}",f"{api}/entity/{k['id']}/file"]).decode()
        subprocess.call(["curl","-sS","--retry","4","-o",f,u])
    df = pd.read_csv(f, usecols=keep, skipinitialspace=True)
    for c in keep[4:]: df[c] = df[c].astype("float32")
    df.to_parquet(out); os.remove(f)
    print(name, len(df), flush=True)
print("DONE")
