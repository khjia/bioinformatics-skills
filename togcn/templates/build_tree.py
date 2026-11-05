"""Rebuild the BFS tree (with parents) of a TO-GCN run, so that the figure can show
TF inner layers and non-TF target genes attached to their parent TF."""
import sys
import numpy as np

run = sys.argv[1] if len(sys.argv) > 1 else "run95_all_tpm"

names = []
with open("data/all_tpm.tsv") as f:
    for line in f:
        names.append(line.split("\t", 1)[0])
n = len(names)
idx = {g: i for i, g in enumerate(names)}

eu, ev = [], []
with open(f"{run}/Node_relation.csv") as f:
    next(f)
    for line in f:
        a, b, _ = line.rstrip("\n").split(",")
        eu.append(idx[a]); ev.append(idx[b])
eu = np.asarray(eu, dtype=np.int32)
ev = np.asarray(ev, dtype=np.int32)
print("edges:", len(eu))

src = np.concatenate([eu, ev])
dst = np.concatenate([ev, eu])
order = np.argsort(src, kind="stable")
src_s = src[order]
adj = dst[order]
ptr = np.searchsorted(src_s, np.arange(n + 1), side="left")

level = np.full(n, -1, dtype=np.int16)
parent = np.full(n, -1, dtype=np.int32)
frontier = np.array([idx[l.strip()] for l in open("data/seeds.txt") if l.strip()], dtype=np.int32)
level[frontier] = 0
lv = 0
while frontier.size:
    lv += 1
    nxt = []
    for u in frontier.tolist():
        for p in range(ptr[u], ptr[u + 1]):
            v = adj[p]
            if level[v] < 0:
                level[v] = lv
                parent[v] = u
                nxt.append(v)
    frontier = np.array(nxt, dtype=np.int32)
print("reached:", int((level >= 0).sum()), "max level:", int(level.max()))

tf = set()
for line in open("/tmp/itak_pea/tf_classification.txt"):
    tf.add(line.split("\t")[0].split(".")[0])
istf = np.array([g in tf for g in names])

ok = level >= 0
with open("figs/tree_all_raw.tsv", "w") as f:
    f.write("gene\tlevel\tparent\tis_TF\tparent_level\tparent_is_TF\n")
    for i in np.nonzero(ok)[0].tolist():
        p = parent[i]
        f.write(f"{names[i]}\t{level[i]}\t{names[p] if p >= 0 else 'NA'}\t{int(istf[i])}\t"
                f"{level[p] if p >= 0 else -1}\t{int(istf[p]) if p >= 0 else 0}\n")

print("TFs with level:", int(istf[ok].sum()), "non-TF:", int((~istf[ok]).sum()))
