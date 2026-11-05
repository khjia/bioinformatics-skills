"""Merge TO-GCN layers under BOTH a lower and an upper size bound.

Rules
  (1) merge: any group with < min_n genes is merged into the immediate neighbour
      whose mean temporal profile is more similar (Pearson r);
  (2) no group may exceed max_n genes after a merge; if a merge would break the
      cap the other neighbour is used (if both break it, the smaller result wins);
  (3) a native layer that is ALREADY larger than max_n is split into
      ceil(n / max_n) consecutive parts ordered by the gene-level centre-of-mass
      time, so the temporal ordering inside the layer is preserved.

Layer order is never permuted.
"""
import argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--run", default="run95_all_tpm")
ap.add_argument("--tag", default="all_raw")
ap.add_argument("--min_n", type=int, default=100)
ap.add_argument("--max_n", type=int, default=2000)
ap.add_argument("--split_large", action="store_true",
                help="split native layers bigger than max_n (default: keep them intact)")
ap.add_argument("--verbose", action="store_true")
a = ap.parse_args()

times = np.array([0, 30, 60, 90, 120])
rows = []
for line in open("data/all_tpm.tsv"):
    p = line.rstrip("\n").split("\t")
    rows.append((p[0], np.log2(np.array([float(x) for x in p[1:]]) + 1)))
names = [r[0] for r in rows]
M = np.vstack([r[1] for r in rows])
t = np.repeat(times, 3)
mu = np.column_stack([M[:, t == tt].mean(1) for tt in times])
gz = (mu - mu.mean(1, keepdims=True)) / mu.std(1, keepdims=True)
idx = {g: i for i, g in enumerate(names)}

lay = {}
for i, line in enumerate(open(f"{a.run}/Node_level.csv")):
    if i == 0:
        continue
    g, l = line.strip().split(",")
    lay[g] = int(l)


def prof_of(genes):
    return gz[[idx[g] for g in genes]].mean(0)


def com(genes):
    X = gz[[idx[g] for g in genes]]
    w = np.clip(X, 0, None)
    s = w.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(s > 0, (w * times).sum(1) / np.where(s > 0, s, 1), times[X.argmax(1)])
    return c


groups = []
for l in sorted(set(lay.values())):
    gs = [g for g, v in lay.items() if v == l and g in idx]
    if gs:
        groups.append({"disp": f"L{l}", "layers": [l], "genes": gs, "prof": prof_of(gs)})

# ---------- (3) optional split of native layers larger than max_n ----------
# OFF by default: big native layers are kept intact, max_n only limits merges.
new = []
for g in groups:
    n = len(g["genes"])
    k = int(np.ceil(n / a.max_n)) if a.split_large else 1
    if k <= 1:
        new.append(g)
        continue
    c = com(g["genes"])
    order = np.argsort(c, kind="stable")
    base = g["disp"]
    for j in range(k):
        sl = [g["genes"][x] for x in order[j * n // k:(j + 1) * n // k]]
        new.append({"disp": f"{base}{'abcdefgh'[j]}", "layers": g["layers"],
                    "genes": sl, "prof": prof_of(sl)})
    if a.verbose:
        print(f"[split] {base} (n={n}) -> {k} parts of ~{n // k}")
groups = new

# ---------- (1)+(2) merge small groups, never above the cap ----------
while True:
    small = [(len(g["genes"]), i) for i, g in enumerate(groups) if len(g["genes"]) < a.min_n]
    if not small or len(groups) == 1:
        break
    n_i, i = min(small)
    cands = []
    if i > 0:
        rl = np.corrcoef(groups[i]["prof"], groups[i - 1]["prof"])[0, 1]
        cands.append((rl, i - 1, len(groups[i - 1]["genes"]) + n_i))
    if i < len(groups) - 1:
        rr = np.corrcoef(groups[i]["prof"], groups[i + 1]["prof"])[0, 1]
        cands.append((rr, i, len(groups[i + 1]["genes"]) + n_i))
    ok = [c for c in cands if c[2] <= a.max_n]
    _, j, cap_n = max(ok, key=lambda c: c[0]) if ok else min(cands, key=lambda c: c[2])
    A, B = groups[j], groups[j + 1]
    genes = A["genes"] + B["genes"]
    layers = sorted(set(A["layers"] + B["layers"]))
    disp = f"L{layers[0]}" if len(layers) == 1 else f"L{layers[0]}-{layers[-1]}"
    if a.verbose:
        print(f"[merge] {groups[i]['disp']} (n={n_i}) -> {A['disp']}"
              f"   result n={cap_n}{'   (cap broken)' if cap_n > a.max_n else ''}")
    groups = groups[:j] + [{"disp": disp, "layers": layers, "genes": genes,
                            "prof": prof_of(genes)}] + groups[j + 2:]

print(f"\n== {a.tag}: {len(groups)} groups (min_n={a.min_n}, max_n={a.max_n}) ==")
print(f"{'group':>9s} {'n':>6s}  profile(z)                  peak  adj_r")
out = []
for i, g in enumerate(groups):
    r = f"{np.corrcoef(g['prof'], groups[i+1]['prof'])[0,1]:.3f}" if i < len(groups) - 1 else "-"
    print(f"{g['disp']:>9s} {len(g['genes']):6d}  " + " ".join(f"{x:5.2f}" for x in g["prof"]) +
          f"  {times[int(np.argmax(g['prof']))]:4d}  {r}")
    for x in g["genes"]:
        out.append((x, g["disp"], i, lay[x]))

out.sort(key=lambda r: (r[2], r[3], r[0]))
with open(f"figs/groups_{a.tag}.csv", "w") as f:
    f.write("gene,group,group_order,layer\n")
    for row in out:
        f.write(",".join(map(str, row)) + "\n")
print("written figs/groups_%s.csv" % a.tag)
