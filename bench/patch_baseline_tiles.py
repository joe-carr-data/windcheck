"""Apples to apples: cut the corpus into patch-sized tiles and count.

A crossing is inside a tile only when BOTH quads land in the SAME tile,
which is exactly the visibility condition for a patch. Tiling handles
clustering automatically -- a defect that produces thousands of pairs
still only spoils the tiles it actually occupies.
"""
import csv, json
from pathlib import Path
import numpy as np, tifffile

W = 24
CERTS = Path("results/certificates")   # override with --certs
DIRS = ["out/crossing", "out/crossing_0139", "out/crossing_0814",
        "out/crossing_1667", "out/crossing_s1"]

tiles_tot = tiles_bad = 0
for cert in sorted(CERTS.glob("*_certificate.json")):
    stem = cert.name[:-len("_certificate.json")]
    mesh = Path(json.loads(cert.read_text())["mesh"]["path"])
    csvs = [p for d in DIRS for p in Path(d).glob(f"{stem}_d*.csv")]
    if not mesh.is_dir() or not csvs:
        continue
    x = tifffile.imread(mesh / "x.tif"); y = tifffile.imread(mesh / "y.tif")
    z = tifffile.imread(mesh / "z.tif")
    valid = ~((x == -1) & (y == -1) & (z == -1))
    valid &= np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    m = mesh / "mask.tif"
    if m.exists():
        mk = np.asarray(tifffile.imread(m))
        if mk.shape == valid.shape:
            valid &= mk.astype(bool)

    V, U = valid.shape
    nV, nU = (V + W - 1) // W, (U + W - 1) // W
    # a tile counts only if it holds enough valid cells to be a real patch
    occ = np.zeros((nV, nU), np.int64)
    for i in range(nV):
        for j in range(nU):
            occ[i, j] = valid[i*W:(i+1)*W, j*W:(j+1)*W].sum()
    live = occ >= 164          # the smallest patch in the pilot sample
    bad = np.zeros((nV, nU), bool)
    for c in csvs:
        with c.open() as fh:
            for r in csv.DictReader(fh):
                if r["verdict"] != "transverse":
                    continue
                v1, u1, v2, u2 = (int(r["v1"]), int(r["u1"]),
                                  int(r["v2"]), int(r["u2"]))
                if v1 // W == v2 // W and u1 // W == u2 // W:
                    bad[v1 // W, u1 // W] = True
    tiles_tot += int(live.sum())
    tiles_bad += int((live & bad).sum())

print(f"patch-sized tiles ({W}x{W}) in the corpus with >=164 valid cells: {tiles_tot:,}")
print(f"tiles containing a self-intersection:                            {tiles_bad:,}")
p = tiles_bad / tiles_tot
print(f"fraction:                                                        {100*p:.4f}%")
print()
n = 500
print(f"if 500 verified patches behaved like corpus tiles, expected with crossings: {p*n:.2f}")
print(f"probability all 500 come out clean at that rate:  {(1-p)**n:.3e}")
json.dump({"window": W, "tiles": tiles_tot, "tiles_with_crossing": tiles_bad,
           "fraction": p}, open("out/patches/corpus_tile_rate.json", "w"), indent=1)
