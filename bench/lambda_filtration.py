"""Connectivity filtration for inter-component crossing events.

For each event whose branches are inter-component at maxedge=60, find
lambda(e): the smallest quad-edge threshold at which the branches connect,
by adding four-valid quads in ascending max-edge order to a union-find and
checking pending events at checkpoints (bracket width recorded). Events
still disconnected with every four-valid quad are mask-disconnected.

    uv run python bench/lambda_filtration.py --json out/spectrum_full_d0.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
from windcheck import tifxyz                              # noqa: E402
from windcheck.check import load_pairs                    # noqa: E402
from windcheck.intrinsic import (SurfaceGraph, oriented_events,  # noqa: E402
                                 retained_quads)

CORPORA = {"Scroll 1": ("data/scroll1_tifxyz", "20230205180739", "out/crossing_s1"),
           "Scroll 5": ("data/scroll5_tifxyz", "20241024131839", "out/crossing"),
           "PHerc0139": ("data/PHerc0139_tifxyz", "20250728140407", "out/crossing_0139"),
           "PHerc0814": ("data/PHerc0814_tifxyz", "20250804134230", "out/crossing_0814"),
           "PHerc1667": ("data/PHerc1667_tifxyz", "20231117161658", "out/crossing_1667")}
CHECK_EVERY = 500


class DSU:
    def __init__(self, n):
        self.p = np.arange(n, dtype=np.int64)

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", type=Path, default=Path("out/spectrum_full_d0.json"))
    ap.add_argument("--out", type=Path, default=Path("out/lambda_filtration.json"))
    a = ap.parse_args()
    rows = json.loads(a.json.read_text())
    out = []
    for r in rows:
        if not r["events_inter_component"]:
            continue
        root, volume, work = CORPORA[r["corpus"]]
        d = Path(root) / r["segment"]
        s = tifxyz.read(sorted(d.glob(f"mesh/*{volume}*.tifxyz"))[0])
        rec = load_pairs(Path(work) / f"{r['segment'][:40]}_d0.csv")
        nv, nu = s.points.shape[:2]
        rec = rec[(rec["v1"] < nv - 1) & (rec["v2"] < nv - 1)
                  & (rec["u1"] < nu - 1) & (rec["u2"] < nu - 1)]
        g60 = SurfaceGraph(s.points, s.valid, 0)          # maxedge=60 reference
        evs = [e for e in oriented_events(rec) if not e["ambiguous"]]

        def comps(regs, graph):
            return {int(graph.comp[c]) for q in regs
                    for c in graph.quad_corners(*q)}

        inter = [e for e in evs if not (comps(e["region_a"], g60)
                                        & comps(e["region_b"], g60))]
        if not inter:
            continue
        # all four-valid quads, ranked by their own max edge
        P, V = s.points, s.valid
        Q4 = V[:-1, :-1] & V[1:, :-1] & V[:-1, 1:] & V[1:, 1:]
        p00, p10, p01, p11 = P[:-1, :-1], P[1:, :-1], P[:-1, 1:], P[1:, 1:]
        e = np.zeros(Q4.shape)
        for x, y in ((p00, p01), (p01, p11), (p11, p10), (p00, p10),
                     (p00, p11), (p10, p01)):
            e = np.maximum(e, np.linalg.norm(x - y, axis=-1))
        qv, qu = np.nonzero(Q4)
        order = np.argsort(e[qv, qu])
        qv, qu, qe = qv[order], qu[order], e[qv, qu][order]
        dsu = DSU(nv * nu)
        pend = [{"rows": len(ev["rows"]),
                 "ra": [v * nu + u for q in ev["region_a"]
                        for v, u in [q, (q[0]+1, q[1])]],
                 "rb": [v * nu + u for q in ev["region_b"]
                        for v, u in [q, (q[0]+1, q[1])]]} for ev in inter]
        lam = {}
        for i in range(len(qv)):
            base = qv[i] * nu + qu[i]
            for o in (1, nu, nu + 1):
                dsu.union(base, base + o)
            if i % CHECK_EVERY == 0 or i == len(qv) - 1:
                for j, ev in enumerate(pend):
                    if j in lam:
                        continue
                    if {dsu.find(x) for x in ev["ra"]} & \
                       {dsu.find(x) for x in ev["rb"]}:
                        lam[j] = float(qe[i])
                if len(lam) == len(pend):
                    break
        for j, ev in enumerate(pend):
            out.append({"corpus": r["corpus"], "segment": r["segment"],
                        "n_pairs": ev["rows"],
                        "lambda_vx": lam.get(j),          # None = mask-disc.
                        "voxel_um": r["voxel_um"]})
        got = [v for v in lam.values()]
        print(f"{r['corpus']:10s} {r['segment'][:34]:36s} inter {len(inter):4d} "
              f"mask-disc {len(pend)-len(lam):3d} "
              f"lambda med {np.median(got) if got else float('nan'):8.1f} vx",
              flush=True)
    a.out.write_text(json.dumps(out, indent=1))
    lams = [o["lambda_vx"] for o in out if o["lambda_vx"] is not None]
    print(f"\n{len(out)} inter events: {len(out)-len(lams)} mask-disconnected, "
          f"{len(lams)} connect; lambda median {np.median(lams):.0f} vx, "
          f"p90 {np.percentile(lams, 90):.0f}, max {max(lams):.0f}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
