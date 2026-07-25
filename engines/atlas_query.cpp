// atlas_query — nearest-surface winding lookup over a reference atlas.
//
// Given a set of reference surfaces, each tagged with a winding number, and a
// set of query points, report for every query point:
//
//     w1, d1   winding of the nearest reference surface, and the distance to it
//     w2, d2   the same for the nearest surface of a *different* winding
//
// d1 is a true point-to-surface distance (nearest point on the bilinear quad,
// via its two triangles), not a distance to the nearest grid sample. That
// matters: the reference grids are sampled every ~20 voxels while adjacent
// sheets sit ~30 voxels apart, so a nearest-sample lookup has only a ~1.5x
// margin and cannot certify anything. Nearest-point drives d1 to ~0 on the
// correct sheet while d2 stays at the sheet spacing.
//
// The ratio d1/d2 is the confidence: an assignment is only used when the point
// is unambiguously on one sheet. A violation reported at low confidence is not
// a violation, it is a coin flip.
//
// Acceleration: uniform grid over triangle AABBs, searched in expanding shells
// so the walk stops as soon as no unvisited shell can beat the current best.
//
// Build:
//   clang++ -O3 -std=c++17 -pthread -o atlas_query engines/atlas_query.cpp

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

struct Vec3 {
    float x = 0, y = 0, z = 0;
};

inline Vec3 sub(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// Squared distance from p to triangle (a,b,c). Ericson, Real-Time Collision
// Detection §5.1.5: classify against the vertex/edge Voronoi regions, then
// fall through to the interior barycentric case.
float dist2_point_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 ab = sub(b, a), ac = sub(c, a), ap = sub(p, a);
    const float d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0.0f && d2 <= 0.0f) return dot(ap, ap);

    const Vec3 bp = sub(p, b);
    const float d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0.0f && d4 <= d3) return dot(bp, bp);

    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        const float v = d1 / (d1 - d3);
        const Vec3 q{a.x + v * ab.x, a.y + v * ab.y, a.z + v * ab.z};
        const Vec3 pq = sub(p, q);
        return dot(pq, pq);
    }

    const Vec3 cp = sub(p, c);
    const float d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0.0f && d5 <= d6) return dot(cp, cp);

    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        const float w = d2 / (d2 - d6);
        const Vec3 q{a.x + w * ac.x, a.y + w * ac.y, a.z + w * ac.z};
        const Vec3 pq = sub(p, q);
        return dot(pq, pq);
    }

    const float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
        const float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const Vec3 q{b.x + w * (c.x - b.x), b.y + w * (c.y - b.y), b.z + w * (c.z - b.z)};
        const Vec3 pq = sub(p, q);
        return dot(pq, pq);
    }

    const float denom = 1.0f / (va + vb + vc);
    const float v = vb * denom, w = vc * denom;
    const Vec3 q{a.x + ab.x * v + ac.x * w, a.y + ab.y * v + ac.y * w, a.z + ab.z * v + ac.z * w};
    const Vec3 pq = sub(p, q);
    return dot(pq, pq);
}

struct Tri {
    Vec3 a, b, c;
    int32_t winding;  // group tag: a winding number, or (self-gap mode) a u-index
};

// Self-gap mode. Instead of "which reference sheet is nearest", ask "how far is
// the nearest OTHER PART OF THIS SAME TRACE". Tag every triangle with its
// u-index along the trace and refuse any triangle within `exclude_u` of the
// query point's own u. What survives is the neighbouring wrap.
//
// This is the whole reason the mode exists: it never needs a scroll axis or a
// radius. Fitting radius-from-an-axis fails its own positive control on the
// labelled sheets (residual sd 79.8 vx against a 12-17 vx sheet spacing),
// because a Herculaneum scroll is crushed, not round. Comparing a trace to
// itself is invariant to that deformation.
//
// A correct trace reads one sheet spacing here. Two spacings means it skipped
// a wrap; near zero means it doubled back onto one it had already traced.
struct QueryPoint {
    Vec3 p;
    int32_t group = 0;
};

// ---------------------------------------------------------------- uniform grid

struct Grid {
    float ox = 0, oy = 0, oz = 0;  // origin
    float cell = 32.0f;
    int32_t nx = 0, ny = 0, nz = 0;
    std::vector<uint32_t> start;  // size nx*ny*nz + 1
    std::vector<uint32_t> items;  // triangle indices, bucketed

    inline int64_t index(int32_t i, int32_t j, int32_t k) const {
        return (static_cast<int64_t>(k) * ny + j) * nx + i;
    }
    inline int32_t ci(float v, float o) const {
        return static_cast<int32_t>(std::floor((v - o) / cell));
    }
};

Grid build_grid(const std::vector<Tri>& tris, float cell) {
    Grid g;
    g.cell = cell;
    float lo[3] = {1e30f, 1e30f, 1e30f}, hi[3] = {-1e30f, -1e30f, -1e30f};
    for (const Tri& t : tris) {
        for (const Vec3* v : {&t.a, &t.b, &t.c}) {
            lo[0] = std::min(lo[0], v->x); hi[0] = std::max(hi[0], v->x);
            lo[1] = std::min(lo[1], v->y); hi[1] = std::max(hi[1], v->y);
            lo[2] = std::min(lo[2], v->z); hi[2] = std::max(hi[2], v->z);
        }
    }
    g.ox = lo[0] - cell; g.oy = lo[1] - cell; g.oz = lo[2] - cell;
    g.nx = static_cast<int32_t>((hi[0] - g.ox) / cell) + 2;
    g.ny = static_cast<int32_t>((hi[1] - g.oy) / cell) + 2;
    g.nz = static_cast<int32_t>((hi[2] - g.oz) / cell) + 2;

    const int64_t ncells = static_cast<int64_t>(g.nx) * g.ny * g.nz;
    std::fprintf(stderr, "grid %dx%dx%d = %lld cells, cell=%.1f\n", g.nx, g.ny, g.nz,
                 static_cast<long long>(ncells), cell);

    // two passes: count, then fill
    std::vector<uint32_t> counts(ncells + 1, 0);
    auto for_each_cell = [&](const Tri& t, auto&& fn) {
        const float tlo[3] = {std::min({t.a.x, t.b.x, t.c.x}), std::min({t.a.y, t.b.y, t.c.y}),
                              std::min({t.a.z, t.b.z, t.c.z})};
        const float thi[3] = {std::max({t.a.x, t.b.x, t.c.x}), std::max({t.a.y, t.b.y, t.c.y}),
                              std::max({t.a.z, t.b.z, t.c.z})};
        const int32_t i0 = g.ci(tlo[0], g.ox), i1 = g.ci(thi[0], g.ox);
        const int32_t j0 = g.ci(tlo[1], g.oy), j1 = g.ci(thi[1], g.oy);
        const int32_t k0 = g.ci(tlo[2], g.oz), k1 = g.ci(thi[2], g.oz);
        for (int32_t k = k0; k <= k1; ++k)
            for (int32_t j = j0; j <= j1; ++j)
                for (int32_t i = i0; i <= i1; ++i) fn(g.index(i, j, k));
    };

    for (const Tri& t : tris) for_each_cell(t, [&](int64_t c) { counts[c + 1]++; });
    for (int64_t c = 0; c < ncells; ++c) counts[c + 1] += counts[c];
    g.start = counts;
    g.items.resize(counts[ncells]);
    std::vector<uint32_t> cursor(g.start.begin(), g.start.end() - 1);
    for (uint32_t ti = 0; ti < tris.size(); ++ti)
        for_each_cell(tris[ti], [&](int64_t c) { g.items[cursor[c]++] = ti; });

    std::fprintf(stderr, "grid insertions: %zu (%.2f per triangle)\n", g.items.size(),
                 static_cast<double>(g.items.size()) / static_cast<double>(tris.size()));
    return g;
}

struct Result {
    int32_t w1 = -1;
    float d1 = std::numeric_limits<float>::infinity();
    int32_t w2 = -1;
    float d2 = std::numeric_limits<float>::infinity();
};

// Expanding-shell nearest search. Stops when the closest possible point in the
// next shell is farther than the best found so far -- correctness depends on
// that bound, so it is computed from the shell's inner face, not its centre.
Result query(const QueryPoint& qp, const std::vector<Tri>& tris, const Grid& g,
             int32_t max_shell, int32_t exclude_u, float max_dist) {
    Result r;
    const Vec3& p = qp.p;
    const int32_t i0 = g.ci(p.x, g.ox), j0 = g.ci(p.y, g.oy), k0 = g.ci(p.z, g.oz);

    for (int32_t s = 0; s <= max_shell; ++s) {
        // Nothing in shell s can be closer than (s-1)*cell.
        //
        // Two exits. The second matters: keying the exit solely on d2 means a
        // query that never finds a second winding -- a single-winding atlas, or
        // a point out beyond the edge of the reference set -- walks every shell
        // in the grid. `max_shell` bounds that. Results at the bound are
        // reported as +inf rather than as a large finite distance, because
        // "no surface within the search radius" is a different statement from
        // "the nearest surface is exactly this far", and conflating them would
        // let an unanswered query masquerade as a confident one.
        if (s > 0) {
            const float floor_dist = static_cast<float>(s - 1) * g.cell;
            if (r.w2 >= 0 && floor_dist * floor_dist > r.d2) break;
        }
        for (int32_t k = k0 - s; k <= k0 + s; ++k) {
            if (k < 0 || k >= g.nz) continue;
            for (int32_t j = j0 - s; j <= j0 + s; ++j) {
                if (j < 0 || j >= g.ny) continue;
                for (int32_t i = i0 - s; i <= i0 + s; ++i) {
                    if (i < 0 || i >= g.nx) continue;
                    // shell only: skip the interior, already visited
                    const bool on_shell = (std::abs(i - i0) == s || std::abs(j - j0) == s ||
                                           std::abs(k - k0) == s);
                    if (!on_shell) continue;
                    const int64_t c = g.index(i, j, k);
                    for (uint32_t n = g.start[c]; n < g.start[c + 1]; ++n) {
                        const Tri& t = tris[g.items[n]];
                        // self-gap mode: ignore the query's own neighbourhood
                        if (exclude_u > 0 && std::abs(t.winding - qp.group) < exclude_u) continue;
                        const float d2v = dist2_point_triangle(p, t.a, t.b, t.c);
                        if (t.winding == r.w1) {
                            if (d2v < r.d1) r.d1 = d2v;
                        } else if (d2v < r.d1) {
                            r.w2 = r.w1; r.d2 = r.d1;
                            r.w1 = t.winding; r.d1 = d2v;
                        } else if (t.winding != r.w2 && d2v < r.d2) {
                            r.w2 = t.winding; r.d2 = d2v;
                        } else if (t.winding == r.w2 && d2v < r.d2) {
                            r.d2 = d2v;
                        }
                    }
                }
            }
        }
    }
    r.d1 = std::sqrt(r.d1);
    r.d2 = std::sqrt(r.d2);

    // A result is only certified global if the search actually covered a radius
    // larger than it. Exiting by exhausting `max_shell` does NOT certify
    // anything beyond `max_dist`, and returning such a value as a plain float
    // makes an uncertified answer indistinguishable from a certified one.
    //
    // This was a real bug, caught by adversarial review: the null control on a
    // single-wrap sheet reported "0.5% finite, median 331 vx" under
    // max_dist=256, and it was dismissed as far-field noise. It was this.
    // Nearby hits are unaffected -- a d1 well inside max_dist was certified by
    // the shell bound before the loop ended -- so no sub-6 vx result changes.
    const float certified = max_dist - g.cell;
    if (r.d1 > certified) { r.d1 = std::numeric_limits<float>::infinity(); r.w1 = -1; }
    if (r.d2 > certified) { r.d2 = std::numeric_limits<float>::infinity(); r.w2 = -1; }
    return r;
}

// ------------------------------------------------------------------------ io

template <typename T>
bool read_exact(std::FILE* f, T* dst, size_t n) {
    return std::fread(dst, sizeof(T), n, f) == n;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr,
                     "usage: %s <atlas.bin> <query.bin> <out.bin> [threads] [cell] [max_dist] [exclude_u]\n"
                     "  atlas.bin : 'WCAT' u32 version, u32 n_surf,\n"
                     "              per surface: i32 winding, u32 rows, u32 cols,\n"
                     "                           rows*cols*3 f32 xyz, rows*cols u8 valid\n"
                     "  query.bin : 'WCQP' u32 version, u32 n, n*3 f32\n"
                     "  out.bin   : per query: i32 w1, f32 d1, i32 w2, f32 d2\n",
                     argv[0]);
        return 2;
    }
    const int threads = argc > 4 ? std::atoi(argv[4])
                                 : static_cast<int>(std::thread::hardware_concurrency());
    const float cell = argc > 5 ? static_cast<float>(std::atof(argv[5])) : 32.0f;
    const int32_t exclude_u = argc > 7 ? std::atoi(argv[7]) : 0;
    const bool self_gap = exclude_u > 0;

    // ---- load atlas, tessellate into triangles
    std::FILE* fa = std::fopen(argv[1], "rb");
    if (!fa) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    char magic[4];
    uint32_t version = 0, n_surf = 0;
    if (!read_exact(fa, magic, 4) || std::memcmp(magic, "WCAT", 4) != 0 ||
        !read_exact(fa, &version, 1) || !read_exact(fa, &n_surf, 1)) {
        std::fprintf(stderr, "bad atlas header\n"); return 1;
    }

    std::vector<Tri> tris;
    for (uint32_t s = 0; s < n_surf; ++s) {
        int32_t winding = 0;
        uint32_t rows = 0, cols = 0;
        if (!read_exact(fa, &winding, 1) || !read_exact(fa, &rows, 1) || !read_exact(fa, &cols, 1)) {
            std::fprintf(stderr, "bad surface header %u\n", s); return 1;
        }

        // rows == 0 marks a triangle soup rather than a grid, so an unstructured
        // mesh (OBJ) can be checked with the same kernel. `cols` then holds the
        // triangle count, and each triangle carries its own tag -- for a grid
        // the tag is derived from the column index, here it comes from the
        // mesh's own vt parametrisation.
        if (rows == 0) {
            const uint32_t n_tri = cols;
            std::vector<float> buf(9);
            int32_t tag = 0;
            for (uint32_t t = 0; t < n_tri; ++t) {
                if (!read_exact(fa, buf.data(), 9) || !read_exact(fa, &tag, 1)) {
                    std::fprintf(stderr, "truncated triangle soup %u\n", s); return 1;
                }
                tris.push_back({{buf[0], buf[1], buf[2]},
                                {buf[3], buf[4], buf[5]},
                                {buf[6], buf[7], buf[8]},
                                self_gap ? tag : winding});
            }
            continue;
        }
        std::vector<float> xyz(static_cast<size_t>(rows) * cols * 3);
        std::vector<uint8_t> valid(static_cast<size_t>(rows) * cols);
        if (!read_exact(fa, xyz.data(), xyz.size()) ||
            !read_exact(fa, valid.data(), valid.size())) {
            std::fprintf(stderr, "truncated surface %u\n", s); return 1;
        }
        auto at = [&](uint32_t r, uint32_t c) {
            const size_t i = (static_cast<size_t>(r) * cols + c) * 3;
            return Vec3{xyz[i], xyz[i + 1], xyz[i + 2]};
        };
        auto ok = [&](uint32_t r, uint32_t c) {
            return valid[static_cast<size_t>(r) * cols + c] != 0;
        };
        // Each fully-valid 2x2 cell becomes two triangles. Cells with a missing
        // corner are dropped rather than patched: a hole in the reference atlas
        // is honest, an invented surface is not.
        for (uint32_t r = 0; r + 1 < rows; ++r) {
            for (uint32_t c = 0; c + 1 < cols; ++c) {
                if (!ok(r, c) || !ok(r, c + 1) || !ok(r + 1, c) || !ok(r + 1, c + 1)) continue;
                const Vec3 p00 = at(r, c), p01 = at(r, c + 1);
                const Vec3 p10 = at(r + 1, c), p11 = at(r + 1, c + 1);
                // In self-gap mode the tag is the column index, so the query can
                // exclude its own neighbourhood along the trace. Tagging by
                // column (not by a coarser band) avoids cutting seams into the
                // surface, which would show up as spurious gaps.
                const int32_t tag = self_gap ? static_cast<int32_t>(c) : winding;
                tris.push_back({p00, p01, p11, tag});
                tris.push_back({p00, p11, p10, tag});
            }
        }
    }
    std::fclose(fa);
    std::fprintf(stderr, "atlas: %u surfaces, %zu triangles\n", n_surf, tris.size());
    if (tris.empty()) { std::fprintf(stderr, "empty atlas\n"); return 1; }

    const Grid grid = build_grid(tris, cell);

    // ---- load queries
    std::FILE* fq = std::fopen(argv[2], "rb");
    if (!fq) { std::fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
    uint32_t nq = 0;
    if (!read_exact(fq, magic, 4) || !read_exact(fq, &version, 1) || !read_exact(fq, &nq, 1)) {
        std::fprintf(stderr, "bad query header\n"); return 1;
    }
    const bool has_groups = std::memcmp(magic, "WCQ2", 4) == 0;
    if (!has_groups && std::memcmp(magic, "WCQP", 4) != 0) {
        std::fprintf(stderr, "bad query magic\n"); return 1;
    }
    if (self_gap && !has_groups) {
        std::fprintf(stderr, "self-gap mode needs a WCQ2 query file (with u-indices)\n");
        return 1;
    }
    std::vector<QueryPoint> pts(nq);
    {
        std::vector<float> xyz(static_cast<size_t>(nq) * 3);
        if (!read_exact(fq, xyz.data(), xyz.size())) {
            std::fprintf(stderr, "truncated queries\n"); return 1;
        }
        std::vector<int32_t> grp(nq, 0);
        if (has_groups && !read_exact(fq, grp.data(), nq)) {
            std::fprintf(stderr, "truncated query groups\n"); return 1;
        }
        for (uint32_t i = 0; i < nq; ++i)
            pts[i] = {{xyz[i * 3], xyz[i * 3 + 1], xyz[i * 3 + 2]}, grp[i]};
    }
    std::fclose(fq);
    std::fprintf(stderr, "queries: %u, threads: %d\n", nq, threads);

    // ---- run
    // Search radius. A point hundreds of voxels from any reference sheet
    // carries no winding information, so there is nothing to gain by finding
    // its nearest surface exactly -- and an unbounded walk over a sparse atlas
    // costs the whole grid per query.
    const float max_dist = argc > 6 ? static_cast<float>(std::atof(argv[6])) : 256.0f;
    const int32_t max_shell = static_cast<int32_t>(std::ceil(max_dist / cell)) + 1;
    std::fprintf(stderr, "max_dist %.0f -> max_shell %d\n", max_dist, max_shell);
    std::vector<Result> out(nq);
    std::atomic<uint32_t> next{0};
    const uint32_t chunk = 4096;
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&] {
            for (;;) {
                const uint32_t lo = next.fetch_add(chunk);
                if (lo >= nq) return;
                const uint32_t hi = std::min(lo + chunk, nq);
                for (uint32_t i = lo; i < hi; ++i)
                    out[i] = query(pts[i], tris, grid, max_shell, exclude_u, max_dist);
            }
        });
    }
    for (std::thread& th : pool) th.join();

    // ---- write
    std::FILE* fo = std::fopen(argv[3], "wb");
    if (!fo) { std::fprintf(stderr, "cannot open %s\n", argv[3]); return 1; }
    for (const Result& r : out) {
        std::fwrite(&r.w1, sizeof(int32_t), 1, fo);
        std::fwrite(&r.d1, sizeof(float), 1, fo);
        std::fwrite(&r.w2, sizeof(int32_t), 1, fo);
        std::fwrite(&r.d2, sizeof(float), 1, fo);
    }
    std::fclose(fo);
    std::fprintf(stderr, "wrote %u results to %s\n", nq, argv[3]);
    return 0;
}
