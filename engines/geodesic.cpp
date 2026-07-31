// Intrinsic-separation spectrum engine: C++ port of src/windcheck/intrinsic.py.
//
// A transverse self-intersection has two preimages in the mesh grid that occupy
// one place in space. This engine measures, per crossing EVENT, the distance
// between them ALONG the censused surface: crossing pairs are grouped in
// product space with an orientation parity (conflicts flagged ambiguous, never
// measured), each pair's triangle-triangle intersection segment is recomputed
// and its endpoints enter the graph as virtual seeds tied barycentrically to
// their OWN triangle's corners, and a multi-source early-exit Dijkstra returns
// the exact locus-to-locus walk on the retained-triangle edge graph.
//
// THIS FILE IS A PORT, NOT A REDESIGN. The Python module is the reference; its
// review-locked semantics (rounds 11-14, notes/DECISIONS.md) are replicated
// operation for operation:
//   - retained quad: four valid corners AND six pairwise corner distances
//     (four sides + both diagonals) <= maxedge; maxedge 0 disables.
//   - graph: vertices = grid points in >= 1 retained quad; edges = the four
//     sides of every retained quad (each emitted once) + the chosen diagonal.
//   - float discipline: tifxyz points are float32. numpy computes the edge
//     weights as SEQUENTIAL float32 arithmetic (subtract, square, sum, sqrt
//     all in f32) and only widens the stored weight to double inside the
//     Dijkstra relaxation. This port does exactly that. All other geometry
//     (tri-tri segments, seed offsets) is double from float32 inputs, with
//     numpy's left-to-right accumulation order. Build with -ffp-contract=off:
//     a fused multiply-add would round differently from numpy and break the
//     bit-level regression contract.
//   - the exact answer is the contract; the Python pop-budget fallback is an
//     implementation detail this port does not need (its Dijkstra always runs
//     to the exact early exit), so distance_method is always "early_exit" and
//     distance_exact is always true.
//
// Build (plain):
//   clang++ -O3 -std=c++17 -pthread -ffp-contract=off -Wall -Wextra \
//       -o engines/geodesic engines/geodesic.cpp
// Build (jemalloc, optional; helps the many small allocations in grouping):
//   clang++ -O3 -std=c++17 -pthread -ffp-contract=off -Wall -Wextra \
//       -L$(brew --prefix jemalloc)/lib -ljemalloc \
//       -o engines/geodesic engines/geodesic.cpp
//
// Usage:
//   ./geodesic <atlas.bin> <pairs.csv> <out.json> [threads=0(auto)]
//              [diagonal=0] [maxedge=60] [voxel_um=7.91]
//
// The atlas is the same WCAT format selfcross consumes (windcheck.atlas
// .write_atlas); the pairs CSV is a selfcross output CSV, of which only the
// rows classified "transverse" are used, bounds-filtered exactly as
// bench/spectrum_census.py filters them. Output is a JSON array of event
// objects carrying the same field names and values as segment_spectrum, plus
// region_a/region_b (the parity-normalised quad sets) so a regression harness
// can match events by region signature rather than by order.
//
// Debug modes (regression harness / analytic fixtures):
//   ./geodesic --selftest
//   ./geodesic --test <atlas.bin> <diagonal> <maxedge> dist v1 u1 v2 u2
//   ./geodesic --test <atlas.bin> <diagonal> <maxedge> idx v u
//   ./geodesic --test <atlas.bin> <diagonal> <maxedge> ncomp
//   ./geodesic --test <atlas.bin> <diagonal> <maxedge> seeded nA (v u off)*nA
//                                                            nB (v u off)*nB

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <queue>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr double INF = std::numeric_limits<double>::infinity();

// ------------------------------------------------------------------ vectors

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

// numpy accumulates a 3-term dot left to right; keep that order exactly.
inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// np.linalg.norm on float32 rows: subtraction, squares, left-to-right sum and
// sqrt ALL in float32. The stored edge weight is this f32 value; it is only
// widened to double where Python widens it (weights[..].tolist() in the
// Dijkstra loop, np.maximum(float64, ...) in the maxedge test).
inline float f32dist(const float* a, const float* b) {
    float dx = a[0] - b[0];
    float dy = a[1] - b[1];
    float dz = a[2] - b[2];
    float s = dx * dx;
    s = s + dy * dy;
    s = s + dz * dz;
    return std::sqrt(s);
}

// --------------------------------------------------------------------- atlas

struct Surface {
    int64_t rows = 0, cols = 0;
    std::vector<float> pts;      // rows*cols*3
    std::vector<uint8_t> valid;  // rows*cols
    inline bool ok(int64_t v, int64_t u) const {
        return valid[(size_t)(v * cols + u)] != 0;
    }
    inline const float* at(int64_t v, int64_t u) const {
        return &pts[(size_t)(v * cols + u) * 3];
    }
};

bool read_atlas(const char* path, Surface& S) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    char magic[4];
    f.read(magic, 4);
    if (std::memcmp(magic, "WCAT", 4) != 0) return false;
    uint32_t version = 0, count = 0;
    f.read((char*)&version, 4);
    f.read((char*)&count, 4);
    if (count < 1) return false;
    int32_t winding = 0;
    uint32_t rows = 0, cols = 0;
    f.read((char*)&winding, 4);
    f.read((char*)&rows, 4);
    f.read((char*)&cols, 4);
    if (rows == 0) { std::fprintf(stderr, "geodesic: OBJ soup not supported\n"); return false; }
    S.rows = rows; S.cols = cols;
    S.pts.resize((size_t)rows * cols * 3);
    S.valid.resize((size_t)rows * cols);
    f.read((char*)S.pts.data(), (std::streamsize)S.pts.size() * 4);
    f.read((char*)S.valid.data(), (std::streamsize)S.valid.size());
    return (bool)f;
}

// ---------------------------------------------------------------- pairs CSV

// windcheck.check.load_pairs: only rows whose 5th column is "transverse";
// columns v1,u1,v2,u2,verdict, then penetration (the intersection-segment
// length) and angle when the header carries them.
struct PairRow {
    int32_t v1, u1, v2, u2;
    double pen;
};

bool load_pairs(const char* path, std::vector<PairRow>& out) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f, line)) return true;          // empty file: no rows
    const bool has_margin = line.find("penetration") != std::string::npos;
    std::vector<std::string> p;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        p.clear();
        size_t start = 0;
        while (true) {
            size_t c = line.find(',', start);
            if (c == std::string::npos) { p.push_back(line.substr(start)); break; }
            p.push_back(line.substr(start, c - start));
            start = c + 1;
        }
        if (p.size() < 5 || p[4] != "transverse") continue;
        PairRow r;
        r.v1 = (int32_t)std::strtol(p[0].c_str(), nullptr, 10);
        r.u1 = (int32_t)std::strtol(p[1].c_str(), nullptr, 10);
        r.v2 = (int32_t)std::strtol(p[2].c_str(), nullptr, 10);
        r.u2 = (int32_t)std::strtol(p[3].c_str(), nullptr, 10);
        r.pen = (has_margin && p.size() > 5)
                    ? std::strtod(p[5].c_str(), nullptr)
                    : std::numeric_limits<double>::quiet_NaN();
        out.push_back(r);
    }
    return true;
}

// --------------------------------------------------------------- edge graph

// SurfaceGraph: the censused triangulation's edge graph, stride 1.
// Vertices are numbered in row-major grid order (numpy nonzero order); the
// numbering matters only internally but is kept identical to Python anyway.
struct Graph {
    int64_t nv = 0, nu = 0;             // grid shape
    int64_t n = 0;                      // vertex count
    int64_t ncomp = 0;
    int diagonal = 0;
    double maxedge = 60.0;
    std::vector<int64_t> idx;           // nv*nu -> vertex id or -1
    std::vector<uint8_t> Q;             // (nv-1)*(nu-1) retained-quad mask
    std::vector<size_t> indptr;         // CSR offsets, n+1
    std::vector<int64_t> indices;       // CSR neighbours
    std::vector<float> weights;         // CSR weights, float32 by contract
    std::vector<int64_t> comp;          // component label per vertex
    std::vector<double> X;              // n*3 vertex positions, double

    inline bool quad(int64_t v, int64_t u) const {
        return Q[(size_t)(v * (nu - 1) + u)] != 0;
    }
    inline int64_t vid(int64_t v, int64_t u) const {
        return idx[(size_t)(v * nu + u)];
    }
    inline Vec3 pos(int64_t i) const {
        const double* p = &X[(size_t)i * 3];
        return {p[0], p[1], p[2]};
    }
};

struct DSU {                             // plain union-find for components
    std::vector<int64_t> p;
    explicit DSU(int64_t n) : p(n) { for (int64_t i = 0; i < n; ++i) p[i] = i; }
    int64_t find(int64_t x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; }
    void unite(int64_t a, int64_t b) { a = find(a); b = find(b); if (a != b) p[b] = a; }
};

void build_graph(const Surface& S, int diagonal, double maxedge, Graph& g) {
    const int64_t nv = S.rows, nu = S.cols;
    g.nv = nv; g.nu = nu; g.diagonal = diagonal; g.maxedge = maxedge;
    assert(nv >= 1 && nu >= 1);

    // retained quads: four valid corners, six f32 corner distances <= maxedge
    g.Q.assign((size_t)std::max<int64_t>(0, (nv - 1) * (nu - 1)), 0);
    for (int64_t v = 0; v + 1 < nv; ++v) {
        for (int64_t u = 0; u + 1 < nu; ++u) {
            if (!S.ok(v, u) || !S.ok(v + 1, u) || !S.ok(v, u + 1) || !S.ok(v + 1, u + 1))
                continue;
            if (maxedge > 0) {
                const float* p00 = S.at(v, u);
                const float* p10 = S.at(v + 1, u);
                const float* p01 = S.at(v, u + 1);
                const float* p11 = S.at(v + 1, u + 1);
                // the six distances in retained_quads order; each is an f32
                // value compared against maxedge after widening, exactly as
                // np.maximum(float64_zeros, f32_norm) <= maxedge does
                double e = 0.0;
                e = std::max(e, (double)f32dist(p00, p01));
                e = std::max(e, (double)f32dist(p01, p11));
                e = std::max(e, (double)f32dist(p11, p10));
                e = std::max(e, (double)f32dist(p00, p10));
                e = std::max(e, (double)f32dist(p00, p11));
                e = std::max(e, (double)f32dist(p10, p01));
                if (!(e <= maxedge)) continue;
            }
            g.Q[(size_t)(v * (nu - 1) + u)] = 1;
        }
    }

    // vertices: grid points in >= 1 retained quad, numbered row-major
    g.idx.assign((size_t)(nv * nu), -1);
    auto used_by_quad = [&](int64_t v, int64_t u) {   // any retained quad touching (v,u)
        for (int64_t dv = -1; dv <= 0; ++dv)
            for (int64_t du = -1; du <= 0; ++du) {
                int64_t qv = v + dv, qu = u + du;
                if (qv >= 0 && qv < nv - 1 && qu >= 0 && qu < nu - 1 && g.quad(qv, qu))
                    return true;
            }
        return false;
    };
    int64_t next_id = 0;
    for (int64_t v = 0; v < nv; ++v)
        for (int64_t u = 0; u < nu; ++u)
            if (used_by_quad(v, u)) g.idx[(size_t)(v * nu + u)] = next_id++;
    g.n = next_id;

    // vertex positions widened to double once (Python's X array)
    g.X.assign((size_t)g.n * 3, 0.0);
    for (int64_t v = 0; v < nv; ++v)
        for (int64_t u = 0; u < nu; ++u) {
            int64_t i = g.vid(v, u);
            if (i < 0) continue;
            const float* p = S.at(v, u);
            g.X[(size_t)i * 3 + 0] = p[0];
            g.X[(size_t)i * 3 + 1] = p[1];
            g.X[(size_t)i * 3 + 2] = p[2];
        }

    // edges, each emitted once: quad sides + chosen diagonal
    struct E { int64_t a, b; float w; };
    std::vector<E> edges;
    edges.reserve((size_t)std::max<int64_t>(0, 3 * (nv - 1) * (nu - 1) / 2));
    auto quad_at = [&](int64_t v, int64_t u) {
        return v >= 0 && v < nv - 1 && u >= 0 && u < nu - 1 && g.quad(v, u);
    };
    for (int64_t v = 0; v < nv; ++v)                     // (v,u)-(v,u+1)
        for (int64_t u = 0; u + 1 < nu; ++u)
            if (quad_at(v, u) || quad_at(v - 1, u))
                edges.push_back({g.vid(v, u), g.vid(v, u + 1),
                                 f32dist(S.at(v, u), S.at(v, u + 1))});
    for (int64_t v = 0; v + 1 < nv; ++v)                 // (v,u)-(v+1,u)
        for (int64_t u = 0; u < nu; ++u)
            if (quad_at(v, u) || quad_at(v, u - 1))
                edges.push_back({g.vid(v, u), g.vid(v + 1, u),
                                 f32dist(S.at(v, u), S.at(v + 1, u))});
    if (diagonal == 0) {
        for (int64_t v = 0; v + 1 < nv; ++v)
            for (int64_t u = 0; u + 1 < nu; ++u)
                if (g.quad(v, u))
                    edges.push_back({g.vid(v, u), g.vid(v + 1, u + 1),
                                     f32dist(S.at(v, u), S.at(v + 1, u + 1))});
    } else if (diagonal == 1) {
        for (int64_t v = 0; v + 1 < nv; ++v)
            for (int64_t u = 0; u + 1 < nu; ++u)
                if (g.quad(v, u))
                    edges.push_back({g.vid(v + 1, u), g.vid(v, u + 1),
                                     f32dist(S.at(v + 1, u), S.at(v, u + 1))});
    }
    // diagonal == -1 (perimeter only) exists in Python for deformation fields;
    // the spectrum never uses it, so this port does not offer it.

    // CSR, both directions
    g.indptr.assign((size_t)g.n + 1, 0);
    for (const E& e : edges) {
        assert(e.a >= 0 && e.a < g.n && e.b >= 0 && e.b < g.n);
        ++g.indptr[(size_t)e.a + 1];
        ++g.indptr[(size_t)e.b + 1];
    }
    for (size_t i = 1; i < g.indptr.size(); ++i) g.indptr[i] += g.indptr[i - 1];
    g.indices.assign(g.indptr.back(), 0);
    g.weights.assign(g.indptr.back(), 0.f);
    std::vector<size_t> cur(g.indptr.begin(), g.indptr.end() - 1);
    for (const E& e : edges) {
        size_t pa = cur[(size_t)e.a]++;
        size_t pb = cur[(size_t)e.b]++;
        g.indices[pa] = e.b; g.weights[pa] = e.w;
        g.indices[pb] = e.a; g.weights[pb] = e.w;
    }

    // connected components
    DSU dsu(g.n);
    for (const E& e : edges) dsu.unite(e.a, e.b);
    g.comp.assign((size_t)g.n, -1);
    std::unordered_map<int64_t, int64_t> label;
    label.reserve(64);
    for (int64_t i = 0; i < g.n; ++i) {
        int64_t r = dsu.find(i);
        auto it = label.find(r);
        if (it == label.end()) it = label.emplace(r, (int64_t)label.size()).first;
        g.comp[(size_t)i] = it->second;
    }
    g.ncomp = (int64_t)label.size();
}

// quad_triangles, selfcross's corner order (intrinsic.SurfaceGraph)
inline void quad_triangles(const Graph& g, int64_t v, int64_t u, int64_t t[2][3]) {
    assert(v >= 0 && v + 1 < g.nv && u >= 0 && u + 1 < g.nu);
    if (g.diagonal == 0) {
        t[0][0] = g.vid(v, u);     t[0][1] = g.vid(v, u + 1);     t[0][2] = g.vid(v + 1, u + 1);
        t[1][0] = g.vid(v, u);     t[1][1] = g.vid(v + 1, u + 1); t[1][2] = g.vid(v + 1, u);
    } else {
        t[0][0] = g.vid(v, u);     t[0][1] = g.vid(v, u + 1);     t[0][2] = g.vid(v + 1, u);
        t[1][0] = g.vid(v, u + 1); t[1][1] = g.vid(v + 1, u + 1); t[1][2] = g.vid(v + 1, u);
    }
}

// ----------------------------------------------------------- event grouping

// oriented_events: union-find over crossing-pair rows with orientation parity.
// The iteration order (rows ascending, q1 before q2, dv/du in -1,0,1 order,
// the j <= i skip and the per-row seen set) matches Python line for line so
// identical unions happen in the identical order.
struct Event {
    std::vector<int32_t> rows;                 // member row indices, ascending
    std::vector<uint8_t> flip;                 // parity per member, same order
    std::set<std::pair<int32_t, int32_t>> A, B;
    bool ambiguous = false;
    bool self_touching = false;
};

inline uint64_t qkey(int64_t v, int64_t u) {
    return ((uint64_t)(uint32_t)(int32_t)v << 32) | (uint32_t)(int32_t)u;
}
inline bool adj(int32_t av, int32_t au, int32_t bv, int32_t bu) {
    return std::abs(av - bv) <= 1 && std::abs(au - bu) <= 1;
}

std::vector<Event> oriented_events(const std::vector<PairRow>& rec) {
    const int32_t n = (int32_t)rec.size();
    std::vector<Event> out;
    if (n == 0) return out;

    std::vector<int32_t> parent(n), parity(n, 0);
    std::vector<uint8_t> conflict(n, 0);
    for (int32_t i = 0; i < n; ++i) parent[i] = i;

    std::vector<int32_t> chain;
    auto find = [&](int32_t x, int32_t& p_out) -> int32_t {
        chain.clear();
        while (parent[x] != x) { chain.push_back(x); x = parent[x]; }
        int32_t p = 0;
        for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
            p ^= parity[*it];
            parent[*it] = x;
            parity[*it] = p;
        }
        p_out = chain.empty() ? 0 : parity[chain.front()];
        return x;
    };
    auto unite = [&](int32_t a, int32_t b, int32_t rel) {
        int32_t pa, pb;
        int32_t ra = find(a, pa);
        int32_t rb = find(b, pb);
        if (ra == rb) {
            if ((pa ^ pb) != rel) conflict[ra] = 1;
            return;
        }
        parent[rb] = ra;
        parity[rb] = pa ^ pb ^ rel;
        conflict[ra] |= conflict[rb];
    };

    std::unordered_map<uint64_t, std::vector<int32_t>> cell;
    cell.reserve((size_t)n * 2);
    for (int32_t i = 0; i < n; ++i) {
        cell[qkey(rec[i].v1, rec[i].u1)].push_back(i);
        cell[qkey(rec[i].v2, rec[i].u2)].push_back(i);
    }

    std::vector<int32_t> seen_stamp(n, -1);
    for (int32_t i = 0; i < n; ++i) {
        const int32_t qs[2][2] = {{rec[i].v1, rec[i].u1}, {rec[i].v2, rec[i].u2}};
        for (int side = 0; side < 2; ++side) {
            const int32_t qv = qs[side][0], qu = qs[side][1];
            for (int32_t dv = -1; dv <= 1; ++dv) {
                for (int32_t du = -1; du <= 1; ++du) {
                    auto it = cell.find(qkey(qv + dv, qu + du));
                    if (it == cell.end()) continue;
                    for (int32_t j : it->second) {
                        if (j <= i || seen_stamp[j] == i) continue;
                        seen_stamp[j] = i;
                        const bool straight =
                            adj(rec[i].v1, rec[i].u1, rec[j].v1, rec[j].u1) &&
                            adj(rec[i].v2, rec[i].u2, rec[j].v2, rec[j].u2);
                        const bool swapped =
                            adj(rec[i].v1, rec[i].u1, rec[j].v2, rec[j].u2) &&
                            adj(rec[i].v2, rec[i].u2, rec[j].v1, rec[j].u1);
                        if (straight && swapped) {
                            unite(i, j, 0);
                            unite(i, j, 1);            // forces the conflict flag
                        } else if (straight) {
                            unite(i, j, 0);
                        } else if (swapped) {
                            unite(i, j, 1);
                        }
                    }
                }
            }
        }
    }

    // groups in first-encounter order (Python dict insertion order)
    std::unordered_map<int32_t, size_t> gidx;
    gidx.reserve((size_t)n);
    std::vector<int32_t> roots;
    for (int32_t i = 0; i < n; ++i) {
        int32_t p;
        int32_t root = find(i, p);
        auto it = gidx.find(root);
        if (it == gidx.end()) {
            it = gidx.emplace(root, out.size()).first;
            out.emplace_back();
            roots.push_back(root);
        }
        Event& ev = out[it->second];
        ev.rows.push_back(i);
        ev.flip.push_back((uint8_t)p);
        auto q1 = std::make_pair(rec[i].v1, rec[i].u1);
        auto q2 = std::make_pair(rec[i].v2, rec[i].u2);
        (p == 0 ? ev.A : ev.B).insert(q1);
        (p == 0 ? ev.B : ev.A).insert(q2);
    }
    for (size_t k = 0; k < out.size(); ++k) {
        Event& ev = out[k];
        ev.ambiguous = conflict[roots[k]] != 0;
        bool touching = false;
        for (const auto& a : ev.A) {
            for (int32_t dv = -1; dv <= 1 && !touching; ++dv)
                for (int32_t du = -1; du <= 1 && !touching; ++du)
                    if (ev.B.count({a.first + dv, a.second + du})) touching = true;
            if (touching) break;
        }
        ev.self_touching = touching;
    }
    return out;
}

// -------------------------------------------------- barycentric endpoints

// _tri_tri_segment, operation for operation. Interval method: each triangle
// clipped by the other's plane, both clipped segments projected on the
// normalised cross-product direction, overlap of the parameter intervals.
// Returns false for near-coplanar input, a degenerate clip, an empty overlap
// or a sub-eps first-segment span -- exactly the Python None cases.
bool tri_tri_segment(const Vec3 A[3], const Vec3 B[3], Vec3& p_lo, Vec3& p_hi) {
    const double eps = 1e-12;

    Vec3 n1 = cross(A[1] - A[0], A[2] - A[0]);
    Vec3 n2 = cross(B[1] - B[0], B[2] - B[0]);
    Vec3 d = cross(n1, n2);
    const double L = norm(d);
    if (L < eps * std::max(norm(n1), 1.0) * std::max(norm(n2), 1.0))
        return false;                            // near-coplanar: fallback
    d = {d.x / L, d.y / L, d.z / L};             // numpy divides, not *(1/L)

    struct TP { double t; Vec3 p; };
    auto clip = [](const Vec3 T[3], const Vec3& nrm, double d0,
                   Vec3 pts[2]) -> int {
        double s[3];
        for (int i = 0; i < 3; ++i) s[i] = dot(T[i], nrm) + d0;
        int np = 0;
        for (int i = 0; i < 3; ++i) {
            const int j = (i + 1) % 3;
            const double si = s[i], sj = s[j];
            if (si == 0.0) {
                if (np < 2) pts[np] = T[i];
                ++np;
            }
            if ((si > 0) != (sj > 0) && si != 0.0 && sj != 0.0) {
                const double t = si / (si - sj);
                if (np < 2) pts[np] = T[i] + (T[j] - T[i]) * t;
                ++np;
            }
        }
        return np;
    };
    // NOTE Python appends T[i] + t*(T[j]-T[i]) with t a scalar broadcast:
    // per component a + t*(b-a); Vec3 (b-a)*t + a reproduces it exactly.

    Vec3 segA[2], segB[2];
    if (clip(A, n2, -dot(n2, B[0]), segA) < 2) return false;
    if (clip(B, n1, -dot(n1, A[0]), segB) < 2) return false;

    auto lex_less = [](const Vec3& a, const Vec3& b) {
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return a.z < b.z;
    };
    auto order = [&](Vec3 s[2], TP o[2]) {       // Python sorted((t, tuple(p)))
        o[0] = {dot(s[0], d), s[0]};
        o[1] = {dot(s[1], d), s[1]};
        if (o[1].t < o[0].t ||
            (o[1].t == o[0].t && lex_less(o[1].p, o[0].p)))
            std::swap(o[0], o[1]);
    };
    TP ta[2], tb[2];
    order(segA, ta);
    order(segB, tb);

    const double lo = std::max(ta[0].t, tb[0].t);
    const double hi = std::min(ta[1].t, tb[1].t);
    if (hi <= lo) return false;
    const Vec3 a0 = ta[0].p, a1 = ta[1].p;
    const double span = ta[1].t - ta[0].t;
    if (span < eps) return false;
    const double s_lo = (lo - ta[0].t) / span;
    const double s_hi = (hi - ta[0].t) / span;
    p_lo = a0 + (a1 - a0) * s_lo;
    p_hi = a0 + (a1 - a0) * s_hi;
    return true;
}

// np.linalg.norm(X[c] - x) on double 3-vectors: left-to-right, double
inline double dnorm3(const Vec3& a, const Vec3& b) {
    const Vec3 v = a - b;
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

using Seed = std::pair<int64_t, double>;         // (vertex, offset)

// _pair_seeds: barycentric seeds for one crossing pair; returns endpoint_exact
bool pair_seeds(const Graph& g, int64_t v1, int64_t u1, int64_t v2, int64_t u2,
                std::vector<Seed>& sa, std::vector<Seed>& sb) {
    int64_t t1[2][3], t2[2][3];
    quad_triangles(g, v1, u1, t1);
    quad_triangles(g, v2, u2, t2);
    struct Combo { const int64_t* t1; const int64_t* t2; Vec3 p, q; };
    Combo combos[4];
    int nc = 0;
    for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
            bool bad = false;
            for (int k = 0; k < 3; ++k)
                if (t1[a][k] < 0 || t2[b][k] < 0) { bad = true; break; }
            if (bad) continue;
            Vec3 TA[3], TB[3];
            for (int k = 0; k < 3; ++k) { TA[k] = g.pos(t1[a][k]); TB[k] = g.pos(t2[b][k]); }
            Vec3 p, q;
            if (tri_tri_segment(TA, TB, p, q))
                combos[nc++] = {t1[a], t2[b], p, q};
        }
    }
    if (nc == 0) {
        // corner fallback, zero offsets, quad_corners order (v,u),(v+1,u),
        // (v,u+1),(v+1,u+1), only corners present in the graph
        const int64_t cs1[4][2] = {{v1, u1}, {v1 + 1, u1}, {v1, u1 + 1}, {v1 + 1, u1 + 1}};
        const int64_t cs2[4][2] = {{v2, u2}, {v2 + 1, u2}, {v2, u2 + 1}, {v2 + 1, u2 + 1}};
        for (int k = 0; k < 4; ++k) {
            int64_t c = g.vid(cs1[k][0], cs1[k][1]);
            if (c >= 0) sa.emplace_back(c, 0.0);
        }
        for (int k = 0; k < 4; ++k) {
            int64_t c = g.vid(cs2[k][0], cs2[k][1]);
            if (c >= 0) sb.emplace_back(c, 0.0);
        }
        return false;
    }
    for (int ci = 0; ci < nc; ++ci) {
        const Combo& c = combos[ci];
        for (int e = 0; e < 2; ++e) {
            const Vec3& x = e == 0 ? c.p : c.q;
            for (int k = 0; k < 3; ++k)
                sa.emplace_back(c.t1[k], dnorm3(g.pos(c.t1[k]), x));
            for (int k = 0; k < 3; ++k)
                sb.emplace_back(c.t2[k], dnorm3(g.pos(c.t2[k]), x));
        }
    }
    return true;
}

// ------------------------------------------------------------ seeded search

// Per-thread Dijkstra state: flat double array with a touched-list reset.
struct Dijk {
    std::vector<double> dist;
    std::vector<int64_t> touched;
    std::priority_queue<std::pair<double, int64_t>,
                        std::vector<std::pair<double, int64_t>>,
                        std::greater<std::pair<double, int64_t>>> heap;
    explicit Dijk(int64_t n) : dist((size_t)n, INF) {}
    void reset() {
        for (int64_t t : touched) dist[(size_t)t] = INF;
        touched.clear();
        while (!heap.empty()) heap.pop();
    }
};

// seeded_distance: min over virtual endpoints of (in-A weight + path + in-B
// weight). Multi-source Dijkstra from the A seeds; every settled vertex that
// carries a B offset updates the best answer; the search stops once the heap
// minimum can no longer beat it. Exact by construction -- the Python pop
// budget and super-source fallback are not needed here.
double seeded_distance(const Graph& g, Dijk& D,
                       const std::vector<Seed>& seeds_a,
                       const std::vector<Seed>& seeds_b) {
    std::unordered_map<int64_t, double> b_off;
    b_off.reserve(seeds_b.size() * 2);
    for (const Seed& s : seeds_b) {
        auto it = b_off.find(s.first);
        if (it == b_off.end()) b_off.emplace(s.first, s.second);
        else it->second = std::min(it->second, s.second);
    }
    std::unordered_set<int64_t> comp_a;
    for (const Seed& s : seeds_a) comp_a.insert(g.comp[(size_t)s.first]);
    bool reachable = false;
    for (const auto& kv : b_off)
        if (comp_a.count(g.comp[(size_t)kv.first])) { reachable = true; break; }
    if (!reachable) return INF;

    D.reset();
    for (const Seed& s : seeds_a) {
        if (s.second < D.dist[(size_t)s.first]) {
            if (D.dist[(size_t)s.first] == INF) D.touched.push_back(s.first);
            D.dist[(size_t)s.first] = s.second;
            D.heap.emplace(s.second, s.first);
        }
    }
    double best = INF;
    while (!D.heap.empty()) {
        auto [d, nd] = D.heap.top();
        D.heap.pop();
        if (d >= best) break;
        if (d > D.dist[(size_t)nd]) continue;
        auto bit = b_off.find(nd);
        if (bit != b_off.end()) best = std::min(best, d + bit->second);
        const size_t lo = g.indptr[(size_t)nd], hi = g.indptr[(size_t)nd + 1];
        for (size_t e = lo; e < hi; ++e) {
            const int64_t m = g.indices[e];
            const double cand = d + (double)g.weights[e];  // f32 weight widened
            if (cand < D.dist[(size_t)m]) {
                if (D.dist[(size_t)m] == INF) D.touched.push_back(m);
                D.dist[(size_t)m] = cand;
                D.heap.emplace(cand, m);
            }
        }
    }
    return best;
}

// --------------------------------------------------------------- spectrum

struct EventResult {
    int64_t n_pairs = 0;
    bool ambiguous = false;
    bool self_touching = false;
    int64_t du_max = 0;
    double median_len = 0.0;             // nanmedian of pen; may be NaN
    bool has_separation = false;
    double separation_mm = 0.0;
    int same_component = -1;             // -1 null, 0 false, 1 true
    int endpoint_exact = -1;             // -1 null, 0 false, 1 true
    int distance_exact = -1;             // -1 null, 1 true
    bool has_method = false;             // distance_method: "early_exit"|null
    const Event* ev = nullptr;
};

// np.nanmedian: ignore NaNs; all-NaN -> NaN; even count -> mean of middles
double nanmedian(std::vector<double>& vals) {
    std::vector<double> v;
    v.reserve(vals.size());
    for (double x : vals)
        if (!std::isnan(x)) v.push_back(x);
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const size_t m = v.size() / 2;
    if (v.size() % 2 == 1) return v[m];
    return (v[m - 1] + v[m]) / 2.0;
}

void measure_event(const Graph& g, const std::vector<PairRow>& rec,
                   const Event& ev, double mm, Dijk& D, EventResult& R) {
    R.ev = &ev;
    R.n_pairs = (int64_t)ev.rows.size();
    R.ambiguous = ev.ambiguous;
    R.self_touching = ev.self_touching;
    R.du_max = 0;
    std::vector<double> pens;
    pens.reserve(ev.rows.size());
    for (int32_t i : ev.rows) {
        R.du_max = std::max<int64_t>(R.du_max, std::abs((int64_t)rec[i].u1 - rec[i].u2));
        pens.push_back(rec[i].pen);
    }
    R.median_len = nanmedian(pens);
    if (ev.ambiguous) return;                    // reported, never measured

    // event_separation
    std::vector<Seed> seeds_a, seeds_b, sa, sb;
    bool exact = true;
    for (size_t k = 0; k < ev.rows.size(); ++k) {
        const PairRow& r = rec[ev.rows[k]];
        sa.clear(); sb.clear();
        exact &= pair_seeds(g, r.v1, r.u1, r.v2, r.u2, sa, sb);
        if (ev.flip[k]) std::swap(sa, sb);
        seeds_a.insert(seeds_a.end(), sa.begin(), sa.end());
        seeds_b.insert(seeds_b.end(), sb.begin(), sb.end());
    }
    if (seeds_a.empty() || seeds_b.empty()) {
        // Python's early return here lacks the distance_exact key and would
        // raise in segment_spectrum; it cannot occur on censused input. Keep
        // the flags null and warn so a regression run surfaces it.
        R.endpoint_exact = 0;
        std::fprintf(stderr, "geodesic: WARNING event with empty seed set\n");
        return;
    }
    std::unordered_set<int64_t> comps;
    for (const Seed& s : seeds_a) comps.insert(g.comp[(size_t)s.first]);
    for (const Seed& s : seeds_b) comps.insert(g.comp[(size_t)s.first]);
    const bool same = comps.size() == 1;
    const double sep = seeded_distance(g, D, seeds_a, seeds_b);
    const bool finite = std::isfinite(sep);
    R.has_separation = finite;
    if (finite) R.separation_mm = sep * mm;      // rounded at print time
    R.same_component = (same && finite) ? 1 : 0;
    R.endpoint_exact = exact ? 1 : 0;
    R.distance_exact = 1;
    R.has_method = true;
}

// ------------------------------------------------------------------ output

void write_json(std::FILE* fp, const std::vector<EventResult>& results) {
    std::fprintf(fp, "[");
    char buf[64];
    for (size_t i = 0; i < results.size(); ++i) {
        const EventResult& R = results[i];
        std::fprintf(fp, i ? ",\n {" : "\n {");
        std::fprintf(fp, "\"n_pairs\": %lld", (long long)R.n_pairs);
        std::fprintf(fp, ", \"ambiguous\": %s", R.ambiguous ? "true" : "false");
        std::fprintf(fp, ", \"self_touching\": %s", R.self_touching ? "true" : "false");
        std::fprintf(fp, ", \"du_max\": %lld", (long long)R.du_max);
        if (std::isnan(R.median_len))            // json.dumps writes NaN bare
            std::fprintf(fp, ", \"median_intersection_length_vx\": NaN");
        else
            std::fprintf(fp, ", \"median_intersection_length_vx\": %.17g", R.median_len);
        if (R.ambiguous || !R.has_separation) {
            std::fprintf(fp, ", \"separation_mm\": null");
        } else {
            // Python: round(sep*mm, 4). %.4f is the same correctly-rounded
            // 4-decimal value; a JSON parser recovers the identical double.
            std::snprintf(buf, sizeof buf, "%.4f", R.separation_mm);
            std::fprintf(fp, ", \"separation_mm\": %s", buf);
        }
        auto tri = [&](const char* name, int v) {
            if (v < 0) std::fprintf(fp, ", \"%s\": null", name);
            else std::fprintf(fp, ", \"%s\": %s", name, v ? "true" : "false");
        };
        tri("same_component", R.ambiguous ? -1 : R.same_component);
        tri("endpoint_exact", R.ambiguous ? -1 : R.endpoint_exact);
        tri("distance_exact", R.ambiguous ? -1 : R.distance_exact);
        if (R.has_method) std::fprintf(fp, ", \"distance_method\": \"early_exit\"");
        else std::fprintf(fp, ", \"distance_method\": null");
        auto region = [&](const char* name,
                          const std::set<std::pair<int32_t, int32_t>>& S) {
            std::fprintf(fp, ", \"%s\": [", name);
            bool first = true;
            for (const auto& q : S) {
                std::fprintf(fp, "%s[%d, %d]", first ? "" : ", ", q.first, q.second);
                first = false;
            }
            std::fprintf(fp, "]");
        };
        region("region_a", R.ev->A);
        region("region_b", R.ev->B);
        std::fprintf(fp, "}");
    }
    std::fprintf(fp, "\n]\n");
}

// ---------------------------------------------------------------- selftest

int selftest() {
    int fails = 0;
    auto check = [&](bool ok, const char* what) {
        std::fprintf(stderr, "  %-52s %s\n", what, ok ? "ok" : "FAIL");
        if (!ok) ++fails;
    };
    // tests/test_intrinsic.py::test_tri_tri_segment_exact
    {
        const Vec3 T1[3] = {{0, 0, 0}, {4, 0, 0}, {0, 4, 0}};
        const Vec3 T2[3] = {{0, 1, -1}, {4, 1, -1}, {0, 1, 3}};
        Vec3 p, q;
        bool got = tri_tri_segment(T1, T2, p, q);
        check(got, "tri_tri: flat vs wall intersects");
        if (got) {
            if (p.x > q.x) std::swap(p, q);
            auto near3 = [](const Vec3& a, double x, double y, double z) {
                return std::fabs(a.x - x) < 1e-9 && std::fabs(a.y - y) < 1e-9
                    && std::fabs(a.z - z) < 1e-9;
            };
            check(near3(p, 0, 1, 0) && near3(q, 3, 1, 0),
                  "tri_tri: segment is y=1, z=0, x in [0,3]");
        }
        const Vec3 T3[3] = {{0, 0, 5}, {4, 0, 5}, {0, 4, 5}};
        Vec3 pp, qq;
        check(!tri_tri_segment(T1, T3, pp, qq), "tri_tri: parallel planes -> none");
    }
    // nanmedian semantics
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> a{3.0, 1.0, 2.0};
        std::vector<double> b{4.0, nan, 1.0, 2.0};
        std::vector<double> c{nan, nan};
        check(nanmedian(a) == 2.0, "nanmedian: odd count");
        check(nanmedian(b) == 2.0, "nanmedian: NaN ignored");
        check(std::isnan(nanmedian(c)), "nanmedian: all-NaN -> NaN");
    }
    std::fprintf(stderr, "geodesic selftest: %s\n", fails ? "FAILED" : "passed");
    return fails ? 1 : 0;
}

// --------------------------------------------------------------- test mode

int test_mode(int argc, char** argv) {
    // geodesic --test <atlas> <diagonal> <maxedge> <cmd> [args...]
    if (argc < 6) { std::fprintf(stderr, "geodesic --test: missing args\n"); return 2; }
    Surface S;
    if (!read_atlas(argv[2], S)) { std::fprintf(stderr, "geodesic: bad atlas %s\n", argv[2]); return 1; }
    const int diagonal = std::atoi(argv[3]);
    const double maxedge = std::atof(argv[4]);
    Graph g;
    build_graph(S, diagonal, maxedge, g);
    Dijk D(g.n);
    const std::string cmd = argv[5];
    if (cmd == "ncomp") {
        std::printf("%lld\n", (long long)g.ncomp);
        return 0;
    }
    if (cmd == "idx" && argc >= 8) {
        int64_t v = std::atoll(argv[6]), u = std::atoll(argv[7]);
        if (v < 0 || v >= g.nv || u < 0 || u >= g.nu) { std::printf("-1\n"); return 0; }
        std::printf("%lld\n", (long long)g.vid(v, u));
        return 0;
    }
    if (cmd == "dist" && argc >= 10) {
        int64_t va = g.vid(std::atoll(argv[6]), std::atoll(argv[7]));
        int64_t vb = g.vid(std::atoll(argv[8]), std::atoll(argv[9]));
        if (va < 0 || vb < 0) { std::printf("novertex\n"); return 0; }
        double d = seeded_distance(g, D, {{va, 0.0}}, {{vb, 0.0}});
        if (std::isfinite(d)) std::printf("%.12g\n", d);
        else std::printf("inf\n");
        return 0;
    }
    if (cmd == "seeded" && argc >= 7) {
        int a = 6;
        auto read_seeds = [&](std::vector<Seed>& s) -> bool {
            if (a >= argc) return false;
            int k = std::atoi(argv[a++]);
            for (int i = 0; i < k; ++i) {
                if (a + 2 >= argc) return false;
                int64_t v = std::atoll(argv[a]), u = std::atoll(argv[a + 1]);
                double off = std::atof(argv[a + 2]);
                a += 3;
                int64_t id = g.vid(v, u);
                if (id < 0) return false;
                s.emplace_back(id, off);
            }
            return true;
        };
        std::vector<Seed> sa, sb;
        if (!read_seeds(sa) || !read_seeds(sb)) {
            std::fprintf(stderr, "geodesic --test seeded: bad seed list\n");
            return 2;
        }
        double d = seeded_distance(g, D, sa, sb);
        if (std::isfinite(d)) std::printf("%.12g\n", d);
        else std::printf("inf\n");
        return 0;
    }
    std::fprintf(stderr, "geodesic --test: unknown command %s\n", cmd.c_str());
    return 2;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc >= 2 && std::strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc >= 2 && std::strcmp(argv[1], "--test") == 0) return test_mode(argc, argv);
    if (argc < 4) {
        std::fprintf(stderr,
            "usage: %s <atlas.bin> <pairs.csv> <out.json> [threads=0(auto)] "
            "[diagonal=0] [maxedge=60] [voxel_um=7.91]\n"
            "       %s --selftest\n"
            "       %s --test <atlas.bin> <diagonal> <maxedge> <cmd> [args...]\n",
            argv[0], argv[0], argv[0]);
        return 2;
    }
    const char* atlas_path = argv[1];
    const char* pairs_path = argv[2];
    const char* out_path = argv[3];
    int nthreads = argc > 4 ? std::atoi(argv[4]) : 0;
    const int diagonal = argc > 5 ? std::atoi(argv[5]) : 0;
    const double maxedge = argc > 6 ? std::atof(argv[6]) : 60.0;
    const double voxel_um = argc > 7 ? std::atof(argv[7]) : 7.91;
    if (diagonal != 0 && diagonal != 1) {
        std::fprintf(stderr, "geodesic: diagonal must be 0 or 1\n");
        return 2;
    }
    if (nthreads <= 0) nthreads = (int)std::thread::hardware_concurrency();
    if (nthreads <= 0) nthreads = 4;

    Surface S;
    if (!read_atlas(atlas_path, S)) {
        std::fprintf(stderr, "geodesic: bad atlas %s\n", atlas_path);
        return 1;
    }
    // A missing CSV is an empty census, exactly as windcheck.check.load_pairs
    // treats it (segments with zero crossings have no CSV at all).
    std::vector<PairRow> raw;
    if (!load_pairs(pairs_path, raw))
        std::fprintf(stderr, "geodesic: no pairs csv %s (treating as empty)\n",
                     pairs_path);
    // the census bounds filter (bench/spectrum_census.py): quad indices must
    // address a quad of THIS grid
    std::vector<PairRow> rec;
    rec.reserve(raw.size());
    for (const PairRow& r : raw)
        if (r.v1 >= 0 && r.v1 < S.rows - 1 && r.v2 >= 0 && r.v2 < S.rows - 1 &&
            r.u1 >= 0 && r.u1 < S.cols - 1 && r.u2 >= 0 && r.u2 < S.cols - 1)
            rec.push_back(r);

    const auto t0 = std::chrono::steady_clock::now();
    if (rec.empty()) {
        std::FILE* fp = std::fopen(out_path, "w");
        if (!fp) { std::fprintf(stderr, "geodesic: cannot write %s\n", out_path); return 1; }
        std::fprintf(fp, "[]\n");
        std::fclose(fp);
        std::fprintf(stderr, "geodesic: grid %lldx%lld | pairs 0 | events 0 | nothing to do\n",
                     (long long)S.rows, (long long)S.cols);
        std::printf("{\"pairs\":0,\"events\":0,\"ambiguous\":0,\"diagonal\":%d,"
                    "\"maxedge\":%.1f}\n", diagonal, maxedge);
        return 0;
    }

    Graph g;
    build_graph(S, diagonal, maxedge, g);

    std::vector<Event> events = oriented_events(rec);
    // largest-first, stable on grouping order (Python's sorted is stable)
    std::stable_sort(events.begin(), events.end(),
                     [](const Event& a, const Event& b) {
                         return a.rows.size() > b.rows.size();
                     });

    const double mm = voxel_um / 1000.0;
    std::vector<EventResult> results(events.size());
    std::atomic<size_t> next{0};
    auto worker = [&]() {
        Dijk D(g.n);
        for (;;) {
            const size_t i = next.fetch_add(1);
            if (i >= events.size()) break;
            measure_event(g, rec, events[i], mm, D, results[i]);
        }
    };
    std::vector<std::thread> pool;
    const int nt = std::min<int>(nthreads, (int)std::max<size_t>(1, events.size()));
    for (int t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();

    std::FILE* fp = std::fopen(out_path, "w");
    if (!fp) { std::fprintf(stderr, "geodesic: cannot write %s\n", out_path); return 1; }
    write_json(fp, results);
    std::fclose(fp);

    size_t n_amb = 0, n_inter = 0;
    for (const EventResult& R : results) {
        if (R.ambiguous) ++n_amb;
        else if (R.same_component == 0) ++n_inter;
    }
    const double secs = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::fprintf(stderr,
        "geodesic: grid %lldx%lld, %lld vertices, %zu edges, %lld comps | "
        "pairs %zu (%zu dropped by bounds) | events %zu (%zu ambiguous, "
        "%zu inter-component) | diagonal %d, maxedge %.1f, %d threads | %.1fs\n",
        (long long)S.rows, (long long)S.cols, (long long)g.n,
        g.indices.size() / 2, (long long)g.ncomp, rec.size(),
        raw.size() - rec.size(), events.size(), n_amb, n_inter,
        diagonal, maxedge, nt, secs);
    std::printf("{\"pairs\":%zu,\"events\":%zu,\"ambiguous\":%zu,"
                "\"inter_component\":%zu,\"vertices\":%lld,\"edges\":%zu,"
                "\"diagonal\":%d,\"maxedge\":%.1f,\"seconds\":%.2f}\n",
                rec.size(), events.size(), n_amb, n_inter, (long long)g.n,
                g.indices.size() / 2, diagonal, maxedge, secs);
    return 0;
}
