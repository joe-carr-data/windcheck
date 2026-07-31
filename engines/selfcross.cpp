// Exact nonlocal self-intersection census for a tifxyz quad mesh.
//
// Everything windcheck has measured so far is PROXIMITY: "these two parts of the
// trace are within 6 voxels". That question is irreducibly ambiguous, because
// wraps in a crushed scroll genuinely lie very close together, so a small
// distance can never by itself mean the trace is wrong.
//
// A transverse crossing is different in kind. An embedded surface cannot pass
// through itself, no matter how tightly it is packed. There is no threshold to
// argue about and no interpretation to get wrong. So this tool answers a binary
// question about published data: do triangles of one trace actually cross each
// other, and where.
//
// Classification is deliberately finer than yes/no, because the three cases have
// different meanings:
//
//   TRANSVERSE  the two triangle interiors genuinely pass through each other.
//               This is the defect-shaped signal.
//   COPLANAR    the triangles are (near) coplanar and their projections overlap.
//               Two sheets pressed flat against each other look like this, so it
//               is reported, never treated as a crossing.
//   GRAZING     a vertex or edge lies within EPS of the other plane, so the sign
//               data that the interval test relies on is not trustworthy at
//               double precision. Reported separately rather than guessed at.
//
// Adjacency is excluded by grid index, not by geometry: two triangles whose quad
// origins are within `--exclude` cells in BOTH v and u are neighbours on the
// sheet, and their sharing space is what a surface does. Anything beyond that
// window is nonlocal, and a nonlocal crossing has no innocent reading.
//
// Build:
//   clang++ -O3 -std=c++17 -pthread -o engines/selfcross engines/selfcross.cpp
//
// Usage:
//   ./selfcross <atlas.bin> <out.csv> [threads] [cell] [exclude] [diagonal] [maxedge]
//
// `maxedge` drops any triangle with an edge longer than that. This is not
// cosmetic. A tifxyz grid can contain discontinuities where two adjacent valid
// cells sit far apart in 3D; the triangle built across such a gap spans many
// wraps and crosses everything it passes through, purely because the mesh has a
// hole. Observed edges reach 1848 voxels against a 20-voxel grid pitch, and
// they occur in the multi-wrap traces while the labelled single-wrap segments
// top out at 37, so leaving them in would manufacture exactly the result this
// census is meant to test. 0 disables the filter.
//
// `diagonal` selects which way each quad is split. The choice is not cosmetic: a
// twisted quad can cross under one diagonal and not the other, so running both
// and comparing separates real crossings from tessellation artifacts.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <mutex>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

// Plane-distance tolerance. A fixed 1e-6 voxels was wrong: the input
// coordinates are float32 and run to ~1.3e4 voxels, where one ULP is already
// ~1e-3 voxels. A tolerance a thousand times finer than the data's own
// representation cannot separate signal from representation noise, so the
// tolerance is derived from the operands instead.
//
// The distance is dot(unit_normal, p - q), so its error scales with the
// magnitude of the coordinates involved. FLT_EPSILON is the float32 relative
// step; the factor of 16 covers the subtraction, the three-term dot product and
// the cross product that produced the normal.
constexpr double FLT_EPS = 1.1920929e-7;
constexpr double EPS_SCALE = 16.0;

inline double plane_eps(double mag) {
    return std::max(1e-9, EPS_SCALE * FLT_EPS * mag);
}

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

struct Tri {
    Vec3 a, b, c;
    int32_t v = 0, u = 0;      // grid origin of the owning quad
    int32_t t = 0;             // local triangle index within the quad (0/1),
                               // in the documented order for the diagonal
};

enum Verdict { NONE = 0, TRANSVERSE = 1, COPLANAR = 2, GRAZING = 3 };

// A crossing is not a boolean. Two triangles whose intersection intervals merely
// touch at one endpoint share a point without either interior passing through
// the other, and calling that a crossing overstates the finding. So the test
// returns how far they actually penetrate, and the caller decides.
struct Result2 {
    Verdict verdict = NONE;
    double penetration = 0.0;   // length of interval overlap, voxels
    double angle_deg = 0.0;     // between the two triangle planes
};

// ---------------------------------------------------------------- coplanar 2D

// Project to the dominant plane and run an exact-ish 2D triangle overlap.
bool tri2d_overlap(const Vec3 tri1[3], const Vec3 tri2[3], const Vec3& n) {
    int drop = 0;
    double ax = std::fabs(n.x), ay = std::fabs(n.y), az = std::fabs(n.z);
    if (ay > ax && ay > az) drop = 1;
    else if (az > ax && az > ay) drop = 2;

    auto proj = [drop](const Vec3& p) -> std::pair<double, double> {
        if (drop == 0) return {p.y, p.z};
        if (drop == 1) return {p.x, p.z};
        return {p.x, p.y};
    };
    std::pair<double, double> P[3], Q[3];
    for (int i = 0; i < 3; ++i) { P[i] = proj(tri1[i]); Q[i] = proj(tri2[i]); }

    auto cross2 = [](const std::pair<double, double>& o,
                     const std::pair<double, double>& a,
                     const std::pair<double, double>& b) {
        return (a.first - o.first) * (b.second - o.second)
             - (a.second - o.second) * (b.first - o.first);
    };
    // separating-axis over both triangles' edges
    auto separated = [&](const std::pair<double, double> A[3],
                         const std::pair<double, double> B[3]) {
        for (int i = 0; i < 3; ++i) {
            int j = (i + 1) % 3;
            double ex = A[j].first - A[i].first;
            double ey = A[j].second - A[i].second;
            // outward normal (-ey, ex); check if all of B is strictly outside
            double refs[3];
            for (int k = 0; k < 3; ++k)
                refs[k] = -ey * (B[k].first - A[i].first)
                        + ex * (B[k].second - A[i].second);
            double selfSide = -ey * (A[(i + 2) % 3].first - A[i].first)
                            + ex * (A[(i + 2) % 3].second - A[i].second);
            bool allOpposite = true;
            for (int k = 0; k < 3; ++k)
                if (refs[k] * selfSide > -1e-12) { allOpposite = false; break; }
            if (allOpposite) return true;
        }
        return false;
    };
    (void)cross2;
    if (separated(P, Q)) return false;
    if (separated(Q, P)) return false;
    return true;
}

// ------------------------------------------------------- triangle vs triangle

// Moller's interval-overlap test, with the degenerate branches surfaced instead
// of collapsed into "no intersection".
Result2 tri_tri(const Tri& T1, const Tri& T2, double touch_tol) {
    const Vec3 p1[3] = {T1.a, T1.b, T1.c};
    const Vec3 p2[3] = {T2.a, T2.b, T2.c};

    // Bounding-box separation, first and cheapest.
    //
    // Without this, two triangles hundreds of voxels apart could still be
    // classified GRAZING whenever one triangle's vertex happened to lie near the
    // other's INFINITE plane -- which says nothing about whether they touch. The
    // symptom was that the grazing count tracked the broad-phase cell size
    // (24,013 at cell 20 against 266,337 at cell 80) while the transverse count
    // stayed put. Disjoint boxes cannot intersect, so this is both a speedup and
    // a correctness fix.
    {
        double lo1[3] = {std::min({p1[0].x, p1[1].x, p1[2].x}),
                         std::min({p1[0].y, p1[1].y, p1[2].y}),
                         std::min({p1[0].z, p1[1].z, p1[2].z})};
        double hi1[3] = {std::max({p1[0].x, p1[1].x, p1[2].x}),
                         std::max({p1[0].y, p1[1].y, p1[2].y}),
                         std::max({p1[0].z, p1[1].z, p1[2].z})};
        double lo2[3] = {std::min({p2[0].x, p2[1].x, p2[2].x}),
                         std::min({p2[0].y, p2[1].y, p2[2].y}),
                         std::min({p2[0].z, p2[1].z, p2[2].z})};
        double hi2[3] = {std::max({p2[0].x, p2[1].x, p2[2].x}),
                         std::max({p2[0].y, p2[1].y, p2[2].y}),
                         std::max({p2[0].z, p2[1].z, p2[2].z})};
        for (int k = 0; k < 3; ++k)
            if (hi1[k] < lo2[k] - touch_tol || hi2[k] < lo1[k] - touch_tol)
                return {NONE, 0.0, 0.0};
    }

    Vec3 n2 = cross(p2[1] - p2[0], p2[2] - p2[0]);
    double l2 = norm(n2);
    Vec3 n1 = cross(p1[1] - p1[0], p1[2] - p1[0]);
    double l1 = norm(n1);
    if (l1 < 1e-12 || l2 < 1e-12) return {GRAZING, 0.0, 0.0};   // degenerate
    n1 = n1 * (1.0 / l1);
    n2 = n2 * (1.0 / l2);

    double d1[3], d2[3];
    for (int i = 0; i < 3; ++i) d1[i] = dot(n2, p1[i] - p2[0]);
    for (int i = 0; i < 3; ++i) d2[i] = dot(n1, p2[i] - p1[0]);

    // Tolerance from the actual operand magnitudes, not a global constant.
    double mag = 0.0;
    for (const Vec3& p : {p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]})
        mag = std::max(mag, std::max(std::fabs(p.x),
                       std::max(std::fabs(p.y), std::fabs(p.z))));
    const double EPS = plane_eps(mag);

    bool graze = false;
    for (int i = 0; i < 3; ++i) {
        if (std::fabs(d1[i]) < EPS) graze = true;
        if (std::fabs(d2[i]) < EPS) graze = true;
    }

    // fully on one side of the other plane -> disjoint
    auto allSameSign = [EPS](const double d[3]) {
        return (d[0] > EPS && d[1] > EPS && d[2] > EPS)
            || (d[0] < -EPS && d[1] < -EPS && d[2] < -EPS);
    };
    if (allSameSign(d1) || allSameSign(d2)) return {NONE, 0.0, 0.0};

    const double ang = std::acos(std::min(1.0, std::fabs(dot(n1, n2))))
                     * 180.0 / 3.14159265358979323846;

    bool coplanar = std::fabs(d1[0]) < EPS && std::fabs(d1[1]) < EPS
                 && std::fabs(d1[2]) < EPS;
    if (coplanar)
        return tri2d_overlap(p1, p2, n1) ? Result2{COPLANAR, 0.0, ang}
                                         : Result2{NONE, 0.0, ang};
    if (graze) return {GRAZING, 0.0, ang};

    // intersection line direction, and the two intervals along it
    Vec3 D = cross(n1, n2);
    double ax = std::fabs(D.x), ay = std::fabs(D.y), az = std::fabs(D.z);
    int idx = (ay > ax && ay > az) ? 1 : ((az > ax && az > ay) ? 2 : 0);
    auto comp = [idx](const Vec3& p) { return idx == 0 ? p.x : (idx == 1 ? p.y : p.z); };

    // Compute the parameter interval of one triangle on the intersection line.
    // The vertex that sits alone on its side of the other plane is the apex.
    auto interval = [&](const Vec3 p[3], const double d[3], double& t0, double& t1) {
        int apex = -1;
        for (int i = 0; i < 3; ++i) {
            int j = (i + 1) % 3, k = (i + 2) % 3;
            if ((d[i] > 0 && d[j] <= 0 && d[k] <= 0)
             || (d[i] < 0 && d[j] >= 0 && d[k] >= 0)) { apex = i; break; }
        }
        if (apex < 0) return false;
        int j = (apex + 1) % 3, k = (apex + 2) % 3;
        double pa = comp(p[apex]);
        t0 = pa + (comp(p[j]) - pa) * (d[apex] / (d[apex] - d[j]));
        t1 = pa + (comp(p[k]) - pa) * (d[apex] / (d[apex] - d[k]));
        if (t0 > t1) std::swap(t0, t1);
        return true;
    };

    double a0, a1, b0, b1;
    if (!interval(p1, d1, a0, a1)) return {GRAZING, 0.0, ang};
    if (!interval(p2, d2, b0, b1)) return {GRAZING, 0.0, ang};

    // Strict positive overlap. Intervals [0,1] and [1,2] meet at a point, which
    // is a touch, not a transverse crossing.
    //
    // The intervals are parameterised by one coordinate of the intersection
    // line, so their overlap is a PROJECTION, shorter than the true segment by
    // the direction cosine. Dividing by it gives the Euclidean length of the
    // shared segment, which is what "penetration" should mean and what a reader
    // will assume it means. The factor lies between 1 and sqrt(3) because idx is
    // the largest component of D.
    const double dlen = norm(D);
    const double dcos = dlen > 1e-12
        ? std::fabs((idx == 0 ? D.x : (idx == 1 ? D.y : D.z)) / dlen) : 1.0;
    const double proj = std::min(a1, b1) - std::max(a0, b0);
    if (proj <= 0.0) return {NONE, 0.0, ang};
    const double pen = proj / std::max(dcos, 1e-9);   // Euclidean segment length
    if (pen <= touch_tol) return {GRAZING, pen, ang};
    return {TRANSVERSE, pen, ang};
}

// --------------------------------------------------------------------- IO

struct Surface {
    int32_t rows = 0, cols = 0;
    std::vector<float> pts;      // rows*cols*3
    std::vector<uint8_t> valid;  // rows*cols
    inline bool ok(int v, int u) const { return valid[(size_t)v * cols + u] != 0; }
    inline Vec3 at(int v, int u) const {
        size_t o = ((size_t)v * cols + u) * 3;
        return {pts[o], pts[o + 1], pts[o + 2]};
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
    if (rows == 0) { std::fprintf(stderr, "selfcross: OBJ soup not supported\n"); return false; }
    S.rows = (int)rows; S.cols = (int)cols;
    S.pts.resize((size_t)rows * cols * 3);
    S.valid.resize((size_t)rows * cols);
    f.read((char*)S.pts.data(), (std::streamsize)S.pts.size() * 4);
    f.read((char*)S.valid.data(), (std::streamsize)S.valid.size());
    return (bool)f;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
            "usage: %s <atlas.bin> <out.csv> [threads] [cell] [exclude] [diagonal]\n",
            argv[0]);
        return 2;
    }
    const char* in = argv[1];
    const char* out = argv[2];
    int nthreads = argc > 3 ? std::atoi(argv[3]) : 0;
    double cell = argc > 4 ? std::atof(argv[4]) : 40.0;
    int exclude = argc > 5 ? std::atoi(argv[5]) : 1;
    int diagonal = argc > 6 ? std::atoi(argv[6]) : 0;
    double maxedge = argc > 7 ? std::atof(argv[7]) : 60.0;
    // Touch tolerance on the interval overlap. float32 coordinates near 5,000
    // have a spacing of ~5e-4 voxels, so anything finer than that is below the
    // input's own resolution and cannot be a meaningful margin.
    double touch_tol = argc > 8 ? std::atof(argv[8]) : 1e-3;
    if (nthreads <= 0) nthreads = (int)std::thread::hardware_concurrency();
    if (nthreads <= 0) nthreads = 4;

    Surface S;
    if (!read_atlas(in, S)) { std::fprintf(stderr, "selfcross: bad atlas %s\n", in); return 1; }

    // ------------------------------------------------------ build triangles
    std::vector<Tri> tris;
    size_t dropped = 0;
    tris.reserve((size_t)S.rows * S.cols);
    for (int v = 0; v + 1 < S.rows; ++v) {
        for (int u = 0; u + 1 < S.cols; ++u) {
            if (!S.ok(v, u) || !S.ok(v + 1, u) || !S.ok(v, u + 1) || !S.ok(v + 1, u + 1))
                continue;
            Vec3 p00 = S.at(v, u), p10 = S.at(v + 1, u);
            Vec3 p01 = S.at(v, u + 1), p11 = S.at(v + 1, u + 1);
            if (maxedge > 0) {
                double e = std::max(std::max(norm(p01 - p00), norm(p11 - p01)),
                          std::max(std::max(norm(p10 - p11), norm(p00 - p10)),
                                   std::max(norm(p11 - p00), norm(p10 - p01))));
                if (e > maxedge) { ++dropped; continue; }
            }
            if (diagonal == 0) {
                tris.push_back({p00, p01, p11, v, u, 0});
                tris.push_back({p00, p11, p10, v, u, 1});
            } else {
                tris.push_back({p00, p01, p10, v, u, 0});
                tris.push_back({p01, p11, p10, v, u, 1});
            }
        }
    }
    std::fprintf(stderr, "selfcross: grid %dx%d, %zu triangles (%zu quads dropped "
                 "for edge > %.0f), cell %.1f, exclude %d, diagonal %d, %d threads\n",
                 S.rows, S.cols, tris.size(), dropped, maxedge, cell, exclude,
                 diagonal, nthreads);
    if (tris.empty()) { std::FILE* fp = std::fopen(out, "w");
        if (fp) { std::fprintf(fp,
            "v1,u1,v2,u2,verdict,penetration,angle_deg,tri1,tri2\n");
            std::fclose(fp); }
        std::fprintf(stderr, "selfcross: no triangles\n"); return 0; }

    // ------------------------------------------------------- uniform grid
    Vec3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
    for (const Tri& t : tris) {
        for (const Vec3& p : {t.a, t.b, t.c}) {
            lo.x = std::min(lo.x, p.x); hi.x = std::max(hi.x, p.x);
            lo.y = std::min(lo.y, p.y); hi.y = std::max(hi.y, p.y);
            lo.z = std::min(lo.z, p.z); hi.z = std::max(hi.z, p.z);
        }
    }
    const int nx = std::max(1, (int)((hi.x - lo.x) / cell) + 1);
    const int ny = std::max(1, (int)((hi.y - lo.y) / cell) + 1);
    const int nz = std::max(1, (int)((hi.z - lo.z) / cell) + 1);
    auto cellIdx = [&](int i, int j, int k) {
        return ((size_t)k * ny + j) * nx + i;
    };
    std::unordered_map<size_t, std::vector<uint32_t>> buckets;
    buckets.reserve(tris.size());
    for (uint32_t ti = 0; ti < tris.size(); ++ti) {
        const Tri& t = tris[ti];
        double tx0 = std::min({t.a.x, t.b.x, t.c.x}), tx1 = std::max({t.a.x, t.b.x, t.c.x});
        double ty0 = std::min({t.a.y, t.b.y, t.c.y}), ty1 = std::max({t.a.y, t.b.y, t.c.y});
        double tz0 = std::min({t.a.z, t.b.z, t.c.z}), tz1 = std::max({t.a.z, t.b.z, t.c.z});
        int i0 = (int)((tx0 - lo.x) / cell), i1 = (int)((tx1 - lo.x) / cell);
        int j0 = (int)((ty0 - lo.y) / cell), j1 = (int)((ty1 - lo.y) / cell);
        int k0 = (int)((tz0 - lo.z) / cell), k1 = (int)((tz1 - lo.z) / cell);
        for (int k = k0; k <= k1; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i)
                    buckets[cellIdx(i, j, k)].push_back(ti);
    }
    std::vector<std::pair<size_t, const std::vector<uint32_t>*>> cells;
    cells.reserve(buckets.size());
    for (auto& kv : buckets)
        if (kv.second.size() > 1) cells.emplace_back(kv.first, &kv.second);
    std::sort(cells.begin(), cells.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::fprintf(stderr, "selfcross: %zu occupied cells with >1 triangle\n", cells.size());

    // ---------------------------------------------------------- narrow phase
    struct Hit { int32_t v1, u1, v2, u2; int verdict; float pen, ang;
                 int32_t t1, t2; };
    std::vector<std::vector<Hit>> perThread(nthreads);
    std::atomic<size_t> next{0};
    std::atomic<size_t> tested{0};

    // DETERMINISTIC OWNERSHIP.
    //
    // A triangle pair can share several grid cells. The previous version kept a
    // per-thread `seen` set, so such a pair was tested once per thread that
    // happened to receive one of those cells: the same trace returned 11,211
    // crossings on one thread and 21,213 on twelve, and the number moved with
    // scheduling. Every count it produced was meaningless.
    //
    // Instead, each pair is owned by exactly one cell -- the one containing the
    // minimum corner of the two AABBs' overlap -- and is tested only there. No
    // shared state, no locks, and the answer no longer depends on thread count
    // or on the order cells are handed out.
    auto ownerCell = [&](const Tri& A, const Tri& B) {
        double ox = std::max(std::min({A.a.x, A.b.x, A.c.x}),
                             std::min({B.a.x, B.b.x, B.c.x}));
        double oy = std::max(std::min({A.a.y, A.b.y, A.c.y}),
                             std::min({B.a.y, B.b.y, B.c.y}));
        double oz = std::max(std::min({A.a.z, A.b.z, A.c.z}),
                             std::min({B.a.z, B.b.z, B.c.z}));
        int i = (int)((ox - lo.x) / cell);
        int j = (int)((oy - lo.y) / cell);
        int k = (int)((oz - lo.z) / cell);
        return cellIdx(i, j, k);
    };

    auto worker = [&](int id) {
        auto& outv = perThread[id];
        size_t local = 0;
        for (;;) {
            size_t ci = next.fetch_add(1);
            if (ci >= cells.size()) break;
            const size_t key = cells[ci].first;
            const std::vector<uint32_t>& c = *cells[ci].second;
            for (size_t a = 0; a < c.size(); ++a) {
                for (size_t b = a + 1; b < c.size(); ++b) {
                    const Tri& T1 = tris[c[a]];
                    const Tri& T2 = tris[c[b]];
                    // Adjacency. Quads within Chebyshev distance 1 in the grid
                    // share at least one vertex, so exclude = 1 is exactly
                    // shared-vertex/shared-edge exclusion rather than an
                    // arbitrary window.
                    if (std::abs(T1.v - T2.v) <= exclude
                        && std::abs(T1.u - T2.u) <= exclude) continue;
                    if (ownerCell(T1, T2) != key) continue;   // tested elsewhere
                    ++local;
                    Result2 r = tri_tri(T1, T2, touch_tol);
                    if (r.verdict != NONE) {
                        // Canonical endpoint order in the ENGINE (round 19):
                        // the lexicographically smaller (v,u,t) is side 1, so
                        // every reader sees one identity per geometric pair.
                        bool swap = std::make_tuple(T1.v, T1.u, T1.t)
                                  > std::make_tuple(T2.v, T2.u, T2.t);
                        const Tri& A_ = swap ? T2 : T1;
                        const Tri& B_ = swap ? T1 : T2;
                        outv.push_back({A_.v, A_.u, B_.v, B_.u, (int)r.verdict,
                                        (float)r.penetration,
                                        (float)r.angle_deg, A_.t, B_.t});
                    }
                }
            }
        }
        tested.fetch_add(local);
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < nthreads; ++t) pool.emplace_back(worker, t);
    for (auto& t : pool) t.join();

    size_t nT = 0, nC = 0, nG = 0;
    std::FILE* fp = std::fopen(out, "w");
    if (!fp) { std::fprintf(stderr, "selfcross: cannot write %s\n", out); return 1; }
    // schema v2: +tri1,tri2 (local triangle index per side). Legacy readers
    // that index columns 0-6 are unaffected; strict transactional acceptance
    // requires these columns (round 19).
    std::fprintf(fp, "v1,u1,v2,u2,verdict,penetration,angle_deg,tri1,tri2\n");
    // Rows were written in per-thread completion order, so identical meshes
    // produced identical row SETS in scheduling-dependent order -- and
    // different CSV bytes/hashes, and order-dependent downstream tie-breaks.
    // (v1,u1,tri1,v2,u2,tri2) is unique per tested pair: a total order.
    std::vector<Hit> allHits;
    { size_t tot = 0;
      for (auto& vec : perThread) tot += vec.size();
      allHits.reserve(tot);
      for (auto& vec : perThread)
          allHits.insert(allHits.end(), vec.begin(), vec.end()); }
    std::sort(allHits.begin(), allHits.end(),
              [](const Hit& x, const Hit& y) {
        return std::tie(x.v1, x.u1, x.t1, x.v2, x.u2, x.t2)
             < std::tie(y.v1, y.u1, y.t1, y.v2, y.u2, y.t2);
    });
    for (const Hit& h : allHits) {
            const char* name = h.verdict == TRANSVERSE ? "transverse"
                             : h.verdict == COPLANAR ? "coplanar" : "grazing";
            if (h.verdict == TRANSVERSE) ++nT;
            else if (h.verdict == COPLANAR) ++nC;
            else ++nG;
            std::fprintf(fp, "%d,%d,%d,%d,%s,%.6g,%.2f,%d,%d\n",
                         h.v1, h.u1, h.v2, h.u2, name, h.pen, h.ang,
                         h.t1, h.t2);
        }
    std::fclose(fp);

    std::fprintf(stderr,
        "selfcross: pairs tested %zu | transverse %zu | coplanar %zu | grazing %zu\n",
        tested.load(), nT, nC, nG);
    std::printf("{\"triangles\":%zu,\"quads_dropped\":%zu,\"pairs_tested\":%zu,"
                "\"transverse\":%zu,\"coplanar\":%zu,\"grazing\":%zu,"
                "\"exclude\":%d,\"diagonal\":%d,\"maxedge\":%.1f}\n",
                tris.size(), dropped, tested.load(), nT, nC, nG, exclude,
                diagonal, maxedge);
    return 0;
}
