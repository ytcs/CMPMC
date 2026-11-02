#pragma once

#include "cmpmc/core/core.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

// Electric fields interpolated on finite-element meshes: piecewise-linear
// potentials on 3D tetrahedral meshes (TriLinearMesh3D) and on 2D triangular
// meshes, planar or axisymmetric (BiLinearMesh2D). Positions in metres,
// potentials in volts, fields in V/m.

namespace cmpmc::geometry {

// ============================================================================
// TriLinearMesh3D: 3D tetrahedral FEM mesh, electric field and potential
//
//    Physics & Numerical Method:
//    - Given a 3D tetrahedral mesh with nodal coordinates X[n] (m) and
//      electrostatic potentials V[n] (V):
//    - Each tetrahedron t = (v0, v1, v2, v3) has edge matrix
//        T = [ X[v0] - X[v3],  X[v1] - X[v3],  X[v2] - X[v3] ]
//      We precompute T^{-1} and the constant interior electric field
//        E_tet = -Grad(V) = -(T^{-1})^T * [ V[v0]-V[v3], V[v1]-V[v3], V[v2]-V[v3] ]^T
//      at mesh construction time (eliminating per-step matrix inversions).
//    - Point location uses the barycentric facet-walking algorithm:
//        [lam0, lam1, lam2]^T = T^{-1} * (p - X[v3]),  lam3 = 1 - lam0 - lam1 - lam2
//      If min(lam_i) >= -eps, `p` is inside tetrahedron `t`; otherwise we step
//      across the facet opposite to argmin(lam_i) via `neighbors_[t][argmin]`.
// ============================================================================
class TriLinearMesh3D {
public:
    using Tet = std::array<int, 4>;

    TriLinearMesh3D() = default;

    TriLinearMesh3D(std::vector<Vec3> vertices,
                    std::vector<double> potentials,
                    std::vector<Tet> tetrahedra,
                    double V_min_clamp = -1e30,
                    double V_max_clamp =  1e30)
        : X_(std::move(vertices)),
          V_(std::move(potentials)),
          Tetra_(std::move(tetrahedra)) {
        for (double& v : V_) {
            v = std::clamp(v, V_min_clamp, V_max_clamp);
        }
        fill_neighbors();
        fill_tinv_and_grad();
    }

    bool empty() const noexcept { return Tetra_.empty(); }
    size_t num_vertices() const noexcept { return X_.size(); }
    size_t num_tetrahedra() const noexcept { return Tetra_.size(); }

    // Walk tetrahedral mesh starting from `tet_hint` (updated in-place).
    // Returns true if `p` is inside the mesh, false if walk hits an exterior hull facet.
    bool find_tetrahedron(const Vec3& p, int& tet_hint, std::array<double, 4>& bary) const noexcept {
        if (Tetra_.empty()) return false;
        if (tet_hint < 0 || tet_hint >= static_cast<int>(Tetra_.size())) {
            tet_hint = static_cast<int>(Tetra_.size() / 2);
        }

        constexpr double kCartTolerance = -1e-10;
        const int max_steps = static_cast<int>(Tetra_.size()) + 4;

        for (int step = 0; step < max_steps; ++step) {
            bary_coord(p, tet_hint, bary);

            double min_b = bary[0];
            int min_idx  = 0;
            for (int i = 1; i < 4; ++i) {
                if (bary[i] < min_b) {
                    min_b   = bary[i];
                    min_idx = i;
                }
            }

            if (min_b >= kCartTolerance) {
                return true;
            }

            const int next_tet = Neighbors_[tet_hint][min_idx];
            if (next_tet < 0) {
                return false; // Reached convex/exterior boundary of mesh
            }
            tet_hint = next_tet;
        }
        return false;
    }

    // Interpolated electrostatic potential V(p) [Volts]
    double potential_at(const Vec3& p, int& tet_hint) const noexcept {
        std::array<double, 4> bary{};
        if (!find_tetrahedron(p, tet_hint, bary)) return 0.0;
        const Tet& tet = Tetra_[tet_hint];
        return bary[0] * V_[tet[0]] +
               bary[1] * V_[tet[1]] +
               bary[2] * V_[tet[2]] +
               bary[3] * V_[tet[3]];
    }

    // Electric field E(p) = -Grad V(p) [V/m]
    Vec3 field_at(const Vec3& p, int& tet_hint) const noexcept {
        std::array<double, 4> bary{};
        if (!find_tetrahedron(p, tet_hint, bary)) return {0.0, 0.0, 0.0};
        return Grad_[tet_hint];
    }

    // Stateless convenience overload (uses thread_local tetrahedron cache for O(1) tracking)
    Vec3 at(const Vec3& p) const noexcept {
        thread_local int cached_tet = 0;
        return field_at(p, cached_tet);
    }

    double potential_at(const Vec3& p) const noexcept {
        thread_local int cached_tet = 0;
        return potential_at(p, cached_tet);
    }

    // Helper: Build a 3D Kuhn-triangulated Cartesian box mesh (6 tetrahedra per hex cell)
    // over [-half_x, +half_x] x [-half_y, +half_y] x [-half_z, +half_z]
    template <typename PotentialFn>
    static TriLinearMesh3D build_cartesian_grid(double half_x, double half_y, double half_z,
                                                int nx, int ny, int nz,
                                                PotentialFn&& potential_fn) {
        nx = std::max(1, nx);
        ny = std::max(1, ny);
        nz = std::max(1, nz);

        const int vx = nx + 1;
        const int vy = ny + 1;
        const int vz = nz + 1;
        auto v_idx = [vy, vz](int i, int j, int k) noexcept {
            return (i * vy + j) * vz + k;
        };

        std::vector<Vec3> verts(vx * vy * vz);
        std::vector<double> pots(vx * vy * vz);

        for (int i = 0; i <= nx; ++i) {
            const double x = -half_x + 2.0 * half_x * (static_cast<double>(i) / nx);
            for (int j = 0; j <= ny; ++j) {
                const double y = -half_y + 2.0 * half_y * (static_cast<double>(j) / ny);
                for (int k = 0; k <= nz; ++k) {
                    const double z = -half_z + 2.0 * half_z * (static_cast<double>(k) / nz);
                    const int id = v_idx(i, j, k);
                    verts[id] = Vec3{x, y, z};
                    pots[id]  = potential_fn(Vec3{x, y, z});
                }
            }
        }

        // Subdivide each hexahedral voxel into 6 Kuhn tetrahedra sharing the main diagonal (v000 -> v111)
        std::vector<Tet> tets;
        tets.reserve(static_cast<size_t>(nx * ny * nz * 6));

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    const int v000 = v_idx(i,     j,     k);
                    const int v100 = v_idx(i + 1, j,     k);
                    const int v010 = v_idx(i,     j + 1, k);
                    const int v110 = v_idx(i + 1, j + 1, k);
                    const int v001 = v_idx(i,     j,     k + 1);
                    const int v101 = v_idx(i + 1, j,     k + 1);
                    const int v011 = v_idx(i,     j + 1, k + 1);
                    const int v111 = v_idx(i + 1, j + 1, k + 1);

                    tets.push_back({v000, v100, v110, v111});
                    tets.push_back({v000, v100, v101, v111});
                    tets.push_back({v000, v010, v110, v111});
                    tets.push_back({v000, v010, v011, v111});
                    tets.push_back({v000, v001, v101, v111});
                    tets.push_back({v000, v001, v011, v111});
                }
            }
        }

        return TriLinearMesh3D(std::move(verts), std::move(pots), std::move(tets));
    }

    // Load a mesh from two whitespace-separated ASCII files (the layout also used by
    // G4CMP mesh fields): `coords_file` has one node per line "x y z V";
    // `tets_file` has one tetrahedron per line "i0 i1 i2 i3" (0-based node indices).
    // Lines starting with '#' or '%' are skipped. Coordinates and potentials are
    // multiplied by coord_unit_m and potential_unit_V.
    static TriLinearMesh3D load_from_files(const std::string& coords_file,
                                           const std::string& tets_file,
                                           double coord_unit_m = 1.0,
                                           double potential_unit_V = 1.0) {
        std::vector<Vec3> verts;
        std::vector<double> pots;
        {
            std::ifstream in(coords_file);
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#' || line[0] == '%') continue;
                std::istringstream iss(line);
                double x, y, z, v;
                if (iss >> x >> y >> z >> v) {
                    verts.push_back(Vec3{x * coord_unit_m, y * coord_unit_m, z * coord_unit_m});
                    pots.push_back(v * potential_unit_V);
                }
            }
        }
        std::vector<Tet> tets;
        {
            std::ifstream in(tets_file);
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#' || line[0] == '%') continue;
                std::istringstream iss(line);
                int a, b, c, d;
                if (iss >> a >> b >> c >> d) {
                    Tet t{a, b, c, d};
                    std::sort(t.begin(), t.end());
                    tets.push_back(t);
                }
            }
        }
        return TriLinearMesh3D(std::move(verts), std::move(pots), std::move(tets));
    }

private:
    std::vector<Vec3> X_;
    std::vector<double> V_;
    std::vector<Tet> Tetra_;
    std::vector<Tet> Neighbors_;
    std::vector<Mat3> Tinv_;
    std::vector<Vec3> Grad_; // Stored as E = -Grad(V) [V/m]

    inline void bary_coord(const Vec3& p, int tet_idx, std::array<double, 4>& bary) const noexcept {
        const Vec3 dx = p - X_[Tetra_[tet_idx][3]];
        const Vec3 lam012 = Tinv_[tet_idx] * dx;
        bary[0] = lam012.x;
        bary[1] = lam012.y;
        bary[2] = lam012.z;
        bary[3] = 1.0 - lam012.x - lam012.y - lam012.z;
    }

    struct FacetEntry {
        std::array<int, 3> nodes;
        int tet_id;
        int opp_slot;

        bool operator<(const FacetEntry& o) const noexcept {
            return nodes < o.nodes;
        }
    };

    void fill_neighbors() {
        const size_t n_tets = Tetra_.size();
        Neighbors_.assign(n_tets, {-1, -1, -1, -1});

        for (Tet& t : Tetra_) {
            std::sort(t.begin(), t.end());
        }

        std::vector<FacetEntry> facets;
        facets.reserve(n_tets * 4);
        for (size_t it = 0; it < n_tets; ++it) {
            const Tet& t = Tetra_[it];
            facets.push_back({{t[1], t[2], t[3]}, static_cast<int>(it), 0});
            facets.push_back({{t[0], t[2], t[3]}, static_cast<int>(it), 1});
            facets.push_back({{t[0], t[1], t[3]}, static_cast<int>(it), 2});
            facets.push_back({{t[0], t[1], t[2]}, static_cast<int>(it), 3});
        }

        std::sort(facets.begin(), facets.end());

        for (size_t i = 0; i + 1 < facets.size(); ++i) {
            if (facets[i].nodes == facets[i + 1].nodes) {
                Neighbors_[facets[i].tet_id][facets[i].opp_slot]         = facets[i + 1].tet_id;
                Neighbors_[facets[i + 1].tet_id][facets[i + 1].opp_slot] = facets[i].tet_id;
                ++i;
            }
        }
    }

    void fill_tinv_and_grad() {
        const size_t n_tets = Tetra_.size();
        Tinv_.resize(n_tets);
        Grad_.resize(n_tets);

        for (size_t it = 0; it < n_tets; ++it) {
            const Tet& tet = Tetra_[it];
            const Vec3& x3 = X_[tet[3]];
            const Vec3 c0  = X_[tet[0]] - x3;
            const Vec3 c1  = X_[tet[1]] - x3;
            const Vec3 c2  = X_[tet[2]] - x3;

            // Matrix T has columns c0, c1, c2
            const double a = c0.x, b = c1.x, c = c2.x;
            const double d = c0.y, e = c1.y, f = c2.y;
            const double g = c0.z, h = c1.z, k = c2.z;

            const double det = a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
            if (std::abs(det) <= 1e-30) {
                Tinv_[it] = Mat3{};
                Grad_[it] = {0.0, 0.0, 0.0};
                continue;
            }
            const double inv_det = 1.0 / det;

            Mat3 inv;
            inv.m[0][0] =  (e * k - f * h) * inv_det;
            inv.m[0][1] = -(b * k - c * h) * inv_det;
            inv.m[0][2] =  (b * f - c * e) * inv_det;
            inv.m[1][0] = -(d * k - f * g) * inv_det;
            inv.m[1][1] =  (a * k - c * g) * inv_det;
            inv.m[1][2] = -(a * f - c * d) * inv_det;
            inv.m[2][0] =  (d * h - e * g) * inv_det;
            inv.m[2][1] = -(a * h - b * g) * inv_det;
            inv.m[2][2] =  (a * e - b * d) * inv_det;

            Tinv_[it] = inv;

            // E = -Grad(V) = -(T^{-1})^T * dV
            const Vec3 dV{
                V_[tet[0]] - V_[tet[3]],
                V_[tet[1]] - V_[tet[3]],
                V_[tet[2]] - V_[tet[3]]
            };
            Grad_[it] = inv.mul_transpose(dV) * (-1.0);
        }
    }
};

// ============================================================================
// BiLinearMesh2D: 2D triangular FEM mesh, electric field and potential
//
//    Supports either:
//    - Axisymmetric cylindrical (r, z) meshes (`axisymmetric_rz = true`), where
//      3D position (x, y, z) is projected to (r = hypot(x, y), z) and the radial
//      electric field E_r is rotated back into (E_x, E_y, E_z)
//    - Planar 2D (x, y) meshes (`axisymmetric_rz = false`).
// ============================================================================
class BiLinearMesh2D {
public:
    using Tri = std::array<int, 3>;

    struct Point2D {
        double u = 0.0; // r (or x)
        double v = 0.0; // z (or y)
    };

    BiLinearMesh2D() = default;

    BiLinearMesh2D(std::vector<Point2D> vertices,
                   std::vector<double> potentials,
                   std::vector<Tri> triangles,
                   bool axisymmetric_rz = true)
        : X_(std::move(vertices)),
          V_(std::move(potentials)),
          Tri_(std::move(triangles)),
          axisymmetric_rz_(axisymmetric_rz) {
        fill_neighbors();
        fill_tinv_and_grad();
    }

    bool empty() const noexcept { return Tri_.empty(); }
    size_t num_vertices() const noexcept { return X_.size(); }
    size_t num_triangles() const noexcept { return Tri_.size(); }

    bool find_triangle(Point2D p, int& tri_hint, std::array<double, 3>& bary) const noexcept {
        if (Tri_.empty()) return false;
        if (tri_hint < 0 || tri_hint >= static_cast<int>(Tri_.size())) {
            tri_hint = static_cast<int>(Tri_.size() / 2);
        }

        constexpr double kCartTolerance = -1e-10;
        const int max_steps = static_cast<int>(Tri_.size()) + 4;

        for (int step = 0; step < max_steps; ++step) {
            bary_coord(p, tri_hint, bary);

            double min_b = bary[0];
            int min_idx  = 0;
            for (int i = 1; i < 3; ++i) {
                if (bary[i] < min_b) {
                    min_b   = bary[i];
                    min_idx = i;
                }
            }

            if (min_b >= kCartTolerance) {
                return true;
            }

            const int next_tri = Neighbors_[tri_hint][min_idx];
            if (next_tri < 0) return false;
            tri_hint = next_tri;
        }
        return false;
    }

    double potential_at(const Vec3& p3d, int& tri_hint) const noexcept {
        const Point2D p2d = project(p3d);
        std::array<double, 3> bary{};
        if (!find_triangle(p2d, tri_hint, bary)) return 0.0;
        const Tri& t = Tri_[tri_hint];
        return bary[0] * V_[t[0]] + bary[1] * V_[t[1]] + bary[2] * V_[t[2]];
    }

    Vec3 field_at(const Vec3& p3d, int& tri_hint) const noexcept {
        const Point2D p2d = project(p3d);
        std::array<double, 3> bary{};
        if (!find_triangle(p2d, tri_hint, bary)) return {0.0, 0.0, 0.0};

        const Point2D g = Grad_[tri_hint]; // (E_u, E_v) = (-dV/du, -dV/dv)
        if (!axisymmetric_rz_) {
            return {g.u, g.v, 0.0};
        }
        const double r = p2d.u;
        if (r > 1e-15) {
            return {g.u * (p3d.x / r), g.u * (p3d.y / r), g.v};
        }
        return {0.0, 0.0, g.v};
    }

    Vec3 at(const Vec3& p3d) const noexcept {
        thread_local int cached_tri = 0;
        return field_at(p3d, cached_tri);
    }

    double potential_at(const Vec3& p3d) const noexcept {
        thread_local int cached_tri = 0;
        return potential_at(p3d, cached_tri);
    }

    // Helper: Build a 2D triangulated grid on [u_min, u_max] x [v_min, v_max]
    template <typename PotentialFn2D>
    static BiLinearMesh2D build_grid(double u_min, double u_max,
                                     double v_min, double v_max,
                                     int nu, int nv,
                                     PotentialFn2D&& potential_fn,
                                     bool axisymmetric_rz = true) {
        nu = std::max(1, nu);
        nv = std::max(1, nv);
        const int vu = nu + 1;
        const int vv = nv + 1;
        auto v_idx = [vv](int i, int j) noexcept { return i * vv + j; };

        std::vector<Point2D> verts(vu * vv);
        std::vector<double> pots(vu * vv);

        for (int i = 0; i <= nu; ++i) {
            const double u = u_min + (u_max - u_min) * (static_cast<double>(i) / nu);
            for (int j = 0; j <= nv; ++j) {
                const double v = v_min + (v_max - v_min) * (static_cast<double>(j) / nv);
                const int id = v_idx(i, j);
                verts[id] = {u, v};
                pots[id]  = potential_fn(u, v);
            }
        }

        std::vector<Tri> tris;
        tris.reserve(static_cast<size_t>(nu * nv * 2));
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nv; ++j) {
                const int v00 = v_idx(i,     j);
                const int v10 = v_idx(i + 1, j);
                const int v01 = v_idx(i,     j + 1);
                const int v11 = v_idx(i + 1, j + 1);
                tris.push_back({v00, v10, v11});
                tris.push_back({v00, v01, v11});
            }
        }

        return BiLinearMesh2D(std::move(verts), std::move(pots), std::move(tris), axisymmetric_rz);
    }

private:
    std::vector<Point2D> X_;
    std::vector<double> V_;
    std::vector<Tri> Tri_;
    std::vector<Tri> Neighbors_;
    std::vector<std::array<double, 4>> Tinv_; // 2x2 row-major inverse
    std::vector<Point2D> Grad_;               // E = -Grad(V)
    bool axisymmetric_rz_ = true;

    inline Point2D project(const Vec3& p) const noexcept {
        if (axisymmetric_rz_) {
            return {std::hypot(p.x, p.y), p.z};
        }
        return {p.x, p.y};
    }

    inline void bary_coord(Point2D p, int tri_idx, std::array<double, 3>& bary) const noexcept {
        const Point2D& x2 = X_[Tri_[tri_idx][2]];
        const double du = p.u - x2.u;
        const double dv = p.v - x2.v;
        const auto& inv = Tinv_[tri_idx];
        bary[0] = inv[0] * du + inv[1] * dv;
        bary[1] = inv[2] * du + inv[3] * dv;
        bary[2] = 1.0 - bary[0] - bary[1];
    }

    struct EdgeEntry {
        std::array<int, 2> nodes;
        int tri_id;
        int opp_slot;

        bool operator<(const EdgeEntry& o) const noexcept {
            return nodes < o.nodes;
        }
    };

    void fill_neighbors() {
        const size_t n_tris = Tri_.size();
        Neighbors_.assign(n_tris, {-1, -1, -1});

        for (Tri& t : Tri_) {
            std::sort(t.begin(), t.end());
        }

        std::vector<EdgeEntry> edges;
        edges.reserve(n_tris * 3);
        for (size_t it = 0; it < n_tris; ++it) {
            const Tri& t = Tri_[it];
            edges.push_back({{t[1], t[2]}, static_cast<int>(it), 0});
            edges.push_back({{t[0], t[2]}, static_cast<int>(it), 1});
            edges.push_back({{t[0], t[1]}, static_cast<int>(it), 2});
        }

        std::sort(edges.begin(), edges.end());
        for (size_t i = 0; i + 1 < edges.size(); ++i) {
            if (edges[i].nodes == edges[i + 1].nodes) {
                Neighbors_[edges[i].tri_id][edges[i].opp_slot]         = edges[i + 1].tri_id;
                Neighbors_[edges[i + 1].tri_id][edges[i + 1].opp_slot] = edges[i].tri_id;
                ++i;
            }
        }
    }

    void fill_tinv_and_grad() {
        const size_t n_tris = Tri_.size();
        Tinv_.resize(n_tris);
        Grad_.resize(n_tris);

        for (size_t it = 0; it < n_tris; ++it) {
            const Tri& t = Tri_[it];
            const Point2D& x2 = X_[t[2]];
            const double a = X_[t[0]].u - x2.u;
            const double b = X_[t[1]].u - x2.u;
            const double c = X_[t[0]].v - x2.v;
            const double d = X_[t[1]].v - x2.v;

            const double det = a * d - b * c;
            if (std::abs(det) <= 1e-30) {
                Tinv_[it] = {0.0, 0.0, 0.0, 0.0};
                Grad_[it] = {0.0, 0.0};
                continue;
            }
            const double inv_det = 1.0 / det;
            const std::array<double, 4> inv = {
                 d * inv_det, -b * inv_det,
                -c * inv_det,  a * inv_det
            };
            Tinv_[it] = inv;

            const double dv0 = V_[t[0]] - V_[t[2]];
            const double dv1 = V_[t[1]] - V_[t[2]];
            Grad_[it] = {
                -(inv[0] * dv0 + inv[2] * dv1),
                -(inv[1] * dv0 + inv[3] * dv1)
            };
        }
    }
};

} // namespace cmpmc::geometry
