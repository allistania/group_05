#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <array>
#include <limits>
#include <tuple>
#include <cstdio>      // for sprintf
#include <mpi.h>

using namespace std;

const double GAMMA = 1.4;
const double CFL = 0.8;
const double FINAL_TIME = 1.0;
const double OUTPUT_INTERVAL = 0.1;

const double INLET_RHO = 1.0;
const double INLET_U   = 1.0;
const double INLET_V   = 0.0;
const double INLET_P   = 1.0;

const double OUTLET_P = 1.0;

int mpi_rank, mpi_size;

double WaveFunk(double p, double pk, double rho_k) {
    if (p > pk) {
        double Ak = 2.0 / ((GAMMA + 1) * rho_k);
        double Bk = (GAMMA - 1) / (GAMMA + 1) * pk;
        return (p - pk) * sqrt(Ak / (p + Bk));
    } else {
        return 2.0 * sqrt(GAMMA * pk / rho_k) / (GAMMA - 1) *
               (pow(p / pk, (GAMMA - 1) / (2 * GAMMA)) - 1.0);
    }
}

double ProisvWaveFunk(double p, double pk, double rho_k) {
    if (p > pk) {
        double Ak = 2.0 / ((GAMMA + 1) * rho_k);
        double Bk = (GAMMA - 1) / (GAMMA + 1) * pk;
        return sqrt(Ak / (p + Bk)) * (1.0 - (p - pk) / (2.0 * (p + Bk)));
    } else {
        return (1.0 / (rho_k * sqrt(GAMMA * pk / rho_k))) *
               pow(p / pk, -(GAMMA + 1) / (2 * GAMMA));
    }
}

double solve_p_star(double rhoL, double uL, double pL,
                    double rhoR, double uR, double pR) {
    double tol = 1e-8;
    int max_iter = 100;
    double p0 = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rhoL + rhoR) *
                (sqrt(GAMMA * pL / rhoL) + sqrt(GAMMA * pR / rhoR));
    p0 = max(tol, p0);
    double p_star = p0;

    for (int i = 0; i < max_iter; i++) {
        double fL = WaveFunk(p_star, pL, rhoL);
        double fR = WaveFunk(p_star, pR, rhoR);
        double f = fL + fR + uR - uL;

        double dfL = ProisvWaveFunk(p_star, pL, rhoL);
        double dfR = ProisvWaveFunk(p_star, pR, rhoR);
        double df = dfL + dfR;

        if (abs(df) < 1e-12) break;

        double dp = -f / df;
        p_star += dp;

        if (p_star < tol) p_star = tol;
        if (abs(dp) < tol * p_star) break;
    }
    return p_star;
}

double u_star(double p_star, double rhoL, double uL, double pL,
              double rhoR, double uR, double pR) {
    double fL = WaveFunk(p_star, pL, rhoL);
    double fR = WaveFunk(p_star, pR, rhoR);
    return 0.5 * (uL + uR) + 0.5 * (fR - fL);
}

tuple<double, double, double> sample(
    double p_star, double u_star_val, double rhoL, double uL, double pL,
    double rhoR, double uR, double pR, double x_over_t) {

    double cL = sqrt(GAMMA * pL / rhoL);
    double cR = sqrt(GAMMA * pR / rhoR);

    if (x_over_t < u_star_val) {
        if (p_star > pL) {
            double SL = uL - sqrt(((GAMMA + 1) / (2 * GAMMA)) * (p_star / pL) +
                                 (GAMMA - 1) / (2 * GAMMA)) * cL;
            if (x_over_t < SL)
                return {rhoL, uL, pL};
            else {
                double rho_star = rhoL * ((p_star / pL) + (GAMMA - 1) / (GAMMA + 1)) /
                                 (((GAMMA - 1) / (GAMMA + 1)) * (p_star / pL) + 1);
                return {rho_star, u_star_val, p_star};
            }
        } else {
            double SHL = uL - cL;
            double STL = u_star_val - sqrt(GAMMA * p_star /
                                          (rhoL * pow(p_star / pL, 1.0 / GAMMA)));
            if (x_over_t < SHL)
                return {rhoL, uL, pL};
            else if (x_over_t < STL) {
                double u = 2.0 / (GAMMA + 1) * (cL + (GAMMA - 1) / 2 * uL + x_over_t);
                double c = 2.0 / (GAMMA + 1) * (cL + (GAMMA - 1) / 2 * (uL - x_over_t));
                double rho = rhoL * pow(c / cL, 2.0 / (GAMMA - 1));
                double p = pL * pow(rho / rhoL, GAMMA);
                return {rho, u, p};
            } else {
                double rho_star = rhoL * pow(p_star / pL, 1.0 / GAMMA);
                return {rho_star, u_star_val, p_star};
            }
        }
    } else {
        if (p_star > pR) {
            double SR = uR + sqrt(((GAMMA + 1) / (2 * GAMMA)) * (p_star / pR) +
                                 (GAMMA - 1) / (2 * GAMMA)) * cR;
            if (x_over_t > SR)
                return {rhoR, uR, pR};
            else {
                double rho_star = rhoR * ((p_star / pR) + (GAMMA - 1) / (GAMMA + 1)) /
                                 (((GAMMA - 1) / (GAMMA + 1)) * (p_star / pR) + 1);
                return {rho_star, u_star_val, p_star};
            }
        } else {
            double SHR = uR + cR;
            double STR = u_star_val + sqrt(GAMMA * p_star /
                                          (rhoR * pow(p_star / pR, 1.0 / GAMMA)));
            if (x_over_t > SHR)
                return {rhoR, uR, pR};
            else if (x_over_t > STR) {
                double u = 2.0 / (GAMMA + 1) * (-cR + (GAMMA - 1) / 2 * uR + x_over_t);
                double c = 2.0 / (GAMMA + 1) * (cR - (GAMMA - 1) / 2 * (uR - x_over_t));
                double rho = rhoR * pow(c / cR, 2.0 / (GAMMA - 1));
                double p = pR * pow(rho / rhoR, GAMMA);
                return {rho, u, p};
            } else {
                double rho_star = rhoR * pow(p_star / pR, 1.0 / GAMMA);
                return {rho_star, u_star_val, p_star};
            }
        }
    }
}

tuple<double, double, double> exactRiemannSolver(double rhoL, double uL, double pL,
                                                  double rhoR, double uR, double pR,
                                                  double xi) {
    double p_star = solve_p_star(rhoL, uL, pL, rhoR, uR, pR);
    double u_star_val = u_star(p_star, rhoL, uL, pL, rhoR, uR, pR);
    return sample(p_star, u_star_val, rhoL, uL, pL, rhoR, uR, pR, xi);
}

struct Point {
    double x, y, z;
};

struct Triangle {
    int nodes[3];
    int edges[3];
    double area;
    double centroid[2];
};

struct Edge {
    int cellA, cellB;
    int edgeIdxA, edgeIdxB;
    double length;
    double nx, ny;
    int bcTag;
};

struct BoundaryLine {
    int tag;
    int n1, n2;
};

struct PairHash {
    size_t operator()(const pair<int,int>& p) const {
        return (static_cast<size_t>(p.first) << 32) | static_cast<size_t>(p.second);
    }
};

struct LocalEdge {
    int localA;
    int localB;
    double nx, ny;
    double length;
    int bcTag;
    bool is_boundary;
    bool ownA, ownB;
};

struct NeighborExchange {
    int rank;
    vector<int> send_ids;
    vector<int> recv_ids;
    vector<int> send_local;
    vector<int> recv_local;
    MPI_Request send_req, recv_req;
};

bool readGmshMesh(const string& filename,
                  vector<Point>& points,
                  vector<Triangle>& triangles,
                  vector<BoundaryLine>& boundaryLines) {
    ifstream in(filename);
    if (!in) {
        cerr << "Cannot open input file: " << filename << "\n";
        return false;
    }

    string line;
    while (getline(in, line)) if (line == "$MeshFormat") break;
    if (in.eof()) { cerr << "No $MeshFormat found.\n"; return false; }
    getline(in, line);

    while (getline(in, line)) if (line == "$Nodes") break;
    if (in.eof()) { cerr << "No $Nodes found.\n"; return false; }

    int numEntityBlocks, numNodes, minNodeTag, maxNodeTag;
    in >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;
    points.resize(numNodes);
    vector<int> nodeMap(maxNodeTag - minNodeTag + 1, -1);
    int idx = 0;
    for (int blk = 0; blk < numEntityBlocks; ++blk) {
        int dim, tag, parametric, numNodesInBlock;
        in >> dim >> tag >> parametric >> numNodesInBlock;
        vector<int> tags(numNodesInBlock);
        for (int i = 0; i < numNodesInBlock; ++i) in >> tags[i];
        for (int i = 0; i < numNodesInBlock; ++i) {
            double x, y, z;
            in >> x >> y >> z;
            int nodeTag = tags[i];
            int localIdx = nodeTag - minNodeTag;
            nodeMap[localIdx] = idx;
            points[idx].x = x; points[idx].y = y; points[idx].z = z;
            ++idx;
        }
    }

    while (getline(in, line)) if (line == "$Elements") break;
    if (in.eof()) { cerr << "No $Elements found.\n"; return false; }

    int numElemBlocks, numElements, minElemTag, maxElemTag;
    in >> numElemBlocks >> numElements >> minElemTag >> maxElemTag;

    vector<array<int,3>> triangleNodes;
    boundaryLines.clear();

    for (int blk = 0; blk < numElemBlocks; ++blk) {
        int dim, tag, elemType, numElemsInBlock;
        in >> dim >> tag >> elemType >> numElemsInBlock;
        if (elemType == 2) {
            for (int i = 0; i < numElemsInBlock; ++i) {
                int id, n1, n2, n3;
                in >> id >> n1 >> n2 >> n3;
                int i1 = nodeMap[n1 - minNodeTag];
                int i2 = nodeMap[n2 - minNodeTag];
                int i3 = nodeMap[n3 - minNodeTag];
                if (i1 >= 0 && i2 >= 0 && i3 >= 0)
                    triangleNodes.push_back({i1, i2, i3});
            }
        } else if (elemType == 1) {
            for (int i = 0; i < numElemsInBlock; ++i) {
                int id, n1, n2;
                in >> id >> n1 >> n2;
                int i1 = nodeMap[n1 - minNodeTag];
                int i2 = nodeMap[n2 - minNodeTag];
                boundaryLines.push_back({tag, i1, i2});
            }
        } else {
            int numNodesInElem = (elemType == 1 ? 2 : elemType == 2 ? 3 : elemType == 3 ? 4 : 0);
            for (int i = 0; i < numElemsInBlock; ++i) {
                int id; in >> id;
                for (int j = 0; j < numNodesInElem; ++j) { int dummy; in >> dummy; }
            }
        }
    }
    in.close();

    triangles.clear();
    for (auto& tri : triangleNodes) {
        int i0 = tri[0], i1 = tri[1], i2 = tri[2];
        double x1 = points[i1].x - points[i0].x, y1 = points[i1].y - points[i0].y;
        double x2 = points[i2].x - points[i0].x, y2 = points[i2].y - points[i0].y;
        double area = 0.5 * (x1*y2 - x2*y1);
        if (area < 0) {
            swap(i1, i2);
            area = -area;
        }
        Triangle t;
        t.nodes[0] = i0; t.nodes[1] = i1; t.nodes[2] = i2;
        t.area = area;
        t.centroid[0] = (points[i0].x + points[i1].x + points[i2].x) / 3.0;
        t.centroid[1] = (points[i0].y + points[i1].y + points[i2].y) / 3.0;
        triangles.push_back(t);
    }
    return true;
}

void buildEdges(const vector<Point>& points,
                vector<Triangle>& triangles,
                const vector<BoundaryLine>& boundaryLines,
                vector<Edge>& edges) {
    edges.clear();
    unordered_map<pair<int,int>, int, PairHash> edgeMap;

    for (size_t i = 0; i < triangles.size(); ++i) {
        Triangle& t = triangles[i];
        for (int e = 0; e < 3; ++e) {
            int n1 = t.nodes[e];
            int n2 = t.nodes[(e+1)%3];
            if (n1 > n2) swap(n1, n2);
            auto key = make_pair(n1, n2);
            auto it = edgeMap.find(key);
            if (it == edgeMap.end()) {
                Edge edge;
                edge.cellA = i;
                edge.cellB = -1;
                edge.edgeIdxA = e;
                edge.edgeIdxB = -1;
                edge.length = 0.0;
                edge.nx = edge.ny = 0.0;
                edge.bcTag = 0;
                int idx = edges.size();
                edges.push_back(edge);
                edgeMap[key] = idx;
                t.edges[e] = idx;
            } else {
                int idx = it->second;
                Edge& edge = edges[idx];
                edge.cellB = i;
                edge.edgeIdxB = e;
                t.edges[e] = idx;
            }
        }
    }

    unordered_map<pair<int,int>, int, PairHash> lineTagMap;
    for (const auto& line : boundaryLines) {
        int n1 = line.n1, n2 = line.n2;
        if (n1 > n2) swap(n1, n2);
        lineTagMap[make_pair(n1, n2)] = line.tag;
    }

    for (auto& edge : edges) {
        if (edge.cellB == -1) {
            const Triangle& t = triangles[edge.cellA];
            int n1 = t.nodes[edge.edgeIdxA];
            int n2 = t.nodes[(edge.edgeIdxA+1)%3];
            if (n1 > n2) swap(n1, n2);
            auto it = lineTagMap.find(make_pair(n1, n2));
            if (it != lineTagMap.end()) {
                edge.bcTag = it->second;
            } else {
                edge.bcTag = 1;
            }
        }
    }

    for (auto& edge : edges) {
        const Triangle& tA = triangles[edge.cellA];
        int i1 = tA.nodes[edge.edgeIdxA];
        int i2 = tA.nodes[(edge.edgeIdxA+1)%3];
        double dx = points[i2].x - points[i1].x;
        double dy = points[i2].y - points[i1].y;
        edge.length = sqrt(dx*dx + dy*dy);
        double nx = dy / edge.length;
        double ny = -dx / edge.length;
        double mx = (points[i1].x + points[i2].x) * 0.5;
        double my = (points[i1].y + points[i2].y) * 0.5;
        double cx = tA.centroid[0], cy = tA.centroid[1];
        double dot = (cx - mx) * nx + (cy - my) * ny;
        if (dot < 0) { nx = -nx; ny = -ny; }
        edge.nx = nx;
        edge.ny = ny;
    }
}

void primitiveToConservative(double rho, double u, double v, double p,
                             array<double,4>& U) {
    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = p/(GAMMA-1.0) + 0.5*rho*(u*u + v*v);
}

void conservativeToPrimitive(const array<double,4>& U,
                             double& rho, double& u, double& v, double& p) {
    rho = U[0];
    u = U[1] / rho;
    v = U[2] / rho;
    double E = U[3];
    p = (GAMMA-1.0) * (E - 0.5*rho*(u*u + v*v));
}

void computeFluxExact(const array<double,4>& UL, const array<double,4>& UR,
                      double nx, double ny, array<double,4>& flux) {
    double rhoL, uL, vL, pL;
    double rhoR, uR, vR, pR;
    conservativeToPrimitive(UL, rhoL, uL, vL, pL);
    conservativeToPrimitive(UR, rhoR, uR, vR, pR);

    double uNormL = uL*nx + vL*ny;
    double uNormR = uR*nx + vR*ny;
    double uTangL = -uL*ny + vL*nx;
    double uTangR = -uR*ny + vR*nx;

    auto [rhoFace, uNormFace, pFace] = exactRiemannSolver(rhoL, uNormL, pL,
                                                           rhoR, uNormR, pR,
                                                           0.0);
    double uTangFace = (uNormFace > 0.0) ? uTangL : uTangR;
    double uFace = uNormFace * nx - uTangFace * ny;
    double vFace = uNormFace * ny + uTangFace * nx;

    flux[0] = rhoFace * uNormFace;
    flux[1] = rhoFace * uNormFace * uFace + pFace * nx;
    flux[2] = rhoFace * uNormFace * vFace + pFace * ny;
    double EFace = pFace/(GAMMA-1.0) + 0.5*rhoFace*(uFace*uFace + vFace*vFace);
    flux[3] = (EFace + pFace) * uNormFace;
}

void computeWallFlux(const array<double,4>& UL, double nx, double ny,
                     array<double,4>& flux) {
    double rho, u, v, p;
    conservativeToPrimitive(UL, rho, u, v, p);
    double uNorm = u*nx + v*ny;
    double uTang = -u*ny + v*nx;
    double uNormWall = -uNorm;
    double uFace = uNormWall * nx - uTang * ny;
    double vFace = uNormWall * ny + uTang * nx;
    flux[0] = 0.0;
    flux[1] = p * nx;
    flux[2] = p * ny;
    flux[3] = 0.0;
}

void computeInletFlux(const array<double,4>& UL, double nx, double ny,
                      array<double,4>& flux) {
    double rhoL, uL, vL, pL;
    conservativeToPrimitive(UL, rhoL, uL, vL, pL);
    double rhoR = INLET_RHO;
    double uR   = INLET_U;
    double vR   = INLET_V;
    double pR   = INLET_P;

    double uNormL = uL*nx + vL*ny;
    double uNormR = uR*nx + vR*ny;
    double uTangL = -uL*ny + vL*nx;
    double uTangR = -uR*ny + vR*nx;

    auto [rhoFace, uNormFace, pFace] = exactRiemannSolver(rhoL, uNormL, pL,
                                                           rhoR, uNormR, pR,
                                                           0.0);
    double uTangFace = uTangR;
    double uFace = uNormFace * nx - uTangFace * ny;
    double vFace = uNormFace * ny + uTangFace * nx;

    flux[0] = rhoFace * uNormFace;
    flux[1] = rhoFace * uNormFace * uFace + pFace * nx;
    flux[2] = rhoFace * uNormFace * vFace + pFace * ny;
    double EFace = pFace/(GAMMA-1.0) + 0.5*rhoFace*(uFace*uFace + vFace*vFace);
    flux[3] = (EFace + pFace) * uNormFace;
}

void computeOutletFlux(const array<double,4>& UL, double nx, double ny,
                       array<double,4>& flux) {
    double rhoL, uL, vL, pL;
    conservativeToPrimitive(UL, rhoL, uL, vL, pL);
    double rhoR = rhoL;
    double uR   = uL;
    double vR   = vL;
    double pR   = OUTLET_P;

    double uNormL = uL*nx + vL*ny;
    double uNormR = uR*nx + vR*ny;
    double uTangL = -uL*ny + vL*nx;
    double uTangR = -uR*ny + vR*nx;

    auto [rhoFace, uNormFace, pFace] = exactRiemannSolver(rhoL, uNormL, pL,
                                                           rhoR, uNormR, pR,
                                                           0.0);
    double uTangFace = (uNormFace > 0.0) ? uTangL : uTangR;
    double uFace = uNormFace * nx - uTangFace * ny;
    double vFace = uNormFace * ny + uTangFace * nx;

    flux[0] = rhoFace * uNormFace;
    flux[1] = rhoFace * uNormFace * uFace + pFace * nx;
    flux[2] = rhoFace * uNormFace * vFace + pFace * ny;
    double EFace = pFace/(GAMMA-1.0) + 0.5*rhoFace*(uFace*uFace + vFace*vFace);
    flux[3] = (EFace + pFace) * uNormFace;
}

static void rcb_split(const vector<Triangle>& triangles, vector<int>& part,
                      vector<int>& ids, int rank_start, int proc_count) {
    if (proc_count == 1) {
        for (int id : ids) part[id] = rank_start;
        return;
    }
    double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
    for (int id : ids) {
        xmin = min(xmin, triangles[id].centroid[0]);
        xmax = max(xmax, triangles[id].centroid[0]);
        ymin = min(ymin, triangles[id].centroid[1]);
        ymax = max(ymax, triangles[id].centroid[1]);
    }
    bool split_x = (xmax - xmin) >= (ymax - ymin);
    sort(ids.begin(), ids.end(),
        [&](int a, int b) {
            return split_x ? triangles[a].centroid[0] < triangles[b].centroid[0]
                           : triangles[a].centroid[1] < triangles[b].centroid[1];
        });
    int left_procs = proc_count / 2;
    int split_idx = ids.size() * left_procs / proc_count;
    if (split_idx == 0) split_idx = 1;
    if (split_idx == (int)ids.size()) split_idx = ids.size() - 1;
    vector<int> left_ids(ids.begin(), ids.begin() + split_idx);
    vector<int> right_ids(ids.begin() + split_idx, ids.end());
    rcb_split(triangles, part, left_ids, rank_start, left_procs);
    rcb_split(triangles, part, right_ids, rank_start + left_procs, proc_count - left_procs);
}

void rcb_partition(const vector<Triangle>& triangles, vector<int>& part, int num_procs) {
    int n = triangles.size();
    part.resize(n);
    vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;
    rcb_split(triangles, part, indices, 0, num_procs);
}

void buildLocalStructures(const vector<Triangle>& triangles,
                          const vector<Edge>& edges,
                          const vector<int>& part,
                          vector<int>& my_cells,
                          vector<int>& ghost_cells,
                          vector<NeighborExchange>& neighbors,
                          vector<int>& local_to_global,
                          vector<LocalEdge>& local_edges) {
    int n_cells = triangles.size();
    my_cells.clear();
    for (int i = 0; i < n_cells; ++i)
        if (part[i] == mpi_rank)
            my_cells.push_back(i);
    int n_own = my_cells.size();

    unordered_map<int, int> global_to_local;
    for (int i = 0; i < n_own; ++i)
        global_to_local[my_cells[i]] = i;

    vector<int> ghost_set;
    for (const Edge& e : edges) {
        if (e.cellB == -1) continue;
        int cA = e.cellA, cB = e.cellB;
        if (part[cA] == mpi_rank && part[cB] != mpi_rank)
            ghost_set.push_back(cB);
        if (part[cB] == mpi_rank && part[cA] != mpi_rank)
            ghost_set.push_back(cA);
    }
    sort(ghost_set.begin(), ghost_set.end());
    ghost_set.erase(unique(ghost_set.begin(), ghost_set.end()), ghost_set.end());
    ghost_cells = ghost_set;
    int n_ghost = ghost_cells.size();

    local_to_global.resize(n_own + n_ghost);
    for (int i = 0; i < n_own; ++i) local_to_global[i] = my_cells[i];
    for (int i = 0; i < n_ghost; ++i) {
        local_to_global[n_own + i] = ghost_cells[i];
        global_to_local[ghost_cells[i]] = n_own + i;
    }

    struct TempNeighbor {
        int rank;
        vector<int> send_ids;
        vector<int> recv_ids;
    };
    vector<TempNeighbor> temp_neighbors;
    for (const Edge& e : edges) {
        if (e.cellB == -1) continue;
        int cA = e.cellA, cB = e.cellB;
        if (part[cA] == mpi_rank && part[cB] != mpi_rank) {
            int neigh = part[cB];
            bool found = false;
            for (auto& tn : temp_neighbors) {
                if (tn.rank == neigh) {
                    tn.send_ids.push_back(cA);
                    tn.recv_ids.push_back(cB);
                    found = true;
                    break;
                }
            }
            if (!found) {
                TempNeighbor tn;
                tn.rank = neigh;
                tn.send_ids.push_back(cA);
                tn.recv_ids.push_back(cB);
                temp_neighbors.push_back(tn);
            }
        }
        if (part[cB] == mpi_rank && part[cA] != mpi_rank) {
            int neigh = part[cA];
            bool found = false;
            for (auto& tn : temp_neighbors) {
                if (tn.rank == neigh) {
                    tn.send_ids.push_back(cB);
                    tn.recv_ids.push_back(cA);
                    found = true;
                    break;
                }
            }
            if (!found) {
                TempNeighbor tn;
                tn.rank = neigh;
                tn.send_ids.push_back(cB);
                tn.recv_ids.push_back(cA);
                temp_neighbors.push_back(tn);
            }
        }
    }

    neighbors.clear();
    for (auto& tn : temp_neighbors) {
        sort(tn.send_ids.begin(), tn.send_ids.end());
        tn.send_ids.erase(unique(tn.send_ids.begin(), tn.send_ids.end()), tn.send_ids.end());
        sort(tn.recv_ids.begin(), tn.recv_ids.end());
        tn.recv_ids.erase(unique(tn.recv_ids.begin(), tn.recv_ids.end()), tn.recv_ids.end());

        NeighborExchange nb;
        nb.rank = tn.rank;
        nb.send_ids = tn.send_ids;
        nb.recv_ids = tn.recv_ids;
        nb.send_local.resize(nb.send_ids.size());
        for (size_t i = 0; i < nb.send_ids.size(); ++i)
            nb.send_local[i] = global_to_local[nb.send_ids[i]];
        nb.recv_local.resize(nb.recv_ids.size());
        for (size_t i = 0; i < nb.recv_ids.size(); ++i)
            nb.recv_local[i] = global_to_local[nb.recv_ids[i]];
        neighbors.push_back(nb);
    }

    local_edges.clear();
    for (const Edge& e : edges) {
        int cA = e.cellA, cB = e.cellB;
        bool ownA = (part[cA] == mpi_rank);
        bool ownB = (cB != -1 && part[cB] == mpi_rank);
        if (!ownA && !ownB) continue;
        LocalEdge le;
        le.nx = e.nx; le.ny = e.ny; le.length = e.length;
        le.bcTag = e.bcTag;
        le.is_boundary = (cB == -1);
        le.localA = global_to_local[cA];
        le.ownA = ownA;
        if (cB != -1) {
            le.localB = global_to_local[cB];
            le.ownB = ownB;
        } else {
            le.localB = -1;
            le.ownB = false;
        }
        local_edges.push_back(le);
    }
}

void initializeSolutionLocal(const vector<Triangle>& triangles,
                             const vector<int>& my_cells,
                             vector<array<double,4>>& U_local) {
    int n_own = my_cells.size();
    for (int i = 0; i < n_own; ++i) {
        int gid = my_cells[i];
        double xc = triangles[gid].centroid[0];
        double rho, u, v, p;
        if (xc < -0.5) {
            rho = 1.0;
            u   = 1.0;
            v   = 0.0;
            p   = 1.0;
        } else {
            rho = 0.125;
            u   = 0.0;
            v   = 0.0;
            p   = 0.1;
        }
        primitiveToConservative(rho, u, v, p, U_local[i]);
    }
}

void exchangeHalo(const vector<NeighborExchange>& neighbors,
                  vector<array<double,4>>& U_local) {
    int n_neigh = neighbors.size();
    vector<vector<double>> send_bufs(n_neigh);
    vector<vector<double>> recv_bufs(n_neigh);
    vector<MPI_Request> reqs(2 * n_neigh);

    for (int i = 0; i < n_neigh; ++i) {
        const auto& nb = neighbors[i];
        int n_recv = nb.recv_local.size();
        recv_bufs[i].resize(n_recv * 4);
        MPI_Irecv(recv_bufs[i].data(), n_recv * 4, MPI_DOUBLE, nb.rank, 0,
                  MPI_COMM_WORLD, &reqs[2*i]);
    }
    for (int i = 0; i < n_neigh; ++i) {
        const auto& nb = neighbors[i];
        int n_send = nb.send_local.size();
        send_bufs[i].resize(n_send * 4);
        for (int j = 0; j < n_send; ++j) {
            int local_idx = nb.send_local[j];
            for (int k = 0; k < 4; ++k)
                send_bufs[i][j*4 + k] = U_local[local_idx][k];
        }
        MPI_Isend(send_bufs[i].data(), n_send * 4, MPI_DOUBLE, nb.rank, 0,
                  MPI_COMM_WORLD, &reqs[2*i + 1]);
    }
    MPI_Waitall(2 * n_neigh, reqs.data(), MPI_STATUSES_IGNORE);
    for (int i = 0; i < n_neigh; ++i) {
        const auto& nb = neighbors[i];
        for (size_t j = 0; j < nb.recv_local.size(); ++j) {
            int local_idx = nb.recv_local[j];
            for (int k = 0; k < 4; ++k)
                U_local[local_idx][k] = recv_bufs[i][j*4 + k];
        }
    }
}

void writeVTK(const string& filename,
              const vector<Point>& points,
              const vector<Triangle>& triangles,
              const vector<array<double,4>>& U,
              const vector<int>& part,
              double time) {
    if (mpi_rank != 0) return;
    ofstream out(filename);
    if (!out) { cerr << "Cannot create output file: " << filename << "\n"; return; }

    out << "# vtk DataFile Version 3.0\n";
    out << "Mesh from Gmsh, time = " << time << "\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << points.size() << " double\n";
    for (const auto& p : points) out << p.x << " " << p.y << " " << p.z << "\n";

    int numCells = triangles.size();
    int totalIndices = 0;
    for (const auto& t : triangles) totalIndices += 3 + 1;
    out << "CELLS " << numCells << " " << totalIndices << "\n";
    for (const auto& t : triangles) out << "3 " << t.nodes[0] << " " << t.nodes[1] << " " << t.nodes[2] << "\n";

    out << "CELL_TYPES " << numCells << "\n";
    for (size_t i = 0; i < triangles.size(); ++i) out << "5\n";

    out << "CELL_DATA " << numCells << "\n";
    out << "SCALARS process_rank int 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < triangles.size(); ++i) out << part[i] << "\n";
    out << "SCALARS density double 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < triangles.size(); ++i) out << U[i][0] << "\n";
    out << "SCALARS pressure double 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < triangles.size(); ++i) {
        double rho, u, v, p;
        conservativeToPrimitive(U[i], rho, u, v, p);
        out << p << "\n";
    }
    out << "VECTORS velocity double\n";
    for (size_t i = 0; i < triangles.size(); ++i) {
        double rho, u, v, p;
        conservativeToPrimitive(U[i], rho, u, v, p);
        out << u << " " << v << " 0.0\n";
    }
    out.close();
    if (mpi_rank == 0)
        cout << "VTK file written: " << filename << " at time " << time << "\n";
}

void gatherSolution(const vector<array<double,4>>& U_local,
                    const vector<int>& my_cells,
                    const vector<Triangle>& triangles,
                    vector<array<double,4>>& global_U) {
    int n_cells = triangles.size();
    global_U.resize(n_cells);
    int n_own = my_cells.size();

    vector<int> recv_counts(mpi_size), displs(mpi_size);
    MPI_Gather(&n_own, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> all_indices;
    if (mpi_rank == 0) {
        all_indices.resize(n_cells);
        int offset = 0;
        for (int i = 0; i < mpi_size; ++i) {
            displs[i] = offset;
            offset += recv_counts[i];
        }
    }
    vector<int> send_idx = my_cells;
    MPI_Gatherv(send_idx.data(), n_own, MPI_INT,
                all_indices.data(), recv_counts.data(), displs.data(), MPI_INT,
                0, MPI_COMM_WORLD);

    vector<double> send_buf(n_own * 4);
    for (int i = 0; i < n_own; ++i) {
        for (int k = 0; k < 4; ++k)
            send_buf[i*4 + k] = U_local[i][k];
    }

    vector<double> recv_buf;
    if (mpi_rank == 0) {
        recv_buf.resize(n_cells * 4);
    }

    vector<int> recv_counts_data(mpi_size), displs_data(mpi_size);
    if (mpi_rank == 0) {
        for (int i = 0; i < mpi_size; ++i) {
            recv_counts_data[i] = recv_counts[i] * 4;
            displs_data[i] = displs[i] * 4;
        }
    }

    MPI_Gatherv(send_buf.data(), n_own * 4, MPI_DOUBLE,
                recv_buf.data(), recv_counts_data.data(), displs_data.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        for (int p = 0; p < mpi_size; ++p) {
            int start = displs[p];
            int count = recv_counts[p];
            for (int j = 0; j < count; ++j) {
                int gid = all_indices[start + j];
                for (int k = 0; k < 4; ++k)
                    global_U[gid][k] = recv_buf[(start + j) * 4 + k];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    double start_time = MPI_Wtime();

    if (argc < 2 || argc > 3) {
        if (mpi_rank == 0)
            cout << "Usage: " << argv[0] << " file.msh [output_prefix]\n";
        MPI_Finalize();
        return 1;
    }

    string meshFile = argv[1];
    string outputPrefix;
    if (argc >= 3) {
        outputPrefix = argv[2];
    } else {
        outputPrefix = meshFile;
        size_t dot = outputPrefix.rfind('.');
        if (dot != string::npos && outputPrefix.substr(dot) == ".msh")
            outputPrefix = outputPrefix.substr(0, dot);
    }

    vector<Point> points;
    vector<Triangle> triangles;
    vector<BoundaryLine> boundaryLines;
    if (!readGmshMesh(meshFile, points, triangles, boundaryLines)) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    vector<Edge> edges;
    buildEdges(points, triangles, boundaryLines, edges);

    vector<int> part;
    rcb_partition(triangles, part, mpi_size);

    vector<int> my_cells, ghost_cells, local_to_global;
    vector<NeighborExchange> neighbors;
    vector<LocalEdge> local_edges;
    buildLocalStructures(triangles, edges, part, my_cells, ghost_cells,
                         neighbors, local_to_global, local_edges);

    int n_own = my_cells.size();
    int n_ghost = ghost_cells.size();
    int n_local = n_own + n_ghost;

    vector<array<double,4>> U_local(n_local);
    initializeSolutionLocal(triangles, my_cells, U_local);

    vector<double> local_area(n_local);
    for (int i = 0; i < n_local; ++i) {
        int gid = local_to_global[i];
        local_area[i] = triangles[gid].area;
    }

    double t = 0.0;
    int iter = 0;
    double nextOutput = OUTPUT_INTERVAL;
    int fileCounter = 0;

    if (mpi_rank == 0) {
        cout << "Starting parallel Godunov solver on " << mpi_size << " processes.\n";
        cout << "Number of cells: " << triangles.size() << "\n";
        cout << "My cells: " << n_own << ", ghost: " << n_ghost << "\n";
    }

    while (t < FINAL_TIME) {
        exchangeHalo(neighbors, U_local);

        double dt_local = 1e30;
        double max_speed_local = 0.0;
        for (int i = 0; i < n_own; ++i) {
            double area = local_area[i];
            double rho, u, v, p;
            conservativeToPrimitive(U_local[i], rho, u, v, p);
            double c = sqrt(GAMMA * p / rho);
            double speed = sqrt(u*u + v*v) + c;
            if (speed > max_speed_local) max_speed_local = speed;
            if (area < dt_local) dt_local = area;
        }
        if (max_speed_local > 1e-12)
            dt_local = CFL * dt_local / max_speed_local;
        else
            dt_local = 1e-6;

        double dt_global;
        MPI_Allreduce(&dt_local, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (t + dt_global > FINAL_TIME) dt_global = FINAL_TIME - t;

        vector<array<double,4>> flux_local(n_local, {0.0,0.0,0.0,0.0});
        for (const LocalEdge& le : local_edges) {
            const array<double,4>& UL = U_local[le.localA];
            array<double,4> flux;
            if (le.is_boundary) {
                if (le.bcTag == 3) computeWallFlux(UL, le.nx, le.ny, flux);
                else if (le.bcTag == 4) computeInletFlux(UL, le.nx, le.ny, flux);
                else if (le.bcTag == 11) computeOutletFlux(UL, le.nx, le.ny, flux);
                else computeWallFlux(UL, le.nx, le.ny, flux);
            } else {
                const array<double,4>& UR = U_local[le.localB];
                computeFluxExact(UL, UR, le.nx, le.ny, flux);
            }
            for (int k = 0; k < 4; ++k) flux[k] *= le.length;
            for (int k = 0; k < 4; ++k) flux_local[le.localA][k] -= flux[k];
            if (!le.is_boundary && le.ownB) {
                for (int k = 0; k < 4; ++k) flux_local[le.localB][k] += flux[k];
            }
        }
        for (int i = 0; i < n_own; ++i) {
            double area = local_area[i];
            for (int k = 0; k < 4; ++k) {
                U_local[i][k] += (dt_global / area) * flux_local[i][k];
            }
            double rho = U_local[i][0];
            double E = U_local[i][3];
            double u = U_local[i][1] / rho;
            double v = U_local[i][2] / rho;
            double p = (GAMMA-1.0) * (E - 0.5*rho*(u*u + v*v));
            if (rho < 1e-12) rho = 1e-12;
            if (p < 1e-12) p = 1e-12;
            primitiveToConservative(rho, u, v, p, U_local[i]);
        }

        t += dt_global;
        iter++;

        if (t >= nextOutput) {
            vector<array<double,4>> global_U;
            gatherSolution(U_local, my_cells, triangles, global_U);
            if (mpi_rank == 0) {
                char outName[256];
                sprintf(outName, "%s_%06d.vtk", outputPrefix.c_str(), fileCounter);
                writeVTK(outName, points, triangles, global_U, part, t);
                fileCounter++;
            }
            nextOutput += OUTPUT_INTERVAL;
        }

        if (iter % 100 == 0 && mpi_rank == 0) {
            cout << "Iteration " << iter << ", time = " << t << ", dt = " << dt_global << "\n";
        }
    }

    vector<array<double,4>> global_U;
    gatherSolution(U_local, my_cells, triangles, global_U);
    if (mpi_rank == 0) {
        string finalVTK = outputPrefix + "_final.vtk";
        writeVTK(finalVTK.c_str(), points, triangles, global_U, part, t);
        cout << "Simulation finished. Final time = " << t << "\n";
    }

    double end_time = MPI_Wtime();
    if (mpi_rank == 0) {
        cout << "Total wall time: " << (end_time - start_time) << " seconds" << endl;
    }
    MPI_Finalize();
    return 0;
}