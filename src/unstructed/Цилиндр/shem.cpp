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

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path) mkdir(path, 0777)
#endif

using namespace std;

const double GAMMA = 1.4;
const double CFL = 0.8;
const double FINAL_TIME = 4.0;
const double OUTPUT_INTERVAL = 0.05;

const double INLET_RHO = 1.0;
const double INLET_U   = 1.0;
const double INLET_V   = 0.0;
const double INLET_P   = 1.0;

const double OUTLET_P = 1.0;

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
                                                  double xi = 0.0) {
    double p_star = solve_p_star(rhoL, uL, pL, rhoR, uR, pR);
    double u_star_val = u_star(p_star, rhoL, uL, pL, rhoR, uR, pR);
    return sample(p_star, u_star_val, rhoL, uL, pL, rhoR, uR, pR, xi);
}
// CHECK: MESH_STRUCTS
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
// CHECK: BUILD_FACES
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

double computeTimeStep(const vector<Triangle>& triangles,
                       const vector<Edge>& edges,
                       const vector<array<double,4>>& U,
                       double cfl) {
    double minArea = numeric_limits<double>::max();
    double maxSpeed = 0.0;
    for (size_t i = 0; i < triangles.size(); ++i) {
        double rho, u, v, p;
        conservativeToPrimitive(U[i], rho, u, v, p);
        double c = sqrt(GAMMA * p / rho);
        double speed = sqrt(u*u + v*v) + c;
        if (triangles[i].area < minArea) minArea = triangles[i].area;
        if (speed > maxSpeed) maxSpeed = speed;
    }
    double maxEdgeLen = 0.0;
    for (const auto& e : edges) if (e.length > maxEdgeLen) maxEdgeLen = e.length;
    return cfl * minArea / (maxSpeed * maxEdgeLen);
}
// CHECK: UNSTRUCT_SCHEMES
void updateSolution(vector<array<double,4>>& U,
                    const vector<Triangle>& triangles,
                    const vector<Edge>& edges,
                    double dt) {
    vector<array<double,4>> Unew = U;

    for (const auto& edge : edges) {
        const array<double,4>& UL = U[edge.cellA];
        array<double,4> flux;
        bool isBoundary = (edge.cellB == -1);
// CHECK: GHOST_STATE
        if (isBoundary) {
            if (edge.bcTag == 3) {
                computeWallFlux(UL, edge.nx, edge.ny, flux);
            } else if (edge.bcTag == 4) {
                computeInletFlux(UL, edge.nx, edge.ny, flux);
            } else if (edge.bcTag == 11) {
                computeOutletFlux(UL, edge.nx, edge.ny, flux);
            } else {
                computeWallFlux(UL, edge.nx, edge.ny, flux);
            }
        } else {
            const array<double,4>& UR = U[edge.cellB];
            computeFluxExact(UL, UR, edge.nx, edge.ny, flux);
        }

        for (int i = 0; i < 4; ++i) flux[i] *= edge.length;

        for (int i = 0; i < 4; ++i)
            Unew[edge.cellA][i] -= (dt / triangles[edge.cellA].area) * flux[i];

        if (!isBoundary) {
            for (int i = 0; i < 4; ++i)
                Unew[edge.cellB][i] += (dt / triangles[edge.cellB].area) * flux[i];
        }
    }

    for (size_t i = 0; i < U.size(); ++i) {
        double rho, u, v, p;
        conservativeToPrimitive(Unew[i], rho, u, v, p);
        if (rho < 1e-12) rho = 1e-12;
        if (p < 1e-12) p = 1e-12;
        primitiveToConservative(rho, u, v, p, U[i]);
    }
}

void writeVTK(const string& filename,
              const vector<Point>& points,
              const vector<Triangle>& triangles,
              const vector<array<double,4>>& U,
              double time) {
    //              results,           
    MKDIR("results");

    //                             (        )
    size_t lastSlash = filename.find_last_of("/\\");
    string basename = (lastSlash == string::npos) ? filename : filename.substr(lastSlash + 1);

    string outName = "results/" + basename;
    //                        .vtk
    size_t dot = outName.rfind('.');
    if (dot != string::npos) {
        outName = outName.substr(0, dot) + ".vtk";
    } else {
        outName += ".vtk";
    }

    ofstream out(outName);
    if (!out) { cerr << "Cannot create output file.\n"; return; }

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
    cout << "VTK file written: " << outName << " at time " << time << "\n";
}

void initializeSolution(const vector<Triangle>& triangles,
                        vector<array<double,4>>& U) {
    U.resize(triangles.size());
    for (size_t i = 0; i < triangles.size(); ++i) {
        double rho = 1.0, u = 0.0, v = 0.0, p = 1.0;
        primitiveToConservative(rho, u, v, p, U[i]);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " file.msh\n";
        return 1;
    }

    vector<Point> points;
    vector<Triangle> triangles;
    vector<BoundaryLine> boundaryLines;

    if (!readGmshMesh(argv[1], points, triangles, boundaryLines)) return 1;

    cout << "Points: " << points.size() << "\n";
    cout << "Triangles: " << triangles.size() << "\n";
    cout << "Boundary lines: " << boundaryLines.size() << "\n";

    vector<Edge> edges;
    buildEdges(points, triangles, boundaryLines, edges);
    cout << "Edges: " << edges.size() << "\n";

    vector<array<double,4>> U;
    initializeSolution(triangles, U);
    string baseName = argv[1];
    size_t dot = baseName.rfind('.');
    if (dot != string::npos && baseName.substr(dot) == ".msh")
        baseName = baseName.substr(0, dot);
    double t = 0.0;
    int iter = 0;
    double nextOutput = OUTPUT_INTERVAL;
    int fileCounter = 0;
    cout << "Starting Godunov time integration with exact Riemann solver...\n";
    while (t < FINAL_TIME) {
        double dt = computeTimeStep(triangles, edges, U, CFL);
        if (t + dt > FINAL_TIME) dt = FINAL_TIME - t;
        updateSolution(U, triangles, edges, dt);
        t += dt;
        iter++;
        if (t >= nextOutput) {
            char outName[256];
            sprintf(outName, "%s_%06d.vtk", baseName.c_str(), fileCounter);
            writeVTK(outName, points, triangles, U, t);
            fileCounter++;
            nextOutput += OUTPUT_INTERVAL;
        }
        if (iter % 100 == 0) {
            cout << "Iteration " << iter << ", time = " << t << ", dt = " << dt << "\n";
        }
    }

    writeVTK(argv[1], points, triangles, U, t);
    cout << "Simulation finished. Final time = " << t << "\n";
    return 0;
}