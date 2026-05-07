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
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _WIN32
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
#endif
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;

const double GAMMA = 1.4;
const double CFL = 0.3;
const double FINAL_TIME = 1.0;        // увеличено для выхода на стационар
const double OUTPUT_INTERVAL = 0.1;    // реже выводим результаты
const double ALPHA_DEG = 0.0;           // угол атаки (градусы)

// Параметры невозмущённого потока
const double MACH_INF = 0.8;
const double CHORD = 1.0;
const double INLET_P = 1.0;
const double INLET_RHO = 1.0;
const double OUTLET_P = 1.0;

// Скорость звука и компоненты скорости в связанной с профилем системе координат
double SOUND_SPEED, INLET_U, INLET_V;

// ------------------------------------------------------------
// Структуры сетки и функции чтения (без изменений)
// ------------------------------------------------------------
struct Point { double x, y, z; };
struct Triangle {
    int nodes[3], edges[3];
    double area, centroid[2];
};
struct Edge {
    int cellA, cellB, edgeIdxA, edgeIdxB;
    double length, nx, ny;
    int bcTag;
};
struct BoundaryLine { int tag, n1, n2; };
struct PairHash {
    size_t operator()(const pair<int,int>& p) const {
        return (static_cast<size_t>(p.first) << 32) | p.second;
    }
};

bool readGmshMesh(const string& filename, vector<Point>& points,
                  vector<Triangle>& triangles, vector<BoundaryLine>& boundaryLines) {
    ifstream in(filename);
    if (!in) { cerr << "Cannot open file: " << filename << "\n"; return false; }
    string line;
    while (getline(in, line)) if (line == "$MeshFormat") break;
    if (in.eof()) { cerr << "No $MeshFormat\n"; return false; }
    getline(in, line);
    while (getline(in, line)) if (line == "$Nodes") break;
    if (in.eof()) { cerr << "No $Nodes\n"; return false; }
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
            double x, y, z; in >> x >> y >> z;
            int nodeTag = tags[i];
            int localIdx = nodeTag - minNodeTag;
            nodeMap[localIdx] = idx;
            points[idx].x = x; points[idx].y = y; points[idx].z = z;
            ++idx;
        }
    }
    while (getline(in, line)) if (line == "$Elements") break;
    if (in.eof()) { cerr << "No $Elements\n"; return false; }
    int numElemBlocks, numElements, minElemTag, maxElemTag;
    in >> numElemBlocks >> numElements >> minElemTag >> maxElemTag;
    vector<array<int,3>> triangleNodes;
    boundaryLines.clear();
    for (int blk = 0; blk < numElemBlocks; ++blk) {
        int dim, tag, elemType, numElemsInBlock;
        in >> dim >> tag >> elemType >> numElemsInBlock;
        if (elemType == 2) {  // треугольники
            for (int i = 0; i < numElemsInBlock; ++i) {
                int id, n1, n2, n3; in >> id >> n1 >> n2 >> n3;
                int i1 = nodeMap[n1 - minNodeTag];
                int i2 = nodeMap[n2 - minNodeTag];
                int i3 = nodeMap[n3 - minNodeTag];
                if (i1 >= 0 && i2 >= 0 && i3 >= 0)
                    triangleNodes.push_back({i1, i2, i3});
            }
        } else if (elemType == 1) { // линии границ
            for (int i = 0; i < numElemsInBlock; ++i) {
                int id, n1, n2; in >> id >> n1 >> n2;
                int i1 = nodeMap[n1 - minNodeTag];
                int i2 = nodeMap[n2 - minNodeTag];
                boundaryLines.push_back({tag, i1, i2});
            }
        } else {
            int nNode = (elemType == 1 ? 2 : elemType == 2 ? 3 : elemType == 3 ? 4 : 0);
            for (int i = 0; i < numElemsInBlock; ++i) {
                int id; in >> id;
                for (int j = 0; j < nNode; ++j) { int dummy; in >> dummy; }
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
        if (area < 0) { swap(i1, i2); area = -area; }
        Triangle t;
        t.nodes[0] = i0; t.nodes[1] = i1; t.nodes[2] = i2;
        t.area = area;
        t.centroid[0] = (points[i0].x + points[i1].x + points[i2].x) / 3.0;
        t.centroid[1] = (points[i0].y + points[i1].y + points[i2].y) / 3.0;
        triangles.push_back(t);
    }
    return true;
}

void buildEdges(const vector<Point>& points, vector<Triangle>& triangles,
                const vector<BoundaryLine>& boundaryLines, vector<Edge>& edges) {
    edges.clear();
    unordered_map<pair<int,int>, int, PairHash> edgeMap;
    for (size_t i = 0; i < triangles.size(); ++i) {
        Triangle& t = triangles[i];
        for (int e = 0; e < 3; ++e) {
            int n1 = t.nodes[e], n2 = t.nodes[(e+1)%3];
            if (n1 > n2) swap(n1, n2);
            auto key = make_pair(n1, n2);
            auto it = edgeMap.find(key);
            if (it == edgeMap.end()) {
                Edge edge;
                edge.cellA = i; edge.cellB = -1;
                edge.edgeIdxA = e; edge.edgeIdxB = -1;
                edge.length = 0.0; edge.nx = edge.ny = 0.0; edge.bcTag = 0;
                int idx = edges.size();
                edges.push_back(edge);
                edgeMap[key] = idx;
                t.edges[e] = idx;
            } else {
                int idx = it->second;
                Edge& edge = edges[idx];
                edge.cellB = i; edge.edgeIdxB = e;
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
    for (Edge& edge : edges) {
        if (edge.cellB == -1) {
            const Triangle& t = triangles[edge.cellA];
            int n1 = t.nodes[edge.edgeIdxA], n2 = t.nodes[(edge.edgeIdxA+1)%3];
            if (n1 > n2) swap(n1, n2);
            auto it = lineTagMap.find(make_pair(n1, n2));
            edge.bcTag = (it != lineTagMap.end()) ? it->second : 1;
        }
    }
    // Вычисление нормалей с явным заданием для границ
    for (Edge& edge : edges) {
        const Triangle& tA = triangles[edge.cellA];
        int i1 = tA.nodes[edge.edgeIdxA], i2 = tA.nodes[(edge.edgeIdxA+1)%3];
        double dx = points[i2].x - points[i1].x, dy = points[i2].y - points[i1].y;
        edge.length = sqrt(dx*dx + dy*dy);
        double nx = dy / edge.length, ny = -dx / edge.length;
        double mx = (points[i1].x + points[i2].x) * 0.5;
        double my = (points[i1].y + points[i2].y) * 0.5;
        double cx = tA.centroid[0], cy = tA.centroid[1];
        if (edge.cellB != -1) { // внутреннее ребро
            if ((cx - mx)*nx + (cy - my)*ny < 0) { nx = -nx; ny = -ny; }
            edge.nx = nx; edge.ny = ny;
        } else { // граничное ребро
            // Явно задаём нормали по геометрии
            if (fabs(mx + 3.0) < 1e-5) {         // левая граница
                edge.nx = -1.0; edge.ny = 0.0;
            } else if (fabs(mx - 8.0) < 1e-5) {  // правая граница
                edge.nx = 1.0; edge.ny = 0.0;
            } else if (fabs(my - 2.0) < 1e-5) {  // верхняя граница
                edge.nx = 0.0; edge.ny = 1.0;
            } else if (fabs(my + 2.0) < 1e-5) {  // нижняя граница
                edge.nx = 0.0; edge.ny = -1.0;
            } else {
                // Профиль – оставляем нормаль, ориентированную наружу из ячейки (внутрь твёрдого тела)
                if ((cx - mx)*nx + (cy - my)*ny < 0) { nx = -nx; ny = -ny; }
                edge.nx = nx; edge.ny = ny;
            }
        }
    }
}

// ------------------------------------------------------------
// Решатель Римана (точный)
// ------------------------------------------------------------
inline double WaveFunk(double p, double pk, double rho_k) {
    if (p > pk) {
        double Ak = 2.0 / ((GAMMA + 1) * rho_k);
        double Bk = (GAMMA - 1) / (GAMMA + 1) * pk;
        return (p - pk) * sqrt(Ak / (p + Bk));
    } else {
        return 2.0 * sqrt(GAMMA * pk / rho_k) / (GAMMA - 1) *
               (pow(p / pk, (GAMMA - 1) / (2 * GAMMA)) - 1.0);
    }
}
inline double ProisvWaveFunk(double p, double pk, double rho_k) {
    if (p > pk) {
        double Ak = 2.0 / ((GAMMA + 1) * rho_k);
        double Bk = (GAMMA - 1) / (GAMMA + 1) * pk;
        return sqrt(Ak / (p + Bk)) * (1.0 - (p - pk) / (2.0 * (p + Bk)));
    } else {
        return (1.0 / (rho_k * sqrt(GAMMA * pk / rho_k))) *
               pow(p / pk, -(GAMMA + 1) / (2 * GAMMA));
    }
}
double solve_p_star(double rhoL, double uL, double pL, double rhoR, double uR, double pR) {
    double tol = 1e-8;
    int max_iter = 100;
    double p0 = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rhoL + rhoR) *
                (sqrt(GAMMA * pL / rhoL) + sqrt(GAMMA * pR / rhoR));
    p0 = max(tol, p0);
    double p_star = p0;
    for (int i = 0; i < max_iter; ++i) {
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
tuple<double, double, double> sample(double p_star, double u_star_val,
    double rhoL, double uL, double pL, double rhoR, double uR, double pR, double x_over_t) {
    double cL = sqrt(GAMMA * pL / rhoL), cR = sqrt(GAMMA * pR / rhoR);
    if (x_over_t < u_star_val) {
        if (p_star > pL) {
            double SL = uL - sqrt(((GAMMA+1)/(2*GAMMA))*(p_star/pL) + (GAMMA-1)/(2*GAMMA)) * cL;
            if (x_over_t < SL) return {rhoL, uL, pL};
            else {
                double rho_star = rhoL * ((p_star/pL) + (GAMMA-1)/(GAMMA+1)) /
                                 (((GAMMA-1)/(GAMMA+1))*(p_star/pL) + 1);
                return {rho_star, u_star_val, p_star};
            }
        } else {
            double SHL = uL - cL, STL = u_star_val - sqrt(GAMMA*p_star/(rhoL*pow(p_star/pL,1/GAMMA)));
            if (x_over_t < SHL) return {rhoL, uL, pL};
            else if (x_over_t < STL) {
                double u = 2.0/(GAMMA+1)*(cL + (GAMMA-1)/2*uL + x_over_t);
                double c = 2.0/(GAMMA+1)*(cL + (GAMMA-1)/2*(uL - x_over_t));
                double rho = rhoL * pow(c/cL, 2.0/(GAMMA-1));
                double p = pL * pow(rho/rhoL, GAMMA);
                return {rho, u, p};
            } else {
                double rho_star = rhoL * pow(p_star/pL, 1.0/GAMMA);
                return {rho_star, u_star_val, p_star};
            }
        }
    } else {
        if (p_star > pR) {
            double SR = uR + sqrt(((GAMMA+1)/(2*GAMMA))*(p_star/pR) + (GAMMA-1)/(2*GAMMA)) * cR;
            if (x_over_t > SR) return {rhoR, uR, pR};
            else {
                double rho_star = rhoR * ((p_star/pR) + (GAMMA-1)/(GAMMA+1)) /
                                 (((GAMMA-1)/(GAMMA+1))*(p_star/pR) + 1);
                return {rho_star, u_star_val, p_star};
            }
        } else {
            double SHR = uR + cR, STR = u_star_val + sqrt(GAMMA*p_star/(rhoR*pow(p_star/pR,1/GAMMA)));
            if (x_over_t > SHR) return {rhoR, uR, pR};
            else if (x_over_t > STR) {
                double u = 2.0/(GAMMA+1)*(-cR + (GAMMA-1)/2*uR + x_over_t);
                double c = 2.0/(GAMMA+1)*(cR - (GAMMA-1)/2*(uR - x_over_t));
                double rho = rhoR * pow(c/cR, 2.0/(GAMMA-1));
                double p = pR * pow(rho/rhoR, GAMMA);
                return {rho, u, p};
            } else {
                double rho_star = rhoR * pow(p_star/pR, 1.0/GAMMA);
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

// ------------------------------------------------------------
// Вспомогательные функции для перехода между переменными
// ------------------------------------------------------------
void primitiveToConservative(double rho, double u, double v, double p, array<double,4>& U) {
    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = p/(GAMMA-1.0) + 0.5*rho*(u*u + v*v);
}
void conservativeToPrimitive(const array<double,4>& U, double& rho, double& u, double& v, double& p) {
    rho = U[0];
    u = U[1] / rho;
    v = U[2] / rho;
    double E = U[3];
    p = (GAMMA-1.0) * (E - 0.5*rho*(u*u + v*v));
}

// ------------------------------------------------------------
// Вычисление потока через грань (точный решатель Римана)
// ------------------------------------------------------------
void computeFluxExact(const array<double,4>& UL, const array<double,4>& UR,
                      double nx, double ny, array<double,4>& flux) {
    double rhoL, uL, vL, pL, rhoR, uR, vR, pR;
    conservativeToPrimitive(UL, rhoL, uL, vL, pL);
    conservativeToPrimitive(UR, rhoR, uR, vR, pR);
    double uNormL = uL*nx + vL*ny, uNormR = uR*nx + vR*ny;
    double uTangL = -uL*ny + vL*nx, uTangR = -uR*ny + vR*nx;
    auto [rhoFace, uNormFace, pFace] = exactRiemannSolver(rhoL, uNormL, pL, rhoR, uNormR, pR, 0.0);
    double uTangFace = (uNormFace > 0.0) ? uTangL : uTangR;
    double uFace = uNormFace * nx - uTangFace * ny;
    double vFace = uNormFace * ny + uTangFace * nx;
    flux[0] = rhoFace * uNormFace;
    flux[1] = rhoFace * uNormFace * uFace + pFace * nx;
    flux[2] = rhoFace * uNormFace * vFace + pFace * ny;
    double EFace = pFace/(GAMMA-1.0) + 0.5*rhoFace*(uFace*uFace + vFace*vFace);
    flux[3] = (EFace + pFace) * uNormFace;
}
void computeInletFluxSimple(const array<double,4>& UL, double nx, double ny, array<double,4>& flux) {
    double alpha_rad = ALPHA_DEG * M_PI / 180.0;
    double uR = MACH_INF * SOUND_SPEED * cos(alpha_rad);
    double vR = MACH_INF * SOUND_SPEED * sin(alpha_rad);
    double rhoR = INLET_RHO;
    double pR = INLET_P;
    double u_n = uR * nx + vR * ny;
    flux[0] = rhoR * u_n;
    flux[1] = rhoR * u_n * uR + pR * nx;
    flux[2] = rhoR * u_n * vR + pR * ny;
    double ER = pR/(GAMMA-1.0) + 0.5 * rhoR * (uR*uR + vR*vR);
    flux[3] = (ER + pR) * u_n;
}

void computeOutletFluxSimple(const array<double,4>& UL, double nx, double ny, array<double,4>& flux) {
    double rhoL, uL, vL, pL;
    conservativeToPrimitive(UL, rhoL, uL, vL, pL);
    double u_n = uL * nx + vL * ny;
    flux[0] = rhoL * u_n;
    flux[1] = rhoL * u_n * uL + pL * nx;
    flux[2] = rhoL * u_n * vL + pL * ny;
    double EL = pL/(GAMMA-1.0) + 0.5 * rhoL * (uL*uL + vL*vL);
    flux[3] = (EL + pL) * u_n;
}
// ------------------------------------------------------------
// Граничные условия
// ------------------------------------------------------------
void computeWallFlux(const array<double,4>& UL, double nx, double ny, array<double,4>& flux) {
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

void computeInletFlux(const array<double,4>& UL, double nx, double ny, array<double,4>& flux) {
    double rhoL, uL, vL, pL;
    conservativeToPrimitive(UL, rhoL, uL, vL, pL);
    double alpha_rad = ALPHA_DEG * M_PI / 180.0;
    double uR = MACH_INF * SOUND_SPEED * cos(alpha_rad);
    double vR = MACH_INF * SOUND_SPEED * sin(alpha_rad);
    double rhoR = INLET_RHO;
    double pR   = INLET_P;
    double uNormL = uL*nx + vL*ny, uNormR = uR*nx + vR*ny;
    double uTangL = -uL*ny + vL*nx, uTangR = -uR*ny + vR*nx;
    auto [rhoFace, uNormFace, pFace] = exactRiemannSolver(rhoL, uNormL, pL, rhoR, uNormR, pR, 0.0);
    double uTangFace = uTangR;
    double uFace = uNormFace * nx - uTangFace * ny;
    double vFace = uNormFace * ny + uTangFace * nx;
    flux[0] = rhoFace * uNormFace;
    flux[1] = rhoFace * uNormFace * uFace + pFace * nx;
    flux[2] = rhoFace * uNormFace * vFace + pFace * ny;
    double EFace = pFace/(GAMMA-1.0) + 0.5*rhoFace*(uFace*uFace + vFace*vFace);
    flux[3] = (EFace + pFace) * uNormFace;
}

void computeOutletFlux(const array<double,4>& UL, double nx, double ny, array<double,4>& flux) {
    double rhoL, uL, vL, pL;
    conservativeToPrimitive(UL, rhoL, uL, vL, pL);
    double rhoR = rhoL, uR = uL, vR = vL, pR = OUTLET_P;
    double uNormL = uL*nx + vL*ny, uNormR = uR*nx + vR*ny;
    double uTangL = -uL*ny + vL*nx, uTangR = -uR*ny + vR*nx;
    auto [rhoFace, uNormFace, pFace] = exactRiemannSolver(rhoL, uNormL, pL, rhoR, uNormR, pR, 0.0);
    double uTangFace = (uNormFace > 0.0) ? uTangL : uTangR;
    double uFace = uNormFace * nx - uTangFace * ny;
    double vFace = uNormFace * ny + uTangFace * nx;
    flux[0] = rhoFace * uNormFace;
    flux[1] = rhoFace * uNormFace * uFace + pFace * nx;
    flux[2] = rhoFace * uNormFace * vFace + pFace * ny;
    double EFace = pFace/(GAMMA-1.0) + 0.5*rhoFace*(uFace*uFace + vFace*vFace);
    flux[3] = (EFace + pFace) * uNormFace;
}

// ------------------------------------------------------------
// Вычисление глобального шага по времени (исправлено!)
// ------------------------------------------------------------
double computeTimeStep(const vector<Triangle>& triangles, const vector<Edge>& edges,
                       const vector<array<double,4>>& U, double cfl) {
    double minArea = 1e30, maxSpeed = 0.0;
    for (size_t i = 0; i < triangles.size(); ++i) {
        double rho, u, v, p;
        conservativeToPrimitive(U[i], rho, u, v, p);
        double c = sqrt(GAMMA * p / rho);
        double speed = sqrt(u*u + v*v) + c;
        double area = triangles[i].area;
        if (area < minArea) minArea = area;
        if (speed > maxSpeed) maxSpeed = speed;
    }
    double maxEdgeLen = 0.0;
    for (const auto& e : edges) if (e.length > maxEdgeLen) maxEdgeLen = e.length;
    // Исправленный возврат
    return cfl * minArea / (maxSpeed * maxEdgeLen);
}

// ------------------------------------------------------------
// Обновление решения (схема Годунова)
// ------------------------------------------------------------
void updateSolution(vector<array<double,4>>& U, const vector<Triangle>& triangles,
                    const vector<Edge>& edges, double dt) {
    vector<array<double,4>> Unew = U;
    for (const auto& edge : edges) {
        const array<double,4>& UL = U[edge.cellA];
        array<double,4> flux;
        bool isBoundary = (edge.cellB == -1);
        if (isBoundary) {
            /*if (edge.bcTag == 1)
                computeWallFlux(UL, edge.nx, edge.ny, flux);
            else if (edge.bcTag == 2)
                computeInletFluxSimple(UL, edge.nx, edge.ny, flux);
            else if (edge.bcTag == 3)
                computeOutletFluxSimple(UL, edge.nx, edge.ny, flux);
            else
                computeWallFlux(UL, edge.nx, edge.ny, flux);*/
            if (edge.bcTag == 1)           // стенка (тег 1 в вашей сетке – wall, но в файле тег 1 для стенки)
                computeWallFlux(UL, edge.nx, edge.ny, flux);
            else if (edge.bcTag == 2)       // вход
                computeInletFlux(UL, edge.nx, edge.ny, flux);
            else if (edge.bcTag == 3)      // выход
                computeOutletFlux(UL, edge.nx, edge.ny, flux);
            else
                computeWallFlux(UL, edge.nx, edge.ny, flux);
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
    // Применяем обновление и обеспечиваем положительность
    for (size_t i = 0; i < U.size(); ++i) {
        double rho, u, v, p;
        conservativeToPrimitive(Unew[i], rho, u, v, p);
        if (rho < 1e-12) rho = 1e-12;
        if (p < 1e-12) p = 1e-12;
        primitiveToConservative(rho, u, v, p, U[i]);
    }
}

// ------------------------------------------------------------
// Запись VTK (добавлено поле числа Маха)
// ------------------------------------------------------------
void writeVTK(const string& baseName, const string& outputDir,
              const vector<Point>& points, const vector<Triangle>& triangles,
              const vector<array<double,4>>& U, double time, int fileCounter) {
    char fileName[512];
    sprintf(fileName, "%s/simulation/%s_%06d.vtk", outputDir.c_str(), baseName.c_str(), fileCounter);
    ofstream out(fileName);
    if (!out) { cerr << "Cannot create VTK file.\n"; return; }
    out << "# vtk DataFile Version 3.0\nMesh from Gmsh, time = " << time << "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << points.size() << " double\n";
    for (const auto& p : points) out << p.x << " " << p.y << " " << p.z << "\n";
    int numCells = triangles.size(), totalIndices = 0;
    for (const auto& t : triangles) totalIndices += 4;
    out << "CELLS " << numCells << " " << totalIndices << "\n";
    for (const auto& t : triangles) out << "3 " << t.nodes[0] << " " << t.nodes[1] << " " << t.nodes[2] << "\n";
    out << "CELL_TYPES " << numCells << "\n";
    for (size_t i = 0; i < triangles.size(); ++i) out << "5\n";
    out << "CELL_DATA " << numCells << "\n";
    // Плотность
    out << "SCALARS density double 1\nLOOKUP_TABLE default\n";
    for (const auto& Ucell : U) out << Ucell[0] << "\n";
    // Давление
    out << "SCALARS pressure double 1\nLOOKUP_TABLE default\n";
    for (const auto& Ucell : U) {
        double rho, u, v, p; conservativeToPrimitive(Ucell, rho, u, v, p);
        out << p << "\n";
    }
    // Скорость
    out << "VECTORS velocity double\n";
    for (const auto& Ucell : U) {
        double rho, u, v, p; conservativeToPrimitive(Ucell, rho, u, v, p);
        out << u << " " << v << " 0.0\n";
    }
    // Число Маха (добавлено)
    out << "SCALARS Mach double 1\nLOOKUP_TABLE default\n";
    for (const auto& Ucell : U) {
        double rho, u, v, p; conservativeToPrimitive(Ucell, rho, u, v, p);
        double c = sqrt(GAMMA * p / rho);
        double Mach = sqrt(u*u + v*v) / c;
        out << Mach << "\n";
    }
    out.close();
    cout << "VTK written: " << fileName << " at t=" << time << "\n";
}

// ------------------------------------------------------------
// Инициализация решения (исправлена: однородный поток с углом атаки)
// ------------------------------------------------------------
void initializeSolution(const vector<Triangle>& triangles, vector<array<double,4>>& U) {
    U.resize(triangles.size());
    double alpha_rad = ALPHA_DEG * M_PI / 180.0;
    double u_inf = MACH_INF * SOUND_SPEED * cos(alpha_rad);
    double v_inf = MACH_INF * SOUND_SPEED * sin(alpha_rad);
    double rho_inf = INLET_RHO;
    double p_inf   = INLET_P;
    for (size_t i = 0; i < triangles.size(); ++i) {
        primitiveToConservative(rho_inf, u_inf, v_inf, p_inf, U[i]);
    }
}
// ------------------------------------------------------------
// Сохранение Cp вдоль поверхности профиля
// ------------------------------------------------------------
void saveCpProfile(const vector<Point>& points,
                   const vector<Triangle>& triangles,
                   const vector<Edge>& edges,
                   const vector<array<double,4>>& U,
                   double rho_inf, double u_inf, double chord,
                   const string& filename) {
    ofstream out(filename);
    if (!out) {
        cerr << "Cannot create " << filename << endl;
        return;
    }
    out << "# x/c    y/c    Cp\n";
    double alpha_rad = ALPHA_DEG * M_PI / 180.0;
    double q_inf = 0.5 * rho_inf * (u_inf*u_inf + (MACH_INF*SOUND_SPEED*sin(alpha_rad))*(MACH_INF*SOUND_SPEED*sin(alpha_rad)));
    for (const auto& edge : edges) {
        if (edge.bcTag != 1) continue;   // только стенка профиля
        const Triangle& cell = triangles[edge.cellA];
        int i1 = cell.nodes[edge.edgeIdxA];
        int i2 = cell.nodes[(edge.edgeIdxA+1)%3];
        double mx = 0.5 * (points[i1].x + points[i2].x);
        double my = 0.5 * (points[i1].y + points[i2].y);
        double rho, u, v, p;
        conservativeToPrimitive(U[edge.cellA], rho, u, v, p);
        double Cp = (p - INLET_P) / q_inf;
        out << mx/chord << " " << my/chord << " " << Cp << "\n";
    }
    out.close();
    cout << "Cp profile saved to " << filename << endl;
}
// ------------------------------------------------------------
// Вычисление аэродинамических коэффициентов (исправлена фильтрация)
// ------------------------------------------------------------
void computeAerodynamicCoefficients(const vector<Point>& points,
                                    const vector<Triangle>& triangles,
                                    const vector<Edge>& edges,
                                    const vector<array<double,4>>& U,
                                    double rho_inf, double u_inf, double chord,
                                    double& CL, double& CD) {
    double Fx = 0.0, Fy = 0.0;
    for (const auto& edge : edges) {
        // Физическая группа 1 соответствует стенке (профилю)
        if (edge.bcTag != 1) continue;
        const Triangle& cell = triangles[edge.cellA];
        int i1 = cell.nodes[edge.edgeIdxA];
        int i2 = cell.nodes[(edge.edgeIdxA + 1) % 3];
        double mx = 0.5 * (points[i1].x + points[i2].x);
        double my = 0.5 * (points[i1].y + points[i2].y);
        // Дополнительная фильтрация не требуется, оставим только для контроля
        // (можно удалить, но оставим для уверенности, что не учитываются рёбра вне профиля)
        if (mx < -0.2 || mx > 1.2 || fabs(my) > 0.3) continue;
        double rho, u, v, p;
        conservativeToPrimitive(U[edge.cellA], rho, u, v, p);
        double L = edge.length;
        Fx += -p * edge.nx * L;
        Fy += -p * edge.ny * L;
    }
    double alpha_rad = ALPHA_DEG * M_PI / 180.0;
    double q = 0.5 * rho_inf * (u_inf*u_inf + (MACH_INF*SOUND_SPEED*sin(alpha_rad))*(MACH_INF*SOUND_SPEED*sin(alpha_rad))) * chord;
    CL = (Fy * cos(alpha_rad) - Fx * sin(alpha_rad)) / q;
    CD = (Fx * cos(alpha_rad) + Fy * sin(alpha_rad)) / q;
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " file.msh\n";
        return 1;
    }
    // Вычисление скорости звука и компонент невозмущённого потока
    SOUND_SPEED = sqrt(GAMMA * INLET_P / INLET_RHO);
    double alpha_rad = ALPHA_DEG * M_PI / 180.0;
    INLET_U = MACH_INF * SOUND_SPEED * cos(alpha_rad);
    INLET_V = MACH_INF * SOUND_SPEED * sin(alpha_rad);

    vector<Point> points;
    vector<Triangle> triangles;
    vector<BoundaryLine> boundaryLines;
    if (!readGmshMesh(argv[1], points, triangles, boundaryLines)) return 1;
    cout << "Points: " << points.size() << ", Triangles: " << triangles.size() << "\n";

    vector<Edge> edges;
    buildEdges(points, triangles, boundaryLines, edges);
    cout << "Edges: " << edges.size() << "\n";

    string inputPath = argv[1];
    size_t lastSlash = inputPath.find_last_of("/\\");
    string fileName = (lastSlash == string::npos) ? inputPath : inputPath.substr(lastSlash+1);
    size_t dot = fileName.rfind('.');
    string baseName = (dot != string::npos) ? fileName.substr(0, dot) : fileName;

    string resultsDir = "results";
    mkdir(resultsDir.c_str(), 0755);
    mkdir((resultsDir + "/simulation").c_str(), 0755);
    mkdir((resultsDir + "/mesh_quality").c_str(), 0755);

    // Инициализация решения (однородный поток)
    vector<array<double,4>> U;
    initializeSolution(triangles, U);

    double t = 0.0;
    int iter = 0;
    int fileCounter = 0;
    double nextOutput = OUTPUT_INTERVAL;
    cout << "Starting time integration with exact Riemann solver...\n";
    while (t < FINAL_TIME) {
        double dt = computeTimeStep(triangles, edges, U, CFL);
        if (t + dt > FINAL_TIME) dt = FINAL_TIME - t;
        updateSolution(U, triangles, edges, dt);
        t += dt;
        ++iter;
        if (t >= nextOutput) {
            writeVTK(baseName, resultsDir, points, triangles, U, t, fileCounter);
            ++fileCounter;
            nextOutput += OUTPUT_INTERVAL;
        }
        if (iter % 100 == 0) {
            cout << "Iter " << iter << ", t = " << t << ", dt = " << dt << "\n";
// После каждого шага (или каждые 100 итераций) выводите:
double rho_min = 1e30, rho_max = -1e30;
for (const auto& Ucell : U) {
    if (Ucell[0] < rho_min) rho_min = Ucell[0];
    if (Ucell[0] > rho_max) rho_max = Ucell[0];
}
cout << "rho: [" << rho_min << ", " << rho_max << "]" << endl;

        }
    }
    // Финальный вывод
    writeVTK(baseName, resultsDir, points, triangles, U, t, fileCounter);

    double CL, CD;
    computeAerodynamicCoefficients(points, triangles, edges, U,
                                   INLET_RHO, INLET_U, CHORD, CL, CD);
    cout << "\n=== AERODYNAMIC COEFFICIENTS ===\n";
    cout << "CL = " << CL << "\n";
    cout << "CD = " << CD << "\n";

    string coeffFile = resultsDir + "/" + baseName + "_coefficients.dat";
    ofstream fcoeff(coeffFile);
    fcoeff << "# alpha = " << ALPHA_DEG << " deg\n";
    fcoeff << "# CL = " << CL << "\n";
    fcoeff << "# CD = " << CD << "\n";
    fcoeff.close();
    cout << "Coefficients saved to " << coeffFile << "\n";
// В конце main():
string cpFile = resultsDir + "/" + baseName + "_cp_profile.dat";
saveCpProfile(points, triangles, edges, U, INLET_RHO, INLET_U, CHORD, cpFile);
    return 0;
}