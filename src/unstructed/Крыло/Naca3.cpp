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
#include <filesystem>

namespace fs = std::filesystem;
using namespace std;

const double GAMMA = 1.4;
const double CFL = 0.8;
const double OUTPUT_INTERVAL = 0.1;
const double CHORD = 1.0;

// Глобальные параметры набегающего потока (устанавливаются в main)
double INLET_U = 0.0;
double INLET_V = 0.0;
double INLET_P = 1.0;
double INLET_RHO = 1.0;
double OUTLET_P = 1.0;
double REYNOLDS = 1.0e6;
double ALPHA_DEG = 0.0;

// ------------------------------------------------------------
// Утилиты для работы с директориями
// ------------------------------------------------------------
void cleanDirectory(const string& dirPath) {
    struct stat info;
    if (stat(dirPath.c_str(), &info) != 0) return;
    try {
        for (const auto& entry : fs::directory_iterator(dirPath))
            fs::remove_all(entry.path());
        cout << "Cleaned directory: " << dirPath << "\n";
    } catch (const fs::filesystem_error& e) {
        cerr << "Warning: Could not clean directory " << dirPath << ": " << e.what() << "\n";
    }
}

bool createDirectory(const string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
#ifdef _WIN32
        int ret = _mkdir(path.c_str());
#else
        int ret = mkdir(path.c_str(), 0755);
#endif
        if (ret != 0) {
            cerr << "Failed to create directory: " << path << "\n";
            return false;
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        cerr << path << " exists but is not a directory.\n";
        return false;
    }
    return true;
}

// ------------------------------------------------------------
// Точный решатель Римана
// ------------------------------------------------------------
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
        // x_over_t >= u_star_val
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

// ------------------------------------------------------------
// Структуры сетки
// ------------------------------------------------------------
struct Point { double x, y, z; };
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

// ------------------------------------------------------------
// Чтение сетки Gmsh (поддержка MSH2 и старого формата)
// ------------------------------------------------------------
bool readGmshMesh(const string& filename,
                  vector<Point>& points,
                  vector<Triangle>& triangles,
                  vector<BoundaryLine>& boundaryLines) {
    ifstream in(filename);
    if (!in) { cerr << "Cannot open input file: " << filename << "\n"; return false; }
    string line;

    auto findSection = [&](const string& section) -> bool {
        while (getline(in, line)) {
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            if (line.find(section) != string::npos) return true;
        }
        return false;
    };
    auto parseInts = [](const string& s) {
        vector<int> v; stringstream ss(s); int x;
        while (ss >> x) v.push_back(x);
        return v;
    };

    // -- $MeshFormat --
    if (!findSection("$MeshFormat")) { cerr << "No $MeshFormat\n"; return false; }
    getline(in, line); // версия

    // -- $PhysicalNames (пропускаем) --
    if (findSection("$PhysicalNames")) {
        getline(in, line); int n; stringstream(line) >> n;
        for (int i = 0; i < n; ++i) getline(in, line);
    }

    // -- $Nodes --
    if (!findSection("$Nodes")) { cerr << "No $Nodes\n"; return false; }
    getline(in, line);
    vector<int> hdr = parseInts(line);
    int numNodes, minNodeTag, maxNodeTag;
    bool oldNodes = (hdr.size() == 1);
    if (oldNodes) {
        numNodes = hdr[0];
        minNodeTag = numeric_limits<int>::max();
        maxNodeTag = 0;
    } else if (hdr.size() == 4) {
        numNodes = hdr[1];
        minNodeTag = hdr[2];
        maxNodeTag = hdr[3];
    } else { cerr << "Bad $Nodes header\n"; return false; }

    points.resize(numNodes);
    vector<int> nodeMap, nodeTags(numNodes);

    if (!oldNodes) {
        nodeMap.assign(maxNodeTag - minNodeTag + 1, -1);
        int idx = 0, numBlocks = hdr[0];
        for (int b = 0; b < numBlocks; ++b) {
            int dim, tag, param, cnt; in >> dim >> tag >> param >> cnt;
            vector<int> t(cnt); for (int i=0;i<cnt;++i) in>>t[i];
            for (int i=0;i<cnt;++i) {
                double x,y,z; in>>x>>y>>z;
                int local = t[i] - minNodeTag;
                nodeMap[local] = idx;
                points[idx++] = {x,y,z};
            }
        }
    } else {
        // старый формат: сразу тег x y z
        for (int i=0;i<numNodes;++i) {
            int t; double x,y,z; in>>t>>x>>y>>z;
            nodeTags[i] = t;
            points[i] = {x,y,z};
            if (t < minNodeTag) minNodeTag = t;
            if (t > maxNodeTag) maxNodeTag = t;
        }
        nodeMap.assign(maxNodeTag - minNodeTag + 1, -1);
        for (int i=0;i<numNodes;++i) nodeMap[nodeTags[i] - minNodeTag] = i;
    }

    // -- $Elements --
    if (!findSection("$Elements")) { cerr << "No $Elements\n"; return false; }
    getline(in, line); hdr = parseInts(line);
    int numElements, minElemTag, maxElemTag;
    bool oldElems = (hdr.size() == 1);
    vector<array<int,3>> triNodes;
    boundaryLines.clear();

    if (!oldElems) {
        int numBlocks = hdr[0]; numElements = hdr[1]; minElemTag = hdr[2]; maxElemTag = hdr[3];
        for (int b=0;b<numBlocks;++b) {
            int dim, physTag, etype, cnt; in>>dim>>physTag>>etype>>cnt;
            if (etype == 2) {
                for (int i=0;i<cnt;++i) {
                    int id, n1,n2,n3; in>>id>>n1>>n2>>n3;
                    int i1=nodeMap[n1-minNodeTag], i2=nodeMap[n2-minNodeTag], i3=nodeMap[n3-minNodeTag];
                    if (i1>=0&&i2>=0&&i3>=0) triNodes.push_back({i1,i2,i3});
                }
            } else if (etype == 1) {
                for (int i=0;i<cnt;++i) {
                    int id, n1,n2; in>>id>>n1>>n2;
                    int i1=nodeMap[n1-minNodeTag], i2=nodeMap[n2-minNodeTag];
                    boundaryLines.push_back({physTag, i1, i2});
                }
            } else { in.ignore(100000, '\n'); }
        }
    } else {
        numElements = hdr[0];
        for (int i=0;i<numElements;++i) {
            int id, etype; in>>id>>etype;
            if (etype == 2) {
                int n1,n2,n3; in>>n1>>n2>>n3;
                int i1=nodeMap[n1-minNodeTag], i2=nodeMap[n2-minNodeTag], i3=nodeMap[n3-minNodeTag];
                if (i1>=0&&i2>=0&&i3>=0) triNodes.push_back({i1,i2,i3});
            } else if (etype == 1) {
                int n1,n2; in>>n1>>n2;
                int i1=nodeMap[n1-minNodeTag], i2=nodeMap[n2-minNodeTag];
                // физическая группа неизвестна – примем 1 (стенка), но здесь нужно особое внимание,
                // потому что без физических групп все граничные линии будут стенкой.
                // Однако в вашем файле есть $PhysicalNames, но старый формат не привязывает теги к линиям.
                // Поэтому старый формат НЕ подходит, если нужны корректные bcTag.
                boundaryLines.push_back({1, i1, i2}); // <-- это ошибка!
            }
        }
    }

    in.close();

    // треугольники
    triangles.clear();
    for (auto& t : triNodes) {
        double x1=points[t[1]].x-points[t[0]].x, y1=points[t[1]].y-points[t[0]].y;
        double x2=points[t[2]].x-points[t[0]].x, y2=points[t[2]].y-points[t[0]].y;
        double area = 0.5*(x1*y2 - x2*y1);
        if (area < 0) { swap(t[1], t[2]); area = -area; }
        Triangle tr;
        tr.nodes[0]=t[0]; tr.nodes[1]=t[1]; tr.nodes[2]=t[2];
        tr.area = area;
        tr.centroid[0] = (points[t[0]].x+points[t[1]].x+points[t[2]].x)/3.0;
        tr.centroid[1] = (points[t[0]].y+points[t[1]].y+points[t[2]].y)/3.0;
        triangles.push_back(tr);
    }
    return true;
}// ------------------------------------------------------------
// Построение рёбер
// ------------------------------------------------------------
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
            if (it != lineTagMap.end()) edge.bcTag = it->second;
            else edge.bcTag = 1;
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

// ------------------------------------------------------------
// Преобразования консервативных переменных
// ------------------------------------------------------------
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

// ------------------------------------------------------------
// Вычисление потока через грань (точный решатель Римана)
// ------------------------------------------------------------
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

// ------------------------------------------------------------
// Граничные условия
// ------------------------------------------------------------
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

void computeFarfieldFlux(const array<double,4>& UL, double nx, double ny,
                         double rho_inf, double u_inf, double v_inf, double p_inf,
                         array<double,4>& flux) {
    double rhoL, uL, vL, pL;
    conservativeToPrimitive(UL, rhoL, uL, vL, pL);
    double uNormL = uL*nx + vL*ny;
    double uNormR = u_inf*nx + v_inf*ny;
    double uTangL = -uL*ny + vL*nx;
    double uTangR = -u_inf*ny + v_inf*nx;
    auto [rhoFace, uNormFace, pFace] = exactRiemannSolver(rhoL, uNormL, pL,
                                                          rho_inf, uNormR, p_inf,
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

// ------------------------------------------------------------
// Расчёт шага по времени
// ------------------------------------------------------------
double computeTimeStep(const vector<Triangle>& triangles,
                       const vector<Edge>& edges,
                       const vector<array<double,4>>& U,
                       double cfl) {
    double maxSpeed = 0.0;
    for (size_t i = 0; i < U.size(); ++i) {
        double rho, u, v, p;
        conservativeToPrimitive(U[i], rho, u, v, p);
        if (rho <= 0 || p <= 0) continue; // игнорируем нефизичные ячейки
        double c = sqrt(GAMMA * p / rho);
        double speed = sqrt(u*u + v*v) + c;
        if (speed > maxSpeed) maxSpeed = speed;
    }
    if (maxSpeed <= 0) return 0.0;

    double minArea = numeric_limits<double>::max();
    for (const auto& tri : triangles) {
        if (tri.area > 1e-12 && tri.area < minArea) minArea = tri.area;
    }
    if (minArea <= 0) return 0.0;

    // Характерный размер ячейки – корень из площади (консервативная оценка)
    double minCellSize = sqrt(minArea);
    double dt = cfl * minCellSize / maxSpeed;
    return dt;
}
// ------------------------------------------------------------
// Обновление решения по схеме Годунова
// ------------------------------------------------------------
void updateSolution(vector<array<double,4>>& U,
                    const vector<Triangle>& triangles,
                    const vector<Edge>& edges,
                    double dt,
                    double rho_inf, double u_inf, double v_inf, double p_inf) {
    vector<array<double,4>> Unew = U;
    for (const auto& edge : edges) {
        const array<double,4>& UL = U[edge.cellA];
        array<double,4> flux;
        bool isBoundary = (edge.cellB == -1);
        if (isBoundary) {
            // bcTag: 1 – стенка (профиль), 2 – дальнее поле (внешние границы)
            if (edge.bcTag == 1) {
                computeWallFlux(UL, edge.nx, edge.ny, flux);
            } else {
                computeFarfieldFlux(UL, edge.nx, edge.ny,
                                    rho_inf, u_inf, v_inf, p_inf, flux);
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

// ------------------------------------------------------------
// Запись VTK (добавлено число Маха)
// ------------------------------------------------------------
void writeVTK(const string& baseName, const string& outputDir,
              const vector<Point>& points,
              const vector<Triangle>& triangles,
              const vector<array<double,4>>& U,
              double time, int fileCounter) {
    char fileName[512];
    sprintf(fileName, "%s/simulation/%s_%06d.vtk", outputDir.c_str(), baseName.c_str(), fileCounter);
    ofstream out(fileName);
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

    out << "SCALARS Mach double 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < triangles.size(); ++i) {
        double rho, u, v, p;
        conservativeToPrimitive(U[i], rho, u, v, p);
        double c = sqrt(GAMMA * p / rho);
        out << sqrt(u*u + v*v) / c << "\n";
    }

    out.close();
    cout << "VTK file written: " << fileName << " at time " << time << "\n";
}

// ------------------------------------------------------------
// Вычисление аэродинамических коэффициентов
// ------------------------------------------------------------
void computeAerodynamicCoefficients(
    const vector<Point>& points,
    const vector<Triangle>& triangles,
    const vector<Edge>& edges,
    const vector<array<double,4>>& U,
    double alpha_deg,
    double rho_inf,
    double u_inf,
    double chord,
    double& CL,
    double& CD) {
    double Fx = 0.0, Fy = 0.0;
    double q_inf = 0.5 * rho_inf * u_inf * u_inf;
    for (const auto& edge : edges) {
        if (edge.bcTag != 1) continue;   // только стенка (профиль)
        const Triangle& cell = triangles[edge.cellA];
        int i1 = cell.nodes[edge.edgeIdxA];
        int i2 = cell.nodes[(edge.edgeIdxA+1)%3];
        double mx = 0.5 * (points[i1].x + points[i2].x);
        double my = 0.5 * (points[i1].y + points[i2].y);
        // Фильтр по координатам: профиль примерно от 0 до 1 по x, y около 0
        if (mx < -0.2 || mx > 1.2 || fabs(my) > 0.3) continue;
        double rho, u, v, p;
        conservativeToPrimitive(U[edge.cellA], rho, u, v, p);
        Fx += -p * edge.nx * edge.length;
        Fy += -p * edge.ny * edge.length;
    }
    double alpha = alpha_deg * M_PI / 180.0;
    double q = q_inf * chord;
    CL = (Fy * cos(alpha) - Fx * sin(alpha)) / q;
    CD = (Fx * cos(alpha) + Fy * sin(alpha)) / q;
}

// ------------------------------------------------------------
// Сбор данных о распределении Cp по поверхности
// ------------------------------------------------------------
void collectSurfacePressure(const vector<Point>& points,
                            const vector<Triangle>& triangles,
                            const vector<Edge>& edges,
                            const vector<array<double,4>>& U,
                            double rho_inf, double u_inf, double p_inf,
                            const string& caseName,
                            double alpha_deg,
                            const string& outputDir) {
    double q_inf = 0.5 * rho_inf * u_inf * u_inf;
    string cpFileName = outputDir + "/cp_" + caseName + "_alpha=" + to_string(alpha_deg) + ".dat";
    ofstream cpFile(cpFileName);
    if (!cpFile) { cerr << "Cannot create " << cpFileName << "\n"; return; }
    cpFile << "# x y Cp\n";
    for (const auto& edge : edges) {
        if (edge.bcTag != 1) continue;   // только профиль
        const Triangle& cell = triangles[edge.cellA];
        int i1 = cell.nodes[edge.edgeIdxA];
        int i2 = cell.nodes[(edge.edgeIdxA+1)%3];
        double mx = 0.5 * (points[i1].x + points[i2].x);
        double my = 0.5 * (points[i1].y + points[i2].y);
        double rho, u, v, p;
        conservativeToPrimitive(U[edge.cellA], rho, u, v, p);
        double cp = (p - p_inf) / q_inf;
        cpFile << mx << " " << my << " " << cp << "\n";
    }
    cpFile.close();
    cout << "Surface pressure data saved to " << cpFileName << "\n";
}

// ------------------------------------------------------------
// Инициализация решения (дальнее поле)
// ------------------------------------------------------------
void initializeSolution(const vector<Triangle>& triangles,
                        vector<array<double,4>>& U,
                        double rho_inf, double u_inf, double v_inf, double p_inf) {
    U.resize(triangles.size());
    for (size_t i = 0; i < triangles.size(); ++i)
        primitiveToConservative(rho_inf, u_inf, v_inf, p_inf, U[i]);
}

// ------------------------------------------------------------
// Главная функция
// ------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " mesh.msh test_case alpha [final_time]\n";
        cout << "  test_case: naca0012 or cylinder\n";
        cout << "  alpha: angle of attack in degrees\n";
        return 1;
    }
    string meshFile = argv[1];
    string testCase = argv[2];
    double alpha_deg = atof(argv[3]);
    double finalTime = (argc >= 5) ? atof(argv[4]) : 10.0;

    // Установка параметров набегающего потока
    double M_inf, rho_inf, p_inf, u_inf, v_inf;
    if (testCase == "naca0012") {
        M_inf = 0.8;
        rho_inf = 1.0;
        p_inf = 1.0;                       // давление 1.0 -> c = sqrt(1.4) ? 1.183
        u_inf = M_inf * cos(alpha_deg * M_PI / 180.0);
        v_inf = M_inf * sin(alpha_deg * M_PI / 180.0);
    } else if (testCase == "cylinder") {
        M_inf = 0.3;
        rho_inf = 1.0;
        p_inf = 1.0 / GAMMA;               // для цилиндра оставим историческую нормировку (c=1)
        u_inf = M_inf;
        v_inf = 0.0;
        alpha_deg = 0.0;
    } else {
        cerr << "Unknown test case: " << testCase << "\n";
        return 1;
    }

    // Сохраняем глобальные переменные (если потребуется)
    INLET_U = u_inf;
    INLET_V = v_inf;
    INLET_P = p_inf;
    INLET_RHO = rho_inf;
    ALPHA_DEG = alpha_deg;

    // Чтение сетки
    vector<Point> points;
    vector<Triangle> triangles;
    vector<BoundaryLine> boundaryLines;
    if (!readGmshMesh(meshFile, points, triangles, boundaryLines)) return 1;
    cout << "Points: " << points.size() << ", Triangles: " << triangles.size() << "\n";

    vector<Edge> edges;
    buildEdges(points, triangles, boundaryLines, edges);
    cout << "Edges: " << edges.size() << "\n";

    // Создание директорий
    string resultsDir = "results";
    createDirectory(resultsDir);
    createDirectory(resultsDir + "/simulation");
    cleanDirectory(resultsDir + "/simulation");

    // Инициализация решения
    vector<array<double,4>> U;
    initializeSolution(triangles, U, rho_inf, u_inf, v_inf, p_inf);

    double t = 0.0;
    int iter = 0;
    double nextOutput = OUTPUT_INTERVAL;
    int fileCounter = 0;

    cout << "Starting Godunov time integration...\n";
    cout << "Case: " << testCase << ", M=" << M_inf << ", alpha=" << alpha_deg << " deg\n";
    while (t < finalTime) {
        double dt = computeTimeStep(triangles, edges, U, CFL);
        if (t + dt > finalTime) dt = finalTime - t;
        updateSolution(U, triangles, edges, dt, rho_inf, u_inf, v_inf, p_inf);
        t += dt;
        iter++;
        if (t >= nextOutput) {
            string baseName = testCase + "_alpha=" + to_string(alpha_deg);
            writeVTK(baseName, resultsDir, points, triangles, U, t, fileCounter);
            fileCounter++;
            nextOutput += OUTPUT_INTERVAL;
        }
        if (iter % 100 == 0)
            cout << "Iteration " << iter << ", time = " << t << ", dt = " << dt << "\n";
    }
    // Финальный VTK
    string baseName = testCase + "_alpha=" + to_string(alpha_deg);
    writeVTK(baseName, resultsDir, points, triangles, U, t, fileCounter);

    // Расчёт коэффициентов
    double CL, CD;
    computeAerodynamicCoefficients(points, triangles, edges, U,
                                   alpha_deg, rho_inf, u_inf, CHORD, CL, CD);
    cout << "\n=== AERODYNAMIC COEFFICIENTS ===\n";
    cout << "CL = " << CL << "\n";
    cout << "CD = " << CD << "\n";

    string coeffFile = resultsDir + "/coefficients_" + testCase + "_alpha=" + to_string(alpha_deg) + ".dat";
    ofstream fcoeff(coeffFile);
    fcoeff << "# Test case: " << testCase << "\n";
    fcoeff << "# M_inf = " << M_inf << ", alpha = " << alpha_deg << " deg\n";
    fcoeff << "# CL = " << CL << "\n";
    fcoeff << "# CD = " << CD << "\n";
    fcoeff.close();

    // Сохранение распределения Cp
    collectSurfacePressure(points, triangles, edges, U,
                           rho_inf, u_inf, p_inf, testCase, alpha_deg, resultsDir);

    cout << "Done.\n";
    return 0;
}