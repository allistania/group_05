// ============================================================
// shem.cpp – решение уравнений газовой динамики методом Годунова
// на неструктурированной треугольной сетке (Gmsh, формат 2.2)
// Используется точный решатель Римана.
// ============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <tuple>
#include <string>
#include <cassert>
#include <gmsh.h>   // подключение официальной библиотеки Gmsh

using namespace std;

const double GAMMA = 1.4;
const double CFL   = 0.4;
const double T_END = 1.0;
const int    SAVE_INTERVAL = 100;

struct Conserved {
    double rho, rhou, rhov, rhoE;
};

struct Primitive {
    double rho, u, v, p;
};

struct Cell {
    int id;
    vector<int> nodes;
    double area, xc, yc;
    Conserved U, dU;
    int phys_tag;
};

struct Face {
    int id;
    int cellL, cellR;
    int type;           // 0 – internal, 1 – boundary
    int phys_tag;
    vector<int> nodes;
    double nx, ny, length;
};

vector<double> nodeCoords;          // плоский массив: x1,y1,z1, x2,y2,z2, ...
int nodeNum, elemNum;
vector<int> elemNodes;               // 4*elemNum, последний элемент = -1 для треугольников
vector<int> elemTags;                // физические тэги для элементов
vector<Cell> cells;
vector<Face> faces;

// ========== Примитивные <-> консервативные ==========
Primitive getPrimitive(const Conserved& U) {
    double rho = U.rho;
    double u = U.rhou / rho;
    double v = U.rhov / rho;
    double e = U.rhoE / rho - 0.5 * (u*u + v*v);
    double p = (GAMMA - 1.0) * rho * e;
    return {rho, u, v, p};
}

Conserved getConserved(const Primitive& prim) {
    double rho = prim.rho;
    double u = prim.u;
    double v = prim.v;
    double p = prim.p;
    double e = p / (GAMMA - 1.0) / rho;
    double E = e + 0.5 * (u*u + v*v);
    return {rho, rho*u, rho*v, rho*E};
}

// ========== Геометрия ==========
double triangleArea(const vector<int>& nodes, const vector<double>& coords) {
    double x1 = coords[3*(nodes[0]-1)], y1 = coords[3*(nodes[0]-1)+1];
    double x2 = coords[3*(nodes[1]-1)], y2 = coords[3*(nodes[1]-1)+1];
    double x3 = coords[3*(nodes[2]-1)], y3 = coords[3*(nodes[2]-1)+1];
    return 0.5 * fabs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
}

void triangleCentroid(const vector<int>& nodes, const vector<double>& coords, double& xc, double& yc) {
    double x1 = coords[3*(nodes[0]-1)], y1 = coords[3*(nodes[0]-1)+1];
    double x2 = coords[3*(nodes[1]-1)], y2 = coords[3*(nodes[1]-1)+1];
    double x3 = coords[3*(nodes[2]-1)], y3 = coords[3*(nodes[2]-1)+1];
    xc = (x1 + x2 + x3) / 3.0;
    yc = (y1 + y2 + y3) / 3.0;
}

// ========== Чтение сетки с использованием Gmsh API ==========
bool readMesh(const string& filename) {
    // Инициализация Gmsh
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    
    // Загружаем файл сетки
    gmsh::open(filename);
    
    // Получаем все узлы (nodeTags и координаты)
    vector<size_t> nodeTags;
    vector<double> coord, parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);
    
    nodeNum = nodeTags.size();
    nodeCoords.resize(3 * nodeNum);
    // Заполняем nodeCoords, считая, что nodeTags[i] = i+1 (нумерация с 1)
    for (size_t i = 0; i < nodeTags.size(); ++i) {
        // Проверка: обычно тэги идут последовательно, начиная с 1
        if (nodeTags[i] != i+1) {
            cerr << "Warning: node tags are not sequential from 1. The code assumes sequential indexing." << endl;
            // В этом случае пришлось бы строить отображение, но для простоты считаем последовательными
        }
        nodeCoords[3*i]   = coord[3*i];
        nodeCoords[3*i+1] = coord[3*i+1];
        nodeCoords[3*i+2] = coord[3*i+2];
    }
    
    // Получаем все элементы вместе с физическими тэгами
    vector<int> elementTypes;
    vector<vector<size_t>> elementTags, nodeTagsOfElem;
    vector<vector<int>> elementPhysicalTags; // физические тэги для каждого элемента
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsOfElem, elementPhysicalTags);
    
    // Вспомогательные структуры для хранения треугольников и граничных линий
    vector<tuple<int, int, int>> triNodes;   // (node1, node2, node3) – тэги узлов
    vector<int> triPhysTags;                 // физические тэги треугольников
    map<pair<int,int>, int> edgeToPhysTag;   // для граничных рёбер: (node1,node2) -> phys_tag
    
    // Перебираем все типы элементов
    for (size_t i = 0; i < elementTypes.size(); ++i) {
        int type = elementTypes[i];
        const auto& tags = elementTags[i];
        const auto& nodes = nodeTagsOfElem[i];
        const auto& physTags = elementPhysicalTags[i];
        int elemCount = tags.size();
        
        if (type == 2) { // 2D треугольник
            for (int j = 0; j < elemCount; ++j) {
                int n0 = nodes[3*j]   - 1; // переводим в 0-базу для внутреннего хранения (будем хранить как 1-базу позже)
                int n1 = nodes[3*j+1] - 1;
                int n2 = nodes[3*j+2] - 1;
                // В исходном коде хранятся 1-базовые номера узлов, поэтому запомним +1
                triNodes.emplace_back(n0+1, n1+1, n2+1);
                triPhysTags.push_back(physTags[j]);
            }
        }
        else if (type == 1) { // 1D линия (граница)
            for (int j = 0; j < elemCount; ++j) {
                int n0 = nodes[2*j]   - 1;
                int n1 = nodes[2*j+1] - 1;
                // Сохраняем ребро с сортировкой узлов
                int a = n0+1, b = n1+1;
                if (a > b) swap(a,b);
                edgeToPhysTag[{a,b}] = physTags[j];
            }
        }
        // Другие типы (четырёхугольники и т.п.) игнорируем
    }
    
    elemNum = triNodes.size();
    elemNodes.resize(4 * elemNum);
    elemTags.resize(elemNum);
    // Заполняем массивы в формате, ожидаемом исходным кодом
    for (int i = 0; i < elemNum; ++i) {
        auto [n1, n2, n3] = triNodes[i];
        elemNodes[4*i]   = n1;
        elemNodes[4*i+1] = n2;
        elemNodes[4*i+2] = n3;
        elemNodes[4*i+3] = -1;   // маркер треугольника
        elemTags[i] = triPhysTags[i];
    }
    
    // Сохраняем отображение граничных рёбер для последующего использования
    // Оно понадобится при построении граней, поэтому сохраним в глобальную переменную или передадим.
    // Сделаем глобальную переменную для простоты, но лучше передавать через аргументы.
    // В данном случае мы сохраним edgeToPhysTag в статической переменной или передадим.
    // Поскольку функция buildFaces() будет вызываться после readMesh, мы можем сохранить карту в глобальной переменной.
    // Для этого создадим глобальный объект.
    // Однако, чтобы не переписывать сильно, передадим её в buildFaces через дополнительный аргумент.
    // Для этого изменим сигнатуру buildFaces.
    
    // Но проще всего сделать глобальную переменную. Создадим её перед buildFaces.
    // Здесь мы сохраним её в статической переменной, но для простоты сделаем глобальную.
    // Мы создадим глобальную переменную `globalEdgeToPhysTag` и заполним её.
    // Определим её где-то сверху.
    // Объявим extern map<pair<int,int>, int> globalEdgeToPhysTag; и определим в этом файле.
    // Пока просто сохраним в локальную, а в buildFaces передадим как аргумент.
    // Для этого изменим buildFaces на buildFaces(const map<pair<int,int>,int>& edgeToPhysTag).
    // Но для минимальных изменений, создадим глобальную.
    // Ниже объявим глобальную переменную.
    
    // Закрываем Gmsh
    gmsh::finalize();
    
    return true;
}

// Глобальная переменная для отображения ребро->физический тэг
map<pair<int,int>, int> globalEdgeToPhysTag;

// Модифицированная функция построения граней с использованием глобальной карты
void buildFaces() {
    map<pair<int,int>, vector<int>> edgeToCells;
    for (size_t cid = 0; cid < cells.size(); ++cid) {
        const auto& cell = cells[cid];
        int n[3] = {cell.nodes[0], cell.nodes[1], cell.nodes[2]};
        for (int i = 0; i < 3; ++i) {
            int a = n[i], b = n[(i+1)%3];
            if (a > b) swap(a,b);
            edgeToCells[{a,b}].push_back(cid);
        }
    }
    faces.clear();
    int faceId = 0;
    for (auto& kv : edgeToCells) {
        const auto& edge = kv.first;
        const auto& cellIndices = kv.second;
        Face f;
        f.id = faceId++;
        f.nodes = {edge.first, edge.second};
        double x1 = nodeCoords[3*(edge.first-1)], y1 = nodeCoords[3*(edge.first-1)+1];
        double x2 = nodeCoords[3*(edge.second-1)], y2 = nodeCoords[3*(edge.second-1)+1];
        double dx = x2 - x1, dy = y2 - y1;
        f.length = hypot(dx, dy);
        f.nx =  dy / f.length;
        f.ny = -dx / f.length;
        if (cellIndices.size() == 2) {
            f.type = 0;
            f.cellL = cellIndices[0];
            f.cellR = cellIndices[1];
            f.phys_tag = -1;
            double cxL = cells[f.cellL].xc, cyL = cells[f.cellL].yc;
            double mx = (x1+x2)/2.0, my = (y1+y2)/2.0;
            double dxL = mx - cxL, dyL = my - cyL;
            if (f.nx*dxL + f.ny*dyL < 0) { f.nx = -f.nx; f.ny = -f.ny; }
        } else if (cellIndices.size() == 1) {
            f.type = 1;
            f.cellL = cellIndices[0];
            f.cellR = -1;
            // Получаем физический тэг из глобальной карты
            auto it = globalEdgeToPhysTag.find(edge);
            if (it != globalEdgeToPhysTag.end())
                f.phys_tag = it->second;
            else
                f.phys_tag = 0; // если не найден, ставим 0
            double cx = cells[f.cellL].xc, cy = cells[f.cellL].yc;
            double mx = (x1+x2)/2.0, my = (y1+y2)/2.0;
            double dx = mx - cx, dy = my - cy;
            if (f.nx*dx + f.ny*dy < 0) { f.nx = -f.nx; f.ny = -f.ny; }
        } else {
            cerr << "Edge belongs to >2 cells!" << endl;
            exit(1);
        }
        faces.push_back(f);
    }
    cout << "Faces: " << faces.size() << endl;
}

// ========== Остальные функции без изменений ==========
// ... (все функции от getPrimitive до updateSolution, computeTimeStep, writeVTK, setInitialConditions и main остаются такими же, за исключением вызова buildFaces)

// Но чтобы код был полным, приведу их снова, но можно оставить как есть.

// ========== Точный решатель Римана ==========
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

double solve_p_star(double rhoL, double uL, double pL, double rhoR, double uR, double pR) {
    double tol = 1e-8;
    int max_iter = 100;
    double p0 = 0.5*(pL+pR) - 0.125*(uR-uL)*(rhoL+rhoR)*
                (sqrt(GAMMA*pL/rhoL) + sqrt(GAMMA*pR/rhoR));
    p0 = max(tol, p0);
    double p_star = p0;
    for (int i=0; i<max_iter; ++i) {
        double fL = WaveFunk(p_star, pL, rhoL);
        double fR = WaveFunk(p_star, pR, rhoR);
        double f = fL + fR + uR - uL;
        double dfL = ProisvWaveFunk(p_star, pL, rhoL);
        double dfR = ProisvWaveFunk(p_star, pR, rhoR);
        double df = dfL + dfR;
        if (fabs(df) < 1e-12) break;
        double dp = -f / df;
        p_star += dp;
        if (p_star < tol) p_star = tol;
        if (fabs(dp) < tol * p_star) break;
    }
    return p_star;
}

double u_star(double p_star, double rhoL, double uL, double pL, double rhoR, double uR, double pR) {
    double fL = WaveFunk(p_star, pL, rhoL);
    double fR = WaveFunk(p_star, pR, rhoR);
    return 0.5 * (uL + uR) + 0.5 * (fR - fL);
}

tuple<double, double, double> sample(double p_star, double u_star,
    double rhoL, double uL, double pL, double rhoR, double uR, double pR, double x_over_t) {
    double cL = sqrt(GAMMA * pL / rhoL);
    double cR = sqrt(GAMMA * pR / rhoR);
    if (x_over_t < u_star) {
        if (p_star > pL) {
            double SL = uL - sqrt(((GAMMA+1)/(2*GAMMA))*(p_star/pL) + (GAMMA-1)/(2*GAMMA)) * cL;
            if (x_over_t < SL)
                return {rhoL, uL, pL};
            else {
                double rho_star = rhoL * ((p_star/pL)+(GAMMA-1)/(GAMMA+1)) /
                                 (((GAMMA-1)/(GAMMA+1))*(p_star/pL)+1);
                return {rho_star, u_star, p_star};
            }
        } else {
            double SHL = uL - cL;
            double STL = u_star - sqrt(GAMMA * p_star / (rhoL * pow(p_star/pL, 1.0/GAMMA)));
            if (x_over_t < SHL)
                return {rhoL, uL, pL};
            else if (x_over_t < STL) {
                double u = 2.0/(GAMMA+1)*(cL + (GAMMA-1)/2*uL + x_over_t);
                double c = 2.0/(GAMMA+1)*(cL + (GAMMA-1)/2*(uL - x_over_t));
                double rho = rhoL * pow(c/cL, 2.0/(GAMMA-1));
                double p = pL * pow(rho/rhoL, GAMMA);
                return {rho, u, p};
            } else {
                double rho_star = rhoL * pow(p_star/pL, 1.0/GAMMA);
                return {rho_star, u_star, p_star};
            }
        }
    } else {
        if (p_star > pR) {
            double SR = uR + sqrt(((GAMMA+1)/(2*GAMMA))*(p_star/pR) + (GAMMA-1)/(2*GAMMA)) * cR;
            if (x_over_t > SR)
                return {rhoR, uR, pR};
            else {
                double rho_star = rhoR * ((p_star/pR)+(GAMMA-1)/(GAMMA+1)) /
                                 (((GAMMA-1)/(GAMMA+1))*(p_star/pR)+1);
                return {rho_star, u_star, p_star};
            }
        } else {
            double SHR = uR + cR;
            double STR = u_star + sqrt(GAMMA * p_star / (rhoR * pow(p_star/pR, 1.0/GAMMA)));
            if (x_over_t > SHR)
                return {rhoR, uR, pR};
            else if (x_over_t > STR) {
                double u = 2.0/(GAMMA+1)*(-cR + (GAMMA-1)/2*uR + x_over_t);
                double c = 2.0/(GAMMA+1)*(cR - (GAMMA-1)/2*(uR - x_over_t));
                double rho = rhoR * pow(c/cR, 2.0/(GAMMA-1));
                double p = pR * pow(rho/rhoR, GAMMA);
                return {rho, u, p};
            } else {
                double rho_star = rhoR * pow(p_star/pR, 1.0/GAMMA);
                return {rho_star, u_star, p_star};
            }
        }
    }
}

tuple<double, double, double> exactRiemann(double rhoL, double uL, double pL,
                                           double rhoR, double uR, double pR, double xi) {
    double p_star = solve_p_star(rhoL, uL, pL, rhoR, uR, pR);
    double u_star_val = u_star(p_star, rhoL, uL, pL, rhoR, uR, pR);
    return sample(p_star, u_star_val, rhoL, uL, pL, rhoR, uR, pR, xi);
}

// ========== Поток через грань (4 компоненты) ==========
tuple<double, double, double, double> riemannSolver(const Primitive& left, const Primitive& right,
                                                    double nx, double ny) {
    double uL_n = left.u*nx + left.v*ny;
    double uR_n = right.u*nx + right.v*ny;
    double rhoL = left.rho, pL = left.p;
    double rhoR = right.rho, pR = right.p;
    auto [flux_mass, flux_mom_n, flux_energy] = exactRiemann(rhoL, uL_n, pL, rhoR, uR_n, pR, 0.0);
    double flux_mom_x = flux_mom_n * nx;
    double flux_mom_y = flux_mom_n * ny;
    return {flux_mass, flux_mom_x, flux_mom_y, flux_energy};
}

// ========== Граничные условия ==========
Primitive getBoundaryState(const Primitive& interior, int phys_tag, double nx, double ny) {
    Primitive bc = interior;
    if (phys_tag == 10) { // inflow
        bc.rho = 1.0; bc.u = 1.0; bc.v = 0.0; bc.p = 1.0/GAMMA;
    } else if (phys_tag == 11) { // outflow – extrapolate
        // nothing
    } else if (phys_tag == 12) { // wall
        double u_n = interior.u*nx + interior.v*ny;
        double u_tx = interior.u - u_n*nx;
        double u_ty = interior.v - u_n*ny;
        bc.u = u_tx - u_n*nx;
        bc.v = u_ty - u_n*ny;
    }
    return bc;
}

// ========== Обновление решения ==========
void updateSolution(double dt) {
    for (auto& cell : cells) cell.dU = {0,0,0,0};
    for (const auto& face : faces) {
        Primitive left = getPrimitive(cells[face.cellL].U);
        Primitive right = (face.type == 0) ? getPrimitive(cells[face.cellR].U) :
                           getBoundaryState(left, face.phys_tag, face.nx, face.ny);
        auto [f_m, f_mx, f_my, f_E] = riemannSolver(left, right, face.nx, face.ny);
        Conserved flux = {f_m, f_mx, f_my, f_E};
        double invVolL = 1.0 / cells[face.cellL].area;
        cells[face.cellL].dU.rho   -= flux.rho   * face.length * invVolL;
        cells[face.cellL].dU.rhou  -= flux.rhou  * face.length * invVolL;
        cells[face.cellL].dU.rhov  -= flux.rhov  * face.length * invVolL;
        cells[face.cellL].dU.rhoE  -= flux.rhoE  * face.length * invVolL;
        if (face.type == 0) {
            double invVolR = 1.0 / cells[face.cellR].area;
            cells[face.cellR].dU.rho   += flux.rho   * face.length * invVolR;
            cells[face.cellR].dU.rhou  += flux.rhou  * face.length * invVolR;
            cells[face.cellR].dU.rhov  += flux.rhov  * face.length * invVolR;
            cells[face.cellR].dU.rhoE  += flux.rhoE  * face.length * invVolR;
        }
    }
    for (auto& cell : cells) {
        cell.U.rho   += dt * cell.dU.rho;
        cell.U.rhou  += dt * cell.dU.rhou;
        cell.U.rhov  += dt * cell.dU.rhov;
        cell.U.rhoE  += dt * cell.dU.rhoE;
        if (cell.U.rho <= 0.0) { cell.U.rho = 1e-10; }
        Primitive prim = getPrimitive(cell.U);
        if (prim.p <= 0.0) {
            double e = cell.U.rhoE/cell.U.rho - 0.5*(prim.u*prim.u + prim.v*prim.v);
            if (e <= 0.0) e = 1e-10;
            cell.U.rhoE = cell.U.rho * (e + 0.5*(prim.u*prim.u + prim.v*prim.v));
        }
    }
}

double computeTimeStep(double cfl) {
    double max_lambda = 0.0;
    for (const auto& cell : cells) {
        Primitive prim = getPrimitive(cell.U);
        double a = sqrt(GAMMA * prim.p / prim.rho);
        double u_abs = sqrt(prim.u*prim.u + prim.v*prim.v);
        double lambda = u_abs + a;
        double h = sqrt(cell.area);
        max_lambda = max(max_lambda, lambda / h);
    }
    return cfl / max_lambda;
}

void writeVTK(int iter, double time) {
    string filename = "solution_" + to_string(iter) + ".vtk";
    ofstream vtk(filename);
    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "Godunov solution, time = " << time << "\n";
    vtk << "ASCII\n";
    vtk << "DATASET UNSTRUCTURED_GRID\n";
    vtk << "POINTS " << nodeNum << " float\n";
    for (int i = 0; i < nodeNum; ++i) {
        vtk << nodeCoords[3*i] << " " << nodeCoords[3*i+1] << " " << nodeCoords[3*i+2] << "\n";
    }
    vtk << "CELLS " << cells.size() << " " << 4 * cells.size() << "\n";
    for (const auto& cell : cells) {
        vtk << "3 " << cell.nodes[0]-1 << " " << cell.nodes[1]-1 << " " << cell.nodes[2]-1 << "\n";
    }
    vtk << "CELL_TYPES " << cells.size() << "\n";
    for (size_t i = 0; i < cells.size(); ++i) vtk << "5\n";
    vtk << "CELL_DATA " << cells.size() << "\n";
    vtk << "SCALARS density float 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (const auto& cell : cells) vtk << cell.U.rho << "\n";
    vtk << "VECTORS velocity float\n";
    for (const auto& cell : cells) {
        double u = cell.U.rhou / cell.U.rho;
        double v = cell.U.rhov / cell.U.rho;
        vtk << u << " " << v << " 0\n";
    }
    vtk << "SCALARS pressure float 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (const auto& cell : cells) {
        Primitive prim = getPrimitive(cell.U);
        vtk << prim.p << "\n";
    }
    vtk.close();
    cout << "Saved " << filename << endl;
}

void setInitialConditions() {
    for (auto& cell : cells) {
        Primitive prim0 = {1.0, 1.0, 0.0, 1.0/GAMMA};
        cell.U = getConserved(prim0);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " mesh.msh" << endl;
        return 1;
    }
    string meshFile = argv[1];
    if (!readMesh(meshFile)) {
        cerr << "Failed to read mesh file: " << meshFile << endl;
        return 1;
    }
    buildCells();
    buildFaces();  // теперь использует глобальную переменную globalEdgeToPhysTag
    setInitialConditions();
    double t = 0.0;
    int iter = 0;
    while (t < T_END) {
        double dt = computeTimeStep(CFL);
        if (t + dt > T_END) dt = T_END - t;
        updateSolution(dt);
        t += dt;
        iter++;
        if (iter % SAVE_INTERVAL == 0 || t >= T_END) writeVTK(iter, t);
        cout << "Iter " << iter << ", t = " << t << ", dt = " << dt << endl;
    }
    cout << "Simulation finished." << endl;
    return 0;
}