#include "parameters.h"
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cmath>

// функция для создания теста
void writeTest() {
    std::string test;
    std::cout << "Название теста ";
    std::cin >> test;
    std::ofstream file(test);
    std::string rho_L, u_L, p_L, rho_R, u_R, p_R, t_end;
    std::cout << "rho_L= ";
    std::cin >> rho_L;
    std::cout << "u_L= ";
    std::cin >> u_L;
    std::cout << "p_L= ";
    std::cin >> p_L;
    std::cout << "rho_R= ";
    std::cin >> rho_R;
    std::cout << "u_R= ";
    std::cin >> u_R;
    std::cout << "p_R= ";
    std::cin >> p_R;
    std::cout << "t_end= ";
    std::cin >> t_end;

    if (file.is_open()) {
        file << "rho_L= " << rho_L << "\n"; 
        file << "u_L= " << u_L << "\n";
        file << "p_L= " << p_L << "\n";
        file << "rho_R= " << rho_R << "\n"; 
        file << "u_R= " << u_R << "\n";
        file << "p_R= " << p_R << "\n";
        file << "t_end= " << t_end << "\n";
        file << "\n";
        file << "Руководство к тесту:\n";
        file << "В тестах содержится 7 параметров для тестов\n";
        file << "1: rho_L и rho_R - это начальные плотности газа слева и справа на бесконечности\n";
        file << "2: p_L и p_R - это начальные давления газа слева и справа на бесконечности\n";
        file << "3: u_L и u_R - это начальные скорости газа слева и справа на бесконечности\n"; 
        file << "4: t_end - время расчёта";

        file.close();
        std::cout << "Данные записаны: " << test << std::endl;
    } else {
        std::cerr << "Ошибка не удалось записать в файл " << test << std::endl;
    }
}

//инициализация параметров 
void start() {
    std::string start;
    std::cout << "Name of init file: ";
    std::cin >> start;
    std::ofstream file(start);
    
    std::string tau, h, cells, cfl,solver_name, left_bc, right_bc, test_name, count, t_output, key;
    std::vector<std::string> available_solvers = {
        "godunov", "godunov-kolgan", "godunov-kolgan-rodionov", "eno", "weno" , "Fletcher"
    };
    
    std::cout << "Time step: ";
    std::cin >> tau;
    std::cout << "Coordinate step: ";
    std::cin >> h;
    std::cout << "Nubmer of cells: ";
    std::cin >> cells;
    std::cout << "CFK: ";
    std::cin >> cfl;
    std::cout << "Left boundary condition (wall/flow/periodic): ";
    std::cin >> left_bc;
    std::cout << "Right boundary condition (wall/flow/periodic): ";
    std::cin >> right_bc;
    std::cout << "Sod test (sod1/sod2/sod3/sod4/sod5/custom): ";
    std::cin >> test_name;
    std::string t_end_user;
    std::cout << "End time (enter 'default' or specific value): ";
    std::cin >> t_end_user; 
    std::cout << "count: ";
    std::cin >> count;
    std::cout << "Output interval for gif (time between frames): ";
    std::cin >> t_output;
    std::cout << "true - if time steps and talse if count steps: ";
    std::cin >> key;
    std::cout << "\nAvailable solvers: ";
    for (size_t i = 0; i< available_solvers.size(); i++) {
        std::cout << available_solvers[i];
        if (i< available_solvers.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    int num_solvers;
    std::cout << "Numbers of solvers";
    std::cin >> num_solvers;
    
    std::vector<std::string> selected_solvers;
    for (int i = 0; i< num_solvers; i++){
        std::string solver;
        std::cout << "Enter solver number " << i+1<< ": ";
        std::cin >> solver; 
        selected_solvers.push_back(solver);
    }

    std::string solver_names_str; 
    for (size_t i = 0; i<selected_solvers.size(); i++){
        solver_names_str +=selected_solvers[i];
        if (i< selected_solvers.size() -1) solver_names_str+= ",";
    }

    if (file.is_open()) {
        file << "tau = " << tau << "\n"; 
        file << "h = " << h << "\n"; 
        file << "cells = " << cells << "\n";
        file << "cfl = " << cfl << "\n";
        file << "left_bc = " << left_bc << "\n";
        file << "right_bc = " << right_bc << "\n";
        file << "test_name = " << test_name << "\n";
        file << "t_end = " << t_end_user << "\n";
        file << "count = " << count << "\n";
        file << "t_output = " << t_output << "\n";
        file << "key = " << key << "\n";
        file << "solver_names = " << solver_names_str << "\n";
        file << "\n";
        file << "Руководство к файлу инициализации:\n";
        file << "В файле содержится 7 параметров для инициализации\n";
        file << "1: tau и h - это начальные шаги по времени соответственно\n";
        file << "2: cells - это количество расчётных ячеек\n";
        file << "3: CFL - число Куранта — Фридрихса — Леви\n";
        file << "4: left_dc и right_bc - это, соответственно, левое и правое граничное условие. Реализованы трех видов\n"; 
        file << "wall - отражение от стенки \n";
        file << "flow - поток \n";
        file << "periodic - периодическое ГУ \n";
        file << "5: test_name - имя теста. Если нужен готовый тест Сода, то пишем sod{номер теста}. Например, sod1\n";
        file << "6: разностная схема: Годунова, Годунова-Колгана, Годунова-Колгана-Родионова, ЕНО или ВЕнО \n";
        file << "7: параметр count отвечает за запись файла для гифки. Например, если count=100, то каждую сотую операцию будет запись файла \n";
        file << "8: в отличие от параметра count, t_output отвечает за запись файла для гифки не через count-ый шаг, а через t_output-ое время\n";
        file << "9: True если по времени запись и False если по шагам:\n";
        file.close();
        std::cout << "Данные записаны: " << start << std::endl;
    } else {
        std::cerr << "Ошибка не удалось записать в файл " << start << std::endl;
    }
}

bool readTest(const std::string& filename, double& rho_L, double& u_L, double& p_L,
              double& rho_R, double& u_R, double& p_R, double& t_end) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: failed opening of file " << filename << std::endl;
        return false;
    }

    std::string line;
    std::map<std::string, double*> paramMap = {
        {"rho_L", &rho_L},
        {"u_L", &u_L},
        {"p_L", &p_L},
        {"rho_R", &rho_R},
        {"u_R", &u_R},
        {"p_R", &p_R},
        {"t_end", &t_end}
    };
    int paramsFound = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        size_t equalPos = line.find('=');
        if (equalPos == std::string::npos) continue;
        std::string paramName = line.substr(0, equalPos);
        size_t start = paramName.find_first_not_of(" \t");
        size_t end = paramName.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            paramName = paramName.substr(start, end - start + 1);
        }

        std::string valueStr = line.substr(equalPos + 1);
        try {
            double value = std::stod(valueStr);
            if (paramMap.find(paramName) != paramMap.end()) {
                *paramMap[paramName] = value;
                paramsFound++;
                std::cout << "Read parameter: " << paramName << " = " << value << std::endl;
            } else {
                std::cout << "Предупреждение: неизвестный параметр '" << paramName << "'" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Ошибка преобразования значения для параметра '" << paramName << "'" << std::endl;
        }
    }    
    file.close();
    if (paramsFound == 7) {
        return true;
    } else {
        std::cerr << "Предупреждение: прочитано только " << paramsFound << " из 7 параметров" << std::endl;
        return paramsFound > 0;
    }
}

bool readInit(const std::string& filename, double& tau, double& h, 
              int& cells, double& cfl, std::string& left_bc, 
              std::string& right_bc, std::string& test_name, int& count, std::vector<std::string>& solver_names,
              std::string& t_end_user, double& t_output, std::string& key) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: failed opening of file " << filename << std::endl;
        return false;
    }

    std::string line;
    std::string solver_names_str;
    std::map<std::string, void*> paramMap = {
        {"tau", &tau},
        {"h", &h},
        {"cells", &cells},
        {"cfl", &cfl},
        {"left_bc", &left_bc},
        {"right_bc", &right_bc},
        {"test_name", &test_name},
        {"count", &count},
        {"solver_names", &solver_names_str},
        {"t_end", &t_end_user}, 
        {"t_output", &t_output}, 
        {"key", &key}
    };
    int paramsFound = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line.empty() || line.find("Руководство") != std::string::npos) {
        break;  // прекращаем чтение при пустой строке или начале руководства
        }
        size_t equalPos = line.find('=');
        if (equalPos == std::string::npos) continue;
        std::string paramName = line.substr(0, equalPos);
        size_t start = paramName.find_first_not_of(" \t");
        size_t end = paramName.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            paramName = paramName.substr(start, end - start + 1);
        }

        std::string valueStr = line.substr(equalPos + 1);
        // Убираем пробелы вокруг значения
        start = valueStr.find_first_not_of(" \t");
        end = valueStr.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            valueStr = valueStr.substr(start, end - start + 1);
        }

        try {
            if (paramMap.find(paramName) != paramMap.end()) {
                if (paramName == "tau" || paramName == "h" || paramName == "cfl" || paramName == "t_output") {
                    *static_cast<double*>(paramMap[paramName]) = std::stod(valueStr);
                } else if (paramName == "cells" || paramName == "count") {
                    *static_cast<int*>(paramMap[paramName]) = std::stoi(valueStr);
                } else {
                    *static_cast<std::string*>(paramMap[paramName]) = valueStr;
                }
                paramsFound++;
                std::cout << "Read parameter: " << paramName << " = " << valueStr << std::endl;
            } else {
                std::cout << "Warning: unknown parameter '" << paramName << "'" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Ошибка преобразования значения для параметра '" << paramName << "'" << std::endl;
        }
    }    
    file.close();


        if (!solver_names_str.empty()) {
        size_t start = 0;
        size_t end = solver_names_str.find(',');
        while (end != std::string::npos) {
            std::string solver = solver_names_str.substr(start, end - start);
            // Убираем пробелы вокруг имени схемы
            size_t solver_start = solver.find_first_not_of(" \t");
            size_t solver_end = solver.find_last_not_of(" \t");
            if (solver_start != std::string::npos && solver_end != std::string::npos) {
                solver = solver.substr(solver_start, solver_end - solver_start + 1);
            }
            solver_names.push_back(solver);
            start = end + 1;
            end = solver_names_str.find(',', start);
        }
        // Добавляем последнюю схему
        std::string last_solver = solver_names_str.substr(start);
        size_t solver_start = last_solver.find_first_not_of(" \t");
        size_t solver_end = last_solver.find_last_not_of(" \t");
        if (solver_start != std::string::npos && solver_end != std::string::npos) {
            last_solver = last_solver.substr(solver_start, solver_end - solver_start + 1);
        }
        solver_names.push_back(last_solver);
        
        std::cout << "Found " << solver_names.size() << " solvers: ";
        for (size_t i = 0; i < solver_names.size(); i++) {
            std::cout << solver_names[i];
            if (i < solver_names.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    
    if (paramsFound >= 7) {
        return true;
    } else {
        std::cerr << "Warning! Read only " << paramsFound << " 7 params" << std::endl;
        return paramsFound > 0;
    }
}

//для получения фикт ячеек

int getGhostCellsBySolver (const std::string& solver_name){
    if (solver_name == "godunov-kolgan" || 
        solver_name == "godunov-kolgan-rodionov" || 
        solver_name == "eno" ||
        solver_name == "weno"){return 2;}
    if (solver_name == "godunov"){return 1;}

    return 1;
}

// Функция для установки готовых тестов Сода
void setupSODTest(const std::string& test_name, 
                  double& rho_L, double& u_L, double& p_L,
                  double& rho_R, double& u_R, double& p_R, 
                  double& t_end,
                  const std::string& t_end_user) {
    
    if (test_name == "sod1") {
        rho_L = 1.0; u_L = 0.0; p_L = 1.0;
        rho_R = 0.125; u_R = 0.0; p_R = 0.1;
        t_end = 0.25;
    }
    else if (test_name == "sod2") {
        rho_L = 1.0; u_L = -2.0; p_L = 0.4;
        rho_R = 1.0; u_R = 2.0; p_R = 0.4;
        t_end = 0.15;
    }
    else if (test_name == "sod3") {
        rho_L = 1.0; u_L = 0.0; p_L = 1000.0;
        rho_R = 1.0; u_R = 0.0; p_R = 0.01;
        t_end = 0.012;
    }
    else if (test_name == "sod4") {
        rho_L = 5.99924; u_L = 19.5975; p_L = 460.894;
        rho_R = 5.99242; u_R = -6.19633; p_R = 46.0950;
        t_end = 0.035;
    }
    else if (test_name == "sod5") {
        rho_L = 1.0; u_L = 0.0; p_L = 0.01;
        rho_R = 1.0; u_R = 0.0; p_R = 100.0;
        t_end = 0.035;
    }
    else {
        // По умолчанию sod1
        rho_L = 1.0; u_L = 0.0; p_L = 1.0;
        rho_R = 0.125; u_R = 0.0; p_R = 0.1;
        t_end = 0.2;
    }
    
    // ПРОВЕРЯЕМ ПОЛЬЗОВАТЕЛЬСКОЕ ВРЕМЯ
    if (t_end_user != "default" && !t_end_user.empty()) {
        try {
            double user_time = std::stod(t_end_user);
            t_end = user_time;  // используем пользовательское время
            std::cout << "Using user-specified time: " << t_end << std::endl;
        } catch (const std::exception& e) {
            std::cout << "Invalid time format, using default: " << t_end << std::endl;
        }
    } else {
        std::cout << "Using default time: " << t_end << std::endl;
    }
    
    std::cout << "Test: " << test_name << std::endl;
    std::cout << "rho_L=" << rho_L << " u_L=" << u_L << " p_L=" << p_L << std::endl;
    std::cout << "rho_R=" << rho_R << " u_R=" << u_R << " p_R=" << p_R << std::endl;
    std::cout << "t_end=" << t_end << std::endl;
}

// Функция для применения граничных условий
void applyBoundaryConditions(double* rho, double* u, double* p,
                            int cells, int ghost_cells,
                            const std::string& left_bc,
                            const std::string& right_bc) {
    
    // Индексация:
    // Внутренние ячейки: [ghost_cells, ghost_cells + cells - 1]
    // Левые фиктивные:   [0, ghost_cells - 1]
    // Правые фиктивные:  [ghost_cells + cells, ghost_cells + cells + ghost_cells - 1]

    // --- Левая граница ---
    if (left_bc == "wall") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[i] = rho[ghost_cells];
            u[i]   = -u[ghost_cells];
            p[i]   = p[ghost_cells];
        }
    }
    else if (left_bc == "flow") {   // экстраполяция нулевого порядка
        for (int i = 0; i < ghost_cells; ++i) {
            rho[i] = rho[ghost_cells];
            u[i]   = u[ghost_cells];
            p[i]   = p[ghost_cells];
        }
    }
    else if (left_bc == "periodic") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[i] = rho[cells + ghost_cells - 1 - i];
            u[i]   = u[cells + ghost_cells - 1 - i];
            p[i]   = p[cells + ghost_cells - 1 - i];
        }
    }
    else if (left_bc == "inlet") {
        // Фиксированные параметры входа (можно вынести в аргументы или константы)
        const double rho_in = 1.0;
        const double u_in   = 1;   
        const double p_in   = 1.0;
        for (int i = 0; i <= ghost_cells; ++i) {
            rho[i] = rho_in;
            u[i]   = u_in;
            p[i]   = p_in;
        }
    }
    else if (left_bc == "outlet") {
        // Выход: экстраполяция нулевого порядка из первой внутренней ячейки
        for (int i = 0; i < ghost_cells; ++i) {
            rho[i] = rho[ghost_cells];
            u[i]   = u[ghost_cells];
            p[i]   = p[ghost_cells];
        }
    }

    // --- Правая граница ---
    if (right_bc == "wall") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho[cells + ghost_cells - 1];
            u[cells + ghost_cells + i]   = -u[cells + ghost_cells - 1];
            p[cells + ghost_cells + i]   = p[cells + ghost_cells - 1];
        }
    }
    else if (right_bc == "flow") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho[cells + ghost_cells - 1];
            u[cells + ghost_cells + i]   = u[cells + ghost_cells - 1];
            p[cells + ghost_cells + i]   = p[cells + ghost_cells - 1];
        }
    }
    else if (right_bc == "periodic") {
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho[ghost_cells + i];
            u[cells + ghost_cells + i]   = u[ghost_cells + i];
            p[cells + ghost_cells + i]   = p[ghost_cells + i];
        }
    }
    else if (right_bc == "inlet") {
        // Вход справа: скорость отрицательная (поток влево)
        const double rho_in = 1.0;
        const double u_in   = -0.5;   // отрицательная скорость
        const double p_in   = 1.0;
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho_in;
            u[cells + ghost_cells + i]   = u_in;
            p[cells + ghost_cells + i]   = p_in;
        }
    }
else if (right_bc == "outlet") {
    const double GAMMA = 1.4;
    const double P_OUT = 1.0;   // заданное давление на выходе

    int last = ghost_cells + cells - 1;      // последняя внутренняя ячейка
    int last2 = last - 1;                     // предпоследняя (для экстраполяции энтропии)

    double rhoL = rho[last];
    double uL   = u[last];
    double pL   = p[last];
    double cL   = std::sqrt(GAMMA * pL / rhoL);

    // Экстраполируем энтропию (s = p / rho^gamma) изнутри
    double sL = pL / pow(rhoL, GAMMA);
    // Можно взять среднее или экстраполяцию нулевого порядка
    double s_boundary = sL;   // либо использовать last2 для линейной экстраполяции

    // Проверяем тип течения на границе (по числу Маха)
    double M = uL / cL;

    if (M >= 1.0) {
        // Сверхзвуковой отток — все величины экстраполируются
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rhoL;
            u[cells + ghost_cells + i]   = uL;
            p[cells + ghost_cells + i]   = pL;
        }
    }
    else if (M > 0.0) {
        // Дозвуковой отток: одна характеристика входит (давление задано),
        // остальные выходят (экстраполируются)
        // Вычисляем плотность и скорость на границе из заданного давления и экстраполированной энтропии
        double rho_b = pow(P_OUT / s_boundary, 1.0 / GAMMA);
        double c_b   = std::sqrt(GAMMA * P_OUT / rho_b);
        // Инвариант Римана, приходящий изнутри: J_minus = u - 2c/(gamma-1)
        double J_minus = uL - 2.0 * cL / (GAMMA - 1.0);
        // Используем его для определения скорости на границе
        double u_b = J_minus + 2.0 * c_b / (GAMMA - 1.0);

        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rho_b;
            u[cells + ghost_cells + i]   = u_b;
            p[cells + ghost_cells + i]   = P_OUT;
        }
    }
    else {
        // Поток направлен внутрь — это уже не выход, а вход. Обрабатываем как вход (если такое возможно)
        // Например, экстраполяция (аварийная ситуация)
        for (int i = 0; i < ghost_cells; ++i) {
            rho[cells + ghost_cells + i] = rhoL;
            u[cells + ghost_cells + i]   = uL;
            p[cells + ghost_cells + i]   = pL;
        }
    }
}}