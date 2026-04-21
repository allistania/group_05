#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <iostream>
#include <cmath>
// «апись таблицы данных дл€ фиксированного момента времени
void outputCSV(std::string name, 
               const std::vector<double>& grid_P,
               const std::vector<double>& u,
               const std::vector<double>& P,
               const std::vector<double>& rho,
               const std::vector<double>& I,
               double t);

// ‘ункци€ дл€ создани€ записи с несколькими моментами времени
std::function<bool(double, const std::vector<double>&, const std::vector<double>&,
                   const std::vector<double>&, const std::vector<double>&,
                   const std::vector<double>&)>
createTimeSeriesWriter(double output_interval, std::string base_name);
void clearDirectory(const std::string& directoryPath);

std::tuple<double, double, double> hll(double rhoL, double uL, double pL,
                                       double rhoR, double uR, double pR,
                                       double xi, double gamma = 1.4);

std::tuple<double, double, double> hllc(double rhoL, double uL, double pL,
                                        double rhoR, double uR, double pR,
                                        double xi, double gamma = 1.4);
std::tuple<double, double, double> rusanov(double rhoL, double uL, double pL,
                                          double rhoR, double uR, double pR,
                                          double xi);
#endif