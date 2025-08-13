#include <iostream>

void initial_memory(int mem_size);

int get_regions(double* y_coord1, double* z_coord1, double* y_coord2, double* z_coord2, int N1, int N2, double* surface_norm, bool &error_flag);

double ecalcgen(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
                 double ym, double yp, double sig, double cene, double n, double c3now, int caseind);

double ecalcc2(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
               double ym, double yp, double sig, double cene, double n, double c3now, int caseind);

double ecalcc1(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
               double ym, double yp, double sig, double cene, double n, double c3now, int caseind);

double calc_gen(double sig, double lc, double atre1, double atre2, double repe, double catr1, double catr2,double crep, const std::vector<std::vector<double>>& region, double Epsilon, double* surface_norm);

double surface_energy(double sig, double lc, double atre1, double atre2, double repe, double catr1, double catr2,double crep, double Epsilon,
                      double* y_coord1, double* z_coord1, double* y_coord2, double* z_coord2, int N1, int N2, double* surface_norm);
