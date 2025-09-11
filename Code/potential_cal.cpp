#include "PolyUtils.h"
#include "PolygonClipping.h"
#include <iostream>
#include <cmath>
#include <vector>

static std::vector<std::vector<double>> region_y;
static std::vector<std::vector<double>> region_z;
static std::vector<std::vector<double>> region_x;
static std::vector<double> tempy;
static std::vector<double> tempz;
static std::vector<std::vector<double>> region_coord;

double mini_v(double a, double b) {
    return (b < a) ? b : a;
}

double maxi_v(double a, double b) {
    return (b > a) ? b : a;
}

void initial_memory(int mem_size) {
    region_x.resize(mem_size, std::vector<double>(4));
    region_y.resize(mem_size, std::vector<double>(4));
    region_z.resize(mem_size, std::vector<double>(4));
    tempy.resize(mem_size);
    tempz.resize(mem_size);
    region_coord.resize(4, std::vector<double>(2));
}

int get_regions(double* y_coord1, double* z_coord1, double* y_coord2, double* z_coord2, int N1, int N2, double* surface_norm, bool &error_flag) {
    int i;
    int j;
    int k;
    int m;
    int n;
    int flag = 1;
    double temp;
    double y1;
    double z1;
    double y2;
    double z2;
    int Num_regions = 0;
    int count_num = 0;

    error_flag = false; // Initial state
    //std::cout << "Initial error_flag in get_regions: " << error_flag << std::endl;

    std::vector<PolyClip::Point2d> vertices1;
    for (i = 0; i < N1; i++) {
        vertices1.push_back(PolyClip::Point2d(*(y_coord1 + i), *(z_coord1 + i)));
    }
    PolyClip::Polygon polygon1(vertices1);

    std::vector<PolyClip::Point2d> vertices2;
    for (i = 0; i < N2; i++) {
        vertices2.push_back(PolyClip::Point2d(*(y_coord2 + i), *(z_coord2 + i)));
    }
    PolyClip::Polygon polygon2(vertices2);

    PolyClip::PloygonOpration::DetectIntersection(polygon1, polygon2, error_flag);
    //std::cout << "error_flag_get_regions " << error_flag << std::endl;
    if (error_flag) {
        return 0;
    }

    std::vector<std::vector<PolyClip::Point2d>> possible_result;
    if (PolyClip::PloygonOpration::Mark(polygon1, polygon2, possible_result, PolyClip::MarkIntersection)) {
        std::vector<std::vector<PolyClip::Point2d>> results = PolyClip::PloygonOpration::ExtractIntersectionResults(polygon1);
        for (i = 0; i < results.size(); ++i) {
            for (auto p : results[i]) {
                tempy[count_num] = p.x_;
                tempz[count_num] = p.y_;

                if (count_num > 0 && std::abs(tempy[count_num] - tempy[count_num - 1]) < 1e-5 && std::abs(tempz[count_num] - tempz[count_num - 1]) < 1e-5) {
                    continue;
                }
                //printf("%lf\n",tempy[count_num]);
                count_num++;
            }
            //std::cout << "intersection area is\n";
            //std::cout << PolyClip::Utils::CalculatePolygonArea(results[0]) << "\n";
        }
    } else {
        return 0;
    }
    count_num--;

    if (std::abs(tempy[count_num - 1] - tempy[0]) < 1e-5 && std::abs(tempz[count_num - 1] - tempz[0]) < 1e-5) {
        count_num--;
    }
    /*printf("Input coordinates (polygon1)(y,z)\n");
    for(i=0;i<N1;i++){
        printf("%lf     %lf\n",y_coord1[i],z_coord1[i]);
        printf("\n");
    }
    printf("Input coordinates (polygon2)(y,z)\n");
    for(i=0;i<N2;i++){
        printf("%lf     %lf\n",y_coord2[i],z_coord2[i]);
        printf("\n");
    }
    printf("Output coordinates(y,z)\n");
    for(i=0;i<count_num;i++){
        printf("%lf     %lf\n",tempy[i],tempz[i]);
        printf("\n");
    }

	printf("%d\n",count_num);*/

    for (i = 0; i < count_num; i++) {
        if (count_num <= 2) { break; }
        for (j = 0; j < count_num; j++) {
            flag = 1;
            m = (j - 1 + count_num) % count_num;
            n = (j + 1) % count_num;
            y1 = tempy[j] - tempy[m];
            y2 = tempy[n] - tempy[j];
            z1 = tempz[j] - tempz[m];
            z2 = tempz[n] - tempz[j];
            //std::cout << "tempy[m] tempz[m] " << tempy[m] << " " << tempz[m] << std::endl;
            for (k = 0; k < count_num; k++) {
                if (count_num <= 3) { flag = -1; break; }
                if (k == m || k == n || k == j) { continue; }
                if ((tempy[j] - tempy[m]) * (tempz[k] - tempz[m]) - (tempy[k] - tempy[m]) * (tempz[j] - tempz[m]) < 0.0) { flag = -1; break; }
                if ((tempy[n] - tempy[j]) * (tempz[k] - tempz[j]) - (tempy[k] - tempy[j]) * (tempz[n] - tempz[j]) < 0.0) { flag = -1; break; }
                if ((tempy[m] - tempy[n]) * (tempz[k] - tempz[n]) - (tempy[k] - tempy[n]) * (tempz[m] - tempz[n]) < 0.0) { flag = -1; break; }
            }

            if (flag > 0) { continue; }
            region_y[Num_regions][0] = tempy[m];
            region_y[Num_regions][1] = tempy[j];
            region_y[Num_regions][2] = tempy[n];

            region_z[Num_regions][0] = tempz[m];
            region_z[Num_regions][1] = tempz[j];
            region_z[Num_regions][2] = tempz[n];
            //std::cout << "Num_regions " << Num_regions << std::endl;
            //std::cout << "region_x[Num_regions] " << region_x[Num_regions][0] << " " << region_x[Num_regions][1] << " " << region_x[Num_regions][2] << std::endl;
            region_x[Num_regions][0] = -(*(surface_norm + 1) * (tempy[m]) + *(surface_norm + 2) * (tempz[m]) + *(surface_norm + 3)) / (*surface_norm);
            //std::cout << "region_x[Num_regions][0] " << region_x[Num_regions][0] << std::endl;
            bool check1 = (region_x[Num_regions][0] - 1.0) < -(1e-5); // region_x[Num_regions][0] < 1.0, -1e-5 is selected to avoid numeric errors
            if (check1) { return -1; }
            region_x[Num_regions][1] = -(*(surface_norm + 1) * (tempy[j]) + *(surface_norm + 2) * (tempz[j]) + *(surface_norm + 3)) / (*surface_norm);
            //std::cout << "region_x[Num_regions][1] " << region_x[Num_regions][1] << std::endl;
            bool check2 = (region_x[Num_regions][1] - 1.0) < -(1e-5); // region_x[Num_regions][1] < 1.0, -1e-5 is selected to avoid numeric errors
            if (check2) { return -1; }
            region_x[Num_regions][2] = -(*(surface_norm + 1) * (tempy[n]) + *(surface_norm + 2) * (tempz[n]) + *(surface_norm + 3)) / (*surface_norm);
            //std::cout << "region_x[Num_regions][2] " << region_x[Num_regions][2] << std::endl;
            bool check3 = (region_x[Num_regions][2] - 1.0) < -(1e-5); // region_x[Num_regions][2] < 1.0, -1e-5 is selected to avoid numeric errors
            if (check3) { return -1; }
            Num_regions++;

            for (k = j; k < count_num - 1; k++) {
                tempy[k] = tempy[k + 1];
                tempz[k] = tempz[k + 1];
            }
            count_num--;
            i--;
            break;
        }
    }
/*
    for(i=0;i<Num_regions;i++){
        printf("%lf     %lf      %lf      %lf\n",*(*(region_x+i)),*(*(region_x+i)+1),*(*(region_x+i)+2),*(*(region_x+i)+3));
        printf("%lf     %lf      %lf      %lf\n",*(*(region_y+i)),*(*(region_y+i)+1),*(*(region_y+i)+2),*(*(region_y+i)+3));
        printf("%lf     %lf      %lf      %lf\n",*(*(region_z+i)),*(*(region_z+i)+1),*(*(region_z+i)+2),*(*(region_z+i)+3));
        printf("\n");
    }
*/

    for (i = 0; i < Num_regions; i++) {
        for (j = 0; j < 3; j++) {
            if (region_y[i][0] <= region_y[i][1] && region_y[i][1] <= region_y[i][2]) { break; }
            if (region_y[i][0] >= region_y[i][1] && region_y[i][1] >= region_y[i][2]) { break; }
            double temp_y = region_y[i][2];
            region_y[i][2] = region_y[i][1];
            region_y[i][1] = region_y[i][0];
            region_y[i][0] = temp_y;

            double temp_x = region_x[i][2];
            region_x[i][2] = region_x[i][1];
            region_x[i][1] = region_x[i][0];
            region_x[i][0] = temp_x;

            double temp_z = region_z[i][2];
            region_z[i][2] = region_z[i][1];
            region_z[i][1] = region_z[i][0];
            region_z[i][0] = temp_z;
        }
        region_y[i][3] = region_y[i][2];
        region_z[i][3] = region_z[i][2];
        region_x[i][3] = region_x[i][2];

        y1 = (region_z[i][3] - region_z[i][0]) / (region_y[i][3] - region_y[i][0]);
        if (std::abs((region_y[i][3] - region_y[i][0])) < 1e-6) { continue; }
        z1 = y1 * (-region_y[i][0]) + region_z[i][0];

        z2 = y1 * (region_y[i][1]) + z1;
        region_y[i][2] = region_y[i][1];
        region_z[i][2] = z2;
        region_x[i][2] = -(*(surface_norm + 1) * (region_y[i][2]) + *(surface_norm + 2) * (z2) + *(surface_norm + 3)) / (*surface_norm);
    }

    return Num_regions;
}

double ecalcgen(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
                 double ym, double yp, double sig, double cene, double n, double c3now, int caseind) {

    double u1, u2, energy;
    sig = 1.0; // all distances are written in terms of sigma

    if (caseind == 1) {
        u1 = (yp - ym) * pow((a1 * c2 + c3now), 1.0 - n);
        u2 = (ym - yp) * pow((a2 * c2 + c3now), 1.0 - n);
        energy = cene * (u1 + u2) * pow(sig, n - 2.0) / (c2 * (n - 1.0) * (n - 2.0));
    } else if (caseind == 2) {
        u1 = (yp - ym) * pow((a1 * c2 + c3now), 1.0 - n);
        u2 = pow((a2 * c2 + c3now + (c1 + b2 * c2) * yp), 2.0 - n) - pow((a2 * c2 + c3now + (c1 + b2 * c2) * ym), 2.0 - n);
        u2 /= (c1 + b2 * c2);
        energy = cene * (u1 + u2) * pow(sig, n - 2.0) / (c2 * (n - 1.0) * (n - 2.0));
    } else if (caseind == 3) {
        u1 = pow((a1 * c2 + c3now + (c1 + b1 * c2) * ym), 2.0 - n) -
             pow((a1 * c2 + c3now + (c1 + b1 * c2) * yp), 2.0 - n);
        u1 /= (c1 + b1 * c2);
        u2 = (ym - yp) * pow((a2 * c2 + c3now), 1.0 - n);
        energy = cene * (u1 + u2) * pow(sig, n - 2.0) / (c2 * (n - 1.0) * (n - 2.0));
    } else if (caseind == 4) {
        u1 = pow((a1 * c2 + c3now + (c1 + b1 * c2) * ym), 2.0 - n) -
             pow((a1 * c2 + c3now + (c1 + b1 * c2) * yp), 2.0 - n);
        u1 /= (c1 + b1 * c2);
        u2 = pow((a2 * c2 + c3now + (c1 + b2 * c2) * yp), 2.0 - n) -
             pow((a2 * c2 + c3now + (c1 + b2 * c2) * ym), 2.0 - n);
        u2 /= (c1 + b2 * c2);
        energy = cene * (u1 + u2) * pow(sig, n - 2.0) / (c2 * (n - 1.0) * (n - 2.0));
    }

    return energy;
}

// Computes the analytical expressions for the energies when dx=dx(y)
double ecalcc2(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
               double ym, double yp, double sig, double cene, double n, double c3now, int caseind) {

    double energy;
    sig = 1.0; // all distances are written in terms of sigma

    double u1 = pow((c1 * yp + c3now), 1.0 - n) * ((a1 * c1 - a2 * c1) * (n - 2.0) + (b1 - b2) * (c3now + c1 * (n - 1.0) * yp));
    double u2 = -pow((c1 * ym + c3now), 1.0 - n) * ((a1 * c1 - a2 * c1) * (n - 2.0) + (b1 - b2) * (c3now + c1 * (n - 1.0) * ym));

    energy = cene * pow(sig, n - 2.0) * (u1 + u2) / (pow(c1, 2.0) * (n - 2.0) * (n - 1.0));
    return energy;
}

// Computes the analytical expressions for the energies when dx=dx(z)
double ecalcc1(double a1, double a2, double b1, double b2, double c1, double c2, double c3,
               double ym, double yp, double sig, double cene, double n, double c3now, int caseind) {

    double u1, u2, energy;
    sig = 1.0; // all distances are written in terms of sigma

    if (caseind == 1) {
        energy = cene * pow(sig, n - 2.0) * ((yp - ym) * pow((a2 * c2 + c3now), 1.0 - n) + (ym - yp) * pow((a1 * c2 + c3now), 1.0 - n)) / (c2 * (1.0 - n));
    } else if (caseind == 2) {
        u1 = (yp - ym) * pow((a2 * c2 + c3now), 1.0 - n);
        u2 = pow((a1 * c2 + c3now + b1 * c2 * ym), 2.0 - n) - pow((a1 * c2 + c3now + b1 * c2 * yp), 2.0 - n);
        u2 /= (b1 * c2 * (2.0 - n));
        energy = cene * pow(sig, n - 2.0) * (u1 + u2) / (c2 * (1.0 - n));
    } else if (caseind == 3) {
        u1 = (ym - yp) * pow((a1 * c2 + c3now), 1.0 - n);
        u2 = pow((a2 * c2 + c3now + b2 * c2 * yp), 2.0 - n) - pow((a2 * c2 + c3now + b2 * c2 * ym), 2.0 - n);
        u2 /= (b2 * c2 * (2.0 - n));
        energy = cene * pow(sig, n - 2.0) * (u1 + u2) / (c2 * (1.0 - n));
    } else if (caseind == 4) {
        u1 = -pow((a2 * c2 + c3now + b2 * c2 * ym), 2.0 - n) + pow((a2 * c2 + c3now + b2 * c2 * yp), 2.0 - n);
        u1 /= (2.0 * b2 * c2 - b2 * c2 * n);
        u2 = pow((a1 * c2 + c3now + b1 * c2 * ym), 2.0 - n) - pow((a1 * c2 + c3now + b1 * c2 * yp), 2.0 - n);
        u2 /= (2.0 * b1 * c2 - b1 * c2 * n);
        energy = cene * (u1 + u2) * pow(sig, n - 2.0) / (c2 * (1.0 - n));
    }

    return energy;
}

double calc_gen(double sig, double lc, double atre1, double atre2, double repe, double catr1, double catr2,double crep, const std::vector<std::vector<double>>& region, double Epsilon, double* surface_norm) {
    double uatr, urep, utot;
    utot = 0.0;
    sig = sig/10.0; // make A nm since this formula is written in terms of nm
    double dslim = 2.0;
    double slopetol1 = pow(10,-2);
    double slopetol2 = pow(10,-6);
    double yp = maxi_v(region[1][0],region[0][0]);
    double ym = mini_v(region[1][0],region[0][0]);
    double b2 = (region[1][1] - region[0][1]) / (region[1][0] - region[0][0]);
    double a2 = region[0][1] - b2 * region[0][0];
    double b1 = (region[3][1] - region[2][1]) / (region[3][0] - region[2][0]);
    double a1 = region[3][1] - b1 * region[3][0];

    double c1 = -*(surface_norm+1)/(*surface_norm);
    double c2 = -*(surface_norm+2)/(*surface_norm);
    double c3 = -*(surface_norm+3)/(*surface_norm);
    double n;

    //printf("c1:%lf    c2:%lf    c3:%lf\n",c1,c2,c3);
    //printf("b2:%lf    a2:%lf\n",b2,a2);
    //printf("b1:%lf    a1:%lf\n",b1,a1);
    //printf("yp:%lf    ym:%lf    ym-yp:%lf\n",yp,ym,ym-yp);
    //printf("area:%lf\n",abs(((region[0][1] - region[2][1]) + (region[1][1] - region[3][1])) * (region[1][0] - region[0][0]) / 2.0));

    if (abs(ym-yp)<1e-6){return 0.0;}
    if (abs(b1-b2)<1e-8 && abs(a1-a2)<1e-8){return 0.0;}
    //if (abs(region[1][0] - region[0][0])<1e-6){return 0.0;}
    //if (abs(region[3][0] - region[2][0])<1e-6){return 0.0;}

    if (c3 > dslim) {
        if (abs(c2) < slopetol1 && abs(c1) < slopetol1) {
            double area = abs(((region[0][1] - region[2][1]) + (region[1][1] - region[3][1])) * (region[1][0] - region[0][0])) / 2.0;
            n = atre2;
            uatr = catr2 * area * pow(c3, -n);
            urep = 0.0;
        } else if (abs(c1) < slopetol1) {
            int caseindc1;
            if (abs(b1) < slopetol2 && abs(b2) < slopetol2) {
                caseindc1 = 1;
            } else if (abs(b2) < slopetol2) {
                caseindc1 = 2;
            } else if (abs(b1) < slopetol2) {
                caseindc1 = 3;
            } else {
                caseindc1 = 4;
            }
            //cout << "Parameters" << endl;
            //cout << a1 << a2 << b1 << b2 << c1 << c2 << c3 << ym << yp << sig << catr2 << atre2 << c3 << caseindc1 << endl;
            uatr = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3, caseindc1);
            urep = 0.0;
        } else if (abs(c2) < slopetol1) {
            n = atre2;
            uatr = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3, 1);
            urep = 0.0;
        } else {
            double check1 = abs(c1 + b1 * c2);
            double check2 = abs(c1 + b2 * c2);
            int caseindgen;
            if (check1 < slopetol2 && check2 < slopetol2) {
                caseindgen = 1;
            } else if (check1 < slopetol2) {
                caseindgen = 2;
            } else if (check2 < slopetol2) {
                caseindgen = 3;
            } else {
                caseindgen = 4;
            }
            uatr = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3, caseindgen);
            urep = 0.0;
        }
    } else {
        double c3t = dslim;
        if (abs(c2) < slopetol1 && abs(c1) < slopetol1) {
            double area = abs(((region[0][1] - region[2][1]) + (region[1][1] - region[3][1])) * (region[1][0] - region[0][0])) / 2.0;
            n = atre1;
            uatr = catr1 * area * pow(c3, -n);
            double nt2 = atre2;
            double nt1 = atre1;
            double uatrt2 = catr2 * area * pow(c3t, -nt2);
            double uatrt1 = catr1 * area * pow(c3t, -nt1);
            uatr = uatr + uatrt2 - uatrt1;
            n = repe;
            urep = crep * area * pow(c3, -n);
            double ureptemp = crep * area * pow(c3t, -n);
            urep = urep - ureptemp;
        } else if (abs(c1) < slopetol1) {
            int caseindc1;
            if (abs(b1) < slopetol2 && abs(b2) < slopetol2) {
                caseindc1 = 1;
            } else if (abs(b2) < slopetol2) {
                caseindc1 = 2;
            } else if (abs(b1) < slopetol2) {
                caseindc1 = 3;
            } else {
                caseindc1 = 4;
            }
            uatr = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3, caseindc1);
            double uatrt1 = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3t, caseindc1);
            double uatrt2 = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3t, caseindc1);
            uatr = uatr + uatrt2 - uatrt1;
            n = repe;
            urep = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3, caseindc1);
            double ureptemp = ecalcc1(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3t, caseindc1);
            urep = urep - ureptemp;
        } else if (abs(c2) < slopetol1) {
            uatr = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3, 1);
            double uatrt1 = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3t, 1);
            double uatrt2 = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3t, 1);
            uatr = uatr + uatrt2 - uatrt1;
            n = repe;
            urep = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3, 1);
            double ureptemp = ecalcc2(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3t, 1);
            urep = urep - ureptemp;
        } else {
            double check1 = abs(c1 + b1 * c2);
            double check2 = abs(c1 + b2 * c2);
            int caseindgen;
            if (check1 < slopetol2 && check2 < slopetol2) {
                caseindgen = 1;
            } else if (check1 < slopetol2) {
                caseindgen = 2;
            } else if (check2 < slopetol2) {
                caseindgen = 3;
            } else {
                caseindgen = 4;
            }
            uatr = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3, caseindgen);
            double uatrt1 = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr1, atre1, c3t, caseindgen);
            double uatrt2 = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, catr2, atre2, c3t, caseindgen);
            uatr = uatr + uatrt2 - uatrt1;
            n = repe;
            urep = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3, caseindgen);
            double ureptemp = ecalcgen(a1, a2, b1, b2, c1, c2, c3, ym, yp, sig, crep, repe, c3t, caseindgen);
            urep = urep - ureptemp;
        }
    }

    uatr = uatr * Epsilon;
    urep = urep * Epsilon;
    utot = uatr + urep;
    //printf("total energy:%lf\n\n",utot);
    return utot;
}

double surface_energy(double sig, double lc, double atre1, double atre2, double repe, double catr1, double catr2, double crep, double Epsilon,
                      double* y_coord1, double* z_coord1, double* y_coord2, double* z_coord2, int N1, int N2, double* surface_norm) {
    double energy = 0.0;
    bool error_flag = false;

    if (std::abs(surface_norm[0]) < 1e-5) {
        return 0.0;
    }

    //std::cout << "y_coord1 " << y_coord1[0] << " " << y_coord1[1] << " " << y_coord1[2] << " " << y_coord1[3] << std::endl;
    //std::cout << "y_coord2 " << y_coord2[0] << " " << y_coord2[1] << " " << y_coord2[2] << " " << y_coord2[3] << std::endl;
    //std::cout << "z_coord1 " << z_coord1[0] << " " << z_coord1[1] << " " << z_coord1[2] << " " << z_coord1[3] << std::endl;
    //std::cout << "z_coord2 " << z_coord2[0] << " " << z_coord2[1] << " " << z_coord2[2] << " " << z_coord2[3] << std::endl;
    //std::cout << "surface_norm " << surface_norm[0] << " " << surface_norm[1] << " " << surface_norm[2] << " " << surface_norm[3] << std::endl;
    //std::cout << "N1 N2 " << N1 << " " << N2 << std::endl;

    int Num_region = get_regions(y_coord1, z_coord1, y_coord2, z_coord2, N1, N2, surface_norm, error_flag);
    //std::cout << "error_flag_surfaceenergy " << error_flag << std::endl;
    //std::cout << "Num_region " << Num_region << std::endl;

    if (Num_region == 0) { return 0.0; }
    if (error_flag) { return 1e10; } // Representing infinity
    if (Num_region < 0) { return 1e10; } // Representing infinity

    for (int i = 0; i < Num_region; i++) {
        // Check the order of region_z to assign region_coord correctly
        if (region_z[i][2] <= region_z[i][1]) {
            region_coord[0] = { region_y[i][0], region_z[i][0] };
            region_coord[1] = { region_y[i][1], region_z[i][1] };
            region_coord[2] = { region_y[i][2], region_z[i][2] };
            region_coord[3] = { region_y[i][0], region_z[i][0] };
        } else {
            region_coord[0] = { region_y[i][0], region_z[i][0] };
            region_coord[1] = { region_y[i][2], region_z[i][2] };
            region_coord[2] = { region_y[i][1], region_z[i][1] };
            region_coord[3] = { region_y[i][0], region_z[i][0] };
        }

        energy += calc_gen(sig, lc, atre1, atre2, repe, catr1, catr2, crep, region_coord, Epsilon, surface_norm);

        if (region_z[i][2] <= region_z[i][1]) {
            region_coord[0] = { region_y[i][3], region_z[i][3] };
            region_coord[1] = { region_y[i][1], region_z[i][1] };
            region_coord[2] = { region_y[i][2], region_z[i][2] };
            region_coord[3] = { region_y[i][3], region_z[i][3] };
        } else {
            region_coord[0] = { region_y[i][3], region_z[i][3] };
            region_coord[1] = { region_y[i][2], region_z[i][2] };
            region_coord[2] = { region_y[i][1], region_z[i][1] };
            region_coord[3] = { region_y[i][3], region_z[i][3] };
        }

        energy += calc_gen(sig, lc, atre1, atre2, repe, catr1, catr2, crep, region_coord, Epsilon, surface_norm);
    }

    return energy;
}
