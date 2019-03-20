// Copyright Kiselev Denis 2019

#include <cfloat>
#include <iostream>
#include <random>
#include <omp.h>

#define DEFAULT_NPOINTS 100000

// Ellipsoid options
#define ELLPS_A 1.0
#define ELLPS_B 2.0
#define ELLPS_C 3.0

// Limits of integration
#define X1 -1.0
#define X2 1.0
#define Y1 -2.0
#define Y2 2.0
#define Z1 -3.0
#define Z2 3.0

// This is a simple elipsoid function, which we
// will integrate using the Monte Carlo method. 
double ellipsoid(double x, double y, double z) {
    return - 1.0 + pow(x / ELLPS_A, 2)
                 + pow(y / ELLPS_B, 2)
                 + pow(z / ELLPS_C, 2);
}

int main(int argc, char *argv[]) {
    // Creating a generator and distributions
    // for random double numbers in the ranges
    std::default_random_engine generator;
    std::uniform_real_distribution<> xAxisDistr(X1, X2);
    std::uniform_real_distribution<> yAxisDistr(Y1, Y2);
    std::uniform_real_distribution<> zAxisDistr(Z1, Z2);
    int nPoints = (argc == 1) ? DEFAULT_NPOINTS : atoi(argv[1]);



    // Sequential code
    double t1 = omp_get_wtime();
    int nPointsInEllipsoid = 0;
    generator.seed(0);
    for (int i = 0; i < nPoints; i++) {
        double value = ellipsoid(xAxisDistr(generator),
                                 yAxisDistr(generator),
                                 zAxisDistr(generator));
        if (value <= 0) nPointsInEllipsoid++;
    }
    double avgValueOfHits = static_cast<double>(nPointsInEllipsoid) / nPoints;
    double result = (X2 - X1) * (Y2 - Y1) * (Z2 - Z1) * avgValueOfHits;
    double seqTime = omp_get_wtime() - t1;
    std::cout.precision(8);
    std::cout << std::fixed;
    std::cout << "Sequential alg:" << std::endl
              << "\tResult: " << result << std::endl
              << "\tTime: " << seqTime << " sec" << std::endl << std::endl;



    // Parallel code
    t1 = omp_get_wtime();
    nPointsInEllipsoid = 0;
    #pragma omp parallel private(generator) reduction(+: nPointsInEllipsoid)
    {
        generator.seed(omp_get_thread_num());
        #pragma omp for schedule(static)
        for (int i = 0; i < nPoints; i++) {
            double value = ellipsoid(xAxisDistr(generator),
                                     yAxisDistr(generator),
                                     zAxisDistr(generator));
            if (value <= 0) nPointsInEllipsoid++;
        }
    }
    avgValueOfHits = static_cast<double>(nPointsInEllipsoid) / nPoints;
    result = (X2 - X1) * (Y2 - Y1) * (Z2 - Z1) * avgValueOfHits;
    double parTime = omp_get_wtime() - t1;
    std::cout << "Parallel alg:" << std::endl
              << "\tResult: " << result << std::endl
              << "\tTime: " << parTime << " sec" << std::endl << std::endl;



    double realRes = 4.0 / 3.0 * std::acos(-1) * ELLPS_A * ELLPS_B * ELLPS_C;
    double acceleration = seqTime / parTime;
    std::cout << "Real result: " << realRes << std::endl
              << "Acceleration: " << acceleration << std::endl;
              
    return 0;
}
