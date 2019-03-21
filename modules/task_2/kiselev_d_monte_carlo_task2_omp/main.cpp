// Copyright Kiselev Denis 2019

#include <cfloat>
#include <functional>
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

double ellipsoid(std::vector<double> args) {
    return -1.0 + pow(args[0] / ELLPS_A, 2)
                + pow(args[1] / ELLPS_B, 2)
                + pow(args[2] / ELLPS_C, 2);
}

double integrateByMonteCarlo(
        std::function<double(std::vector<double>)> func,
        std::vector<std::pair<double, double> > limits, int nPoints) {
    int dimension = limits.size();
    std::vector<std::uniform_real_distribution<> > distrs(dimension);

    double measure = 1.0;
    for (size_t i = 0; i < dimension; i++) {
        measure *= limits[i].second - limits[i].first;
        distrs[i].param(std::uniform_real_distribution<>::param_type(
            limits[i].first, limits[i].second));
    }

    int nPointsInEllipsoid = 0;
    std::default_random_engine generator();
    std::vector<double> args(dimension);
    for (int i = 0; i < nPoints; i++) {
        for (size_t j = 0; j < dimension; j++)
            args[j] = distrs[j](generator);
        double value = func(args);
        if (value <= 0) nPointsInEllipsoid++;
    }

    double avgValueOfHits = static_cast<double>(nPointsInEllipsoid) / nPoints;
    return measure * avgValueOfHits;
}

double integrateByMonteCarloParallel(
        std::function<double(std::vector<double>)> func,
		std::vector<std::pair<double, double> > limits, int nPoints) {
    std::default_random_engine generator;
    std::vector<std::uniform_real_distribution<> > distrs(limits.size());
    for (size_t i = 0; i < limits.size(); i++) {
        distrs[i].param(std::uniform_real_distribution<>::param_type(
            limits[i].first, limits[i].second));
    }

    int nPointsInEllipsoid = 0;
    std::vector<double> args(limits.size());
    #pragma omp parallel firstprivate(args) private(generator) reduction(+: nPointsInEllipsoid)
    {
        generator.seed(omp_get_thread_num());
        #pragma omp for schedule(static)
        for (int i = 0; i < nPoints; i++) {
            for (size_t j = 0; j < args.size(); j++)
                args[j] = distrs[j](generator);
            double value = func(args);
            if (value <= 0) nPointsInEllipsoid++;
        }
    }
    double avgValueOfHits = static_cast<double>(nPointsInEllipsoid) / nPoints;
    double measure = 1.0;
    for (auto limit : limits)
        measure *= limit.second - limit.first;
    return measure * avgValueOfHits;
}

int main(int argc, char *argv[]) {
    int nPoints = (argc == 1) ? DEFAULT_NPOINTS : atoi(argv[1]);
    int maxThreads = omp_get_max_threads();

    // Sequential code
    double t1 = omp_get_wtime();
    omp_set_num_threads(1);
    double result = inte(ellipsoid, { { X1, X2 }, { Y1, Y2 }, { Z1, Z2 } }, nPoints);
    double seqTime = omp_get_wtime() - t1;
    std::cout.precision(8);
    std::cout << std::fixed;
    std::cout << "Sequential alg:" << std::endl
              << "\tResult: " << result << std::endl
              << "\tTime: " << seqTime << " sec" << std::endl << std::endl;

    // Parallel code
    t1 = omp_get_wtime();
    omp_set_num_threads(maxThreads);
    result = integrate(ellipsoid, { { X1, X2 }, { Y1, Y2 }, { Z1, Z2 } }, nPoints);
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
