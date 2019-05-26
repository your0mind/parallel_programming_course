/// Copyright 2019 Volkov Daniil
#include<tbb/tbb.h>
#include<iostream>
#include<functional>


double MainFunction1(double x, double y) {
    return -x * x - y * y;
}
double TopFunction(double x) {
    return -x * x;
}
double LowerFunction(double x) {
    return -x * x;
}
double oneDimensionalIntegral(double x, double func(double, double), double a, double b, int n);
double TBBIntegration(
    double ax, double bx, int n) {
    double step = (bx - ax) / n;
    double res = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(1, n-1, 100), 0.0,
        [&](const tbb::blocked_range<size_t>& range, double result) -> double {
        size_t begin = range.begin(), end = range.end();
        for (size_t i = begin; i < end; i++) {
            result += oneDimensionalIntegral(ax + i * step, MainFunction1,
            TopFunction(ax + i * step), LowerFunction(ax + i * step), n);
        }
        return result;
    }, std::plus<double>());
    res *= 2;
    res += oneDimensionalIntegral(ax, MainFunction1, TopFunction(ax), LowerFunction(ax), n) +
        oneDimensionalIntegral(bx, MainFunction1, TopFunction(bx), LowerFunction(bx), n);
    res *= (step / 2);
    return res;
}

double Trapezoid(double func(double, double), double ax, double bx,
    double TopLimit(double), double LowerLimit(double), int n) {
    double step = (bx - ax) / n;
    double result = 0;

    for (int i = 1; i < n - 1; i++) {
        result += oneDimensionalIntegral(ax + i * step, func,
         TopLimit(ax + i * step), LowerLimit(ax + i * step), n);
    }
    result *= 2;
    result += oneDimensionalIntegral(ax, func, TopLimit(ax), LowerLimit(ax), n) +
        oneDimensionalIntegral(bx, func, TopLimit(bx), LowerLimit(bx), n);
    result *= (step / 2);
    return result;
}

double oneDimensionalIntegral(double x, double func(double, double),
double b, double a, int n) {
    double step = (b - a) / n;
    double result = 0;
    for (int i = 1; i <= n - 1; i++) {
        result += func(x, a + i * step);
    }
    result *= 2;
    result += func(x, a) + func(x, b);
    result *= (step / 2);
    return result;
}



int main(int argc, char *argv[]) {
    double  SerialTime, ParallelTime;
    int numThreads, n;
    numThreads = 12;
    n = 10000;

    tbb::tick_count start = tbb::tick_count::now();
    std::cout << "Serial result: " << Trapezoid(MainFunction1, 0, 10, TopFunction,
        LowerFunction, n);
    std::cout << "\nSerial Time: " <<
        (SerialTime = ((tbb::tick_count::now()) - start).seconds());
    start = tbb::tick_count::now();
    std::cout << "\nParallel result: " << TBBIntegration(0, 10, 10000);
    std::cout << "\nParallel Time: " <<
        (ParallelTime = ((tbb::tick_count::now()) - start).seconds());
    std::cout << "\nBoost: " << (SerialTime / ParallelTime) << "\nEfficiency: ";
    std::cout << (SerialTime / ParallelTime) / numThreads;
}
