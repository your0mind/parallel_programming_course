// Copyright Kiselev Denis 2019

#include <iostream>
#include <random>
#include <ctime>

/* This is a simple parabola function, which we
   will integrate using the Monte Carlo method. */
double parabola(double x) {
    return x * x;
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    // Getting program call arguments
    int nPoints = atoi(argv[1]);
    double leftBorder = atof(argv[2]);
    double rightBorder = atof(argv[3]);

    // Creating a generator and distribution
    // for random double numbers in the range
    std::default_random_engine generator(time(NULL));
    std::uniform_real_distribution<double> distribution(leftBorder,
                                                        rightBorder);

    // The sum of the random values of our parabola
    double sum = 0.0;
    for (int i = 0; i < nPoints; i++) {
        double randX = distribution(generator);
        sum += parabola(randX);
    }

    // Calculation of the result according to the Monte Carlo method
    double avgValue = sum / nPoints;
    double result = (rightBorder - leftBorder) * avgValue;

    std::cout.precision(8);
    std::cout << "Result: " << std::fixed << result << std::endl;
    return 0;
}
