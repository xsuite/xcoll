// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SIMPSON_H
#define XCOLL_GEOM_SIMPSON_H

// Simpson's Rule function for numerical integration
double simpson(double (*func)(double,void*), double a, double b, int n, void *params) {
    if (n % 2 == 1) {
        n++;                                         // Requires an even number of subintervals
    }
    double h = (b - a) / n;                          // Step size
    double sum = func(a, params) + func(b, params);  // f(a) + f(b)

    // Sum the middle terms
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 2 == 0) {
            sum += 2.0 * func(x, params);            // Even indices (except the endpoints) are multiplied by 2
        } else {
            sum += 4.0 * func(x, params);            // Odd indices are multiplied by 4
        }
    }
    sum *= h / 3;                                    // Multiply by step size / 3
    return sum;
}

#endif /* XCOLL_GEOM_SIMPSON_H */