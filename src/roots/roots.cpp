#include <iostream>
#include <iomanip>
#include "roots.hpp"
#include <cmath>



bool bisection(std::function<double(double)> f,double a, double b, double *root) {
    double lower, upper, midpoint;
    if (f(a)*f(b) > 0) {
        return false;
    }

    if (a < b) {
        lower = a;
        upper = b;
    }
    else {
        lower = b;
        upper = a; 
    }
   
        int iteration = 0;
        int max_iter = 1e6;
        double epsilon = 1e-6;
        midpoint = (lower + upper)/2;
        double f_lower = f(lower);
        double f_upper = f(upper);

    while (std::fabs(f(midpoint)) > epsilon) { 
        iteration ++;
        if (iteration >= max_iter) return false;
        midpoint = (lower + upper)/2;
        double f_midpoint = f(midpoint);

        /*
        std::cout << "Iteration: " << iteration << std::endl;
        std::cout << "Lower: " << lower << std::endl;
        std::cout << "Upper: " << upper << std::endl;
        std::cout << "Midpoint: " << midpoint << std::endl;
        std::cout << "f_Midpoint: " << f_midpoint << "\n\n";
        */

        if (f_lower * f_midpoint < 0) {
            // root on left interval
            upper = midpoint;
            f_upper = f_midpoint;
        }
        else { 
            // root on right interval
            lower = midpoint;
            f_lower = f_midpoint;
        }
    } 
    *root = midpoint;
    return true;
}

bool regula_falsi(std::function<double(double)> f, double a, double b, double *root) {
    double lower, upper, midpoint;
    if (f(a)*f(b) > 0) {
        return false;
    }

    if (a < b) {
        lower = a;
        upper = b;
    }
    else {
        lower = b;
        upper = a; 
    }
   
        int iteration = 0;
        int max_iter = 1e6;
        double epsilon = 1e-6;
        double f_lower = f(lower);
        double f_upper = f(upper);
        midpoint = (lower*f_upper-upper*f_lower)/(f_upper-f_lower);


    while (std::fabs(f(midpoint)) > epsilon) { 

        iteration ++;
        if (iteration >= max_iter) return false;
        midpoint = (lower*f_upper-upper*f_lower)/(f_upper-f_lower);
        double f_midpoint = f(midpoint);

       /*
        std::cout << "Iteration: " << iteration << std::endl;
        std::cout << "Lower: " << lower << std::endl;
        std::cout << "Upper: " << upper << std::endl;
        std::cout << "Midpoint: " << midpoint << std::endl;
        std::cout << "f_Midpoint: " << f_midpoint << "\n\n";
        */

        if (f_lower * f_midpoint < 0) {
            // root on left interval
            upper = midpoint;
            f_upper = f_midpoint;
            f_lower /= 2; //illinois method to avoid gettinig stuck
        }
        else { 
            // root on right interval
            lower = midpoint;
            f_lower = f_midpoint;
            f_upper /= 2;
        }
    } 
    *root = midpoint;
    return true;
}

bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root) {
    double x0 = c;
    double x1;
    double diff = 1;

    int iteration = 0;
    int max_iter = 1e6;
    double epsilon = 1e-6;
    double upsilon = 1e-12;

    while (std::fabs(diff) > epsilon) {
        if (iteration >= max_iter) return false;
        iteration ++;
        if (std::fabs(g(x0)) < upsilon) return false;
        x1 = x0 - f(x0)/g(x0);
        diff = x1 - x0;
        if (x1 < a || x1 > b) return false;
        /*
        std::cout << "Iteration: " << iteration << std::endl;
        std::cout << "Xn: " << x0 << std::endl;
        std::cout << "Xn+1: " << x1 << "\n\n";
        */
        x0 = x1;
    }
    *root = x0;
    return true;
}

bool secant(std::function<double(double)> f, double a, double b, double c, double *root) {
    double x0 = c;
    double x1 = c+1e-3;
    double x2;
    double diff = 1;

    int iteration = 0;
    int max_iter = 1e6;
    double epsilon = 1e-6;
    double upsilon = 1e-12;

    while (std::fabs(diff) > epsilon) {
        if (iteration >= max_iter) return false;
        iteration ++;
        double f_x1 = f(x1);
        double f_x0 = f(x0);
        if (std::fabs(f_x1-f_x0) < upsilon) return false;
        x2 = x1 - f_x1*(x1-x0)/(f_x1-f_x0);
        diff = x2 - x1;
        if (x2 < a || x2 > b) return false;
        
        /*
        std::cout << "Iteration: " << iteration << std::endl;
        std::cout << "Xn-1: " << x0 << std::endl;
        std::cout << "Xn: " << x1 << std::endl;
        std::cout << "Xn+1: " << x2 << "\n\n";
        */
       
        x0 = x1;
        x1 = x2;
    }
    *root = x1;
    return true;
}

