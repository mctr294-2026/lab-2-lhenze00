#include <iostream>
#include <iomanip>
#include "roots.hpp"
#include <cmath>


double poly1(double x)
{
    // initial bracket [-200, 300]
    return x * x * x - x * x + 2;
}

double poly1_deriv(double x)
{
    return 3 * x * x - 2 * x;
}

double poly2(double x)
{
    // initial bracket [-1,1]
    return 2 * x * x * x - 4 * x * x + 3 * x;
}

double poly2_deriv(double x)
{
    return 6 * x * x - 8 * x + 3;
}


double bisection(std::function<double(double)> f,double a, double b) {
    double lower, upper, midpoint;
    if (f(a)*f(b) > 0) {
        std::cerr << "The function is not bracketed" << std::endl;
        return -1;
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
        double epsilon = std::pow(10,-6);
        midpoint = (lower + upper)/2;
        double f_lower = f(lower);
        double f_upper = f(upper);

    while (std::fabs(f(midpoint)) > epsilon) { 
        iteration ++;
        midpoint = (lower + upper)/2;
        double f_midpoint = f(midpoint);

        bool debug = true;
        if (debug) {
            std::cout << "Iteration: " << iteration << std::endl;
            std::cout << "Lower: " << lower << std::endl;
            std::cout << "Upper: " << upper << std::endl;
            std::cout << "Midpoint: " << midpoint << std::endl;
            std::cout << "f_Midpoint: " << f_midpoint << "\n\n";
        }

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
    return midpoint;
}

double regula_falsi(std::function<double(double)> f, double a, double b) {
    double lower, upper, midpoint;
    if (f(a)*f(b) > 0) {
        std::cerr << "The function is not bracketed" << std::endl;
        return -1;
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
        double epsilon = std::pow(10,-6);
        double f_lower = f(lower);
        double f_upper = f(upper);
        midpoint = (lower*f_upper-upper*f_lower)/(f_upper-f_lower);


    while (std::fabs(f(midpoint)) > epsilon) { 
        iteration ++;
        midpoint = (lower*f_upper-upper*f_lower)/(f_upper-f_lower);
        double f_midpoint = f(midpoint);

        bool debug = true;
        if (debug) {
            std::cout << "Iteration: " << iteration << std::endl;
            std::cout << "Lower: " << lower << std::endl;
            std::cout << "Upper: " << upper << std::endl;
            std::cout << "Midpoint: " << midpoint << std::endl;
            std::cout << "f_Midpoint: " << f_midpoint << "\n\n";
        }

        if (f_lower * f_midpoint < 0) {
            // root on left interval
            upper = midpoint;
            f_upper = f_midpoint;
            f_lower /= 2;
        }
        else { 
            // root on right interval
            lower = midpoint;
            f_lower = f_midpoint;
            f_upper /= 2;
        }
    } 
    return midpoint;
}


double newtons_method(std::function<double(double)> f, std::function<double(double)> g, double a) {
    double x0 = a;
    double x1;
    double diff = 1;

    int iteration = 0;
    int max_iter = 1000;
    double epsilon = std::pow(10,-6);
    double upsilon = 1e-12;

    while (std::fabs(diff) > epsilon && iteration < max_iter) {
        iteration ++;
         if (std::fabs(g(x0)) < upsilon) {
            throw std::runtime_error ("Newton's method failed: Divide by zero error");
        }
        x1 = x0 - f(x0)/g(x0);
        diff = x1 - x0;
        std::cout << "Iteration: " << iteration << std::endl;
        std::cout << "Xn: " << x0 << std::endl;
        std::cout << "Xn+1: " << x1 << "\n\n";
        x0 = x1;
    }
    return x1;
}


double secant_method(std::function<double(double)> f, double a, double b, double c) {
    double x0 = c;
    double x1 = c+1e-3;
    double x2;
    double diff = 1;

    int iteration = 0;
    int max_iters = 1000;
    double epsilon = 1e-6;
    double upsilon = 1e-12;

    while (std::fabs(diff) > epsilon && iteration < max_iters) {
        iteration ++;
        double f_x1 = f(x1);
        double f_x0 = f(x0);
        if (std::fabs(f_x1-f_x0) < upsilon) {
            throw std::runtime_error ("Secant method failed: Divide by zero error");
        }
        x2 = x1 - f_x1*(x1-x0)/(f_x1-f_x0);
        diff = x2 - x1;
        std::cout << "Iteration: " << iteration << std::endl;
        std::cout << "Xn-1: " << x0 << std::endl;
        std::cout << "Xn: " << x1 << std::endl;
        std::cout << "Xn+1: " << x2 << "\n\n";
        
        x0 = x1;
        x1 = x2;
    }
    return x1;
}


double poly3(double x) {
    return x*x-4;
}





int main() {
    //std::cout << "Running Bisection" << std::endl;
    //std::cout << bisection(poly1,-200,300) << std::endl;

    //std::cout << "Running Regula-Falsi" << std::endl;
    //std::cout << regula_falsi(poly1,-200,300) << std::endl;

    //std::cout << "Running Newton-Raphson" << std::endl;
    //std::cout << newtons_method(poly1,poly1_deriv,300) << std::endl;

    std::cout << "Running Secant" << std::endl;
    std::cout << secant_method(poly3, 0.0,3.0,1.0) << std::endl;

    return 0;
}





/*
bool secant(std::function<double(double)> f, double a, double b, double c) {
    double x0 = a;
    double x1 = c;
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
        
        std::cout << "Iteration: " << iteration << std::endl;
        std::cout << "Xn-1: " << x0 << std::endl;
        std::cout << "Xn: " << x1 << std::endl;
        std::cout << "Xn+1: " << x2 << "\n\n";
        
        x0 = x1;
        x1 = x2;
    }
    return true;
}

int main() {
    std::cout << "Running secant: " << secant(poly3,0.0,3.0,1.0) << "\n";

    return 0;
}

*/