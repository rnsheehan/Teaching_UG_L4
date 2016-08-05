#ifndef INTEGRATION_LIB_H
#define INTEGRATION_LIB_H

// Functions for computing approximate values of definite integrals
// it is assumed that the nodes are evenly spaced
// R. Sheehan 20 - 1 - 2014

namespace integration_funcs{

	double rect_rule(double (*f)(double x), double xl, double xu, int n_nodes); // rectangular rule
	double trap_rule(double (*f)(double x), double xl, double xu, int n_nodes); // trapezoidal rule
	double simp_rule(double (*f)(double x), double xl, double xu, int n_nodes); // Simpson's rule
	double gauss_quad(double (*f)(double x), double xl, double xu); // Gaussian quadrature

	// Function for computing an approximate value of the integral of a data set
	// The data is integrated as provided
	double data_integrate_trap_rule(double *x_vals, double *y_vals, int n_data); // trapezoidal rule integral for a data set

}

#endif