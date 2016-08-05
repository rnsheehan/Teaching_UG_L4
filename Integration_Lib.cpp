#ifndef ATTACH_H
#include "Attach.h"
#endif

// Functions for computing approximate values of definite integrals
// it is assumed that the nodes are evenly spaced
// R. Sheehan 20 - 1 - 2014

double rect_rule(double (*f)(double x), double xl, double xu, int n_nodes)
{
	// rectangular rule approximation to the integral of f between xl and xu
	// the function is evaluated at the left endpoint of each subinterval 
	// double (*f)(double x) is a pointer to a function of type double that has a single double input parameter

	// You cannot test for floating point equality with xu == xl
	// this is why we test for fabs(xu - xl) > EPS
	// Floating point numbers are only "equal" when their difference is less than some small number

	bool c1 = (n_nodes > 1 ? true : false); // test to see if n_nodes > 1
	bool c2 = ( fabs(xu - xl) > EPS ? true : false ); // test to see if xl != xu

	if(c1 && c2){

		// Only attempt to compute the integral if more than 1 nodes are used
		// Also require that xu != xl
		// if c1 and c2 are both true the code will proceed

		int N = n_nodes-1; // number of terms in the sum

		double lower_bound = min(xu,xl); // lower edge of the integration region
		double upper_bound = max(xu,xl); // upper edge of the integration region
		double h = (upper_bound - lower_bound) / (static_cast<double>(N)); // compute node spacing
		double eval_pos; // f will be evaluated at this position
		double t; // use this variable to store the value of the terms in the sum as they are computed
		double s; // use this variable to store the value of the sum which approximates the integral

		eval_pos = xl; // initialise the value at which the function will be evaluated
	
		s = 0.0; // initialise the sum value

		for(int i = 1; i <= N; i++){
	
			t = (*f)(eval_pos); // evalute f at the i^{th} node position

			s += t; // compute the value of the sum

			eval_pos += h; // increment the node position
		}

		s *= h; // multiply the value of the sum by the node spacing

		return s; 
	}
	else{
		// the integral is not evaluated

		if(!c1){
			cout<<"\nIntegral was not evaluated because the number of nodes was not greater than one\n";
		}

		if(!c2){
			cout<<"\nIntegral was not evaluated because the endpoints of the region of integration are equal\n";
		}

		return 0.0; 

	}
}

double trap_rule(double (*f)(double x), double xl, double xu, int n_nodes)
{
	// trapezoidal rule approximation to the integral of f between xl and xu
	// double (*f)(double x) is a pointer to a function of type double that has a single double input parameter

	// You cannot test for floating point equality with xu == xl
	// this is why we test for fabs(xu - xl) > EPS
	// Floating point numbers are only "equal" when their difference is less than some small number

	bool c1 = (n_nodes > 1 ? true : false); // test to see if n_nodes > 1
	bool c2 = ( fabs(xu - xl) > EPS ? true : false ); // test to see if xl != xu

	if(c1 && c2){

		// Only attempt to compute the integral if more than 1 nodes are used
		// Also require that xu != xl
		// if c1 and c2 are both true the code will proceed

		int N = n_nodes-1; // number of terms in the sum

		double lower_bound = min(xu,xl); // lower edge of the integration region
		double upper_bound = max(xu,xl); // upper edge of the integration region
		double h = (upper_bound - lower_bound) / (static_cast<double>(N)); // compute node spacing
		double half_h = 0.5*h; // half the inter-node spacing
		double eval_pos; // f will be evaluated at this position
		double t; // use this variable to store the value of the terms in the sum as they are computed
		double s; // use this variable to store the value of the sum which approximates the integral

		eval_pos = xl + h; // initialise the value at which the function will be evaluated
	
		s = 0.0; // initialise the sum value

		for(int i = 2; i <= N; i++){
	
			t = (*f)(eval_pos); // evalute f at the i^{th} node position

			s += t; // compute the value of the sum

			eval_pos += h; // increment the node position
		}

		s *= 2.0; // Multiply the sum value by two 

		s += ( (*f)(xl) + (*f)(xu) ); // add function values at the first and last node positions

		s *= half_h; // multiply the value of the sum by half the node spacing

		return s; 
	}
	else{
		
		// the integral is not evaluated

		if(!c1){
			cout<<"\nIntegral was not evaluated because the number of nodes was not greater than one\n";
		}

		if(!c2){
			cout<<"\nIntegral was not evaluated because the endpoints of the region of integration are equal\n";
		}

		return 0.0; 

	}
}

double simp_rule(double (*f)(double x), double xl, double xu, int n_nodes)
{
	// Simpson's rule approximation to the integral of f between xl and xu
	// double (*f)(double x) is a pointer to a function of type double that has a single double input parameter

	// You cannot test for floating point equality with xu == xl
	// this is why we test for fabs(xu - xl) > EPS
	// Floating point numbers are only "equal" when their difference is less than some small number

	bool c1 = (n_nodes > 1 ? true : false); // test to see if n_nodes > 1
	bool c2 = ( fabs(xu - xl) > EPS ? true : false ); // test to see if xl != xu
	bool c3 = ( (n_nodes % 2) != 0 ? true : false ); // Simpson's rule also requires the n_nodes be odd, so we require a non-zero remainder upon division by 2

	if(c1 && c2 && c3){

		// Only attempt to compute the integral if more than 1 nodes are used
		// Also require that xu != xl
		// if c1 and c2 and c3 are all true the code will proceed

		int N = n_nodes-1; // number of terms in the sum

		double lower_bound = min(xu,xl); // lower edge of the integration region
		double upper_bound = max(xu,xl); // upper edge of the integration region
		double h = (upper_bound - lower_bound) / (static_cast<double>(N)); // compute node spacing
		double third_h = (h / 3.0); // one-third the inter-node spacing
		double eval_pos; // f will be evaluated at this position
		double t; // use this variable to store the value of the terms in the sum as they are computed
		double even_sum; // use this variable to store the result of the sum over the function values at the even numbered nodes 
		double odd_sum; // use this variable to store the result of the sum over the function values at the odd numbered nodes
		double s; // use this variable to store the value of the sum which approximates the integral

		eval_pos = xl + h; // initialise the value at which the function will be evaluated
	
		// initialise the sum values
		s = 0.0; 
		even_sum = 0.0; 
		odd_sum = 0.0; 

		for(int i = 2; i <= N; i++){
	
			t = (*f)(eval_pos); // evalute f at the i^{th} node position

			if(i%2 == 0){
				// even numbered node
				even_sum += t; 
			}
			else{
				// odd numbered node
				odd_sum += t; 
			}

			eval_pos += h; // increment the node position
		}

		even_sum *= 4.0; // multiply the sum over the even numbered nodes by 4

		odd_sum *= 2.0; // multiply the sum over the odd numbered nodes by 2

		s = ( (*f)(xl) + (*f)(xu) + even_sum + odd_sum ); // add function values at the first and last node positions to the even and odd sums

		s *= third_h; // multiply the total value of the sum by one-third the node spacing

		return s; 
	}
	else{
		
		// the integral is not evaluated
		
		if(!c1){
			cout<<"\nIntegral was not evaluated because the number of nodes was not greater than one\n";
		}

		if(!c2){
			cout<<"\nIntegral was not evaluated because the endpoints of the region of integration are equal\n";
		}

		if(!c3){
			cout<<"\nIntegral was not evaluated because the number of nodes was even\n";
		}

		return 0.0; 
	}
}

double gauss_quad(double (*f)(double x), double xl, double xu)
{
	// Gaussian quadrature approximation to the integral of f between xl and xu
	// 20-point Gauss-Legendre quadrature is accurate to at least 14 decimal places
	// double (*f)(double x) is a pointer to a function of type double that has a single double input parameter

	// You cannot test for floating point equality with xu == xl
	// this is why we test for fabs(xu - xl) > EPS
	// Floating point numbers are only "equal" when their difference is less than some small number

	bool c2 = ( fabs(xu - xl) > EPS ? true : false ); // test to see if xl != xu

	if(c2){

		// Require that xu != xl
		// if c2 is true the code will proceed

		static const int N = 10; // number of nodes used to perform Gaussian quadrature

		double lower_bound = min(xu,xl); // lower edge of the integration region
		double upper_bound = max(xu,xl); // upper edge of the integration region
		double xr,xm,dx,s; // variables used to transform the region of integration

		// node list used to perform Gaussian quadrature
		static double x[]={0.0, 0.076526521133497333755, 0.227785851141645078080, 0.373706088715419560673, 0.510867001950827098004,
			0.636053680726515025453, 0.746331906460150792614, 0.839116971822218823395, 0.912234428251325905868,
			0.963971927277913791268, 0.993128599185094924786}; // Nodes taken from Abramowitz and Stegun 

		// weight list used to perform Gaussian quadrature
		static double w[]={0.0, 0.152753387130725850698, 0.149172986472603746788, 0.142096109318382051329, 0.131688638449176626898,
			0.118194531961518417312, 0.101930119817240435037, 0.083276741576704748725, 0.062672048334109063570,
			0.040601429800386941331, 0.017614007139152118312}; // Weights taken from Abramowitz and Stegun

		xm = 0.5*(upper_bound + lower_bound);
		
		xr = 0.5*(upper_bound - lower_bound);
		
		s = 0.0; // initialise the sum values

		for(int j = 1; j <= N ; j++){

			dx = xr*x[j];

			s += w[j]*( f(xm + dx) + f(xm - dx) );
		}

		s *= xr; // Scale the answer back to the range of integration

		return s; 
	}
	else{
		
		// the integral is not evaluated

		if(!c2){
			cout<<"\nIntegral was not evaluated because the endpoints of the region of integration are equal\n";
		}

		return 0.0; 

	}
}

double data_integrate_trap_rule(double *x_vals, double *y_vals, int n_data)
{
	// trapezoidal rule integral for a data set of length n_data
	// it is assumed that the data in x_vals is ordered

	bool c1 = (n_data > 1 ? true : false); // test to see if n_data > 1

	if(c1){
	
		int N = n_data-1; // number of terms in the sum

		double hi, t, s; // variables used to compute the sum

		s = 0.0; // initialise the value of the sum

		for(int i=1; i<= N; i++){
			
			hi = ( x_vals[i+1] - x_vals[i] ); // spacing between node i+1 and node i

			t = ( y_vals[i] + y_vals[i+1] ); // sum of f(x_{i}) and f(x_{i+1})

			s += ( hi * t ); // compute the sum
		}

		s *= 0.5; // multiply the sum by half

		return s; 
	}
	else{

		// the integral is not evaluated

		if(!c1){
			cout<<"\nIntegral was not evaluated because the number of nodes was not greater than one\n";
		}

		return 0.0; 
	}
}