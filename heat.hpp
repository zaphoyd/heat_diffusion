#ifndef HEAT_HPP
#define HEAT_HPP

#include <cmath>
#include <vector>

enum smode {
	JSON = 0,
	BINARY = 1
};

enum boundary_style {
	CONSTANT = 0,	// constant
	PERIODIC = 1	// wraps around
};

enum initial_condition {
	FLAT = 0,
	GAUSSIAN = 1,
	GAUSSIAN_NOISE = 2,
	SPIKE = 3
};

enum disc_method {
	FTCS = 0,
	CRANK_NICHOLSON = 1,
	JACOBI = 2,
	GAUSS_SEIDEL = 3,
	SOR = 4,
	MULTIGRID = 5
};

enum relax_scheme {
	RS_JACOBI = 0,
	RS_GAUSS_SEIDEL = 1,
	RS_RB_GAUSS_SEIDEL = 2
};

/// return a gaussian distribution
/* Returns an vector populated with a gaussian distribution scaled to between
 * zero and one with a specified resolution.
 *
 * @param nx number of values in output vector
 * @return filled vector of doubles
 */
std::vector<double> gaussian(size_t nx) {
	std::vector<double> g(nx);
		
	double sx = 1.0 / (nx+1);
	for (size_t i = 0; i < nx; i++) {
		g[i] = exp(-1*(pow((5*(double(i)/double(nx))-2.5),2)));
	}
	
	return g;
}

#endif // HEAT_HPP
