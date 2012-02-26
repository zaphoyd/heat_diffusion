#include "../heat.hpp"
#include "../object1d.hpp"
#include "../object2d.hpp"
#include "../object3d.hpp"

#include <iostream>
#include <string>

bool print1d(const object1d& o,size_t ts) {
	std::string sep;
	
	std::cout << "[";
	
	for (size_t i = 0; i < o.nx(); i++) {
		std::cout << sep << o[i];
		sep = ",";
	}
	
	std::cout << "]" << std::endl;
	
	return true;
}

bool print2d(const object2d& o,size_t ts) {
	std::string sep;
	
	std::cout << "[" << std::endl;
	
	for (size_t i = 0; i < o.nx(); i++) {
		std::cout << "[";
		
		sep = "";
		
		for (size_t j = 0; j < o.ny(); j++) {
			std::cout << sep << o[i][j];
			sep = ",";
		}
		std::cout << "]" << std::endl;
	}
	
	std::cout << "]" << std::endl;
	
	return true;
}

bool print3d(const object3d& o,size_t ts) {
	size_t zslice = o.nz()/2;
	
	std::string sep;
	
	std::cout << "[" << std::endl;
	
	for (size_t i = 0; i < o.nx(); i++) {
		std::cout << "[";
		
		sep = "";
		
		for (size_t j = 0; j < o.ny(); j++) {
			std::cout << sep << o[i][j][zslice];
			sep = ",";
		}
		std::cout << "]" << std::endl;
	}
	
	std::cout << "]" << std::endl;
	
	return true;
}

int main() {
	double lx = 1.0;
	double ly = 1.0;
	double lz = 1.0;
	
	size_t nx = 3;
	size_t ny = 3;
	size_t nz = 3;
	
	double alpha = 0.0005;
	double dt = 0.001;
	
	object1d o1(lx,nx,alpha);
    object1d s1(lx,nx,0.0);
    o1.init(GAUSSIAN);
    s1.init(FLAT);
	
	object2d o2(lx,ly,nx,ny,alpha);
    object2d s2(lx,ly,nx,ny,0.0);
    o2.init(GAUSSIAN);
    s2.init(FLAT);
	
	object3d o3(lx,ly,lz,nx,ny,nz,alpha);
    object3d s3(lx,ly,lz,nx,ny,nz,0.0);
    o3.init(GAUSSIAN);
    s3.init(FLAT);

	print1d(o1,0);
	print1d(s1,0);
	std::cout << "diff: " << o1.compute_mean_abs_diff(s1) << std::endl;

	/*o2.ftcs(
		5, dt,			// timesteps, dt
		CONSTANT, 0.0,	// boundaries
		s2,				// source term
		&print2d,		// callback function
		5				// callback interval
	);*/
	
	/*o2.crank_nicholson(
		3, dt,			// timesteps, dt
		CONSTANT, 0.0,	// boundaries
		s2,				// source term
		&print2d,		// callback function
		5				// callback interval
	);*/

	o1.jacobi (
		5, dt,
		CONSTANT, 0.0,
		s1,
		&print1d,
		1
	);
	
	/*o1.crank_nicholson(
		5, dt,			// timesteps, dt
		CONSTANT, 0.0,	// boundaries
		s1,				// source term
		&print1d,		// callback function
		5				// callback interval
	);*/
	
	return 0;
}
