#include "../heat.hpp"
#include "../object1d.hpp"
#include "../object2d.hpp"
#include "../object3d.hpp"

#include <chrono>
#include <map>
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

// use this if you don't want sub-steps printed.
bool empty(const object3d& o,size_t ts) {
	return true;
}

int main() {
	double lx = 1.0;
	double ly = 1.0;
	double lz = 1.0;
	
	size_t nx = 70;
	size_t ny = 70;
	size_t nz = 70;
	
	double alpha = 0.0005;
	
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
	
	/*std::cout << nx << "x" << ny << "x" << nz << " FTCS simulation:" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	o3.ftcs (10000, 0.001, // timesteps, dt
		CONSTANT, 0.0, // boundary conditions
		s3, // source term
		&empty, 10 // callback settings
	);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = end-start;
	std::cout << (duration.count()/1000000000.0) << " seconds" << std::endl;*/
	
	std::map<disc_method,std::string> sims;
	
	sims[JACOBI] = "Jacobi";
	sims[GAUSS_SEIDEL] = "Gauss-Seidel";
	sims[SOR] = "Successive Over-relaxation";
	
	for (auto &i : sims) {
		std::cout << nx << "x" << ny << "x" << nz << " " << i.second << " simulation:" << std::endl;
		auto start = std::chrono::high_resolution_clock::now();
		o3.iterative (i.first, 1.65,
			10, 1.0, // timesteps, dt
			CONSTANT, 0.0, // boundary conditions
			s3, // source term
			&empty, 10 // callback settings
		);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = end-start;
		std::cout << (duration.count()/1000000000.0) << " seconds" << std::endl;
	}
	
	
	return 0;
}
