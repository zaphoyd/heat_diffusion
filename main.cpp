#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <vector>


class object1d {
public:
	typedef std::vector<double> data_type;

	object1d(double len,size_t steps, double alpha) : m_data(steps),m_x_len(len),m_alpha(alpha) {}
	
	void init() {
		double s = m_x_len / (m_data.size()+1);
		for (size_t i = 0; i < m_data.size(); i++) {
			m_data[i] = exp(-1*pow((5*(i+1)*s-2.5),2));
		}
	}

	double x_len() const {
		return m_x_len;
	}

	size_t x_steps() const {
		return m_data.size();
	}

	double alpha() const {
		return m_alpha;
	}
	
	const std::vector<double> data() const {
		return m_data;
	}
	
	double& operator[] (const size_t index) {
		return m_data[index];
	}

	void set(const std::vector<double>& new_data) {
		m_data = new_data;
	}

	std::string str() const {
		std::stringstream o;
		std::string sep;

		for (const double& v : m_data) {
			o << sep << v;
			sep = " ";
		}

		return o.str();
	}
private:
	data_type m_data;				// vector representing evenly spaced grid points
	double m_x_len;					// length of x dimension (m)
	double m_alpha;					// thermal diffusivity (m^2/s)
};

enum boundary_style {
	CONSTANT,	// constant
	PERIODIC	// wraps around
};

/// runs a Forward-Time Central-Space discritized simulation of heat diffusion of object o
/*
 *
 * @param dt Length of time step (s)
 * @param steps Number of time steps to simulate
 */
void ftcs1d(object1d& o, double dt, int steps, std::function<void(const object1d&)> callback) {
	// generated parameters
	double dx = o.x_len()/o.x_steps();			// distance between grid points (m)
	double C = o.alpha()*dt/(pow(dx,2));		// C! (unitless)
	size_t xsize = o.x_steps()+2;

	// Set up two buffers to store the current and previous timestep information. These will
	// be swapped rather than copied. Each dimension of the buffer will be the size of the 
	// object plus two for the boundary conditions;
	std::vector<object1d> buf;
	buf.push_back(object1d(o.x_len(),xsize,o.alpha()));
	buf.push_back(object1d(o.x_len(),xsize,o.alpha()));
	
	boundary_style	bstyle = CONSTANT;
	double			bvalue = 0.0;

	for (int i = 0; i < xsize; i++) {
		// check for boundary conditions
		if (i == 0) {
			if (bstyle == CONSTANT) {
				buf[0][0] = bvalue;
			} else if (bstyle == PERIODIC) {
				buf[0][0] = o[o.x_steps()-1];
			}
		} else if (i == xsize-1) {
			if (bstyle == CONSTANT) {
				buf[0][xsize-1] = bvalue;
			} else if (bstyle == PERIODIC) {
				buf[0][xsize-1] = o[0];
			}
		} else {
			buf[0][i] = o[i-1];
		}
	}

	std::cout << "alpha: " << o.alpha() << " dt: " << dt << " dx " << dx << std::endl;
	
	callback(buf[0]);

	for (int t = 0; t < steps; t++) {
		for (size_t x = 1; x < xsize-1; x++) {
			buf[(t+1)%2][x] = buf[t%2][x] + C*(buf[t%2][x+1]-2*buf[t%2][x]+buf[t%2][x-1]);
		}
		
		// simulate
		callback(buf[(t+1)%2]);
	}
}

void print(const object1d& o) {
	std::cout << o.str() << std::endl;
}

int main() {
	object1d foo(1.0,10,0.001);
	
	foo.init();
	
	ftcs1d(foo,0.0005,20,&print);
	
	return 0;
}
