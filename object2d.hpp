#ifndef HEAT_OBJECT2D_HPP
#define HEAT_OBJECT2D_HPP

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class object2d {
public:
	typedef std::vector<std::vector<double>> data_type;

	object2d(double lx, double ly, size_t nx, size_t ny, double alpha)
	 : m_data(nx),
	   m_lx(lx),m_ly(ly),
	   m_nx(nx),m_ny(ny),
	   m_alpha(alpha)
	{
		for (size_t x = 0; x < m_nx; x++) {
			m_data[x].resize(m_ny);
		}
	}
	
	/// Initialize object using an initialization strategy
	/* Valid strategies:
	 *   FLAT: all cells initialized to `value`
	 *   GAUSSIAN: smooth gaussian distribution
	 *   GAUSSIAN_NOISE: gaussian distribution with noise
	 *
	 * @param i Initialization strategy to use
	 * @param value Value used by some initialization strategies. Default is 0.0
	 */
	void init(initial_condition i, double value = 0.0) {
		switch (i) {
			case FLAT:
				for (size_t x = 0; x < m_nx; x++) {
					for (size_t y = 0; y < m_ny; y++) {
						m_data[x][y] = value;
					}
				}
				break;
			case GAUSSIAN:
				set_gaussian();	
				break;
			case GAUSSIAN_NOISE:
				set_gaussian();	
				// TODO: noise
				break;
			default:
				throw std::invalid_argument("Invalid initial condition");
		}
	}
	
	void set_gaussian() {
		auto gx = gaussian(m_nx);
		auto gy = gaussian(m_ny);
		
		for (size_t x = 0; x < m_nx; x++) {
			for (size_t y = 0; y < m_ny; y++) {
				m_data[x][y] = gx[x]*gy[y];
			}
		}
	}

	double lx() const {return m_lx;}
	double ly() const {return m_ly;}

	size_t nx() const {return m_nx;}
	size_t ny() const {return m_ny;}

	double alpha() const {return m_alpha;}
	
	std::vector<double>& operator[] (const size_t index) {return m_data[index];}
	const std::vector<double>& operator[] (const size_t index) const {return m_data[index];}
	
	/// output data as packed binary array
	/* Output m_data as a packed array of 1 byte integers suitable for graphing
	 * or drawing. Output is in row major order.
	 *
	 * @param edges should the edges of the data set be included? Default is yes
	 * @return binary array as std::string
	 */
	std::string binary(bool edges = true) const {
		std::string data;
		
		size_t sx = (edges ? m_nx : m_nx-2);
		size_t sy = (edges ? m_ny : m_ny-2);
		size_t offset = (edges ? 0 : 1);
		
		data.resize(sx*sy);
		
		for (size_t x = 0; x < sx; x++) {
			for (size_t y = 0; y < sy; y++) {
				data[x*sx+y] = uint8_t(m_data[x+offset][y+offset]*255);
			}
		}
		
		return data;
	}
	
	/// output data as JSON array
	/* Output m_data as a JSON array of arrays.
	 *
	 * @param edges should the edges of the data set be included? Default is yes
	 * @return json array as std::string
	 */
	std::string json(bool edges = true) const {
		std::stringstream data;
		
		size_t sx = (edges ? m_nx : m_nx-2);
		size_t sy = (edges ? m_ny : m_ny-2);
		size_t offset = (edges ? 0 : 1);
		
		data << "[";
		
		std::string xsep;
		
		for (size_t x = 0; x < sx; x++) {
			for (size_t y = 0; y < sy; y++) {
				data << xsep << int(m_data[x+offset][y+offset]*255);
				xsep = ",";
			}
		}
		
		data << "]";
		
		return data.str();
	}
	
	/// runs a FTCS discritized simulation of heat diffusion of object o
	/* Run a Forward-Time Central-Space simulation of heat diffusion of object o
	 *
	 * @param ts Number of timesteps to simulate
	 * @param dt Duration of each time step (s)
	 * @param b Boundary handling style. Options are:
	 *              CONSTANT: boundary cells are set to `v` and never changed
	 *              PERIODIC: boundary cells wrap around the object to take on
	 *                        the temperature value on the opposite side.
	 * @param v Value to use for CONSTANT boundary handling style
	 * @param S Time independent source term to be added at each timestep
	 * @param callback function to call periodically during the simulation to
	 *                 provide feedback to the caller and test whether to halt
	 *                 the simulation early.
	 * @param callback_interval Number of timesteps between each callback
	 */
	void ftcs(size_t ts,
			  double dt,
			  boundary_style b, 
			  double v,
			  const object2d& S,
			  std::function<bool(const object2d&,size_t ts)> callback,
			  size_t callback_interval) const
	{
		// generated parameters
		double dx = m_lx / m_nx;			// distance between grid points (m)
		double dy = m_ly / m_ny;			// distance between grid points (m)
		double C = m_alpha*dt/(pow(dx,2));	// C! (unitless)
		size_t nx = m_nx+2;
		size_t ny = m_ny+2;
				
		// Set up two buffers to store the current and previous timestep information. These will
		// be swapped rather than copied. Each dimension of the buffer will be the size of the 
		// object plus two for the boundary conditions;
		std::vector<object2d> buf; // optimization: buf should only be cleared if necessary
		buf.push_back(object2d(m_lx,m_ly,nx,ny,m_alpha));
		buf.push_back(object2d(m_lx,m_ly,nx,ny,m_alpha));
		
		// fill the buffer array with data from o
		for (int x = 1; x < nx-1; x++) {
			for (int y = 1; y < ny-1; y++) {
				buf[0][x][y] = m_data[x-1][y-1];
			}
		}
		
		// compute boundary cell
		if (b == CONSTANT) {
			compute_constant_boundaries(buf[0],v);
		} else {
			compute_periodic_boundaries(buf[0]);
		}
		
		// run simulation
		size_t t;
		for (t = 0; t < ts; t++) {
			if (t%callback_interval == 0) {
				if (!callback(buf[(t+1)%2],t)) {
					break;
				}
			}
			
			// simulate non-boundary squares
			for (size_t x = 1; x < nx-1; x++) {
				for (size_t y = 1; y < ny-1; y++) {
					buf[(t+1)%2][x][y] = buf[t%2][x][y] + S[x-1][y-1] + 
										C*(
											buf[t%2][x][y-1] +
											buf[t%2][x][y+1] +
											buf[t%2][x-1][y] +
											buf[t%2][x+1][y] -
											4*buf[t%2][x][y]);
				}
			}
			
			// boundary conditions
			// The above simulation loop doesn't change the edges. If they were
			// constant this is correct, if they were periodic they need to be 
			// re-filled.
			if (b == PERIODIC) {
				compute_periodic_boundaries(buf[(t+1)%2]);
			}
		}
		callback(buf[(t+1)%2],t);
	}
	
	/// Sets boundary cells of an object to a constant value.
	/* 
	 * @param o Object to write to
	 * @param val Value to write
	 */
	void compute_constant_boundaries(object2d& o, double val) const {
    	for (size_t x = 0; x < o.nx(); x++) {
			o[x][0] = val;
			o[x][o.ny()-1] = val;
		}
		for (size_t y = 1; y < o.ny()-1; y++) {
			o[0][y] = val;
			o[o.nx()-1][y] = val;
		}
    }
    
    /// Sets boundary cells of an object to wrap around.
	/* 
	 * @param o Object to write to
	 */
    void compute_periodic_boundaries(object2d& o) const {
    	// corner cases. litterally. not used in this simulation, but could
		// be useful with approximation algorithms that used diagonals.
		//o[0][0] = o[1][1];
		//o[o.nx()-1][0] = o[o.nx()-2][1];
		//o[0][ny-1] = o[1][o.ny()-2];
		//o[nx-1][ny-1] = o[o.nx()-2][o.ny()-2];
		for (size_t x = 1; x < o.nx()-1; x++) {
			o[x][0] = o[x][o.ny()-2];
			o[x][o.ny()-1] = o[x][1];
		}
		for (size_t y = 1; y < o.ny()-1; y++) {
			o[0][y] = o[o.nx()-2][y];
			o[o.nx()-1][y] = o[1][y-1];
		}
    }
	
private:
	data_type m_data;				// vector representing evenly spaced grid points
	double m_lx;					// length of x dimension (m)
	double m_ly;					// length of y dimension (m)
	size_t m_nx;					// grid spaces in x dimension
	size_t m_ny;					// grid spaces in y dimension
	double m_alpha;					// thermal diffusivity (m^2/s)
};

#endif // HEAT_OBJECT2D_HPP
