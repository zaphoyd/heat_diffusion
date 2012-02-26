#ifndef HEAT_OBJECT3D_HPP
#define HEAT_OBJECT3D_HPP

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class object3d {
public:
    typedef std::vector<std::vector<std::vector<double>>> data_type;

    object3d(double lx, double ly, double lz, size_t nx, size_t ny, size_t nz, double alpha)
     : m_data(nx),
       m_lx(lx),m_ly(ly),m_lz(lz),
       m_nx(nx),m_ny(ny),m_nz(nz),
       m_alpha(alpha)
    {
        for (size_t x = 0; x < m_nx; x++) {
            m_data[x].resize(m_ny);
            
            for (size_t y = 0; y < m_ny; y++) {
                m_data[x][y].resize(m_nz);
            }
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
                        for (size_t z = 0; z < m_nz; z++) {
                            m_data[x][y][z] = value;
                        }
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
        std::vector<double> gx = gaussian(m_nx);
        std::vector<double> gy = gaussian(m_ny);
        std::vector<double> gz = gaussian(m_nz);
        
        for (size_t x = 0; x < m_nx; x++) {
            for (size_t y = 0; y < m_ny; y++) {
                for (size_t z = 0; z < m_nz; z++) {
                    m_data[x][y][z] = gx[x]*gy[y]*gz[z];
                }
            }
        }
    }

    double lx() const {return m_lx;}
    double ly() const {return m_ly;}
    double lz() const {return m_lz;}

    size_t nx() const {return m_nx;}
    size_t ny() const {return m_ny;}
    size_t nz() const {return m_nz;}

    double alpha() const {return m_alpha;}
    
    std::vector<std::vector<double>>& operator[] (const size_t index) {
        return m_data[index];
    }
    
    const std::vector<std::vector<double>>& operator[] (const size_t index) const {
        return m_data[index];
    }
    
    /// output data as packed binary array
    /* Output m_data as a packed array of 1 byte integers suitable for graphing
     * or drawing. Output is row-column-depth order (Row major)
     *
     * @param edges should the edges of the data set be included? Default is yes
     * @return binary array as std::string
     */
    std::string binary(bool edges = true) const {
        std::string data;
        
        size_t sx = (edges ? m_nx : m_nx-2);
        size_t sy = (edges ? m_ny : m_ny-2);
        size_t sz = (edges ? m_nz : m_nz-2);
        size_t offset = (edges ? 0 : 1);
        
        data.resize(sx*sy*sz);
        
        for (size_t x = 0; x < sx; x++) {
            for (size_t y = 0; y < sy; y++) {
                for (size_t z = 0; z < sz; z++) {
                    data[x*sx+y*sy+z] = uint8_t(m_data[x+offset][y+offset][z+offset]*255);
                }
            }
        }
        
        return data;
    }
    
    /// output a 2D slice of data as packed binary array
    /* Output a 2D slice of m_data as a packed array of 1 byte integers suitable
     * for graphing or drawing.
     *
     * @param z z-index to slice on
     * @param edges should the edges of the data set be included? Default is yes
     * @return binary array as std::string
     */
    std::string binary_slice(size_t z, bool edges = true) const {
        std::string data;
        
        size_t sx = (edges ? m_nx : m_nx-2);
        size_t sy = (edges ? m_ny : m_ny-2);
        size_t sz = (edges ? z : z-1);
        size_t offset = (edges ? 0 : 1);
        
        data.resize(sx*sy);
        
        for (size_t x = 0; x < sx; x++) {
            for (size_t y = 0; y < sy; y++) {
                data[x*sx+y] = int(m_data[x+offset][y+offset][sz]*255);
            }
        }
        
        return data;
    }
    
    /// output a 2D slice of data as a JSON array
    /* Output a 2D slice of m_data as a JSON array of arrays of numbers
     *
     * @param z z-index to slice on
     * @param edges should the edges of the data set be included? Default is yes
     * @return JSON array as std::string
     */
    std::string json_slice(size_t z, bool edges = true) const {
        std::stringstream data;
        
        size_t sx = (edges ? m_nx : m_nx-2);
        size_t sy = (edges ? m_ny : m_ny-2);
        size_t sz = (edges ? m_nz : m_nz-2);
        size_t offset = (edges ? 0 : 1);
        
        data << "[";
        
        std::string xsep;
        std::string ysep;
        
        for (size_t x = 0; x < sx; x++) {
            data << xsep << "[";
            xsep = ",";
            ysep = "";
            for (size_t y = 0; y < sy; y++) {
                data << ysep 
                     << uint8_t(m_data[x+offset][y+offset][z+offset]*255);
                ysep = ",";
            }
            data << "]";
        }
        
        data << "]";
        
        return data.str();
    }

    // TODO ftcs 3d

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
              const object3d& S,
              std::function<bool(const object3d&,size_t ts)> callback,
              size_t callback_interval) const
    {
        // generated parameters
        double dx = m_lx / m_nx;            // distance between grid points (m)
        double dy = m_ly / m_ny;            // distance between grid points (m)
        double dz = m_lz / m_nz;            // distance between grid points (m)
        double C = m_alpha*dt/(pow(dx,2));  // C! (unitless)
        size_t nx = m_nx+2;
        size_t ny = m_ny+2;
        size_t nz = m_nz+2;
                
        // Set up two buffers to store the current and previous timestep information. These will
        // be swapped rather than copied. Each dimension of the buffer will be the size of the 
        // object plus two for the boundary conditions;
        std::vector<object3d> buf; // optimization: buf should only be cleared if necessary
        buf.push_back(object3d(m_lx,m_ly,m_lz,nx,ny,nz,m_alpha));
        buf.push_back(object3d(m_lx,m_ly,m_lz,nx,ny,nz,m_alpha));
        
        // fill the buffer array with data from o
        for (int x = 1; x < nx-1; x++) {
            for (int y = 1; y < ny-1; y++) {
                for (int z = 1; z < nz-1; z++) {
                    buf[0][x][y][z] = m_data[x-1][y-1][z-1];
                }
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
                    for (size_t z = 1; z < nz-1; z++) {
                        buf[(t+1)%2][x][y][z] = buf[t%2][x][y][z] + dt*S[x-1][y-1][z-1] + 
                                        C*(
                                            buf[t%2][x][y-1][z] +
                                            buf[t%2][x][y+1][z] +
                                            buf[t%2][x-1][y][z] +
                                            buf[t%2][x+1][y][z] +
                                            buf[t%2][x][y][z-1] +
                                            buf[t%2][x][y][z+1] -
                                            6*buf[t%2][x][y][z]);
                    }
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
   
		/// runs a simulation of heat diffusion of object o using an iterative method
    /* Run a simulation of heat diffusion of object o using the specified iterative
	 * method. 
     *
	 * @param method method to use for performing the iterative calculation. Options:
	 *               JACOBI: Jacobi Iteration
	 *               GAUSS_SEIDEL: Gauss Seidel iteration
	 *               SOR: Successive Over-relaxation
	 * @param w Relaxation factor for use when method is SOR
     * @param ts Number of timesteps to simulate
     * @param dt Duration of each time step (s)
     * @param bs Boundary handling style. Options are:
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
    void iterative(disc_method method,
			    double w,
			    size_t ts,
                double dt,
                boundary_style bs, 
                double v,
                const object3d& S,
                std::function<bool(const object3d&,size_t ts)> callback,
                size_t callback_interval) const
    {
        double dx = m_lx / m_nx;
		double dy = m_ly / m_ny;
		double dz = m_lz / m_nz;
        double C = m_alpha*dt/(pow(dx,2));
        size_t nx = m_nx+2;
		size_t ny = m_ny+2;
		size_t nz = m_nz+2;
		
		object3d xold(m_lx,m_ly,m_lz,nx,ny,nz,m_alpha);
		object3d xcur(m_lx,m_ly,m_lz,nx,ny,nz,m_alpha);
		object3d xnew(m_lx,m_ly,m_lz,nx,ny,nz,m_alpha);
		
		// fill the buffer array with data from o
        for (size_t x = 1; x < nx-1; x++) {
            for (size_t y = 1; y < ny-1; y++) {
				for (size_t z = 1; z < nz-1; z++) {
                	xold[x][y][z] = m_data[x-1][y-1][z-1];
				}
            }
        }
        
        // compute boundary cell
        if (bs == CONSTANT) {
            compute_constant_boundaries(xold,v);
        } else {
            compute_periodic_boundaries(xold);
        }
		
		xcur = xold;

        size_t t;
		size_t MAX_ITER = 1000;
        double EPSILON = 1e-6;

		double C2 = (method == SOR ? (w*C)/(6*C+1) : C/(6*C+1));
		double C3 = (method == SOR ? w/(6*C+1) : 1/(6*C+1));

		int iter = 0;

		for (t = 0; t < ts; t++) {
			if (t%callback_interval == 0) {
                if (!callback(xcur,t)) {
                    break;
                }
            }
			
			size_t i = 0;
			for (i = 0; i < MAX_ITER; i++) {
				if (method == JACOBI) {
					for (size_t x = 1; x < nx-1; x++) {
						for (size_t y = 1; y < ny-1; y++) {
							for (size_t z = 1; z < nz-1; z++) {
								xnew[x][y][z] = C2*(xcur[x][y-1][z] + 
												    xcur[x-1][y][z] + 
												    xcur[x+1][y][z] +
												    xcur[x][y+1][z] +
												    xcur[x][y][z-1] +
												    xcur[x][y][z+1]) +
											    C3*xold[x][y][z];
							}
						}
					}
					if (xcur.mean_abs_diff(xnew,false) < EPSILON) {
						break;
					}
					xcur = xnew;
				} else if (method == GAUSS_SEIDEL) {
					xnew = xcur;
					for (size_t x = 1; x < nx-1; x++) {
						for (size_t y = 1; y < ny-1; y++) {
							for (size_t z = 1; z < nz-1; z++) {
								xcur[x][y][z] = C2*(xcur[x][y-1][z] + 
												    xcur[x-1][y][z] + 
												    xcur[x+1][y][z] +
												    xcur[x][y+1][z] +
												    xcur[x][y][z-1] +
												    xcur[x][y][z+1]) +
											    C3*xold[x][y][z];
							}
						}
					}
					if (xcur.mean_abs_diff(xnew,false) < EPSILON) {
						break;
					}
				} else {
					xnew = xcur;
					for (size_t x = 1; x < nx-1; x++) {
						for (size_t y = 1; y < ny-1; y++) {
							for (size_t z = 1; z < nz-1; z++) {
								xcur[x][y][z] = (1-w)*xcur[x][y][z] + 
											    C2*(xcur[x][y-1][z] + 
											        xcur[x-1][y][z] + 
											        xcur[x+1][y][z] +
											        xcur[x][y+1][z] +
											        xcur[x][y][z-1] +
											        xcur[x][y][z+1]) +
										        C3*xold[x][y][z];
							}
						}
					}
					if (xcur.mean_abs_diff(xnew,false) < EPSILON) {
						break;
					}
				}
			}
			iter += i;
			
			for (size_t x = 1; x < nx-1; x++) {
				for (size_t y = 1; y < ny-1; y++) {
					for (size_t z = 1; z < nz-1; z++) {
						xcur[x][y][z] += S[x-1][y-1][z-1]*dt;
					}
				}
			}
			
			// boundary conditions
            // The above simulation loop doesn't change the edges. If they were
            // constant this is correct, if they were periodic they need to be 
            // re-filled.
            if (bs == PERIODIC) {
                compute_periodic_boundaries(xcur);
            }

			xold = xcur;
        }

		std::cout << "average iterations: " << (iter/t) << std::endl;
        callback(xcur,t);
	}

    /// Sets boundary cells of an object to a constant value.
    /* 
     * @param o Object to write to
     * @param val Value to write
     */
    void compute_constant_boundaries(object3d& o, double val) const {
        for (size_t x = 0; x < o.nx(); x++) {
            o[x][0][0] = val;
            o[x][0][o.nz()-1] = val;
            o[x][o.ny()-1][0] = val;
            o[x][o.ny()-1][o.nz()-1] = val;
        }
        for (size_t y = 1; y < o.ny()-1; y++) {
            o[0][y][0] = val;
            o[0][y][o.nz()-1] = val;
            o[o.nx()-1][y][0] = val;
            o[o.nx()-1][y][o.nz()-1] = val;
        }
        for (size_t z = 1; z < o.nz()-1; z++) {
            o[0][0][z] = val;
            o[0][o.ny()-1][z] = val;
            o[o.nx()-1][0][z] = val;
            o[o.nx()-1][o.ny()-1][z] = val;
        }
    }
    
    /// Sets boundary cells of an object to wrap around.
    /* 
     * @param o Object to write to
     */
    void compute_periodic_boundaries(object3d& o) const {
        // The idea here is to swap faces rather than edges. Corners in this case represent
        // both literal corners and edges. 
        //
        // TODO: corner/edge cases
        
        // Swap x/y faces
        for (size_t x = 1; x < o.nx()-1; x++) {
        for (size_t y = 1; y < o.ny()-1; y++) {
            o[x][y][0] = o[x][y][o.nz()-2];
            o[x][y][o.nz()-1] = o[x][y][1];
        }
        }

        // Swap x/z faces
        for (size_t x = 1; x < o.nx()-1; x++) {
        for (size_t z = 1; z < o.nz()-1; z++) {
            o[x][0][z] = o[x][o.ny()-2][z];
            o[x][o.ny()-1][z] = o[x][1][z];
        }
        }

        // Swap y/z faces
        for (size_t y = 1; y < o.ny()-1; y++) {
        for (size_t z = 1; z < o.nz()-1; z++) {
            o[0][y][z] = o[o.nx()-2][y][z];
            o[o.nx()-1][y][z] = o[1][y][z];
        }
        }
    }
	
	/// Compute the mean of the absolute value of the difference
	/*
	 *
	 * @param o Object to subtract
	 * @param edges Whether or not to include the edge rows in the calculation. Default: true
	 */
	double mean_abs_diff(object3d& o, bool edges = true) const {
		double val = 0;
		
		if (o.nx() != nx() || o.ny() != ny() || o.nz() != nz()) {
			throw std::invalid_argument("objects must be the same size to use mean_abs_diff");
		}
		
		size_t start = (edges ? 0 : 1);
		size_t xend = (edges ? m_nx : m_nx-1);
		size_t yend = (edges ? m_ny : m_ny-1);
		size_t zend = (edges ? m_nz : m_nz-1);

		for (size_t x = start; x < xend; x++) {
			for (size_t y = start; y < yend; y++) {
				for (size_t z = start; z < zend; z++) {
					val += fabs(m_data[x][y][z] - o[x][y][z]);
				}
			}
		}
		return val / (m_nx*m_ny);
	}
private:
    data_type m_data;               // vector representing evenly spaced grid points
    double m_lx;                    // length of x dimension (m)
    double m_ly;                    // length of y dimension (m)
    double m_lz;                    // length of z dimension (m)
    size_t m_nx;                    // grid spaces in x dimension
    size_t m_ny;                    // grid spaces in y dimension
    size_t m_nz;                    // grid spaces in z dimension
    double m_alpha;                 // thermal diffusivity (m^2/s)
};

#endif // HEAT_OBJECT3D_HPP
