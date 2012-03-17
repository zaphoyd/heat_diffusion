#ifndef HEAT_OBJECT3D_HPP
#define HEAT_OBJECT3D_HPP

#include <cassert>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

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
    
	object3d(size_t nx, size_t ny, size_t nz, double alpha)
	 : m_data(nx),
	   m_lx(1.0),m_ly(1.0),m_lz(1.0),
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
	
	object3d() {}

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
			case SPIKE:
				for (size_t x = 0; x < m_nx; x++) {                 
                    for (size_t y = 0; y < m_ny; y++) {                     
                        for (size_t z = 0; z < m_nz; z++) {
                            m_data[x][y][z] = 0.0;
                        }
                    }
                }
				m_data[m_nx/2][m_ny/2][m_nz/2] = value;
				break;
            default:
                throw std::invalid_argument("Invalid initial condition");
        }
    }
    
	object3d& operator+=(const object3d& rhs) {
		for (size_t x = 0; x < m_nx; x++) {
            for (size_t y = 0; y < m_ny; y++) {
                for (size_t z = 0; z < m_nz; z++) {
                    m_data[x][y][z] += rhs[x][y][z];
                }
            }
        }
		return *this;
	}

    void set_gaussian() {
        std::vector<double> gx = gaussian(m_nx);
        std::vector<double> gy = gaussian(m_ny);
        std::vector<double> gz = gaussian(m_nz);
        
        for (size_t x = 0; x < m_nx; x++) {
            for (size_t y = 0; y < m_ny; y++) {
                for (size_t z = 0; z < m_nz; z++) {
                    m_data[x][y][z] = gx[x]*gy[y]*gz[z];
                    //m_data[x][y][z] = 0.5;
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
	
		double max = 1.0;
        for (size_t x = 0; x < sx; x++) {
            for (size_t y = 0; y < sy; y++) {
                for (size_t z = 0; z < sz; z++) {
					if (m_data[x][y][z] > max) {
						max = m_data[x][y][z];
					}
                }
            }
        }
		
		//std::cout << "max: " << max << std::endl;

		for (size_t x = 0; x < sx; x++) {
            for (size_t y = 0; y < sy; y++) {
                for (size_t z = 0; z < sz; z++) {
                    data[x*sx+y*sy+z] = uint8_t((m_data[x+offset][y+offset][z+offset]/max)*255);
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
        
		double max = 0.0;
        for (size_t x = 0; x < sx; x++) {
            for (size_t y = 0; y < sy; y++) {
				if (m_data[x+offset][y+offset][sz] > max) {
					max = m_data[x+offset][y+offset][sz];
				}
            }
        }
		max = 1.0;
		//std::cout << "max: " << max << std::endl;
 
        for (size_t x = 0; x < sx; x++) {
            for (size_t y = 0; y < sy; y++) {
                data[x*sx+y] = int((m_data[x+offset][y+offset][sz]/max)*255);
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
			
			/*if (nx > 16) {
				std::cout << "at: " << t << ": " 
						  << xcur[nx/2+1][nx/2+1][nx/2+1] << ", " 
						  << xcur[nx/2+2][nx/2+2][nx/2+2] << ", " 
						  << xcur[nx/2+11][nx/2+11][nx/2+11] << std::endl;
			}*/
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
	
	/* Full Multigrid Algorithm for solution of the steady state heat
	 * equation with forcing.  On input u[1..n][1..n] contains the
	 * right-hand side ρ, while on output it returns the solution.  The
	 * dimension n must be of the form 2j + 1 for some integer j. (j is
	 * actually the number of grid levels used in the solution, called ng
	 * below.) ncycle is the number of V-cycles to be used at each level.
	 */
	void multigrid(relax_scheme scheme,
			       //double w
				   int ncycle,
			       size_t ts,
                   double dt,
                   //boundary_style bs, 
                   //double v,
                   const object3d& S,
                   std::function<bool(const object3d&,size_t ts)> callback,
                   size_t callback_interval) {
		object3d solution[NGMAX+1]; // solution at each grid level
		object3d rhs[NGMAX+1]; // rhs at each grid level
		object3d res[NGMAX+1]; // residual at each grid level
		object3d scratch[NGMAX+1]; // rhs during intial solution of FMG
		object3d source[NGMAX+1];   // time dependent source term
		
		size_t n = nx(); 		// Size of object
		
		// grid levels including the solution level
		// the coarsest grid is level 0
		// the finest (ie input grid) is level (ng-1)
		size_t ng = std::log2(n); 
		
		// input checks
		if (n != ny() || n != nz()) {
			throw std::invalid_argument("mglin currently only supports cubic domains");
		}
		if (n != 1+(1L << ng)) {
			throw std::invalid_argument("n-1 must be a power of 2 in mglin.");
		}
		if (ng > NGMAX) {
			throw std::invalid_argument("Object size exceeds maximum multi-grid depth. Increase NGMAX in mglin.");
		}
		
		size_t nn = n; // nn represents the size of our working grid.
		
		// the initial conditions start in solution and source terms in source
		solution[ng-1] = *this;
		source[ng-1] = S;
		scratch[ng-1] = object3d(nn,nn,nn,m_alpha);
		
		// copy restricted versions of the input down the level stack
		for (int i = ng-2; i >= 0; i--) {
			nn = nn/2+1;
			solution[i] = object3d(nn,nn,nn,m_alpha);
			source[i] = object3d(nn,nn,nn,m_alpha);
			scratch[i] = object3d(nn,nn,nn,m_alpha);
			rstrct(solution[i+1],solution[i]);
			rstrct(source[i+1],source[i]);

			//callback(source[i+1],0);
			//sleep(1);
		}
		
		size_t t;
		
		//callback(rho[ng-4],0);
		
		//solve_small(rho[0],solution[0]);
		
		/*nn = n;
		rhs[ng-1] = solution[ng-1];
		for (t = 0; t < ts; t++) {
			for (int i = 0; i < 25; i++) {
				relax(dt,solution[ng-1],rhs[ng-1],solution[ng-1]);
			}
			
			for (size_t x = 0; x < nn; x++) {
				for (size_t y = 0; y < nn; y++) {
					for (size_t z = 0; z < nn; z++) {
						solution[ng-1][x][y][z] += source[ng-1][x][y][z]*dt;
					}
				}
			}

			rhs[ng-1] = solution[ng-1];
			std::cout << solution[ng-1][64][64][64] << ", " << solution[ng-1][65][65][65] << std::endl;
		}

		callback(solution[ng-1],t);
		return;*/

		// allocate memory for the intermediate levels
		nn = 3;
		//solution[0] = object3d(nn,nn,nn,m_alpha);
		rhs[0] = object3d(nn,nn,nn,m_alpha);

		for (size_t level = 1; level < ng; level++) {
			nn = 2*nn-1;
			
			rhs[level] = object3d(nn,nn,nn,m_alpha);
			res[level] = object3d(nn,nn,nn,m_alpha);

		}
		
		for (t = 0; t < ts; t++) {
			if (t%callback_interval == 0) {
                if (!callback(solution[ng-1],t)) {
                    break;
                }
            }
            
            // At each timestep we will run the multigrid algorithm with the 
            // output of the previous timestep as the rhs. This involves
            
            
            // Set up data
            
			// copy the solution from the previous timestep into rhs
			rhs[0] = solution[0];
	
			// solve level 0 directly
			solve_small(dt,solution[0],rhs[0]);
	
			// solve levels 1 through ng-1
			for (size_t level = 1; level < ng; level++) {				
				// copy the solution from the previous timestep to rhs
				rhs[level] = solution[level];
	
				// copy lhs from previous level of this timestep by interpolating.
				trilinear_interpolate(solution[level-1],solution[level]);
				
				// perform `ncycle` vcycles starting at the current level
				for (int jcycle = 0; jcycle < ncycle; jcycle++) {
					// start at current level and go down until the second coarsest
					for (size_t cycle_level = level; cycle_level > 0; cycle_level--) {
						for (int j = 0; j < NPRE; j++) {
							// solution sweeps
							if (scheme == RS_JACOBI) {
								relax(dt, solution[cycle_level], rhs[cycle_level], scratch[cycle_level]);
								solution[cycle_level] = scratch[cycle_level];
							} else if (scheme == RS_GAUSS_SEIDEL) {
								// identical to JACOBI except solution is passed
								// as both the input and output parameter
								relax(dt, solution[cycle_level], rhs[cycle_level], solution[cycle_level]);
							} else if (scheme == RS_RB_GAUSS_SEIDEL) {
								relax_red_black(dt, solution[cycle_level], rhs[cycle_level], solution[cycle_level]);
							}
						}
						
						// compute residual, store in res
						residual(dt,solution[cycle_level],rhs[cycle_level],res[cycle_level]);
											
						// restrict residuals as rhs of next coarsest scale
						rstrct(res[cycle_level],rhs[cycle_level-1]);
	
						solution[cycle_level-1].init(FLAT);
					}
					
					// solve coarsest directly
					solve_small(dt,solution[0],rhs[0]);
					
					// cycle back up to current level
					for (size_t cycle_level = 1; cycle_level < level; cycle_level++) {
						// interpolate error and add to previous solution guess
						trilinear_interpolate(solution[cycle_level-1],res[cycle_level]);
						
						solution[cycle_level] += res[cycle_level];
						
						// do NPOST sweeps
						for (int j = 0; j < NPOST; j++) {
							if (scheme == RS_JACOBI) {
								relax(dt, solution[cycle_level], rhs[cycle_level], scratch[cycle_level]);
								solution[cycle_level] = scratch[cycle_level];
							} else if (scheme == RS_GAUSS_SEIDEL) {
								// identical to JACOBI except solution is passed
								// as both the input and output parameter
								relax(dt, solution[cycle_level], rhs[cycle_level], solution[cycle_level]);
							} else if (scheme == RS_RB_GAUSS_SEIDEL) {
								relax_red_black(dt, solution[cycle_level], rhs[cycle_level], solution[cycle_level]);
							}
						}
					}
				}
	
			}
			// at this point the solution after this timestep will be stored in 
			// solution[ng-1]
			
			// Add time dependent source term
			for (size_t x = 0; x < nn; x++) {
				for (size_t y = 0; y < nn; y++) {
					for (size_t z = 0; z < nn; z++) {
						solution[ng-1][x][y][z] += source[ng-1][x][y][z]*dt;
					}
				}
			}
			
			/*if (n > 16) {
				std::cout << "at: " << t << ": " 
			    	      << solution[ng-1][n/2][n/2][n/2] << ", " 
			        	  << solution[ng-1][n/2+1][n/2+1][n/2+1] << ", " 
			        	  << solution[ng-1][n/2+10][n/2+10][n/2+10] << std::endl;
			}*/
		}
		
		// callback with solution at final timestep
		callback(solution[ng-1],t);
	}
	
	// calculate minus the residual and store in res
	void residual(double dt,object3d& u, object3d& rhs, object3d& res) {
		assert(u.nx() == rhs.nx());
		
		size_t n = u.nx();
		double dx = 1.0/(n-1);
		double dx2i = 1.0/(dx*dx);
		
		double C = m_alpha*dt/(dx*dx);

		// Interior Points
		for (size_t x = 1; x < n-1; x++) {
			for (size_t y = 1; y < n-1; y++) {
				for (size_t z = 1; z < n-1; z++) {
					res[x][y][z] = -1 * C*(u[x+1][y][z]+u[x-1][y][z]
					                    +u[x][y+1][z]+u[x][y-1][z]
					                    +u[x][y][z+1]+u[x][y][z-1]
					                    -6.0*u[x][y][z])
					               +rhs[x][y][z];
				}	
			}	
		}
		
		// Boundary Points
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				res[i][j][0] = res[i][j][n-1] = 0.0;
				res[i][0][j] = res[i][n-1][j] = 0.0;
				res[0][i][j] = res[n-1][i][j] = 0.0;
			}
		}
	}
	
	// Half-weighting restriction. source represents the fine-grid solution input
	// the course-grid solution is stored in target.
	void rstrct(object3d& fine, object3d& coarse) {
		size_t nfine = fine.nx();
		size_t ncoarse = coarse.nx();

		assert(nfine == 2*ncoarse-1);
		
		size_t xc,yc,zc; // coarse counters
		size_t xf,yf,zf; // fine counters
		double s = 1/12;  // weight of non center points
		
		// Interior Points
		for (xc = 1,xf = 2; xc < ncoarse-1; xc++,xf+=2) {
			for (yc = 1,yf = 2; yc < ncoarse-1; yc++,yf+=2) {
				for (zc = 1,zf = 2; zc < ncoarse-1; zc++,zf+=2) {
					coarse[xc][yc][zc] = 0.5*fine[xf][yf][zf]+
										 s*(fine[xf+1][yf][zf]+
										    fine[xf-1][yf][zf]+
										    fine[xf][yf+1][zf]+
										    fine[xf][yf-1][zf]+
										    fine[xf][yf][zf+1]+
										    fine[xf][yf][zf-1]);
				}
			}
		}
		
		// Boundary Points
		// x/y faces
		for (xc = 0,xf = 0; xc < ncoarse; xc++,xf+=2) {
			for (yc = 0,yf = 0; yc < ncoarse; yc++,yf+=2) {
				coarse[xc][yc][0] = fine[xf][yf][0];
				coarse[xc][yc][ncoarse-1] = fine[xf][yf][nfine-1];
			}
		}

		// x/z faces
		for (xc = 0,xf = 0; xc < ncoarse; xc++,xf+=2) {
			for (zc = 0,zf = 0; zc < ncoarse; zc++,zf+=2) {
				coarse[xc][0][zc] = fine[xf][0][zf];
				coarse[xc][ncoarse-1][zc] = fine[xf][nfine-1][zf];
			}
		}
		
		// y/z faces
		for (yc = 0,yf = 0; yc < ncoarse; yc++,yf+=2) {
			for (zc = 0,zf = 0; zc < ncoarse; zc++,zf+=2) {
				coarse[0][yc][zc] = fine[0][yf][zf];
				coarse[ncoarse-1][yc][zc] = fine[nfine-1][yf][zf];
			}
		}
		
		// x/y plane columns
		for (yc = 0,yf = 0; yc < ncoarse; yc++,yf+=2) {
			for (zc = 0,zf = 0; zc < ncoarse; zc++,zf+=2) {
				coarse[0][yc][zc] = fine[0][yf][zf];
				coarse[ncoarse-1][yc][zc] = fine[nfine-1][yf][zf];
			}
		}
	}
	
	// Fills in `fine` based on information interpolated from coarse
	void trilinear_interpolate(object3d& coarse, object3d& fine) {
		size_t nfine = fine.nx();
		size_t ncoarse = coarse.nx();

		assert(nfine == 2*ncoarse-1);
		
		size_t xc,yc,zc; // coarse counters
		size_t xf,yf,zf; // fine counters
		
		// exact values
		for (xc = 0,xf = 0; xc < ncoarse; xc++,xf+=2) {
			for (yc = 0,yf = 0; yc < ncoarse; yc++,yf+=2) {
				for (zc = 0,zf = 0; zc < ncoarse; zc++,zf+=2) {
					fine[xf][yf][zf] = coarse[xc][yc][zc];
				}
			}
		}
		
		// interpolate along x axis
		for (xf = 1; xf < nfine; xf+=2) {
			for (yf = 0; yf < nfine; yf+=2) {
				for (zf = 0; zf < nfine; zf+=2) {
					fine[xf][yf][zf] = 0.5*(fine[xf+1][yf][zf]+fine[xf-1][yf][zf]);
				}
			}
		}
		
		// interpolate along y axis
		for (xf = 0; xf < nfine; xf++) {
			for (yf = 1; yf < nfine; yf+=2) {
				for (zf = 0; zf < nfine; zf+=2) {
					fine[xf][yf][zf] = 0.5*(fine[xf][yf+1][zf]+fine[xf][yf-1][zf]);
				}
			}
		}
		
		// interpolate along z axis
		for (xf = 0; xf < nfine; xf++) {
			for (yf = 0; yf < nfine; yf++) {
				for (zf = 1; zf < nfine; zf+=2) {
					fine[xf][yf][zf] = 0.5*(fine[xf][yf][zf+1]+fine[xf][yf][zf-1]);
				}
			}
		}
	}
	
	// performs a single step of the Jacobi, Gauss-Seidel
	// if lhs != output then algorithm is Jacobi
	// if lhs == output then algorithm is GS
	void relax(double dt, object3d& lhs, object3d& rhs, object3d& output) {
		assert(output.nx() == lhs.nx());
		assert(output.nx() == rhs.nx());
		
		size_t n = output.nx();
		double dx = 1.0/n;

		double C = m_alpha*dt/(dx*dx);
		double C2 = C/(6*C+1.0);
		double C3 = 1.0/(6*C+1.0);
		
		// plain Gauss-Seidel / Jacobi
		
		for (size_t x = 1; x < n-1; x++) {
			for (size_t y = 1; y < n-1; y++) {
				for (size_t z = 1; z < n-1; z++) {
					output[x][y][z] = C2*(lhs[x+1][y][z]
										  +lhs[x-1][y][z]
										  +lhs[x][y+1][z]
										  +lhs[x][y-1][z]
										  +lhs[x][y][z+1]
										  +lhs[x][y][z-1])
									 +C3*rhs[x][y][z];
				}	
			}	
		}
	}

	// performs a single step of the Red/Black G-S alg
	void relax_red_black(double dt, object3d& lhs, object3d& rhs, object3d& output) {
		assert(output.nx() == lhs.nx());
		assert(output.nx() == rhs.nx());
		
		size_t n = output.nx();
		double dx = 1.0/n;

		double C = m_alpha*dt/(dx*dx);
		double C2 = C/(6*C+1.0);
		double C3 = 1.0/(6*C+1.0);
		
		// Red/Black Gauss-Seidel
		size_t offset = 1;
		for (int pass = 0; pass < 2; pass++) {
			for (size_t x = 1; x < n-1; x++,offset=3-offset) {
				for (size_t y = 1; y < n-1; y++,offset=3-offset) {
					for (size_t z = offset; z < n-1; z+=2) {
						output[x][y][z] = C2*(lhs[x+1][y][z]
											  +lhs[x-1][y][z]
											  +lhs[x][y+1][z]
											  +lhs[x][y-1][z]
											  +lhs[x][y][z+1]
											  +lhs[x][y][z-1])
										 +C3*rhs[x][y][z];
					}	
				}	
			}
		}
	}

	void solve_small(double dt, object3d& solution, object3d& rhs) {
		assert(solution.nx() == 3);
		assert(rhs.nx() == 3);
		
		size_t n = rhs.nx();
		double dx = 1.0/n;

		double C = m_alpha*dt/(dx*dx);
		double C3 = 1.0/(6*C+1.0);

		solution[1][1][1] = C3*rhs[1][1][1];
	}
private:
	static const unsigned int NGMAX = 15; // Maximum grid levels
	static const unsigned int NPRE = 2; // number of GS sweeps before vcycles
	static const unsigned int NPOST = 2; // number of GS sweeps after vcycles
	
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
