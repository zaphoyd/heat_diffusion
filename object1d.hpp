#ifndef HEAT_OBJECT1D_HPP
#define HEAT_OBJECT1D_HPP

#include "matrix.hpp"

#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class object1d {
public:
    typedef std::vector<double> data_type;

    object1d(double lx, size_t nx, double alpha)
     : m_data(nx),
       m_lx(lx),
       m_nx(nx),
       m_alpha(alpha) {}
    
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
                for (size_t i = 0; i < m_nx; i++) {
                    m_data[i] = value;
                }
                break;
            case GAUSSIAN:
                m_data = gaussian(m_nx);
                break;
            case GAUSSIAN_NOISE:
                m_data = gaussian(m_nx);
                // TODO: noise
                break;
            default:
                throw std::invalid_argument("Invalid initial condition");
        }
    }

    double lx() const {return m_lx;}
    size_t nx() const {return m_nx;}

    double alpha() const {return m_alpha;}
    
    double& operator[] (const size_t index) {
        return m_data[index];
    }
    
    const double& operator[] (const size_t index) const {
        return m_data[index];
    }

    /// output data as packed binary array
    /* Output m_data as a packed array of 1 byte integers suitable for graphing
     * or drawing. Output is in row major order.
     *
     * @param edges should the edges of the data set be included? Default is yes
     */
    std::string binary(bool edges = true) const {
        std::string data;
        
        size_t s = (edges ? m_nx : m_nx-2);
        size_t offset = (edges ? 0 : 1);
        data.resize(s);
    
        for (size_t x = 0; x < s; x++) {
            data[x] = uint8_t(m_data[x+offset]*255);
        }
        
        return data;
    }
    
    /// output data as JSON array
    /* Output m_data as a JSON array.
     *
     * @param edges should the edges of the data set be included? Default is yes
     * @return json array as std::string
     */
    std::string json(bool edges = true) const {
        std::stringstream data;
        
        size_t sx = (edges ? m_nx : m_nx-2);
        size_t offset = (edges ? 0 : 1);
        
        data << "[";
        
        std::string xsep;
        
        for (size_t x = 0; x < sx; x++) {
            data << xsep << int(m_data[x+offset]*255);
            xsep = ",";
        }
        
        data << "]";
        
        return data.str();
    }
    
    /// runs a FTCS discritized simulation of heat diffusion of object o
    /* Run a Forward-Time Central-Space simulation of heat diffusion of object o
     *
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
    void ftcs(size_t ts,
              double dt,
              boundary_style bs, 
              double v,
              const object1d& S,
              std::function<bool(const object1d&,size_t ts)> callback,
              size_t callback_interval) const
    {
        // generated parameters
        double dx = m_lx / m_nx;            // distance between grid points (m)
        double C = m_alpha*dt/(pow(dx,2));  // C! (unitless)
        size_t nx = m_nx+2;
                
        // Set up two buffers to store the current and previous timestep information. These will
        // be swapped rather than copied. Each dimension of the buffer will be the size of the 
        // object plus two for the boundary conditions;
        std::vector<object1d> buf; // optimization: buf should only be cleared if necessary
        buf.push_back(object1d(m_lx,nx,m_alpha));
        buf.push_back(object1d(m_lx,nx,m_alpha));
        
        // fill the buffer array with data from o
        for (int x = 1; x < nx-1; x++) {
            buf[0][x] = m_data[x-1];
        }
        
        // compute boundary cell
        if (bs == CONSTANT) {
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
                buf[(t+1)%2][x] = buf[t%2][x] + dt*S[x] + 
                    C*( buf[t%2][x-1] + buf[t%2][x+1] - 2*buf[t%2][x] );
            }
            
            // boundary conditions
            // The above simulation loop doesn't change the edges. If they were
            // constant this is correct, if they were periodic they need to be 
            // re-filled.
            if (bs == PERIODIC) {
                compute_periodic_boundaries(buf[(t+1)%2]);
            }
        }
        callback(buf[(t+1)%2],t);
    }
    
    /// runs a Crank Nichsolson discritized simulation of heat diffusion of object o
    /* Run a Crank Nicholson simulation of heat diffusion of object o
     *
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
    void crank_nicholson(size_t ts,
                         double dt,
                         boundary_style bs, 
                         double v,
                         const object1d& S,
                         std::function<bool(const object1d&,size_t ts)> callback,
                         size_t callback_interval) const
    {
        double dx = m_lx / m_nx;
        double C = m_alpha*dt/(pow(dx,2));
        size_t nx = m_nx;
        
        if (nx < 2) {
            throw std::invalid_argument("object must have size at least 2");
        }

        object1d buf(m_lx,nx,m_alpha);

        matrix<double> A(m_nx,m_nx);    // Coefficient matrix
        std::vector<double> b(m_nx);    // 
        std::vector<double> x(m_nx);    // tnew
        
        // load initial conditions
        for (size_t i = 0; i < nx; i++) {
            x[i] = m_data[i];
        }

        size_t t;
        for (t = 0; t < ts; t++) {
            if (t%callback_interval == 0) {
                // this should be optimized with a object1d copy constructor from vector or a 
                // better callback that writes vectors rather than objects back to the wire.
                for (size_t i = 0; i < nx; i++) {
                    buf[i] = x[i];
                }

                if (!callback(buf,t)) {
                    break;
                }
            }
            // Fill in A
            A[0][0] = (1+2*C);
            A[0][1] = -1*C;
            //A[0][nx-1] = 1;   // in boundary condition wrap?
            A[nx-1][nx-2] = -1*C;
            A[nx-1][nx-1] = (1+2*C);
            //A[nx-1][0] = 1;   // in boundary condition wrap?
            for (size_t i = 1; i < nx-1; i++) {
                A[i][i-1] = -1*C;
                A[i][i] = (1+2*C);
                A[i][i+1] = -1*C;
            }
            

            // Fill in b
            if (bs == CONSTANT) {
                b[0] = x[0] + (C/2.0)*(v+x[1]-2*x[0]);
                b[nx-1] = x[nx-1] + (C/2)*(x[nx-2]+v-2*x[nx-1]);
            } else {
                b[0] = x[0] + (C/2.0)*(x[nx-1]+x[1]-2*x[0]);
                b[nx-1] = x[nx-1] + (C/2)*(x[nx-2]+x[0]-2*x[nx-1]);
            }
            for (size_t i = 1; i < nx-1; i++) {
                b[i] = x[i] + (C/2.0)*(x[i-1]+x[i+1]-2*x[i]);
            }
            
            std::cout << "ts: " << t << std::endl;
            std::cout << "initial A: " << std::endl << A << std::endl;
            std::cout << "initial b: " << b << std::endl;
            
            upper_triangulate(A,b);
            
            std::cout << "solved A: " << std::endl << A << std::endl;
            std::cout << "solved b: " << b << std::endl;
            
            back_sub(A,b,x);
            
            std::cout << "solved x: " << x << std::endl;

            for (size_t i = 0; i < nx; i++) {
                x[i] += S[i]*dt;
            }
        }

        // this should be optimized with a object1d copy constructor from vector or a 
        // better callback that writes vectors rather than objects back to the wire.
        for (size_t i = 0; i < nx; i++) {
            buf[i] = x[i];
        }
        callback(buf,t);
    }
    
	/// runs a simulation of heat diffusion of object o using Jacobi iteration
    /* Run a simulation of heat diffusion of object o using Jacobi iteration and
	 * backwards euler discretization.
     *
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
    void jacobi(size_t ts,
                double dt,
                boundary_style bs, 
                double v,
                const object1d& S,
                std::function<bool(const object1d&,size_t ts)> callback,
                size_t callback_interval) const
    {
        double dx = m_lx / m_nx;
        double C = m_alpha*dt/(pow(dx,2));
        size_t nx = m_nx;
		
		object1d xold(m_lx,nx,m_alpha);
		object1d xcur(m_lx,nx,m_alpha);
		object1d xnew(m_lx,nx,m_alpha);
		
		// load initial values
		xold = *this;
		xcur = xold;

        size_t t;
		size_t MAX_ITER = 1000;
        double EPSILON = 1e-6;

		double C2 = C/(2*C+1);
		double C3 = 1/(2*C+1);

		int iter = 0;

		for (t = 0; t < ts; t++) {
			if (t%callback_interval == 0) {
                if (!callback(xcur,t)) {
                    break;
                }
            }
			
			size_t i = 0;
			for (i = 0; i < MAX_ITER; i++) {
				xnew[0] = C2*((bs == CONSTANT ? v : xcur[nx-1]) + xcur[1]) + C3*xold[0];
				for (size_t x = 1; x < nx-1; x++) {
					xnew[x] = C2*(xcur[x-1] + xcur[x+1]) + C3*xold[x];
				}
				xnew[nx-1] = C2*(xcur[nx-1] + (bs == CONSTANT ? v : xcur[0])) + C3*xold[nx-1];
				
				if (xcur.mean_abs_diff(xnew) < EPSILON) {
					break;
				}
				xcur = xnew;
			}
			iter += i;
			
			for (size_t x = 0; x < nx; x++) {
				xcur[x] += S[x]*dt;
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
    void compute_constant_boundaries(object1d& o, double val) const {
        o[0] = val;
        o[o.nx()-1] = val;
    }
    
    /// Sets boundary cells of an object to wrap around.
    /* 
     * @param o Object to write to
     */
    void compute_periodic_boundaries(object1d& o) const {
        o[0] = o[o.nx()-2];
        o[o.nx()-1] = o[1];
    }

	/// Compute the mean of the absolute value of the difference
	double mean_abs_diff(object1d& o) const {
		double val = 0;
		
		if (o.nx() != nx()) {
			throw std::invalid_argument("objects must be the same size to use mean_abs_diff");
		}
		
		for (size_t x = 0; x < m_nx; x++) {
			val += fabs(m_data[x] - o[x]);
		}
		return val / m_nx;
	}
private:
    data_type m_data;               // data vector
    double m_lx;                    // length of x dimension (m)
    size_t m_nx;                    // grid spaces in x dimension
    double m_alpha;                 // thermal diffusivity (m^2/s)
};


#endif // HEAT_OBJECT1D_HPP
