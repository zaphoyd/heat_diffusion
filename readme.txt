3D Heat Diffusion solver.

Overview of components:

1. Core solver
heat.hpp
matrix.hpp
object1d.hpp
object2d.hpp
object3d.hpp

Features of the solver:
Objects can be created with any length and subdivided into discrete partitions. Length (lx,ly,lz) and partitioning (nx,ny,nz) are independent for each dimension. Objects can be initialized using a FLAT or GAUSSIAN method.

Two solvers are (mostly) implimented. Forward-Time Central-Space is implimented for all three dimensions. Crank Nicholson is implimented in 1D and 2D. The matrix and object storage/serialization code was designed before I understood exactly how Crank-Nicholson actually worked. As a result it is heavily optimized for FTCS (dimensionally independent lookups) and embarassingly slow at CN which needs memory contiguous matricies and vectors to really work at all. I have not implimented 3D Crank Nicholson because without rewriting my matrix and object data structures it would be too slow to even test for correctness.

Both solvers support adding an arbitrary time independent source term and both flat (with the option to specify a non-zero boundary value) and periodic where the boundary conditions wrap around.

2. heat_test
This is a simple command line utility that bootstraps the solver code and allows running the solver in small input tests. It has no external dependencies beyond the C++11 STL.

3. heat_server
This is a more complicated simulation server that uses WebSocket to communicate with a command and visualization client running in a browser. The server accepts input in the form of a simple command language that allows specifying the simulation input parameters. Visualization output is performed either as a 2D plot (1D problems only) or as a 2D grayscale image. White represents hot, black cold.

The command language works as follows:

command:arg1=val1;arg2=val2;

trailing punctuation is always retained.

The commands recognized are `simulate` and `cancel`. Simulate will start a simulation on one of the processing threads. Cancel will signal the simulation thread to cancel the simulation the next time that processing thread checks in. Processing threads check in to test for canceling and deliver a simulation snapshot+timestep at a specifiable interval. Simulator I/O is handled in one thread and there are a configurable number of simulation processing threads. Simulations are assigned to processing threads in FIFO order.

Arguments for the `simulate` command:
timesteps - how many iterations to run. Allowed values [1-1000000] default: 10000

dimensions - how many dimensions to use for the simulated object. Allowed values [1,2,3] default: 1

lx - length of x dimension in m. Required. Allowed values [1-1000], default: 1
ly - length of y dimension in m. Required for 2D and 3D. Allowed values [1-1000], default: 0
lz - length of z dimension in m. Required for 3D. Allowed values [1-1000], default: 0

nx - partitions in x dimension in m. Required. Allowed values [1-1000], default: 400
ny - partitions in y dimension in m. Required for 2D and 3D. Allowed values [1-1000], default: 0
nz - partitions in z dimension in m. Required for 3D. Allowed values [1-1000], default: 0

initial - initialization strategy. Allowed values [0,1], default: 1; 0=FLAT,1=GAUSSIAN
boundary - boundary handling strategy. Allowed values [0,1], default: 0; 0=CONSTANT,1=PERIODIC

method - which solver to use. Allowed values [0,1], default: 0; 0=FTCS,1=CRANK_NICHOLSON
zslice - which z slice to display in 3D simulations. Allowed values [0,nz-1]. default: 0

callback_interval - timestep interval used to check for cancel and return intermediate results.
                    Allowed values [1,timesteps], default: 100
smode - simulation snapshot format. Allowed values [0,1] default: 0; 0=JSON,1=BINARY.
        heat_client automatically chooses and sends this value based on browser capabilities

Examples
simulate:
cancel:

The client will operate on any modern browser with WebSocket support. This includes:
Safari 5.0.1+
iOS 4.2+
FireFox 6+
Chrome 6+
IE10+

On browsers with full RFC6455 support it will use binary websocket for drastically improved visualization performance. These browsers are:
Chrome 16+
FireFox 11+

Building the server requires Boost 1.47+ and WebSocket++ (github.com/zaphoyd/websocketpp) both built in C++11 mode. I have an example server running that you can play with at heat.zaphoyd.net.

Example make statement:

make CXX=g++4.6 DEBUG=0 WEBSOCKETPP_PATH=/path/to/websocketpp/headers WEBSOCKETPP=/path/to/websocketpplib.a BOOST_INCLUDE_PATH=/path/to/boost/headers BOOST_LIB_PATH=/path/to/boost/libs

