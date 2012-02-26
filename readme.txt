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
simulate:timesteps=100000;
simulate:timesteps=1000;dimensions=2;nx=100;ny=100;ly=1.0;
simulate:timesteps=100000;dimensions=2;nx=100;ny=100;ly=1.0;callback_interval=1000;

The client will operate on any modern browser with WebSocket support. This includes:
Safari 5.0.1+
iOS 4.2+
FireFox 6+
Chrome 6+
IE10+

On browsers with full RFC6455 support it will use binary websocket for drastically improved visualization performance. These browsers are:
Chrome 16+
FireFox 11+

Building the server requires Boost and WebSocket++ (github.com/zaphoyd/websocketpp). If you have trouble building it let me know. I also have an example server running that you can play with at heat.zaphoyd.net

Example make statement:

make CXX=g++ DEBUG=0 WEBSOCKETPP_PATH=/path/to/websocketpp/headers WEBSOCKETPP=/path/to/websocketpplib.a BOOST_INCLUDE_PATH=/path/to/boost/headers BOOST_LIB_PATH=/path/to/boost/libs

I have successfully built heat_server on the linux cluster using the following:
git clone git://github.com/zaphoyd/websocketpp.git
git clone git://github.com/zaphoyd/heat_diffusion.git
cd websocketpp
make
cd ../heat_diffusion/heat_server
make CXX=g++ CSPP_LINUX=1 WEBSOCKETPP=../../websocketpp/libwebsocketpp.a WEBSOCKETPP_PATH=../../websocketpp/src

The server runs but there appears to be firewalls in place that prevent anything but localhost from actually talking to it. You should be able to use the local heat_client.html file to connect to a locally running server at ws://localhost:9002



Performance Analysis:

FTCS 3D 1000 time steps
10x10x10 cube: 00:00:00.143543
20x20x20 cube: 00:00:01.131008
30x30x30 cube: 00:00:03.829194
40x40x40 cube: 00:00:09.083705
50x50x50 cube: 00:00:17.804704
60x60x60 cube: 00:00:30.732263
70x70x70 cube: 00:00:48.869378
80x80x80 cube: 00:01:12.005312
90x90x90 cube: 00:01:44.246772
100x100x100 cube: 00:02:22.738739

size      µs          µs/size  µs/timestep
1,000     143543      143.543  143.543        
8,000     1131008     141.376  1131.008
27,000    3829194     141.822  3829.194
64,000    9083705     141.932  9083.705
125,000	  17804704	  142.437  17804.704
216,000	  30732263	  142.278  30732.263
343,000	  48869378	  142.476  48869.378
512,000	  72005312	  140.635  72005.312
729,000	  104246772	  142.999  104246.772
1,000,000  142738739  142.738  142738.739


2012-02-26 Update
- Added JACOBI, GAUSS-SEIDEL, and SOR (successive over-relaxation) solvers for
  1D, 2D, and 3D.

Examples:
3D simulation long enough to show a nice average iterations spread:
simulate:timesteps=25;method=4;nx=100;ny=100;nz=100;ly=1;lz=1;dimensions=3;callback_interval=5;zslice=50;dt=1.0;

Nice, big, well performing 2D simulation:
simulate:timesteps=30;method=2;nx=400;ny=400;ly=1;dimensions=2;callback_interval=1;dt=1.0;

Performance Analysis of iterative methods

2D
300x300 square simulated for 100 seconds with dt=1.0s (except FTCS)
FTCS: 00:00:49.067901 (dt=0.001s)
JACOBI: 00:00:23.413229 (avg iterations: 325)
HAUS-SEIDEL: 00:00:20.495972 (avg iterations: 221)
SOR: 00:00:07.686357 (avg iterations: 75)

3D
100x100x100 cube simulated for 10 seconds with dt=0.1s (except FTCS)
FTCS: 00:00:10.090005 (dt=0.01s)
JACOBI: 00:00:39.082192, average iterations: 28
GAUSS-SEIDEL: 00:00:23.739189, average iterations: 17
SOR: 00:00:10.210998, average iterations: 6

100x100x100 cube simulated for 25 seconds with dt=1.0s (except FTCS)
FTCS: 00:00:24.048602 (dt=0.01s)
JACOBI: 00:01:16.475670 (avg iterations: 238)
GAUSS-SEIDEL: 00:00:38.909941 (avg iterations: 131)
SOR: 00:00:10.885491 (avg iterations: 33)

At higher dt FTCS blows up and is wildly inaccurate. At shorter simulation
periods, FTCS is competitive. At longer time periods SOR pulls ahead drasticly
as FTCS is stuck on smaller timesteps.
