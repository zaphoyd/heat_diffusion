# 3D Heat Diffusion solver

## Overview of components

### Core solver

Objects can be created with any length and subdivided into discrete partitions. Length (lx,ly,lz) and partitioning (nx,ny,nz) are independent for each dimension. Objects can be initialized using a FLAT or GAUSSIAN method.

Two solvers are (mostly) implimented. Forward-Time Central-Space is implimented for all three dimensions. Crank Nicholson is implimented in 1D and 2D. The matrix and object storage/serialization code was designed before I understood exactly how Crank-Nicholson actually worked. As a result it is heavily optimized for FTCS (dimensionally independent lookups) and embarassingly slow at CN which needs memory contiguous matricies and vectors to really work at all. I have not implimented 3D Crank Nicholson because without rewriting my matrix and object data structures it would be too slow to even test for correctness.

Both solvers support adding an arbitrary time independent source term and both flat (with the option to specify a non-zero boundary value) and periodic where the boundary conditions wrap around.

## heat_test
This is a simple command line utility that bootstraps the solver code and allows running the solver in small input tests. It has no external dependencies beyond the C++11 STL.

## heat_server
This is a more complicated simulation server that uses WebSocket to communicate with a command and visualization client running in a browser. The server accepts input in the form of a simple command language that allows specifying the simulation input parameters. Visualization output is performed either as a 2D plot (1D problems only) or as a 2D grayscale image. White represents hot, black cold.

The command language works as follows:

`command:arg1=val1;arg2=val2;`

trailing punctuation is always retained.

The commands recognized are `simulate` and `cancel`. Simulate will start a simulation on one of the processing threads. Cancel will signal the simulation thread to cancel the simulation the next time that processing thread checks in. Processing threads check in to test for canceling and deliver a simulation snapshot+timestep at a specifiable interval. Simulator I/O is handled in one thread and there are a configurable number of simulation processing threads. Simulations are assigned to processing threads in FIFO order.

Arguments for the `simulate` command:

| Argument | Effect | Required? | Allowed Values | Default |
| --- | --- | --- | --- | --- |
| `timesteps` | how many iterations to run. | No |  1-1000000 | 10000 |
| `dimensions` | how many dimensions to use for the simulated object. | No | 1,2,3 | 1 |
| `lx` | length of x dimension in m | Yes | 1-1000 | 1 |
| `ly` | length of y dimension in m | 2D and 3D only | 1-1000 | 0 |
| `lz` | length of z dimension in m | 3D only | 1-1000 | 0 |
| `nx` | partitions in x dimension in m | Yes | 1-1000 | 400 |
| `ny` | partitions in y dimension in m | 2D and 3D only | 1-1000 | 0 |
| `nz` | partitions in z dimension in m | 3D only | 1-1000 | 0 |
| `initial` | initialization strategy, `0=FLAT,1=GAUSSIAN` | No | 0,1 | 1 | 
| `boundary` | boundary handling strategy, `0=CONSTANT,1=PERIODIC` | No | 0,1 | 0| 
| `method` | which solver to use, `0=FTCS,1=CRANK_NICHOLSON,2=JACOBI,3=GAUSS-SEIDEL,4=SOR` | No | 0,1,2,3,4 | 0 | 
| `zslice` | which z slice to display in 3D simulations | No | 0,nz-1 | 0 |
| `callback_interval` | timestep interval used to check for cancel and return intermediate results. | No | 1,`timesteps` | 100 |
| `smode` | simulation snapshot format, `0=JSON,1=BINARY`, `heat_client` automatically chooses and sends this value based on browser capabilities | No | 0,1 | 0 | 
        

# Examples
| Command | Effect |
| --- | --- |
| `simulate:` | Run simulation with default parameters |
| `cancel:` | Cancel the currently running simulation (requires server to be started in one of the multithreaded modes) |
| `simulate:timesteps=100000;` | Run default simulation but adjust timesteps |
| `simulate:timesteps=1000;dimensions=2;nx=100;ny=100;ly=1.0;` | Run 2D simulation with default solver and custom problem size |
| `simulate:timesteps=100000;dimensions=2;nx=100; ny=100;ly=1.0;callback_interval=1000;` | Run a longer simulation and ask for progress update renderings every 1000 timesteps |
| `simulate:timesteps=25;method=4;nx=100;ny=100;nz=100; ly=1;lz=1;dimensions=3;callback_interval=5;zslice=50;dt=1.0;` | 3D simulation long enough to show a nice average iterations spread |
| `simulate:timesteps=30;method=2;nx=400;ny=400; ly=1;dimensions=2;callback_interval=1;dt=1.0;` | Nice, big, well performing 2D simulation |

# Usage

## Browser Support

The client will operate on any modern browser with WebSocket support. This includes:
Safari 5.0.1+
iOS 4.2+
FireFox 6+
Chrome 6+
IE10+

On browsers with full RFC6455 support it will use binary websocket for drastically improved visualization performance. These browsers are:
Chrome 16+
FireFox 11+

## Building the server

Building the server requires Boost and WebSocket++ (github.com/zaphoyd/websocketpp). If you have trouble building it let me know.

### Example make statement:

`make CXX=g++ DEBUG=0 WEBSOCKETPP_PATH=/path/to/websocketpp/headers WEBSOCKETPP=/path/to/websocketpplib.a BOOST_INCLUDE_PATH=/path/to/boost/headers BOOST_LIB_PATH=/path/to/boost/libs`

# Performance Analysis:

## FTCS 3D 1000 time steps

| Size | Time |
| --- | --- |
| 10x10x10 cube | 00:00:00.143543 |
| 20x20x20 cube | 00:00:01.131008 |
| 30x30x30 cube | 00:00:03.829194 |
| 40x40x40 cube | 00:00:09.083705 |
| 50x50x50 cube | 00:00:17.804704 |
| 60x60x60 cube | 00:00:30.732263 |
| 70x70x70 cube | 00:00:48.869378 |
| 80x80x80 cube | 00:01:12.005312 |
| 90x90x90 cube | 00:01:44.246772 |
| 100x100x100 cube | 00:02:22.738739 |

| size | µs | µs/size | µs/timestep |
| --- | --- | --- | --- |
|1,000 | 143543 | 143.543 | 143.543 |
|8,000 | 1131008 | 141.376 | 1131.008 |
|27,000 | 3829194 | 141.822 | 3829.194 |
|64,000 | 9083705 | 141.932 | 9083.705 |
|125,000 | 17804704 | 142.437 | 17804.704 |
|216,000 | 30732263 | 142.278 | 30732.263 |
|343,000 | 48869378 | 142.476 | 48869.378 |
|512,000 | 72005312 | 140.635 | 72005.312 |
|729,000 | 104246772 | 142.999 | 104246.772 |
|1,000,000 | 142738739 | 142.738 | 142738.739 |

## Performance Analysis of iterative methods

### 2D
300x300 square simulated for 100 seconds
| Method | Time | dt | avg iterations |
| --- | --- | --- | --- |
| FTCS | 00:00:49.067901 | 0.001s | N/A |
| JACOBI | 00:00:23.413229 | 1.0s | 325 |
| GAUSS-SEIDEL | 00:00:20.495972 | 1.0s | 221 |
| SOR | 00:00:07.686357 | 1.0s | 75 |

### 3D
100x100x100 cube simulated for 10 seconds
| Method | Time | dt | avg iterations |
| --- | --- | --- | --- |
| FTCS | 00:00:10.090005 | 0.01s | N/A |
| JACOBI | 00:00:39.082192 | 0.1s | 28 |
| GAUSS-SEIDEL | 00:00:23.739189 | 0.1s | 17 |
| SOR | 00:00:10.210998 | 0.1s | 6 |

100x100x100 cube simulated for 25 seconds
| Method | Time | dt | avg iterations |
| --- | --- | --- | --- |
|FTCS | 00:00:24.048602 | 0.01s | N/A |
|JACOBI | 00:01:16.475670 | 1.0s | 238 |
|GAUSS-SEIDEL | 00:00:38.909941 | 1.0s | 131 |
|SOR | 00:00:10.885491 | 1.0s | 33 |

At higher dt FTCS blows up and is wildly inaccurate. At shorter simulation
periods, FTCS is competitive. At longer time periods SOR pulls ahead drastically
as FTCS is stuck on smaller timesteps.
