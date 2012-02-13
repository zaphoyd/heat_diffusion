/*
 * Copyright (c) 2011, Peter Thorson. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the WebSocket++ Project nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL PETER THORSON BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

// Heat Library
#include "../heat.hpp"
#include "../object1d.hpp"
#include "../object2d.hpp"
#include "../object3d.hpp"

// WebSocket++ Library
#include "websocketpp.hpp"
#include "wscmd.hpp"

// Boost Library
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

// Standard Library
#include <memory>
#include <vector>
#include <functional>
#include <cstring>
#include <sstream>

using websocketpp::server;

// A request encapsulates all of the information necesssary to perform a request
// the coordinator will fill in this information from the websocket connection 
// and add it to the processing queue. Sleeping in this example is a placeholder
// for any long serial task.
struct request {
	// connection related
    server::handler::connection_ptr con;
    
    // simulation parameters
    int					dimensions;
    int					lx;
    int					ly;
    int					lz;
    size_t				nx;
    size_t				ny;
    size_t				nz;
    initial_condition	initial;
    boundary_style		boundaries;
    double				bvalue;
    size_t				timesteps;
    double				dt;
    double 				alpha;
    
    // simulation state
    bool				stop;
	size_t				zslice;
	smode				mode;
	size_t				callback_interval;
    
    boost::posix_time::ptime start_time;
    boost::posix_time::ptime end_time;
    
    request()
     : dimensions(1),
       lx(1), ly(0), lz(0),
       nx(400), ny(0), nz(0),
       initial(GAUSSIAN),
       boundaries(CONSTANT), bvalue(0.0),
       timesteps(100000), dt(0.001),
       alpha(0.0005),
       stop(false),
	   zslice(0),
	   mode(JSON),
	   callback_interval(100) {}
    
    /// Simulation callback hook
    /* print is called from within the simulation loops to provide an update on
     * the status of the simulation. Included is a reference to the current 
     * simulation data and the timestep that this represents.
     * 
     * The return value of print controls whether or not the simulation should
     * continue or halt immediately. Returning `false` will halt the simulation
     */
    template <typename T>
    bool print(const T& o,size_t ts) {
    	if (mode == JSON) {
    		std::string foo = "{\"type\":\"data\",\"value\":"+o.json(false)+"}";
    		//std::cout << foo << std::endl;
    		con->send(foo);
    	} else if (mode == BINARY) {
    		con->send(o.binary(false), websocketpp::frame::opcode::BINARY);
    	}
		
		return !stop;
	}
    
	bool print(const object3d& o, size_t ts) {
		if (mode == JSON) {
    		//con->send("{\"type\":\"data\",\"value\":"+o.json_slice(zslice,false)+"}");
    	} else if (mode == BINARY) {
    		con->send(o.binary(false), websocketpp::frame::opcode::BINARY);
    	}

		return !stop;
	}

    void process() {
        send_text("Starting Simulation");
        start_time = boost::posix_time::microsec_clock::local_time();
        
        if (dimensions == 1) {
        	object1d object(lx,nx,alpha);
        	object1d source(lx,nx,0.0);
			object.init(initial);
			source.init(FLAT);
			
			stop = false;
			
			object.ftcs(
				timesteps, dt,
				boundaries, bvalue,
				source,
				std::bind(
					&request::print<object1d>,
					this,
					std::placeholders::_1,
					std::placeholders::_2
				),
				callback_interval
			);
        } else if (dimensions == 2) {
			object2d object(lx,ly,nx,ny,alpha);
        	object2d source(lx,ly,nx,ny,0.0);
			object.init(initial);
			source.init(FLAT);
			
			stop = false;
			
			object.ftcs(
				timesteps, dt,
				boundaries, bvalue,
				source,
				std::bind(
					&request::print<object2d>,
					this,
					std::placeholders::_1,
					std::placeholders::_2
				),
				callback_interval
			);
        } else if (dimensions == 3) {
        	send_text("3D simulation not implimented yet");
        } else {
        	send_text("Invalid Dimensions, shouldn't be here");
        }
        
        end_time = boost::posix_time::microsec_clock::local_time();
        boost::posix_time::time_period len(start_time,end_time);
        
        std::stringstream o;
        o << "Simulation Completed in " << len.length();
        
        send_text(o.str());
    }
    
    void send_text(const std::string& msg) {
    	con->send("{\"type\":\"message\",\"value\":\""+msg+"\"}");
    }
    
    void cancel() {
    	stop = true;
    }
};

typedef std::shared_ptr<request> request_ptr;

// The coordinator is a simple wrapper around an STL queue. add_request inserts
// a new request. get_request returns the next available request and blocks
// (using condition variables) in the case that the queue is empty.
class request_coordinator {
public:
    void add_request(request_ptr r) {
        boost::unique_lock<boost::mutex> lock(m_lock);
        m_requests.push(r);
        lock.unlock();
        m_cond.notify_one();
    }
    
    request_ptr get_request() {
        boost::unique_lock<boost::mutex> lock(m_lock);
        
        while (m_requests.empty()) {
            m_cond.wait(lock);
        }
        
        request_ptr r = m_requests.front();
        m_requests.pop();
        return r;
    }
private:
    std::queue<request_ptr>		m_requests;
    boost::mutex                m_lock;
    boost::condition_variable   m_cond;
};

class server_handler : public server::handler {
public:    
    server_handler(request_coordinator& c) : m_coordinator(c) {}
    
    void on_message(connection_ptr con,message_ptr msg) {
        wscmd::cmd command = wscmd::parse(msg->get_payload());
        
        if (command.command == "simulate") {
        	// simulate
        	// dimensions=[1,2,3]
        	// lx = [1-1000]
        	// ly = [1-1000]
        	// lz = [1-1000]
        	// nx = [1-1000]
        	// ny = [1-1000]
        	// nz = [1-1000]
        	// initial = [0]
        	// timesteps = [1-10000000]
        	// boundary = [0,1]
        	// dt = [float]
        	// alpha = [float]
        	
        	request_ptr r(new request());
        	
        	r->con = con;
        	
        	if (command.args["dimensions"] != "") {
        		if (std::atoi(command.args["dimensions"].c_str()) >= 1 && std::atoi(command.args["dimensions"].c_str()) <= 3) {
        			r->dimensions = std::atoi(command.args["dimensions"].c_str());
        		} else {
        			send_text(con,"invalid dimensions");
            		return;
        		}
        	}
        	
        	if (command.args["lx"] != "") {
        		if (std::atoi(command.args["lx"].c_str()) >= 1 && std::atoi(command.args["lx"].c_str()) <= 1000) {
        			r->lx = std::atoi(command.args["lx"].c_str());
        		} else {
        			send_text(con,"invalid lx");
            		return;
        		}
        	}
        	
        	if (command.args["ly"] != "") {
        		if (std::atoi(command.args["ly"].c_str()) >= 1 && std::atoi(command.args["ly"].c_str()) <= 1000 && r->dimensions >= 2) {
        			r->ly = std::atoi(command.args["ly"].c_str());
        		} else {
        			send_text(con,"invalid ly");
            		return;
        		}
        	}
        	
        	if (command.args["lz"] != "") {
        		if (std::atoi(command.args["lz"].c_str()) >= 1 && std::atoi(command.args["lz"].c_str()) <= 1000 && r->dimensions == 3) {
        			r->lz = std::atoi(command.args["lz"].c_str());
        		} else {
        			send_text(con,"invalid lz");
            		return;
        		}
        	}
        	
        	if (command.args["nx"] != "") {
        		if (std::atoi(command.args["nx"].c_str()) >= 1 && std::atoi(command.args["nx"].c_str()) <= 10000) {
        			r->nx = std::atoi(command.args["nx"].c_str());
        		} else {
        			send_text(con,"invalid nx");
            		return;
        		}
        	}
        	
        	if (command.args["ny"] != "") {
        		if (std::atoi(command.args["ny"].c_str()) >= 1 && std::atoi(command.args["ny"].c_str()) <= 10000 && r->dimensions >= 2) {
        			r->ny = std::atoi(command.args["ny"].c_str());
        		} else {
        			send_text(con,"invalid ny");
            		return;
        		}
        	}
        	
        	if (command.args["nz"] != "") {
        		if (std::atoi(command.args["nz"].c_str()) >= 1 && std::atoi(command.args["nz"].c_str()) <= 10000 && r->dimensions == 3) {
        			r->nz = std::atoi(command.args["nz"].c_str());
        		} else {
        			send_text(con,"invalid nz");
            		return;
        		}
        	}
        	
        	if (command.args["timesteps"] != "") {
        		if (std::atoi(command.args["timesteps"].c_str()) >= 1 && std::atoi(command.args["timesteps"].c_str()) <= 100000000) {
        			r->timesteps = std::atoi(command.args["timesteps"].c_str());
        		} else {
        			send_text(con,"invalid timesteps");
            		return;
        		}
        	}
        	
        	if (command.args["initial"] != "") {
        		if (std::atoi(command.args["initial"].c_str()) == 0) {
        			r->initial = initial_condition(std::atoi(command.args["initial"].c_str()));
        		} else {
        			send_text(con,"invalid initial");
            		return;
        		}
        	}
        	
        	if (command.args["boundary"] != "") {
        		if (std::atoi(command.args["boundary"].c_str()) >= 0 && std::atoi(command.args["boundary"].c_str()) <= 1 ) {
        			r->boundaries = boundary_style(std::atoi(command.args["boundary"].c_str()));
        		} else {
        			send_text(con,"invalid boundary");
            		return;
        		}
        	}
        	
        	if (command.args["callback_interval"] != "") {
        		if (std::atoi(command.args["callback_interval"].c_str()) > 0 && std::atoi(command.args["callback_interval"].c_str()) <= r->timesteps ) {
        			r->callback_interval = std::atoi(command.args["callback_interval"].c_str());
        		} else {
        			send_text(con,"invalid callback_interval");
            		return;
        		}
        	}
        	
        	if (command.args["smode"] != "") {
        		if (command.args["smode"] == "0") {
        			r->mode = JSON;
        		} else if (command.args["smode"] == "1") {
        			r->mode = BINARY;
        		} else {
        			send_text(con,"invalid smode");
            		return;
        		}
        	}
        	
        	m_last_request = r;
        	
        	m_coordinator.add_request(r);
        } else if (command.command == "cancel") {
        	if (m_last_request) {
        		send_text(con,"Attempting to cancel simulation...");
        		m_last_request->cancel();
        	} else {
        		send_text(con,"No simulation to cancel.");
        	}
        } else {
        	send_text(con,"unrecognized message");
            return;
        }
    }
    
    void send_text(connection_ptr con,const std::string& msg) {
    	con->send("{\"type\":\"message\",\"value\":\""+msg+"\"}");
    }
private:
    request_coordinator&	m_coordinator;
    request_ptr				m_last_request;
};

// process_requests is the body function for a processing thread. It loops 
// forever reading requests, processing them serially, then reading another 
// request. 
void process_requests(request_coordinator* coordinator);

void process_requests(request_coordinator* coordinator) {
    request_ptr r;
    
    while (1) {
        r = coordinator->get_request();
        
        r->process();
    }
}

// concurrent server takes two arguments. A port to bind to and a number of 
// worker threads to create. The thread count must be an integer greater than
// or equal to zero.
// 
// num_threads=0 Standard non-threaded WebSocket++ mode. Handlers will block
//               i/o operations and other handlers.
// num_threads=1 One thread processes requests serially the other handles i/o
//               This allows new connections and requests to be made while the
//               processing thread is busy, but does allow long jobs to 
//               monopolize the processor increasing request latency.
// num_threads>1 Multiple processing threads will work on the single queue of
//               requests provided by the i/o thread. This enables out of order
//               completion of requests. The number of threads can be tuned 
//               based on hardware concurrency available and expected load and
//               job length.
int main(int argc, char* argv[]) {
    unsigned short port = 9002;
    unsigned short num_threads = 2;
    
    try {
        if (argc == 2) {
            std::stringstream buffer(argv[1]);
            buffer >> port;
        }
        
        if (argc == 3) {
            std::stringstream buffer(argv[2]);
            buffer >> num_threads;
        }
            
        request_coordinator rc;
        
        server::handler::ptr h;
        if (num_threads == 0) {
            std::cout << "this is a bad idea" << std::endl;
            return 0;
        } else {
            h = server::handler::ptr(new server_handler(rc));
        }
        
        server echo_endpoint(h);
        
        echo_endpoint.alog().unset_level(websocketpp::log::alevel::ALL);
        echo_endpoint.elog().unset_level(websocketpp::log::elevel::ALL);
        
        echo_endpoint.elog().set_level(websocketpp::log::elevel::ERROR);
        echo_endpoint.elog().set_level(websocketpp::log::elevel::FATAL);
        
        std::list<boost::shared_ptr<boost::thread> > threads;
        
        for (int i = 0; i < num_threads; i++) {
            threads.push_back(boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&process_requests, &rc))));
        }
        
        std::cout << "Starting WebSocket sleep server on port " << port << " with " << num_threads << " processing threads." << std::endl;
        echo_endpoint.listen(port);
    } catch (std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    
    return 0;
}
