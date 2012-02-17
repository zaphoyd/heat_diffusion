BOOST_LIB_PATH		= /Users/zaphoyd/Documents/boost_1_48_0_libcpp/stage/lib
BOOST_INCLUDE_PATH	= /Users/zaphoyd/Documents/boost_1_48_0_libcpp
WEBSOCKETPP			= /Users/zaphoyd/Documents/websocketpp/libwebsocketpp_cpp11.a
WEBSOCKETPP_PATH	= /Users/zaphoyd/Documents/websocketpp/src
CPP11               = -std=c++0x
DEBUG				= 0

CFLAGS = $(CPP11) -I$(BOOST_INCLUDE_PATH) -I$(WEBSOCKETPP_PATH)
LDFLAGS = -L$(BOOST_LIB_PATH)

CXX		= clang++

ifeq ($(CXX), clang++)
	CFLAGS += -stdlib=libc++
endif

ifeq ($(DEBUG), 0)
	CFLAGS += -O2 -DNDEBUG
else
	CFLAGS += -g -O0
endif

LDFLAGS += $(WEBSOCKETPP) $(BOOST_LIB_PATH)/libboost_system.a $(BOOST_LIB_PATH)/libboost_date_time.a $(BOOST_LIB_PATH)/libboost_regex.a $(BOOST_LIB_PATH)/libboost_thread.a -lpthread

heat_server: heat_server.o
	$(CXX) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $^

# cleanup by removing generated files
#
.PHONY:		clean
clean:
		rm -f *.o heat_server
