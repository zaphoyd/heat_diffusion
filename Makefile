CFLAGS = -O2 -std=c++0x
LDFLAGS = 

CXX		?= clang++
LDFLAGS := $(LDFLAGS) 

ifeq ($(CXX), clang++)
	CFLAGS := -stdlib=libc++
endif

heat_diffusion: main.o
	$(CXX) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $^

# cleanup by removing generated files
#
.PHONY:		clean
clean:
		rm -f *.o heat_diffusion
