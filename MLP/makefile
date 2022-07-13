
src = $(wildcard src/*.cpp)
obj = $(src:.cpp=.o)

CXXFLAGS = -std=c++11 -I include/ -O3

MLP: $(obj)
	$(CXX) -o $@ $^ $(CXXFLAGS)
clean:
	rm -rf *o MLP