CXX = g++
CXXFLAGS = -O3 -Wall
BIN = introDist Dxyz

all: $(BIN)

introDist: introDist.cpp introDist.h generalUtils.cpp generalUtils.h
	$(CXX) $(CXXFLAGs) introDist.cpp generalUtils.cpp -o introDist
Dxyz: Dxyz.cpp Dxyz.h
	$(CXX) $(CXXFLAGS) Dxyz.cpp -o Dxyz -lz -lboost_iostreams

clean:
	rm -f $(BIN) *.o *.d

.PHONY: clean all
