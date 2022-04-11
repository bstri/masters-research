CXXFLAGS = -O0 -std=c++11 -g -I /opt/ohpc/pub/libs/spack/opt/spack/linux-centos7-sandybridge/gcc-4.8.5/gmp-6.1.2-hc7piarcvdeoqjtu3ikpsrapvkbivdbm/include/ -L /opt/ohpc/pub/libs/spack/opt/spack/linux-centos7-sandybridge/gcc-4.8.5/gmp-6.1.2-hc7piarcvdeoqjtu3ikpsrapvkbivdbm/lib/
LIBS = -lgmpxx -lgmp
DEMO = demo
BASIS4 = basis4
BASIS6 = basis6
BASIS7 = basis7
BASIS8 = basis8
MAIN = main
TM = traceMatrixSolve
F = Fadlev
IMPL = mpImpl
TRALGS = TraceAlgs
HEADERS = Matrix.h MPList.h
OBJS = $(F).o $(IMPL).o
EXE = $(DEMO) $(MAIN) $(TM)

.PHONY: all
all: $(EXE)

$(BASIS4): $(OBJS) $(BASIS4).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(BASIS6): $(OBJS) $(BASIS6).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(BASIS7): $(OBJS) $(BASIS7).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(BASIS8): $(OBJS) $(BASIS8).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(DEMO): $(OBJS) $(DEMO).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(MAIN): $(OBJS) $(TRALGS).o $(MAIN).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(TM): $(OBJS) $(TRALGS).o $(TM).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(BASIS4).o: $(HEADERS)

$(BASIS6).o: $(HEADERS)

$(BASIS7).o: $(HEADERS)

$(BASIS8).o: $(HEADERS)

$(DEMO).o: $(F).h $(HEADERS)

$(MAIN).o: $(F).h $(HEADERS)

$(TM).o: $(F).h $(HEADERS)

$(F).o: $(HEADERS)

$(TRALGS).o: Matrix.h	

%.o: %.cpp %.h
	$(CXX) -c $< $(CXXFLAGS)

submit:
	zip FaddeevCppProject.zip -r9 *.cpp *.h *.txt Makefile make_*

clean:
	$(RM) *.o $(EXE)