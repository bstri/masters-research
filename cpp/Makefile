CXXFLAGS = -Wall -std=c++14 -O0 -g
LIBS = -lgmpxx -lgmp
DEMO = demo
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

$(DEMO): $(OBJS) $(DEMO).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(MAIN): $(OBJS) $(TRALGS).o $(MAIN).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(TM): $(OBJS) $(TRALGS).o $(TM).o
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

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