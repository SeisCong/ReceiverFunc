BIN=trace_decon
DATADIR=schemas/css3.0.ext
DATA=decon
SUBDIR=/contrib

PF=trace_decon.pf

ldlibs=-L$(XMOTIFLIB) -L$(BOOSTLIB) -lseisw -lXm -lXt \
	-lgclgrid -lseispp -lperf -L/usr/local/lib -lgsl -lgslcblas -lseispp $(TRLIBS) $(DBLIBS)  \
       $(F77LIBS) $(X11LIBS) -lseisw

ANTELOPEMAKELOCAL = $(ANTELOPE)/local/include/antelopemake.local
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)
CXX=g++
CXXFLAGS += -I/usr/local/include -I$(BOOSTINCLUDE) -fopenmp -O3
LDLIBS += -L/sw/lib/gcc5/lib 
LDFLAGS += -Wl,-rpath,/usr/local/lib/ 

OBJS= trace_decon.o ComplexArray.o DeconOperator.o dpss.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
