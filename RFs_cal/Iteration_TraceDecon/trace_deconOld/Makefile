BIN=trace_decon
DATADIR=schemas/css3.0.ext
DATA=decon
SUBDIR=/contrib

PF=trace_decon.pf

ldlibs= -L$(XMOTIFLIB) -L$(BOOSTLIB) -lseisw -lXm -lXt \
	-lgclgrid -lseispp -lperf  -L/usr/local/lib -lgsl -lgslcblas  $(TRLIBS) $(DBLIBS)  \
       $(F77LIBS) $(X11LIBS)

ANTELOPEMAKELOCAL = /opt/antelope/5.7/local/include/antelopemake.local
#$(ANTELOPE)
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)
 CXX=/usr/local/bin/clang
 CXXFLAGS += -I/usr/local/include -I$(BOOSTINCLUDE)  -fopenmp -O1
#-I/usr/local/opt/llvm/include
# LDLIBS += -L/sw2/lib/gcc7/lib -L$(ANTELOPE)/contrib/include
 LDFLAGS += -Wl,-rpath,/usr/local/lib/
#-L/usr/local/opt/llvm/lib
# 
# 
 OBJS= trace_decon.o ComplexArray.o DeconOperator.o dpss.o
 $(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -g -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)