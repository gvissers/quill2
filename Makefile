FEATURES = 
FEATURES += DEBUG

CXX = ccache g++
OPT = -std=c++0x -Wall -W -ggdb
#OPT = -std=c++0x -Wall -W -O2

VPATH = io gaussint

INCLUDES = -I . -I /home/ge/Programs/lithium/include -I eigen
CXXFLAGS = $(OPT) $(INCLUDES) \
	$(foreach FEATURE, $(FEATURES), -D$(FEATURE))
LDFLAGS = -L /home/ge/Programs/lithium/lib -lli_base \
	-Wl,-rpath=/home/ge/Programs/lithium/lib

OBJS = main.o Basis.o BasisSet.o boys.o CGTO.o CGTOPair.o CGTOQuad.o \
	CommentFilter.o Dispatcher.o Element.o Geometry.o gto_elec_rep.o \
	gto_kinetic.o gto_one_elec.o gto_nuc_attr.o gto_overlap.o \
	IndentFilter.o JobFilter.o LineIStream.o manipulators.o \
	PeriodicTable.o support.o XYZMatrix.o ZMatrix.o

.PHONY: doc

all: main

main: $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	\rm -f *.o *.d main

doc:
	doxygen

%.d: %.cc *.hh
	set -e;\
	$(CXX) -MM $(CXXFLAGS) $< | sed 's/\($*\).o[ :]*/\1.o $@: /g' > $@

ifneq "$(MAKECMDGOALS)" "clean"
-include $(OBJS:.o=.d)
endif
