FEATURES = 

CXX = ccache g++
#OPT = -std=c++0x -Wall -W -ggdb -DDEBUG
OPT = -std=c++0x -Wall -W -O2 -DNDEBUG
#OPT = -std=c++0x -Wall -W -O3 -ffast-math -DNDEBUG

#CXX = ccache clang++
#OPT = -std=c++0x -stdlib=libc++ -O2 -DNDEBUG

VPATH = io

INCLUDES = -I . -I /home/ge/Programs/lithium/include -I eigen
CXXFLAGS = $(OPT) $(INCLUDES) \
	-DEIGEN_ARRAYBASE_PLUGIN=\"eigen_addons.hh\" \
	-DEIGEN_FUNCTORS_PLUGIN=\"quill_functors.hh\" \
	$(foreach FEATURE, $(FEATURES), -D$(FEATURE))
LDFLAGS = -L /home/ge/Programs/lithium/lib -lli_base \
	-Wl,-rpath=/home/ge/Programs/lithium/lib

OBJS = main.o Basis.o BasisSet.o boys.o CGTO.o CGTOPair.o CGTOQuad.o \
	CGTOShell.o CGTOShellList.o CGTOShellQuad.o CGTOSpecPair.o \
	CGTOSpecQuad_aaaa.o CGTOSpecQuad_aacc.o CGTOSpecQuad_aacd.o \
	CGTOSpecQuad_abcc.o CGTOSpecQuad_abcd.o CommentFilter.o DIIS.o \
	Dispatcher.o Element.o Geometry.o HartreeFock.o IndentFilter.o \
	JobFilter.o LineIStream.o manipulators.o PeriodicTable.o support.o \
	XYZMatrix.o ZMatrix.o

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
