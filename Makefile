INCLUDES = -I /home/ge/Programs/lithium/include -I eigen
#CXXFLAGS = -Wall -W -Weffc++ -ggdb $(INCLUDES)
CXXFLAGS = -Wall -W -ggdb $(INCLUDES)
LDFLAGS = -L /home/ge/Programs/lithium/lib -lli_base \
	-Wl,-rpath=/home/ge/Programs/lithium/lib

OBJS = main.o BasisSet.o Element.o Geometry.o Indenter.o \
	JobFilter.o LineGetter.o LineIStream.o PeriodicTable.o support.o

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
