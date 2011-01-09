INCLUDES = -I /home/ge/Programs/lithium/include -I eigen
#CXXFLAGS = -Wall -W -Weffc++ -ggdb $(INCLUDES)
CXXFLAGS = -Wall -W -ggdb $(INCLUDES)
LDFLAGS = -L /home/ge/Programs/lithium/lib -lli_base \
	-Wl,-rpath=/home/ge/Programs/lithium/lib

OBJS = main.o BasisSet.o Element.o LineGetter.o PeriodicTable.o support.o

all: main

main: $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	\rm -f *.o *.d main

%.d: %.cc *.hh
	set -e;\
	$(CXX) -MM -MG $(CXXFLAGS) $< | sed 's/\($*\).o[ :]*/\1.o $@: /g' > $@

ifneq "$(MAKECMDGOALS)" "clean"
-include $(OBJS:.o=.d)
endif
