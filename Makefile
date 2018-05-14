CFLAGS	         = -g -w2
LIBS=-L/home/b216449/soft/boost_1_67_0
ODIR=obj
SRCDIR=src

include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules

_DEPS = geom.hpp
DEPS  = $(patsubst %,$(SRCDIR)/%,$(_DEPS))

_OBJ  = pirt.o geom.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SRCDIR)/%.cxx $(DEPS)
	@mkdir -p $(@D)
	${PETSC_CXXCOMPILE} -c -o $@ $< $(CFLAGS) ${PETSC_KSP_LIB} ${PETSC_CC_INCLUDES} ${LIBS}

all: pirt

pirt: $(OBJ)
	@echo $(OBJ)
	-${CLINKER} -o $@ $^ $(CFLAGS) ${PETSC_KSP_LIB} ${LIBCS}

.PHONY: clean

clean:
	rm -f $(OBJ)
	rm -f pirt
