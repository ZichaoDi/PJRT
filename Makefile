CFLAGS	         =

ODIR=obj
SRCDIR=.

include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules

_DEPS =
DEPS  = $(patsubst %,$(SRCDIR)/%,$(_DEPS))

_OBJ  =
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	@mkdir -p $(@D)
	${PETSC_COMPILE} -c -o $@ $< $(CFLAGS) ${PETSC_KSP_LIB} ${PETSC_CC_INCLUDES}

all: prt

pjrt: $(ODIR)/prt.o $(OBJ)
	-${CLINKER} -o $@ $^ $(CFLAGS) ${PETSC_KSP_LIB}

.PHONY: clean

clean:
	rm -f $(OBJ)
	rm -f pjt
