include /kaskade7/Makefile.Local
include $(KASKADE7)/Makefile.Rules

TARGET := emiModel
SOURCES := EMI_model.cpp
OBJECTS := $(SOURCES:.cpp=.o)
DEPS := $(SOURCES:.cpp=.d)

# MPI is opt-in for this executable. The MPI compiler wrapper supplies the
# MPI headers and libraries; KASKADE_HAVE_MPI guards MPI-only source code.
MPI ?= 0
MPICXX ?= mpicxx
ifeq ($(MPI),1)
  CXX := $(MPICXX)
  # MUMPS ships a compatibility mpi.h in /usr/include/mumps.  Put the
  # selected MPI implementation first so <mpi.h> resolves to OpenMPI.
  MPI_INCLUDE_FLAGS ?= -I/usr/lib64/mpi/gcc/openmpi4/include
  MPI_FLAGS := $(MPI_INCLUDE_FLAGS) -DKASKADE_HAVE_MPI
  FLAGS := $(MPI_FLAGS) $(FLAGS)
else
  MPI_FLAGS :=
endif

.PHONY: all kaskade clean help

all: $(TARGET)

kaskade:
	$(MAKE) -C $(KASKADE7) kasklib

$(KASKADE7)/libs/libkaskade.a:
	$(MAKE) -C $(KASKADE7) kasklib

$(TARGET): $(OBJECTS) $(KASKADE7)/libs/libkaskade.a
	$(CXX) $(FLAGS) $(OBJECTS) $(KASKADE7)/libs/*.o $(KASKADELIB) \
	$(DUNELIB) \
	$(UGLIB) \
	$(BOOSTLIB) \
	$(DIRECTSOLVERLIB) \
	$(ITSOLLIB) \
	$(HYPRELIB) \
	$(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@

clean:
	rm -f $(TARGET) $(OBJECTS) $(DEPS) gccerr.txt *.d.extended
	rm -f *.vtu *.vtk *.m *.txt
	rm -rf output matlab_dir

help:
	@echo "Targets:"
	@echo "  make          Build $(TARGET) from $(SOURCES)"
	@echo "  make kaskade  Build the Kaskade library in /kaskade7"
	@echo "  make clean    Remove scratch build/output files"
	@echo "  make MPI=1 MPICXX=mpicxx  Build the MPI-enabled smoke-test executable"
	@echo ""
	@echo "Run inside Podman with:"
	@echo "  podman run -it --rm -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 -v /data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC:/BDDC -w /BDDC/emiModelBddcSdc_new_scratch localhost/kaskade7-build-environment:latest"

-include $(DEPS)
