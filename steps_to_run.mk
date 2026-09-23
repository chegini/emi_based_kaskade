
TARGET := emiModel
NTHREADS ?= 2
MPI_RANKS ?= 2
MPICXX ?= mpicxx
THREAD_ARGS := --nThreads $(NTHREADS)
PROFILE_ARGS ?=
VTK ?= 1
BLAS_ENV := OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
COMPRESSION_BITS ?= 16

# Compression scenarios. The report separates logical raw traffic from the
# transmitted payload and reports the shared codebook estimate separately.
BDDCC_REPORT_ARGS := --bddcCompressionReport 1
BDDCC_QUANTIZATION_ARGS := --bddcCompression 1 --bddcCompressionBits $(COMPRESSION_BITS) \
  --bddcGraphLifting 0 --bddcHuffman 0 --bddcBitlength 0 $(BDDCC_REPORT_ARGS)
BDDCC_HUFFMAN_ARGS := --bddcCompression 1 --bddcCompressionBits $(COMPRESSION_BITS) \
  --bddcGraphLifting 0 --bddcHuffman 1 --bddcBitlength 0 $(BDDCC_REPORT_ARGS)
BDDCC_FULL_ARGS := --bddcCompression 1 --bddcCompressionBits $(COMPRESSION_BITS) \
  --bddcGraphLifting 0 --bddcHuffman 1 --bddcBitlength 1 $(BDDCC_REPORT_ARGS)
BDDCC_GRAPH_ARGS := --bddcCompression 1 --bddcCompressionBits $(COMPRESSION_BITS) \
  --bddcGraphLifting 1 --bddcHuffman 1 --bddcBitlength 1 $(BDDCC_REPORT_ARGS)

COMMON_ARGS := \
  --input ../input_emi_mesh/10Cells3d_10extra_mesh.vtu \
  --extra_set ../input_emi_mesh/10Cells3d_10extra_list_extracellular.txt \
  --intra_set ../input_emi_mesh/10Cells3d_10extra_list_intracellular.txt \
  --excited ../input_emi_mesh/10Cells3d_10extra_early_excited.txt

TIME_ARGS := \
  --finalTime 1 \
  --dt 0.01 \
  --maximumNumberOfTimeSteps 1

SDC_ARGS := \
  --minimumSdcSweeps 2 \
  --maximumSdcSweeps 4

VTK_TEST_TIME_ARGS := \
  --finalTime 0.02 \
  --dt 0.01 \
  --maximumNumberOfTimeSteps 2

BDDC_ARGS := \
  --bddcIterations 100 \
  --bddcVerbose 0

AA_ARGS := \
  --algebraicAdaptivity 1 \
  --algebraicAdaptivityTolerance 1e-4

.PHONY: help build build-mpi clean run-mpi-smoke \
  run-emi run-sdc run-sdc-aa run-bddc run-bddc-quantization \
  run-bddc-huffman run-bddc-full run-bddc-graph run-bddc-compression \
  run-sdc-bddc run-sdc-bddc-aa run-sdc-bddc-compression \
  run-sdc-bddc-aa-compression run-sdc-bddc-aa-vtk-off-test run-all

help:
	@echo "Build:"
	@echo "  make -f steps_to_run.mk build"
	@echo "  make -f steps_to_run.mk build-mpi MPICXX=mpicxx"
	@echo "  make -f steps_to_run.mk run-mpi-smoke MPI_RANKS=2"
	@echo ""
	@echo "Run scenarios:"
	@echo "  make -f steps_to_run.mk run-emi"
	@echo "  make -f steps_to_run.mk run-sdc"
	@echo "  make -f steps_to_run.mk run-sdc-aa"
	@echo "  make -f steps_to_run.mk run-bddc"
	@echo "  make -f steps_to_run.mk run-bddc-quantization"
	@echo "  make -f steps_to_run.mk run-bddc-huffman"
	@echo "  make -f steps_to_run.mk run-bddc-full"
	@echo "  make -f steps_to_run.mk run-bddc-graph"
	@echo "  make -f steps_to_run.mk run-bddc-compression"
	@echo "  make -f steps_to_run.mk run-sdc-bddc"
	@echo "  make -f steps_to_run.mk run-sdc-bddc-aa"
	@echo "  make -f steps_to_run.mk run-sdc-bddc-compression"
	@echo "  make -f steps_to_run.mk run-sdc-bddc-aa-compression"
	@echo "  make -f steps_to_run.mk run-sdc-bddc-aa-vtk-off-test"
	@echo "  make -f steps_to_run.mk run-all"

build:
	$(MAKE) -f Makefile

build-mpi:
	$(MAKE) -f Makefile MPI=1 MPICXX=$(MPICXX)

run-mpi-smoke: build-mpi
	$(BLAS_ENV) $(MPICXX) -show >/dev/null
	$(BLAS_ENV) mpirun -np $(MPI_RANKS) ./$(TARGET) \
	  --mpi 1 --bddc 0 --vtk 0 --finalTime 0.01 --dt 0.01 --maximumNumberOfTimeSteps 1 \
	  $(THREAD_ARGS) $(COMMON_ARGS) --dir output/mpi_smoke

clean:
	$(MAKE) -f Makefile clean

run-emi: build
	$(BLAS_ENV) ./$(TARGET) $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) --dir output/emi --vtk 1

run-sdc: build
	$(BLAS_ENV) ./$(TARGET) --sdc 1 $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) --dir output/sdc --vtk 1

run-sdc-aa: build
	$(BLAS_ENV) ./$(TARGET) --sdc 1 $(THREAD_ARGS) $(AA_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) --dir output/sdc_aa --vtk 1

run-bddc: build
	$(BLAS_ENV) ./$(TARGET) --bddc 1 --bddcCompression 0 --bddcCompressionReport 1 \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(BDDC_ARGS) --dir output/bddc --vtk 1

run-bddc-quantization: build
	mkdir -p output/compression/bddc_quantization
	$(BLAS_ENV) ./$(TARGET) $(BDDCC_QUANTIZATION_ARGS) $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(BDDC_ARGS) \
	  --dir output/compression/bddc_quantization --vtk $(VTK)

run-bddc-huffman: build
	mkdir -p output/compression/bddc_huffman
	$(BLAS_ENV) ./$(TARGET) $(BDDCC_HUFFMAN_ARGS) $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(BDDC_ARGS) \
	  --dir output/compression/bddc_huffman --vtk $(VTK)

run-bddc-full: build
	mkdir -p output/compression/bddc_full
	$(BLAS_ENV) ./$(TARGET) $(BDDCC_FULL_ARGS) $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(BDDC_ARGS) \
	  --dir output/compression/bddc_full --vtk $(VTK)

run-bddc-graph: build
	mkdir -p output/compression/bddc_graph
	$(BLAS_ENV) ./$(TARGET) $(BDDCC_GRAPH_ARGS) $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(BDDC_ARGS) \
	  --dir output/compression/bddc_graph --vtk $(VTK)

run-bddc-compression: build
	$(MAKE) -f steps_to_run.mk run-bddc-graph VTK=1

run-sdc-bddc: build
	$(BLAS_ENV) ./$(TARGET) --sdc 1 --bddc 1 --bddcCompression 0 --bddcCompressionReport 1 \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) $(PROFILE_ARGS) \
	  --dir output/sdc_bddc --vtk $(VTK)

run-sdc-bddc-aa: build
	$(BLAS_ENV) ./$(TARGET) --sdc 1 --bddc 1 --bddcCompression 0 --bddcCompressionReport 1 $(AA_ARGS) \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) $(PROFILE_ARGS) \
	  --dir output/sdc_bddc_aa --vtk $(VTK)

run-sdc-bddc-compression: build
	mkdir -p output/compression/sdc_bddc
	$(BLAS_ENV) ./$(TARGET) --sdc 1 --bddc 1 $(BDDCC_GRAPH_ARGS) \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) $(PROFILE_ARGS) \
	  --dir output/compression/sdc_bddc --vtk $(VTK)

run-sdc-bddc-aa-compression: build
	mkdir -p output/compression/sdc_bddc_aa
	$(BLAS_ENV) ./$(TARGET) --sdc 1 --bddc 1 $(BDDCC_GRAPH_ARGS) $(AA_ARGS) \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) $(PROFILE_ARGS) \
	  --dir output/compression/sdc_bddc_aa --vtk $(VTK)

run-sdc-bddc-aa-vtk-off-test: build
	./$(TARGET) --sdc 1 --bddc 1 $(AA_ARGS) \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(VTK_TEST_TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) \
	  --dir output/vtk_disabled_after_rebuild --vtk 0

run-all: run-emi run-sdc run-sdc-aa run-bddc run-bddc-quantization \
         run-bddc-huffman run-bddc-full run-bddc-graph \
         run-sdc-bddc run-sdc-bddc-aa run-sdc-bddc-compression \
         run-sdc-bddc-aa-compression
