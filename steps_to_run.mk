
TARGET := emiModel
NTHREADS ?= 2
THREAD_ARGS := --nThreads $(NTHREADS)

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

BDDC_ARGS := \
  --bddcIterations 100 \
  --bddcVerbose 0

AA_ARGS := \
  --algebraicAdaptivity 1 \
  --algebraicAdaptivityTolerance 1e-4

.PHONY: help build clean \
  run-emi run-sdc run-sdc-aa run-bddc run-bddc-compression \
  run-sdc-bddc run-sdc-bddc-aa run-sdc-bddc-compression \
  run-sdc-bddc-aa-compression run-all

help:
	@echo "Build:"
	@echo "  make -f steps_to_run.mk build"
	@echo ""
	@echo "Run scenarios:"
	@echo "  make -f steps_to_run.mk run-emi"
	@echo "  make -f steps_to_run.mk run-sdc"
	@echo "  make -f steps_to_run.mk run-sdc-aa"
	@echo "  make -f steps_to_run.mk run-bddc"
	@echo "  make -f steps_to_run.mk run-bddc-compression"
	@echo "  make -f steps_to_run.mk run-sdc-bddc"
	@echo "  make -f steps_to_run.mk run-sdc-bddc-aa"
	@echo "  make -f steps_to_run.mk run-sdc-bddc-compression"
	@echo "  make -f steps_to_run.mk run-sdc-bddc-aa-compression"
	@echo "  make -f steps_to_run.mk run-all"

build:
	$(MAKE) -f Makefile

clean:
	$(MAKE) -f Makefile clean

run-emi: build
	./$(TARGET) $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) --dir output/emi --vtk 1

run-sdc: build
	./$(TARGET) --sdc 1 $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) --dir output/sdc --vtk 1

run-sdc-aa: build
	./$(TARGET) --sdc 1 $(THREAD_ARGS) $(AA_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) --dir output/sdc_aa --vtk 1

run-bddc: build
	./$(TARGET) --bddc 1 --bddcCompression 0 $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(BDDC_ARGS) --dir output/bddc --vtk 1

run-bddc-compression: build
	mkdir -p output/compression/bddc
	./$(TARGET) --bddc 1 --bddcCompression 1 --bddcCompressionBits 16 \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(BDDC_ARGS) --dir output/compression/bddc --vtk 1

run-sdc-bddc: build
	./$(TARGET) --sdc 1 --bddc 1 --bddcCompression 0 \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) \
	  --dir output/sdc_bddc --vtk 1

run-sdc-bddc-aa: build
	./$(TARGET) --sdc 1 --bddc 1 --bddcCompression 0 $(AA_ARGS) \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) \
	  --dir output/sdc_bddc_aa --vtk 1

run-sdc-bddc-compression: build
	mkdir -p output/compression/sdc_bddc
	./$(TARGET) --sdc 1 --bddc 1 --bddcCompression 1 --bddcCompressionBits 16 \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) \
	  --dir output/compression/sdc_bddc --vtk 1

run-sdc-bddc-aa-compression: build
	mkdir -p output/compression/sdc_bddc_aa
	./$(TARGET) --sdc 1 --bddc 1 --bddcCompression 1 --bddcCompressionBits 16 $(AA_ARGS) \
	  $(THREAD_ARGS) $(COMMON_ARGS) $(TIME_ARGS) $(SDC_ARGS) $(BDDC_ARGS) \
	  --dir output/compression/sdc_bddc_aa --vtk 1

run-all: run-emi run-sdc run-sdc-aa run-bddc run-bddc-compression \
         run-sdc-bddc run-sdc-bddc-aa run-sdc-bddc-compression \
         run-sdc-bddc-aa-compression
