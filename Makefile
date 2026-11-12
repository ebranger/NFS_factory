# Linux build of NFS_factory (the Windows build uses the Visual Studio solutions, see "Compiling the code.txt").
#
#   make           the CPU build, build/NFS_factory, and the utilities
#   make gpu       the GPU build, build/NFS_factory_gpu, with the CUDA Toolkit
#   make utils     the utilities of Utilities/: build/RelationsToBinaryFormat and build/RelationFilterSplit
#   make both      all of them
#   make clean
#
# Needs MPIR built and installed with its C++ interface (see "Compiling the code.txt"), linked statically as in the Windows build.
# Linux/include/intrin.h stands in for MSVC's intrinsics. Utilities/RelationFilterSplit.cpp is saved as UTF-16, which GCC cannot
# read, so it is compiled from a UTF-8 copy made by iconv in build/utils/.
#
# Settings, on the command line (make gpu CUDA_ARCH="75 86" NVCC_CCBIN=g++-13):
#   MPIR_DIR     where MPIR is installed, with include/mpir.h and lib/libmpir.a (default /opt/mpir)
#   OPT          optimization flags of the host code (default -O2 -march=native: the program runs on the machine it is built on)
#   CUDA_ARCH    compute capabilities to build for, "native" (the GPUs of this machine, the default) or a list such as
#                "75 80 86 89 90 120", which also adds PTX for the last one
#   NVCC         the CUDA compiler (default nvcc)
#   NVCC_CCBIN   the host compiler for nvcc, if the default one is newer than the CUDA version supports (CUDA 12.4: GCC 13 at most)
#   CUDA_LIBDIR  the folder of libcudart_static.a, if it is not in a default library folder

CC = gcc
CXX = g++
NVCC ?= nvcc
NVCC_CCBIN ?=
CUDA_ARCH ?= native
MPIR_DIR ?= /opt/mpir
CUDA_LIBDIR ?= $(dir $(shell command -v $(NVCC)))../lib64
OPT ?= -O2 -march=native

INCLUDES = -ILinux/include -I$(MPIR_DIR)/include -Iint_ntt
CPPFLAGS = $(INCLUDES) -MMD -MP
CFLAGS = $(OPT) -fopenmp
CXXFLAGS = -std=c++17 $(OPT) -fopenmp
LIBS = $(MPIR_DIR)/lib/libmpirxx.a $(MPIR_DIR)/lib/libmpir.a

ifeq ($(CUDA_ARCH),native)
NVCC_ARCH = -arch=native
else
NVCC_ARCH = $(foreach a,$(CUDA_ARCH),-gencode arch=compute_$(a),code=sm_$(a)) \
            -gencode arch=compute_$(lastword $(CUDA_ARCH)),code=compute_$(lastword $(CUDA_ARCH))
endif
NVCCFLAGS = -std=c++17 -O3 --default-stream per-thread $(NVCC_ARCH) $(if $(NVCC_CCBIN),-ccbin $(NVCC_CCBIN)) \
            -Xcompiler "$(OPT) -fopenmp"
GPU_LIBS = -L$(CUDA_LIBDIR) -lcudart_static -ldl -lrt -lpthread

SOURCES = NFS_factory/Batch_smooth.cpp NFS_factory/Fft_product.cpp NFS_factory/Main.cpp NFS_factory/Polynomial.cpp \
          NFS_factory/Prime_product.cpp NFS_factory/Profile.cpp NFS_factory/Smooth_factor.cpp NFS_factory/Tree_backend.cpp \
          NFS_factory/cofactorize_siqs.c int_ntt/int_ntt.cpp int_ntt/ntt_core.cpp
CUDA_SOURCES = NFS_factory/Gpu_tree_backend.cu int_ntt/ntt_gpu.cu

objects = $(patsubst %,build/$(1)/%.o,$(basename $(2)))
CPU_OBJECTS = $(call objects,cpu,$(SOURCES))
GPU_OBJECTS = $(call objects,gpu,$(SOURCES) $(CUDA_SOURCES))
UTILS = build/RelationsToBinaryFormat build/RelationFilterSplit
UTIL_OBJECTS = build/utils/RelationsToBinaryFormat.o build/utils/RelationFilterSplit.o

.PHONY: all cpu gpu utils both clean
all: cpu utils
cpu: build/NFS_factory
gpu: build/NFS_factory_gpu
utils: $(UTILS)
both: cpu gpu utils

build/NFS_factory: $(CPU_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

build/NFS_factory_gpu: $(GPU_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(GPU_LIBS)

build/RelationsToBinaryFormat: build/utils/RelationsToBinaryFormat.o
	$(CXX) $(CXXFLAGS) -o $@ $^

build/RelationFilterSplit: build/utils/RelationFilterSplit.o build/cpu/NFS_factory/Polynomial.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

build/utils/RelationsToBinaryFormat.o: Utilities/RelationsToBinaryFormat.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

build/utils/RelationFilterSplit.cpp: Utilities/RelationFilterSplit.cpp
	@mkdir -p $(@D)
	iconv -f UTF-16 -t UTF-8 $< > $@

build/utils/RelationFilterSplit.o: build/utils/RelationFilterSplit.cpp
	$(CXX) $(CPPFLAGS) -INFS_factory $(CXXFLAGS) -c -o $@ $<

build/cpu/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

build/cpu/%.o: %.c
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

build/gpu/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -DNFS_GPU $(CXXFLAGS) -c -o $@ $<

build/gpu/%.o: %.c
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) -DNFS_GPU $(CFLAGS) -c -o $@ $<

build/gpu/%.o: %.cu
	@mkdir -p $(@D)
	$(NVCC) $(INCLUDES) -DNFS_GPU $(NVCCFLAGS) -MMD -MP -c -o $@ $<

clean:
	rm -rf build

-include $(CPU_OBJECTS:.o=.d) $(GPU_OBJECTS:.o=.d) $(UTIL_OBJECTS:.o=.d)
