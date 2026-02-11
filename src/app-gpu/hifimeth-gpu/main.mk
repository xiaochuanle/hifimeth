ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ./$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ./$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := hifimeth-gpu
SOURCES  := \
	5mc_call_gpu.cpp \
	eval_kmer_features.cpp \
	main.cpp \
	mod_options.cpp \
	program_info.cpp \

SRC_INCDIRS  := . ../../../3rdparty/libtorch/include ../../../3rdparty/libtorch/include/torch/csrc/api/include

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhifimeth
TGT_PREREQS := libhifimeth.a

SUBMAKEFILES :=
