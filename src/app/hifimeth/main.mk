ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ./$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ./$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := hifimeth
SOURCES  := \
	cov_to_bed.cpp \
	eval_kmer_features.cpp \
	eval.cpp \
	main.cpp \
	mod_batch.cpp \
	mod_main.cpp \
	mod_options.cpp \
	pileup_correlation.cpp \
	pileup.cpp \
	program_info.cpp \
	subsample_bam.cpp

SRC_INCDIRS  :=

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhifimeth
TGT_PREREQS := libhifimeth.a

SUBMAKEFILES :=
