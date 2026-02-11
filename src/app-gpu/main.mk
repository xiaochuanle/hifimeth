ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ./$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ./$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libhifimeth.a

SOURCES      := \
	../corelib/5mc_context.cpp \
	../corelib/5mc_motif_finder.cpp \
	../corelib/bam_info.cpp \
	../corelib/bam_mod_parser.cpp \
	../corelib/get_core_count.cpp \
	../corelib/getRSS.c \
	../corelib/getMemorySize.c \
	../corelib/build_mod_bam.cpp \
	../corelib/hbn_aux.cpp \
	../corelib/hbn_seqdb.cpp \
	../corelib/line_reader.cpp \
	../corelib/spookyhash.c

SRC_INCDIRS  := 

SUBMAKEFILES := ./hifimeth-gpu/main.mk
