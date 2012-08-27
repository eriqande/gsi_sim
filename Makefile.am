AM_CFLAGS =  -I$(top_srcdir)/eca-shared/ecalibs \
	 -I$(top_srcdir)/eca-shared/ranlib/src

bin_PROGRAMS = gsisim

gsisim_SOURCES = $(top_srcdir)/eca-shared/ranlib/linpack/linpack.c \
    $(top_srcdir)/eca-shared/ranlib/src/com.c \
    $(top_srcdir)/eca-shared/ranlib/src/ranlib.c \
    $(top_srcdir)/eca-shared/ecalibs/ECA_MemAlloc.c \
    $(top_srcdir)/eca-shared/ecalibs/ECA_Opt3.c \
    $(top_srcdir)/eca-shared/ecalibs/ECA_Opt3.h \
    $(top_srcdir)/eca-shared/ecalibs/ECA_ReadData.c \
    $(top_srcdir)/eca-shared/ecalibs/ECA_utilities.c \
    $(top_srcdir)/eca-shared/ecalibs/MathStatRand.c \
    $(top_srcdir)/eca-shared/ecalibs/MCTypesEtc.c \
    $(top_srcdir)/src/gsi_sim.c


AM_LDFLAGS = -lm
