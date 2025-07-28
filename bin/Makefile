# Disable all of make's built-in rules (similar to Fortran's implicit none)
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
##############################
# MAKE FOR spline_ryab
##############################

.KEEP_STATE:

# Do you want to compile with debug options? 1=Yes, 0=No
IS_DEBUG = 0

#sys := $(shell uname)

# SYSTYPE ?= "ifort"
# SYSTYPE := "lf"
SYSTYPE ?= gfortran

MD := mkdir
MV := mv
DEL = rm -f


########################################################################
# Paths

ifeq ($(SPLHOME),)
 PWD       := $(abspath $(shell pwd))
#  PWD       = .
 SPLHOME   := $(PWD)/..
#  SPLHOME := $(abspath $(shell pwd))
 export SPLHOME
endif
SRCDIR    := $(SPLHOME)/src
BINDIR    := $(SPLHOME)/bin
RUNDIR    := $(SPLHOME)/test

LIBDIR    := $(SPLHOME)/lib


vpath %.f90 $(SRCDIR) $(SPLHOME)/test
vpath %.o   $(BINDIR)

VPATH := $(SRCDIR) $(SRCDIR)/test $(BINDIR) $(LIBDIR)
#INCL_DIR := -I$(SRCDIR)/common  -I$(SRCDIR)/stl -I$(SRCDIR)/util -I$(SRCDIR)/vladsf 
INCL_DIR := -I$(BINDIR)  -I$(LIBDIR)


########################################
## For Intel fortran
ifeq ($(SYSTYPE),"ifort")
  #  IS_STL_HDF5 = 1
  FC = ifort  # Intel fortran
  #   FFLAGS = -c -save -zero -O3 -mcmodel=medium -shared-intel
  #   FFLAGS = -c -save -zero -O3 -fp-model precise
  #   FFLAGS = -c -save -zero -O3 -fp-model strict -CB -g -traceback # -- for debug
  #   FFLAGS = -c -save -zero -fpe0 -no-ftz  -check bounds -g -traceback # -- for debug
  #   FFLAGS = -c -save -zero -O3 -fp-model strict
  #   FFLAGS = -c -O3 -save -zero -tpp7   -ip # -xW # -- optimize for Pentium
  #   FFLAGS =   -c -save -static -zero -O2 -fpe0 -traceback # -CB -traceback
  FFLAGS_FIX = -c -132 -save -zero -O3 -fp-model strict 
  # FFLAGS_FIX += -fpp   # preprocessor
  FFLAGS_FIX += -cpp   # preprocessor
  #  FFLAGS_FIX = -c -132 -save -zero -fp-model strict -CB -g -traceback -debug inline-debug-info -fpe0 -DD # -- for debug
  #  FFLAGS_FIX = -c -132 -save -zero -fp-model strict -CB -g -fpe0 -traceback -debug # inline-debug-info # -- for debug
  FFLAGS := -free
  FFLAGS := $(FFLAGS_FIX) $(FFLAGS)
  #define
  FFLAGS := $(FFLAGS) -D_INTEL

  
  LIBS= 
endif


########################################
## For gfortran fortran
ifeq ($(SYSTYPE),gfortran)
  FC = gfortran
  #  FFLAGS_FIX = -g -fbounds-check  -fno-automatic
  #  FFLAGS_FIX = -g  -fno-automatic
  #  FFLAGS_FIX = -O2  -ffixed-line-length-132 -fno-automatic -Wall
  ifeq ($(IS_DEBUG),1)
    FFLAGS_FIX += -g -O0 -fbounds-check -Wuninitialized -Wall
    FFLAGS_FIX += -ffpe-trap=invalid,zero,overflow -finit-real=nan
    # FFLAGS_FIX += -g  #  debugging
    # FFLAGS_FIX += -ffpe-trap=invalid,zero,overflow  #  debugging
    # FFLAGS_FIX += -finit-real=nan  #  debugging
    FFLAGS_FIX += -Wall  #  debugging
    FFLAGS_FIX += -fcheck=all  #  debugging
    # FFLAGS_FIX += -Warray-temporaries
    FFLAGS_FIX += -Wextra   -Wconversion #  debugging
    FFLAGS_FIX += -pedantic #  debugging
    FFLAGS_FIX += -fbacktrace
    FFLAGS_FIX += -frange-check
    FFLAGS_FIX += -fdump-core
    FFLAGS_FIX += -Wno-do-subscript  -Wno-surprising #  tmp
    FFLAGS_FIX += -Wno-unused-dummy-argument #  tmp
    FFLAGS_FIX += -Wno-maybe-uninitialized # gcc bug 
  else
    # FFLAGS_FIX += -O2  # optimize 
    FFLAGS_FIX += -O3  # optimize 
  endif
  # FFLAGS_FIX =  -fbounds-check   #  -fno-range-check 
  FFLAGS_FIX += -ffixed-line-length-none
  FFLAGS_FIX += -ffree-line-length-none #-ffixed-line-length-132 
  # FFLAGS_FIX += -fdefault-real-8 
  FFLAGS_FIX += -cpp   # preprocessor
  # FFLAGS_FIX += -O2   # optimize 
  # FFLAGS_FIX += -O0   # valgrind --track-origins=yes --keep-stacktraces=alloc-and-fre --leak-check=full 
  
  FFLAGS_FIX += -finit-local-zero
  # FFLAGS_FIX += -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan  #  debugging
  # FFLAGS_FIX += -Wno-function-elimination   # tmp https://gcc.gnu.org/onlinedocs/gfortran/Evaluation-of-logical-expressions.html
	FFLAGS_FIX += -ftrapping-math
  FFLAGS_FIX += -D_GFORTRAN
  # FFLAGS_FIX += -fopenmp  

  # FFLAGS_FIX += -fmax-stack-var-size=4500000
  
  # FFLAGS_FIX += -std=f2008   #  debugging
  # FFLAGS := -ffree-form
  # FFLAGS := $(FFLAGS_FIX) $(FFLAGS)
  FFLAGS := $(FFLAGS_FIX) $(FFLAGS) -ffree-form

  LDFLAGS := -o 
  LIBS := 
  
  
endif
########################################


#--------------------------
INCL =  

OBJDIR = obj

#--------------------------------------------
FILES_ay3d = kinds.f90 ryabenkiy.f90  main.f90
TEMP := $(patsubst %.f90,%.o, $(FILES_ay3d))
OBJ_ay3d :=  $(TEMP:.f=.o)
PROG_ay3d = xtest_ay3d.exe 

#--------------------------------------------
FILES_loaddata =  test_loaddata.f90
TEMP := $(patsubst %.f90,%.o, $(FILES_loaddata))
OBJloaddata :=  $(TEMP:.f=.o)
PROGloaddata = xtest_loaddata.exe 

#--------------------------------------------
FILES_ryabenky_spline =  ryabenky_spline.f90 test_ryabenky_nonuniform.f90
TEMP := $(patsubst %.f90,%.o, $(FILES_ryabenky_spline))
OBJryabenky_spline :=  $(TEMP:.f=.o)
PROGryabenky_spline = xtest_ryabenky_spline.exe 



.IGNORE:

###################################################


%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@  $(INCL_DIR)

%.o: %.mod ## if make wants m2c thinking that mod is from Modula language

%.mod : %.o
	@if [! -f $@ ]; then \
	@rm $< \
	$(MAKE) $< \
	fi

.PHONY:

all: ay3d
# 	@echo $(LIBS)


ay3d: $(OBJ_ay3d)  ## compile spline by AV Yudin for 3d
	$(FC) $(LDFLAGS) $(PROG_ay3d) $(OBJ_ay3d) $(LIBS)
	$(MV) $(PROG_ay3d)  $(RUNDIR)     
	@ls -l $(RUNDIR)/$(PROG_ay3d)



spline: $(OBJryabenky_spline)  ## compile spline lib and test
	$(FC) $(LDFLAGS) $(PROGryabenky_spline) $(OBJryabenky_spline) $(LIBS)
	$(MV) $(PROGryabenky_spline)  $(RUNDIR)     
	@ls -l $(RUNDIR)/$(PROGryabenky_spline)


load: $(OBJloaddata)  ## load dump data 
	$(FC) $(LDFLAGS) $(PROGloaddata) $(OBJloaddata) $(LIBS)
	$(MV) $(PROGloaddata)  $(RUNDIR)     
	@ls -l $(RUNDIR)/$(PROGloaddata)


help: ## Show this dynamic help
	$(info Available targets)
	@awk '/^[a-zA-Z\-_0-9]+:/ {                                   \
      nb = sub( /^## /, "", helpMsg );                             \
      if(nb == 0) {                                                \
        helpMsg = $$0;                                             \
        nb = sub( /^[^:]*:.* ## /, "", helpMsg );                  \
      }                                                            \
      if (nb)                                                      \
        printf "\033[1;31m%-" width "s\033[0m %s\n", $$1, helpMsg; \
                                                                   \
      helpMsg = $$0 }'                                             \
    width=$$(grep -o '^[a-zA-Z_0-9]\+:' $(MAKEFILE_LIST) |  wc -L)  \
    $(MAKEFILE_LIST)


.PHONY:

clean:  ## rm all temporry files
	$(DEL) $(BINDIR)/*.o 
	$(DEL) $(SRCDIR)/*.mod $(SRCDIR)/*.o
	find $(BINDIR)/  -maxdepth 1 \( -not -type l -name "*.mod" \)  -print0 | xargs -0  -I {} rm {}

cleandata:  ## rm ../run/strad/nohup.out ../run/strad/core* ../run/strad/fort.*
	$(DEL) $(RUNDIR)/*.dump

cleanall:  clean  ## rm all temporry files and executable files
	$(DEL) $(RUNDIR)/$(PROGloaddata) $(RUNDIR)/$(PROGryabenky_spline)  $(RUNDIR)/$(PROG_ay3d)
