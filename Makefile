#----------------------------------------------------------------
#   Parent Makefile for spline_ryab code
#   This file is just a wrapper for the various sub-makes
#
#   See bin/Makefile for the main Makefile
#
#   (c) 2025 bakl
#----------------------------------------------------------------
# ifeq ($(STLHOME),)
# # STLHOME := $(abspath $(shell pwd))
#  STLHOME := .
#  export STLHOME
# endif

.PHONY: all

all:
	@cd bin; ${MAKE} ${MAKECMDGOALS}

%::
	@cd bin; ${MAKE} "${MAKECMDGOALS}"

clean:
	@cd bin; ${MAKE} clean