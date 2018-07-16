all: archdefs

export ARCHDEFS = $(CURDIR)/archdefs
export SIFDECODE = $(CURDIR)/sifdecode
export CUTEST = $(CURDIR)/cutest
export PYCUTEST_CACHE = $(CURDIR)/cache

export PATH := $(PATH):$(SIFDECODE)/bin:$(CUTEST)/bin
export MANPATH := $(MANPATH):$(SIFDECODE)/man:$(CUTEST)/man
export MASTSIF = $(CURDIR)/sif
export MYARCH = pc.lnx.gfo

export PYCUTEST = $(CUTEST)/src/tools
export PYCUTEST_CACHE = $(CURDIR)/cache
export PYTHONPATH = $(PYCUTEST_CACHE):$(CUTEST)/src/python

PYTHON=python


#
# Patch
#

patch: build/cutest-patch-stamp

build/cutest-patch-stamp:
	for repo in "$(CUTEST)" "$(ARCHDEFS)" "$(SIFDECODE)"; do \
		git -C "$$repo" clean -f -d -x; \
		git -C "$$repo" reset --hard; \
	done
	install -d build
	install -d $(CUTEST)/src/python
	cp -f pycutestitf.c $(CUTEST)/src/tools
	cp -f pycutestmgr.py $(CUTEST)/src/python
	touch "$@"

#
# Build
#

build: build/cutest-build-stamp

build/cutest-build-stamp: build/cutest-patch-stamp
	install -d build
	cd cutest && bash ../archdefs/install_optrove
	touch "$@"

#
# Run
# 

run: build
	install -d $(PYCUTEST_CACHE)
	OPT="-O2 -w" $(PYTHON) $(SCRIPT)


.PHONY: patch build

