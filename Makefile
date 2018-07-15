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


#
# Patch
#

patch: build/cutest-patch-stamp

build/cutest-patch-stamp:
	git -C $(CUTEST) clean -f -d -x
	git -C $(CUTEST) reset --hard
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
	python


.PHONY: patch build

