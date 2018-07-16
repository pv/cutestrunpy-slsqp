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


all: diff

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

#
# Virtualenv
#

export NPY_NUM_BUILD_JOBS=4

env:
	$(PYTHON) -mvirtualenv env
	./env/bin/python -mpip install numpy scipy psutil

env-dev:
	$(PYTHON) -mvirtualenv env-dev
	./env-dev/bin/python -mpip install numpy Cython Tempita psutil
	./env-dev/bin/python -mpip install git+https://github.com/scipy/scipy@refs/pull/8986/head

run-installed: build env
	make -s run "PYTHON=$(CURDIR)/env/bin/python" "SCRIPT=cutest_slsqp.py" PYCUTEST_CACHE="$(CURDIR)/cache/installed" 2>&1|tee run-installed.log

run-dev: build env-dev
	make -s run "PYTHON=$(CURDIR)/env-dev/bin/python" "SCRIPT=cutest_slsqp.py" PYCUTEST_CACHE="$(CURDIR)/cache/dev" 2>&1|tee run-dev.log


diff: run-installed run-dev
	diff -u run-installed.log run-dev.log || true


.PHONY: patch build run run-installed run-dev all diff

