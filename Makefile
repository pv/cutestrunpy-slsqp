all: archdefs

export ARCHDEFS="$(CURDIR)/archdefs"
export SIFDECODE="$(CURDIR)/sifdecode"
export CUTEST="$(CURDIR)/cutest"


cutest-build-stamp: cutest archdefs
	cd cutest && bash ../archdefs/install_optrove
	touch cutest-build-stamp

#
# Patch
#

cutest-build-stamp: cutest
	install -d $(CUTEST)/src/python
	cp -f pycutestitf.c $(CUTEST)/src/tools
	cp -f pycutestmgr.py $(CUTEST)/src/python
