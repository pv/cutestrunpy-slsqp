# (C)2011 Arpad Buermen
# Licensed under LGPL V2.1

"""PyCUTEst problem manager

Currently works only under Linux.

Requres CUTEst installation for double precision built with gcc and gfortran. 
CUTEst library must be built with -fPIC (as position independent code). 

SIF decoder binary (sifdecode) and gfortran must be accessible through $PATH. 

PYCUTEST_CACHE environmental variable specifies the location of problem 
interface cache. If it is not specified the current directory (.) is used 
for caching. Problems are cached in $PYCUTEST_CACHE/pycutest. 

Once a problem is built only the following files in 
$PYCUTEST_CACHE/pycutest/PROBLEM_NAME are needed for the problem to work: 
  _pycutestitf.so	-- CUTEst problem interface module
  __init__.py       -- module initialization and wrappers
  OUTSDIF.d         -- problem information

PYTHONPATH must include $PYCUTEST_CACHE and $CUTEST/src/python. 

A cached problem can be imported with 'from pycutest import PROBLEM_NAME'. One 
can also use importProblem() which returns a reference to the problem module. 

Available functions
clearCache     -- remove Python interface to problem from cache
prepareProblem -- decode problem and build Python interface
importProblem  -- import problem interface module from cache
                  (prepareProblem must be called first)
"""

import os, shutil, sys
from os import environ, getcwd
import subprocess
from glob import glob

__all__ = [ 'clearCache', 'prepareProblem', 'importProblem', 'isCached' ]

#
# The sutup.py script with a placeholder for platform-dependent part. 
#
setupScript="""#!/usr/bin/env python
# (C)2011 Arpad Buermen
# Licensed under LGPL V2.1

#
# Do not edit. This is a computer-generated file. 
#

from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc
import os
from subprocess import call
from glob import glob
import numpy, numpy.distutils

#
# OS specific
#

%s

#
# End of OS specific
#

interface_source_file=os.path.join(os.environ['PYCUTEST'], 'pycutestitf.c')

# Settings
setup(name='PyCuter automatic test function interface builder',
	version='1.0',
	description='Builds a CUTEst test function interface for Python.',
	long_description='Builds a CUTEst test function interface for Python.', 
	author='Arpad Buermen',
	author_email='arpadb@fides.fe.uni-lj.si',
	url='',
	platforms='Linux', 
	license='LGPL V2.1', 
	packages=[],
	ext_modules=[
		Extension(
			'_pycutestitf', 
			[interface_source_file], 
			include_dirs=include_dirs, 
			define_macros=define_macros, 
			extra_objects=objFileList, 
			libraries=libraries,
		)
	]
)
"""

#
# Linux-specific part of setup.py
#
setupScriptLinux="""
define_macros=[('LINUX', None)]
include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs()
include_dirs.append(os.path.join(os.environ['CUTEST'], 'include'))
include_dirs.append(os.path.join(include_dirs[0], 'numpy'))
objFileList=glob('*.o')
libraries=['gfortran', 'lapack']
objFileList.append(os.path.join(os.environ['CUTEST'], 'objects', os.environ['MYARCH'], 'double', 'libcutest.a'))
"""

#
# Problem interface module initialization file with placeholders for 
# efirst, lfirst, and nvfirst. This also defines the wrapper functions. 
# A placeholder is included for problem name and ordering. 
#
initScript="""# PyCUTEst problem interface module intialization file
# (C)2011 Arpad Buermen
# Licensed under LGPL V2.1

\"\"\"Interface module for CUTEst problem %s with ordering 
  efirst=%s, lfirst=%s, nvfirst=%s
sifdecode parameters : %s
sifdecode options    : %s

Available functions
getinfo    -- get problem information
varnames   -- get names of problem's variables
connames   -- get names of problem's constraints
objcons    -- objective and constraints
obj        -- objective and objective gradient
cons       -- constraints and constraints gradients/Jacobian
lagjac     -- gradient of objective/Lagrangian and constraints Jacobian 
jprod      -- product of constraints Jacobian with a vector
hess       -- Hessian of objective/Lagrangian
ihess      -- Hessian of objective/constraint
hprod      -- product of Hessian of objective/Lagrangian with a vector
gradhess   -- gradient and Hessian of objective (unconstrained problems) or
              gradient of objective/Lagrangian, Jacobian of constraints and 
              Hessian of Lagrangian (constrained problems)
scons      -- constraints and sparse Jacobian of constraints
slagjac    -- gradient of objective/Lagrangian and sparse Jacobian
sphess     -- sparse Hessian of objective/Lagrangian
isphess    -- sparse Hessian of objective/constraint
gradsphess -- gradient and sparse Hessian of objective (unconstrained probl.) 
              or gradient of objective/Lagrangian, sparse Jacobian of 
              constraints and sparse Hessian of Lagrangian (constrained probl.)
report     -- get usage statistics
\"\"\"

from _pycutestitf import *
import _pycutestitf
import os
from scipy.sparse import coo_matrix
from numpy import zeros

# Get the directory where the bunary module (and OUTSDIF.d) are found. 
(_directory, _module)=os.path.split(_pycutestitf.__file__)

# Problem info structure and dimension
info=None
n=None
m=None

# Constraints and variable ordering
efirst=%s
lfirst=%s
nvfirst=%s

# Remember current directory and go to module directory where OUTSDIF.d is located
fromDir=os.getcwd()
os.chdir(_directory)
	
# Get problem dimension
(n, m)=_pycutestitf._dims()
	
# Set up the problem and get basic information
info=_pycutestitf._setup(efirst, lfirst, nvfirst)
	
# Store constraint and variable ordering information 
if m>0:
	info['efirst']=efirst
	info['lfirst']=lfirst
info['nvfirst']=nvfirst 

# Store sifdecode parameters and options
info['sifparams']=%s
info['sifoptions']=%s

# Go back to initial directory
os.chdir(fromDir)

# Return problem info
def getinfo():
	\"\"\"
	Return the problem info dictionary.
	
	info=geinfo()
	
	Output 
	info -- dictionary with the summary of test function's properties
	
	The dictionary has the following members:
	name       -- problem name
	n          -- number of variables
	m          -- number of constraints (excluding bounds)
	x          -- initial point (1D array of length n)
	bl         -- 1D array of length n with lower bounds on variables 
	bu         -- 1D array of length n with upper bounds on variables
	nnzh       -- number of nonzero elements in the diagonal and upper triangle of
	              sparse Hessian
	vartype    -- 1D integer array of length n storing variable type
	              0=real,  1=boolean (0 or 1), 2=integer
	nvfirst    -- boolean flag indicating that nonlinear variables were placed
	              before linear variables
	sifparams  -- parameters passed to sifdecode with the -param option 
	              None if no parameters were given
	sifoptions -- additional options passed to sifdecode
	              None if no additional options were given. 
			   
	For constrained problems the following additional members are available
	nnzj    -- number of nonzero elements in sparse Jacobian of constraints
	v       -- 1D array of length m with initial values of Lagrange multipliers
	cl      -- 1D array of length m with lower bounds on constraint functions
	cu      -- 1D array of length m with upper bounds on constraint functions
	equatn  -- 1D boolean array of length m indicating whether a constraint
	           is an equation constraint
	linear  -- 1D boolean array of length m indicating whether a constraint
	           is a linear constraint
	efirst  -- boolean flag indicating that equation constraints were places
	           before inequation constraints
	lfirst  -- boolean flag indicating that linear constraints were placed 
	           before nonlinear constraints
	\"\"\"
	return info

def varnames():
	\"\"\"
	Return the names of problem's variables.
	
	nameList=varnames()
	
	nameList -- a list of strings representing the names of problem's variables.
	            The variabels are ordered according to nvfirst flag. 
	\"\"\"
	return _pycutestitf._varnames()
	
def connames():
	\"\"\"
	Return the names of problem's constraints.
	
	nameList=connames()
	
	nameList -- a list of strings representing the names of problem constraints. 
	            The constraints are ordered according to efirst and lfirst flags. 
	\"\"\"
	return _pycutestitf._connames()
	
# Sparse tool wrappers (return scipy.sparse.coo_matrix matrices)
# _scons() wrapper
def scons(x, i=None):
	\"\"\"Returns the value of constraints and 
	the sparse Jacobian of constraints at x.
	
	(c, J)=_scons(x)      -- Jacobian of constraints
	(ci, gi)=_scons(x, i) -- i-th constraint and its gradient
	
	Input
	x -- 1D array of length n with the values of variables
	i -- integer index of constraint (between 0 and m-1)
	
	Output
	c  -- 1D array of length m holding the values of constraints at x
	J  -- a scipy.sparse.coo_matrix of size m-by-n holding the Jacobian at x
	ci -- 1D array of length 1 holding the value of i-th constraint at x
	gi -- a scipy.sparse.coo_matrix of size 1-by-n holding the gradient of i-th constraint at x
	
	This function is a wrapper for _scons().
	\"\"\"
	
	if i is None:
		(c, Ji, Jif, Jv)=_pycutestitf._scons(x)
		return (c, coo_matrix((Jv, (Jif, Ji)), shape=(m, n)))
	else:
		(c, gi, gv)=_pycutestitf._scons(x, i)
		return (c, coo_matrix((gv, (zeros(n), gi)), shape=(1, n)))

# _slagjac() wrapper
def slagjac(x, v=None):
	\"\"\"Returns the sparse gradient of objective at x or Lagrangian at (x, v),
	and the sparse Jacobian of constraints at x.
	
	(g, J)=_slagjac(x)    -- objective gradient and Jacobian
	(g, J)=_slagjac(x, v) -- Lagrangian gradient and Jacobian
	
	Input
	x -- 1D array of length n with the values of variables
	v -- 1D array of length m with the values of Lagrange multipliers
	
	Output
	g -- a scipy.sparse.coo_matrix of size 1-by-n holding the gradient of objective at x or
	     the gradient of Lagrangian at (x, v)
	J -- a scipy.sparse.coo_matrix of size m-by-n holding the sparse Jacobian
	     of constraints at x
	
	This function is a wrapper for _slagjac().
	\"\"\"
	
	if v is None:
		(gi, gv, Ji, Jfi, Jv)=_pycutestitf._slagjac(x)
	else:
		(gi, gv, Ji, Jfi, Jv)=_pycutestitf._slagjac(x, v)
	return (
		coo_matrix((gv, (zeros(n), gi)), shape=(1, n)), 
		coo_matrix((Jv, (Jfi, Ji)), shape=(m, n))
	)

# _sphess() wrapper
def sphess(x, v=None):
	\"\"\"Returns the sparse Hessian of the objective at x (unconstrained problems) 
	or the sparse Hessian of the Lagrangian (constrained problems) at (x, v).
	
	H=_sphess(x)    -- Hessian of objective (unconstrained problems)
	H=_sphess(x, v) -- Hessian of Lagrangian (constrained problems)
	
	Input
	x -- 1D array of length n with the values of variables
	v -- 1D array of length m with the values of Lagrange multipliers
	
	Output
	H -- a scipy.sparse.coo_matrix of size n-by-n holding the sparse Hessian
	     of objective at x or the sparse Hessian of the Lagrangian at (x, v)
	
	This function is a wrapper for _sphess().
	\"\"\"
	
	if v is None:
		(Hi, Hj, Hv)=_pycutestitf._sphess(x)
	else:
		(Hi, Hj, Hv)=_pycutestitf._sphess(x, v)
	return coo_matrix((Hv, (Hi, Hj)), shape=(n, n))

# _isphess() wrapper
def isphess(x, i=None):
	\"\"\"Returns the sparse Hessian of the objective or the sparse Hessian of i-th
	constraint at x.
	
	H=_isphess(x)    -- Hessian of objective 
	H=_isphess(x, i) -- Hessian of i-th constraint
	
	Input
	x -- 1D array of length n with the values of variables
	i -- integer holding the index of constraint (between 0 and m-1)
	
	Output
	H -- a scipy.sparse.coo_matrix of size n-by-n holding the sparse Hessian
	     of objective or the sparse Hessian i-th constraint at x 
	
	This function is a wrapper for _isphess().
	\"\"\"
	
	if i is None:
		(Hi, Hj, Hv)=_pycutestitf._isphess(x)
	else:
		(Hi, Hj, Hv)=_pycutestitf._isphess(x, i)
	return coo_matrix((Hv, (Hi, Hj)), shape=(n, n))

# _gradsphess() wrapper
def gradsphess(x, v=None, lagrFlag=False):
	\"\"\"Returns the sparse Hessian of the Lagrangian, the sparse Jacobian of
	constraints, and the gradient of the objective or Lagrangian.
	
	(g, H)=gradsphess(x)              -- unconstrained problems
	(g, J, H)=gradsphess(x, v, gradl) -- constrained problems
	
	Input
	x     -- 1D array of length n with the values of variables
	v     -- 1D array of length m with the values of Lagrange multipliers
	gradl -- boolean flag. If False the gradient of the objective is returned, 
	         if True the gradient of the Lagrangian is returned. 
	         Default is False
			 
	Output
	g -- a scipy.sparse.coo_matrix of size 1-by-n holding the gradient of objective at x or 
	     the gradient of Lagrangian at (x, v)
	J -- a scipy.sparse.coo_matrix of size m-by-n holding the sparse Jacobian
	     of constraints at x
	H -- a scipy.sparse.coo_matrix of size n-by-n holding the sparse Hessian
	     of objective at x or the sparse Hessian of the Lagrangian at (x, v) 
	
	This function is a wrapper for _gradsphess().
	\"\"\"
	
	if v is None:
		(g, Hi, Hj, Hv)=_pycutestitf._gradsphess(x)
		return (coo_matrix(g), coo_matrix((Hv, (Hi, Hj)), shape=(n, n)))
	else:
		(gi, gv, Ji, Jfi, Jv, Hi, Hj, Hv)=_pycutestitf._gradsphess(x, v, lagrFlag)
		return (
			coo_matrix((gv, (zeros(n), gi)), shape=(1, n)), 
			coo_matrix((Jv, (Jfi, Ji)), shape=(m, n)), 
			coo_matrix((Hv, (Hi, Hj)), shape=(n, n))
		)

# Clean up
del os, fromDir, efirst, lfirst, nvfirst
""" 

def _cachePath():
	"""Return the path to PyCUTEst cache (PYCUTEST_CACHE environmental variable). 
	If PYCUTEST_CACHE is not set, return the full path to current work directory. 
	"""
	
	if 'PYCUTEST_CACHE' in environ:
		return environ['PYCUTEST_CACHE']
	else:
		return os.getcwd()

def isCached(cachedName):
	"""
	Return True if a problem is in cache.
	
	Keyword arguments:
	cachedName -- cache entry name
	"""
	
	# The problem's cache entry 
	problemDir=os.path.join(cachePath, 'pycutest', cachedName)
	
	# See if a directory with problem's name exists
	return os.path.isdir(problemDir)
	
def clearCache(cachedName):
	"""
	Removes a cache entry from cache. 
	
	Keyword arguments:
	cachedName -- cache entry name
	""" 
	
	# The problem's cache entry 
	problemDir=os.path.join(cachePath, 'pycutest', cachedName)
	
	# See if a directory with problem's name exists
	if os.path.isdir(problemDir):
		# It exists, delete it. 
		shutil.rmtree(problemDir, True)
	elif os.path.isfile(problemDir):
		# It is a file, delete it. 
		os.remove(problemDir)
	
def prepareCache(cachedName):
	"""
	Prepares a cache entry. 
	If an entry already exists it is deleted first. 
	
	Keyword arguments:
	cachedName -- cache entry name
	"""
	
	# The directory with test function entries
	pycutestDir=os.path.join(cachePath, 'pycutest')
	
	# The problem's cache entry 
	problemDir=os.path.join(pycutestDir, cachedName)
	
	# See if a folder named pycutest exists in the cache path. 
	if not os.path.isdir(pycutestDir):
		# Create it. If this fails, give up. The user should delete manualy the
		# offending file which prevents the creation of a directory. 
		os.mkdir(pycutestDir)
	
	# See in pycutestDir if there is an __init__.py file. 
	initfile=os.path.join(pycutestDir, '__init__.py')
	if not os.path.isfile(initfile):
		# Create it
		f=open(initfile, 'w+')
		f.write("#PyCUTEst cache initialization file\n")
		f.close()
	
	# Remove old entry
	clearCache(cachedName)
	
	# Create folder with problem's name
	os.mkdir(problemDir)
	
def decodeAndCompileProblem(problemName, destination=None, sifParams=None, sifOptions=None, quiet=True):
	"""
	Call sifdecode on given problem and compile the resulting .f files. 
	Use gfortran with -fPIC option for compiling. 
	Collect the resulting object file names and return them. 
	This function is OS dependent. Currently works only for Linux. 
	
	Keyword arguments:
	problemName -- CUTEst problem name
	destination -- the name under which the compiled problem interface is stored in the cache
	               If not given problemName is used. 
	sifParams   -- parameters passed to sifdecode using the -param command line option 
	               given in the form of a dictionary with parameter name as key. Values 
				   are converted to strings using str() and every parameter contributes
				     -param key=str(value)
				   to the sifdecode's command line options. 
	sifOptions  -- additional options passed to sifdecode given in the form of a list of strings. 
	quiet       -- supress output (default True)
	
	*destination* must not contain dots because it is a part of a Python module name. 
	"""
	
	# Default destination
	if destination is None:
		destination=problemName
	
	# The problem's cache entry 
	problemDir=os.path.join(cachePath, 'pycutest', destination)
	
	# Remember current work directory and go to cache
	fromDir=os.getcwd()
	os.chdir(problemDir)
	
	# Additional args
	args=[]
	
	# Handle params
	if sifParams is not None:
		for (key, value) in sifParams.iteritems():
			#if type(key) is not str:
			#	raise Exception, "sifParams keys must be strings"
			args+=['-param', key+"="+str(value)]
	
	# Handle options
	if sifOptions is not None: 
		for opt in sifOptions:
			#if type(opt) is not str:
			#	raise Exception, "sifOptions must consist of strings"
			args+=[str(opt)]
	
	# Call sifdecode (assume it is in PATH)
	spawnOK=True
	p=None
	try:
		# Start sifdecode
		p=subprocess.Popen(
				['sifdecoder']+args+[problemName], 
				universal_newlines=True, 
				stdout=subprocess.PIPE, stderr=subprocess.STDOUT
			)
		
		# Collect output
		messages=p.stdout.read()
		
		# Now wait for the process to finish. If we don't wait p might get garbage-collected before the
		# actual process finishes which can result in a crash of the interpreter. 
		retcode=p.wait()
		
		# Check return code. Nonzero return code means that something has gone bad. 
		if retcode!=0:
			spawnOK=False
	except:
		spawnOK=False
	
	if not spawnOK or not quiet:
		print(messages)
		
	# Collect all .f files
	filelist=glob('*.f')
	
	# Compile FORTRAN files 
	for filename in filelist:
		cmd=['gfortran', '-fPIC', '-c', filename]
		if not quiet:
			for s in cmd:
				print(s,) 
			print()
		if subprocess.call(cmd)!=0:
			raise Exception#, "gfortran call failed for "+filename
		
	# Collect list of all object files (.o)
	objFileList=glob('*.o')
	
	# Go back to original work directory
	os.chdir(fromDir)
	
	return objFileList
	
def compileAndInstallInterface(problemName, objFileList, destination=None, sifParams=None, sifOptions=None, 
								efirst=False, lfirst=False, nvfirst=False, quiet=True):
	"""
	Compiles and installs the binary interface module. 
	Uses distutils to achieve this. 
	Assumes decodeAndCompile() successfully completed its task. 
	This function is OS dependent. Currently works only for Linux. 
	
	Keyword arguments:
	problemName -- CUTEst problem name
	destination -- the name under which the compiled problem interface is stored in the cache
	               If not given problemName is used. 
	objFileList -- list of object files that were generated using gfortran
	sifParams   -- parameters passed to sifdecode using the -param command line option 
	               given in the form of a dictionary with parameter name as key. Values 
				   are converted to strings using str() and every parameter contributes
				     -param key=str(value)
				   to the sifdecode's command line options. 
	sifOptions  -- additional options passed to sifdecode given in the form of a list of strings. 
	efirst      -- order equation constraints first (default True)
	lfirst      -- order linear constraints first (default True)
	nvfirst     -- order nonlinear variables before linear variables 
	               (default False)
	quiet       -- supress output (default True)
	
	*destination* must not contain dots because it is a part of a Python module name. 
	"""
	
	# Default destination
	if destination is None:
		destination=problemName
	
	# The problem's cache entry 
	problemDir=os.path.join(cachePath, 'pycutest', destination)
	
	# Remember current work directory and go to cache
	fromDir=os.getcwd()
	os.chdir(problemDir)
	
	# Prepare a script file 
	f=open('setup.py', 'w+')
	f.write(setupScript % setupScriptLinux)
	f.close()
	
	# Convert sifParams to a string
	sifParamsStr=""
	if sifParams is not None:
		for (key, value) in sifParams.iteritems():
			sifParamsStr+="%s=%s " % (str(key), str(value))
	
	# Convert sifOptions to a string
	sifOptionsStr=""
	if sifOptions is not None:
		for opt in sifOptions:
			sifOptionsStr+=str(opt)+" "
			
	# Prepare -q option for setup.py
	if quiet:
		quietopt=['-q']
	else:
		quietopt=[]
	
	# Call 'python setup.py build'
	if subprocess.call([sys.executable, 'setup.py']+quietopt+['build'])!=0:
		raise Exception#, "Failed to build the Python interface module"
	
	# Call 'python setup.py install --install-lib .'
	if subprocess.call([sys.executable, 'setup.py']+quietopt+['install', '--install-lib', '.'])!=0:
		raise Exception#, "Failed to install the Python interface module"
		
	# Create __init__.py
	f=open('__init__.py', 'w+')
	f.write(initScript % (
			problemName, 
			str(bool(efirst)), str(bool(lfirst)), str(bool(nvfirst)), 
			sifParamsStr, sifOptionsStr, 
			str(bool(efirst)), str(bool(lfirst)), str(bool(nvfirst)), 
			str(sifParams), str(sifOptions)
		)
	)
	f.close()
	
	# Go back to original work directory
	os.chdir(fromDir)

def prepareProblem(problemName, destination=None, sifParams=None, sifOptions=None, 
					efirst=False, lfirst=False, nvfirst=False, quiet=True):
	"""
	Prepares a problem interface module, imports and initiazlizes it, 
	and returns a reference to the imported module. 
	
	Keyword arguments:
	problemName -- CUTEst problem name
	destination -- the name under which the compiled problem interface is stored in the cache
	               If not given problemName is used. 
	sifParams   -- parameters passed to sifdecode using the -param command line option 
	               given in the form of a dictionary with parameter name as key. Values 
				   are converted to strings using str() and every parameter contributes
				     -param key=str(value)
				   to the sifdecode's command line options. 
	sifOptions  -- additional options passed to sifdecode given in the form of a list of strings. 
	efirst      -- order equation constraints first (default True)
	lfirst      -- order linear constraints first (default True)
	nvfirst     -- order nonlinear variables before linear variables 
	               (default False)
	quiet       -- supress output (default True)
	
	*destination* must not contain dots because it is a part of a Python module name. 
	"""
	
	# Default destination
	if destination is None:
		destination=problemName
	
	# Build it
	prepareCache(destination)
	objList=decodeAndCompileProblem(problemName, destination, sifParams, sifOptions, quiet)
	compileAndInstallInterface(problemName, objList, destination, sifParams, sifOptions, efirst, lfirst, nvfirst, quiet)
	
	# Import interface module. Initialization is done by __init__.py. 
	return importProblem(destination)
	
def importProblem(cachedName):
	"""
	Imports and initializes a problem module with CUTEst interface functions. 
	The module must be available in cache (see prepareProblem()). 
	
	Keyword arguments:
	cachedName -- name under which the problem is stored in cache
	"""
	
	# Import interface module. Initialization is done by __init__.py. 
	return __import__('pycutest.'+cachedName, globals(), locals(), [cachedName]) 
	
	
# Initialization (performed at first import)
# Save full path to PyCUTEst cache in cachePath. 
cachePath=_cachePath()
