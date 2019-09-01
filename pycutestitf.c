/* CUTEst2 interface to Python and NumPy */
/* (c)2011 Arpad Buermen */
/* Licensed under LGPL V2.1 */

/* Note that in Windows we do not use Debug compile because we don't have the debug version
   of Python libraries and interpreter. We use Release version instead where optimizations
   are disabled. Such a Release version can be debugged. 
 */
 
/* Unused CUTEst tools - sparse finite element matrices and banded matrices
     cdimse
	 ceh
	 csgreh
	 ubandh
	 udimse
	 ueh
     ugreh
	 
   CUTEst tools that are not used because they duplicate functionality or are obsolete
     cnames		... used pbname, varnames, connames
	 udimen		... used cdimen 
	 ufn		... used uofg
	 ugr		... used uofg
	 unames		... used pbname, varnames
     cscfg		... obsolete
	 cscifg		... obsolete
*/
 
#include "Python.h"
#include "cutest.h"
#include "arrayobject.h"
#include <math.h>
#include <stdio.h>

/* Debug switch - uncomment to enable debug messages */
/* #undefine PYDEBUG */

/* Debug file */
#define df stdout

#ifndef WIN32
#define __declspec(a) 
#endif

/* Safeguard against C++ symbol mangling */
#ifdef __cplusplus
extern "C" {
#endif

/* NumPy type for CUTEst integer 
   INTEGER type in FORTRAN is by default 4 bytes wide (tested with gfortran 
   on AMD64 Linux). The C prototypes of FORTRAN-based CUTEst tools use 
   integer C type which is typedefed to long int in cutest.h. 
   On IA32 Linux long int is 4 bytes wide. Unfortunately on AMD64 it is 8 
   bytes wide. This causes problems if we are passing arguments to
   FORTRAN functions. A serious problem appears when an input or output 
   argument is an array because the elements of the array are twice the size
   the FORTRAN expects them to be. 
   Until this bug is fixed in CUTEst the solution is to use npy_integer for
   variables and npy_integer_type_num for creating NumPy arrays. To get rid 
   of warnings integer arguments are explicitly converted to (integer *) when 
   they are passed to FORTRAN functions. This is harmless since FORTRAN 
   arguments are passed by reference (as pointers) and the type conversion does 
   not affect the pointer's value. 
   The same goes for logical typedef of CUTEst. This interface uses npy_logical 
   which is actually int, but type casts npy_logical* pointers to logical*.
*/
static int npy_integer_type_num=PyArray_INT;
typedef int npy_integer; 
typedef int npy_logical; 


/* Prototypes */
static PyObject *cutest__dims(PyObject *self, PyObject *args);
static PyObject *cutest__setup(PyObject *self, PyObject *args);
static PyObject *cutest__varnames(PyObject *self, PyObject *args);
static PyObject *cutest__connames(PyObject *self, PyObject *args);
static PyObject *cutest_objcons(PyObject *self, PyObject *args);
static PyObject *cutest_obj(PyObject *self, PyObject *args);
static PyObject *cutest_cons(PyObject *self, PyObject *args);
static PyObject *cutest_lagjac(PyObject *self, PyObject *args);
static PyObject *cutest_jprod(PyObject *self, PyObject *args);
static PyObject *cutest_hess(PyObject *self, PyObject *args);
static PyObject *cutest_ihess(PyObject *self, PyObject *args);
static PyObject *cutest_hprod(PyObject *self, PyObject *args);
static PyObject *cutest_gradhess(PyObject *self, PyObject *args);
static PyObject *cutest__scons(PyObject *self, PyObject *args); 
static PyObject *cutest__slagjac(PyObject *self, PyObject *args); 
static PyObject *cutest__sphess(PyObject *self, PyObject *args); 
static PyObject *cutest__isphess(PyObject *self, PyObject *args); 
// static PyObject *cutest__gradsphess(PyObject *self, PyObject *args); 
// static PyObject *cutest_report(PyObject *self, PyObject *args); 

/* Persistent data */
#define STR_LEN 10
static npy_integer CUTEst_nvar = 0;		/* number of variables */
static npy_integer CUTEst_ncon = 0;		/* number of constraints */
static npy_integer CUTEst_nnzj = 0;		/* nnz in Jacobian */
static npy_integer CUTEst_nnzh = 0;		/* nnz in upper triangular Hessian */
static char CUTEst_probName[STR_LEN+1];		/* problem name */
static char setupCalled = 0;			/* Flag to indicate if setup was called */
static char dataFileOpen = 0;			/* Flag to indicate if OUTSDIf is open */

static npy_integer funit = 42;			/* FORTRAN unit number for OUTSDIF.d */
static npy_integer iout = 6;			/* FORTRAN unit number for error output */
static npy_integer io_buffer = 11;		/* FORTRAN unit internal input/output */
npy_integer status;				/* Exit flag from CUTEst tools */
static char  fName[] = "OUTSDIF.d"; 	/* Data file name */
/* Logical constants for FORTRAN calls */
static logical somethingFalse = FALSE_, somethingTrue = TRUE_;


/* Helper functions */

/* Open data file, return 0 on error. */
int open_datafile(void) {
	npy_integer  ioErr;					/* Exit flag from OPEN and CLOSE */
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: opening data file\n");
#endif
	ioErr = 0;
	if (! dataFileOpen) 
		FORTRAN_open((integer *)&funit, fName, (integer *)&ioErr);
	if (ioErr) {
		PyErr_SetString(PyExc_Exception, "Failed to open data file");
		return 0;
	}
	dataFileOpen = 1;
	return 1;
}

/* Close data file, return 0 on error. */
int close_datafile(void) {
	npy_integer ioErr;					/* Exit flag from OPEN and CLOSE */
	ioErr = 0;
	FORTRAN_close((integer *)&funit, (integer *)&ioErr);
	if (ioErr) {
		PyErr_SetString(PyExc_Exception, "Error closing data file");
		return 0;
	}
	dataFileOpen = 0;
	return 1; 
}

/* Check if the problem is set up, return 0 if it is not. */
int check_setup(void) {
	if (!setupCalled) {
		PyErr_SetString(PyExc_Exception, "Problem is not set up");
		return 0;
	}
	return 1;
}

/* Trim trailing spaces from a string starting at index n. */
void trim_string(char *s, int n) {
	int i;
	
	for(i=n;i>=0;i--) {
		if (s[i]!=' ')
			break;
	}
	s[i+1]=0;
}

/* Decrese reference counf for newly created dictionary members */
PyObject *decRefDict(PyObject *dict) {
	PyObject *key, *value;
	Py_ssize_t pos; 
	pos=0;
	while (PyDict_Next(dict, &pos, &key, &value)) {
		Py_XDECREF(value);
	}
	return dict;
}

/* Decrease reference count for newly created tuple members */
PyObject *decRefTuple(PyObject *tuple) {
	Py_ssize_t pos; 
	for(pos=0;pos<PyTuple_Size(tuple);pos++) {
		Py_XDECREF(PyTuple_GetItem(tuple, pos));
	}
	return tuple;
}

/* Extract sparse gradient and Jacobian in form of NumPy arrays */
void extract_sparse_gradient_jacobian(npy_integer nnzjplusno, npy_integer *sji, npy_integer *sjfi, npy_double *sjv, 
		PyObject **Mgi, PyObject **Mgv, PyObject **MJi, PyObject **MJfi, PyObject **MJv) {
	npy_double *gv, *Jv;
	npy_integer *gi, *Ji, *Jfi, nnzg, i, jg, jj; 
	npy_intp dims[1]; 
	
	/* Get number of nonzeros in gradient vector */
	nnzg=0;
	for(i=0;i<nnzjplusno;i++) {
		if (sjfi[i]==0)
			nnzg++;
	}
	
	/* Alocate and fill objective/Lagrangian gradient data and Jacobian data, 
	   convert indices from FORTRAN to C. */
	dims[0]=nnzg; 
	*Mgi=PyArray_SimpleNew(1, dims, npy_integer_type_num); 
	*Mgv=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
	gi=(npy_integer *)PyArray_DATA(*Mgi);
	gv=(npy_double *)PyArray_DATA(*Mgv);
	dims[0]=nnzjplusno-nnzg; 
	*MJi=PyArray_SimpleNew(1, dims, npy_integer_type_num); 
	*MJfi=PyArray_SimpleNew(1, dims, npy_integer_type_num); 
	*MJv=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
	Ji=(npy_integer *)PyArray_DATA(*MJi);
	Jfi=(npy_integer *)PyArray_DATA(*MJfi);
	Jv=(npy_double *)PyArray_DATA(*MJv);
	jg=0;
	jj=0;
	for(i=0;i<nnzjplusno;i++) {
		if (sjfi[i]==0) {
			gi[jg]=sji[i]-1;
			gv[jg]=sjv[i];
			jg++;
		} else {
			Ji[jj]=sji[i]-1;
			Jfi[jj]=sjfi[i]-1;
			Jv[jj]=sjv[i];
			jj++;
		}
	}
}

/* Extract sparse Hessian in form of NumPy arrays 
   from sparse ijv format of Hessian's upper triangle + diagonal. 
   Add elements to lower triangle. */
void extract_sparse_hessian(npy_integer nnzho, npy_integer *si, npy_integer *sj, npy_double *sv, 
					PyObject **MHi, PyObject **MHj, PyObject **MHv) {
	npy_integer *Hi, *Hj, nnzdiag, i, j;
	npy_double *Hv; 
	npy_intp dims[1];
	
	
	/* Get number of nonzeros on the diagonal */
	nnzdiag=0;
	for(i=0;i<nnzho;i++) {
		if (si[i]==sj[i])
			nnzdiag++;
	}
	
	/* Alocate and fill objective/Lagrangian gradient data and Jacobian data, 
	   convert indices from FORTRAN to C, fill lower triangle. */
	dims[0]=2*nnzho-nnzdiag; /* Do not duplicate diagonal elements */
	*MHi=PyArray_SimpleNew(1, dims, npy_integer_type_num); 
	*MHj=PyArray_SimpleNew(1, dims, npy_integer_type_num); 
	*MHv=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
	Hi=(npy_integer *)PyArray_DATA(*MHi);
	Hj=(npy_integer *)PyArray_DATA(*MHj);
	Hv=(npy_double *)PyArray_DATA(*MHv);
	j=0;
	for(i=0;i<nnzho;i++) {
		Hi[j]=si[i]-1;
		Hj[j]=sj[i]-1;
		Hv[j]=sv[i];
		j++;
		if (si[i]!=sj[i]) {
			Hi[j]=sj[i]-1;
			Hj[j]=si[i]-1;
			Hv[j]=sv[i];
			j++;
		}
	}
}


/* Functions */

static char cutest__dims_doc[]=
"Returns the dimension of the problem and the number of constraints.\n"
"\n"
"(n, m)=_dims()\n"
"\n"
"Output\n"
"n -- number of variables\n"
"m -- number of constraints\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"__init__.py script when the test function interface is loaded.\n"
"If you decide to call it anyway, the working directory at the time of call\n"
"must be the one where the file OUTSIF.d can be found.\n"
"\n"
"CUTEst tools used: CDIMEN\n";

static PyObject *cutest__dims(PyObject *self, PyObject *args) {
	if (PyObject_Length(args)!=0)
		PyErr_SetString(PyExc_Exception, "_dims() takes no arguments");

	if (!open_datafile())
		return NULL;
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: Calling CDIMEN\n");
#endif
	CUTEST_cdimen((integer *)&status, (integer *)&funit, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon);
#ifdef PYDEBUG
		fprintf(df, "PyCUTEst:   n = %-d, m = %-d\n", CUTEst_nvar, CUTEst_ncon);
#endif

	return decRefTuple(PyTuple_Pack(2, PyLong_FromLong((long)CUTEst_nvar), PyLong_FromLong((long)CUTEst_ncon))); 
}


static char cutest__setup_doc[]=
"Sets up the problem.\n"
"\n"
"data=_setup(efirst, lfirst, nvfirst)\n"
"\n"
"Input\n"
"efirst  -- if True, equation constraints are ordered before inequations.\n"
"           Defaults to False.\n"
"lfirst  -- if True, linear constraints are ordered before nonlinear ones.\n"
"           Defaults to False.\n"
"nvfirst -- if True, nonlinear variables are ordered before linear ones.\n"
"           Defaults to False.\n"
"\n"
"Setting both efirst and lfirst to True results in the following ordering:\n"
"linear equations, followed by linear inequations, nonlinear equations,\n"
"and finally nonlinear inequations.\n"
"\n"
"Output\n"
"data -- dictionary with the summary of test function's properties\n"
"\n"
"The problem data dictionary has the following members:\n"
"name    -- problem name\n"
"n       -- number of variables\n"
"m       -- number of constraints (excluding bounds)\n"
"x       -- initial point (1D array of length n)\n"
"bl      -- vector of lower bounds on variables (1D array of length n)\n"
"bu      -- vector of upper bounds on variables (1D array of length n)\n"
"nnzh    -- number of nonzero elements in the diagonal and upper triangle of\n"
"           sparse Hessian\n"
"vartype -- 1D integer array of length n storing variable type\n"
"           0=real,  1=boolean (0 or 1), 2=integer\n"
"\n"
"For constrained problems the following additional members are available\n"
"nnzj    -- number of nonzero elements in sparse Jacobian of constraints\n"
"v       -- initial value of Lagrange multipliers (1D array of length m)\n"
"cl      -- lower bounds on constraint functions (1D array of length m)\n"
"cu      -- upper bounds on constraint functions (1D array of length m)\n"
"equatn  -- 1D boolean array of length m indicating whether a constraint\n"
"           is an equation constraint\n"
"linear  -- 1D boolean array of length m indicating whether a constraint\n"
"           is a linear constraint\n"
"\n"
"-1e+20 and 1e+20 in bl, bu, cl, and cu stand for -infinity and +infinity.\n"
"\n"
"This function must be called before any other CUTEst function is called.\n"
"The only exception is the _dims() function.\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"__init__.py script when the test function interface is loaded.\n"
"If you decide to call it anyway, the working directory at the time of call\n"
"must be the one where the file OUTSIF.d can be found.\n"
"\n"
"CUTEst tools used: CDIMEN, CSETUP, USETUP, CVARTY, UVARTY, \n"
"                   CDIMSH, UDIMSH, CDIMSJ, PBNAME\n";

static PyObject *cutest__setup(PyObject *self, PyObject *args) {
	npy_logical  efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_;
	int eFirst, lFirst, nvFirst;
	int e_order = 1, l_order = 1, v_order = 0;
	PyObject *dict;
	PyObject *Mx, *Mbl, *Mbu, *Mv=NULL, *Mcl=NULL, *Mcu=NULL, *Meq=NULL, *Mlinear=NULL;
	PyObject *Mvt; 
	doublereal *x, *bl, *bu, *v=NULL, *cl=NULL, *cu=NULL;
	npy_integer *vartypes; 
	npy_logical *equatn=NULL, *linear=NULL;
	npy_intp dims[1];
	int i;
	
	if (PyObject_Length(args)!=0 && PyObject_Length(args)!=3) {
		PyErr_SetString(PyExc_Exception, "_setup() takes 0 or 3 arguments");
		return NULL; 
	}
	
	if (PyObject_Length(args)==3) {
		if (!PyArg_ParseTuple(args, "iii", &eFirst, &lFirst, &nvFirst)) {
			return NULL; 
		}

		efirst = eFirst  ? TRUE_ : FALSE_;
		lfirst = lFirst  ? TRUE_ : FALSE_;
		nvfrst = nvFirst ? TRUE_ : FALSE_;
	}
	
	if (!open_datafile())
		return NULL; 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: Calling CDIMEN\n");
#endif
	CUTEST_cdimen((integer *)&status, (integer *)&funit, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon);
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst:   n = %-d, m = %-d\n", CUTEst_nvar, CUTEst_ncon);
#endif

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: Allocating space\n");
#endif
	/* Numpy arrays */
	dims[0]=CUTEst_nvar;
	Mx=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	Mbl=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	Mbu=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	Mvt=PyArray_SimpleNew(1, dims, PyArray_INT); 
	if (CUTEst_ncon>0) {
		dims[0]=CUTEst_ncon;
		Mv=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
		Mcl=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
		Mcu=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
		Meq=PyArray_SimpleNew(1, dims, PyArray_BOOL);
		Mlinear=PyArray_SimpleNew(1, dims, PyArray_BOOL);
	}
	
	/* Get internal data buffers */
	/* Assume that npy_double is equivalent to double and npy_int is equivalent to integer */
	x = (npy_double *)PyArray_DATA(Mx); 
	bl = (npy_double *)PyArray_DATA(Mbl); 
	bu = (npy_double *)PyArray_DATA(Mbu); 
	if (CUTEst_ncon>0) {
		v = (npy_double *)PyArray_DATA(Mv); 
		cl = (npy_double *)PyArray_DATA(Mcl); 
		cu = (npy_double *)PyArray_DATA(Mcu); 
	
		/* Create temporary CUTEst logical arrays */
		equatn = (npy_logical *)malloc(CUTEst_ncon*sizeof(npy_logical));
		linear = (npy_logical *)malloc(CUTEst_ncon*sizeof(npy_logical));
	}
	vartypes=(npy_integer *)malloc(CUTEst_nvar*sizeof(npy_integer)); 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: Calling [UC]SETUP\n");
#endif
	if (CUTEst_ncon > 0)
		CUTEST_csetup((integer *)&status, (integer *)&funit, (integer *)&iout, (integer *)&io_buffer,
			     (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, bl, bu,
				v, cl, cu, (logical *)equatn, (logical *)linear,
				&e_order, &l_order, &v_order);
	else
		CUTEST_usetup((integer *)&status, (integer *)&funit, (integer *)&iout, (integer *)&io_buffer,
			      (integer *)&CUTEst_nvar, x, bl, bu);
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst:   n = %-d, m = %-d\n", CUTEst_nvar, CUTEst_ncon);
#endif

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: Calling [UC]VARTY\n");
#endif
	if (CUTEst_ncon > 0)
		CUTEST_cvartype((integer *)&status, (integer *)&CUTEst_nvar, (integer *)vartypes);
	else
		CUTEST_uvartype((integer *)&status, (integer *)&CUTEst_nvar, (integer *)vartypes);
		
	/* Copy logical values to NumPy bool arrays and free temporary storage */
	if (CUTEst_ncon > 0) {
		for(i=0; i<CUTEst_ncon; i++) {
			*((npy_bool*)(PyArray_GETPTR1(Meq, i)))=(equatn[i]==TRUE_ ? NPY_TRUE : NPY_FALSE);
			*((npy_bool*)(PyArray_GETPTR1(Mlinear, i)))=(linear[i]==TRUE_ ? NPY_TRUE : NPY_FALSE);
		}
		free(equatn);
		free(linear);
	} 
	
	/* Copy variable types to NumPy integer arrays and free temporary storage */
	for(i=0; i<CUTEst_nvar; i++) {
		*((npy_int*)PyArray_GETPTR1(Mvt, i))=vartypes[i];
	}
	free(vartypes);
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: Calling [CU]DIMSH\n");
#endif
	if (CUTEst_ncon>0)
		CUTEST_cdimsh((integer *)&status, (integer *)&CUTEst_nnzh);
	else 
		CUTEST_udimsh((integer *)&status, (integer *)&CUTEst_nnzh);
		
	if (CUTEst_ncon > 0) {
#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: Calling CDIMSJ\n");
#endif
		CUTEST_cdimsj((integer *)&status, (integer *)&CUTEst_nnzj);
		CUTEst_nnzj -= CUTEst_nvar;
	} 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst:   nnzh = %-d, nnzj = %-d\n", CUTEst_nnzh, CUTEst_nnzj);
	fprintf(df, "PyCUTEst: Finding out problem name\n");
#endif
	CUTEST_probname((integer *)&status, CUTEst_probName);
	trim_string(CUTEst_probName, STR_LEN-1);

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst:   %-s\n", CUTEst_probName);
	fprintf(df, "PyCUTEst: Closing data file\n");
#endif
	close_datafile();

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: Building structure\n");
#endif
	dict=PyDict_New();
	PyDict_SetItemString(dict, "n", PyLong_FromLong((long)CUTEst_nvar)); 
	PyDict_SetItemString(dict, "m", PyLong_FromLong((long)CUTEst_ncon)); 
	PyDict_SetItemString(dict, "nnzh", PyLong_FromLong((long)CUTEst_nnzh)); 
	PyDict_SetItemString(dict, "x", Mx); 
	PyDict_SetItemString(dict, "bl", Mbl); 
	PyDict_SetItemString(dict, "bu", Mbu); 
	PyDict_SetItemString(dict, "name", PyUnicode_FromString(CUTEst_probName)); 
	PyDict_SetItemString(dict, "vartype", Mvt); 
	if (CUTEst_ncon > 0) {
		PyDict_SetItemString(dict, "nnzj", PyLong_FromLong((long)CUTEst_nnzj)); 
		PyDict_SetItemString(dict, "v", (PyObject*)Mv); 
		PyDict_SetItemString(dict, "cl", (PyObject*)Mcl); 
		PyDict_SetItemString(dict, "cu", (PyObject*)Mcu); 
		PyDict_SetItemString(dict, "equatn", (PyObject*)Meq); 
		PyDict_SetItemString(dict, "linear", (PyObject*)Mlinear); 
	}
	
	setupCalled = 1;
	
	return decRefDict(dict); 
}


static char cutest__varnames_doc[]=
"Returns the names of variables in the problem.\n"
"\n"
"namelist=_varnames()\n"
"\n"
"Output\n"
"namelist -- list of length n holding strings holding names of variables\n"
"\n"
"The list reflects the ordering imposed by the nvfirst argument to _setup().\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"__init__.py script when the test function interface is loaded.\n"
"\n"
"CUTEst tools used: VARNAMES\n";

static PyObject *cutest__varnames(PyObject *self, PyObject *args) {
	char *Fvnames, Fvname[STR_LEN+1], *ptr;
	PyObject *list; 
	int i, j; 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
		
	if (PyObject_Length(args)!=0) {
		PyErr_SetString(PyExc_Exception, "_varnames() takes 0 arguments");
		return NULL; 
	}

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: allocating space\n");
#endif
	Fvnames=(char *)malloc(CUTEst_nvar*STR_LEN*sizeof(char));
	list=PyList_New(0);
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling VARNAMES\n");
#endif
	CUTEST_varnames((integer *)&status, (integer *)&CUTEst_nvar, Fvnames);

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: building results\n");
#endif
	for(i=0;i<CUTEst_nvar;i++) {
		ptr=Fvnames+i*STR_LEN; 
		for(j=0;j<STR_LEN;j++) {
			Fvname[j]=*ptr;
			ptr++;
		}
		trim_string(Fvname, STR_LEN-1);
		PyList_Append(list, PyUnicode_FromString(Fvname)); 
	}
	
	free(Fvnames); 
	
	return list;
}


static char cutest__connames_doc[]=
"Returns the names of constraints in the problem.\n"
"\n"
"namelist=_connames()\n"
"\n"
"Output\n"
"namelist -- list of length m holding strings holding names of constraints\n"
"\n"
"The list is ordered in the way specified by efirst and lfirst arguments to\n"
"_setup().\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"__init__.py script when the test function interface is loaded.\n"
"\n"
"CUTEst tools used: CONNAMES\n";

static PyObject *cutest__connames(PyObject *self, PyObject *args) {
	char *Fcnames, Fcname[STR_LEN+1], *ptr;
	PyObject *list; 
	int i, j; 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
		
	if (PyObject_Length(args)!=0) {
		PyErr_SetString(PyExc_Exception, "_connames() takes 0 arguments");
		return NULL; 
	}
	
	list=PyList_New(0);
	
	if (CUTEst_ncon>0) {

#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: allocating space\n");
#endif
		Fcnames=(char *)malloc(CUTEst_ncon*STR_LEN*sizeof(char));
	
#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: calling CONNAMES\n");
#endif
		CUTEST_connames((integer *)&status, (integer *)&CUTEst_ncon, Fcnames);

#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: building results\n");
#endif
		for(i=0;i<CUTEst_ncon;i++) {
			ptr=Fcnames+i*STR_LEN; 
			for(j=0;j<STR_LEN;j++) {
				Fcname[j]=*ptr;
				ptr++;
			}
			trim_string(Fcname, STR_LEN-1);
			PyList_Append(list, PyUnicode_FromString(Fcname)); 
		}
	
		free(Fcnames); 
	}
	
	return list;
}


static char cutest_objcons_doc[]=
"Returns the value of objective and constraints at x.\n"
"\n"
"(f, c)=objcons(x)\n"
"\n"
"Input\n"
"x -- 1D array of length n with the values of variables\n"
"\n"
"Output\n"
"f -- 1D array of length 1 holding the value of the function at x\n"
"c -- 1D array of length m holding the values of constraints at x\n"
"\n"
"CUTEst tools used: CFN\n";

static PyObject *cutest_objcons(PyObject *self, PyObject *args) {
	PyObject *arg1, *Mf, *Mc; 
	doublereal *x, *f, *c; 
	npy_intp dims[1]; 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	if (!PyArg_ParseTuple(args, "O", &arg1)) 
		return NULL; 
	
	/* Check if x is double and of correct length and shape */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1)&& PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	dims[0]=1;
	Mf=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	f=(npy_double *)PyArray_DATA(Mf);
	dims[0]=CUTEst_ncon;
	Mc=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	c=(npy_double *)PyArray_DATA(Mc);
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling CFN\n");
#endif
	CUTEST_cfn((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, f, c);
	
	return decRefTuple(PyTuple_Pack(2, Mf, Mc)); 
}


static char cutest_obj_doc[]=
"Returns the value of objective and its gradient at x.\n"
"\n"
"f=obj(x)\n"
"(f, g)=obj(x, gradFlag)\n"
"\n"
"Input\n"
"x        -- 1D array of length n with the values of variables\n"
"gradFlag -- if given the function returns f and g; can be anything\n"
"\n"
"Output\n"
"f -- 1D array of length 1 holding the value of the function at x\n"
"g -- 1D array of length n holding the value of the gradient of f at x\n"
"\n"
"CUTEst tools used: UOFG, COFG\n";

static PyObject *cutest_obj(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2;
	PyObject *Mf, *Mg=NULL; 
	doublereal *x, *f, *g=NULL; 
	npy_intp dims[1];

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	if (!PyArg_ParseTuple(args, "O|O", &arg1, &arg2)) 
		return NULL; 
	
	/* Check if x is double and of correct length and shape */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	dims[0]=1;
	Mf=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	f=(npy_double *)PyArray_DATA(Mf);
	if (PyObject_Length(args)>1) {
		dims[0]=CUTEst_nvar;
		Mg=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
		g=(npy_double *)PyArray_DATA(Mg);
	}

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling [UC]OFG\n");
#endif
	if (CUTEst_ncon == 0) {
		if (PyObject_Length(args)==1) {
			CUTEST_uofg((integer *)&status, (integer *)&CUTEst_nvar, x, f, NULL, &somethingFalse);
			return (PyObject *)Mf; 
		} else {
			CUTEST_uofg((integer *)&status, (integer *)&CUTEst_nvar, x, f, g, &somethingTrue);
			return decRefTuple(PyTuple_Pack(2, Mf, Mg)); 
		}
	} else {
		if (PyObject_Length(args)==1) {
			CUTEST_cofg((integer *)&status, (integer *)&CUTEst_nvar, x, f, NULL, &somethingFalse);
			return (PyObject *)Mf; 
		} else {
			CUTEST_cofg((integer *)&status, (integer *)&CUTEst_nvar, x, f, g, &somethingTrue);
			return decRefTuple(PyTuple_Pack(2, Mf, Mg)); 
		}
	}
}


static char cutest_cons_doc[]=
"Returns the value of constraints and the Jacobian of constraints at x.\n"
"\n"
"c=cons(x)                 -- constraints\n"
"ci=cons(x, False, i)      -- i-th constraint\n"
"(c, J)=cons(x, True)      -- Jacobian of constraints\n"
"(ci, Ji)=cons(x, True, i) -- i-th constraint and its gradient\n"
"\n"
"Input\n"
"x -- 1D array of length n with the values of variables\n"
"i -- integer index of constraint (between 0 and m-1)\n"
"\n"
"Output\n"
"c  -- 1D array of length m holding the values of constraints at x\n"
"ci -- 1D array of length 1 holding the value of i-th constraint at x\n"
"J  -- 2D array with m rows of n columns holding Jacobian of constraints at x\n"
"Ji -- 1D array of length n holding the gradient of i-th constraintn"
"      (also equal to the i-th row of Jacobian)\n"
"\n"
"CUTEst tools used: CCFG, CCIFG\n";

static PyObject *cutest_cons(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *Mc, *MJ; 
	doublereal *x, *c, *J; 
	int derivs, index, wantSingle;
	npy_integer icon; 
	npy_integer zero = 0;
	npy_intp dims[2];
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg2=NULL;
	if (!PyArg_ParseTuple(args, "O|Oi", &arg1, &arg2, &index)) 
		return NULL; 
	
	/* Check if x is double and of correct length and shape */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	/* Do we want derivatives */
	if (arg2!=NULL && arg2==Py_True)
		derivs=1;
	else
		derivs=0;
	
	/* Do we want a particular derivative */
	if (PyObject_Length(args)==3) {
		/* Check index */
		if (index<0 || index>=CUTEst_ncon) {
			PyErr_SetString(PyExc_Exception, "Argument 3 must be an integer between 0 and ncon-1");
			return NULL; 
		}
		icon=index+1;
		wantSingle=1;
	} else {
		wantSingle=0;
	}
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	if (!wantSingle) {
		dims[0]=CUTEst_ncon;
		Mc=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
		c=(npy_double *)PyArray_DATA(Mc);
		if (derivs) {
			dims[0]=CUTEst_ncon;
			dims[1]=CUTEst_nvar;
			/* Create a FORTRAN style array (first index stride is 1) */
			MJ=PyArray_New(&PyArray_Type, 2, dims, PyArray_DOUBLE, NULL, NULL, 0, NPY_F_CONTIGUOUS, NULL);
			J=(npy_double *)PyArray_DATA(MJ);
		}
	} else {
		dims[0]=1;
		Mc=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
		c=(npy_double *)PyArray_DATA(Mc);
		if (derivs) {
			dims[0]=CUTEst_nvar;
			MJ=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
			J=(npy_double *)PyArray_DATA(MJ);
		}
	}
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling CCFG/CCIFG\n");
#endif
	if (!wantSingle) {
		if (!derivs) {
			CUTEST_ccfg((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, c,
					&somethingFalse, (integer *)&zero, (integer *)&zero, NULL, &somethingFalse);
			return Mc;
		} else {
			CUTEST_ccfg((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, c,
					&somethingFalse, (integer *)&CUTEst_ncon, (integer *)&CUTEst_nvar, J,
					&somethingTrue);
			return decRefTuple(PyTuple_Pack(2, Mc, MJ)); 
		}
	} else {
		if (!derivs) {
			CUTEST_ccifg((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&icon, x, c, NULL, &somethingFalse);
			return Mc; 
		} else {
			CUTEST_ccifg((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&icon, x, c, J, &somethingTrue);
			return decRefTuple(PyTuple_Pack(2, Mc, MJ)); 
		}
	}
}


static char cutest_lagjac_doc[]=
"Returns the gradient of the objective or Lagrangian, and the Jacobian of\n"
"constraints at x. The gradient is the gradient with respect to the problem's\n"
"variables (has n components).\n"
"\n"
"(g, J)=lagjac(x)    -- objective gradient and the Jacobian of constraints\n"
"(g, J)=lagjac(x, v) -- Lagrangian gradient and the Jacobian of constraints\n"
"\n"
"Input\n"
"x -- 1D array of length n with the values of variables\n"
"v -- 1D array of length m with the Lagrange multipliers\n"
"\n"
"Output\n"
"g  -- 1D array of length n holding the gradient at x\n"
"J  -- 2D array with m rows of n columns holding Jacobian of constraints at x\n"
"\n"
"CUTEst tools used: CGR\n";

static PyObject *cutest_lagjac(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *Mg, *MJ; 
	doublereal *x, *v=NULL, *g, *J; 
	int lagrangian; 
	npy_intp dims[2];
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg2=NULL;
	if (!PyArg_ParseTuple(args, "O|O", &arg1, &arg2)) 
		return NULL; 
	
	/* Check if x is double and of correct length and shape */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	/* Check if v is double and of correct length and shape. */
	if (arg2!=NULL) {
		if (!(PyArray_Check(arg2) && PyArray_ISFLOAT(arg2) && PyArray_TYPE(arg2)==NPY_DOUBLE && PyArray_NDIM(arg2)==1 && PyArray_DIM(arg2, 0)==CUTEst_ncon)) {
			PyErr_SetString(PyExc_Exception, "Argument 2 must be a 1D double array of length ncon");
			return NULL; 
		}
		lagrangian=1;
	} else {
		lagrangian=0;
	}

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	if (lagrangian) 
		v=(npy_double *)PyArray_DATA(arg2);
	dims[0]=CUTEst_nvar;
	Mg=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	g=(npy_double *)PyArray_DATA(Mg);
	dims[0]=CUTEst_ncon;
	dims[1]=CUTEst_nvar;
	/* Create a FORTRAN style array (first index stride is 1) */
	MJ=PyArray_New(&PyArray_Type, 2, dims, PyArray_DOUBLE, NULL, NULL, 0, NPY_F_CONTIGUOUS, NULL);
	J=(npy_double *)PyArray_DATA(MJ);

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling CGR\n");
#endif
	if (!lagrangian) {
		CUTEST_cgr((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, NULL,
			&somethingFalse, g, &somethingFalse, (integer *)&CUTEst_ncon, (integer *)&CUTEst_nvar, J);
	} else {
		CUTEST_cgr((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, v,
			&somethingTrue, g, &somethingFalse, (integer *)&CUTEst_ncon, (integer *)&CUTEst_nvar, J);
	}
	
	return decRefTuple(PyTuple_Pack(2, Mg, MJ)); 
}


static char cutest_jprod_doc[]=
"Returns the product of constraints Jacobian at x with vector p\n"
"\n"
"r=jprod(transpose, p, x) -- computes Jacobian at x before product calculation\n"
"r=jprod(transpose, p)    -- uses last computed Jacobian\n"
"\n"
"Input\n"
"transpose -- boolean flag indicating that the Jacobian should be trasposed\n"
"             before the product is calculated\n"
"p         -- the vector that will be multiplied with the Jacobian\n"
"             1D array of length n (m) if transpose if False (True)\n"
"x         -- 1D array of length n holding the values of variables used in the\n"
"             evaluation of the constraints Jacobian\n"
"\n"
"Output\n"
"r  -- 1D array of length m if transpose=False (or n if transpose=True)\n"
"      with the result\n"
"\n"
"CUTEst tools used: CJPROD\n";

static PyObject *cutest_jprod(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *arg3, *Mr; 
	doublereal *p, *x, *r; 
	int transpose; 
	npy_intp dims[1];

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg3=NULL;
	if (!PyArg_ParseTuple(args, "OO|O", &arg1, &arg2, &arg3)) 
		return NULL; 
	
	/* Check if arg1 is True */
	if (arg1==Py_True) 
		transpose=1;
	else
		transpose=0; 
	
	/* Check if p is double and of correct dimension */
	if (!(PyArray_Check(arg2) && PyArray_ISFLOAT(arg2) && PyArray_TYPE(arg2)==NPY_DOUBLE && PyArray_NDIM(arg2)==1)) {
		PyErr_SetString(PyExc_Exception, "Argument 2 must be a 1D double array");
		return NULL; 
	}
	
	/* Check length of p when J is not transposed */
	if (!transpose && !(PyArray_DIM(arg2, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 2 must be of length nvar (J is not transposed)");
		return NULL;
	}
	
	/* Check length of p when J is transposed */
	if (transpose && !(PyArray_DIM(arg2, 0)==CUTEst_ncon)) {
		PyErr_SetString(PyExc_Exception, "Argument 2 must be of length ncon (J is transposed)");
		return NULL; 
	}
	
	/* Check if x is double and of correct length and shape. */
	if (arg3!=NULL) {
		if (!(arg3!=NULL && PyArray_Check(arg3) && PyArray_ISFLOAT(arg3) && PyArray_TYPE(arg3)==NPY_DOUBLE && PyArray_NDIM(arg3)==1 && PyArray_DIM(arg3, 0)==CUTEst_nvar)) {
			PyErr_SetString(PyExc_Exception, "Argument 3 must be a 1D double array of length nvar");
			return NULL; 
		}
	}

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	p=(npy_double *)PyArray_DATA(arg2);
	if (arg3!=NULL) 
		x=(npy_double *)PyArray_DATA(arg3);
	else
		x=NULL;
	if (!transpose) {
		dims[0]=CUTEst_ncon;
	} else {
		dims[0]=CUTEst_nvar;
	}
	Mr=PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
	r=(npy_double *)PyArray_DATA(Mr);

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling CJPROD\n");
#endif
	if (!transpose) {
		if (arg3==NULL) {
			CUTEST_cjprod((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, &somethingTrue,
					&somethingFalse, NULL, p, (integer *)&CUTEst_nvar, r, (integer *)&CUTEst_ncon);
		} else {
			CUTEST_cjprod((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, &somethingFalse,
					&somethingFalse, x, p, (integer *)&CUTEst_nvar, r, (integer *)&CUTEst_ncon);
		}
	} else {
		if (arg3==NULL) {
			CUTEST_cjprod((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, &somethingTrue,
					&somethingTrue, NULL, p, (integer *)&CUTEst_ncon, r, (integer *)&CUTEst_nvar);
		} else {
			CUTEST_cjprod((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, &somethingFalse,
					&somethingTrue, x, p, (integer *)&CUTEst_ncon, r, (integer *)&CUTEst_nvar);
		}
	}
	
	return Mr;
}


static char cutest_hess_doc[]=
"Returns the Hessian of the objective (for unconstrained problems) or the\n"
"Hessian of the Lagrangian (for constrained problems) at x.\n" 
"\n"
"H=hess(x)    -- Hessian of objective at x for unconstrained problems\n"
"H=hess(x, v) -- Hessian of Lagrangian at (x, v) for constrained problems\n"
"\n"
"The first form can only be used for unconstrained problems. The second one\n"
"can only be used for constrained problems. For obtaining the Hessian of the\n"
"objective in case of a constrained problem use ihess().\n"
"\n"
"The Hessian is meant with respect to problem variables (has dimension n).\n"
"\n"
"Input\n"
"x         -- 1D array of length n holding the values of variables\n"
"v         -- 1D array of length m holding the values of Lagrange multipliers\n"
"\n"
"Output\n"
"H  -- 2D array with n rows of n columns holding the Hessian at x (or (x, v))\n"
"\n"
"CUTEst tools used: CDH, UDH\n";

static PyObject *cutest_hess(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *MH; 
	doublereal *x, *v=NULL, *H; 
	npy_intp dims[2];
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg2=NULL;
	if (!PyArg_ParseTuple(args, "O|O", &arg1, &arg2)) 
		return NULL; 
	
	/* Check if x is double and of correct dimension */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	if (CUTEst_ncon>0) {
		/* Check if v is double and of correct dimension */
		if (arg2!=NULL) {
			if (!(PyArray_Check(arg2) && PyArray_ISFLOAT(arg2) && PyArray_TYPE(arg2)==NPY_DOUBLE && PyArray_NDIM(arg2)==1 && PyArray_DIM(arg2, 0)==CUTEst_ncon)) {
				PyErr_SetString(PyExc_Exception, "Argument 2 must be a 1D double array of length ncon");
				return NULL; 
			}
		} else {
			PyErr_SetString(PyExc_Exception, "Argument 2 must be specified for constrained problems. Use ihess().");
			return NULL; 
		}
	}
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	if (CUTEst_ncon>0) 
		v=(npy_double *)PyArray_DATA(arg2);
	dims[0]=CUTEst_nvar;
	dims[1]=CUTEst_nvar;
	/* Create a FORTRAN style array (first index stride is 1) */
	MH=PyArray_New(&PyArray_Type, 2, dims, PyArray_DOUBLE, NULL, NULL, 0, NPY_F_CONTIGUOUS, NULL);
	H=(npy_double *)PyArray_DATA(MH);

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling [CU]DH\n");
#endif
	if (CUTEst_ncon>0) {
		CUTEST_cdh((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, v, (integer *)&CUTEst_nvar, H);
	} else {
		CUTEST_udh((integer *)&status, (integer *)&CUTEst_nvar, x, (integer *)&CUTEst_nvar, H);
	}
	
	return MH; 
}


static char cutest_ihess_doc[]=
"Returns the Hessian of the objective or the Hessian of i-th constraint at x.\n"
"\n"
"H=ihess(x)    -- Hessian of the objective\n"
"H=ihess(x, i) -- Hessian of i-th constraint\n"
"\n"
"The first form can only be used for unconstrained problems. The second one\n"
"can only be used for constrained problems. For obtaining the Hessian of the\n"
"objective in case of a constrained problem use ihess().\n"
"\n"
"The Hessian is meant with respect to problem variables (has dimension n).\n"
"\n"
"Input\n"
"x -- 1D array of length n holding the values of variables\n"
"i -- integer holding the index of the constraint (between 0 and m-1)\n"
"\n"
"Output\n"
"H  -- 2D array with n rows of n columns holding the Hessian at x\n"
"\n"
"CUTEst tools used: CIDH, UDH\n";

static PyObject *cutest_ihess(PyObject *self, PyObject *args) {
	PyObject *arg1, *MH; 
	doublereal *x, *H; 
	npy_intp dims[2];
	int i; 
	npy_integer icon; 

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	if (!PyArg_ParseTuple(args, "O|i", &arg1, &i)) 
		return NULL; 
	
	if (PyObject_Length(args)>1) {
		icon=i+1;
		if (i<0 || i>=CUTEst_ncon) {
			PyErr_SetString(PyExc_Exception, "Argument 2 must be between 0 and ncon-1");
			return NULL; 
		}
	} else {
		icon=0; 
	}
	
	/* Check if x is double and of correct dimension */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	dims[0]=CUTEst_nvar;
	dims[1]=CUTEst_nvar;
	/* Create a FORTRAN style array (first index stride is 1) */
	MH=PyArray_New(&PyArray_Type, 2, dims, PyArray_DOUBLE, NULL, NULL, 0, NPY_F_CONTIGUOUS, NULL);
	H=(npy_double *)PyArray_DATA(MH);	
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling CIDH/UDH\n");
#endif
	if (CUTEst_ncon>0) {
		CUTEST_cidh((integer *)&status, (integer *)&CUTEst_nvar, x, (integer *)&icon, (integer *)&CUTEst_nvar, H);
	} else {
		CUTEST_udh((integer *)&status, (integer *)&CUTEst_nvar, x, (integer *)&CUTEst_nvar, H);
	}
	
	return MH;
}


static char cutest_hprod_doc[]=
"Returns the product of Hessian at x and vector p.\n"
"The Hessian is either the Hessian of objective or the Hessian of Lagrangian.\n" 
"\n"
"r=hprod(p, x, v) -- use Hessian of Lagrangian at x (constrained problem)\n"
"r=hprod(p, x)    -- use Hessian of objective at x (unconstrained problem)\n"
"r=hprod(p)       -- use last computed Hessian\n"
"\n"
"The first form can only be used for constrained problems. The second one\n"
"can only be used for unconstrained problems.\n"
"\n"
"The Hessian is meant with respect to problem variables (has dimension n).\n"
"\n"
"Input\n"
"p -- 1D array of length n holding the components of the vector\n"
"x -- 1D array of length n holding the values of variables\n"
"v -- 1D array of length m holding the values of Lagrange multipliers\n"
"\n"
"Output\n"
"r  -- 1D array of length n holding the result\n"
"\n"
"CUTEst tools used: CPROD, UPROD\n";

static PyObject *cutest_hprod(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *arg3, *Mr; 
	doublereal *p, *x=NULL, *v=NULL, *r; 
	npy_intp dims[1];

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg2=arg3=NULL; 
	if (!PyArg_ParseTuple(args, "O|OO", &arg1, &arg2, &arg3)) 
		return NULL; 
	
	if (CUTEst_ncon>0) {
		if (PyObject_Length(args)==2) {
			PyErr_SetString(PyExc_Exception, "Need 1 or 3 arguments for constrained problems");
			return NULL; 
		}
	} else {
		if (PyObject_Length(args)==3) {
			PyErr_SetString(PyExc_Exception, "Need 1 or 2 arguments for unconstrained problems");
			return NULL; 
		}
	}
	
	/* Check if p is double and of correct dimension */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	/* Check if x is double and of correct dimension */
	if (arg2!=NULL && !(PyArray_Check(arg2) && PyArray_ISFLOAT(arg2) && PyArray_TYPE(arg2)==NPY_DOUBLE && PyArray_NDIM(arg2)==1 && PyArray_DIM(arg2, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 2 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	/* Check if v is double and of correct dimension */
	if (arg3!=NULL && !(PyArray_Check(arg3) && PyArray_ISFLOAT(arg3) && PyArray_TYPE(arg3)==NPY_DOUBLE && PyArray_NDIM(arg3)==1 && PyArray_DIM(arg3, 0)==CUTEst_ncon)) {
		PyErr_SetString(PyExc_Exception, "Argument 3 must be a 1D double array of length ncon");
		return NULL; 
	}
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	p=(npy_double *)PyArray_DATA(arg1);
	if (arg2!=NULL)
		x=(npy_double *)PyArray_DATA(arg2);
	if (arg3!=NULL)
		v=(npy_double *)PyArray_DATA(arg3);
	dims[0]=CUTEst_nvar;
	Mr=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
	r=(npy_double *)PyArray_DATA(Mr);	

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling [CU]PROD\n");
#endif
	if (CUTEst_ncon>0) {
		if (arg2==NULL)
			CUTEST_chprod((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, &somethingTrue, NULL, NULL, p, r);
		else 
			CUTEST_chprod((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, &somethingFalse, x, v, p, r);
	} else {
		if (arg2==NULL)
			CUTEST_uhprod((integer *)&status, (integer *)&CUTEst_nvar, &somethingTrue, NULL, p, r);
		else 
			CUTEST_uhprod((integer *)&status, (integer *)&CUTEst_nvar, &somethingFalse, x, p, r);
	}
	
	return Mr;
}


static char cutest_gradhess_doc[]=
"Returns the Hessian of the Lagrangian, the Jacobian of constraints, and the\n"
"gradient of the objective or the gradient of the Lagrangian at.\n" 
"\n"
"(g, H)=gradhess(x)       -- for unconstrained problems\n"
"(g, J, H)=gradhess(x, v, gradl) -- for constrained problems\n"
"\n"
"The first form can only be used for unconstrained problems. The second one\n"
"can only be used for constrained problems.\n"
"\n"
"The Hessian is meant with respect to problem variables (has dimension n).\n"
"\n"
"Input\n"
"x     -- 1D array of length n holding the values of variables\n"
"v     -- 1D array of length m holding the values of Lagrange multipliers\n"
"gradl -- boolean flag. If False the gradient of the objective is returned, \n"
"         if True the gradient of the Lagrangian is returned.\n"
"         Default is False.\n"
"\n"
"Output\n"
"g  -- 1D array of length n holding\n"
"      the gradient of objective at x (for unconstrained problems) or\n"
"      the gradient of Lagrangian at (x, v) (for constrained problems)\n"
"J  -- 2D array with m rows and n columns holding the Jacobian of constraints\n"
"H  -- 2D array with n rows and n columns holding\n"
"      the Hessian of the objective (for unconstrained problems) or\n"
"      the Hessian of the Lagrangian (for constrained problems)\n"
"\n"
"CUTEst tools used: CGRDH, UGRDH\n";

static PyObject *cutest_gradhess(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *arg3, *Mg, *MH, *MJ; 
	doublereal *x, *v=NULL, *g, *H, *J; 
	npy_logical grlagf; 
	npy_intp dims[2];
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg2=arg3=NULL; 
	if (!PyArg_ParseTuple(args, "O|OO", &arg1, &arg2, &arg3)) 
		return NULL; 
	
	if (CUTEst_ncon>0) {
		if (PyObject_Length(args)<2) {
			PyErr_SetString(PyExc_Exception, "Need at least 2 arguments for constrained problems");
			return NULL; 
		}
	} else {
		if (PyObject_Length(args)!=1) {
			PyErr_SetString(PyExc_Exception, "Need 1 argument for unconstrained problems");
			return NULL; 
		}
	}
	
	/* Check if x is double and of correct dimension */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	/* Check if v is double and of correct dimension */
	if (arg2!=NULL && !(PyArray_Check(arg2) && PyArray_ISFLOAT(arg2) && PyArray_TYPE(arg2)==NPY_DOUBLE && PyArray_NDIM(arg2)==1 && PyArray_DIM(arg2, 0)==CUTEst_ncon)) {
		PyErr_SetString(PyExc_Exception, "Argument 2 must be a 1D double array of length ncon");
		return NULL; 
	}
	
	/* Are we computing the gradient of the Lagrangian */
	if (arg3!=NULL && arg3==Py_True) {
		grlagf=TRUE_;
	} else {
		grlagf=FALSE_; 
	}

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	if (arg2!=NULL)
		v=(npy_double *)PyArray_DATA(arg2);
	dims[0]=CUTEst_nvar;
	Mg=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
	g=(npy_double *)PyArray_DATA(Mg);	
	dims[0]=CUTEst_nvar;
	dims[1]=CUTEst_nvar;
	/* Create a FORTRAN style array (first index stride is 1) */
	MH=PyArray_New(&PyArray_Type, 2, dims, PyArray_DOUBLE, NULL, NULL, 0, NPY_F_CONTIGUOUS, NULL);
	H=(npy_double *)PyArray_DATA(MH);
	dims[0]=CUTEst_ncon;
	dims[1]=CUTEst_nvar;
	/* Create a FORTRAN style array (first index stride is 1) */
	MJ=PyArray_New(&PyArray_Type, 2, dims, PyArray_DOUBLE, NULL, NULL, 0, NPY_F_CONTIGUOUS, NULL);
	J=(npy_double *)PyArray_DATA(MJ);

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling [CU]GRDH\n");
#endif
	if (CUTEst_ncon>0) {
		CUTEST_cgrdh((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, v,
			     (logical *)&grlagf, g, &somethingFalse, (integer *)&CUTEst_ncon,
			     (integer *)&CUTEst_nvar, J, (integer *)&CUTEst_nvar, H);
		return decRefTuple(PyTuple_Pack(3, Mg, MJ, MH)); 
	} else {
		CUTEST_ugrdh((integer *)&status, (integer *)&CUTEst_nvar, x, g, (integer *)&CUTEst_nvar, H);
		return decRefTuple(PyTuple_Pack(2, Mg, MH)); 
	}
}


static char cutest__scons_doc[]=
"Returns the value of constraints and the sparse Jacobian of constraints at x.\n"
"\n"
"(c, Jvi, Jfi, Jv)=_scons(x) -- Jacobian of constraints\n"
"(ci, gi, gv)=_scons(x, i)   -- i-th constraint and its gradient\n"
"\n"
"Input\n"
"x -- 1D array of length n with the values of variables\n"
"i -- integer index of constraint (between 0 and m-1)\n"
"\n"
"Output\n"
"c   -- 1D array of length m holding the values of constraints at x\n"
"Jvi -- 1D array of integers holding the column indices (0 .. n-1)\n"
"       of nozero elements in sparse Jacobian of constraints\n"
"Jfi -- 1D array of integers holding the row indices (0 .. m-1)\n"
"       of nozero elements in sparse Jacobian of constraints\n"
"Jv  -- 1D array holding the values of nonzero elements in the sparse Jacobian\n"
"       of constraints at x. Has the same length as Jvi and Jfi.\n"
"ci  -- 1D array of length 1 with the value of i-th constraint at x\n"
"gi  -- 1D array of integers holding the indices (0 .. n-1) of nonzero\n"
"       elements in the sparse gradient vector of i-th constraint\n"
"gv  -- 1D array holding the values of nonzero elements in the sparse gradient\n"
"       vector representing the gradient of i-th constraint at x.\n"
"       Has the same length as gi. gi and gv corespond to the i-th row of\n"
"       constraints Jacobian at x.\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"wrapper function scons().\n"
"\n"
"CUTEst tools used: CCFSG, CCIFSG\n";

static PyObject *cutest__scons(PyObject *self, PyObject *args) {
	PyObject *arg1, *Mc, *MJi, *MJfi, *MJv, *Mgi, *Mgv; 
	doublereal *c, *Jv, *gv, *x, *sv; 
	npy_integer *Ji, *Jfi, *gi, *si;  
	npy_integer index, nnzsgc; 
	int i; 
	npy_intp dims[1]; 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	if (!PyArg_ParseTuple(args, "O|i", &arg1, &i)) 
		return NULL; 
		
	/* Check if x is double and of correct dimension */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	if (PyObject_Length(args)==2) {
		if (i<0 || i>=CUTEst_ncon) {
			PyErr_SetString(PyExc_Exception, "Argument 2 must be an integer between 0 and ncon-1");
			return NULL; 
		}
		index=i+1;
	}
	
	if (PyObject_Length(args)==1) {
#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
		x=(npy_double *)PyArray_DATA(arg1);
		dims[0]=CUTEst_nnzj;
		MJi=PyArray_SimpleNew(1, dims, npy_integer_type_num); 
		MJfi=PyArray_SimpleNew(1, dims, npy_integer_type_num); 
		MJv=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
		Ji=(npy_integer *)PyArray_DATA(MJi);	
		Jfi=(npy_integer *)PyArray_DATA(MJfi);	
		Jv=(npy_double *)PyArray_DATA(MJv);	
		dims[0]=CUTEst_ncon; 
		Mc=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
		c=(npy_double *)PyArray_DATA(Mc);	
	
#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: calling CCFSG\n");
#endif
		CUTEST_ccfsg((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, c, (integer *)&CUTEst_nnzj,
              (integer *)&CUTEst_nnzj, Jv, (integer *)Ji, (integer *)Jfi, &somethingTrue);
		
		/* Convert FORTRAN indices to C indices */
		for(i=0;i<CUTEst_nnzj;i++) {
			Ji[i]--;
			Jfi[i]--;
		}
		
		return decRefTuple(PyTuple_Pack(4, Mc, MJi, MJfi, MJv));
	} else {
#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
		x=(npy_double *)PyArray_DATA(arg1);
		si=(npy_integer *)malloc(CUTEst_nvar*sizeof(npy_integer));
		sv=(npy_double *)malloc(CUTEst_nvar*sizeof(npy_double));
		dims[0]=1;
		Mc=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
		c=(npy_double *)PyArray_DATA(Mc);	

#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: calling CCIFSG\n");
#endif
		CUTEST_ccifsg((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&index, x, c, (integer *)&nnzsgc, (integer *)&CUTEst_nvar, sv, (integer *)si, &somethingTrue);
		
		/* Allocate and copy results, convert indices from FORTRAN to C, free storage */
		dims[0]=nnzsgc; 
		Mgi=PyArray_SimpleNew(1, dims, npy_integer_type_num); 
		Mgv=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
		gi=(npy_integer *)PyArray_DATA(Mgi);
		gv=(npy_double *)PyArray_DATA(Mgv);	
		for (i=0;i<nnzsgc;i++) {
			gi[i]=si[i]-1;
			gv[i]=sv[i];
		}
		free(si);
		free(sv);
		
		return decRefTuple(PyTuple_Pack(3, Mc, Mgi, Mgv)); 
	}
}


static char cutest__slagjac_doc[]=
"Returns the sparse gradient of objective at x or Lagrangian at (x, v), \n"
"and the sparse Jacobian of constraints at x.\n"
"\n"
"(gi, gv, Jvi, Jfi, Jv)=_slagjac(x)    -- objective gradient and Jacobian\n"
"(gi, gv, Jvi, Jfi, Jv)=_slagjac(x, v) -- Lagrangian gradient and Jacobian\n"
"\n"
"Input\n"
"x -- 1D array of length n with the values of variables\n"
"v -- 1D array of length m with the values of Lagrange multipliers\n"
"\n"
"Output\n"
"gi  -- 1D array of integers holding the indices (0 .. n-1) of nonzero\n"
"       elements in the sparse gradient vector\n"
"gv  -- 1D array holding the values of nonzero elements in the sparse gradient\n"
"       vector. Has the same length as gi.\n"
"Jvi -- 1D array of integers holding the column indices (0 .. n-1)\n"
"       of nozero elements in sparse Jacobian of constraints\n"
"Jfi -- 1D array of integers holding the row indices (0 .. m-1)\n"
"       of nozero elements in sparse Jacobian of constraints\n"
"Jv  -- 1D array holding the values of nonzero elements in the sparse Jacobian\n"
"       of constraints at x. Has the same length as Jvi and Jfi.\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"wrapper function slagjac().\n"
"\n"
"CUTEst tools used: CSGR\n";

static PyObject *cutest__slagjac(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *Mgi, *Mgv, *MJi, *MJfi, *MJv; 
	doublereal *x, *v=NULL, *sv; 
	npy_integer *si, *sfi;
	npy_integer nnzjplusn, nnzjplusno; 
	int lagrangian; 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg2=NULL;
	if (!PyArg_ParseTuple(args, "O|O", &arg1, &arg2)) 
		return NULL; 
	
	/* Check if x is double and of correct length and shape */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	/* Check if v is double and of correct length and shape. */
	if (arg2!=NULL) {
		if (!(PyArray_Check(arg2) && PyArray_ISFLOAT(arg2) && PyArray_TYPE(arg2)==NPY_DOUBLE && PyArray_NDIM(arg2)==1 && PyArray_DIM(arg2, 0)==CUTEst_ncon)) {
			PyErr_SetString(PyExc_Exception, "Argument 2 must be a 1D double array of length ncon");
			return NULL; 
		}
		lagrangian=1;
	} else {
		lagrangian=0;
	}

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	if (lagrangian) 
		v=(npy_double *)PyArray_DATA(arg2);
	nnzjplusn=CUTEst_nnzj+CUTEst_nvar; 
	si=(npy_integer *)malloc(nnzjplusn*sizeof(npy_integer));
	sfi=(npy_integer *)malloc(nnzjplusn*sizeof(npy_integer));
	sv=(npy_double *)malloc(nnzjplusn*sizeof(npy_double)); 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling CSGR\n");
#endif
	/* Must use different variable for output NNZJ and input LCJAC */
	if (!lagrangian) {
		CUTEST_csgr((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, NULL,
			   &somethingFalse, (integer *)&nnzjplusno, (integer *)&nnzjplusn, sv, (integer *)si, (integer *)sfi);
	} else {
		CUTEST_csgr((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, v,
			   &somethingFalse, (integer *)&nnzjplusno, (integer *)&nnzjplusn, sv, (integer *)si, (integer *)sfi);
	}
	
	extract_sparse_gradient_jacobian(nnzjplusno, si, sfi, sv, &Mgi, &Mgv, &MJi, &MJfi, &MJv); 
	
	/* Free temporary storage */
	free(si);
	free(sfi);
	free(sv);
	
	return decRefTuple(PyTuple_Pack(5, Mgi, Mgv, MJi, MJfi, MJv)); 
}


static char cutest__sphess_doc[]=
"Returns the sparse Hessian of the objective at x (unconstrained problems) or\n"
"the sparse Hessian of the Lagrangian (constrained problems) at (x, v).\n"
"\n"
"(Hi, Hj, Hv)=_sphess(x)    -- Hessian of objective (unconstrained problems)\n"
"(Hi, Hj, Hv)=_sphess(x, v) -- Hessian of Lagrangian (constrained problems)\n"
"\n"
"Input\n"
"x -- 1D array of length n with the values of variables\n"
"v -- 1D array of length m with the values of Lagrange multipliers\n"
"\n"
"Output\n"
"Hi -- 1D array of integers holding the row indices (0 .. n-1)\n"
"      of nozero elements in sparse Hessian\n"
"Hj -- 1D array of integers holding the column indices (0 .. n-1)\n"
"      of nozero elements in sparse Hessian\n"
"Hv -- 1D array holding the values of nonzero elements in the sparse Hessian\n"
"      Has the same length as Hi and Hj.\n"
"\n"
"Hi, Hj, and Hv represent the full Hessian and not only the diagonal and the\n"
"upper triangle. To obtain the Hessian of the objective of constrained\n"
"problems use _isphess().\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"wrapper function sphess().\n"
"\n"
"CUTEst tools used: CSH, USH\n";

static PyObject *cutest__sphess(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *MHi, *MHj, *MHv; 
	doublereal *x, *v=NULL, *sv; 
	npy_integer *si, *sj, nnzho; 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg2=NULL;
	if (!PyArg_ParseTuple(args, "O|O", &arg1, &arg2)) 
		return NULL; 
	
	/* Check if x is double and of correct dimension */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	if (CUTEst_ncon>0) {
		/* Check if v is double and of correct dimension */
		if (arg2!=NULL) {
			if (!(PyArray_Check(arg2) && PyArray_ISFLOAT(arg2) && PyArray_TYPE(arg2)==NPY_DOUBLE && PyArray_NDIM(arg2)==1 && PyArray_DIM(arg2, 0)==CUTEst_ncon)) {
				PyErr_SetString(PyExc_Exception, "Argument 2 must be a 1D double array of length ncon");
				return NULL; 
			}
		} else {
			PyErr_SetString(PyExc_Exception, "Argument 2 must be specified for constrained problems. Use _isphess().");
			return NULL; 
		}
	}
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	if (CUTEst_ncon>0) 
		v=(npy_double *)PyArray_DATA(arg2);
	si=(npy_integer *)malloc(CUTEst_nnzh*sizeof(npy_integer));
	sj=(npy_integer *)malloc(CUTEst_nnzh*sizeof(npy_integer));
	sv=(npy_double *)malloc(CUTEst_nnzh*sizeof(npy_double)); 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling [CU]SH\n");
#endif
	if (CUTEst_ncon>0) {
		CUTEST_csh((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, v,
			   (integer *)&nnzho, (integer *)&CUTEst_nnzh, sv, (integer *)si, (integer *)sj);
	} else {
		CUTEST_ush((integer *)&status, (integer *)&CUTEst_nvar, x,
			   (integer *)&nnzho, (integer *)&CUTEst_nnzh, sv, (integer *)si, (integer *)sj);
	}
	
	extract_sparse_hessian(nnzho, si, sj, sv, &MHi, &MHj, &MHv); 
	
	/* Free temporary storage */
	free(si);
	free(sj);
	free(sv);
	
	return decRefTuple(PyTuple_Pack(3, MHi, MHj, MHv)); 
}


static char cutest__isphess_doc[]=
"Returns the sparse Hessian of the objective or the sparse Hessian of i-th\n"
"constraint at x.\n"
"\n"
"(Hi, Hj, Hv)=_isphess(x)    -- Hessian of objective\n"
"(Hi, Hj, Hv)=_isphess(x, i) -- Hessian of i-th constraint\n"
"\n"
"Input\n"
"x -- 1D array of length n with the values of variables\n"
"i -- integer holding the index of constraint (between 0 and m-1)\n"
"\n"
"Output\n"
"Hi -- 1D array of integers holding the row indices (0 .. n-1)\n"
"      of nozero elements in sparse Hessian\n"
"Hj -- 1D array of integers holding the column indices (0 .. n-1)\n"
"      of nozero elements in sparse Hessian\n"
"Hv -- 1D array holding the values of nonzero elements in the sparse Hessian\n"
"      Has the same length as Hi and Hj.\n"
"\n"
"Hi, Hj, and Hv represent the full Hessian and not only the diagonal and the\n"
"upper triangle.\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"wrapper function isphess().\n"
"\n"
"CUTEst tools used: CISH, USH\n";

static PyObject *cutest__isphess(PyObject *self, PyObject *args) {
	PyObject *arg1, *MHi, *MHj, *MHv; 
	doublereal *x, *sv; 
	npy_integer *si, *sj, nnzho, i; 
	npy_integer icon; 

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	if (!PyArg_ParseTuple(args, "O|i", &arg1, &i)) 
		return NULL; 
	
	if (PyObject_Length(args)>1) {
		icon=i+1;
		if (i<0 || i>=CUTEst_ncon) {
			PyErr_SetString(PyExc_Exception, "Argument 2 must be between 0 and ncon-1");
			return NULL; 
		}
	} else {
		icon=0; 
	}
	
	/* Check if x is double and of correct dimension */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	si=(npy_integer *)malloc(CUTEst_nnzh*sizeof(npy_integer));
	sj=(npy_integer *)malloc(CUTEst_nnzh*sizeof(npy_integer));
	sv=(npy_double *)malloc(CUTEst_nnzh*sizeof(npy_double)); 
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: calling CISH/USH\n");
#endif
	if (CUTEst_ncon>0) {
		CUTEST_cish((integer *)&status, (integer *)&CUTEst_nvar, x, (integer *)&icon,
			    (integer *)&nnzho, (integer *)&CUTEst_nnzh, sv, (integer *)si, (integer *)sj);
	} else {
		CUTEST_ush((integer *)&status, (integer *)&CUTEst_nvar, x,
			   (integer *)&nnzho, (integer *)&CUTEst_nnzh, sv, (integer *)si, (integer *)sj);
	}
	
	extract_sparse_hessian(nnzho, si, sj, sv, &MHi, &MHj, &MHv); 
	
	/* Free temporary storage */
	free(si);
	free(sj);
	free(sv);
	
	return decRefTuple(PyTuple_Pack(3, MHi, MHj, MHv)); 
}


static char cutest__gradsphess_doc[]=
"Returns the sparse Hessian of the Lagrangian, the sparse Jacobian of\n"
"constraints, and the gradient of the objective or Lagrangian.\n" 
"\n"
"(g, Hi, Hj, Hv)=_gradsphess(x) -- unconstrained problems\n"
"(gi, gv, Jvi, Jfi, Jv, Hi, Hj, Hv)=_gradsphess(x, v, gradl)\n"
"                               -- constrained problems\n"
"\n"
"Input\n"
"x     -- 1D array of length n with the values of variables\n"
"v     -- 1D array of length m holding the values of Lagrange multipliers\n"
"gradl -- boolean flag. If False the gradient of the objective is returned, \n"
"         if True the gradient of the Lagrangian is returned. Default is False.\n"
"\n"
"Output\n"
"g   -- 1D array of length n with the gradient of objective or Lagrangian\n"
"Hi  -- 1D array of integers holding the row indices (0 .. n-1)\n"
"       of nozero elements in sparse Hessian\n"
"Hj  -- 1D array of integers holding the column indices (0 .. n-1)\n"
"       of nozero elements in sparse Hessian\n"
"Hv  -- 1D array holding the values of nonzero elements in the sparse Hessian\n"
"       Has the same length as Hi and Hj.\n"
"gi  -- 1D array of integers holding the indices (0 .. n-1) of nonzero\n"
"       elements in the sparse gradient vector\n"
"gv  -- 1D array holding the values of nonzero elements in the sparse gradient\n"
"       vector. Has the same length as gi.\n"
"Jvi -- 1D array of integers holding the column indices (0 .. n-1)\n"
"       of nozero elements in sparse Jacobian of constraints\n"
"Jfi -- 1D array of integers holding the row indices (0 .. m-1)\n"
"       of nozero elements in sparse Jacobian of constraints\n"
"Jv  -- 1D array holding the values of nonzero elements in the sparse Jacobian\n"
"       of constraints at x. Has the same length as Jvi and Jfi.\n"
"\n"
"For constrained problems the gradient is returned in sparse format.\n"
"\n"
"Hi, Hj, and Hv represent the full Hessian and not only the diagonal and the\n"
"upper triangle.\n"
"\n"
"This function is not supposed to be called by the user. It is called by the\n"
"wrapper function gradsphess().\n"
"\n"
"CUTEst tools used: CSGRSH, UGRSH\n";

static PyObject *cutest__gradsphess(PyObject *self, PyObject *args) {
	PyObject *arg1, *arg2, *arg3, *Mg=NULL, *Mgi, *Mgv, *MJi, *MJfi, *MJv, *MHi, *MHj, *MHv; 
	doublereal *x, *v, *g, *sv, *sjv; 
	npy_integer lagrangian;
	npy_integer *si, *sj, *sji, *sjfi, nnzho, nnzjplusn, nnzjplusno; 
	npy_intp dims[1]; 

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
	
	arg2=NULL;
	arg3=NULL;
	if (!PyArg_ParseTuple(args, "O|OO", &arg1, &arg2, &arg3)) 
		return NULL; 	
	
	/* Check bool argument */
	if (arg3!=NULL && arg3==Py_True)
		lagrangian=1;
	else
		lagrangian=0;
	
	/* Check if x is double and of correct dimension */
	if (!(PyArray_Check(arg1) && PyArray_ISFLOAT(arg1) && PyArray_TYPE(arg1)==NPY_DOUBLE && PyArray_NDIM(arg1)==1 && PyArray_DIM(arg1, 0)==CUTEst_nvar)) {
		PyErr_SetString(PyExc_Exception, "Argument 1 must be a 1D double array of length nvar");
		return NULL; 
	}
	
	if (CUTEst_ncon>0) {
		/* Check if v is double and of correct dimension */
		if (arg2!=NULL) {
			/* Check if v is double and of correct dimension */
			if (!(PyArray_Check(arg2) && PyArray_ISFLOAT(arg2) && PyArray_TYPE(arg2)==NPY_DOUBLE && PyArray_NDIM(arg2)==1 && PyArray_DIM(arg2, 0)==CUTEst_ncon)) {
				PyErr_SetString(PyExc_Exception, "Argument 2 must be a 1D double array of length ncon");
				return NULL; 
			}
		} else {
			PyErr_SetString(PyExc_Exception, "Argument 2 must be specified for constrained problems.");
			return NULL; 
		}
	}

#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: preparing for evaluation\n");
#endif
	x=(npy_double *)PyArray_DATA(arg1);
	si=(npy_integer *)malloc(CUTEst_nnzh*sizeof(npy_integer));
	sj=(npy_integer *)malloc(CUTEst_nnzh*sizeof(npy_integer));
	sv=(npy_double *)malloc(CUTEst_nnzh*sizeof(npy_double)); 
		
	if (CUTEst_ncon>0) {
		v=(npy_double *)PyArray_DATA(arg2);
		nnzjplusn=CUTEst_nnzj+CUTEst_nvar;
		sji=(npy_integer *)malloc(nnzjplusn*sizeof(npy_integer));
		sjfi=(npy_integer *)malloc(nnzjplusn*sizeof(npy_integer));
		sjv=(npy_double *)malloc(nnzjplusn*sizeof(npy_double)); 
		
#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: calling CSGRSH\n");
#endif
		if (lagrangian) {
			CUTEST_csgrsh((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, v, &somethingTrue,
				 	(integer *)&nnzjplusno, (integer *)&nnzjplusn, sjv, (integer *)sji, (integer *)sjfi, 
					(integer *)&nnzho, (integer *)&CUTEst_nnzh, sv, (integer *)si, (integer *)sj);     
		} else { 
			CUTEST_csgrsh((integer *)&status, (integer *)&CUTEst_nvar, (integer *)&CUTEst_ncon, x, v, &somethingFalse,
					(integer *)&nnzjplusno, (integer *)&nnzjplusn, sjv, (integer *)sji, (integer *)sjfi, 
					(integer *)&nnzho, (integer *)&CUTEst_nnzh, sv, (integer *)si, (integer *)sj);
		}
		
		extract_sparse_gradient_jacobian(nnzjplusno, sji, sjfi, sjv, &Mgi, &Mgv, &MJi, &MJfi, &MJv); 
		
		/* Free temporary storage - Jacobian */
		free(sji);
		free(sjfi);
		free(sjv);
	} else {
		dims[0]=CUTEst_nvar; 
		Mg=PyArray_SimpleNew(1, dims, PyArray_DOUBLE); 
		g=(npy_double *)PyArray_DATA(Mg);
		
#ifdef PYDEBUG
		fprintf(df, "PyCUTEst: calling UGRSH\n");
#endif
		CUTEST_ugrsh((integer *)&status, (integer *)&CUTEst_nvar, x, g, 
			     (integer *)&nnzho, (integer *)&CUTEst_nnzh, sv, (integer *)si, (integer *)sj);
	}
	
	extract_sparse_hessian(nnzho, si, sj, sv, &MHi, &MHj, &MHv); 
	
	/* Free temporary storage - Hessian */
	free(si);
	free(sj);
	free(sv);
		
	if (CUTEst_ncon>0) {
		return decRefTuple(PyTuple_Pack(8, Mgi, Mgv, MJi, MJfi, MJv, MHi, MHj, MHv)); 
	} else {
		return decRefTuple(PyTuple_Pack(4, Mg, MHi, MHj, MHv)); 
	}
}


static char cutest_report_doc[]=
"Reports usage statistics.\n"
"\n"
"stat=report()\n"
"\n"
"Output\n"
"stat -- dictionary with the usage statistics\n"
"\n"
"The usage statistics dictionary has the following members:\n"
"f      -- number of objective evaluations\n"
"g      -- number of objective gradient evaluations\n"
"H      -- number of objective Hessian evaluations\n"
"Hprod  -- number of Hessian multiplications with a vector\n"
"tsetup -- CPU time used in setup\n"
"trun   -- CPU time used in run\n"
"\n"
"For constrained problems the following additional members are available\n"
"c      -- number of constraint evaluations\n"
"cg     -- number of constraint gradient evaluations\n"
"cH     -- number of constraint Hessian evaluations\n"
"\n"
"CUTEst tools used: CREPRT, UREPRT\n";

static PyObject *cutest_report(PyObject *self, PyObject *args) {
	doublereal calls[7], time[2];
	PyObject *dict;
	
#ifdef PYDEBUG
	fprintf(df, "PyCUTEst: checking arguments\n");
#endif
	if (!check_setup())
		return NULL; 
		
	if (PyObject_Length(args)!=0) {
		PyErr_SetString(PyExc_Exception, "report() takes no arguments");
		return NULL; 
	}
	
	if (CUTEst_ncon>0) 
		CUTEST_creport((integer *)&status, calls, time); 
	else
		CUTEST_ureport((integer *)&status, calls, time); 
		
	dict=PyDict_New();
	PyDict_SetItemString(dict, "f", PyLong_FromLong((long)(calls[0]))); 
	PyDict_SetItemString(dict, "g", PyLong_FromLong((long)(calls[1]))); 
	PyDict_SetItemString(dict, "H", PyLong_FromLong((long)(calls[2]))); 
	PyDict_SetItemString(dict, "Hprod", PyLong_FromLong((long)(calls[3]))); 
	if (CUTEst_ncon>0) {
		PyDict_SetItemString(dict, "c", PyLong_FromLong((long)(calls[4]))); 
		PyDict_SetItemString(dict, "cg", PyLong_FromLong((long)(calls[5]))); 
		PyDict_SetItemString(dict, "cH", PyLong_FromLong((long)(calls[6]))); 
	}
	PyDict_SetItemString(dict, "tsetup", PyFloat_FromDouble((long)(time[0]))); 
	PyDict_SetItemString(dict, "trun", PyFloat_FromDouble((long)(time[1]))); 
	
	return decRefDict(dict); 
}

/* Methods table */
static PyMethodDef _methods[] = {
	{"_dims", cutest__dims, METH_VARARGS, cutest__dims_doc},
	{"_setup", cutest__setup, METH_VARARGS, cutest__setup_doc},
	{"_varnames", cutest__varnames, METH_VARARGS, cutest__varnames_doc},
	{"_connames", cutest__connames, METH_VARARGS, cutest__connames_doc},
	{"objcons", cutest_objcons, METH_VARARGS, cutest_objcons_doc},
	{"obj", cutest_obj, METH_VARARGS, cutest_obj_doc},
	{"cons", cutest_cons, METH_VARARGS, cutest_cons_doc},
	{"lagjac", cutest_lagjac, METH_VARARGS, cutest_lagjac_doc},
	{"jprod", cutest_jprod, METH_VARARGS, cutest_jprod_doc},
	{"hess", cutest_hess, METH_VARARGS, cutest_hess_doc},
	{"ihess", cutest_ihess, METH_VARARGS, cutest_ihess_doc},
	{"hprod", cutest_hprod, METH_VARARGS, cutest_hprod_doc},
	{"gradhess", cutest_gradhess, METH_VARARGS, cutest_gradhess_doc},
	{"_scons", cutest__scons, METH_VARARGS, cutest__scons_doc},
	{"_slagjac", cutest__slagjac, METH_VARARGS, cutest__slagjac_doc},
	{"_sphess", cutest__sphess, METH_VARARGS, cutest__sphess_doc},
	{"_isphess", cutest__isphess, METH_VARARGS, cutest__isphess_doc},
	{"_gradsphess", cutest__gradsphess, METH_VARARGS, cutest__gradsphess_doc},
	{"report", cutest_report, METH_VARARGS, cutest_report_doc}, 
	{NULL, NULL}     /* Marks the end of this structure */
};

/* Module initialization 
   Module name must be _rawfile in compile and link */

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_pycutestitf",
    NULL,
    -1,
    _methods,
    NULL,
    NULL,
    NULL,
    NULL
};

__declspec(dllexport)
PyObject *PyInit__pycutestitf(void)  {
        PyObject *m, *d;
        m = PyModule_Create(&moduledef);
	import_array();  /* Must be present for NumPy.  Called first after above line. */
        return m;
}

#ifdef __cplusplus
}
#endif
