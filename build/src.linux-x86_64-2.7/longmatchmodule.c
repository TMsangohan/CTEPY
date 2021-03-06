/* File: longmatchmodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * See http://cens.ioc.ee/projects/f2py2e/
 * Generation date: Thu Jan 19 11:16:04 2017
 * $Revision:$
 * $Date:$
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
/*need_includes0*/

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *longmatch_error;
static PyObject *longmatch_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
#define PRINTPYOBJERR(obj)\
  fprintf(stderr,"longmatch.error is related to ");\
  PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
  fprintf(stderr,"\n");

#define pyobj_from_double1(v) (PyFloat_FromDouble(v))
#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
  PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
  fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) _F2PYWRAP##F
#else
#define F_WRAPPEDFUNC(f,F) _f2pywrap##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) _F2PYWRAP##F##_
#else
#define F_WRAPPEDFUNC(f,F) _f2pywrap##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) F2PYWRAP##F
#else
#define F_WRAPPEDFUNC(f,F) f2pywrap##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) F2PYWRAP##F##_
#else
#define F_WRAPPEDFUNC(f,F) f2pywrap##f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_WRAPPEDFUNC_US(f,F) F_WRAPPEDFUNC(f##_,F##_)
#else
#define F_WRAPPEDFUNC_US(f,F) F_WRAPPEDFUNC(f,F)
#endif

/* New SciPy */
#define TRYPYARRAYTEMPLATECHAR case NPY_STRING: *(char *)(PyArray_DATA(arr))=*v; break;
#define TRYPYARRAYTEMPLATELONG case NPY_LONG: *(long *)(PyArray_DATA(arr))=*v; break;
#define TRYPYARRAYTEMPLATEOBJECT case NPY_OBJECT: (PyArray_DESCR(arr)->f->setitem)(pyobj_from_ ## ctype ## 1(*v),PyArray_DATA(arr)); break;

#define TRYPYARRAYTEMPLATE(ctype,typecode) \
        PyArrayObject *arr = NULL;\
        if (!obj) return -2;\
        if (!PyArray_Check(obj)) return -1;\
        if (!(arr=(PyArrayObject *)obj)) {fprintf(stderr,"TRYPYARRAYTEMPLATE:");PRINTPYOBJERR(obj);return 0;}\
        if (PyArray_DESCR(arr)->type==typecode)  {*(ctype *)(PyArray_DATA(arr))=*v; return 1;}\
        switch (PyArray_TYPE(arr)) {\
                case NPY_DOUBLE: *(double *)(PyArray_DATA(arr))=*v; break;\
                case NPY_INT: *(int *)(PyArray_DATA(arr))=*v; break;\
                case NPY_LONG: *(long *)(PyArray_DATA(arr))=*v; break;\
                case NPY_FLOAT: *(float *)(PyArray_DATA(arr))=*v; break;\
                case NPY_CDOUBLE: *(double *)(PyArray_DATA(arr))=*v; break;\
                case NPY_CFLOAT: *(float *)(PyArray_DATA(arr))=*v; break;\
                case NPY_BOOL: *(npy_bool *)(PyArray_DATA(arr))=(*v!=0); break;\
                case NPY_UBYTE: *(unsigned char *)(PyArray_DATA(arr))=*v; break;\
                case NPY_BYTE: *(signed char *)(PyArray_DATA(arr))=*v; break;\
                case NPY_SHORT: *(short *)(PyArray_DATA(arr))=*v; break;\
                case NPY_USHORT: *(npy_ushort *)(PyArray_DATA(arr))=*v; break;\
                case NPY_UINT: *(npy_uint *)(PyArray_DATA(arr))=*v; break;\
                case NPY_ULONG: *(npy_ulong *)(PyArray_DATA(arr))=*v; break;\
                case NPY_LONGLONG: *(npy_longlong *)(PyArray_DATA(arr))=*v; break;\
                case NPY_ULONGLONG: *(npy_ulonglong *)(PyArray_DATA(arr))=*v; break;\
                case NPY_LONGDOUBLE: *(npy_longdouble *)(PyArray_DATA(arr))=*v; break;\
                case NPY_CLONGDOUBLE: *(npy_longdouble *)(PyArray_DATA(arr))=*v; break;\
                case NPY_OBJECT: (PyArray_DESCR(arr)->f->setitem)(pyobj_from_ ## ctype ## 1(*v),PyArray_DATA(arr), arr); break;\
        default: return -2;\
        };\
        return 1


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int try_pyarr_from_double(PyObject* obj,double* v) {
  TRYPYARRAYTEMPLATE(double,'d');
}

static int double_from_pyobj(double* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyFloat_Check(obj)) {
#ifdef __sgi
    *v = PyFloat_AsDouble(obj);
#else
    *v = PyFloat_AS_DOUBLE(obj);
#endif
    return 1;
  }
  tmp = PyNumber_Float(obj);
  if (tmp) {
#ifdef __sgi
    *v = PyFloat_AsDouble(tmp);
#else
    *v = PyFloat_AS_DOUBLE(tmp);
#endif
    Py_DECREF(tmp);
    return 1;
  }
  if (PyComplex_Check(obj))
    tmp = PyObject_GetAttrString(obj,"real");
  else if (PyString_Check(obj) || PyUnicode_Check(obj))
    /*pass*/;
  else if (PySequence_Check(obj))
    tmp = PySequence_GetItem(obj,0);
  if (tmp) {
    PyErr_Clear();
    if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = longmatch_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}

static int int_from_pyobj(int* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyInt_Check(obj)) {
    *v = (int)PyInt_AS_LONG(obj);
    return 1;
  }
  tmp = PyNumber_Int(obj);
  if (tmp) {
    *v = PyInt_AS_LONG(tmp);
    Py_DECREF(tmp);
    return 1;
  }
  if (PyComplex_Check(obj))
    tmp = PyObject_GetAttrString(obj,"real");
  else if (PyString_Check(obj) || PyUnicode_Check(obj))
    /*pass*/;
  else if (PySequence_Check(obj))
    tmp = PySequence_GetItem(obj,0);
  if (tmp) {
    PyErr_Clear();
    if (int_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = longmatch_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC(samplelongmatched,SAMPLELONGMATCHED)(int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);
extern void F_WRAPPEDFUNC(ham,HAM)(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);
extern void F_WRAPPEDFUNC(hammax,HAMMAX)(double*,double*,double*,double*,double*,double*,double*,double*);
extern void F_WRAPPEDFUNC(hamzero,HAMZERO)(double*,double*,double*,double*,double*,double*,double*,double*);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/***************************** samplelongmatched *****************************/
static char doc_f2py_rout_longmatch_samplelongmatched[] = "\
t,pt,sigs = samplelongmatched(np,fnharm,fnharm2,ham1sig,tauhat,ptmax,v00,omega0,hammax,clight,tcoeff,vrf1,vrf2,iseed)\n\nWrapper for ``samplelongmatched``.\
\n\nParameters\n----------\n"
"np : input int\n"
"fnharm : input float\n"
"fnharm2 : input float\n"
"ham1sig : input float\n"
"tauhat : input float\n"
"ptmax : input float\n"
"v00 : input float\n"
"omega0 : input float\n"
"hammax : input float\n"
"clight : input float\n"
"tcoeff : input float\n"
"vrf1 : input float\n"
"vrf2 : input float\n"
"iseed : in/output rank-0 array(float,'d')\n"
"\nReturns\n-------\n"
"t : rank-1 array('d') with bounds (np)\n"
"pt : rank-1 array('d') with bounds (np)\n"
"sigs : float";
/* extern void F_FUNC(samplelongmatched,SAMPLELONGMATCHED)(int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*); */
static PyObject *f2py_rout_longmatch_samplelongmatched(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  int np = 0;
  PyObject *np_capi = Py_None;
  double *t = NULL;
  npy_intp t_Dims[1] = {-1};
  const int t_Rank = 1;
  PyArrayObject *capi_t_tmp = NULL;
  int capi_t_intent = 0;
  double *pt = NULL;
  npy_intp pt_Dims[1] = {-1};
  const int pt_Rank = 1;
  PyArrayObject *capi_pt_tmp = NULL;
  int capi_pt_intent = 0;
  double fnharm = 0;
  PyObject *fnharm_capi = Py_None;
  double fnharm2 = 0;
  PyObject *fnharm2_capi = Py_None;
  double ham1sig = 0;
  PyObject *ham1sig_capi = Py_None;
  double sigs = 0;
  double tauhat = 0;
  PyObject *tauhat_capi = Py_None;
  double ptmax = 0;
  PyObject *ptmax_capi = Py_None;
  double v00 = 0;
  PyObject *v00_capi = Py_None;
  double omega0 = 0;
  PyObject *omega0_capi = Py_None;
  double hammax = 0;
  PyObject *hammax_capi = Py_None;
  double clight = 0;
  PyObject *clight_capi = Py_None;
  double tcoeff = 0;
  PyObject *tcoeff_capi = Py_None;
  double vrf1 = 0;
  PyObject *vrf1_capi = Py_None;
  double vrf2 = 0;
  PyObject *vrf2_capi = Py_None;
  double iseed = 0;
  PyObject *iseed_capi = Py_None;
  static char *capi_kwlist[] = {"np","fnharm","fnharm2","ham1sig","tauhat","ptmax","v00","omega0","hammax","clight","tcoeff","vrf1","vrf2","iseed",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOOOOOOOOO:longmatch.samplelongmatched",\
    capi_kwlist,&np_capi,&fnharm_capi,&fnharm2_capi,&ham1sig_capi,&tauhat_capi,&ptmax_capi,&v00_capi,&omega0_capi,&hammax_capi,&clight_capi,&tcoeff_capi,&vrf1_capi,&vrf2_capi,&iseed_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable v00 */
    f2py_success = double_from_pyobj(&v00,v00_capi,"longmatch.samplelongmatched() 7th argument (v00) can't be converted to double");
  if (f2py_success) {
  /* Processing variable hammax */
    f2py_success = double_from_pyobj(&hammax,hammax_capi,"longmatch.samplelongmatched() 9th argument (hammax) can't be converted to double");
  if (f2py_success) {
  /* Processing variable clight */
    f2py_success = double_from_pyobj(&clight,clight_capi,"longmatch.samplelongmatched() 10th argument (clight) can't be converted to double");
  if (f2py_success) {
  /* Processing variable fnharm */
    f2py_success = double_from_pyobj(&fnharm,fnharm_capi,"longmatch.samplelongmatched() 2nd argument (fnharm) can't be converted to double");
  if (f2py_success) {
  /* Processing variable omega0 */
    f2py_success = double_from_pyobj(&omega0,omega0_capi,"longmatch.samplelongmatched() 8th argument (omega0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable np */
    f2py_success = int_from_pyobj(&np,np_capi,"longmatch.samplelongmatched() 1st argument (np) can't be converted to int");
  if (f2py_success) {
  /* Processing variable sigs */
  /* Processing variable iseed */
    f2py_success = double_from_pyobj(&iseed,iseed_capi,"longmatch.samplelongmatched() 14th argument (iseed) can't be converted to double");
  if (f2py_success) {
  /* Processing variable ptmax */
    f2py_success = double_from_pyobj(&ptmax,ptmax_capi,"longmatch.samplelongmatched() 6th argument (ptmax) can't be converted to double");
  if (f2py_success) {
  /* Processing variable fnharm2 */
    f2py_success = double_from_pyobj(&fnharm2,fnharm2_capi,"longmatch.samplelongmatched() 3rd argument (fnharm2) can't be converted to double");
  if (f2py_success) {
  /* Processing variable tauhat */
    f2py_success = double_from_pyobj(&tauhat,tauhat_capi,"longmatch.samplelongmatched() 5th argument (tauhat) can't be converted to double");
  if (f2py_success) {
  /* Processing variable vrf2 */
    f2py_success = double_from_pyobj(&vrf2,vrf2_capi,"longmatch.samplelongmatched() 13rd argument (vrf2) can't be converted to double");
  if (f2py_success) {
  /* Processing variable ham1sig */
    f2py_success = double_from_pyobj(&ham1sig,ham1sig_capi,"longmatch.samplelongmatched() 4th argument (ham1sig) can't be converted to double");
  if (f2py_success) {
  /* Processing variable vrf1 */
    f2py_success = double_from_pyobj(&vrf1,vrf1_capi,"longmatch.samplelongmatched() 12nd argument (vrf1) can't be converted to double");
  if (f2py_success) {
  /* Processing variable tcoeff */
    f2py_success = double_from_pyobj(&tcoeff,tcoeff_capi,"longmatch.samplelongmatched() 11st argument (tcoeff) can't be converted to double");
  if (f2py_success) {
  /* Processing variable t */
  t_Dims[0]=np;
  capi_t_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_t_tmp = array_from_pyobj(NPY_DOUBLE,t_Dims,t_Rank,capi_t_intent,Py_None);
  if (capi_t_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(longmatch_error,"failed in converting hidden `t' of longmatch.samplelongmatched to C/Fortran array" );
  } else {
    t = (double *)(PyArray_DATA(capi_t_tmp));

  /* Processing variable pt */
  pt_Dims[0]=np;
  capi_pt_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_pt_tmp = array_from_pyobj(NPY_DOUBLE,pt_Dims,pt_Rank,capi_pt_intent,Py_None);
  if (capi_pt_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(longmatch_error,"failed in converting hidden `pt' of longmatch.samplelongmatched to C/Fortran array" );
  } else {
    pt = (double *)(PyArray_DATA(capi_pt_tmp));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(&np,t,pt,&fnharm,&fnharm2,&ham1sig,&sigs,&tauhat,&ptmax,&v00,&omega0,&hammax,&clight,&tcoeff,&vrf1,&vrf2,&iseed);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
  f2py_success = try_pyarr_from_double(iseed_capi,&iseed);
  if (f2py_success) {
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("NNd",capi_t_tmp,capi_pt_tmp,sigs);
/*closepyobjfrom*/
  } /*if (f2py_success) of iseed pyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  }  /*if (capi_pt_tmp == NULL) ... else of pt*/
  /* End of cleaning variable pt */
  }  /*if (capi_t_tmp == NULL) ... else of t*/
  /* End of cleaning variable t */
  } /*if (f2py_success) of tcoeff*/
  /* End of cleaning variable tcoeff */
  } /*if (f2py_success) of vrf1*/
  /* End of cleaning variable vrf1 */
  } /*if (f2py_success) of ham1sig*/
  /* End of cleaning variable ham1sig */
  } /*if (f2py_success) of vrf2*/
  /* End of cleaning variable vrf2 */
  } /*if (f2py_success) of tauhat*/
  /* End of cleaning variable tauhat */
  } /*if (f2py_success) of fnharm2*/
  /* End of cleaning variable fnharm2 */
  } /*if (f2py_success) of ptmax*/
  /* End of cleaning variable ptmax */
  } /*if (f2py_success) of iseed*/
  /* End of cleaning variable iseed */
  /* End of cleaning variable sigs */
  } /*if (f2py_success) of np*/
  /* End of cleaning variable np */
  } /*if (f2py_success) of omega0*/
  /* End of cleaning variable omega0 */
  } /*if (f2py_success) of fnharm*/
  /* End of cleaning variable fnharm */
  } /*if (f2py_success) of clight*/
  /* End of cleaning variable clight */
  } /*if (f2py_success) of hammax*/
  /* End of cleaning variable hammax */
  } /*if (f2py_success) of v00*/
  /* End of cleaning variable v00 */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/************************** end of samplelongmatched **************************/

/************************************ ham ************************************/
static char doc_f2py_rout_longmatch_ham[] = "\
ham = ham(t,pt,tcoeff,omega0,v00,fnharm,fnharm2,vrf1,vrf2)\n\nWrapper for ``ham``.\
\n\nParameters\n----------\n"
"t : input float\n"
"pt : input float\n"
"tcoeff : input float\n"
"omega0 : input float\n"
"v00 : input float\n"
"fnharm : input float\n"
"fnharm2 : input float\n"
"vrf1 : input float\n"
"vrf2 : input float\n"
"\nReturns\n-------\n"
"ham : float";
/* extern void F_WRAPPEDFUNC(ham,HAM)(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*); */
static PyObject *f2py_rout_longmatch_ham(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double ham = 0;
  double t = 0;
  PyObject *t_capi = Py_None;
  double pt = 0;
  PyObject *pt_capi = Py_None;
  double tcoeff = 0;
  PyObject *tcoeff_capi = Py_None;
  double omega0 = 0;
  PyObject *omega0_capi = Py_None;
  double v00 = 0;
  PyObject *v00_capi = Py_None;
  double fnharm = 0;
  PyObject *fnharm_capi = Py_None;
  double fnharm2 = 0;
  PyObject *fnharm2_capi = Py_None;
  double vrf1 = 0;
  PyObject *vrf1_capi = Py_None;
  double vrf2 = 0;
  PyObject *vrf2_capi = Py_None;
  static char *capi_kwlist[] = {"t","pt","tcoeff","omega0","v00","fnharm","fnharm2","vrf1","vrf2",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOOOO:longmatch.ham",\
    capi_kwlist,&t_capi,&pt_capi,&tcoeff_capi,&omega0_capi,&v00_capi,&fnharm_capi,&fnharm2_capi,&vrf1_capi,&vrf2_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable v00 */
    f2py_success = double_from_pyobj(&v00,v00_capi,"longmatch.ham() 5th argument (v00) can't be converted to double");
  if (f2py_success) {
  /* Processing variable pt */
    f2py_success = double_from_pyobj(&pt,pt_capi,"longmatch.ham() 2nd argument (pt) can't be converted to double");
  if (f2py_success) {
  /* Processing variable fnharm2 */
    f2py_success = double_from_pyobj(&fnharm2,fnharm2_capi,"longmatch.ham() 7th argument (fnharm2) can't be converted to double");
  if (f2py_success) {
  /* Processing variable fnharm */
    f2py_success = double_from_pyobj(&fnharm,fnharm_capi,"longmatch.ham() 6th argument (fnharm) can't be converted to double");
  if (f2py_success) {
  /* Processing variable tcoeff */
    f2py_success = double_from_pyobj(&tcoeff,tcoeff_capi,"longmatch.ham() 3rd argument (tcoeff) can't be converted to double");
  if (f2py_success) {
  /* Processing variable vrf2 */
    f2py_success = double_from_pyobj(&vrf2,vrf2_capi,"longmatch.ham() 9th argument (vrf2) can't be converted to double");
  if (f2py_success) {
  /* Processing variable t */
    f2py_success = double_from_pyobj(&t,t_capi,"longmatch.ham() 1st argument (t) can't be converted to double");
  if (f2py_success) {
  /* Processing variable vrf1 */
    f2py_success = double_from_pyobj(&vrf1,vrf1_capi,"longmatch.ham() 8th argument (vrf1) can't be converted to double");
  if (f2py_success) {
  /* Processing variable omega0 */
    f2py_success = double_from_pyobj(&omega0,omega0_capi,"longmatch.ham() 4th argument (omega0) can't be converted to double");
  if (f2py_success) {
  /* Processing variable ham */
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
  (*f2py_func)(&ham,&t,&pt,&tcoeff,&omega0,&v00,&fnharm,&fnharm2,&vrf1,&vrf2);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("d",ham);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  /* End of cleaning variable ham */
  } /*if (f2py_success) of omega0*/
  /* End of cleaning variable omega0 */
  } /*if (f2py_success) of vrf1*/
  /* End of cleaning variable vrf1 */
  } /*if (f2py_success) of t*/
  /* End of cleaning variable t */
  } /*if (f2py_success) of vrf2*/
  /* End of cleaning variable vrf2 */
  } /*if (f2py_success) of tcoeff*/
  /* End of cleaning variable tcoeff */
  } /*if (f2py_success) of fnharm*/
  /* End of cleaning variable fnharm */
  } /*if (f2py_success) of fnharm2*/
  /* End of cleaning variable fnharm2 */
  } /*if (f2py_success) of pt*/
  /* End of cleaning variable pt */
  } /*if (f2py_success) of v00*/
  /* End of cleaning variable v00 */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/********************************* end of ham *********************************/

/*********************************** hammax ***********************************/
static char doc_f2py_rout_longmatch_hammax[] = "\
hammax = hammax(tauhat,omega0,v00,fnharm,fnharm2,vrf1,vrf2)\n\nWrapper for ``hammax``.\
\n\nParameters\n----------\n"
"tauhat : input float\n"
"omega0 : input float\n"
"v00 : input float\n"
"fnharm : input float\n"
"fnharm2 : input float\n"
"vrf1 : input float\n"
"vrf2 : input float\n"
"\nReturns\n-------\n"
"hammax : float";
/* extern void F_WRAPPEDFUNC(hammax,HAMMAX)(double*,double*,double*,double*,double*,double*,double*,double*); */
static PyObject *f2py_rout_longmatch_hammax(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(double*,double*,double*,double*,double*,double*,double*,double*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double hammax = 0;
  double tauhat = 0;
  PyObject *tauhat_capi = Py_None;
  double omega0 = 0;
  PyObject *omega0_capi = Py_None;
  double v00 = 0;
  PyObject *v00_capi = Py_None;
  double fnharm = 0;
  PyObject *fnharm_capi = Py_None;
  double fnharm2 = 0;
  PyObject *fnharm2_capi = Py_None;
  double vrf1 = 0;
  PyObject *vrf1_capi = Py_None;
  double vrf2 = 0;
  PyObject *vrf2_capi = Py_None;
  static char *capi_kwlist[] = {"tauhat","omega0","v00","fnharm","fnharm2","vrf1","vrf2",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOO:longmatch.hammax",\
    capi_kwlist,&tauhat_capi,&omega0_capi,&v00_capi,&fnharm_capi,&fnharm2_capi,&vrf1_capi,&vrf2_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable v00 */
    f2py_success = double_from_pyobj(&v00,v00_capi,"longmatch.hammax() 3rd argument (v00) can't be converted to double");
  if (f2py_success) {
  /* Processing variable hammax */
  /* Processing variable fnharm2 */
    f2py_success = double_from_pyobj(&fnharm2,fnharm2_capi,"longmatch.hammax() 5th argument (fnharm2) can't be converted to double");
  if (f2py_success) {
  /* Processing variable fnharm */
    f2py_success = double_from_pyobj(&fnharm,fnharm_capi,"longmatch.hammax() 4th argument (fnharm) can't be converted to double");
  if (f2py_success) {
  /* Processing variable tauhat */
    f2py_success = double_from_pyobj(&tauhat,tauhat_capi,"longmatch.hammax() 1st argument (tauhat) can't be converted to double");
  if (f2py_success) {
  /* Processing variable vrf2 */
    f2py_success = double_from_pyobj(&vrf2,vrf2_capi,"longmatch.hammax() 7th argument (vrf2) can't be converted to double");
  if (f2py_success) {
  /* Processing variable vrf1 */
    f2py_success = double_from_pyobj(&vrf1,vrf1_capi,"longmatch.hammax() 6th argument (vrf1) can't be converted to double");
  if (f2py_success) {
  /* Processing variable omega0 */
    f2py_success = double_from_pyobj(&omega0,omega0_capi,"longmatch.hammax() 2nd argument (omega0) can't be converted to double");
  if (f2py_success) {
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
  (*f2py_func)(&hammax,&tauhat,&omega0,&v00,&fnharm,&fnharm2,&vrf1,&vrf2);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("d",hammax);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  } /*if (f2py_success) of omega0*/
  /* End of cleaning variable omega0 */
  } /*if (f2py_success) of vrf1*/
  /* End of cleaning variable vrf1 */
  } /*if (f2py_success) of vrf2*/
  /* End of cleaning variable vrf2 */
  } /*if (f2py_success) of tauhat*/
  /* End of cleaning variable tauhat */
  } /*if (f2py_success) of fnharm*/
  /* End of cleaning variable fnharm */
  } /*if (f2py_success) of fnharm2*/
  /* End of cleaning variable fnharm2 */
  /* End of cleaning variable hammax */
  } /*if (f2py_success) of v00*/
  /* End of cleaning variable v00 */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/******************************* end of hammax *******************************/

/********************************** hamzero **********************************/
static char doc_f2py_rout_longmatch_hamzero[] = "\
hamzero = hamzero(phizero,omega0,v00,fnharm,fnharm2,vrf1,vrf2)\n\nWrapper for ``hamzero``.\
\n\nParameters\n----------\n"
"phizero : input float\n"
"omega0 : input float\n"
"v00 : input float\n"
"fnharm : input float\n"
"fnharm2 : input float\n"
"vrf1 : input float\n"
"vrf2 : input float\n"
"\nReturns\n-------\n"
"hamzero : float";
/* extern void F_WRAPPEDFUNC(hamzero,HAMZERO)(double*,double*,double*,double*,double*,double*,double*,double*); */
static PyObject *f2py_rout_longmatch_hamzero(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(double*,double*,double*,double*,double*,double*,double*,double*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  double hamzero = 0;
  double phizero = 0;
  PyObject *phizero_capi = Py_None;
  double omega0 = 0;
  PyObject *omega0_capi = Py_None;
  double v00 = 0;
  PyObject *v00_capi = Py_None;
  double fnharm = 0;
  PyObject *fnharm_capi = Py_None;
  double fnharm2 = 0;
  PyObject *fnharm2_capi = Py_None;
  double vrf1 = 0;
  PyObject *vrf1_capi = Py_None;
  double vrf2 = 0;
  PyObject *vrf2_capi = Py_None;
  static char *capi_kwlist[] = {"phizero","omega0","v00","fnharm","fnharm2","vrf1","vrf2",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOO:longmatch.hamzero",\
    capi_kwlist,&phizero_capi,&omega0_capi,&v00_capi,&fnharm_capi,&fnharm2_capi,&vrf1_capi,&vrf2_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable phizero */
    f2py_success = double_from_pyobj(&phizero,phizero_capi,"longmatch.hamzero() 1st argument (phizero) can't be converted to double");
  if (f2py_success) {
  /* Processing variable v00 */
    f2py_success = double_from_pyobj(&v00,v00_capi,"longmatch.hamzero() 3rd argument (v00) can't be converted to double");
  if (f2py_success) {
  /* Processing variable hamzero */
  /* Processing variable fnharm2 */
    f2py_success = double_from_pyobj(&fnharm2,fnharm2_capi,"longmatch.hamzero() 5th argument (fnharm2) can't be converted to double");
  if (f2py_success) {
  /* Processing variable fnharm */
    f2py_success = double_from_pyobj(&fnharm,fnharm_capi,"longmatch.hamzero() 4th argument (fnharm) can't be converted to double");
  if (f2py_success) {
  /* Processing variable vrf2 */
    f2py_success = double_from_pyobj(&vrf2,vrf2_capi,"longmatch.hamzero() 7th argument (vrf2) can't be converted to double");
  if (f2py_success) {
  /* Processing variable vrf1 */
    f2py_success = double_from_pyobj(&vrf1,vrf1_capi,"longmatch.hamzero() 6th argument (vrf1) can't be converted to double");
  if (f2py_success) {
  /* Processing variable omega0 */
    f2py_success = double_from_pyobj(&omega0,omega0_capi,"longmatch.hamzero() 2nd argument (omega0) can't be converted to double");
  if (f2py_success) {
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
  (*f2py_func)(&hamzero,&phizero,&omega0,&v00,&fnharm,&fnharm2,&vrf1,&vrf2);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("d",hamzero);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  } /*if (f2py_success) of omega0*/
  /* End of cleaning variable omega0 */
  } /*if (f2py_success) of vrf1*/
  /* End of cleaning variable vrf1 */
  } /*if (f2py_success) of vrf2*/
  /* End of cleaning variable vrf2 */
  } /*if (f2py_success) of fnharm*/
  /* End of cleaning variable fnharm */
  } /*if (f2py_success) of fnharm2*/
  /* End of cleaning variable fnharm2 */
  /* End of cleaning variable hamzero */
  } /*if (f2py_success) of v00*/
  /* End of cleaning variable v00 */
  } /*if (f2py_success) of phizero*/
  /* End of cleaning variable phizero */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/******************************* end of hamzero *******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"samplelongmatched",-1,{{-1}},0,(char *)F_FUNC(samplelongmatched,SAMPLELONGMATCHED),(f2py_init_func)f2py_rout_longmatch_samplelongmatched,doc_f2py_rout_longmatch_samplelongmatched},
  {"ham",-1,{{-1}},0,(char *)F_WRAPPEDFUNC(ham,HAM),(f2py_init_func)f2py_rout_longmatch_ham,doc_f2py_rout_longmatch_ham},
  {"hammax",-1,{{-1}},0,(char *)F_WRAPPEDFUNC(hammax,HAMMAX),(f2py_init_func)f2py_rout_longmatch_hammax,doc_f2py_rout_longmatch_hammax},
  {"hamzero",-1,{{-1}},0,(char *)F_WRAPPEDFUNC(hamzero,HAMZERO),(f2py_init_func)f2py_rout_longmatch_hamzero,doc_f2py_rout_longmatch_hamzero},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "longmatch",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};
#endif

#if PY_VERSION_HEX >= 0x03000000
#define RETVAL m
PyMODINIT_FUNC PyInit_longmatch(void) {
#else
#define RETVAL
PyMODINIT_FUNC initlongmatch(void) {
#endif
  int i;
  PyObject *m,*d, *s;
#if PY_VERSION_HEX >= 0x03000000
  m = longmatch_module = PyModule_Create(&moduledef);
#else
  m = longmatch_module = Py_InitModule("longmatch", f2py_module_methods);
#endif
  Py_TYPE(&PyFortran_Type) = &PyType_Type;
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module longmatch (failed to import numpy)"); return RETVAL;}
  d = PyModule_GetDict(m);
  s = PyString_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
#if PY_VERSION_HEX >= 0x03000000
  s = PyUnicode_FromString(
#else
  s = PyString_FromString(
#endif
    "This module 'longmatch' is auto-generated with f2py (version:2).\nFunctions:\n"
"  t,pt,sigs = samplelongmatched(np,fnharm,fnharm2,ham1sig,tauhat,ptmax,v00,omega0,hammax,clight,tcoeff,vrf1,vrf2,iseed)\n"
"  ham = ham(t,pt,tcoeff,omega0,v00,fnharm,fnharm2,vrf1,vrf2)\n"
"  hammax = hammax(tauhat,omega0,v00,fnharm,fnharm2,vrf1,vrf2)\n"
"  hamzero = hamzero(phizero,omega0,v00,fnharm,fnharm2,vrf1,vrf2)\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  longmatch_error = PyErr_NewException ("longmatch.error", NULL, NULL);
  Py_DECREF(s);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++)
    PyDict_SetItemString(d, f2py_routine_defs[i].name,PyFortranObject_NewAsAttr(&f2py_routine_defs[i]));


    {
      extern double F_FUNC(ham,HAM)(void);
      PyObject* o = PyDict_GetItemString(d,"ham");
      PyObject_SetAttrString(o,"_cpointer", F2PyCapsule_FromVoidPtr((void*)F_FUNC(ham,HAM),NULL));
#if PY_VERSION_HEX >= 0x03000000
      PyObject_SetAttrString(o,"__name__", PyUnicode_FromString("ham"));
#else
      PyObject_SetAttrString(o,"__name__", PyString_FromString("ham"));
#endif
    }
    

    {
      extern double F_FUNC(hammax,HAMMAX)(void);
      PyObject* o = PyDict_GetItemString(d,"hammax");
      PyObject_SetAttrString(o,"_cpointer", F2PyCapsule_FromVoidPtr((void*)F_FUNC(hammax,HAMMAX),NULL));
#if PY_VERSION_HEX >= 0x03000000
      PyObject_SetAttrString(o,"__name__", PyUnicode_FromString("hammax"));
#else
      PyObject_SetAttrString(o,"__name__", PyString_FromString("hammax"));
#endif
    }
    

    {
      extern double F_FUNC(hamzero,HAMZERO)(void);
      PyObject* o = PyDict_GetItemString(d,"hamzero");
      PyObject_SetAttrString(o,"_cpointer", F2PyCapsule_FromVoidPtr((void*)F_FUNC(hamzero,HAMZERO),NULL));
#if PY_VERSION_HEX >= 0x03000000
      PyObject_SetAttrString(o,"__name__", PyUnicode_FromString("hamzero"));
#else
      PyObject_SetAttrString(o,"__name__", PyString_FromString("hamzero"));
#endif
    }
    
/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"longmatch");
#endif

  return RETVAL;
}
#ifdef __cplusplus
}
#endif
