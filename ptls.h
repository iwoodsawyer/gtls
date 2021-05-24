/* Starting from version 7.8, MATLAB BLAS expects ptrdiff_t arguments for integers */
#if MATLAB_VERSION >= 0x0708
#include <stddef.h>
#include <stdlib.h>
#endif 
#include <string.h>
        
/* Define MX_HAS_INTERLEAVED_COMPLEX for version <9.4 */
#ifndef MX_HAS_INTERLEAVED_COMPLEX
#define MX_HAS_INTERLEAVED_COMPLEX 0
#endif

/* Starting from version 7.6, MATLAB BLAS is seperated */
#if MATLAB_VERSION >= 0x0705
#include <blas.h>
#else
#define dcabs1 FORTRAN_WRAPPER(dcabs1)
extern doublereal dcabs1(
        doublereal *z
        );
#endif
#include <lapack.h>
#include "f2c.h"

#define ptls FORTRAN_WRAPPER(ptls)
int ptls(doublereal *, integer *, integer *, integer *,
    integer *, integer *, doublereal *, doublereal *,
    integer *, doublereal *, logical *, doublereal *,
    integer *, logical *, doublereal *, doublereal *,
    integer *, integer *)

#define dtls FORTRAN_WRAPPER(dtls)
int dtls(doublereal *, integer *, integer *, integer*,
    integer *, doublereal *, doublereal *, integer *,
    doublereal *, integer *, doublereal *, doublereal *,
    char *, integer *, integer *)