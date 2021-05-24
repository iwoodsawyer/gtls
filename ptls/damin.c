/* damin.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
    on Microsoft Windows system, link with libf2c.lib;
    on Linux or Unix systems, link with .../path/to/libf2c.a -lm
    or, if you install libf2c.a in a standard place, with -lf2c -lm
    -- in that order, at the end of the command line, as in
        cc *.o -lf2c -lm
    Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

        http://www.netlib.org/f2c/libf2c.zip
*/

#include "ptls.h"

doublereal damin(doublereal *x, integer *nx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;
    static doublereal dx;
    static integer ix, nincx;


/*     PURPOSE: */

/*     The function DAMIN computes the absolute minimal value of NX */
/*     elements in an array. */
/*     The function returns the value zero if NX < 1. */

/*     ARGUMENT LIST: */

/*     X - DOUBLE PRECISION array of DIMENSION (NX x INCX). */
/*         X is the one-dimensional array of which the absolute minimal */
/*         value of the elements is to be computed. */
/*     NX - INTEGER. */
/*         NX is the number of elements in X to be examined. */
/*     INCX - INTEGER. */
/*         INCX is the increment to be taken in the array X, defining */
/*         the distance between two consecutive elements. */
/*         INCX = 1, if all elements are contiguous in memory. */

/*     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    if (*nx < 1) {
    return ret_val;
    }

    ret_val = abs(x[1]);
    if (*nx == 1) {
    return ret_val;
    }

    ix = *incx + 1;
    nincx = *nx * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = ix; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    dx = (d__1 = x[i__], abs(d__1));
    if (dx < ret_val) {
        ret_val = dx;
    }
/* L1: */
    }
    return ret_val;
/* *** Last line of DAMIN ********************************************** */
} /* damin_ */

