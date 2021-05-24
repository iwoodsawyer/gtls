/* nsingv.f -- translated by f2c (version 20100827).
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

integer nsingv(doublereal *q, doublereal *e, integer *k, doublereal *theta, 
    doublereal *tol1, doublereal *tol2)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal r__, t;
    static integer numeig;


/*     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --e;
    --q;

    /* Function Body */
    if (*theta < 0.) {
    ret_val = 0;
    return ret_val;
    }
    t = -(*theta) - *tol1;
    numeig = *k;
    if (abs(q[1]) <= *tol2) {
    r__ = t;
    } else {
    r__ = t - q[1] * (q[1] / t);
    if (r__ > 0.) {
        --numeig;
    }
    }

    i__1 = *k;
    for (j = 2; j <= i__1; ++j) {
    if ((d__1 = e[j], abs(d__1)) <= *tol2) {
        r__ = t;
    } else {
        r__ = t - e[j] * (e[j] / r__);
        if (r__ > 0.) {
        --numeig;
        }
    }
    if ((d__1 = q[j], abs(d__1)) <= *tol2) {
        r__ = t;
    } else {
        r__ = t - q[j] * (q[j] / r__);
        if (r__ > 0.) {
        --numeig;
        }
    }
/* L1: */
    }

    ret_val = numeig;
    return ret_val;
/* *** Last line of NSINGV ********************************************* */
} /* nsingv_ */

