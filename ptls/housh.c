/* housh.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int housh(doublereal *dummy, integer *k, integer *j, 
    doublereal *tol, logical *zero, doublereal *s)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal dum1, alfa;


/*     PURPOSE: */

/*     The subroutine HOUSH computes a Householder transformation H, */
/*     H = I - S x UU', that 'mirrors' a vector DUMMY(1,...,K) to the */
/*     J-th unit vector. */

/*     ARGUMENT LIST : */

/*     DUMMY - DOUBLE PRECISION array of DIMENSION (K). */
/*         A row or column vector of a matrix that has to be mirrored */
/*         to the corresponding unit vector EJ = (0,...,1,0,...,0). */
/*         On return, DUMMY contains the U-vector of the transformation */
/*         matrix H = I - S x UU'. */
/*     K - INTEGER. */
/*         The dimension of DUMMY. */
/*     J - INTEGER. */
/*         The transformation preserves the J-th element of DUMMY to */
/*         become zero. All the other elements are transformed to zero. */
/*     TOL - DOUBLE PRECISION. */
/*         If on entry, norm(DUMMY) < TOL, ZERO is put equal to .TRUE. */
/*     ZERO - LOGICAL. */
/*         See the description of TOL. */
/*     S - DOUBLE PRECISION. */
/*         On return, S contains the scalar S of the transformation */
/*         matrix H. */

/*     REFERENCES: */

/*     [1] A. Emami-Naeine and P. Van Dooren, Computation of Zeros of */
/*         Linear Multivariable Systems. */
/*         Automatica, 18,No.4 (1982), 415-430. */

/*     CONTRIBUTOR: P. Van Dooren (Philips Res. Laboratory, Brussels). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --dummy;

    /* Function Body */
    *zero = TRUE_;
    *s = 0.;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
    d__1 = dummy[i__];
    *s += d__1 * d__1;
/* L10: */
    }
    alfa = sqrt(*s);
    if (alfa <= *tol) {
    return 0;
    }
    *zero = FALSE_;
    dum1 = dummy[*j];
    if (dum1 > 0.) {
    alfa = -alfa;
    }
    dummy[*j] = dum1 - alfa;
    *s = 1. / (*s - alfa * dum1);
    return 0;
/* *** Last line of HOUSH ********************************************** */
} /* housh_ */

