/* init.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int init(doublereal *x, integer *ldx, integer *m, integer *n)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;


/*     PURPOSE: */

/*     The subroutine INIT initializes an M by N matrix X with the M by N */
/*     identity matrix, characterized by unit diagonal entries and zero */
/*     off-diagonal elements. */

/*     ARGUMENT LIST: */

/*     X - DOUBLE PRECISION array of DIMENSION (LDX,N) */
/*         On return, X contains the M by N identity matrix. */
/*     LDX - INTEGER */
/*         LDX is the leading dimension of the array X (LDX >= M). */
/*     M - INTEGER */
/*         M is the number of rows of the matrix X. */
/*     N - INTEGER */
/*         N is the number of columns of the matrix X. */

/*     CONTRIBUTOR: S. Van Huffel, (ESAT Laboratory, KU Leuven). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    i__2 = *m;
    for (i__ = 1; i__ <= i__2; ++i__) {
        x[i__ + j * x_dim1] = 0.;
/* L10: */
    }
    x[j + j * x_dim1] = 1.;
/* L20: */
    }
    return 0;
/* *** Last line of INIT *********************************************** */
} /* init_ */

