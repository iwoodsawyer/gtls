/* tr2.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tr2(doublereal *a, integer *lda, doublereal *u, 
    doublereal *s, integer *i1, integer *i2, integer *j1, integer *j2)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal y, inprod;


/*     PURPOSE: */
/*     The subroutine TR2 performs the Householder transformation */
/*     H = I - S x UU' on the columns J1+1 to J1+J2  of matrix A, this */
/*     from rows I1 to I2. */

/*     PARAMETERS: */

/*     A - DOUBLE PRECISION array of DIMENSION (LDA,J1+J2). */
/*         Contains the submatrix of A onto which the Householder */
/*         transformation H is applied. */
/*     LDA - INTEGER. */
/*         The leading dimension of the array A (LDA >= I2). */
/*     U - DOUBLE PRECISION array of DIMENSION (J2). */
/*         Contains the transformation vector of the transformation */
/*         matrix H. */
/*     S - DOUBLE PRECISION. */
/*         Contains the scalar S of the transformation matrix H. */
/*     I1 - INTEGER. */
/*         Contains the first row index of A (see purpose). */
/*     I2 - INTEGER. */
/*         Contains the last row index of A (see purpose, I2 >= I1). */
/*     J1 - INTEGER. */
/*         Contains the first column index of A (see purpose). */
/*     J2 - INTEGER. */
/*         Contains the last column index of A (see purpose). */

/*     CONTRIBUTOR: P. Van Dooren (Philips Res. Laboratory, Brussels). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --u;

    /* Function Body */
    i__1 = *i2;
    for (i__ = *i1; i__ <= i__1; ++i__) {
    inprod = 0.;
    i__2 = *j2;
    for (j = 1; j <= i__2; ++j) {
        inprod += u[j] * a[i__ + (*j1 + j) * a_dim1];
/* L10: */
    }
    y = inprod * *s;
    i__2 = *j2;
    for (j = 1; j <= i__2; ++j) {
        a[i__ + (*j1 + j) * a_dim1] -= u[j] * y;
/* L20: */
    }
    }
    return 0;
/* *** Last line of TR2 ************************************************ */
} /* tr2_ */

