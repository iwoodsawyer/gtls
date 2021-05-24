/* bidiag.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int bidiag(doublereal *x, integer *ldx, integer *n, integer 
    *p, doublereal *q, doublereal *e, doublereal *wrk)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, l, m;
    static doublereal t;
    static integer lu, lp1, pp1, nct, npp, nrt;
    extern doublereal ddot(integer *, doublereal *, integer *, doublereal *, 
        integer *), dnrm2(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal(integer *, doublereal *, doublereal *, 
        integer *), daxpy(integer *, doublereal *, doublereal *, integer 
        *, doublereal *, integer *);


/*     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. External Subroutines/Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

/*     Reduce X to bidiagonal form, storing the diagonal elements in Q */
/*     and the superdiagonal elements in E. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --wrk;
    --q;
    --e;

    /* Function Body */
    npp = *n + *p;
    pp1 = *p + 1;
/* Computing MIN */
    i__1 = *n - 1;
    nct = min(i__1,*p);
/* Computing MAX */
/* Computing MIN */
    i__3 = *p - 2;
    i__1 = 0, i__2 = min(i__3,*n);
    nrt = max(i__1,i__2);
    lu = max(nct,nrt);
    i__1 = lu;
    for (l = 1; l <= i__1; ++l) {
    lp1 = l + 1;
    if (l <= nct) {

/*           Compute the transformation for the L-th column and place the */
/*           L-th diagonal in Q(L). */

        i__2 = *n - l + 1;
        q[l] = dnrm2(&i__2, &x[l + l * x_dim1], &c__1);
        if (q[l] != 0.) {
        if (x[l + l * x_dim1] != 0.) {
            q[l] = d_sign(&q[l], &x[l + l * x_dim1]);
        }
        i__2 = *n - l + 1;
        d__1 = 1. / q[l];
        dscal(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
        x[l + l * x_dim1] += 1.;
        }
        q[l] = -q[l];
    }

    i__2 = *p;
    for (j = lp1; j <= i__2; ++j) {
        if (l > nct) {
        goto L1;
        }
        if (q[l] == 0.) {
        goto L1;
        }

/*              Apply the transformation. */

        i__3 = *n - l + 1;
        t = -ddot(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
            c__1) / x[l + l * x_dim1];
        i__3 = *n - l + 1;
        daxpy(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
            c__1);
L1:

/*           Place the L-th row of X into WRK for the subsequent */
/*           calculation of the row transformation. */

        wrk[j] = x[l + j * x_dim1];
/* L2: */
    }
    if (l <= nrt) {

/*           Compute the L-th row transformation and place the L-th */
/*           superdiagonal in E(L). */

        i__2 = *p - l;
        wrk[l] = dnrm2(&i__2, &wrk[lp1], &c__1);
        if (wrk[l] != 0.) {
        if (wrk[lp1] != 0.) {
            wrk[l] = d_sign(&wrk[l], &wrk[lp1]);
        }
        i__2 = *p - l;
        d__1 = 1. / wrk[l];
        dscal(&i__2, &d__1, &wrk[lp1], &c__1);
        wrk[lp1] += 1.;
        }
        wrk[l] = -wrk[l];
        e[lp1] = wrk[l];
        if (lp1 <= *n && wrk[l] != 0.) {

/*              Apply the transformation. */

        i__2 = npp;
        for (i__ = pp1; i__ <= i__2; ++i__) {
            wrk[i__] = 0.;
/* L3: */
        }
        i__2 = *p;
        for (j = lp1; j <= i__2; ++j) {
            i__3 = *n - l;
            daxpy(&i__3, &wrk[j], &x[lp1 + j * x_dim1], &c__1, &wrk[
                pp1], &c__1);
/* L4: */
        }
        i__2 = *p;
        for (j = lp1; j <= i__2; ++j) {
            i__3 = *n - l;
            d__1 = -wrk[j] / wrk[lp1];
            daxpy(&i__3, &d__1, &wrk[pp1], &c__1, &x[lp1 + j * 
                x_dim1], &c__1);
/* L5: */
        }
        }
        i__2 = *p;
        for (i__ = lp1; i__ <= i__2; ++i__) {
        x[l + i__ * x_dim1] = wrk[i__];
/* L6: */
        }
    }
/* L7: */
    }

/*     Set up the final bidiagonal matrix elements, if necessary. */

    e[1] = 0.;
/* Computing MIN */
    i__1 = *p - 1;
    m = min(i__1,*n);
    if (nct < *p) {
    q[*n] = x[*n + *n * x_dim1];
    }
    if (nrt < m) {
    e[*p] = x[*p - 1 + *p * x_dim1];
    }
    return 0;
/* *** Last line of BIDIAG ********************************************* */
} /* bidiag_ */

