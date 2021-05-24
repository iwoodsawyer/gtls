/* estim.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int estim(doublereal *q, doublereal *e, integer *n, integer 
    *l, doublereal *theta, doublereal *tol1, doublereal *tol2, integer *iwarn)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal y, z__, h1, h2, th;
    static integer num, numz;
    static doublereal sumz;
    extern doublereal damin(doublereal *, integer *, integer *);
    extern integer nsingv(doublereal *, doublereal *, integer *, doublereal *
        , doublereal *, doublereal *);


/*     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. External Subroutines/Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --e;
    --q;

    /* Function Body */
    *iwarn = 0;
    if (*l < 0 || *l > *n) {
    return 0;
    }

/*     Step 1: Initialization of THETA. */
/*             ----------------------- */
    if (*l == 0) {
    *theta = 0.;
    }
    if (*theta < 0.) {
    if (*l == 1) {

/*           An upper bound which is close if S(N-1) >> S(N): */

        *theta = damin(&q[1], n, &c__1);
    } else {

/*           An experimentally established estimate which is good if */
/*           S(N-L) >> S(N-L+1): */

        *theta = (d__1 = q[*n - *l + 1], abs(d__1));
    }
    }

/*     Step 2: Check quality initial estimate THETA. */
/*             ------------------------------------ */
    num = nsingv(&q[1], &e[1], n, theta, tol1, tol2);
    if (num == *l) {
    return 0;
    }

/*     Step 3: Initialization starting values for bisection method. */
/*             --------------------------------------------------- */
/*     Let S(i), i=1,...,N, be the singular values of J in decreasing */
/*     order. Then, the computed Y and Z will be such that */
/*        (number of S(i) <= Y) < L < (number of S(i) <= Z). */

    if (num < *l) {
    y = *theta;
    th = abs(q[1]);
    z__ = th;
    numz = *n;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
        h1 = (d__1 = e[i__], abs(d__1));
        h2 = (d__1 = q[i__], abs(d__1));
/* Computing MAX */
        d__1 = th + h1, d__2 = h2 + h1;
        sumz = max(d__1,d__2);
        if (sumz > z__) {
        z__ = sumz;
        }
        th = h2;
/* L1: */
    }
    } else {
    z__ = *theta;
    y = 0.;
    numz = num;
    }

/*     Step 4: Bisection method for finding the upper bound on the L */
/*             smallest singular values of the bidiagonal. */
/*             ------------------------------------------ */
/*     A sequence of subintervals [Y,Z] is produced such that */
/*            (number of S(i) <= Y) < L < (number of S(i) <= Z). */
/*     NUM : number of S(i) <= TH, */
/*     NUMZ: number of S(i) <= Z. */

/*     WHILE ((NUM .NE. L) .AND. (Z-Y) .GT. TOL1) DO */
L2:
    if (num != *l && z__ - y > *tol1) {
    th = (y + z__) / 2.;
    num = nsingv(&q[1], &e[1], n, &th, tol1, tol2);
    if (num < *l) {
        y = th;
    } else {
        z__ = th;
        numz = num;
    }
    goto L2;
    }
/*     END WHILE 2 */

/*     If (Z - Y) <= TOL1, then at least two singular values of J lie in */
/*     the interval [Y,Z] within a distance < TOL1 from each other. */
/*     S(N-L) ans S(N-L+1) are then assumed to coincide. */
/*     L is increased, and a warning is given. */

    if (num != *l) {
    *theta = z__;
    *l = numz;
    *iwarn = 1;
    } else {
    *theta = th;
    }
    return 0;
/* *** Last line of ESTIM ********************************************** */
} /* estim_ */

