/* dtls.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int dtls(doublereal *c__, integer *ldc, integer *m, integer 
    *n, integer *l, doublereal *s, doublereal *x, integer *ldx, 
    doublereal *wrk, integer *rank, doublereal *tol1, doublereal *tol2, 
    char *comprt, integer *ierr, integer *iwarn)
{
    /* System generated locals */
    integer c_dim1, c_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, p, j1, n1, r1, mc, nj, nl, nl1;
    extern /* Subroutine */ int tr2(doublereal *, integer *, doublereal *, 
        doublereal *, integer *, integer *, integer *, integer *);
    extern doublereal ddot(integer *, doublereal *, integer *, doublereal *, 
        integer *);
    static logical ctol;
    static doublereal temp, smax;
    static logical zero;
    static doublereal smax2;
    extern /* Subroutine */ int dqrdc(doublereal *, integer *, integer *, 
        integer *, doublereal *, integer *, doublereal *, integer *);
    static logical crank;
    extern /* Subroutine */ int dsvdc(doublereal *, integer *, integer *, 
        integer *, doublereal *, doublereal *, doublereal *, integer *, 
        doublereal *, integer *, doublereal *, integer *, integer *), 
        dcopy(integer *, doublereal *, integer *, doublereal *, integer *
        ), daxpy(integer *, doublereal *, doublereal *, integer *, 
        doublereal *, integer *), housh(doublereal *, integer *, integer 
        *, doublereal *, logical *, doublereal *);
    static integer dumar1[1];
    static doublereal dumar2[1] /* was [1][1] */;


/*     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. External Subroutines/Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. Parameters .. */
/*     .. Executable Statements .. */

/*     Determine whether RANK and/or TOL1 is to be computed. */

    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --wrk;
    --s;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    crank = TRUE_;
    ctol = TRUE_;
    if (*(unsigned char *)comprt == 'N' || *(unsigned char *)comprt == 'n') {
    crank = FALSE_;
    ctol = FALSE_;
    } else {
    if (*(unsigned char *)comprt == 'R' || *(unsigned char *)comprt == 
        'r') {
        ctol = FALSE_;
    }
    if (*(unsigned char *)comprt == 'T' || *(unsigned char *)comprt == 
        't') {
        crank = FALSE_;
    }
    }

    *ierr = 0;
    *iwarn = 0;
    nl = *n + *l;
    k = max(*m,nl);
    p = min(*m,*n);
    if (*m < 1) {
    *ierr = 1;
    }
    if (*n < 1) {
    *ierr = 2;
    }
    if (*l < 1) {
    *ierr = 3;
    }
    if (*ldc < k) {
    *ierr = 4;
    }
    if (*ldx < *n) {
    *ierr = 5;
    }
    if (! crank) {
    if (*rank > p) {
        *ierr = 6;
    }
    }
    *tol1 = max(*tol1,1e-16);
    *tol2 = max(*tol2,1e-16);
    if (*ierr != 0) {
    return 0;
    }

/*     Initialize the solution matrix X. */

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
        x[i__ + j * x_dim1] = 0.;
/* L1: */
    }
/* L2: */
    }

/*     Subroutine DTLS solves a set of linear equations by a Total Least */
/*     Squares Approximation. */

/*     Step 1: */
/*     1.a. : if M .GE. 5*NL/3, then transform [A   ;B   ] into upper */
/*                                               M,N  M,L */
/*            triangular form R by Householder transformations. */

    mc = *m;
    if (*m * 3 >= nl * 5) {
    dqrdc(&c__[c_offset], ldc, m, &nl, &s[1], dumar1, &wrk[1], &c__0);
    i__1 = nl;
    for (j = 1; j <= i__1; ++j) {
        j1 = j + 1;
        i__2 = nl;
        for (i__ = j1; i__ <= i__2; ++i__) {
        c__[i__ + j * c_dim1] = 0.;
/* L3: */
        }
/* L4: */
    }
    mc = nl;
    }
/*                                     T */
/*     1.b. : compute the SVD of  U S V  of [A;B] (or R). */

    dsvdc(&c__[c_offset], ldc, &mc, &nl, &s[1], &wrk[1], dumar2, &c__1, &c__[
        c_offset], ldc, &wrk[nl + 1], &c__1, ierr);
    if (*ierr != 0) {
    *ierr += 1000;
    return 0;
    }

/*     Step 2: Compute the rank of the approximation [A+DA;B+DB]. */

    smax = *tol1;
    if (ctol) {
    smax = sqrt(k * 2.) * smax;
    }
/* Computing 2nd power */
    d__1 = smax;
    smax2 = d__1 * d__1;
    if (crank) {
    *rank = p;
/*        WHILE (RANK .GT. 0) .AND. (S(RANK) .LE. SMAX) DO */
L5:
    if (*rank > 0) {
        if (s[*rank] <= smax) {
        --(*rank);
        goto L5;
        }
    }
/*        END WHILE */
    }

/*     Step 3: Compute the Householder matrix Q and matrices F and Y */
/*     such that F is nonsingular. */

/*     REPEAT */

/*        Adjust the rank if S(RANK) has multiplicity > 1. */

L6:
    r1 = *rank + 1;
/*        WHILE (RANK .GT. 0) .AND. (S(RANK)**2 - S(R1)**2 .LE. SMAX2) DO */
L7:
    if (*rank > 0) {
/* Computing 2nd power */
    d__1 = s[*rank];
/* Computing 2nd power */
    d__2 = s[r1];
    if (d__1 * d__1 - d__2 * d__2 <= smax2) {
        --(*rank);
        *iwarn = 1;
        goto L7;
    }
    }
/*        END WHILE */
    if (*rank == 0) {
    return 0;
    }
    r1 = *rank + 1;

/*        Compute the Householder matrix Q and matrices F and Y. */

    nl1 = max(*n,r1) + 1;
    zero = FALSE_;
    i__ = nl;
/*        WHILE ((.NOT.ZERO) .AND. (I .GE. NL1)) DO */
L8:
    if (! zero && i__ >= nl1) {
    k = i__ - *rank;
    dcopy(&k, &c__[i__ + r1 * c_dim1], ldc, &wrk[1], &c__1);
    housh(&wrk[1], &k, &k, tol2, &zero, &temp);
    if (! zero) {
        tr2(&c__[c_offset], ldc, &wrk[1], &temp, &c__1, &i__, rank, &k);
    }
    --i__;
    goto L8;
    }
/*        END WHILE */
    n1 = *n + 1;
    if (zero || (d__1 = c__[n1 + n1 * c_dim1], abs(d__1)) <= *tol2) {
    --(*rank);
    *iwarn = 2;
    goto L6;
    }
/*     UNTIL ((.NOT.ZERO) .AND. (ABS(C(N1,N1) .GT. TOL2)) */

/*     Step 4: Solve X F = -Y by forward elimination, */
/*             (F is upper triangular). */

    d__1 = -1. / c__[n1 + n1 * c_dim1];
    daxpy(n, &d__1, &c__[n1 * c_dim1 + 1], &c__1, &x[x_offset], &c__1);
    i__1 = *l;
    for (j = 2; j <= i__1; ++j) {
    nj = *n + j;
    temp = c__[nj + nj * c_dim1];
    j1 = j - 1;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
        x[i__ + j * x_dim1] = -(c__[i__ + nj * c_dim1] + ddot(&j1, &c__[
            n1 + nj * c_dim1], &c__1, &x[i__ + x_dim1], ldx)) / temp;
/* L10: */
    }
/* L11: */
    }
    return 0;
/* *** Last line of DTLS *********************************************** */
} /* dtls_ */

