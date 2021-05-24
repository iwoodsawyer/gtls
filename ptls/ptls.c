/* ptls.f -- translated by f2c (version 20100827).
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
static logical c_false = FALSE_;
static logical c_true = TRUE_;

/* Subroutine */ int ptls(doublereal *c__, integer *ldc, integer *m, integer 
    *n, integer *l, integer *rank, doublereal *theta, doublereal *x, 
    integer *ldx, doublereal *q, logical *inul, doublereal *wrk, integer *
    iwrk, logical *lwrk, doublereal *tol1, doublereal *tol2, integer *
    ierr, integer *iwarn)
{
    /* System generated locals */
    integer c_dim1, c_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, p, i1, j1, n1, mc;
    static doublereal hh;
    static integer ii, nj, nl, pp1, mnl1;
    extern doublereal ddot(integer *, doublereal *, integer *, doublereal *, 
        integer *);
    extern /* Subroutine */ int init(doublereal *, integer *, integer *, 
        integer *);
    static doublereal temp;
    extern /* Subroutine */ int qrql(doublereal *, integer *, doublereal *, 
        integer *, integer *, integer *, integer *, doublereal *, 
        doublereal *, doublereal *, logical *, doublereal *, doublereal *,
         logical *, logical *, integer *, integer *);
    static logical zero;
    extern /* Subroutine */ int dqrdc(doublereal *, integer *, integer *, 
        integer *, doublereal *, integer *, doublereal *, integer *), 
        dcopy(integer *, doublereal *, integer *, doublereal *, integer *
        ), daxpy(integer *, doublereal *, doublereal *, integer *, 
        doublereal *, integer *), housh(doublereal *, integer *, integer 
        *, doublereal *, logical *, doublereal *);
    static integer dumar1[1];
    static doublereal dumar2[1] /* was [1][1] */;
    extern /* Subroutine */ int bidiag(doublereal *, integer *, integer *, 
        integer *, doublereal *, doublereal *, doublereal *), cancel(
        doublereal *, integer *, doublereal *, integer *, doublereal *, 
        doublereal *, integer *, integer *, integer *, integer *, 
        doublereal *, logical *, logical *);
    static doublereal inprod;


/*     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. External Subroutines/Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --lwrk;
    --iwrk;
    --inul;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --q;
    --wrk;

    /* Function Body */
    *ierr = 0;
    *iwarn = 0;
    if (*m < 1) {
    *ierr = 1;
    }
    if (*n < 1) {
    *ierr = 2;
    }
    if (*l < 1) {
    *ierr = 3;
    }
/* Computing MAX */
    i__1 = *m, i__2 = *n + *l;
    if (*ldc < max(i__1,i__2)) {
    *ierr = 4;
    }
    if (*ldx < *n) {
    *ierr = 5;
    }
    if (*rank > min(*m,*n)) {
    *ierr = 6;
    }
    if (*rank < 0 && *theta < 0.) {
    *ierr = 7;
    }
    if (*tol1 < 0.) {
    *ierr = 8;
    }
    if (*tol2 < 0.) {
    *ierr = 9;
    }
    if (*ierr != 0) {
    return 0;
    }

/*     Initializations. */

    nl = *n + *l;
    mnl1 = *m + nl + 1;
    p = min(*m,nl);
    i__1 = p;
    for (i__ = 1; i__ <= i__1; ++i__) {
    inul[i__] = FALSE_;
    lwrk[i__] = FALSE_;
/* L1: */
    }
    pp1 = p + 1;
    i__1 = nl;
    for (i__ = pp1; i__ <= i__1; ++i__) {
    inul[i__] = TRUE_;
    lwrk[i__] = FALSE_;
/* L2: */
    }

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
        x[i__ + j * x_dim1] = 0.;
/* L3: */
    }
    }

/*     Subroutine PTLS solves a set of linear equations by a Total Least */
/*     Squares Approximation, based on the Partial SVD. */

/*     Step 1: Bidiagonalization phase */
/*             ----------------------- */
/*     1.a): If M >= 5*(N+L)/3, transform C into upper triangular form R */
/*           by Householder transformations. */

    mc = *m;
    if (*m * 3 >= nl * 5) {
    dqrdc(&c__[c_offset], ldc, m, &nl, &q[1], dumar1, &wrk[1], &c__0);
    i__2 = nl;
    for (j = 1; j <= i__2; ++j) {
        j1 = j + 1;
        i__1 = nl;
        for (i__ = j1; i__ <= i__1; ++i__) {
        c__[i__ + j * c_dim1] = 0.;
/* L4: */
        }
    }
    mc = nl;
    }

/*     1.b): Transform C (or R) into bidiagonal form Q using Householder */
/*           transformations. */

    bidiag(&c__[c_offset], ldc, &mc, &nl, &q[1], &q[p + 1], &wrk[1]);

/*     Store the Householder transformations performed onto the rows of C */
/*     in the last storage locations of the work array WRK. */

/* Computing MIN */
    i__1 = nl - 2;
    mc = min(i__1,*m);
    if (mc > 0) {
    k = mnl1;
    i__1 = mc;
    for (ii = 1; ii <= i__1; ++ii) {
        j = mc - ii + 1;
        nj = nl - j;
        dcopy(&nj, &c__[j + (j + 1) * c_dim1], ldc, &wrk[k], &c__1);
        k += nj;
/* L5: */
    }
    }

/*     1.c): Initialize the right singular base matrix V with the identi- */
/*           ty matrix (V overwrites C). */

    init(&c__[c_offset], ldc, &nl, &nl);

/*     1.d): If M < N+L, bring the bidiagonal Q to M by M by cancelling */
/*           its last superdiagonal element using Givens rotations. */

    pp1 = p;
    if (*m < nl) {
    pp1 = p + 1;
    cancel(dumar2, &c__1, &c__[c_offset], ldc, &q[1], &q[pp1], &p, &pp1, 
        &pp1, &pp1, tol2, &c_false, &c_true);
    }

/*     REPEAT */

/*        Compute the Householder matrix Q and matrices F and Y such that */
/*        F is nonsingular. */
/*        Step 2: Partial diagonalization phase. */
/*                ----------------------------- */
/*        Diagonalize the bidiagonal Q partially until convergence to */
/*        the desired right singular subspace. */

L6:
    qrql(dumar2, &c__1, &c__[c_offset], ldc, &p, &pp1, rank, theta, &q[1], &
        q[p + 1], &inul[1], tol1, tol2, &c_false, &c_true, ierr, iwarn);

    if (*ierr != 0) {
    return 0;
    }

/*        Step 3: Back transformation phase. */
/*                ------------------------- */
/*        Apply the Householder transformations (stored in WRK) perfor- */
/*        med onto the rows of C during the bidiagonalization phase, to */
/*        the selected base vectors in the right singular base matrix */
/*        of C. */

    i__1 = nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
    if (inul[i__] && ! lwrk[i__]) {
        k = mnl1;
        i__2 = mc;
        for (ii = 1; ii <= i__2; ++ii) {
        j = mc - ii + 1;
        nj = nl - j;
        j1 = j + 1;
        if ((d__1 = wrk[k], abs(d__1)) > *tol2) {
            temp = -ddot(&nj, &wrk[k], &c__1, &c__[j1 + i__ * c_dim1]
                , &c__1) / wrk[k];
            daxpy(&nj, &temp, &wrk[k], &c__1, &c__[j1 + i__ * c_dim1]
                , &c__1);
            k += nj;
        }
/* L7: */
        }
        lwrk[i__] = TRUE_;
    }
/* L8: */
    }
    if (*rank <= 0) {
    return 0;
    }

/*        Step 4: Compute matrices F and Y using Householder transf. Q. */
/*                ------------------------ */
    k = 0;
    i__1 = nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
    if (inul[i__]) {
        ++k;
        iwrk[k] = i__;
    }
/* L9: */
    }

    if (k < *l) {

/*           Rank TLS approximation is larger than min(M,N). */

    *ierr = 3;
    return 0;
    }
    n1 = *n + 1;
    zero = FALSE_;
    i__ = nl;

/*        WHILE ((K > 1) .AND. (I > N) .AND. (.NOT.ZERO)) DO */
L10:
    if (k > 1 && i__ > *n && ! zero) {
    i__1 = k;
    for (j = 1; j <= i__1; ++j) {
        wrk[j] = c__[i__ + iwrk[j] * c_dim1];
/* L11: */
    }

/*           Compute Householder transformation. */

    housh(&wrk[1], &k, &k, tol2, &zero, &temp);
    if (! zero) {

/*              Apply Householder transformation onto the selected base */
/*              vectors. */

        i__1 = i__;
        for (i1 = 1; i1 <= i__1; ++i1) {
        inprod = 0.;
        i__2 = k;
        for (j = 1; j <= i__2; ++j) {
            inprod += wrk[j] * c__[i1 + iwrk[j] * c_dim1];
/* L12: */
        }
        hh = inprod * temp;
        i__2 = k;
        for (j = 1; j <= i__2; ++j) {
            j1 = iwrk[j];
            c__[i1 + j1 * c_dim1] -= wrk[j] * hh;
/* L13: */
        }
        }

        --k;
    }
    --i__;
    goto L10;
    }
/*        END WHILE 10 */

    if (! zero) {
    k = n1 - *rank;
    }

/*        If F singular, lower the rank of the TLS approximation . */

    if ((d__1 = c__[n1 + iwrk[k] * c_dim1], abs(d__1)) <= *tol2) {
    --(*rank);
    *iwarn = 2;
    *theta = -1.;
    goto L6;
    }
/*     UNTIL ((.NOT.ZERO) .AND. (F nonsingular)) */

/*     Step 5: Compute TLS solution */
/*             -------------------- */
/*     Solve X F = -Y  by forward elimination  (F is upper triangular). */
    nj = iwrk[k];
    d__1 = -1. / c__[n1 + nj * c_dim1];
    daxpy(n, &d__1, &c__[nj * c_dim1 + 1], &c__1, &x[x_offset], &c__1);
    i__2 = *l;
    for (j = 2; j <= i__2; ++j) {
    j1 = j - 1;
    nj = iwrk[k + j1];
    temp = c__[*n + j + nj * c_dim1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        x[i__ + j * x_dim1] = -(c__[i__ + nj * c_dim1] + ddot(&j1, &c__[
            n1 + nj * c_dim1], &c__1, &x[i__ + x_dim1], ldx)) / temp;
/* L14: */
    }
/* L15: */
    }
    return 0;
/* *** Last line of PTLS *********************************************** */
} /* ptls_ */

