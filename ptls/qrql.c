/* qrql.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int qrql(doublereal *u, integer *ldu, doublereal *v, 
    integer *ldv, integer *m, integer *n, integer *rank, doublereal *
    theta, doublereal *q, doublereal *e, logical *inul, doublereal *tol1, 
    doublereal *tol2, logical *wantu, logical *wantv, integer *ierr, 
    integer *iwarn)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, p, r__;
    static doublereal x;
    static integer i1;
    static logical noc12;
    static integer iter;
    static doublereal shift;
    extern /* Subroutine */ int estim(doublereal *, doublereal *, integer *, 
        integer *, doublereal *, doublereal *, doublereal *, integer *);
    static integer maxit;
    extern /* Subroutine */ int cancel(doublereal *, integer *, doublereal *,
         integer *, doublereal *, doublereal *, integer *, integer *, 
        integer *, integer *, doublereal *, logical *, logical *);
    static integer numeig;
    extern integer nsingv(doublereal *, doublereal *, integer *, doublereal *
        , doublereal *, doublereal *);
    extern /* Subroutine */ int qlstep(doublereal *, integer *, doublereal *,
         integer *, doublereal *, doublereal *, integer *, integer *, 
        integer *, integer *, doublereal *, logical *, logical *), 
        qrstep(doublereal *, integer *, doublereal *, integer *, 
        doublereal *, doublereal *, integer *, integer *, integer *, 
        integer *, doublereal *, logical *, logical *);


/*     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven). */

/*     REVISIONS: 1988, February 15. */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. External Subroutines/Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --q;
    --e;
    --inul;

    /* Function Body */
    *ierr = 0;
    p = min(*m,*n);

/*     Estimate THETA (if not fixed by the user), and set R. */

    if (*rank >= 0) {
    j = p - *rank;
    estim(&q[1], &e[1], &p, &j, theta, tol1, tol2, iwarn);
    if (j <= 0) {
        return 0;
    }
    r__ = p - j;
    } else {
    r__ = 0;
    }

    maxit = 50;
    *rank = p;
    i__1 = p;
    for (i__ = 1; i__ <= i__1; ++i__) {
    if (inul[i__]) {
        --(*rank);
    }
/* L1: */
    }
    e[1] = 0.;

/*     From now K is the smallest known index such that the subbidia- */
/*     gonals with indices > K belong to C1 or C2. */
/*     RANK = P - SUM(dimensions of known elements of C2). */

    k = p;
/*     WHILE (C3 NOT EMPTY) DO */
L2:
    if (*rank > r__ && k > 0) {
/*        WHILE (K .GT. 0 .AND INUL(K)) DO */

/*        Search for the rightmost index of a subbidiagonal of C1 or C3. */

L3:
    if (k > 0) {
        if (inul[k]) {
        --k;
        goto L3;
        }
    }
/*        END WHILE 3 */

    if (k == 0) {
        return 0;
    }

    iter = 0;
    noc12 = TRUE_;
/*        WHILE ((ITER < MAXIT) .AND. (No element of C1 or C2 found)) DO */
L4:
    if (iter < maxit && noc12) {

/*           Search for negligible Q(I) or E(I). */

        i__ = k;
        x = (d__1 = q[i__], abs(d__1));
        shift = x;
/*           WHILE (ABS(Q(I)) > TOL2 .AND. ABS(E(I)) > TOL2) DO */
L5:
        if (x > *tol2 && (d__1 = e[i__], abs(d__1)) > *tol2) {
        --i__;
        x = (d__1 = q[i__], abs(d__1));
        if (x < shift) {
            shift = x;
        }
        goto L5;
        }
/*           END WHILE 5 */

/*           Classify the subbidiagonal found. */

        j = k - i__ + 1;
        if (x <= *tol2 || k == i__) {
        noc12 = FALSE_;
        } else {
        numeig = nsingv(&q[i__], &e[i__], &j, theta, tol1, tol2);
        if (numeig >= j || numeig <= 0) {
            noc12 = FALSE_;
        }
        }
        if (noc12) {
        if (shift >= *theta) {
            shift = 0.;
        }
        if ((d__1 = q[k], abs(d__1)) <= (d__2 = q[i__], abs(d__2))) {
            qrstep(&u[u_offset], ldu, &v[v_offset], ldv, &q[1], &e[1]
                , m, n, &i__, &k, &shift, wantu, wantv);
        } else {
            qlstep(&u[u_offset], ldu, &v[v_offset], ldv, &q[1], &e[1]
                , m, n, &i__, &k, &shift, wantu, wantv);
        }
        ++iter;
        }
        goto L4;
    }
/*        END WHILE 4 */

    if (iter == maxit) {
        *ierr = 10;
        return 0;
    }

    if (x <= *tol2) {

/*           Split at negligible diagonal element abs(Q(I)) <= TOL2. */

        cancel(&u[u_offset], ldu, &v[v_offset], ldv, &q[1], &e[1], m, n, 
            &i__, &k, tol2, wantu, wantv);
        inul[i__] = TRUE_;
        --(*rank);
    } else {

/*           A negligible superdiagonal element abs(E(I)) <= TOL2 has */
/*           been found, the corresponding subbidiagonal belongs to */
/*           C1 or C2. Treat this subbidiagonal. */

        if (j >= 2) {
        if (numeig == j) {
            i__1 = k;
            for (i1 = i__; i1 <= i__1; ++i1) {
            inul[i1] = TRUE_;
/* L6: */
            }
            *rank -= j;
            k -= j;
        } else {
            k = i__ - 1;
        }
        } else {
        if (x <= *theta + *tol1) {
            inul[i__] = TRUE_;
            --(*rank);
        }
        --k;
        }
    }
    goto L2;
    }
/*     END WHILE 2 */
    return 0;
/* *** Last line of QRQL *********************************************** */
} /* qrql_ */

