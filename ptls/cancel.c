/* cancel.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int cancel(doublereal *u, integer *ldu, doublereal *v, 
    integer *ldv, doublereal *q, doublereal *e, integer *m, integer *n, 
    integer *i__, integer *k, doublereal *tol, logical *wantu, logical *
    wantv)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__, f, g, h__;
    static integer l;
    static doublereal s;
    static integer i1, l1;
    extern /* Subroutine */ int drot(integer *, doublereal *, integer *, 
        doublereal *, integer *, doublereal *, doublereal *);


/*     PURPOSE: */

/*     Either, subroutine CANCEL separates a zero singular value of a */
/*     subbidiagonal matrix of order k, k <= p, of the bidiagonal */

/*               !Q(1) E(2)  0    ...   0  ! */
/*               ! 0   Q(2) E(3)        .  ! */
/*           J = ! .                    .  ! */
/*               ! .                   E(p)! */
/*               ! 0   ...             Q(p)! */

/*     with p = min(M,N), by annihilating one or two superdiagonal */
/*     elements E(i) and/or E(i+1). */
/*     Or, CANCEL annihilates the element E(M+1) of the bidiagonal matrix */

/*               !Q(1) E(2)  0    ...   0     0   ! */
/*               ! 0   Q(2) E(3)        .     .   ! */
/*           J = ! .                    .     .   ! */
/*               ! .                   E(M)   .   ! */
/*               ! 0   ...             Q(M) E(M+1)! */

/*     ARGUMENT LIST: */

/*     U - DOUBLE PRECISION array of DIMENSION (LDU,p). */
/*         On entry, U contains the M by p (p = min(M,N)) left trans- */
/*         formation matrix. */
/*         On return, the Givens rotations S on the left, annihilating */
/*         E(i+1), have been postmultiplied into U. */
/*         NOTE: U is not referenced if WANTU = .FALSE. . */
/*     LDU - INTEGER. */
/*         LDU is the leading dimension of the array U (LDU >= M). */
/*     V - DOUBLE PRECISION array of DIMENSION (LDV,s). */
/*         On entry, V contains the N by s (s = min(M+1,N)) right trans- */
/*         formation matrix. */
/*         On return, the Givens rotations T on the right, annihilating */
/*         E(i), have been postmultiplied into V. */
/*         NOTE: V is not referenced if WANTV = .FALSE. . */
/*     LDV - INTEGER. */
/*         LDV is the leading dimension of the array V (LDV >= N). */
/*     Q - DOUBLE PRECISION array of DIMENSION (p). */
/*         On entry, Q contains the diagonal entries of the bidiagonal J. */
/*         p = min(M,N). */
/*         On return, Q contains the diagonal elements of the transformed */
/*         bidiagonal S' J T. */
/*     E - DOUBLE PRECISION array of DIMENSION (s). */
/*         On entry, E(i), i=2,...,s, contain the superdiagonal entries */
/*         of the bidiagonal J. s = min(M+1,N), E(1) = 0.0D0. */
/*         On return, E contains the superdiagonal elements of the trans- */
/*         formed bidiagonal S' J T. */
/*     M - INTEGER. */
/*         M is the number of rows of the matrix U. */
/*     N - INTEGER. */
/*         N is the number of rows of the matrix V. */
/*     I - INTEGER. */
/*         Either, I is the index of the negligible diagonal entry Q(I) */
/*         of the bidiagonal J, i,e. abs(Q(I)) <= TOL, I <= p. */
/*         Or, I = M + 1 if E(M+1) is to be annihilated. */
/*     K - INTEGER. */
/*         Either, K is the index of the last diagonal entry of the con- */
/*         sidered subbidiagonal of J, i.e. abs(E(K+1)) <= TOL, K <= p. */
/*         Or, K = M + 1 if E(M+1) is to be annihilated. */
/*     TOL - DOUBLE PRECISION. */
/*         Specifies that matrix elements Q(i), which are <= TOL in */
/*         absolute value, are considered to be zero. */
/*     WANTU - LOGICAL. */
/*         Logical indicating the need for postmultiplying the Givens */
/*         rotations S on the left into U. */
/*     WANTV - LOGICAL. */
/*         Logical indicating the need for postmultiplying the Givens */
/*         rotations T on the right into V. */

/*     EXTERNAL SUBROUTINES and FUNCTIONS: */

/*     DROT from BLAS. */

/*     METHOD DESCRIPTION: */

/*     Let the considered subbidiagonal be */

/*               !Q(1) E(2)  0                    ...   0  ! */
/*               ! 0   Q(2) E(3)                  ...      ! */
/*               ! .                              ...      ! */
/*               !             Q(i-1) E(i)              .  ! */
/*          Jk = !                    Q(i) E(i+1)       .  ! */
/*               !                         Q(i+1) .        ! */
/*               ! .                              ..       ! */
/*               ! .                                   E(k)! */
/*               ! 0    ...                       ...  Q(k)! */

/*     A zero singular value of Jk manifests itself by a zero diagonal */
/*     entry Q(i) or in practice, a negligible value of Q(i). */
/*     We call Q(i) negligible if abs(Q(i)) <= TOL. */
/*     When such a negligible diagonal element Q(i) in Jk is present, */
/*     the subbidiagonal Jk is splitted by the routine CANCEL into 2 or */
/*     3 unreduced subbidiagonals by annihilating E(i+1) (if i<k) using */
/*     Givens rotations S on the left and by annihilating E(i) (if i>1) */
/*     using Givens rotations T on the right until Jk is reduced to the */
/*     form : */

/*               !Q(1) E(2)  0                ...   0  ! */
/*               ! 0         .                ...      ! */
/*               ! .                          ...      ! */
/*               !         Q(i-1) 0                 .  ! */
/*     S' Jk T = !                0   0             .  ! */
/*               !                   Q(i+1)   .        ! */
/*               ! .                          ..       ! */
/*               ! .                               E(k)! */
/*               ! 0    ...                   ...  Q(k)! */

/*     For more details, see [1, pp.11.12-11.14]. */
/*     The case of the annihilation of E(M+1) can be treated by the same */
/*     process. This may be seen by augmenting the matrix J with an extra */
/*     row of zeros, i.e. by introducing Q(M+1) = 0. */

/*     REFERENCES: */

/*     [1] J.J. Dongarra, J.R. Bunch, C.B. Moler and G.W. Stewart, */
/*         LINPACK User's Guide. SIAM, Philadelphia (1979). */

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

    /* Function Body */
    if (*i__ <= *m) {
    q[*i__] = 0.;
    }

/*     Annihilate E(I+1) (if I < K). */

    if (*i__ < *k) {
    c__ = 0.;
    s = 1.;
    i1 = *i__ + 1;
    i__1 = *k;
    for (l = i1; l <= i__1; ++l) {
        g = e[l];
        f = s * g;
        e[l] = c__ * g;
        if (abs(f) <= *tol) {
        goto L2;
        }
        g = q[l];
/* Computing 2nd power */
        d__1 = f;
/* Computing 2nd power */
        d__2 = g;
        h__ = sqrt(d__1 * d__1 + d__2 * d__2);
        q[l] = h__;
        c__ = g / h__;
        s = -f / h__;
        if (*wantu) {
        drot(m, &u[*i__ * u_dim1 + 1], &c__1, &u[l * u_dim1 + 1], &
            c__1, &c__, &s);
        }
/* L1: */
    }
    }

/*     Annihilate E(I) (if I > 1). */

L2:
    if (*i__ > 1) {
    i1 = *i__ - 1;
    f = e[*i__];
    e[*i__] = 0.;
    i__1 = i1;
    for (l1 = 1; l1 <= i__1; ++l1) {
        if (abs(f) <= *tol) {
        return 0;
        }
        l = *i__ - l1;
        g = q[l];
        if (abs(g) <= *tol) {
        g = 0.;
        h__ = abs(f);
        } else {
/* Computing 2nd power */
        d__1 = f;
/* Computing 2nd power */
        d__2 = g;
        h__ = sqrt(d__1 * d__1 + d__2 * d__2);
        }
        q[l] = h__;
        c__ = g / h__;
        s = -f / h__;
        g = e[l];
        f = s * g;
        e[l] = c__ * g;
        if (*wantv) {
        drot(n, &v[*i__ * v_dim1 + 1], &c__1, &v[l * v_dim1 + 1], &
            c__1, &c__, &s);
        }
/* L3: */
    }
    e[1] = 0.;
    }
    return 0;
/* *** Last line of CANCEL ********************************************* */
} /* cancel_ */

