/* qrstep.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int qrstep(doublereal *u, integer *ldu, doublereal *v, 
    integer *ldv, doublereal *q, doublereal *e, integer *m, integer *n, 
    integer *i__, integer *k, doublereal *shift, logical *wantu, logical *
    wantv)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__, f, g, h__;
    static integer j;
    static doublereal s, x, y, z__;
    static integer i1, ll;
    extern /* Subroutine */ int drot(integer *, doublereal *, integer *, 
        doublereal *, integer *, doublereal *, doublereal *);


/*     PURPOSE: */

/*     The subroutine QRSTEP performs one QR iteration step onto the */
/*     unreduced subbidiagonal Jk: */

/*              !Q(i) E(i+1)  0  ...    0  ! */
/*              ! 0   Q(i+1) E(i+2)     .  ! */
/*         Jk = ! .                     .  ! */
/*              ! .                        ! */
/*              ! .                    E(k)! */
/*              ! 0   ...              Q(k)! */

/*     with k <= p and i >= 1, p = min(M,N), of the bidiagonal J: */

/*              !Q(1) E(2)  0    ...   0  ! */
/*              ! 0   Q(2) E(3)        .  ! */
/*          J = ! .                    .  ! */
/*              ! .                   E(p)! */
/*              ! 0   ...             Q(p)! */

/*     Hereby, Jk is transformed to  S'Jk T with S and T products of */
/*     Givens rotations. These Givens rotations S (resp.,T) will be post- */
/*     multiplied into U (resp.,V), if WANTU (resp.,WANTV) = .TRUE. */

/*     ARGUMENT LIST: */

/*     U - DOUBLE PRECISION array of DIMENSION (LDU,min(M,N)). */
/*         On entry, U may contain the M by p (p=min(M,N)) left transfor- */
/*         mation matrix. */
/*         On return, if WANTU = .TRUE., the Givens rotations S on the */
/*         left have been postmultiplied into U. */
/*         NOTE: U is not referenced if WANTU = .FALSE. */
/*     LDU - INTEGER. */
/*         LDU is the leading dimension of the array U (LDU >= M). */
/*     V - DOUBLE PRECISION array of DIMENSION (LDV,min(M,N)). */
/*         On entry, V may contain the N by p (p=min(M,N)) right trans- */
/*         formation matrix. */
/*         On return, if WANTV = .TRUE., the Givens rotations T on the */
/*         right have been postmultiplied into V. */
/*         NOTE: V is not referenced if WANTV is .false. */
/*     LDV - INTEGER. */
/*         LDV is the leading dimension of the array V (LDV >= N). */
/*     Q - DOUBLE PRECISION array of DIMENSION (min(M,N)). */
/*         On entry, Q contains the diagonal entries of the bidiagonal J. */
/*         On return, Q contains the diagonal entries of the transformed */
/*         matrix S' J T. */
/*     E - DOUBLE PRECISION array of DIMENSION (min(M,N)). */
/*         On entry, E contains the superdiagonal entries of J. */
/*         On return, E contains the superdiagonal entries of the trans- */
/*         formed matrix S' J T. E(i) = 0. */
/*     M - INTEGER. */
/*         M is the number of rows of the matrix U. */
/*     N - INTEGER. */
/*         N is the number of rows of the matrix V. */
/*     I - INTEGER. */
/*         I is the index of the first diagonal entry of the considered */
/*         unreduced subbidiagonal Jk of J. */
/*     K - INTEGER. */
/*         K is the index of the last diagonal entry of the considered */
/*         unreduced subbidiagonal Jk of J. */
/*     SHIFT - DOUBLE PRECISION. */
/*         Value of the shift used in the QR iteration step. */
/*     WANTU - LOGICAL. */
/*         WANTU = .TRUE. if the Givens rotations S must be postmulti- */
/*         plied on the left into U, else .FALSE. */
/*     WANTV - LOGICAL. */
/*         WANTV = .TRUE. if the Givens rotations T must be postmulti- */
/*         plied  on the left into V, else .FALSE. */

/*     EXTERNAL SUBROUTINES AND FUNCTIONS: */

/*     DROT from BLAS. */

/*     METHOD DESCRIPTION: */

/*     QR iterations diagonalize the bidiagonal by zeroing the super- */
/*     diagonal elements of Jk from bottom to top. */
/*     The routine QRSTEP overwrites Jk with the bidiagonal matrix */
/*     S' Jk T where S and T are Givens rotations. */
/*     T is essentially the orthogonal matrix that would be obtained by */
/*     applying one implicit symmetric shift QR step onto the matrix */
/*     Jk'Jk. This step factors the matrix (Jk'Jk - shift*I) into a */
/*     product of an orthogonal matrix T and an upper triangular matrix. */
/*     See [1,Sec.8.2-8.3] for more details. */

/*     REFERENCES: */
/*     [1] G.H. Golub and C.F. Van Loan, Matrix Computations. The Johns */
/*         Hopkins University Press, Baltimore,Maryland (1983). */

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
    x = q[*i__];
    g = *shift;
    f = (x - g) * (x + g) / x;
    c__ = 1.;
    s = 1.;
    i1 = *i__ + 1;
    i__1 = *k;
    for (j = i1; j <= i__1; ++j) {
    ll = j - 1;
    g = e[j];
    y = q[j];
    h__ = s * g;
    g = c__ * g;
/* Computing 2nd power */
    d__1 = f;
/* Computing 2nd power */
    d__2 = h__;
    z__ = sqrt(d__1 * d__1 + d__2 * d__2);
    e[ll] = z__;
    c__ = f / z__;
    s = h__ / z__;
    f = x * c__ + g * s;
    g = -x * s + g * c__;
    h__ = y * s;
    y *= c__;
    if (*wantv) {
        drot(n, &v[ll * v_dim1 + 1], &c__1, &v[j * v_dim1 + 1], &c__1, &
            c__, &s);
    }
/* Computing 2nd power */
    d__1 = f;
/* Computing 2nd power */
    d__2 = h__;
    z__ = sqrt(d__1 * d__1 + d__2 * d__2);
    q[ll] = z__;
    c__ = f / z__;
    s = h__ / z__;
    f = c__ * g + s * y;
    x = -s * g + c__ * y;
    if (*wantu) {
        drot(m, &u[ll * u_dim1 + 1], &c__1, &u[j * u_dim1 + 1], &c__1, &
            c__, &s);
    }
/* L10: */
    }
    e[*k] = f;
    q[*k] = x;
    e[*i__] = 0.;
    return 0;
/* *** Last line of QRSTEP ********************************************* */
} /* qrstep_ */

