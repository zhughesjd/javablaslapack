package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackCL extends Library
{

	public static LapackCL instance = (LapackCL) Native.loadLibrary("liblapack",LapackCL.class);

/**
*> \brief \b CLABRD reduces the first nb rows and columns of a general matrix to a bidiagonal form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLABRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clabrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clabrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clabrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y,<br>
*                          LDY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, LDX, LDY, M, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * )<br>
*       COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ), X( LDX, * ),<br>
*      $                   Y( LDY, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLABRD reduces the first NB rows and columns of a complex general<br>
*> m by n matrix A to upper or lower real bidiagonal form by a unitary<br>
*> transformation Q**H * A * P, and returns the matrices X and Y which<br>
*> are needed to apply the transformation to the unreduced part of A.<br>
*><br>
*> If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower<br>
*> bidiagonal form.<br>
*><br>
*> This is an auxiliary routine called by CGEBRD<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows in the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns in the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The number of leading rows and columns of A to be reduced.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the m by n general matrix to be reduced.<br>
*>          On exit, the first NB rows and columns of the matrix are<br>
*>          overwritten; the rest of the array is unchanged.<br>
*>          If m >= n, elements on and below the diagonal in the first NB<br>
*>            columns, with the array TAUQ, represent the unitary<br>
*>            matrix Q as a product of elementary reflectors; and<br>
*>            elements above the diagonal in the first NB rows, with the<br>
*>            array TAUP, represent the unitary matrix P as a product<br>
*>            of elementary reflectors.<br>
*>          If m < n, elements below the diagonal in the first NB<br>
*>            columns, with the array TAUQ, represent the unitary<br>
*>            matrix Q as a product of elementary reflectors, and<br>
*>            elements on and above the diagonal in the first NB rows,<br>
*>            with the array TAUP, represent the unitary matrix P as<br>
*>            a product of elementary reflectors.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (NB)<br>
*>          The diagonal elements of the first NB rows and columns of<br>
*>          the reduced matrix.  D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (NB)<br>
*>          The off-diagonal elements of the first NB rows and columns of<br>
*>          the reduced matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ<br>
*> \verbatim<br>
*>          TAUQ is COMPLEX array dimension (NB)<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the unitary matrix Q. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP<br>
*> \verbatim<br>
*>          TAUP is COMPLEX array, dimension (NB)<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the unitary matrix P. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (LDX,NB)<br>
*>          The m-by-nb matrix X required to update the unreduced part<br>
*>          of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X. LDX >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (LDY,NB)<br>
*>          The n-by-nb matrix Y required to update the unreduced part<br>
*>          of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDY<br>
*> \verbatim<br>
*>          LDY is INTEGER<br>
*>          The leading dimension of the array Y. LDY >= max(1,N).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrices Q and P are represented as products of elementary<br>
*>  reflectors:<br>
*><br>
*>     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)<br>
*><br>
*>  Each H(i) and G(i) has the form:<br>
*><br>
*>     H(i) = I - tauq * v * v**H  and G(i) = I - taup * u * u**H<br>
*><br>
*>  where tauq and taup are complex scalars, and v and u are complex<br>
*>  vectors.<br>
*><br>
*>  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in<br>
*>  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in<br>
*>  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).<br>
*><br>
*>  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in<br>
*>  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in<br>
*>  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).<br>
*><br>
*>  The elements of the vectors v and u together form the m-by-nb matrix<br>
*>  V and the nb-by-n matrix U**H which are needed, with X and Y, to apply<br>
*>  the transformation to the unreduced part of the matrix, using a block<br>
*>  update of the form:  A := A - V*Y**H - X*U**H.<br>
*><br>
*>  The contents of A on exit are illustrated by the following examples<br>
*>  with nb = 2:<br>
*><br>
*>  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):<br>
*><br>
*>    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )<br>
*>    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )<br>
*>    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )<br>
*>    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )<br>
*>    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )<br>
*>    (  v1  v2  a   a   a  )<br>
*><br>
*>  where a denotes an element of the original matrix which is unchanged,<br>
*>  vi denotes an element of the vector defining H(i), and ui an element<br>
*>  of the vector defining G(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clabrd_(INTEGER M,INTEGER N,INTEGER NB,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAUQ,float[] TAUP,float[] X,INTEGER LDX,float[] Y,INTEGER LDY);
/**
*> \brief \b CLACGV conjugates a complex vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLACGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacgv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacgv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacgv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLACGV( N, X, INCX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLACGV conjugates a complex vector of length N.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The length of the vector X.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension<br>
*>                         (1+(N-1)*abs(INCX))<br>
*>          On entry, the vector of length N to be conjugated.<br>
*>          On exit, X is overwritten with conjg(X).<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The spacing between successive elements of X.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clacgv_(INTEGER N,float[] X,INTEGER INCX);
/**
*> \brief \b CLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLACN2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacn2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacn2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacn2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLACN2( N, V, X, EST, KASE, ISAVE )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            KASE, N<br>
*       REAL               EST<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISAVE( 3 )<br>
*       COMPLEX            V( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLACN2 estimates the 1-norm of a square, complex matrix A.<br>
*> Reverse communication is used for evaluating matrix-vector products.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The order of the matrix.  N >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension (N)<br>
*>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)<br>
*>         (W is not returned).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>         On an intermediate return, X should be overwritten by<br>
*>               A * X,   if KASE=1,<br>
*>               A**H * X,  if KASE=2,<br>
*>         where A**H is the conjugate transpose of A, and CLACN2 must be<br>
*>         re-called with all the other parameters unchanged.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EST<br>
*> \verbatim<br>
*>          EST is REAL<br>
*>         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be<br>
*>         unchanged from the previous call to CLACN2.<br>
*>         On exit, EST is an estimate (a lower bound) for norm(A). <br>
*> \endverbatim<br>
*><br>
*> \param[in,out] KASE<br>
*> \verbatim<br>
*>          KASE is INTEGER<br>
*>         On the initial call to CLACN2, KASE should be 0.<br>
*>         On an intermediate return, KASE will be 1 or 2, indicating<br>
*>         whether X should be overwritten by A * X  or A**H * X.<br>
*>         On the final return from CLACN2, KASE will again be 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ISAVE<br>
*> \verbatim<br>
*>          ISAVE is INTEGER array, dimension (3)<br>
*>         ISAVE is used to save variables between calls to SLACN2<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Originally named CONEST, dated March 16, 1988.<br>
*><br>
*>  Last modified:  April, 1999<br>
*><br>
*>  This is a thread safe version of CLACON, which uses the array ISAVE<br>
*>  in place of a SAVE statement, as follows:<br>
*><br>
*>     CLACON     CLACN2<br>
*>      JUMP     ISAVE(1)<br>
*>      J        ISAVE(2)<br>
*>      ITER     ISAVE(3)<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Nick Higham, University of Manchester<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  N.J. Higham, "FORTRAN codes for estimating the one-norm of<br>
*>  a real or complex matrix, with applications to condition estimation",<br>
*>  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.<br>
*><br>
*  =====================================================================<br>
*/
	public void clacn2_(INTEGER N,float[] V,float[] X,REAL EST,INTEGER KASE,int[] ISAVE);
/**
*> \brief \b CLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLACON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLACON( N, V, X, EST, KASE )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            KASE, N<br>
*       REAL               EST<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            V( N ), X( N )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLACON estimates the 1-norm of a square, complex matrix A.<br>
*> Reverse communication is used for evaluating matrix-vector products.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The order of the matrix.  N >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension (N)<br>
*>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)<br>
*>         (W is not returned).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>         On an intermediate return, X should be overwritten by<br>
*>               A * X,   if KASE=1,<br>
*>               A**H * X,  if KASE=2,<br>
*>         where A**H is the conjugate transpose of A, and CLACON must be<br>
*>         re-called with all the other parameters unchanged.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EST<br>
*> \verbatim<br>
*>          EST is REAL<br>
*>         On entry with KASE = 1 or 2 and JUMP = 3, EST should be<br>
*>         unchanged from the previous call to CLACON.<br>
*>         On exit, EST is an estimate (a lower bound) for norm(A). <br>
*> \endverbatim<br>
*><br>
*> \param[in,out] KASE<br>
*> \verbatim<br>
*>          KASE is INTEGER<br>
*>         On the initial call to CLACON, KASE should be 0.<br>
*>         On an intermediate return, KASE will be 1 or 2, indicating<br>
*>         whether X should be overwritten by A * X  or A**H * X.<br>
*>         On the final return from CLACON, KASE will again be 0.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*>  Originally named CONEST, dated March 16, 1988. \n<br>
*>  Last modified:  April, 1999<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Nick Higham, University of Manchester<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  N.J. Higham, "FORTRAN codes for estimating the one-norm of<br>
*>  a real or complex matrix, with applications to condition estimation",<br>
*>  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.<br>
*><br>
*  =====================================================================<br>
*/
	public void clacon_(INTEGER N,float[] V,float[] X,REAL EST,INTEGER KASE);
/**
*> \brief \b CLACP2 copies all or part of a real two-dimensional array to a complex array.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLACP2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacp2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacp2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacp2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLACP2( UPLO, M, N, A, LDA, B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            LDA, LDB, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * )<br>
*       COMPLEX            B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLACP2 copies all or part of a real two-dimensional matrix A to a<br>
*> complex matrix B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies the part of the matrix A to be copied to B.<br>
*>          = 'U':      Upper triangular part<br>
*>          = 'L':      Lower triangular part<br>
*>          Otherwise:  All of the matrix A<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
*>          The m by n matrix A.  If UPLO = 'U', only the upper trapezium<br>
*>          is accessed; if UPLO = 'L', only the lower trapezium is<br>
*>          accessed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,N)<br>
*>          On exit, B = A in the locations specified by UPLO.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,M).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clacp2_(CHARACTER UPLO,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB);
/**
*> \brief \b CLACPY copies all or part of one two-dimensional array to another.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLACPY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacpy.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacpy.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacpy.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLACPY( UPLO, M, N, A, LDA, B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            LDA, LDB, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLACPY copies all or part of a two-dimensional matrix A to another<br>
*> matrix B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies the part of the matrix A to be copied to B.<br>
*>          = 'U':      Upper triangular part<br>
*>          = 'L':      Lower triangular part<br>
*>          Otherwise:  All of the matrix A<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The m by n matrix A.  If UPLO = 'U', only the upper trapezium<br>
*>          is accessed; if UPLO = 'L', only the lower trapezium is<br>
*>          accessed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,N)<br>
*>          On exit, B = A in the locations specified by UPLO.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,M).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clacpy_(CHARACTER UPLO,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB);
/**
*> \brief \b CLACRM multiplies a complex matrix by a square real matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLACRM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacrm.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacrm.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacrm.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, LDB, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               B( LDB, * ), RWORK( * )<br>
*       COMPLEX            A( LDA, * ), C( LDC, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLACRM performs a very simple matrix-matrix multiplication:<br>
*>          C := A * B,<br>
*> where A is M by N and complex; B is N by N and real;<br>
*> C is M by N and complex.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A and of the matrix C.<br>
*>          M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns and rows of the matrix B and<br>
*>          the number of columns of the matrix C.<br>
*>          N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA, N)<br>
*>          A contains the M by N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >=max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB, N)<br>
*>          B contains the N by N matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >=max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is COMPLEX array, dimension (LDC, N)<br>
*>          C contains the M by N matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDC >=max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (2*M*N)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clacrm_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] C,INTEGER LDC,float[] RWORK);
/**
*> \brief \b CLACRT performs a linear transformation of a pair of complex vectors.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLACRT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacrt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacrt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacrt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLACRT( N, CX, INCX, CY, INCY, C, S )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, INCY, N<br>
*       COMPLEX            C, S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            CX( * ), CY( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLACRT performs the operation<br>
*><br>
*>    (  c  s )( x )  ==> ( x )<br>
*>    ( -s  c )( y )      ( y )<br>
*><br>
*> where c and s are complex and the vectors x and y are complex.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of elements in the vectors CX and CY.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CX<br>
*> \verbatim<br>
*>          CX is COMPLEX array, dimension (N)<br>
*>          On input, the vector x.<br>
*>          On output, CX is overwritten with c*x + s*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between successive values of CX.  INCX <> 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CY<br>
*> \verbatim<br>
*>          CY is COMPLEX array, dimension (N)<br>
*>          On input, the vector y.<br>
*>          On output, CY is overwritten with -s*x + c*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>          The increment between successive values of CY.  INCY <> 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is COMPLEX<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is COMPLEX<br>
*>          C and S define the matrix<br>
*>             [  C   S  ].<br>
*>             [ -S   C  ]<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clacrt_(INTEGER N,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY,float[] C,float[] S);
/**
*> \brief \b CLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLADIV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cladiv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cladiv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cladiv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       COMPLEX FUNCTION CLADIV( X, Y )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX            X, Y<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLADIV := X / Y, where X and Y are complex.  The computation of X / Y<br>
*> will not overflow on an intermediary step unless the results<br>
*> overflows.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX<br>
*>          The complex scalars X and Y.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public COMPLEX8 cladiv_(float[] X,float[] Y);
/**
*> \brief \b CLAED0 used by sstedc. Computes all eigenvalues and corresponding eigenvectors of an unreduced symmetric tridiagonal matrix using the divide and conquer method.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAED0 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claed0.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claed0.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claed0.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDQ, LDQS, N, QSIZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               D( * ), E( * ), RWORK( * )<br>
*       COMPLEX            Q( LDQ, * ), QSTORE( LDQS, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Using the divide and conquer method, CLAED0 computes all eigenvalues<br>
*> of a symmetric tridiagonal matrix which is one diagonal block of<br>
*> those from reducing a dense or band Hermitian matrix and<br>
*> corresponding eigenvectors of the dense or band matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] QSIZ<br>
*> \verbatim<br>
*>          QSIZ is INTEGER<br>
*>         The dimension of the unitary matrix used to reduce<br>
*>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The dimension of the symmetric tridiagonal matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry, the diagonal elements of the tridiagonal matrix.<br>
*>         On exit, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>         On entry, the off-diagonal elements of the tridiagonal matrix.<br>
*>         On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX array, dimension (LDQ,N)<br>
*>         On entry, Q must contain an QSIZ x N matrix whose columns<br>
*>         unitarily orthonormal. It is a part of the unitary matrix<br>
*>         that reduces the full dense Hermitian matrix to a<br>
*>         (reducible) symmetric tridiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>         The leading dimension of the array Q.  LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array,<br>
*>         the dimension of IWORK must be at least<br>
*>                      6 + 6*N + 5*N*lg N<br>
*>                      ( lg( N ) = smallest integer k<br>
*>                                  such that 2^k >= N )<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array,<br>
*>                               dimension (1 + 3*N + 2*N*lg N + 3*N**2)<br>
*>                        ( lg( N ) = smallest integer k<br>
*>                                    such that 2^k >= N )<br>
*> \endverbatim<br>
*><br>
*> \param[out] QSTORE<br>
*> \verbatim<br>
*>          QSTORE is COMPLEX array, dimension (LDQS, N)<br>
*>         Used to store parts of<br>
*>         the eigenvector matrix when the updating matrix multiplies<br>
*>         take place.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQS<br>
*> \verbatim<br>
*>          LDQS is INTEGER<br>
*>         The leading dimension of the array QSTORE.<br>
*>         LDQS >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  The algorithm failed to compute an eigenvalue while<br>
*>                working on the submatrix lying in rows and columns<br>
*>                INFO/(N+1) through mod(INFO,N+1).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void claed0_(INTEGER QSIZ,INTEGER N,float[] D,float[] E,float[] Q,INTEGER LDQ,float[] QSTORE,INTEGER LDQS,float[] RWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b CLAED7 used by sstedc. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is dense.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAED7 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claed7.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claed7.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claed7.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,<br>
*                          LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM,<br>
*                          GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ,<br>
*      $                   TLVLS<br>
*       REAL               RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),<br>
*      $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )<br>
*       REAL               D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * )<br>
*       COMPLEX            Q( LDQ, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAED7 computes the updated eigensystem of a diagonal<br>
*> matrix after modification by a rank-one symmetric matrix. This<br>
*> routine is used only for the eigenproblem which requires all<br>
*> eigenvalues and optionally eigenvectors of a dense or banded<br>
*> Hermitian matrix that has been reduced to tridiagonal form.<br>
*><br>
*>   T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out)<br>
*><br>
*>   where Z = Q**Hu, u is a vector of length N with ones in the<br>
*>   CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.<br>
*><br>
*>    The eigenvectors of the original matrix are stored in Q, and the<br>
*>    eigenvalues are in D.  The algorithm consists of three stages:<br>
*><br>
*>       The first stage consists of deflating the size of the problem<br>
*>       when there are multiple eigenvalues or if there is a zero in<br>
*>       the Z vector.  For each such occurrence the dimension of the<br>
*>       secular equation problem is reduced by one.  This stage is<br>
*>       performed by the routine SLAED2.<br>
*><br>
*>       The second stage consists of calculating the updated<br>
*>       eigenvalues. This is done by finding the roots of the secular<br>
*>       equation via the routine SLAED4 (as called by SLAED3).<br>
*>       This routine also calculates the eigenvectors of the current<br>
*>       problem.<br>
*><br>
*>       The final stage consists of computing the updated eigenvectors<br>
*>       directly using the updated eigenvalues.  The eigenvectors for<br>
*>       the current problem are multiplied with the eigenvectors from<br>
*>       the overall problem.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The dimension of the symmetric tridiagonal matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CUTPNT<br>
*> \verbatim<br>
*>          CUTPNT is INTEGER<br>
*>         Contains the location of the last eigenvalue in the leading<br>
*>         sub-matrix.  min(1,N) <= CUTPNT <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] QSIZ<br>
*> \verbatim<br>
*>          QSIZ is INTEGER<br>
*>         The dimension of the unitary matrix used to reduce<br>
*>         the full matrix to tridiagonal form.  QSIZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TLVLS<br>
*> \verbatim<br>
*>          TLVLS is INTEGER<br>
*>         The total number of merging levels in the overall divide and<br>
*>         conquer tree.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CURLVL<br>
*> \verbatim<br>
*>          CURLVL is INTEGER<br>
*>         The current level in the overall merge routine,<br>
*>         0 <= curlvl <= tlvls.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CURPBM<br>
*> \verbatim<br>
*>          CURPBM is INTEGER<br>
*>         The current problem in the current level in the overall<br>
*>         merge routine (counting from upper left to lower right).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry, the eigenvalues of the rank-1-perturbed matrix.<br>
*>         On exit, the eigenvalues of the repaired matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX array, dimension (LDQ,N)<br>
*>         On entry, the eigenvectors of the rank-1-perturbed matrix.<br>
*>         On exit, the eigenvectors of the repaired tridiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>         The leading dimension of the array Q.  LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         Contains the subdiagonal element used to create the rank-1<br>
*>         modification.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDXQ<br>
*> \verbatim<br>
*>          INDXQ is INTEGER array, dimension (N)<br>
*>         This contains the permutation which will reintegrate the<br>
*>         subproblem just solved back into sorted order,<br>
*>         ie. D( INDXQ( I = 1, N ) ) will be in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array,<br>
*>                                 dimension (3*N+2*QSIZ*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (QSIZ*N)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] QSTORE<br>
*> \verbatim<br>
*>          QSTORE is REAL array, dimension (N**2+1)<br>
*>         Stores eigenvectors of submatrices encountered during<br>
*>         divide and conquer, packed together. QPTR points to<br>
*>         beginning of the submatrices.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] QPTR<br>
*> \verbatim<br>
*>          QPTR is INTEGER array, dimension (N+2)<br>
*>         List of indices pointing to beginning of submatrices stored<br>
*>         in QSTORE. The submatrices are numbered starting at the<br>
*>         bottom left of the divide and conquer tree, from left to<br>
*>         right and bottom to top.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PRMPTR<br>
*> \verbatim<br>
*>          PRMPTR is INTEGER array, dimension (N lg N)<br>
*>         Contains a list of pointers which indicate where in PERM a<br>
*>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)<br>
*>         indicates the size of the permutation and also the size of<br>
*>         the full, non-deflated problem.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PERM<br>
*> \verbatim<br>
*>          PERM is INTEGER array, dimension (N lg N)<br>
*>         Contains the permutations (from deflation and sorting) to be<br>
*>         applied to each eigenblock.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVPTR<br>
*> \verbatim<br>
*>          GIVPTR is INTEGER array, dimension (N lg N)<br>
*>         Contains a list of pointers which indicate where in GIVCOL a<br>
*>         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)<br>
*>         indicates the number of Givens rotations.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVCOL<br>
*> \verbatim<br>
*>          GIVCOL is INTEGER array, dimension (2, N lg N)<br>
*>         Each pair of numbers indicates a pair of columns to take place<br>
*>         in a Givens rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVNUM<br>
*> \verbatim<br>
*>          GIVNUM is REAL array, dimension (2, N lg N)<br>
*>         Each number indicates the S value to be used in the<br>
*>         corresponding Givens rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = 1, an eigenvalue did not converge<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void claed7_(INTEGER N,INTEGER CUTPNT,INTEGER QSIZ,INTEGER TLVLS,INTEGER CURLVL,INTEGER CURPBM,float[] D,float[] Q,INTEGER LDQ,REAL RHO,int[] INDXQ,float[] QSTORE,int[] QPTR,int[] PRMPTR,int[] PERM,int[] GIVPTR,int[] GIVCOL,float[] GIVNUM,float[] WORK,float[] RWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b CLAED8 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original matrix is dense.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAED8 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claed8.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claed8.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claed8.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMDA,<br>
*                          Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR,<br>
*                          GIVCOL, GIVNUM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ<br>
*       REAL               RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),<br>
*      $                   INDXQ( * ), PERM( * )<br>
*       REAL               D( * ), DLAMDA( * ), GIVNUM( 2, * ), W( * ),<br>
*      $                   Z( * )<br>
*       COMPLEX            Q( LDQ, * ), Q2( LDQ2, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAED8 merges the two sets of eigenvalues together into a single<br>
*> sorted set.  Then it tries to deflate the size of the problem.<br>
*> There are two ways in which deflation can occur:  when two or more<br>
*> eigenvalues are close together or if there is a tiny element in the<br>
*> Z vector.  For each such occurrence the order of the related secular<br>
*> equation problem is reduced by one.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[out] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>         Contains the number of non-deflated eigenvalues.<br>
*>         This is the order of the related secular equation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The dimension of the symmetric tridiagonal matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] QSIZ<br>
*> \verbatim<br>
*>          QSIZ is INTEGER<br>
*>         The dimension of the unitary matrix used to reduce<br>
*>         the dense or band matrix to tridiagonal form.<br>
*>         QSIZ >= N if ICOMPQ = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX array, dimension (LDQ,N)<br>
*>         On entry, Q contains the eigenvectors of the partially solved<br>
*>         system which has been previously updated in matrix<br>
*>         multiplies with other partially solved eigensystems.<br>
*>         On exit, Q contains the trailing (N-K) updated eigenvectors<br>
*>         (those which were deflated) in its last N-K columns.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>         The leading dimension of the array Q.  LDQ >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry, D contains the eigenvalues of the two submatrices to<br>
*>         be combined.  On exit, D contains the trailing (N-K) updated<br>
*>         eigenvalues (those which were deflated) sorted into increasing<br>
*>         order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         Contains the off diagonal element associated with the rank-1<br>
*>         cut which originally split the two submatrices which are now<br>
*>         being recombined. RHO is modified during the computation to<br>
*>         the value required by SLAED3.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CUTPNT<br>
*> \verbatim<br>
*>          CUTPNT is INTEGER<br>
*>         Contains the location of the last eigenvalue in the leading<br>
*>         sub-matrix.  MIN(1,N) <= CUTPNT <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (N)<br>
*>         On input this vector contains the updating vector (the last<br>
*>         row of the first sub-eigenvector matrix and the first row of<br>
*>         the second sub-eigenvector matrix).  The contents of Z are<br>
*>         destroyed during the updating process.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DLAMDA<br>
*> \verbatim<br>
*>          DLAMDA is REAL array, dimension (N)<br>
*>         Contains a copy of the first K eigenvalues which will be used<br>
*>         by SLAED3 to form the secular equation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q2<br>
*> \verbatim<br>
*>          Q2 is COMPLEX array, dimension (LDQ2,N)<br>
*>         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,<br>
*>         Contains a copy of the first K eigenvectors which will be used<br>
*>         by SLAED7 in a matrix multiply (SGEMM) to update the new<br>
*>         eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ2<br>
*> \verbatim<br>
*>          LDQ2 is INTEGER<br>
*>         The leading dimension of the array Q2.  LDQ2 >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>         This will hold the first k values of the final<br>
*>         deflation-altered z-vector and will be passed to SLAED3.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDXP<br>
*> \verbatim<br>
*>          INDXP is INTEGER array, dimension (N)<br>
*>         This will contain the permutation used to place deflated<br>
*>         values of D at the end of the array. On output INDXP(1:K)<br>
*>         points to the nondeflated D-values and INDXP(K+1:N)<br>
*>         points to the deflated eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDX<br>
*> \verbatim<br>
*>          INDX is INTEGER array, dimension (N)<br>
*>         This will contain the permutation used to sort the contents of<br>
*>         D into ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INDXQ<br>
*> \verbatim<br>
*>          INDXQ is INTEGER array, dimension (N)<br>
*>         This contains the permutation which separately sorts the two<br>
*>         sub-problems in D into ascending order.  Note that elements in<br>
*>         the second half of this permutation must first have CUTPNT<br>
*>         added to their values in order to be accurate.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PERM<br>
*> \verbatim<br>
*>          PERM is INTEGER array, dimension (N)<br>
*>         Contains the permutations (from deflation and sorting) to be<br>
*>         applied to each eigenblock.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVPTR<br>
*> \verbatim<br>
*>          GIVPTR is INTEGER<br>
*>         Contains the number of Givens rotations which took place in<br>
*>         this subproblem.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVCOL<br>
*> \verbatim<br>
*>          GIVCOL is INTEGER array, dimension (2, N)<br>
*>         Each pair of numbers indicates a pair of columns to take place<br>
*>         in a Givens rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVNUM<br>
*> \verbatim<br>
*>          GIVNUM is REAL array, dimension (2, N)<br>
*>         Each number indicates the S value to be used in the<br>
*>         corresponding Givens rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void claed8_(INTEGER K,INTEGER N,INTEGER QSIZ,float[] Q,INTEGER LDQ,float[] D,REAL RHO,INTEGER CUTPNT,float[] Z,float[] DLAMDA,float[] Q2,INTEGER LDQ2,float[] W,int[] INDXP,int[] INDX,int[] INDXQ,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,float[] GIVNUM,INTEGER INFO);
/**
*> \brief \b CLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAEIN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claein.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claein.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claein.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK,<br>
*                          EPS3, SMLNUM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            NOINIT, RIGHTV<br>
*       INTEGER            INFO, LDB, LDH, N<br>
*       REAL               EPS3, SMLNUM<br>
*       COMPLEX            W<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               RWORK( * )<br>
*       COMPLEX            B( LDB, * ), H( LDH, * ), V( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAEIN uses inverse iteration to find a right or left eigenvector<br>
*> corresponding to the eigenvalue W of a complex upper Hessenberg<br>
*> matrix H.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] RIGHTV<br>
*> \verbatim<br>
*>          RIGHTV is LOGICAL<br>
*>          = .TRUE. : compute right eigenvector;<br>
*>          = .FALSE.: compute left eigenvector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NOINIT<br>
*> \verbatim<br>
*>          NOINIT is LOGICAL<br>
*>          = .TRUE. : no initial vector supplied in V<br>
*>          = .FALSE.: initial vector supplied in V.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix H.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] H<br>
*> \verbatim<br>
*>          H is COMPLEX array, dimension (LDH,N)<br>
*>          The upper Hessenberg matrix H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>          The leading dimension of the array H.  LDH >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] W<br>
*> \verbatim<br>
*>          W is COMPLEX<br>
*>          The eigenvalue of H whose corresponding right or left<br>
*>          eigenvector is to be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension (N)<br>
*>          On entry, if NOINIT = .FALSE., V must contain a starting<br>
*>          vector for inverse iteration; otherwise V need not be set.<br>
*>          On exit, V contains the computed eigenvector, normalized so<br>
*>          that the component of largest magnitude has magnitude 1; here<br>
*>          the magnitude of a complex number (x,y) is taken to be<br>
*>          |x| + |y|.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,N)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[in] EPS3<br>
*> \verbatim<br>
*>          EPS3 is REAL<br>
*>          A small machine-dependent value which is used to perturb<br>
*>          close eigenvalues, and to replace zero pivots.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SMLNUM<br>
*> \verbatim<br>
*>          SMLNUM is REAL<br>
*>          A machine-dependent value close to the underflow threshold.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          = 1:  inverse iteration did not converge; V is set to the<br>
*>                last iterate.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claein_(LOGICAL RIGHTV,LOGICAL NOINIT,INTEGER N,float[] H,INTEGER LDH,float[] W,float[] V,float[] B,INTEGER LDB,float[] RWORK,REAL EPS3,REAL SMLNUM,INTEGER INFO);
/**
*> \brief \b CLAESY computes the eigenvalues and eigenvectors of a 2-by-2 complex symmetric matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAESY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claesy.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claesy.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claesy.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX            A, B, C, CS1, EVSCAL, RT1, RT2, SN1<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix<br>
*>    ( ( A, B );( B, C ) )<br>
*> provided the norm of the matrix of eigenvectors is larger than<br>
*> some threshold value.<br>
*><br>
*> RT1 is the eigenvalue of larger absolute value, and RT2 of<br>
*> smaller absolute value.  If the eigenvectors are computed, then<br>
*> on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence<br>
*><br>
*> [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ]<br>
*> [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ]<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX<br>
*>          The ( 1, 1 ) element of input matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX<br>
*>          The ( 1, 2 ) element of input matrix.  The ( 2, 1 ) element<br>
*>          is also given by B, since the 2-by-2 matrix is symmetric.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is COMPLEX<br>
*>          The ( 2, 2 ) element of input matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT1<br>
*> \verbatim<br>
*>          RT1 is COMPLEX<br>
*>          The eigenvalue of larger modulus.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT2<br>
*> \verbatim<br>
*>          RT2 is COMPLEX<br>
*>          The eigenvalue of smaller modulus.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EVSCAL<br>
*> \verbatim<br>
*>          EVSCAL is COMPLEX<br>
*>          The complex value by which the eigenvector matrix was scaled<br>
*>          to make it orthonormal.  If EVSCAL is zero, the eigenvectors<br>
*>          were not computed.  This means one of two things:  the 2-by-2<br>
*>          matrix could not be diagonalized, or the norm of the matrix<br>
*>          of eigenvectors before scaling was larger than the threshold<br>
*>          value THRESH (set below).<br>
*> \endverbatim<br>
*><br>
*> \param[out] CS1<br>
*> \verbatim<br>
*>          CS1 is COMPLEX<br>
*> \endverbatim<br>
*><br>
*> \param[out] SN1<br>
*> \verbatim<br>
*>          SN1 is COMPLEX<br>
*>          If EVSCAL .NE. 0,  ( CS1, SN1 ) is the unit right eigenvector<br>
*>          for RT1.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexSYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claesy_(float[] A,float[] B,float[] C,float[] RT1,float[] RT2,float[] EVSCAL,float[] CS1,float[] SN1);
/**
*> \brief \b CLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAEV2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claev2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claev2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claev2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAEV2( A, B, C, RT1, RT2, CS1, SN1 )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               CS1, RT1, RT2<br>
*       COMPLEX            A, B, C, SN1<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAEV2 computes the eigendecomposition of a 2-by-2 Hermitian matrix<br>
*>    [  A         B  ]<br>
*>    [  CONJG(B)  C  ].<br>
*> On return, RT1 is the eigenvalue of larger absolute value, RT2 is the<br>
*> eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right<br>
*> eigenvector for RT1, giving the decomposition<br>
*><br>
*> [ CS1  CONJG(SN1) ] [    A     B ] [ CS1 -CONJG(SN1) ] = [ RT1  0  ]<br>
*> [-SN1     CS1     ] [ CONJG(B) C ] [ SN1     CS1     ]   [  0  RT2 ].<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX<br>
*>         The (1,1) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX<br>
*>         The (1,2) element and the conjugate of the (2,1) element of<br>
*>         the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is COMPLEX<br>
*>         The (2,2) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT1<br>
*> \verbatim<br>
*>          RT1 is REAL<br>
*>         The eigenvalue of larger absolute value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT2<br>
*> \verbatim<br>
*>          RT2 is REAL<br>
*>         The eigenvalue of smaller absolute value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CS1<br>
*> \verbatim<br>
*>          CS1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] SN1<br>
*> \verbatim<br>
*>          SN1 is COMPLEX<br>
*>         The vector (CS1, SN1) is a unit right eigenvector for RT1.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  RT1 is accurate to a few ulps barring over/underflow.<br>
*><br>
*>  RT2 may be inaccurate if there is massive cancellation in the<br>
*>  determinant A*C-B*B; higher precision or correctly rounded or<br>
*>  correctly truncated arithmetic would be needed to compute RT2<br>
*>  accurately in all cases.<br>
*><br>
*>  CS1 and SN1 are accurate to a few ulps barring over/underflow.<br>
*><br>
*>  Overflow is possible only if RT1 is within a factor of 5 of overflow.<br>
*>  Underflow is harmless if the input data is 0 or exceeds<br>
*>     underflow_threshold / macheps.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void claev2_(float[] A,float[] B,float[] C,REAL RT1,REAL RT2,REAL CS1,float[] SN1);
/**
*> \brief \b CLAG2Z converts a complex single precision matrix to a complex double precision matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAG2Z + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clag2z.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clag2z.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clag2z.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAG2Z( M, N, SA, LDSA, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDSA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            SA( LDSA, * )<br>
*       COMPLEX*16         A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAG2Z converts a COMPLEX matrix, SA, to a COMPLEX*16 matrix, A.<br>
*><br>
*> Note that while it is possible to overflow while converting<br>
*> from double to single, it is not possible to overflow when<br>
*> converting from single to double.<br>
*><br>
*> This is an auxiliary routine so there is no argument checking.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of lines of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SA<br>
*> \verbatim<br>
*>          SA is COMPLEX array, dimension (LDSA,N)<br>
*>          On entry, the M-by-N coefficient matrix SA.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDSA<br>
*> \verbatim<br>
*>          LDSA is INTEGER<br>
*>          The leading dimension of the array SA.  LDSA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On exit, the M-by-N coefficient matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complex16OTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clag2z_(INTEGER M,INTEGER N,float[] SA,INTEGER LDSA,double[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b CLAGS2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAGS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clags2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clags2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clags2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV,<br>
*                          SNV, CSQ, SNQ )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            UPPER<br>
*       REAL               A1, A3, B1, B3, CSQ, CSU, CSV<br>
*       COMPLEX            A2, B2, SNQ, SNU, SNV<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAGS2 computes 2-by-2 unitary matrices U, V and Q, such<br>
*> that if ( UPPER ) then<br>
*><br>
*>           U**H *A*Q = U**H *( A1 A2 )*Q = ( x  0  )<br>
*>                             ( 0  A3 )     ( x  x  )<br>
*> and<br>
*>           V**H*B*Q = V**H *( B1 B2 )*Q = ( x  0  )<br>
*>                            ( 0  B3 )     ( x  x  )<br>
*><br>
*> or if ( .NOT.UPPER ) then<br>
*><br>
*>           U**H *A*Q = U**H *( A1 0  )*Q = ( x  x  )<br>
*>                             ( A2 A3 )     ( 0  x  )<br>
*> and<br>
*>           V**H *B*Q = V**H *( B1 0  )*Q = ( x  x  )<br>
*>                             ( B2 B3 )     ( 0  x  )<br>
*> where<br>
*><br>
*>   U = (   CSU    SNU ), V = (  CSV    SNV ),<br>
*>       ( -SNU**H  CSU )      ( -SNV**H CSV )<br>
*><br>
*>   Q = (   CSQ    SNQ )<br>
*>       ( -SNQ**H  CSQ )<br>
*><br>
*> The rows of the transformed A and B are parallel. Moreover, if the<br>
*> input 2-by-2 matrix A is not zero, then the transformed (1,1) entry<br>
*> of A is not zero. If the input matrices A and B are both not zero,<br>
*> then the transformed (2,2) element of B is not zero, except when the<br>
*> first rows of input A and B are parallel and the second rows are<br>
*> zero.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPPER<br>
*> \verbatim<br>
*>          UPPER is LOGICAL<br>
*>          = .TRUE.: the input matrices A and B are upper triangular.<br>
*>          = .FALSE.: the input matrices A and B are lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A1<br>
*> \verbatim<br>
*>          A1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] A2<br>
*> \verbatim<br>
*>          A2 is COMPLEX<br>
*> \endverbatim<br>
*><br>
*> \param[in] A3<br>
*> \verbatim<br>
*>          A3 is REAL<br>
*>          On entry, A1, A2 and A3 are elements of the input 2-by-2<br>
*>          upper (lower) triangular matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B1<br>
*> \verbatim<br>
*>          B1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] B2<br>
*> \verbatim<br>
*>          B2 is COMPLEX<br>
*> \endverbatim<br>
*><br>
*> \param[in] B3<br>
*> \verbatim<br>
*>          B3 is REAL<br>
*>          On entry, B1, B2 and B3 are elements of the input 2-by-2<br>
*>          upper (lower) triangular matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSU<br>
*> \verbatim<br>
*>          CSU is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNU<br>
*> \verbatim<br>
*>          SNU is COMPLEX<br>
*>          The desired unitary matrix U.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSV<br>
*> \verbatim<br>
*>          CSV is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNV<br>
*> \verbatim<br>
*>          SNV is COMPLEX<br>
*>          The desired unitary matrix V.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSQ<br>
*> \verbatim<br>
*>          CSQ is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNQ<br>
*> \verbatim<br>
*>          SNQ is COMPLEX<br>
*>          The desired unitary matrix Q.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2011<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clags2_(LOGICAL UPPER,REAL A1,float[] A2,REAL A3,REAL B1,float[] B2,REAL B3,REAL CSU,float[] SNU,REAL CSV,float[] SNV,REAL CSQ,float[] SNQ);
/**
*> \brief \b CLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matrix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAGTM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clagtm.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clagtm.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clagtm.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA,<br>
*                          B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            LDB, LDX, N, NRHS<br>
*       REAL               ALPHA, BETA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAGTM performs a matrix-vector product of the form<br>
*><br>
*>    B := alpha * A * X + beta * B<br>
*><br>
*> where A is a tridiagonal matrix of order N, B and X are N by NRHS<br>
*> matrices, and alpha and beta are real scalars, each of which may be<br>
*> 0., 1., or -1.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the operation applied to A.<br>
*>          = 'N':  No transpose, B := alpha * A * X + beta * B<br>
*>          = 'T':  Transpose,    B := alpha * A**T * X + beta * B<br>
*>          = 'C':  Conjugate transpose, B := alpha * A**H * X + beta * B<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrices X and B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise,<br>
*>          it is assumed to be 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DL<br>
*> \verbatim<br>
*>          DL is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) sub-diagonal elements of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is COMPLEX array, dimension (N)<br>
*>          The diagonal elements of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) super-diagonal elements of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>          The N by NRHS matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(N,1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,<br>
*>          it is assumed to be 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>          On entry, the N by NRHS matrix B.<br>
*>          On exit, B is overwritten by the matrix expression<br>
*>          B := alpha * A * X + beta * B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(N,1).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clagtm_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,REAL ALPHA,float[] DL,float[] D,float[] DU,float[] X,INTEGER LDX,REAL BETA,float[] B,INTEGER LDB);
/**
*> \brief \b CLAHEF computes a partial factorization of a complex Hermitian indefinite matrix using the Bunch-Kaufman diagonal pivoting method (blocked algorithm, calling Level 3 BLAS).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download CLAHEF + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAHEF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KB, LDA, LDW, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), W( LDW, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAHEF computes a partial factorization of a complex Hermitian<br>
*> matrix A using the Bunch-Kaufman diagonal pivoting method. The<br>
*> partial factorization has the form:<br>
*><br>
*> A  =  ( I  U12 ) ( A11  0  ) (  I      0     )  if UPLO = 'U', or:<br>
*>       ( 0  U22 ) (  0   D  ) ( U12**H U22**H )<br>
*><br>
*> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L'<br>
*>       ( L21  I ) (  0  A22 ) (  0      I     )<br>
*><br>
*> where the order of D is at most NB. The actual order is returned in<br>
*> the argument KB, and is either NB or NB-1, or N if N <= NB.<br>
*> Note that U**H denotes the conjugate transpose of U.<br>
*><br>
*> CLAHEF is an auxiliary routine called by CHETRF. It uses blocked code<br>
*> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or<br>
*> A22 (if UPLO = 'L').<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          Hermitian matrix A is stored:<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The maximum number of columns of the matrix A that should be<br>
*>          factored.  NB should be at least 2 to allow for 2-by-2 pivot<br>
*>          blocks.<br>
*> \endverbatim<br>
*><br>
*> \param[out] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of columns of A that were actually factored.<br>
*>          KB is either NB-1 or NB, or N if N <= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit, A contains details of the partial factorization.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D.<br>
*><br>
*>          If UPLO = 'U':<br>
*>             Only the last KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns<br>
*>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)<br>
*>             is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>             Only the first KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns<br>
*>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)<br>
*>             is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (LDW,NB)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDW<br>
*> \verbatim<br>
*>          LDW is INTEGER<br>
*>          The leading dimension of the array W.  LDW >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complexHEcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>  November 2013,  Igor Kozachenko,<br>
*>                  Computer Science Division,<br>
*>                  University of California, Berkeley<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void clahef_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] IPIV,float[] W,INTEGER LDW,INTEGER INFO);
/**
* \brief \b CLAHEF_ROOK computes a partial factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (blocked algorithm, calling Level 3 BLAS).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download CLAHEF_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAHEF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KB, LDA, LDW, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), W( LDW, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAHEF_ROOK computes a partial factorization of a complex Hermitian<br>
*> matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting<br>
*> method. The partial factorization has the form:<br>
*><br>
*> A  =  ( I  U12 ) ( A11  0  ) (  I      0     )  if UPLO = 'U', or:<br>
*>       ( 0  U22 ) (  0   D  ) ( U12**H U22**H )<br>
*><br>
*> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L'<br>
*>       ( L21  I ) (  0  A22 ) (  0      I     )<br>
*><br>
*> where the order of D is at most NB. The actual order is returned in<br>
*> the argument KB, and is either NB or NB-1, or N if N <= NB.<br>
*> Note that U**H denotes the conjugate transpose of U.<br>
*><br>
*> CLAHEF_ROOK is an auxiliary routine called by CHETRF_ROOK. It uses<br>
*> blocked code (calling Level 3 BLAS) to update the submatrix<br>
*> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          Hermitian matrix A is stored:<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The maximum number of columns of the matrix A that should be<br>
*>          factored.  NB should be at least 2 to allow for 2-by-2 pivot<br>
*>          blocks.<br>
*> \endverbatim<br>
*><br>
*> \param[out] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of columns of A that were actually factored.<br>
*>          KB is either NB-1 or NB, or N if N <= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit, A contains details of the partial factorization.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D.<br>
*><br>
*>          If UPLO = 'U':<br>
*>             Only the last KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k-1 and -IPIV(k-1) were inerchaged,<br>
*>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>             Only the first KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>             were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k+1 and -IPIV(k+1) were inerchaged,<br>
*>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (LDW,NB)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDW<br>
*> \verbatim<br>
*>          LDW is INTEGER<br>
*>          The leading dimension of the array W.  LDW >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complexHEcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>  November 2013, Igor Kozachenko,<br>
*>                  Computer Science Division,<br>
*>                  University of California, Berkeley<br>
*><br>
*>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,<br>
*>                  School of Mathematics,<br>
*>                  University of Manchester<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void clahef_rook_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] IPIV,float[] W,INTEGER LDW,INTEGER INFO);
/**
*> \brief \b CLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAHQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,<br>
*                          IHIZ, Z, LDZ, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLAHQR is an auxiliary routine called by CHSEQR to update the<br>
*>    eigenvalues and Schur decomposition already computed by CHSEQR, by<br>
*>    dealing with the Hessenberg submatrix in rows and columns ILO to<br>
*>    IHI.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTT<br>
*> \verbatim<br>
*>          WANTT is LOGICAL<br>
*>          = .TRUE. : the full Schur form T is required;<br>
*>          = .FALSE.: only eigenvalues are required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          = .TRUE. : the matrix of Schur vectors Z is required;<br>
*>          = .FALSE.: Schur vectors are not required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix H.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILO<br>
*> \verbatim<br>
*>          ILO is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHI<br>
*> \verbatim<br>
*>          IHI is INTEGER<br>
*>          It is assumed that H is already upper triangular in rows and<br>
*>          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).<br>
*>          CLAHQR works primarily with the Hessenberg submatrix in rows<br>
*>          and columns ILO to IHI, but applies transformations to all of<br>
*>          H if WANTT is .TRUE..<br>
*>          1 <= ILO <= max(1,IHI); IHI <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is COMPLEX array, dimension (LDH,N)<br>
*>          On entry, the upper Hessenberg matrix H.<br>
*>          On exit, if INFO is zero and if WANTT is .TRUE., then H<br>
*>          is upper triangular in rows and columns ILO:IHI.  If INFO<br>
*>          is zero and if WANTT is .FALSE., then the contents of H<br>
*>          are unspecified on exit.  The output state of H in case<br>
*>          INF is positive is below under the description of INFO.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>          The leading dimension of the array H. LDH >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (N)<br>
*>          The computed eigenvalues ILO to IHI are stored in the<br>
*>          corresponding elements of W. If WANTT is .TRUE., the<br>
*>          eigenvalues are stored in the same order as on the diagonal<br>
*>          of the Schur form returned in H, with W(i) = H(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILOZ<br>
*> \verbatim<br>
*>          ILOZ is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHIZ<br>
*> \verbatim<br>
*>          IHIZ is INTEGER<br>
*>          Specify the rows of Z to which transformations must be<br>
*>          applied if WANTZ is .TRUE..<br>
*>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (LDZ,N)<br>
*>          If WANTZ is .TRUE., on entry Z must contain the current<br>
*>          matrix Z of transformations accumulated by CHSEQR, and on<br>
*>          exit Z has been updated; transformations are applied only to<br>
*>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).<br>
*>          If WANTZ is .FALSE., Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z. LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           =   0: successful exit<br>
*>          .GT. 0: if INFO = i, CLAHQR failed to compute all the<br>
*>                  eigenvalues ILO to IHI in a total of 30 iterations<br>
*>                  per eigenvalue; elements i+1:ihi of W contain<br>
*>                  those eigenvalues which have been successfully<br>
*>                  computed.<br>
*><br>
*>                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,<br>
*>                  the remaining unconverged eigenvalues are the<br>
*>                  eigenvalues of the upper Hessenberg matrix<br>
*>                  rows and columns ILO thorugh INFO of the final,<br>
*>                  output value of H.<br>
*><br>
*>                  If INFO .GT. 0 and WANTT is .TRUE., then on exit<br>
*>          (*)       (initial value of H)*U  = U*(final value of H)<br>
*>                  where U is an orthognal matrix.    The final<br>
*>                  value of H is upper Hessenberg and triangular in<br>
*>                  rows and columns INFO+1 through IHI.<br>
*><br>
*>                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit<br>
*>                      (final value of Z)  = (initial value of Z)*U<br>
*>                  where U is the orthogonal matrix in (*)<br>
*>                  (regardless of the value of WANTT.)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2015<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>     02-96 Based on modifications by<br>
*>     David Day, Sandia National Laboratory, USA<br>
*><br>
*>     12-04 Further modifications by<br>
*>     Ralph Byers, University of Kansas, USA<br>
*>     This is a modified version of CLAHQR from LAPACK version 3.0.<br>
*>     It is (1) more robust against overflow and underflow and<br>
*>     (2) adopts the more conservative Ahues & Tisseur stopping<br>
*>     criterion (LAWN 122, 1997).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clahqr_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] W,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER INFO);
/**
*> \brief \b CLAHR2 reduces the specified number of first columns of a general rectangular matrix A so that elements below the specified subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformation to the unreduced part of A.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAHR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            K, LDA, LDT, LDY, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), T( LDT, NB ), TAU( NB ),<br>
*      $                   Y( LDY, NB )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAHR2 reduces the first NB columns of A complex general n-BY-(n-k+1)<br>
*> matrix A so that elements below the k-th subdiagonal are zero. The<br>
*> reduction is performed by an unitary similarity transformation<br>
*> Q**H * A * Q. The routine returns the matrices V and T which determine<br>
*> Q as a block reflector I - V*T*v**H, and also the matrix Y = A * V * T.<br>
*><br>
*> This is an auxiliary routine called by CGEHRD.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The offset for the reduction. Elements below the k-th<br>
*>          subdiagonal in the first NB columns are reduced to zero.<br>
*>          K < N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The number of columns to be reduced.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N-K+1)<br>
*>          On entry, the n-by-(n-k+1) general matrix A.<br>
*>          On exit, the elements on and above the k-th subdiagonal in<br>
*>          the first NB columns are overwritten with the corresponding<br>
*>          elements of the reduced matrix; the elements below the k-th<br>
*>          subdiagonal, with the array TAU, represent the matrix Q as a<br>
*>          product of elementary reflectors. The other columns of A are<br>
*>          unchanged. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX array, dimension (NB)<br>
*>          The scalar factors of the elementary reflectors. See Further<br>
*>          Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is COMPLEX array, dimension (LDT,NB)<br>
*>          The upper triangular matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (LDY,NB)<br>
*>          The n-by-nb matrix Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDY<br>
*> \verbatim<br>
*>          LDY is INTEGER<br>
*>          The leading dimension of the array Y. LDY >= N.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of nb elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(nb).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in<br>
*>  A(i+k+1:n,i), and tau in TAU(i).<br>
*><br>
*>  The elements of the vectors v together form the (n-k+1)-by-nb matrix<br>
*>  V which is needed, with T and Y, to apply the transformation to the<br>
*>  unreduced part of the matrix, using an update of the form:<br>
*>  A := (I - V*T*V**H) * (A - Y*V**H).<br>
*><br>
*>  The contents of A on exit are illustrated by the following example<br>
*>  with n = 7, k = 3 and nb = 2:<br>
*><br>
*>     ( a   a   a   a   a )<br>
*>     ( a   a   a   a   a )<br>
*>     ( a   a   a   a   a )<br>
*>     ( h   h   a   a   a )<br>
*>     ( v1  h   a   a   a )<br>
*>     ( v1  v2  a   a   a )<br>
*>     ( v1  v2  a   a   a )<br>
*><br>
*>  where a denotes an element of the original matrix A, h denotes a<br>
*>  modified element of the upper Hessenberg matrix H, and vi denotes an<br>
*>  element of the vector defining H(i).<br>
*><br>
*>  This subroutine is a slight modification of LAPACK-3.0's DLAHRD<br>
*>  incorporating improvements proposed by Quintana-Orti and Van de<br>
*>  Gejin. Note that the entries of A(1:K,2:NB) differ from those<br>
*>  returned by the original LAPACK-3.0's DLAHRD routine. (This<br>
*>  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.)<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  Gregorio Quintana-Orti and Robert van de Geijn, "Improving the<br>
*>  performance of reduction to Hessenberg form," ACM Transactions on<br>
*>  Mathematical Software, 32(2):180-194, June 2006.<br>
*><br>
*  =====================================================================<br>
*/
	public void clahr2_(INTEGER N,INTEGER K,INTEGER NB,float[] A,INTEGER LDA,float[] TAU,float[] T,INTEGER LDT,float[] Y,INTEGER LDY);
/**
*> \brief \b CLAIC1 applies one step of incremental condition estimation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAIC1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claic1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claic1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claic1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            J, JOB<br>
*       REAL               SEST, SESTPR<br>
*       COMPLEX            C, GAMMA, S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            W( J ), X( J )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAIC1 applies one step of incremental condition estimation in<br>
*> its simplest version:<br>
*><br>
*> Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j<br>
*> lower triangular matrix L, such that<br>
*>          twonorm(L*x) = sest<br>
*> Then CLAIC1 computes sestpr, s, c such that<br>
*> the vector<br>
*>                 [ s*x ]<br>
*>          xhat = [  c  ]<br>
*> is an approximate singular vector of<br>
*>                 [ L      0  ]<br>
*>          Lhat = [ w**H gamma ]<br>
*> in the sense that<br>
*>          twonorm(Lhat*xhat) = sestpr.<br>
*><br>
*> Depending on JOB, an estimate for the largest or smallest singular<br>
*> value is computed.<br>
*><br>
*> Note that [s c]**H and sestpr**2 is an eigenpair of the system<br>
*><br>
*>     diag(sest*sest, 0) + [alpha  gamma] * [ conjg(alpha) ]<br>
*>                                           [ conjg(gamma) ]<br>
*><br>
*> where  alpha =  x**H*w.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is INTEGER<br>
*>          = 1: an estimate for the largest singular value is computed.<br>
*>          = 2: an estimate for the smallest singular value is computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] J<br>
*> \verbatim<br>
*>          J is INTEGER<br>
*>          Length of X and W<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (J)<br>
*>          The j-vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SEST<br>
*> \verbatim<br>
*>          SEST is REAL<br>
*>          Estimated singular value of j by j matrix L<br>
*> \endverbatim<br>
*><br>
*> \param[in] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (J)<br>
*>          The j-vector w.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GAMMA<br>
*> \verbatim<br>
*>          GAMMA is COMPLEX<br>
*>          The diagonal element gamma.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SESTPR<br>
*> \verbatim<br>
*>          SESTPR is REAL<br>
*>          Estimated singular value of (j+1) by (j+1) matrix Lhat.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is COMPLEX<br>
*>          Sine needed in forming xhat.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is COMPLEX<br>
*>          Cosine needed in forming xhat.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claic1_(INTEGER JOB,INTEGER J,float[] X,REAL SEST,float[] W,float[] GAMMA,REAL SESTPR,float[] S,float[] C);
/**
*> \brief \b CLALS0 applies back multiplying factors in solving the least squares problem using divide and conquer SVD approach. Used by sgelsd.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLALS0 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clals0.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clals0.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clals0.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX,<br>
*                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,<br>
*                          POLES, DIFL, DIFR, Z, K, C, S, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL,<br>
*      $                   LDGNUM, NL, NR, NRHS, SQRE<br>
*       REAL               C, S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( LDGCOL, * ), PERM( * )<br>
*       REAL               DIFL( * ), DIFR( LDGNUM, * ),<br>
*      $                   GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ),<br>
*      $                   RWORK( * ), Z( * )<br>
*       COMPLEX            B( LDB, * ), BX( LDBX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLALS0 applies back the multiplying factors of either the left or the<br>
*> right singular vector matrix of a diagonal matrix appended by a row<br>
*> to the right hand side matrix B in solving the least squares problem<br>
*> using the divide-and-conquer SVD approach.<br>
*><br>
*> For the left singular vector matrix, three types of orthogonal<br>
*> matrices are involved:<br>
*><br>
*> (1L) Givens rotations: the number of such rotations is GIVPTR; the<br>
*>      pairs of columns/rows they were applied to are stored in GIVCOL;<br>
*>      and the C- and S-values of these rotations are stored in GIVNUM.<br>
*><br>
*> (2L) Permutation. The (NL+1)-st row of B is to be moved to the first<br>
*>      row, and for J=2:N, PERM(J)-th row of B is to be moved to the<br>
*>      J-th row.<br>
*><br>
*> (3L) The left singular vector matrix of the remaining matrix.<br>
*><br>
*> For the right singular vector matrix, four types of orthogonal<br>
*> matrices are involved:<br>
*><br>
*> (1R) The right singular vector matrix of the remaining matrix.<br>
*><br>
*> (2R) If SQRE = 1, one extra Givens rotation to generate the right<br>
*>      null space.<br>
*><br>
*> (3R) The inverse transformation of (2L).<br>
*><br>
*> (4R) The inverse transformation of (1L).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ICOMPQ<br>
*> \verbatim<br>
*>          ICOMPQ is INTEGER<br>
*>         Specifies whether singular vectors are to be computed in<br>
*>         factored form:<br>
*>         = 0: Left singular vector matrix.<br>
*>         = 1: Right singular vector matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NL<br>
*> \verbatim<br>
*>          NL is INTEGER<br>
*>         The row dimension of the upper block. NL >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NR<br>
*> \verbatim<br>
*>          NR is INTEGER<br>
*>         The row dimension of the lower block. NR >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SQRE<br>
*> \verbatim<br>
*>          SQRE is INTEGER<br>
*>         = 0: the lower block is an NR-by-NR square matrix.<br>
*>         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.<br>
*><br>
*>         The bidiagonal matrix has row dimension N = NL + NR + 1,<br>
*>         and column dimension M = N + SQRE.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>         The number of columns of B and BX. NRHS must be at least 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension ( LDB, NRHS )<br>
*>         On input, B contains the right hand sides of the least<br>
*>         squares problem in rows 1 through M. On output, B contains<br>
*>         the solution X in rows 1 through N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>         The leading dimension of B. LDB must be at least<br>
*>         max(1,MAX( M, N ) ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BX<br>
*> \verbatim<br>
*>          BX is COMPLEX array, dimension ( LDBX, NRHS )<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDBX<br>
*> \verbatim<br>
*>          LDBX is INTEGER<br>
*>         The leading dimension of BX.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PERM<br>
*> \verbatim<br>
*>          PERM is INTEGER array, dimension ( N )<br>
*>         The permutations (from deflation and sorting) applied<br>
*>         to the two blocks.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVPTR<br>
*> \verbatim<br>
*>          GIVPTR is INTEGER<br>
*>         The number of Givens rotations which took place in this<br>
*>         subproblem.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVCOL<br>
*> \verbatim<br>
*>          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 )<br>
*>         Each pair of numbers indicates a pair of rows/columns<br>
*>         involved in a Givens rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDGCOL<br>
*> \verbatim<br>
*>          LDGCOL is INTEGER<br>
*>         The leading dimension of GIVCOL, must be at least N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVNUM<br>
*> \verbatim<br>
*>          GIVNUM is REAL array, dimension ( LDGNUM, 2 )<br>
*>         Each number indicates the C or S value used in the<br>
*>         corresponding Givens rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDGNUM<br>
*> \verbatim<br>
*>          LDGNUM is INTEGER<br>
*>         The leading dimension of arrays DIFR, POLES and<br>
*>         GIVNUM, must be at least K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] POLES<br>
*> \verbatim<br>
*>          POLES is REAL array, dimension ( LDGNUM, 2 )<br>
*>         On entry, POLES(1:K, 1) contains the new singular<br>
*>         values obtained from solving the secular equation, and<br>
*>         POLES(1:K, 2) is an array containing the poles in the secular<br>
*>         equation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIFL<br>
*> \verbatim<br>
*>          DIFL is REAL array, dimension ( K ).<br>
*>         On entry, DIFL(I) is the distance between I-th updated<br>
*>         (undeflated) singular value and the I-th (undeflated) old<br>
*>         singular value.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIFR<br>
*> \verbatim<br>
*>          DIFR is REAL array, dimension ( LDGNUM, 2 ).<br>
*>         On entry, DIFR(I, 1) contains the distances between I-th<br>
*>         updated (undeflated) singular value and the I+1-th<br>
*>         (undeflated) old singular value. And DIFR(I, 2) is the<br>
*>         normalizing factor for the I-th right singular vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( K )<br>
*>         Contain the components of the deflation-adjusted updating row<br>
*>         vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>         Contains the dimension of the non-deflated matrix,<br>
*>         This is the order of the related secular equation. 1 <= K <=N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL<br>
*>         C contains garbage if SQRE =0 and the C-value of a Givens<br>
*>         rotation related to the right null space if SQRE = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL<br>
*>         S contains garbage if SQRE =0 and the S-value of a Givens<br>
*>         rotation related to the right null space if SQRE = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension<br>
*>         ( K*(1+NRHS) + 2*NRHS )<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2015<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Ren-Cang Li, Computer Science Division, University of<br>
*>       California at Berkeley, USA \n<br>
*>     Osni Marques, LBNL/NERSC, USA \n<br>
*<br>
*  =====================================================================<br>
*/
	public void clals0_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER NRHS,float[] B,INTEGER LDB,float[] BX,INTEGER LDBX,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,float[] GIVNUM,INTEGER LDGNUM,float[] POLES,float[] DIFL,float[] DIFR,float[] Z,INTEGER K,REAL C,REAL S,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLALSA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clalsa.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clalsa.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clalsa.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U,<br>
*                          LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR,<br>
*                          GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS,<br>
*      $                   SMLSIZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),<br>
*      $                   K( * ), PERM( LDGCOL, * )<br>
*       REAL               C( * ), DIFL( LDU, * ), DIFR( LDU, * ),<br>
*      $                   GIVNUM( LDU, * ), POLES( LDU, * ), RWORK( * ),<br>
*      $                   S( * ), U( LDU, * ), VT( LDU, * ), Z( LDU, * )<br>
*       COMPLEX            B( LDB, * ), BX( LDBX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLALSA is an itermediate step in solving the least squares problem<br>
*> by computing the SVD of the coefficient matrix in compact form (The<br>
*> singular vectors are computed as products of simple orthorgonal<br>
*> matrices.).<br>
*><br>
*> If ICOMPQ = 0, CLALSA applies the inverse of the left singular vector<br>
*> matrix of an upper bidiagonal matrix to the right hand side; and if<br>
*> ICOMPQ = 1, CLALSA applies the right singular vector matrix to the<br>
*> right hand side. The singular vector matrices were generated in<br>
*> compact form by CLALSA.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ICOMPQ<br>
*> \verbatim<br>
*>          ICOMPQ is INTEGER<br>
*>         Specifies whether the left or the right singular vector<br>
*>         matrix is involved.<br>
*>         = 0: Left singular vector matrix<br>
*>         = 1: Right singular vector matrix<br>
*> \endverbatim<br>
*><br>
*> \param[in] SMLSIZ<br>
*> \verbatim<br>
*>          SMLSIZ is INTEGER<br>
*>         The maximum size of the subproblems at the bottom of the<br>
*>         computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The row and column dimensions of the upper bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>         The number of columns of B and BX. NRHS must be at least 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension ( LDB, NRHS )<br>
*>         On input, B contains the right hand sides of the least<br>
*>         squares problem in rows 1 through M.<br>
*>         On output, B contains the solution X in rows 1 through N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>         The leading dimension of B in the calling subprogram.<br>
*>         LDB must be at least max(1,MAX( M, N ) ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BX<br>
*> \verbatim<br>
*>          BX is COMPLEX array, dimension ( LDBX, NRHS )<br>
*>         On exit, the result of applying the left or right singular<br>
*>         vector matrix to B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDBX<br>
*> \verbatim<br>
*>          LDBX is INTEGER<br>
*>         The leading dimension of BX.<br>
*> \endverbatim<br>
*><br>
*> \param[in] U<br>
*> \verbatim<br>
*>          U is REAL array, dimension ( LDU, SMLSIZ ).<br>
*>         On entry, U contains the left singular vector matrices of all<br>
*>         subproblems at the bottom level.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER, LDU = > N.<br>
*>         The leading dimension of arrays U, VT, DIFL, DIFR,<br>
*>         POLES, GIVNUM, and Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VT<br>
*> \verbatim<br>
*>          VT is REAL array, dimension ( LDU, SMLSIZ+1 ).<br>
*>         On entry, VT**H contains the right singular vector matrices of<br>
*>         all subproblems at the bottom level.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER array, dimension ( N ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIFL<br>
*> \verbatim<br>
*>          DIFL is REAL array, dimension ( LDU, NLVL ).<br>
*>         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIFR<br>
*> \verbatim<br>
*>          DIFR is REAL array, dimension ( LDU, 2 * NLVL ).<br>
*>         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record<br>
*>         distances between singular values on the I-th level and<br>
*>         singular values on the (I -1)-th level, and DIFR(*, 2 * I)<br>
*>         record the normalizing factors of the right singular vectors<br>
*>         matrices of subproblems on I-th level.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( LDU, NLVL ).<br>
*>         On entry, Z(1, I) contains the components of the deflation-<br>
*>         adjusted updating row vector for subproblems on the I-th<br>
*>         level.<br>
*> \endverbatim<br>
*><br>
*> \param[in] POLES<br>
*> \verbatim<br>
*>          POLES is REAL array, dimension ( LDU, 2 * NLVL ).<br>
*>         On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old<br>
*>         singular values involved in the secular equations on the I-th<br>
*>         level.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVPTR<br>
*> \verbatim<br>
*>          GIVPTR is INTEGER array, dimension ( N ).<br>
*>         On entry, GIVPTR( I ) records the number of Givens<br>
*>         rotations performed on the I-th problem on the computation<br>
*>         tree.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVCOL<br>
*> \verbatim<br>
*>          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 * NLVL ).<br>
*>         On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the<br>
*>         locations of Givens rotations performed on the I-th level on<br>
*>         the computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDGCOL<br>
*> \verbatim<br>
*>          LDGCOL is INTEGER, LDGCOL = > N.<br>
*>         The leading dimension of arrays GIVCOL and PERM.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PERM<br>
*> \verbatim<br>
*>          PERM is INTEGER array, dimension ( LDGCOL, NLVL ).<br>
*>         On entry, PERM(*, I) records permutations done on the I-th<br>
*>         level of the computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GIVNUM<br>
*> \verbatim<br>
*>          GIVNUM is REAL array, dimension ( LDU, 2 * NLVL ).<br>
*>         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S-<br>
*>         values of Givens rotations performed on the I-th level on the<br>
*>         computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension ( N ).<br>
*>         On entry, if the I-th subproblem is not square,<br>
*>         C( I ) contains the C-value of a Givens rotation related to<br>
*>         the right null space of the I-th subproblem.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension ( N ).<br>
*>         On entry, if the I-th subproblem is not square,<br>
*>         S( I ) contains the S-value of a Givens rotation related to<br>
*>         the right null space of the I-th subproblem.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension at least<br>
*>         MAX( (SMLSZ+1)*NRHS*3, N*(1+NRHS) + 2*NRHS ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array.<br>
*>         The dimension must be at least 3 * N<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Ren-Cang Li, Computer Science Division, University of<br>
*>       California at Berkeley, USA \n<br>
*>     Osni Marques, LBNL/NERSC, USA \n<br>
*<br>
*  =====================================================================<br>
*/
	public void clalsa_(INTEGER ICOMPQ,INTEGER SMLSIZ,INTEGER N,INTEGER NRHS,float[] B,INTEGER LDB,float[] BX,INTEGER LDBX,float[] U,INTEGER LDU,float[] VT,int[] K,float[] DIFL,float[] DIFR,float[] Z,float[] POLES,int[] GIVPTR,int[] GIVCOL,INTEGER LDGCOL,int[] PERM,float[] GIVNUM,float[] C,float[] S,float[] RWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b CLALSD uses the singular value decomposition of A to solve the least squares problem.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLALSD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clalsd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clalsd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clalsd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND,<br>
*                          RANK, WORK, RWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               D( * ), E( * ), RWORK( * )<br>
*       COMPLEX            B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLALSD uses the singular value decomposition of A to solve the least<br>
*> squares problem of finding X to minimize the Euclidean norm of each<br>
*> column of A*X-B, where A is N-by-N upper bidiagonal, and X and B<br>
*> are N-by-NRHS. The solution X overwrites B.<br>
*><br>
*> The singular values of A smaller than RCOND times the largest<br>
*> singular value are treated as zero in solving the least squares<br>
*> problem; in this case a minimum norm solution is returned.<br>
*> The actual singular values are returned in D in ascending order.<br>
*><br>
*> This code makes very mild assumptions about floating point<br>
*> arithmetic. It will work on machines with a guard digit in<br>
*> add/subtract, or on those binary machines without guard digits<br>
*> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.<br>
*> It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>         = 'U': D and E define an upper bidiagonal matrix.<br>
*>         = 'L': D and E define a  lower bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SMLSIZ<br>
*> \verbatim<br>
*>          SMLSIZ is INTEGER<br>
*>         The maximum size of the subproblems at the bottom of the<br>
*>         computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The dimension of the  bidiagonal matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>         The number of columns of B. NRHS must be at least 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry D contains the main diagonal of the bidiagonal<br>
*>         matrix. On exit, if INFO = 0, D contains its singular values.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>         Contains the super-diagonal entries of the bidiagonal matrix.<br>
*>         On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>         On input, B contains the right hand sides of the least<br>
*>         squares problem. On output, B contains the solution X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>         The leading dimension of B in the calling subprogram.<br>
*>         LDB must be at least max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>         The singular values of A less than or equal to RCOND times<br>
*>         the largest singular value are treated as zero in solving<br>
*>         the least squares problem. If RCOND is negative,<br>
*>         machine precision is used instead.<br>
*>         For example, if diag(S)*X=B were the least squares problem,<br>
*>         where diag(S) is a diagonal matrix of singular values, the<br>
*>         solution would be X(i) = B(i) / S(i) if S(i) is greater than<br>
*>         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to<br>
*>         RCOND*max(S).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RANK<br>
*> \verbatim<br>
*>          RANK is INTEGER<br>
*>         The number of singular values of A greater than RCOND times<br>
*>         the largest singular value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (N * NRHS).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension at least<br>
*>         (9*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS +<br>
*>         MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ),<br>
*>         where<br>
*>         NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (3*N*NLVL + 11*N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>         = 0:  successful exit.<br>
*>         < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>         > 0:  The algorithm failed to compute a singular value while<br>
*>               working on the submatrix lying in rows and columns<br>
*>               INFO/(N+1) through MOD(INFO,N+1).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Ren-Cang Li, Computer Science Division, University of<br>
*>       California at Berkeley, USA \n<br>
*>     Osni Marques, LBNL/NERSC, USA \n<br>
*<br>
*  =====================================================================<br>
*/
	public void clalsd_(CHARACTER UPLO,INTEGER SMLSIZ,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,REAL RCOND,INTEGER RANK,float[] WORK,float[] RWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b CLANGB returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of general band matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANGB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clangb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clangb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clangb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANGB( NORM, N, KL, KU, AB, LDAB,<br>
*                        WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            KL, KU, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANGB  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the element of  largest absolute value  of an<br>
*> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.<br>
*> \endverbatim<br>
*><br>
*> \return CLANGB<br>
*> \verbatim<br>
*><br>
*>    CLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANGB as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANGB is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>          The number of sub-diagonals of the matrix A.  KL >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>          The number of super-diagonals of the matrix A.  KU >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th<br>
*>          column of A is stored in the j-th column of the array AB as<br>
*>          follows:<br>
*>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I'; otherwise, WORK is not<br>
*>          referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGBauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clangb_(CHARACTER NORM,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] WORK);
/**
*> \brief \b CLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANGE + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clange.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clange.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clange.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANGE( NORM, M, N, A, LDA, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANGE  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> complex matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return CLANGE<br>
*> \verbatim<br>
*><br>
*>    CLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANGE as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.  When M = 0,<br>
*>          CLANGE is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.  When N = 0,<br>
*>          CLANGE is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(M,1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not<br>
*>          referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGEauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clange_(CHARACTER NORM,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b CLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general tridiagonal matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANGT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clangt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clangt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clangt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANGT( NORM, N, DL, D, DU )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            D( * ), DL( * ), DU( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANGT  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> complex tridiagonal matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return CLANGT<br>
*> \verbatim<br>
*><br>
*>    CLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANGT as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANGT is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DL<br>
*> \verbatim<br>
*>          DL is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) sub-diagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is COMPLEX array, dimension (N)<br>
*>          The diagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) super-diagonal elements of A.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clangt_(CHARACTER NORM,INTEGER N,float[] DL,float[] D,float[] DU);
/**
*> \brief \b CLANHB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a Hermitian band matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANHB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANHB( NORM, UPLO, N, K, AB, LDAB,<br>
*                        WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            K, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANHB  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the element of  largest absolute value  of an<br>
*> n by n hermitian band matrix A,  with k super-diagonals.<br>
*> \endverbatim<br>
*><br>
*> \return CLANHB<br>
*> \verbatim<br>
*><br>
*>    CLANHB = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANHB as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          band matrix A is supplied.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANHB is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of super-diagonals or sub-diagonals of the<br>
*>          band matrix A.  K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The upper or lower triangle of the hermitian band matrix A,<br>
*>          stored in the first K+1 rows of AB.  The j-th column of A is<br>
*>          stored in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).<br>
*>          Note that the imaginary parts of the diagonal elements need<br>
*>          not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= K+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,<br>
*>          WORK is not referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clanhb_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,INTEGER K,float[] AB,INTEGER LDAB,float[] WORK);
/**
*> \brief \b CLANHE returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a complex Hermitian matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANHE + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhe.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhe.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhe.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANHE( NORM, UPLO, N, A, LDA, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANHE  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> complex hermitian matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return CLANHE<br>
*> \verbatim<br>
*><br>
*>    CLANHE = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANHE as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          hermitian matrix A is to be referenced.<br>
*>          = 'U':  Upper triangular part of A is referenced<br>
*>          = 'L':  Lower triangular part of A is referenced<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANHE is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The hermitian matrix A.  If UPLO = 'U', the leading n by n<br>
*>          upper triangular part of A contains the upper triangular part<br>
*>          of the matrix A, and the strictly lower triangular part of A<br>
*>          is not referenced.  If UPLO = 'L', the leading n by n lower<br>
*>          triangular part of A contains the lower triangular part of<br>
*>          the matrix A, and the strictly upper triangular part of A is<br>
*>          not referenced. Note that the imaginary parts of the diagonal<br>
*>          elements need not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(N,1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,<br>
*>          WORK is not referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexHEauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clanhe_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b CLANHF returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a Hermitian matrix in RFP format.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANHF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLANHF( NORM, TRANSR, UPLO, N, A, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, TRANSR, UPLO<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( 0: * )<br>
*       COMPLEX            A( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANHF  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> complex Hermitian matrix A in RFP format.<br>
*> \endverbatim<br>
*><br>
*> \return CLANHF<br>
*> \verbatim<br>
*><br>
*>    CLANHF = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER<br>
*>            Specifies the value to be returned in CLANHF as described<br>
*>            above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER<br>
*>            Specifies whether the RFP format of A is normal or<br>
*>            conjugate-transposed format.<br>
*>            = 'N':  RFP format is Normal<br>
*>            = 'C':  RFP format is Conjugate-transposed<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER<br>
*>            On entry, UPLO specifies whether the RFP matrix A came from<br>
*>            an upper or lower triangular matrix as follows:<br>
*><br>
*>            UPLO = 'U' or 'u' RFP A came from an upper triangular<br>
*>            matrix<br>
*><br>
*>            UPLO = 'L' or 'l' RFP A came from a  lower triangular<br>
*>            matrix<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>            The order of the matrix A.  N >= 0.  When N = 0, CLANHF is<br>
*>            set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension ( N*(N+1)/2 );<br>
*>            On entry, the matrix A in RFP Format.<br>
*>            RFP Format is described by TRANSR, UPLO and N as follows:<br>
*>            If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;<br>
*>            K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If<br>
*>            TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A<br>
*>            as defined when TRANSR = 'N'. The contents of RFP A are<br>
*>            defined by UPLO as follows: If UPLO = 'U' the RFP A<br>
*>            contains the ( N*(N+1)/2 ) elements of upper packed A<br>
*>            either in normal or conjugate-transpose Format. If<br>
*>            UPLO = 'L' the RFP A contains the ( N*(N+1) /2 ) elements<br>
*>            of lower packed A either in normal or conjugate-transpose<br>
*>            Format. The LDA of RFP A is (N+1)/2 when TRANSR = 'C'. When<br>
*>            TRANSR is 'N' the LDA is N+1 when N is even and is N when<br>
*>            is odd. See the Note below for more details.<br>
*>            Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (LWORK),<br>
*>            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,<br>
*>            WORK is not referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2015<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  We first consider Standard Packed Format when N is even.<br>
*>  We give an example where N = 6.<br>
*><br>
*>      AP is Upper             AP is Lower<br>
*><br>
*>   00 01 02 03 04 05       00<br>
*>      11 12 13 14 15       10 11<br>
*>         22 23 24 25       20 21 22<br>
*>            33 34 35       30 31 32 33<br>
*>               44 45       40 41 42 43 44<br>
*>                  55       50 51 52 53 54 55<br>
*><br>
*><br>
*>  Let TRANSR = 'N'. RFP holds AP as follows:<br>
*>  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last<br>
*>  three columns of AP upper. The lower triangle A(4:6,0:2) consists of<br>
*>  conjugate-transpose of the first three columns of AP upper.<br>
*>  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first<br>
*>  three columns of AP lower. The upper triangle A(0:2,0:2) consists of<br>
*>  conjugate-transpose of the last three columns of AP lower.<br>
*>  To denote conjugate we place -- above the element. This covers the<br>
*>  case N even and TRANSR = 'N'.<br>
*><br>
*>         RFP A                   RFP A<br>
*><br>
*>                                -- -- --<br>
*>        03 04 05                33 43 53<br>
*>                                   -- --<br>
*>        13 14 15                00 44 54<br>
*>                                      --<br>
*>        23 24 25                10 11 55<br>
*><br>
*>        33 34 35                20 21 22<br>
*>        --<br>
*>        00 44 45                30 31 32<br>
*>        -- --<br>
*>        01 11 55                40 41 42<br>
*>        -- -- --<br>
*>        02 12 22                50 51 52<br>
*><br>
*>  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-<br>
*>  transpose of RFP A above. One therefore gets:<br>
*><br>
*><br>
*>           RFP A                   RFP A<br>
*><br>
*>     -- -- -- --                -- -- -- -- -- --<br>
*>     03 13 23 33 00 01 02    33 00 10 20 30 40 50<br>
*>     -- -- -- -- --                -- -- -- -- --<br>
*>     04 14 24 34 44 11 12    43 44 11 21 31 41 51<br>
*>     -- -- -- -- -- --                -- -- -- --<br>
*>     05 15 25 35 45 55 22    53 54 55 22 32 42 52<br>
*><br>
*><br>
*>  We next  consider Standard Packed Format when N is odd.<br>
*>  We give an example where N = 5.<br>
*><br>
*>     AP is Upper                 AP is Lower<br>
*><br>
*>   00 01 02 03 04              00<br>
*>      11 12 13 14              10 11<br>
*>         22 23 24              20 21 22<br>
*>            33 34              30 31 32 33<br>
*>               44              40 41 42 43 44<br>
*><br>
*><br>
*>  Let TRANSR = 'N'. RFP holds AP as follows:<br>
*>  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last<br>
*>  three columns of AP upper. The lower triangle A(3:4,0:1) consists of<br>
*>  conjugate-transpose of the first two   columns of AP upper.<br>
*>  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first<br>
*>  three columns of AP lower. The upper triangle A(0:1,1:2) consists of<br>
*>  conjugate-transpose of the last two   columns of AP lower.<br>
*>  To denote conjugate we place -- above the element. This covers the<br>
*>  case N odd  and TRANSR = 'N'.<br>
*><br>
*>         RFP A                   RFP A<br>
*><br>
*>                                   -- --<br>
*>        02 03 04                00 33 43<br>
*>                                      --<br>
*>        12 13 14                10 11 44<br>
*><br>
*>        22 23 24                20 21 22<br>
*>        --<br>
*>        00 33 34                30 31 32<br>
*>        -- --<br>
*>        01 11 44                40 41 42<br>
*><br>
*>  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-<br>
*>  transpose of RFP A above. One therefore gets:<br>
*><br>
*><br>
*>           RFP A                   RFP A<br>
*><br>
*>     -- -- --                   -- -- -- -- -- --<br>
*>     02 12 22 00 01             00 10 20 30 40 50<br>
*>     -- -- -- --                   -- -- -- -- --<br>
*>     03 13 23 33 11             33 11 21 31 41 51<br>
*>     -- -- -- -- --                   -- -- -- --<br>
*>     04 14 24 34 44             43 44 22 32 42 52<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public float clanhf_(CHARACTER NORM,CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,float[] WORK);
/**
*> \brief \b CLANHP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a complex Hermitian matrix supplied in packed form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANHP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANHP( NORM, UPLO, N, AP, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANHP  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> complex hermitian matrix A,  supplied in packed form.<br>
*> \endverbatim<br>
*><br>
*> \return CLANHP<br>
*> \verbatim<br>
*><br>
*>    CLANHP = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANHP as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          hermitian matrix A is supplied.<br>
*>          = 'U':  Upper triangular part of A is supplied<br>
*>          = 'L':  Lower triangular part of A is supplied<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANHP is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangle of the hermitian matrix A, packed<br>
*>          columnwise in a linear array.  The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          Note that the  imaginary parts of the diagonal elements need<br>
*>          not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,<br>
*>          WORK is not referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clanhp_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,float[] AP,float[] WORK);
/**
*> \brief \b CLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANHS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANHS( NORM, N, A, LDA, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANHS  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> Hessenberg matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return CLANHS<br>
*> \verbatim<br>
*><br>
*>    CLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANHS as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANHS is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The n by n upper Hessenberg matrix A; the part of A below the<br>
*>          first sub-diagonal is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(N,1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I'; otherwise, WORK is not<br>
*>          referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clanhs_(CHARACTER NORM,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b CLANHT returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a complex Hermitian tridiagonal matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANHT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanht.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanht.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanht.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANHT( NORM, N, D, E )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * )<br>
*       COMPLEX            E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANHT  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> complex Hermitian tridiagonal matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return CLANHT<br>
*> \verbatim<br>
*><br>
*>    CLANHT = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANHT as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANHT is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The diagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) sub-diagonal or super-diagonal elements of A.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clanht_(CHARACTER NORM,INTEGER N,float[] D,float[] E);
/**
*> \brief \b CLANSB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric band matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANSB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANSB( NORM, UPLO, N, K, AB, LDAB,<br>
*                        WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            K, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANSB  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the element of  largest absolute value  of an<br>
*> n by n symmetric band matrix A,  with k super-diagonals.<br>
*> \endverbatim<br>
*><br>
*> \return CLANSB<br>
*> \verbatim<br>
*><br>
*>    CLANSB = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANSB as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          band matrix A is supplied.<br>
*>          = 'U':  Upper triangular part is supplied<br>
*>          = 'L':  Lower triangular part is supplied<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANSB is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of super-diagonals or sub-diagonals of the<br>
*>          band matrix A.  K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The upper or lower triangle of the symmetric band matrix A,<br>
*>          stored in the first K+1 rows of AB.  The j-th column of A is<br>
*>          stored in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= K+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,<br>
*>          WORK is not referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clansb_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,INTEGER K,float[] AB,INTEGER LDAB,float[] WORK);
/**
*> \brief \b CLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric matrix supplied in packed form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANSP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANSP( NORM, UPLO, N, AP, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANSP  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> complex symmetric matrix A,  supplied in packed form.<br>
*> \endverbatim<br>
*><br>
*> \return CLANSP<br>
*> \verbatim<br>
*><br>
*>    CLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANSP as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is supplied.<br>
*>          = 'U':  Upper triangular part of A is supplied<br>
*>          = 'L':  Lower triangular part of A is supplied<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANSP is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangle of the symmetric matrix A, packed<br>
*>          columnwise in a linear array.  The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,<br>
*>          WORK is not referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clansp_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,float[] AP,float[] WORK);
/**
*> \brief \b CLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a complex symmetric matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANSY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clansy.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clansy.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clansy.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANSY( NORM, UPLO, N, A, LDA, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANSY  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> complex symmetric matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return CLANSY<br>
*> \verbatim<br>
*><br>
*>    CLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANSY as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is to be referenced.<br>
*>          = 'U':  Upper triangular part of A is referenced<br>
*>          = 'L':  Lower triangular part of A is referenced<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANSY is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The symmetric matrix A.  If UPLO = 'U', the leading n by n<br>
*>          upper triangular part of A contains the upper triangular part<br>
*>          of the matrix A, and the strictly lower triangular part of A<br>
*>          is not referenced.  If UPLO = 'L', the leading n by n lower<br>
*>          triangular part of A contains the lower triangular part of<br>
*>          the matrix A, and the strictly upper triangular part of A is<br>
*>          not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(N,1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,<br>
*>          WORK is not referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2015<br>
*<br>
*> \ingroup complexSYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clansy_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b CLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a triangular band matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANTB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clantb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clantb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clantb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANTB( NORM, UPLO, DIAG, N, K, AB,<br>
*                        LDAB, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            K, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANTB  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the element of  largest absolute value  of an<br>
*> n by n triangular band matrix A,  with ( k + 1 ) diagonals.<br>
*> \endverbatim<br>
*><br>
*> \return CLANTB<br>
*> \verbatim<br>
*><br>
*>    CLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANTB as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the matrix A is upper or lower triangular.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          Specifies whether or not the matrix A is unit triangular.<br>
*>          = 'N':  Non-unit triangular<br>
*>          = 'U':  Unit triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANTB is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of super-diagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of sub-diagonals of the matrix A if UPLO = 'L'.<br>
*>          K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The upper or lower triangular band matrix A, stored in the<br>
*>          first k+1 rows of AB.  The j-th column of A is stored<br>
*>          in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).<br>
*>          Note that when DIAG = 'U', the elements of the array AB<br>
*>          corresponding to the diagonal elements of the matrix A are<br>
*>          not referenced, but are assumed to be one.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= K+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I'; otherwise, WORK is not<br>
*>          referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clantb_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,INTEGER K,float[] AB,INTEGER LDAB,float[] WORK);
/**
*> \brief \b CLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a triangular matrix supplied in packed form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANTP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clantp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clantp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clantp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANTP( NORM, UPLO, DIAG, N, AP, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANTP  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> triangular matrix A, supplied in packed form.<br>
*> \endverbatim<br>
*><br>
*> \return CLANTP<br>
*> \verbatim<br>
*><br>
*>    CLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANTP as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the matrix A is upper or lower triangular.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          Specifies whether or not the matrix A is unit triangular.<br>
*>          = 'N':  Non-unit triangular<br>
*>          = 'U':  Unit triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, CLANTP is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangular matrix A, packed columnwise in<br>
*>          a linear array.  The j-th column of A is stored in the array<br>
*>          AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          Note that when DIAG = 'U', the elements of the array AP<br>
*>          corresponding to the diagonal elements of the matrix A are<br>
*>          not referenced, but are assumed to be one.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= N when NORM = 'I'; otherwise, WORK is not<br>
*>          referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clantp_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] AP,float[] WORK);
/**
*> \brief \b CLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a trapezoidal or triangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLANTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clantr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clantr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clantr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION CLANTR( NORM, UPLO, DIAG, M, N, A, LDA,<br>
*                        WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               WORK( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLANTR  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> trapezoidal or triangular matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return CLANTR<br>
*> \verbatim<br>
*><br>
*>    CLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
*>             (<br>
*>             ( norm1(A),         NORM = '1', 'O' or 'o'<br>
*>             (<br>
*>             ( normI(A),         NORM = 'I' or 'i'<br>
*>             (<br>
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'<br>
*><br>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),<br>
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and<br>
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of<br>
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in CLANTR as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the matrix A is upper or lower trapezoidal.<br>
*>          = 'U':  Upper trapezoidal<br>
*>          = 'L':  Lower trapezoidal<br>
*>          Note that A is triangular instead of trapezoidal if M = N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          Specifies whether or not the matrix A has unit diagonal.<br>
*>          = 'N':  Non-unit diagonal<br>
*>          = 'U':  Unit diagonal<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0, and if<br>
*>          UPLO = 'U', M <= N.  When M = 0, CLANTR is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0, and if<br>
*>          UPLO = 'L', N <= M.  When N = 0, CLANTR is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The trapezoidal matrix A (A is triangular if M = N).<br>
*>          If UPLO = 'U', the leading m by n upper trapezoidal part of<br>
*>          the array A contains the upper trapezoidal matrix, and the<br>
*>          strictly lower triangular part of A is not referenced.<br>
*>          If UPLO = 'L', the leading m by n lower trapezoidal part of<br>
*>          the array A contains the lower trapezoidal matrix, and the<br>
*>          strictly upper triangular part of A is not referenced.  Note<br>
*>          that when DIAG = 'U', the diagonal elements of A are not<br>
*>          referenced and are assumed to be one.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(M,1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)),<br>
*>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not<br>
*>          referenced.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float clantr_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b CLAPLL measures the linear dependence of two vectors.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAPLL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clapll.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clapll.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clapll.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAPLL( N, X, INCX, Y, INCY, SSMIN )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, INCY, N<br>
*       REAL               SSMIN<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Given two column vectors X and Y, let<br>
*><br>
*>                      A = ( X Y ).<br>
*><br>
*> The subroutine first computes the QR factorization of A = Q*R,<br>
*> and then computes the SVD of the 2-by-2 upper triangular matrix R.<br>
*> The smaller singular value of R is returned in SSMIN, which is used<br>
*> as the measurement of the linear dependency of the vectors X and Y.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The length of the vectors X and Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (1+(N-1)*INCX)<br>
*>          On entry, X contains the N-vector X.<br>
*>          On exit, X is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between successive elements of X. INCX > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (1+(N-1)*INCY)<br>
*>          On entry, Y contains the N-vector Y.<br>
*>          On exit, Y is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>          The increment between successive elements of Y. INCY > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SSMIN<br>
*> \verbatim<br>
*>          SSMIN is REAL<br>
*>          The smallest singular value of the N-by-2 matrix A = ( X Y ).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clapll_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,REAL SSMIN);
/**
*> \brief \b CLAPMR rearranges rows of a matrix as specified by a permutation vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAPMR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clapmr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clapmr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clapmr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAPMR( FORWRD, M, N, X, LDX, K )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            FORWRD<br>
*       INTEGER            LDX, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            K( * )<br>
*       COMPLEX            X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAPMR rearranges the rows of the M by N matrix X as specified<br>
*> by the permutation K(1),K(2),...,K(M) of the integers 1,...,M.<br>
*> If FORWRD = .TRUE.,  forward permutation:<br>
*><br>
*>      X(K(I),*) is moved X(I,*) for I = 1,2,...,M.<br>
*><br>
*> If FORWRD = .FALSE., backward permutation:<br>
*><br>
*>      X(I,*) is moved to X(K(I),*) for I = 1,2,...,M.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] FORWRD<br>
*> \verbatim<br>
*>          FORWRD is LOGICAL<br>
*>          = .TRUE., forward permutation<br>
*>          = .FALSE., backward permutation<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix X. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix X. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (LDX,N)<br>
*>          On entry, the M by N matrix X.<br>
*>          On exit, X contains the permuted matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X, LDX >= MAX(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] K<br>
*> \verbatim<br>
*>          K is INTEGER array, dimension (M)<br>
*>          On entry, K contains the permutation vector. K is used as<br>
*>          internal workspace, but reset to its original value on<br>
*>          output.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clapmr_(LOGICAL FORWRD,INTEGER M,INTEGER N,float[] X,INTEGER LDX,int[] K);
/**
*> \brief \b CLAPMT performs a forward or backward permutation of the columns of a matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAPMT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clapmt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clapmt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clapmt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAPMT( FORWRD, M, N, X, LDX, K )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            FORWRD<br>
*       INTEGER            LDX, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            K( * )<br>
*       COMPLEX            X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAPMT rearranges the columns of the M by N matrix X as specified<br>
*> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.<br>
*> If FORWRD = .TRUE.,  forward permutation:<br>
*><br>
*>      X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.<br>
*><br>
*> If FORWRD = .FALSE., backward permutation:<br>
*><br>
*>      X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] FORWRD<br>
*> \verbatim<br>
*>          FORWRD is LOGICAL<br>
*>          = .TRUE., forward permutation<br>
*>          = .FALSE., backward permutation<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix X. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix X. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (LDX,N)<br>
*>          On entry, the M by N matrix X.<br>
*>          On exit, X contains the permuted matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X, LDX >= MAX(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] K<br>
*> \verbatim<br>
*>          K is INTEGER array, dimension (N)<br>
*>          On entry, K contains the permutation vector. K is used as<br>
*>          internal workspace, but reset to its original value on<br>
*>          output.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clapmt_(LOGICAL FORWRD,INTEGER M,INTEGER N,float[] X,INTEGER LDX,int[] K);
/**
*> \brief \b CLAQGB scales a general band matrix, using row and column scaling factors computed by sgbequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQGB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqgb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqgb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqgb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,<br>
*                          AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED<br>
*       INTEGER            KL, KU, LDAB, M, N<br>
*       REAL               AMAX, COLCND, ROWCND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * ), R( * )<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQGB equilibrates a general M by N band matrix A with KL<br>
*> subdiagonals and KU superdiagonals using the row and scaling factors<br>
*> in the vectors R and C.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>          The number of subdiagonals within the band of A.  KL >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>          The number of superdiagonals within the band of A.  KU >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.<br>
*>          The j-th column of A is stored in the j-th column of the<br>
*>          array AB as follows:<br>
*>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)<br>
*><br>
*>          On exit, the equilibrated matrix, in the same storage format<br>
*>          as A.  See EQUED for the form of the equilibrated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDA >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] R<br>
*> \verbatim<br>
*>          R is REAL array, dimension (M)<br>
*>          The row scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>          The column scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ROWCND<br>
*> \verbatim<br>
*>          ROWCND is REAL<br>
*>          Ratio of the smallest R(i) to the largest R(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] COLCND<br>
*> \verbatim<br>
*>          COLCND is REAL<br>
*>          Ratio of the smallest C(i) to the largest C(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix entry.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies the form of equilibration that was done.<br>
*>          = 'N':  No equilibration<br>
*>          = 'R':  Row equilibration, i.e., A has been premultiplied by<br>
*>                  diag(R).<br>
*>          = 'C':  Column equilibration, i.e., A has been postmultiplied<br>
*>                  by diag(C).<br>
*>          = 'B':  Both row and column equilibration, i.e., A has been<br>
*>                  replaced by diag(R) * A * diag(C).<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  THRESH is a threshold value used to decide if row or column scaling<br>
*>  should be done based on the ratio of the row or column scaling<br>
*>  factors.  If ROWCND < THRESH, row scaling is done, and if<br>
*>  COLCND < THRESH, column scaling is done.<br>
*><br>
*>  LARGE and SMALL are threshold values used to decide if row scaling<br>
*>  should be done based on the absolute size of the largest matrix<br>
*>  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGBauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claqgb_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b CLAQGE scales a general rectangular matrix, using row and column scaling factors computed by sgeequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQGE + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqge.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqge.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqge.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,<br>
*                          EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED<br>
*       INTEGER            LDA, M, N<br>
*       REAL               AMAX, COLCND, ROWCND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * ), R( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQGE equilibrates a general M by N matrix A using the row and<br>
*> column scaling factors in the vectors R and C.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the M by N matrix A.<br>
*>          On exit, the equilibrated matrix.  See EQUED for the form of<br>
*>          the equilibrated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(M,1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] R<br>
*> \verbatim<br>
*>          R is REAL array, dimension (M)<br>
*>          The row scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>          The column scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ROWCND<br>
*> \verbatim<br>
*>          ROWCND is REAL<br>
*>          Ratio of the smallest R(i) to the largest R(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] COLCND<br>
*> \verbatim<br>
*>          COLCND is REAL<br>
*>          Ratio of the smallest C(i) to the largest C(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix entry.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies the form of equilibration that was done.<br>
*>          = 'N':  No equilibration<br>
*>          = 'R':  Row equilibration, i.e., A has been premultiplied by<br>
*>                  diag(R).<br>
*>          = 'C':  Column equilibration, i.e., A has been postmultiplied<br>
*>                  by diag(C).<br>
*>          = 'B':  Both row and column equilibration, i.e., A has been<br>
*>                  replaced by diag(R) * A * diag(C).<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  THRESH is a threshold value used to decide if row or column scaling<br>
*>  should be done based on the ratio of the row or column scaling<br>
*>  factors.  If ROWCND < THRESH, row scaling is done, and if<br>
*>  COLCND < THRESH, column scaling is done.<br>
*><br>
*>  LARGE and SMALL are threshold values used to decide if row scaling<br>
*>  should be done based on the absolute size of the largest matrix<br>
*>  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGEauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claqge_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b CLAQHB scales a Hermitian band matrix, using scaling factors computed by cpbequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQHB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqhb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqhb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqhb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQHB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            KD, LDAB, N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               S( * )<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQHB equilibrates an Hermitian band matrix A using the scaling<br>
*> factors in the vector S.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is stored.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KD<br>
*> \verbatim<br>
*>          KD is INTEGER<br>
*>          The number of super-diagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*><br>
*>          On exit, if INFO = 0, the triangular factor U or L from the<br>
*>          Cholesky factorization A = U**H *U or A = L*L**H of the band<br>
*>          matrix A, in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          Ratio of the smallest S(i) to the largest S(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix entry.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies whether or not equilibration was done.<br>
*>          = 'N':  No equilibration.<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  THRESH is a threshold value used to decide if scaling should be done<br>
*>  based on the ratio of the scaling factors.  If SCOND < THRESH,<br>
*>  scaling is done.<br>
*><br>
*>  LARGE and SMALL are threshold values used to decide if scaling should<br>
*>  be done based on the absolute size of the largest matrix element.<br>
*>  If AMAX > LARGE or AMAX < SMALL, scaling is done.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claqhb_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b CLAQHE scales a Hermitian matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQHE + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqhe.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqhe.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqhe.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            LDA, N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               S( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQHE equilibrates a Hermitian matrix A using the scaling factors<br>
*> in the vector S.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          Hermitian matrix A is stored.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          n by n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n by n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if EQUED = 'Y', the equilibrated matrix:<br>
*>          diag(S) * A * diag(S).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(N,1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          Ratio of the smallest S(i) to the largest S(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix entry.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies whether or not equilibration was done.<br>
*>          = 'N':  No equilibration.<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  THRESH is a threshold value used to decide if scaling should be done<br>
*>  based on the ratio of the scaling factors.  If SCOND < THRESH,<br>
*>  scaling is done.<br>
*><br>
*>  LARGE and SMALL are threshold values used to decide if scaling should<br>
*>  be done based on the absolute size of the largest matrix element.<br>
*>  If AMAX > LARGE or AMAX < SMALL, scaling is done.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexHEauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claqhe_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b CLAQHP scales a Hermitian matrix stored in packed form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQHP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqhp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqhp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqhp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               S( * )<br>
*       COMPLEX            AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQHP equilibrates a Hermitian matrix A using the scaling factors<br>
*> in the vector S.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          Hermitian matrix A is stored.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the equilibrated matrix:  diag(S) * A * diag(S), in<br>
*>          the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          Ratio of the smallest S(i) to the largest S(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix entry.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies whether or not equilibration was done.<br>
*>          = 'N':  No equilibration.<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  THRESH is a threshold value used to decide if scaling should be done<br>
*>  based on the ratio of the scaling factors.  If SCOND < THRESH,<br>
*>  scaling is done.<br>
*><br>
*>  LARGE and SMALL are threshold values used to decide if scaling should<br>
*>  be done based on the absolute size of the largest matrix element.<br>
*>  If AMAX > LARGE or AMAX < SMALL, scaling is done.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claqhp_(CHARACTER UPLO,INTEGER N,float[] AP,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b CLAQP2 computes a QR factorization with column pivoting of the matrix block.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQP2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqp2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqp2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqp2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,<br>
*                          WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, M, N, OFFSET<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            JPVT( * )<br>
*       REAL               VN1( * ), VN2( * )<br>
*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQP2 computes a QR factorization with column pivoting of<br>
*> the block A(OFFSET+1:M,1:N).<br>
*> The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] OFFSET<br>
*> \verbatim<br>
*>          OFFSET is INTEGER<br>
*>          The number of rows of the matrix A that must be pivoted<br>
*>          but no factorized. OFFSET >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is <br>
*>          the triangular factor obtained; the elements in block<br>
*>          A(OFFSET+1:M,1:N) below the diagonal, together with the<br>
*>          array TAU, represent the orthogonal matrix Q as a product of<br>
*>          elementary reflectors. Block A(1:OFFSET,1:N) has been<br>
*>          accordingly pivoted, but no factorized.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] JPVT<br>
*> \verbatim<br>
*>          JPVT is INTEGER array, dimension (N)<br>
*>          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted<br>
*>          to the front of A*P (a leading column); if JPVT(i) = 0,<br>
*>          the i-th column of A is a free column.<br>
*>          On exit, if JPVT(i) = k, then the i-th column of A*P<br>
*>          was the k-th column of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VN1<br>
*> \verbatim<br>
*>          VN1 is REAL array, dimension (N)<br>
*>          The vector with the partial column norms.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VN2<br>
*> \verbatim<br>
*>          VN2 is REAL array, dimension (N)<br>
*>          The vector with the exact column norms.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (N)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain<br>
*>    X. Sun, Computer Science Dept., Duke University, USA<br>
*> \n<br>
*>  Partial column norm updating strategy modified on April 2011<br>
*>    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,<br>
*>    University of Zagreb, Croatia.<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*> LAPACK Working Note 176<br>
*<br>
*> \htmlonly<br>
*> <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a> <br>
*> \endhtmlonly <br>
*<br>
*  =====================================================================<br>
*/
	public void claqp2_(INTEGER M,INTEGER N,INTEGER OFFSET,float[] A,INTEGER LDA,int[] JPVT,float[] TAU,float[] VN1,float[] VN2,float[] WORK);
/**
*> \brief \b CLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by using BLAS level 3.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQPS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqps.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqps.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqps.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1,<br>
*                          VN2, AUXV, F, LDF )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            KB, LDA, LDF, M, N, NB, OFFSET<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            JPVT( * )<br>
*       REAL               VN1( * ), VN2( * )<br>
*       COMPLEX            A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQPS computes a step of QR factorization with column pivoting<br>
*> of a complex M-by-N matrix A by using Blas-3.  It tries to factorize<br>
*> NB columns from A starting from the row OFFSET+1, and updates all<br>
*> of the matrix with Blas-3 xGEMM.<br>
*><br>
*> In some cases, due to catastrophic cancellations, it cannot<br>
*> factorize NB columns.  Hence, the actual number of factorized<br>
*> columns is returned in KB.<br>
*><br>
*> Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A. N >= 0<br>
*> \endverbatim<br>
*><br>
*> \param[in] OFFSET<br>
*> \verbatim<br>
*>          OFFSET is INTEGER<br>
*>          The number of rows of A that have been factorized in<br>
*>          previous steps.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The number of columns to factorize.<br>
*> \endverbatim<br>
*><br>
*> \param[out] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of columns actually factorized.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, block A(OFFSET+1:M,1:KB) is the triangular<br>
*>          factor obtained and block A(1:OFFSET,1:N) has been<br>
*>          accordingly pivoted, but no factorized.<br>
*>          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has<br>
*>          been updated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] JPVT<br>
*> \verbatim<br>
*>          JPVT is INTEGER array, dimension (N)<br>
*>          JPVT(I) = K <==> Column K of the full matrix A has been<br>
*>          permuted into position I in AP.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX array, dimension (KB)<br>
*>          The scalar factors of the elementary reflectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VN1<br>
*> \verbatim<br>
*>          VN1 is REAL array, dimension (N)<br>
*>          The vector with the partial column norms.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VN2<br>
*> \verbatim<br>
*>          VN2 is REAL array, dimension (N)<br>
*>          The vector with the exact column norms.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AUXV<br>
*> \verbatim<br>
*>          AUXV is COMPLEX array, dimension (NB)<br>
*>          Auxiliar vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] F<br>
*> \verbatim<br>
*>          F is COMPLEX array, dimension (LDF,NB)<br>
*>          Matrix  F**H = L * Y**H * A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDF<br>
*> \verbatim<br>
*>          LDF is INTEGER<br>
*>          The leading dimension of the array F. LDF >= max(1,N).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain<br>
*>    X. Sun, Computer Science Dept., Duke University, USA<br>
*><br>
*> \n<br>
*>  Partial column norm updating strategy modified on April 2011<br>
*>    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,<br>
*>    University of Zagreb, Croatia.<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*> LAPACK Working Note 176<br>
*<br>
*> \htmlonly<br>
*> <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a> <br>
*> \endhtmlonly <br>
*<br>
*  =====================================================================<br>
*/
	public void claqps_(INTEGER M,INTEGER N,INTEGER OFFSET,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] JPVT,float[] TAU,float[] VN1,float[] VN2,float[] AUXV,float[] F,INTEGER LDF);
/**
*> \brief \b CLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQR0 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr0.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr0.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr0.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,<br>
*                          IHIZ, Z, LDZ, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLAQR0 computes the eigenvalues of a Hessenberg matrix H<br>
*>    and, optionally, the matrices T and Z from the Schur decomposition<br>
*>    H = Z T Z**H, where T is an upper triangular matrix (the<br>
*>    Schur form), and Z is the unitary matrix of Schur vectors.<br>
*><br>
*>    Optionally Z may be postmultiplied into an input unitary<br>
*>    matrix Q so that this routine can give the Schur factorization<br>
*>    of a matrix A which has been reduced to the Hessenberg form H<br>
*>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTT<br>
*> \verbatim<br>
*>          WANTT is LOGICAL<br>
*>          = .TRUE. : the full Schur form T is required;<br>
*>          = .FALSE.: only eigenvalues are required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          = .TRUE. : the matrix of Schur vectors Z is required;<br>
*>          = .FALSE.: Schur vectors are not required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           The order of the matrix H.  N .GE. 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILO<br>
*> \verbatim<br>
*>          ILO is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHI<br>
*> \verbatim<br>
*>          IHI is INTEGER<br>
*>           It is assumed that H is already upper triangular in rows<br>
*>           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,<br>
*>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a<br>
*>           previous call to CGEBAL, and then passed to CGEHRD when the<br>
*>           matrix output by CGEBAL is reduced to Hessenberg form.<br>
*>           Otherwise, ILO and IHI should be set to 1 and N,<br>
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.<br>
*>           If N = 0, then ILO = 1 and IHI = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is COMPLEX array, dimension (LDH,N)<br>
*>           On entry, the upper Hessenberg matrix H.<br>
*>           On exit, if INFO = 0 and WANTT is .TRUE., then H<br>
*>           contains the upper triangular matrix T from the Schur<br>
*>           decomposition (the Schur form). If INFO = 0 and WANT is<br>
*>           .FALSE., then the contents of H are unspecified on exit.<br>
*>           (The output value of H when INFO.GT.0 is given under the<br>
*>           description of INFO below.)<br>
*><br>
*>           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and<br>
*>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>           The leading dimension of the array H. LDH .GE. max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (N)<br>
*>           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored<br>
*>           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are<br>
*>           stored in the same order as on the diagonal of the Schur<br>
*>           form returned in H, with W(i) = H(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILOZ<br>
*> \verbatim<br>
*>          ILOZ is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHIZ<br>
*> \verbatim<br>
*>          IHIZ is INTEGER<br>
*>           Specify the rows of Z to which transformations must be<br>
*>           applied if WANTZ is .TRUE..<br>
*>           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (LDZ,IHI)<br>
*>           If WANTZ is .FALSE., then Z is not referenced.<br>
*>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is<br>
*>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the<br>
*>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).<br>
*>           (The output value of Z when INFO.GT.0 is given under<br>
*>           the description of INFO below.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>           The leading dimension of the array Z.  if WANTZ is .TRUE.<br>
*>           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension LWORK<br>
*>           On exit, if LWORK = -1, WORK(1) returns an estimate of<br>
*>           the optimal value for LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK.  LWORK .GE. max(1,N)<br>
*>           is sufficient, but LWORK typically as large as 6*N may<br>
*>           be required for optimal performance.  A workspace query<br>
*>           to determine the optimal workspace size is recommended.<br>
*><br>
*>           If LWORK = -1, then CLAQR0 does a workspace query.<br>
*>           In this case, CLAQR0 checks the input parameters and<br>
*>           estimates the optimal workspace size for the given<br>
*>           values of N, ILO and IHI.  The estimate is returned<br>
*>           in WORK(1).  No error message related to LWORK is<br>
*>           issued by XERBLA.  Neither H nor Z are accessed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>             =  0:  successful exit<br>
*>           .GT. 0:  if INFO = i, CLAQR0 failed to compute all of<br>
*>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR<br>
*>                and WI contain those eigenvalues which have been<br>
*>                successfully computed.  (Failures are rare.)<br>
*><br>
*>                If INFO .GT. 0 and WANT is .FALSE., then on exit,<br>
*>                the remaining unconverged eigenvalues are the eigen-<br>
*>                values of the upper Hessenberg matrix rows and<br>
*>                columns ILO through INFO of the final, output<br>
*>                value of H.<br>
*><br>
*>                If INFO .GT. 0 and WANTT is .TRUE., then on exit<br>
*><br>
*>           (*)  (initial value of H)*U  = U*(final value of H)<br>
*><br>
*>                where U is a unitary matrix.  The final<br>
*>                value of  H is upper Hessenberg and triangular in<br>
*>                rows and columns INFO+1 through IHI.<br>
*><br>
*>                If INFO .GT. 0 and WANTZ is .TRUE., then on exit<br>
*><br>
*>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)<br>
*>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U<br>
*><br>
*>                where U is the unitary matrix in (*) (regard-<br>
*>                less of the value of WANTT.)<br>
*><br>
*>                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not<br>
*>                accessed.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3<br>
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages<br>
*>       929--947, 2002.<br>
*> \n<br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal<br>
*>       of Matrix Analysis, volume 23, pages 948--973, 2002.<br>
*><br>
*  =====================================================================<br>
*/
	public void claqr0_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] W,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b CLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQR1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQR1( N, H, LDH, S1, S2, V )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX            S1, S2<br>
*       INTEGER            LDH, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            H( LDH, * ), V( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>      Given a 2-by-2 or 3-by-3 matrix H, CLAQR1 sets v to a<br>
*>      scalar multiple of the first column of the product<br>
*><br>
*>      (*)  K = (H - s1*I)*(H - s2*I)<br>
*><br>
*>      scaling to avoid overflows and most underflows.<br>
*><br>
*>      This is useful for starting double implicit shift bulges<br>
*>      in the QR algorithm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is integer<br>
*>              Order of the matrix H. N must be either 2 or 3.<br>
*> \endverbatim<br>
*><br>
*> \param[in] H<br>
*> \verbatim<br>
*>          H is COMPLEX array of dimension (LDH,N)<br>
*>              The 2-by-2 or 3-by-3 matrix H in (*).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is integer<br>
*>              The leading dimension of H as declared in<br>
*>              the calling procedure.  LDH.GE.N<br>
*> \endverbatim<br>
*><br>
*> \param[in] S1<br>
*> \verbatim<br>
*>          S1 is COMPLEX<br>
*> \endverbatim<br>
*><br>
*> \param[in] S2<br>
*> \verbatim<br>
*>          S2 is COMPLEX<br>
*><br>
*>          S1 and S2 are the shifts defining K in (*) above.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is COMPLEX array of dimension N<br>
*>              A scalar multiple of the first column of the<br>
*>              matrix K in (*).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void claqr1_(INTEGER N,float[] H,INTEGER LDH,float[] S1,float[] S2,float[] V);
/**
*> \brief \b CLAQR2 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,<br>
*                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,<br>
*                          NV, WV, LDWV, WORK, LWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,<br>
*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),<br>
*      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLAQR2 is identical to CLAQR3 except that it avoids<br>
*>    recursion by calling CLAHQR instead of CLAQR4.<br>
*><br>
*>    Aggressive early deflation:<br>
*><br>
*>    This subroutine accepts as input an upper Hessenberg matrix<br>
*>    H and performs an unitary similarity transformation<br>
*>    designed to detect and deflate fully converged eigenvalues from<br>
*>    a trailing principal submatrix.  On output H has been over-<br>
*>    written by a new Hessenberg matrix that is a perturbation of<br>
*>    an unitary similarity transformation of H.  It is to be<br>
*>    hoped that the final version of H has many zero subdiagonal<br>
*>    entries.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTT<br>
*> \verbatim<br>
*>          WANTT is LOGICAL<br>
*>          If .TRUE., then the Hessenberg matrix H is fully updated<br>
*>          so that the triangular Schur factor may be<br>
*>          computed (in cooperation with the calling subroutine).<br>
*>          If .FALSE., then only enough of H is updated to preserve<br>
*>          the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          If .TRUE., then the unitary matrix Z is updated so<br>
*>          so that the unitary Schur factor may be computed<br>
*>          (in cooperation with the calling subroutine).<br>
*>          If .FALSE., then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix H and (if WANTZ is .TRUE.) the<br>
*>          order of the unitary matrix Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KTOP<br>
*> \verbatim<br>
*>          KTOP is INTEGER<br>
*>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.<br>
*>          KBOT and KTOP together determine an isolated block<br>
*>          along the diagonal of the Hessenberg matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KBOT<br>
*> \verbatim<br>
*>          KBOT is INTEGER<br>
*>          It is assumed without a check that either<br>
*>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together<br>
*>          determine an isolated block along the diagonal of the<br>
*>          Hessenberg matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NW<br>
*> \verbatim<br>
*>          NW is INTEGER<br>
*>          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is COMPLEX array, dimension (LDH,N)<br>
*>          On input the initial N-by-N section of H stores the<br>
*>          Hessenberg matrix undergoing aggressive early deflation.<br>
*>          On output H has been transformed by a unitary<br>
*>          similarity transformation, perturbed, and the returned<br>
*>          to Hessenberg form that (it is to be hoped) has some<br>
*>          zero subdiagonal entries.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is integer<br>
*>          Leading dimension of H just as declared in the calling<br>
*>          subroutine.  N .LE. LDH<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILOZ<br>
*> \verbatim<br>
*>          ILOZ is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHIZ<br>
*> \verbatim<br>
*>          IHIZ is INTEGER<br>
*>          Specify the rows of Z to which transformations must be<br>
*>          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (LDZ,N)<br>
*>          IF WANTZ is .TRUE., then on output, the unitary<br>
*>          similarity transformation mentioned above has been<br>
*>          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.<br>
*>          If WANTZ is .FALSE., then Z is unreferenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is integer<br>
*>          The leading dimension of Z just as declared in the<br>
*>          calling subroutine.  1 .LE. LDZ.<br>
*> \endverbatim<br>
*><br>
*> \param[out] NS<br>
*> \verbatim<br>
*>          NS is integer<br>
*>          The number of unconverged (ie approximate) eigenvalues<br>
*>          returned in SR and SI that may be used as shifts by the<br>
*>          calling subroutine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ND<br>
*> \verbatim<br>
*>          ND is integer<br>
*>          The number of converged eigenvalues uncovered by this<br>
*>          subroutine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SH<br>
*> \verbatim<br>
*>          SH is COMPLEX array, dimension KBOT<br>
*>          On output, approximate eigenvalues that may<br>
*>          be used for shifts are stored in SH(KBOT-ND-NS+1)<br>
*>          through SR(KBOT-ND).  Converged eigenvalues are<br>
*>          stored in SH(KBOT-ND+1) through SH(KBOT).<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension (LDV,NW)<br>
*>          An NW-by-NW work array.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is integer scalar<br>
*>          The leading dimension of V just as declared in the<br>
*>          calling subroutine.  NW .LE. LDV<br>
*> \endverbatim<br>
*><br>
*> \param[in] NH<br>
*> \verbatim<br>
*>          NH is integer scalar<br>
*>          The number of columns of T.  NH.GE.NW.<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is COMPLEX array, dimension (LDT,NW)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is integer<br>
*>          The leading dimension of T just as declared in the<br>
*>          calling subroutine.  NW .LE. LDT<br>
*> \endverbatim<br>
*><br>
*> \param[in] NV<br>
*> \verbatim<br>
*>          NV is integer<br>
*>          The number of rows of work array WV available for<br>
*>          workspace.  NV.GE.NW.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WV<br>
*> \verbatim<br>
*>          WV is COMPLEX array, dimension (LDWV,NW)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDWV<br>
*> \verbatim<br>
*>          LDWV is integer<br>
*>          The leading dimension of W just as declared in the<br>
*>          calling subroutine.  NW .LE. LDV<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension LWORK.<br>
*>          On exit, WORK(1) is set to an estimate of the optimal value<br>
*>          of LWORK for the given values of N, NW, KTOP and KBOT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is integer<br>
*>          The dimension of the work array WORK.  LWORK = 2*NW<br>
*>          suffices, but greater efficiency may result from larger<br>
*>          values of LWORK.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; CLAQR2<br>
*>          only estimates the optimal workspace size for the given<br>
*>          values of N, NW, KTOP and KBOT.  The estimate is returned<br>
*>          in WORK(1).  No error message related to LWORK is issued<br>
*>          by XERBLA.  Neither H nor Z are accessed.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void claqr2_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NW,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER NS,INTEGER ND,float[] SH,float[] V,INTEGER LDV,INTEGER NH,float[] T,INTEGER LDT,INTEGER NV,float[] WV,INTEGER LDWV,float[] WORK,INTEGER LWORK);
/**
*> \brief \b CLAQR3 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQR3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,<br>
*                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,<br>
*                          NV, WV, LDWV, WORK, LWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,<br>
*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),<br>
*      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    Aggressive early deflation:<br>
*><br>
*>    CLAQR3 accepts as input an upper Hessenberg matrix<br>
*>    H and performs an unitary similarity transformation<br>
*>    designed to detect and deflate fully converged eigenvalues from<br>
*>    a trailing principal submatrix.  On output H has been over-<br>
*>    written by a new Hessenberg matrix that is a perturbation of<br>
*>    an unitary similarity transformation of H.  It is to be<br>
*>    hoped that the final version of H has many zero subdiagonal<br>
*>    entries.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTT<br>
*> \verbatim<br>
*>          WANTT is LOGICAL<br>
*>          If .TRUE., then the Hessenberg matrix H is fully updated<br>
*>          so that the triangular Schur factor may be<br>
*>          computed (in cooperation with the calling subroutine).<br>
*>          If .FALSE., then only enough of H is updated to preserve<br>
*>          the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          If .TRUE., then the unitary matrix Z is updated so<br>
*>          so that the unitary Schur factor may be computed<br>
*>          (in cooperation with the calling subroutine).<br>
*>          If .FALSE., then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix H and (if WANTZ is .TRUE.) the<br>
*>          order of the unitary matrix Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KTOP<br>
*> \verbatim<br>
*>          KTOP is INTEGER<br>
*>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.<br>
*>          KBOT and KTOP together determine an isolated block<br>
*>          along the diagonal of the Hessenberg matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KBOT<br>
*> \verbatim<br>
*>          KBOT is INTEGER<br>
*>          It is assumed without a check that either<br>
*>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together<br>
*>          determine an isolated block along the diagonal of the<br>
*>          Hessenberg matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NW<br>
*> \verbatim<br>
*>          NW is INTEGER<br>
*>          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is COMPLEX array, dimension (LDH,N)<br>
*>          On input the initial N-by-N section of H stores the<br>
*>          Hessenberg matrix undergoing aggressive early deflation.<br>
*>          On output H has been transformed by a unitary<br>
*>          similarity transformation, perturbed, and the returned<br>
*>          to Hessenberg form that (it is to be hoped) has some<br>
*>          zero subdiagonal entries.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is integer<br>
*>          Leading dimension of H just as declared in the calling<br>
*>          subroutine.  N .LE. LDH<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILOZ<br>
*> \verbatim<br>
*>          ILOZ is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHIZ<br>
*> \verbatim<br>
*>          IHIZ is INTEGER<br>
*>          Specify the rows of Z to which transformations must be<br>
*>          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (LDZ,N)<br>
*>          IF WANTZ is .TRUE., then on output, the unitary<br>
*>          similarity transformation mentioned above has been<br>
*>          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.<br>
*>          If WANTZ is .FALSE., then Z is unreferenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is integer<br>
*>          The leading dimension of Z just as declared in the<br>
*>          calling subroutine.  1 .LE. LDZ.<br>
*> \endverbatim<br>
*><br>
*> \param[out] NS<br>
*> \verbatim<br>
*>          NS is integer<br>
*>          The number of unconverged (ie approximate) eigenvalues<br>
*>          returned in SR and SI that may be used as shifts by the<br>
*>          calling subroutine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ND<br>
*> \verbatim<br>
*>          ND is integer<br>
*>          The number of converged eigenvalues uncovered by this<br>
*>          subroutine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SH<br>
*> \verbatim<br>
*>          SH is COMPLEX array, dimension KBOT<br>
*>          On output, approximate eigenvalues that may<br>
*>          be used for shifts are stored in SH(KBOT-ND-NS+1)<br>
*>          through SR(KBOT-ND).  Converged eigenvalues are<br>
*>          stored in SH(KBOT-ND+1) through SH(KBOT).<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension (LDV,NW)<br>
*>          An NW-by-NW work array.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is integer scalar<br>
*>          The leading dimension of V just as declared in the<br>
*>          calling subroutine.  NW .LE. LDV<br>
*> \endverbatim<br>
*><br>
*> \param[in] NH<br>
*> \verbatim<br>
*>          NH is integer scalar<br>
*>          The number of columns of T.  NH.GE.NW.<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is COMPLEX array, dimension (LDT,NW)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is integer<br>
*>          The leading dimension of T just as declared in the<br>
*>          calling subroutine.  NW .LE. LDT<br>
*> \endverbatim<br>
*><br>
*> \param[in] NV<br>
*> \verbatim<br>
*>          NV is integer<br>
*>          The number of rows of work array WV available for<br>
*>          workspace.  NV.GE.NW.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WV<br>
*> \verbatim<br>
*>          WV is COMPLEX array, dimension (LDWV,NW)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDWV<br>
*> \verbatim<br>
*>          LDWV is integer<br>
*>          The leading dimension of W just as declared in the<br>
*>          calling subroutine.  NW .LE. LDV<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension LWORK.<br>
*>          On exit, WORK(1) is set to an estimate of the optimal value<br>
*>          of LWORK for the given values of N, NW, KTOP and KBOT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is integer<br>
*>          The dimension of the work array WORK.  LWORK = 2*NW<br>
*>          suffices, but greater efficiency may result from larger<br>
*>          values of LWORK.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; CLAQR3<br>
*>          only estimates the optimal workspace size for the given<br>
*>          values of N, NW, KTOP and KBOT.  The estimate is returned<br>
*>          in WORK(1).  No error message related to LWORK is issued<br>
*>          by XERBLA.  Neither H nor Z are accessed.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void claqr3_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NW,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER NS,INTEGER ND,float[] SH,float[] V,INTEGER LDV,INTEGER NH,float[] T,INTEGER LDT,INTEGER NV,float[] WV,INTEGER LDWV,float[] WORK,INTEGER LWORK);
/**
*> \brief \b CLAQR4 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQR4 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr4.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr4.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr4.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,<br>
*                          IHIZ, Z, LDZ, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLAQR4 implements one level of recursion for CLAQR0.<br>
*>    It is a complete implementation of the small bulge multi-shift<br>
*>    QR algorithm.  It may be called by CLAQR0 and, for large enough<br>
*>    deflation window size, it may be called by CLAQR3.  This<br>
*>    subroutine is identical to CLAQR0 except that it calls CLAQR2<br>
*>    instead of CLAQR3.<br>
*><br>
*>    CLAQR4 computes the eigenvalues of a Hessenberg matrix H<br>
*>    and, optionally, the matrices T and Z from the Schur decomposition<br>
*>    H = Z T Z**H, where T is an upper triangular matrix (the<br>
*>    Schur form), and Z is the unitary matrix of Schur vectors.<br>
*><br>
*>    Optionally Z may be postmultiplied into an input unitary<br>
*>    matrix Q so that this routine can give the Schur factorization<br>
*>    of a matrix A which has been reduced to the Hessenberg form H<br>
*>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTT<br>
*> \verbatim<br>
*>          WANTT is LOGICAL<br>
*>          = .TRUE. : the full Schur form T is required;<br>
*>          = .FALSE.: only eigenvalues are required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          = .TRUE. : the matrix of Schur vectors Z is required;<br>
*>          = .FALSE.: Schur vectors are not required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           The order of the matrix H.  N .GE. 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILO<br>
*> \verbatim<br>
*>          ILO is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHI<br>
*> \verbatim<br>
*>          IHI is INTEGER<br>
*>           It is assumed that H is already upper triangular in rows<br>
*>           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,<br>
*>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a<br>
*>           previous call to CGEBAL, and then passed to CGEHRD when the<br>
*>           matrix output by CGEBAL is reduced to Hessenberg form.<br>
*>           Otherwise, ILO and IHI should be set to 1 and N,<br>
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.<br>
*>           If N = 0, then ILO = 1 and IHI = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is COMPLEX array, dimension (LDH,N)<br>
*>           On entry, the upper Hessenberg matrix H.<br>
*>           On exit, if INFO = 0 and WANTT is .TRUE., then H<br>
*>           contains the upper triangular matrix T from the Schur<br>
*>           decomposition (the Schur form). If INFO = 0 and WANT is<br>
*>           .FALSE., then the contents of H are unspecified on exit.<br>
*>           (The output value of H when INFO.GT.0 is given under the<br>
*>           description of INFO below.)<br>
*><br>
*>           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and<br>
*>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>           The leading dimension of the array H. LDH .GE. max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (N)<br>
*>           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored<br>
*>           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are<br>
*>           stored in the same order as on the diagonal of the Schur<br>
*>           form returned in H, with W(i) = H(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILOZ<br>
*> \verbatim<br>
*>          ILOZ is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHIZ<br>
*> \verbatim<br>
*>          IHIZ is INTEGER<br>
*>           Specify the rows of Z to which transformations must be<br>
*>           applied if WANTZ is .TRUE..<br>
*>           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (LDZ,IHI)<br>
*>           If WANTZ is .FALSE., then Z is not referenced.<br>
*>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is<br>
*>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the<br>
*>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).<br>
*>           (The output value of Z when INFO.GT.0 is given under<br>
*>           the description of INFO below.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>           The leading dimension of the array Z.  if WANTZ is .TRUE.<br>
*>           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension LWORK<br>
*>           On exit, if LWORK = -1, WORK(1) returns an estimate of<br>
*>           the optimal value for LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK.  LWORK .GE. max(1,N)<br>
*>           is sufficient, but LWORK typically as large as 6*N may<br>
*>           be required for optimal performance.  A workspace query<br>
*>           to determine the optimal workspace size is recommended.<br>
*><br>
*>           If LWORK = -1, then CLAQR4 does a workspace query.<br>
*>           In this case, CLAQR4 checks the input parameters and<br>
*>           estimates the optimal workspace size for the given<br>
*>           values of N, ILO and IHI.  The estimate is returned<br>
*>           in WORK(1).  No error message related to LWORK is<br>
*>           issued by XERBLA.  Neither H nor Z are accessed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>             =  0:  successful exit<br>
*>           .GT. 0:  if INFO = i, CLAQR4 failed to compute all of<br>
*>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR<br>
*>                and WI contain those eigenvalues which have been<br>
*>                successfully computed.  (Failures are rare.)<br>
*><br>
*>                If INFO .GT. 0 and WANT is .FALSE., then on exit,<br>
*>                the remaining unconverged eigenvalues are the eigen-<br>
*>                values of the upper Hessenberg matrix rows and<br>
*>                columns ILO through INFO of the final, output<br>
*>                value of H.<br>
*><br>
*>                If INFO .GT. 0 and WANTT is .TRUE., then on exit<br>
*><br>
*>           (*)  (initial value of H)*U  = U*(final value of H)<br>
*><br>
*>                where U is a unitary matrix.  The final<br>
*>                value of  H is upper Hessenberg and triangular in<br>
*>                rows and columns INFO+1 through IHI.<br>
*><br>
*>                If INFO .GT. 0 and WANTZ is .TRUE., then on exit<br>
*><br>
*>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)<br>
*>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U<br>
*><br>
*>                where U is the unitary matrix in (*) (regard-<br>
*>                less of the value of WANTT.)<br>
*><br>
*>                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not<br>
*>                accessed.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3<br>
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages<br>
*>       929--947, 2002.<br>
*> \n<br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal<br>
*>       of Matrix Analysis, volume 23, pages 948--973, 2002.<br>
*><br>
*  =====================================================================<br>
*/
	public void claqr4_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] W,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b CLAQR5 performs a single small-bulge multi-shift QR sweep.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQR5 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr5.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr5.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr5.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,<br>
*                          H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,<br>
*                          WV, LDWV, NH, WH, LDWH )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,<br>
*      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ),<br>
*      $                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLAQR5 called by CLAQR0 performs a<br>
*>    single small-bulge multi-shift QR sweep.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTT<br>
*> \verbatim<br>
*>          WANTT is logical scalar<br>
*>             WANTT = .true. if the triangular Schur factor<br>
*>             is being computed.  WANTT is set to .false. otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is logical scalar<br>
*>             WANTZ = .true. if the unitary Schur factor is being<br>
*>             computed.  WANTZ is set to .false. otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KACC22<br>
*> \verbatim<br>
*>          KACC22 is integer with value 0, 1, or 2.<br>
*>             Specifies the computation mode of far-from-diagonal<br>
*>             orthogonal updates.<br>
*>        = 0: CLAQR5 does not accumulate reflections and does not<br>
*>             use matrix-matrix multiply to update far-from-diagonal<br>
*>             matrix entries.<br>
*>        = 1: CLAQR5 accumulates reflections and uses matrix-matrix<br>
*>             multiply to update the far-from-diagonal matrix entries.<br>
*>        = 2: CLAQR5 accumulates reflections, uses matrix-matrix<br>
*>             multiply to update the far-from-diagonal matrix entries,<br>
*>             and takes advantage of 2-by-2 block structure during<br>
*>             matrix multiplies.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is integer scalar<br>
*>             N is the order of the Hessenberg matrix H upon which this<br>
*>             subroutine operates.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KTOP<br>
*> \verbatim<br>
*>          KTOP is integer scalar<br>
*> \endverbatim<br>
*><br>
*> \param[in] KBOT<br>
*> \verbatim<br>
*>          KBOT is integer scalar<br>
*>             These are the first and last rows and columns of an<br>
*>             isolated diagonal block upon which the QR sweep is to be<br>
*>             applied. It is assumed without a check that<br>
*>                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0<br>
*>             and<br>
*>                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NSHFTS<br>
*> \verbatim<br>
*>          NSHFTS is integer scalar<br>
*>             NSHFTS gives the number of simultaneous shifts.  NSHFTS<br>
*>             must be positive and even.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] S<br>
*> \verbatim<br>
*>          S is COMPLEX array of size (NSHFTS)<br>
*>             S contains the shifts of origin that define the multi-<br>
*>             shift QR sweep.  On output S may be reordered.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is COMPLEX array of size (LDH,N)<br>
*>             On input H contains a Hessenberg matrix.  On output a<br>
*>             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied<br>
*>             to the isolated diagonal block in rows and columns KTOP<br>
*>             through KBOT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is integer scalar<br>
*>             LDH is the leading dimension of H just as declared in the<br>
*>             calling procedure.  LDH.GE.MAX(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILOZ<br>
*> \verbatim<br>
*>          ILOZ is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHIZ<br>
*> \verbatim<br>
*>          IHIZ is INTEGER<br>
*>             Specify the rows of Z to which transformations must be<br>
*>             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array of size (LDZ,IHIZ)<br>
*>             If WANTZ = .TRUE., then the QR Sweep unitary<br>
*>             similarity transformation is accumulated into<br>
*>             Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.<br>
*>             If WANTZ = .FALSE., then Z is unreferenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is integer scalar<br>
*>             LDA is the leading dimension of Z just as declared in<br>
*>             the calling procedure. LDZ.GE.N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is COMPLEX array of size (LDV,NSHFTS/2)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is integer scalar<br>
*>             LDV is the leading dimension of V as declared in the<br>
*>             calling procedure.  LDV.GE.3.<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is COMPLEX array of size<br>
*>             (LDU,3*NSHFTS-3)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is integer scalar<br>
*>             LDU is the leading dimension of U just as declared in the<br>
*>             in the calling subroutine.  LDU.GE.3*NSHFTS-3.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NH<br>
*> \verbatim<br>
*>          NH is integer scalar<br>
*>             NH is the number of columns in array WH available for<br>
*>             workspace. NH.GE.1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WH<br>
*> \verbatim<br>
*>          WH is COMPLEX array of size (LDWH,NH)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDWH<br>
*> \verbatim<br>
*>          LDWH is integer scalar<br>
*>             Leading dimension of WH just as declared in the<br>
*>             calling procedure.  LDWH.GE.3*NSHFTS-3.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NV<br>
*> \verbatim<br>
*>          NV is integer scalar<br>
*>             NV is the number of rows in WV agailable for workspace.<br>
*>             NV.GE.1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WV<br>
*> \verbatim<br>
*>          WV is COMPLEX array of size<br>
*>             (LDWV,3*NSHFTS-3)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDWV<br>
*> \verbatim<br>
*>          LDWV is integer scalar<br>
*>             LDWV is the leading dimension of WV as declared in the<br>
*>             in the calling subroutine.  LDWV.GE.NV.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3<br>
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages<br>
*>       929--947, 2002.<br>
*><br>
*  =====================================================================<br>
*/
	public void claqr5_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER KACC22,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NSHFTS,float[] S,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] V,INTEGER LDV,float[] U,INTEGER LDU,INTEGER NV,float[] WV,INTEGER LDWV,INTEGER NH,float[] WH,INTEGER LDWH);
/**
*> \brief \b CLAQSB scales a symmetric/Hermitian band matrix, using scaling factors computed by spbequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQSB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqsb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqsb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqsb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQSB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            KD, LDAB, N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               S( * )<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQSB equilibrates a symmetric band matrix A using the scaling<br>
*> factors in the vector S.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is stored.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KD<br>
*> \verbatim<br>
*>          KD is INTEGER<br>
*>          The number of super-diagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*><br>
*>          On exit, if INFO = 0, the triangular factor U or L from the<br>
*>          Cholesky factorization A = U**H *U or A = L*L**H of the band<br>
*>          matrix A, in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          Ratio of the smallest S(i) to the largest S(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix entry.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies whether or not equilibration was done.<br>
*>          = 'N':  No equilibration.<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  THRESH is a threshold value used to decide if scaling should be done<br>
*>  based on the ratio of the scaling factors.  If SCOND < THRESH,<br>
*>  scaling is done.<br>
*><br>
*>  LARGE and SMALL are threshold values used to decide if scaling should<br>
*>  be done based on the absolute size of the largest matrix element.<br>
*>  If AMAX > LARGE or AMAX < SMALL, scaling is done.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claqsb_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b CLAQSP scales a symmetric/Hermitian matrix in packed storage, using scaling factors computed by sppequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQSP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqsp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqsp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqsp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQSP( UPLO, N, AP, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               S( * )<br>
*       COMPLEX            AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQSP equilibrates a symmetric matrix A using the scaling factors<br>
*> in the vector S.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is stored.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the equilibrated matrix:  diag(S) * A * diag(S), in<br>
*>          the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          Ratio of the smallest S(i) to the largest S(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix entry.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies whether or not equilibration was done.<br>
*>          = 'N':  No equilibration.<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  THRESH is a threshold value used to decide if scaling should be done<br>
*>  based on the ratio of the scaling factors.  If SCOND < THRESH,<br>
*>  scaling is done.<br>
*><br>
*>  LARGE and SMALL are threshold values used to decide if scaling should<br>
*>  be done based on the absolute size of the largest matrix element.<br>
*>  If AMAX > LARGE or AMAX < SMALL, scaling is done.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claqsp_(CHARACTER UPLO,INTEGER N,float[] AP,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b CLAQSY scales a symmetric/Hermitian matrix, using scaling factors computed by spoequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAQSY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqsy.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqsy.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqsy.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            LDA, N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               S( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAQSY equilibrates a symmetric matrix A using the scaling factors<br>
*> in the vector S.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is stored.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          n by n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n by n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if EQUED = 'Y', the equilibrated matrix:<br>
*>          diag(S) * A * diag(S).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(N,1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          Ratio of the smallest S(i) to the largest S(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix entry.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies whether or not equilibration was done.<br>
*>          = 'N':  No equilibration.<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  THRESH is a threshold value used to decide if scaling should be done<br>
*>  based on the ratio of the scaling factors.  If SCOND < THRESH,<br>
*>  scaling is done.<br>
*><br>
*>  LARGE and SMALL are threshold values used to decide if scaling should<br>
*>  be done based on the absolute size of the largest matrix element.<br>
*>  If AMAX > LARGE or AMAX < SMALL, scaling is done.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexSYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claqsy_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b CLAR1V computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn of the tridiagonal matrix LDLT - λI.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAR1V + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clar1v.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clar1v.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clar1v.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD,<br>
*                  PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA,<br>
*                  R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            WANTNC<br>
*       INTEGER   B1, BN, N, NEGCNT, R<br>
*       REAL               GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID,<br>
*      $                   RQCORR, ZTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISUPPZ( * )<br>
*       REAL               D( * ), L( * ), LD( * ), LLD( * ),<br>
*      $                  WORK( * )<br>
*       COMPLEX          Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAR1V computes the (scaled) r-th column of the inverse of<br>
*> the sumbmatrix in rows B1 through BN of the tridiagonal matrix<br>
*> L D L**T - sigma I. When sigma is close to an eigenvalue, the<br>
*> computed vector is an accurate eigenvector. Usually, r corresponds<br>
*> to the index where the eigenvector is largest in magnitude.<br>
*> The following steps accomplish this computation :<br>
*> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,<br>
*> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,<br>
*> (c) Computation of the diagonal elements of the inverse of<br>
*>     L D L**T - sigma I by combining the above transforms, and choosing<br>
*>     r as the index where the diagonal of the inverse is (one of the)<br>
*>     largest in magnitude.<br>
*> (d) Computation of the (scaled) r-th column of the inverse using the<br>
*>     twisted factorization obtained by combining the top part of the<br>
*>     the stationary and the bottom part of the progressive transform.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           The order of the matrix L D L**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B1<br>
*> \verbatim<br>
*>          B1 is INTEGER<br>
*>           First index of the submatrix of L D L**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BN<br>
*> \verbatim<br>
*>          BN is INTEGER<br>
*>           Last index of the submatrix of L D L**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LAMBDA<br>
*> \verbatim<br>
*>          LAMBDA is REAL<br>
*>           The shift. In order to compute an accurate eigenvector,<br>
*>           LAMBDA should be a good approximation to an eigenvalue<br>
*>           of L D L**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is REAL array, dimension (N-1)<br>
*>           The (n-1) subdiagonal elements of the unit bidiagonal matrix<br>
*>           L, in elements 1 to N-1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>           The n diagonal elements of the diagonal matrix D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LD<br>
*> \verbatim<br>
*>          LD is REAL array, dimension (N-1)<br>
*>           The n-1 elements L(i)*D(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LLD<br>
*> \verbatim<br>
*>          LLD is REAL array, dimension (N-1)<br>
*>           The n-1 elements L(i)*L(i)*D(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>           The minimum pivot in the Sturm sequence.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GAPTOL<br>
*> \verbatim<br>
*>          GAPTOL is REAL<br>
*>           Tolerance that indicates when eigenvector entries are negligible<br>
*>           w.r.t. their contribution to the residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (N)<br>
*>           On input, all entries of Z must be set to 0.<br>
*>           On output, Z contains the (scaled) r-th column of the<br>
*>           inverse. The scaling is such that Z(R) equals 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTNC<br>
*> \verbatim<br>
*>          WANTNC is LOGICAL<br>
*>           Specifies whether NEGCNT has to be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] NEGCNT<br>
*> \verbatim<br>
*>          NEGCNT is INTEGER<br>
*>           If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin<br>
*>           in the  matrix factorization L D L**T, and NEGCNT = -1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ZTZ<br>
*> \verbatim<br>
*>          ZTZ is REAL<br>
*>           The square of the 2-norm of Z.<br>
*> \endverbatim<br>
*><br>
*> \param[out] MINGMA<br>
*> \verbatim<br>
*>          MINGMA is REAL<br>
*>           The reciprocal of the largest (in magnitude) diagonal<br>
*>           element of the inverse of L D L**T - sigma I.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] R<br>
*> \verbatim<br>
*>          R is INTEGER<br>
*>           The twist index for the twisted factorization used to<br>
*>           compute Z.<br>
*>           On input, 0 <= R <= N. If R is input as 0, R is set to<br>
*>           the index where (L D L**T - sigma I)^{-1} is largest<br>
*>           in magnitude. If 1 <= R <= N, R is unchanged.<br>
*>           On output, R contains the twist index used to compute Z.<br>
*>           Ideally, R designates the position of the maximum entry in the<br>
*>           eigenvector.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISUPPZ<br>
*> \verbatim<br>
*>          ISUPPZ is INTEGER array, dimension (2)<br>
*>           The support of the vector in Z, i.e., the vector Z is<br>
*>           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] NRMINV<br>
*> \verbatim<br>
*>          NRMINV is REAL<br>
*>           NRMINV = 1/SQRT( ZTZ )<br>
*> \endverbatim<br>
*><br>
*> \param[out] RESID<br>
*> \verbatim<br>
*>          RESID is REAL<br>
*>           The residual of the FP vector.<br>
*>           RESID = ABS( MINGMA )/SQRT( ZTZ )<br>
*> \endverbatim<br>
*><br>
*> \param[out] RQCORR<br>
*> \verbatim<br>
*>          RQCORR is REAL<br>
*>           The Rayleigh Quotient correction to LAMBDA.<br>
*>           RQCORR = MINGMA*TMP<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (4*N)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Beresford Parlett, University of California, Berkeley, USA \n<br>
*> Jim Demmel, University of California, Berkeley, USA \n<br>
*> Inderjit Dhillon, University of Texas, Austin, USA \n<br>
*> Osni Marques, LBNL/NERSC, USA \n<br>
*> Christof Voemel, University of California, Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void clar1v_(INTEGER N,INTEGER B1,INTEGER BN,REAL LAMBDA,float[] D,float[] L,float[] LD,float[] LLD,REAL PIVMIN,REAL GAPTOL,float[] Z,LOGICAL WANTNC,INTEGER NEGCNT,REAL ZTZ,REAL MINGMA,INTEGER R,int[] ISUPPZ,REAL NRMINV,REAL RESID,REAL RQCORR,float[] WORK);
/**
*> \brief \b CLAR2V applies a vector of plane rotations with real cosines and complex sines from both sides to a sequence of 2-by-2 symmetric/Hermitian matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAR2V + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clar2v.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clar2v.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clar2v.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAR2V( N, X, Y, Z, INCX, C, S, INCC )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCC, INCX, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * )<br>
*       COMPLEX            S( * ), X( * ), Y( * ), Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAR2V applies a vector of complex plane rotations with real cosines<br>
*> from both sides to a sequence of 2-by-2 complex Hermitian matrices,<br>
*> defined by the elements of the vectors x, y and z. For i = 1,2,...,n<br>
*><br>
*>    (       x(i)  z(i) ) :=<br>
*>    ( conjg(z(i)) y(i) )<br>
*><br>
*>      (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )<br>
*>      ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of plane rotations to be applied.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (1+(N-1)*INCX)<br>
*>          The vector x; the elements of x are assumed to be real.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (1+(N-1)*INCX)<br>
*>          The vector y; the elements of y are assumed to be real.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (1+(N-1)*INCX)<br>
*>          The vector z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between elements of X, Y and Z. INCX > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (1+(N-1)*INCC)<br>
*>          The cosines of the plane rotations.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is COMPLEX array, dimension (1+(N-1)*INCC)<br>
*>          The sines of the plane rotations.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCC<br>
*> \verbatim<br>
*>          INCC is INTEGER<br>
*>          The increment between elements of C and S. INCC > 0.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clar2v_(INTEGER N,float[] X,float[] Y,float[] Z,INTEGER INCX,float[] C,float[] S,INTEGER INCC);
/**
*> \brief \b CLARCM copies all or part of a real two-dimensional array to a complex array.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARCM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarcm.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarcm.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarcm.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, LDB, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), RWORK( * )<br>
*       COMPLEX            B( LDB, * ), C( LDC, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARCM performs a very simple matrix-matrix multiplication:<br>
*>          C := A * B,<br>
*> where A is M by M and real; B is M by N and complex;<br>
*> C is M by N and complex.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A and of the matrix C.<br>
*>          M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns and rows of the matrix B and<br>
*>          the number of columns of the matrix C.<br>
*>          N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA, M)<br>
*>          A contains the M by M matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >=max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB, N)<br>
*>          B contains the M by N matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >=max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is COMPLEX array, dimension (LDC, N)<br>
*>          C contains the M by N matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDC >=max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (2*M*N)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clarcm_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] C,INTEGER LDC,float[] RWORK);
/**
*> \brief \b CLARF applies an elementary reflector to a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE<br>
*       INTEGER            INCV, LDC, M, N<br>
*       COMPLEX            TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            C( LDC, * ), V( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARF applies a complex elementary reflector H to a complex M-by-N<br>
*> matrix C, from either the left or the right. H is represented in the<br>
*> form<br>
*><br>
*>       H = I - tau * v * v**H<br>
*><br>
*> where tau is a complex scalar and v is a complex vector.<br>
*><br>
*> If tau = 0, then H is taken to be the unit matrix.<br>
*><br>
*> To apply H**H (the conjugate transpose of H), supply conjg(tau) instead<br>
*> tau.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': form  H * C<br>
*>          = 'R': form  C * H<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension<br>
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'<br>
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'<br>
*>          The vector v in the representation of H. V is not used if<br>
*>          TAU = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCV<br>
*> \verbatim<br>
*>          INCV is INTEGER<br>
*>          The increment between elements of v. INCV <> 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX<br>
*>          The value tau in the representation of H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',<br>
*>          or C * H if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDC >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension<br>
*>                         (N) if SIDE = 'L'<br>
*>                      or (M) if SIDE = 'R'<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clarf_(CHARACTER SIDE,INTEGER M,INTEGER N,float[] V,INTEGER INCV,float[] TAU,float[] C,INTEGER LDC,float[] WORK);
/**
*> \brief \b CLARFB applies a block reflector or its conjugate-transpose to a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARFB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,<br>
*                          T, LDT, C, LDC, WORK, LDWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, SIDE, STOREV, TRANS<br>
*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ),<br>
*      $                   WORK( LDWORK, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARFB applies a complex block reflector H or its transpose H**H to a<br>
*> complex M-by-N matrix C, from either the left or the right.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply H or H**H from the Left<br>
*>          = 'R': apply H or H**H from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply H (No transpose)<br>
*>          = 'C': apply H**H (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIRECT<br>
*> \verbatim<br>
*>          DIRECT is CHARACTER*1<br>
*>          Indicates how H is formed from a product of elementary<br>
*>          reflectors<br>
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)<br>
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)<br>
*> \endverbatim<br>
*><br>
*> \param[in] STOREV<br>
*> \verbatim<br>
*>          STOREV is CHARACTER*1<br>
*>          Indicates how the vectors which define the elementary<br>
*>          reflectors are stored:<br>
*>          = 'C': Columnwise<br>
*>          = 'R': Rowwise<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The order of the matrix T (= the number of elementary<br>
*>          reflectors whose product defines the block reflector).<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension<br>
*>                                (LDV,K) if STOREV = 'C'<br>
*>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'<br>
*>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'<br>
*>          The matrix V. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V.<br>
*>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);<br>
*>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);<br>
*>          if STOREV = 'R', LDV >= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] T<br>
*> \verbatim<br>
*>          T is COMPLEX array, dimension (LDT,K)<br>
*>          The triangular K-by-K matrix T in the representation of the<br>
*>          block reflector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. LDT >= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDC >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (LDWORK,K)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDWORK<br>
*> \verbatim<br>
*>          LDWORK is INTEGER<br>
*>          The leading dimension of the array WORK.<br>
*>          If SIDE = 'L', LDWORK >= max(1,N);<br>
*>          if SIDE = 'R', LDWORK >= max(1,M).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2013<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The shape of the matrix V and the storage of the vectors which define<br>
*>  the H(i) is best illustrated by the following example with n = 5 and<br>
*>  k = 3. The elements equal to 1 are not stored; the corresponding<br>
*>  array elements are modified but restored on exit. The rest of the<br>
*>  array is not used.<br>
*><br>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':<br>
*><br>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )<br>
*>                   ( v1  1    )                     (     1 v2 v2 v2 )<br>
*>                   ( v1 v2  1 )                     (        1 v3 v3 )<br>
*>                   ( v1 v2 v3 )<br>
*>                   ( v1 v2 v3 )<br>
*><br>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':<br>
*><br>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )<br>
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )<br>
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )<br>
*>                   (     1 v3 )<br>
*>                   (        1 )<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clarfb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] C,INTEGER LDC,float[] WORK,INTEGER LDWORK);
/**
*> \brief \b CLARFG generates an elementary reflector (Householder matrix).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARFG + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfg.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfg.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfg.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       COMPLEX            ALPHA, TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARFG generates a complex elementary reflector H of order n, such<br>
*> that<br>
*><br>
*>       H**H * ( alpha ) = ( beta ),   H**H * H = I.<br>
*>              (   x   )   (   0  )<br>
*><br>
*> where alpha and beta are scalars, with beta real, and x is an<br>
*> (n-1)-element complex vector. H is represented in the form<br>
*><br>
*>       H = I - tau * ( 1 ) * ( 1 v**H ) ,<br>
*>                     ( v )<br>
*><br>
*> where tau is a complex scalar and v is a complex (n-1)-element<br>
*> vector. Note that H is not hermitian.<br>
*><br>
*> If the elements of x are all zero and alpha is real, then tau = 0<br>
*> and H is taken to be the unit matrix.<br>
*><br>
*> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the elementary reflector.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>          On entry, the value alpha.<br>
*>          On exit, it is overwritten with the value beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension<br>
*>                         (1+(N-2)*abs(INCX))<br>
*>          On entry, the vector x.<br>
*>          On exit, it is overwritten with the vector v.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between elements of X. INCX > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX<br>
*>          The value tau.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clarfg_(INTEGER N,float[] ALPHA,float[] X,INTEGER INCX,float[] TAU);
/**
*> \brief \b CLARFGP generates an elementary reflector (Householder matrix) with non-negative beta.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARFGP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfgp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfgp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfgp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARFGP( N, ALPHA, X, INCX, TAU )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       COMPLEX            ALPHA, TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARFGP generates a complex elementary reflector H of order n, such<br>
*> that<br>
*><br>
*>       H**H * ( alpha ) = ( beta ),   H**H * H = I.<br>
*>              (   x   )   (   0  )<br>
*><br>
*> where alpha and beta are scalars, beta is real and non-negative, and<br>
*> x is an (n-1)-element complex vector.  H is represented in the form<br>
*><br>
*>       H = I - tau * ( 1 ) * ( 1 v**H ) ,<br>
*>                     ( v )<br>
*><br>
*> where tau is a complex scalar and v is a complex (n-1)-element<br>
*> vector. Note that H is not hermitian.<br>
*><br>
*> If the elements of x are all zero and alpha is real, then tau = 0<br>
*> and H is taken to be the unit matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the elementary reflector.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>          On entry, the value alpha.<br>
*>          On exit, it is overwritten with the value beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension<br>
*>                         (1+(N-2)*abs(INCX))<br>
*>          On entry, the vector x.<br>
*>          On exit, it is overwritten with the vector v.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between elements of X. INCX > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX<br>
*>          The value tau.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2015<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clarfgp_(INTEGER N,float[] ALPHA,float[] X,INTEGER INCX,float[] TAU);
/**
*> \brief \b CLARFT forms the triangular factor T of a block reflector H = I - vtvH<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARFT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarft.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarft.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarft.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, STOREV<br>
*       INTEGER            K, LDT, LDV, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARFT forms the triangular factor T of a complex block reflector H<br>
*> of order n, which is defined as a product of k elementary reflectors.<br>
*><br>
*> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;<br>
*><br>
*> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.<br>
*><br>
*> If STOREV = 'C', the vector which defines the elementary reflector<br>
*> H(i) is stored in the i-th column of the array V, and<br>
*><br>
*>    H  =  I - V * T * V**H<br>
*><br>
*> If STOREV = 'R', the vector which defines the elementary reflector<br>
*> H(i) is stored in the i-th row of the array V, and<br>
*><br>
*>    H  =  I - V**H * T * V<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] DIRECT<br>
*> \verbatim<br>
*>          DIRECT is CHARACTER*1<br>
*>          Specifies the order in which the elementary reflectors are<br>
*>          multiplied to form the block reflector:<br>
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)<br>
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)<br>
*> \endverbatim<br>
*><br>
*> \param[in] STOREV<br>
*> \verbatim<br>
*>          STOREV is CHARACTER*1<br>
*>          Specifies how the vectors which define the elementary<br>
*>          reflectors are stored (see also Further Details):<br>
*>          = 'C': columnwise<br>
*>          = 'R': rowwise<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the block reflector H. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The order of the triangular factor T (= the number of<br>
*>          elementary reflectors). K >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension<br>
*>                               (LDV,K) if STOREV = 'C'<br>
*>                               (LDV,N) if STOREV = 'R'<br>
*>          The matrix V. See further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V.<br>
*>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is COMPLEX array, dimension (LDT,K)<br>
*>          The k by k triangular factor T of the block reflector.<br>
*>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is<br>
*>          lower triangular. The rest of the array is not used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. LDT >= K.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The shape of the matrix V and the storage of the vectors which define<br>
*>  the H(i) is best illustrated by the following example with n = 5 and<br>
*>  k = 3. The elements equal to 1 are not stored.<br>
*><br>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':<br>
*><br>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )<br>
*>                   ( v1  1    )                     (     1 v2 v2 v2 )<br>
*>                   ( v1 v2  1 )                     (        1 v3 v3 )<br>
*>                   ( v1 v2 v3 )<br>
*>                   ( v1 v2 v3 )<br>
*><br>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':<br>
*><br>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )<br>
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )<br>
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )<br>
*>                   (     1 v3 )<br>
*>                   (        1 )<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clarft_(CHARACTER DIRECT,CHARACTER STOREV,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] TAU,float[] T,INTEGER LDT);
/**
*> \brief \b CLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling when the reflector has order ≤ 10.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARFX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE<br>
*       INTEGER            LDC, M, N<br>
*       COMPLEX            TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            C( LDC, * ), V( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARFX applies a complex elementary reflector H to a complex m by n<br>
*> matrix C, from either the left or the right. H is represented in the<br>
*> form<br>
*><br>
*>       H = I - tau * v * v**H<br>
*><br>
*> where tau is a complex scalar and v is a complex vector.<br>
*><br>
*> If tau = 0, then H is taken to be the unit matrix<br>
*><br>
*> This version uses inline code if H has order < 11.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': form  H * C<br>
*>          = 'R': form  C * H<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension (M) if SIDE = 'L'<br>
*>                                        or (N) if SIDE = 'R'<br>
*>          The vector v in the representation of H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX<br>
*>          The value tau in the representation of H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array, dimension (LDC,N)<br>
*>          On entry, the m by n matrix C.<br>
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',<br>
*>          or C * H if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (N) if SIDE = 'L'<br>
*>                                            or (M) if SIDE = 'R'<br>
*>          WORK is not referenced if H has order < 11.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clarfx_(CHARACTER SIDE,INTEGER M,INTEGER N,float[] V,float[] TAU,float[] C,INTEGER LDC,float[] WORK);
/**
*> \brief \b CLARGV generates a vector of plane rotations with real cosines and complex sines.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clargv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clargv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clargv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARGV( N, X, INCX, Y, INCY, C, INCC )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCC, INCX, INCY, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * )<br>
*       COMPLEX            X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARGV generates a vector of complex plane rotations with real<br>
*> cosines, determined by elements of the complex vectors x and y.<br>
*> For i = 1,2,...,n<br>
*><br>
*>    (        c(i)   s(i) ) ( x(i) ) = ( r(i) )<br>
*>    ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )<br>
*><br>
*>    where c(i)**2 + ABS(s(i))**2 = 1<br>
*><br>
*> The following conventions are used (these are the same as in CLARTG,<br>
*> but differ from the BLAS1 routine CROTG):<br>
*>    If y(i)=0, then c(i)=1 and s(i)=0.<br>
*>    If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of plane rotations to be generated.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (1+(N-1)*INCX)<br>
*>          On entry, the vector x.<br>
*>          On exit, x(i) is overwritten by r(i), for i = 1,...,n.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between elements of X. INCX > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (1+(N-1)*INCY)<br>
*>          On entry, the vector y.<br>
*>          On exit, the sines of the plane rotations.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>          The increment between elements of Y. INCY > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (1+(N-1)*INCC)<br>
*>          The cosines of the plane rotations.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCC<br>
*> \verbatim<br>
*>          INCC is INTEGER<br>
*>          The increment between elements of C. INCC > 0.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  6-6-96 - Modified with a new algorithm by W. Kahan and J. Demmel<br>
*><br>
*>  This version has a few statements commented out for thread safety<br>
*>  (machine parameters are computed on each entry). 10 feb 03, SJH.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clargv_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] C,INTEGER INCC);
/**
*> \brief \b CLARNV returns a vector of random numbers from a uniform or normal distribution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARNV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarnv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarnv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarnv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARNV( IDIST, ISEED, N, X )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IDIST, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISEED( 4 )<br>
*       COMPLEX            X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARNV returns a vector of n random complex numbers from a uniform or<br>
*> normal distribution.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] IDIST<br>
*> \verbatim<br>
*>          IDIST is INTEGER<br>
*>          Specifies the distribution of the random numbers:<br>
*>          = 1:  real and imaginary parts each uniform (0,1)<br>
*>          = 2:  real and imaginary parts each uniform (-1,1)<br>
*>          = 3:  real and imaginary parts each normal (0,1)<br>
*>          = 4:  uniformly distributed on the disc abs(z) < 1<br>
*>          = 5:  uniformly distributed on the circle abs(z) = 1<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ISEED<br>
*> \verbatim<br>
*>          ISEED is INTEGER array, dimension (4)<br>
*>          On entry, the seed of the random number generator; the array<br>
*>          elements must be between 0 and 4095, and ISEED(4) must be<br>
*>          odd.<br>
*>          On exit, the seed is updated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of random numbers to be generated.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>          The generated random numbers.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  This routine calls the auxiliary routine SLARUV to generate random<br>
*>  real numbers from a uniform (0,1) distribution, in batches of up to<br>
*>  128 using vectorisable code. The Box-Muller method is used to<br>
*>  transform numbers from a uniform to a normal distribution.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clarnv_(INTEGER IDIST,int[] ISEED,INTEGER N,float[] X);
/**
*> \brief \b CLARRV computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenvalues of L D LT.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARRV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarrv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarrv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarrv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARRV( N, VL, VU, D, L, PIVMIN,<br>
*                          ISPLIT, M, DOL, DOU, MINRGP,<br>
*                          RTOL1, RTOL2, W, WERR, WGAP,<br>
*                          IBLOCK, INDEXW, GERS, Z, LDZ, ISUPPZ,<br>
*                          WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            DOL, DOU, INFO, LDZ, M, N<br>
*       REAL               MINRGP, PIVMIN, RTOL1, RTOL2, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IBLOCK( * ), INDEXW( * ), ISPLIT( * ),<br>
*      $                   ISUPPZ( * ), IWORK( * )<br>
*       REAL               D( * ), GERS( * ), L( * ), W( * ), WERR( * ),<br>
*      $                   WGAP( * ), WORK( * )<br>
*       COMPLEX           Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARRV computes the eigenvectors of the tridiagonal matrix<br>
*> T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.<br>
*> The input eigenvalues should have been computed by SLARRE.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*>          Lower bound of the interval that contains the desired<br>
*>          eigenvalues. VL < VU. Needed to compute gaps on the left or right<br>
*>          end of the extremal eigenvalues in the desired RANGE.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
*>          Upper bound of the interval that contains the desired<br>
*>          eigenvalues. VL < VU. Needed to compute gaps on the left or right<br>
*>          end of the extremal eigenvalues in the desired RANGE.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the N diagonal elements of the diagonal matrix D.<br>
*>          On exit, D may be overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] L<br>
*> \verbatim<br>
*>          L is REAL array, dimension (N)<br>
*>          On entry, the (N-1) subdiagonal elements of the unit<br>
*>          bidiagonal matrix L are in elements 1 to N-1 of L<br>
*>          (if the matrix is not split.) At the end of each block<br>
*>          is stored the corresponding shift as given by SLARRE.<br>
*>          On exit, L is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot allowed in the Sturm sequence.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ISPLIT<br>
*> \verbatim<br>
*>          ISPLIT is INTEGER array, dimension (N)<br>
*>          The splitting points, at which T breaks up into blocks.<br>
*>          The first block consists of rows/columns 1 to<br>
*>          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1<br>
*>          through ISPLIT( 2 ), etc.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The total number of input eigenvalues.  0 <= M <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DOL<br>
*> \verbatim<br>
*>          DOL is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] DOU<br>
*> \verbatim<br>
*>          DOU is INTEGER<br>
*>          If the user wants to compute only selected eigenvectors from all<br>
*>          the eigenvalues supplied, he can specify an index range DOL:DOU.<br>
*>          Or else the setting DOL=1, DOU=M should be applied.<br>
*>          Note that DOL and DOU refer to the order in which the eigenvalues<br>
*>          are stored in W.<br>
*>          If the user wants to compute only selected eigenpairs, then<br>
*>          the columns DOL-1 to DOU+1 of the eigenvector space Z contain the<br>
*>          computed eigenvectors. All other columns of Z are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MINRGP<br>
*> \verbatim<br>
*>          MINRGP is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] RTOL1<br>
*> \verbatim<br>
*>          RTOL1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] RTOL2<br>
*> \verbatim<br>
*>          RTOL2 is REAL<br>
*>           Parameters for bisection.<br>
*>           An interval [LEFT,RIGHT] has converged if<br>
*>           RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>          The first M elements of W contain the APPROXIMATE eigenvalues for<br>
*>          which eigenvectors are to be computed.  The eigenvalues<br>
*>          should be grouped by split-off block and ordered from<br>
*>          smallest to largest within the block ( The output array<br>
*>          W from SLARRE is expected here ). Furthermore, they are with<br>
*>          respect to the shift of the corresponding root representation<br>
*>          for their block. On exit, W holds the eigenvalues of the<br>
*>          UNshifted matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] WERR<br>
*> \verbatim<br>
*>          WERR is REAL array, dimension (N)<br>
*>          The first M elements contain the semiwidth of the uncertainty<br>
*>          interval of the corresponding eigenvalue in W<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] WGAP<br>
*> \verbatim<br>
*>          WGAP is REAL array, dimension (N)<br>
*>          The separation from the right neighbor eigenvalue in W.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IBLOCK<br>
*> \verbatim<br>
*>          IBLOCK is INTEGER array, dimension (N)<br>
*>          The indices of the blocks (submatrices) associated with the<br>
*>          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue<br>
*>          W(i) belongs to the first block from the top, =2 if W(i)<br>
*>          belongs to the second block, etc.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INDEXW<br>
*> \verbatim<br>
*>          INDEXW is INTEGER array, dimension (N)<br>
*>          The indices of the eigenvalues within each block (submatrix);<br>
*>          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the<br>
*>          i-th eigenvalue W(i) is the 10-th eigenvalue in the second block.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GERS<br>
*> \verbatim<br>
*>          GERS is REAL array, dimension (2*N)<br>
*>          The N Gerschgorin intervals (the i-th Gerschgorin interval<br>
*>          is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals should<br>
*>          be computed from the original UNshifted matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is array, dimension (LDZ, max(1,M) )<br>
*>          If INFO = 0, the first M columns of Z contain the<br>
*>          orthonormal eigenvectors of the matrix T<br>
*>          corresponding to the input eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISUPPZ<br>
*> \verbatim<br>
*>          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) )<br>
*>          The support of the eigenvectors in Z, i.e., the indices<br>
*>          indicating the nonzero elements in Z. The I-th eigenvector<br>
*>          is nonzero only in elements ISUPPZ( 2*I-1 ) through<br>
*>          ISUPPZ( 2*I ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (12*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (7*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*><br>
*>          > 0:  A problem occurred in CLARRV.<br>
*>          < 0:  One of the called subroutines signaled an internal problem.<br>
*>                Needs inspection of the corresponding parameter IINFO<br>
*>                for further information.<br>
*><br>
*>          =-1:  Problem in SLARRB when refining a child's eigenvalues.<br>
*>          =-2:  Problem in SLARRF when computing the RRR of a child.<br>
*>                When a child is inside a tight cluster, it can be difficult<br>
*>                to find an RRR. A partial remedy from the user's point of<br>
*>                view is to make the parameter MINRGP smaller and recompile.<br>
*>                However, as the orthogonality of the computed vectors is<br>
*>                proportional to 1/MINRGP, the user should be aware that<br>
*>                he might be trading in precision when he decreases MINRGP.<br>
*>          =-3:  Problem in SLARRB when refining a single eigenvalue<br>
*>                after the Rayleigh correction was rejected.<br>
*>          = 5:  The Rayleigh Quotient Iteration failed to converge to<br>
*>                full accuracy in MAXITR steps.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Beresford Parlett, University of California, Berkeley, USA \n<br>
*> Jim Demmel, University of California, Berkeley, USA \n<br>
*> Inderjit Dhillon, University of Texas, Austin, USA \n<br>
*> Osni Marques, LBNL/NERSC, USA \n<br>
*> Christof Voemel, University of California, Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void clarrv_(INTEGER N,REAL VL,REAL VU,float[] D,float[] L,REAL PIVMIN,int[] ISPLIT,INTEGER M,INTEGER DOL,INTEGER DOU,REAL MINRGP,REAL RTOL1,REAL RTOL2,float[] W,float[] WERR,float[] WGAP,int[] IBLOCK,int[] INDEXW,float[] GERS,float[] Z,INTEGER LDZ,int[] ISUPPZ,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b CLARSCL2 performs reciprocal diagonal scaling on a vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARSCL2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarscl2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarscl2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarscl2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARSCL2 ( M, N, D, X, LDX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDX<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            X( LDX, * )<br>
*       REAL               D( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARSCL2 performs a reciprocal diagonal scaling on an vector:<br>
*>   x <-- inv(D) * x<br>
*> where the REAL diagonal matrix D is stored as a vector.<br>
*><br>
*> Eventually to be replaced by BLAS_cge_diag_scale in the new BLAS<br>
*> standard.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>     The number of rows of D and X. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of columns of X. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, length M<br>
*>     Diagonal matrix D, stored as a vector of length M.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (LDX,N)<br>
*>     On entry, the vector X to be scaled by D.<br>
*>     On exit, the scaled vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>     The leading dimension of the vector X. LDX >= M.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void clarscl2_(INTEGER M,INTEGER N,float[] D,float[] X,INTEGER LDX);
/**
*> \brief \b CLARTG generates a plane rotation with real cosine and complex sine.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARTG + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clartg.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clartg.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clartg.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARTG( F, G, CS, SN, R )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               CS<br>
*       COMPLEX            F, G, R, SN<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARTG generates a plane rotation so that<br>
*><br>
*>    [  CS  SN  ]     [ F ]     [ R ]<br>
*>    [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.<br>
*>    [ -SN  CS  ]     [ G ]     [ 0 ]<br>
*><br>
*> This is a faster version of the BLAS1 routine CROTG, except for<br>
*> the following differences:<br>
*>    F and G are unchanged on return.<br>
*>    If G=0, then CS=1 and SN=0.<br>
*>    If F=0, then CS=0 and SN is chosen so that R is real.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] F<br>
*> \verbatim<br>
*>          F is COMPLEX<br>
*>          The first component of vector to be rotated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] G<br>
*> \verbatim<br>
*>          G is COMPLEX<br>
*>          The second component of vector to be rotated.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CS<br>
*> \verbatim<br>
*>          CS is REAL<br>
*>          The cosine of the rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SN<br>
*> \verbatim<br>
*>          SN is COMPLEX<br>
*>          The sine of the rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] R<br>
*> \verbatim<br>
*>          R is COMPLEX<br>
*>          The nonzero component of the rotated vector.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel<br>
*><br>
*>  This version has a few statements commented out for thread safety<br>
*>  (machine parameters are computed on each entry). 10 feb 03, SJH.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clartg_(float[] F,float[] G,REAL CS,float[] SN,float[] R);
/**
*> \brief \b CLARTV applies a vector of plane rotations with real cosines and complex sines to the elements of a pair of vectors.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARTV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clartv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clartv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clartv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARTV( N, X, INCX, Y, INCY, C, S, INCC )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCC, INCX, INCY, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * )<br>
*       COMPLEX            S( * ), X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARTV applies a vector of complex plane rotations with real cosines<br>
*> to elements of the complex vectors x and y. For i = 1,2,...,n<br>
*><br>
*>    ( x(i) ) := (        c(i)   s(i) ) ( x(i) )<br>
*>    ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of plane rotations to be applied.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (1+(N-1)*INCX)<br>
*>          The vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between elements of X. INCX > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (1+(N-1)*INCY)<br>
*>          The vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>          The increment between elements of Y. INCY > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (1+(N-1)*INCC)<br>
*>          The cosines of the plane rotations.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is COMPLEX array, dimension (1+(N-1)*INCC)<br>
*>          The sines of the plane rotations.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCC<br>
*> \verbatim<br>
*>          INCC is INTEGER<br>
*>          The increment between elements of C and S. INCC > 0.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clartv_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] C,float[] S,INTEGER INCC);
/**
*> \brief \b CLARZ applies an elementary reflector (as returned by stzrzf) to a general matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE<br>
*       INTEGER            INCV, L, LDC, M, N<br>
*       COMPLEX            TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            C( LDC, * ), V( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARZ applies a complex elementary reflector H to a complex<br>
*> M-by-N matrix C, from either the left or the right. H is represented<br>
*> in the form<br>
*><br>
*>       H = I - tau * v * v**H<br>
*><br>
*> where tau is a complex scalar and v is a complex vector.<br>
*><br>
*> If tau = 0, then H is taken to be the unit matrix.<br>
*><br>
*> To apply H**H (the conjugate transpose of H), supply conjg(tau) instead<br>
*> tau.<br>
*><br>
*> H is a product of k elementary reflectors as returned by CTZRZF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': form  H * C<br>
*>          = 'R': form  C * H<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The number of entries of the vector V containing<br>
*>          the meaningful part of the Householder vectors.<br>
*>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension (1+(L-1)*abs(INCV))<br>
*>          The vector v in the representation of H as returned by<br>
*>          CTZRZF. V is not used if TAU = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCV<br>
*> \verbatim<br>
*>          INCV is INTEGER<br>
*>          The increment between elements of v. INCV <> 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX<br>
*>          The value tau in the representation of H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',<br>
*>          or C * H if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDC >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension<br>
*>                         (N) if SIDE = 'L'<br>
*>                      or (M) if SIDE = 'R'<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clarz_(CHARACTER SIDE,INTEGER M,INTEGER N,INTEGER L,float[] V,INTEGER INCV,float[] TAU,float[] C,INTEGER LDC,float[] WORK);
/**
*> \brief \b CLARZB applies a block reflector or its conjugate-transpose to a general matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARZB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarzb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarzb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarzb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V,<br>
*                          LDV, T, LDT, C, LDC, WORK, LDWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, SIDE, STOREV, TRANS<br>
*       INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ),<br>
*      $                   WORK( LDWORK, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARZB applies a complex block reflector H or its transpose H**H<br>
*> to a complex distributed M-by-N  C from the left or the right.<br>
*><br>
*> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply H or H**H from the Left<br>
*>          = 'R': apply H or H**H from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply H (No transpose)<br>
*>          = 'C': apply H**H (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIRECT<br>
*> \verbatim<br>
*>          DIRECT is CHARACTER*1<br>
*>          Indicates how H is formed from a product of elementary<br>
*>          reflectors<br>
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)<br>
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)<br>
*> \endverbatim<br>
*><br>
*> \param[in] STOREV<br>
*> \verbatim<br>
*>          STOREV is CHARACTER*1<br>
*>          Indicates how the vectors which define the elementary<br>
*>          reflectors are stored:<br>
*>          = 'C': Columnwise                        (not supported yet)<br>
*>          = 'R': Rowwise<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The order of the matrix T (= the number of elementary<br>
*>          reflectors whose product defines the block reflector).<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The number of columns of the matrix V containing the<br>
*>          meaningful part of the Householder reflectors.<br>
*>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension (LDV,NV).<br>
*>          If STOREV = 'C', NV = K; if STOREV = 'R', NV = L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V.<br>
*>          If STOREV = 'C', LDV >= L; if STOREV = 'R', LDV >= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] T<br>
*> \verbatim<br>
*>          T is COMPLEX array, dimension (LDT,K)<br>
*>          The triangular K-by-K matrix T in the representation of the<br>
*>          block reflector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. LDT >= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDC >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (LDWORK,K)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDWORK<br>
*> \verbatim<br>
*>          LDWORK is INTEGER<br>
*>          The leading dimension of the array WORK.<br>
*>          If SIDE = 'L', LDWORK >= max(1,N);<br>
*>          if SIDE = 'R', LDWORK >= max(1,M).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clarzb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,INTEGER L,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] C,INTEGER LDC,float[] WORK,INTEGER LDWORK);
/**
*> \brief \b CLARZT forms the triangular factor T of a block reflector H = I - vtvH.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLARZT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarzt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarzt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarzt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, STOREV<br>
*       INTEGER            K, LDT, LDV, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLARZT forms the triangular factor T of a complex block reflector<br>
*> H of order > n, which is defined as a product of k elementary<br>
*> reflectors.<br>
*><br>
*> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;<br>
*><br>
*> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.<br>
*><br>
*> If STOREV = 'C', the vector which defines the elementary reflector<br>
*> H(i) is stored in the i-th column of the array V, and<br>
*><br>
*>    H  =  I - V * T * V**H<br>
*><br>
*> If STOREV = 'R', the vector which defines the elementary reflector<br>
*> H(i) is stored in the i-th row of the array V, and<br>
*><br>
*>    H  =  I - V**H * T * V<br>
*><br>
*> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] DIRECT<br>
*> \verbatim<br>
*>          DIRECT is CHARACTER*1<br>
*>          Specifies the order in which the elementary reflectors are<br>
*>          multiplied to form the block reflector:<br>
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)<br>
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)<br>
*> \endverbatim<br>
*><br>
*> \param[in] STOREV<br>
*> \verbatim<br>
*>          STOREV is CHARACTER*1<br>
*>          Specifies how the vectors which define the elementary<br>
*>          reflectors are stored (see also Further Details):<br>
*>          = 'C': columnwise                        (not supported yet)<br>
*>          = 'R': rowwise<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the block reflector H. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The order of the triangular factor T (= the number of<br>
*>          elementary reflectors). K >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V<br>
*> \verbatim<br>
*>          V is COMPLEX array, dimension<br>
*>                               (LDV,K) if STOREV = 'C'<br>
*>                               (LDV,N) if STOREV = 'R'<br>
*>          The matrix V. See further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V.<br>
*>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is COMPLEX array, dimension (LDT,K)<br>
*>          The k by k triangular factor T of the block reflector.<br>
*>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is<br>
*>          lower triangular. The rest of the array is not used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. LDT >= K.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The shape of the matrix V and the storage of the vectors which define<br>
*>  the H(i) is best illustrated by the following example with n = 5 and<br>
*>  k = 3. The elements equal to 1 are not stored; the corresponding<br>
*>  array elements are modified but restored on exit. The rest of the<br>
*>  array is not used.<br>
*><br>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':<br>
*><br>
*>                                              ______V_____<br>
*>         ( v1 v2 v3 )                        /            \<br>
*>         ( v1 v2 v3 )                      ( v1 v1 v1 v1 v1 . . . . 1 )<br>
*>     V = ( v1 v2 v3 )                      ( v2 v2 v2 v2 v2 . . . 1   )<br>
*>         ( v1 v2 v3 )                      ( v3 v3 v3 v3 v3 . . 1     )<br>
*>         ( v1 v2 v3 )<br>
*>            .  .  .<br>
*>            .  .  .<br>
*>            1  .  .<br>
*>               1  .<br>
*>                  1<br>
*><br>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':<br>
*><br>
*>                                                        ______V_____<br>
*>            1                                          /            \<br>
*>            .  1                           ( 1 . . . . v1 v1 v1 v1 v1 )<br>
*>            .  .  1                        ( . 1 . . . v2 v2 v2 v2 v2 )<br>
*>            .  .  .                        ( . . 1 . . v3 v3 v3 v3 v3 )<br>
*>            .  .  .<br>
*>         ( v1 v2 v3 )<br>
*>         ( v1 v2 v3 )<br>
*>     V = ( v1 v2 v3 )<br>
*>         ( v1 v2 v3 )<br>
*>         ( v1 v2 v3 )<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clarzt_(CHARACTER DIRECT,CHARACTER STOREV,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] TAU,float[] T,INTEGER LDT);
/**
*> \brief \b CLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLASCL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clascl.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clascl.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clascl.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TYPE<br>
*       INTEGER            INFO, KL, KU, LDA, M, N<br>
*       REAL               CFROM, CTO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLASCL multiplies the M by N complex matrix A by the real scalar<br>
*> CTO/CFROM.  This is done without over/underflow as long as the final<br>
*> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that<br>
*> A may be full, upper triangular, lower triangular, upper Hessenberg,<br>
*> or banded.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TYPE<br>
*> \verbatim<br>
*>          TYPE is CHARACTER*1<br>
*>          TYPE indices the storage type of the input matrix.<br>
*>          = 'G':  A is a full matrix.<br>
*>          = 'L':  A is a lower triangular matrix.<br>
*>          = 'U':  A is an upper triangular matrix.<br>
*>          = 'H':  A is an upper Hessenberg matrix.<br>
*>          = 'B':  A is a symmetric band matrix with lower bandwidth KL<br>
*>                  and upper bandwidth KU and with the only the lower<br>
*>                  half stored.<br>
*>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL<br>
*>                  and upper bandwidth KU and with the only the upper<br>
*>                  half stored.<br>
*>          = 'Z':  A is a band matrix with lower bandwidth KL and upper<br>
*>                  bandwidth KU. See CGBTRF for storage details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>          The lower bandwidth of A.  Referenced only if TYPE = 'B',<br>
*>          'Q' or 'Z'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>          The upper bandwidth of A.  Referenced only if TYPE = 'B',<br>
*>          'Q' or 'Z'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CFROM<br>
*> \verbatim<br>
*>          CFROM is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] CTO<br>
*> \verbatim<br>
*>          CTO is REAL<br>
*><br>
*>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed<br>
*>          without over/underflow if the final result CTO*A(I,J)/CFROM<br>
*>          can be represented without over/underflow.  CFROM must be<br>
*>          nonzero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the<br>
*>          storage type.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*>          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M);<br>
*>             TYPE = 'B', LDA >= KL+1;<br>
*>             TYPE = 'Q', LDA >= KU+1;<br>
*>             TYPE = 'Z', LDA >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          0  - successful exit<br>
*>          <0 - if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clascl_(CHARACTER TYPE,INTEGER KL,INTEGER KU,REAL CFROM,REAL CTO,INTEGER M,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b CLASCL2 performs diagonal scaling on a vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLASCL2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clascl2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clascl2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clascl2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLASCL2 ( M, N, D, X, LDX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDX<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * )<br>
*       COMPLEX            X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLASCL2 performs a diagonal scaling on a vector:<br>
*>   x <-- D * x<br>
*> where the diagonal REAL matrix D is stored as a vector.<br>
*><br>
*> Eventually to be replaced by BLAS_cge_diag_scale in the new BLAS<br>
*> standard.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>     The number of rows of D and X. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of columns of X. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, length M<br>
*>     Diagonal matrix D, stored as a vector of length M.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (LDX,N)<br>
*>     On entry, the vector X to be scaled by D.<br>
*>     On exit, the scaled vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>     The leading dimension of the vector X. LDX >= M.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void clascl2_(INTEGER M,INTEGER N,float[] D,float[] X,INTEGER LDX);
/**
*> \brief \b CLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLASET + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claset.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claset.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claset.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLASET( UPLO, M, N, ALPHA, BETA, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            LDA, M, N<br>
*       COMPLEX            ALPHA, BETA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLASET initializes a 2-D array A to BETA on the diagonal and<br>
*> ALPHA on the offdiagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies the part of the matrix A to be set.<br>
*>          = 'U':      Upper triangular part is set. The lower triangle<br>
*>                      is unchanged.<br>
*>          = 'L':      Lower triangular part is set. The upper triangle<br>
*>                      is unchanged.<br>
*>          Otherwise:  All of the matrix A is set.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          On entry, M specifies the number of rows of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          On entry, N specifies the number of columns of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>          All the offdiagonal array elements are set to ALPHA.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>          All the diagonal array elements are set to BETA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the m by n matrix A.<br>
*>          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;<br>
*>                   A(i,i) = BETA , 1 <= i <= min(m,n)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2015<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void claset_(CHARACTER UPLO,INTEGER M,INTEGER N,float[] ALPHA,float[] BETA,float[] A,INTEGER LDA);
/**
*> \brief \b CLASR applies a sequence of plane rotations to a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLASR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, PIVOT, SIDE<br>
*       INTEGER            LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * ), S( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLASR applies a sequence of real plane rotations to a complex matrix<br>
*> A, from either the left or the right.<br>
*><br>
*> When SIDE = 'L', the transformation takes the form<br>
*><br>
*>    A := P*A<br>
*><br>
*> and when SIDE = 'R', the transformation takes the form<br>
*><br>
*>    A := A*P**T<br>
*><br>
*> where P is an orthogonal matrix consisting of a sequence of z plane<br>
*> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',<br>
*> and P**T is the transpose of P.<br>
*> <br>
*> When DIRECT = 'F' (Forward sequence), then<br>
*> <br>
*>    P = P(z-1) * ... * P(2) * P(1)<br>
*> <br>
*> and when DIRECT = 'B' (Backward sequence), then<br>
*> <br>
*>    P = P(1) * P(2) * ... * P(z-1)<br>
*> <br>
*> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation<br>
*> <br>
*>    R(k) = (  c(k)  s(k) )<br>
*>         = ( -s(k)  c(k) ).<br>
*> <br>
*> When PIVOT = 'V' (Variable pivot), the rotation is performed<br>
*> for the plane (k,k+1), i.e., P(k) has the form<br>
*> <br>
*>    P(k) = (  1                                            )<br>
*>           (       ...                                     )<br>
*>           (              1                                )<br>
*>           (                   c(k)  s(k)                  )<br>
*>           (                  -s(k)  c(k)                  )<br>
*>           (                                1              )<br>
*>           (                                     ...       )<br>
*>           (                                            1  )<br>
*> <br>
*> where R(k) appears as a rank-2 modification to the identity matrix in<br>
*> rows and columns k and k+1.<br>
*> <br>
*> When PIVOT = 'T' (Top pivot), the rotation is performed for the<br>
*> plane (1,k+1), so P(k) has the form<br>
*> <br>
*>    P(k) = (  c(k)                    s(k)                 )<br>
*>           (         1                                     )<br>
*>           (              ...                              )<br>
*>           (                     1                         )<br>
*>           ( -s(k)                    c(k)                 )<br>
*>           (                                 1             )<br>
*>           (                                      ...      )<br>
*>           (                                             1 )<br>
*> <br>
*> where R(k) appears in rows and columns 1 and k+1.<br>
*> <br>
*> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is<br>
*> performed for the plane (k,z), giving P(k) the form<br>
*> <br>
*>    P(k) = ( 1                                             )<br>
*>           (      ...                                      )<br>
*>           (             1                                 )<br>
*>           (                  c(k)                    s(k) )<br>
*>           (                         1                     )<br>
*>           (                              ...              )<br>
*>           (                                     1         )<br>
*>           (                 -s(k)                    c(k) )<br>
*> <br>
*> where R(k) appears in rows and columns k and z.  The rotations are<br>
*> performed without ever forming P(k) explicitly.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          Specifies whether the plane rotation matrix P is applied to<br>
*>          A on the left or the right.<br>
*>          = 'L':  Left, compute A := P*A<br>
*>          = 'R':  Right, compute A:= A*P**T<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVOT<br>
*> \verbatim<br>
*>          PIVOT is CHARACTER*1<br>
*>          Specifies the plane for which P(k) is a plane rotation<br>
*>          matrix.<br>
*>          = 'V':  Variable pivot, the plane (k,k+1)<br>
*>          = 'T':  Top pivot, the plane (1,k+1)<br>
*>          = 'B':  Bottom pivot, the plane (k,z)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIRECT<br>
*> \verbatim<br>
*>          DIRECT is CHARACTER*1<br>
*>          Specifies whether P is a forward or backward sequence of<br>
*>          plane rotations.<br>
*>          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)<br>
*>          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  If m <= 1, an immediate<br>
*>          return is effected.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  If n <= 1, an<br>
*>          immediate return is effected.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension<br>
*>                  (M-1) if SIDE = 'L'<br>
*>                  (N-1) if SIDE = 'R'<br>
*>          The cosines c(k) of the plane rotations.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension<br>
*>                  (M-1) if SIDE = 'L'<br>
*>                  (N-1) if SIDE = 'R'<br>
*>          The sines s(k) of the plane rotations.  The 2-by-2 plane<br>
*>          rotation part of the matrix P(k), R(k), has the form<br>
*>          R(k) = (  c(k)  s(k) )<br>
*>                 ( -s(k)  c(k) ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The M-by-N matrix A.  On exit, A is overwritten by P*A if<br>
*>          SIDE = 'R' or by A*P**T if SIDE = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clasr_(CHARACTER SIDE,CHARACTER PIVOT,CHARACTER DIRECT,INTEGER M,INTEGER N,float[] C,float[] S,float[] A,INTEGER LDA);
/**
*> \brief \b CLASSQ updates a sum of squares represented in scaled form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLASSQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/classq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/classq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/classq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       REAL               SCALE, SUMSQ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLASSQ returns the values scl and ssq such that<br>
*><br>
*>    ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,<br>
*><br>
*> where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is<br>
*> assumed to be at least unity and the value of ssq will then satisfy<br>
*><br>
*>    1.0 .le. ssq .le. ( sumsq + 2*n ).<br>
*><br>
*> scale is assumed to be non-negative and scl returns the value<br>
*><br>
*>    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),<br>
*>           i<br>
*><br>
*> scale and sumsq must be supplied in SCALE and SUMSQ respectively.<br>
*> SCALE and SUMSQ are overwritten by scl and ssq respectively.<br>
*><br>
*> The routine makes only one pass through the vector X.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of elements to be used from the vector X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>          The vector x as described above.<br>
*>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between successive values of the vector X.<br>
*>          INCX > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          On entry, the value  scale  in the equation above.<br>
*>          On exit, SCALE is overwritten with the value  scl .<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SUMSQ<br>
*> \verbatim<br>
*>          SUMSQ is REAL<br>
*>          On entry, the value  sumsq  in the equation above.<br>
*>          On exit, SUMSQ is overwritten with the value  ssq .<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void classq_(INTEGER N,float[] X,INTEGER INCX,REAL SCALE,REAL SUMSQ);
/**
*> \brief \b CLASWP performs a series of row interchanges on a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLASWP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claswp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claswp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claswp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLASWP( N, A, LDA, K1, K2, IPIV, INCX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, K1, K2, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLASWP performs a series of row interchanges on the matrix A.<br>
*> One row interchange is initiated for each of rows K1 through K2 of A.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the matrix of column dimension N to which the row<br>
*>          interchanges will be applied.<br>
*>          On exit, the permuted matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K1<br>
*> \verbatim<br>
*>          K1 is INTEGER<br>
*>          The first element of IPIV for which a row interchange will<br>
*>          be done.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K2<br>
*> \verbatim<br>
*>          K2 is INTEGER<br>
*>          The last element of IPIV for which a row interchange will<br>
*>          be done.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (K2*abs(INCX))<br>
*>          The vector of pivot indices.  Only the elements in positions<br>
*>          K1 through K2 of IPIV are accessed.<br>
*>          IPIV(K) = L implies rows K and L are to be interchanged.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between successive values of IPIV.  If IPIV<br>
*>          is negative, the pivots are applied in reverse order.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Modified by<br>
*>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void claswp_(INTEGER N,float[] A,INTEGER LDA,INTEGER K1,INTEGER K2,int[] IPIV,INTEGER INCX);
/**
*> \brief \b CLASYF computes a partial factorization of a complex symmetric matrix using the Bunch-Kaufman diagonal pivoting method.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download CLASYF + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasyf.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasyf.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasyf.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KB, LDA, LDW, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), W( LDW, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLASYF computes a partial factorization of a complex symmetric matrix<br>
*> A using the Bunch-Kaufman diagonal pivoting method. The partial<br>
*> factorization has the form:<br>
*><br>
*> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:<br>
*>       ( 0  U22 ) (  0   D  ) ( U12**T U22**T )<br>
*><br>
*> A  =  ( L11  0 ) ( D    0  ) ( L11**T L21**T )  if UPLO = 'L'<br>
*>       ( L21  I ) ( 0   A22 ) (  0       I    )<br>
*><br>
*> where the order of D is at most NB. The actual order is returned in<br>
*> the argument KB, and is either NB or NB-1, or N if N <= NB.<br>
*> Note that U**T denotes the transpose of U.<br>
*><br>
*> CLASYF is an auxiliary routine called by CSYTRF. It uses blocked code<br>
*> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or<br>
*> A22 (if UPLO = 'L').<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is stored:<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The maximum number of columns of the matrix A that should be<br>
*>          factored.  NB should be at least 2 to allow for 2-by-2 pivot<br>
*>          blocks.<br>
*> \endverbatim<br>
*><br>
*> \param[out] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of columns of A that were actually factored.<br>
*>          KB is either NB-1 or NB, or N if N <= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit, A contains details of the partial factorization.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D.<br>
*><br>
*>          If UPLO = 'U':<br>
*>             Only the last KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns<br>
*>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)<br>
*>             is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>             Only the first KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns<br>
*>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)<br>
*>             is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (LDW,NB)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDW<br>
*> \verbatim<br>
*>          LDW is INTEGER<br>
*>          The leading dimension of the array W.  LDW >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complexSYcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>  November 2013,  Igor Kozachenko,<br>
*>                  Computer Science Division,<br>
*>                  University of California, Berkeley<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void clasyf_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] IPIV,float[] W,INTEGER LDW,INTEGER INFO);
/**
*> \brief \b CLASYF_ROOK computes a partial factorization of a complex symmetric matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download CLASYF_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasyf_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasyf_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasyf_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLASYF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KB, LDA, LDW, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), W( LDW, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLASYF_ROOK computes a partial factorization of a complex symmetric<br>
*> matrix A using the bounded Bunch-Kaufman ("rook") diagonal<br>
*> pivoting method. The partial factorization has the form:<br>
*><br>
*> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:<br>
*>       ( 0  U22 ) (  0   D  ) ( U12**T U22**T )<br>
*><br>
*> A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L'<br>
*>       ( L21  I ) (  0  A22 ) (  0       I    )<br>
*><br>
*> where the order of D is at most NB. The actual order is returned in<br>
*> the argument KB, and is either NB or NB-1, or N if N <= NB.<br>
*><br>
*> CLASYF_ROOK is an auxiliary routine called by CSYTRF_ROOK. It uses<br>
*> blocked code (calling Level 3 BLAS) to update the submatrix<br>
*> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is stored:<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The maximum number of columns of the matrix A that should be<br>
*>          factored.  NB should be at least 2 to allow for 2-by-2 pivot<br>
*>          blocks.<br>
*> \endverbatim<br>
*><br>
*> \param[out] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of columns of A that were actually factored.<br>
*>          KB is either NB-1 or NB, or N if N <= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit, A contains details of the partial factorization.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D.<br>
*><br>
*>          If UPLO = 'U':<br>
*>             Only the last KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k-1 and -IPIV(k-1) were inerchaged,<br>
*>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>             Only the first KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>             were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k+1 and -IPIV(k+1) were inerchaged,<br>
*>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (LDW,NB)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDW<br>
*> \verbatim<br>
*>          LDW is INTEGER<br>
*>          The leading dimension of the array W.  LDW >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complexSYcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>  November 2013,     Igor Kozachenko,<br>
*>                  Computer Science Division,<br>
*>                  University of California, Berkeley<br>
*><br>
*>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,<br>
*>                  School of Mathematics,<br>
*>                  University of Manchester<br>
*><br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void clasyf_rook_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] IPIV,float[] W,INTEGER LDW,INTEGER INFO);
/**
*> \brief \b CLATBS solves a triangular banded system of equations.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLATBS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatbs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatbs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatbs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X,<br>
*                          SCALE, CNORM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
*       REAL               SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               CNORM( * )<br>
*       COMPLEX            AB( LDAB, * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLATBS solves one of the triangular systems<br>
*><br>
*>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,<br>
*><br>
*> with scaling to prevent overflow, where A is an upper or lower<br>
*> triangular band matrix.  Here A**T denotes the transpose of A, x and b<br>
*> are n-element vectors, and s is a scaling factor, usually less than<br>
*> or equal to 1, chosen so that the components of x will be less than<br>
*> the overflow threshold.  If the unscaled problem will not cause<br>
*> overflow, the Level 2 BLAS routine CTBSV is called.  If the matrix A<br>
*> is singular (A(j,j) = 0 for some j), then s is set to 0 and a<br>
*> non-trivial solution to A*x = 0 is returned.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the matrix A is upper or lower triangular.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the operation applied to A.<br>
*>          = 'N':  Solve A * x = s*b     (No transpose)<br>
*>          = 'T':  Solve A**T * x = s*b  (Transpose)<br>
*>          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          Specifies whether or not the matrix A is unit triangular.<br>
*>          = 'N':  Non-unit triangular<br>
*>          = 'U':  Unit triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] NORMIN<br>
*> \verbatim<br>
*>          NORMIN is CHARACTER*1<br>
*>          Specifies whether CNORM has been set or not.<br>
*>          = 'Y':  CNORM contains the column norms on entry<br>
*>          = 'N':  CNORM is not set on entry.  On exit, the norms will<br>
*>                  be computed and stored in CNORM.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KD<br>
*> \verbatim<br>
*>          KD is INTEGER<br>
*>          The number of subdiagonals or superdiagonals in the<br>
*>          triangular matrix A.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The upper or lower triangular band matrix A, stored in the<br>
*>          first KD+1 rows of the array. The j-th column of A is stored<br>
*>          in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>          On entry, the right hand side b of the triangular system.<br>
*>          On exit, X is overwritten by the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          The scaling factor s for the triangular system<br>
*>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.<br>
*>          If SCALE = 0, the matrix A is singular or badly scaled, and<br>
*>          the vector x is an exact or approximate solution to A*x = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CNORM<br>
*> \verbatim<br>
*>          CNORM is REAL array, dimension (N)<br>
*><br>
*>          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)<br>
*>          contains the norm of the off-diagonal part of the j-th column<br>
*>          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal<br>
*>          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)<br>
*>          must be greater than or equal to the 1-norm.<br>
*><br>
*>          If NORMIN = 'N', CNORM is an output argument and CNORM(j)<br>
*>          returns the 1-norm of the offdiagonal part of the j-th column<br>
*>          of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -k, the k-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  A rough bound on x is computed; if that is less than overflow, CTBSV<br>
*>  is called, otherwise, specific code is used which checks for possible<br>
*>  overflow or divide-by-zero at every operation.<br>
*><br>
*>  A columnwise scheme is used for solving A*x = b.  The basic algorithm<br>
*>  if A is lower triangular is<br>
*><br>
*>       x[1:n] := b[1:n]<br>
*>       for j = 1, ..., n<br>
*>            x(j) := x(j) / A(j,j)<br>
*>            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]<br>
*>       end<br>
*><br>
*>  Define bounds on the components of x after j iterations of the loop:<br>
*>     M(j) = bound on x[1:j]<br>
*>     G(j) = bound on x[j+1:n]<br>
*>  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.<br>
*><br>
*>  Then for iteration j+1 we have<br>
*>     M(j+1) <= G(j) / | A(j+1,j+1) |<br>
*>     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |<br>
*>            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )<br>
*><br>
*>  where CNORM(j+1) is greater than or equal to the infinity-norm of<br>
*>  column j+1 of A, not counting the diagonal.  Hence<br>
*><br>
*>     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )<br>
*>                  1<=i<=j<br>
*>  and<br>
*><br>
*>     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )<br>
*>                                   1<=i< j<br>
*><br>
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTBSV if the<br>
*>  reciprocal of the largest M(j), j=1,..,n, is larger than<br>
*>  max(underflow, 1/overflow).<br>
*><br>
*>  The bound on x(j) is also used to determine when a step in the<br>
*>  columnwise method can be performed without fear of overflow.  If<br>
*>  the computed bound is greater than a large constant, x is scaled to<br>
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to<br>
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.<br>
*><br>
*>  Similarly, a row-wise scheme is used to solve A**T *x = b  or<br>
*>  A**H *x = b.  The basic algorithm for A upper triangular is<br>
*><br>
*>       for j = 1, ..., n<br>
*>            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)<br>
*>       end<br>
*><br>
*>  We simultaneously compute two bounds<br>
*>       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j<br>
*>       M(j) = bound on x(i), 1<=i<=j<br>
*><br>
*>  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we<br>
*>  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.<br>
*>  Then the bound on x(j) is<br>
*><br>
*>       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |<br>
*><br>
*>            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )<br>
*>                      1<=i<=j<br>
*><br>
*>  and we can safely call CTBSV if 1/M(n) and 1/G(n) are both greater<br>
*>  than max(underflow, 1/overflow).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clatbs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
/**
*> \brief \b CLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contribution to the reciprocal Dif-estimate.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLATDF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatdf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatdf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatdf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV,<br>
*                          JPIV )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IJOB, LDZ, N<br>
*       REAL               RDSCAL, RDSUM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), JPIV( * )<br>
*       COMPLEX            RHS( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLATDF computes the contribution to the reciprocal Dif-estimate<br>
*> by solving for x in Z * x = b, where b is chosen such that the norm<br>
*> of x is as large as possible. It is assumed that LU decomposition<br>
*> of Z has been computed by CGETC2. On entry RHS = f holds the<br>
*> contribution from earlier solved sub-systems, and on return RHS = x.<br>
*><br>
*> The factorization of Z returned by CGETC2 has the form<br>
*> Z = P * L * U * Q, where P and Q are permutation matrices. L is lower<br>
*> triangular with unit diagonal elements and U is upper triangular.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] IJOB<br>
*> \verbatim<br>
*>          IJOB is INTEGER<br>
*>          IJOB = 2: First compute an approximative null-vector e<br>
*>              of Z using CGECON, e is normalized and solve for<br>
*>              Zx = +-e - f with the sign giving the greater value of<br>
*>              2-norm(x).  About 5 times as expensive as Default.<br>
*>          IJOB .ne. 2: Local look ahead strategy where<br>
*>              all entries of the r.h.s. b is chosen as either +1 or<br>
*>              -1.  Default.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (LDZ, N)<br>
*>          On entry, the LU part of the factorization of the n-by-n<br>
*>          matrix Z computed by CGETC2:  Z = P * L * U * Q<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDA >= max(1, N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RHS<br>
*> \verbatim<br>
*>          RHS is COMPLEX array, dimension (N).<br>
*>          On entry, RHS contains contributions from other subsystems.<br>
*>          On exit, RHS contains the solution of the subsystem with<br>
*>          entries according to the value of IJOB (see above).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RDSUM<br>
*> \verbatim<br>
*>          RDSUM is REAL<br>
*>          On entry, the sum of squares of computed contributions to<br>
*>          the Dif-estimate under computation by CTGSYL, where the<br>
*>          scaling factor RDSCAL (see below) has been factored out.<br>
*>          On exit, the corresponding sum of squares updated with the<br>
*>          contributions from the current sub-system.<br>
*>          If TRANS = 'T' RDSUM is not touched.<br>
*>          NOTE: RDSUM only makes sense when CTGSY2 is called by CTGSYL.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RDSCAL<br>
*> \verbatim<br>
*>          RDSCAL is REAL<br>
*>          On entry, scaling factor used to prevent overflow in RDSUM.<br>
*>          On exit, RDSCAL is updated w.r.t. the current contributions<br>
*>          in RDSUM.<br>
*>          If TRANS = 'T', RDSCAL is not touched.<br>
*>          NOTE: RDSCAL only makes sense when CTGSY2 is called by<br>
*>          CTGSYL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N).<br>
*>          The pivot indices; for 1 <= i <= N, row i of the<br>
*>          matrix has been interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] JPIV<br>
*> \verbatim<br>
*>          JPIV is INTEGER array, dimension (N).<br>
*>          The pivot indices; for 1 <= j <= N, column j of the<br>
*>          matrix has been interchanged with column JPIV(j).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*>  This routine is a further developed implementation of algorithm<br>
*>  BSOLVE in [1] using complete pivoting in the LU factorization.<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,<br>
*>     Umea University, S-901 87 Umea, Sweden.<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>   [1]   Bo Kagstrom and Lars Westin,<br>
*>         Generalized Schur Methods with Condition Estimators for<br>
*>         Solving the Generalized Sylvester Equation, IEEE Transactions<br>
*>         on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751.<br>
*><br>
*>   [2]   Peter Poromaa,<br>
*>         On Efficient and Robust Estimators for the Separation<br>
*>         between two Regular Matrix Pairs with Applications in<br>
*>         Condition Estimation. Report UMINF-95.05, Department of<br>
*>         Computing Science, Umea University, S-901 87 Umea, Sweden,<br>
*>         1995.<br>
*<br>
*  =====================================================================<br>
*/
	public void clatdf_(INTEGER IJOB,INTEGER N,float[] Z,INTEGER LDZ,float[] RHS,REAL RDSUM,REAL RDSCAL,int[] IPIV,int[] JPIV);
/**
*> \brief \b CLATPS solves a triangular system of equations with the matrix held in packed storage.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLATPS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatps.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatps.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatps.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE,<br>
*                          CNORM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO<br>
*       INTEGER            INFO, N<br>
*       REAL               SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               CNORM( * )<br>
*       COMPLEX            AP( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLATPS solves one of the triangular systems<br>
*><br>
*>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,<br>
*><br>
*> with scaling to prevent overflow, where A is an upper or lower<br>
*> triangular matrix stored in packed form.  Here A**T denotes the<br>
*> transpose of A, A**H denotes the conjugate transpose of A, x and b<br>
*> are n-element vectors, and s is a scaling factor, usually less than<br>
*> or equal to 1, chosen so that the components of x will be less than<br>
*> the overflow threshold.  If the unscaled problem will not cause<br>
*> overflow, the Level 2 BLAS routine CTPSV is called. If the matrix A<br>
*> is singular (A(j,j) = 0 for some j), then s is set to 0 and a<br>
*> non-trivial solution to A*x = 0 is returned.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the matrix A is upper or lower triangular.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the operation applied to A.<br>
*>          = 'N':  Solve A * x = s*b     (No transpose)<br>
*>          = 'T':  Solve A**T * x = s*b  (Transpose)<br>
*>          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          Specifies whether or not the matrix A is unit triangular.<br>
*>          = 'N':  Non-unit triangular<br>
*>          = 'U':  Unit triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] NORMIN<br>
*> \verbatim<br>
*>          NORMIN is CHARACTER*1<br>
*>          Specifies whether CNORM has been set or not.<br>
*>          = 'Y':  CNORM contains the column norms on entry<br>
*>          = 'N':  CNORM is not set on entry.  On exit, the norms will<br>
*>                  be computed and stored in CNORM.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangular matrix A, packed columnwise in<br>
*>          a linear array.  The j-th column of A is stored in the array<br>
*>          AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>          On entry, the right hand side b of the triangular system.<br>
*>          On exit, X is overwritten by the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          The scaling factor s for the triangular system<br>
*>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.<br>
*>          If SCALE = 0, the matrix A is singular or badly scaled, and<br>
*>          the vector x is an exact or approximate solution to A*x = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CNORM<br>
*> \verbatim<br>
*>          CNORM is REAL array, dimension (N)<br>
*><br>
*>          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)<br>
*>          contains the norm of the off-diagonal part of the j-th column<br>
*>          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal<br>
*>          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)<br>
*>          must be greater than or equal to the 1-norm.<br>
*><br>
*>          If NORMIN = 'N', CNORM is an output argument and CNORM(j)<br>
*>          returns the 1-norm of the offdiagonal part of the j-th column<br>
*>          of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -k, the k-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  A rough bound on x is computed; if that is less than overflow, CTPSV<br>
*>  is called, otherwise, specific code is used which checks for possible<br>
*>  overflow or divide-by-zero at every operation.<br>
*><br>
*>  A columnwise scheme is used for solving A*x = b.  The basic algorithm<br>
*>  if A is lower triangular is<br>
*><br>
*>       x[1:n] := b[1:n]<br>
*>       for j = 1, ..., n<br>
*>            x(j) := x(j) / A(j,j)<br>
*>            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]<br>
*>       end<br>
*><br>
*>  Define bounds on the components of x after j iterations of the loop:<br>
*>     M(j) = bound on x[1:j]<br>
*>     G(j) = bound on x[j+1:n]<br>
*>  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.<br>
*><br>
*>  Then for iteration j+1 we have<br>
*>     M(j+1) <= G(j) / | A(j+1,j+1) |<br>
*>     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |<br>
*>            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )<br>
*><br>
*>  where CNORM(j+1) is greater than or equal to the infinity-norm of<br>
*>  column j+1 of A, not counting the diagonal.  Hence<br>
*><br>
*>     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )<br>
*>                  1<=i<=j<br>
*>  and<br>
*><br>
*>     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )<br>
*>                                   1<=i< j<br>
*><br>
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTPSV if the<br>
*>  reciprocal of the largest M(j), j=1,..,n, is larger than<br>
*>  max(underflow, 1/overflow).<br>
*><br>
*>  The bound on x(j) is also used to determine when a step in the<br>
*>  columnwise method can be performed without fear of overflow.  If<br>
*>  the computed bound is greater than a large constant, x is scaled to<br>
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to<br>
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.<br>
*><br>
*>  Similarly, a row-wise scheme is used to solve A**T *x = b  or<br>
*>  A**H *x = b.  The basic algorithm for A upper triangular is<br>
*><br>
*>       for j = 1, ..., n<br>
*>            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)<br>
*>       end<br>
*><br>
*>  We simultaneously compute two bounds<br>
*>       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j<br>
*>       M(j) = bound on x(i), 1<=i<=j<br>
*><br>
*>  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we<br>
*>  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.<br>
*>  Then the bound on x(j) is<br>
*><br>
*>       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |<br>
*><br>
*>            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )<br>
*>                      1<=i<=j<br>
*><br>
*>  and we can safely call CTPSV if 1/M(n) and 1/G(n) are both greater<br>
*>  than max(underflow, 1/overflow).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clatps_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,float[] AP,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
/**
*> \brief \b CLATRD reduces the first nb rows and columns of a symmetric/Hermitian matrix A to real tridiagonal form by an unitary similarity transformation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLATRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            LDA, LDW, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               E( * )<br>
*       COMPLEX            A( LDA, * ), TAU( * ), W( LDW, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLATRD reduces NB rows and columns of a complex Hermitian matrix A to<br>
*> Hermitian tridiagonal form by a unitary similarity<br>
*> transformation Q**H * A * Q, and returns the matrices V and W which are<br>
*> needed to apply the transformation to the unreduced part of A.<br>
*><br>
*> If UPLO = 'U', CLATRD reduces the last NB rows and columns of a<br>
*> matrix, of which the upper triangle is supplied;<br>
*> if UPLO = 'L', CLATRD reduces the first NB rows and columns of a<br>
*> matrix, of which the lower triangle is supplied.<br>
*><br>
*> This is an auxiliary routine called by CHETRD.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          Hermitian matrix A is stored:<br>
*>          = 'U': Upper triangular<br>
*>          = 'L': Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The number of rows and columns to be reduced.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit:<br>
*>          if UPLO = 'U', the last NB columns have been reduced to<br>
*>            tridiagonal form, with the diagonal elements overwriting<br>
*>            the diagonal elements of A; the elements above the diagonal<br>
*>            with the array TAU, represent the unitary matrix Q as a<br>
*>            product of elementary reflectors;<br>
*>          if UPLO = 'L', the first NB columns have been reduced to<br>
*>            tridiagonal form, with the diagonal elements overwriting<br>
*>            the diagonal elements of A; the elements below the diagonal<br>
*>            with the array TAU, represent the  unitary matrix Q as a<br>
*>            product of elementary reflectors.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal<br>
*>          elements of the last NB columns of the reduced matrix;<br>
*>          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of<br>
*>          the first NB columns of the reduced matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX array, dimension (N-1)<br>
*>          The scalar factors of the elementary reflectors, stored in<br>
*>          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (LDW,NB)<br>
*>          The n-by-nb matrix W required to update the unreduced part<br>
*>          of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDW<br>
*> \verbatim<br>
*>          LDW is INTEGER<br>
*>          The leading dimension of the array W. LDW >= max(1,N).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', the matrix Q is represented as a product of elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(n) H(n-1) . . . H(n-nb+1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),<br>
*>  and tau in TAU(i-1).<br>
*><br>
*>  If UPLO = 'L', the matrix Q is represented as a product of elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(nb).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),<br>
*>  and tau in TAU(i).<br>
*><br>
*>  The elements of the vectors v together form the n-by-nb matrix V<br>
*>  which is needed, with W, to apply the transformation to the unreduced<br>
*>  part of the matrix, using a Hermitian rank-2k update of the form:<br>
*>  A := A - V*W**H - W*V**H.<br>
*><br>
*>  The contents of A on exit are illustrated by the following examples<br>
*>  with n = 5 and nb = 2:<br>
*><br>
*>  if UPLO = 'U':                       if UPLO = 'L':<br>
*><br>
*>    (  a   a   a   v4  v5 )              (  d                  )<br>
*>    (      a   a   v4  v5 )              (  1   d              )<br>
*>    (          a   1   v5 )              (  v1  1   a          )<br>
*>    (              d   1  )              (  v1  v2  a   a      )<br>
*>    (                  d  )              (  v1  v2  a   a   a  )<br>
*><br>
*>  where d denotes a diagonal element of the reduced matrix, a denotes<br>
*>  an element of the original matrix that is unchanged, and vi denotes<br>
*>  an element of the vector defining H(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clatrd_(CHARACTER UPLO,INTEGER N,INTEGER NB,float[] A,INTEGER LDA,float[] E,float[] TAU,float[] W,INTEGER LDW);
/**
*> \brief \b CLATRS solves a triangular system of equations with the scale factor set to prevent overflow.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLATRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,<br>
*                          CNORM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       REAL               SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               CNORM( * )<br>
*       COMPLEX            A( LDA, * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLATRS solves one of the triangular systems<br>
*><br>
*>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,<br>
*><br>
*> with scaling to prevent overflow.  Here A is an upper or lower<br>
*> triangular matrix, A**T denotes the transpose of A, A**H denotes the<br>
*> conjugate transpose of A, x and b are n-element vectors, and s is a<br>
*> scaling factor, usually less than or equal to 1, chosen so that the<br>
*> components of x will be less than the overflow threshold.  If the<br>
*> unscaled problem will not cause overflow, the Level 2 BLAS routine<br>
*> CTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),<br>
*> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the matrix A is upper or lower triangular.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the operation applied to A.<br>
*>          = 'N':  Solve A * x = s*b     (No transpose)<br>
*>          = 'T':  Solve A**T * x = s*b  (Transpose)<br>
*>          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          Specifies whether or not the matrix A is unit triangular.<br>
*>          = 'N':  Non-unit triangular<br>
*>          = 'U':  Unit triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] NORMIN<br>
*> \verbatim<br>
*>          NORMIN is CHARACTER*1<br>
*>          Specifies whether CNORM has been set or not.<br>
*>          = 'Y':  CNORM contains the column norms on entry<br>
*>          = 'N':  CNORM is not set on entry.  On exit, the norms will<br>
*>                  be computed and stored in CNORM.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The triangular matrix A.  If UPLO = 'U', the leading n by n<br>
*>          upper triangular part of the array A contains the upper<br>
*>          triangular matrix, and the strictly lower triangular part of<br>
*>          A is not referenced.  If UPLO = 'L', the leading n by n lower<br>
*>          triangular part of the array A contains the lower triangular<br>
*>          matrix, and the strictly upper triangular part of A is not<br>
*>          referenced.  If DIAG = 'U', the diagonal elements of A are<br>
*>          also not referenced and are assumed to be 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max (1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>          On entry, the right hand side b of the triangular system.<br>
*>          On exit, X is overwritten by the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          The scaling factor s for the triangular system<br>
*>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.<br>
*>          If SCALE = 0, the matrix A is singular or badly scaled, and<br>
*>          the vector x is an exact or approximate solution to A*x = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CNORM<br>
*> \verbatim<br>
*>          CNORM is REAL array, dimension (N)<br>
*><br>
*>          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)<br>
*>          contains the norm of the off-diagonal part of the j-th column<br>
*>          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal<br>
*>          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)<br>
*>          must be greater than or equal to the 1-norm.<br>
*><br>
*>          If NORMIN = 'N', CNORM is an output argument and CNORM(j)<br>
*>          returns the 1-norm of the offdiagonal part of the j-th column<br>
*>          of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -k, the k-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  A rough bound on x is computed; if that is less than overflow, CTRSV<br>
*>  is called, otherwise, specific code is used which checks for possible<br>
*>  overflow or divide-by-zero at every operation.<br>
*><br>
*>  A columnwise scheme is used for solving A*x = b.  The basic algorithm<br>
*>  if A is lower triangular is<br>
*><br>
*>       x[1:n] := b[1:n]<br>
*>       for j = 1, ..., n<br>
*>            x(j) := x(j) / A(j,j)<br>
*>            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]<br>
*>       end<br>
*><br>
*>  Define bounds on the components of x after j iterations of the loop:<br>
*>     M(j) = bound on x[1:j]<br>
*>     G(j) = bound on x[j+1:n]<br>
*>  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.<br>
*><br>
*>  Then for iteration j+1 we have<br>
*>     M(j+1) <= G(j) / | A(j+1,j+1) |<br>
*>     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |<br>
*>            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )<br>
*><br>
*>  where CNORM(j+1) is greater than or equal to the infinity-norm of<br>
*>  column j+1 of A, not counting the diagonal.  Hence<br>
*><br>
*>     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )<br>
*>                  1<=i<=j<br>
*>  and<br>
*><br>
*>     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )<br>
*>                                   1<=i< j<br>
*><br>
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTRSV if the<br>
*>  reciprocal of the largest M(j), j=1,..,n, is larger than<br>
*>  max(underflow, 1/overflow).<br>
*><br>
*>  The bound on x(j) is also used to determine when a step in the<br>
*>  columnwise method can be performed without fear of overflow.  If<br>
*>  the computed bound is greater than a large constant, x is scaled to<br>
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to<br>
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.<br>
*><br>
*>  Similarly, a row-wise scheme is used to solve A**T *x = b  or<br>
*>  A**H *x = b.  The basic algorithm for A upper triangular is<br>
*><br>
*>       for j = 1, ..., n<br>
*>            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)<br>
*>       end<br>
*><br>
*>  We simultaneously compute two bounds<br>
*>       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j<br>
*>       M(j) = bound on x(i), 1<=i<=j<br>
*><br>
*>  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we<br>
*>  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.<br>
*>  Then the bound on x(j) is<br>
*><br>
*>       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |<br>
*><br>
*>            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )<br>
*>                      1<=i<=j<br>
*><br>
*>  and we can safely call CTRSV if 1/M(n) and 1/G(n) are both greater<br>
*>  than max(underflow, 1/overflow).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clatrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,float[] A,INTEGER LDA,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
/**
*> \brief \b CLATRZ factors an upper trapezoidal matrix by means of unitary transformations.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLATRZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatrz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatrz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatrz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLATRZ( M, N, L, A, LDA, TAU, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            L, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLATRZ factors the M-by-(M+L) complex upper trapezoidal matrix<br>
*> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z by means<br>
*> of unitary transformations, where  Z is an (M+L)-by-(M+L) unitary<br>
*> matrix and, R and A1 are M-by-M upper triangular matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The number of columns of the matrix A containing the<br>
*>          meaningful part of the Householder vectors. N-M >= L >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the leading M-by-N upper trapezoidal part of the<br>
*>          array A must contain the matrix to be factorized.<br>
*>          On exit, the leading M-by-M upper triangular part of A<br>
*>          contains the upper triangular matrix R, and elements N-L+1 to<br>
*>          N of the first M rows of A, with the array TAU, represent the<br>
*>          unitary matrix Z as a product of M elementary reflectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX array, dimension (M)<br>
*>          The scalar factors of the elementary reflectors.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (M)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The factorization is obtained by Householder's method.  The kth<br>
*>  transformation matrix, Z( k ), which is used to introduce zeros into<br>
*>  the ( m - k + 1 )th row of A, is given in the form<br>
*><br>
*>     Z( k ) = ( I     0   ),<br>
*>              ( 0  T( k ) )<br>
*><br>
*>  where<br>
*><br>
*>     T( k ) = I - tau*u( k )*u( k )**H,   u( k ) = (   1    ),<br>
*>                                                 (   0    )<br>
*>                                                 ( z( k ) )<br>
*><br>
*>  tau is a scalar and z( k ) is an l element vector. tau and z( k )<br>
*>  are chosen to annihilate the elements of the kth row of A2.<br>
*><br>
*>  The scalar tau is returned in the kth element of TAU and the vector<br>
*>  u( k ) in the kth row of A2, such that the elements of z( k ) are<br>
*>  in  a( k, l + 1 ), ..., a( k, n ). The elements of R are returned in<br>
*>  the upper triangular part of A1.<br>
*><br>
*>  Z is given by<br>
*><br>
*>     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void clatrz_(INTEGER M,INTEGER N,INTEGER L,float[] A,INTEGER LDA,float[] TAU,float[] WORK);
/**
*> \brief \b CLAUU2 computes the product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAUU2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clauu2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clauu2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clauu2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAUU2( UPLO, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAUU2 computes the product U * U**H or L**H * L, where the triangular<br>
*> factor U or L is stored in the upper or lower triangular part of<br>
*> the array A.<br>
*><br>
*> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,<br>
*> overwriting the factor U in A.<br>
*> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,<br>
*> overwriting the factor L in A.<br>
*><br>
*> This is the unblocked form of the algorithm, calling Level 2 BLAS.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the triangular factor stored in the array A<br>
*>          is upper or lower triangular:<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the triangular factor U or L.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the triangular factor U or L.<br>
*>          On exit, if UPLO = 'U', the upper triangle of A is<br>
*>          overwritten with the upper triangle of the product U * U**H;<br>
*>          if UPLO = 'L', the lower triangle of A is overwritten with<br>
*>          the lower triangle of the product L**H * L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -k, the k-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clauu2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b CLAUUM computes the product UUH or LHL, where U and L are upper or lower triangular matrices (blocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLAUUM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clauum.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clauum.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clauum.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLAUUM( UPLO, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLAUUM computes the product U * U**H or L**H * L, where the triangular<br>
*> factor U or L is stored in the upper or lower triangular part of<br>
*> the array A.<br>
*><br>
*> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,<br>
*> overwriting the factor U in A.<br>
*> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,<br>
*> overwriting the factor L in A.<br>
*><br>
*> This is the blocked form of the algorithm, calling Level 3 BLAS.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the triangular factor stored in the array A<br>
*>          is upper or lower triangular:<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the triangular factor U or L.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the triangular factor U or L.<br>
*>          On exit, if UPLO = 'U', the upper triangle of A is<br>
*>          overwritten with the upper triangle of the product U * U**H;<br>
*>          if UPLO = 'L', the lower triangle of A is overwritten with<br>
*>          the lower triangle of the product L**H * L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -k, the k-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void clauum_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b CLA_GBAMV performs a matrix-vector operation to calculate error bounds.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GBAMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbamv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbamv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbamv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X,<br>
*                             INCX, BETA, Y, INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               ALPHA, BETA<br>
*       INTEGER            INCX, INCY, LDAB, M, N, KL, KU, TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AB( LDAB, * ), X( * )<br>
*       REAL               Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_GBAMV  performs one of the matrix-vector operations<br>
*><br>
*>         y := alpha*abs(A)*abs(x) + beta*abs(y),<br>
*>    or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n matrix.<br>
*><br>
*> This function is primarily used in calculating error bounds.<br>
*> To protect against underflow during evaluation, components in<br>
*> the resulting vector are perturbed away from zero by (N+1)<br>
*> times the underflow threshold.  To prevent unnecessarily large<br>
*> errors for block-structure embedded in general matrices,<br>
*> "symbolically" zero components are not perturbed.  A zero<br>
*> entry is considered "symbolic" if all multiplications involved<br>
*> in computing that entry have at least one zero multiplicand.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is INTEGER<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)<br>
*>             BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y)<br>
*>             BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y)<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>           The number of subdiagonals within the band of A.  KL >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>           The number of superdiagonals within the band of A.  KU >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,n)<br>
*>           Before entry, the leading m by n part of the array AB must<br>
*>           contain the matrix of coefficients.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>           On entry, LDAB specifies the first dimension of AB as declared<br>
*>           in the calling (sub) program. LDAB must be at least<br>
*>           max( 1, m ).<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array, dimension<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry with BETA non-zero, the incremented array Y<br>
*>           must contain the vector y. On exit, Y is overwritten by the<br>
*>           updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*>           Unchanged on exit.<br>
*><br>
*>  Level 2 Blas routine.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_gbamv_(INTEGER TRANS,INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,REAL ALPHA,float[] AB,INTEGER LDAB,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
/**
*> \brief \b CLA_GBRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for general banded matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GBRCOND_C + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbrcond_c.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbrcond_c.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbrcond_c.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_GBRCOND_C( TRANS, N, KL, KU, AB, LDAB, AFB,<br>
*                                    LDAFB, IPIV, C, CAPPLY, INFO, WORK,<br>
*                                    RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       LOGICAL            CAPPLY<br>
*       INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), WORK( * )<br>
*       REAL               C( * ), RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_GBRCOND_C Computes the infinity norm condition number of<br>
*>    op(A) * inv(diag(C)) where C is a REAL vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>     Specifies the form of the system of equations:<br>
*>       = 'N':  A * X = B     (No transpose)<br>
*>       = 'T':  A**T * X = B  (Transpose)<br>
*>       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>     The number of subdiagonals within the band of A.  KL >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>     The number of superdiagonals within the band of A.  KU >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.<br>
*>     The j-th column of A is stored in the j-th column of the<br>
*>     array AB as follows:<br>
*>     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>     The leading dimension of the array AB.  LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AFB<br>
*> \verbatim<br>
*>          AFB is COMPLEX array, dimension (LDAFB,N)<br>
*>     Details of the LU factorization of the band matrix A, as<br>
*>     computed by CGBTRF.  U is stored as an upper triangular<br>
*>     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,<br>
*>     and the multipliers used during the factorization are stored<br>
*>     in rows KL+KU+2 to 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     The pivot indices from the factorization A = P*L*U<br>
*>     as computed by CGBTRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * inv(diag(C)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CAPPLY<br>
*> \verbatim<br>
*>          CAPPLY is LOGICAL<br>
*>     If .TRUE. then access the vector C in the formula above.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_gbrcond_c_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,float[] C,LOGICAL CAPPLY,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_GBRCOND_X computes the infinity norm condition number of op(A)*diag(x) for general banded matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GBRCOND_X + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbrcond_x.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbrcond_x.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbrcond_x.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_GBRCOND_X( TRANS, N, KL, KU, AB, LDAB, AFB,<br>
*                                    LDAFB, IPIV, X, INFO, WORK, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ),<br>
*      $                   X( * )<br>
*       REAL               RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_GBRCOND_X Computes the infinity norm condition number of<br>
*>    op(A) * diag(X) where X is a COMPLEX vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>     Specifies the form of the system of equations:<br>
*>       = 'N':  A * X = B     (No transpose)<br>
*>       = 'T':  A**T * X = B  (Transpose)<br>
*>       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>     The number of subdiagonals within the band of A.  KL >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>     The number of superdiagonals within the band of A.  KU >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.<br>
*>     The j-th column of A is stored in the j-th column of the<br>
*>     array AB as follows:<br>
*>     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>     The leading dimension of the array AB.  LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AFB<br>
*> \verbatim<br>
*>          AFB is COMPLEX array, dimension (LDAFB,N)<br>
*>     Details of the LU factorization of the band matrix A, as<br>
*>     computed by CGBTRF.  U is stored as an upper triangular<br>
*>     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,<br>
*>     and the multipliers used during the factorization are stored<br>
*>     in rows KL+KU+2 to 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     The pivot indices from the factorization A = P*L*U<br>
*>     as computed by CGBTRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>     The vector X in the formula op(A) * diag(X).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_gbrcond_x_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,float[] X,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_GBRFSX_EXTENDED improves the computed solution to a system of linear equations for general banded matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GBRFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbrfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbrfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbrfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_GBRFSX_EXTENDED ( PREC_TYPE, TRANS_TYPE, N, KL, KU,<br>
*                                       NRHS, AB, LDAB, AFB, LDAFB, IPIV,<br>
*                                       COLEQU, C, B, LDB, Y, LDY,<br>
*                                       BERR_OUT, N_NORMS, ERR_BNDS_NORM,<br>
*                                       ERR_BNDS_COMP, RES, AYB, DY,<br>
*                                       Y_TAIL, RCOND, ITHRESH, RTHRESH,<br>
*                                       DZ_UB, IGNORE_CWISE, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDAB, LDAFB, LDB, LDY, N, KL, KU, NRHS,<br>
*      $                   PREC_TYPE, TRANS_TYPE, N_NORMS, ITHRESH<br>
*       LOGICAL            COLEQU, IGNORE_CWISE<br>
*       REAL               RTHRESH, DZ_UB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),<br>
*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )<br>
*       REAL               C( * ), AYB(*), RCOND, BERR_OUT( * ),<br>
*      $                   ERR_BNDS_NORM( NRHS, * ),<br>
*      $                   ERR_BNDS_COMP( NRHS, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_GBRFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by CGBRFSX to perform iterative refinement.<br>
*> In addition to normwise error bound, the code provides maximum<br>
*> componentwise error bound if possible. See comments for ERR_BNDS_NORM<br>
*> and ERR_BNDS_COMP for details of the error bounds. Note that this<br>
*> subroutine is only resonsible for setting the second fields of<br>
*> ERR_BNDS_NORM and ERR_BNDS_COMP.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] PREC_TYPE<br>
*> \verbatim<br>
*>          PREC_TYPE is INTEGER<br>
*>     Specifies the intermediate precision to be used in refinement.<br>
*>     The value is defined by ILAPREC(P) where P is a CHARACTER and<br>
*>     P    = 'S':  Single<br>
*>          = 'D':  Double<br>
*>          = 'I':  Indigenous<br>
*>          = 'X', 'E':  Extra<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS_TYPE<br>
*> \verbatim<br>
*>          TRANS_TYPE is INTEGER<br>
*>     Specifies the transposition operation on A.<br>
*>     The value is defined by ILATRANS(T) where T is a CHARACTER and<br>
*>     T    = 'N':  No transpose<br>
*>          = 'T':  Transpose<br>
*>          = 'C':  Conjugate transpose<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>     The number of subdiagonals within the band of A.  KL >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>     The number of superdiagonals within the band of A.  KU >= 0<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right-hand-sides, i.e., the number of columns of the<br>
*>     matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>     On entry, the N-by-N matrix AB.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>     The leading dimension of the array AB.  LDAB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AFB<br>
*> \verbatim<br>
*>          AFB is COMPLEX array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by CGBTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     The pivot indices from the factorization A = P*L*U<br>
*>     as computed by CGBTRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] COLEQU<br>
*> \verbatim<br>
*>          COLEQU is LOGICAL<br>
*>     If .TRUE. then column equilibration was done to A before calling<br>
*>     this routine. This is needed to compute the solution and error<br>
*>     bounds correctly.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The column scale factors for A. If COLEQU = .FALSE., C<br>
*>     is not accessed. If C is input, each element of C should be a power<br>
*>     of the radix to ensure a reliable solution and error estimates.<br>
*>     Scaling by powers of the radix does not cause rounding errors unless<br>
*>     the result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>     The right-hand-side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by CGBTRS.<br>
*>     On exit, the improved solution matrix Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDY<br>
*> \verbatim<br>
*>          LDY is INTEGER<br>
*>     The leading dimension of the array Y.  LDY >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR_OUT<br>
*> \verbatim<br>
*>          BERR_OUT is REAL array, dimension (NRHS)<br>
*>     On exit, BERR_OUT(j) contains the componentwise relative backward<br>
*>     error for right-hand-side j from the formula<br>
*>         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )<br>
*>     where abs(Z) is the componentwise absolute value of the matrix<br>
*>     or vector Z. This is computed by CLA_LIN_BERR.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_NORMS<br>
*> \verbatim<br>
*>          N_NORMS is INTEGER<br>
*>     Determines which error bounds to return (see ERR_BNDS_NORM<br>
*>     and ERR_BNDS_COMP).<br>
*>     If N_NORMS >= 1 return normwise error bounds.<br>
*>     If N_NORMS >= 2 return componentwise error bounds.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERR_BNDS_NORM<br>
*> \verbatim<br>
*>          ERR_BNDS_NORM is REAL array, dimension<br>
*>                    (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     normwise relative error, which is defined as follows:<br>
*><br>
*>     Normwise relative error in the ith solution vector:<br>
*>             max_j (abs(XTRUE(j,i) - X(j,i)))<br>
*>            ------------------------------<br>
*>                  max_j abs(X(j,i))<br>
*><br>
*>     The array is indexed by the type of error information as described<br>
*>     below. There currently are up to three pieces of information<br>
*>     returned.<br>
*><br>
*>     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERR_BNDS_NORM(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERR_BNDS_COMP<br>
*> \verbatim<br>
*>          ERR_BNDS_COMP is REAL array, dimension<br>
*>                    (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     componentwise relative error, which is defined as follows:<br>
*><br>
*>     Componentwise relative error in the ith solution vector:<br>
*>                    abs(XTRUE(j,i) - X(j,i))<br>
*>             max_j ----------------------<br>
*>                         abs(X(j,i))<br>
*><br>
*>     The array is indexed by the right-hand side i (on which the<br>
*>     componentwise relative error depends), and the type of error<br>
*>     information as described below. There currently are up to three<br>
*>     pieces of information returned for each right-hand side. If<br>
*>     componentwise accuracy is not requested (PARAMS(3) = 0.0), then<br>
*>     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most<br>
*>     the first (:,N_ERR_BNDS) entries are returned.<br>
*><br>
*>     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERR_BNDS_COMP(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RES<br>
*> \verbatim<br>
*>          RES is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the trailing bits of the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ITHRESH<br>
*> \verbatim<br>
*>          ITHRESH is INTEGER<br>
*>     The maximum number of residual computations allowed for<br>
*>     refinement. The default is 10. For 'aggressive' set to 100 to<br>
*>     permit convergence using approximate factorizations or<br>
*>     factorizations other than LU. If the factorization uses a<br>
*>     technique other than Gaussian elimination, the guarantees in<br>
*>     ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RTHRESH<br>
*> \verbatim<br>
*>          RTHRESH is REAL<br>
*>     Determines when to stop refinement if the error estimate stops<br>
*>     decreasing. Refinement will stop when the next solution no longer<br>
*>     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is<br>
*>     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The<br>
*>     default value is 0.5. For 'aggressive' set to 0.9 to permit<br>
*>     convergence on extremely ill-conditioned matrices. See LAWN 165<br>
*>     for more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DZ_UB<br>
*> \verbatim<br>
*>          DZ_UB is REAL<br>
*>     Determines when to start considering componentwise convergence.<br>
*>     Componentwise convergence is only considered after each component<br>
*>     of the solution Y is stable, which we definte as the relative<br>
*>     change in each component being less than DZ_UB. The default value<br>
*>     is 0.25, requiring the first bit to be stable. See LAWN 165 for<br>
*>     more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IGNORE_CWISE<br>
*> \verbatim<br>
*>          IGNORE_CWISE is LOGICAL<br>
*>     If .TRUE. then ignore componentwise convergence. Default value<br>
*>     is .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>       < 0:  if INFO = -i, the ith argument to CGBTRS had an illegal<br>
*>             value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_gbrfsx_extended_(INTEGER PREC_TYPE,INTEGER TRANS_TYPE,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,REAL SLAMCH,CHARACTER CHLA_TRANSTYPE);
/**
*> \brief \b CLA_GBRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a general banded matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GBRPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbrpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbrpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbrpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB,<br>
*                                   LDAFB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N, KL, KU, NCOLS, LDAB, LDAFB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_GBRPVGRW computes the reciprocal pivot growth factor<br>
*> norm(A)/norm(U). The "max absolute element" norm is used. If this is<br>
*> much less than 1, the stability of the LU factorization of the<br>
*> (equilibrated) matrix A could be poor. This also means that the<br>
*> solution X, estimated condition numbers, and error bounds could be<br>
*> unreliable.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>     The number of subdiagonals within the band of A.  KL >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>     The number of superdiagonals within the band of A.  KU >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NCOLS<br>
*> \verbatim<br>
*>          NCOLS is INTEGER<br>
*>     The number of columns of the matrix A.  NCOLS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.<br>
*>     The j-th column of A is stored in the j-th column of the<br>
*>     array AB as follows:<br>
*>     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>     The leading dimension of the array AB.  LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AFB<br>
*> \verbatim<br>
*>          AFB is COMPLEX array, dimension (LDAFB,N)<br>
*>     Details of the LU factorization of the band matrix A, as<br>
*>     computed by CGBTRF.  U is stored as an upper triangular<br>
*>     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,<br>
*>     and the multipliers used during the factorization are stored<br>
*>     in rows KL+KU+2 to 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_gbrpvgrw_(INTEGER N,INTEGER KL,INTEGER KU,INTEGER NCOLS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB);
/**
*> \brief \b CLA_GEAMV computes a matrix-vector product using a general matrix to calculate error bounds.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GEAMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_geamv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_geamv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_geamv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA,<br>
*                              Y, INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               ALPHA, BETA<br>
*       INTEGER            INCX, INCY, LDA, M, N<br>
*       INTEGER            TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), X( * )<br>
*       REAL               Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_GEAMV  performs one of the matrix-vector operations<br>
*><br>
*>         y := alpha*abs(A)*abs(x) + beta*abs(y),<br>
*>    or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n matrix.<br>
*><br>
*> This function is primarily used in calculating error bounds.<br>
*> To protect against underflow during evaluation, components in<br>
*> the resulting vector are perturbed away from zero by (N+1)<br>
*> times the underflow threshold.  To prevent unnecessarily large<br>
*> errors for block-structure embedded in general matrices,<br>
*> "symbolically" zero components are not perturbed.  A zero<br>
*> entry is considered "symbolic" if all multiplications involved<br>
*> in computing that entry have at least one zero multiplicand.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is INTEGER<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)<br>
*>             BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y)<br>
*>             BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y)<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,n)<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array, dimension<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry with BETA non-zero, the incremented array Y<br>
*>           must contain the vector y. On exit, Y is overwritten by the<br>
*>           updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*>           Unchanged on exit.<br>
*><br>
*>  Level 2 Blas routine.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_geamv_(INTEGER TRANS,INTEGER M,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
/**
*> \brief \b CLA_GERCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for general matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GERCOND_C + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gercond_c.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gercond_c.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gercond_c.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_GERCOND_C( TRANS, N, A, LDA, AF, LDAF, IPIV, C,<br>
*                                    CAPPLY, INFO, WORK, RWORK )<br>
* <br>
*       .. Scalar Aguments ..<br>
*       CHARACTER          TRANS<br>
*       LOGICAL            CAPPLY<br>
*       INTEGER            N, LDA, LDAF, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * )<br>
*       REAL               C( * ), RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*>    CLA_GERCOND_C computes the infinity norm condition number of<br>
*>    op(A) * inv(diag(C)) where C is a REAL vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>     Specifies the form of the system of equations:<br>
*>       = 'N':  A * X = B     (No transpose)<br>
*>       = 'T':  A**T * X = B  (Transpose)<br>
*>       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by CGETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     The pivot indices from the factorization A = P*L*U<br>
*>     as computed by CGETRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * inv(diag(C)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CAPPLY<br>
*> \verbatim<br>
*>          CAPPLY is LOGICAL<br>
*>     If .TRUE. then access the vector C in the formula above.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_gercond_c_(CHARACTER TRANS,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] C,LOGICAL CAPPLY,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_GERCOND_X computes the infinity norm condition number of op(A)*diag(x) for general matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GERCOND_X + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gercond_x.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gercond_x.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gercond_x.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_GERCOND_X( TRANS, N, A, LDA, AF, LDAF, IPIV, X,<br>
*                                    INFO, WORK, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            N, LDA, LDAF, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )<br>
*       REAL               RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*>    CLA_GERCOND_X computes the infinity norm condition number of<br>
*>    op(A) * diag(X) where X is a COMPLEX vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>     Specifies the form of the system of equations:<br>
*>       = 'N':  A * X = B     (No transpose)<br>
*>       = 'T':  A**T * X = B  (Transpose)<br>
*>       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by CGETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     The pivot indices from the factorization A = P*L*U<br>
*>     as computed by CGETRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>     The vector X in the formula op(A) * diag(X).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_gercond_x_(CHARACTER TRANS,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] X,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_GERFSX_EXTENDED<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GERFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gerfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gerfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gerfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A,<br>
*                                       LDA, AF, LDAF, IPIV, COLEQU, C, B,<br>
*                                       LDB, Y, LDY, BERR_OUT, N_NORMS,<br>
*                                       ERRS_N, ERRS_C, RES, AYB, DY,<br>
*                                       Y_TAIL, RCOND, ITHRESH, RTHRESH,<br>
*                                       DZ_UB, IGNORE_CWISE, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE,<br>
*      $                   TRANS_TYPE, N_NORMS<br>
*       LOGICAL            COLEQU, IGNORE_CWISE<br>
*       INTEGER            ITHRESH<br>
*       REAL               RTHRESH, DZ_UB<br>
*       ..<br>
*       .. Array Arguments<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )<br>
*       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ),<br>
*      $                   ERRS_N( NRHS, * ), ERRS_C( NRHS, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*> CLA_GERFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by CGERFSX to perform iterative refinement.<br>
*> In addition to normwise error bound, the code provides maximum<br>
*> componentwise error bound if possible. See comments for ERRS_N<br>
*> and ERRS_C for details of the error bounds. Note that this<br>
*> subroutine is only resonsible for setting the second fields of<br>
*> ERRS_N and ERRS_C.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] PREC_TYPE<br>
*> \verbatim<br>
*>          PREC_TYPE is INTEGER<br>
*>     Specifies the intermediate precision to be used in refinement.<br>
*>     The value is defined by ILAPREC(P) where P is a CHARACTER and<br>
*>     P    = 'S':  Single<br>
*>          = 'D':  Double<br>
*>          = 'I':  Indigenous<br>
*>          = 'X', 'E':  Extra<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS_TYPE<br>
*> \verbatim<br>
*>          TRANS_TYPE is INTEGER<br>
*>     Specifies the transposition operation on A.<br>
*>     The value is defined by ILATRANS(T) where T is a CHARACTER and<br>
*>     T    = 'N':  No transpose<br>
*>          = 'T':  Transpose<br>
*>          = 'C':  Conjugate transpose<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right-hand-sides, i.e., the number of columns of the<br>
*>     matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by CGETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     The pivot indices from the factorization A = P*L*U<br>
*>     as computed by CGETRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] COLEQU<br>
*> \verbatim<br>
*>          COLEQU is LOGICAL<br>
*>     If .TRUE. then column equilibration was done to A before calling<br>
*>     this routine. This is needed to compute the solution and error<br>
*>     bounds correctly.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The column scale factors for A. If COLEQU = .FALSE., C<br>
*>     is not accessed. If C is input, each element of C should be a power<br>
*>     of the radix to ensure a reliable solution and error estimates.<br>
*>     Scaling by powers of the radix does not cause rounding errors unless<br>
*>     the result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>     The right-hand-side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by CGETRS.<br>
*>     On exit, the improved solution matrix Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDY<br>
*> \verbatim<br>
*>          LDY is INTEGER<br>
*>     The leading dimension of the array Y.  LDY >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR_OUT<br>
*> \verbatim<br>
*>          BERR_OUT is REAL array, dimension (NRHS)<br>
*>     On exit, BERR_OUT(j) contains the componentwise relative backward<br>
*>     error for right-hand-side j from the formula<br>
*>         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )<br>
*>     where abs(Z) is the componentwise absolute value of the matrix<br>
*>     or vector Z. This is computed by CLA_LIN_BERR.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_NORMS<br>
*> \verbatim<br>
*>          N_NORMS is INTEGER<br>
*>     Determines which error bounds to return (see ERRS_N<br>
*>     and ERRS_C).<br>
*>     If N_NORMS >= 1 return normwise error bounds.<br>
*>     If N_NORMS >= 2 return componentwise error bounds.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERRS_N<br>
*> \verbatim<br>
*>          ERRS_N is REAL array, dimension (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     normwise relative error, which is defined as follows:<br>
*><br>
*>     Normwise relative error in the ith solution vector:<br>
*>             max_j (abs(XTRUE(j,i) - X(j,i)))<br>
*>            ------------------------------<br>
*>                  max_j abs(X(j,i))<br>
*><br>
*>     The array is indexed by the type of error information as described<br>
*>     below. There currently are up to three pieces of information<br>
*>     returned.<br>
*><br>
*>     The first index in ERRS_N(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERRS_N(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERRS_C<br>
*> \verbatim<br>
*>          ERRS_C is REAL array, dimension (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     componentwise relative error, which is defined as follows:<br>
*><br>
*>     Componentwise relative error in the ith solution vector:<br>
*>                    abs(XTRUE(j,i) - X(j,i))<br>
*>             max_j ----------------------<br>
*>                         abs(X(j,i))<br>
*><br>
*>     The array is indexed by the right-hand side i (on which the<br>
*>     componentwise relative error depends), and the type of error<br>
*>     information as described below. There currently are up to three<br>
*>     pieces of information returned for each right-hand side. If<br>
*>     componentwise accuracy is not requested (PARAMS(3) = 0.0), then<br>
*>     ERRS_C is not accessed.  If N_ERR_BNDS .LT. 3, then at most<br>
*>     the first (:,N_ERR_BNDS) entries are returned.<br>
*><br>
*>     The first index in ERRS_C(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERRS_C(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RES<br>
*> \verbatim<br>
*>          RES is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the trailing bits of the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ITHRESH<br>
*> \verbatim<br>
*>          ITHRESH is INTEGER<br>
*>     The maximum number of residual computations allowed for<br>
*>     refinement. The default is 10. For 'aggressive' set to 100 to<br>
*>     permit convergence using approximate factorizations or<br>
*>     factorizations other than LU. If the factorization uses a<br>
*>     technique other than Gaussian elimination, the guarantees in<br>
*>     ERRS_N and ERRS_C may no longer be trustworthy.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RTHRESH<br>
*> \verbatim<br>
*>          RTHRESH is REAL<br>
*>     Determines when to stop refinement if the error estimate stops<br>
*>     decreasing. Refinement will stop when the next solution no longer<br>
*>     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is<br>
*>     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The<br>
*>     default value is 0.5. For 'aggressive' set to 0.9 to permit<br>
*>     convergence on extremely ill-conditioned matrices. See LAWN 165<br>
*>     for more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DZ_UB<br>
*> \verbatim<br>
*>          DZ_UB is REAL<br>
*>     Determines when to start considering componentwise convergence.<br>
*>     Componentwise convergence is only considered after each component<br>
*>     of the solution Y is stable, which we definte as the relative<br>
*>     change in each component being less than DZ_UB. The default value<br>
*>     is 0.25, requiring the first bit to be stable. See LAWN 165 for<br>
*>     more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IGNORE_CWISE<br>
*> \verbatim<br>
*>          IGNORE_CWISE is LOGICAL<br>
*>     If .TRUE. then ignore componentwise convergence. Default value<br>
*>     is .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>       < 0:  if INFO = -i, the ith argument to CGETRS had an illegal<br>
*>             value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2011<br>
*<br>
*> \ingroup complexGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_gerfsx_extended_(INTEGER PREC_TYPE,INTEGER TRANS_TYPE,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERRS_N,float[] ERRS_C,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,REAL SLAMCH,CHARACTER CHLA_TRANSTYPE);
/**
*> \brief \b CLA_GERPVGRW multiplies a square real matrix by a complex matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_GERPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gerpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gerpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gerpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N, NCOLS, LDA, LDAF<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*> CLA_GERPVGRW computes the reciprocal pivot growth factor<br>
*> norm(A)/norm(U). The "max absolute element" norm is used. If this is<br>
*> much less than 1, the stability of the LU factorization of the<br>
*> (equilibrated) matrix A could be poor. This also means that the<br>
*> solution X, estimated condition numbers, and error bounds could be<br>
*> unreliable.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NCOLS<br>
*> \verbatim<br>
*>          NCOLS is INTEGER<br>
*>     The number of columns of the matrix A. NCOLS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by CGETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_gerpvgrw_(INTEGER N,INTEGER NCOLS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF);
/**
*> \brief \b CLA_HEAMV computes a matrix-vector product using a Hermitian indefinite matrix to calculate error bounds.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_HEAMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_heamv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_heamv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_heamv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_HEAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,<br>
*                             INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               ALPHA, BETA<br>
*       INTEGER            INCX, INCY, LDA, N, UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), X( * )<br>
*       REAL               Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_SYAMV  performs the matrix-vector operation<br>
*><br>
*>         y := alpha*abs(A)*abs(x) + beta*abs(y),<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> n by n symmetric matrix.<br>
*><br>
*> This function is primarily used in calculating error bounds.<br>
*> To protect against underflow during evaluation, components in<br>
*> the resulting vector are perturbed away from zero by (N+1)<br>
*> times the underflow threshold.  To prevent unnecessarily large<br>
*> errors for block-structure embedded in general matrices,<br>
*> "symbolically" zero components are not perturbed.  A zero<br>
*> entry is considered "symbolic" if all multiplications involved<br>
*> in computing that entry have at least one zero multiplicand.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is INTEGER<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the array A is to be referenced as<br>
*>           follows:<br>
*><br>
*>              UPLO = BLAS_UPPER   Only the upper triangular part of A<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = BLAS_LOWER   Only the lower triangular part of A<br>
*>                                  is to be referenced.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL .<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) )<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL .<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array, dimension<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) )<br>
*>           Before entry with BETA non-zero, the incremented array Y<br>
*>           must contain the vector y. On exit, Y is overwritten by the<br>
*>           updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexHEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*>  -- Modified for the absolute-value product, April 2006<br>
*>     Jason Riedy, UC Berkeley<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cla_heamv_(INTEGER UPLO,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
/**
*> \brief \b CLA_HERCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for Hermitian indefinite matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_HERCOND_C + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_hercond_c.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_hercond_c.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_hercond_c.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_HERCOND_C( UPLO, N, A, LDA, AF, LDAF, IPIV, C,<br>
*                                    CAPPLY, INFO, WORK, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       LOGICAL            CAPPLY<br>
*       INTEGER            N, LDA, LDAF, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * )<br>
*       REAL               C ( * ), RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_HERCOND_C computes the infinity norm condition number of<br>
*>    op(A) * inv(diag(C)) where C is a REAL vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by CHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     Details of the interchanges and the block structure of D<br>
*>     as determined by CHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * inv(diag(C)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CAPPLY<br>
*> \verbatim<br>
*>          CAPPLY is LOGICAL<br>
*>     If .TRUE. then access the vector C in the formula above.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexHEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_hercond_c_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] C,LOGICAL CAPPLY,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_HERCOND_X computes the infinity norm condition number of op(A)*diag(x) for Hermitian indefinite matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_HERCOND_X + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_hercond_x.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_hercond_x.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_hercond_x.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_HERCOND_X( UPLO, N, A, LDA, AF, LDAF, IPIV, X,<br>
*                                    INFO, WORK, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            N, LDA, LDAF, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )<br>
*       REAL               RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_HERCOND_X computes the infinity norm condition number of<br>
*>    op(A) * diag(X) where X is a COMPLEX vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by CHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     Details of the interchanges and the block structure of D<br>
*>     as determined by CHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>     The vector X in the formula op(A) * diag(X).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexHEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_hercond_x_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] X,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_HERFSX_EXTENDED improves the computed solution to a system of linear equations for Hermitian indefinite matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_HERFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_herfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_herfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_herfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_HERFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA,<br>
*                                       AF, LDAF, IPIV, COLEQU, C, B, LDB,<br>
*                                       Y, LDY, BERR_OUT, N_NORMS,<br>
*                                       ERR_BNDS_NORM, ERR_BNDS_COMP, RES,<br>
*                                       AYB, DY, Y_TAIL, RCOND, ITHRESH,<br>
*                                       RTHRESH, DZ_UB, IGNORE_CWISE,<br>
*                                       INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE,<br>
*      $                   N_NORMS, ITHRESH<br>
*       CHARACTER          UPLO<br>
*       LOGICAL            COLEQU, IGNORE_CWISE<br>
*       REAL               RTHRESH, DZ_UB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )<br>
*       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ),<br>
*      $                   ERR_BNDS_NORM( NRHS, * ),<br>
*      $                   ERR_BNDS_COMP( NRHS, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_HERFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by CHERFSX to perform iterative refinement.<br>
*> In addition to normwise error bound, the code provides maximum<br>
*> componentwise error bound if possible. See comments for ERR_BNDS_NORM<br>
*> and ERR_BNDS_COMP for details of the error bounds. Note that this<br>
*> subroutine is only resonsible for setting the second fields of<br>
*> ERR_BNDS_NORM and ERR_BNDS_COMP.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] PREC_TYPE<br>
*> \verbatim<br>
*>          PREC_TYPE is INTEGER<br>
*>     Specifies the intermediate precision to be used in refinement.<br>
*>     The value is defined by ILAPREC(P) where P is a CHARACTER and<br>
*>     P    = 'S':  Single<br>
*>          = 'D':  Double<br>
*>          = 'I':  Indigenous<br>
*>          = 'X', 'E':  Extra<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right-hand-sides, i.e., the number of columns of the<br>
*>     matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by CHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     Details of the interchanges and the block structure of D<br>
*>     as determined by CHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COLEQU<br>
*> \verbatim<br>
*>          COLEQU is LOGICAL<br>
*>     If .TRUE. then column equilibration was done to A before calling<br>
*>     this routine. This is needed to compute the solution and error<br>
*>     bounds correctly.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The column scale factors for A. If COLEQU = .FALSE., C<br>
*>     is not accessed. If C is input, each element of C should be a power<br>
*>     of the radix to ensure a reliable solution and error estimates.<br>
*>     Scaling by powers of the radix does not cause rounding errors unless<br>
*>     the result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>     The right-hand-side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension<br>
*>                    (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by CHETRS.<br>
*>     On exit, the improved solution matrix Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDY<br>
*> \verbatim<br>
*>          LDY is INTEGER<br>
*>     The leading dimension of the array Y.  LDY >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR_OUT<br>
*> \verbatim<br>
*>          BERR_OUT is REAL array, dimension (NRHS)<br>
*>     On exit, BERR_OUT(j) contains the componentwise relative backward<br>
*>     error for right-hand-side j from the formula<br>
*>         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )<br>
*>     where abs(Z) is the componentwise absolute value of the matrix<br>
*>     or vector Z. This is computed by CLA_LIN_BERR.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_NORMS<br>
*> \verbatim<br>
*>          N_NORMS is INTEGER<br>
*>     Determines which error bounds to return (see ERR_BNDS_NORM<br>
*>     and ERR_BNDS_COMP).<br>
*>     If N_NORMS >= 1 return normwise error bounds.<br>
*>     If N_NORMS >= 2 return componentwise error bounds.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERR_BNDS_NORM<br>
*> \verbatim<br>
*>          ERR_BNDS_NORM is REAL array, dimension<br>
*>                    (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     normwise relative error, which is defined as follows:<br>
*><br>
*>     Normwise relative error in the ith solution vector:<br>
*>             max_j (abs(XTRUE(j,i) - X(j,i)))<br>
*>            ------------------------------<br>
*>                  max_j abs(X(j,i))<br>
*><br>
*>     The array is indexed by the type of error information as described<br>
*>     below. There currently are up to three pieces of information<br>
*>     returned.<br>
*><br>
*>     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERR_BNDS_NORM(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERR_BNDS_COMP<br>
*> \verbatim<br>
*>          ERR_BNDS_COMP is REAL array, dimension<br>
*>                    (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     componentwise relative error, which is defined as follows:<br>
*><br>
*>     Componentwise relative error in the ith solution vector:<br>
*>                    abs(XTRUE(j,i) - X(j,i))<br>
*>             max_j ----------------------<br>
*>                         abs(X(j,i))<br>
*><br>
*>     The array is indexed by the right-hand side i (on which the<br>
*>     componentwise relative error depends), and the type of error<br>
*>     information as described below. There currently are up to three<br>
*>     pieces of information returned for each right-hand side. If<br>
*>     componentwise accuracy is not requested (PARAMS(3) = 0.0), then<br>
*>     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most<br>
*>     the first (:,N_ERR_BNDS) entries are returned.<br>
*><br>
*>     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERR_BNDS_COMP(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RES<br>
*> \verbatim<br>
*>          RES is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the trailing bits of the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ITHRESH<br>
*> \verbatim<br>
*>          ITHRESH is INTEGER<br>
*>     The maximum number of residual computations allowed for<br>
*>     refinement. The default is 10. For 'aggressive' set to 100 to<br>
*>     permit convergence using approximate factorizations or<br>
*>     factorizations other than LU. If the factorization uses a<br>
*>     technique other than Gaussian elimination, the guarantees in<br>
*>     ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RTHRESH<br>
*> \verbatim<br>
*>          RTHRESH is REAL<br>
*>     Determines when to stop refinement if the error estimate stops<br>
*>     decreasing. Refinement will stop when the next solution no longer<br>
*>     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is<br>
*>     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The<br>
*>     default value is 0.5. For 'aggressive' set to 0.9 to permit<br>
*>     convergence on extremely ill-conditioned matrices. See LAWN 165<br>
*>     for more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DZ_UB<br>
*> \verbatim<br>
*>          DZ_UB is REAL<br>
*>     Determines when to start considering componentwise convergence.<br>
*>     Componentwise convergence is only considered after each component<br>
*>     of the solution Y is stable, which we definte as the relative<br>
*>     change in each component being less than DZ_UB. The default value<br>
*>     is 0.25, requiring the first bit to be stable. See LAWN 165 for<br>
*>     more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IGNORE_CWISE<br>
*> \verbatim<br>
*>          IGNORE_CWISE is LOGICAL<br>
*>     If .TRUE. then ignore componentwise convergence. Default value<br>
*>     is .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>       < 0:  if INFO = -i, the ith argument to CLA_HERFSX_EXTENDED had an illegal<br>
*>             value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexHEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_herfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
/**
*> \brief \b CLA_HERPVGRW<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_HERPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_herpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_herpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_herpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_HERPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV,<br>
*                                   WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*1        UPLO<br>
*       INTEGER            N, INFO, LDA, LDAF<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * )<br>
*       REAL               WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*> CLA_HERPVGRW computes the reciprocal pivot growth factor<br>
*> norm(A)/norm(U). The "max absolute element" norm is used. If this is<br>
*> much less than 1, the stability of the LU factorization of the<br>
*> (equilibrated) matrix A could be poor. This also means that the<br>
*> solution X, estimated condition numbers, and error bounds could be<br>
*> unreliable.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>     The value of INFO returned from SSYTRF, .i.e., the pivot in<br>
*>     column INFO is exactly 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by CHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     Details of the interchanges and the block structure of D<br>
*>     as determined by CHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexHEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_herpvgrw_(byte[] UPLO,INTEGER N,INTEGER INFO,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] WORK);
/**
*> \brief \b CLA_LIN_BERR computes a component-wise relative backward error.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_LIN_BERR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_lin_berr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_lin_berr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_lin_berr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N, NZ, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AYB( N, NRHS ), BERR( NRHS )<br>
*       COMPLEX            RES( N, NRHS )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_LIN_BERR computes componentwise relative backward error from<br>
*>    the formula<br>
*>        max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )<br>
*>    where abs(Z) is the componentwise absolute value of the matrix<br>
*>    or vector Z.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NZ<br>
*> \verbatim<br>
*>          NZ is INTEGER<br>
*>     We add (NZ+1)*SLAMCH( 'Safe minimum' ) to R(i) in the numerator to<br>
*>     guard against spuriously zero residuals. Default value is N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right hand sides, i.e., the number of columns<br>
*>     of the matrices AYB, RES, and BERR.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RES<br>
*> \verbatim<br>
*>          RES is COMPLEX array, dimension (N,NRHS)<br>
*>     The residual matrix, i.e., the matrix R in the relative backward<br>
*>     error formula above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N, NRHS)<br>
*>     The denominator in the relative backward error formula above, i.e.,<br>
*>     the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B<br>
*>     are from iterative refinement (see cla_gerfsx_extended.f).<br>
*> \endverbatim<br>
*>     <br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is REAL array, dimension (NRHS)<br>
*>     The componentwise relative backward error from the formula above.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_lin_berr_(INTEGER N,INTEGER NZ,INTEGER NRHS,float[] RES,float[] AYB,float[] BERR);
/**
*> \brief \b CLA_PORCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for Hermitian positive-definite matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_PORCOND_C + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_porcond_c.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_porcond_c.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_porcond_c.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_PORCOND_C( UPLO, N, A, LDA, AF, LDAF, C, CAPPLY,<br>
*                                    INFO, WORK, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       LOGICAL            CAPPLY<br>
*       INTEGER            N, LDA, LDAF, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * )<br>
*       REAL               C( * ), RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_PORCOND_C Computes the infinity norm condition number of<br>
*>    op(A) * inv(diag(C)) where C is a REAL vector<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The triangular factor U or L from the Cholesky factorization<br>
*>     A = U**H*U or A = L*L**H, as computed by CPOTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * inv(diag(C)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CAPPLY<br>
*> \verbatim<br>
*>          CAPPLY is LOGICAL<br>
*>     If .TRUE. then access the vector C in the formula above.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_porcond_c_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] C,LOGICAL CAPPLY,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_PORCOND_X computes the infinity norm condition number of op(A)*diag(x) for Hermitian positive-definite matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_PORCOND_X + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_porcond_x.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_porcond_x.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_porcond_x.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_PORCOND_X( UPLO, N, A, LDA, AF, LDAF, X, INFO,<br>
*                                    WORK, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            N, LDA, LDAF, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )<br>
*       REAL               RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_PORCOND_X Computes the infinity norm condition number of<br>
*>    op(A) * diag(X) where X is a COMPLEX vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The triangular factor U or L from the Cholesky factorization<br>
*>     A = U**H*U or A = L*L**H, as computed by CPOTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>     The vector X in the formula op(A) * diag(X).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_porcond_x_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] X,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_PORFSX_EXTENDED improves the computed solution to a system of linear equations for symmetric or Hermitian positive-definite matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_PORFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_porfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_porfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_porfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_PORFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA,<br>
*                                       AF, LDAF, COLEQU, C, B, LDB, Y,<br>
*                                       LDY, BERR_OUT, N_NORMS,<br>
*                                       ERR_BNDS_NORM, ERR_BNDS_COMP, RES,<br>
*                                       AYB, DY, Y_TAIL, RCOND, ITHRESH,<br>
*                                       RTHRESH, DZ_UB, IGNORE_CWISE,<br>
*                                       INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE,<br>
*      $                   N_NORMS, ITHRESH<br>
*       CHARACTER          UPLO<br>
*       LOGICAL            COLEQU, IGNORE_CWISE<br>
*       REAL               RTHRESH, DZ_UB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )<br>
*       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ),<br>
*      $                   ERR_BNDS_NORM( NRHS, * ),<br>
*      $                   ERR_BNDS_COMP( NRHS, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_PORFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by CPORFSX to perform iterative refinement.<br>
*> In addition to normwise error bound, the code provides maximum<br>
*> componentwise error bound if possible. See comments for ERR_BNDS_NORM<br>
*> and ERR_BNDS_COMP for details of the error bounds. Note that this<br>
*> subroutine is only resonsible for setting the second fields of<br>
*> ERR_BNDS_NORM and ERR_BNDS_COMP.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] PREC_TYPE<br>
*> \verbatim<br>
*>          PREC_TYPE is INTEGER<br>
*>     Specifies the intermediate precision to be used in refinement.<br>
*>     The value is defined by ILAPREC(P) where P is a CHARACTER and<br>
*>     P    = 'S':  Single<br>
*>          = 'D':  Double<br>
*>          = 'I':  Indigenous<br>
*>          = 'X', 'E':  Extra<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right-hand-sides, i.e., the number of columns of the<br>
*>     matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The triangular factor U or L from the Cholesky factorization<br>
*>     A = U**T*U or A = L*L**T, as computed by CPOTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] COLEQU<br>
*> \verbatim<br>
*>          COLEQU is LOGICAL<br>
*>     If .TRUE. then column equilibration was done to A before calling<br>
*>     this routine. This is needed to compute the solution and error<br>
*>     bounds correctly.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The column scale factors for A. If COLEQU = .FALSE., C<br>
*>     is not accessed. If C is input, each element of C should be a power<br>
*>     of the radix to ensure a reliable solution and error estimates.<br>
*>     Scaling by powers of the radix does not cause rounding errors unless<br>
*>     the result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>     The right-hand-side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension<br>
*>                    (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by CPOTRS.<br>
*>     On exit, the improved solution matrix Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDY<br>
*> \verbatim<br>
*>          LDY is INTEGER<br>
*>     The leading dimension of the array Y.  LDY >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR_OUT<br>
*> \verbatim<br>
*>          BERR_OUT is REAL array, dimension (NRHS)<br>
*>     On exit, BERR_OUT(j) contains the componentwise relative backward<br>
*>     error for right-hand-side j from the formula<br>
*>         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )<br>
*>     where abs(Z) is the componentwise absolute value of the matrix<br>
*>     or vector Z. This is computed by CLA_LIN_BERR.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_NORMS<br>
*> \verbatim<br>
*>          N_NORMS is INTEGER<br>
*>     Determines which error bounds to return (see ERR_BNDS_NORM<br>
*>     and ERR_BNDS_COMP).<br>
*>     If N_NORMS >= 1 return normwise error bounds.<br>
*>     If N_NORMS >= 2 return componentwise error bounds.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERR_BNDS_NORM<br>
*> \verbatim<br>
*>          ERR_BNDS_NORM is REAL array, dimension<br>
*>                    (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     normwise relative error, which is defined as follows:<br>
*><br>
*>     Normwise relative error in the ith solution vector:<br>
*>             max_j (abs(XTRUE(j,i) - X(j,i)))<br>
*>            ------------------------------<br>
*>                  max_j abs(X(j,i))<br>
*><br>
*>     The array is indexed by the type of error information as described<br>
*>     below. There currently are up to three pieces of information<br>
*>     returned.<br>
*><br>
*>     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERR_BNDS_NORM(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERR_BNDS_COMP<br>
*> \verbatim<br>
*>          ERR_BNDS_COMP is REAL array, dimension<br>
*>                    (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     componentwise relative error, which is defined as follows:<br>
*><br>
*>     Componentwise relative error in the ith solution vector:<br>
*>                    abs(XTRUE(j,i) - X(j,i))<br>
*>             max_j ----------------------<br>
*>                         abs(X(j,i))<br>
*><br>
*>     The array is indexed by the right-hand side i (on which the<br>
*>     componentwise relative error depends), and the type of error<br>
*>     information as described below. There currently are up to three<br>
*>     pieces of information returned for each right-hand side. If<br>
*>     componentwise accuracy is not requested (PARAMS(3) = 0.0), then<br>
*>     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most<br>
*>     the first (:,N_ERR_BNDS) entries are returned.<br>
*><br>
*>     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERR_BNDS_COMP(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RES<br>
*> \verbatim<br>
*>          RES is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the trailing bits of the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ITHRESH<br>
*> \verbatim<br>
*>          ITHRESH is INTEGER<br>
*>     The maximum number of residual computations allowed for<br>
*>     refinement. The default is 10. For 'aggressive' set to 100 to<br>
*>     permit convergence using approximate factorizations or<br>
*>     factorizations other than LU. If the factorization uses a<br>
*>     technique other than Gaussian elimination, the guarantees in<br>
*>     ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RTHRESH<br>
*> \verbatim<br>
*>          RTHRESH is REAL<br>
*>     Determines when to stop refinement if the error estimate stops<br>
*>     decreasing. Refinement will stop when the next solution no longer<br>
*>     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is<br>
*>     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The<br>
*>     default value is 0.5. For 'aggressive' set to 0.9 to permit<br>
*>     convergence on extremely ill-conditioned matrices. See LAWN 165<br>
*>     for more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DZ_UB<br>
*> \verbatim<br>
*>          DZ_UB is REAL<br>
*>     Determines when to start considering componentwise convergence.<br>
*>     Componentwise convergence is only considered after each component<br>
*>     of the solution Y is stable, which we definte as the relative<br>
*>     change in each component being less than DZ_UB. The default value<br>
*>     is 0.25, requiring the first bit to be stable. See LAWN 165 for<br>
*>     more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IGNORE_CWISE<br>
*> \verbatim<br>
*>          IGNORE_CWISE is LOGICAL<br>
*>     If .TRUE. then ignore componentwise convergence. Default value<br>
*>     is .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>       < 0:  if INFO = -i, the ith argument to CPOTRS had an illegal<br>
*>             value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_porfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
/**
*> \brief \b CLA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Hermitian positive-definite matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_PORPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_porpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_porpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_porpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, LDAF, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*1        UPLO<br>
*       INTEGER            NCOLS, LDA, LDAF<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * )<br>
*       REAL               WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*> CLA_PORPVGRW computes the reciprocal pivot growth factor<br>
*> norm(A)/norm(U). The "max absolute element" norm is used. If this is<br>
*> much less than 1, the stability of the LU factorization of the<br>
*> (equilibrated) matrix A could be poor. This also means that the<br>
*> solution X, estimated condition numbers, and error bounds could be<br>
*> unreliable.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NCOLS<br>
*> \verbatim<br>
*>          NCOLS is INTEGER<br>
*>     The number of columns of the matrix A. NCOLS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The triangular factor U or L from the Cholesky factorization<br>
*>     A = U**T*U or A = L*L**T, as computed by CPOTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_porpvgrw_(byte[] UPLO,INTEGER NCOLS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] WORK);
/**
*> \brief \b CLA_SYAMV computes a matrix-vector product using a symmetric indefinite matrix to calculate error bounds.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_SYAMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_syamv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_syamv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_syamv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,<br>
*                             INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               ALPHA, BETA<br>
*       INTEGER            INCX, INCY, LDA, N<br>
*       INTEGER            UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), X( * )<br>
*       REAL               Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_SYAMV  performs the matrix-vector operation<br>
*><br>
*>         y := alpha*abs(A)*abs(x) + beta*abs(y),<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> n by n symmetric matrix.<br>
*><br>
*> This function is primarily used in calculating error bounds.<br>
*> To protect against underflow during evaluation, components in<br>
*> the resulting vector are perturbed away from zero by (N+1)<br>
*> times the underflow threshold.  To prevent unnecessarily large<br>
*> errors for block-structure embedded in general matrices,<br>
*> "symbolically" zero components are not perturbed.  A zero<br>
*> entry is considered "symbolic" if all multiplications involved<br>
*> in computing that entry have at least one zero multiplicand.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is INTEGER<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the array A is to be referenced as<br>
*>           follows:<br>
*><br>
*>              UPLO = BLAS_UPPER   Only the upper triangular part of A<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = BLAS_LOWER   Only the lower triangular part of A<br>
*>                                  is to be referenced.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL .<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) )<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL .<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array, dimension<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) )<br>
*>           Before entry with BETA non-zero, the incremented array Y<br>
*>           must contain the vector y. On exit, Y is overwritten by the<br>
*>           updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexSYcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*>  -- Modified for the absolute-value product, April 2006<br>
*>     Jason Riedy, UC Berkeley<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cla_syamv_(INTEGER UPLO,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
/**
*> \brief \b CLA_SYRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for symmetric indefinite matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_SYRCOND_C + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_syrcond_c.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_syrcond_c.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_syrcond_c.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_SYRCOND_C( UPLO, N, A, LDA, AF, LDAF, IPIV, C,<br>
*                                    CAPPLY, INFO, WORK, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       LOGICAL            CAPPLY<br>
*       INTEGER            N, LDA, LDAF, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * )<br>
*       REAL               C( * ), RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_SYRCOND_C Computes the infinity norm condition number of<br>
*>    op(A) * inv(diag(C)) where C is a REAL vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by CSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     Details of the interchanges and the block structure of D<br>
*>     as determined by CSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * inv(diag(C)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CAPPLY<br>
*> \verbatim<br>
*>          CAPPLY is LOGICAL<br>
*>     If .TRUE. then access the vector C in the formula above.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_syrcond_c_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] C,LOGICAL CAPPLY,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_SYRCOND_X computes the infinity norm condition number of op(A)*diag(x) for symmetric indefinite matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_SYRCOND_X + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_syrcond_x.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_syrcond_x.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_syrcond_x.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_SYRCOND_X( UPLO, N, A, LDA, AF, LDAF, IPIV, X,<br>
*                                    INFO, WORK, RWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            N, LDA, LDAF, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )<br>
*       REAL               RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_SYRCOND_X Computes the infinity norm condition number of<br>
*>    op(A) * diag(X) where X is a COMPLEX vector.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by CSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     Details of the interchanges and the block structure of D<br>
*>     as determined by CSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>     The vector X in the formula op(A) * diag(X).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>     i > 0:  The ith argument is invalid.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_syrcond_x_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] X,INTEGER INFO,float[] WORK,float[] RWORK);
/**
*> \brief \b CLA_SYRFSX_EXTENDED improves the computed solution to a system of linear equations for symmetric indefinite matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_SYRFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_syrfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_syrfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_syrfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_SYRFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA,<br>
*                                       AF, LDAF, IPIV, COLEQU, C, B, LDB,<br>
*                                       Y, LDY, BERR_OUT, N_NORMS,<br>
*                                       ERR_BNDS_NORM, ERR_BNDS_COMP, RES,<br>
*                                       AYB, DY, Y_TAIL, RCOND, ITHRESH,<br>
*                                       RTHRESH, DZ_UB, IGNORE_CWISE,<br>
*                                       INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE,<br>
*      $                   N_NORMS, ITHRESH<br>
*       CHARACTER          UPLO<br>
*       LOGICAL            COLEQU, IGNORE_CWISE<br>
*       REAL               RTHRESH, DZ_UB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )<br>
*       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ),<br>
*      $                   ERR_BNDS_NORM( NRHS, * ),<br>
*      $                   ERR_BNDS_COMP( NRHS, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CLA_SYRFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by CSYRFSX to perform iterative refinement.<br>
*> In addition to normwise error bound, the code provides maximum<br>
*> componentwise error bound if possible. See comments for ERR_BNDS_NORM<br>
*> and ERR_BNDS_COMP for details of the error bounds. Note that this<br>
*> subroutine is only resonsible for setting the second fields of<br>
*> ERR_BNDS_NORM and ERR_BNDS_COMP.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] PREC_TYPE<br>
*> \verbatim<br>
*>          PREC_TYPE is INTEGER<br>
*>     Specifies the intermediate precision to be used in refinement.<br>
*>     The value is defined by ILAPREC(P) where P is a CHARACTER and<br>
*>     P    = 'S':  Single<br>
*>          = 'D':  Double<br>
*>          = 'I':  Indigenous<br>
*>          = 'X', 'E':  Extra<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right-hand-sides, i.e., the number of columns of the<br>
*>     matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by CSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     Details of the interchanges and the block structure of D<br>
*>     as determined by CSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COLEQU<br>
*> \verbatim<br>
*>          COLEQU is LOGICAL<br>
*>     If .TRUE. then column equilibration was done to A before calling<br>
*>     this routine. This is needed to compute the solution and error<br>
*>     bounds correctly.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The column scale factors for A. If COLEQU = .FALSE., C<br>
*>     is not accessed. If C is input, each element of C should be a power<br>
*>     of the radix to ensure a reliable solution and error estimates.<br>
*>     Scaling by powers of the radix does not cause rounding errors unless<br>
*>     the result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>     The right-hand-side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension<br>
*>                    (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by CSYTRS.<br>
*>     On exit, the improved solution matrix Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDY<br>
*> \verbatim<br>
*>          LDY is INTEGER<br>
*>     The leading dimension of the array Y.  LDY >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR_OUT<br>
*> \verbatim<br>
*>          BERR_OUT is REAL array, dimension (NRHS)<br>
*>     On exit, BERR_OUT(j) contains the componentwise relative backward<br>
*>     error for right-hand-side j from the formula<br>
*>         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )<br>
*>     where abs(Z) is the componentwise absolute value of the matrix<br>
*>     or vector Z. This is computed by CLA_LIN_BERR.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_NORMS<br>
*> \verbatim<br>
*>          N_NORMS is INTEGER<br>
*>     Determines which error bounds to return (see ERR_BNDS_NORM<br>
*>     and ERR_BNDS_COMP).<br>
*>     If N_NORMS >= 1 return normwise error bounds.<br>
*>     If N_NORMS >= 2 return componentwise error bounds.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERR_BNDS_NORM<br>
*> \verbatim<br>
*>          ERR_BNDS_NORM is REAL array, dimension<br>
*>                    (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     normwise relative error, which is defined as follows:<br>
*><br>
*>     Normwise relative error in the ith solution vector:<br>
*>             max_j (abs(XTRUE(j,i) - X(j,i)))<br>
*>            ------------------------------<br>
*>                  max_j abs(X(j,i))<br>
*><br>
*>     The array is indexed by the type of error information as described<br>
*>     below. There currently are up to three pieces of information<br>
*>     returned.<br>
*><br>
*>     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERR_BNDS_NORM(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ERR_BNDS_COMP<br>
*> \verbatim<br>
*>          ERR_BNDS_COMP is REAL array, dimension<br>
*>                    (NRHS, N_ERR_BNDS)<br>
*>     For each right-hand side, this array contains information about<br>
*>     various error bounds and condition numbers corresponding to the<br>
*>     componentwise relative error, which is defined as follows:<br>
*><br>
*>     Componentwise relative error in the ith solution vector:<br>
*>                    abs(XTRUE(j,i) - X(j,i))<br>
*>             max_j ----------------------<br>
*>                         abs(X(j,i))<br>
*><br>
*>     The array is indexed by the right-hand side i (on which the<br>
*>     componentwise relative error depends), and the type of error<br>
*>     information as described below. There currently are up to three<br>
*>     pieces of information returned for each right-hand side. If<br>
*>     componentwise accuracy is not requested (PARAMS(3) = 0.0), then<br>
*>     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most<br>
*>     the first (:,N_ERR_BNDS) entries are returned.<br>
*><br>
*>     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith<br>
*>     right-hand side.<br>
*><br>
*>     The second index in ERR_BNDS_COMP(:,err) contains the following<br>
*>     three fields:<br>
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the<br>
*>              reciprocal condition number is less than the threshold<br>
*>              sqrt(n) * slamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * slamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * slamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     This subroutine is only responsible for setting the second field<br>
*>     above.<br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RES<br>
*> \verbatim<br>
*>          RES is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is COMPLEX array, dimension (N)<br>
*>     Workspace to hold the trailing bits of the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ITHRESH<br>
*> \verbatim<br>
*>          ITHRESH is INTEGER<br>
*>     The maximum number of residual computations allowed for<br>
*>     refinement. The default is 10. For 'aggressive' set to 100 to<br>
*>     permit convergence using approximate factorizations or<br>
*>     factorizations other than LU. If the factorization uses a<br>
*>     technique other than Gaussian elimination, the guarantees in<br>
*>     ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RTHRESH<br>
*> \verbatim<br>
*>          RTHRESH is REAL<br>
*>     Determines when to stop refinement if the error estimate stops<br>
*>     decreasing. Refinement will stop when the next solution no longer<br>
*>     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is<br>
*>     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The<br>
*>     default value is 0.5. For 'aggressive' set to 0.9 to permit<br>
*>     convergence on extremely ill-conditioned matrices. See LAWN 165<br>
*>     for more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DZ_UB<br>
*> \verbatim<br>
*>          DZ_UB is REAL<br>
*>     Determines when to start considering componentwise convergence.<br>
*>     Componentwise convergence is only considered after each component<br>
*>     of the solution Y is stable, which we definte as the relative<br>
*>     change in each component being less than DZ_UB. The default value<br>
*>     is 0.25, requiring the first bit to be stable. See LAWN 165 for<br>
*>     more details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IGNORE_CWISE<br>
*> \verbatim<br>
*>          IGNORE_CWISE is LOGICAL<br>
*>     If .TRUE. then ignore componentwise convergence. Default value<br>
*>     is .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit.<br>
*>       < 0:  if INFO = -i, the ith argument to CLA_SYRFSX_EXTENDED had an illegal<br>
*>             value<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_syrfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
/**
*> \brief \b CLA_SYRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefinite matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_SYRPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_syrpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_syrpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_syrpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION CLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV,<br>
*                                   WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*1        UPLO<br>
*       INTEGER            N, INFO, LDA, LDAF<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * )<br>
*       REAL               WORK( * )<br>
*       INTEGER            IPIV( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*> CLA_SYRPVGRW computes the reciprocal pivot growth factor<br>
*> norm(A)/norm(U). The "max absolute element" norm is used. If this is<br>
*> much less than 1, the stability of the LU factorization of the<br>
*> (equilibrated) matrix A could be poor. This also means that the<br>
*> solution X, estimated condition numbers, and error bounds could be<br>
*> unreliable.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The number of linear equations, i.e., the order of the<br>
*>     matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>     The value of INFO returned from CSYTRF, .i.e., the pivot in<br>
*>     column INFO is exactly 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by CSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     Details of the interchanges and the block structure of D<br>
*>     as determined by CSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2015<br>
*<br>
*> \ingroup complexSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float cla_syrpvgrw_(byte[] UPLO,INTEGER N,INTEGER INFO,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] WORK);
/**
*> \brief \b CLA_WWADDW adds a vector into a doubled-single vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CLA_WWADDW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_wwaddw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_wwaddw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_wwaddw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CLA_WWADDW( N, X, Y, W )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            X( * ), Y( * ), W( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CLA_WWADDW adds a vector W into a doubled-single vector (X, Y).<br>
*><br>
*>    This works for all extant IBM's hex and binary floating point<br>
*>    arithmetics, but not for decimal.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>            The length of vectors X, Y, and W.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array, dimension (N)<br>
*>            The first part of the doubled-single accumulation vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array, dimension (N)<br>
*>            The second part of the doubled-single accumulation vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] W<br>
*> \verbatim<br>
*>          W is COMPLEX array, dimension (N)<br>
*>            The vector to be added.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date September 2012<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cla_wwaddw_(INTEGER N,float[] X,float[] Y,float[] W);

}