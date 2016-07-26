package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSL extends Library
{

	public static LapackSL instance = (LapackSL) Native.loadLibrary("liblapack",LapackSL.class);

/**
*> \brief \b SLABAD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLABAD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slabad.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slabad.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slabad.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLABAD( SMALL, LARGE )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               LARGE, SMALL<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLABAD takes as input the values computed by SLAMCH for underflow and<br>
*> overflow, and returns the square root of each of these values if the<br>
*> log of LARGE is sufficiently large.  This subroutine is intended to<br>
*> identify machines with a large exponent range, such as the Crays, and<br>
*> redefine the underflow and overflow limits to be the square roots of<br>
*> the values computed by SLAMCH.  This subroutine is needed because<br>
*> SLAMCH does not compensate for poor arithmetic in the upper half of<br>
*> the exponent range, as is found on a Cray.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in,out] SMALL<br>
*> \verbatim<br>
*>          SMALL is REAL<br>
*>          On entry, the underflow threshold as computed by SLAMCH.<br>
*>          On exit, if LOG10(LARGE) is sufficiently large, the square<br>
*>          root of SMALL, otherwise unchanged.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] LARGE<br>
*> \verbatim<br>
*>          LARGE is REAL<br>
*>          On entry, the overflow threshold as computed by SLAMCH.<br>
*>          On exit, if LOG10(LARGE) is sufficiently large, the square<br>
*>          root of LARGE, otherwise unchanged.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slabad_(REAL SMALL,REAL LARGE);
/**
*> \brief \b SLABRD reduces the first nb rows and columns of a general matrix to a bidiagonal form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLABRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slabrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slabrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slabrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y,<br>
*                          LDY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, LDX, LDY, M, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), D( * ), E( * ), TAUP( * ),<br>
*      $                   TAUQ( * ), X( LDX, * ), Y( LDY, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLABRD reduces the first NB rows and columns of a real general<br>
*> m by n matrix A to upper or lower bidiagonal form by an orthogonal<br>
*> transformation Q**T * A * P, and returns the matrices X and Y which<br>
*> are needed to apply the transformation to the unreduced part of A.<br>
*><br>
*> If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower<br>
*> bidiagonal form.<br>
*><br>
*> This is an auxiliary routine called by SGEBRD<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the m by n general matrix to be reduced.<br>
*>          On exit, the first NB rows and columns of the matrix are<br>
*>          overwritten; the rest of the array is unchanged.<br>
*>          If m >= n, elements on and below the diagonal in the first NB<br>
*>            columns, with the array TAUQ, represent the orthogonal<br>
*>            matrix Q as a product of elementary reflectors; and<br>
*>            elements above the diagonal in the first NB rows, with the<br>
*>            array TAUP, represent the orthogonal matrix P as a product<br>
*>            of elementary reflectors.<br>
*>          If m < n, elements below the diagonal in the first NB<br>
*>            columns, with the array TAUQ, represent the orthogonal<br>
*>            matrix Q as a product of elementary reflectors, and<br>
*>            elements on and above the diagonal in the first NB rows,<br>
*>            with the array TAUP, represent the orthogonal matrix P as<br>
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
*>          TAUQ is REAL array dimension (NB)<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix Q. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP<br>
*> \verbatim<br>
*>          TAUP is REAL array, dimension (NB)<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix P. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (LDX,NB)<br>
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
*>          Y is REAL array, dimension (LDY,NB)<br>
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
*> \ingroup realOTHERauxiliary<br>
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
*>     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T<br>
*><br>
*>  where tauq and taup are real scalars, and v and u are real vectors.<br>
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
*>  V and the nb-by-n matrix U**T which are needed, with X and Y, to apply<br>
*>  the transformation to the unreduced part of the matrix, using a block<br>
*>  update of the form:  A := A - V*Y**T - X*U**T.<br>
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
	public void slabrd_(INTEGER M,INTEGER N,INTEGER NB,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAUQ,float[] TAUP,float[] X,INTEGER LDX,float[] Y,INTEGER LDY);
/**
*> \brief \b SLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLACN2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slacn2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slacn2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slacn2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLACN2( N, V, X, ISGN, EST, KASE, ISAVE )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            KASE, N<br>
*       REAL               EST<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISGN( * ), ISAVE( 3 )<br>
*       REAL               V( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLACN2 estimates the 1-norm of a square, real matrix A.<br>
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
*>          V is REAL array, dimension (N)<br>
*>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)<br>
*>         (W is not returned).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (N)<br>
*>         On an intermediate return, X should be overwritten by<br>
*>               A * X,   if KASE=1,<br>
*>               A**T * X,  if KASE=2,<br>
*>         and SLACN2 must be re-called with all the other parameters<br>
*>         unchanged.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISGN<br>
*> \verbatim<br>
*>          ISGN is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EST<br>
*> \verbatim<br>
*>          EST is REAL<br>
*>         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be<br>
*>         unchanged from the previous call to SLACN2.<br>
*>         On exit, EST is an estimate (a lower bound) for norm(A). <br>
*> \endverbatim<br>
*><br>
*> \param[in,out] KASE<br>
*> \verbatim<br>
*>          KASE is INTEGER<br>
*>         On the initial call to SLACN2, KASE should be 0.<br>
*>         On an intermediate return, KASE will be 1 or 2, indicating<br>
*>         whether X should be overwritten by A * X  or A**T * X.<br>
*>         On the final return from SLACN2, KASE will again be 0.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Originally named SONEST, dated March 16, 1988.<br>
*><br>
*>  This is a thread safe version of SLACON, which uses the array ISAVE<br>
*>  in place of a SAVE statement, as follows:<br>
*><br>
*>     SLACON     SLACN2<br>
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
	public void slacn2_(INTEGER N,float[] V,float[] X,int[] ISGN,REAL EST,INTEGER KASE,int[] ISAVE);
/**
*> \brief \b SLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLACON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slacon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slacon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slacon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLACON( N, V, X, ISGN, EST, KASE )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            KASE, N<br>
*       REAL               EST<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISGN( * )<br>
*       REAL               V( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLACON estimates the 1-norm of a square, real matrix A.<br>
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
*>          V is REAL array, dimension (N)<br>
*>         On the final return, V = A*W,  where  EST = norm(V)/norm(W)<br>
*>         (W is not returned).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (N)<br>
*>         On an intermediate return, X should be overwritten by<br>
*>               A * X,   if KASE=1,<br>
*>               A**T * X,  if KASE=2,<br>
*>         and SLACON must be re-called with all the other parameters<br>
*>         unchanged.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISGN<br>
*> \verbatim<br>
*>          ISGN is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EST<br>
*> \verbatim<br>
*>          EST is REAL<br>
*>         On entry with KASE = 1 or 2 and JUMP = 3, EST should be<br>
*>         unchanged from the previous call to SLACON.<br>
*>         On exit, EST is an estimate (a lower bound) for norm(A). <br>
*> \endverbatim<br>
*><br>
*> \param[in,out] KASE<br>
*> \verbatim<br>
*>          KASE is INTEGER<br>
*>         On the initial call to SLACON, KASE should be 0.<br>
*>         On an intermediate return, KASE will be 1 or 2, indicating<br>
*>         whether X should be overwritten by A * X  or A**T * X.<br>
*>         On the final return from SLACON, KASE will again be 0.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>  Nick Higham, University of Manchester. \n<br>
*>  Originally named SONEST, dated March 16, 1988.<br>
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
	public void slacon_(INTEGER N,float[] V,float[] X,int[] ISGN,REAL EST,INTEGER KASE);
/**
*> \brief \b SLACPY copies all or part of one two-dimensional array to another.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLACPY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slacpy.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slacpy.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slacpy.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLACPY( UPLO, M, N, A, LDA, B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            LDA, LDB, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLACPY copies all or part of a two-dimensional matrix A to another<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          The m by n matrix A.  If UPLO = 'U', only the upper triangle<br>
*>          or trapezoid is accessed; if UPLO = 'L', only the lower<br>
*>          triangle or trapezoid is accessed.<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slacpy_(CHARACTER UPLO,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB);
/**
*> \brief \b SLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLADIV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sladiv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sladiv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sladiv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLADIV( A, B, C, D, P, Q )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               A, B, C, D, P, Q<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLADIV performs complex division in  real arithmetic<br>
*><br>
*>                       a + i*b<br>
*>            p + i*q = ---------<br>
*>                       c + i*d<br>
*><br>
*> The algorithm is due to Michael Baudin and Robert L. Smith<br>
*> and can be found in the paper<br>
*> "A Robust Complex Division in Scilab"<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL<br>
*>          The scalars a, b, c, and d in the above expression.<br>
*> \endverbatim<br>
*><br>
*> \param[out] P<br>
*> \verbatim<br>
*>          P is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is REAL<br>
*>          The scalars p and q in the above expression.<br>
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
*> \date January 2013<br>
*<br>
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void sladiv_(REAL A,REAL B,REAL C,REAL D,REAL P,REAL Q);
	public void sladiv1_(REAL A,REAL B,REAL C,REAL D,REAL P,REAL Q);
	public float sladiv2_(REAL A,REAL B,REAL C,REAL D,REAL R,REAL T);
/**
*> \brief \b SLAE2 computes the eigenvalues of a 2-by-2 symmetric matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAE2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slae2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slae2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slae2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAE2( A, B, C, RT1, RT2 )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               A, B, C, RT1, RT2<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix<br>
*>    [  A   B  ]<br>
*>    [  B   C  ].<br>
*> On return, RT1 is the eigenvalue of larger absolute value, and RT2<br>
*> is the eigenvalue of smaller absolute value.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL<br>
*>          The (1,1) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL<br>
*>          The (1,2) and (2,1) elements of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL<br>
*>          The (2,2) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT1<br>
*> \verbatim<br>
*>          RT1 is REAL<br>
*>          The eigenvalue of larger absolute value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT2<br>
*> \verbatim<br>
*>          RT2 is REAL<br>
*>          The eigenvalue of smaller absolute value.<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
*>  Overflow is possible only if RT1 is within a factor of 5 of overflow.<br>
*>  Underflow is harmless if the input data is 0 or exceeds<br>
*>     underflow_threshold / macheps.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slae2_(REAL A,REAL B,REAL C,REAL RT1,REAL RT2);
/**
*> \brief \b SLAEBZ computes the number of eigenvalues of a real symmetric tridiagonal matrix which are less than or equal to a given value, and performs other tasks required by the routine sstebz.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAEBZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaebz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaebz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaebz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL,<br>
*                          RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT,<br>
*                          NAB, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX<br>
*       REAL               ABSTOL, PIVMIN, RELTOL<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * )<br>
*       REAL               AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAEBZ contains the iteration loops which compute and use the<br>
*> function N(w), which is the count of eigenvalues of a symmetric<br>
*> tridiagonal matrix T less than or equal to its argument  w.  It<br>
*> performs a choice of two types of loops:<br>
*><br>
*> IJOB=1, followed by<br>
*> IJOB=2: It takes as input a list of intervals and returns a list of<br>
*>         sufficiently small intervals whose union contains the same<br>
*>         eigenvalues as the union of the original intervals.<br>
*>         The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.<br>
*>         The output interval (AB(j,1),AB(j,2)] will contain<br>
*>         eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.<br>
*><br>
*> IJOB=3: It performs a binary search in each input interval<br>
*>         (AB(j,1),AB(j,2)] for a point  w(j)  such that<br>
*>         N(w(j))=NVAL(j), and uses  C(j)  as the starting point of<br>
*>         the search.  If such a w(j) is found, then on output<br>
*>         AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output<br>
*>         (AB(j,1),AB(j,2)] will be a small interval containing the<br>
*>         point where N(w) jumps through NVAL(j), unless that point<br>
*>         lies outside the initial interval.<br>
*><br>
*> Note that the intervals are in all cases half-open intervals,<br>
*> i.e., of the form  (a,b] , which includes  b  but not  a .<br>
*><br>
*> To avoid underflow, the matrix should be scaled so that its largest<br>
*> element is no greater than  overflow**(1/2) * underflow**(1/4)<br>
*> in absolute value.  To assure the most accurate computation<br>
*> of small eigenvalues, the matrix should be scaled to be<br>
*> not much smaller than that, either.<br>
*><br>
*> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal<br>
*> Matrix", Report CS41, Computer Science Dept., Stanford<br>
*> University, July 21, 1966<br>
*><br>
*> Note: the arguments are, in general, *not* checked for unreasonable<br>
*> values.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] IJOB<br>
*> \verbatim<br>
*>          IJOB is INTEGER<br>
*>          Specifies what is to be done:<br>
*>          = 1:  Compute NAB for the initial intervals.<br>
*>          = 2:  Perform bisection iteration to find eigenvalues of T.<br>
*>          = 3:  Perform bisection iteration to invert N(w), i.e.,<br>
*>                to find a point which has a specified number of<br>
*>                eigenvalues of T to its left.<br>
*>          Other values will cause SLAEBZ to return with INFO=-1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NITMAX<br>
*> \verbatim<br>
*>          NITMAX is INTEGER<br>
*>          The maximum number of "levels" of bisection to be<br>
*>          performed, i.e., an interval of width W will not be made<br>
*>          smaller than 2^(-NITMAX) * W.  If not all intervals<br>
*>          have converged after NITMAX iterations, then INFO is set<br>
*>          to the number of non-converged intervals.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The dimension n of the tridiagonal matrix T.  It must be at<br>
*>          least 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MMAX<br>
*> \verbatim<br>
*>          MMAX is INTEGER<br>
*>          The maximum number of intervals.  If more than MMAX intervals<br>
*>          are generated, then SLAEBZ will quit with INFO=MMAX+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MINP<br>
*> \verbatim<br>
*>          MINP is INTEGER<br>
*>          The initial number of intervals.  It may not be greater than<br>
*>          MMAX.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NBMIN<br>
*> \verbatim<br>
*>          NBMIN is INTEGER<br>
*>          The smallest number of intervals that should be processed<br>
*>          using a vector loop.  If zero, then only the scalar loop<br>
*>          will be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is REAL<br>
*>          The minimum (absolute) width of an interval.  When an<br>
*>          interval is narrower than ABSTOL, or than RELTOL times the<br>
*>          larger (in magnitude) endpoint, then it is considered to be<br>
*>          sufficiently small, i.e., converged.  This must be at least<br>
*>          zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RELTOL<br>
*> \verbatim<br>
*>          RELTOL is REAL<br>
*>          The minimum relative width of an interval.  When an interval<br>
*>          is narrower than ABSTOL, or than RELTOL times the larger (in<br>
*>          magnitude) endpoint, then it is considered to be<br>
*>          sufficiently small, i.e., converged.  Note: this should<br>
*>          always be at least radix*machine epsilon.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum absolute value of a "pivot" in the Sturm<br>
*>          sequence loop.<br>
*>          This must be at least  max |e(j)**2|*safe_min  and at<br>
*>          least safe_min, where safe_min is at least<br>
*>          the smallest number that can divide one without overflow.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
*>          The offdiagonal elements of the tridiagonal matrix T in<br>
*>          positions 1 through N-1.  E(N) is arbitrary.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E2<br>
*> \verbatim<br>
*>          E2 is REAL array, dimension (N)<br>
*>          The squares of the offdiagonal elements of the tridiagonal<br>
*>          matrix T.  E2(N) is ignored.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] NVAL<br>
*> \verbatim<br>
*>          NVAL is INTEGER array, dimension (MINP)<br>
*>          If IJOB=1 or 2, not referenced.<br>
*>          If IJOB=3, the desired values of N(w).  The elements of NVAL<br>
*>          will be reordered to correspond with the intervals in AB.<br>
*>          Thus, NVAL(j) on output will not, in general be the same as<br>
*>          NVAL(j) on input, but it will correspond with the interval<br>
*>          (AB(j,1),AB(j,2)] on output.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is REAL array, dimension (MMAX,2)<br>
*>          The endpoints of the intervals.  AB(j,1) is  a(j), the left<br>
*>          endpoint of the j-th interval, and AB(j,2) is b(j), the<br>
*>          right endpoint of the j-th interval.  The input intervals<br>
*>          will, in general, be modified, split, and reordered by the<br>
*>          calculation.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (MMAX)<br>
*>          If IJOB=1, ignored.<br>
*>          If IJOB=2, workspace.<br>
*>          If IJOB=3, then on input C(j) should be initialized to the<br>
*>          first search point in the binary search.<br>
*> \endverbatim<br>
*><br>
*> \param[out] MOUT<br>
*> \verbatim<br>
*>          MOUT is INTEGER<br>
*>          If IJOB=1, the number of eigenvalues in the intervals.<br>
*>          If IJOB=2 or 3, the number of intervals output.<br>
*>          If IJOB=3, MOUT will equal MINP.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] NAB<br>
*> \verbatim<br>
*>          NAB is INTEGER array, dimension (MMAX,2)<br>
*>          If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)).<br>
*>          If IJOB=2, then on input, NAB(i,j) should be set.  It must<br>
*>             satisfy the condition:<br>
*>             N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)),<br>
*>             which means that in interval i only eigenvalues<br>
*>             NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,<br>
*>             NAB(i,j)=N(AB(i,j)), from a previous call to SLAEBZ with<br>
*>             IJOB=1.<br>
*>             On output, NAB(i,j) will contain<br>
*>             max(na(k),min(nb(k),N(AB(i,j)))), where k is the index of<br>
*>             the input interval that the output interval<br>
*>             (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the<br>
*>             the input values of NAB(k,1) and NAB(k,2).<br>
*>          If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),<br>
*>             unless N(w) > NVAL(i) for all search points  w , in which<br>
*>             case NAB(i,1) will not be modified, i.e., the output<br>
*>             value will be the same as the input value (modulo<br>
*>             reorderings -- see NVAL and AB), or unless N(w) < NVAL(i)<br>
*>             for all search points  w , in which case NAB(i,2) will<br>
*>             not be modified.  Normally, NAB should be set to some<br>
*>             distinctive value(s) before SLAEBZ is called.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MMAX)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MMAX)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:       All intervals converged.<br>
*>          = 1--MMAX: The last INFO intervals did not converge.<br>
*>          = MMAX+1:  More than MMAX intervals were generated.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>      This routine is intended to be called only by other LAPACK<br>
*>  routines, thus the interface is less user-friendly.  It is intended<br>
*>  for two purposes:<br>
*><br>
*>  (a) finding eigenvalues.  In this case, SLAEBZ should have one or<br>
*>      more initial intervals set up in AB, and SLAEBZ should be called<br>
*>      with IJOB=1.  This sets up NAB, and also counts the eigenvalues.<br>
*>      Intervals with no eigenvalues would usually be thrown out at<br>
*>      this point.  Also, if not all the eigenvalues in an interval i<br>
*>      are desired, NAB(i,1) can be increased or NAB(i,2) decreased.<br>
*>      For example, set NAB(i,1)=NAB(i,2)-1 to get the largest<br>
*>      eigenvalue.  SLAEBZ is then called with IJOB=2 and MMAX<br>
*>      no smaller than the value of MOUT returned by the call with<br>
*>      IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1<br>
*>      through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the<br>
*>      tolerance specified by ABSTOL and RELTOL.<br>
*><br>
*>  (b) finding an interval (a',b'] containing eigenvalues w(f),...,w(l).<br>
*>      In this case, start with a Gershgorin interval  (a,b).  Set up<br>
*>      AB to contain 2 search intervals, both initially (a,b).  One<br>
*>      NVAL element should contain  f-1  and the other should contain  l<br>
*>      , while C should contain a and b, resp.  NAB(i,1) should be -1<br>
*>      and NAB(i,2) should be N+1, to flag an error if the desired<br>
*>      interval does not lie in (a,b).  SLAEBZ is then called with<br>
*>      IJOB=3.  On exit, if w(f-1) < w(f), then one of the intervals --<br>
*>      j -- will have AB(j,1)=AB(j,2) and NAB(j,1)=NAB(j,2)=f-1, while<br>
*>      if, to the specified tolerance, w(f-k)=...=w(f+r), k > 0 and r<br>
*>      >= 0, then the interval will have  N(AB(j,1))=NAB(j,1)=f-k and<br>
*>      N(AB(j,2))=NAB(j,2)=f+r.  The cases w(l) < w(l+1) and<br>
*>      w(l-r)=...=w(l+k) are handled similarly.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slaebz_(INTEGER IJOB,INTEGER NITMAX,INTEGER N,INTEGER MMAX,INTEGER MINP,INTEGER NBMIN,REAL ABSTOL,REAL RELTOL,REAL PIVMIN,float[] D,float[] E,float[] E2,int[] NVAL,float[] AB,float[] C,INTEGER MOUT,int[] NAB,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLAED0 used by sstedc. Computes all eigenvalues and corresponding eigenvectors of an unreduced symmetric tridiagonal matrix using the divide and conquer method.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED0 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed0.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed0.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed0.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS,<br>
*                          WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAED0 computes all eigenvalues and corresponding eigenvectors of a<br>
*> symmetric tridiagonal matrix using the divide and conquer method.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ICOMPQ<br>
*> \verbatim<br>
*>          ICOMPQ is INTEGER<br>
*>          = 0:  Compute eigenvalues only.<br>
*>          = 1:  Compute eigenvectors of original dense symmetric matrix<br>
*>                also.  On entry, Q contains the orthogonal matrix used<br>
*>                to reduce the original matrix to tridiagonal form.<br>
*>          = 2:  Compute eigenvalues and eigenvectors of tridiagonal<br>
*>                matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] QSIZ<br>
*> \verbatim<br>
*>          QSIZ is INTEGER<br>
*>         The dimension of the orthogonal matrix used to reduce<br>
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
*>         On entry, the main diagonal of the tridiagonal matrix.<br>
*>         On exit, its eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>         The off-diagonal elements of the tridiagonal matrix.<br>
*>         On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ, N)<br>
*>         On entry, Q must contain an N-by-N orthogonal matrix.<br>
*>         If ICOMPQ = 0    Q is not referenced.<br>
*>         If ICOMPQ = 1    On entry, Q is a subset of the columns of the<br>
*>                          orthogonal matrix used to reduce the full<br>
*>                          matrix to tridiagonal form corresponding to<br>
*>                          the subset of the full matrix which is being<br>
*>                          decomposed at this time.<br>
*>         If ICOMPQ = 2    On entry, Q will be the identity matrix.<br>
*>                          On exit, Q contains the eigenvectors of the<br>
*>                          tridiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>         The leading dimension of the array Q.  If eigenvectors are<br>
*>         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] QSTORE<br>
*> \verbatim<br>
*>          QSTORE is REAL array, dimension (LDQS, N)<br>
*>         Referenced only when ICOMPQ = 1.  Used to store parts of<br>
*>         the eigenvector matrix when the updating matrix multiplies<br>
*>         take place.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQS<br>
*> \verbatim<br>
*>          LDQS is INTEGER<br>
*>         The leading dimension of the array QSTORE.  If ICOMPQ = 1,<br>
*>         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array,<br>
*>         If ICOMPQ = 0 or 1, the dimension of WORK must be at least<br>
*>                     1 + 3*N + 2*N*lg N + 3*N**2<br>
*>                     ( lg( N ) = smallest integer k<br>
*>                                 such that 2^k >= N )<br>
*>         If ICOMPQ = 2, the dimension of WORK must be at least<br>
*>                     4*N + N**2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array,<br>
*>         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least<br>
*>                        6 + 6*N + 5*N*lg N.<br>
*>                        ( lg( N ) = smallest integer k<br>
*>                                    such that 2^k >= N )<br>
*>         If ICOMPQ = 2, the dimension of IWORK must be at least<br>
*>                        3 + 5*N.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void slaed0_(INTEGER ICOMPQ,INTEGER QSIZ,INTEGER N,float[] D,float[] E,float[] Q,INTEGER LDQ,float[] QSTORE,INTEGER LDQS,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLAED1 used by sstedc. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is tridiagonal.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            CUTPNT, INFO, LDQ, N<br>
*       REAL               RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            INDXQ( * ), IWORK( * )<br>
*       REAL               D( * ), Q( LDQ, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAED1 computes the updated eigensystem of a diagonal<br>
*> matrix after modification by a rank-one symmetric matrix.  This<br>
*> routine is used only for the eigenproblem which requires all<br>
*> eigenvalues and eigenvectors of a tridiagonal matrix.  SLAED7 handles<br>
*> the case in which eigenvalues only or eigenvalues and eigenvectors<br>
*> of a full symmetric matrix (which was reduced to tridiagonal form)<br>
*> are desired.<br>
*><br>
*>   T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)<br>
*><br>
*>    where Z = Q**T*u, u is a vector of length N with ones in the<br>
*>    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.<br>
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
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry, the eigenvalues of the rank-1-perturbed matrix.<br>
*>         On exit, the eigenvalues of the repaired matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ,N)<br>
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
*> \param[in,out] INDXQ<br>
*> \verbatim<br>
*>          INDXQ is INTEGER array, dimension (N)<br>
*>         On entry, the permutation which separately sorts the two<br>
*>         subproblems in D into ascending order.<br>
*>         On exit, the permutation which will reintegrate the<br>
*>         subproblems back into sorted order,<br>
*>         i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         The subdiagonal entry used to create the rank-1 modification.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CUTPNT<br>
*> \verbatim<br>
*>          CUTPNT is INTEGER<br>
*>         The location of the last eigenvalue in the leading sub-matrix.<br>
*>         min(1,N) <= CUTPNT <= N/2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (4*N + N**2)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (4*N)<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA \n<br>
*>  Modified by Francoise Tisseur, University of Tennessee<br>
*><br>
*  =====================================================================<br>
*/
	public void slaed1_(INTEGER N,float[] D,float[] Q,INTEGER LDQ,int[] INDXQ,REAL RHO,INTEGER CUTPNT,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLAED2 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original matrix is tridiagonal.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W,<br>
*                          Q2, INDX, INDXC, INDXP, COLTYP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDQ, N, N1<br>
*       REAL               RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ),<br>
*      $                   INDXQ( * )<br>
*       REAL               D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),<br>
*      $                   W( * ), Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAED2 merges the two sets of eigenvalues together into a single<br>
*> sorted set.  Then it tries to deflate the size of the problem.<br>
*> There are two ways in which deflation can occur:  when two or more<br>
*> eigenvalues are close together or if there is a tiny entry in the<br>
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
*>         The number of non-deflated eigenvalues, and the order of the<br>
*>         related secular equation. 0 <= K <=N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The dimension of the symmetric tridiagonal matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N1<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*>         The location of the last eigenvalue in the leading sub-matrix.<br>
*>         min(1,N) <= N1 <= N/2.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry, D contains the eigenvalues of the two submatrices to<br>
*>         be combined.<br>
*>         On exit, D contains the trailing (N-K) updated eigenvalues<br>
*>         (those which were deflated) sorted into increasing order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ, N)<br>
*>         On entry, Q contains the eigenvectors of two submatrices in<br>
*>         the two square blocks with corners at (1,1), (N1,N1)<br>
*>         and (N1+1, N1+1), (N,N).<br>
*>         On exit, Q contains the trailing (N-K) updated eigenvectors<br>
*>         (those which were deflated) in its last N-K columns.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>         The leading dimension of the array Q.  LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] INDXQ<br>
*> \verbatim<br>
*>          INDXQ is INTEGER array, dimension (N)<br>
*>         The permutation which separately sorts the two sub-problems<br>
*>         in D into ascending order.  Note that elements in the second<br>
*>         half of this permutation must first have N1 added to their<br>
*>         values. Destroyed on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         On entry, the off-diagonal element associated with the rank-1<br>
*>         cut which originally split the two submatrices which are now<br>
*>         being recombined.<br>
*>         On exit, RHO has been modified to the value required by<br>
*>         SLAED3.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (N)<br>
*>         On entry, Z contains the updating vector (the last<br>
*>         row of the first sub-eigenvector matrix and the first row of<br>
*>         the second sub-eigenvector matrix).<br>
*>         On exit, the contents of Z have been destroyed by the updating<br>
*>         process.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DLAMDA<br>
*> \verbatim<br>
*>          DLAMDA is REAL array, dimension (N)<br>
*>         A copy of the first K eigenvalues which will be used by<br>
*>         SLAED3 to form the secular equation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>         The first k values of the final deflation-altered z-vector<br>
*>         which will be passed to SLAED3.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q2<br>
*> \verbatim<br>
*>          Q2 is REAL array, dimension (N1**2+(N-N1)**2)<br>
*>         A copy of the first K eigenvectors which will be used by<br>
*>         SLAED3 in a matrix multiply (SGEMM) to solve for the new<br>
*>         eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDX<br>
*> \verbatim<br>
*>          INDX is INTEGER array, dimension (N)<br>
*>         The permutation used to sort the contents of DLAMDA into<br>
*>         ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDXC<br>
*> \verbatim<br>
*>          INDXC is INTEGER array, dimension (N)<br>
*>         The permutation used to arrange the columns of the deflated<br>
*>         Q matrix into three groups:  the first group contains non-zero<br>
*>         elements only at and above N1, the second contains<br>
*>         non-zero elements only below N1, and the third is dense.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDXP<br>
*> \verbatim<br>
*>          INDXP is INTEGER array, dimension (N)<br>
*>         The permutation used to place deflated values of D at the end<br>
*>         of the array.  INDXP(1:K) points to the nondeflated D-values<br>
*>         and INDXP(K+1:N) points to the deflated eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[out] COLTYP<br>
*> \verbatim<br>
*>          COLTYP is INTEGER array, dimension (N)<br>
*>         During execution, a label which will indicate which of the<br>
*>         following types a column in the Q2 matrix is:<br>
*>         1 : non-zero in the upper half only;<br>
*>         2 : dense;<br>
*>         3 : non-zero in the lower half only;<br>
*>         4 : deflated.<br>
*>         On exit, COLTYP(i) is the number of columns of type i,<br>
*>         for i=1 to 4 only.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA \n<br>
*>  Modified by Francoise Tisseur, University of Tennessee<br>
*><br>
*  =====================================================================<br>
*/
	public void slaed2_(INTEGER K,INTEGER N,INTEGER N1,float[] D,float[] Q,INTEGER LDQ,int[] INDXQ,REAL RHO,float[] Z,float[] DLAMDA,float[] W,float[] Q2,int[] INDX,int[] INDXC,int[] INDXP,int[] COLTYP,INTEGER INFO);
/**
*> \brief \b SLAED3 used by sstedc. Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is tridiagonal.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX,<br>
*                          CTOT, W, S, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDQ, N, N1<br>
*       REAL               RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            CTOT( * ), INDX( * )<br>
*       REAL               D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),<br>
*      $                   S( * ), W( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAED3 finds the roots of the secular equation, as defined by the<br>
*> values in D, W, and RHO, between 1 and K.  It makes the<br>
*> appropriate calls to SLAED4 and then updates the eigenvectors by<br>
*> multiplying the matrix of eigenvectors of the pair of eigensystems<br>
*> being combined by the matrix of eigenvectors of the K-by-K system<br>
*> which is solved here.<br>
*><br>
*> This code makes very mild assumptions about floating point<br>
*> arithmetic. It will work on machines with a guard digit in<br>
*> add/subtract, or on those binary machines without guard digits<br>
*> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.<br>
*> It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of terms in the rational function to be solved by<br>
*>          SLAED4.  K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of rows and columns in the Q matrix.<br>
*>          N >= K (deflation may result in N>K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] N1<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*>          The location of the last eigenvalue in the leading submatrix.<br>
*>          min(1,N) <= N1 <= N/2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          D(I) contains the updated eigenvalues for<br>
*>          1 <= I <= K.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ,N)<br>
*>          Initially the first K columns are used as workspace.<br>
*>          On output the columns 1 to K contain<br>
*>          the updated eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>          The value of the parameter in the rank one update equation.<br>
*>          RHO >= 0 required.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DLAMDA<br>
*> \verbatim<br>
*>          DLAMDA is REAL array, dimension (K)<br>
*>          The first K elements of this array contain the old roots<br>
*>          of the deflated updating problem.  These are the poles<br>
*>          of the secular equation. May be changed on output by<br>
*>          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,<br>
*>          Cray-2, or Cray C-90, as described above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q2<br>
*> \verbatim<br>
*>          Q2 is REAL array, dimension (LDQ2, N)<br>
*>          The first K columns of this matrix contain the non-deflated<br>
*>          eigenvectors for the split problem.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INDX<br>
*> \verbatim<br>
*>          INDX is INTEGER array, dimension (N)<br>
*>          The permutation used to arrange the columns of the deflated<br>
*>          Q matrix into three groups (see SLAED2).<br>
*>          The rows of the eigenvectors found by SLAED4 must be likewise<br>
*>          permuted before the matrix multiply can take place.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CTOT<br>
*> \verbatim<br>
*>          CTOT is INTEGER array, dimension (4)<br>
*>          A count of the total number of the various types of columns<br>
*>          in Q, as described in INDX.  The fourth column type is any<br>
*>          column which has been deflated.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (K)<br>
*>          The first K elements of this array contain the components<br>
*>          of the deflation-adjusted updating vector. Destroyed on<br>
*>          output.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N1 + 1)*K<br>
*>          Will contain the eigenvectors of the repaired matrix which<br>
*>          will be multiplied by the previously accumulated eigenvectors<br>
*>          to update the system.<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA \n<br>
*>  Modified by Francoise Tisseur, University of Tennessee<br>
*><br>
*  =====================================================================<br>
*/
	public void slaed3_(INTEGER K,INTEGER N,INTEGER N1,float[] D,float[] Q,INTEGER LDQ,REAL RHO,float[] DLAMDA,float[] Q2,int[] INDX,int[] CTOT,float[] W,float[] S,INTEGER INFO);
/**
*> \brief \b SLAED4 used by sstedc. Finds a single root of the secular equation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED4 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed4.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed4.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed4.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            I, INFO, N<br>
*       REAL               DLAM, RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), DELTA( * ), Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This subroutine computes the I-th updated eigenvalue of a symmetric<br>
*> rank-one modification to a diagonal matrix whose elements are<br>
*> given in the array d, and that<br>
*><br>
*>            D(i) < D(j)  for  i < j<br>
*><br>
*> and that RHO > 0.  This is arranged by the calling routine, and is<br>
*> no loss in generality.  The rank-one modified system is thus<br>
*><br>
*>            diag( D )  +  RHO * Z * Z_transpose.<br>
*><br>
*> where we assume the Euclidean norm of Z is 1.<br>
*><br>
*> The method consists of approximating the rational functions in the<br>
*> secular equation by simpler interpolating rational functions.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The length of all arrays.<br>
*> \endverbatim<br>
*><br>
*> \param[in] I<br>
*> \verbatim<br>
*>          I is INTEGER<br>
*>         The index of the eigenvalue to be computed.  1 <= I <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         The original eigenvalues.  It is assumed that they are in<br>
*>         order, D(I) < D(J)  for I < J.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (N)<br>
*>         The components of the updating vector.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DELTA<br>
*> \verbatim<br>
*>          DELTA is REAL array, dimension (N)<br>
*>         If N .GT. 2, DELTA contains (D(j) - lambda_I) in its  j-th<br>
*>         component.  If N = 1, then DELTA(1) = 1. If N = 2, see SLAED5<br>
*>         for detail. The vector DELTA contains the information necessary<br>
*>         to construct the eigenvectors by SLAED3 and SLAED9.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         The scalar in the symmetric updating formula.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DLAM<br>
*> \verbatim<br>
*>          DLAM is REAL<br>
*>         The computed lambda_I, the I-th updated eigenvalue.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>         = 0:  successful exit<br>
*>         > 0:  if INFO = 1, the updating process failed.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  Logical variable ORGATI (origin-at-i?) is used for distinguishing<br>
*>  whether D(i) or D(i+1) is treated as the origin.<br>
*><br>
*>            ORGATI = .true.    origin at i<br>
*>            ORGATI = .false.   origin at i+1<br>
*><br>
*>   Logical variable SWTCH3 (switch-for-3-poles?) is for noting<br>
*>   if we are working with THREE poles!<br>
*><br>
*>   MAXIT is the maximum number of iterations allowed for each<br>
*>   eigenvalue.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ren-Cang Li, Computer Science Division, University of California<br>
*>     at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slaed4_(INTEGER N,INTEGER I,float[] D,float[] Z,float[] DELTA,REAL RHO,REAL DLAM,INTEGER INFO);
/**
*> \brief \b SLAED5 used by sstedc. Solves the 2-by-2 secular equation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED5 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed5.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed5.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed5.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED5( I, D, Z, DELTA, RHO, DLAM )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            I<br>
*       REAL               DLAM, RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( 2 ), DELTA( 2 ), Z( 2 )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This subroutine computes the I-th eigenvalue of a symmetric rank-one<br>
*> modification of a 2-by-2 diagonal matrix<br>
*><br>
*>            diag( D )  +  RHO * Z * transpose(Z) .<br>
*><br>
*> The diagonal elements in the array D are assumed to satisfy<br>
*><br>
*>            D(i) < D(j)  for  i < j .<br>
*><br>
*> We also assume RHO > 0 and that the Euclidean norm of the vector<br>
*> Z is one.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] I<br>
*> \verbatim<br>
*>          I is INTEGER<br>
*>         The index of the eigenvalue to be computed.  I = 1 or I = 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (2)<br>
*>         The original eigenvalues.  We assume D(1) < D(2).<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (2)<br>
*>         The components of the updating vector.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DELTA<br>
*> \verbatim<br>
*>          DELTA is REAL array, dimension (2)<br>
*>         The vector DELTA contains the information necessary<br>
*>         to construct the eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         The scalar in the symmetric updating formula.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DLAM<br>
*> \verbatim<br>
*>          DLAM is REAL<br>
*>         The computed lambda_I, the I-th updated eigenvalue.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ren-Cang Li, Computer Science Division, University of California<br>
*>     at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slaed5_(INTEGER I,float[] D,float[] Z,float[] DELTA,REAL RHO,REAL DLAM);
/**
*> \brief \b SLAED6 used by sstedc. Computes one Newton step in solution of the secular equation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED6 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed6.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed6.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed6.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            ORGATI<br>
*       INTEGER            INFO, KNITER<br>
*       REAL               FINIT, RHO, TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( 3 ), Z( 3 )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAED6 computes the positive or negative root (closest to the origin)<br>
*> of<br>
*>                  z(1)        z(2)        z(3)<br>
*> f(x) =   rho + --------- + ---------- + ---------<br>
*>                 d(1)-x      d(2)-x      d(3)-x<br>
*><br>
*> It is assumed that<br>
*><br>
*>       if ORGATI = .true. the root is between d(2) and d(3);<br>
*>       otherwise it is between d(1) and d(2)<br>
*><br>
*> This routine will be called by SLAED4 when necessary. In most cases,<br>
*> the root sought is the smallest in magnitude, though it might not be<br>
*> in some extremely rare situations.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] KNITER<br>
*> \verbatim<br>
*>          KNITER is INTEGER<br>
*>               Refer to SLAED4 for its significance.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ORGATI<br>
*> \verbatim<br>
*>          ORGATI is LOGICAL<br>
*>               If ORGATI is true, the needed root is between d(2) and<br>
*>               d(3); otherwise it is between d(1) and d(2).  See<br>
*>               SLAED4 for further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>               Refer to the equation f(x) above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (3)<br>
*>               D satisfies d(1) < d(2) < d(3).<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (3)<br>
*>               Each of the elements in z must be positive.<br>
*> \endverbatim<br>
*><br>
*> \param[in] FINIT<br>
*> \verbatim<br>
*>          FINIT is REAL<br>
*>               The value of f at 0. It is more accurate than the one<br>
*>               evaluated inside this routine (if someone wants to do<br>
*>               so).<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is REAL<br>
*>               The root of the equation f(x).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>               = 0: successful exit<br>
*>               > 0: if INFO = 1, failure to converge<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  10/02/03: This version has a few statements commented out for thread<br>
*>  safety (machine parameters are computed on each entry). SJH.<br>
*><br>
*>  05/10/06: Modified from a new version of Ren-Cang Li, use<br>
*>     Gragg-Thornton-Warner cubic convergent scheme for better stability.<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ren-Cang Li, Computer Science Division, University of California<br>
*>     at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slaed6_(INTEGER KNITER,LOGICAL ORGATI,REAL RHO,float[] D,float[] Z,REAL FINIT,REAL TAU,INTEGER INFO);
/**
*> \brief \b SLAED7 used by sstedc. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is dense.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED7 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed7.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed7.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed7.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,<br>
*                          LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR,<br>
*                          PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N,<br>
*      $                   QSIZ, TLVLS<br>
*       REAL               RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),<br>
*      $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )<br>
*       REAL               D( * ), GIVNUM( 2, * ), Q( LDQ, * ),<br>
*      $                   QSTORE( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAED7 computes the updated eigensystem of a diagonal<br>
*> matrix after modification by a rank-one symmetric matrix. This<br>
*> routine is used only for the eigenproblem which requires all<br>
*> eigenvalues and optionally eigenvectors of a dense symmetric matrix<br>
*> that has been reduced to tridiagonal form.  SLAED1 handles<br>
*> the case in which all eigenvalues and eigenvectors of a symmetric<br>
*> tridiagonal matrix are desired.<br>
*><br>
*>   T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)<br>
*><br>
*>    where Z = Q**Tu, u is a vector of length N with ones in the<br>
*>    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.<br>
*><br>
*>    The eigenvectors of the original matrix are stored in Q, and the<br>
*>    eigenvalues are in D.  The algorithm consists of three stages:<br>
*><br>
*>       The first stage consists of deflating the size of the problem<br>
*>       when there are multiple eigenvalues or if there is a zero in<br>
*>       the Z vector.  For each such occurrence the dimension of the<br>
*>       secular equation problem is reduced by one.  This stage is<br>
*>       performed by the routine SLAED8.<br>
*><br>
*>       The second stage consists of calculating the updated<br>
*>       eigenvalues. This is done by finding the roots of the secular<br>
*>       equation via the routine SLAED4 (as called by SLAED9).<br>
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
*> \param[in] ICOMPQ<br>
*> \verbatim<br>
*>          ICOMPQ is INTEGER<br>
*>          = 0:  Compute eigenvalues only.<br>
*>          = 1:  Compute eigenvectors of original dense symmetric matrix<br>
*>                also.  On entry, Q contains the orthogonal matrix used<br>
*>                to reduce the original matrix to tridiagonal form.<br>
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
*>         The dimension of the orthogonal matrix used to reduce<br>
*>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.<br>
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
*>         0 <= CURLVL <= TLVLS.<br>
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
*>          Q is REAL array, dimension (LDQ, N)<br>
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
*> \param[out] INDXQ<br>
*> \verbatim<br>
*>          INDXQ is INTEGER array, dimension (N)<br>
*>         The permutation which will reintegrate the subproblem just<br>
*>         solved back into sorted order, i.e., D( INDXQ( I = 1, N ) )<br>
*>         will be in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         The subdiagonal element used to create the rank-1<br>
*>         modification.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CUTPNT<br>
*> \verbatim<br>
*>          CUTPNT is INTEGER<br>
*>         Contains the location of the last eigenvalue in the leading<br>
*>         sub-matrix.  min(1,N) <= CUTPNT <= N.<br>
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
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (3*N+2*QSIZ*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (4*N)<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void slaed7_(INTEGER ICOMPQ,INTEGER N,INTEGER QSIZ,INTEGER TLVLS,INTEGER CURLVL,INTEGER CURPBM,float[] D,float[] Q,INTEGER LDQ,int[] INDXQ,REAL RHO,INTEGER CUTPNT,float[] QSTORE,int[] QPTR,int[] PRMPTR,int[] PERM,int[] GIVPTR,int[] GIVCOL,float[] GIVNUM,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLAED8 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original matrix is dense.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED8 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed8.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed8.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed8.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO,<br>
*                          CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR,<br>
*                          GIVCOL, GIVNUM, INDXP, INDX, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N,<br>
*      $                   QSIZ<br>
*       REAL               RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),<br>
*      $                   INDXQ( * ), PERM( * )<br>
*       REAL               D( * ), DLAMDA( * ), GIVNUM( 2, * ),<br>
*      $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAED8 merges the two sets of eigenvalues together into a single<br>
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
*> \param[in] ICOMPQ<br>
*> \verbatim<br>
*>          ICOMPQ is INTEGER<br>
*>          = 0:  Compute eigenvalues only.<br>
*>          = 1:  Compute eigenvectors of original dense symmetric matrix<br>
*>                also.  On entry, Q contains the orthogonal matrix used<br>
*>                to reduce the original matrix to tridiagonal form.<br>
*> \endverbatim<br>
*><br>
*> \param[out] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>         The number of non-deflated eigenvalues, and the order of the<br>
*>         related secular equation.<br>
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
*>         The dimension of the orthogonal matrix used to reduce<br>
*>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry, the eigenvalues of the two submatrices to be<br>
*>         combined.  On exit, the trailing (N-K) updated eigenvalues<br>
*>         (those which were deflated) sorted into increasing order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ,N)<br>
*>         If ICOMPQ = 0, Q is not referenced.  Otherwise,<br>
*>         on entry, Q contains the eigenvectors of the partially solved<br>
*>         system which has been previously updated in matrix<br>
*>         multiplies with other partially solved eigensystems.<br>
*>         On exit, Q contains the trailing (N-K) updated eigenvectors<br>
*>         (those which were deflated) in its last N-K columns.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>         The leading dimension of the array Q.  LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] INDXQ<br>
*> \verbatim<br>
*>          INDXQ is INTEGER array, dimension (N)<br>
*>         The permutation which separately sorts the two sub-problems<br>
*>         in D into ascending order.  Note that elements in the second<br>
*>         half of this permutation must first have CUTPNT added to<br>
*>         their values in order to be accurate.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         On entry, the off-diagonal element associated with the rank-1<br>
*>         cut which originally split the two submatrices which are now<br>
*>         being recombined.<br>
*>         On exit, RHO has been modified to the value required by<br>
*>         SLAED3.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CUTPNT<br>
*> \verbatim<br>
*>          CUTPNT is INTEGER<br>
*>         The location of the last eigenvalue in the leading<br>
*>         sub-matrix.  min(1,N) <= CUTPNT <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (N)<br>
*>         On entry, Z contains the updating vector (the last row of<br>
*>         the first sub-eigenvector matrix and the first row of the<br>
*>         second sub-eigenvector matrix).<br>
*>         On exit, the contents of Z are destroyed by the updating<br>
*>         process.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DLAMDA<br>
*> \verbatim<br>
*>          DLAMDA is REAL array, dimension (N)<br>
*>         A copy of the first K eigenvalues which will be used by<br>
*>         SLAED3 to form the secular equation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q2<br>
*> \verbatim<br>
*>          Q2 is REAL array, dimension (LDQ2,N)<br>
*>         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,<br>
*>         a copy of the first K eigenvectors which will be used by<br>
*>         SLAED7 in a matrix multiply (SGEMM) to update the new<br>
*>         eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ2<br>
*> \verbatim<br>
*>          LDQ2 is INTEGER<br>
*>         The leading dimension of the array Q2.  LDQ2 >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>         The first k values of the final deflation-altered z-vector and<br>
*>         will be passed to SLAED3.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PERM<br>
*> \verbatim<br>
*>          PERM is INTEGER array, dimension (N)<br>
*>         The permutations (from deflation and sorting) to be applied<br>
*>         to each eigenblock.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVPTR<br>
*> \verbatim<br>
*>          GIVPTR is INTEGER<br>
*>         The number of Givens rotations which took place in this<br>
*>         subproblem.<br>
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
*> \param[out] INDXP<br>
*> \verbatim<br>
*>          INDXP is INTEGER array, dimension (N)<br>
*>         The permutation used to place deflated values of D at the end<br>
*>         of the array.  INDXP(1:K) points to the nondeflated D-values<br>
*>         and INDXP(K+1:N) points to the deflated eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDX<br>
*> \verbatim<br>
*>          INDX is INTEGER array, dimension (N)<br>
*>         The permutation used to sort the contents of D into ascending<br>
*>         order.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void slaed8_(INTEGER ICOMPQ,INTEGER K,INTEGER N,INTEGER QSIZ,float[] D,float[] Q,INTEGER LDQ,int[] INDXQ,REAL RHO,INTEGER CUTPNT,float[] Z,float[] DLAMDA,float[] Q2,INTEGER LDQ2,float[] W,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,float[] GIVNUM,int[] INDXP,int[] INDX,INTEGER INFO);
/**
*> \brief \b SLAED9 used by sstedc. Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is dense.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAED9 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed9.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed9.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed9.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W,<br>
*                          S, LDS, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, KSTART, KSTOP, LDQ, LDS, N<br>
*       REAL               RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ),<br>
*      $                   W( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAED9 finds the roots of the secular equation, as defined by the<br>
*> values in D, Z, and RHO, between KSTART and KSTOP.  It makes the<br>
*> appropriate calls to SLAED4 and then stores the new matrix of<br>
*> eigenvectors for use in calculating the next level of Z vectors.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of terms in the rational function to be solved by<br>
*>          SLAED4.  K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KSTART<br>
*> \verbatim<br>
*>          KSTART is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] KSTOP<br>
*> \verbatim<br>
*>          KSTOP is INTEGER<br>
*>          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP<br>
*>          are to be computed.  1 <= KSTART <= KSTOP <= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of rows and columns in the Q matrix.<br>
*>          N >= K (delation may result in N > K).<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          D(I) contains the updated eigenvalues<br>
*>          for KSTART <= I <= KSTOP.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ,N)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  LDQ >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>          The value of the parameter in the rank one update equation.<br>
*>          RHO >= 0 required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DLAMDA<br>
*> \verbatim<br>
*>          DLAMDA is REAL array, dimension (K)<br>
*>          The first K elements of this array contain the old roots<br>
*>          of the deflated updating problem.  These are the poles<br>
*>          of the secular equation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (K)<br>
*>          The first K elements of this array contain the components<br>
*>          of the deflation-adjusted updating vector.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (LDS, K)<br>
*>          Will contain the eigenvectors of the repaired matrix which<br>
*>          will be stored for subsequent Z vector calculation and<br>
*>          multiplied by the previously accumulated eigenvectors<br>
*>          to update the system.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDS<br>
*> \verbatim<br>
*>          LDS is INTEGER<br>
*>          The leading dimension of S.  LDS >= max( 1, K ).<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void slaed9_(INTEGER K,INTEGER KSTART,INTEGER KSTOP,INTEGER N,float[] D,float[] Q,INTEGER LDQ,REAL RHO,float[] DLAMDA,float[] W,float[] S,INTEGER LDS,INTEGER INFO);
/**
*> \brief \b SLAEDA used by sstedc. Computes the Z vector determining the rank-one modification of the diagonal matrix. Used when the original matrix is dense.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAEDA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaeda.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaeda.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaeda.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,<br>
*                          GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            CURLVL, CURPBM, INFO, N, TLVLS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ),<br>
*      $                   PRMPTR( * ), QPTR( * )<br>
*       REAL               GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAEDA computes the Z vector corresponding to the merge step in the<br>
*> CURLVLth step of the merge process with TLVLS steps for the CURPBMth<br>
*> problem.<br>
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
*> \param[in] PRMPTR<br>
*> \verbatim<br>
*>          PRMPTR is INTEGER array, dimension (N lg N)<br>
*>         Contains a list of pointers which indicate where in PERM a<br>
*>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)<br>
*>         indicates the size of the permutation and incidentally the<br>
*>         size of the full, non-deflated problem.<br>
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
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (N**2)<br>
*>         Contains the square eigenblocks from previous levels, the<br>
*>         starting positions for blocks are given by QPTR.<br>
*> \endverbatim<br>
*><br>
*> \param[in] QPTR<br>
*> \verbatim<br>
*>          QPTR is INTEGER array, dimension (N+2)<br>
*>         Contains a list of pointers which indicate where in Q an<br>
*>         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates<br>
*>         the size of the block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (N)<br>
*>         On output this vector contains the updating vector (the last<br>
*>         row of the first sub-eigenvector matrix and the first row of<br>
*>         the second sub-eigenvector matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ZTEMP<br>
*> \verbatim<br>
*>          ZTEMP is REAL array, dimension (N)<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void slaeda_(INTEGER N,INTEGER TLVLS,INTEGER CURLVL,INTEGER CURPBM,int[] PRMPTR,int[] PERM,int[] GIVPTR,int[] GIVCOL,float[] GIVNUM,float[] Q,int[] QPTR,float[] Z,float[] ZTEMP,INTEGER INFO);
/**
*> \brief \b SLAEIN computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAEIN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaein.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaein.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaein.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAEIN( RIGHTV, NOINIT, N, H, LDH, WR, WI, VR, VI, B,<br>
*                          LDB, WORK, EPS3, SMLNUM, BIGNUM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            NOINIT, RIGHTV<br>
*       INTEGER            INFO, LDB, LDH, N<br>
*       REAL               BIGNUM, EPS3, SMLNUM, WI, WR<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               B( LDB, * ), H( LDH, * ), VI( * ), VR( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAEIN uses inverse iteration to find a right or left eigenvector<br>
*> corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg<br>
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
*>          = .TRUE. : no initial vector supplied in (VR,VI).<br>
*>          = .FALSE.: initial vector supplied in (VR,VI).<br>
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
*>          H is REAL array, dimension (LDH,N)<br>
*>          The upper Hessenberg matrix H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>          The leading dimension of the array H.  LDH >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] WR<br>
*> \verbatim<br>
*>          WR is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] WI<br>
*> \verbatim<br>
*>          WI is REAL<br>
*>          The real and imaginary parts of the eigenvalue of H whose<br>
*>          corresponding right or left eigenvector is to be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VR<br>
*> \verbatim<br>
*>          VR is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VI<br>
*> \verbatim<br>
*>          VI is REAL array, dimension (N)<br>
*>          On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain<br>
*>          a real starting vector for inverse iteration using the real<br>
*>          eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI<br>
*>          must contain the real and imaginary parts of a complex<br>
*>          starting vector for inverse iteration using the complex<br>
*>          eigenvalue (WR,WI); otherwise VR and VI need not be set.<br>
*>          On exit, if WI = 0.0 (real eigenvalue), VR contains the<br>
*>          computed real eigenvector; if WI.ne.0.0 (complex eigenvalue),<br>
*>          VR and VI contain the real and imaginary parts of the<br>
*>          computed complex eigenvector. The eigenvector is normalized<br>
*>          so that the component of largest magnitude has magnitude 1;<br>
*>          here the magnitude of a complex number (x,y) is taken to be<br>
*>          |x| + |y|.<br>
*>          VI is not referenced if WI = 0.0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,N)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= N+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (N)<br>
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
*> \param[in] BIGNUM<br>
*> \verbatim<br>
*>          BIGNUM is REAL<br>
*>          A machine-dependent value close to the overflow threshold.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          = 1:  inverse iteration did not converge; VR is set to the<br>
*>                last iterate, and so is VI if WI.ne.0.0.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaein_(LOGICAL RIGHTV,LOGICAL NOINIT,INTEGER N,float[] H,INTEGER LDH,REAL WR,REAL WI,float[] VR,float[] VI,float[] B,INTEGER LDB,float[] WORK,REAL EPS3,REAL SMLNUM,REAL BIGNUM,INTEGER INFO);
/**
*> \brief \b SLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAEV2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaev2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaev2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaev2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAEV2( A, B, C, RT1, RT2, CS1, SN1 )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               A, B, C, CS1, RT1, RT2, SN1<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix<br>
*>    [  A   B  ]<br>
*>    [  B   C  ].<br>
*> On return, RT1 is the eigenvalue of larger absolute value, RT2 is the<br>
*> eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right<br>
*> eigenvector for RT1, giving the decomposition<br>
*><br>
*>    [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]<br>
*>    [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL<br>
*>          The (1,1) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL<br>
*>          The (1,2) element and the conjugate of the (2,1) element of<br>
*>          the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL<br>
*>          The (2,2) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT1<br>
*> \verbatim<br>
*>          RT1 is REAL<br>
*>          The eigenvalue of larger absolute value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT2<br>
*> \verbatim<br>
*>          RT2 is REAL<br>
*>          The eigenvalue of smaller absolute value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CS1<br>
*> \verbatim<br>
*>          CS1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] SN1<br>
*> \verbatim<br>
*>          SN1 is REAL<br>
*>          The vector (CS1, SN1) is a unit right eigenvector for RT1.<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
	public void slaev2_(REAL A,REAL B,REAL C,REAL RT1,REAL RT2,REAL CS1,REAL SN1);
/**
*> \brief \b SLAEXC swaps adjacent diagonal blocks of a real upper quasi-triangular matrix in Schur canonical form, by an orthogonal similarity transformation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAEXC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaexc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaexc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaexc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            WANTQ<br>
*       INTEGER            INFO, J1, LDQ, LDT, N, N1, N2<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               Q( LDQ, * ), T( LDT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in<br>
*> an upper quasi-triangular matrix T by an orthogonal similarity<br>
*> transformation.<br>
*><br>
*> T must be in Schur canonical form, that is, block upper triangular<br>
*> with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block<br>
*> has its diagonal elemnts equal and its off-diagonal elements of<br>
*> opposite sign.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTQ<br>
*> \verbatim<br>
*>          WANTQ is LOGICAL<br>
*>          = .TRUE. : accumulate the transformation in the matrix Q;<br>
*>          = .FALSE.: do not accumulate the transformation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix T. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] T<br>
*> \verbatim<br>
*>          T is REAL array, dimension (LDT,N)<br>
*>          On entry, the upper quasi-triangular matrix T, in Schur<br>
*>          canonical form.<br>
*>          On exit, the updated matrix T, again in Schur canonical form.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. LDT >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ,N)<br>
*>          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.<br>
*>          On exit, if WANTQ is .TRUE., the updated matrix Q.<br>
*>          If WANTQ is .FALSE., Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.<br>
*>          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] J1<br>
*> \verbatim<br>
*>          J1 is INTEGER<br>
*>          The index of the first row of the first block T11.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N1<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*>          The order of the first block T11. N1 = 0, 1 or 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N2<br>
*> \verbatim<br>
*>          N2 is INTEGER<br>
*>          The order of the second block T22. N2 = 0, 1 or 2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          = 1: the transformed matrix T would be too far from Schur<br>
*>               form; the blocks are not swapped and T and Q are<br>
*>               unchanged.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaexc_(LOGICAL WANTQ,INTEGER N,float[] T,INTEGER LDT,float[] Q,INTEGER LDQ,INTEGER J1,INTEGER N1,INTEGER N2,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLAG2 computes the eigenvalues of a 2-by-2 generalized eigenvalue problem, with scaling as necessary to avoid over-/underflow.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAG2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slag2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slag2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slag2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1,<br>
*                         WR2, WI )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, LDB<br>
*       REAL               SAFMIN, SCALE1, SCALE2, WI, WR1, WR2<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue<br>
*> problem  A - w B, with scaling as necessary to avoid over-/underflow.<br>
*><br>
*> The scaling factor "s" results in a modified eigenvalue equation<br>
*><br>
*>     s A - w B<br>
*><br>
*> where  s  is a non-negative scaling factor chosen so that  w,  w B,<br>
*> and  s A  do not overflow and, if possible, do not underflow, either.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA, 2)<br>
*>          On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm<br>
*>          is less than 1/SAFMIN.  Entries less than<br>
*>          sqrt(SAFMIN)*norm(A) are subject to being treated as zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB, 2)<br>
*>          On entry, the 2 x 2 upper triangular matrix B.  It is<br>
*>          assumed that the one-norm of B is less than 1/SAFMIN.  The<br>
*>          diagonals should be at least sqrt(SAFMIN) times the largest<br>
*>          element of B (in absolute value); if a diagonal is smaller<br>
*>          than that, then  +/- sqrt(SAFMIN) will be used instead of<br>
*>          that diagonal.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SAFMIN<br>
*> \verbatim<br>
*>          SAFMIN is REAL<br>
*>          The smallest positive number s.t. 1/SAFMIN does not<br>
*>          overflow.  (This should always be SLAMCH('S') -- it is an<br>
*>          argument in order to avoid having to call SLAMCH frequently.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE1<br>
*> \verbatim<br>
*>          SCALE1 is REAL<br>
*>          A scaling factor used to avoid over-/underflow in the<br>
*>          eigenvalue equation which defines the first eigenvalue.  If<br>
*>          the eigenvalues are complex, then the eigenvalues are<br>
*>          ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the<br>
*>          exponent range of the machine), SCALE1=SCALE2, and SCALE1<br>
*>          will always be positive.  If the eigenvalues are real, then<br>
*>          the first (real) eigenvalue is  WR1 / SCALE1 , but this may<br>
*>          overflow or underflow, and in fact, SCALE1 may be zero or<br>
*>          less than the underflow threshold if the exact eigenvalue<br>
*>          is sufficiently large.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE2<br>
*> \verbatim<br>
*>          SCALE2 is REAL<br>
*>          A scaling factor used to avoid over-/underflow in the<br>
*>          eigenvalue equation which defines the second eigenvalue.  If<br>
*>          the eigenvalues are complex, then SCALE2=SCALE1.  If the<br>
*>          eigenvalues are real, then the second (real) eigenvalue is<br>
*>          WR2 / SCALE2 , but this may overflow or underflow, and in<br>
*>          fact, SCALE2 may be zero or less than the underflow<br>
*>          threshold if the exact eigenvalue is sufficiently large.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WR1<br>
*> \verbatim<br>
*>          WR1 is REAL<br>
*>          If the eigenvalue is real, then WR1 is SCALE1 times the<br>
*>          eigenvalue closest to the (2,2) element of A B**(-1).  If the<br>
*>          eigenvalue is complex, then WR1=WR2 is SCALE1 times the real<br>
*>          part of the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WR2<br>
*> \verbatim<br>
*>          WR2 is REAL<br>
*>          If the eigenvalue is real, then WR2 is SCALE2 times the<br>
*>          other eigenvalue.  If the eigenvalue is complex, then<br>
*>          WR1=WR2 is SCALE1 times the real part of the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is REAL<br>
*>          If the eigenvalue is real, then WI is zero.  If the<br>
*>          eigenvalue is complex, then WI is SCALE1 times the imaginary<br>
*>          part of the eigenvalues.  WI will always be non-negative.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slag2_(float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL SAFMIN,REAL SCALE1,REAL SCALE2,REAL WR1,REAL WR2,REAL WI);
/**
*> \brief \b SLAG2D converts a single precision matrix to a double precision matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAG2D + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slag2d.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slag2d.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slag2d.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAG2D( M, N, SA, LDSA, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDSA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               SA( LDSA, * )<br>
*       DOUBLE PRECISION   A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAG2D converts a SINGLE PRECISION matrix, SA, to a DOUBLE<br>
*> PRECISION matrix, A.<br>
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
*>          SA is REAL array, dimension (LDSA,N)<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slag2d_(INTEGER M,INTEGER N,float[] SA,INTEGER LDSA,double[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b SLAGS2 computes 2-by-2 orthogonal matrices U, V, and Q, and applies them to matrices A and B such that the rows of the transformed A and B are parallel.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAGS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slags2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slags2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slags2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV,<br>
*                          SNV, CSQ, SNQ )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            UPPER<br>
*       REAL               A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, SNQ,<br>
*      $                   SNU, SNV<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such<br>
*> that if ( UPPER ) then<br>
*><br>
*>           U**T *A*Q = U**T *( A1 A2 )*Q = ( x  0  )<br>
*>                             ( 0  A3 )     ( x  x  )<br>
*> and<br>
*>           V**T*B*Q = V**T *( B1 B2 )*Q = ( x  0  )<br>
*>                            ( 0  B3 )     ( x  x  )<br>
*><br>
*> or if ( .NOT.UPPER ) then<br>
*><br>
*>           U**T *A*Q = U**T *( A1 0  )*Q = ( x  x  )<br>
*>                             ( A2 A3 )     ( 0  x  )<br>
*> and<br>
*>           V**T*B*Q = V**T*( B1 0  )*Q = ( x  x  )<br>
*>                           ( B2 B3 )     ( 0  x  )<br>
*><br>
*> The rows of the transformed A and B are parallel, where<br>
*><br>
*>   U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ )<br>
*>       ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ )<br>
*><br>
*> Z**T denotes the transpose of Z.<br>
*><br>
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
*>          A2 is REAL<br>
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
*>          B2 is REAL<br>
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
*>          SNU is REAL<br>
*>          The desired orthogonal matrix U.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSV<br>
*> \verbatim<br>
*>          CSV is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNV<br>
*> \verbatim<br>
*>          SNV is REAL<br>
*>          The desired orthogonal matrix V.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSQ<br>
*> \verbatim<br>
*>          CSQ is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNQ<br>
*> \verbatim<br>
*>          SNQ is REAL<br>
*>          The desired orthogonal matrix Q.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slags2_(LOGICAL UPPER,REAL A1,REAL A2,REAL A3,REAL B1,REAL B2,REAL B3,REAL CSU,REAL SNU,REAL CSV,REAL SNV,REAL CSQ,REAL SNQ);
/**
*> \brief \b SLAGTF computes an LU factorization of a matrix T-λI, where T is a general tridiagonal matrix, and λ a scalar, using partial pivoting with row interchanges.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAGTF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagtf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagtf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagtf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, N<br>
*       REAL               LAMBDA, TOL<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IN( * )<br>
*       REAL               A( * ), B( * ), C( * ), D( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAGTF factorizes the matrix (T - lambda*I), where T is an n by n<br>
*> tridiagonal matrix and lambda is a scalar, as<br>
*><br>
*>    T - lambda*I = PLU,<br>
*><br>
*> where P is a permutation matrix, L is a unit lower tridiagonal matrix<br>
*> with at most one non-zero sub-diagonal elements per column and U is<br>
*> an upper triangular matrix with at most two non-zero super-diagonal<br>
*> elements per column.<br>
*><br>
*> The factorization is obtained by Gaussian elimination with partial<br>
*> pivoting and implicit row scaling.<br>
*><br>
*> The parameter LAMBDA is included in the routine so that SLAGTF may<br>
*> be used, in conjunction with SLAGTS, to obtain eigenvectors of T by<br>
*> inverse iteration.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (N)<br>
*>          On entry, A must contain the diagonal elements of T.<br>
*><br>
*>          On exit, A is overwritten by the n diagonal elements of the<br>
*>          upper triangular matrix U of the factorization of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LAMBDA<br>
*> \verbatim<br>
*>          LAMBDA is REAL<br>
*>          On entry, the scalar lambda.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (N-1)<br>
*>          On entry, B must contain the (n-1) super-diagonal elements of<br>
*>          T.<br>
*><br>
*>          On exit, B is overwritten by the (n-1) super-diagonal<br>
*>          elements of the matrix U of the factorization of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N-1)<br>
*>          On entry, C must contain the (n-1) sub-diagonal elements of<br>
*>          T.<br>
*><br>
*>          On exit, C is overwritten by the (n-1) sub-diagonal elements<br>
*>          of the matrix L of the factorization of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TOL<br>
*> \verbatim<br>
*>          TOL is REAL<br>
*>          On entry, a relative tolerance used to indicate whether or<br>
*>          not the matrix (T - lambda*I) is nearly singular. TOL should<br>
*>          normally be chose as approximately the largest relative error<br>
*>          in the elements of T. For example, if the elements of T are<br>
*>          correct to about 4 significant figures, then TOL should be<br>
*>          set to about 5*10**(-4). If TOL is supplied as less than eps,<br>
*>          where eps is the relative machine precision, then the value<br>
*>          eps is used in place of TOL.<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N-2)<br>
*>          On exit, D is overwritten by the (n-2) second super-diagonal<br>
*>          elements of the matrix U of the factorization of T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IN<br>
*> \verbatim<br>
*>          IN is INTEGER array, dimension (N)<br>
*>          On exit, IN contains details of the permutation matrix P. If<br>
*>          an interchange occurred at the kth step of the elimination,<br>
*>          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)<br>
*>          returns the smallest positive integer j such that<br>
*><br>
*>             abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,<br>
*><br>
*>          where norm( A(j) ) denotes the sum of the absolute values of<br>
*>          the jth row of the matrix A. If no such j exists then IN(n)<br>
*>          is returned as zero. If IN(n) is returned as positive, then a<br>
*>          diagonal element of U is small, indicating that<br>
*>          (T - lambda*I) is singular or nearly singular,<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0   : successful exit<br>
*>          .lt. 0: if INFO = -k, the kth argument had an illegal value<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slagtf_(INTEGER N,float[] A,REAL LAMBDA,float[] B,float[] C,REAL TOL,float[] D,int[] IN,INTEGER INFO);
/**
*> \brief \b SLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matrix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAGTM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagtm.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagtm.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagtm.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA,<br>
*                          B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            LDB, LDX, N, NRHS<br>
*       REAL               ALPHA, BETA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               B( LDB, * ), D( * ), DL( * ), DU( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAGTM performs a matrix-vector product of the form<br>
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
*>          = 'T':  Transpose,    B := alpha * A'* X + beta * B<br>
*>          = 'C':  Conjugate transpose = Transpose<br>
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
*>          DL is REAL array, dimension (N-1)<br>
*>          The (n-1) sub-diagonal elements of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The diagonal elements of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is REAL array, dimension (N-1)<br>
*>          The (n-1) super-diagonal elements of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (LDX,NRHS)<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slagtm_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,REAL ALPHA,float[] DL,float[] D,float[] DU,float[] X,INTEGER LDX,REAL BETA,float[] B,INTEGER LDB);
/**
*> \brief \b SLAGTS solves the system of equations (T-λI)x = y or (T-λI)Tx = y,where T is a general tridiagonal matrix and λ a scalar, using the LU factorization computed by slagtf.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAGTS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagts.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagts.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagts.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, JOB, N<br>
*       REAL               TOL<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IN( * )<br>
*       REAL               A( * ), B( * ), C( * ), D( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAGTS may be used to solve one of the systems of equations<br>
*><br>
*>    (T - lambda*I)*x = y   or   (T - lambda*I)**T*x = y,<br>
*><br>
*> where T is an n by n tridiagonal matrix, for x, following the<br>
*> factorization of (T - lambda*I) as<br>
*><br>
*>    (T - lambda*I) = P*L*U ,<br>
*><br>
*> by routine SLAGTF. The choice of equation to be solved is<br>
*> controlled by the argument JOB, and in each case there is an option<br>
*> to perturb zero or very small diagonal elements of U, this option<br>
*> being intended for use in applications such as inverse iteration.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is INTEGER<br>
*>          Specifies the job to be performed by SLAGTS as follows:<br>
*>          =  1: The equations  (T - lambda*I)x = y  are to be solved,<br>
*>                but diagonal elements of U are not to be perturbed.<br>
*>          = -1: The equations  (T - lambda*I)x = y  are to be solved<br>
*>                and, if overflow would otherwise occur, the diagonal<br>
*>                elements of U are to be perturbed. See argument TOL<br>
*>                below.<br>
*>          =  2: The equations  (T - lambda*I)**Tx = y  are to be solved,<br>
*>                but diagonal elements of U are not to be perturbed.<br>
*>          = -2: The equations  (T - lambda*I)**Tx = y  are to be solved<br>
*>                and, if overflow would otherwise occur, the diagonal<br>
*>                elements of U are to be perturbed. See argument TOL<br>
*>                below.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (N)<br>
*>          On entry, A must contain the diagonal elements of U as<br>
*>          returned from SLAGTF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (N-1)<br>
*>          On entry, B must contain the first super-diagonal elements of<br>
*>          U as returned from SLAGTF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N-1)<br>
*>          On entry, C must contain the sub-diagonal elements of L as<br>
*>          returned from SLAGTF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N-2)<br>
*>          On entry, D must contain the second super-diagonal elements<br>
*>          of U as returned from SLAGTF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IN<br>
*> \verbatim<br>
*>          IN is INTEGER array, dimension (N)<br>
*>          On entry, IN must contain details of the matrix P as returned<br>
*>          from SLAGTF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array, dimension (N)<br>
*>          On entry, the right hand side vector y.<br>
*>          On exit, Y is overwritten by the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] TOL<br>
*> \verbatim<br>
*>          TOL is REAL<br>
*>          On entry, with  JOB .lt. 0, TOL should be the minimum<br>
*>          perturbation to be made to very small diagonal elements of U.<br>
*>          TOL should normally be chosen as about eps*norm(U), where eps<br>
*>          is the relative machine precision, but if TOL is supplied as<br>
*>          non-positive, then it is reset to eps*max( abs( u(i,j) ) ).<br>
*>          If  JOB .gt. 0  then TOL is not referenced.<br>
*><br>
*>          On exit, TOL is changed as described above, only if TOL is<br>
*>          non-positive on entry. Otherwise TOL is unchanged.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0   : successful exit<br>
*>          .lt. 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          .gt. 0: overflow would occur when computing the INFO(th)<br>
*>                  element of the solution vector x. This can only occur<br>
*>                  when JOB is supplied as positive and either means<br>
*>                  that a diagonal element of U is very small, or that<br>
*>                  the elements of the right-hand side vector y are very<br>
*>                  large.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slagts_(INTEGER JOB,INTEGER N,float[] A,float[] B,float[] C,float[] D,int[] IN,float[] Y,REAL TOL,INTEGER INFO);
/**
*> \brief \b SLAGV2 computes the Generalized Schur factorization of a real 2-by-2 matrix pencil (A,B) where B is upper triangular.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAGV2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagv2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagv2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagv2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAGV2( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL,<br>
*                          CSR, SNR )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, LDB<br>
*       REAL               CSL, CSR, SNL, SNR<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), ALPHAI( 2 ), ALPHAR( 2 ),<br>
*      $                   B( LDB, * ), BETA( 2 )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAGV2 computes the Generalized Schur factorization of a real 2-by-2<br>
*> matrix pencil (A,B) where B is upper triangular. This routine<br>
*> computes orthogonal (rotation) matrices given by CSL, SNL and CSR,<br>
*> SNR such that<br>
*><br>
*> 1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0<br>
*>    types), then<br>
*><br>
*>    [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]<br>
*>    [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]<br>
*><br>
*>    [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]<br>
*>    [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],<br>
*><br>
*> 2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,<br>
*>    then<br>
*><br>
*>    [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]<br>
*>    [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]<br>
*><br>
*>    [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]<br>
*>    [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]<br>
*><br>
*>    where b11 >= b22 > 0.<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA, 2)<br>
*>          On entry, the 2 x 2 matrix A.<br>
*>          On exit, A is overwritten by the ``A-part'' of the<br>
*>          generalized Schur form.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          THe leading dimension of the array A.  LDA >= 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB, 2)<br>
*>          On entry, the upper triangular 2 x 2 matrix B.<br>
*>          On exit, B is overwritten by the ``B-part'' of the<br>
*>          generalized Schur form.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          THe leading dimension of the array B.  LDB >= 2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is REAL array, dimension (2)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is REAL array, dimension (2)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is REAL array, dimension (2)<br>
*>          (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the<br>
*>          pencil (A,B), k=1,2, i = sqrt(-1).  Note that BETA(k) may<br>
*>          be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSL<br>
*> \verbatim<br>
*>          CSL is REAL<br>
*>          The cosine of the left rotation matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNL<br>
*> \verbatim<br>
*>          SNL is REAL<br>
*>          The sine of the left rotation matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSR<br>
*> \verbatim<br>
*>          CSR is REAL<br>
*>          The cosine of the right rotation matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNR<br>
*> \verbatim<br>
*>          SNR is REAL<br>
*>          The sine of the right rotation matrix.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void slagv2_(float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHAR,float[] ALPHAI,float[] BETA,REAL CSL,REAL SNL,REAL CSR,REAL SNR);
/**
*> \brief \b SLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAHQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slahqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slahqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slahqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,<br>
*                          ILOZ, IHIZ, Z, LDZ, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLAHQR is an auxiliary routine called by SHSEQR to update the<br>
*>    eigenvalues and Schur decomposition already computed by SHSEQR, by<br>
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
*>          It is assumed that H is already upper quasi-triangular in<br>
*>          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless<br>
*>          ILO = 1). SLAHQR works primarily with the Hessenberg<br>
*>          submatrix in rows and columns ILO to IHI, but applies<br>
*>          transformations to all of H if WANTT is .TRUE..<br>
*>          1 <= ILO <= max(1,IHI); IHI <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is REAL array, dimension (LDH,N)<br>
*>          On entry, the upper Hessenberg matrix H.<br>
*>          On exit, if INFO is zero and if WANTT is .TRUE., H is upper<br>
*>          quasi-triangular in rows and columns ILO:IHI, with any<br>
*>          2-by-2 diagonal blocks in standard form. If INFO is zero<br>
*>          and WANTT is .FALSE., the contents of H are unspecified on<br>
*>          exit.  The output state of H if INFO is nonzero is given<br>
*>          below under the description of INFO.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>          The leading dimension of the array H. LDH >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is REAL array, dimension (N)<br>
*>          The real and imaginary parts, respectively, of the computed<br>
*>          eigenvalues ILO to IHI are stored in the corresponding<br>
*>          elements of WR and WI. If two eigenvalues are computed as a<br>
*>          complex conjugate pair, they are stored in consecutive<br>
*>          elements of WR and WI, say the i-th and (i+1)th, with<br>
*>          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the<br>
*>          eigenvalues are stored in the same order as on the diagonal<br>
*>          of the Schur form returned in H, with WR(i) = H(i,i), and, if<br>
*>          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,<br>
*>          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).<br>
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
*>          Z is REAL array, dimension (LDZ,N)<br>
*>          If WANTZ is .TRUE., on entry Z must contain the current<br>
*>          matrix Z of transformations accumulated by SHSEQR, and on<br>
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
*>          .GT. 0: If INFO = i, SLAHQR failed to compute all the<br>
*>                  eigenvalues ILO to IHI in a total of 30 iterations<br>
*>                  per eigenvalue; elements i+1:ihi of WR and WI<br>
*>                  contain those eigenvalues which have been<br>
*>                  successfully computed.<br>
*><br>
*>                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,<br>
*>                  the remaining unconverged eigenvalues are the<br>
*>                  eigenvalues of the upper Hessenberg matrix rows<br>
*>                  and columns ILO thorugh INFO of the final, output<br>
*>                  value of H.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     02-96 Based on modifications by<br>
*>     David Day, Sandia National Laboratory, USA<br>
*><br>
*>     12-04 Further modifications by<br>
*>     Ralph Byers, University of Kansas, USA<br>
*>     This is a modified version of SLAHQR from LAPACK version 3.0.<br>
*>     It is (1) more robust against overflow and underflow and<br>
*>     (2) adopts the more conservative Ahues & Tisseur stopping<br>
*>     criterion (LAWN 122, 1997).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slahqr_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] WR,float[] WI,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER INFO);
/**
*> \brief \b SLAHR2 reduces the specified number of first columns of a general rectangular matrix A so that elements below the specified subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformation to the unreduced part of A.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAHR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slahr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slahr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slahr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            K, LDA, LDT, LDY, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL              A( LDA, * ), T( LDT, NB ), TAU( NB ),<br>
*      $                   Y( LDY, NB )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAHR2 reduces the first NB columns of A real general n-BY-(n-k+1)<br>
*> matrix A so that elements below the k-th subdiagonal are zero. The<br>
*> reduction is performed by an orthogonal similarity transformation<br>
*> Q**T * A * Q. The routine returns the matrices V and T which determine<br>
*> Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.<br>
*><br>
*> This is an auxiliary routine called by SGEHRD.<br>
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
*>          A is REAL array, dimension (LDA,N-K+1)<br>
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
*>          TAU is REAL array, dimension (NB)<br>
*>          The scalar factors of the elementary reflectors. See Further<br>
*>          Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is REAL array, dimension (LDT,NB)<br>
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
*>          Y is REAL array, dimension (LDY,NB)<br>
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
*> \ingroup realOTHERauxiliary<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in<br>
*>  A(i+k+1:n,i), and tau in TAU(i).<br>
*><br>
*>  The elements of the vectors v together form the (n-k+1)-by-nb matrix<br>
*>  V which is needed, with T and Y, to apply the transformation to the<br>
*>  unreduced part of the matrix, using an update of the form:<br>
*>  A := (I - V*T*V**T) * (A - Y*V**T).<br>
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
	public void slahr2_(INTEGER N,INTEGER K,INTEGER NB,float[] A,INTEGER LDA,float[] TAU,float[] T,INTEGER LDT,float[] Y,INTEGER LDY);
/**
*> \brief \b SLAIC1 applies one step of incremental condition estimation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAIC1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaic1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaic1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaic1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            J, JOB<br>
*       REAL               C, GAMMA, S, SEST, SESTPR<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               W( J ), X( J )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAIC1 applies one step of incremental condition estimation in<br>
*> its simplest version:<br>
*><br>
*> Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j<br>
*> lower triangular matrix L, such that<br>
*>          twonorm(L*x) = sest<br>
*> Then SLAIC1 computes sestpr, s, c such that<br>
*> the vector<br>
*>                 [ s*x ]<br>
*>          xhat = [  c  ]<br>
*> is an approximate singular vector of<br>
*>                 [ L      0  ]<br>
*>          Lhat = [ w**T gamma ]<br>
*> in the sense that<br>
*>          twonorm(Lhat*xhat) = sestpr.<br>
*><br>
*> Depending on JOB, an estimate for the largest or smallest singular<br>
*> value is computed.<br>
*><br>
*> Note that [s c]**T and sestpr**2 is an eigenpair of the system<br>
*><br>
*>     diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]<br>
*>                                           [ gamma ]<br>
*><br>
*> where  alpha =  x**T*w.<br>
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
*>          X is REAL array, dimension (J)<br>
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
*>          W is REAL array, dimension (J)<br>
*>          The j-vector w.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GAMMA<br>
*> \verbatim<br>
*>          GAMMA is REAL<br>
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
*>          S is REAL<br>
*>          Sine needed in forming xhat.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is REAL<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaic1_(INTEGER JOB,INTEGER J,float[] X,REAL SEST,float[] W,REAL GAMMA,REAL SESTPR,REAL S,REAL C);
/**
*> \brief \b SLAISNAN tests input for NaN by comparing two arguments for inequality.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAISNAN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaisnan.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaisnan.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaisnan.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       LOGICAL FUNCTION SLAISNAN( SIN1, SIN2 )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               SIN1, SIN2<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This routine is not for general use.  It exists solely to avoid<br>
*> over-optimization in SISNAN.<br>
*><br>
*> SLAISNAN checks for NaNs by comparing its two arguments for<br>
*> inequality.  NaN is the only floating-point value where NaN != NaN<br>
*> returns .TRUE.  To check for NaNs, pass the same variable as both<br>
*> arguments.<br>
*><br>
*> A compiler must assume that the two arguments are<br>
*> not the same variable, and the test will not be optimized away.<br>
*> Interprocedural or whole-program optimization may delete this<br>
*> test.  The ISNAN functions will be replaced by the correct<br>
*> Fortran 03 intrinsic once the intrinsic is widely available.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIN1<br>
*> \verbatim<br>
*>          SIN1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIN2<br>
*> \verbatim<br>
*>          SIN2 is REAL<br>
*>          Two numbers to compare for inequality.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public boolean slaisnan_(REAL SIN1,REAL SIN2);
/**
*> \brief \b SLALN2 solves a 1-by-1 or 2-by-2 linear system of equations of the specified form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLALN2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaln2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaln2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaln2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,<br>
*                          LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            LTRANS<br>
*       INTEGER            INFO, LDA, LDB, LDX, NA, NW<br>
*       REAL               CA, D1, D2, SCALE, SMIN, WI, WR, XNORM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), B( LDB, * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLALN2 solves a system of the form  (ca A - w D ) X = s B<br>
*> or (ca A**T - w D) X = s B   with possible scaling ("s") and<br>
*> perturbation of A.  (A**T means A-transpose.)<br>
*><br>
*> A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA<br>
*> real diagonal matrix, w is a real or complex value, and X and B are<br>
*> NA x 1 matrices -- real if w is real, complex if w is complex.  NA<br>
*> may be 1 or 2.<br>
*><br>
*> If w is complex, X and B are represented as NA x 2 matrices,<br>
*> the first column of each being the real part and the second<br>
*> being the imaginary part.<br>
*><br>
*> "s" is a scaling factor (.LE. 1), computed by SLALN2, which is<br>
*> so chosen that X can be computed without overflow.  X is further<br>
*> scaled if necessary to assure that norm(ca A - w D)*norm(X) is less<br>
*> than overflow.<br>
*><br>
*> If both singular values of (ca A - w D) are less than SMIN,<br>
*> SMIN*identity will be used instead of (ca A - w D).  If only one<br>
*> singular value is less than SMIN, one element of (ca A - w D) will be<br>
*> perturbed enough to make the smallest singular value roughly SMIN.<br>
*> If both singular values are at least SMIN, (ca A - w D) will not be<br>
*> perturbed.  In any case, the perturbation will be at most some small<br>
*> multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values<br>
*> are computed by infinity-norm approximations, and thus will only be<br>
*> correct to a factor of 2 or so.<br>
*><br>
*> Note: all input quantities are assumed to be smaller than overflow<br>
*> by a reasonable factor.  (See BIGNUM.)<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] LTRANS<br>
*> \verbatim<br>
*>          LTRANS is LOGICAL<br>
*>          =.TRUE.:  A-transpose will be used.<br>
*>          =.FALSE.: A will be used (not transposed.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] NA<br>
*> \verbatim<br>
*>          NA is INTEGER<br>
*>          The size of the matrix A.  It may (only) be 1 or 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NW<br>
*> \verbatim<br>
*>          NW is INTEGER<br>
*>          1 if "w" is real, 2 if "w" is complex.  It may only be 1<br>
*>          or 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SMIN<br>
*> \verbatim<br>
*>          SMIN is REAL<br>
*>          The desired lower bound on the singular values of A.  This<br>
*>          should be a safe distance away from underflow or overflow,<br>
*>          say, between (underflow/machine precision) and  (machine<br>
*>          precision * overflow ).  (See BIGNUM and ULP.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] CA<br>
*> \verbatim<br>
*>          CA is REAL<br>
*>          The coefficient c, which A is multiplied by.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,NA)<br>
*>          The NA x NA matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of A.  It must be at least NA.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D1<br>
*> \verbatim<br>
*>          D1 is REAL<br>
*>          The 1,1 element in the diagonal matrix D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D2<br>
*> \verbatim<br>
*>          D2 is REAL<br>
*>          The 2,2 element in the diagonal matrix D.  Not used if NW=1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NW)<br>
*>          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is<br>
*>          complex), column 1 contains the real part of B and column 2<br>
*>          contains the imaginary part.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of B.  It must be at least NA.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WR<br>
*> \verbatim<br>
*>          WR is REAL<br>
*>          The real part of the scalar "w".<br>
*> \endverbatim<br>
*><br>
*> \param[in] WI<br>
*> \verbatim<br>
*>          WI is REAL<br>
*>          The imaginary part of the scalar "w".  Not used if NW=1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (LDX,NW)<br>
*>          The NA x NW matrix X (unknowns), as computed by SLALN2.<br>
*>          If NW=2 ("w" is complex), on exit, column 1 will contain<br>
*>          the real part of X and column 2 will contain the imaginary<br>
*>          part.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of X.  It must be at least NA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          The scale factor that B must be multiplied by to insure<br>
*>          that overflow does not occur when computing X.  Thus,<br>
*>          (ca A - w D) X  will be SCALE*B, not B (ignoring<br>
*>          perturbations of A.)  It will be at most 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] XNORM<br>
*> \verbatim<br>
*>          XNORM is REAL<br>
*>          The infinity-norm of X, when X is regarded as an NA x NW<br>
*>          real matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          An error flag.  It will be set to zero if no error occurs,<br>
*>          a negative number if an argument is in error, or a positive<br>
*>          number if  ca A - w D  had to be perturbed.<br>
*>          The possible values are:<br>
*>          = 0: No error occurred, and (ca A - w D) did not have to be<br>
*>                 perturbed.<br>
*>          = 1: (ca A - w D) had to be perturbed to make its smallest<br>
*>               (or only) singular value greater than SMIN.<br>
*>          NOTE: In the interests of speed, this routine does not<br>
*>                check the inputs for errors.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaln2_(LOGICAL LTRANS,INTEGER NA,INTEGER NW,REAL SMIN,REAL CA,float[] A,INTEGER LDA,REAL D1,REAL D2,float[] B,INTEGER LDB,REAL WR,REAL WI,float[] X,INTEGER LDX,REAL SCALE,REAL XNORM,INTEGER INFO);
/**
*> \brief \b SLALS0 applies back multiplying factors in solving the least squares problem using divide and conquer SVD approach. Used by sgelsd.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLALS0 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slals0.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slals0.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slals0.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX,<br>
*                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,<br>
*                          POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL,<br>
*      $                   LDGNUM, NL, NR, NRHS, SQRE<br>
*       REAL               C, S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( LDGCOL, * ), PERM( * )<br>
*       REAL               B( LDB, * ), BX( LDBX, * ), DIFL( * ),<br>
*      $                   DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ),<br>
*      $                   POLES( LDGNUM, * ), WORK( * ), Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLALS0 applies back the multiplying factors of either the left or the<br>
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
*>          B is REAL array, dimension ( LDB, NRHS )<br>
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
*>          BX is REAL array, dimension ( LDBX, NRHS )<br>
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
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension ( K )<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void slals0_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER NRHS,float[] B,INTEGER LDB,float[] BX,INTEGER LDBX,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,float[] GIVNUM,INTEGER LDGNUM,float[] POLES,float[] DIFL,float[] DIFR,float[] Z,INTEGER K,REAL C,REAL S,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLALSA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slalsa.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slalsa.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slalsa.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U,<br>
*                          LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR,<br>
*                          GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS,<br>
*      $                   SMLSIZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),<br>
*      $                   K( * ), PERM( LDGCOL, * )<br>
*       REAL               B( LDB, * ), BX( LDBX, * ), C( * ),<br>
*      $                   DIFL( LDU, * ), DIFR( LDU, * ),<br>
*      $                   GIVNUM( LDU, * ), POLES( LDU, * ), S( * ),<br>
*      $                   U( LDU, * ), VT( LDU, * ), WORK( * ),<br>
*      $                   Z( LDU, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLALSA is an itermediate step in solving the least squares problem<br>
*> by computing the SVD of the coefficient matrix in compact form (The<br>
*> singular vectors are computed as products of simple orthorgonal<br>
*> matrices.).<br>
*><br>
*> If ICOMPQ = 0, SLALSA applies the inverse of the left singular vector<br>
*> matrix of an upper bidiagonal matrix to the right hand side; and if<br>
*> ICOMPQ = 1, SLALSA applies the right singular vector matrix to the<br>
*> right hand side. The singular vector matrices were generated in<br>
*> compact form by SLALSA.<br>
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
*>          B is REAL array, dimension ( LDB, NRHS )<br>
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
*>          BX is REAL array, dimension ( LDBX, NRHS )<br>
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
*>         On entry, VT**T contains the right singular vector matrices of<br>
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
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array.<br>
*>         The dimension must be at least N.<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void slalsa_(INTEGER ICOMPQ,INTEGER SMLSIZ,INTEGER N,INTEGER NRHS,float[] B,INTEGER LDB,float[] BX,INTEGER LDBX,float[] U,INTEGER LDU,float[] VT,int[] K,float[] DIFL,float[] DIFR,float[] Z,float[] POLES,int[] GIVPTR,int[] GIVCOL,INTEGER LDGCOL,int[] PERM,float[] GIVNUM,float[] C,float[] S,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLALSD uses the singular value decomposition of A to solve the least squares problem.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLALSD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slalsd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slalsd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slalsd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND,<br>
*                          RANK, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               B( LDB, * ), D( * ), E( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLALSD uses the singular value decomposition of A to solve the least<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          WORK is REAL array, dimension at least<br>
*>         (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2),<br>
*>         where NLVL = max(0, INT(log_2 (N/(SMLSIZ+1))) + 1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension at least<br>
*>         (3*N*NLVL + 11*N)<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void slalsd_(CHARACTER UPLO,INTEGER SMLSIZ,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,REAL RCOND,INTEGER RANK,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLAMRG creates a permutation list to merge the entries of two independently sorted sets into a single set sorted in ascending order.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAMRG + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slamrg.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slamrg.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slamrg.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAMRG( N1, N2, A, STRD1, STRD2, INDEX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N1, N2, STRD1, STRD2<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            INDEX( * )<br>
*       REAL               A( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAMRG will create a permutation list which will merge the elements<br>
*> of A (which is composed of two independently sorted sets) into a<br>
*> single set which is sorted in ascending order.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N1<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] N2<br>
*> \verbatim<br>
*>          N2 is INTEGER<br>
*>         These arguments contain the respective lengths of the two<br>
*>         sorted lists to be merged.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (N1+N2)<br>
*>         The first N1 elements of A contain a list of numbers which<br>
*>         are sorted in either ascending or descending order.  Likewise<br>
*>         for the final N2 elements.<br>
*> \endverbatim<br>
*><br>
*> \param[in] STRD1<br>
*> \verbatim<br>
*>          STRD1 is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] STRD2<br>
*> \verbatim<br>
*>          STRD2 is INTEGER<br>
*>         These are the strides to be taken through the array A.<br>
*>         Allowable strides are 1 and -1.  They indicate whether a<br>
*>         subset of A is sorted in ascending (STRDx = 1) or descending<br>
*>         (STRDx = -1) order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDEX<br>
*> \verbatim<br>
*>          INDEX is INTEGER array, dimension (N1+N2)<br>
*>         On exit this array will contain a permutation such that<br>
*>         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be<br>
*>         sorted in ascending order.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slamrg_(INTEGER N1,INTEGER N2,float[] A,INTEGER STRD1,INTEGER STRD2,int[] INDEX);
/**
*> \brief \b SLANEG computes the Sturm count.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANEG + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaneg.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaneg.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaneg.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION SLANEG( N, D, LLD, SIGMA, PIVMIN, R )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N, R<br>
*       REAL               PIVMIN, SIGMA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), LLD( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANEG computes the Sturm count, the number of negative pivots<br>
*> encountered while factoring tridiagonal T - sigma I = L D L^T.<br>
*> This implementation works directly on the factors without forming<br>
*> the tridiagonal matrix T.  The Sturm count is also the number of<br>
*> eigenvalues of T less than sigma.<br>
*><br>
*> This routine is called from SLARRB.<br>
*><br>
*> The current routine does not use the PIVMIN parameter but rather<br>
*> requires IEEE-754 propagation of Infinities and NaNs.  This<br>
*> routine also has no input range restrictions but does require<br>
*> default exception handling such that x/0 produces Inf when x is<br>
*> non-zero, and Inf/Inf produces NaN.  For more information, see:<br>
*><br>
*>   Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in<br>
*>   Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on<br>
*>   Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624<br>
*>   (Tech report version in LAWN 172 with the same title.)<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The N diagonal elements of the diagonal matrix D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LLD<br>
*> \verbatim<br>
*>          LLD is REAL array, dimension (N-1)<br>
*>          The (N-1) elements L(i)*L(i)*D(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIGMA<br>
*> \verbatim<br>
*>          SIGMA is REAL<br>
*>          Shift amount in T - sigma I = L D L^T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot in the Sturm sequence.  May be used<br>
*>          when zero pivots are encountered on non-IEEE-754<br>
*>          architectures.<br>
*> \endverbatim<br>
*><br>
*> \param[in] R<br>
*> \verbatim<br>
*>          R is INTEGER<br>
*>          The twist index for the twisted factorization that is used<br>
*>          for the negcount.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Osni Marques, LBNL/NERSC, USA \n<br>
*>     Christof Voemel, University of California, Berkeley, USA \n<br>
*>     Jason Riedy, University of California, Berkeley, USA \n<br>
*><br>
*  =====================================================================<br>
*/
	public int slaneg_(INTEGER N,float[] D,float[] LLD,REAL SIGMA,REAL PIVMIN,INTEGER R);
/**
*> \brief \b SLANGB returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of general band matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANGB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slangb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slangb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slangb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANGB( NORM, N, KL, KU, AB, LDAB,<br>
*                        WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            KL, KU, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANGB  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the element of  largest absolute value  of an<br>
*> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.<br>
*> \endverbatim<br>
*><br>
*> \return SLANGB<br>
*> \verbatim<br>
*><br>
*>    SLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANGB as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANGB is<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*> \ingroup realGBauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slangb_(CHARACTER NORM,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] WORK);
/**
*> \brief \b SLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANGE + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slange.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slange.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slange.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANGE  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> real matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return SLANGE<br>
*> \verbatim<br>
*><br>
*>    SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANGE as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.  When M = 0,<br>
*>          SLANGE is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0.  When N = 0,<br>
*>          SLANGE is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realGEauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slange_(CHARACTER NORM,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b SLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general tridiagonal matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANGT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slangt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slangt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slangt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANGT( NORM, N, DL, D, DU )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), DL( * ), DU( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANGT  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> real tridiagonal matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return SLANGT<br>
*> \verbatim<br>
*><br>
*>    SLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANGT as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANGT is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DL<br>
*> \verbatim<br>
*>          DL is REAL array, dimension (N-1)<br>
*>          The (n-1) sub-diagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The diagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is REAL array, dimension (N-1)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slangt_(CHARACTER NORM,INTEGER N,float[] DL,float[] D,float[] DU);
/**
*> \brief \b SLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANHS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slanhs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slanhs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slanhs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANHS( NORM, N, A, LDA, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANHS  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> Hessenberg matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return SLANHS<br>
*> \verbatim<br>
*><br>
*>    SLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANHS as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANHS is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slanhs_(CHARACTER NORM,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b SLANSB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric band matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANSB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANSB( NORM, UPLO, N, K, AB, LDAB,<br>
*                        WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            K, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANSB  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the element of  largest absolute value  of an<br>
*> n by n symmetric band matrix A,  with k super-diagonals.<br>
*> \endverbatim<br>
*><br>
*> \return SLANSB<br>
*> \verbatim<br>
*><br>
*>    SLANSB = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANSB as described<br>
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
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANSB is<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slansb_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,INTEGER K,float[] AB,INTEGER LDAB,float[] WORK);
/**
*> \brief \b SLANSF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANSF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLANSF( NORM, TRANSR, UPLO, N, A, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, TRANSR, UPLO<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( 0: * ), WORK( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANSF returns the value of the one norm, or the Frobenius norm, or<br>
*> the infinity norm, or the element of largest absolute value of a<br>
*> real symmetric matrix A in RFP format.<br>
*> \endverbatim<br>
*><br>
*> \return SLANSF<br>
*> \verbatim<br>
*><br>
*>    SLANSF = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          NORM is CHARACTER*1<br>
*>          Specifies the value to be returned in SLANSF as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          Specifies whether the RFP format of A is normal or<br>
*>          transposed format.<br>
*>          = 'N':  RFP format is Normal;<br>
*>          = 'T':  RFP format is Transpose.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the RFP matrix A came from<br>
*>           an upper or lower triangular matrix as follows:<br>
*>           = 'U': RFP A came from an upper triangular matrix;<br>
*>           = 'L': RFP A came from a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A. N >= 0. When N = 0, SLANSF is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension ( N*(N+1)/2 );<br>
*>          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')<br>
*>          part of the symmetric matrix A stored in RFP format. See the<br>
*>          "Notes" below for more details.<br>
*>          Unchanged on exit.<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  We first consider Rectangular Full Packed (RFP) Format when N is<br>
*>  even. We give an example where N = 6.<br>
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
*>  the transpose of the first three columns of AP upper.<br>
*>  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first<br>
*>  three columns of AP lower. The upper triangle A(0:2,0:2) consists of<br>
*>  the transpose of the last three columns of AP lower.<br>
*>  This covers the case N even and TRANSR = 'N'.<br>
*><br>
*>         RFP A                   RFP A<br>
*><br>
*>        03 04 05                33 43 53<br>
*>        13 14 15                00 44 54<br>
*>        23 24 25                10 11 55<br>
*>        33 34 35                20 21 22<br>
*>        00 44 45                30 31 32<br>
*>        01 11 55                40 41 42<br>
*>        02 12 22                50 51 52<br>
*><br>
*>  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the<br>
*>  transpose of RFP A above. One therefore gets:<br>
*><br>
*><br>
*>           RFP A                   RFP A<br>
*><br>
*>     03 13 23 33 00 01 02    33 00 10 20 30 40 50<br>
*>     04 14 24 34 44 11 12    43 44 11 21 31 41 51<br>
*>     05 15 25 35 45 55 22    53 54 55 22 32 42 52<br>
*><br>
*><br>
*>  We then consider Rectangular Full Packed (RFP) Format when N is<br>
*>  odd. We give an example where N = 5.<br>
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
*>  the transpose of the first two columns of AP upper.<br>
*>  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first<br>
*>  three columns of AP lower. The upper triangle A(0:1,1:2) consists of<br>
*>  the transpose of the last two columns of AP lower.<br>
*>  This covers the case N odd and TRANSR = 'N'.<br>
*><br>
*>         RFP A                   RFP A<br>
*><br>
*>        02 03 04                00 33 43<br>
*>        12 13 14                10 11 44<br>
*>        22 23 24                20 21 22<br>
*>        00 33 34                30 31 32<br>
*>        01 11 44                40 41 42<br>
*><br>
*>  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the<br>
*>  transpose of RFP A above. One therefore gets:<br>
*><br>
*>           RFP A                   RFP A<br>
*><br>
*>     02 12 22 00 01             00 10 20 30 40 50<br>
*>     03 13 23 33 11             33 11 21 31 41 51<br>
*>     04 14 24 34 44             43 44 22 32 42 52<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public float slansf_(CHARACTER NORM,CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,float[] WORK);
/**
*> \brief \b SLANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric matrix supplied in packed form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANSP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANSP( NORM, UPLO, N, AP, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANSP  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> real symmetric matrix A,  supplied in packed form.<br>
*> \endverbatim<br>
*><br>
*> \return SLANSP<br>
*> \verbatim<br>
*><br>
*>    SLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANSP as described<br>
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
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANSP is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slansp_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,float[] AP,float[] WORK);
/**
*> \brief \b SLANST returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a real symmetric tridiagonal matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANST + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slanst.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slanst.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slanst.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANST( NORM, N, D, E )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANST  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> real symmetric tridiagonal matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return SLANST<br>
*> \verbatim<br>
*><br>
*>    SLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANST as described<br>
*>          above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANST is<br>
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
*>          E is REAL array, dimension (N-1)<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slanst_(CHARACTER NORM,INTEGER N,float[] D,float[] E);
/**
*> \brief \b SLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a real symmetric matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANSY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansy.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansy.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansy.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANSY( NORM, UPLO, N, A, LDA, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM, UPLO<br>
*       INTEGER            LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANSY  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> real symmetric matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return SLANSY<br>
*> \verbatim<br>
*><br>
*>    SLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANSY as described<br>
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
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANSY is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realSYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slansy_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b SLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a triangular band matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANTB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANTB( NORM, UPLO, DIAG, N, K, AB,<br>
*                        LDAB, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            K, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANTB  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the element of  largest absolute value  of an<br>
*> n by n triangular band matrix A,  with ( k + 1 ) diagonals.<br>
*> \endverbatim<br>
*><br>
*> \return SLANTB<br>
*> \verbatim<br>
*><br>
*>    SLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANTB as described<br>
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
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANTB is<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slantb_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,INTEGER K,float[] AB,INTEGER LDAB,float[] WORK);
/**
*> \brief \b SLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a triangular matrix supplied in packed form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANTP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANTP( NORM, UPLO, DIAG, N, AP, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANTP  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> triangular matrix A, supplied in packed form.<br>
*> \endverbatim<br>
*><br>
*> \return SLANTP<br>
*> \verbatim<br>
*><br>
*>    SLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANTP as described<br>
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
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANTP is<br>
*>          set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slantp_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] AP,float[] WORK);
/**
*> \brief \b SLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a trapezoidal or triangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLANTR( NORM, UPLO, DIAG, M, N, A, LDA,<br>
*                        WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANTR  returns the value of the one norm,  or the Frobenius norm, or<br>
*> the  infinity norm,  or the  element of  largest absolute value  of a<br>
*> trapezoidal or triangular matrix A.<br>
*> \endverbatim<br>
*><br>
*> \return SLANTR<br>
*> \verbatim<br>
*><br>
*>    SLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'<br>
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
*>          Specifies the value to be returned in SLANTR as described<br>
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
*>          UPLO = 'U', M <= N.  When M = 0, SLANTR is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.  N >= 0, and if<br>
*>          UPLO = 'L', N <= M.  When N = 0, SLANTR is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slantr_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
/**
*> \brief \b SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric matrix in standard form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLANV2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slanv2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slanv2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slanv2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric<br>
*> matrix in standard form:<br>
*><br>
*>      [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]<br>
*>      [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]<br>
*><br>
*> where either<br>
*> 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or<br>
*> 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex<br>
*> conjugate eigenvalues.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL<br>
*>          On entry, the elements of the input matrix.<br>
*>          On exit, they are overwritten by the elements of the<br>
*>          standardised Schur form.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT1R<br>
*> \verbatim<br>
*>          RT1R is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT1I<br>
*> \verbatim<br>
*>          RT1I is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT2R<br>
*> \verbatim<br>
*>          RT2R is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] RT2I<br>
*> \verbatim<br>
*>          RT2I is REAL<br>
*>          The real and imaginary parts of the eigenvalues. If the<br>
*>          eigenvalues are a complex conjugate pair, RT1I > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] CS<br>
*> \verbatim<br>
*>          CS is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] SN<br>
*> \verbatim<br>
*>          SN is REAL<br>
*>          Parameters of the rotation matrix.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Modified by V. Sima, Research Institute for Informatics, Bucharest,<br>
*>  Romania, to reduce the risk of cancellation errors,<br>
*>  when computing real eigenvalues, and to ensure, if possible, that<br>
*>  abs(RT1R) >= abs(RT2R).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slanv2_(REAL A,REAL B,REAL C,REAL D,REAL RT1R,REAL RT1I,REAL RT2R,REAL RT2I,REAL CS,REAL SN);
/**
*> \brief \b SLAPLL measures the linear dependence of two vectors.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAPLL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapll.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapll.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapll.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAPLL( N, X, INCX, Y, INCY, SSMIN )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, INCY, N<br>
*       REAL               SSMIN<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               X( * ), Y( * )<br>
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
*>          X is REAL array,<br>
*>                         dimension (1+(N-1)*INCX)<br>
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
*>          Y is REAL array,<br>
*>                         dimension (1+(N-1)*INCY)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slapll_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,REAL SSMIN);
/**
*> \brief \b SLAPMR rearranges rows of a matrix as specified by a permutation vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAPMR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapmr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapmr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapmr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAPMR( FORWRD, M, N, X, LDX, K )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            FORWRD<br>
*       INTEGER            LDX, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            K( * )<br>
*       REAL               X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAPMR rearranges the rows of the M by N matrix X as specified<br>
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
*>          X is REAL array, dimension (LDX,N)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slapmr_(LOGICAL FORWRD,INTEGER M,INTEGER N,float[] X,INTEGER LDX,int[] K);
/**
*> \brief \b SLAPMT performs a forward or backward permutation of the columns of a matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAPMT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapmt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapmt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapmt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAPMT( FORWRD, M, N, X, LDX, K )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            FORWRD<br>
*       INTEGER            LDX, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            K( * )<br>
*       REAL               X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAPMT rearranges the columns of the M by N matrix X as specified<br>
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
*>          X is REAL array, dimension (LDX,N)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slapmt_(LOGICAL FORWRD,INTEGER M,INTEGER N,float[] X,INTEGER LDX,int[] K);
/**
*> \brief \b SLAPY2 returns sqrt(x2+y2).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAPY2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLAPY2( X, Y )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               X, Y<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary<br>
*> overflow.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is REAL<br>
*>          X and Y specify the values x and y.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slapy2_(REAL X,REAL Y);
/**
*> \brief \b SLAPY3 returns sqrt(x2+y2+z2).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAPY3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL             FUNCTION SLAPY3( X, Y, Z )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               X, Y, Z<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause<br>
*> unnecessary overflow.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL<br>
*>          X, Y and Z specify the values x, y and z.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public float slapy3_(REAL X,REAL Y,REAL Z);
/**
*> \brief \b SLAQGB scales a general band matrix, using row and column scaling factors computed by sgbequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQGB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqgb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqgb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqgb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,<br>
*                          AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED<br>
*       INTEGER            KL, KU, LDAB, M, N<br>
*       REAL               AMAX, COLCND, ROWCND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), C( * ), R( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAQGB equilibrates a general M by N band matrix A with KL<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*> \ingroup realGBauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaqgb_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b SLAQGE scales a general rectangular matrix, using row and column scaling factors computed by sgeequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQGE + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqge.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqge.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqge.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,<br>
*                          EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED<br>
*       INTEGER            LDA, M, N<br>
*       REAL               AMAX, COLCND, ROWCND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( * ), R( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAQGE equilibrates a general M by N matrix A using the row and<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realGEauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaqge_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b SLAQP2 computes a QR factorization with column pivoting of the matrix block.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQP2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqp2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqp2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqp2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,<br>
*                          WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, M, N, OFFSET<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            JPVT( * )<br>
*       REAL               A( LDA, * ), TAU( * ), VN1( * ), VN2( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAQP2 computes a QR factorization with column pivoting of<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is <br>
*>          the triangular factor obtained; the elements in block <br>
*>          A(OFFSET+1:M,1:N) below the diagonal, together with the <br>
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
*>          TAU is REAL array, dimension (min(M,N))<br>
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
*>          WORK is REAL array, dimension (N)<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slaqp2_(INTEGER M,INTEGER N,INTEGER OFFSET,float[] A,INTEGER LDA,int[] JPVT,float[] TAU,float[] VN1,float[] VN2,float[] WORK);
/**
*> \brief \b SLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by using BLAS level 3.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQPS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqps.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqps.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqps.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1,<br>
*                          VN2, AUXV, F, LDF )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            KB, LDA, LDF, M, N, NB, OFFSET<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            JPVT( * )<br>
*       REAL               A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ),<br>
*      $                   VN1( * ), VN2( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAQPS computes a step of QR factorization with column pivoting<br>
*> of a real M-by-N matrix A by using Blas-3.  It tries to factorize<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          TAU is REAL array, dimension (KB)<br>
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
*>          AUXV is REAL array, dimension (NB)<br>
*>          Auxiliar vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] F<br>
*> \verbatim<br>
*>          F is REAL array, dimension (LDF,NB)<br>
*>          Matrix F**T = L*Y**T*A.<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slaqps_(INTEGER M,INTEGER N,INTEGER OFFSET,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] JPVT,float[] TAU,float[] VN1,float[] VN2,float[] AUXV,float[] F,INTEGER LDF);
/**
*> \brief \b SLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQR0 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr0.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr0.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr0.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,<br>
*                          ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               H( LDH, * ), WI( * ), WORK( * ), WR( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLAQR0 computes the eigenvalues of a Hessenberg matrix H<br>
*>    and, optionally, the matrices T and Z from the Schur decomposition<br>
*>    H = Z T Z**T, where T is an upper quasi-triangular matrix (the<br>
*>    Schur form), and Z is the orthogonal matrix of Schur vectors.<br>
*><br>
*>    Optionally Z may be postmultiplied into an input orthogonal<br>
*>    matrix Q so that this routine can give the Schur factorization<br>
*>    of a matrix A which has been reduced to the Hessenberg form H<br>
*>    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.<br>
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
*>           previous call to SGEBAL, and then passed to SGEHRD when the<br>
*>           matrix output by SGEBAL is reduced to Hessenberg form.<br>
*>           Otherwise, ILO and IHI should be set to 1 and N,<br>
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.<br>
*>           If N = 0, then ILO = 1 and IHI = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is REAL array, dimension (LDH,N)<br>
*>           On entry, the upper Hessenberg matrix H.<br>
*>           On exit, if INFO = 0 and WANTT is .TRUE., then H contains<br>
*>           the upper quasi-triangular matrix T from the Schur<br>
*>           decomposition (the Schur form); 2-by-2 diagonal blocks<br>
*>           (corresponding to complex conjugate pairs of eigenvalues)<br>
*>           are returned in standard form, with H(i,i) = H(i+1,i+1)<br>
*>           and H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and WANTT is<br>
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
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is REAL array, dimension (IHI)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is REAL array, dimension (IHI)<br>
*>           The real and imaginary parts, respectively, of the computed<br>
*>           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)<br>
*>           and WI(ILO:IHI). If two eigenvalues are computed as a<br>
*>           complex conjugate pair, they are stored in consecutive<br>
*>           elements of WR and WI, say the i-th and (i+1)th, with<br>
*>           WI(i) .GT. 0 and WI(i+1) .LT. 0. If WANTT is .TRUE., then<br>
*>           the eigenvalues are stored in the same order as on the<br>
*>           diagonal of the Schur form returned in H, with<br>
*>           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal<br>
*>           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and<br>
*>           WI(i+1) = -WI(i).<br>
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
*>          Z is REAL array, dimension (LDZ,IHI)<br>
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
*>          WORK is REAL array, dimension LWORK<br>
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
*>           If LWORK = -1, then SLAQR0 does a workspace query.<br>
*>           In this case, SLAQR0 checks the input parameters and<br>
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
*>           .GT. 0:  if INFO = i, SLAQR0 failed to compute all of<br>
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
*>                where U is an orthogonal matrix.  The final<br>
*>                value of H is upper Hessenberg and quasi-triangular<br>
*>                in rows and columns INFO+1 through IHI.<br>
*><br>
*>                If INFO .GT. 0 and WANTZ is .TRUE., then on exit<br>
*><br>
*>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)<br>
*>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U<br>
*><br>
*>                where U is the orthogonal matrix in (*) (regard-<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slaqr0_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] WR,float[] WI,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQR1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               SI1, SI2, SR1, SR2<br>
*       INTEGER            LDH, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               H( LDH, * ), V( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>      Given a 2-by-2 or 3-by-3 matrix H, SLAQR1 sets v to a<br>
*>      scalar multiple of the first column of the product<br>
*><br>
*>      (*)  K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I)<br>
*><br>
*>      scaling to avoid overflows and most underflows. It<br>
*>      is assumed that either<br>
*><br>
*>              1) sr1 = sr2 and si1 = -si2<br>
*>          or<br>
*>              2) si1 = si2 = 0.<br>
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
*>          H is REAL array of dimension (LDH,N)<br>
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
*> \param[in] SR1<br>
*> \verbatim<br>
*>          SR1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] SI1<br>
*> \verbatim<br>
*>          SI1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] SR2<br>
*> \verbatim<br>
*>          SR2 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] SI2<br>
*> \verbatim<br>
*>          SI2 is REAL<br>
*>              The shifts in (*).<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is REAL array of dimension N<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slaqr1_(INTEGER N,float[] H,INTEGER LDH,REAL SR1,REAL SI1,REAL SR2,REAL SI2,float[] V);
/**
*> \brief \b SLAQR2 performs the orthogonal similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,<br>
*                          IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,<br>
*                          LDT, NV, WV, LDWV, WORK, LWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,<br>
*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),<br>
*      $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLAQR2 is identical to SLAQR3 except that it avoids<br>
*>    recursion by calling SLAHQR instead of SLAQR4.<br>
*><br>
*>    Aggressive early deflation:<br>
*><br>
*>    This subroutine accepts as input an upper Hessenberg matrix<br>
*>    H and performs an orthogonal similarity transformation<br>
*>    designed to detect and deflate fully converged eigenvalues from<br>
*>    a trailing principal submatrix.  On output H has been over-<br>
*>    written by a new Hessenberg matrix that is a perturbation of<br>
*>    an orthogonal similarity transformation of H.  It is to be<br>
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
*>          so that the quasi-triangular Schur factor may be<br>
*>          computed (in cooperation with the calling subroutine).<br>
*>          If .FALSE., then only enough of H is updated to preserve<br>
*>          the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          If .TRUE., then the orthogonal matrix Z is updated so<br>
*>          so that the orthogonal Schur factor may be computed<br>
*>          (in cooperation with the calling subroutine).<br>
*>          If .FALSE., then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix H and (if WANTZ is .TRUE.) the<br>
*>          order of the orthogonal matrix Z.<br>
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
*>          H is REAL array, dimension (LDH,N)<br>
*>          On input the initial N-by-N section of H stores the<br>
*>          Hessenberg matrix undergoing aggressive early deflation.<br>
*>          On output H has been transformed by an orthogonal<br>
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
*>          Z is REAL array, dimension (LDZ,N)<br>
*>          IF WANTZ is .TRUE., then on output, the orthogonal<br>
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
*> \param[out] SR<br>
*> \verbatim<br>
*>          SR is REAL array, dimension KBOT<br>
*> \endverbatim<br>
*><br>
*> \param[out] SI<br>
*> \verbatim<br>
*>          SI is REAL array, dimension KBOT<br>
*>          On output, the real and imaginary parts of approximate<br>
*>          eigenvalues that may be used for shifts are stored in<br>
*>          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and<br>
*>          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.<br>
*>          The real and imaginary parts of converged eigenvalues<br>
*>          are stored in SR(KBOT-ND+1) through SR(KBOT) and<br>
*>          SI(KBOT-ND+1) through SI(KBOT), respectively.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is REAL array, dimension (LDV,NW)<br>
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
*>          T is REAL array, dimension (LDT,NW)<br>
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
*>          WV is REAL array, dimension (LDWV,NW)<br>
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
*>          WORK is REAL array, dimension LWORK.<br>
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
*>          If LWORK = -1, then a workspace query is assumed; SLAQR2<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slaqr2_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NW,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER NS,INTEGER ND,float[] SR,float[] SI,float[] V,INTEGER LDV,INTEGER NH,float[] T,INTEGER LDT,INTEGER NV,float[] WV,INTEGER LDWV,float[] WORK,INTEGER LWORK);
/**
*> \brief \b SLAQR3 performs the orthogonal similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQR3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,<br>
*                          IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,<br>
*                          LDT, NV, WV, LDWV, WORK, LWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,<br>
*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),<br>
*      $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),<br>
*      $                   Z( LDZ, * )<br>
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
*>    SLAQR3 accepts as input an upper Hessenberg matrix<br>
*>    H and performs an orthogonal similarity transformation<br>
*>    designed to detect and deflate fully converged eigenvalues from<br>
*>    a trailing principal submatrix.  On output H has been over-<br>
*>    written by a new Hessenberg matrix that is a perturbation of<br>
*>    an orthogonal similarity transformation of H.  It is to be<br>
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
*>          so that the quasi-triangular Schur factor may be<br>
*>          computed (in cooperation with the calling subroutine).<br>
*>          If .FALSE., then only enough of H is updated to preserve<br>
*>          the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          If .TRUE., then the orthogonal matrix Z is updated so<br>
*>          so that the orthogonal Schur factor may be computed<br>
*>          (in cooperation with the calling subroutine).<br>
*>          If .FALSE., then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix H and (if WANTZ is .TRUE.) the<br>
*>          order of the orthogonal matrix Z.<br>
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
*>          H is REAL array, dimension (LDH,N)<br>
*>          On input the initial N-by-N section of H stores the<br>
*>          Hessenberg matrix undergoing aggressive early deflation.<br>
*>          On output H has been transformed by an orthogonal<br>
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
*>          Z is REAL array, dimension (LDZ,N)<br>
*>          IF WANTZ is .TRUE., then on output, the orthogonal<br>
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
*> \param[out] SR<br>
*> \verbatim<br>
*>          SR is REAL array, dimension KBOT<br>
*> \endverbatim<br>
*><br>
*> \param[out] SI<br>
*> \verbatim<br>
*>          SI is REAL array, dimension KBOT<br>
*>          On output, the real and imaginary parts of approximate<br>
*>          eigenvalues that may be used for shifts are stored in<br>
*>          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and<br>
*>          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.<br>
*>          The real and imaginary parts of converged eigenvalues<br>
*>          are stored in SR(KBOT-ND+1) through SR(KBOT) and<br>
*>          SI(KBOT-ND+1) through SI(KBOT), respectively.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is REAL array, dimension (LDV,NW)<br>
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
*>          T is REAL array, dimension (LDT,NW)<br>
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
*>          WV is REAL array, dimension (LDWV,NW)<br>
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
*>          WORK is REAL array, dimension LWORK.<br>
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
*>          If LWORK = -1, then a workspace query is assumed; SLAQR3<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slaqr3_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NW,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER NS,INTEGER ND,float[] SR,float[] SI,float[] V,INTEGER LDV,INTEGER NH,float[] T,INTEGER LDT,INTEGER NV,float[] WV,INTEGER LDWV,float[] WORK,INTEGER LWORK);
/**
*> \brief \b SLAQR4 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQR4 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr4.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr4.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr4.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,<br>
*                          ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               H( LDH, * ), WI( * ), WORK( * ), WR( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLAQR4 implements one level of recursion for SLAQR0.<br>
*>    It is a complete implementation of the small bulge multi-shift<br>
*>    QR algorithm.  It may be called by SLAQR0 and, for large enough<br>
*>    deflation window size, it may be called by SLAQR3.  This<br>
*>    subroutine is identical to SLAQR0 except that it calls SLAQR2<br>
*>    instead of SLAQR3.<br>
*><br>
*>    SLAQR4 computes the eigenvalues of a Hessenberg matrix H<br>
*>    and, optionally, the matrices T and Z from the Schur decomposition<br>
*>    H = Z T Z**T, where T is an upper quasi-triangular matrix (the<br>
*>    Schur form), and Z is the orthogonal matrix of Schur vectors.<br>
*><br>
*>    Optionally Z may be postmultiplied into an input orthogonal<br>
*>    matrix Q so that this routine can give the Schur factorization<br>
*>    of a matrix A which has been reduced to the Hessenberg form H<br>
*>    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.<br>
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
*>           previous call to SGEBAL, and then passed to SGEHRD when the<br>
*>           matrix output by SGEBAL is reduced to Hessenberg form.<br>
*>           Otherwise, ILO and IHI should be set to 1 and N,<br>
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.<br>
*>           If N = 0, then ILO = 1 and IHI = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is REAL array, dimension (LDH,N)<br>
*>           On entry, the upper Hessenberg matrix H.<br>
*>           On exit, if INFO = 0 and WANTT is .TRUE., then H contains<br>
*>           the upper quasi-triangular matrix T from the Schur<br>
*>           decomposition (the Schur form); 2-by-2 diagonal blocks<br>
*>           (corresponding to complex conjugate pairs of eigenvalues)<br>
*>           are returned in standard form, with H(i,i) = H(i+1,i+1)<br>
*>           and H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and WANTT is<br>
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
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is REAL array, dimension (IHI)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is REAL array, dimension (IHI)<br>
*>           The real and imaginary parts, respectively, of the computed<br>
*>           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)<br>
*>           and WI(ILO:IHI). If two eigenvalues are computed as a<br>
*>           complex conjugate pair, they are stored in consecutive<br>
*>           elements of WR and WI, say the i-th and (i+1)th, with<br>
*>           WI(i) .GT. 0 and WI(i+1) .LT. 0. If WANTT is .TRUE., then<br>
*>           the eigenvalues are stored in the same order as on the<br>
*>           diagonal of the Schur form returned in H, with<br>
*>           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal<br>
*>           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and<br>
*>           WI(i+1) = -WI(i).<br>
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
*>          Z is REAL array, dimension (LDZ,IHI)<br>
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
*>          WORK is REAL array, dimension LWORK<br>
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
*>           If LWORK = -1, then SLAQR4 does a workspace query.<br>
*>           In this case, SLAQR4 checks the input parameters and<br>
*>           estimates the optimal workspace size for the given<br>
*>           values of N, ILO and IHI.  The estimate is returned<br>
*>           in WORK(1).  No error message related to LWORK is<br>
*>           issued by XERBLA.  Neither H nor Z are accessed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>             =  0:  successful exit<br>
*>           .GT. 0:  if INFO = i, SLAQR4 failed to compute all of<br>
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
*>                where U is a orthogonal matrix.  The final<br>
*>                value of  H is upper Hessenberg and triangular in<br>
*>                rows and columns INFO+1 through IHI.<br>
*><br>
*>                If INFO .GT. 0 and WANTZ is .TRUE., then on exit<br>
*><br>
*>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)<br>
*>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U<br>
*><br>
*>                where U is the orthogonal matrix in (*) (regard-<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slaqr4_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] WR,float[] WI,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SLAQR5 performs a single small-bulge multi-shift QR sweep.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQR5 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr5.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr5.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr5.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS,<br>
*                          SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U,<br>
*                          LDU, NV, WV, LDWV, NH, WH, LDWH )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,<br>
*      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV<br>
*       LOGICAL            WANTT, WANTZ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               H( LDH, * ), SI( * ), SR( * ), U( LDU, * ),<br>
*      $                   V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLAQR5, called by SLAQR0, performs a<br>
*>    single small-bulge multi-shift QR sweep.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTT<br>
*> \verbatim<br>
*>          WANTT is logical scalar<br>
*>             WANTT = .true. if the quasi-triangular Schur factor<br>
*>             is being computed.  WANTT is set to .false. otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is logical scalar<br>
*>             WANTZ = .true. if the orthogonal Schur factor is being<br>
*>             computed.  WANTZ is set to .false. otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KACC22<br>
*> \verbatim<br>
*>          KACC22 is integer with value 0, 1, or 2.<br>
*>             Specifies the computation mode of far-from-diagonal<br>
*>             orthogonal updates.<br>
*>        = 0: SLAQR5 does not accumulate reflections and does not<br>
*>             use matrix-matrix multiply to update far-from-diagonal<br>
*>             matrix entries.<br>
*>        = 1: SLAQR5 accumulates reflections and uses matrix-matrix<br>
*>             multiply to update the far-from-diagonal matrix entries.<br>
*>        = 2: SLAQR5 accumulates reflections, uses matrix-matrix<br>
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
*> \param[in,out] SR<br>
*> \verbatim<br>
*>          SR is REAL array of size (NSHFTS)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SI<br>
*> \verbatim<br>
*>          SI is REAL array of size (NSHFTS)<br>
*>             SR contains the real parts and SI contains the imaginary<br>
*>             parts of the NSHFTS shifts of origin that define the<br>
*>             multi-shift QR sweep.  On output SR and SI may be<br>
*>             reordered.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is REAL array of size (LDH,N)<br>
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
*>          Z is REAL array of size (LDZ,IHIZ)<br>
*>             If WANTZ = .TRUE., then the QR Sweep orthogonal<br>
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
*>          V is REAL array of size (LDV,NSHFTS/2)<br>
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
*>          U is REAL array of size<br>
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
*>          WH is REAL array of size (LDWH,NH)<br>
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
*>          WV is REAL array of size<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slaqr5_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER KACC22,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NSHFTS,float[] SR,float[] SI,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] V,INTEGER LDV,float[] U,INTEGER LDU,INTEGER NV,float[] WV,INTEGER LDWV,INTEGER NH,float[] WH,INTEGER LDWH);
/**
*> \brief \b SLAQSB scales a symmetric/Hermitian band matrix, using scaling factors computed by spbequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQSB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqsb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqsb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqsb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQSB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            KD, LDAB, N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), S( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAQSB equilibrates a symmetric band matrix A using the scaling<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*><br>
*>          On exit, if INFO = 0, the triangular factor U or L from the<br>
*>          Cholesky factorization A = U**T*U or A = L*L**T of the band<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaqsb_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b SLAQSP scales a symmetric/Hermitian matrix in packed storage, using scaling factors computed by sppequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQSP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqsp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqsp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqsp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQSP( UPLO, N, AP, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), S( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAQSP equilibrates a symmetric matrix A using the scaling factors<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaqsp_(CHARACTER UPLO,INTEGER N,float[] AP,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b SLAQSY scales a symmetric/Hermitian matrix, using scaling factors computed by spoequ.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQSY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqsy.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqsy.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqsy.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, UPLO<br>
*       INTEGER            LDA, N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), S( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAQSY equilibrates a symmetric matrix A using the scaling factors<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realSYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaqsy_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
/**
*> \brief \b SLAQTR solves a real quasi-triangular system of equations, or a complex quasi-triangular system of special form, in real arithmetic.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAQTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            LREAL, LTRAN<br>
*       INTEGER            INFO, LDT, N<br>
*       REAL               SCALE, W<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               B( * ), T( LDT, * ), WORK( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAQTR solves the real quasi-triangular system<br>
*><br>
*>              op(T)*p = scale*c,               if LREAL = .TRUE.<br>
*><br>
*> or the complex quasi-triangular systems<br>
*><br>
*>            op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.<br>
*><br>
*> in real arithmetic, where T is upper quasi-triangular.<br>
*> If LREAL = .FALSE., then the first diagonal block of T must be<br>
*> 1 by 1, B is the specially structured matrix<br>
*><br>
*>                B = [ b(1) b(2) ... b(n) ]<br>
*>                    [       w            ]<br>
*>                    [           w        ]<br>
*>                    [              .     ]<br>
*>                    [                 w  ]<br>
*><br>
*> op(A) = A or A**T, A**T denotes the transpose of<br>
*> matrix A.<br>
*><br>
*> On input, X = [ c ].  On output, X = [ p ].<br>
*>               [ d ]                  [ q ]<br>
*><br>
*> This subroutine is designed for the condition number estimation<br>
*> in routine STRSNA.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] LTRAN<br>
*> \verbatim<br>
*>          LTRAN is LOGICAL<br>
*>          On entry, LTRAN specifies the option of conjugate transpose:<br>
*>             = .FALSE.,    op(T+i*B) = T+i*B,<br>
*>             = .TRUE.,     op(T+i*B) = (T+i*B)**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LREAL<br>
*> \verbatim<br>
*>          LREAL is LOGICAL<br>
*>          On entry, LREAL specifies the input matrix structure:<br>
*>             = .FALSE.,    the input is complex<br>
*>             = .TRUE.,     the input is real<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          On entry, N specifies the order of T+i*B. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] T<br>
*> \verbatim<br>
*>          T is REAL array, dimension (LDT,N)<br>
*>          On entry, T contains a matrix in Schur canonical form.<br>
*>          If LREAL = .FALSE., then the first diagonal block of T must<br>
*>          be 1 by 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the matrix T. LDT >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (N)<br>
*>          On entry, B contains the elements to form the matrix<br>
*>          B as described above.<br>
*>          If LREAL = .TRUE., B is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] W<br>
*> \verbatim<br>
*>          W is REAL<br>
*>          On entry, W is the diagonal element of the matrix B.<br>
*>          If LREAL = .TRUE., W is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          On exit, SCALE is the scale factor.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (2*N)<br>
*>          On entry, X contains the right hand side of the system.<br>
*>          On exit, X is overwritten by the solution.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          On exit, INFO is set to<br>
*>             0: successful exit.<br>
*>               1: the some diagonal 1 by 1 block has been perturbed by<br>
*>                  a small number SMIN to keep nonsingularity.<br>
*>               2: the some diagonal 2 by 2 block has been perturbed by<br>
*>                  a small number in SLALN2 to keep nonsingularity.<br>
*>          NOTE: In the interests of speed, this routine does not<br>
*>                check the inputs for errors.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaqtr_(LOGICAL LTRAN,LOGICAL LREAL,INTEGER N,float[] T,INTEGER LDT,float[] B,REAL W,REAL SCALE,float[] X,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLAR1V computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn of the tridiagonal matrix LDLT - λI.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAR1V + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slar1v.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slar1v.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slar1v.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD,<br>
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
*       REAL             Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAR1V computes the (scaled) r-th column of the inverse of<br>
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
*>          Z is REAL array, dimension (N)<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slar1v_(INTEGER N,INTEGER B1,INTEGER BN,REAL LAMBDA,float[] D,float[] L,float[] LD,float[] LLD,REAL PIVMIN,REAL GAPTOL,float[] Z,LOGICAL WANTNC,INTEGER NEGCNT,REAL ZTZ,REAL MINGMA,INTEGER R,int[] ISUPPZ,REAL NRMINV,REAL RESID,REAL RQCORR,float[] WORK);
/**
*> \brief \b SLAR2V applies a vector of plane rotations with real cosines and real sines from both sides to a sequence of 2-by-2 symmetric/Hermitian matrices.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAR2V + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slar2v.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slar2v.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slar2v.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAR2V( N, X, Y, Z, INCX, C, S, INCC )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCC, INCX, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * ), S( * ), X( * ), Y( * ), Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAR2V applies a vector of real plane rotations from both sides to<br>
*> a sequence of 2-by-2 real symmetric matrices, defined by the elements<br>
*> of the vectors x, y and z. For i = 1,2,...,n<br>
*><br>
*>    ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )<br>
*>    ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )<br>
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
*>          X is REAL array,<br>
*>                         dimension (1+(N-1)*INCX)<br>
*>          The vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array,<br>
*>                         dimension (1+(N-1)*INCX)<br>
*>          The vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is REAL array,<br>
*>                         dimension (1+(N-1)*INCX)<br>
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
*>          S is REAL array, dimension (1+(N-1)*INCC)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slar2v_(INTEGER N,float[] X,float[] Y,float[] Z,INTEGER INCX,float[] C,float[] S,INTEGER INCC);
/**
*> \brief \b SLARF applies an elementary reflector to a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE<br>
*       INTEGER            INCV, LDC, M, N<br>
*       REAL               TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( LDC, * ), V( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARF applies a real elementary reflector H to a real m by n matrix<br>
*> C, from either the left or the right. H is represented in the form<br>
*><br>
*>       H = I - tau * v * v**T<br>
*><br>
*> where tau is a real scalar and v is a real vector.<br>
*><br>
*> If tau = 0, then H is taken to be the unit matrix.<br>
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
*>          V is REAL array, dimension<br>
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
*>          TAU is REAL<br>
*>          The value tau in the representation of H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
*>          On entry, the m by n matrix C.<br>
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
*>          WORK is REAL array, dimension<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slarf_(CHARACTER SIDE,INTEGER M,INTEGER N,float[] V,INTEGER INCV,REAL TAU,float[] C,INTEGER LDC,float[] WORK);
/**
*> \brief \b SLARFB applies a block reflector or its transpose to a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARFB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,<br>
*                          T, LDT, C, LDC, WORK, LDWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, SIDE, STOREV, TRANS<br>
*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ),<br>
*      $                   WORK( LDWORK, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARFB applies a real block reflector H or its transpose H**T to a<br>
*> real m by n matrix C, from either the left or the right.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply H or H**T from the Left<br>
*>          = 'R': apply H or H**T from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply H (No transpose)<br>
*>          = 'T': apply H**T (Transpose)<br>
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
*>          V is REAL array, dimension<br>
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
*>          T is REAL array, dimension (LDT,K)<br>
*>          The triangular k by k matrix T in the representation of the<br>
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
*>          C is REAL array, dimension (LDC,N)<br>
*>          On entry, the m by n matrix C.<br>
*>          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.<br>
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
*>          WORK is REAL array, dimension (LDWORK,K)<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slarfb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] C,INTEGER LDC,float[] WORK,INTEGER LDWORK);
/**
*> \brief \b SLARFG generates an elementary reflector (Householder matrix).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARFG + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfg.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfg.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfg.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       REAL               ALPHA, TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARFG generates a real elementary reflector H of order n, such<br>
*> that<br>
*><br>
*>       H * ( alpha ) = ( beta ),   H**T * H = I.<br>
*>           (   x   )   (   0  )<br>
*><br>
*> where alpha and beta are scalars, and x is an (n-1)-element real<br>
*> vector. H is represented in the form<br>
*><br>
*>       H = I - tau * ( 1 ) * ( 1 v**T ) ,<br>
*>                     ( v )<br>
*><br>
*> where tau is a real scalar and v is a real (n-1)-element<br>
*> vector.<br>
*><br>
*> If the elements of x are all zero, then tau = 0 and H is taken to be<br>
*> the unit matrix.<br>
*><br>
*> Otherwise  1 <= tau <= 2.<br>
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
*>          ALPHA is REAL<br>
*>          On entry, the value alpha.<br>
*>          On exit, it is overwritten with the value beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension<br>
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
*>          TAU is REAL<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slarfg_(INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,REAL TAU);
/**
*> \brief \b SLARFGP generates an elementary reflector (Householder matrix) with non-negative beta.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARFGP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfgp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfgp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfgp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARFGP( N, ALPHA, X, INCX, TAU )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       REAL               ALPHA, TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARFGP generates a real elementary reflector H of order n, such<br>
*> that<br>
*><br>
*>       H * ( alpha ) = ( beta ),   H**T * H = I.<br>
*>           (   x   )   (   0  )<br>
*><br>
*> where alpha and beta are scalars, beta is non-negative, and x is<br>
*> an (n-1)-element real vector.  H is represented in the form<br>
*><br>
*>       H = I - tau * ( 1 ) * ( 1 v**T ) ,<br>
*>                     ( v )<br>
*><br>
*> where tau is a real scalar and v is a real (n-1)-element<br>
*> vector.<br>
*><br>
*> If the elements of x are all zero, then tau = 0 and H is taken to be<br>
*> the unit matrix.<br>
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
*>          ALPHA is REAL<br>
*>          On entry, the value alpha.<br>
*>          On exit, it is overwritten with the value beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension<br>
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
*>          TAU is REAL<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slarfgp_(INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,REAL TAU);
/**
*> \brief \b SLARFT forms the triangular factor T of a block reflector H = I - vtvH<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARFT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarft.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarft.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarft.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, STOREV<br>
*       INTEGER            K, LDT, LDV, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               T( LDT, * ), TAU( * ), V( LDV, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARFT forms the triangular factor T of a real block reflector H<br>
*> of order n, which is defined as a product of k elementary reflectors.<br>
*><br>
*> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;<br>
*><br>
*> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.<br>
*><br>
*> If STOREV = 'C', the vector which defines the elementary reflector<br>
*> H(i) is stored in the i-th column of the array V, and<br>
*><br>
*>    H  =  I - V * T * V**T<br>
*><br>
*> If STOREV = 'R', the vector which defines the elementary reflector<br>
*> H(i) is stored in the i-th row of the array V, and<br>
*><br>
*>    H  =  I - V**T * T * V<br>
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
*>          V is REAL array, dimension<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is REAL array, dimension (LDT,K)<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slarft_(CHARACTER DIRECT,CHARACTER STOREV,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] TAU,float[] T,INTEGER LDT);
/**
*> \brief \b SLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling when the reflector has order ≤ 10.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARFX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE<br>
*       INTEGER            LDC, M, N<br>
*       REAL               TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( LDC, * ), V( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARFX applies a real elementary reflector H to a real m by n<br>
*> matrix C, from either the left or the right. H is represented in the<br>
*> form<br>
*><br>
*>       H = I - tau * v * v**T<br>
*><br>
*> where tau is a real scalar and v is a real vector.<br>
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
*>          V is REAL array, dimension (M) if SIDE = 'L'<br>
*>                                     or (N) if SIDE = 'R'<br>
*>          The vector v in the representation of H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is REAL<br>
*>          The value tau in the representation of H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
*>          On entry, the m by n matrix C.<br>
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',<br>
*>          or C * H if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDA >= (1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension<br>
*>                      (N) if SIDE = 'L'<br>
*>                      or (M) if SIDE = 'R'<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slarfx_(CHARACTER SIDE,INTEGER M,INTEGER N,float[] V,REAL TAU,float[] C,INTEGER LDC,float[] WORK);
/**
*> \brief \b SLARGV generates a vector of plane rotations with real cosines and real sines.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slargv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slargv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slargv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARGV( N, X, INCX, Y, INCY, C, INCC )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCC, INCX, INCY, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * ), X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARGV generates a vector of real plane rotations, determined by<br>
*> elements of the real vectors x and y. For i = 1,2,...,n<br>
*><br>
*>    (  c(i)  s(i) ) ( x(i) ) = ( a(i) )<br>
*>    ( -s(i)  c(i) ) ( y(i) ) = (   0  )<br>
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
*>          X is REAL array,<br>
*>                         dimension (1+(N-1)*INCX)<br>
*>          On entry, the vector x.<br>
*>          On exit, x(i) is overwritten by a(i), for i = 1,...,n.<br>
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
*>          Y is REAL array,<br>
*>                         dimension (1+(N-1)*INCY)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slargv_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] C,INTEGER INCC);
/**
*> \brief \b SLARNV returns a vector of random numbers from a uniform or normal distribution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARNV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarnv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarnv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarnv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARNV( IDIST, ISEED, N, X )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IDIST, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISEED( 4 )<br>
*       REAL               X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARNV returns a vector of n random real numbers from a uniform or<br>
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
*>          = 1:  uniform (0,1)<br>
*>          = 2:  uniform (-1,1)<br>
*>          = 3:  normal (0,1)<br>
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
*>          X is REAL array, dimension (N)<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
	public void slarnv_(INTEGER IDIST,int[] ISEED,INTEGER N,float[] X);
/**
*> \brief \b SLARRA computes the splitting points with the specified threshold.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarra.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarra.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarra.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRA( N, D, E, E2, SPLTOL, TNRM,<br>
*                           NSPLIT, ISPLIT, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, N, NSPLIT<br>
*       REAL                SPLTOL, TNRM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISPLIT( * )<br>
*       REAL               D( * ), E( * ), E2( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Compute the splitting points with threshold SPLTOL.<br>
*> SLARRA sets any "small" off-diagonal elements to zero.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix. N > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the N diagonal elements of the tridiagonal<br>
*>          matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
*>          On entry, the first (N-1) entries contain the subdiagonal<br>
*>          elements of the tridiagonal matrix T; E(N) need not be set.<br>
*>          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT,<br>
*>          are set to zero, the other entries of E are untouched.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E2<br>
*> \verbatim<br>
*>          E2 is REAL array, dimension (N)<br>
*>          On entry, the first (N-1) entries contain the SQUARES of the<br>
*>          subdiagonal elements of the tridiagonal matrix T;<br>
*>          E2(N) need not be set.<br>
*>          On exit, the entries E2( ISPLIT( I ) ),<br>
*>          1 <= I <= NSPLIT, have been set to zero<br>
*> \endverbatim<br>
*><br>
*> \param[in] SPLTOL<br>
*> \verbatim<br>
*>          SPLTOL is REAL<br>
*>          The threshold for splitting. Two criteria can be used:<br>
*>          SPLTOL<0 : criterion based on absolute off-diagonal value<br>
*>          SPLTOL>0 : criterion that preserves relative accuracy<br>
*> \endverbatim<br>
*><br>
*> \param[in] TNRM<br>
*> \verbatim<br>
*>          TNRM is REAL<br>
*>          The norm of the matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] NSPLIT<br>
*> \verbatim<br>
*>          NSPLIT is INTEGER<br>
*>          The number of blocks T splits into. 1 <= NSPLIT <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISPLIT<br>
*> \verbatim<br>
*>          ISPLIT is INTEGER array, dimension (N)<br>
*>          The splitting points, at which T breaks up into blocks.<br>
*>          The first block consists of rows/columns 1 to ISPLIT(1),<br>
*>          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),<br>
*>          etc., and the NSPLIT-th consists of rows/columns<br>
*>          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
	public void slarra_(INTEGER N,float[] D,float[] E,float[] E2,REAL SPLTOL,REAL TNRM,INTEGER NSPLIT,int[] ISPLIT,INTEGER INFO);
/**
*> \brief \b SLARRB provides limited bisection to locate eigenvalues for more accuracy.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRB( N, D, LLD, IFIRST, ILAST, RTOL1,<br>
*                          RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK,<br>
*                          PIVMIN, SPDIAM, TWIST, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IFIRST, ILAST, INFO, N, OFFSET, TWIST<br>
*       REAL               PIVMIN, RTOL1, RTOL2, SPDIAM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               D( * ), LLD( * ), W( * ),<br>
*      $                   WERR( * ), WGAP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Given the relatively robust representation(RRR) L D L^T, SLARRB<br>
*> does "limited" bisection to refine the eigenvalues of L D L^T,<br>
*> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial<br>
*> guesses for these eigenvalues are input in W, the corresponding estimate<br>
*> of the error in these guesses and their gaps are input in WERR<br>
*> and WGAP, respectively. During bisection, intervals<br>
*> [left, right] are maintained by storing their mid-points and<br>
*> semi-widths in the arrays W and WERR respectively.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The N diagonal elements of the diagonal matrix D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LLD<br>
*> \verbatim<br>
*>          LLD is REAL array, dimension (N-1)<br>
*>          The (N-1) elements L(i)*L(i)*D(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IFIRST<br>
*> \verbatim<br>
*>          IFIRST is INTEGER<br>
*>          The index of the first eigenvalue to be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILAST<br>
*> \verbatim<br>
*>          ILAST is INTEGER<br>
*>          The index of the last eigenvalue to be computed.<br>
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
*>          Tolerance for the convergence of the bisection intervals.<br>
*>          An interval [LEFT,RIGHT] has converged if<br>
*>          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )<br>
*>          where GAP is the (estimated) distance to the nearest<br>
*>          eigenvalue.<br>
*> \endverbatim<br>
*><br>
*> \param[in] OFFSET<br>
*> \verbatim<br>
*>          OFFSET is INTEGER<br>
*>          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET<br>
*>          through ILAST-OFFSET elements of these arrays are to be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are<br>
*>          estimates of the eigenvalues of L D L^T indexed IFIRST throug<br>
*>          ILAST.<br>
*>          On output, these estimates are refined.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] WGAP<br>
*> \verbatim<br>
*>          WGAP is REAL array, dimension (N-1)<br>
*>          On input, the (estimated) gaps between consecutive<br>
*>          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between<br>
*>          eigenvalues I and I+1. Note that if IFIRST.EQ.ILAST<br>
*>          then WGAP(IFIRST-OFFSET) must be set to ZERO.<br>
*>          On output, these gaps are refined.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] WERR<br>
*> \verbatim<br>
*>          WERR is REAL array, dimension (N)<br>
*>          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are<br>
*>          the errors in the estimates of the corresponding elements in W.<br>
*>          On output, these errors are refined.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (2*N)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot in the Sturm sequence.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SPDIAM<br>
*> \verbatim<br>
*>          SPDIAM is REAL<br>
*>          The spectral diameter of the matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TWIST<br>
*> \verbatim<br>
*>          TWIST is INTEGER<br>
*>          The twist index for the twisted factorization that is used<br>
*>          for the negcount.<br>
*>          TWIST = N: Compute negcount from L D L^T - LAMBDA I = L+ D+ L+^T<br>
*>          TWIST = 1: Compute negcount from L D L^T - LAMBDA I = U- D- U-^T<br>
*>          TWIST = R: Compute negcount from L D L^T - LAMBDA I = N(r) D(r) N(r)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          Error flag.<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
	public void slarrb_(INTEGER N,float[] D,float[] LLD,INTEGER IFIRST,INTEGER ILAST,REAL RTOL1,REAL RTOL2,INTEGER OFFSET,float[] W,float[] WGAP,float[] WERR,float[] WORK,int[] IWORK,REAL PIVMIN,REAL SPDIAM,INTEGER TWIST,INTEGER INFO);
/**
*> \brief \b SLARRC computes the number of eigenvalues of the symmetric tridiagonal matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRC( JOBT, N, VL, VU, D, E, PIVMIN,<br>
*                                   EIGCNT, LCNT, RCNT, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBT<br>
*       INTEGER            EIGCNT, INFO, LCNT, N, RCNT<br>
*       REAL               PIVMIN, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Find the number of eigenvalues of the symmetric tridiagonal matrix T<br>
*> that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T<br>
*> if JOBT = 'L'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBT<br>
*> \verbatim<br>
*>          JOBT is CHARACTER*1<br>
*>          = 'T':  Compute Sturm count for matrix T.<br>
*>          = 'L':  Compute Sturm count for matrix L D L^T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix. N > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*>          The lower bound for the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
*>          The upper bound for the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T.<br>
*>          JOBT = 'L': The N diagonal elements of the diagonal matrix D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
*>          JOBT = 'T': The N-1 offdiagonal elements of the matrix T.<br>
*>          JOBT = 'L': The N-1 offdiagonal elements of the matrix L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot in the Sturm sequence for T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] EIGCNT<br>
*> \verbatim<br>
*>          EIGCNT is INTEGER<br>
*>          The number of eigenvalues of the symmetric tridiagonal matrix T<br>
*>          that are in the interval (VL,VU]<br>
*> \endverbatim<br>
*><br>
*> \param[out] LCNT<br>
*> \verbatim<br>
*>          LCNT is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCNT<br>
*> \verbatim<br>
*>          RCNT is INTEGER<br>
*>          The left and right negcounts of the interval.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
	public void slarrc_(CHARACTER JOBT,INTEGER N,REAL VL,REAL VU,float[] D,float[] E,REAL PIVMIN,INTEGER EIGCNT,INTEGER LCNT,INTEGER RCNT,INTEGER INFO);
/**
*> \brief \b SLARRD computes the eigenvalues of a symmetric tridiagonal matrix to suitable accuracy.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRD( RANGE, ORDER, N, VL, VU, IL, IU, GERS,<br>
*                           RELTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT,<br>
*                           M, W, WERR, WL, WU, IBLOCK, INDEXW,<br>
*                           WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          ORDER, RANGE<br>
*       INTEGER            IL, INFO, IU, M, N, NSPLIT<br>
*       REAL                PIVMIN, RELTOL, VL, VU, WL, WU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IBLOCK( * ), INDEXW( * ),<br>
*      $                   ISPLIT( * ), IWORK( * )<br>
*       REAL               D( * ), E( * ), E2( * ),<br>
*      $                   GERS( * ), W( * ), WERR( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARRD computes the eigenvalues of a symmetric tridiagonal<br>
*> matrix T to suitable accuracy. This is an auxiliary code to be<br>
*> called from SSTEMR.<br>
*> The user may ask for all eigenvalues, all eigenvalues<br>
*> in the half-open interval (VL, VU], or the IL-th through IU-th<br>
*> eigenvalues.<br>
*><br>
*> To avoid overflow, the matrix must be scaled so that its<br>
*> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest<br>
*> accuracy, it should not be much smaller than that.<br>
*><br>
*> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal<br>
*> Matrix", Report CS41, Computer Science Dept., Stanford<br>
*> University, July 21, 1966.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] RANGE<br>
*> \verbatim<br>
*>          RANGE is CHARACTER*1<br>
*>          = 'A': ("All")   all eigenvalues will be found.<br>
*>          = 'V': ("Value") all eigenvalues in the half-open interval<br>
*>                           (VL, VU] will be found.<br>
*>          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the<br>
*>                           entire matrix) will be found.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ORDER<br>
*> \verbatim<br>
*>          ORDER is CHARACTER*1<br>
*>          = 'B': ("By Block") the eigenvalues will be grouped by<br>
*>                              split-off block (see IBLOCK, ISPLIT) and<br>
*>                              ordered from smallest to largest within<br>
*>                              the block.<br>
*>          = 'E': ("Entire matrix")<br>
*>                              the eigenvalues for the entire matrix<br>
*>                              will be ordered from smallest to<br>
*>                              largest.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the tridiagonal matrix T.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues.  Eigenvalues less than or equal<br>
*>          to VL, or greater than VU, will not be returned.  VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for eigenvalues.  Eigenvalues less than or equal<br>
*>          to VL, or greater than VU, will not be returned.  VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          smallest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GERS<br>
*> \verbatim<br>
*>          GERS is REAL array, dimension (2*N)<br>
*>          The N Gerschgorin intervals (the i-th Gerschgorin interval<br>
*>          is (GERS(2*i-1), GERS(2*i)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] RELTOL<br>
*> \verbatim<br>
*>          RELTOL is REAL<br>
*>          The minimum relative width of an interval.  When an interval<br>
*>          is narrower than RELTOL times the larger (in<br>
*>          magnitude) endpoint, then it is considered to be<br>
*>          sufficiently small, i.e., converged.  Note: this should<br>
*>          always be at least radix*machine epsilon.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The n diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          The (n-1) off-diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E2<br>
*> \verbatim<br>
*>          E2 is REAL array, dimension (N-1)<br>
*>          The (n-1) squared off-diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot allowed in the Sturm sequence for T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NSPLIT<br>
*> \verbatim<br>
*>          NSPLIT is INTEGER<br>
*>          The number of diagonal blocks in the matrix T.<br>
*>          1 <= NSPLIT <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ISPLIT<br>
*> \verbatim<br>
*>          ISPLIT is INTEGER array, dimension (N)<br>
*>          The splitting points, at which T breaks up into submatrices.<br>
*>          The first submatrix consists of rows/columns 1 to ISPLIT(1),<br>
*>          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),<br>
*>          etc., and the NSPLIT-th consists of rows/columns<br>
*>          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.<br>
*>          (Only the first NSPLIT elements will actually be used, but<br>
*>          since the user cannot know a priori what value NSPLIT will<br>
*>          have, N words must be reserved for ISPLIT.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The actual number of eigenvalues found. 0 <= M <= N.<br>
*>          (See also the description of INFO=2,3.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>          On exit, the first M elements of W will contain the<br>
*>          eigenvalue approximations. SLARRD computes an interval<br>
*>          I_j = (a_j, b_j] that includes eigenvalue j. The eigenvalue<br>
*>          approximation is given as the interval midpoint<br>
*>          W(j)= ( a_j + b_j)/2. The corresponding error is bounded by<br>
*>          WERR(j) = abs( a_j - b_j)/2<br>
*> \endverbatim<br>
*><br>
*> \param[out] WERR<br>
*> \verbatim<br>
*>          WERR is REAL array, dimension (N)<br>
*>          The error bound on the corresponding eigenvalue approximation<br>
*>          in W.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WL<br>
*> \verbatim<br>
*>          WL is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] WU<br>
*> \verbatim<br>
*>          WU is REAL<br>
*>          The interval (WL, WU] contains all the wanted eigenvalues.<br>
*>          If RANGE='V', then WL=VL and WU=VU.<br>
*>          If RANGE='A', then WL and WU are the global Gerschgorin bounds<br>
*>                        on the spectrum.<br>
*>          If RANGE='I', then WL and WU are computed by SLAEBZ from the<br>
*>                        index range specified.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IBLOCK<br>
*> \verbatim<br>
*>          IBLOCK is INTEGER array, dimension (N)<br>
*>          At each row/column j where E(j) is zero or small, the<br>
*>          matrix T is considered to split into a block diagonal<br>
*>          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which<br>
*>          block (from 1 to the number of blocks) the eigenvalue W(i)<br>
*>          belongs.  (SLARRD may use the remaining N-M elements as<br>
*>          workspace.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDEXW<br>
*> \verbatim<br>
*>          INDEXW is INTEGER array, dimension (N)<br>
*>          The indices of the eigenvalues within each block (submatrix);<br>
*>          for example, INDEXW(i)= j and IBLOCK(i)=k imply that the<br>
*>          i-th eigenvalue W(i) is the j-th eigenvalue in block k.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  some or all of the eigenvalues failed to converge or<br>
*>                were not computed:<br>
*>                =1 or 3: Bisection failed to converge for some<br>
*>                        eigenvalues; these eigenvalues are flagged by a<br>
*>                        negative block number.  The effect is that the<br>
*>                        eigenvalues may not be as accurate as the<br>
*>                        absolute and relative tolerances.  This is<br>
*>                        generally caused by unexpectedly inaccurate<br>
*>                        arithmetic.<br>
*>                =2 or 3: RANGE='I' only: Not all of the eigenvalues<br>
*>                        IL:IU were found.<br>
*>                        Effect: M < IU+1-IL<br>
*>                        Cause:  non-monotonic arithmetic, causing the<br>
*>                                Sturm sequence to be non-monotonic.<br>
*>                        Cure:   recalculate, using RANGE='A', and pick<br>
*>                                out eigenvalues IL:IU.  In some cases,<br>
*>                                increasing the PARAMETER "FUDGE" may<br>
*>                                make things work.<br>
*>                = 4:    RANGE='I', and the Gershgorin interval<br>
*>                        initially used was too small.  No eigenvalues<br>
*>                        were computed.<br>
*>                        Probable cause: your machine has sloppy<br>
*>                                        floating-point arithmetic.<br>
*>                        Cure: Increase the PARAMETER "FUDGE",<br>
*>                              recompile, and try again.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  FUDGE   REAL, default = 2<br>
*>          A "fudge factor" to widen the Gershgorin intervals.  Ideally,<br>
*>          a value of 1 should work, but on machines with sloppy<br>
*>          arithmetic, this needs to be larger.  The default for<br>
*>          publicly released versions should be large enough to handle<br>
*>          the worst machine around.  Note that this has no effect<br>
*>          on accuracy of the solution.<br>
*> \endverbatim<br>
*><br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     W. Kahan, University of California, Berkeley, USA \n<br>
*>     Beresford Parlett, University of California, Berkeley, USA \n<br>
*>     Jim Demmel, University of California, Berkeley, USA \n<br>
*>     Inderjit Dhillon, University of Texas, Austin, USA \n<br>
*>     Osni Marques, LBNL/NERSC, USA \n<br>
*>     Christof Voemel, University of California, Berkeley, USA \n<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slarrd_(CHARACTER RANGE,CHARACTER ORDER,INTEGER N,REAL VL,REAL VU,INTEGER IL,INTEGER IU,float[] GERS,REAL RELTOL,float[] D,float[] E,float[] E2,REAL PIVMIN,INTEGER NSPLIT,int[] ISPLIT,INTEGER M,float[] W,float[] WERR,REAL WL,REAL WU,int[] IBLOCK,int[] INDEXW,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLARRE given the tridiagonal matrix T, sets small off-diagonal elements to zero and for each unreduced block Ti, finds base representations and eigenvalues.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRE + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarre.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarre.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarre.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRE( RANGE, N, VL, VU, IL, IU, D, E, E2,<br>
*                           RTOL1, RTOL2, SPLTOL, NSPLIT, ISPLIT, M,<br>
*                           W, WERR, WGAP, IBLOCK, INDEXW, GERS, PIVMIN,<br>
*                           WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          RANGE<br>
*       INTEGER            IL, INFO, IU, M, N, NSPLIT<br>
*       REAL               PIVMIN, RTOL1, RTOL2, SPLTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * ),<br>
*      $                   INDEXW( * )<br>
*       REAL               D( * ), E( * ), E2( * ), GERS( * ),<br>
*      $                   W( * ),WERR( * ), WGAP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> To find the desired eigenvalues of a given real symmetric<br>
*> tridiagonal matrix T, SLARRE sets any "small" off-diagonal<br>
*> elements to zero, and for each unreduced block T_i, it finds<br>
*> (a) a suitable shift at one end of the block's spectrum,<br>
*> (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and<br>
*> (c) eigenvalues of each L_i D_i L_i^T.<br>
*> The representations and eigenvalues found are then used by<br>
*> SSTEMR to compute the eigenvectors of T.<br>
*> The accuracy varies depending on whether bisection is used to<br>
*> find a few eigenvalues or the dqds algorithm (subroutine SLASQ2) to<br>
*> conpute all and then discard any unwanted one.<br>
*> As an added benefit, SLARRE also outputs the n<br>
*> Gerschgorin intervals for the matrices L_i D_i L_i^T.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] RANGE<br>
*> \verbatim<br>
*>          RANGE is CHARACTER*1<br>
*>          = 'A': ("All")   all eigenvalues will be found.<br>
*>          = 'V': ("Value") all eigenvalues in the half-open interval<br>
*>                           (VL, VU] will be found.<br>
*>          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the<br>
*>                           entire matrix) will be found.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix. N > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound for the eigenvalues.<br>
*>          Eigenvalues less than or equal to VL, or greater than VU,<br>
*>          will not be returned.  VL < VU.<br>
*>          If RANGE='I' or ='A', SLARRE computes bounds on the desired<br>
*>          part of the spectrum.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
*>          If RANGE='V', the upper bound for the eigenvalues.<br>
*>          Eigenvalues less than or equal to VL, or greater than VU,<br>
*>          will not be returned.  VL < VU.<br>
*>          If RANGE='I' or ='A', SLARRE computes bounds on the desired<br>
*>          part of the spectrum.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          smallest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the N diagonal elements of the tridiagonal<br>
*>          matrix T.<br>
*>          On exit, the N diagonal elements of the diagonal<br>
*>          matrices D_i.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
*>          On entry, the first (N-1) entries contain the subdiagonal<br>
*>          elements of the tridiagonal matrix T; E(N) need not be set.<br>
*>          On exit, E contains the subdiagonal elements of the unit<br>
*>          bidiagonal matrices L_i. The entries E( ISPLIT( I ) ),<br>
*>          1 <= I <= NSPLIT, contain the base points sigma_i on output.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E2<br>
*> \verbatim<br>
*>          E2 is REAL array, dimension (N)<br>
*>          On entry, the first (N-1) entries contain the SQUARES of the<br>
*>          subdiagonal elements of the tridiagonal matrix T;<br>
*>          E2(N) need not be set.<br>
*>          On exit, the entries E2( ISPLIT( I ) ),<br>
*>          1 <= I <= NSPLIT, have been set to zero<br>
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
*> \param[in] SPLTOL<br>
*> \verbatim<br>
*>          SPLTOL is REAL<br>
*>          The threshold for splitting.<br>
*> \endverbatim<br>
*><br>
*> \param[out] NSPLIT<br>
*> \verbatim<br>
*>          NSPLIT is INTEGER<br>
*>          The number of blocks T splits into. 1 <= NSPLIT <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISPLIT<br>
*> \verbatim<br>
*>          ISPLIT is INTEGER array, dimension (N)<br>
*>          The splitting points, at which T breaks up into blocks.<br>
*>          The first block consists of rows/columns 1 to ISPLIT(1),<br>
*>          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),<br>
*>          etc., and the NSPLIT-th consists of rows/columns<br>
*>          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The total number of eigenvalues (of all L_i D_i L_i^T)<br>
*>          found.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>          The first M elements contain the eigenvalues. The<br>
*>          eigenvalues of each of the blocks, L_i D_i L_i^T, are<br>
*>          sorted in ascending order ( SLARRE may use the<br>
*>          remaining N-M elements as workspace).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WERR<br>
*> \verbatim<br>
*>          WERR is REAL array, dimension (N)<br>
*>          The error bound on the corresponding eigenvalue in W.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WGAP<br>
*> \verbatim<br>
*>          WGAP is REAL array, dimension (N)<br>
*>          The separation from the right neighbor eigenvalue in W.<br>
*>          The gap is only with respect to the eigenvalues of the same block<br>
*>          as each block has its own representation tree.<br>
*>          Exception: at the right end of a block we store the left gap<br>
*> \endverbatim<br>
*><br>
*> \param[out] IBLOCK<br>
*> \verbatim<br>
*>          IBLOCK is INTEGER array, dimension (N)<br>
*>          The indices of the blocks (submatrices) associated with the<br>
*>          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue<br>
*>          W(i) belongs to the first block from the top, =2 if W(i)<br>
*>          belongs to the second block, etc.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INDEXW<br>
*> \verbatim<br>
*>          INDEXW is INTEGER array, dimension (N)<br>
*>          The indices of the eigenvalues within each block (submatrix);<br>
*>          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the<br>
*>          i-th eigenvalue W(i) is the 10-th eigenvalue in block 2<br>
*> \endverbatim<br>
*><br>
*> \param[out] GERS<br>
*> \verbatim<br>
*>          GERS is REAL array, dimension (2*N)<br>
*>          The N Gerschgorin intervals (the i-th Gerschgorin interval<br>
*>          is (GERS(2*i-1), GERS(2*i)).<br>
*> \endverbatim<br>
*><br>
*> \param[out] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot in the Sturm sequence for T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (6*N)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (5*N)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          > 0:  A problem occurred in SLARRE.<br>
*>          < 0:  One of the called subroutines signaled an internal problem.<br>
*>                Needs inspection of the corresponding parameter IINFO<br>
*>                for further information.<br>
*><br>
*>          =-1:  Problem in SLARRD.<br>
*>          = 2:  No base representation could be found in MAXTRY iterations.<br>
*>                Increasing MAXTRY and recompilation might be a remedy.<br>
*>          =-3:  Problem in SLARRB when computing the refined root<br>
*>                representation for SLASQ2.<br>
*>          =-4:  Problem in SLARRB when preforming bisection on the<br>
*>                desired part of the spectrum.<br>
*>          =-5:  Problem in SLASQ2.<br>
*>          =-6:  Problem in SLASQ2.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The base representations are required to suffer very little<br>
*>  element growth and consequently define all their eigenvalues to<br>
*>  high relative accuracy.<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Beresford Parlett, University of California, Berkeley, USA \n<br>
*>     Jim Demmel, University of California, Berkeley, USA \n<br>
*>     Inderjit Dhillon, University of Texas, Austin, USA \n<br>
*>     Osni Marques, LBNL/NERSC, USA \n<br>
*>     Christof Voemel, University of California, Berkeley, USA \n<br>
*><br>
*  =====================================================================<br>
*/
	public void slarre_(CHARACTER RANGE,INTEGER N,REAL VL,REAL VU,INTEGER IL,INTEGER IU,float[] D,float[] E,float[] E2,REAL RTOL1,REAL RTOL2,REAL SPLTOL,INTEGER NSPLIT,int[] ISPLIT,INTEGER M,float[] W,float[] WERR,float[] WGAP,int[] IBLOCK,int[] INDEXW,float[] GERS,REAL PIVMIN,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLARRF finds a new relatively robust representation such that at least one of the eigenvalues is relatively isolated.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRF( N, D, L, LD, CLSTRT, CLEND,<br>
*                          W, WGAP, WERR,<br>
*                          SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA,<br>
*                          DPLUS, LPLUS, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            CLSTRT, CLEND, INFO, N<br>
*       REAL               CLGAPL, CLGAPR, PIVMIN, SIGMA, SPDIAM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), DPLUS( * ), L( * ), LD( * ),<br>
*      $          LPLUS( * ), W( * ), WGAP( * ), WERR( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Given the initial representation L D L^T and its cluster of close<br>
*> eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ...<br>
*> W( CLEND ), SLARRF finds a new relatively robust representation<br>
*> L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the<br>
*> eigenvalues of L(+) D(+) L(+)^T is relatively isolated.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix (subblock, if the matrix split).<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The N diagonal elements of the diagonal matrix D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is REAL array, dimension (N-1)<br>
*>          The (N-1) subdiagonal elements of the unit bidiagonal<br>
*>          matrix L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LD<br>
*> \verbatim<br>
*>          LD is REAL array, dimension (N-1)<br>
*>          The (N-1) elements L(i)*D(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CLSTRT<br>
*> \verbatim<br>
*>          CLSTRT is INTEGER<br>
*>          The index of the first eigenvalue in the cluster.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CLEND<br>
*> \verbatim<br>
*>          CLEND is INTEGER<br>
*>          The index of the last eigenvalue in the cluster.<br>
*> \endverbatim<br>
*><br>
*> \param[in] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension<br>
*>          dimension is >=  (CLEND-CLSTRT+1)<br>
*>          The eigenvalue APPROXIMATIONS of L D L^T in ascending order.<br>
*>          W( CLSTRT ) through W( CLEND ) form the cluster of relatively<br>
*>          close eigenalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] WGAP<br>
*> \verbatim<br>
*>          WGAP is REAL array, dimension<br>
*>          dimension is >=  (CLEND-CLSTRT+1)<br>
*>          The separation from the right neighbor eigenvalue in W.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WERR<br>
*> \verbatim<br>
*>          WERR is REAL array, dimension<br>
*>          dimension is >=  (CLEND-CLSTRT+1)<br>
*>          WERR contain the semiwidth of the uncertainty<br>
*>          interval of the corresponding eigenvalue APPROXIMATION in W<br>
*> \endverbatim<br>
*><br>
*> \param[in] SPDIAM<br>
*> \verbatim<br>
*>          SPDIAM is REAL<br>
*>          estimate of the spectral diameter obtained from the<br>
*>          Gerschgorin intervals<br>
*> \endverbatim<br>
*><br>
*> \param[in] CLGAPL<br>
*> \verbatim<br>
*>          CLGAPL is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] CLGAPR<br>
*> \verbatim<br>
*>          CLGAPR is REAL<br>
*>          absolute gap on each end of the cluster.<br>
*>          Set by the calling routine to protect against shifts too close<br>
*>          to eigenvalues outside the cluster.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot allowed in the Sturm sequence.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SIGMA<br>
*> \verbatim<br>
*>          SIGMA is REAL<br>
*>          The shift used to form L(+) D(+) L(+)^T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DPLUS<br>
*> \verbatim<br>
*>          DPLUS is REAL array, dimension (N)<br>
*>          The N diagonal elements of the diagonal matrix D(+).<br>
*> \endverbatim<br>
*><br>
*> \param[out] LPLUS<br>
*> \verbatim<br>
*>          LPLUS is REAL array, dimension (N-1)<br>
*>          The first (N-1) elements of LPLUS contain the subdiagonal<br>
*>          elements of the unit bidiagonal matrix L(+).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          Signals processing OK (=0) or failure (=1)<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
	public void slarrf_(INTEGER N,float[] D,float[] L,float[] LD,INTEGER CLSTRT,INTEGER CLEND,float[] W,float[] WGAP,float[] WERR,REAL SPDIAM,REAL CLGAPL,REAL CLGAPR,REAL PIVMIN,REAL SIGMA,float[] DPLUS,float[] LPLUS,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLARRJ performs refinement of the initial estimates of the eigenvalues of the matrix T.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRJ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrj.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrj.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrj.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRJ( N, D, E2, IFIRST, ILAST,<br>
*                          RTOL, OFFSET, W, WERR, WORK, IWORK,<br>
*                          PIVMIN, SPDIAM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IFIRST, ILAST, INFO, N, OFFSET<br>
*       REAL               PIVMIN, RTOL, SPDIAM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               D( * ), E2( * ), W( * ),<br>
*      $                   WERR( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Given the initial eigenvalue approximations of T, SLARRJ<br>
*> does  bisection to refine the eigenvalues of T,<br>
*> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial<br>
*> guesses for these eigenvalues are input in W, the corresponding estimate<br>
*> of the error in these guesses in WERR. During bisection, intervals<br>
*> [left, right] are maintained by storing their mid-points and<br>
*> semi-widths in the arrays W and WERR respectively.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The N diagonal elements of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E2<br>
*> \verbatim<br>
*>          E2 is REAL array, dimension (N-1)<br>
*>          The Squares of the (N-1) subdiagonal elements of T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IFIRST<br>
*> \verbatim<br>
*>          IFIRST is INTEGER<br>
*>          The index of the first eigenvalue to be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILAST<br>
*> \verbatim<br>
*>          ILAST is INTEGER<br>
*>          The index of the last eigenvalue to be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RTOL<br>
*> \verbatim<br>
*>          RTOL is REAL<br>
*>          Tolerance for the convergence of the bisection intervals.<br>
*>          An interval [LEFT,RIGHT] has converged if<br>
*>          RIGHT-LEFT.LT.RTOL*MAX(|LEFT|,|RIGHT|).<br>
*> \endverbatim<br>
*><br>
*> \param[in] OFFSET<br>
*> \verbatim<br>
*>          OFFSET is INTEGER<br>
*>          Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET<br>
*>          through ILAST-OFFSET elements of these arrays are to be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are<br>
*>          estimates of the eigenvalues of L D L^T indexed IFIRST through<br>
*>          ILAST.<br>
*>          On output, these estimates are refined.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] WERR<br>
*> \verbatim<br>
*>          WERR is REAL array, dimension (N)<br>
*>          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are<br>
*>          the errors in the estimates of the corresponding elements in W.<br>
*>          On output, these errors are refined.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (2*N)<br>
*>          Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot in the Sturm sequence for T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SPDIAM<br>
*> \verbatim<br>
*>          SPDIAM is REAL<br>
*>          The spectral diameter of T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          Error flag.<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
	public void slarrj_(INTEGER N,float[] D,float[] E2,INTEGER IFIRST,INTEGER ILAST,REAL RTOL,INTEGER OFFSET,float[] W,float[] WERR,float[] WORK,int[] IWORK,REAL PIVMIN,REAL SPDIAM,INTEGER INFO);
/**
*> \brief \b SLARRK computes one eigenvalue of a symmetric tridiagonal matrix T to suitable accuracy.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrk.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrk.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrk.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRK( N, IW, GL, GU,<br>
*                           D, E2, PIVMIN, RELTOL, W, WERR, INFO)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER   INFO, IW, N<br>
*       REAL                PIVMIN, RELTOL, GL, GU, W, WERR<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E2( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARRK computes one eigenvalue of a symmetric tridiagonal<br>
*> matrix T to suitable accuracy. This is an auxiliary code to be<br>
*> called from SSTEMR.<br>
*><br>
*> To avoid overflow, the matrix must be scaled so that its<br>
*> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest<br>
*> accuracy, it should not be much smaller than that.<br>
*><br>
*> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal<br>
*> Matrix", Report CS41, Computer Science Dept., Stanford<br>
*> University, July 21, 1966.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the tridiagonal matrix T.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IW<br>
*> \verbatim<br>
*>          IW is INTEGER<br>
*>          The index of the eigenvalues to be returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] GL<br>
*> \verbatim<br>
*>          GL is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] GU<br>
*> \verbatim<br>
*>          GU is REAL<br>
*>          An upper and a lower bound on the eigenvalue.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The n diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E2<br>
*> \verbatim<br>
*>          E2 is REAL array, dimension (N-1)<br>
*>          The (n-1) squared off-diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PIVMIN<br>
*> \verbatim<br>
*>          PIVMIN is REAL<br>
*>          The minimum pivot allowed in the Sturm sequence for T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RELTOL<br>
*> \verbatim<br>
*>          RELTOL is REAL<br>
*>          The minimum relative width of an interval.  When an interval<br>
*>          is narrower than RELTOL times the larger (in<br>
*>          magnitude) endpoint, then it is considered to be<br>
*>          sufficiently small, i.e., converged.  Note: this should<br>
*>          always be at least radix*machine epsilon.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] WERR<br>
*> \verbatim<br>
*>          WERR is REAL<br>
*>          The error bound on the corresponding eigenvalue approximation<br>
*>          in W.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:       Eigenvalue converged<br>
*>          = -1:      Eigenvalue did NOT converge<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  FUDGE   REAL            , default = 2<br>
*>          A "fudge factor" to widen the Gershgorin intervals.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slarrk_(INTEGER N,INTEGER IW,REAL GL,REAL GU,float[] D,float[] E2,REAL PIVMIN,REAL RELTOL,REAL W,REAL WERR,INTEGER INFO);
/**
*> \brief \b SLARRR performs tests to decide whether the symmetric tridiagonal matrix T warrants expensive computations which guarantee high relative accuracy in the eigenvalues.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRR( N, D, E, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * )<br>
*       ..<br>
*  <br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Perform tests to decide whether the symmetric tridiagonal matrix T<br>
*> warrants expensive computations which guarantee high relative accuracy<br>
*> in the eigenvalues.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix. N > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The N diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
*>          On entry, the first (N-1) entries contain the subdiagonal<br>
*>          elements of the tridiagonal matrix T; E(N) is set to ZERO.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          INFO = 0(default) : the matrix warrants computations preserving<br>
*>                              relative accuracy.<br>
*>          INFO = 1          : the matrix warrants computations guaranteeing<br>
*>                              only absolute accuracy.<br>
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
*> \ingroup auxOTHERauxiliary<br>
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
	public void slarrr_(INTEGER N,float[] D,float[] E,INTEGER INFO);
/**
*> \brief \b SLARRV computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenvalues of L D LT.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARRV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARRV( N, VL, VU, D, L, PIVMIN,<br>
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
*       REAL              Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARRV computes the eigenvectors of the tridiagonal matrix<br>
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
*>          Z is REAL array, dimension (LDZ, max(1,M) )<br>
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
*>          > 0:  A problem occurred in SLARRV.<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slarrv_(INTEGER N,REAL VL,REAL VU,float[] D,float[] L,REAL PIVMIN,int[] ISPLIT,INTEGER M,INTEGER DOL,INTEGER DOU,REAL MINRGP,REAL RTOL1,REAL RTOL2,float[] W,float[] WERR,float[] WGAP,int[] IBLOCK,int[] INDEXW,float[] GERS,float[] Z,INTEGER LDZ,int[] ISUPPZ,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLARSCL2 performs reciprocal diagonal scaling on a vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARSCL2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarscl2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarscl2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarscl2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARSCL2 ( M, N, D, X, LDX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDX<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARSCL2 performs a reciprocal diagonal scaling on an vector:<br>
*>   x <-- inv(D) * x<br>
*> where the diagonal matrix D is stored as a vector.<br>
*><br>
*> Eventually to be replaced by BLAS_sge_diag_scale in the new BLAS<br>
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
*>          X is REAL array, dimension (LDX,N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slarscl2_(INTEGER M,INTEGER N,float[] D,float[] X,INTEGER LDX);
/**
*> \brief \b SLARTG generates a plane rotation with real cosine and real sine.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARTG + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartg.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartg.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartg.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARTG( F, G, CS, SN, R )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               CS, F, G, R, SN<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARTG generate a plane rotation so that<br>
*><br>
*>    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.<br>
*>    [ -SN  CS  ]     [ G ]     [ 0 ]<br>
*><br>
*> This is a slower, more accurate version of the BLAS1 routine SROTG,<br>
*> with the following other differences:<br>
*>    F and G are unchanged on return.<br>
*>    If G=0, then CS=1 and SN=0.<br>
*>    If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any<br>
*>       floating point operations (saves work in SBDSQR when<br>
*>       there are zeros on the diagonal).<br>
*><br>
*> If F exceeds G in magnitude, CS will be positive.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] F<br>
*> \verbatim<br>
*>          F is REAL<br>
*>          The first component of vector to be rotated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] G<br>
*> \verbatim<br>
*>          G is REAL<br>
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
*>          SN is REAL<br>
*>          The sine of the rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] R<br>
*> \verbatim<br>
*>          R is REAL<br>
*>          The nonzero component of the rotated vector.<br>
*><br>
*>  This version has a few statements commented out for thread safety<br>
*>  (machine parameters are computed on each entry). 10 feb 03, SJH.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slartg_(REAL F,REAL G,REAL CS,REAL SN,REAL R);
/**
*> \brief \b SLARTGP generates a plane rotation so that the diagonal is nonnegative.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARTGP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartgp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartgp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartgp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARTGP( F, G, CS, SN, R )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               CS, F, G, R, SN<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARTGP generates a plane rotation so that<br>
*><br>
*>    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.<br>
*>    [ -SN  CS  ]     [ G ]     [ 0 ]<br>
*><br>
*> This is a slower, more accurate version of the Level 1 BLAS routine SROTG,<br>
*> with the following other differences:<br>
*>    F and G are unchanged on return.<br>
*>    If G=0, then CS=(+/-)1 and SN=0.<br>
*>    If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.<br>
*><br>
*> The sign is chosen so that R >= 0.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] F<br>
*> \verbatim<br>
*>          F is REAL<br>
*>          The first component of vector to be rotated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] G<br>
*> \verbatim<br>
*>          G is REAL<br>
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
*>          SN is REAL<br>
*>          The sine of the rotation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] R<br>
*> \verbatim<br>
*>          R is REAL<br>
*>          The nonzero component of the rotated vector.<br>
*><br>
*>  This version has a few statements commented out for thread safety<br>
*>  (machine parameters are computed on each entry). 10 feb 03, SJH.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slartgp_(REAL F,REAL G,REAL CS,REAL SN,REAL R);
/**
*> \brief \b SLARTGS generates a plane rotation designed to introduce a bulge in implicit QR iteration for the bidiagonal SVD problem.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARTGS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartgs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartgs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartgs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARTGS( X, Y, SIGMA, CS, SN )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL                    CS, SIGMA, SN, X, Y<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARTGS generates a plane rotation designed to introduce a bulge in<br>
*> Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD<br>
*> problem. X and Y are the top-row entries, and SIGMA is the shift.<br>
*> The computed CS and SN define a plane rotation satisfying<br>
*><br>
*>    [  CS  SN  ]  .  [ X^2 - SIGMA ]  =  [ R ],<br>
*>    [ -SN  CS  ]     [    X * Y    ]     [ 0 ]<br>
*><br>
*> with R nonnegative.  If X^2 - SIGMA and X * Y are 0, then the<br>
*> rotation is by PI/2.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL<br>
*>          The (1,1) entry of an upper bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is REAL<br>
*>          The (1,2) entry of an upper bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIGMA<br>
*> \verbatim<br>
*>          SIGMA is REAL<br>
*>          The shift.<br>
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
*>          SN is REAL<br>
*>          The sine of the rotation.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slartgs_(REAL X,REAL Y,REAL SIGMA,REAL CS,REAL SN);
/**
*> \brief \b SLARTV applies a vector of plane rotations with real cosines and real sines to the elements of a pair of vectors.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARTV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARTV( N, X, INCX, Y, INCY, C, S, INCC )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCC, INCX, INCY, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( * ), S( * ), X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARTV applies a vector of real plane rotations to elements of the<br>
*> real vectors x and y. For i = 1,2,...,n<br>
*><br>
*>    ( x(i) ) := (  c(i)  s(i) ) ( x(i) )<br>
*>    ( y(i) )    ( -s(i)  c(i) ) ( y(i) )<br>
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
*>          X is REAL array,<br>
*>                         dimension (1+(N-1)*INCX)<br>
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
*>          Y is REAL array,<br>
*>                         dimension (1+(N-1)*INCY)<br>
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
*>          S is REAL array, dimension (1+(N-1)*INCC)<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slartv_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] C,float[] S,INTEGER INCC);
/**
*> \brief \b SLARUV returns a vector of n random real numbers from a uniform distribution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARUV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaruv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaruv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaruv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARUV( ISEED, N, X )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISEED( 4 )<br>
*       REAL               X( N )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARUV returns a vector of n random real numbers from a uniform (0,1)<br>
*> distribution (n <= 128).<br>
*><br>
*> This is an auxiliary routine called by SLARNV and CLARNV.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
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
*>          The number of random numbers to be generated. N <= 128.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (N)<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  This routine uses a multiplicative congruential method with modulus<br>
*>  2**48 and multiplier 33952834046453 (see G.S.Fishman,<br>
*>  'Multiplicative congruential random number generators with modulus<br>
*>  2**b: an exhaustive analysis for b = 32 and a partial analysis for<br>
*>  b = 48', Math. Comp. 189, pp 331-344, 1990).<br>
*><br>
*>  48-bit integers are stored in 4 integer array elements with 12 bits<br>
*>  per element. Hence the routine is portable across machines with<br>
*>  integers of 32 bits or more.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slaruv_(int[] ISEED,INTEGER N,float[] X);
/**
*> \brief \b SLARZ applies an elementary reflector (as returned by stzrzf) to a general matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE<br>
*       INTEGER            INCV, L, LDC, M, N<br>
*       REAL               TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( LDC, * ), V( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARZ applies a real elementary reflector H to a real M-by-N<br>
*> matrix C, from either the left or the right. H is represented in the<br>
*> form<br>
*><br>
*>       H = I - tau * v * v**T<br>
*><br>
*> where tau is a real scalar and v is a real vector.<br>
*><br>
*> If tau = 0, then H is taken to be the unit matrix.<br>
*><br>
*><br>
*> H is a product of k elementary reflectors as returned by STZRZF.<br>
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
*>          V is REAL array, dimension (1+(L-1)*abs(INCV))<br>
*>          The vector v in the representation of H as returned by<br>
*>          STZRZF. V is not used if TAU = 0.<br>
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
*>          TAU is REAL<br>
*>          The value tau in the representation of H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void slarz_(CHARACTER SIDE,INTEGER M,INTEGER N,INTEGER L,float[] V,INTEGER INCV,REAL TAU,float[] C,INTEGER LDC,float[] WORK);
/**
*> \brief \b SLARZB applies a block reflector or its transpose to a general matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARZB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarzb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarzb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarzb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V,<br>
*                          LDV, T, LDT, C, LDC, WORK, LDWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, SIDE, STOREV, TRANS<br>
*       INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ),<br>
*      $                   WORK( LDWORK, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARZB applies a real block reflector H or its transpose H**T to<br>
*> a real distributed M-by-N  C from the left or the right.<br>
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
*>          = 'L': apply H or H**T from the Left<br>
*>          = 'R': apply H or H**T from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply H (No transpose)<br>
*>          = 'C': apply H**T (Transpose)<br>
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
*>          V is REAL array, dimension (LDV,NV).<br>
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
*>          T is REAL array, dimension (LDT,K)<br>
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
*>          C is REAL array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.<br>
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
*>          WORK is REAL array, dimension (LDWORK,K)<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void slarzb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,INTEGER L,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] C,INTEGER LDC,float[] WORK,INTEGER LDWORK);
/**
*> \brief \b SLARZT forms the triangular factor T of a block reflector H = I - vtvH.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLARZT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarzt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarzt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarzt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, STOREV<br>
*       INTEGER            K, LDT, LDV, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               T( LDT, * ), TAU( * ), V( LDV, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLARZT forms the triangular factor T of a real block reflector<br>
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
*>    H  =  I - V * T * V**T<br>
*><br>
*> If STOREV = 'R', the vector which defines the elementary reflector<br>
*> H(i) is stored in the i-th row of the array V, and<br>
*><br>
*>    H  =  I - V**T * T * V<br>
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
*>          V is REAL array, dimension<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is REAL array, dimension (LDT,K)<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void slarzt_(CHARACTER DIRECT,CHARACTER STOREV,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] TAU,float[] T,INTEGER LDT);
/**
*> \brief \b SLAS2 computes singular values of a 2-by-2 triangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slas2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slas2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slas2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAS2( F, G, H, SSMIN, SSMAX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               F, G, H, SSMAX, SSMIN<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAS2  computes the singular values of the 2-by-2 matrix<br>
*>    [  F   G  ]<br>
*>    [  0   H  ].<br>
*> On return, SSMIN is the smaller singular value and SSMAX is the<br>
*> larger singular value.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] F<br>
*> \verbatim<br>
*>          F is REAL<br>
*>          The (1,1) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] G<br>
*> \verbatim<br>
*>          G is REAL<br>
*>          The (1,2) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] H<br>
*> \verbatim<br>
*>          H is REAL<br>
*>          The (2,2) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SSMIN<br>
*> \verbatim<br>
*>          SSMIN is REAL<br>
*>          The smaller singular value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SSMAX<br>
*> \verbatim<br>
*>          SSMAX is REAL<br>
*>          The larger singular value.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Barring over/underflow, all output quantities are correct to within<br>
*>  a few units in the last place (ulps), even in the absence of a guard<br>
*>  digit in addition/subtraction.<br>
*><br>
*>  In IEEE arithmetic, the code works correctly if one matrix element is<br>
*>  infinite.<br>
*><br>
*>  Overflow will not occur unless the largest singular value itself<br>
*>  overflows, or is within a few ulps of overflow. (On machines with<br>
*>  partial overflow, like the Cray, overflow may occur if the largest<br>
*>  singular value is within a factor of 2 of overflow.)<br>
*><br>
*>  Underflow is harmless if underflow is gradual. Otherwise, results<br>
*>  may correspond to a matrix modified by perturbations of size near<br>
*>  the underflow threshold.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slas2_(REAL F,REAL G,REAL H,REAL SSMIN,REAL SSMAX);
/**
*> \brief \b SLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASCL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slascl.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slascl.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slascl.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TYPE<br>
*       INTEGER            INFO, KL, KU, LDA, M, N<br>
*       REAL               CFROM, CTO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASCL multiplies the M by N real matrix A by the real scalar<br>
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
*>                  bandwidth KU. See SGBTRF for storage details.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slascl_(CHARACTER TYPE,INTEGER KL,INTEGER KU,REAL CFROM,REAL CTO,INTEGER M,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b SLASCL2 performs diagonal scaling on a vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASCL2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slascl2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slascl2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slascl2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASCL2 ( M, N, D, X, LDX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDX<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASCL2 performs a diagonal scaling on a vector:<br>
*>   x <-- D * x<br>
*> where the diagonal matrix D is stored as a vector.<br>
*><br>
*> Eventually to be replaced by BLAS_sge_diag_scale in the new BLAS<br>
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
*>          X is REAL array, dimension (LDX,N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slascl2_(INTEGER M,INTEGER N,float[] D,float[] X,INTEGER LDX);
/**
*> \brief \b SLASD0 computes the singular values of a real upper bidiagonal n-by-m matrix B with diagonal d and off-diagonal e. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASD0 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd0.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd0.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd0.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDU, LDVT, N, SMLSIZ, SQRE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               D( * ), E( * ), U( LDU, * ), VT( LDVT, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Using a divide and conquer approach, SLASD0 computes the singular<br>
*> value decomposition (SVD) of a real upper bidiagonal N-by-M<br>
*> matrix B with diagonal D and offdiagonal E, where M = N + SQRE.<br>
*> The algorithm computes orthogonal matrices U and VT such that<br>
*> B = U * S * VT. The singular values S are overwritten on D.<br>
*><br>
*> A related subroutine, SLASDA, computes only the singular values,<br>
*> and optionally, the singular vectors in compact form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         On entry, the row dimension of the upper bidiagonal matrix.<br>
*>         This is also the dimension of the main diagonal array D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SQRE<br>
*> \verbatim<br>
*>          SQRE is INTEGER<br>
*>         Specifies the column dimension of the bidiagonal matrix.<br>
*>         = 0: The bidiagonal matrix has column dimension M = N;<br>
*>         = 1: The bidiagonal matrix has column dimension M = N+1;<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry D contains the main diagonal of the bidiagonal<br>
*>         matrix.<br>
*>         On exit D, if INFO = 0, contains its singular values.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (M-1)<br>
*>         Contains the subdiagonal entries of the bidiagonal matrix.<br>
*>         On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is REAL array, dimension at least (LDQ, N)<br>
*>         On exit, U contains the left singular vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>         On entry, leading dimension of U.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VT<br>
*> \verbatim<br>
*>          VT is REAL array, dimension at least (LDVT, M)<br>
*>         On exit, VT**T contains the right singular vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>         On entry, leading dimension of VT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SMLSIZ<br>
*> \verbatim<br>
*>          SMLSIZ is INTEGER<br>
*>         On entry, maximum size of the subproblems at the<br>
*>         bottom of the computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (8*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (3*M**2+2*M)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = 1, a singular value did not converge<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd0_(INTEGER N,INTEGER SQRE,float[] D,float[] E,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,INTEGER SMLSIZ,int[] IWORK,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLASD1 computes the SVD of an upper bidiagonal matrix B of the specified size. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASD1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD1( NL, NR, SQRE, D, ALPHA, BETA, U, LDU, VT, LDVT,<br>
*                          IDXQ, IWORK, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDU, LDVT, NL, NR, SQRE<br>
*       REAL               ALPHA, BETA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IDXQ( * ), IWORK( * )<br>
*       REAL               D( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASD1 computes the SVD of an upper bidiagonal N-by-M matrix B,<br>
*> where N = NL + NR + 1 and M = N + SQRE. SLASD1 is called from SLASD0.<br>
*><br>
*> A related subroutine SLASD7 handles the case in which the singular<br>
*> values (and the singular vectors in factored form) are desired.<br>
*><br>
*> SLASD1 computes the SVD as follows:<br>
*><br>
*>               ( D1(in)    0    0       0 )<br>
*>   B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)<br>
*>               (   0       0   D2(in)   0 )<br>
*><br>
*>     = U(out) * ( D(out) 0) * VT(out)<br>
*><br>
*> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M<br>
*> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros<br>
*> elsewhere; and the entry b is empty if SQRE = 0.<br>
*><br>
*> The left singular vectors of the original matrix are stored in U, and<br>
*> the transpose of the right singular vectors are stored in VT, and the<br>
*> singular values are in D.  The algorithm consists of three stages:<br>
*><br>
*>    The first stage consists of deflating the size of the problem<br>
*>    when there are multiple singular values or when there are zeros in<br>
*>    the Z vector.  For each such occurrence the dimension of the<br>
*>    secular equation problem is reduced by one.  This stage is<br>
*>    performed by the routine SLASD2.<br>
*><br>
*>    The second stage consists of calculating the updated<br>
*>    singular values. This is done by finding the square roots of the<br>
*>    roots of the secular equation via the routine SLASD4 (as called<br>
*>    by SLASD3). This routine also calculates the singular vectors of<br>
*>    the current problem.<br>
*><br>
*>    The final stage consists of computing the updated singular vectors<br>
*>    directly using the updated singular values.  The singular vectors<br>
*>    for the current problem are multiplied with the singular vectors<br>
*>    from the overall problem.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NL<br>
*> \verbatim<br>
*>          NL is INTEGER<br>
*>         The row dimension of the upper block.  NL >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NR<br>
*> \verbatim<br>
*>          NR is INTEGER<br>
*>         The row dimension of the lower block.  NR >= 1.<br>
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
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (NL+NR+1).<br>
*>         N = NL+NR+1<br>
*>         On entry D(1:NL,1:NL) contains the singular values of the<br>
*>         upper block; and D(NL+2:N) contains the singular values of<br>
*>         the lower block. On exit D(1:N) contains the singular values<br>
*>         of the modified matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>         Contains the diagonal element associated with the added row.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>         Contains the off-diagonal element associated with the added<br>
*>         row.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] U<br>
*> \verbatim<br>
*>          U is REAL array, dimension (LDU,N)<br>
*>         On entry U(1:NL, 1:NL) contains the left singular vectors of<br>
*>         the upper block; U(NL+2:N, NL+2:N) contains the left singular<br>
*>         vectors of the lower block. On exit U contains the left<br>
*>         singular vectors of the bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>         The leading dimension of the array U.  LDU >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VT<br>
*> \verbatim<br>
*>          VT is REAL array, dimension (LDVT,M)<br>
*>         where M = N + SQRE.<br>
*>         On entry VT(1:NL+1, 1:NL+1)**T contains the right singular<br>
*>         vectors of the upper block; VT(NL+2:M, NL+2:M)**T contains<br>
*>         the right singular vectors of the lower block. On exit<br>
*>         VT**T contains the right singular vectors of the<br>
*>         bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>         The leading dimension of the array VT.  LDVT >= max( 1, M ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IDXQ<br>
*> \verbatim<br>
*>          IDXQ is INTEGER array, dimension (N)<br>
*>         This contains the permutation which will reintegrate the<br>
*>         subproblem just solved back into sorted order, i.e.<br>
*>         D( IDXQ( I = 1, N ) ) will be in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (3*M**2+2*M)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = 1, a singular value did not converge<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd1_(INTEGER NL,INTEGER NR,INTEGER SQRE,float[] D,REAL ALPHA,REAL BETA,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,int[] IDXQ,int[] IWORK,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLASD2 merges the two sets of singular values together into a single sorted set. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASD2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD2( NL, NR, SQRE, K, D, Z, ALPHA, BETA, U, LDU, VT,<br>
*                          LDVT, DSIGMA, U2, LDU2, VT2, LDVT2, IDXP, IDX,<br>
*                          IDXC, IDXQ, COLTYP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDU, LDU2, LDVT, LDVT2, NL, NR, SQRE<br>
*       REAL               ALPHA, BETA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            COLTYP( * ), IDX( * ), IDXC( * ), IDXP( * ),<br>
*      $                   IDXQ( * )<br>
*       REAL               D( * ), DSIGMA( * ), U( LDU, * ),<br>
*      $                   U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ),<br>
*      $                   Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASD2 merges the two sets of singular values together into a single<br>
*> sorted set.  Then it tries to deflate the size of the problem.<br>
*> There are two ways in which deflation can occur:  when two or more<br>
*> singular values are close together or if there is a tiny entry in the<br>
*> Z vector.  For each such occurrence the order of the related secular<br>
*> equation problem is reduced by one.<br>
*><br>
*> SLASD2 is called from SLASD1.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NL<br>
*> \verbatim<br>
*>          NL is INTEGER<br>
*>         The row dimension of the upper block.  NL >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NR<br>
*> \verbatim<br>
*>          NR is INTEGER<br>
*>         The row dimension of the lower block.  NR >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SQRE<br>
*> \verbatim<br>
*>          SQRE is INTEGER<br>
*>         = 0: the lower block is an NR-by-NR square matrix.<br>
*>         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.<br>
*><br>
*>         The bidiagonal matrix has N = NL + NR + 1 rows and<br>
*>         M = N + SQRE >= N columns.<br>
*> \endverbatim<br>
*><br>
*> \param[out] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>         Contains the dimension of the non-deflated matrix,<br>
*>         This is the order of the related secular equation. 1 <= K <=N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>         On entry D contains the singular values of the two submatrices<br>
*>         to be combined.  On exit D contains the trailing (N-K) updated<br>
*>         singular values (those which were deflated) sorted into<br>
*>         increasing order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (N)<br>
*>         On exit Z contains the updating row vector in the secular<br>
*>         equation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>         Contains the diagonal element associated with the added row.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>         Contains the off-diagonal element associated with the added<br>
*>         row.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] U<br>
*> \verbatim<br>
*>          U is REAL array, dimension (LDU,N)<br>
*>         On entry U contains the left singular vectors of two<br>
*>         submatrices in the two square blocks with corners at (1,1),<br>
*>         (NL, NL), and (NL+2, NL+2), (N,N).<br>
*>         On exit U contains the trailing (N-K) updated left singular<br>
*>         vectors (those which were deflated) in its last N-K columns.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>         The leading dimension of the array U.  LDU >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VT<br>
*> \verbatim<br>
*>          VT is REAL array, dimension (LDVT,M)<br>
*>         On entry VT**T contains the right singular vectors of two<br>
*>         submatrices in the two square blocks with corners at (1,1),<br>
*>         (NL+1, NL+1), and (NL+2, NL+2), (M,M).<br>
*>         On exit VT**T contains the trailing (N-K) updated right singular<br>
*>         vectors (those which were deflated) in its last N-K columns.<br>
*>         In case SQRE =1, the last row of VT spans the right null<br>
*>         space.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>         The leading dimension of the array VT.  LDVT >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DSIGMA<br>
*> \verbatim<br>
*>          DSIGMA is REAL array, dimension (N)<br>
*>         Contains a copy of the diagonal elements (K-1 singular values<br>
*>         and one zero) in the secular equation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] U2<br>
*> \verbatim<br>
*>          U2 is REAL array, dimension (LDU2,N)<br>
*>         Contains a copy of the first K-1 left singular vectors which<br>
*>         will be used by SLASD3 in a matrix multiply (SGEMM) to solve<br>
*>         for the new left singular vectors. U2 is arranged into four<br>
*>         blocks. The first block contains a column with 1 at NL+1 and<br>
*>         zero everywhere else; the second block contains non-zero<br>
*>         entries only at and above NL; the third contains non-zero<br>
*>         entries only below NL+1; and the fourth is dense.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU2<br>
*> \verbatim<br>
*>          LDU2 is INTEGER<br>
*>         The leading dimension of the array U2.  LDU2 >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VT2<br>
*> \verbatim<br>
*>          VT2 is REAL array, dimension (LDVT2,N)<br>
*>         VT2**T contains a copy of the first K right singular vectors<br>
*>         which will be used by SLASD3 in a matrix multiply (SGEMM) to<br>
*>         solve for the new right singular vectors. VT2 is arranged into<br>
*>         three blocks. The first block contains a row that corresponds<br>
*>         to the special 0 diagonal element in SIGMA; the second block<br>
*>         contains non-zeros only at and before NL +1; the third block<br>
*>         contains non-zeros only at and after  NL +2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT2<br>
*> \verbatim<br>
*>          LDVT2 is INTEGER<br>
*>         The leading dimension of the array VT2.  LDVT2 >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IDXP<br>
*> \verbatim<br>
*>          IDXP is INTEGER array, dimension (N)<br>
*>         This will contain the permutation used to place deflated<br>
*>         values of D at the end of the array. On output IDXP(2:K)<br>
*>         points to the nondeflated D-values and IDXP(K+1:N)<br>
*>         points to the deflated singular values.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IDX<br>
*> \verbatim<br>
*>          IDX is INTEGER array, dimension (N)<br>
*>         This will contain the permutation used to sort the contents of<br>
*>         D into ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IDXC<br>
*> \verbatim<br>
*>          IDXC is INTEGER array, dimension (N)<br>
*>         This will contain the permutation used to arrange the columns<br>
*>         of the deflated U matrix into three groups:  the first group<br>
*>         contains non-zero entries only at and above NL, the second<br>
*>         contains non-zero entries only below NL+2, and the third is<br>
*>         dense.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IDXQ<br>
*> \verbatim<br>
*>          IDXQ is INTEGER array, dimension (N)<br>
*>         This contains the permutation which separately sorts the two<br>
*>         sub-problems in D into ascending order.  Note that entries in<br>
*>         the first hlaf of this permutation must first be moved one<br>
*>         position backward; and entries in the second half<br>
*>         must first have NL+1 added to their values.<br>
*> \endverbatim<br>
*><br>
*> \param[out] COLTYP<br>
*> \verbatim<br>
*>          COLTYP is INTEGER array, dimension (N)<br>
*>         As workspace, this will contain a label which will indicate<br>
*>         which of the following types a column in the U2 matrix or a<br>
*>         row in the VT2 matrix is:<br>
*>         1 : non-zero in the upper half only<br>
*>         2 : non-zero in the lower half only<br>
*>         3 : dense<br>
*>         4 : deflated<br>
*><br>
*>         On exit, it is an array of dimension 4, with COLTYP(I) being<br>
*>         the dimension of the I-th type columns.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd2_(INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,float[] D,float[] Z,REAL ALPHA,REAL BETA,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] DSIGMA,float[] U2,INTEGER LDU2,float[] VT2,INTEGER LDVT2,int[] IDXP,int[] IDX,int[] IDXC,int[] IDXQ,int[] COLTYP,INTEGER INFO);
/**
*> \brief \b SLASD3 finds all square roots of the roots of the secular equation, as defined by the values in D and Z, and then updates the singular vectors by matrix multiplication. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASD3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2,<br>
*                          LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDQ, LDU, LDU2, LDVT, LDVT2, NL, NR,<br>
*      $                   SQRE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            CTOT( * ), IDXC( * )<br>
*       REAL               D( * ), DSIGMA( * ), Q( LDQ, * ), U( LDU, * ),<br>
*      $                   U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ),<br>
*      $                   Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASD3 finds all the square roots of the roots of the secular<br>
*> equation, as defined by the values in D and Z.  It makes the<br>
*> appropriate calls to SLASD4 and then updates the singular<br>
*> vectors by matrix multiplication.<br>
*><br>
*> This code makes very mild assumptions about floating point<br>
*> arithmetic. It will work on machines with a guard digit in<br>
*> add/subtract, or on those binary machines without guard digits<br>
*> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.<br>
*> It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.<br>
*><br>
*> SLASD3 is called from SLASD1.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NL<br>
*> \verbatim<br>
*>          NL is INTEGER<br>
*>         The row dimension of the upper block.  NL >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NR<br>
*> \verbatim<br>
*>          NR is INTEGER<br>
*>         The row dimension of the lower block.  NR >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SQRE<br>
*> \verbatim<br>
*>          SQRE is INTEGER<br>
*>         = 0: the lower block is an NR-by-NR square matrix.<br>
*>         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.<br>
*><br>
*>         The bidiagonal matrix has N = NL + NR + 1 rows and<br>
*>         M = N + SQRE >= N columns.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>         The size of the secular equation, 1 =< K = < N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension(K)<br>
*>         On exit the square roots of the roots of the secular equation,<br>
*>         in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is REAL array,<br>
*>                     dimension at least (LDQ,K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>         The leading dimension of the array Q.  LDQ >= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DSIGMA<br>
*> \verbatim<br>
*>          DSIGMA is REAL array, dimension(K)<br>
*>         The first K elements of this array contain the old roots<br>
*>         of the deflated updating problem.  These are the poles<br>
*>         of the secular equation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is REAL array, dimension (LDU, N)<br>
*>         The last N - K columns of this matrix contain the deflated<br>
*>         left singular vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>         The leading dimension of the array U.  LDU >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] U2<br>
*> \verbatim<br>
*>          U2 is REAL array, dimension (LDU2, N)<br>
*>         The first K columns of this matrix contain the non-deflated<br>
*>         left singular vectors for the split problem.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU2<br>
*> \verbatim<br>
*>          LDU2 is INTEGER<br>
*>         The leading dimension of the array U2.  LDU2 >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VT<br>
*> \verbatim<br>
*>          VT is REAL array, dimension (LDVT, M)<br>
*>         The last M - K columns of VT**T contain the deflated<br>
*>         right singular vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>         The leading dimension of the array VT.  LDVT >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VT2<br>
*> \verbatim<br>
*>          VT2 is REAL array, dimension (LDVT2, N)<br>
*>         The first K columns of VT2**T contain the non-deflated<br>
*>         right singular vectors for the split problem.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT2<br>
*> \verbatim<br>
*>          LDVT2 is INTEGER<br>
*>         The leading dimension of the array VT2.  LDVT2 >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IDXC<br>
*> \verbatim<br>
*>          IDXC is INTEGER array, dimension (N)<br>
*>         The permutation used to arrange the columns of U (and rows of<br>
*>         VT) into three groups:  the first group contains non-zero<br>
*>         entries only at and above (or before) NL +1; the second<br>
*>         contains non-zero entries only at and below (or after) NL+2;<br>
*>         and the third is dense. The first column of U and the row of<br>
*>         VT are treated separately, however.<br>
*><br>
*>         The rows of the singular vectors found by SLASD4<br>
*>         must be likewise permuted before the matrix multiplies can<br>
*>         take place.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CTOT<br>
*> \verbatim<br>
*>          CTOT is INTEGER array, dimension (4)<br>
*>         A count of the total number of the various types of columns<br>
*>         in U (or rows in VT), as described in IDXC. The fourth column<br>
*>         type is any column which has been deflated.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (K)<br>
*>         The first K elements of this array contain the components<br>
*>         of the deflation-adjusted updating row vector.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>         = 0:  successful exit.<br>
*>         < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>         > 0:  if INFO = 1, a singular value did not converge<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd3_(INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,float[] D,float[] Q,INTEGER LDQ,float[] DSIGMA,float[] U,INTEGER LDU,float[] U2,INTEGER LDU2,float[] VT,INTEGER LDVT,float[] VT2,INTEGER LDVT2,int[] IDXC,int[] CTOT,float[] Z,INTEGER INFO);
/**
*> \brief \b SLASD4 computes the square root of the i-th updated eigenvalue of a positive symmetric rank-one modification to a positive diagonal matrix. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download SLASD4 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd4.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd4.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd4.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       INTEGER            I, INFO, N<br>
*       REAL               RHO, SIGMA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), DELTA( * ), WORK( * ), Z( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This subroutine computes the square root of the I-th updated<br>
*> eigenvalue of a positive symmetric rank-one modification to<br>
*> a positive diagonal matrix whose entries are given as the squares<br>
*> of the corresponding entries in the array d, and that<br>
*><br>
*>        0 <= D(i) < D(j)  for  i < j<br>
*><br>
*> and that RHO > 0. This is arranged by the calling routine, and is<br>
*> no loss in generality.  The rank-one modified system is thus<br>
*><br>
*>        diag( D ) * diag( D ) +  RHO * Z * Z_transpose.<br>
*><br>
*> where we assume the Euclidean norm of Z is 1.<br>
*><br>
*> The method consists of approximating the rational functions in the<br>
*> secular equation by simpler interpolating rational functions.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The length of all arrays.<br>
*> \endverbatim<br>
*><br>
*> \param[in] I<br>
*> \verbatim<br>
*>          I is INTEGER<br>
*>         The index of the eigenvalue to be computed.  1 <= I <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension ( N )<br>
*>         The original eigenvalues.  It is assumed that they are in<br>
*>         order, 0 <= D(I) < D(J)  for I < J.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( N )<br>
*>         The components of the updating vector.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DELTA<br>
*> \verbatim<br>
*>          DELTA is REAL array, dimension ( N )<br>
*>         If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th<br>
*>         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA<br>
*>         contains the information necessary to construct the<br>
*>         (singular) eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         The scalar in the symmetric updating formula.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SIGMA<br>
*> \verbatim<br>
*>          SIGMA is REAL<br>
*>         The computed sigma_I, the I-th updated eigenvalue.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension ( N )<br>
*>         If N .ne. 1, WORK contains (D(j) + sigma_I) in its  j-th<br>
*>         component.  If N = 1, then WORK( 1 ) = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>         = 0:  successful exit<br>
*>         > 0:  if INFO = 1, the updating process failed.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  Logical variable ORGATI (origin-at-i?) is used for distinguishing<br>
*>  whether D(i) or D(i+1) is treated as the origin.<br>
*><br>
*>            ORGATI = .true.    origin at i<br>
*>            ORGATI = .false.   origin at i+1<br>
*><br>
*>  Logical variable SWTCH3 (switch-for-3-poles?) is for noting<br>
*>  if we are working with THREE poles!<br>
*><br>
*>  MAXIT is the maximum number of iterations allowed for each<br>
*>  eigenvalue.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ren-Cang Li, Computer Science Division, University of California<br>
*>     at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd4_(INTEGER N,INTEGER I,float[] D,float[] Z,float[] DELTA,REAL RHO,REAL SIGMA,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLASD5 computes the square root of the i-th eigenvalue of a positive symmetric rank-one modification of a 2-by-2 diagonal matrix. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASD5 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd5.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd5.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd5.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD5( I, D, Z, DELTA, RHO, DSIGMA, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            I<br>
*       REAL               DSIGMA, RHO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( 2 ), DELTA( 2 ), WORK( 2 ), Z( 2 )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This subroutine computes the square root of the I-th eigenvalue<br>
*> of a positive symmetric rank-one modification of a 2-by-2 diagonal<br>
*> matrix<br>
*><br>
*>            diag( D ) * diag( D ) +  RHO * Z * transpose(Z) .<br>
*><br>
*> The diagonal entries in the array D are assumed to satisfy<br>
*><br>
*>            0 <= D(i) < D(j)  for  i < j .<br>
*><br>
*> We also assume RHO > 0 and that the Euclidean norm of the vector<br>
*> Z is one.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] I<br>
*> \verbatim<br>
*>          I is INTEGER<br>
*>         The index of the eigenvalue to be computed.  I = 1 or I = 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (2)<br>
*>         The original eigenvalues.  We assume 0 <= D(1) < D(2).<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (2)<br>
*>         The components of the updating vector.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DELTA<br>
*> \verbatim<br>
*>          DELTA is REAL array, dimension (2)<br>
*>         Contains (D(j) - sigma_I) in its  j-th component.<br>
*>         The vector DELTA contains the information necessary<br>
*>         to construct the eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RHO<br>
*> \verbatim<br>
*>          RHO is REAL<br>
*>         The scalar in the symmetric updating formula.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DSIGMA<br>
*> \verbatim<br>
*>          DSIGMA is REAL<br>
*>         The computed sigma_I, the I-th updated eigenvalue.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2)<br>
*>         WORK contains (D(j) + sigma_I) in its  j-th component.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ren-Cang Li, Computer Science Division, University of California<br>
*>     at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd5_(INTEGER I,float[] D,float[] Z,float[] DELTA,REAL RHO,REAL DSIGMA,float[] WORK);
/**
*> \brief \b SLASD6 computes the SVD of an updated upper bidiagonal matrix obtained by merging two smaller ones by appending a row. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASD6 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd6.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd6.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd6.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD6( ICOMPQ, NL, NR, SQRE, D, VF, VL, ALPHA, BETA,<br>
*                          IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM,<br>
*                          LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL,<br>
*      $                   NR, SQRE<br>
*       REAL               ALPHA, BETA, C, S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( LDGCOL, * ), IDXQ( * ), IWORK( * ),<br>
*      $                   PERM( * )<br>
*       REAL               D( * ), DIFL( * ), DIFR( * ),<br>
*      $                   GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ),<br>
*      $                   VF( * ), VL( * ), WORK( * ), Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASD6 computes the SVD of an updated upper bidiagonal matrix B<br>
*> obtained by merging two smaller ones by appending a row. This<br>
*> routine is used only for the problem which requires all singular<br>
*> values and optionally singular vector matrices in factored form.<br>
*> B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.<br>
*> A related subroutine, SLASD1, handles the case in which all singular<br>
*> values and singular vectors of the bidiagonal matrix are desired.<br>
*><br>
*> SLASD6 computes the SVD as follows:<br>
*><br>
*>               ( D1(in)    0    0       0 )<br>
*>   B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)<br>
*>               (   0       0   D2(in)   0 )<br>
*><br>
*>     = U(out) * ( D(out) 0) * VT(out)<br>
*><br>
*> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M<br>
*> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros<br>
*> elsewhere; and the entry b is empty if SQRE = 0.<br>
*><br>
*> The singular values of B can be computed using D1, D2, the first<br>
*> components of all the right singular vectors of the lower block, and<br>
*> the last components of all the right singular vectors of the upper<br>
*> block. These components are stored and updated in VF and VL,<br>
*> respectively, in SLASD6. Hence U and VT are not explicitly<br>
*> referenced.<br>
*><br>
*> The singular values are stored in D. The algorithm consists of two<br>
*> stages:<br>
*><br>
*>       The first stage consists of deflating the size of the problem<br>
*>       when there are multiple singular values or if there is a zero<br>
*>       in the Z vector. For each such occurrence the dimension of the<br>
*>       secular equation problem is reduced by one. This stage is<br>
*>       performed by the routine SLASD7.<br>
*><br>
*>       The second stage consists of calculating the updated<br>
*>       singular values. This is done by finding the roots of the<br>
*>       secular equation via the routine SLASD4 (as called by SLASD8).<br>
*>       This routine also updates VF and VL and computes the distances<br>
*>       between the updated singular values and the old singular<br>
*>       values.<br>
*><br>
*> SLASD6 is called from SLASDA.<br>
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
*>         = 0: Compute singular values only.<br>
*>         = 1: Compute singular vectors in factored form as well.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NL<br>
*> \verbatim<br>
*>          NL is INTEGER<br>
*>         The row dimension of the upper block.  NL >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NR<br>
*> \verbatim<br>
*>          NR is INTEGER<br>
*>         The row dimension of the lower block.  NR >= 1.<br>
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
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (NL+NR+1).<br>
*>         On entry D(1:NL,1:NL) contains the singular values of the<br>
*>         upper block, and D(NL+2:N) contains the singular values<br>
*>         of the lower block. On exit D(1:N) contains the singular<br>
*>         values of the modified matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VF<br>
*> \verbatim<br>
*>          VF is REAL array, dimension (M)<br>
*>         On entry, VF(1:NL+1) contains the first components of all<br>
*>         right singular vectors of the upper block; and VF(NL+2:M)<br>
*>         contains the first components of all right singular vectors<br>
*>         of the lower block. On exit, VF contains the first components<br>
*>         of all right singular vectors of the bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is REAL array, dimension (M)<br>
*>         On entry, VL(1:NL+1) contains the  last components of all<br>
*>         right singular vectors of the upper block; and VL(NL+2:M)<br>
*>         contains the last components of all right singular vectors of<br>
*>         the lower block. On exit, VL contains the last components of<br>
*>         all right singular vectors of the bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>         Contains the diagonal element associated with the added row.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>         Contains the off-diagonal element associated with the added<br>
*>         row.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IDXQ<br>
*> \verbatim<br>
*>          IDXQ is INTEGER array, dimension (N)<br>
*>         This contains the permutation which will reintegrate the<br>
*>         subproblem just solved back into sorted order, i.e.<br>
*>         D( IDXQ( I = 1, N ) ) will be in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PERM<br>
*> \verbatim<br>
*>          PERM is INTEGER array, dimension ( N )<br>
*>         The permutations (from deflation and sorting) to be applied<br>
*>         to each block. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVPTR<br>
*> \verbatim<br>
*>          GIVPTR is INTEGER<br>
*>         The number of Givens rotations which took place in this<br>
*>         subproblem. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVCOL<br>
*> \verbatim<br>
*>          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 )<br>
*>         Each pair of numbers indicates a pair of columns to take place<br>
*>         in a Givens rotation. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDGCOL<br>
*> \verbatim<br>
*>          LDGCOL is INTEGER<br>
*>         leading dimension of GIVCOL, must be at least N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVNUM<br>
*> \verbatim<br>
*>          GIVNUM is REAL array, dimension ( LDGNUM, 2 )<br>
*>         Each number indicates the C or S value to be used in the<br>
*>         corresponding Givens rotation. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDGNUM<br>
*> \verbatim<br>
*>          LDGNUM is INTEGER<br>
*>         The leading dimension of GIVNUM and POLES, must be at least N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] POLES<br>
*> \verbatim<br>
*>          POLES is REAL array, dimension ( LDGNUM, 2 )<br>
*>         On exit, POLES(1,*) is an array containing the new singular<br>
*>         values obtained from solving the secular equation, and<br>
*>         POLES(2,*) is an array containing the poles in the secular<br>
*>         equation. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIFL<br>
*> \verbatim<br>
*>          DIFL is REAL array, dimension ( N )<br>
*>         On exit, DIFL(I) is the distance between I-th updated<br>
*>         (undeflated) singular value and the I-th (undeflated) old<br>
*>         singular value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIFR<br>
*> \verbatim<br>
*>          DIFR is REAL array,<br>
*>                   dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and<br>
*>                   dimension ( K ) if ICOMPQ = 0.<br>
*>          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not<br>
*>          defined and will not be referenced.<br>
*><br>
*>          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the<br>
*>          normalizing factors for the right singular vector matrix.<br>
*><br>
*>         See SLASD8 for details on DIFL and DIFR.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( M )<br>
*>         The first elements of this array contain the components<br>
*>         of the deflation-adjusted updating row vector.<br>
*> \endverbatim<br>
*><br>
*> \param[out] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>         Contains the dimension of the non-deflated matrix,<br>
*>         This is the order of the related secular equation. 1 <= K <=N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is REAL<br>
*>         C contains garbage if SQRE =0 and the C-value of a Givens<br>
*>         rotation related to the right null space if SQRE = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL<br>
*>         S contains garbage if SQRE =0 and the S-value of a Givens<br>
*>         rotation related to the right null space if SQRE = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension ( 4 * M )<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension ( 3 * N )<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = 1, a singular value did not converge<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd6_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,float[] D,float[] VF,float[] VL,REAL ALPHA,REAL BETA,int[] IDXQ,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,float[] GIVNUM,INTEGER LDGNUM,float[] POLES,float[] DIFL,float[] DIFR,float[] Z,INTEGER K,REAL C,REAL S,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLASD7 merges the two sets of singular values together into a single sorted set. Then it tries to deflate the size of the problem. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASD7 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd7.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd7.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd7.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL,<br>
*                          VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ,<br>
*                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,<br>
*                          C, S, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL,<br>
*      $                   NR, SQRE<br>
*       REAL               ALPHA, BETA, C, S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( LDGCOL, * ), IDX( * ), IDXP( * ),<br>
*      $                   IDXQ( * ), PERM( * )<br>
*       REAL               D( * ), DSIGMA( * ), GIVNUM( LDGNUM, * ),<br>
*      $                   VF( * ), VFW( * ), VL( * ), VLW( * ), Z( * ),<br>
*      $                   ZW( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASD7 merges the two sets of singular values together into a single<br>
*> sorted set. Then it tries to deflate the size of the problem. There<br>
*> are two ways in which deflation can occur:  when two or more singular<br>
*> values are close together or if there is a tiny entry in the Z<br>
*> vector. For each such occurrence the order of the related<br>
*> secular equation problem is reduced by one.<br>
*><br>
*> SLASD7 is called from SLASD6.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ICOMPQ<br>
*> \verbatim<br>
*>          ICOMPQ is INTEGER<br>
*>          Specifies whether singular vectors are to be computed<br>
*>          in compact form, as follows:<br>
*>          = 0: Compute singular values only.<br>
*>          = 1: Compute singular vectors of upper<br>
*>               bidiagonal matrix in compact form.<br>
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
*>         The bidiagonal matrix has<br>
*>         N = NL + NR + 1 rows and<br>
*>         M = N + SQRE >= N columns.<br>
*> \endverbatim<br>
*><br>
*> \param[out] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>         Contains the dimension of the non-deflated matrix, this is<br>
*>         the order of the related secular equation. 1 <= K <=N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension ( N )<br>
*>         On entry D contains the singular values of the two submatrices<br>
*>         to be combined. On exit D contains the trailing (N-K) updated<br>
*>         singular values (those which were deflated) sorted into<br>
*>         increasing order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( M )<br>
*>         On exit Z contains the updating row vector in the secular<br>
*>         equation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ZW<br>
*> \verbatim<br>
*>          ZW is REAL array, dimension ( M )<br>
*>         Workspace for Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VF<br>
*> \verbatim<br>
*>          VF is REAL array, dimension ( M )<br>
*>         On entry, VF(1:NL+1) contains the first components of all<br>
*>         right singular vectors of the upper block; and VF(NL+2:M)<br>
*>         contains the first components of all right singular vectors<br>
*>         of the lower block. On exit, VF contains the first components<br>
*>         of all right singular vectors of the bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VFW<br>
*> \verbatim<br>
*>          VFW is REAL array, dimension ( M )<br>
*>         Workspace for VF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is REAL array, dimension ( M )<br>
*>         On entry, VL(1:NL+1) contains the  last components of all<br>
*>         right singular vectors of the upper block; and VL(NL+2:M)<br>
*>         contains the last components of all right singular vectors<br>
*>         of the lower block. On exit, VL contains the last components<br>
*>         of all right singular vectors of the bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VLW<br>
*> \verbatim<br>
*>          VLW is REAL array, dimension ( M )<br>
*>         Workspace for VL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>         Contains the diagonal element associated with the added row.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>         Contains the off-diagonal element associated with the added<br>
*>         row.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DSIGMA<br>
*> \verbatim<br>
*>          DSIGMA is REAL array, dimension ( N )<br>
*>         Contains a copy of the diagonal elements (K-1 singular values<br>
*>         and one zero) in the secular equation.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IDX<br>
*> \verbatim<br>
*>          IDX is INTEGER array, dimension ( N )<br>
*>         This will contain the permutation used to sort the contents of<br>
*>         D into ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IDXP<br>
*> \verbatim<br>
*>          IDXP is INTEGER array, dimension ( N )<br>
*>         This will contain the permutation used to place deflated<br>
*>         values of D at the end of the array. On output IDXP(2:K)<br>
*>         points to the nondeflated D-values and IDXP(K+1:N)<br>
*>         points to the deflated singular values.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IDXQ<br>
*> \verbatim<br>
*>          IDXQ is INTEGER array, dimension ( N )<br>
*>         This contains the permutation which separately sorts the two<br>
*>         sub-problems in D into ascending order.  Note that entries in<br>
*>         the first half of this permutation must first be moved one<br>
*>         position backward; and entries in the second half<br>
*>         must first have NL+1 added to their values.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PERM<br>
*> \verbatim<br>
*>          PERM is INTEGER array, dimension ( N )<br>
*>         The permutations (from deflation and sorting) to be applied<br>
*>         to each singular block. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVPTR<br>
*> \verbatim<br>
*>          GIVPTR is INTEGER<br>
*>         The number of Givens rotations which took place in this<br>
*>         subproblem. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVCOL<br>
*> \verbatim<br>
*>          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 )<br>
*>         Each pair of numbers indicates a pair of columns to take place<br>
*>         in a Givens rotation. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDGCOL<br>
*> \verbatim<br>
*>          LDGCOL is INTEGER<br>
*>         The leading dimension of GIVCOL, must be at least N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVNUM<br>
*> \verbatim<br>
*>          GIVNUM is REAL array, dimension ( LDGNUM, 2 )<br>
*>         Each number indicates the C or S value to be used in the<br>
*>         corresponding Givens rotation. Not referenced if ICOMPQ = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDGNUM<br>
*> \verbatim<br>
*>          LDGNUM is INTEGER<br>
*>         The leading dimension of GIVNUM, must be at least N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is REAL<br>
*>         C contains garbage if SQRE =0 and the C-value of a Givens<br>
*>         rotation related to the right null space if SQRE = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL<br>
*>         S contains garbage if SQRE =0 and the S-value of a Givens<br>
*>         rotation related to the right null space if SQRE = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>         = 0:  successful exit.<br>
*>         < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd7_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,float[] D,float[] Z,float[] ZW,float[] VF,float[] VFW,float[] VL,float[] VLW,REAL ALPHA,REAL BETA,float[] DSIGMA,int[] IDX,int[] IDXP,int[] IDXQ,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,float[] GIVNUM,INTEGER LDGNUM,REAL C,REAL S,INTEGER INFO);
/**
*> \brief \b SLASD8 finds the square roots of the roots of the secular equation, and stores, for each element in D, the distance to its two nearest poles. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASD8 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd8.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd8.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd8.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR,<br>
*                          DSIGMA, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            ICOMPQ, INFO, K, LDDIFR<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), DIFL( * ), DIFR( LDDIFR, * ),<br>
*      $                   DSIGMA( * ), VF( * ), VL( * ), WORK( * ),<br>
*      $                   Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASD8 finds the square roots of the roots of the secular equation,<br>
*> as defined by the values in DSIGMA and Z. It makes the appropriate<br>
*> calls to SLASD4, and stores, for each  element in D, the distance<br>
*> to its two nearest poles (elements in DSIGMA). It also updates<br>
*> the arrays VF and VL, the first and last components of all the<br>
*> right singular vectors of the original bidiagonal matrix.<br>
*><br>
*> SLASD8 is called from SLASD6.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ICOMPQ<br>
*> \verbatim<br>
*>          ICOMPQ is INTEGER<br>
*>          Specifies whether singular vectors are to be computed in<br>
*>          factored form in the calling routine:<br>
*>          = 0: Compute singular values only.<br>
*>          = 1: Compute singular vectors in factored form as well.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of terms in the rational function to be solved<br>
*>          by SLASD4.  K >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension ( K )<br>
*>          On output, D contains the updated singular values.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( K )<br>
*>          On entry, the first K elements of this array contain the<br>
*>          components of the deflation-adjusted updating row vector.<br>
*>          On exit, Z is updated.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VF<br>
*> \verbatim<br>
*>          VF is REAL array, dimension ( K )<br>
*>          On entry, VF contains  information passed through DBEDE8.<br>
*>          On exit, VF contains the first K components of the first<br>
*>          components of all right singular vectors of the bidiagonal<br>
*>          matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is REAL array, dimension ( K )<br>
*>          On entry, VL contains  information passed through DBEDE8.<br>
*>          On exit, VL contains the first K components of the last<br>
*>          components of all right singular vectors of the bidiagonal<br>
*>          matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIFL<br>
*> \verbatim<br>
*>          DIFL is REAL array, dimension ( K )<br>
*>          On exit, DIFL(I) = D(I) - DSIGMA(I).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIFR<br>
*> \verbatim<br>
*>          DIFR is REAL array,<br>
*>                   dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and<br>
*>                   dimension ( K ) if ICOMPQ = 0.<br>
*>          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not<br>
*>          defined and will not be referenced.<br>
*><br>
*>          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the<br>
*>          normalizing factors for the right singular vector matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDDIFR<br>
*> \verbatim<br>
*>          LDDIFR is INTEGER<br>
*>          The leading dimension of DIFR, must be at least K.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DSIGMA<br>
*> \verbatim<br>
*>          DSIGMA is REAL array, dimension ( K )<br>
*>          On entry, the first K elements of this array contain the old<br>
*>          roots of the deflated updating problem.  These are the poles<br>
*>          of the secular equation.<br>
*>          On exit, the elements of DSIGMA may be very slightly altered<br>
*>          in value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension at least 3 * K<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = 1, a singular value did not converge<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasd8_(INTEGER ICOMPQ,INTEGER K,float[] D,float[] Z,float[] VF,float[] VL,float[] DIFL,float[] DIFR,INTEGER LDDIFR,float[] DSIGMA,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLASDA computes the singular value decomposition (SVD) of a real upper bidiagonal matrix with diagonal d and off-diagonal e. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASDA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasda.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasda.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasda.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K,<br>
*                          DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL,<br>
*                          PERM, GIVNUM, C, S, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),<br>
*      $                   K( * ), PERM( LDGCOL, * )<br>
*       REAL               C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ),<br>
*      $                   E( * ), GIVNUM( LDU, * ), POLES( LDU, * ),<br>
*      $                   S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ),<br>
*      $                   Z( LDU, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Using a divide and conquer approach, SLASDA computes the singular<br>
*> value decomposition (SVD) of a real upper bidiagonal N-by-M matrix<br>
*> B with diagonal D and offdiagonal E, where M = N + SQRE. The<br>
*> algorithm computes the singular values in the SVD B = U * S * VT.<br>
*> The orthogonal matrices U and VT are optionally computed in<br>
*> compact form.<br>
*><br>
*> A related subroutine, SLASD0, computes the singular values and<br>
*> the singular vectors in explicit form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ICOMPQ<br>
*> \verbatim<br>
*>          ICOMPQ is INTEGER<br>
*>         Specifies whether singular vectors are to be computed<br>
*>         in compact form, as follows<br>
*>         = 0: Compute singular values only.<br>
*>         = 1: Compute singular vectors of upper bidiagonal<br>
*>              matrix in compact form.<br>
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
*>         The row dimension of the upper bidiagonal matrix. This is<br>
*>         also the dimension of the main diagonal array D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SQRE<br>
*> \verbatim<br>
*>          SQRE is INTEGER<br>
*>         Specifies the column dimension of the bidiagonal matrix.<br>
*>         = 0: The bidiagonal matrix has column dimension M = N;<br>
*>         = 1: The bidiagonal matrix has column dimension M = N + 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension ( N )<br>
*>         On entry D contains the main diagonal of the bidiagonal<br>
*>         matrix. On exit D, if INFO = 0, contains its singular values.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension ( M-1 )<br>
*>         Contains the subdiagonal entries of the bidiagonal matrix.<br>
*>         On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is REAL array,<br>
*>         dimension ( LDU, SMLSIZ ) if ICOMPQ = 1, and not referenced<br>
*>         if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left<br>
*>         singular vector matrices of all subproblems at the bottom<br>
*>         level.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER, LDU = > N.<br>
*>         The leading dimension of arrays U, VT, DIFL, DIFR, POLES,<br>
*>         GIVNUM, and Z.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VT<br>
*> \verbatim<br>
*>          VT is REAL array,<br>
*>         dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced<br>
*>         if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT**T contains the right<br>
*>         singular vector matrices of all subproblems at the bottom<br>
*>         level.<br>
*> \endverbatim<br>
*><br>
*> \param[out] K<br>
*> \verbatim<br>
*>          K is INTEGER array, dimension ( N )<br>
*>         if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0.<br>
*>         If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th<br>
*>         secular equation on the computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIFL<br>
*> \verbatim<br>
*>          DIFL is REAL array, dimension ( LDU, NLVL ),<br>
*>         where NLVL = floor(log_2 (N/SMLSIZ))).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIFR<br>
*> \verbatim<br>
*>          DIFR is REAL array,<br>
*>                  dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and<br>
*>                  dimension ( N ) if ICOMPQ = 0.<br>
*>         If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1)<br>
*>         record distances between singular values on the I-th<br>
*>         level and singular values on the (I -1)-th level, and<br>
*>         DIFR(1:N, 2 * I ) contains the normalizing factors for<br>
*>         the right singular vector matrix. See SLASD8 for details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array,<br>
*>                  dimension ( LDU, NLVL ) if ICOMPQ = 1 and<br>
*>                  dimension ( N ) if ICOMPQ = 0.<br>
*>         The first K elements of Z(1, I) contain the components of<br>
*>         the deflation-adjusted updating row vector for subproblems<br>
*>         on the I-th level.<br>
*> \endverbatim<br>
*><br>
*> \param[out] POLES<br>
*> \verbatim<br>
*>          POLES is REAL array,<br>
*>         dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not referenced<br>
*>         if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and<br>
*>         POLES(1, 2*I) contain  the new and old singular values<br>
*>         involved in the secular equations on the I-th level.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVPTR<br>
*> \verbatim<br>
*>          GIVPTR is INTEGER array,<br>
*>         dimension ( N ) if ICOMPQ = 1, and not referenced if<br>
*>         ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR( I ) records<br>
*>         the number of Givens rotations performed on the I-th<br>
*>         problem on the computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVCOL<br>
*> \verbatim<br>
*>          GIVCOL is INTEGER array,<br>
*>         dimension ( LDGCOL, 2 * NLVL ) if ICOMPQ = 1, and not<br>
*>         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,<br>
*>         GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations<br>
*>         of Givens rotations performed on the I-th level on the<br>
*>         computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDGCOL<br>
*> \verbatim<br>
*>          LDGCOL is INTEGER, LDGCOL = > N.<br>
*>         The leading dimension of arrays GIVCOL and PERM.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PERM<br>
*> \verbatim<br>
*>          PERM is INTEGER array, dimension ( LDGCOL, NLVL )<br>
*>         if ICOMPQ = 1, and not referenced<br>
*>         if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records<br>
*>         permutations done on the I-th level of the computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[out] GIVNUM<br>
*> \verbatim<br>
*>          GIVNUM is REAL array,<br>
*>         dimension ( LDU,  2 * NLVL ) if ICOMPQ = 1, and not<br>
*>         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,<br>
*>         GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S-<br>
*>         values of Givens rotations performed on the I-th level on<br>
*>         the computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is REAL array,<br>
*>         dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0.<br>
*>         If ICOMPQ = 1 and the I-th subproblem is not square, on exit,<br>
*>         C( I ) contains the C-value of a Givens rotation related to<br>
*>         the right null space of the I-th subproblem.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension ( N ) if<br>
*>         ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1<br>
*>         and the I-th subproblem is not square, on exit, S( I )<br>
*>         contains the S-value of a Givens rotation related to<br>
*>         the right null space of the I-th subproblem.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension<br>
*>         (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (7*N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = 1, a singular value did not converge<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasda_(INTEGER ICOMPQ,INTEGER SMLSIZ,INTEGER N,INTEGER SQRE,float[] D,float[] E,float[] U,INTEGER LDU,float[] VT,int[] K,float[] DIFL,float[] DIFR,float[] Z,float[] POLES,int[] GIVPTR,int[] GIVCOL,INTEGER LDGCOL,int[] PERM,float[] GIVNUM,float[] C,float[] S,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SLASDQ computes the SVD of a real bidiagonal matrix with diagonal d and off-diagonal e. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASDQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasdq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasdq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasdq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASDQ( UPLO, SQRE, N, NCVT, NRU, NCC, D, E, VT, LDVT,<br>
*                          U, LDU, C, LDC, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU, SQRE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               C( LDC, * ), D( * ), E( * ), U( LDU, * ),<br>
*      $                   VT( LDVT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASDQ computes the singular value decomposition (SVD) of a real<br>
*> (upper or lower) bidiagonal matrix with diagonal D and offdiagonal<br>
*> E, accumulating the transformations if desired. Letting B denote<br>
*> the input bidiagonal matrix, the algorithm computes orthogonal<br>
*> matrices Q and P such that B = Q * S * P**T (P**T denotes the transpose<br>
*> of P). The singular values S are overwritten on D.<br>
*><br>
*> The input matrix U  is changed to U  * Q  if desired.<br>
*> The input matrix VT is changed to P**T * VT if desired.<br>
*> The input matrix C  is changed to Q**T * C  if desired.<br>
*><br>
*> See "Computing  Small Singular Values of Bidiagonal Matrices With<br>
*> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,<br>
*> LAPACK Working Note #3, for a detailed description of the algorithm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>        On entry, UPLO specifies whether the input bidiagonal matrix<br>
*>        is upper or lower bidiagonal, and whether it is square are<br>
*>        not.<br>
*>           UPLO = 'U' or 'u'   B is upper bidiagonal.<br>
*>           UPLO = 'L' or 'l'   B is lower bidiagonal.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SQRE<br>
*> \verbatim<br>
*>          SQRE is INTEGER<br>
*>        = 0: then the input matrix is N-by-N.<br>
*>        = 1: then the input matrix is N-by-(N+1) if UPLU = 'U' and<br>
*>             (N+1)-by-N if UPLU = 'L'.<br>
*><br>
*>        The bidiagonal matrix has<br>
*>        N = NL + NR + 1 rows and<br>
*>        M = N + SQRE >= N columns.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>        On entry, N specifies the number of rows and columns<br>
*>        in the matrix. N must be at least 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NCVT<br>
*> \verbatim<br>
*>          NCVT is INTEGER<br>
*>        On entry, NCVT specifies the number of columns of<br>
*>        the matrix VT. NCVT must be at least 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRU<br>
*> \verbatim<br>
*>          NRU is INTEGER<br>
*>        On entry, NRU specifies the number of rows of<br>
*>        the matrix U. NRU must be at least 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NCC<br>
*> \verbatim<br>
*>          NCC is INTEGER<br>
*>        On entry, NCC specifies the number of columns of<br>
*>        the matrix C. NCC must be at least 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>        On entry, D contains the diagonal entries of the<br>
*>        bidiagonal matrix whose SVD is desired. On normal exit,<br>
*>        D contains the singular values in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array.<br>
*>        dimension is (N-1) if SQRE = 0 and N if SQRE = 1.<br>
*>        On entry, the entries of E contain the offdiagonal entries<br>
*>        of the bidiagonal matrix whose SVD is desired. On normal<br>
*>        exit, E will contain 0. If the algorithm does not converge,<br>
*>        D and E will contain the diagonal and superdiagonal entries<br>
*>        of a bidiagonal matrix orthogonally equivalent to the one<br>
*>        given as input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VT<br>
*> \verbatim<br>
*>          VT is REAL array, dimension (LDVT, NCVT)<br>
*>        On entry, contains a matrix which on exit has been<br>
*>        premultiplied by P**T, dimension N-by-NCVT if SQRE = 0<br>
*>        and (N+1)-by-NCVT if SQRE = 1 (not referenced if NCVT=0).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>        On entry, LDVT specifies the leading dimension of VT as<br>
*>        declared in the calling (sub) program. LDVT must be at<br>
*>        least 1. If NCVT is nonzero LDVT must also be at least N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] U<br>
*> \verbatim<br>
*>          U is REAL array, dimension (LDU, N)<br>
*>        On entry, contains a  matrix which on exit has been<br>
*>        postmultiplied by Q, dimension NRU-by-N if SQRE = 0<br>
*>        and NRU-by-(N+1) if SQRE = 1 (not referenced if NRU=0).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>        On entry, LDU  specifies the leading dimension of U as<br>
*>        declared in the calling (sub) program. LDU must be at<br>
*>        least max( 1, NRU ) .<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC, NCC)<br>
*>        On entry, contains an N-by-NCC matrix which on exit<br>
*>        has been premultiplied by Q**T  dimension N-by-NCC if SQRE = 0<br>
*>        and (N+1)-by-NCC if SQRE = 1 (not referenced if NCC=0).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>        On entry, LDC  specifies the leading dimension of C as<br>
*>        declared in the calling (sub) program. LDC must be at<br>
*>        least 1. If NCC is nonzero, LDC must also be at least N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (4*N)<br>
*>        Workspace. Only referenced if one of NCVT, NRU, or NCC is<br>
*>        nonzero, and if N is at least 2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>        On exit, a value of 0 indicates a successful exit.<br>
*>        If INFO < 0, argument number -INFO is illegal.<br>
*>        If INFO > 0, the algorithm did not converge, and INFO<br>
*>        specifies how many superdiagonals did not converge.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasdq_(CHARACTER UPLO,INTEGER SQRE,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,float[] D,float[] E,float[] VT,INTEGER LDVT,float[] U,INTEGER LDU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLASDT creates a tree of subproblems for bidiagonal divide and conquer. Used by sbdsdc.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASDT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasdt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasdt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasdt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LVL, MSUB, N, ND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            INODE( * ), NDIML( * ), NDIMR( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASDT creates a tree of subproblems for bidiagonal divide and<br>
*> conquer.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          On entry, the number of diagonal elements of the<br>
*>          bidiagonal matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] LVL<br>
*> \verbatim<br>
*>          LVL is INTEGER<br>
*>          On exit, the number of levels on the computation tree.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ND<br>
*> \verbatim<br>
*>          ND is INTEGER<br>
*>          On exit, the number of nodes on the tree.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INODE<br>
*> \verbatim<br>
*>          INODE is INTEGER array, dimension ( N )<br>
*>          On exit, centers of subproblems.<br>
*> \endverbatim<br>
*><br>
*> \param[out] NDIML<br>
*> \verbatim<br>
*>          NDIML is INTEGER array, dimension ( N )<br>
*>          On exit, row dimensions of left children.<br>
*> \endverbatim<br>
*><br>
*> \param[out] NDIMR<br>
*> \verbatim<br>
*>          NDIMR is INTEGER array, dimension ( N )<br>
*>          On exit, row dimensions of right children.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MSUB<br>
*> \verbatim<br>
*>          MSUB is INTEGER<br>
*>          On entry, the maximum row dimension each subproblem at the<br>
*>          bottom of the tree can be of.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void slasdt_(INTEGER N,INTEGER LVL,INTEGER ND,int[] INODE,int[] NDIML,int[] NDIMR,INTEGER MSUB);
/**
*> \brief \b SLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASET + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaset.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaset.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaset.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            LDA, M, N<br>
*       REAL               ALPHA, BETA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASET initializes an m-by-n matrix A to BETA on the diagonal and<br>
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
*>          = 'U':      Upper triangular part is set; the strictly lower<br>
*>                      triangular part of A is not changed.<br>
*>          = 'L':      Lower triangular part is set; the strictly upper<br>
*>                      triangular part of A is not changed.<br>
*>          Otherwise:  All of the matrix A is set.<br>
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
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>          The constant to which the offdiagonal elements are to be set.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>          The constant to which the diagonal elements are to be set.<br>
*> \endverbatim<br>
*><br>
*> \param[out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
*>          On exit, the leading m-by-n submatrix of A is set as follows:<br>
*><br>
*>          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,<br>
*>          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,<br>
*>          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,<br>
*><br>
*>          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slaset_(CHARACTER UPLO,INTEGER M,INTEGER N,REAL ALPHA,REAL BETA,float[] A,INTEGER LDA);
/**
*> \brief \b SLASQ1 computes the singular values of a real square bidiagonal matrix. Used by sbdsqr.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASQ1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASQ1( N, D, E, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASQ1 computes the singular values of a real N-by-N bidiagonal<br>
*> matrix with diagonal D and off-diagonal E. The singular values<br>
*> are computed to high relative accuracy, in the absence of<br>
*> denormalization, underflow and overflow. The algorithm was first<br>
*> presented in<br>
*><br>
*> "Accurate singular values and differential qd algorithms" by K. V.<br>
*> Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,<br>
*> 1994,<br>
*><br>
*> and the present implementation is described in "An implementation of<br>
*> the dqds Algorithm (Positive Case)", LAPACK Working Note.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>        The number of rows and columns in the matrix. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>        On entry, D contains the diagonal elements of the<br>
*>        bidiagonal matrix whose SVD is desired. On normal exit,<br>
*>        D contains the singular values in decreasing order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
*>        On entry, elements E(1:N-1) contain the off-diagonal elements<br>
*>        of the bidiagonal matrix whose SVD is desired.<br>
*>        On exit, E is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>        = 0: successful exit<br>
*>        < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>        > 0: the algorithm failed<br>
*>             = 1, a split was marked by a positive value in E<br>
*>             = 2, current block of Z not diagonalized after 100*N<br>
*>                  iterations (in inner while loop)  On exit D and E<br>
*>                  represent a matrix with the same singular values<br>
*>                  which the calling subroutine could use to finish the<br>
*>                  computation, or even feed back into SLASQ1<br>
*>             = 3, termination criterion of outer while loop not met <br>
*>                  (program created more than N unreduced blocks)<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slasq1_(INTEGER N,float[] D,float[] E,float[] WORK,INTEGER INFO);
/**
*> \brief \b SLASQ2 computes all the eigenvalues of the symmetric positive definite tridiagonal matrix associated with the qd Array Z to high relative accuracy. Used by sbdsqr and sstegr.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASQ2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASQ2( N, Z, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASQ2 computes all the eigenvalues of the symmetric positive <br>
*> definite tridiagonal matrix associated with the qd array Z to high<br>
*> relative accuracy are computed to high relative accuracy, in the<br>
*> absence of denormalization, underflow and overflow.<br>
*><br>
*> To see the relation of Z to the tridiagonal matrix, let L be a<br>
*> unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and<br>
*> let U be an upper bidiagonal matrix with 1's above and diagonal<br>
*> Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the<br>
*> symmetric tridiagonal to which it is similar.<br>
*><br>
*> Note : SLASQ2 defines a logical variable, IEEE, which is true<br>
*> on machines which follow ieee-754 floating-point standard in their<br>
*> handling of infinities and NaNs, and false otherwise. This variable<br>
*> is passed to SLASQ3.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>        The number of rows and columns in the matrix. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( 4*N )<br>
*>        On entry Z holds the qd array. On exit, entries 1 to N hold<br>
*>        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the<br>
*>        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If<br>
*>        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )<br>
*>        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of<br>
*>        shifts that failed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>        = 0: successful exit<br>
*>        < 0: if the i-th argument is a scalar and had an illegal<br>
*>             value, then INFO = -i, if the i-th argument is an<br>
*>             array and the j-entry had an illegal value, then<br>
*>             INFO = -(i*100+j)<br>
*>        > 0: the algorithm failed<br>
*>              = 1, a split was marked by a positive value in E<br>
*>              = 2, current block of Z not diagonalized after 100*N<br>
*>                   iterations (in inner while loop).  On exit Z holds<br>
*>                   a qd array with the same eigenvalues as the given Z.<br>
*>              = 3, termination criterion of outer while loop not met <br>
*>                   (program created more than N unreduced blocks)<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Local Variables: I0:N0 defines a current unreduced segment of Z.<br>
*>  The shifts are accumulated in SIGMA. Iteration count is in ITER.<br>
*>  Ping-pong is controlled by PP (alternates between 0 and 1).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slasq2_(INTEGER N,float[] Z,INTEGER INFO);
/**
*> \brief \b SLASQ3 checks for deflation, computes a shift and calls dqds. Used by sbdsqr.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASQ3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,<br>
*                          ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,<br>
*                          DN2, G, TAU )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            IEEE<br>
*       INTEGER            I0, ITER, N0, NDIV, NFAIL, PP<br>
*       REAL               DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G,<br>
*      $                   QMAX, SIGMA, TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.<br>
*> In case of failure it changes shifts, and tries again until output<br>
*> is positive.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] I0<br>
*> \verbatim<br>
*>          I0 is INTEGER<br>
*>         First index.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] N0<br>
*> \verbatim<br>
*>          N0 is INTEGER<br>
*>         Last index.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( 4*N0 )<br>
*>         Z holds the qd array.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] PP<br>
*> \verbatim<br>
*>          PP is INTEGER<br>
*>         PP=0 for ping, PP=1 for pong.<br>
*>         PP=2 indicates that flipping was applied to the Z array   <br>
*>         and that the initial tests for deflation should not be <br>
*>         performed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DMIN<br>
*> \verbatim<br>
*>          DMIN is REAL<br>
*>         Minimum value of d.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SIGMA<br>
*> \verbatim<br>
*>          SIGMA is REAL<br>
*>         Sum of shifts used in current segment.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DESIG<br>
*> \verbatim<br>
*>          DESIG is REAL<br>
*>         Lower order part of SIGMA<br>
*> \endverbatim<br>
*><br>
*> \param[in] QMAX<br>
*> \verbatim<br>
*>          QMAX is REAL<br>
*>         Maximum value of q.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] NFAIL<br>
*> \verbatim<br>
*>          NFAIL is INTEGER<br>
*>         Increment NFAIL by 1 each time the shift was too big.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ITER<br>
*> \verbatim<br>
*>          ITER is INTEGER<br>
*>         Increment ITER by 1 for each iteration.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] NDIV<br>
*> \verbatim<br>
*>          NDIV is INTEGER<br>
*>         Increment NDIV by 1 for each division.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IEEE<br>
*> \verbatim<br>
*>          IEEE is LOGICAL<br>
*>         Flag for IEEE or non IEEE arithmetic (passed to SLASQ5).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] TTYPE<br>
*> \verbatim<br>
*>          TTYPE is INTEGER<br>
*>         Shift type.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DMIN1<br>
*> \verbatim<br>
*>          DMIN1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DMIN2<br>
*> \verbatim<br>
*>          DMIN2 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DN<br>
*> \verbatim<br>
*>          DN is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DN1<br>
*> \verbatim<br>
*>          DN1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DN2<br>
*> \verbatim<br>
*>          DN2 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] G<br>
*> \verbatim<br>
*>          G is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] TAU<br>
*> \verbatim<br>
*>          TAU is REAL<br>
*><br>
*>         These are passed as arguments in order to save their values<br>
*>         between calls to SLASQ3.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slasq3_(INTEGER I0,INTEGER N0,float[] Z,INTEGER PP,REAL DMIN,REAL SIGMA,REAL DESIG,REAL QMAX,INTEGER NFAIL,INTEGER ITER,INTEGER NDIV,LOGICAL IEEE,INTEGER TTYPE,REAL DMIN1,REAL DMIN2,REAL DN,REAL DN1,REAL DN2,REAL G,REAL TAU);
/**
*> \brief \b SLASQ4 computes an approximation to the smallest eigenvalue using values of d from the previous transform. Used by sbdsqr.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASQ4 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq4.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq4.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq4.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,<br>
*                          DN1, DN2, TAU, TTYPE, G )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            I0, N0, N0IN, PP, TTYPE<br>
*       REAL               DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASQ4 computes an approximation TAU to the smallest eigenvalue<br>
*> using values of d from the previous transform.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] I0<br>
*> \verbatim<br>
*>          I0 is INTEGER<br>
*>        First index.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N0<br>
*> \verbatim<br>
*>          N0 is INTEGER<br>
*>        Last index.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( 4*N0 )<br>
*>        Z holds the qd array.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PP<br>
*> \verbatim<br>
*>          PP is INTEGER<br>
*>        PP=0 for ping, PP=1 for pong.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N0IN<br>
*> \verbatim<br>
*>          N0IN is INTEGER<br>
*>        The value of N0 at start of EIGTEST.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DMIN<br>
*> \verbatim<br>
*>          DMIN is REAL<br>
*>        Minimum value of d.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DMIN1<br>
*> \verbatim<br>
*>          DMIN1 is REAL<br>
*>        Minimum value of d, excluding D( N0 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] DMIN2<br>
*> \verbatim<br>
*>          DMIN2 is REAL<br>
*>        Minimum value of d, excluding D( N0 ) and D( N0-1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] DN<br>
*> \verbatim<br>
*>          DN is REAL<br>
*>        d(N)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DN1<br>
*> \verbatim<br>
*>          DN1 is REAL<br>
*>        d(N-1)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DN2<br>
*> \verbatim<br>
*>          DN2 is REAL<br>
*>        d(N-2)<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is REAL<br>
*>        This is the shift.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TTYPE<br>
*> \verbatim<br>
*>          TTYPE is INTEGER<br>
*>        Shift type.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] G<br>
*> \verbatim<br>
*>          G is REAL<br>
*>        G is passed as an argument in order to save its value between<br>
*>        calls to SLASQ4.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  CNST1 = 9/16<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slasq4_(INTEGER I0,INTEGER N0,float[] Z,INTEGER PP,INTEGER N0IN,REAL DMIN,REAL DMIN1,REAL DMIN2,REAL DN,REAL DN1,REAL DN2,REAL TAU,INTEGER TTYPE,REAL G);
/**
*> \brief \b SLASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASQ5 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq5.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq5.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq5.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN,<br>
*                          DNM1, DNM2, IEEE, EPS )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            IEEE<br>
*       INTEGER            I0, N0, PP<br>
*       REAL               EPS, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, SIGMA, TAU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASQ5 computes one dqds transform in ping-pong form, one<br>
*> version for IEEE machines another for non IEEE machines.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] I0<br>
*> \verbatim<br>
*>          I0 is INTEGER<br>
*>        First index.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N0<br>
*> \verbatim<br>
*>          N0 is INTEGER<br>
*>        Last index.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( 4*N )<br>
*>        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid<br>
*>        an extra argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PP<br>
*> \verbatim<br>
*>          PP is INTEGER<br>
*>        PP=0 for ping, PP=1 for pong.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is REAL<br>
*>        This is the shift.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIGMA<br>
*>          SIGMA is REAL<br>
*>        This is the accumulated shift up to this step.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DMIN<br>
*> \verbatim<br>
*>          DMIN is REAL<br>
*>        Minimum value of d.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DMIN1<br>
*> \verbatim<br>
*>          DMIN1 is REAL<br>
*>        Minimum value of d, excluding D( N0 ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DMIN2<br>
*> \verbatim<br>
*>          DMIN2 is REAL<br>
*>        Minimum value of d, excluding D( N0 ) and D( N0-1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DN<br>
*> \verbatim<br>
*>          DN is REAL<br>
*>        d(N0), the last value of d.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DNM1<br>
*> \verbatim<br>
*>          DNM1 is REAL<br>
*>        d(N0-1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DNM2<br>
*> \verbatim<br>
*>          DNM2 is REAL<br>
*>        d(N0-2).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IEEE<br>
*> \verbatim<br>
*>          IEEE is LOGICAL<br>
*>        Flag for IEEE or non IEEE arithmetic.<br>
*> \endverbatim<br>
*><br>
*> \param[in] EPS<br>
*> \verbatim<br>
*>         EPS is REAL<br>
*>        This is the value of epsilon used.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slasq5_(INTEGER I0,INTEGER N0,float[] Z,INTEGER PP,REAL TAU,REAL SIGMA,REAL DMIN,REAL DMIN1,REAL DMIN2,REAL DN,REAL DNM1,REAL DNM2,LOGICAL IEEE,REAL EPS);
/**
*> \brief \b SLASQ6 computes one dqd transform in ping-pong form. Used by sbdsqr and sstegr.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASQ6 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq6.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq6.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq6.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN,<br>
*                          DNM1, DNM2 )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            I0, N0, PP<br>
*       REAL               DMIN, DMIN1, DMIN2, DN, DNM1, DNM2<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               Z( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASQ6 computes one dqd (shift equal to zero) transform in<br>
*> ping-pong form, with protection against underflow and overflow.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] I0<br>
*> \verbatim<br>
*>          I0 is INTEGER<br>
*>        First index.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N0<br>
*> \verbatim<br>
*>          N0 is INTEGER<br>
*>        Last index.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension ( 4*N )<br>
*>        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid<br>
*>        an extra argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in] PP<br>
*> \verbatim<br>
*>          PP is INTEGER<br>
*>        PP=0 for ping, PP=1 for pong.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DMIN<br>
*> \verbatim<br>
*>          DMIN is REAL<br>
*>        Minimum value of d.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DMIN1<br>
*> \verbatim<br>
*>          DMIN1 is REAL<br>
*>        Minimum value of d, excluding D( N0 ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DMIN2<br>
*> \verbatim<br>
*>          DMIN2 is REAL<br>
*>        Minimum value of d, excluding D( N0 ) and D( N0-1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DN<br>
*> \verbatim<br>
*>          DN is REAL<br>
*>        d(N0), the last value of d.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DNM1<br>
*> \verbatim<br>
*>          DNM1 is REAL<br>
*>        d(N0-1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DNM2<br>
*> \verbatim<br>
*>          DNM2 is REAL<br>
*>        d(N0-2).<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slasq6_(INTEGER I0,INTEGER N0,float[] Z,INTEGER PP,REAL DMIN,REAL DMIN1,REAL DMIN2,REAL DN,REAL DNM1,REAL DNM2);
/**
*> \brief \b SLASR applies a sequence of plane rotations to a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIRECT, PIVOT, SIDE<br>
*       INTEGER            LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( * ), S( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASR applies a sequence of plane rotations to a real matrix A,<br>
*> from either the left or the right.<br>
*> <br>
*> When SIDE = 'L', the transformation takes the form<br>
*> <br>
*>    A := P*A<br>
*> <br>
*> and when SIDE = 'R', the transformation takes the form<br>
*> <br>
*>    A := A*P**T<br>
*> <br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slasr_(CHARACTER SIDE,CHARACTER PIVOT,CHARACTER DIRECT,INTEGER M,INTEGER N,float[] C,float[] S,float[] A,INTEGER LDA);
/**
*> \brief \b SLASRT sorts numbers in increasing or decreasing order.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASRT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasrt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasrt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasrt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASRT( ID, N, D, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          ID<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Sort the numbers in D in increasing order (if ID = 'I') or<br>
*> in decreasing order (if ID = 'D' ).<br>
*><br>
*> Use Quick Sort, reverting to Insertion sort on arrays of<br>
*> size <= 20. Dimension of STACK limits N to about 2**32.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ID<br>
*> \verbatim<br>
*>          ID is CHARACTER*1<br>
*>          = 'I': sort D in increasing order;<br>
*>          = 'D': sort D in decreasing order.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The length of the array D.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the array to be sorted.<br>
*>          On exit, D has been sorted into increasing order<br>
*>          (D(1) <= ... <= D(N) ) or into decreasing order<br>
*>          (D(1) >= ... >= D(N) ), depending on ID.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void slasrt_(CHARACTER ID,INTEGER N,float[] D,INTEGER INFO);
/**
*> \brief \b SLASSQ updates a sum of squares represented in scaled form.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASSQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slassq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slassq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slassq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       REAL               SCALE, SUMSQ<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASSQ  returns the values  scl  and  smsq  such that<br>
*><br>
*>    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,<br>
*><br>
*> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is<br>
*> assumed to be non-negative and  scl  returns the value<br>
*><br>
*>    scl = max( scale, abs( x( i ) ) ).<br>
*><br>
*> scale and sumsq must be supplied in SCALE and SUMSQ and<br>
*> scl and smsq are overwritten on SCALE and SUMSQ respectively.<br>
*><br>
*> The routine makes only one pass through the vector x.<br>
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
*>          X is REAL array, dimension (N)<br>
*>          The vector for which a scaled sum of squares is computed.<br>
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
*>          On exit, SCALE is overwritten with  scl , the scaling factor<br>
*>          for the sum of squares.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SUMSQ<br>
*> \verbatim<br>
*>          SUMSQ is REAL<br>
*>          On entry, the value  sumsq  in the equation above.<br>
*>          On exit, SUMSQ is overwritten with  smsq , the basic sum of<br>
*>          squares from which  scl  has been factored out.<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slassq_(INTEGER N,float[] X,INTEGER INCX,REAL SCALE,REAL SUMSQ);
/**
*> \brief \b SLASV2 computes the singular value decomposition of a 2-by-2 triangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASV2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasv2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasv2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasv2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASV2 computes the singular value decomposition of a 2-by-2<br>
*> triangular matrix<br>
*>    [  F   G  ]<br>
*>    [  0   H  ].<br>
*> On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the<br>
*> smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and<br>
*> right singular vectors for abs(SSMAX), giving the decomposition<br>
*><br>
*>    [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]<br>
*>    [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] F<br>
*> \verbatim<br>
*>          F is REAL<br>
*>          The (1,1) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] G<br>
*> \verbatim<br>
*>          G is REAL<br>
*>          The (1,2) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] H<br>
*> \verbatim<br>
*>          H is REAL<br>
*>          The (2,2) element of the 2-by-2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SSMIN<br>
*> \verbatim<br>
*>          SSMIN is REAL<br>
*>          abs(SSMIN) is the smaller singular value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SSMAX<br>
*> \verbatim<br>
*>          SSMAX is REAL<br>
*>          abs(SSMAX) is the larger singular value.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNL<br>
*> \verbatim<br>
*>          SNL is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSL<br>
*> \verbatim<br>
*>          CSL is REAL<br>
*>          The vector (CSL, SNL) is a unit left singular vector for the<br>
*>          singular value abs(SSMAX).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SNR<br>
*> \verbatim<br>
*>          SNR is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] CSR<br>
*> \verbatim<br>
*>          CSR is REAL<br>
*>          The vector (CSR, SNR) is a unit right singular vector for the<br>
*>          singular value abs(SSMAX).<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Any input parameter may be aliased with any output parameter.<br>
*><br>
*>  Barring over/underflow and assuming a guard digit in subtraction, all<br>
*>  output quantities are correct to within a few units in the last<br>
*>  place (ulps).<br>
*><br>
*>  In IEEE arithmetic, the code works correctly if one matrix element is<br>
*>  infinite.<br>
*><br>
*>  Overflow will not occur unless the largest singular value itself<br>
*>  overflows or is within a few ulps of overflow. (On machines with<br>
*>  partial overflow, like the Cray, overflow may occur if the largest<br>
*>  singular value is within a factor of 2 of overflow.)<br>
*><br>
*>  Underflow is harmless if underflow is gradual. Otherwise, results<br>
*>  may correspond to a matrix modified by perturbations of size near<br>
*>  the underflow threshold.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slasv2_(REAL F,REAL G,REAL H,REAL SSMIN,REAL SSMAX,REAL SNR,REAL CSR,REAL SNL,REAL CSL);
/**
*> \brief \b SLASWP performs a series of row interchanges on a general rectangular matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASWP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaswp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaswp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaswp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, K1, K2, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASWP performs a series of row interchanges on the matrix A.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void slaswp_(INTEGER N,float[] A,INTEGER LDA,INTEGER K1,INTEGER K2,int[] IPIV,INTEGER INCX);
/**
*> \brief \b SLASY2 solves the Sylvester matrix equation where the matrices are of order 1 or 2.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLASY2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasy2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasy2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasy2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,<br>
*                          LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            LTRANL, LTRANR<br>
*       INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2<br>
*       REAL               SCALE, XNORM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in<br>
*><br>
*>        op(TL)*X + ISGN*X*op(TR) = SCALE*B,<br>
*><br>
*> where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or<br>
*> -1.  op(T) = T or T**T, where T**T denotes the transpose of T.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] LTRANL<br>
*> \verbatim<br>
*>          LTRANL is LOGICAL<br>
*>          On entry, LTRANL specifies the op(TL):<br>
*>             = .FALSE., op(TL) = TL,<br>
*>             = .TRUE., op(TL) = TL**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LTRANR<br>
*> \verbatim<br>
*>          LTRANR is LOGICAL<br>
*>          On entry, LTRANR specifies the op(TR):<br>
*>            = .FALSE., op(TR) = TR,<br>
*>            = .TRUE., op(TR) = TR**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ISGN<br>
*> \verbatim<br>
*>          ISGN is INTEGER<br>
*>          On entry, ISGN specifies the sign of the equation<br>
*>          as described before. ISGN may only be 1 or -1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N1<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*>          On entry, N1 specifies the order of matrix TL.<br>
*>          N1 may only be 0, 1 or 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N2<br>
*> \verbatim<br>
*>          N2 is INTEGER<br>
*>          On entry, N2 specifies the order of matrix TR.<br>
*>          N2 may only be 0, 1 or 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TL<br>
*> \verbatim<br>
*>          TL is REAL array, dimension (LDTL,2)<br>
*>          On entry, TL contains an N1 by N1 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDTL<br>
*> \verbatim<br>
*>          LDTL is INTEGER<br>
*>          The leading dimension of the matrix TL. LDTL >= max(1,N1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TR<br>
*> \verbatim<br>
*>          TR is REAL array, dimension (LDTR,2)<br>
*>          On entry, TR contains an N2 by N2 matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDTR<br>
*> \verbatim<br>
*>          LDTR is INTEGER<br>
*>          The leading dimension of the matrix TR. LDTR >= max(1,N2).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,2)<br>
*>          On entry, the N1 by N2 matrix B contains the right-hand<br>
*>          side of the equation.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the matrix B. LDB >= max(1,N1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          On exit, SCALE contains the scale factor. SCALE is chosen<br>
*>          less than or equal to 1 to prevent the solution overflowing.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (LDX,2)<br>
*>          On exit, X contains the N1 by N2 solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the matrix X. LDX >= max(1,N1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] XNORM<br>
*> \verbatim<br>
*>          XNORM is REAL<br>
*>          On exit, XNORM is the infinity-norm of the solution.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          On exit, INFO is set to<br>
*>             0: successful exit.<br>
*>             1: TL and TR have too close eigenvalues, so TL or<br>
*>                TR is perturbed to get a nonsingular equation.<br>
*>          NOTE: In the interests of speed, this routine does not<br>
*>                check the inputs for errors.<br>
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
*> \ingroup realSYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slasy2_(LOGICAL LTRANL,LOGICAL LTRANR,INTEGER ISGN,INTEGER N1,INTEGER N2,float[] TL,INTEGER LDTL,float[] TR,INTEGER LDTR,float[] B,INTEGER LDB,REAL SCALE,float[] X,INTEGER LDX,REAL XNORM,INTEGER INFO);
/**
*> \brief \b SLASYF computes a partial factorization of a real symmetric matrix using the Bunch-Kaufman diagonal pivoting method.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download SLASYF + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KB, LDA, LDW, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), W( LDW, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASYF computes a partial factorization of a real symmetric matrix A<br>
*> using the Bunch-Kaufman diagonal pivoting method. The partial<br>
*> factorization has the form:<br>
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
*> SLASYF is an auxiliary routine called by SSYTRF. It uses blocked code<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          W is REAL array, dimension (LDW,NB)<br>
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
*> \ingroup realSYcomputational<br>
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
	public void slasyf_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] IPIV,float[] W,INTEGER LDW,INTEGER INFO);
/**
*> \brief \b SLASYF_ROOK computes a partial factorization of a real symmetric matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download SLASYF_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLASYF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KB, LDA, LDW, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), W( LDW, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLASYF_ROOK computes a partial factorization of a real symmetric<br>
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
*> SLASYF_ROOK is an auxiliary routine called by SSYTRF_ROOK. It uses<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          W is REAL array, dimension (LDW,NB)<br>
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
*> \ingroup realSYcomputational<br>
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
	public void slasyf_rook_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] IPIV,float[] W,INTEGER LDW,INTEGER INFO);
/**
*> \brief \b SLATBS solves a triangular banded system of equations.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLATBS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatbs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatbs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatbs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X,<br>
*                          SCALE, CNORM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
*       REAL               SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), CNORM( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLATBS solves one of the triangular systems<br>
*><br>
*>    A *x = s*b  or  A**T*x = s*b<br>
*><br>
*> with scaling to prevent overflow, where A is an upper or lower<br>
*> triangular band matrix.  Here A**T denotes the transpose of A, x and b<br>
*> are n-element vectors, and s is a scaling factor, usually less than<br>
*> or equal to 1, chosen so that the components of x will be less than<br>
*> the overflow threshold.  If the unscaled problem will not cause<br>
*> overflow, the Level 2 BLAS routine STBSV is called.  If the matrix A<br>
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
*>          = 'N':  Solve A * x = s*b  (No transpose)<br>
*>          = 'T':  Solve A**T* x = s*b  (Transpose)<br>
*>          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose)<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*>          X is REAL array, dimension (N)<br>
*>          On entry, the right hand side b of the triangular system.<br>
*>          On exit, X is overwritten by the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          The scaling factor s for the triangular system<br>
*>             A * x = s*b  or  A**T* x = s*b.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  A rough bound on x is computed; if that is less than overflow, STBSV<br>
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
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine STBSV if the<br>
*>  reciprocal of the largest M(j), j=1,..,n, is larger than<br>
*>  max(underflow, 1/overflow).<br>
*><br>
*>  The bound on x(j) is also used to determine when a step in the<br>
*>  columnwise method can be performed without fear of overflow.  If<br>
*>  the computed bound is greater than a large constant, x is scaled to<br>
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to<br>
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.<br>
*><br>
*>  Similarly, a row-wise scheme is used to solve A**T*x = b.  The basic<br>
*>  algorithm for A upper triangular is<br>
*><br>
*>       for j = 1, ..., n<br>
*>            x(j) := ( b(j) - A[1:j-1,j]**T * x[1:j-1] ) / A(j,j)<br>
*>       end<br>
*><br>
*>  We simultaneously compute two bounds<br>
*>       G(j) = bound on ( b(i) - A[1:i-1,i]**T * x[1:i-1] ), 1<=i<=j<br>
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
*>  and we can safely call STBSV if 1/M(n) and 1/G(n) are both greater<br>
*>  than max(underflow, 1/overflow).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slatbs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
/**
*> \brief \b SLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contribution to the reciprocal Dif-estimate.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLATDF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatdf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatdf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatdf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV,<br>
*                          JPIV )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IJOB, LDZ, N<br>
*       REAL               RDSCAL, RDSUM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), JPIV( * )<br>
*       REAL               RHS( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLATDF uses the LU factorization of the n-by-n matrix Z computed by<br>
*> SGETC2 and computes a contribution to the reciprocal Dif-estimate<br>
*> by solving Z * x = b for x, and choosing the r.h.s. b such that<br>
*> the norm of x is as large as possible. On entry RHS = b holds the<br>
*> contribution from earlier solved sub-systems, and on return RHS = x.<br>
*><br>
*> The factorization of Z returned by SGETC2 has the form Z = P*L*U*Q,<br>
*> where P and Q are permutation matrices. L is lower triangular with<br>
*> unit diagonal elements and U is upper triangular.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] IJOB<br>
*> \verbatim<br>
*>          IJOB is INTEGER<br>
*>          IJOB = 2: First compute an approximative null-vector e<br>
*>              of Z using SGECON, e is normalized and solve for<br>
*>              Zx = +-e - f with the sign giving the greater value<br>
*>              of 2-norm(x). About 5 times as expensive as Default.<br>
*>          IJOB .ne. 2: Local look ahead strategy where all entries of<br>
*>              the r.h.s. b is chosen as either +1 or -1 (Default).<br>
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
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          On entry, the LU part of the factorization of the n-by-n<br>
*>          matrix Z computed by SGETC2:  Z = P * L * U * Q<br>
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
*>          RHS is REAL array, dimension N.<br>
*>          On entry, RHS contains contributions from other subsystems.<br>
*>          On exit, RHS contains the solution of the subsystem with<br>
*>          entries acoording to the value of IJOB (see above).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RDSUM<br>
*> \verbatim<br>
*>          RDSUM is REAL<br>
*>          On entry, the sum of squares of computed contributions to<br>
*>          the Dif-estimate under computation by STGSYL, where the<br>
*>          scaling factor RDSCAL (see below) has been factored out.<br>
*>          On exit, the corresponding sum of squares updated with the<br>
*>          contributions from the current sub-system.<br>
*>          If TRANS = 'T' RDSUM is not touched.<br>
*>          NOTE: RDSUM only makes sense when STGSY2 is called by STGSYL.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RDSCAL<br>
*> \verbatim<br>
*>          RDSCAL is REAL<br>
*>          On entry, scaling factor used to prevent overflow in RDSUM.<br>
*>          On exit, RDSCAL is updated w.r.t. the current contributions<br>
*>          in RDSUM.<br>
*>          If TRANS = 'T', RDSCAL is not touched.<br>
*>          NOTE: RDSCAL only makes sense when STGSY2 is called by<br>
*>                STGSYL.<br>
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
*> \ingroup realOTHERauxiliary<br>
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
*> \verbatim<br>
*><br>
*><br>
*>  [1] Bo Kagstrom and Lars Westin,<br>
*>      Generalized Schur Methods with Condition Estimators for<br>
*>      Solving the Generalized Sylvester Equation, IEEE Transactions<br>
*>      on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751.<br>
*><br>
*>  [2] Peter Poromaa,<br>
*>      On Efficient and Robust Estimators for the Separation<br>
*>      between two Regular Matrix Pairs with Applications in<br>
*>      Condition Estimation. Report IMINF-95.05, Departement of<br>
*>      Computing Science, Umea University, S-901 87 Umea, Sweden, 1995.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slatdf_(INTEGER IJOB,INTEGER N,float[] Z,INTEGER LDZ,float[] RHS,REAL RDSUM,REAL RDSCAL,int[] IPIV,int[] JPIV);
/**
*> \brief \b SLATPS solves a triangular system of equations with the matrix held in packed storage.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLATPS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatps.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatps.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatps.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE,<br>
*                          CNORM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO<br>
*       INTEGER            INFO, N<br>
*       REAL               SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), CNORM( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLATPS solves one of the triangular systems<br>
*><br>
*>    A *x = s*b  or  A**T*x = s*b<br>
*><br>
*> with scaling to prevent overflow, where A is an upper or lower<br>
*> triangular matrix stored in packed form.  Here A**T denotes the<br>
*> transpose of A, x and b are n-element vectors, and s is a scaling<br>
*> factor, usually less than or equal to 1, chosen so that the<br>
*> components of x will be less than the overflow threshold.  If the<br>
*> unscaled problem will not cause overflow, the Level 2 BLAS routine<br>
*> STPSV is called. If the matrix A is singular (A(j,j) = 0 for some j),<br>
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
*>          = 'N':  Solve A * x = s*b  (No transpose)<br>
*>          = 'T':  Solve A**T* x = s*b  (Transpose)<br>
*>          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose)<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangular matrix A, packed columnwise in<br>
*>          a linear array.  The j-th column of A is stored in the array<br>
*>          AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (N)<br>
*>          On entry, the right hand side b of the triangular system.<br>
*>          On exit, X is overwritten by the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          The scaling factor s for the triangular system<br>
*>             A * x = s*b  or  A**T* x = s*b.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  A rough bound on x is computed; if that is less than overflow, STPSV<br>
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
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine STPSV if the<br>
*>  reciprocal of the largest M(j), j=1,..,n, is larger than<br>
*>  max(underflow, 1/overflow).<br>
*><br>
*>  The bound on x(j) is also used to determine when a step in the<br>
*>  columnwise method can be performed without fear of overflow.  If<br>
*>  the computed bound is greater than a large constant, x is scaled to<br>
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to<br>
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.<br>
*><br>
*>  Similarly, a row-wise scheme is used to solve A**T*x = b.  The basic<br>
*>  algorithm for A upper triangular is<br>
*><br>
*>       for j = 1, ..., n<br>
*>            x(j) := ( b(j) - A[1:j-1,j]**T * x[1:j-1] ) / A(j,j)<br>
*>       end<br>
*><br>
*>  We simultaneously compute two bounds<br>
*>       G(j) = bound on ( b(i) - A[1:i-1,i]**T * x[1:i-1] ), 1<=i<=j<br>
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
*>  and we can safely call STPSV if 1/M(n) and 1/G(n) are both greater<br>
*>  than max(underflow, 1/overflow).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slatps_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,float[] AP,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
/**
*> \brief \b SLATRD reduces the first nb rows and columns of a symmetric/Hermitian matrix A to real tridiagonal form by an orthogonal similarity transformation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLATRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            LDA, LDW, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), E( * ), TAU( * ), W( LDW, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLATRD reduces NB rows and columns of a real symmetric matrix A to<br>
*> symmetric tridiagonal form by an orthogonal similarity<br>
*> transformation Q**T * A * Q, and returns the matrices V and W which are<br>
*> needed to apply the transformation to the unreduced part of A.<br>
*><br>
*> If UPLO = 'U', SLATRD reduces the last NB rows and columns of a<br>
*> matrix, of which the upper triangle is supplied;<br>
*> if UPLO = 'L', SLATRD reduces the first NB rows and columns of a<br>
*> matrix, of which the lower triangle is supplied.<br>
*><br>
*> This is an auxiliary routine called by SSYTRD.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
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
*>            with the array TAU, represent the orthogonal matrix Q as a<br>
*>            product of elementary reflectors;<br>
*>          if UPLO = 'L', the first NB columns have been reduced to<br>
*>            tridiagonal form, with the diagonal elements overwriting<br>
*>            the diagonal elements of A; the elements below the diagonal<br>
*>            with the array TAU, represent the  orthogonal matrix Q as a<br>
*>            product of elementary reflectors.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= (1,N).<br>
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
*>          TAU is REAL array, dimension (N-1)<br>
*>          The scalar factors of the elementary reflectors, stored in<br>
*>          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (LDW,NB)<br>
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
*> \ingroup doubleOTHERauxiliary<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),<br>
*>  and tau in TAU(i).<br>
*><br>
*>  The elements of the vectors v together form the n-by-nb matrix V<br>
*>  which is needed, with W, to apply the transformation to the unreduced<br>
*>  part of the matrix, using a symmetric rank-2k update of the form:<br>
*>  A := A - V*W**T - W*V**T.<br>
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
	public void slatrd_(CHARACTER UPLO,INTEGER N,INTEGER NB,float[] A,INTEGER LDA,float[] E,float[] TAU,float[] W,INTEGER LDW);
/**
*> \brief \b SLATRS solves a triangular system of equations with the scale factor set to prevent overflow.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLATRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,<br>
*                          CNORM, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       REAL               SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), CNORM( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLATRS solves one of the triangular systems<br>
*><br>
*>    A *x = s*b  or  A**T*x = s*b<br>
*><br>
*> with scaling to prevent overflow.  Here A is an upper or lower<br>
*> triangular matrix, A**T denotes the transpose of A, x and b are<br>
*> n-element vectors, and s is a scaling factor, usually less than<br>
*> or equal to 1, chosen so that the components of x will be less than<br>
*> the overflow threshold.  If the unscaled problem will not cause<br>
*> overflow, the Level 2 BLAS routine STRSV is called.  If the matrix A<br>
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
*>          = 'N':  Solve A * x = s*b  (No transpose)<br>
*>          = 'T':  Solve A**T* x = s*b  (Transpose)<br>
*>          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose)<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          X is REAL array, dimension (N)<br>
*>          On entry, the right hand side b of the triangular system.<br>
*>          On exit, X is overwritten by the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          The scaling factor s for the triangular system<br>
*>             A * x = s*b  or  A**T* x = s*b.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  A rough bound on x is computed; if that is less than overflow, STRSV<br>
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
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine STRSV if the<br>
*>  reciprocal of the largest M(j), j=1,..,n, is larger than<br>
*>  max(underflow, 1/overflow).<br>
*><br>
*>  The bound on x(j) is also used to determine when a step in the<br>
*>  columnwise method can be performed without fear of overflow.  If<br>
*>  the computed bound is greater than a large constant, x is scaled to<br>
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to<br>
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.<br>
*><br>
*>  Similarly, a row-wise scheme is used to solve A**T*x = b.  The basic<br>
*>  algorithm for A upper triangular is<br>
*><br>
*>       for j = 1, ..., n<br>
*>            x(j) := ( b(j) - A[1:j-1,j]**T * x[1:j-1] ) / A(j,j)<br>
*>       end<br>
*><br>
*>  We simultaneously compute two bounds<br>
*>       G(j) = bound on ( b(i) - A[1:i-1,i]**T * x[1:i-1] ), 1<=i<=j<br>
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
*>  and we can safely call STRSV if 1/M(n) and 1/G(n) are both greater<br>
*>  than max(underflow, 1/overflow).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void slatrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,float[] A,INTEGER LDA,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
/**
*> \brief \b SLATRZ factors an upper trapezoidal matrix by means of orthogonal transformations.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLATRZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatrz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatrz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatrz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLATRZ( M, N, L, A, LDA, TAU, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            L, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLATRZ factors the M-by-(M+L) real upper trapezoidal matrix<br>
*> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z, by means<br>
*> of orthogonal transformations.  Z is an (M+L)-by-(M+L) orthogonal<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the leading M-by-N upper trapezoidal part of the<br>
*>          array A must contain the matrix to be factorized.<br>
*>          On exit, the leading M-by-M upper triangular part of A<br>
*>          contains the upper triangular matrix R, and elements N-L+1 to<br>
*>          N of the first M rows of A, with the array TAU, represent the<br>
*>          orthogonal matrix Z as a product of M elementary reflectors.<br>
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
*>          TAU is REAL array, dimension (M)<br>
*>          The scalar factors of the elementary reflectors.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (M)<br>
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
*> \ingroup realOTHERcomputational<br>
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
*>     T( k ) = I - tau*u( k )*u( k )**T,   u( k ) = (   1    ),<br>
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
	public void slatrz_(INTEGER M,INTEGER N,INTEGER L,float[] A,INTEGER LDA,float[] TAU,float[] WORK);
/**
*> \brief \b SLAUU2 computes the product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAUU2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slauu2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slauu2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slauu2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAUU2( UPLO, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAUU2 computes the product U * U**T or L**T * L, where the triangular<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the triangular factor U or L.<br>
*>          On exit, if UPLO = 'U', the upper triangle of A is<br>
*>          overwritten with the upper triangle of the product U * U**T;<br>
*>          if UPLO = 'L', the lower triangle of A is overwritten with<br>
*>          the lower triangle of the product L**T * L.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slauu2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b SLAUUM computes the product UUH or LHL, where U and L are upper or lower triangular matrices (blocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLAUUM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slauum.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slauum.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slauum.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLAUUM( UPLO, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLAUUM computes the product U * U**T or L**T * L, where the triangular<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the triangular factor U or L.<br>
*>          On exit, if UPLO = 'U', the upper triangle of A is<br>
*>          overwritten with the upper triangle of the product U * U**T;<br>
*>          if UPLO = 'L', the lower triangle of A is overwritten with<br>
*>          the lower triangle of the product L**T * L.<br>
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
*> \ingroup realOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void slauum_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b SLA_GBAMV performs a matrix-vector operation to calculate error bounds.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_GBAMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gbamv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gbamv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gbamv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X,<br>
*                             INCX, BETA, Y, INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               ALPHA, BETA<br>
*       INTEGER            INCX, INCY, LDAB, M, N, KL, KU, TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLA_GBAMV  performs one of the matrix-vector operations<br>
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
*>          AB is REAL array of DIMENSION ( LDAB, n )<br>
*>           Before entry, the leading m by n part of the array AB must<br>
*>           contain the matrix of coefficients.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>           On entry, LDA specifies the first dimension of AB as declared<br>
*>           in the calling (sub) program. LDAB must be at least<br>
*>           max( 1, m ).<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension<br>
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
*> \ingroup realGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sla_gbamv_(INTEGER TRANS,INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,REAL ALPHA,float[] AB,INTEGER LDAB,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
/**
*> \brief \b SLA_GBRCOND estimates the Skeel condition number for a general banded matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_GBRCOND + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gbrcond.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gbrcond.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gbrcond.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB,<br>
*                                  IPIV, CMODE, C, INFO, WORK, IWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            N, LDAB, LDAFB, INFO, KL, KU, CMODE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * ), IPIV( * )<br>
*       REAL               AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ),<br>
*      $                   C( * )<br>
*      ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLA_GBRCOND Estimates the Skeel condition number of  op(A) * op2(C)<br>
*>    where op2 is determined by CMODE as follows<br>
*>    CMODE =  1    op2(C) = C<br>
*>    CMODE =  0    op2(C) = I<br>
*>    CMODE = -1    op2(C) = inv(C)<br>
*>    The Skeel condition number  cond(A) = norminf( |inv(A)||A| )<br>
*>    is computed by computing scaling factors R such that<br>
*>    diag(R)*A*op2(C) is row equilibrated and computing the standard<br>
*>    infinity-norm condition number.<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*>          AFB is REAL array, dimension (LDAFB,N)<br>
*>     Details of the LU factorization of the band matrix A, as<br>
*>     computed by SGBTRF.  U is stored as an upper triangular<br>
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
*>     as computed by SGBTRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CMODE<br>
*> \verbatim<br>
*>          CMODE is INTEGER<br>
*>     Determines op2(C) in the formula op(A) * op2(C) as follows:<br>
*>     CMODE =  1    op2(C) = C<br>
*>     CMODE =  0    op2(C) = I<br>
*>     CMODE = -1    op2(C) = inv(C)<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * op2(C).<br>
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
*>          WORK is REAL array, dimension (5*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N).<br>
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
*> \ingroup realGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float sla_gbrcond_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,INTEGER CMODE,float[] C,INTEGER INFO,float[] WORK,int[] IWORK);
/**
*> \brief \b SLA_GBRFSX_EXTENDED improves the computed solution to a system of linear equations for general banded matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_GBRFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gbrfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gbrfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gbrfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_GBRFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, KL, KU,<br>
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
*       REAL               AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),<br>
*      $                   Y( LDY, * ), RES(*), DY(*), Y_TAIL(*)<br>
*       REAL               C( * ), AYB(*), RCOND, BERR_OUT(*),<br>
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
*> SLA_GBRFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by SGBRFSX to perform iterative refinement.<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*>          AFB is REAL array, dimension (LDAFB,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by SGBTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>     The leading dimension of the array AF.  LDAFB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     The pivot indices from the factorization A = P*L*U<br>
*>     as computed by SGBTRF; row i of the matrix was interchanged<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          Y is REAL array, dimension (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by SGBTRS.<br>
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
*>     or vector Z. This is computed by SLA_LIN_BERR.<br>
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
*>          RES is REAL array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace. This can be the same workspace passed for Y_TAIL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is REAL array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is REAL array, dimension (N)<br>
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
*>       < 0:  if INFO = -i, the ith argument to SGBTRS had an illegal<br>
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
*> \ingroup realGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sla_gbrfsx_extended_(INTEGER PREC_TYPE,INTEGER TRANS_TYPE,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,REAL SLAMCH,CHARACTER CHLA_TRANSTYPE);
/**
*> \brief \b SLA_GBRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a general banded matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_GBRPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gbrpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gbrpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gbrpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB,<br>
*                                   LDAFB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N, KL, KU, NCOLS, LDAB, LDAFB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), AFB( LDAFB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLA_GBRPVGRW computes the reciprocal pivot growth factor<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*>          AFB is REAL array, dimension (LDAFB,N)<br>
*>     Details of the LU factorization of the band matrix A, as<br>
*>     computed by SGBTRF.  U is stored as an upper triangular<br>
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
*> \ingroup realGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float sla_gbrpvgrw_(INTEGER N,INTEGER KL,INTEGER KU,INTEGER NCOLS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB);
/**
*> \brief \b SLA_GEAMV computes a matrix-vector product using a general matrix to calculate error bounds.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_GEAMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_geamv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_geamv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_geamv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA,<br>
*                              Y, INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               ALPHA, BETA<br>
*       INTEGER            INCX, INCY, LDA, M, N, TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLA_GEAMV  performs one of the matrix-vector operations<br>
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
*>          A is REAL array of DIMENSION ( LDA, n )<br>
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
*>          X is REAL array, dimension<br>
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
*>          Y is REAL<br>
*>           Array of DIMENSION at least<br>
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
*> \ingroup realGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sla_geamv_(INTEGER TRANS,INTEGER M,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
/**
*> \brief \b SLA_GERCOND estimates the Skeel condition number for a general matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_GERCOND + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gercond.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gercond.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gercond.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLA_GERCOND ( TRANS, N, A, LDA, AF, LDAF, IPIV,<br>
*                                   CMODE, C, INFO, WORK, IWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            N, LDA, LDAF, INFO, CMODE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ),<br>
*      $                   C( * )<br>
*      ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLA_GERCOND estimates the Skeel condition number of op(A) * op2(C)<br>
*>    where op2 is determined by CMODE as follows<br>
*>    CMODE =  1    op2(C) = C<br>
*>    CMODE =  0    op2(C) = I<br>
*>    CMODE = -1    op2(C) = inv(C)<br>
*>    The Skeel condition number cond(A) = norminf( |inv(A)||A| )<br>
*>    is computed by computing scaling factors R such that<br>
*>    diag(R)*A*op2(C) is row equilibrated and computing the standard<br>
*>    infinity-norm condition number.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by SGETRF.<br>
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
*>     as computed by SGETRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CMODE<br>
*> \verbatim<br>
*>          CMODE is INTEGER<br>
*>     Determines op2(C) in the formula op(A) * op2(C) as follows:<br>
*>     CMODE =  1    op2(C) = C<br>
*>     CMODE =  0    op2(C) = I<br>
*>     CMODE = -1    op2(C) = inv(C)<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * op2(C).<br>
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
*>          WORK is REAL array, dimension (3*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N).<br>
*>     Workspace.2<br>
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
*> \ingroup realGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float sla_gercond_(CHARACTER TRANS,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,INTEGER CMODE,float[] C,INTEGER INFO,float[] WORK,int[] IWORK);
/**
*> \brief \b SLA_GERFSX_EXTENDED improves the computed solution to a system of linear equations for general matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_GERFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gerfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gerfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gerfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A,<br>
*                                       LDA, AF, LDAF, IPIV, COLEQU, C, B,<br>
*                                       LDB, Y, LDY, BERR_OUT, N_NORMS,<br>
*                                       ERRS_N, ERRS_C, RES,<br>
*                                       AYB, DY, Y_TAIL, RCOND, ITHRESH,<br>
*                                       RTHRESH, DZ_UB, IGNORE_CWISE,<br>
*                                       INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE,<br>
*      $                   TRANS_TYPE, N_NORMS, ITHRESH<br>
*       LOGICAL            COLEQU, IGNORE_CWISE<br>
*       REAL               RTHRESH, DZ_UB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )<br>
*       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ),<br>
*      $                   ERRS_N( NRHS, * ),<br>
*      $                   ERRS_C( NRHS, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLA_GERFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by SGERFSX to perform iterative refinement.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by SGETRF.<br>
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
*>     as computed by SGETRF; row i of the matrix was interchanged<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          Y is REAL array, dimension (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by SGETRS.<br>
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
*>     or vector Z. This is computed by SLA_LIN_BERR.<br>
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
*>          RES is REAL array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace. This can be the same workspace passed for Y_TAIL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is REAL array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is REAL array, dimension (N)<br>
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
*>       < 0:  if INFO = -i, the ith argument to SGETRS had an illegal<br>
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
*> \ingroup realGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sla_gerfsx_extended_(INTEGER PREC_TYPE,INTEGER TRANS_TYPE,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERRS_N,float[] ERRS_C,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,REAL SLAMCH,CHARACTER CHLA_TRANSTYPE);
/**
*> \brief \b SLA_GERPVGRW<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_GERPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gerpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gerpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gerpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N, NCOLS, LDA, LDAF<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), AF( LDAF, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLA_GERPVGRW computes the reciprocal pivot growth factor<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization<br>
*>     A = P*L*U as computed by SGETRF.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup realGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float sla_gerpvgrw_(INTEGER N,INTEGER NCOLS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF);
/**
*> \brief \b SLA_LIN_BERR computes a component-wise relative backward error.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_LIN_BERR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_lin_berr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_lin_berr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_lin_berr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N, NZ, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AYB( N, NRHS ), BERR( NRHS )<br>
*       REAL               RES( N, NRHS )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLA_LIN_BERR computes componentwise relative backward error from<br>
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
*>          RES is REAL array, dimension (N,NRHS)<br>
*>     The residual matrix, i.e., the matrix R in the relative backward<br>
*>     error formula above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N, NRHS)<br>
*>     The denominator in the relative backward error formula above, i.e.,<br>
*>     the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B<br>
*>     are from iterative refinement (see sla_gerfsx_extended.f).<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sla_lin_berr_(INTEGER N,INTEGER NZ,INTEGER NRHS,float[] RES,float[] AYB,float[] BERR);
/**
*> \brief \b SLA_PORCOND estimates the Skeel condition number for a symmetric positive-definite matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_PORCOND + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_porcond.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_porcond.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_porcond.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, CMODE, C,<br>
*                                  INFO, WORK, IWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            N, LDA, LDAF, INFO, CMODE<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ),<br>
*      $                   C( * )<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLA_PORCOND Estimates the Skeel condition number of  op(A) * op2(C)<br>
*>    where op2 is determined by CMODE as follows<br>
*>    CMODE =  1    op2(C) = C<br>
*>    CMODE =  0    op2(C) = I<br>
*>    CMODE = -1    op2(C) = inv(C)<br>
*>    The Skeel condition number  cond(A) = norminf( |inv(A)||A| )<br>
*>    is computed by computing scaling factors R such that<br>
*>    diag(R)*A*op2(C) is row equilibrated and computing the standard<br>
*>    infinity-norm condition number.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The triangular factor U or L from the Cholesky factorization<br>
*>     A = U**T*U or A = L*L**T, as computed by SPOTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] CMODE<br>
*> \verbatim<br>
*>          CMODE is INTEGER<br>
*>     Determines op2(C) in the formula op(A) * op2(C) as follows:<br>
*>     CMODE =  1    op2(C) = C<br>
*>     CMODE =  0    op2(C) = I<br>
*>     CMODE = -1    op2(C) = inv(C)<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * op2(C).<br>
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
*>          WORK is REAL array, dimension (3*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N).<br>
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
*> \ingroup realPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float sla_porcond_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,INTEGER CMODE,float[] C,INTEGER INFO,float[] WORK,int[] IWORK);
/**
*> \brief \b SLA_PORFSX_EXTENDED improves the computed solution to a system of linear equations for symmetric or Hermitian positive-definite matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_PORFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_porfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_porfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_porfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_PORFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA,<br>
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
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
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
*> SLA_PORFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by SPORFSX to perform iterative refinement.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The triangular factor U or L from the Cholesky factorization<br>
*>     A = U**T*U or A = L*L**T, as computed by SPOTRF.<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          Y is REAL array, dimension (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by SPOTRS.<br>
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
*>     or vector Z. This is computed by SLA_LIN_BERR.<br>
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
*>          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>          RES is REAL array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace. This can be the same workspace passed for Y_TAIL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is REAL array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is REAL array, dimension (N)<br>
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
*>       < 0:  if INFO = -i, the ith argument to SPOTRS had an illegal<br>
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
*> \ingroup realPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sla_porfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
/**
*> \brief \b SLA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Hermitian positive-definite matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_PORPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_porpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_porpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_porpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, LDAF, WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*1        UPLO<br>
*       INTEGER            NCOLS, LDA, LDAF<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*> SLA_PORPVGRW computes the reciprocal pivot growth factor<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The triangular factor U or L from the Cholesky factorization<br>
*>     A = U**T*U or A = L*L**T, as computed by SPOTRF.<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup realPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float sla_porpvgrw_(char[] UPLO,INTEGER NCOLS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] WORK);
/**
*> \brief \b SLA_SYAMV computes a matrix-vector product using a symmetric indefinite matrix to calculate error bounds.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_SYAMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syamv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syamv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syamv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,<br>
*                             INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               ALPHA, BETA<br>
*       INTEGER            INCX, INCY, LDA, N, UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SLA_SYAMV  performs the matrix-vector operation<br>
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
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
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
*>          X is REAL array, dimension<br>
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
*> \ingroup realSYcomputational<br>
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
	public void sla_syamv_(INTEGER UPLO,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
/**
*> \brief \b SLA_SYRCOND estimates the Skeel condition number for a symmetric indefinite matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_SYRCOND + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syrcond.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syrcond.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syrcond.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, CMODE,<br>
*                                  C, INFO, WORK, IWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            N, LDA, LDAF, INFO, CMODE<br>
*       ..<br>
*       .. Array Arguments<br>
*       INTEGER            IWORK( * ), IPIV( * )<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C)<br>
*>    where op2 is determined by CMODE as follows<br>
*>    CMODE =  1    op2(C) = C<br>
*>    CMODE =  0    op2(C) = I<br>
*>    CMODE = -1    op2(C) = inv(C)<br>
*>    The Skeel condition number cond(A) = norminf( |inv(A)||A| )<br>
*>    is computed by computing scaling factors R such that<br>
*>    diag(R)*A*op2(C) is row equilibrated and computing the standard<br>
*>    infinity-norm condition number.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by SSYTRF.<br>
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
*>     as determined by SSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CMODE<br>
*> \verbatim<br>
*>          CMODE is INTEGER<br>
*>     Determines op2(C) in the formula op(A) * op2(C) as follows:<br>
*>     CMODE =  1    op2(C) = C<br>
*>     CMODE =  0    op2(C) = I<br>
*>     CMODE = -1    op2(C) = inv(C)<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (N)<br>
*>     The vector C in the formula op(A) * op2(C).<br>
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
*>          WORK is REAL array, dimension (3*N).<br>
*>     Workspace.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N).<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float sla_syrcond_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,INTEGER CMODE,float[] C,INTEGER INFO,float[] WORK,int[] IWORK);
/**
*> \brief \b SLA_SYRFSX_EXTENDED improves the computed solution to a system of linear equations for symmetric indefinite matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_SYRFSX_EXTENDED + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syrfsx_extended.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syrfsx_extended.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syrfsx_extended.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_SYRFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA,<br>
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
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
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
*> <br>
*> SLA_SYRFSX_EXTENDED improves the computed solution to a system of<br>
*> linear equations by performing extra-precise iterative refinement<br>
*> and provides error bounds and backward error estimates for the solution.<br>
*> This subroutine is called by SSYRFSX to perform iterative refinement.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by SSYTRF.<br>
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
*>     as determined by SSYTRF.<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          Y is REAL array, dimension (LDY,NRHS)<br>
*>     On entry, the solution matrix X, as computed by SSYTRS.<br>
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
*>     or vector Z. This is computed by SLA_LIN_BERR.<br>
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
*>          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>          RES is REAL array, dimension (N)<br>
*>     Workspace to hold the intermediate residual.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AYB<br>
*> \verbatim<br>
*>          AYB is REAL array, dimension (N)<br>
*>     Workspace. This can be the same workspace passed for Y_TAIL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY<br>
*> \verbatim<br>
*>          DY is REAL array, dimension (N)<br>
*>     Workspace to hold the intermediate solution.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y_TAIL<br>
*> \verbatim<br>
*>          Y_TAIL is REAL array, dimension (N)<br>
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
*>       < 0:  if INFO = -i, the ith argument to SLA_SYRFSX_EXTENDED had an illegal<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sla_syrfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
/**
*> \brief \b SLA_SYRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefinite matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_SYRPVGRW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syrpvgrw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syrpvgrw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syrpvgrw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV,<br>
*                                   WORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*1        UPLO<br>
*       INTEGER            N, INFO, LDA, LDAF<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> <br>
*> SLA_SYRPVGRW computes the reciprocal pivot growth factor<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The block diagonal matrix D and the multipliers used to<br>
*>     obtain the factor U or L as computed by SSYTRF.<br>
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
*>     as determined by SSYTRF.<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public float sla_syrpvgrw_(char[] UPLO,INTEGER N,INTEGER INFO,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] WORK);
/**
*> \brief \b SLA_WWADDW adds a vector into a doubled-single vector.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SLA_WWADDW + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_wwaddw.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_wwaddw.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_wwaddw.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SLA_WWADDW( N, X, Y, W )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               X( * ), Y( * ), W( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SLA_WWADDW adds a vector W into a doubled-single vector (X, Y).<br>
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
*>          X is REAL array, dimension (N)<br>
*>            The first part of the doubled-single accumulation vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array, dimension (N)<br>
*>            The second part of the doubled-single accumulation vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sla_wwaddw_(INTEGER N,float[] X,float[] Y,float[] W);

}