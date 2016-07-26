package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDO extends Library
{

	public static LapackDO instance = (LapackDO) Native.loadLibrary("liblapack",LapackDO.class);

/**
*> \brief \b DOPGTR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DOPGTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dopgtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dopgtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dopgtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDQ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DOPGTR generates a real orthogonal matrix Q which is defined as the<br>
*> product of n-1 elementary reflectors H(i) of order n, as returned by<br>
*> DSPTRD using packed storage:<br>
*><br>
*> if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),<br>
*><br>
*> if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U': Upper triangular packed storage used in previous<br>
*>                 call to DSPTRD;<br>
*>          = 'L': Lower triangular packed storage used in previous<br>
*>                 call to DSPTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix Q. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)<br>
*>          The vectors which define the elementary reflectors, as<br>
*>          returned by DSPTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DSPTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension (LDQ,N)<br>
*>          The N-by-N orthogonal matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q. LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (N-1)<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dopgtr_(CHARACTER UPLO,INTEGER N,double[] AP,double[] TAU,double[] Q,INTEGER LDQ,double[] WORK,INTEGER INFO);
/**
*> \brief \b DOPMTR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DOPMTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dopmtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dopmtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dopmtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS, UPLO<br>
*       INTEGER            INFO, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   AP( * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DOPMTR overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix of order nq, with nq = m if<br>
*> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of<br>
*> nq-1 elementary reflectors, as returned by DSPTRD using packed<br>
*> storage:<br>
*><br>
*> if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);<br>
*><br>
*> if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U': Upper triangular packed storage used in previous<br>
*>                 call to DSPTRD;<br>
*>          = 'L': Lower triangular packed storage used in previous<br>
*>                 call to DSPTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is DOUBLE PRECISION array, dimension<br>
*>                               (M*(M+1)/2) if SIDE = 'L'<br>
*>                               (N*(N+1)/2) if SIDE = 'R'<br>
*>          The vectors which define the elementary reflectors, as<br>
*>          returned by DSPTRD.  AP is modified by the routine but<br>
*>          restored on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (M-1) if SIDE = 'L'<br>
*>                                     or (N-1) if SIDE = 'R'<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DSPTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension<br>
*>                                   (N) if SIDE = 'L'<br>
*>                                   (M) if SIDE = 'R'<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dopmtr_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANS,INTEGER M,INTEGER N,double[] AP,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORBDB<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORBDB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12,<br>
*                          X21, LDX21, X22, LDX22, THETA, PHI, TAUP1,<br>
*                          TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIGNS, TRANS<br>
*       INTEGER            INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P,<br>
*      $                   Q<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   PHI( * ), THETA( * )<br>
*       DOUBLE PRECISION   TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ),<br>
*      $                   WORK( * ), X11( LDX11, * ), X12( LDX12, * ),<br>
*      $                   X21( LDX21, * ), X22( LDX22, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORBDB simultaneously bidiagonalizes the blocks of an M-by-M<br>
*> partitioned orthogonal matrix X:<br>
*><br>
*>                                 [ B11 | B12 0  0 ]<br>
*>     [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**T<br>
*> X = [-----------] = [---------] [----------------] [---------]   .<br>
*>     [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ]<br>
*>                                 [  0  |  0  0  I ]<br>
*><br>
*> X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is<br>
*> not the case, then X must be transposed and/or permuted. This can be<br>
*> done in constant time using the TRANS and SIGNS options. See DORCSD<br>
*> for details.)<br>
*><br>
*> The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by-<br>
*> (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are<br>
*> represented implicitly by Householder vectors.<br>
*><br>
*> B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented<br>
*> implicitly by angles THETA, PHI.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER<br>
*>          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major<br>
*>                      order;<br>
*>          otherwise:  X, U1, U2, V1T, and V2T are stored in column-<br>
*>                      major order.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIGNS<br>
*> \verbatim<br>
*>          SIGNS is CHARACTER<br>
*>          = 'O':      The lower-left block is made nonpositive (the<br>
*>                      "other" convention);<br>
*>          otherwise:  The upper-right block is made nonpositive (the<br>
*>                      "default" convention).<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows and columns in X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of rows in X11 and X12. 0 <= P <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is INTEGER<br>
*>          The number of columns in X11 and X21. 0 <= Q <=<br>
*>          MIN(P,M-P,M-Q).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X11<br>
*> \verbatim<br>
*>          X11 is DOUBLE PRECISION array, dimension (LDX11,Q)<br>
*>          On entry, the top-left block of the orthogonal matrix to be<br>
*>          reduced. On exit, the form depends on TRANS:<br>
*>          If TRANS = 'N', then<br>
*>             the columns of tril(X11) specify reflectors for P1,<br>
*>             the rows of triu(X11,1) specify reflectors for Q1;<br>
*>          else TRANS = 'T', and<br>
*>             the rows of triu(X11) specify reflectors for P1,<br>
*>             the columns of tril(X11,-1) specify reflectors for Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX11<br>
*> \verbatim<br>
*>          LDX11 is INTEGER<br>
*>          The leading dimension of X11. If TRANS = 'N', then LDX11 >=<br>
*>          P; else LDX11 >= Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X12<br>
*> \verbatim<br>
*>          X12 is DOUBLE PRECISION array, dimension (LDX12,M-Q)<br>
*>          On entry, the top-right block of the orthogonal matrix to<br>
*>          be reduced. On exit, the form depends on TRANS:<br>
*>          If TRANS = 'N', then<br>
*>             the rows of triu(X12) specify the first P reflectors for<br>
*>             Q2;<br>
*>          else TRANS = 'T', and<br>
*>             the columns of tril(X12) specify the first P reflectors<br>
*>             for Q2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX12<br>
*> \verbatim<br>
*>          LDX12 is INTEGER<br>
*>          The leading dimension of X12. If TRANS = 'N', then LDX12 >=<br>
*>          P; else LDX11 >= M-Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X21<br>
*> \verbatim<br>
*>          X21 is DOUBLE PRECISION array, dimension (LDX21,Q)<br>
*>          On entry, the bottom-left block of the orthogonal matrix to<br>
*>          be reduced. On exit, the form depends on TRANS:<br>
*>          If TRANS = 'N', then<br>
*>             the columns of tril(X21) specify reflectors for P2;<br>
*>          else TRANS = 'T', and<br>
*>             the rows of triu(X21) specify reflectors for P2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX21<br>
*> \verbatim<br>
*>          LDX21 is INTEGER<br>
*>          The leading dimension of X21. If TRANS = 'N', then LDX21 >=<br>
*>          M-P; else LDX21 >= Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X22<br>
*> \verbatim<br>
*>          X22 is DOUBLE PRECISION array, dimension (LDX22,M-Q)<br>
*>          On entry, the bottom-right block of the orthogonal matrix to<br>
*>          be reduced. On exit, the form depends on TRANS:<br>
*>          If TRANS = 'N', then<br>
*>             the rows of triu(X22(Q+1:M-P,P+1:M-Q)) specify the last<br>
*>             M-P-Q reflectors for Q2,<br>
*>          else TRANS = 'T', and<br>
*>             the columns of tril(X22(P+1:M-Q,Q+1:M-P)) specify the last<br>
*>             M-P-Q reflectors for P2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX22<br>
*> \verbatim<br>
*>          LDX22 is INTEGER<br>
*>          The leading dimension of X22. If TRANS = 'N', then LDX22 >=<br>
*>          M-P; else LDX22 >= M-Q.<br>
*> \endverbatim<br>
*><br>
*> \param[out] THETA<br>
*> \verbatim<br>
*>          THETA is DOUBLE PRECISION array, dimension (Q)<br>
*>          The entries of the bidiagonal blocks B11, B12, B21, B22 can<br>
*>          be computed from the angles THETA and PHI. See Further<br>
*>          Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is DOUBLE PRECISION array, dimension (Q-1)<br>
*>          The entries of the bidiagonal blocks B11, B12, B21, B22 can<br>
*>          be computed from the angles THETA and PHI. See Further<br>
*>          Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is DOUBLE PRECISION array, dimension (P)<br>
*>          The scalar factors of the elementary reflectors that define<br>
*>          P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is DOUBLE PRECISION array, dimension (M-P)<br>
*>          The scalar factors of the elementary reflectors that define<br>
*>          P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is DOUBLE PRECISION array, dimension (Q)<br>
*>          The scalar factors of the elementary reflectors that define<br>
*>          Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ2<br>
*> \verbatim<br>
*>          TAUQ2 is DOUBLE PRECISION array, dimension (M-Q)<br>
*>          The scalar factors of the elementary reflectors that define<br>
*>          Q2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= M-Q.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The bidiagonal blocks B11, B12, B21, and B22 are represented<br>
*>  implicitly by angles THETA(1), ..., THETA(Q) and PHI(1), ...,<br>
*>  PHI(Q-1). B11 and B21 are upper bidiagonal, while B21 and B22 are<br>
*>  lower bidiagonal. Every entry in each bidiagonal band is a product<br>
*>  of a sine or cosine of a THETA with a sine or cosine of a PHI. See<br>
*>  [1] or DORCSD for details.<br>
*><br>
*>  P1, P2, Q1, and Q2 are represented as products of elementary<br>
*>  reflectors. See DORCSD for details on generating P1, P2, Q1, and Q2<br>
*>  using DORGQR and DORGLQ.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.<br>
*>      Algorithms, 50(1):33-65, 2009.<br>
*><br>
*  =====================================================================<br>
*/
	public void dorbdb_(CHARACTER TRANS,CHARACTER SIGNS,INTEGER M,INTEGER P,INTEGER Q,double[] X11,INTEGER LDX11,double[] X12,INTEGER LDX12,double[] X21,INTEGER LDX21,double[] X22,INTEGER LDX22,double[] THETA,double[] PHI,double[] TAUP1,double[] TAUP2,double[] TAUQ1,double[] TAUQ2,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORBDB1<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORBDB1 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb1.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb1.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb1.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,<br>
*                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   PHI(*), THETA(*)<br>
*       DOUBLE PRECISION   TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),<br>
*      $                   X11(LDX11,*), X21(LDX21,*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> DORBDB1 simultaneously bidiagonalizes the blocks of a tall and skinny<br>
*> matrix X with orthonomal columns:<br>
*><br>
*>                            [ B11 ]<br>
*>      [ X11 ]   [ P1 |    ] [  0  ]<br>
*>      [-----] = [---------] [-----] Q1**T .<br>
*>      [ X21 ]   [    | P2 ] [ B21 ]<br>
*>                            [  0  ]<br>
*><br>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,<br>
*> M-P, or M-Q. Routines DORBDB2, DORBDB3, and DORBDB4 handle cases in<br>
*> which Q is not the minimum dimension.<br>
*><br>
*> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),<br>
*> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by<br>
*> Householder vectors.<br>
*><br>
*> B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by<br>
*> angles THETA, PHI.<br>
*><br>
*>\endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           The number of rows X11 plus the number of rows in X21.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>           The number of rows in X11. 0 <= P <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is INTEGER<br>
*>           The number of columns in X11 and X21. 0 <= Q <=<br>
*>           MIN(P,M-P,M-Q).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X11<br>
*> \verbatim<br>
*>          X11 is DOUBLE PRECISION array, dimension (LDX11,Q)<br>
*>           On entry, the top block of the matrix X to be reduced. On<br>
*>           exit, the columns of tril(X11) specify reflectors for P1 and<br>
*>           the rows of triu(X11,1) specify reflectors for Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX11<br>
*> \verbatim<br>
*>          LDX11 is INTEGER<br>
*>           The leading dimension of X11. LDX11 >= P.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X21<br>
*> \verbatim<br>
*>          X21 is DOUBLE PRECISION array, dimension (LDX21,Q)<br>
*>           On entry, the bottom block of the matrix X to be reduced. On<br>
*>           exit, the columns of tril(X21) specify reflectors for P2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX21<br>
*> \verbatim<br>
*>          LDX21 is INTEGER<br>
*>           The leading dimension of X21. LDX21 >= M-P.<br>
*> \endverbatim<br>
*><br>
*> \param[out] THETA<br>
*> \verbatim<br>
*>          THETA is DOUBLE PRECISION array, dimension (Q)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is DOUBLE PRECISION array, dimension (Q-1)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is DOUBLE PRECISION array, dimension (P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is DOUBLE PRECISION array, dimension (M-P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is DOUBLE PRECISION array, dimension (Q)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK. LWORK >= M-Q.<br>
*> <br>
*>           If LWORK = -1, then a workspace query is assumed; the routine<br>
*>           only calculates the optimal size of the WORK array, returns<br>
*>           this value as the first entry of the WORK array, and no error<br>
*>           message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           = 0:  successful exit.<br>
*>           < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*><br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date July 2012<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by<br>
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry<br>
*>  in each bidiagonal band is a product of a sine or cosine of a THETA<br>
*>  with a sine or cosine of a PHI. See [1] or DORCSD for details.<br>
*><br>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.<br>
*>  See DORCSD2BY1 for details on generating P1, P2, and Q1 using DORGQR<br>
*>  and DORGLQ.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.<br>
*>      Algorithms, 50(1):33-65, 2009.<br>
*><br>
*  =====================================================================<br>
*/
	public void dorbdb1_(INTEGER M,INTEGER P,INTEGER Q,double[] X11,INTEGER LDX11,double[] X21,INTEGER LDX21,double[] THETA,double[] PHI,double[] TAUP1,double[] TAUP2,double[] TAUQ1,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORBDB2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORBDB2 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb2.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb2.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb2.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,<br>
*                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   PHI(*), THETA(*)<br>
*       DOUBLE PRECISION   TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),<br>
*      $                   X11(LDX11,*), X21(LDX21,*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> DORBDB2 simultaneously bidiagonalizes the blocks of a tall and skinny<br>
*> matrix X with orthonomal columns:<br>
*><br>
*>                            [ B11 ]<br>
*>      [ X11 ]   [ P1 |    ] [  0  ]<br>
*>      [-----] = [---------] [-----] Q1**T .<br>
*>      [ X21 ]   [    | P2 ] [ B21 ]<br>
*>                            [  0  ]<br>
*><br>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,<br>
*> Q, or M-Q. Routines DORBDB1, DORBDB3, and DORBDB4 handle cases in<br>
*> which P is not the minimum dimension.<br>
*><br>
*> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),<br>
*> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by<br>
*> Householder vectors.<br>
*><br>
*> B11 and B12 are P-by-P bidiagonal matrices represented implicitly by<br>
*> angles THETA, PHI.<br>
*><br>
*>\endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           The number of rows X11 plus the number of rows in X21.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>           The number of rows in X11. 0 <= P <= min(M-P,Q,M-Q).<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is INTEGER<br>
*>           The number of columns in X11 and X21. 0 <= Q <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X11<br>
*> \verbatim<br>
*>          X11 is DOUBLE PRECISION array, dimension (LDX11,Q)<br>
*>           On entry, the top block of the matrix X to be reduced. On<br>
*>           exit, the columns of tril(X11) specify reflectors for P1 and<br>
*>           the rows of triu(X11,1) specify reflectors for Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX11<br>
*> \verbatim<br>
*>          LDX11 is INTEGER<br>
*>           The leading dimension of X11. LDX11 >= P.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X21<br>
*> \verbatim<br>
*>          X21 is DOUBLE PRECISION array, dimension (LDX21,Q)<br>
*>           On entry, the bottom block of the matrix X to be reduced. On<br>
*>           exit, the columns of tril(X21) specify reflectors for P2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX21<br>
*> \verbatim<br>
*>          LDX21 is INTEGER<br>
*>           The leading dimension of X21. LDX21 >= M-P.<br>
*> \endverbatim<br>
*><br>
*> \param[out] THETA<br>
*> \verbatim<br>
*>          THETA is DOUBLE PRECISION array, dimension (Q)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is DOUBLE PRECISION array, dimension (Q-1)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is DOUBLE PRECISION array, dimension (P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is DOUBLE PRECISION array, dimension (M-P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is DOUBLE PRECISION array, dimension (Q)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK. LWORK >= M-Q.<br>
*> <br>
*>           If LWORK = -1, then a workspace query is assumed; the routine<br>
*>           only calculates the optimal size of the WORK array, returns<br>
*>           this value as the first entry of the WORK array, and no error<br>
*>           message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           = 0:  successful exit.<br>
*>           < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*><br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date July 2012<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by<br>
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry<br>
*>  in each bidiagonal band is a product of a sine or cosine of a THETA<br>
*>  with a sine or cosine of a PHI. See [1] or DORCSD for details.<br>
*><br>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.<br>
*>  See DORCSD2BY1 for details on generating P1, P2, and Q1 using DORGQR<br>
*>  and DORGLQ.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.<br>
*>      Algorithms, 50(1):33-65, 2009.<br>
*><br>
*  =====================================================================<br>
*/
	public void dorbdb2_(INTEGER M,INTEGER P,INTEGER Q,double[] X11,INTEGER LDX11,double[] X21,INTEGER LDX21,double[] THETA,double[] PHI,double[] TAUP1,double[] TAUP2,double[] TAUQ1,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORBDB3<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORBDB3 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb3.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb3.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb3.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORBDB3( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,<br>
*                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   PHI(*), THETA(*)<br>
*       DOUBLE PRECISION   TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),<br>
*      $                   X11(LDX11,*), X21(LDX21,*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> DORBDB3 simultaneously bidiagonalizes the blocks of a tall and skinny<br>
*> matrix X with orthonomal columns:<br>
*><br>
*>                            [ B11 ]<br>
*>      [ X11 ]   [ P1 |    ] [  0  ]<br>
*>      [-----] = [---------] [-----] Q1**T .<br>
*>      [ X21 ]   [    | P2 ] [ B21 ]<br>
*>                            [  0  ]<br>
*><br>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,<br>
*> Q, or M-Q. Routines DORBDB1, DORBDB2, and DORBDB4 handle cases in<br>
*> which M-P is not the minimum dimension.<br>
*><br>
*> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),<br>
*> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by<br>
*> Householder vectors.<br>
*><br>
*> B11 and B12 are (M-P)-by-(M-P) bidiagonal matrices represented<br>
*> implicitly by angles THETA, PHI.<br>
*><br>
*>\endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           The number of rows X11 plus the number of rows in X21.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>           The number of rows in X11. 0 <= P <= M. M-P <= min(P,Q,M-Q).<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is INTEGER<br>
*>           The number of columns in X11 and X21. 0 <= Q <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X11<br>
*> \verbatim<br>
*>          X11 is DOUBLE PRECISION array, dimension (LDX11,Q)<br>
*>           On entry, the top block of the matrix X to be reduced. On<br>
*>           exit, the columns of tril(X11) specify reflectors for P1 and<br>
*>           the rows of triu(X11,1) specify reflectors for Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX11<br>
*> \verbatim<br>
*>          LDX11 is INTEGER<br>
*>           The leading dimension of X11. LDX11 >= P.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X21<br>
*> \verbatim<br>
*>          X21 is DOUBLE PRECISION array, dimension (LDX21,Q)<br>
*>           On entry, the bottom block of the matrix X to be reduced. On<br>
*>           exit, the columns of tril(X21) specify reflectors for P2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX21<br>
*> \verbatim<br>
*>          LDX21 is INTEGER<br>
*>           The leading dimension of X21. LDX21 >= M-P.<br>
*> \endverbatim<br>
*><br>
*> \param[out] THETA<br>
*> \verbatim<br>
*>          THETA is DOUBLE PRECISION array, dimension (Q)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is DOUBLE PRECISION array, dimension (Q-1)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is DOUBLE PRECISION array, dimension (P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is DOUBLE PRECISION array, dimension (M-P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is DOUBLE PRECISION array, dimension (Q)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK. LWORK >= M-Q.<br>
*> <br>
*>           If LWORK = -1, then a workspace query is assumed; the routine<br>
*>           only calculates the optimal size of the WORK array, returns<br>
*>           this value as the first entry of the WORK array, and no error<br>
*>           message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           = 0:  successful exit.<br>
*>           < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date July 2012<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by<br>
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry<br>
*>  in each bidiagonal band is a product of a sine or cosine of a THETA<br>
*>  with a sine or cosine of a PHI. See [1] or DORCSD for details.<br>
*><br>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.<br>
*>  See DORCSD2BY1 for details on generating P1, P2, and Q1 using DORGQR<br>
*>  and DORGLQ.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.<br>
*>      Algorithms, 50(1):33-65, 2009.<br>
*><br>
*  =====================================================================<br>
*/
	public void dorbdb3_(INTEGER M,INTEGER P,INTEGER Q,double[] X11,INTEGER LDX11,double[] X21,INTEGER LDX21,double[] THETA,double[] PHI,double[] TAUP1,double[] TAUP2,double[] TAUQ1,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORBDB4<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORBDB4 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb4.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb4.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb4.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,<br>
*                           TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK,<br>
*                           INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   PHI(*), THETA(*)<br>
*       DOUBLE PRECISION   PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*),<br>
*      $                   WORK(*), X11(LDX11,*), X21(LDX21,*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> DORBDB4 simultaneously bidiagonalizes the blocks of a tall and skinny<br>
*> matrix X with orthonomal columns:<br>
*><br>
*>                            [ B11 ]<br>
*>      [ X11 ]   [ P1 |    ] [  0  ]<br>
*>      [-----] = [---------] [-----] Q1**T .<br>
*>      [ X21 ]   [    | P2 ] [ B21 ]<br>
*>                            [  0  ]<br>
*><br>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,<br>
*> M-P, or Q. Routines DORBDB1, DORBDB2, and DORBDB3 handle cases in<br>
*> which M-Q is not the minimum dimension.<br>
*><br>
*> The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),<br>
*> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by<br>
*> Householder vectors.<br>
*><br>
*> B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented<br>
*> implicitly by angles THETA, PHI.<br>
*><br>
*>\endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           The number of rows X11 plus the number of rows in X21.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>           The number of rows in X11. 0 <= P <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is INTEGER<br>
*>           The number of columns in X11 and X21. 0 <= Q <= M and<br>
*>           M-Q <= min(P,M-P,Q).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X11<br>
*> \verbatim<br>
*>          X11 is DOUBLE PRECISION array, dimension (LDX11,Q)<br>
*>           On entry, the top block of the matrix X to be reduced. On<br>
*>           exit, the columns of tril(X11) specify reflectors for P1 and<br>
*>           the rows of triu(X11,1) specify reflectors for Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX11<br>
*> \verbatim<br>
*>          LDX11 is INTEGER<br>
*>           The leading dimension of X11. LDX11 >= P.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X21<br>
*> \verbatim<br>
*>          X21 is DOUBLE PRECISION array, dimension (LDX21,Q)<br>
*>           On entry, the bottom block of the matrix X to be reduced. On<br>
*>           exit, the columns of tril(X21) specify reflectors for P2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX21<br>
*> \verbatim<br>
*>          LDX21 is INTEGER<br>
*>           The leading dimension of X21. LDX21 >= M-P.<br>
*> \endverbatim<br>
*><br>
*> \param[out] THETA<br>
*> \verbatim<br>
*>          THETA is DOUBLE PRECISION array, dimension (Q)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is DOUBLE PRECISION array, dimension (Q-1)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is DOUBLE PRECISION array, dimension (P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is DOUBLE PRECISION array, dimension (M-P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is DOUBLE PRECISION array, dimension (Q)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHANTOM<br>
*> \verbatim<br>
*>          PHANTOM is DOUBLE PRECISION array, dimension (M)<br>
*>           The routine computes an M-by-1 column vector Y that is<br>
*>           orthogonal to the columns of [ X11; X21 ]. PHANTOM(1:P) and<br>
*>           PHANTOM(P+1:M) contain Householder vectors for Y(1:P) and<br>
*>           Y(P+1:M), respectively.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK. LWORK >= M-Q.<br>
*> <br>
*>           If LWORK = -1, then a workspace query is assumed; the routine<br>
*>           only calculates the optimal size of the WORK array, returns<br>
*>           this value as the first entry of the WORK array, and no error<br>
*>           message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           = 0:  successful exit.<br>
*>           < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date July 2012<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by<br>
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry<br>
*>  in each bidiagonal band is a product of a sine or cosine of a THETA<br>
*>  with a sine or cosine of a PHI. See [1] or DORCSD for details.<br>
*><br>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.<br>
*>  See DORCSD2BY1 for details on generating P1, P2, and Q1 using DORGQR<br>
*>  and DORGLQ.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.<br>
*>      Algorithms, 50(1):33-65, 2009.<br>
*><br>
*  =====================================================================<br>
*/
	public void dorbdb4_(INTEGER M,INTEGER P,INTEGER Q,double[] X11,INTEGER LDX11,double[] X21,INTEGER LDX21,double[] THETA,double[] PHI,double[] TAUP1,double[] TAUP2,double[] TAUQ1,double[] PHANTOM,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORBDB5<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORBDB5 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb5.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb5.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb5.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,<br>
*                           LDQ2, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2,<br>
*      $                   N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> DORBDB5 orthogonalizes the column vector<br>
*>      X = [ X1 ]<br>
*>          [ X2 ]<br>
*> with respect to the columns of<br>
*>      Q = [ Q1 ] .<br>
*>          [ Q2 ]<br>
*> The columns of Q must be orthonormal.<br>
*><br>
*> If the projection is zero according to Kahan's "twice is enough"<br>
*> criterion, then some other vector from the orthogonal complement<br>
*> is returned. This vector is chosen in an arbitrary but deterministic<br>
*> way.<br>
*><br>
*>\endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M1<br>
*> \verbatim<br>
*>          M1 is INTEGER<br>
*>           The dimension of X1 and the number of rows in Q1. 0 <= M1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M2<br>
*> \verbatim<br>
*>          M2 is INTEGER<br>
*>           The dimension of X2 and the number of rows in Q2. 0 <= M2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           The number of columns in Q1 and Q2. 0 <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X1<br>
*> \verbatim<br>
*>          X1 is DOUBLE PRECISION array, dimension (M1)<br>
*>           On entry, the top part of the vector to be orthogonalized.<br>
*>           On exit, the top part of the projected vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX1<br>
*> \verbatim<br>
*>          INCX1 is INTEGER<br>
*>           Increment for entries of X1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X2<br>
*> \verbatim<br>
*>          X2 is DOUBLE PRECISION array, dimension (M2)<br>
*>           On entry, the bottom part of the vector to be<br>
*>           orthogonalized. On exit, the bottom part of the projected<br>
*>           vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX2<br>
*> \verbatim<br>
*>          INCX2 is INTEGER<br>
*>           Increment for entries of X2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q1<br>
*> \verbatim<br>
*>          Q1 is DOUBLE PRECISION array, dimension (LDQ1, N)<br>
*>           The top part of the orthonormal basis matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ1<br>
*> \verbatim<br>
*>          LDQ1 is INTEGER<br>
*>           The leading dimension of Q1. LDQ1 >= M1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q2<br>
*> \verbatim<br>
*>          Q2 is DOUBLE PRECISION array, dimension (LDQ2, N)<br>
*>           The bottom part of the orthonormal basis matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ2<br>
*> \verbatim<br>
*>          LDQ2 is INTEGER<br>
*>           The leading dimension of Q2. LDQ2 >= M2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK. LWORK >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           = 0:  successful exit.<br>
*>           < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date July 2012<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorbdb5_(INTEGER M1,INTEGER M2,INTEGER N,double[] X1,INTEGER INCX1,double[] X2,INTEGER INCX2,double[] Q1,INTEGER LDQ1,double[] Q2,INTEGER LDQ2,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORBDB6<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORBDB6 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb6.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb6.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb6.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,<br>
*                           LDQ2, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2,<br>
*      $                   N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> DORBDB6 orthogonalizes the column vector<br>
*>      X = [ X1 ]<br>
*>          [ X2 ]<br>
*> with respect to the columns of<br>
*>      Q = [ Q1 ] .<br>
*>          [ Q2 ]<br>
*> The columns of Q must be orthonormal.<br>
*><br>
*> If the projection is zero according to Kahan's "twice is enough"<br>
*> criterion, then the zero vector is returned.<br>
*><br>
*>\endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M1<br>
*> \verbatim<br>
*>          M1 is INTEGER<br>
*>           The dimension of X1 and the number of rows in Q1. 0 <= M1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M2<br>
*> \verbatim<br>
*>          M2 is INTEGER<br>
*>           The dimension of X2 and the number of rows in Q2. 0 <= M2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           The number of columns in Q1 and Q2. 0 <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X1<br>
*> \verbatim<br>
*>          X1 is DOUBLE PRECISION array, dimension (M1)<br>
*>           On entry, the top part of the vector to be orthogonalized.<br>
*>           On exit, the top part of the projected vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX1<br>
*> \verbatim<br>
*>          INCX1 is INTEGER<br>
*>           Increment for entries of X1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X2<br>
*> \verbatim<br>
*>          X2 is DOUBLE PRECISION array, dimension (M2)<br>
*>           On entry, the bottom part of the vector to be<br>
*>           orthogonalized. On exit, the bottom part of the projected<br>
*>           vector.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX2<br>
*> \verbatim<br>
*>          INCX2 is INTEGER<br>
*>           Increment for entries of X2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q1<br>
*> \verbatim<br>
*>          Q1 is DOUBLE PRECISION array, dimension (LDQ1, N)<br>
*>           The top part of the orthonormal basis matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ1<br>
*> \verbatim<br>
*>          LDQ1 is INTEGER<br>
*>           The leading dimension of Q1. LDQ1 >= M1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q2<br>
*> \verbatim<br>
*>          Q2 is DOUBLE PRECISION array, dimension (LDQ2, N)<br>
*>           The bottom part of the orthonormal basis matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ2<br>
*> \verbatim<br>
*>          LDQ2 is INTEGER<br>
*>           The leading dimension of Q2. LDQ2 >= M2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK. LWORK >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           = 0:  successful exit.<br>
*>           < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date July 2012<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorbdb6_(INTEGER M1,INTEGER M2,INTEGER N,double[] X1,INTEGER INCX1,double[] X2,INTEGER INCX2,double[] Q1,INTEGER LDQ1,double[] Q2,INTEGER LDQ2,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORCSD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORCSD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorcsd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorcsd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorcsd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       RECURSIVE SUBROUTINE DORCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS,<br>
*                                    SIGNS, M, P, Q, X11, LDX11, X12,<br>
*                                    LDX12, X21, LDX21, X22, LDX22, THETA,<br>
*                                    U1, LDU1, U2, LDU2, V1T, LDV1T, V2T,<br>
*                                    LDV2T, WORK, LWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS<br>
*       INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12,<br>
*      $                   LDX21, LDX22, LWORK, M, P, Q<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   THETA( * )<br>
*       DOUBLE PRECISION   U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),<br>
*      $                   V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ),<br>
*      $                   X12( LDX12, * ), X21( LDX21, * ), X22( LDX22,<br>
*      $                   * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORCSD computes the CS decomposition of an M-by-M partitioned<br>
*> orthogonal matrix X:<br>
*><br>
*>                                 [  I  0  0 |  0  0  0 ]<br>
*>                                 [  0  C  0 |  0 -S  0 ]<br>
*>     [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**T<br>
*> X = [-----------] = [---------] [---------------------] [---------]   .<br>
*>     [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]<br>
*>                                 [  0  S  0 |  0  C  0 ]<br>
*>                                 [  0  0  I |  0  0  0 ]<br>
*><br>
*> X11 is P-by-Q. The orthogonal matrices U1, U2, V1, and V2 are P-by-P,<br>
*> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are<br>
*> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in<br>
*> which R = MIN(P,M-P,Q,M-Q).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBU1<br>
*> \verbatim<br>
*>          JOBU1 is CHARACTER<br>
*>          = 'Y':      U1 is computed;<br>
*>          otherwise:  U1 is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBU2<br>
*> \verbatim<br>
*>          JOBU2 is CHARACTER<br>
*>          = 'Y':      U2 is computed;<br>
*>          otherwise:  U2 is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV1T<br>
*> \verbatim<br>
*>          JOBV1T is CHARACTER<br>
*>          = 'Y':      V1T is computed;<br>
*>          otherwise:  V1T is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV2T<br>
*> \verbatim<br>
*>          JOBV2T is CHARACTER<br>
*>          = 'Y':      V2T is computed;<br>
*>          otherwise:  V2T is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER<br>
*>          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major<br>
*>                      order;<br>
*>          otherwise:  X, U1, U2, V1T, and V2T are stored in column-<br>
*>                      major order.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIGNS<br>
*> \verbatim<br>
*>          SIGNS is CHARACTER<br>
*>          = 'O':      The lower-left block is made nonpositive (the<br>
*>                      "other" convention);<br>
*>          otherwise:  The upper-right block is made nonpositive (the<br>
*>                      "default" convention).<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows and columns in X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of rows in X11 and X12. 0 <= P <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is INTEGER<br>
*>          The number of columns in X11 and X21. 0 <= Q <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X11<br>
*> \verbatim<br>
*>          X11 is DOUBLE PRECISION array, dimension (LDX11,Q)<br>
*>          On entry, part of the orthogonal matrix whose CSD is desired.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX11<br>
*> \verbatim<br>
*>          LDX11 is INTEGER<br>
*>          The leading dimension of X11. LDX11 >= MAX(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X12<br>
*> \verbatim<br>
*>          X12 is DOUBLE PRECISION array, dimension (LDX12,M-Q)<br>
*>          On entry, part of the orthogonal matrix whose CSD is desired.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX12<br>
*> \verbatim<br>
*>          LDX12 is INTEGER<br>
*>          The leading dimension of X12. LDX12 >= MAX(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X21<br>
*> \verbatim<br>
*>          X21 is DOUBLE PRECISION array, dimension (LDX21,Q)<br>
*>          On entry, part of the orthogonal matrix whose CSD is desired.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX21<br>
*> \verbatim<br>
*>          LDX21 is INTEGER<br>
*>          The leading dimension of X11. LDX21 >= MAX(1,M-P).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X22<br>
*> \verbatim<br>
*>          X22 is DOUBLE PRECISION array, dimension (LDX22,M-Q)<br>
*>          On entry, part of the orthogonal matrix whose CSD is desired.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX22<br>
*> \verbatim<br>
*>          LDX22 is INTEGER<br>
*>          The leading dimension of X11. LDX22 >= MAX(1,M-P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] THETA<br>
*> \verbatim<br>
*>          THETA is DOUBLE PRECISION array, dimension (R), in which R =<br>
*>          MIN(P,M-P,Q,M-Q).<br>
*>          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and<br>
*>          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U1<br>
*> \verbatim<br>
*>          U1 is DOUBLE PRECISION array, dimension (P)<br>
*>          If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU1<br>
*> \verbatim<br>
*>          LDU1 is INTEGER<br>
*>          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=<br>
*>          MAX(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U2<br>
*> \verbatim<br>
*>          U2 is DOUBLE PRECISION array, dimension (M-P)<br>
*>          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal<br>
*>          matrix U2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU2<br>
*> \verbatim<br>
*>          LDU2 is INTEGER<br>
*>          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=<br>
*>          MAX(1,M-P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] V1T<br>
*> \verbatim<br>
*>          V1T is DOUBLE PRECISION array, dimension (Q)<br>
*>          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal<br>
*>          matrix V1**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV1T<br>
*> \verbatim<br>
*>          LDV1T is INTEGER<br>
*>          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=<br>
*>          MAX(1,Q).<br>
*> \endverbatim<br>
*><br>
*> \param[out] V2T<br>
*> \verbatim<br>
*>          V2T is DOUBLE PRECISION array, dimension (M-Q)<br>
*>          If JOBV2T = 'Y', V2T contains the (M-Q)-by-(M-Q) orthogonal<br>
*>          matrix V2**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV2T<br>
*> \verbatim<br>
*>          LDV2T is INTEGER<br>
*>          The leading dimension of V2T. If JOBV2T = 'Y', LDV2T >=<br>
*>          MAX(1,M-Q).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*>          If INFO > 0 on exit, WORK(2:R) contains the values PHI(1),<br>
*>          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),<br>
*>          define the matrix in intermediate bidiagonal-block form<br>
*>          remaining after nonconvergence. INFO specifies the number<br>
*>          of nonzero PHI's.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the work array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (M-MIN(P, M-P, Q, M-Q))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  DBBCSD did not converge. See the description of WORK<br>
*>                above for details.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.<br>
*>      Algorithms, 50(1):33-65, 2009.<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorcsd_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,CHARACTER JOBV2T,CHARACTER TRANS,CHARACTER SIGNS,INTEGER M,INTEGER P,INTEGER Q,double[] X11,INTEGER LDX11,double[] X12,INTEGER LDX12,double[] X21,INTEGER LDX21,double[] X22,INTEGER LDX22,double[] THETA,double[] U1,INTEGER LDU1,double[] U2,INTEGER LDU2,double[] V1T,INTEGER LDV1T,double[] V2T,INTEGER LDV2T,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DORCSD2BY1<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORCSD2BY1 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorcsd2by1.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorcsd2by1.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorcsd2by1.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11,<br>
*                              X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T,<br>
*                              LDV1T, WORK, LWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBU1, JOBU2, JOBV1T<br>
*       INTEGER            INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21,<br>
*      $                   M, P, Q<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   THETA(*)<br>
*       DOUBLE PRECISION   U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*),<br>
*      $                   X11(LDX11,*), X21(LDX21,*)<br>
*       INTEGER            IWORK(*)<br>
*       ..<br>
*    <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*> Purpose:<br>
*> ========<br>
*><br>
*> DORCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with<br>
*> orthonormal columns that has been partitioned into a 2-by-1 block<br>
*> structure:<br>
*><br>
*>                                [  I  0  0 ]<br>
*>                                [  0  C  0 ]<br>
*>          [ X11 ]   [ U1 |    ] [  0  0  0 ]<br>
*>      X = [-----] = [---------] [----------] V1**T .<br>
*>          [ X21 ]   [    | U2 ] [  0  0  0 ]<br>
*>                                [  0  S  0 ]<br>
*>                                [  0  0  I ]<br>
*> <br>
*> X11 is P-by-Q. The orthogonal matrices U1, U2, and V1 are P-by-P,<br>
*> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R<br>
*> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which<br>
*> R = MIN(P,M-P,Q,M-Q).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBU1<br>
*> \verbatim<br>
*>          JOBU1 is CHARACTER<br>
*>          = 'Y':      U1 is computed;<br>
*>          otherwise:  U1 is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBU2<br>
*> \verbatim<br>
*>          JOBU2 is CHARACTER<br>
*>          = 'Y':      U2 is computed;<br>
*>          otherwise:  U2 is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV1T<br>
*> \verbatim<br>
*>          JOBV1T is CHARACTER<br>
*>          = 'Y':      V1T is computed;<br>
*>          otherwise:  V1T is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows in X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of rows in X11. 0 <= P <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is INTEGER<br>
*>          The number of columns in X11 and X21. 0 <= Q <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X11<br>
*> \verbatim<br>
*>          X11 is DOUBLE PRECISION array, dimension (LDX11,Q)<br>
*>          On entry, part of the orthogonal matrix whose CSD is desired.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX11<br>
*> \verbatim<br>
*>          LDX11 is INTEGER<br>
*>          The leading dimension of X11. LDX11 >= MAX(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X21<br>
*> \verbatim<br>
*>          X21 is DOUBLE PRECISION array, dimension (LDX21,Q)<br>
*>          On entry, part of the orthogonal matrix whose CSD is desired.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX21<br>
*> \verbatim<br>
*>          LDX21 is INTEGER<br>
*>          The leading dimension of X21. LDX21 >= MAX(1,M-P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] THETA<br>
*> \verbatim<br>
*>          THETA is DOUBLE PRECISION array, dimension (R), in which R =<br>
*>          MIN(P,M-P,Q,M-Q).<br>
*>          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and<br>
*>          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U1<br>
*> \verbatim<br>
*>          U1 is DOUBLE PRECISION array, dimension (P)<br>
*>          If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU1<br>
*> \verbatim<br>
*>          LDU1 is INTEGER<br>
*>          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=<br>
*>          MAX(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U2<br>
*> \verbatim<br>
*>          U2 is DOUBLE PRECISION array, dimension (M-P)<br>
*>          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal<br>
*>          matrix U2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU2<br>
*> \verbatim<br>
*>          LDU2 is INTEGER<br>
*>          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=<br>
*>          MAX(1,M-P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] V1T<br>
*> \verbatim<br>
*>          V1T is DOUBLE PRECISION array, dimension (Q)<br>
*>          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal<br>
*>          matrix V1**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV1T<br>
*> \verbatim<br>
*>          LDV1T is INTEGER<br>
*>          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=<br>
*>          MAX(1,Q).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*>          If INFO > 0 on exit, WORK(2:R) contains the values PHI(1),<br>
*>          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),<br>
*>          define the matrix in intermediate bidiagonal-block form<br>
*>          remaining after nonconvergence. INFO specifies the number<br>
*>          of nonzero PHI's.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the work array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  DBBCSD did not converge. See the description of WORK<br>
*>                above for details.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.<br>
*>      Algorithms, 50(1):33-65, 2009.<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date July 2012<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorcsd2by1_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,INTEGER M,INTEGER P,INTEGER Q,double[] X11,INTEGER LDX11,double[] X21,INTEGER LDX21,double[] THETA,double[] U1,INTEGER LDU1,double[] U2,INTEGER LDU2,double[] V1T,INTEGER LDV1T,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DORG2L generates all or part of the orthogonal matrix Q from a QL factorization determined by sgeqlf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORG2L + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2l.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2l.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2l.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORG2L generates an m by n real matrix Q with orthonormal columns,<br>
*> which is defined as the last n columns of a product of k elementary<br>
*> reflectors of order m<br>
*><br>
*>       Q  =  H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by DGEQLF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q. M >= N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines the<br>
*>          matrix Q. N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the (n-k+i)-th column must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by DGEQLF in the last k columns of its array<br>
*>          argument A.<br>
*>          On exit, the m by n matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The first dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEQLF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument has an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorg2l_(INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORG2R + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2r.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2r.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORG2R generates an m by n real matrix Q with orthonormal columns,<br>
*> which is defined as the first n columns of a product of k elementary<br>
*> reflectors of order m<br>
*><br>
*>       Q  =  H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DGEQRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q. M >= N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines the<br>
*>          matrix Q. N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the i-th column must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by DGEQRF in the first k columns of its array<br>
*>          argument A.<br>
*>          On exit, the m-by-n matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The first dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEQRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument has an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorg2r_(INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORGBR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGBR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgbr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgbr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgbr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          VECT<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGBR generates one of the real orthogonal matrices Q or P**T<br>
*> determined by DGEBRD when reducing a real matrix A to bidiagonal<br>
*> form: A = Q * B * P**T.  Q and P**T are defined as products of<br>
*> elementary reflectors H(i) or G(i) respectively.<br>
*><br>
*> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q<br>
*> is of order M:<br>
*> if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n<br>
*> columns of Q, where m >= n >= k;<br>
*> if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an<br>
*> M-by-M matrix.<br>
*><br>
*> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T<br>
*> is of order N:<br>
*> if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m<br>
*> rows of P**T, where n >= m >= k;<br>
*> if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as<br>
*> an N-by-N matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] VECT<br>
*> \verbatim<br>
*>          VECT is CHARACTER*1<br>
*>          Specifies whether the matrix Q or the matrix P**T is<br>
*>          required, as defined in the transformation applied by DGEBRD:<br>
*>          = 'Q':  generate Q;<br>
*>          = 'P':  generate P**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q or P**T to be returned.<br>
*>          M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q or P**T to be returned.<br>
*>          N >= 0.<br>
*>          If VECT = 'Q', M >= N >= min(M,K);<br>
*>          if VECT = 'P', N >= M >= min(N,K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          If VECT = 'Q', the number of columns in the original M-by-K<br>
*>          matrix reduced by DGEBRD.<br>
*>          If VECT = 'P', the number of rows in the original K-by-N<br>
*>          matrix reduced by DGEBRD.<br>
*>          K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the vectors which define the elementary reflectors,<br>
*>          as returned by DGEBRD.<br>
*>          On exit, the M-by-N matrix Q or P**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension<br>
*>                                (min(M,K)) if VECT = 'Q'<br>
*>                                (min(N,K)) if VECT = 'P'<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i) or G(i), which determines Q or P**T, as<br>
*>          returned by DGEBRD in its array argument TAUQ or TAUP.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,min(M,N)).<br>
*>          For optimum performance LWORK >= min(M,N)*NB, where NB<br>
*>          is the optimal blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date April 2012<br>
*<br>
*> \ingroup doubleGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorgbr_(CHARACTER VECT,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORGHR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGHR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorghr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorghr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorghr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGHR generates a real orthogonal matrix Q which is defined as the<br>
*> product of IHI-ILO elementary reflectors of order N, as returned by<br>
*> DGEHRD:<br>
*><br>
*> Q = H(ilo) H(ilo+1) . . . H(ihi-1).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix Q. N >= 0.<br>
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
*><br>
*>          ILO and IHI must have the same values as in the previous call<br>
*>          of DGEHRD. Q is equal to the unit matrix except in the<br>
*>          submatrix Q(ilo+1:ihi,ilo+1:ihi).<br>
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the vectors which define the elementary reflectors,<br>
*>          as returned by DGEHRD.<br>
*>          On exit, the N-by-N orthogonal matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEHRD.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= IHI-ILO.<br>
*>          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is<br>
*>          the optimal blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorghr_(INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORGL2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGL2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgl2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgl2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgl2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGL2 generates an m by n real matrix Q with orthonormal rows,<br>
*> which is defined as the first m rows of a product of k elementary<br>
*> reflectors of order n<br>
*><br>
*>       Q  =  H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by DGELQF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q. N >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines the<br>
*>          matrix Q. M >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the i-th row must contain the vector which defines<br>
*>          the elementary reflector H(i), for i = 1,2,...,k, as returned<br>
*>          by DGELQF in the first k rows of its array argument A.<br>
*>          On exit, the m-by-n matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The first dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (M)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument has an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorgl2_(INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORGLQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGLQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorglq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorglq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorglq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGLQ generates an M-by-N real matrix Q with orthonormal rows,<br>
*> which is defined as the first M rows of a product of K elementary<br>
*> reflectors of order N<br>
*><br>
*>       Q  =  H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by DGELQF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q. N >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines the<br>
*>          matrix Q. M >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the i-th row must contain the vector which defines<br>
*>          the elementary reflector H(i), for i = 1,2,...,k, as returned<br>
*>          by DGELQF in the first k rows of its array argument A.<br>
*>          On exit, the M-by-N matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The first dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,M).<br>
*>          For optimum performance LWORK >= M*NB, where NB is<br>
*>          the optimal blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument has an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorglq_(INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORGQL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGQL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgql.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgql.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgql.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGQL generates an M-by-N real matrix Q with orthonormal columns,<br>
*> which is defined as the last N columns of a product of K elementary<br>
*> reflectors of order M<br>
*><br>
*>       Q  =  H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by DGEQLF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q. M >= N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines the<br>
*>          matrix Q. N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the (n-k+i)-th column must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by DGEQLF in the last k columns of its array<br>
*>          argument A.<br>
*>          On exit, the M-by-N matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The first dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEQLF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,N).<br>
*>          For optimum performance LWORK >= N*NB, where NB is the<br>
*>          optimal blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument has an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorgql_(INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORGQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGQR generates an M-by-N real matrix Q with orthonormal columns,<br>
*> which is defined as the first N columns of a product of K elementary<br>
*> reflectors of order M<br>
*><br>
*>       Q  =  H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DGEQRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q. M >= N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines the<br>
*>          matrix Q. N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the i-th column must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by DGEQRF in the first k columns of its array<br>
*>          argument A.<br>
*>          On exit, the M-by-N matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The first dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEQRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,N).<br>
*>          For optimum performance LWORK >= N*NB, where NB is the<br>
*>          optimal blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument has an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorgqr_(INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORGR2 generates all or part of the orthogonal matrix Q from an RQ factorization determined by sgerqf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGR2( M, N, K, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGR2 generates an m by n real matrix Q with orthonormal rows,<br>
*> which is defined as the last m rows of a product of k elementary<br>
*> reflectors of order n<br>
*><br>
*>       Q  =  H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DGERQF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q. N >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines the<br>
*>          matrix Q. M >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the (m-k+i)-th row must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by DGERQF in the last k rows of its array argument<br>
*>          A.<br>
*>          On exit, the m by n matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The first dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGERQF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (M)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument has an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorgr2_(INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORGRQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGRQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgrq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgrq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgrq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGRQ generates an M-by-N real matrix Q with orthonormal rows,<br>
*> which is defined as the last M rows of a product of K elementary<br>
*> reflectors of order N<br>
*><br>
*>       Q  =  H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DGERQF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix Q. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix Q. N >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines the<br>
*>          matrix Q. M >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the (m-k+i)-th row must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by DGERQF in the last k rows of its array argument<br>
*>          A.<br>
*>          On exit, the M-by-N matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The first dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGERQF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,M).<br>
*>          For optimum performance LWORK >= M*NB, where NB is the<br>
*>          optimal blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument has an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorgrq_(INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORGTR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORGTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORGTR generates a real orthogonal matrix Q which is defined as the<br>
*> product of n-1 elementary reflectors of order N, as returned by<br>
*> DSYTRD:<br>
*><br>
*> if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),<br>
*><br>
*> if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U': Upper triangle of A contains elementary reflectors<br>
*>                 from DSYTRD;<br>
*>          = 'L': Lower triangle of A contains elementary reflectors<br>
*>                 from DSYTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix Q. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the vectors which define the elementary reflectors,<br>
*>          as returned by DSYTRD.<br>
*>          On exit, the N-by-N orthogonal matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DSYTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,N-1).<br>
*>          For optimum performance LWORK >= (N-1)*NB, where NB is<br>
*>          the optimal blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorgtr_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORM22 multiplies a general matrix by a banded orthogonal matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download DORM22 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm22.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm22.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm22.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*     SUBROUTINE DORM22( SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC,<br>
*    $                   WORK, LWORK, INFO )<br>
*<br>
*     .. Scalar Arguments ..<br>
*     CHARACTER          SIDE, TRANS<br>
*     INTEGER            M, N, N1, N2, LDQ, LDC, LWORK, INFO<br>
*     ..<br>
*     .. Array Arguments ..<br>
*     DOUBLE PRECISION   Q( LDQ, * ), C( LDC, * ), WORK( * )<br>
*     ..<br>
*<br>
*> \par Purpose<br>
*  ============<br>
*><br>
*> \verbatim<br>
*><br>
*><br>
*>  DORM22 overwrites the general real M-by-N matrix C with<br>
*><br>
*>                  SIDE = 'L'     SIDE = 'R'<br>
*>  TRANS = 'N':      Q * C          C * Q<br>
*>  TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*>  where Q is a real orthogonal matrix of order NQ, with NQ = M if<br>
*>  SIDE = 'L' and NQ = N if SIDE = 'R'.<br>
*>  The orthogonal matrix Q processes a 2-by-2 block structure<br>
*><br>
*>         [  Q11  Q12  ]<br>
*>     Q = [            ]<br>
*>         [  Q21  Q22  ],<br>
*><br>
*>  where Q12 is an N1-by-N1 lower triangular matrix and Q21 is an<br>
*>  N2-by-N2 upper triangular matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments<br>
*  =========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  apply Q (No transpose);<br>
*>          = 'C':  apply Q**T (Conjugate transpose).<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N1<br>
*> \param[in] N2<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*>          N2 is INTEGER<br>
*>          The dimension of Q12 and Q21, respectively. N1, N2 >= 0.<br>
*>          The following requirement must be satisfied:<br>
*>          N1 + N2 = M if SIDE = 'L' and N1 + N2 = N if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension<br>
*>                                       (LDQ,M) if SIDE = 'L'<br>
*>                                       (LDQ,N) if SIDE = 'R'<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.<br>
*>          LDQ >= max(1,M) if SIDE = 'L'; LDQ >= max(1,N) if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For optimum performance LWORK >= M*N.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date January 2015<br>
*<br>
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorm22_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER N1,INTEGER N2,double[] Q,INTEGER LDQ,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORM2L multiplies a general matrix by the orthogonal matrix from a QL factorization determined by sgeqlf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORM2L + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2l.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm2l.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm2l.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORM2L overwrites the general real m by n matrix C with<br>
*><br>
*>       Q * C  if SIDE = 'L' and TRANS = 'N', or<br>
*><br>
*>       Q**T * C  if SIDE = 'L' and TRANS = 'T', or<br>
*><br>
*>       C * Q  if SIDE = 'R' and TRANS = 'N', or<br>
*><br>
*>       C * Q**T if SIDE = 'R' and TRANS = 'T',<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by DGEQLF. Q is of order m if SIDE = 'L' and of order n<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left<br>
*>          = 'R': apply Q or Q**T from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply Q  (No transpose)<br>
*>          = 'T': apply Q**T (Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DGEQLF in the last k columns of its array argument A.<br>
*>          A is modified by the routine but restored on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*>          If SIDE = 'L', LDA >= max(1,M);<br>
*>          if SIDE = 'R', LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEQLF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the m by n matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension<br>
*>                                   (N) if SIDE = 'L',<br>
*>                                   (M) if SIDE = 'R'<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorm2l_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORM2R multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORM2R + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2r.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm2r.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm2r.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORM2R overwrites the general real m by n matrix C with<br>
*><br>
*>       Q * C  if SIDE = 'L' and TRANS = 'N', or<br>
*><br>
*>       Q**T* C  if SIDE = 'L' and TRANS = 'T', or<br>
*><br>
*>       C * Q  if SIDE = 'R' and TRANS = 'N', or<br>
*><br>
*>       C * Q**T if SIDE = 'R' and TRANS = 'T',<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left<br>
*>          = 'R': apply Q or Q**T from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply Q  (No transpose)<br>
*>          = 'T': apply Q**T (Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DGEQRF in the first k columns of its array argument A.<br>
*>          A is modified by the routine but restored on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*>          If SIDE = 'L', LDA >= max(1,M);<br>
*>          if SIDE = 'R', LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEQRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the m by n matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension<br>
*>                                   (N) if SIDE = 'L',<br>
*>                                   (M) if SIDE = 'R'<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorm2r_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORMBR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMBR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormbr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormbr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormbr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,<br>
*                          LDC, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS, VECT<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C<br>
*> with<br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C<br>
*> with<br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      P * C          C * P<br>
*> TRANS = 'T':      P**T * C       C * P**T<br>
*><br>
*> Here Q and P**T are the orthogonal matrices determined by DGEBRD when<br>
*> reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and<br>
*> P**T are defined as products of elementary reflectors H(i) and G(i)<br>
*> respectively.<br>
*><br>
*> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the<br>
*> order of the orthogonal matrix Q or P**T that is applied.<br>
*><br>
*> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:<br>
*> if nq >= k, Q = H(1) H(2) . . . H(k);<br>
*> if nq < k, Q = H(1) H(2) . . . H(nq-1).<br>
*><br>
*> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:<br>
*> if k < nq, P = G(1) G(2) . . . G(k);<br>
*> if k >= nq, P = G(1) G(2) . . . G(nq-1).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] VECT<br>
*> \verbatim<br>
*>          VECT is CHARACTER*1<br>
*>          = 'Q': apply Q or Q**T;<br>
*>          = 'P': apply P or P**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q, Q**T, P or P**T from the Left;<br>
*>          = 'R': apply Q, Q**T, P or P**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q  or P;<br>
*>          = 'T':  Transpose, apply Q**T or P**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          If VECT = 'Q', the number of columns in the original<br>
*>          matrix reduced by DGEBRD.<br>
*>          If VECT = 'P', the number of rows in the original<br>
*>          matrix reduced by DGEBRD.<br>
*>          K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                                (LDA,min(nq,K)) if VECT = 'Q'<br>
*>                                (LDA,nq)        if VECT = 'P'<br>
*>          The vectors which define the elementary reflectors H(i) and<br>
*>          G(i), whose products determine the matrices Q and P, as<br>
*>          returned by DGEBRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*>          If VECT = 'Q', LDA >= max(1,nq);<br>
*>          if VECT = 'P', LDA >= max(1,min(nq,K)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (min(nq,K))<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i) or G(i) which determines Q or P, as returned<br>
*>          by DGEBRD in the array argument TAUQ or TAUP.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q<br>
*>          or P*C or P**T*C or C*P or C*P**T.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For optimum performance LWORK >= N*NB if SIDE = 'L', and<br>
*>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal<br>
*>          blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dormbr_(CHARACTER VECT,CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORMHR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMHR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormhr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormhr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormhr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,<br>
*                          LDC, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMHR overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix of order nq, with nq = m if<br>
*> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of<br>
*> IHI-ILO elementary reflectors, as returned by DGEHRD:<br>
*><br>
*> Q = H(ilo) H(ilo+1) . . . H(ihi-1).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
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
*><br>
*>          ILO and IHI must have the same values as in the previous call<br>
*>          of DGEHRD. Q is equal to the unit matrix except in the<br>
*>          submatrix Q(ilo+1:ihi,ilo+1:ihi).<br>
*>          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and<br>
*>          ILO = 1 and IHI = 0, if M = 0;<br>
*>          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and<br>
*>          ILO = 1 and IHI = 0, if N = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                               (LDA,M) if SIDE = 'L'<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The vectors which define the elementary reflectors, as<br>
*>          returned by DGEHRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*>          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension<br>
*>                               (M-1) if SIDE = 'L'<br>
*>                               (N-1) if SIDE = 'R'<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEHRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For optimum performance LWORK >= N*NB if SIDE = 'L', and<br>
*>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal<br>
*>          blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dormhr_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORML2 multiplies a general matrix by the orthogonal matrix from a LQ factorization determined by sgelqf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORML2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorml2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorml2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorml2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORML2 overwrites the general real m by n matrix C with<br>
*><br>
*>       Q * C  if SIDE = 'L' and TRANS = 'N', or<br>
*><br>
*>       Q**T* C  if SIDE = 'L' and TRANS = 'T', or<br>
*><br>
*>       C * Q  if SIDE = 'R' and TRANS = 'N', or<br>
*><br>
*>       C * Q**T if SIDE = 'R' and TRANS = 'T',<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by DGELQF. Q is of order m if SIDE = 'L' and of order n<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left<br>
*>          = 'R': apply Q or Q**T from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply Q  (No transpose)<br>
*>          = 'T': apply Q**T (Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DGELQF in the first k rows of its array argument A.<br>
*>          A is modified by the routine but restored on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the m by n matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension<br>
*>                                   (N) if SIDE = 'L',<br>
*>                                   (M) if SIDE = 'R'<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dorml2_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORMLQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMLQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormlq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormlq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormlq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMLQ overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DGELQF in the first k rows of its array argument A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For good performance, LWORK should generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dormlq_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORMQL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMQL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormql.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormql.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormql.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMQL overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by DGEQLF. Q is of order M if SIDE = 'L' and of order N<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DGEQLF in the last k columns of its array argument A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*>          If SIDE = 'L', LDA >= max(1,M);<br>
*>          if SIDE = 'R', LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEQLF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For good performance, LWORK should generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dormql_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORMQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMQR overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DGEQRF in the first k columns of its array argument A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*>          If SIDE = 'L', LDA >= max(1,M);<br>
*>          if SIDE = 'R', LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGEQRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For good performance, LWORK should generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dormqr_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORMR2 multiplies a general matrix by the orthogonal matrix from a RQ factorization determined by sgerqf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMR2 overwrites the general real m by n matrix C with<br>
*><br>
*>       Q * C  if SIDE = 'L' and TRANS = 'N', or<br>
*><br>
*>       Q**T* C  if SIDE = 'L' and TRANS = 'T', or<br>
*><br>
*>       C * Q  if SIDE = 'R' and TRANS = 'N', or<br>
*><br>
*>       C * Q**T if SIDE = 'R' and TRANS = 'T',<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DGERQF. Q is of order m if SIDE = 'L' and of order n<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left<br>
*>          = 'R': apply Q or Q**T from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply Q  (No transpose)<br>
*>          = 'T': apply Q' (Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DGERQF in the last k rows of its array argument A.<br>
*>          A is modified by the routine but restored on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGERQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the m by n matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension<br>
*>                                   (N) if SIDE = 'L',<br>
*>                                   (M) if SIDE = 'R'<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dormr2_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORMR3 multiplies a general matrix by the orthogonal matrix from a RZ factorization determined by stzrzf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMR3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormr3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormr3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormr3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, L, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMR3 overwrites the general real m by n matrix C with<br>
*><br>
*>       Q * C  if SIDE = 'L' and TRANS = 'N', or<br>
*><br>
*>       Q**T* C  if SIDE = 'L' and TRANS = 'C', or<br>
*><br>
*>       C * Q  if SIDE = 'R' and TRANS = 'N', or<br>
*><br>
*>       C * Q**T if SIDE = 'R' and TRANS = 'C',<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DTZRZF. Q is of order m if SIDE = 'L' and of order n<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left<br>
*>          = 'R': apply Q or Q**T from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply Q  (No transpose)<br>
*>          = 'T': apply Q**T (Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The number of columns of the matrix A containing<br>
*>          the meaningful part of the Householder reflectors.<br>
*>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DTZRZF in the last k rows of its array argument A.<br>
*>          A is modified by the routine but restored on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DTZRZF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the m-by-n matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension<br>
*>                                   (N) if SIDE = 'L',<br>
*>                                   (M) if SIDE = 'R'<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERcomputational<br>
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
	public void dormr3_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER L,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DORMRQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMRQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormrq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormrq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormrq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMRQ overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DGERQF. Q is of order M if SIDE = 'L' and of order N<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DGERQF in the last k rows of its array argument A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DGERQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For good performance, LWORK should generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dormrq_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORMRZ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMRZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormrz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormrz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormrz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, L, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMRZ overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of k<br>
*> elementary reflectors<br>
*><br>
*>       Q = H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by DTZRZF. Q is of order M if SIDE = 'L' and of order N<br>
*> if SIDE = 'R'.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*>          If SIDE = 'L', M >= K >= 0;<br>
*>          if SIDE = 'R', N >= K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The number of columns of the matrix A containing<br>
*>          the meaningful part of the Householder reflectors.<br>
*>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          DTZRZF in the last k rows of its array argument A.<br>
*>          A is modified by the routine but restored on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,K).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DTZRZF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For good performance, LWORK should generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
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
	public void dormrz_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER L,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DORMTR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DORMTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS, UPLO<br>
*       INTEGER            INFO, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DORMTR overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix of order nq, with nq = m if<br>
*> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of<br>
*> nq-1 elementary reflectors, as returned by DSYTRD:<br>
*><br>
*> if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);<br>
*><br>
*> if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**T from the Left;<br>
*>          = 'R': apply Q or Q**T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U': Upper triangle of A contains elementary reflectors<br>
*>                 from DSYTRD;<br>
*>          = 'L': Lower triangle of A contains elementary reflectors<br>
*>                 from DSYTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension<br>
*>                               (LDA,M) if SIDE = 'L'<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The vectors which define the elementary reflectors, as<br>
*>          returned by DSYTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.<br>
*>          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension<br>
*>                               (M-1) if SIDE = 'L'<br>
*>                               (N-1) if SIDE = 'R'<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by DSYTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.<br>
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For optimum performance LWORK >= N*NB if SIDE = 'L', and<br>
*>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal<br>
*>          blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dormtr_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANS,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] C,INTEGER LDC,double[] WORK,INTEGER LWORK,INTEGER INFO);

}