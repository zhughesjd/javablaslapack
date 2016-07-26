package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSO extends Library
{

	public static LapackSO instance = (LapackSO) Native.loadLibrary("liblapack",LapackSO.class);

/**
*> \brief \b SOPGTR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SOPGTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sopgtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sopgtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sopgtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDQ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SOPGTR generates a real orthogonal matrix Q which is defined as the<br>
*> product of n-1 elementary reflectors H(i) of order n, as returned by<br>
*> SSPTRD using packed storage:<br>
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
*>                 call to SSPTRD;<br>
*>          = 'L': Lower triangular packed storage used in previous<br>
*>                 call to SSPTRD.<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          The vectors which define the elementary reflectors, as<br>
*>          returned by SSPTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is REAL array, dimension (N-1)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SSPTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ,N)<br>
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
*>          WORK is REAL array, dimension (N-1)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sopgtr_(CHARACTER UPLO,INTEGER N,float[] AP,float[] TAU,float[] Q,INTEGER LDQ,float[] WORK,INTEGER INFO);
/**
*> \brief \b SOPMTR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SOPMTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sopmtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sopmtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sopmtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS, UPLO<br>
*       INTEGER            INFO, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SOPMTR overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix of order nq, with nq = m if<br>
*> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of<br>
*> nq-1 elementary reflectors, as returned by SSPTRD using packed<br>
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
*>                 call to SSPTRD;<br>
*>          = 'L': Lower triangular packed storage used in previous<br>
*>                 call to SSPTRD.<br>
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
*>          AP is REAL array, dimension<br>
*>                               (M*(M+1)/2) if SIDE = 'L'<br>
*>                               (N*(N+1)/2) if SIDE = 'R'<br>
*>          The vectors which define the elementary reflectors, as<br>
*>          returned by SSPTRD.  AP is modified by the routine but<br>
*>          restored on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TAU<br>
*> \verbatim<br>
*>          TAU is REAL array, dimension (M-1) if SIDE = 'L'<br>
*>                                     or (N-1) if SIDE = 'R'<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SSPTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sopmtr_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANS,INTEGER M,INTEGER N,float[] AP,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORBDB<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORBDB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12,<br>
*                          X21, LDX21, X22, LDX22, THETA, PHI, TAUP1,<br>
*                          TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIGNS, TRANS<br>
*       INTEGER            INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P,<br>
*      $                   Q<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               PHI( * ), THETA( * )<br>
*       REAL               TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ),<br>
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
*> SORBDB simultaneously bidiagonalizes the blocks of an M-by-M<br>
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
*> done in constant time using the TRANS and SIGNS options. See SORCSD<br>
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
*>          X11 is REAL array, dimension (LDX11,Q)<br>
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
*>          X12 is REAL array, dimension (LDX12,M-Q)<br>
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
*>          X21 is REAL array, dimension (LDX21,Q)<br>
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
*>          X22 is REAL array, dimension (LDX22,M-Q)<br>
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
*>          THETA is REAL array, dimension (Q)<br>
*>          The entries of the bidiagonal blocks B11, B12, B21, B22 can<br>
*>          be computed from the angles THETA and PHI. See Further<br>
*>          Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is REAL array, dimension (Q-1)<br>
*>          The entries of the bidiagonal blocks B11, B12, B21, B22 can<br>
*>          be computed from the angles THETA and PHI. See Further<br>
*>          Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is REAL array, dimension (P)<br>
*>          The scalar factors of the elementary reflectors that define<br>
*>          P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is REAL array, dimension (M-P)<br>
*>          The scalar factors of the elementary reflectors that define<br>
*>          P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is REAL array, dimension (Q)<br>
*>          The scalar factors of the elementary reflectors that define<br>
*>          Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ2<br>
*> \verbatim<br>
*>          TAUQ2 is REAL array, dimension (M-Q)<br>
*>          The scalar factors of the elementary reflectors that define<br>
*>          Q2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (LWORK)<br>
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
*> \ingroup realOTHERcomputational<br>
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
*>  [1] or SORCSD for details.<br>
*><br>
*>  P1, P2, Q1, and Q2 are represented as products of elementary<br>
*>  reflectors. See SORCSD for details on generating P1, P2, Q1, and Q2<br>
*>  using SORGQR and SORGLQ.<br>
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
	public void sorbdb_(CHARACTER TRANS,CHARACTER SIGNS,INTEGER M,INTEGER P,INTEGER Q,float[] X11,INTEGER LDX11,float[] X12,INTEGER LDX12,float[] X21,INTEGER LDX21,float[] X22,INTEGER LDX22,float[] THETA,float[] PHI,float[] TAUP1,float[] TAUP2,float[] TAUQ1,float[] TAUQ2,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORBDB1<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORBDB1 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb1.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb1.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb1.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,<br>
*                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               PHI(*), THETA(*)<br>
*       REAL               TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),<br>
*      $                   X11(LDX11,*), X21(LDX21,*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> SORBDB1 simultaneously bidiagonalizes the blocks of a tall and skinny<br>
*> matrix X with orthonomal columns:<br>
*><br>
*>                            [ B11 ]<br>
*>      [ X11 ]   [ P1 |    ] [  0  ]<br>
*>      [-----] = [---------] [-----] Q1**T .<br>
*>      [ X21 ]   [    | P2 ] [ B21 ]<br>
*>                            [  0  ]<br>
*><br>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P,<br>
*> M-P, or M-Q. Routines SORBDB2, SORBDB3, and SORBDB4 handle cases in<br>
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
*>          X11 is REAL array, dimension (LDX11,Q)<br>
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
*>          X21 is REAL array, dimension (LDX21,Q)<br>
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
*>          THETA is REAL array, dimension (Q)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is REAL array, dimension (Q-1)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is REAL array, dimension (P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is REAL array, dimension (M-P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is REAL array, dimension (Q)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (LWORK)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by<br>
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry<br>
*>  in each bidiagonal band is a product of a sine or cosine of a THETA<br>
*>  with a sine or cosine of a PHI. See [1] or SORCSD for details.<br>
*><br>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.<br>
*>  See SORCSD2BY1 for details on generating P1, P2, and Q1 using SORGQR<br>
*>  and SORGLQ.<br>
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
	public void sorbdb1_(INTEGER M,INTEGER P,INTEGER Q,float[] X11,INTEGER LDX11,float[] X21,INTEGER LDX21,float[] THETA,float[] PHI,float[] TAUP1,float[] TAUP2,float[] TAUQ1,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORBDB2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORBDB2 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb2.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb2.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb2.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,<br>
*                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               PHI(*), THETA(*)<br>
*       REAL               TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),<br>
*      $                   X11(LDX11,*), X21(LDX21,*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> SORBDB2 simultaneously bidiagonalizes the blocks of a tall and skinny<br>
*> matrix X with orthonomal columns:<br>
*><br>
*>                            [ B11 ]<br>
*>      [ X11 ]   [ P1 |    ] [  0  ]<br>
*>      [-----] = [---------] [-----] Q1**T .<br>
*>      [ X21 ]   [    | P2 ] [ B21 ]<br>
*>                            [  0  ]<br>
*><br>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,<br>
*> Q, or M-Q. Routines SORBDB1, SORBDB3, and SORBDB4 handle cases in<br>
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
*>          X11 is REAL array, dimension (LDX11,Q)<br>
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
*>          X21 is REAL array, dimension (LDX21,Q)<br>
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
*>          THETA is REAL array, dimension (Q)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is REAL array, dimension (Q-1)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is REAL array, dimension (P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is REAL array, dimension (M-P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is REAL array, dimension (Q)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (LWORK)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by<br>
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry<br>
*>  in each bidiagonal band is a product of a sine or cosine of a THETA<br>
*>  with a sine or cosine of a PHI. See [1] or SORCSD for details.<br>
*><br>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.<br>
*>  See SORCSD2BY1 for details on generating P1, P2, and Q1 using SORGQR<br>
*>  and SORGLQ.<br>
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
	public void sorbdb2_(INTEGER M,INTEGER P,INTEGER Q,float[] X11,INTEGER LDX11,float[] X21,INTEGER LDX21,float[] THETA,float[] PHI,float[] TAUP1,float[] TAUP2,float[] TAUQ1,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORBDB3<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORBDB3 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb3.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb3.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb3.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORBDB3( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,<br>
*                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               PHI(*), THETA(*)<br>
*       REAL               TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),<br>
*      $                   X11(LDX11,*), X21(LDX21,*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> SORBDB3 simultaneously bidiagonalizes the blocks of a tall and skinny<br>
*> matrix X with orthonomal columns:<br>
*><br>
*>                            [ B11 ]<br>
*>      [ X11 ]   [ P1 |    ] [  0  ]<br>
*>      [-----] = [---------] [-----] Q1**T .<br>
*>      [ X21 ]   [    | P2 ] [ B21 ]<br>
*>                            [  0  ]<br>
*><br>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-P must be no larger than P,<br>
*> Q, or M-Q. Routines SORBDB1, SORBDB2, and SORBDB4 handle cases in<br>
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
*>          X11 is REAL array, dimension (LDX11,Q)<br>
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
*>          X21 is REAL array, dimension (LDX21,Q)<br>
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
*>          THETA is REAL array, dimension (Q)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is REAL array, dimension (Q-1)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is REAL array, dimension (P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is REAL array, dimension (M-P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is REAL array, dimension (Q)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (LWORK)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by<br>
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry<br>
*>  in each bidiagonal band is a product of a sine or cosine of a THETA<br>
*>  with a sine or cosine of a PHI. See [1] or SORCSD for details.<br>
*><br>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.<br>
*>  See SORCSD2BY1 for details on generating P1, P2, and Q1 using SORGQR<br>
*>  and SORGLQ.<br>
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
	public void sorbdb3_(INTEGER M,INTEGER P,INTEGER Q,float[] X11,INTEGER LDX11,float[] X21,INTEGER LDX21,float[] THETA,float[] PHI,float[] TAUP1,float[] TAUP2,float[] TAUQ1,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORBDB4<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORBDB4 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb4.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb4.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb4.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,<br>
*                           TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK,<br>
*                           INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               PHI(*), THETA(*)<br>
*       REAL               PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*),<br>
*      $                   WORK(*), X11(LDX11,*), X21(LDX21,*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> SORBDB4 simultaneously bidiagonalizes the blocks of a tall and skinny<br>
*> matrix X with orthonomal columns:<br>
*><br>
*>                            [ B11 ]<br>
*>      [ X11 ]   [ P1 |    ] [  0  ]<br>
*>      [-----] = [---------] [-----] Q1**T .<br>
*>      [ X21 ]   [    | P2 ] [ B21 ]<br>
*>                            [  0  ]<br>
*><br>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P,<br>
*> M-P, or Q. Routines SORBDB1, SORBDB2, and SORBDB3 handle cases in<br>
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
*>          X11 is REAL array, dimension (LDX11,Q)<br>
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
*>          X21 is REAL array, dimension (LDX21,Q)<br>
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
*>          THETA is REAL array, dimension (Q)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHI<br>
*> \verbatim<br>
*>          PHI is REAL array, dimension (Q-1)<br>
*>           The entries of the bidiagonal blocks B11, B21 are defined by<br>
*>           THETA and PHI. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP1<br>
*> \verbatim<br>
*>          TAUP1 is REAL array, dimension (P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP2<br>
*> \verbatim<br>
*>          TAUP2 is REAL array, dimension (M-P)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           P2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ1<br>
*> \verbatim<br>
*>          TAUQ1 is REAL array, dimension (Q)<br>
*>           The scalar factors of the elementary reflectors that define<br>
*>           Q1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PHANTOM<br>
*> \verbatim<br>
*>          PHANTOM is REAL array, dimension (M)<br>
*>           The routine computes an M-by-1 column vector Y that is<br>
*>           orthogonal to the columns of [ X11; X21 ]. PHANTOM(1:P) and<br>
*>           PHANTOM(P+1:M) contain Householder vectors for Y(1:P) and<br>
*>           Y(P+1:M), respectively.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (LWORK)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by<br>
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry<br>
*>  in each bidiagonal band is a product of a sine or cosine of a THETA<br>
*>  with a sine or cosine of a PHI. See [1] or SORCSD for details.<br>
*><br>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.<br>
*>  See SORCSD2BY1 for details on generating P1, P2, and Q1 using SORGQR<br>
*>  and SORGLQ.<br>
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
	public void sorbdb4_(INTEGER M,INTEGER P,INTEGER Q,float[] X11,INTEGER LDX11,float[] X21,INTEGER LDX21,float[] THETA,float[] PHI,float[] TAUP1,float[] TAUP2,float[] TAUQ1,float[] PHANTOM,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORBDB5<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORBDB5 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb5.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb5.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb5.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,<br>
*                           LDQ2, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2,<br>
*      $                   N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> SORBDB5 orthogonalizes the column vector<br>
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
*>          X1 is REAL array, dimension (M1)<br>
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
*>          X2 is REAL array, dimension (M2)<br>
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
*>          Q1 is REAL array, dimension (LDQ1, N)<br>
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
*>          Q2 is REAL array, dimension (LDQ2, N)<br>
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
*>          WORK is REAL array, dimension (LWORK)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorbdb5_(INTEGER M1,INTEGER M2,INTEGER N,float[] X1,INTEGER INCX1,float[] X2,INTEGER INCX2,float[] Q1,INTEGER LDQ1,float[] Q2,INTEGER LDQ2,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORBDB6<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORBDB6 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb6.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb6.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb6.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,<br>
*                           LDQ2, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2,<br>
*      $                   N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)<br>
*       ..<br>
*  <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> SORBDB6 orthogonalizes the column vector<br>
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
*>          X1 is REAL array, dimension (M1)<br>
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
*>          X2 is REAL array, dimension (M2)<br>
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
*>          Q1 is REAL array, dimension (LDQ1, N)<br>
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
*>          Q2 is REAL array, dimension (LDQ2, N)<br>
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
*>          WORK is REAL array, dimension (LWORK)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorbdb6_(INTEGER M1,INTEGER M2,INTEGER N,float[] X1,INTEGER INCX1,float[] X2,INTEGER INCX2,float[] Q1,INTEGER LDQ1,float[] Q2,INTEGER LDQ2,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORCSD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORCSD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorcsd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorcsd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorcsd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       RECURSIVE SUBROUTINE SORCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS,<br>
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
*       REAL               THETA( * )<br>
*       REAL               U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),<br>
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
*> SORCSD computes the CS decomposition of an M-by-M partitioned<br>
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
*>          X11 is REAL array, dimension (LDX11,Q)<br>
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
*>          X12 is REAL array, dimension (LDX12,M-Q)<br>
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
*>          X21 is REAL array, dimension (LDX21,Q)<br>
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
*>          X22 is REAL array, dimension (LDX22,M-Q)<br>
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
*>          THETA is REAL array, dimension (R), in which R =<br>
*>          MIN(P,M-P,Q,M-Q).<br>
*>          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and<br>
*>          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U1<br>
*> \verbatim<br>
*>          U1 is REAL array, dimension (P)<br>
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
*>          U2 is REAL array, dimension (M-P)<br>
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
*>          V1T is REAL array, dimension (Q)<br>
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
*>          V2T is REAL array, dimension (M-Q)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*>          > 0:  SBBCSD did not converge. See the description of WORK<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorcsd_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,CHARACTER JOBV2T,CHARACTER TRANS,CHARACTER SIGNS,INTEGER M,INTEGER P,INTEGER Q,float[] X11,INTEGER LDX11,float[] X12,INTEGER LDX12,float[] X21,INTEGER LDX21,float[] X22,INTEGER LDX22,float[] THETA,float[] U1,INTEGER LDU1,float[] U2,INTEGER LDU2,float[] V1T,INTEGER LDV1T,float[] V2T,INTEGER LDV2T,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SORCSD2BY1<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORCSD2BY1 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorcsd2by1.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorcsd2by1.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorcsd2by1.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11,<br>
*                              X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T,<br>
*                              LDV1T, WORK, LWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBU1, JOBU2, JOBV1T<br>
*       INTEGER            INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21,<br>
*      $                   M, P, Q<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               THETA(*)<br>
*       REAL               U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*),<br>
*      $                   X11(LDX11,*), X21(LDX21,*)<br>
*       INTEGER            IWORK(*)<br>
*       ..<br>
*    <br>
* <br>
*> \par Purpose:<br>
*> =============<br>
*><br>
*>\verbatim<br>
*><br>
*> SORCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with<br>
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
*>          X11 is REAL array, dimension (LDX11,Q)<br>
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
*>          X21 is REAL array, dimension (LDX21,Q)<br>
*>          On entry, part of the orthogonal matrix whose CSD is desired.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX21<br>
*> \verbatim<br>
*>          LDX21 is INTEGER<br>
*>           The leading dimension of X21. LDX21 >= MAX(1,M-P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] THETA<br>
*> \verbatim<br>
*>          THETA is REAL array, dimension (R), in which R =<br>
*>          MIN(P,M-P,Q,M-Q).<br>
*>          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and<br>
*>          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U1<br>
*> \verbatim<br>
*>          U1 is REAL array, dimension (P)<br>
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
*>          U2 is REAL array, dimension (M-P)<br>
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
*>          V1T is REAL array, dimension (Q)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*>          > 0:  SBBCSD did not converge. See the description of WORK<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorcsd2by1_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,INTEGER M,INTEGER P,INTEGER Q,float[] X11,INTEGER LDX11,float[] X21,INTEGER LDX21,float[] THETA,float[] U1,INTEGER LDU1,float[] U2,INTEGER LDU2,float[] V1T,INTEGER LDV1T,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SORG2L generates all or part of the orthogonal matrix Q from a QL factorization determined by sgeqlf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORG2L + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorg2l.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorg2l.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorg2l.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORG2L( M, N, K, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, M, N<br>
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
*> SORG2L generates an m by n real matrix Q with orthonormal columns,<br>
*> which is defined as the last n columns of a product of k elementary<br>
*> reflectors of order m<br>
*><br>
*>       Q  =  H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by SGEQLF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the (n-k+i)-th column must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by SGEQLF in the last k columns of its array<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEQLF.<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorg2l_(INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORG2R + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorg2r.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorg2r.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorg2r.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORG2R( M, N, K, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, M, N<br>
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
*> SORG2R generates an m by n real matrix Q with orthonormal columns,<br>
*> which is defined as the first n columns of a product of k elementary<br>
*> reflectors of order m<br>
*><br>
*>       Q  =  H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by SGEQRF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the i-th column must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by SGEQRF in the first k columns of its array<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEQRF.<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorg2r_(INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORGBR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGBR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgbr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgbr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgbr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          VECT<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
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
*> SORGBR generates one of the real orthogonal matrices Q or P**T<br>
*> determined by SGEBRD when reducing a real matrix A to bidiagonal<br>
*> form: A = Q * B * P**T.  Q and P**T are defined as products of<br>
*> elementary reflectors H(i) or G(i) respectively.<br>
*><br>
*> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q<br>
*> is of order M:<br>
*> if m >= k, Q = H(1) H(2) . . . H(k) and SORGBR returns the first n<br>
*> columns of Q, where m >= n >= k;<br>
*> if m < k, Q = H(1) H(2) . . . H(m-1) and SORGBR returns Q as an<br>
*> M-by-M matrix.<br>
*><br>
*> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T<br>
*> is of order N:<br>
*> if k < n, P**T = G(k) . . . G(2) G(1) and SORGBR returns the first m<br>
*> rows of P**T, where n >= m >= k;<br>
*> if k >= n, P**T = G(n-1) . . . G(2) G(1) and SORGBR returns P**T as<br>
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
*>          required, as defined in the transformation applied by SGEBRD:<br>
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
*>          matrix reduced by SGEBRD.<br>
*>          If VECT = 'P', the number of rows in the original K-by-N<br>
*>          matrix reduced by SGEBRD.<br>
*>          K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the vectors which define the elementary reflectors,<br>
*>          as returned by SGEBRD.<br>
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
*>          TAU is REAL array, dimension<br>
*>                                (min(M,K)) if VECT = 'Q'<br>
*>                                (min(N,K)) if VECT = 'P'<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i) or G(i), which determines Q or P**T, as<br>
*>          returned by SGEBRD in its array argument TAUQ or TAUP.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorgbr_(CHARACTER VECT,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORGHR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGHR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorghr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorghr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorghr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N<br>
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
*> SORGHR generates a real orthogonal matrix Q which is defined as the<br>
*> product of IHI-ILO elementary reflectors of order N, as returned by<br>
*> SGEHRD:<br>
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
*>          of SGEHRD. Q is equal to the unit matrix except in the<br>
*>          submatrix Q(ilo+1:ihi,ilo+1:ihi).<br>
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the vectors which define the elementary reflectors,<br>
*>          as returned by SGEHRD.<br>
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
*>          TAU is REAL array, dimension (N-1)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEHRD.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorghr_(INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORGL2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGL2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgl2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgl2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgl2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGL2( M, N, K, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, M, N<br>
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
*> SORGL2 generates an m by n real matrix Q with orthonormal rows,<br>
*> which is defined as the first m rows of a product of k elementary<br>
*> reflectors of order n<br>
*><br>
*>       Q  =  H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by SGELQF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the i-th row must contain the vector which defines<br>
*>          the elementary reflector H(i), for i = 1,2,...,k, as returned<br>
*>          by SGELQF in the first k rows of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (M)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorgl2_(INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORGLQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGLQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorglq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorglq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorglq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
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
*> SORGLQ generates an M-by-N real matrix Q with orthonormal rows,<br>
*> which is defined as the first M rows of a product of K elementary<br>
*> reflectors of order N<br>
*><br>
*>       Q  =  H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by SGELQF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the i-th row must contain the vector which defines<br>
*>          the elementary reflector H(i), for i = 1,2,...,k, as returned<br>
*>          by SGELQF in the first k rows of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorglq_(INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORGQL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGQL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgql.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgql.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgql.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
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
*> SORGQL generates an M-by-N real matrix Q with orthonormal columns,<br>
*> which is defined as the last N columns of a product of K elementary<br>
*> reflectors of order M<br>
*><br>
*>       Q  =  H(k) . . . H(2) H(1)<br>
*><br>
*> as returned by SGEQLF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the (n-k+i)-th column must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by SGEQLF in the last k columns of its array<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEQLF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorgql_(INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORGQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
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
*> SORGQR generates an M-by-N real matrix Q with orthonormal columns,<br>
*> which is defined as the first N columns of a product of K elementary<br>
*> reflectors of order M<br>
*><br>
*>       Q  =  H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by SGEQRF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the i-th column must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by SGEQRF in the first k columns of its array<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEQRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorgqr_(INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORGR2 generates all or part of the orthogonal matrix Q from an RQ factorization determined by sgerqf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGR2( M, N, K, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, M, N<br>
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
*> SORGR2 generates an m by n real matrix Q with orthonormal rows,<br>
*> which is defined as the last m rows of a product of k elementary<br>
*> reflectors of order n<br>
*><br>
*>       Q  =  H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by SGERQF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the (m-k+i)-th row must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by SGERQF in the last k rows of its array argument<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGERQF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (M)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorgr2_(INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORGRQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGRQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgrq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgrq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgrq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, K, LDA, LWORK, M, N<br>
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
*> SORGRQ generates an M-by-N real matrix Q with orthonormal rows,<br>
*> which is defined as the last M rows of a product of K elementary<br>
*> reflectors of order N<br>
*><br>
*>       Q  =  H(1) H(2) . . . H(k)<br>
*><br>
*> as returned by SGERQF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the (m-k+i)-th row must contain the vector which<br>
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as<br>
*>          returned by SGERQF in the last k rows of its array argument<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGERQF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorgrq_(INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORGTR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORGTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
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
*> SORGTR generates a real orthogonal matrix Q which is defined as the<br>
*> product of n-1 elementary reflectors of order N, as returned by<br>
*> SSYTRD:<br>
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
*>                 from SSYTRD;<br>
*>          = 'L': Lower triangle of A contains elementary reflectors<br>
*>                 from SSYTRD.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the vectors which define the elementary reflectors,<br>
*>          as returned by SSYTRD.<br>
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
*>          TAU is REAL array, dimension (N-1)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SSYTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorgtr_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORM22 multiplies a general matrix by a banded orthogonal matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download SORM22 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorm22.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorm22.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorm22.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*     SUBROUTINE SORM22( SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC,<br>
*    $                   WORK, LWORK, INFO )<br>
*<br>
*     .. Scalar Arguments ..<br>
*     CHARACTER          SIDE, TRANS<br>
*     INTEGER            M, N, N1, N2, LDQ, LDC, LWORK, INFO<br>
*     ..<br>
*     .. Array Arguments ..<br>
*     REAL            Q( LDQ, * ), C( LDC, * ), WORK( * )<br>
*     ..<br>
*<br>
*> \par Purpose<br>
*  ============<br>
*><br>
*> \verbatim<br>
*><br>
*><br>
*>  SORM22 overwrites the general real M-by-N matrix C with<br>
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
*>          Q is REAL array, dimension<br>
*>                              (LDQ,M) if SIDE = 'L'<br>
*>                              (LDQ,N) if SIDE = 'R'<br>
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
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
	public void sorm22_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER N1,INTEGER N2,float[] Q,INTEGER LDQ,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORM2L multiplies a general matrix by the orthogonal matrix from a QL factorization determined by sgeqlf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORM2L + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorm2l.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorm2l.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorm2l.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORM2L overwrites the general real m by n matrix C with<br>
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
*> as returned by SGEQLF. Q is of order m if SIDE = 'L' and of order n<br>
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
*>          A is REAL array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          SGEQLF in the last k columns of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEQLF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorm2l_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORM2R multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORM2R + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorm2r.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorm2r.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorm2r.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORM2R overwrites the general real m by n matrix C with<br>
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
*> as returned by SGEQRF. Q is of order m if SIDE = 'L' and of order n<br>
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
*>          A is REAL array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          SGEQRF in the first k columns of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEQRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorm2r_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORMBR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMBR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormbr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormbr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormbr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,<br>
*                          LDC, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS, VECT<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> If VECT = 'Q', SORMBR overwrites the general real M-by-N matrix C<br>
*> with<br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> If VECT = 'P', SORMBR overwrites the general real M-by-N matrix C<br>
*> with<br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      P * C          C * P<br>
*> TRANS = 'T':      P**T * C       C * P**T<br>
*><br>
*> Here Q and P**T are the orthogonal matrices determined by SGEBRD when<br>
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
*>          matrix reduced by SGEBRD.<br>
*>          If VECT = 'P', the number of rows in the original<br>
*>          matrix reduced by SGEBRD.<br>
*>          K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension<br>
*>                                (LDA,min(nq,K)) if VECT = 'Q'<br>
*>                                (LDA,nq)        if VECT = 'P'<br>
*>          The vectors which define the elementary reflectors H(i) and<br>
*>          G(i), whose products determine the matrices Q and P, as<br>
*>          returned by SGEBRD.<br>
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
*>          TAU is REAL array, dimension (min(nq,K))<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i) or G(i) which determines Q or P, as returned<br>
*>          by SGEBRD in the array argument TAUQ or TAUP.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sormbr_(CHARACTER VECT,CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORMHR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMHR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormhr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormhr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormhr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,<br>
*                          LDC, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMHR overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix of order nq, with nq = m if<br>
*> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of<br>
*> IHI-ILO elementary reflectors, as returned by SGEHRD:<br>
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
*>          of SGEHRD. Q is equal to the unit matrix except in the<br>
*>          submatrix Q(ilo+1:ihi,ilo+1:ihi).<br>
*>          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and<br>
*>          ILO = 1 and IHI = 0, if M = 0;<br>
*>          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and<br>
*>          ILO = 1 and IHI = 0, if N = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension<br>
*>                               (LDA,M) if SIDE = 'L'<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The vectors which define the elementary reflectors, as<br>
*>          returned by SGEHRD.<br>
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
*>          TAU is REAL array, dimension<br>
*>                               (M-1) if SIDE = 'L'<br>
*>                               (N-1) if SIDE = 'R'<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEHRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sormhr_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORML2 multiplies a general matrix by the orthogonal matrix from a LQ factorization determined by sgelqf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORML2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorml2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorml2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorml2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORML2 overwrites the general real m by n matrix C with<br>
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
*> as returned by SGELQF. Q is of order m if SIDE = 'L' and of order n<br>
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
*>          A is REAL array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          SGELQF in the first k rows of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sorml2_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORMLQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMLQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormlq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormlq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormlq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMLQ overwrites the general real M-by-N matrix C with<br>
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
*> as returned by SGELQF. Q is of order M if SIDE = 'L' and of order N<br>
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
*>          A is REAL array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          SGELQF in the first k rows of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If SIDE = 'L', LWORK >= max(1,N);<br>
*>          if SIDE = 'R', LWORK >= max(1,M).<br>
*>          For good performance, LWORK should generally be larger. <br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sormlq_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORMQL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMQL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormql.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormql.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormql.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMQL overwrites the general real M-by-N matrix C with<br>
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
*> as returned by SGEQLF. Q is of order M if SIDE = 'L' and of order N<br>
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
*>          A is REAL array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          SGEQLF in the last k columns of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEQLF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sormql_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORMQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMQR overwrites the general real M-by-N matrix C with<br>
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
*> as returned by SGEQRF. Q is of order M if SIDE = 'L' and of order N<br>
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
*>          A is REAL array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          SGEQRF in the first k columns of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGEQRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sormqr_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORMR2 multiplies a general matrix by the orthogonal matrix from a RQ factorization determined by sgerqf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMR2 overwrites the general real m by n matrix C with<br>
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
*> as returned by SGERQF. Q is of order m if SIDE = 'L' and of order n<br>
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
*>          A is REAL array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          SGERQF in the last k rows of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGERQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sormr2_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORMR3 multiplies a general matrix by the orthogonal matrix from a RZ factorization determined by stzrzf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMR3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormr3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormr3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormr3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, L, LDA, LDC, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMR3 overwrites the general real m by n matrix C with<br>
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
*> as returned by STZRZF. Q is of order m if SIDE = 'L' and of order n<br>
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
*>          A is REAL array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          STZRZF in the last k rows of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by STZRZF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension<br>
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
	public void sormr3_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER L,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
/**
*> \brief \b SORMRQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMRQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormrq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormrq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormrq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMRQ overwrites the general real M-by-N matrix C with<br>
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
*> as returned by SGERQF. Q is of order M if SIDE = 'L' and of order N<br>
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
*>          A is REAL array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          SGERQF in the last k rows of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SGERQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sormrq_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORMRZ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMRZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormrz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormrz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormrz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS<br>
*       INTEGER            INFO, K, L, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMRZ overwrites the general real M-by-N matrix C with<br>
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
*> as returned by STZRZF. Q is of order M if SIDE = 'L' and of order N<br>
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
*>          A is REAL array, dimension<br>
*>                               (LDA,M) if SIDE = 'L',<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The i-th row must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          STZRZF in the last k rows of its array argument A.<br>
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
*>          TAU is REAL array, dimension (K)<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by STZRZF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
	public void sormrz_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER L,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SORMTR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SORMTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormtr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormtr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormtr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          SIDE, TRANS, UPLO<br>
*       INTEGER            INFO, LDA, LDC, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SORMTR overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q * C          C * Q<br>
*> TRANS = 'T':      Q**T * C       C * Q**T<br>
*><br>
*> where Q is a real orthogonal matrix of order nq, with nq = m if<br>
*> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of<br>
*> nq-1 elementary reflectors, as returned by SSYTRD:<br>
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
*>                 from SSYTRD;<br>
*>          = 'L': Lower triangle of A contains elementary reflectors<br>
*>                 from SSYTRD.<br>
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
*>          A is REAL array, dimension<br>
*>                               (LDA,M) if SIDE = 'L'<br>
*>                               (LDA,N) if SIDE = 'R'<br>
*>          The vectors which define the elementary reflectors, as<br>
*>          returned by SSYTRD.<br>
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
*>          TAU is REAL array, dimension<br>
*>                               (M-1) if SIDE = 'L'<br>
*>                               (N-1) if SIDE = 'R'<br>
*>          TAU(i) must contain the scalar factor of the elementary<br>
*>          reflector H(i), as returned by SSYTRD.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sormtr_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANS,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] C,INTEGER LDC,float[] WORK,INTEGER LWORK,INTEGER INFO);

}