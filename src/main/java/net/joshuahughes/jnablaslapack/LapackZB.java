package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZB extends Library
{

	public static LapackZB instance = (LapackZB) Native.loadLibrary("liblapack",LapackZB.class);

/**
*> \brief \b ZBBCSD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZBBCSD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zbbcsd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zbbcsd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zbbcsd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q,<br>
*                          THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T,<br>
*                          V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E,<br>
*                          B22D, B22E, RWORK, LRWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS<br>
*       INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LRWORK, M, P, Q<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   B11D( * ), B11E( * ), B12D( * ), B12E( * ),<br>
*      $                   B21D( * ), B21E( * ), B22D( * ), B22E( * ),<br>
*      $                   PHI( * ), THETA( * ), RWORK( * )<br>
*       COMPLEX*16         U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),<br>
*      $                   V2T( LDV2T, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZBBCSD computes the CS decomposition of a unitary matrix in<br>
*> bidiagonal-block form,<br>
*><br>
*><br>
*>     [ B11 | B12 0  0 ]<br>
*>     [  0  |  0 -I  0 ]<br>
*> X = [----------------]<br>
*>     [ B21 | B22 0  0 ]<br>
*>     [  0  |  0  0  I ]<br>
*><br>
*>                               [  C | -S  0  0 ]<br>
*>                   [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**H<br>
*>                 = [---------] [---------------] [---------]   .<br>
*>                   [    | U2 ] [  S |  C  0  0 ] [    | V2 ]<br>
*>                               [  0 |  0  0  I ]<br>
*><br>
*> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger<br>
*> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be<br>
*> transposed and/or permuted. This can be done in constant time using<br>
*> the TRANS and SIGNS options. See ZUNCSD for details.)<br>
*><br>
*> The bidiagonal matrices B11, B12, B21, and B22 are represented<br>
*> implicitly by angles THETA(1:Q) and PHI(1:Q-1).<br>
*><br>
*> The unitary matrices U1, U2, V1T, and V2T are input/output.<br>
*> The input matrices are pre- or post-multiplied by the appropriate<br>
*> singular vector matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBU1<br>
*> \verbatim<br>
*>          JOBU1 is CHARACTER<br>
*>          = 'Y':      U1 is updated;<br>
*>          otherwise:  U1 is not updated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBU2<br>
*> \verbatim<br>
*>          JOBU2 is CHARACTER<br>
*>          = 'Y':      U2 is updated;<br>
*>          otherwise:  U2 is not updated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV1T<br>
*> \verbatim<br>
*>          JOBV1T is CHARACTER<br>
*>          = 'Y':      V1T is updated;<br>
*>          otherwise:  V1T is not updated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV2T<br>
*> \verbatim<br>
*>          JOBV2T is CHARACTER<br>
*>          = 'Y':      V2T is updated;<br>
*>          otherwise:  V2T is not updated.<br>
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
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows and columns in X, the unitary matrix in<br>
*>          bidiagonal-block form.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of rows in the top-left block of X. 0 <= P <= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Q<br>
*> \verbatim<br>
*>          Q is INTEGER<br>
*>          The number of columns in the top-left block of X.<br>
*>          0 <= Q <= MIN(P,M-P,M-Q).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] THETA<br>
*> \verbatim<br>
*>          THETA is DOUBLE PRECISION array, dimension (Q)<br>
*>          On entry, the angles THETA(1),...,THETA(Q) that, along with<br>
*>          PHI(1), ...,PHI(Q-1), define the matrix in bidiagonal-block<br>
*>          form. On exit, the angles whose cosines and sines define the<br>
*>          diagonal blocks in the CS decomposition.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] PHI<br>
*> \verbatim<br>
*>          PHI is DOUBLE PRECISION array, dimension (Q-1)<br>
*>          The angles PHI(1),...,PHI(Q-1) that, along with THETA(1),...,<br>
*>          THETA(Q), define the matrix in bidiagonal-block form.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] U1<br>
*> \verbatim<br>
*>          U1 is COMPLEX*16 array, dimension (LDU1,P)<br>
*>          On entry, a P-by-P matrix. On exit, U1 is postmultiplied<br>
*>          by the left singular vector matrix common to [ B11 ; 0 ] and<br>
*>          [ B12 0 0 ; 0 -I 0 0 ].<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU1<br>
*> \verbatim<br>
*>          LDU1 is INTEGER<br>
*>          The leading dimension of the array U1, LDU1 >= MAX(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] U2<br>
*> \verbatim<br>
*>          U2 is COMPLEX*16 array, dimension (LDU2,M-P)<br>
*>          On entry, an (M-P)-by-(M-P) matrix. On exit, U2 is<br>
*>          postmultiplied by the left singular vector matrix common to<br>
*>          [ B21 ; 0 ] and [ B22 0 0 ; 0 0 I ].<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU2<br>
*> \verbatim<br>
*>          LDU2 is INTEGER<br>
*>          The leading dimension of the array U2, LDU2 >= MAX(1,M-P).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V1T<br>
*> \verbatim<br>
*>          V1T is COMPLEX*16 array, dimension (LDV1T,Q)<br>
*>          On entry, a Q-by-Q matrix. On exit, V1T is premultiplied<br>
*>          by the conjugate transpose of the right singular vector<br>
*>          matrix common to [ B11 ; 0 ] and [ B21 ; 0 ].<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV1T<br>
*> \verbatim<br>
*>          LDV1T is INTEGER<br>
*>          The leading dimension of the array V1T, LDV1T >= MAX(1,Q).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V2T<br>
*> \verbatim<br>
*>          V2T is COMPLEX*16 array, dimenison (LDV2T,M-Q)<br>
*>          On entry, an (M-Q)-by-(M-Q) matrix. On exit, V2T is<br>
*>          premultiplied by the conjugate transpose of the right<br>
*>          singular vector matrix common to [ B12 0 0 ; 0 -I 0 ] and<br>
*>          [ B22 0 0 ; 0 0 I ].<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV2T<br>
*> \verbatim<br>
*>          LDV2T is INTEGER<br>
*>          The leading dimension of the array V2T, LDV2T >= MAX(1,M-Q).<br>
*> \endverbatim<br>
*><br>
*> \param[out] B11D<br>
*> \verbatim<br>
*>          B11D is DOUBLE PRECISION array, dimension (Q)<br>
*>          When ZBBCSD converges, B11D contains the cosines of THETA(1),<br>
*>          ..., THETA(Q). If ZBBCSD fails to converge, then B11D<br>
*>          contains the diagonal of the partially reduced top-left<br>
*>          block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B11E<br>
*> \verbatim<br>
*>          B11E is DOUBLE PRECISION array, dimension (Q-1)<br>
*>          When ZBBCSD converges, B11E contains zeros. If ZBBCSD fails<br>
*>          to converge, then B11E contains the superdiagonal of the<br>
*>          partially reduced top-left block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B12D<br>
*> \verbatim<br>
*>          B12D is DOUBLE PRECISION array, dimension (Q)<br>
*>          When ZBBCSD converges, B12D contains the negative sines of<br>
*>          THETA(1), ..., THETA(Q). If ZBBCSD fails to converge, then<br>
*>          B12D contains the diagonal of the partially reduced top-right<br>
*>          block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B12E<br>
*> \verbatim<br>
*>          B12E is DOUBLE PRECISION array, dimension (Q-1)<br>
*>          When ZBBCSD converges, B12E contains zeros. If ZBBCSD fails<br>
*>          to converge, then B12E contains the subdiagonal of the<br>
*>          partially reduced top-right block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B21D<br>
*> \verbatim<br>
*>          B21D is DOUBLE PRECISION array, dimension (Q)<br>
*>          When ZBBCSD converges, B21D contains the negative sines of<br>
*>          THETA(1), ..., THETA(Q). If ZBBCSD fails to converge, then<br>
*>          B21D contains the diagonal of the partially reduced bottom-left<br>
*>          block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B21E<br>
*> \verbatim<br>
*>          B21E is DOUBLE PRECISION array, dimension (Q-1)<br>
*>          When ZBBCSD converges, B21E contains zeros. If ZBBCSD fails<br>
*>          to converge, then B21E contains the subdiagonal of the<br>
*>          partially reduced bottom-left block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B22D<br>
*> \verbatim<br>
*>          B22D is DOUBLE PRECISION array, dimension (Q)<br>
*>          When ZBBCSD converges, B22D contains the negative sines of<br>
*>          THETA(1), ..., THETA(Q). If ZBBCSD fails to converge, then<br>
*>          B22D contains the diagonal of the partially reduced bottom-right<br>
*>          block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B22E<br>
*> \verbatim<br>
*>          B22E is DOUBLE PRECISION array, dimension (Q-1)<br>
*>          When ZBBCSD converges, B22E contains zeros. If ZBBCSD fails<br>
*>          to converge, then B22E contains the subdiagonal of the<br>
*>          partially reduced bottom-right block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LRWORK<br>
*> \verbatim<br>
*>          LRWORK is INTEGER<br>
*>          The dimension of the array RWORK. LRWORK >= MAX(1,8*Q).<br>
*><br>
*>          If LRWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal size of the RWORK array,<br>
*>          returns this value as the first entry of the work array, and<br>
*>          no error message related to LRWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if ZBBCSD did not converge, INFO specifies the number<br>
*>                of nonzero entries in PHI, and B11D, B11E, etc.,<br>
*>                contain the partially reduced matrix.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  TOLMUL  DOUBLE PRECISION, default = MAX(10,MIN(100,EPS**(-1/8)))<br>
*>          TOLMUL controls the convergence criterion of the QR loop.<br>
*>          Angles THETA(i), PHI(i) are rounded to 0 or PI/2 when they<br>
*>          are within TOLMUL*EPS of either bound.<br>
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
*> \date June 2016<br>
*<br>
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zbbcsd_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,CHARACTER JOBV2T,CHARACTER TRANS,INTEGER M,INTEGER P,INTEGER Q,double[] THETA,double[] PHI,double[] U1,INTEGER LDU1,double[] U2,INTEGER LDU2,double[] V1T,INTEGER LDV1T,double[] V2T,INTEGER LDV2T,double[] B11D,double[] B11E,double[] B12D,double[] B12E,double[] B21D,double[] B21E,double[] B22D,double[] B22E,double[] RWORK,INTEGER LRWORK,INTEGER INFO);
/**
*> \brief \b ZBDSQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZBDSQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zbdsqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zbdsqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zbdsqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,<br>
*                          LDU, C, LDC, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   D( * ), E( * ), RWORK( * )<br>
*       COMPLEX*16         C( LDC, * ), U( LDU, * ), VT( LDVT, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZBDSQR computes the singular values and, optionally, the right and/or<br>
*> left singular vectors from the singular value decomposition (SVD) of<br>
*> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit<br>
*> zero-shift QR algorithm.  The SVD of B has the form<br>
*> <br>
*>    B = Q * S * P**H<br>
*> <br>
*> where S is the diagonal matrix of singular values, Q is an orthogonal<br>
*> matrix of left singular vectors, and P is an orthogonal matrix of<br>
*> right singular vectors.  If left singular vectors are requested, this<br>
*> subroutine actually returns U*Q instead of Q, and, if right singular<br>
*> vectors are requested, this subroutine returns P**H*VT instead of<br>
*> P**H, for given complex input matrices U and VT.  When U and VT are<br>
*> the unitary matrices that reduce a general matrix A to bidiagonal<br>
*> form: A = U*B*VT, as computed by ZGEBRD, then<br>
*> <br>
*>    A = (U*Q) * S * (P**H*VT)<br>
*> <br>
*> is the SVD of A.  Optionally, the subroutine may also compute Q**H*C<br>
*> for a given complex input matrix C.<br>
*><br>
*> See "Computing  Small Singular Values of Bidiagonal Matrices With<br>
*> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,<br>
*> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,<br>
*> no. 5, pp. 873-912, Sept 1990) and<br>
*> "Accurate singular values and differential qd algorithms," by<br>
*> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics<br>
*> Department, University of California at Berkeley, July 1992<br>
*> for a detailed description of the algorithm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  B is upper bidiagonal;<br>
*>          = 'L':  B is lower bidiagonal.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NCVT<br>
*> \verbatim<br>
*>          NCVT is INTEGER<br>
*>          The number of columns of the matrix VT. NCVT >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRU<br>
*> \verbatim<br>
*>          NRU is INTEGER<br>
*>          The number of rows of the matrix U. NRU >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NCC<br>
*> \verbatim<br>
*>          NCC is INTEGER<br>
*>          The number of columns of the matrix C. NCC >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, the n diagonal elements of the bidiagonal matrix B.<br>
*>          On exit, if INFO=0, the singular values of B in decreasing<br>
*>          order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          On entry, the N-1 offdiagonal elements of the bidiagonal<br>
*>          matrix B.<br>
*>          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E<br>
*>          will contain the diagonal and superdiagonal elements of a<br>
*>          bidiagonal matrix orthogonally equivalent to the one given<br>
*>          as input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VT<br>
*> \verbatim<br>
*>          VT is COMPLEX*16 array, dimension (LDVT, NCVT)<br>
*>          On entry, an N-by-NCVT matrix VT.<br>
*>          On exit, VT is overwritten by P**H * VT.<br>
*>          Not referenced if NCVT = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>          The leading dimension of the array VT.<br>
*>          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] U<br>
*> \verbatim<br>
*>          U is COMPLEX*16 array, dimension (LDU, N)<br>
*>          On entry, an NRU-by-N matrix U.<br>
*>          On exit, U is overwritten by U * Q.<br>
*>          Not referenced if NRU = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>          The leading dimension of the array U.  LDU >= max(1,NRU).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array, dimension (LDC, NCC)<br>
*>          On entry, an N-by-NCC matrix C.<br>
*>          On exit, C is overwritten by Q**H * C.<br>
*>          Not referenced if NCC = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C.<br>
*>          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  If INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  the algorithm did not converge; D and E contain the<br>
*>                elements of a bidiagonal matrix which is orthogonally<br>
*>                similar to the input matrix B;  if INFO = i, i<br>
*>                elements of E have not converged to zero.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))<br>
*>          TOLMUL controls the convergence criterion of the QR loop.<br>
*>          If it is positive, TOLMUL*EPS is the desired relative<br>
*>             precision in the computed singular values.<br>
*>          If it is negative, abs(TOLMUL*EPS*sigma_max) is the<br>
*>             desired absolute accuracy in the computed singular<br>
*>             values (corresponds to relative accuracy<br>
*>             abs(TOLMUL*EPS) in the largest singular value.<br>
*>          abs(TOLMUL) should be between 1 and 1/EPS, and preferably<br>
*>             between 10 (for fast convergence) and .1/EPS<br>
*>             (for there to be some accuracy in the results).<br>
*>          Default is to lose at either one eighth or 2 of the<br>
*>             available decimal digits in each computed singular value<br>
*>             (whichever is smaller).<br>
*><br>
*>  MAXITR  INTEGER, default = 6<br>
*>          MAXITR controls the maximum number of passes of the<br>
*>          algorithm through its inner loop. The algorithms stops<br>
*>          (and so fails to converge) if the number of passes<br>
*>          through the inner loop exceeds MAXITR*N**2.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zbdsqr_(CHARACTER UPLO,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,double[] D,double[] E,double[] VT,INTEGER LDVT,double[] U,INTEGER LDU,double[] C,INTEGER LDC,double[] RWORK,INTEGER INFO);

}