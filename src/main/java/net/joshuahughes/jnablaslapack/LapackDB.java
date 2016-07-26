package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDB extends Library
{

	public static LapackDB instance = (LapackDB) Native.loadLibrary("liblapack",LapackDB.class);

/**
*> \brief \b DBBCSD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DBBCSD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbbcsd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbbcsd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbbcsd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q,<br>
*                          THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T,<br>
*                          V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E,<br>
*                          B22D, B22E, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS<br>
*       INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   B11D( * ), B11E( * ), B12D( * ), B12E( * ),<br>
*      $                   B21D( * ), B21E( * ), B22D( * ), B22E( * ),<br>
*      $                   PHI( * ), THETA( * ), WORK( * )<br>
*       DOUBLE PRECISION   U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),<br>
*      $                   V2T( LDV2T, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DBBCSD computes the CS decomposition of an orthogonal matrix in<br>
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
*>                   [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**T<br>
*>                 = [---------] [---------------] [---------]   .<br>
*>                   [    | U2 ] [  S |  C  0  0 ] [    | V2 ]<br>
*>                               [  0 |  0  0  I ]<br>
*><br>
*> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger<br>
*> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be<br>
*> transposed and/or permuted. This can be done in constant time using<br>
*> the TRANS and SIGNS options. See DORCSD for details.)<br>
*><br>
*> The bidiagonal matrices B11, B12, B21, and B22 are represented<br>
*> implicitly by angles THETA(1:Q) and PHI(1:Q-1).<br>
*><br>
*> The orthogonal matrices U1, U2, V1T, and V2T are input/output.<br>
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
*>          The number of rows and columns in X, the orthogonal matrix in<br>
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
*>          U1 is DOUBLE PRECISION array, dimension (LDU1,P)<br>
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
*>          U2 is DOUBLE PRECISION array, dimension (LDU2,M-P)<br>
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
*>          V1T is DOUBLE PRECISION array, dimension (LDV1T,Q)<br>
*>          On entry, a Q-by-Q matrix. On exit, V1T is premultiplied<br>
*>          by the transpose of the right singular vector<br>
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
*>          V2T is DOUBLE PRECISION array, dimenison (LDV2T,M-Q)<br>
*>          On entry, an (M-Q)-by-(M-Q) matrix. On exit, V2T is<br>
*>          premultiplied by the transpose of the right<br>
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
*>          When DBBCSD converges, B11D contains the cosines of THETA(1),<br>
*>          ..., THETA(Q). If DBBCSD fails to converge, then B11D<br>
*>          contains the diagonal of the partially reduced top-left<br>
*>          block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B11E<br>
*> \verbatim<br>
*>          B11E is DOUBLE PRECISION array, dimension (Q-1)<br>
*>          When DBBCSD converges, B11E contains zeros. If DBBCSD fails<br>
*>          to converge, then B11E contains the superdiagonal of the<br>
*>          partially reduced top-left block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B12D<br>
*> \verbatim<br>
*>          B12D is DOUBLE PRECISION array, dimension (Q)<br>
*>          When DBBCSD converges, B12D contains the negative sines of<br>
*>          THETA(1), ..., THETA(Q). If DBBCSD fails to converge, then<br>
*>          B12D contains the diagonal of the partially reduced top-right<br>
*>          block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B12E<br>
*> \verbatim<br>
*>          B12E is DOUBLE PRECISION array, dimension (Q-1)<br>
*>          When DBBCSD converges, B12E contains zeros. If DBBCSD fails<br>
*>          to converge, then B12E contains the subdiagonal of the<br>
*>          partially reduced top-right block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B21D<br>
*> \verbatim<br>
*>          B21D is DOUBLE PRECISION  array, dimension (Q)<br>
*>          When DBBCSD converges, B21D contains the negative sines of<br>
*>          THETA(1), ..., THETA(Q). If DBBCSD fails to converge, then<br>
*>          B21D contains the diagonal of the partially reduced bottom-left<br>
*>          block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B21E<br>
*> \verbatim<br>
*>          B21E is DOUBLE PRECISION  array, dimension (Q-1)<br>
*>          When DBBCSD converges, B21E contains zeros. If DBBCSD fails<br>
*>          to converge, then B21E contains the subdiagonal of the<br>
*>          partially reduced bottom-left block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B22D<br>
*> \verbatim<br>
*>          B22D is DOUBLE PRECISION  array, dimension (Q)<br>
*>          When DBBCSD converges, B22D contains the negative sines of<br>
*>          THETA(1), ..., THETA(Q). If DBBCSD fails to converge, then<br>
*>          B22D contains the diagonal of the partially reduced bottom-right<br>
*>          block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] B22E<br>
*> \verbatim<br>
*>          B22E is DOUBLE PRECISION  array, dimension (Q-1)<br>
*>          When DBBCSD converges, B22E contains zeros. If DBBCSD fails<br>
*>          to converge, then B22E contains the subdiagonal of the<br>
*>          partially reduced bottom-right block.<br>
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
*>          The dimension of the array WORK. LWORK >= MAX(1,8*Q).<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal size of the WORK array,<br>
*>          returns this value as the first entry of the work array, and<br>
*>          no error message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if DBBCSD did not converge, INFO specifies the number<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dbbcsd_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,CHARACTER JOBV2T,CHARACTER TRANS,INTEGER M,INTEGER P,INTEGER Q,double[] THETA,double[] PHI,double[] U1,INTEGER LDU1,double[] U2,INTEGER LDU2,double[] V1T,INTEGER LDV1T,double[] V2T,INTEGER LDV2T,double[] B11D,double[] B11E,double[] B12D,double[] B12E,double[] B21D,double[] B21E,double[] B22D,double[] B22E,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DBDSDC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DBDSDC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsdc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbdsdc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbdsdc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ,<br>
*                          WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ, UPLO<br>
*       INTEGER            INFO, LDU, LDVT, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IQ( * ), IWORK( * )<br>
*       DOUBLE PRECISION   D( * ), E( * ), Q( * ), U( LDU, * ),<br>
*      $                   VT( LDVT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DBDSDC computes the singular value decomposition (SVD) of a real<br>
*> N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,<br>
*> using a divide and conquer method, where S is a diagonal matrix<br>
*> with non-negative diagonal elements (the singular values of B), and<br>
*> U and VT are orthogonal matrices of left and right singular vectors,<br>
*> respectively. DBDSDC can be used to compute all singular values,<br>
*> and optionally, singular vectors or singular vectors in compact form.<br>
*><br>
*> This code makes very mild assumptions about floating point<br>
*> arithmetic. It will work on machines with a guard digit in<br>
*> add/subtract, or on those binary machines without guard digits<br>
*> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.<br>
*> It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.  See DLASD3 for details.<br>
*><br>
*> The code currently calls DLASDQ if singular values only are desired.<br>
*> However, it can be slightly modified to compute singular values<br>
*> using the divide and conquer method.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  B is upper bidiagonal.<br>
*>          = 'L':  B is lower bidiagonal.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPQ<br>
*> \verbatim<br>
*>          COMPQ is CHARACTER*1<br>
*>          Specifies whether singular vectors are to be computed<br>
*>          as follows:<br>
*>          = 'N':  Compute singular values only;<br>
*>          = 'P':  Compute singular values and compute singular<br>
*>                  vectors in compact form;<br>
*>          = 'I':  Compute singular values and singular vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, the n diagonal elements of the bidiagonal matrix B.<br>
*>          On exit, if INFO=0, the singular values of B.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          On entry, the elements of E contain the offdiagonal<br>
*>          elements of the bidiagonal matrix whose SVD is desired.<br>
*>          On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is DOUBLE PRECISION array, dimension (LDU,N)<br>
*>          If  COMPQ = 'I', then:<br>
*>             On exit, if INFO = 0, U contains the left singular vectors<br>
*>             of the bidiagonal matrix.<br>
*>          For other values of COMPQ, U is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>          The leading dimension of the array U.  LDU >= 1.<br>
*>          If singular vectors are desired, then LDU >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] VT<br>
*> \verbatim<br>
*>          VT is DOUBLE PRECISION array, dimension (LDVT,N)<br>
*>          If  COMPQ = 'I', then:<br>
*>             On exit, if INFO = 0, VT**T contains the right singular<br>
*>             vectors of the bidiagonal matrix.<br>
*>          For other values of COMPQ, VT is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>          The leading dimension of the array VT.  LDVT >= 1.<br>
*>          If singular vectors are desired, then LDVT >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension (LDQ)<br>
*>          If  COMPQ = 'P', then:<br>
*>             On exit, if INFO = 0, Q and IQ contain the left<br>
*>             and right singular vectors in a compact form,<br>
*>             requiring O(N log N) space instead of 2*N**2.<br>
*>             In particular, Q contains all the DOUBLE PRECISION data in<br>
*>             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))<br>
*>             words of memory, where SMLSIZ is returned by ILAENV and<br>
*>             is equal to the maximum size of the subproblems at the<br>
*>             bottom of the computation tree (usually about 25).<br>
*>          For other values of COMPQ, Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IQ<br>
*> \verbatim<br>
*>          IQ is INTEGER array, dimension (LDIQ)<br>
*>          If  COMPQ = 'P', then:<br>
*>             On exit, if INFO = 0, Q and IQ contain the left<br>
*>             and right singular vectors in a compact form,<br>
*>             requiring O(N log N) space instead of 2*N**2.<br>
*>             In particular, IQ contains all INTEGER data in<br>
*>             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))<br>
*>             words of memory, where SMLSIZ is returned by ILAENV and<br>
*>             is equal to the maximum size of the subproblems at the<br>
*>             bottom of the computation tree (usually about 25).<br>
*>          For other values of COMPQ, IQ is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          If COMPQ = 'N' then LWORK >= (4 * N).<br>
*>          If COMPQ = 'P' then LWORK >= (6 * N).<br>
*>          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (8*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  The algorithm failed to compute a singular value.<br>
*>                The update process of divide and conquer failed.<br>
*> \endverbatim<br>
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
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void dbdsdc_(CHARACTER UPLO,CHARACTER COMPQ,INTEGER N,double[] D,double[] E,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] Q,int[] IQ,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DBDSQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DBDSQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbdsqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbdsqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,<br>
*                          LDU, C, LDC, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),<br>
*      $                   VT( LDVT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DBDSQR computes the singular values and, optionally, the right and/or<br>
*> left singular vectors from the singular value decomposition (SVD) of<br>
*> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit<br>
*> zero-shift QR algorithm.  The SVD of B has the form<br>
*> <br>
*>    B = Q * S * P**T<br>
*> <br>
*> where S is the diagonal matrix of singular values, Q is an orthogonal<br>
*> matrix of left singular vectors, and P is an orthogonal matrix of<br>
*> right singular vectors.  If left singular vectors are requested, this<br>
*> subroutine actually returns U*Q instead of Q, and, if right singular<br>
*> vectors are requested, this subroutine returns P**T*VT instead of<br>
*> P**T, for given real input matrices U and VT.  When U and VT are the<br>
*> orthogonal matrices that reduce a general matrix A to bidiagonal<br>
*> form:  A = U*B*VT, as computed by DGEBRD, then<br>
*><br>
*>    A = (U*Q) * S * (P**T*VT)<br>
*><br>
*> is the SVD of A.  Optionally, the subroutine may also compute Q**T*C<br>
*> for a given real input matrix C.<br>
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
*>          matrix B. <br>
*>          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E<br>
*>          will contain the diagonal and superdiagonal elements of a<br>
*>          bidiagonal matrix orthogonally equivalent to the one given<br>
*>          as input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VT<br>
*> \verbatim<br>
*>          VT is DOUBLE PRECISION array, dimension (LDVT, NCVT)<br>
*>          On entry, an N-by-NCVT matrix VT.<br>
*>          On exit, VT is overwritten by P**T * VT.<br>
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
*>          U is DOUBLE PRECISION array, dimension (LDU, N)<br>
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
*>          C is DOUBLE PRECISION array, dimension (LDC, NCC)<br>
*>          On entry, an N-by-NCC matrix C.<br>
*>          On exit, C is overwritten by Q**T * C.<br>
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
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  If INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:<br>
*>             if NCVT = NRU = NCC = 0,<br>
*>                = 1, a split was marked by a positive value in E<br>
*>                = 2, current block of Z not diagonalized after 30*N<br>
*>                     iterations (in inner while loop)<br>
*>                = 3, termination criterion of outer while loop not met <br>
*>                     (program created more than N unreduced blocks)<br>
*>             else NCVT = NRU = NCC = 0,<br>
*>                   the algorithm did not converge; D and E contain the<br>
*>                   elements of a bidiagonal matrix which is orthogonally<br>
*>                   similar to the input matrix B;  if INFO = i, i<br>
*>                   elements of E have not converged to zero.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dbdsqr_(CHARACTER UPLO,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,double[] D,double[] E,double[] VT,INTEGER LDVT,double[] U,INTEGER LDU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DBDSVDX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DBDSVDX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsvdx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbdsvdx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbdsvdx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*     SUBROUTINE DBDSVDX( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU, <br>
*    $                    NS, S, Z, LDZ, WORK, IWORK, INFO )<br>
*<br>
*     .. Scalar Arguments ..<br>
*      CHARACTER          JOBZ, RANGE, UPLO<br>
*      INTEGER            IL, INFO, IU, LDZ, N, NS<br>
*      DOUBLE PRECISION   VL, VU<br>
*     ..<br>
*     .. Array Arguments ..<br>
*      INTEGER            IWORK( * )<br>
*      DOUBLE PRECISION   D( * ), E( * ), S( * ), WORK( * ), <br>
*                         Z( LDZ, * )<br>
*       ..<br>
*  <br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>  DBDSVDX computes the singular value decomposition (SVD) of a real<br>
*>  N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT, <br>
*>  where S is a diagonal matrix with non-negative diagonal elements <br>
*>  (the singular values of B), and U and VT are orthogonal matrices <br>
*>  of left and right singular vectors, respectively.<br>
*><br>
*>  Given an upper bidiagonal B with diagonal D = [ d_1 d_2 ... d_N ] <br>
*>  and superdiagonal E = [ e_1 e_2 ... e_N-1 ], DBDSVDX computes the <br>
*>  singular value decompositon of B through the eigenvalues and <br>
*>  eigenvectors of the N*2-by-N*2 tridiagonal matrix<br>
*>             <br>
*>        |  0  d_1                |      <br>
*>        | d_1  0  e_1            |         <br>
*>  TGK = |     e_1  0  d_2        | <br>
*>        |         d_2  .   .     |      <br>
*>        |              .   .   . |<br>
*><br>
*>  If (s,u,v) is a singular triplet of B with ||u|| = ||v|| = 1, then <br>
*>  (+/-s,q), ||q|| = 1, are eigenpairs of TGK, with q = P * ( u' +/-v' ) / <br>
*>  sqrt(2) = ( v_1 u_1 v_2 u_2 ... v_n u_n ) / sqrt(2), and <br>
*>  P = [ e_{n+1} e_{1} e_{n+2} e_{2} ... ]. <br>
*><br>
*>  Given a TGK matrix, one can either a) compute -s,-v and change signs <br>
*>  so that the singular values (and corresponding vectors) are already in <br>
*>  descending order (as in DGESVD/DGESDD) or b) compute s,v and reorder <br>
*>  the values (and corresponding vectors). DBDSVDX implements a) by <br>
*>  calling DSTEVX (bisection plus inverse iteration, to be replaced <br>
*>  with a version of the Multiple Relative Robust Representation <br>
*>  algorithm. (See P. Willems and B. Lang, A framework for the MR^3 <br>
*>  algorithm: theory and implementation, SIAM J. Sci. Comput., <br>
*>  35:740-766, 2013.)<br>
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
*> \param[in] JOBZ<br>
*> \verbatim<br>
*>          JOBZ is CHARACTER*1<br>
*>          = 'N':  Compute singular values only;<br>
*>          = 'V':  Compute singular values and singular vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RANGE<br>
*> \verbatim<br>
*>          RANGE is CHARACTER*1<br>
*>          = 'A': all singular values will be found.<br>
*>          = 'V': all singular values in the half-open interval [VL,VU)<br>
*>                 will be found.<br>
*>          = 'I': the IL-th through IU-th singular values will be found.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the bidiagonal matrix.  N >= 0.<br>
*> \endverbatim<br>
*> <br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The n diagonal elements of the bidiagonal matrix B.<br>
*> \endverbatim<br>
*> <br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (max(1,N-1))<br>
*>          The (n-1) superdiagonal elements of the bidiagonal matrix<br>
*>          B in elements 1 to N-1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>         VL is DOUBLE PRECISION<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for singular values. VU > VL.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>         VU is DOUBLE PRECISION<br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for singular values. VU > VL.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          smallest singular value to be returned.<br>
*>          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          largest singular value to be returned.<br>
*>          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] NS<br>
*> \verbatim<br>
*>          NS is INTEGER<br>
*>          The total number of singular values found.  0 <= NS <= N.<br>
*>          If RANGE = 'A', NS = N, and if RANGE = 'I', NS = IU-IL+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (N)<br>
*>          The first NS elements contain the selected singular values in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is DOUBLE PRECISION array, dimension (2*N,K) )<br>
*>          If JOBZ = 'V', then if INFO = 0 the first NS columns of Z<br>
*>          contain the singular vectors of the matrix B corresponding to <br>
*>          the selected singular values, with U in rows 1 to N and V<br>
*>          in rows N+1 to N*2, i.e.<br>
*>          Z = [ U ] <br>
*>              [ V ]<br>
*>          If JOBZ = 'N', then Z is not referenced.    <br>
*>          Note: The user must ensure that at least K = NS+1 columns are <br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of <br>
*>          NS is not known in advance and an upper bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z. LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(2,N*2).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (14*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (12*N)<br>
*>          If JOBZ = 'V', then if INFO = 0, the first NS elements of<br>
*>          IWORK are zero. If INFO > 0, then IWORK contains the indices <br>
*>          of the eigenvectors that failed to converge in DSTEVX.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, then i eigenvectors failed to converge<br>
*>                   in DSTEVX. The indices of the eigenvectors<br>
*>                   (as returned by DSTEVX) are stored in the<br>
*>                   array IWORK.<br>
*>                if INFO = N*2 + 1, an internal error occurred.<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHEReigen<br>
*<br>
*  =====================================================================   <br>
*/
	public void dbdsvdx_(CHARACTER UPLO,CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,double[] D,double[] E,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,INTEGER NS,double[] S,double[] Z,INTEGER LDZ,double[] WORK,int[] IWORK,INTEGER INFO);

}