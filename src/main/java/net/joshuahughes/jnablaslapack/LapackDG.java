package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDG extends Library
{

	public static LapackDG instance = (LapackDG) Native.loadLibrary("liblapack",LapackDG.class);

/**
*> \brief \b DGBBRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBBRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbbrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbbrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbbrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q,<br>
*                          LDQ, PT, LDPT, C, LDC, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          VECT<br>
*       INTEGER            INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   AB( LDAB, * ), C( LDC, * ), D( * ), E( * ),<br>
*      $                   PT( LDPT, * ), Q( LDQ, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBBRD reduces a real general m-by-n band matrix A to upper<br>
*> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.<br>
*><br>
*> The routine computes B, and optionally forms Q or P**T, or computes<br>
*> Q**T*C for a given matrix C.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] VECT<br>
*> \verbatim<br>
*>          VECT is CHARACTER*1<br>
*>          Specifies whether or not the matrices Q and P**T are to be<br>
*>          formed.<br>
*>          = 'N': do not form Q or P**T;<br>
*>          = 'Q': form Q only;<br>
*>          = 'P': form P**T only;<br>
*>          = 'B': form both.<br>
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
*> \param[in] NCC<br>
*> \verbatim<br>
*>          NCC is INTEGER<br>
*>          The number of columns of the matrix C.  NCC >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>          The number of subdiagonals of the matrix A. KL >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>          The number of superdiagonals of the matrix A. KU >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          On entry, the m-by-n band matrix A, stored in rows 1 to<br>
*>          KL+KU+1. The j-th column of A is stored in the j-th column of<br>
*>          the array AB as follows:<br>
*>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).<br>
*>          On exit, A is overwritten by values generated during the<br>
*>          reduction.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array A. LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The diagonal elements of the bidiagonal matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (min(M,N)-1)<br>
*>          The superdiagonal elements of the bidiagonal matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension (LDQ,M)<br>
*>          If VECT = 'Q' or 'B', the m-by-m orthogonal matrix Q.<br>
*>          If VECT = 'N' or 'P', the array Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.<br>
*>          LDQ >= max(1,M) if VECT = 'Q' or 'B'; LDQ >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PT<br>
*> \verbatim<br>
*>          PT is DOUBLE PRECISION array, dimension (LDPT,N)<br>
*>          If VECT = 'P' or 'B', the n-by-n orthogonal matrix P'.<br>
*>          If VECT = 'N' or 'Q', the array PT is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDPT<br>
*> \verbatim<br>
*>          LDPT is INTEGER<br>
*>          The leading dimension of the array PT.<br>
*>          LDPT >= max(1,N) if VECT = 'P' or 'B'; LDPT >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,NCC)<br>
*>          On entry, an m-by-ncc matrix C.<br>
*>          On exit, C is overwritten by Q**T*C.<br>
*>          C is not referenced if NCC = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C.<br>
*>          LDC >= max(1,M) if NCC > 0; LDC >= 1 if NCC = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (2*max(M,N))<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgbbrd_(CHARACTER VECT,INTEGER M,INTEGER N,INTEGER NCC,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] D,double[] E,double[] Q,INTEGER LDQ,double[] PT,INTEGER LDPT,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGBCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND,<br>
*                          WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            INFO, KL, KU, LDAB, N<br>
*       DOUBLE PRECISION   ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBCON estimates the reciprocal of the condition number of a real<br>
*> general band matrix A, in either the 1-norm or the infinity-norm,<br>
*> using the LU factorization computed by DGBTRF.<br>
*><br>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the<br>
*> condition number is computed as<br>
*>    RCOND = 1 / ( norm(A) * norm(inv(A)) ).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies whether the 1-norm condition number or the<br>
*>          infinity-norm condition number is required:<br>
*>          = '1' or 'O':  1-norm;<br>
*>          = 'I':         Infinity-norm.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
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
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          Details of the LU factorization of the band matrix A, as<br>
*>          computed by DGBTRF.  U is stored as an upper triangular band<br>
*>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and<br>
*>          the multipliers used during the factorization are stored in<br>
*>          rows KL+KU+2 to 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices; for 1 <= i <= N, row i of the matrix was<br>
*>          interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is DOUBLE PRECISION<br>
*>          If NORM = '1' or 'O', the 1-norm of the original matrix A.<br>
*>          If NORM = 'I', the infinity-norm of the original matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(norm(A) * norm(inv(A))).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgbcon_(CHARACTER NORM,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGBEQU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBEQU + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbequ.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbequ.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbequ.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,<br>
*                          AMAX, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, KL, KU, LDAB, M, N<br>
*       DOUBLE PRECISION   AMAX, COLCND, ROWCND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   AB( LDAB, * ), C( * ), R( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBEQU computes row and column scalings intended to equilibrate an<br>
*> M-by-N band matrix A and reduce its condition number.  R returns the<br>
*> row scale factors and C the column scale factors, chosen to try to<br>
*> make the largest element in each row and column of the matrix B with<br>
*> elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.<br>
*><br>
*> R(i) and C(j) are restricted to be between SMLNUM = smallest safe<br>
*> number and BIGNUM = largest safe number.  Use of these scaling<br>
*> factors is not guaranteed to reduce the condition number of A but<br>
*> works well in practice.<br>
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
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th<br>
*>          column of A is stored in the j-th column of the array AB as<br>
*>          follows:<br>
*>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (M)<br>
*>          If INFO = 0, or INFO > M, R contains the row scale factors<br>
*>          for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, C contains the column scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ROWCND<br>
*> \verbatim<br>
*>          ROWCND is DOUBLE PRECISION<br>
*>          If INFO = 0 or INFO > M, ROWCND contains the ratio of the<br>
*>          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and<br>
*>          AMAX is neither too large nor too small, it is not worth<br>
*>          scaling by R.<br>
*> \endverbatim<br>
*><br>
*> \param[out] COLCND<br>
*> \verbatim<br>
*>          COLCND is DOUBLE PRECISION<br>
*>          If INFO = 0, COLCND contains the ratio of the smallest<br>
*>          C(i) to the largest C(i).  If COLCND >= 0.1, it is not<br>
*>          worth scaling by C.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AMAX<br>
*> \verbatim<br>
*>          AMAX is DOUBLE PRECISION<br>
*>          Absolute value of largest matrix element.  If AMAX is very<br>
*>          close to overflow or very close to underflow, the matrix<br>
*>          should be scaled.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is<br>
*>                <= M:  the i-th row of A is exactly zero<br>
*>                >  M:  the (i-M)-th column of A is exactly zero<br>
*> \endverbatim<br>
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
*> \ingroup doubleGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgbequ_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,INTEGER INFO);
/**
*> \brief \b DGBEQUB<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBEQUB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbequb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbequb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbequb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,<br>
*                           AMAX, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, KL, KU, LDAB, M, N<br>
*       DOUBLE PRECISION   AMAX, COLCND, ROWCND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   AB( LDAB, * ), C( * ), R( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBEQUB computes row and column scalings intended to equilibrate an<br>
*> M-by-N matrix A and reduce its condition number.  R returns the row<br>
*> scale factors and C the column scale factors, chosen to try to make<br>
*> the largest element in each row and column of the matrix B with<br>
*> elements B(i,j)=R(i)*A(i,j)*C(j) have an absolute value of at most<br>
*> the radix.<br>
*><br>
*> R(i) and C(j) are restricted to be a power of the radix between<br>
*> SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use<br>
*> of these scaling factors is not guaranteed to reduce the condition<br>
*> number of A but works well in practice.<br>
*><br>
*> This routine differs from DGEEQU by restricting the scaling factors<br>
*> to a power of the radix.  Baring over- and underflow, scaling by<br>
*> these factors introduces no additional rounding errors.  However, the<br>
*> scaled entries' magnitured are no longer approximately 1 but lie<br>
*> between sqrt(radix) and 1/sqrt(radix).<br>
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
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.<br>
*>          The j-th column of A is stored in the j-th column of the<br>
*>          array AB as follows:<br>
*>          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array A.  LDAB >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (M)<br>
*>          If INFO = 0 or INFO > M, R contains the row scale factors<br>
*>          for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0,  C contains the column scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ROWCND<br>
*> \verbatim<br>
*>          ROWCND is DOUBLE PRECISION<br>
*>          If INFO = 0 or INFO > M, ROWCND contains the ratio of the<br>
*>          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and<br>
*>          AMAX is neither too large nor too small, it is not worth<br>
*>          scaling by R.<br>
*> \endverbatim<br>
*><br>
*> \param[out] COLCND<br>
*> \verbatim<br>
*>          COLCND is DOUBLE PRECISION<br>
*>          If INFO = 0, COLCND contains the ratio of the smallest<br>
*>          C(i) to the largest C(i).  If COLCND >= 0.1, it is not<br>
*>          worth scaling by C.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AMAX<br>
*> \verbatim<br>
*>          AMAX is DOUBLE PRECISION<br>
*>          Absolute value of largest matrix element.  If AMAX is very<br>
*>          close to overflow or very close to underflow, the matrix<br>
*>          should be scaled.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i,  and i is<br>
*>                <= M:  the i-th row of A is exactly zero<br>
*>                >  M:  the (i-M)-th column of A is exactly zero<br>
*> \endverbatim<br>
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
*> \ingroup doubleGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgbequb_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,INTEGER INFO);
/**
*> \brief \b DGBRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB,<br>
*                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),<br>
*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBRFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is banded, and provides<br>
*> error bounds and backward error estimates for the solution.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          The original band matrix A, stored in rows 1 to KL+KU+1.<br>
*>          The j-th column of A is stored in the j-th column of the<br>
*>          array AB as follows:<br>
*>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AFB<br>
*> \verbatim<br>
*>          AFB is DOUBLE PRECISION array, dimension (LDAFB,N)<br>
*>          Details of the LU factorization of the band matrix A, as<br>
*>          computed by DGBTRF.  U is stored as an upper triangular band<br>
*>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and<br>
*>          the multipliers used during the factorization are stored in<br>
*>          rows KL+KU+2 to 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>          The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices from DGBTRF; for 1<=i<=N, row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          The right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by DGBTRS.<br>
*>          On exit, the improved solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The estimated forward error bound for each solution vector<br>
*>          X(j) (the j-th column of the solution matrix X).<br>
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)<br>
*>          is an estimated upper bound for the magnitude of the largest<br>
*>          element in (X(j) - XTRUE) divided by the magnitude of the<br>
*>          largest element in X(j).  The estimate is as reliable as<br>
*>          the estimate for RCOND, and is almost always a slight<br>
*>          overestimate of the true error.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  ITMAX is the maximum number of steps of iterative refinement.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgbrfs_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGBRFSX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBRFSX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbrfsx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbrfsx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbrfsx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBRFSX( TRANS, EQUED, N, KL, KU, NRHS, AB, LDAB, AFB,<br>
*                           LDAFB, IPIV, R, C, B, LDB, X, LDX, RCOND,<br>
*                           BERR, N_ERR_BNDS, ERR_BNDS_NORM,<br>
*                           ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK,<br>
*                           INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS, EQUED<br>
*       INTEGER            INFO, LDAB, LDAFB, LDB, LDX, N, KL, KU, NRHS,<br>
*      $                   NPARAMS, N_ERR_BNDS<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),<br>
*      $                   X( LDX , * ),WORK( * )<br>
*       DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ),<br>
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
*>    DGBRFSX improves the computed solution to a system of linear<br>
*>    equations and provides error bounds and backward error estimates<br>
*>    for the solution.  In addition to normwise error bound, the code<br>
*>    provides maximum componentwise error bound if possible.  See<br>
*>    comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the<br>
*>    error bounds.<br>
*><br>
*>    The original system of linear equations may have been equilibrated<br>
*>    before calling this routine, as described by arguments EQUED, R<br>
*>    and C below. In this case, the solution and error bounds returned<br>
*>    are for the original unequilibrated system.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \verbatim<br>
*>     Some optional parameters are bundled in the PARAMS array.  These<br>
*>     settings determine how refinement is performed, but often the<br>
*>     defaults are acceptable.  If the defaults are acceptable, users<br>
*>     can pass NPARAMS = 0 which prevents the source code from accessing<br>
*>     the PARAMS argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>     Specifies the form of the system of equations:<br>
*>       = 'N':  A * X = B     (No transpose)<br>
*>       = 'T':  A**T * X = B  (Transpose)<br>
*>       = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>     Specifies the form of equilibration that was done to A<br>
*>     before calling this routine. This is needed to compute<br>
*>     the solution and error bounds correctly.<br>
*>       = 'N':  No equilibration<br>
*>       = 'R':  Row equilibration, i.e., A has been premultiplied by<br>
*>               diag(R).<br>
*>       = 'C':  Column equilibration, i.e., A has been postmultiplied<br>
*>               by diag(C).<br>
*>       = 'B':  Both row and column equilibration, i.e., A has been<br>
*>               replaced by diag(R) * A * diag(C).<br>
*>               The right hand side B has been changed accordingly.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The order of the matrix A.  N >= 0.<br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right hand sides, i.e., the number of columns<br>
*>     of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>     The original band matrix A, stored in rows 1 to KL+KU+1.<br>
*>     The j-th column of A is stored in the j-th column of the<br>
*>     array AB as follows:<br>
*>     AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).<br>
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
*>          AFB is DOUBLE PRECISION array, dimension (LDAFB,N)<br>
*>     Details of the LU factorization of the band matrix A, as<br>
*>     computed by DGBTRF.  U is stored as an upper triangular band<br>
*>     matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and<br>
*>     the multipliers used during the factorization are stored in<br>
*>     rows KL+KU+2 to 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>     The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     The pivot indices from DGETRF; for 1<=i<=N, row i of the<br>
*>     matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (N)<br>
*>     The row scale factors for A.  If EQUED = 'R' or 'B', A is<br>
*>     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R<br>
*>     is not accessed.  R is an input argument if FACT = 'F';<br>
*>     otherwise, R is an output argument.  If FACT = 'F' and<br>
*>     EQUED = 'R' or 'B', each element of R must be positive.<br>
*>     If R is output, each element of R is a power of the radix.<br>
*>     If R is input, each element of R should be a power of the radix<br>
*>     to ensure a reliable solution and error estimates. Scaling by<br>
*>     powers of the radix does not cause rounding errors unless the<br>
*>     result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>     The column scale factors for A.  If EQUED = 'C' or 'B', A is<br>
*>     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C<br>
*>     is not accessed.  C is an input argument if FACT = 'F';<br>
*>     otherwise, C is an output argument.  If FACT = 'F' and<br>
*>     EQUED = 'C' or 'B', each element of C must be positive.<br>
*>     If C is output, each element of C is a power of the radix.<br>
*>     If C is input, each element of C should be a power of the radix<br>
*>     to ensure a reliable solution and error estimates. Scaling by<br>
*>     powers of the radix does not cause rounding errors unless the<br>
*>     result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>     The right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>     On entry, the solution matrix X, as computed by DGETRS.<br>
*>     On exit, the improved solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>     The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>     Componentwise relative backward error.  This is the<br>
*>     componentwise relative backward error of each solution vector X(j)<br>
*>     (i.e., the smallest relative change in any element of A or B that<br>
*>     makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_ERR_BNDS<br>
*> \verbatim<br>
*>          N_ERR_BNDS is INTEGER<br>
*>     Number of error bounds to return for each right hand side<br>
*>     and each type (normwise or componentwise).  See ERR_BNDS_NORM and<br>
*>     ERR_BNDS_COMP below.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_NORM<br>
*> \verbatim<br>
*>          ERR_BNDS_NORM is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>              sqrt(n) * dlamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * dlamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * dlamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_COMP<br>
*> \verbatim<br>
*>          ERR_BNDS_COMP is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>              sqrt(n) * dlamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * dlamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * dlamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NPARAMS<br>
*> \verbatim<br>
*>          NPARAMS is INTEGER<br>
*>     Specifies the number of parameters set in PARAMS.  If .LE. 0, the<br>
*>     PARAMS array is never referenced and default values are used.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] PARAMS<br>
*> \verbatim<br>
*>          PARAMS is DOUBLE PRECISION array, dimension (NPARAMS)<br>
*>     Specifies algorithm parameters.  If an entry is .LT. 0.0, then<br>
*>     that entry will be filled with default value used for that<br>
*>     parameter.  Only positions up to NPARAMS are accessed; defaults<br>
*>     are used for higher-numbered parameters.<br>
*><br>
*>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative<br>
*>            refinement or not.<br>
*>         Default: 1.0D+0<br>
*>            = 0.0 : No refinement is performed, and no error bounds are<br>
*>                    computed.<br>
*>            = 1.0 : Use the double-precision refinement algorithm,<br>
*>                    possibly with doubled-single computations if the<br>
*>                    compilation environment does not support DOUBLE<br>
*>                    PRECISION.<br>
*>              (other values are reserved for future use)<br>
*><br>
*>       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual<br>
*>            computations allowed for refinement.<br>
*>         Default: 10<br>
*>         Aggressive: Set to 100 to permit convergence using approximate<br>
*>                     factorizations or factorizations other than LU. If<br>
*>                     the factorization uses a technique other than<br>
*>                     Gaussian elimination, the guarantees in<br>
*>                     err_bnds_norm and err_bnds_comp may no longer be<br>
*>                     trustworthy.<br>
*><br>
*>       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code<br>
*>            will attempt to find a solution with small componentwise<br>
*>            relative error in the double-precision algorithm.  Positive<br>
*>            is true, 0.0 is false.<br>
*>         Default: 1.0 (attempt componentwise convergence)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit. The solution to every right-hand side is<br>
*>         guaranteed.<br>
*>       < 0:  If INFO = -i, the i-th argument had an illegal value<br>
*>       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization<br>
*>         has been completed, but the factor U is exactly singular, so<br>
*>         the solution and error bounds could not be computed. RCOND = 0<br>
*>         is returned.<br>
*>       = N+J: The solution corresponding to the Jth right-hand side is<br>
*>         not guaranteed. The solutions corresponding to other right-<br>
*>         hand sides K with K > J may not be guaranteed as well, but<br>
*>         only the first such right-hand side is reported. If a small<br>
*>         componentwise error is not requested (PARAMS(3) = 0.0) then<br>
*>         the Jth right-hand side is the first with a normwise error<br>
*>         bound that is not guaranteed (the smallest J such<br>
*>         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)<br>
*>         the Jth right-hand side is the first with either a normwise or<br>
*>         componentwise error bound that is not guaranteed (the smallest<br>
*>         J such that either ERR_BNDS_NORM(J,1) = 0.0 or<br>
*>         ERR_BNDS_COMP(J,1) = 0.0). See the definition of<br>
*>         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information<br>
*>         about all of the right-hand sides check ERR_BNDS_NORM or<br>
*>         ERR_BNDS_COMP.<br>
*> \endverbatim<br>
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
	public void dgbrfsx_(CHARACTER TRANS,CHARACTER EQUED,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> DGBSV computes the solution to system of linear equations A * X = B for GB matrices</b> (simple driver)<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBSV computes the solution to a real system of linear equations<br>
*> A * X = B, where A is a band matrix of order N with KL subdiagonals<br>
*> and KU superdiagonals, and X and B are N-by-NRHS matrices.<br>
*><br>
*> The LU decomposition with partial pivoting and row interchanges is<br>
*> used to factor A as A = L * U, where L is a product of permutation<br>
*> and unit lower triangular matrices with KL subdiagonals, and U is<br>
*> upper triangular with KL+KU superdiagonals.  The factored form of A<br>
*> is then used to solve the system of equations A * X = B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of linear equations, i.e., the order of the<br>
*>          matrix A.  N >= 0.<br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          On entry, the matrix A in band storage, in rows KL+1 to<br>
*>          2*KL+KU+1; rows 1 to KL of the array need not be set.<br>
*>          The j-th column of A is stored in the j-th column of the<br>
*>          array AB as follows:<br>
*>          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)<br>
*>          On exit, details of the factorization: U is stored as an<br>
*>          upper triangular band matrix with KL+KU superdiagonals in<br>
*>          rows 1 to KL+KU+1, and the multipliers used during the<br>
*>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.<br>
*>          See below for further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices that define the permutation matrix P;<br>
*>          row i of the matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the N-by-NRHS right hand side matrix B.<br>
*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization<br>
*>                has been completed, but the factor U is exactly<br>
*>                singular, and the solution has not been computed.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGBsolve<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The band storage scheme is illustrated by the following example, when<br>
*>  M = N = 6, KL = 2, KU = 1:<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>      *    *    *    +    +    +       *    *    *   u14  u25  u36<br>
*>      *    *    +    +    +    +       *    *   u13  u24  u35  u46<br>
*>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56<br>
*>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66<br>
*>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *<br>
*>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *<br>
*><br>
*>  Array elements marked * are not used by the routine; elements marked<br>
*>  + need not be set on entry, but are required by the routine to store<br>
*>  elements of U because of fill-in resulting from the row interchanges.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgbsv_(INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> DGBSVX computes the solution to system of linear equations A * X = B for GB matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbsvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbsvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbsvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,<br>
*                          LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX,<br>
*                          RCOND, FERR, BERR, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, TRANS<br>
*       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),<br>
*      $                   BERR( * ), C( * ), FERR( * ), R( * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBSVX uses the LU factorization to compute the solution to a real<br>
*> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,<br>
*> where A is a band matrix of order N with KL subdiagonals and KU<br>
*> superdiagonals, and X and B are N-by-NRHS matrices.<br>
*><br>
*> Error bounds on the solution and a condition estimate are also<br>
*> provided.<br>
*> \endverbatim<br>
*<br>
*> \par Description:<br>
*  =================<br>
*><br>
*> \verbatim<br>
*><br>
*> The following steps are performed by this subroutine:<br>
*><br>
*> 1. If FACT = 'E', real scaling factors are computed to equilibrate<br>
*>    the system:<br>
*>       TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B<br>
*>       TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B<br>
*>       TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B<br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')<br>
*>    or diag(C)*B (if TRANS = 'T' or 'C').<br>
*><br>
*> 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the<br>
*>    matrix A (after equilibration if FACT = 'E') as<br>
*>       A = L * U,<br>
*>    where L is a product of permutation and unit lower triangular<br>
*>    matrices with KL subdiagonals, and U is upper triangular with<br>
*>    KL+KU superdiagonals.<br>
*><br>
*> 3. If some U(i,i)=0, so that U is exactly singular, then the routine<br>
*>    returns with INFO = i. Otherwise, the factored form of A is used<br>
*>    to estimate the condition number of the matrix A.  If the<br>
*>    reciprocal of the condition number is less than machine precision,<br>
*>    INFO = N+1 is returned as a warning, but the routine still goes on<br>
*>    to solve for X and compute error bounds as described below.<br>
*><br>
*> 4. The system of equations is solved for X using the factored form<br>
*>    of A.<br>
*><br>
*> 5. Iterative refinement is applied to improve the computed solution<br>
*>    matrix and calculate error bounds and backward error estimates<br>
*>    for it.<br>
*><br>
*> 6. If equilibration was used, the matrix X is premultiplied by<br>
*>    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so<br>
*>    that it solves the original system before equilibration.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] FACT<br>
*> \verbatim<br>
*>          FACT is CHARACTER*1<br>
*>          Specifies whether or not the factored form of the matrix A is<br>
*>          supplied on entry, and if not, whether the matrix A should be<br>
*>          equilibrated before it is factored.<br>
*>          = 'F':  On entry, AFB and IPIV contain the factored form of<br>
*>                  A.  If EQUED is not 'N', the matrix A has been<br>
*>                  equilibrated with scaling factors given by R and C.<br>
*>                  AB, AFB, and IPIV are not modified.<br>
*>          = 'N':  The matrix A will be copied to AFB and factored.<br>
*>          = 'E':  The matrix A will be equilibrated if necessary, then<br>
*>                  copied to AFB and factored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations.<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of linear equations, i.e., the order of the<br>
*>          matrix A.  N >= 0.<br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.<br>
*>          The j-th column of A is stored in the j-th column of the<br>
*>          array AB as follows:<br>
*>          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)<br>
*><br>
*>          If FACT = 'F' and EQUED is not 'N', then A must have been<br>
*>          equilibrated by the scaling factors in R and/or C.  AB is not<br>
*>          modified if FACT = 'F' or 'N', or if FACT = 'E' and<br>
*>          EQUED = 'N' on exit.<br>
*><br>
*>          On exit, if EQUED .ne. 'N', A is scaled as follows:<br>
*>          EQUED = 'R':  A := diag(R) * A<br>
*>          EQUED = 'C':  A := A * diag(C)<br>
*>          EQUED = 'B':  A := diag(R) * A * diag(C).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AFB<br>
*> \verbatim<br>
*>          AFB is DOUBLE PRECISION array, dimension (LDAFB,N)<br>
*>          If FACT = 'F', then AFB is an input argument and on entry<br>
*>          contains details of the LU factorization of the band matrix<br>
*>          A, as computed by DGBTRF.  U is stored as an upper triangular<br>
*>          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,<br>
*>          and the multipliers used during the factorization are stored<br>
*>          in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is<br>
*>          the factored form of the equilibrated matrix A.<br>
*><br>
*>          If FACT = 'N', then AFB is an output argument and on exit<br>
*>          returns details of the LU factorization of A.<br>
*><br>
*>          If FACT = 'E', then AFB is an output argument and on exit<br>
*>          returns details of the LU factorization of the equilibrated<br>
*>          matrix A (see the description of AB for the form of the<br>
*>          equilibrated matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          If FACT = 'F', then IPIV is an input argument and on entry<br>
*>          contains the pivot indices from the factorization A = L*U<br>
*>          as computed by DGBTRF; row i of the matrix was interchanged<br>
*>          with row IPIV(i).<br>
*><br>
*>          If FACT = 'N', then IPIV is an output argument and on exit<br>
*>          contains the pivot indices from the factorization A = L*U<br>
*>          of the original matrix A.<br>
*><br>
*>          If FACT = 'E', then IPIV is an output argument and on exit<br>
*>          contains the pivot indices from the factorization A = L*U<br>
*>          of the equilibrated matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies the form of equilibration that was done.<br>
*>          = 'N':  No equilibration (always true if FACT = 'N').<br>
*>          = 'R':  Row equilibration, i.e., A has been premultiplied by<br>
*>                  diag(R).<br>
*>          = 'C':  Column equilibration, i.e., A has been postmultiplied<br>
*>                  by diag(C).<br>
*>          = 'B':  Both row and column equilibration, i.e., A has been<br>
*>                  replaced by diag(R) * A * diag(C).<br>
*>          EQUED is an input argument if FACT = 'F'; otherwise, it is an<br>
*>          output argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (N)<br>
*>          The row scale factors for A.  If EQUED = 'R' or 'B', A is<br>
*>          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R<br>
*>          is not accessed.  R is an input argument if FACT = 'F';<br>
*>          otherwise, R is an output argument.  If FACT = 'F' and<br>
*>          EQUED = 'R' or 'B', each element of R must be positive.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>          The column scale factors for A.  If EQUED = 'C' or 'B', A is<br>
*>          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C<br>
*>          is not accessed.  C is an input argument if FACT = 'F';<br>
*>          otherwise, C is an output argument.  If FACT = 'F' and<br>
*>          EQUED = 'C' or 'B', each element of C must be positive.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the right hand side matrix B.<br>
*>          On exit,<br>
*>          if EQUED = 'N', B is not modified;<br>
*>          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by<br>
*>          diag(R)*B;<br>
*>          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is<br>
*>          overwritten by diag(C)*B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X<br>
*>          to the original system of equations.  Note that A and B are<br>
*>          modified on exit if EQUED .ne. 'N', and the solution to the<br>
*>          equilibrated system is inv(diag(C))*X if TRANS = 'N' and<br>
*>          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'<br>
*>          and EQUED = 'R' or 'B'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          The estimate of the reciprocal condition number of the matrix<br>
*>          A after equilibration (if done).  If RCOND is less than the<br>
*>          machine precision (in particular, if RCOND = 0), the matrix<br>
*>          is singular to working precision.  This condition is<br>
*>          indicated by a return code of INFO > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The estimated forward error bound for each solution vector<br>
*>          X(j) (the j-th column of the solution matrix X).<br>
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)<br>
*>          is an estimated upper bound for the magnitude of the largest<br>
*>          element in (X(j) - XTRUE) divided by the magnitude of the<br>
*>          largest element in X(j).  The estimate is as reliable as<br>
*>          the estimate for RCOND, and is almost always a slight<br>
*>          overestimate of the true error.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (3*N)<br>
*>          On exit, WORK(1) contains the reciprocal pivot growth<br>
*>          factor norm(A)/norm(U). The "max absolute element" norm is<br>
*>          used. If WORK(1) is much less than 1, then the stability<br>
*>          of the LU factorization of the (equilibrated) matrix A<br>
*>          could be poor. This also means that the solution X, condition<br>
*>          estimator RCOND, and forward error bound FERR could be<br>
*>          unreliable. If factorization fails with 0<INFO<=N, then<br>
*>          WORK(1) contains the reciprocal pivot growth factor for the<br>
*>          leading INFO columns of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is<br>
*>                <= N:  U(i,i) is exactly zero.  The factorization<br>
*>                       has been completed, but the factor U is exactly<br>
*>                       singular, so the solution and error bounds<br>
*>                       could not be computed. RCOND = 0 is returned.<br>
*>                = N+1: U is nonsingular, but RCOND is less than machine<br>
*>                       precision, meaning that the matrix is singular<br>
*>                       to working precision.  Nevertheless, the<br>
*>                       solution and error bounds are computed because<br>
*>                       there are a number of situations where the<br>
*>                       computed solution can be more accurate than the<br>
*>                       value of RCOND would suggest.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGBsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgbsvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,CHARACTER EQUED,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> DGBSVXX computes the solution to system of linear equations A * X = B for GB matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBSVXX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbsvxx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbsvxx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbsvxx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBSVXX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,<br>
*                           LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX,<br>
*                           RCOND, RPVGRW, BERR, N_ERR_BNDS,<br>
*                           ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS,<br>
*                           WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, TRANS<br>
*       INTEGER            INFO, LDAB, LDAFB, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS, KL, KU<br>
*       DOUBLE PRECISION   RCOND, RPVGRW<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),<br>
*      $                   X( LDX , * ),WORK( * )<br>
*       DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ),<br>
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
*>    DGBSVXX uses the LU factorization to compute the solution to a<br>
*>    double precision system of linear equations  A * X = B,  where A is an<br>
*>    N-by-N matrix and X and B are N-by-NRHS matrices.<br>
*><br>
*>    If requested, both normwise and maximum componentwise error bounds<br>
*>    are returned. DGBSVXX will return a solution with a tiny<br>
*>    guaranteed error (O(eps) where eps is the working machine<br>
*>    precision) unless the matrix is very ill-conditioned, in which<br>
*>    case a warning is returned. Relevant condition numbers also are<br>
*>    calculated and returned.<br>
*><br>
*>    DGBSVXX accepts user-provided factorizations and equilibration<br>
*>    factors; see the definitions of the FACT and EQUED options.<br>
*>    Solving with refinement and using a factorization from a previous<br>
*>    DGBSVXX call will also produce a solution with either O(eps)<br>
*>    errors or warnings, but we cannot make that claim for general<br>
*>    user-provided factorizations and equilibration factors if they<br>
*>    differ from what DGBSVXX would itself produce.<br>
*> \endverbatim<br>
*<br>
*> \par Description:<br>
*  =================<br>
*><br>
*> \verbatim<br>
*><br>
*>    The following steps are performed:<br>
*><br>
*>    1. If FACT = 'E', double precision scaling factors are computed to equilibrate<br>
*>    the system:<br>
*><br>
*>      TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B<br>
*>      TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B<br>
*>      TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B<br>
*><br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')<br>
*>    or diag(C)*B (if TRANS = 'T' or 'C').<br>
*><br>
*>    2. If FACT = 'N' or 'E', the LU decomposition is used to factor<br>
*>    the matrix A (after equilibration if FACT = 'E') as<br>
*><br>
*>      A = P * L * U,<br>
*><br>
*>    where P is a permutation matrix, L is a unit lower triangular<br>
*>    matrix, and U is upper triangular.<br>
*><br>
*>    3. If some U(i,i)=0, so that U is exactly singular, then the<br>
*>    routine returns with INFO = i. Otherwise, the factored form of A<br>
*>    is used to estimate the condition number of the matrix A (see<br>
*>    argument RCOND). If the reciprocal of the condition number is less<br>
*>    than machine precision, the routine still goes on to solve for X<br>
*>    and compute error bounds as described below.<br>
*><br>
*>    4. The system of equations is solved for X using the factored form<br>
*>    of A.<br>
*><br>
*>    5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),<br>
*>    the routine will use iterative refinement to try to get a small<br>
*>    error and error bounds.  Refinement calculates the residual to at<br>
*>    least twice the working precision.<br>
*><br>
*>    6. If equilibration was used, the matrix X is premultiplied by<br>
*>    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so<br>
*>    that it solves the original system before equilibration.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \verbatim<br>
*>     Some optional parameters are bundled in the PARAMS array.  These<br>
*>     settings determine how refinement is performed, but often the<br>
*>     defaults are acceptable.  If the defaults are acceptable, users<br>
*>     can pass NPARAMS = 0 which prevents the source code from accessing<br>
*>     the PARAMS argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in] FACT<br>
*> \verbatim<br>
*>          FACT is CHARACTER*1<br>
*>     Specifies whether or not the factored form of the matrix A is<br>
*>     supplied on entry, and if not, whether the matrix A should be<br>
*>     equilibrated before it is factored.<br>
*>       = 'F':  On entry, AF and IPIV contain the factored form of A.<br>
*>               If EQUED is not 'N', the matrix A has been<br>
*>               equilibrated with scaling factors given by R and C.<br>
*>               A, AF, and IPIV are not modified.<br>
*>       = 'N':  The matrix A will be copied to AF and factored.<br>
*>       = 'E':  The matrix A will be equilibrated if necessary, then<br>
*>               copied to AF and factored.<br>
*> \endverbatim<br>
*><br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right hand sides, i.e., the number of columns<br>
*>     of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.<br>
*>     The j-th column of A is stored in the j-th column of the<br>
*>     array AB as follows:<br>
*>     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)<br>
*><br>
*>     If FACT = 'F' and EQUED is not 'N', then AB must have been<br>
*>     equilibrated by the scaling factors in R and/or C.  AB is not<br>
*>     modified if FACT = 'F' or 'N', or if FACT = 'E' and<br>
*>     EQUED = 'N' on exit.<br>
*><br>
*>     On exit, if EQUED .ne. 'N', A is scaled as follows:<br>
*>     EQUED = 'R':  A := diag(R) * A<br>
*>     EQUED = 'C':  A := A * diag(C)<br>
*>     EQUED = 'B':  A := diag(R) * A * diag(C).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>     The leading dimension of the array AB.  LDAB >= KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AFB<br>
*> \verbatim<br>
*>          AFB is DOUBLE PRECISION array, dimension (LDAFB,N)<br>
*>     If FACT = 'F', then AFB is an input argument and on entry<br>
*>     contains details of the LU factorization of the band matrix<br>
*>     A, as computed by DGBTRF.  U is stored as an upper triangular<br>
*>     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,<br>
*>     and the multipliers used during the factorization are stored<br>
*>     in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is<br>
*>     the factored form of the equilibrated matrix A.<br>
*><br>
*>     If FACT = 'N', then AF is an output argument and on exit<br>
*>     returns the factors L and U from the factorization A = P*L*U<br>
*>     of the original matrix A.<br>
*><br>
*>     If FACT = 'E', then AF is an output argument and on exit<br>
*>     returns the factors L and U from the factorization A = P*L*U<br>
*>     of the equilibrated matrix A (see the description of A for<br>
*>     the form of the equilibrated matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     If FACT = 'F', then IPIV is an input argument and on entry<br>
*>     contains the pivot indices from the factorization A = P*L*U<br>
*>     as computed by DGETRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*><br>
*>     If FACT = 'N', then IPIV is an output argument and on exit<br>
*>     contains the pivot indices from the factorization A = P*L*U<br>
*>     of the original matrix A.<br>
*><br>
*>     If FACT = 'E', then IPIV is an output argument and on exit<br>
*>     contains the pivot indices from the factorization A = P*L*U<br>
*>     of the equilibrated matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>     Specifies the form of equilibration that was done.<br>
*>       = 'N':  No equilibration (always true if FACT = 'N').<br>
*>       = 'R':  Row equilibration, i.e., A has been premultiplied by<br>
*>               diag(R).<br>
*>       = 'C':  Column equilibration, i.e., A has been postmultiplied<br>
*>               by diag(C).<br>
*>       = 'B':  Both row and column equilibration, i.e., A has been<br>
*>               replaced by diag(R) * A * diag(C).<br>
*>     EQUED is an input argument if FACT = 'F'; otherwise, it is an<br>
*>     output argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (N)<br>
*>     The row scale factors for A.  If EQUED = 'R' or 'B', A is<br>
*>     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R<br>
*>     is not accessed.  R is an input argument if FACT = 'F';<br>
*>     otherwise, R is an output argument.  If FACT = 'F' and<br>
*>     EQUED = 'R' or 'B', each element of R must be positive.<br>
*>     If R is output, each element of R is a power of the radix.<br>
*>     If R is input, each element of R should be a power of the radix<br>
*>     to ensure a reliable solution and error estimates. Scaling by<br>
*>     powers of the radix does not cause rounding errors unless the<br>
*>     result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>     The column scale factors for A.  If EQUED = 'C' or 'B', A is<br>
*>     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C<br>
*>     is not accessed.  C is an input argument if FACT = 'F';<br>
*>     otherwise, C is an output argument.  If FACT = 'F' and<br>
*>     EQUED = 'C' or 'B', each element of C must be positive.<br>
*>     If C is output, each element of C is a power of the radix.<br>
*>     If C is input, each element of C should be a power of the radix<br>
*>     to ensure a reliable solution and error estimates. Scaling by<br>
*>     powers of the radix does not cause rounding errors unless the<br>
*>     result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>     On entry, the N-by-NRHS right hand side matrix B.<br>
*>     On exit,<br>
*>     if EQUED = 'N', B is not modified;<br>
*>     if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by<br>
*>        diag(R)*B;<br>
*>     if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is<br>
*>        overwritten by diag(C)*B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>     If INFO = 0, the N-by-NRHS solution matrix X to the original<br>
*>     system of equations.  Note that A and B are modified on exit<br>
*>     if EQUED .ne. 'N', and the solution to the equilibrated system is<br>
*>     inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or 'B', or<br>
*>     inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>     The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RPVGRW<br>
*> \verbatim<br>
*>          RPVGRW is DOUBLE PRECISION<br>
*>     Reciprocal pivot growth.  On exit, this contains the reciprocal<br>
*>     pivot growth factor norm(A)/norm(U). The "max absolute element"<br>
*>     norm is used.  If this is much less than 1, then the stability of<br>
*>     the LU factorization of the (equilibrated) matrix A could be poor.<br>
*>     This also means that the solution X, estimated condition numbers,<br>
*>     and error bounds could be unreliable. If factorization fails with<br>
*>     0<INFO<=N, then this contains the reciprocal pivot growth factor<br>
*>     for the leading INFO columns of A.  In DGESVX, this quantity is<br>
*>     returned in WORK(1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>     Componentwise relative backward error.  This is the<br>
*>     componentwise relative backward error of each solution vector X(j)<br>
*>     (i.e., the smallest relative change in any element of A or B that<br>
*>     makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_ERR_BNDS<br>
*> \verbatim<br>
*>          N_ERR_BNDS is INTEGER<br>
*>     Number of error bounds to return for each right hand side<br>
*>     and each type (normwise or componentwise).  See ERR_BNDS_NORM and<br>
*>     ERR_BNDS_COMP below.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_NORM<br>
*> \verbatim<br>
*>          ERR_BNDS_NORM is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>              sqrt(n) * dlamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * dlamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * dlamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_COMP<br>
*> \verbatim<br>
*>          ERR_BNDS_COMP is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>              sqrt(n) * dlamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * dlamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * dlamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NPARAMS<br>
*> \verbatim<br>
*>          NPARAMS is INTEGER<br>
*>     Specifies the number of parameters set in PARAMS.  If .LE. 0, the<br>
*>     PARAMS array is never referenced and default values are used.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] PARAMS<br>
*> \verbatim<br>
*>          PARAMS is DOUBLE PRECISION array, dimension (NPARAMS)<br>
*>     Specifies algorithm parameters.  If an entry is .LT. 0.0, then<br>
*>     that entry will be filled with default value used for that<br>
*>     parameter.  Only positions up to NPARAMS are accessed; defaults<br>
*>     are used for higher-numbered parameters.<br>
*><br>
*>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative<br>
*>            refinement or not.<br>
*>         Default: 1.0D+0<br>
*>            = 0.0 : No refinement is performed, and no error bounds are<br>
*>                    computed.<br>
*>            = 1.0 : Use the extra-precise refinement algorithm.<br>
*>              (other values are reserved for future use)<br>
*><br>
*>       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual<br>
*>            computations allowed for refinement.<br>
*>         Default: 10<br>
*>         Aggressive: Set to 100 to permit convergence using approximate<br>
*>                     factorizations or factorizations other than LU. If<br>
*>                     the factorization uses a technique other than<br>
*>                     Gaussian elimination, the guarantees in<br>
*>                     err_bnds_norm and err_bnds_comp may no longer be<br>
*>                     trustworthy.<br>
*><br>
*>       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code<br>
*>            will attempt to find a solution with small componentwise<br>
*>            relative error in the double-precision algorithm.  Positive<br>
*>            is true, 0.0 is false.<br>
*>         Default: 1.0 (attempt componentwise convergence)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit. The solution to every right-hand side is<br>
*>         guaranteed.<br>
*>       < 0:  If INFO = -i, the i-th argument had an illegal value<br>
*>       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization<br>
*>         has been completed, but the factor U is exactly singular, so<br>
*>         the solution and error bounds could not be computed. RCOND = 0<br>
*>         is returned.<br>
*>       = N+J: The solution corresponding to the Jth right-hand side is<br>
*>         not guaranteed. The solutions corresponding to other right-<br>
*>         hand sides K with K > J may not be guaranteed as well, but<br>
*>         only the first such right-hand side is reported. If a small<br>
*>         componentwise error is not requested (PARAMS(3) = 0.0) then<br>
*>         the Jth right-hand side is the first with a normwise error<br>
*>         bound that is not guaranteed (the smallest J such<br>
*>         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)<br>
*>         the Jth right-hand side is the first with either a normwise or<br>
*>         componentwise error bound that is not guaranteed (the smallest<br>
*>         J such that either ERR_BNDS_NORM(J,1) = 0.0 or<br>
*>         ERR_BNDS_COMP(J,1) = 0.0). See the definition of<br>
*>         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information<br>
*>         about all of the right-hand sides check ERR_BNDS_NORM or<br>
*>         ERR_BNDS_COMP.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGBsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgbsvxx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,CHARACTER EQUED,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE DLA_GBRPVGRW);
/**
*> \brief \b DGBTF2 computes the LU factorization of a general band matrix using the unblocked version of the algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBTF2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtf2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtf2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtf2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, KL, KU, LDAB, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBTF2 computes an LU factorization of a real m-by-n band matrix A<br>
*> using partial pivoting with row interchanges.<br>
*><br>
*> This is the unblocked version of the algorithm, calling Level 2 BLAS.<br>
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
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          On entry, the matrix A in band storage, in rows KL+1 to<br>
*>          2*KL+KU+1; rows 1 to KL of the array need not be set.<br>
*>          The j-th column of A is stored in the j-th column of the<br>
*>          array AB as follows:<br>
*>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)<br>
*><br>
*>          On exit, details of the factorization: U is stored as an<br>
*>          upper triangular band matrix with KL+KU superdiagonals in<br>
*>          rows 1 to KL+KU+1, and the multipliers used during the<br>
*>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.<br>
*>          See below for further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (min(M,N))<br>
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization<br>
*>               has been completed, but the factor U is exactly<br>
*>               singular, and division by zero will occur if it is used<br>
*>               to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGBcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The band storage scheme is illustrated by the following example, when<br>
*>  M = N = 6, KL = 2, KU = 1:<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>      *    *    *    +    +    +       *    *    *   u14  u25  u36<br>
*>      *    *    +    +    +    +       *    *   u13  u24  u35  u46<br>
*>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56<br>
*>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66<br>
*>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *<br>
*>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *<br>
*><br>
*>  Array elements marked * are not used by the routine; elements marked<br>
*>  + need not be set on entry, but are required by the routine to store<br>
*>  elements of U, because of fill-in resulting from the row<br>
*>  interchanges.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgbtf2_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,int[] IPIV,INTEGER INFO);
/**
*> \brief \b DGBTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, KL, KU, LDAB, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBTRF computes an LU factorization of a real m-by-n band matrix A<br>
*> using partial pivoting with row interchanges.<br>
*><br>
*> This is the blocked version of the algorithm, calling Level 3 BLAS.<br>
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
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          On entry, the matrix A in band storage, in rows KL+1 to<br>
*>          2*KL+KU+1; rows 1 to KL of the array need not be set.<br>
*>          The j-th column of A is stored in the j-th column of the<br>
*>          array AB as follows:<br>
*>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)<br>
*><br>
*>          On exit, details of the factorization: U is stored as an<br>
*>          upper triangular band matrix with KL+KU superdiagonals in<br>
*>          rows 1 to KL+KU+1, and the multipliers used during the<br>
*>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.<br>
*>          See below for further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (min(M,N))<br>
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization<br>
*>               has been completed, but the factor U is exactly<br>
*>               singular, and division by zero will occur if it is used<br>
*>               to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGBcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The band storage scheme is illustrated by the following example, when<br>
*>  M = N = 6, KL = 2, KU = 1:<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>      *    *    *    +    +    +       *    *    *   u14  u25  u36<br>
*>      *    *    +    +    +    +       *    *   u13  u24  u35  u46<br>
*>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56<br>
*>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66<br>
*>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *<br>
*>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *<br>
*><br>
*>  Array elements marked * are not used by the routine; elements marked<br>
*>  + need not be set on entry, but are required by the routine to store<br>
*>  elements of U because of fill-in resulting from the row interchanges.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgbtrf_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,int[] IPIV,INTEGER INFO);
/**
*> \brief \b DGBTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGBTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBTRS solves a system of linear equations<br>
*>    A * X = B  or  A**T * X = B<br>
*> with a general band matrix A using the LU factorization computed<br>
*> by DGBTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations.<br>
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T* X = B  (Transpose)<br>
*>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is DOUBLE PRECISION array, dimension (LDAB,N)<br>
*>          Details of the LU factorization of the band matrix A, as<br>
*>          computed by DGBTRF.  U is stored as an upper triangular band<br>
*>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and<br>
*>          the multipliers used during the factorization are stored in<br>
*>          rows KL+KU+2 to 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices; for 1 <= i <= N, row i of the matrix was<br>
*>          interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the right hand side matrix B.<br>
*>          On exit, the solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleGBcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgbtrs_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b DGEBAK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEBAK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebak.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebak.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebak.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOB, SIDE<br>
*       INTEGER            IHI, ILO, INFO, LDV, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   SCALE( * ), V( LDV, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEBAK forms the right or left eigenvectors of a real general matrix<br>
*> by backward transformation on the computed eigenvectors of the<br>
*> balanced matrix output by DGEBAL.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          Specifies the type of backward transformation required:<br>
*>          = 'N', do nothing, return immediately;<br>
*>          = 'P', do backward transformation for permutation only;<br>
*>          = 'S', do backward transformation for scaling only;<br>
*>          = 'B', do backward transformations for both permutation and<br>
*>                 scaling.<br>
*>          JOB must be the same as the argument JOB supplied to DGEBAL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'R':  V contains right eigenvectors;<br>
*>          = 'L':  V contains left eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of rows of the matrix V.  N >= 0.<br>
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
*>          The integers ILO and IHI determined by DGEBAL.<br>
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SCALE<br>
*> \verbatim<br>
*>          SCALE is DOUBLE PRECISION array, dimension (N)<br>
*>          Details of the permutation and scaling factors, as returned<br>
*>          by DGEBAL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of columns of the matrix V.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension (LDV,M)<br>
*>          On entry, the matrix of right or left eigenvectors to be<br>
*>          transformed, as returned by DHSEIN or DTREVC.<br>
*>          On exit, V is overwritten by the transformed eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V. LDV >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgebak_(CHARACTER JOB,CHARACTER SIDE,INTEGER N,INTEGER ILO,INTEGER IHI,double[] SCALE,INTEGER M,double[] V,INTEGER LDV,INTEGER INFO);
/**
*> \brief \b DGEBAL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEBAL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebal.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebal.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebal.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOB<br>
*       INTEGER            IHI, ILO, INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), SCALE( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEBAL balances a general real matrix A.  This involves, first,<br>
*> permuting A by a similarity transformation to isolate eigenvalues<br>
*> in the first 1 to ILO-1 and last IHI+1 to N elements on the<br>
*> diagonal; and second, applying a diagonal similarity transformation<br>
*> to rows and columns ILO to IHI to make the rows and columns as<br>
*> close in norm as possible.  Both steps are optional.<br>
*><br>
*> Balancing may reduce the 1-norm of the matrix, and improve the<br>
*> accuracy of the computed eigenvalues and/or eigenvectors.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          Specifies the operations to be performed on A:<br>
*>          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0<br>
*>                  for i = 1,...,N;<br>
*>          = 'P':  permute only;<br>
*>          = 'S':  scale only;<br>
*>          = 'B':  both permute and scale.<br>
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
*>          A is DOUBLE array, dimension (LDA,N)<br>
*>          On entry, the input matrix A.<br>
*>          On exit,  A is overwritten by the balanced matrix.<br>
*>          If JOB = 'N', A is not referenced.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ILO<br>
*> \verbatim<br>
*>          ILO is INTEGER<br>
*> \endverbatim<br>
*> \param[out] IHI<br>
*> \verbatim<br>
*>          IHI is INTEGER<br>
*>          ILO and IHI are set to integers such that on exit<br>
*>          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.<br>
*>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is DOUBLE array, dimension (N)<br>
*>          Details of the permutations and scaling factors applied to<br>
*>          A.  If P(j) is the index of the row and column interchanged<br>
*>          with row and column j and D(j) is the scaling factor<br>
*>          applied to row and column j, then<br>
*>          SCALE(j) = P(j)    for j = 1,...,ILO-1<br>
*>                   = D(j)    for j = ILO,...,IHI<br>
*>                   = P(j)    for j = IHI+1,...,N.<br>
*>          The order in which the interchanges are made is N to IHI+1,<br>
*>          then 1 to ILO-1.<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The permutations consist of row and column interchanges which put<br>
*>  the matrix in the form<br>
*><br>
*>             ( T1   X   Y  )<br>
*>     P A P = (  0   B   Z  )<br>
*>             (  0   0   T2 )<br>
*><br>
*>  where T1 and T2 are upper triangular matrices whose eigenvalues lie<br>
*>  along the diagonal.  The column indices ILO and IHI mark the starting<br>
*>  and ending columns of the submatrix B. Balancing consists of applying<br>
*>  a diagonal similarity transformation inv(D) * B * D to make the<br>
*>  1-norms of each row of B and its corresponding column nearly equal.<br>
*>  The output matrix is<br>
*><br>
*>     ( T1     X*D          Y    )<br>
*>     (  0  inv(D)*B*D  inv(D)*Z ).<br>
*>     (  0      0           T2   )<br>
*><br>
*>  Information about the permutations P and the diagonal matrix D is<br>
*>  returned in the vector SCALE.<br>
*><br>
*>  This subroutine is based on the EISPACK routine BALANC.<br>
*><br>
*>  Modified by Tzu-Yi Chen, Computer Science Division, University of<br>
*>    California at Berkeley, USA<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgebal_(CHARACTER JOB,INTEGER N,double[] A,INTEGER LDA,INTEGER ILO,INTEGER IHI,double[] SCALE,INTEGER INFO);
/**
*> \brief \b DGEBD2 reduces a general matrix to bidiagonal form using an unblocked algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEBD2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebd2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebd2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebd2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),<br>
*      $                   TAUQ( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEBD2 reduces a real general m by n matrix A to upper or lower<br>
*> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.<br>
*><br>
*> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows in the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns in the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the m by n general matrix to be reduced.<br>
*>          On exit,<br>
*>          if m >= n, the diagonal and the first superdiagonal are<br>
*>            overwritten with the upper bidiagonal matrix B; the<br>
*>            elements below the diagonal, with the array TAUQ, represent<br>
*>            the orthogonal matrix Q as a product of elementary<br>
*>            reflectors, and the elements above the first superdiagonal,<br>
*>            with the array TAUP, represent the orthogonal matrix P as<br>
*>            a product of elementary reflectors;<br>
*>          if m < n, the diagonal and the first subdiagonal are<br>
*>            overwritten with the lower bidiagonal matrix B; the<br>
*>            elements below the first subdiagonal, with the array TAUQ,<br>
*>            represent the orthogonal matrix Q as a product of<br>
*>            elementary reflectors, and the elements above the diagonal,<br>
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
*>          D is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The diagonal elements of the bidiagonal matrix B:<br>
*>          D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (min(M,N)-1)<br>
*>          The off-diagonal elements of the bidiagonal matrix B:<br>
*>          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;<br>
*>          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ<br>
*> \verbatim<br>
*>          TAUQ is DOUBLE PRECISION array dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix Q. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP<br>
*> \verbatim<br>
*>          TAUP is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix P. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (max(M,N))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit.<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrices Q and P are represented as products of elementary<br>
*>  reflectors:<br>
*><br>
*>  If m >= n,<br>
*><br>
*>     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)<br>
*><br>
*>  Each H(i) and G(i) has the form:<br>
*><br>
*>     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T<br>
*><br>
*>  where tauq and taup are real scalars, and v and u are real vectors;<br>
*>  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);<br>
*>  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);<br>
*>  tauq is stored in TAUQ(i) and taup in TAUP(i).<br>
*><br>
*>  If m < n,<br>
*><br>
*>     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)<br>
*><br>
*>  Each H(i) and G(i) has the form:<br>
*><br>
*>     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T<br>
*><br>
*>  where tauq and taup are real scalars, and v and u are real vectors;<br>
*>  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);<br>
*>  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);<br>
*>  tauq is stored in TAUQ(i) and taup in TAUP(i).<br>
*><br>
*>  The contents of A on exit are illustrated by the following examples:<br>
*><br>
*>  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):<br>
*><br>
*>    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )<br>
*>    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )<br>
*>    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )<br>
*>    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )<br>
*>    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )<br>
*>    (  v1  v2  v3  v4  v5 )<br>
*><br>
*>  where d and e denote diagonal and off-diagonal elements of B, vi<br>
*>  denotes an element of the vector defining H(i), and ui an element of<br>
*>  the vector defining G(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgebd2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] D,double[] E,double[] TAUQ,double[] TAUP,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGEBRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEBRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),<br>
*      $                   TAUQ( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEBRD reduces a general real M-by-N matrix A to upper or lower<br>
*> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.<br>
*><br>
*> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows in the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns in the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N general matrix to be reduced.<br>
*>          On exit,<br>
*>          if m >= n, the diagonal and the first superdiagonal are<br>
*>            overwritten with the upper bidiagonal matrix B; the<br>
*>            elements below the diagonal, with the array TAUQ, represent<br>
*>            the orthogonal matrix Q as a product of elementary<br>
*>            reflectors, and the elements above the first superdiagonal,<br>
*>            with the array TAUP, represent the orthogonal matrix P as<br>
*>            a product of elementary reflectors;<br>
*>          if m < n, the diagonal and the first subdiagonal are<br>
*>            overwritten with the lower bidiagonal matrix B; the<br>
*>            elements below the first subdiagonal, with the array TAUQ,<br>
*>            represent the orthogonal matrix Q as a product of<br>
*>            elementary reflectors, and the elements above the diagonal,<br>
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
*>          D is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The diagonal elements of the bidiagonal matrix B:<br>
*>          D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (min(M,N)-1)<br>
*>          The off-diagonal elements of the bidiagonal matrix B:<br>
*>          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;<br>
*>          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUQ<br>
*> \verbatim<br>
*>          TAUQ is DOUBLE PRECISION array dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix Q. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUP<br>
*> \verbatim<br>
*>          TAUP is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix P. See Further Details.<br>
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
*>          The length of the array WORK.  LWORK >= max(1,M,N).<br>
*>          For optimum performance LWORK >= (M+N)*NB, where NB<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrices Q and P are represented as products of elementary<br>
*>  reflectors:<br>
*><br>
*>  If m >= n,<br>
*><br>
*>     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)<br>
*><br>
*>  Each H(i) and G(i) has the form:<br>
*><br>
*>     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T<br>
*><br>
*>  where tauq and taup are real scalars, and v and u are real vectors;<br>
*>  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);<br>
*>  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);<br>
*>  tauq is stored in TAUQ(i) and taup in TAUP(i).<br>
*><br>
*>  If m < n,<br>
*><br>
*>     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)<br>
*><br>
*>  Each H(i) and G(i) has the form:<br>
*><br>
*>     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T<br>
*><br>
*>  where tauq and taup are real scalars, and v and u are real vectors;<br>
*>  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);<br>
*>  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);<br>
*>  tauq is stored in TAUQ(i) and taup in TAUP(i).<br>
*><br>
*>  The contents of A on exit are illustrated by the following examples:<br>
*><br>
*>  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):<br>
*><br>
*>    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )<br>
*>    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )<br>
*>    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )<br>
*>    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )<br>
*>    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )<br>
*>    (  v1  v2  v3  v4  v5 )<br>
*><br>
*>  where d and e denote diagonal and off-diagonal elements of B, vi<br>
*>  denotes an element of the vector defining H(i), and ui an element of<br>
*>  the vector defining G(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgebrd_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] D,double[] E,double[] TAUQ,double[] TAUP,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGECON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGECON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgecon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgecon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgecon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            INFO, LDA, N<br>
*       DOUBLE PRECISION   ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGECON estimates the reciprocal of the condition number of a general<br>
*> real matrix A, in either the 1-norm or the infinity-norm, using<br>
*> the LU factorization computed by DGETRF.<br>
*><br>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the<br>
*> condition number is computed as<br>
*>    RCOND = 1 / ( norm(A) * norm(inv(A)) ).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies whether the 1-norm condition number or the<br>
*>          infinity-norm condition number is required:<br>
*>          = '1' or 'O':  1-norm;<br>
*>          = 'I':         Infinity-norm.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          The factors L and U from the factorization A = P*L*U<br>
*>          as computed by DGETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is DOUBLE PRECISION<br>
*>          If NORM = '1' or 'O', the 1-norm of the original matrix A.<br>
*>          If NORM = 'I', the infinity-norm of the original matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(norm(A) * norm(inv(A))).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgecon_(CHARACTER NORM,INTEGER N,double[] A,INTEGER LDA,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGEEQU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEEQU + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeequ.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeequ.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeequ.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
*       DOUBLE PRECISION   AMAX, COLCND, ROWCND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEEQU computes row and column scalings intended to equilibrate an<br>
*> M-by-N matrix A and reduce its condition number.  R returns the row<br>
*> scale factors and C the column scale factors, chosen to try to make<br>
*> the largest element in each row and column of the matrix B with<br>
*> elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.<br>
*><br>
*> R(i) and C(j) are restricted to be between SMLNUM = smallest safe<br>
*> number and BIGNUM = largest safe number.  Use of these scaling<br>
*> factors is not guaranteed to reduce the condition number of A but<br>
*> works well in practice.<br>
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
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          The M-by-N matrix whose equilibration factors are<br>
*>          to be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (M)<br>
*>          If INFO = 0 or INFO > M, R contains the row scale factors<br>
*>          for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0,  C contains the column scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ROWCND<br>
*> \verbatim<br>
*>          ROWCND is DOUBLE PRECISION<br>
*>          If INFO = 0 or INFO > M, ROWCND contains the ratio of the<br>
*>          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and<br>
*>          AMAX is neither too large nor too small, it is not worth<br>
*>          scaling by R.<br>
*> \endverbatim<br>
*><br>
*> \param[out] COLCND<br>
*> \verbatim<br>
*>          COLCND is DOUBLE PRECISION<br>
*>          If INFO = 0, COLCND contains the ratio of the smallest<br>
*>          C(i) to the largest C(i).  If COLCND >= 0.1, it is not<br>
*>          worth scaling by C.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AMAX<br>
*> \verbatim<br>
*>          AMAX is DOUBLE PRECISION<br>
*>          Absolute value of largest matrix element.  If AMAX is very<br>
*>          close to overflow or very close to underflow, the matrix<br>
*>          should be scaled.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i,  and i is<br>
*>                <= M:  the i-th row of A is exactly zero<br>
*>                >  M:  the (i-M)-th column of A is exactly zero<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgeequ_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,INTEGER INFO);
/**
*> \brief \b DGEEQUB<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEEQUB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeequb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeequb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeequb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEEQUB( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,<br>
*                           INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
*       DOUBLE PRECISION   AMAX, COLCND, ROWCND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEEQUB computes row and column scalings intended to equilibrate an<br>
*> M-by-N matrix A and reduce its condition number.  R returns the row<br>
*> scale factors and C the column scale factors, chosen to try to make<br>
*> the largest element in each row and column of the matrix B with<br>
*> elements B(i,j)=R(i)*A(i,j)*C(j) have an absolute value of at most<br>
*> the radix.<br>
*><br>
*> R(i) and C(j) are restricted to be a power of the radix between<br>
*> SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use<br>
*> of these scaling factors is not guaranteed to reduce the condition<br>
*> number of A but works well in practice.<br>
*><br>
*> This routine differs from DGEEQU by restricting the scaling factors<br>
*> to a power of the radix.  Baring over- and underflow, scaling by<br>
*> these factors introduces no additional rounding errors.  However, the<br>
*> scaled entries' magnitured are no longer approximately 1 but lie<br>
*> between sqrt(radix) and 1/sqrt(radix).<br>
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
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          The M-by-N matrix whose equilibration factors are<br>
*>          to be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (M)<br>
*>          If INFO = 0 or INFO > M, R contains the row scale factors<br>
*>          for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0,  C contains the column scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ROWCND<br>
*> \verbatim<br>
*>          ROWCND is DOUBLE PRECISION<br>
*>          If INFO = 0 or INFO > M, ROWCND contains the ratio of the<br>
*>          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and<br>
*>          AMAX is neither too large nor too small, it is not worth<br>
*>          scaling by R.<br>
*> \endverbatim<br>
*><br>
*> \param[out] COLCND<br>
*> \verbatim<br>
*>          COLCND is DOUBLE PRECISION<br>
*>          If INFO = 0, COLCND contains the ratio of the smallest<br>
*>          C(i) to the largest C(i).  If COLCND >= 0.1, it is not<br>
*>          worth scaling by C.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AMAX<br>
*> \verbatim<br>
*>          AMAX is DOUBLE PRECISION<br>
*>          Absolute value of largest matrix element.  If AMAX is very<br>
*>          close to overflow or very close to underflow, the matrix<br>
*>          should be scaled.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i,  and i is<br>
*>                <= M:  the i-th row of A is exactly zero<br>
*>                >  M:  the (i-M)-th column of A is exactly zero<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgeequb_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,INTEGER INFO);
/**
*> \brief <b> DGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEES + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgees.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgees.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgees.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,<br>
*                         VS, LDVS, WORK, LWORK, BWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBVS, SORT<br>
*       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            BWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),<br>
*      $                   WR( * )<br>
*       ..<br>
*       .. Function Arguments ..<br>
*       LOGICAL            SELECT<br>
*       EXTERNAL           SELECT<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEES computes for an N-by-N real nonsymmetric matrix A, the<br>
*> eigenvalues, the real Schur form T, and, optionally, the matrix of<br>
*> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).<br>
*><br>
*> Optionally, it also orders the eigenvalues on the diagonal of the<br>
*> real Schur form so that selected eigenvalues are at the top left.<br>
*> The leading columns of Z then form an orthonormal basis for the<br>
*> invariant subspace corresponding to the selected eigenvalues.<br>
*><br>
*> A matrix is in real Schur form if it is upper quasi-triangular with<br>
*> 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the<br>
*> form<br>
*>         [  a  b  ]<br>
*>         [  c  a  ]<br>
*><br>
*> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBVS<br>
*> \verbatim<br>
*>          JOBVS is CHARACTER*1<br>
*>          = 'N': Schur vectors are not computed;<br>
*>          = 'V': Schur vectors are computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SORT<br>
*> \verbatim<br>
*>          SORT is CHARACTER*1<br>
*>          Specifies whether or not to order the eigenvalues on the<br>
*>          diagonal of the Schur form.<br>
*>          = 'N': Eigenvalues are not ordered;<br>
*>          = 'S': Eigenvalues are ordered (see SELECT).<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is a LOGICAL FUNCTION of two DOUBLE PRECISION arguments<br>
*>          SELECT must be declared EXTERNAL in the calling subroutine.<br>
*>          If SORT = 'S', SELECT is used to select eigenvalues to sort<br>
*>          to the top left of the Schur form.<br>
*>          If SORT = 'N', SELECT is not referenced.<br>
*>          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if<br>
*>          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex<br>
*>          conjugate pair of eigenvalues is selected, then both complex<br>
*>          eigenvalues are selected.<br>
*>          Note that a selected complex eigenvalue may no longer<br>
*>          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since<br>
*>          ordering may change the value of complex eigenvalues<br>
*>          (especially if the eigenvalue is ill-conditioned); in this<br>
*>          case INFO is set to N+2 (see INFO below).<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the N-by-N matrix A.<br>
*>          On exit, A has been overwritten by its real Schur form T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SDIM<br>
*> \verbatim<br>
*>          SDIM is INTEGER<br>
*>          If SORT = 'N', SDIM = 0.<br>
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)<br>
*>                         for which SELECT is true. (Complex conjugate<br>
*>                         pairs for which SELECT is true for either<br>
*>                         eigenvalue count as 2.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is DOUBLE PRECISION array, dimension (N)<br>
*>          WR and WI contain the real and imaginary parts,<br>
*>          respectively, of the computed eigenvalues in the same order<br>
*>          that they appear on the diagonal of the output Schur form T.<br>
*>          Complex conjugate pairs of eigenvalues will appear<br>
*>          consecutively with the eigenvalue having the positive<br>
*>          imaginary part first.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VS<br>
*> \verbatim<br>
*>          VS is DOUBLE PRECISION array, dimension (LDVS,N)<br>
*>          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur<br>
*>          vectors.<br>
*>          If JOBVS = 'N', VS is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVS<br>
*> \verbatim<br>
*>          LDVS is INTEGER<br>
*>          The leading dimension of the array VS.  LDVS >= 1; if<br>
*>          JOBVS = 'V', LDVS >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.  LWORK >= max(1,3*N).<br>
*>          For good performance, LWORK must generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BWORK<br>
*> \verbatim<br>
*>          BWORK is LOGICAL array, dimension (N)<br>
*>          Not referenced if SORT = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0: if INFO = i, and i is<br>
*>             <= N: the QR algorithm failed to compute all the<br>
*>                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI<br>
*>                   contain those eigenvalues which have converged; if<br>
*>                   JOBVS = 'V', VS contains the matrix which reduces A<br>
*>                   to its partially converged Schur form.<br>
*>             = N+1: the eigenvalues could not be reordered because some<br>
*>                   eigenvalues were too close to separate (the problem<br>
*>                   is very ill-conditioned);<br>
*>             = N+2: after reordering, roundoff changed values of some<br>
*>                   complex eigenvalues so that leading eigenvalues in<br>
*>                   the Schur form no longer satisfy SELECT=.TRUE.  This<br>
*>                   could also be caused by underflow due to scaling.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dgees_(CHARACTER JOBVS,CHARACTER SORT,LOGICAL SELECT,INTEGER N,double[] A,INTEGER LDA,INTEGER SDIM,double[] WR,double[] WI,double[] VS,INTEGER LDVS,double[] WORK,INTEGER LWORK,boolean[] BWORK,INTEGER INFO);
/**
*> \brief <b> DGEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEESX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeesx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeesx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeesx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM,<br>
*                          WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK,<br>
*                          IWORK, LIWORK, BWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBVS, SENSE, SORT<br>
*       INTEGER            INFO, LDA, LDVS, LIWORK, LWORK, N, SDIM<br>
*       DOUBLE PRECISION   RCONDE, RCONDV<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            BWORK( * )<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),<br>
*      $                   WR( * )<br>
*       ..<br>
*       .. Function Arguments ..<br>
*       LOGICAL            SELECT<br>
*       EXTERNAL           SELECT<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEESX computes for an N-by-N real nonsymmetric matrix A, the<br>
*> eigenvalues, the real Schur form T, and, optionally, the matrix of<br>
*> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).<br>
*><br>
*> Optionally, it also orders the eigenvalues on the diagonal of the<br>
*> real Schur form so that selected eigenvalues are at the top left;<br>
*> computes a reciprocal condition number for the average of the<br>
*> selected eigenvalues (RCONDE); and computes a reciprocal condition<br>
*> number for the right invariant subspace corresponding to the<br>
*> selected eigenvalues (RCONDV).  The leading columns of Z form an<br>
*> orthonormal basis for this invariant subspace.<br>
*><br>
*> For further explanation of the reciprocal condition numbers RCONDE<br>
*> and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where<br>
*> these quantities are called s and sep respectively).<br>
*><br>
*> A real matrix is in real Schur form if it is upper quasi-triangular<br>
*> with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in<br>
*> the form<br>
*>           [  a  b  ]<br>
*>           [  c  a  ]<br>
*><br>
*> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBVS<br>
*> \verbatim<br>
*>          JOBVS is CHARACTER*1<br>
*>          = 'N': Schur vectors are not computed;<br>
*>          = 'V': Schur vectors are computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SORT<br>
*> \verbatim<br>
*>          SORT is CHARACTER*1<br>
*>          Specifies whether or not to order the eigenvalues on the<br>
*>          diagonal of the Schur form.<br>
*>          = 'N': Eigenvalues are not ordered;<br>
*>          = 'S': Eigenvalues are ordered (see SELECT).<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is a LOGICAL FUNCTION of two DOUBLE PRECISION arguments<br>
*>          SELECT must be declared EXTERNAL in the calling subroutine.<br>
*>          If SORT = 'S', SELECT is used to select eigenvalues to sort<br>
*>          to the top left of the Schur form.<br>
*>          If SORT = 'N', SELECT is not referenced.<br>
*>          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if<br>
*>          SELECT(WR(j),WI(j)) is true; i.e., if either one of a<br>
*>          complex conjugate pair of eigenvalues is selected, then both<br>
*>          are.  Note that a selected complex eigenvalue may no longer<br>
*>          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since<br>
*>          ordering may change the value of complex eigenvalues<br>
*>          (especially if the eigenvalue is ill-conditioned); in this<br>
*>          case INFO may be set to N+3 (see INFO below).<br>
*> \endverbatim<br>
*><br>
*> \param[in] SENSE<br>
*> \verbatim<br>
*>          SENSE is CHARACTER*1<br>
*>          Determines which reciprocal condition numbers are computed.<br>
*>          = 'N': None are computed;<br>
*>          = 'E': Computed for average of selected eigenvalues only;<br>
*>          = 'V': Computed for selected right invariant subspace only;<br>
*>          = 'B': Computed for both.<br>
*>          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the N-by-N matrix A.<br>
*>          On exit, A is overwritten by its real Schur form T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SDIM<br>
*> \verbatim<br>
*>          SDIM is INTEGER<br>
*>          If SORT = 'N', SDIM = 0.<br>
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)<br>
*>                         for which SELECT is true. (Complex conjugate<br>
*>                         pairs for which SELECT is true for either<br>
*>                         eigenvalue count as 2.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is DOUBLE PRECISION array, dimension (N)<br>
*>          WR and WI contain the real and imaginary parts, respectively,<br>
*>          of the computed eigenvalues, in the same order that they<br>
*>          appear on the diagonal of the output Schur form T.  Complex<br>
*>          conjugate pairs of eigenvalues appear consecutively with the<br>
*>          eigenvalue having the positive imaginary part first.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VS<br>
*> \verbatim<br>
*>          VS is DOUBLE PRECISION array, dimension (LDVS,N)<br>
*>          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur<br>
*>          vectors.<br>
*>          If JOBVS = 'N', VS is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVS<br>
*> \verbatim<br>
*>          LDVS is INTEGER<br>
*>          The leading dimension of the array VS.  LDVS >= 1, and if<br>
*>          JOBVS = 'V', LDVS >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCONDE<br>
*> \verbatim<br>
*>          RCONDE is DOUBLE PRECISION<br>
*>          If SENSE = 'E' or 'B', RCONDE contains the reciprocal<br>
*>          condition number for the average of the selected eigenvalues.<br>
*>          Not referenced if SENSE = 'N' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCONDV<br>
*> \verbatim<br>
*>          RCONDV is DOUBLE PRECISION<br>
*>          If SENSE = 'V' or 'B', RCONDV contains the reciprocal<br>
*>          condition number for the selected right invariant subspace.<br>
*>          Not referenced if SENSE = 'N' or 'E'.<br>
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
*>          The dimension of the array WORK.  LWORK >= max(1,3*N).<br>
*>          Also, if SENSE = 'E' or 'V' or 'B',<br>
*>          LWORK >= N+2*SDIM*(N-SDIM), where SDIM is the number of<br>
*>          selected eigenvalues computed by this routine.  Note that<br>
*>          N+2*SDIM*(N-SDIM) <= N+N*N/2. Note also that an error is only<br>
*>          returned if LWORK < max(1,3*N), but if SENSE = 'E' or 'V' or<br>
*>          'B' this may not be large enough.<br>
*>          For good performance, LWORK must generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates upper bounds on the optimal sizes of the<br>
*>          arrays WORK and IWORK, returns these values as the first<br>
*>          entries of the WORK and IWORK arrays, and no error messages<br>
*>          related to LWORK or LIWORK are issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))<br>
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK.<br>
*>          LIWORK >= 1; if SENSE = 'V' or 'B', LIWORK >= SDIM*(N-SDIM).<br>
*>          Note that SDIM*(N-SDIM) <= N*N/4. Note also that an error is<br>
*>          only returned if LIWORK < 1, but if SENSE = 'V' or 'B' this<br>
*>          may not be large enough.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates upper bounds on the optimal sizes of<br>
*>          the arrays WORK and IWORK, returns these values as the first<br>
*>          entries of the WORK and IWORK arrays, and no error messages<br>
*>          related to LWORK or LIWORK are issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BWORK<br>
*> \verbatim<br>
*>          BWORK is LOGICAL array, dimension (N)<br>
*>          Not referenced if SORT = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0: if INFO = i, and i is<br>
*>             <= N: the QR algorithm failed to compute all the<br>
*>                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI<br>
*>                   contain those eigenvalues which have converged; if<br>
*>                   JOBVS = 'V', VS contains the transformation which<br>
*>                   reduces A to its partially converged Schur form.<br>
*>             = N+1: the eigenvalues could not be reordered because some<br>
*>                   eigenvalues were too close to separate (the problem<br>
*>                   is very ill-conditioned);<br>
*>             = N+2: after reordering, roundoff changed values of some<br>
*>                   complex eigenvalues so that leading eigenvalues in<br>
*>                   the Schur form no longer satisfy SELECT=.TRUE.  This<br>
*>                   could also be caused by underflow due to scaling.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dgeesx_(CHARACTER JOBVS,CHARACTER SORT,LOGICAL SELECT,CHARACTER SENSE,INTEGER N,double[] A,INTEGER LDA,INTEGER SDIM,double[] WR,double[] WI,double[] VS,INTEGER LDVS,DOUBLE RCONDE,DOUBLE RCONDV,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,boolean[] BWORK,INTEGER INFO);
/**
*> \brief <b> DGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,<br>
*                         LDVR, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBVL, JOBVR<br>
*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   WI( * ), WORK( * ), WR( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEEV computes for an N-by-N real nonsymmetric matrix A, the<br>
*> eigenvalues and, optionally, the left and/or right eigenvectors.<br>
*><br>
*> The right eigenvector v(j) of A satisfies<br>
*>                  A * v(j) = lambda(j) * v(j)<br>
*> where lambda(j) is its eigenvalue.<br>
*> The left eigenvector u(j) of A satisfies<br>
*>               u(j)**H * A = lambda(j) * u(j)**H<br>
*> where u(j)**H denotes the conjugate-transpose of u(j).<br>
*><br>
*> The computed eigenvectors are normalized to have Euclidean norm<br>
*> equal to 1 and largest component real.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBVL<br>
*> \verbatim<br>
*>          JOBVL is CHARACTER*1<br>
*>          = 'N': left eigenvectors of A are not computed;<br>
*>          = 'V': left eigenvectors of A are computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVR<br>
*> \verbatim<br>
*>          JOBVR is CHARACTER*1<br>
*>          = 'N': right eigenvectors of A are not computed;<br>
*>          = 'V': right eigenvectors of A are computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the N-by-N matrix A.<br>
*>          On exit, A has been overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is DOUBLE PRECISION array, dimension (N)<br>
*>          WR and WI contain the real and imaginary parts,<br>
*>          respectively, of the computed eigenvalues.  Complex<br>
*>          conjugate pairs of eigenvalues appear consecutively<br>
*>          with the eigenvalue having the positive imaginary part<br>
*>          first.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION array, dimension (LDVL,N)<br>
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one<br>
*>          after another in the columns of VL, in the same order<br>
*>          as their eigenvalues.<br>
*>          If JOBVL = 'N', VL is not referenced.<br>
*>          If the j-th eigenvalue is real, then u(j) = VL(:,j),<br>
*>          the j-th column of VL.<br>
*>          If the j-th and (j+1)-st eigenvalues form a complex<br>
*>          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and<br>
*>          u(j+1) = VL(:,j) - i*VL(:,j+1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL.  LDVL >= 1; if<br>
*>          JOBVL = 'V', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VR<br>
*> \verbatim<br>
*>          VR is DOUBLE PRECISION array, dimension (LDVR,N)<br>
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one<br>
*>          after another in the columns of VR, in the same order<br>
*>          as their eigenvalues.<br>
*>          If JOBVR = 'N', VR is not referenced.<br>
*>          If the j-th eigenvalue is real, then v(j) = VR(:,j),<br>
*>          the j-th column of VR.<br>
*>          If the j-th and (j+1)-st eigenvalues form a complex<br>
*>          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and<br>
*>          v(j+1) = VR(:,j) - i*VR(:,j+1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR.  LDVR >= 1; if<br>
*>          JOBVR = 'V', LDVR >= N.<br>
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
*>          The dimension of the array WORK.  LWORK >= max(1,3*N), and<br>
*>          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good<br>
*>          performance, LWORK must generally be larger.<br>
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
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = i, the QR algorithm failed to compute all the<br>
*>                eigenvalues, and no eigenvectors have been computed;<br>
*>                elements i+1:N of WR and WI contain eigenvalues which<br>
*>                have converged.<br>
*> \endverbatim<br>
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
*  @precisions fortran d -> s<br>
*<br>
*> \ingroup doubleGEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dgeev_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,double[] A,INTEGER LDA,double[] WR,double[] WI,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,<br>
*                          VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,<br>
*                          RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE<br>
*       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N<br>
*       DOUBLE PRECISION   ABNRM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), RCONDE( * ), RCONDV( * ),<br>
*      $                   SCALE( * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   WI( * ), WORK( * ), WR( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEEVX computes for an N-by-N real nonsymmetric matrix A, the<br>
*> eigenvalues and, optionally, the left and/or right eigenvectors.<br>
*><br>
*> Optionally also, it computes a balancing transformation to improve<br>
*> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,<br>
*> SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues<br>
*> (RCONDE), and reciprocal condition numbers for the right<br>
*> eigenvectors (RCONDV).<br>
*><br>
*> The right eigenvector v(j) of A satisfies<br>
*>                  A * v(j) = lambda(j) * v(j)<br>
*> where lambda(j) is its eigenvalue.<br>
*> The left eigenvector u(j) of A satisfies<br>
*>               u(j)**H * A = lambda(j) * u(j)**H<br>
*> where u(j)**H denotes the conjugate-transpose of u(j).<br>
*><br>
*> The computed eigenvectors are normalized to have Euclidean norm<br>
*> equal to 1 and largest component real.<br>
*><br>
*> Balancing a matrix means permuting the rows and columns to make it<br>
*> more nearly upper triangular, and applying a diagonal similarity<br>
*> transformation D * A * D**(-1), where D is a diagonal matrix, to<br>
*> make its rows and columns closer in norm and the condition numbers<br>
*> of its eigenvalues and eigenvectors smaller.  The computed<br>
*> reciprocal condition numbers correspond to the balanced matrix.<br>
*> Permuting rows and columns will not change the condition numbers<br>
*> (in exact arithmetic) but diagonal scaling will.  For further<br>
*> explanation of balancing, see section 4.10.2 of the LAPACK<br>
*> Users' Guide.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] BALANC<br>
*> \verbatim<br>
*>          BALANC is CHARACTER*1<br>
*>          Indicates how the input matrix should be diagonally scaled<br>
*>          and/or permuted to improve the conditioning of its<br>
*>          eigenvalues.<br>
*>          = 'N': Do not diagonally scale or permute;<br>
*>          = 'P': Perform permutations to make the matrix more nearly<br>
*>                 upper triangular. Do not diagonally scale;<br>
*>          = 'S': Diagonally scale the matrix, i.e. replace A by<br>
*>                 D*A*D**(-1), where D is a diagonal matrix chosen<br>
*>                 to make the rows and columns of A more equal in<br>
*>                 norm. Do not permute;<br>
*>          = 'B': Both diagonally scale and permute A.<br>
*><br>
*>          Computed reciprocal condition numbers will be for the matrix<br>
*>          after balancing and/or permuting. Permuting does not change<br>
*>          condition numbers (in exact arithmetic), but balancing does.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVL<br>
*> \verbatim<br>
*>          JOBVL is CHARACTER*1<br>
*>          = 'N': left eigenvectors of A are not computed;<br>
*>          = 'V': left eigenvectors of A are computed.<br>
*>          If SENSE = 'E' or 'B', JOBVL must = 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVR<br>
*> \verbatim<br>
*>          JOBVR is CHARACTER*1<br>
*>          = 'N': right eigenvectors of A are not computed;<br>
*>          = 'V': right eigenvectors of A are computed.<br>
*>          If SENSE = 'E' or 'B', JOBVR must = 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SENSE<br>
*> \verbatim<br>
*>          SENSE is CHARACTER*1<br>
*>          Determines which reciprocal condition numbers are computed.<br>
*>          = 'N': None are computed;<br>
*>          = 'E': Computed for eigenvalues only;<br>
*>          = 'V': Computed for right eigenvectors only;<br>
*>          = 'B': Computed for eigenvalues and right eigenvectors.<br>
*><br>
*>          If SENSE = 'E' or 'B', both left and right eigenvectors<br>
*>          must also be computed (JOBVL = 'V' and JOBVR = 'V').<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the N-by-N matrix A.<br>
*>          On exit, A has been overwritten.  If JOBVL = 'V' or<br>
*>          JOBVR = 'V', A contains the real Schur form of the balanced<br>
*>          version of the input matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is DOUBLE PRECISION array, dimension (N)<br>
*>          WR and WI contain the real and imaginary parts,<br>
*>          respectively, of the computed eigenvalues.  Complex<br>
*>          conjugate pairs of eigenvalues will appear consecutively<br>
*>          with the eigenvalue having the positive imaginary part<br>
*>          first.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION array, dimension (LDVL,N)<br>
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one<br>
*>          after another in the columns of VL, in the same order<br>
*>          as their eigenvalues.<br>
*>          If JOBVL = 'N', VL is not referenced.<br>
*>          If the j-th eigenvalue is real, then u(j) = VL(:,j),<br>
*>          the j-th column of VL.<br>
*>          If the j-th and (j+1)-st eigenvalues form a complex<br>
*>          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and<br>
*>          u(j+1) = VL(:,j) - i*VL(:,j+1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL.  LDVL >= 1; if<br>
*>          JOBVL = 'V', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VR<br>
*> \verbatim<br>
*>          VR is DOUBLE PRECISION array, dimension (LDVR,N)<br>
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one<br>
*>          after another in the columns of VR, in the same order<br>
*>          as their eigenvalues.<br>
*>          If JOBVR = 'N', VR is not referenced.<br>
*>          If the j-th eigenvalue is real, then v(j) = VR(:,j),<br>
*>          the j-th column of VR.<br>
*>          If the j-th and (j+1)-st eigenvalues form a complex<br>
*>          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and<br>
*>          v(j+1) = VR(:,j) - i*VR(:,j+1).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR.  LDVR >= 1, and if<br>
*>          JOBVR = 'V', LDVR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ILO<br>
*> \verbatim<br>
*>          ILO is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[out] IHI<br>
*> \verbatim<br>
*>          IHI is INTEGER<br>
*>          ILO and IHI are integer values determined when A was<br>
*>          balanced.  The balanced A(i,j) = 0 if I > J and<br>
*>          J = 1,...,ILO-1 or I = IHI+1,...,N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is DOUBLE PRECISION array, dimension (N)<br>
*>          Details of the permutations and scaling factors applied<br>
*>          when balancing A.  If P(j) is the index of the row and column<br>
*>          interchanged with row and column j, and D(j) is the scaling<br>
*>          factor applied to row and column j, then<br>
*>          SCALE(J) = P(J),    for J = 1,...,ILO-1<br>
*>                   = D(J),    for J = ILO,...,IHI<br>
*>                   = P(J)     for J = IHI+1,...,N.<br>
*>          The order in which the interchanges are made is N to IHI+1,<br>
*>          then 1 to ILO-1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ABNRM<br>
*> \verbatim<br>
*>          ABNRM is DOUBLE PRECISION<br>
*>          The one-norm of the balanced matrix (the maximum<br>
*>          of the sum of absolute values of elements of any column).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCONDE<br>
*> \verbatim<br>
*>          RCONDE is DOUBLE PRECISION array, dimension (N)<br>
*>          RCONDE(j) is the reciprocal condition number of the j-th<br>
*>          eigenvalue.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCONDV<br>
*> \verbatim<br>
*>          RCONDV is DOUBLE PRECISION array, dimension (N)<br>
*>          RCONDV(j) is the reciprocal condition number of the j-th<br>
*>          right eigenvector.<br>
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
*>          The dimension of the array WORK.   If SENSE = 'N' or 'E',<br>
*>          LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V',<br>
*>          LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6).<br>
*>          For good performance, LWORK must generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (2*N-2)<br>
*>          If SENSE = 'N' or 'E', not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = i, the QR algorithm failed to compute all the<br>
*>                eigenvalues, and no eigenvectors or condition numbers<br>
*>                have been computed; elements 1:ILO-1 and i+1:N of WR<br>
*>                and WI contain eigenvalues which have converged.<br>
*> \endverbatim<br>
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
*  @precisions fortran d -> s<br>
*<br>
*> \ingroup doubleGEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dgeevx_(CHARACTER BALANC,CHARACTER JOBVL,CHARACTER JOBVR,CHARACTER SENSE,INTEGER N,double[] A,INTEGER LDA,double[] WR,double[] WI,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER ILO,INTEGER IHI,double[] SCALE,DOUBLE ABNRM,double[] RCONDE,double[] RCONDV,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEHD2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgehd2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgehd2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgehd2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, ILO, INFO, LDA, N<br>
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
*> DGEHD2 reduces a real general matrix A to upper Hessenberg form H by<br>
*> an orthogonal similarity transformation:  Q**T * A * Q = H .<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
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
*>          It is assumed that A is already upper triangular in rows<br>
*>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally<br>
*>          set by a previous call to DGEBAL; otherwise they should be<br>
*>          set to 1 and N respectively. See Further Details.<br>
*>          1 <= ILO <= IHI <= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the n by n general matrix to be reduced.<br>
*>          On exit, the upper triangle and the first subdiagonal of A<br>
*>          are overwritten with the upper Hessenberg matrix H, and the<br>
*>          elements below the first subdiagonal, with the array TAU,<br>
*>          represent the orthogonal matrix Q as a product of elementary<br>
*>          reflectors. See Further Details.<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of (ihi-ilo) elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on<br>
*>  exit in A(i+2:ihi,i), and tau in TAU(i).<br>
*><br>
*>  The contents of A are illustrated by the following example, with<br>
*>  n = 7, ilo = 2 and ihi = 6:<br>
*><br>
*>  on entry,                        on exit,<br>
*><br>
*>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )<br>
*>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )<br>
*>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )<br>
*>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )<br>
*>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )<br>
*>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )<br>
*>  (                         a )    (                          a )<br>
*><br>
*>  where a denotes an element of the original matrix A, h denotes a<br>
*>  modified element of the upper Hessenberg matrix H, and vi denotes an<br>
*>  element of the vector defining H(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgehd2_(INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGEHRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEHRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgehrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgehrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgehrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION  A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEHRD reduces a real general matrix A to upper Hessenberg form H by<br>
*> an orthogonal similarity transformation:  Q**T * A * Q = H .<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
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
*>          It is assumed that A is already upper triangular in rows<br>
*>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally<br>
*>          set by a previous call to DGEBAL; otherwise they should be<br>
*>          set to 1 and N respectively. See Further Details.<br>
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the N-by-N general matrix to be reduced.<br>
*>          On exit, the upper triangle and the first subdiagonal of A<br>
*>          are overwritten with the upper Hessenberg matrix H, and the<br>
*>          elements below the first subdiagonal, with the array TAU,<br>
*>          represent the orthogonal matrix Q as a product of elementary<br>
*>          reflectors. See Further Details.<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to<br>
*>          zero.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of the array WORK.  LWORK >= max(1,N).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of (ihi-ilo) elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on<br>
*>  exit in A(i+2:ihi,i), and tau in TAU(i).<br>
*><br>
*>  The contents of A are illustrated by the following example, with<br>
*>  n = 7, ilo = 2 and ihi = 6:<br>
*><br>
*>  on entry,                        on exit,<br>
*><br>
*>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )<br>
*>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )<br>
*>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )<br>
*>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )<br>
*>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )<br>
*>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )<br>
*>  (                         a )    (                          a )<br>
*><br>
*>  where a denotes an element of the original matrix A, h denotes a<br>
*>  modified element of the upper Hessenberg matrix H, and vi denotes an<br>
*>  element of the vector defining H(i).<br>
*><br>
*>  This file is a slight modification of LAPACK-3.0's DGEHRD<br>
*>  subroutine incorporating improvements proposed by Quintana-Orti and<br>
*>  Van de Geijn (2006). (See DLAHR2.)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgehrd_(INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGEJSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEJSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgejsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgejsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgejsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP,<br>
*                          M, N, A, LDA, SVA, U, LDU, V, LDV,<br>
*                          WORK, LWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       IMPLICIT    NONE<br>
*       INTEGER     INFO, LDA, LDU, LDV, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A( LDA, * ), SVA( N ), U( LDU, * ), V( LDV, * ),<br>
*      $            WORK( LWORK )<br>
*       INTEGER     IWORK( * )<br>
*       CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEJSV computes the singular value decomposition (SVD) of a real M-by-N<br>
*> matrix [A], where M >= N. The SVD of [A] is written as<br>
*><br>
*>              [A] = [U] * [SIGMA] * [V]^t,<br>
*><br>
*> where [SIGMA] is an N-by-N (M-by-N) matrix which is zero except for its N<br>
*> diagonal elements, [U] is an M-by-N (or M-by-M) orthonormal matrix, and<br>
*> [V] is an N-by-N orthogonal matrix. The diagonal elements of [SIGMA] are<br>
*> the singular values of [A]. The columns of [U] and [V] are the left and<br>
*> the right singular vectors of [A], respectively. The matrices [U] and [V]<br>
*> are computed and stored in the arrays U and V, respectively. The diagonal<br>
*> of [SIGMA] is computed and stored in the array SVA.<br>
*> DGEJSV can sometimes compute tiny singular values and their singular vectors much<br>
*> more accurately than other SVD routines, see below under Further Details.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBA<br>
*> \verbatim<br>
*>          JOBA is CHARACTER*1<br>
*>        Specifies the level of accuracy:<br>
*>       = 'C': This option works well (high relative accuracy) if A = B * D,<br>
*>             with well-conditioned B and arbitrary diagonal matrix D.<br>
*>             The accuracy cannot be spoiled by COLUMN scaling. The<br>
*>             accuracy of the computed output depends on the condition of<br>
*>             B, and the procedure aims at the best theoretical accuracy.<br>
*>             The relative error max_{i=1:N}|d sigma_i| / sigma_i is<br>
*>             bounded by f(M,N)*epsilon* cond(B), independent of D.<br>
*>             The input matrix is preprocessed with the QRF with column<br>
*>             pivoting. This initial preprocessing and preconditioning by<br>
*>             a rank revealing QR factorization is common for all values of<br>
*>             JOBA. Additional actions are specified as follows:<br>
*>       = 'E': Computation as with 'C' with an additional estimate of the<br>
*>             condition number of B. It provides a realistic error bound.<br>
*>       = 'F': If A = D1 * C * D2 with ill-conditioned diagonal scalings<br>
*>             D1, D2, and well-conditioned matrix C, this option gives<br>
*>             higher accuracy than the 'C' option. If the structure of the<br>
*>             input matrix is not known, and relative accuracy is<br>
*>             desirable, then this option is advisable. The input matrix A<br>
*>             is preprocessed with QR factorization with FULL (row and<br>
*>             column) pivoting.<br>
*>       = 'G'  Computation as with 'F' with an additional estimate of the<br>
*>             condition number of B, where A=D*B. If A has heavily weighted<br>
*>             rows, then using this condition number gives too pessimistic<br>
*>             error bound.<br>
*>       = 'A': Small singular values are the noise and the matrix is treated<br>
*>             as numerically rank defficient. The error in the computed<br>
*>             singular values is bounded by f(m,n)*epsilon*||A||.<br>
*>             The computed SVD A = U * S * V^t restores A up to<br>
*>             f(m,n)*epsilon*||A||.<br>
*>             This gives the procedure the licence to discard (set to zero)<br>
*>             all singular values below N*epsilon*||A||.<br>
*>       = 'R': Similar as in 'A'. Rank revealing property of the initial<br>
*>             QR factorization is used do reveal (using triangular factor)<br>
*>             a gap sigma_{r+1} < epsilon * sigma_r in which case the<br>
*>             numerical RANK is declared to be r. The SVD is computed with<br>
*>             absolute error bounds, but more accurately than with 'A'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBU<br>
*> \verbatim<br>
*>          JOBU is CHARACTER*1<br>
*>        Specifies whether to compute the columns of U:<br>
*>       = 'U': N columns of U are returned in the array U.<br>
*>       = 'F': full set of M left sing. vectors is returned in the array U.<br>
*>       = 'W': U may be used as workspace of length M*N. See the description<br>
*>             of U.<br>
*>       = 'N': U is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV<br>
*> \verbatim<br>
*>          JOBV is CHARACTER*1<br>
*>        Specifies whether to compute the matrix V:<br>
*>       = 'V': N columns of V are returned in the array V; Jacobi rotations<br>
*>             are not explicitly accumulated.<br>
*>       = 'J': N columns of V are returned in the array V, but they are<br>
*>             computed as the product of Jacobi rotations. This option is<br>
*>             allowed only if JOBU .NE. 'N', i.e. in computing the full SVD.<br>
*>       = 'W': V may be used as workspace of length N*N. See the description<br>
*>             of V.<br>
*>       = 'N': V is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBR<br>
*> \verbatim<br>
*>          JOBR is CHARACTER*1<br>
*>        Specifies the RANGE for the singular values. Issues the licence to<br>
*>        set to zero small positive singular values if they are outside<br>
*>        specified range. If A .NE. 0 is scaled so that the largest singular<br>
*>        value of c*A is around DSQRT(BIG), BIG=SLAMCH('O'), then JOBR issues<br>
*>        the licence to kill columns of A whose norm in c*A is less than<br>
*>        DSQRT(SFMIN) (for JOBR.EQ.'R'), or less than SMALL=SFMIN/EPSLN,<br>
*>        where SFMIN=SLAMCH('S'), EPSLN=SLAMCH('E').<br>
*>       = 'N': Do not kill small columns of c*A. This option assumes that<br>
*>             BLAS and QR factorizations and triangular solvers are<br>
*>             implemented to work in that range. If the condition of A<br>
*>             is greater than BIG, use DGESVJ.<br>
*>       = 'R': RESTRICTED range for sigma(c*A) is [DSQRT(SFMIN), DSQRT(BIG)]<br>
*>             (roughly, as described above). This option is recommended.<br>
*>                                            ~~~~~~~~~~~~~~~~~~~~~~~~~~~<br>
*>        For computing the singular values in the FULL range [SFMIN,BIG]<br>
*>        use DGESVJ.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBT<br>
*> \verbatim<br>
*>          JOBT is CHARACTER*1<br>
*>        If the matrix is square then the procedure may determine to use<br>
*>        transposed A if A^t seems to be better with respect to convergence.<br>
*>        If the matrix is not square, JOBT is ignored. This is subject to<br>
*>        changes in the future.<br>
*>        The decision is based on two values of entropy over the adjoint<br>
*>        orbit of A^t * A. See the descriptions of WORK(6) and WORK(7).<br>
*>       = 'T': transpose if entropy test indicates possibly faster<br>
*>        convergence of Jacobi process if A^t is taken as input. If A is<br>
*>        replaced with A^t, then the row pivoting is included automatically.<br>
*>       = 'N': do not speculate.<br>
*>        This option can be used to compute only the singular values, or the<br>
*>        full SVD (U, SIGMA and V). For only one set of singular vectors<br>
*>        (U or V), the caller should provide both U and V, as one of the<br>
*>        matrices is used as workspace if the matrix A is transposed.<br>
*>        The implementer can easily remove this constraint and make the<br>
*>        code more complicated. See the descriptions of U and V.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBP<br>
*> \verbatim<br>
*>          JOBP is CHARACTER*1<br>
*>        Issues the licence to introduce structured perturbations to drown<br>
*>        denormalized numbers. This licence should be active if the<br>
*>        denormals are poorly implemented, causing slow computation,<br>
*>        especially in cases of fast convergence (!). For details see [1,2].<br>
*>        For the sake of simplicity, this perturbations are included only<br>
*>        when the full SVD or only the singular values are requested. The<br>
*>        implementer/user can easily add the perturbation for the cases of<br>
*>        computing one set of singular vectors.<br>
*>       = 'P': introduce perturbation<br>
*>       = 'N': do not perturb<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>         The number of rows of the input matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         The number of columns of the input matrix A. M >= N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SVA<br>
*> \verbatim<br>
*>          SVA is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit,<br>
*>          - For WORK(1)/WORK(2) = ONE: The singular values of A. During the<br>
*>            computation SVA contains Euclidean column norms of the<br>
*>            iterated matrices in the array A.<br>
*>          - For WORK(1) .NE. WORK(2): The singular values of A are<br>
*>            (WORK(1)/WORK(2)) * SVA(1:N). This factored form is used if<br>
*>            sigma_max(A) overflows or if small singular values have been<br>
*>            saved from underflow by scaling the input matrix A.<br>
*>          - If JOBR='R' then some of the singular values may be returned<br>
*>            as exact zeros obtained by "set to zero" because they are<br>
*>            below the numerical rank threshold or are denormalized numbers.<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is DOUBLE PRECISION array, dimension ( LDU, N )<br>
*>          If JOBU = 'U', then U contains on exit the M-by-N matrix of<br>
*>                         the left singular vectors.<br>
*>          If JOBU = 'F', then U contains on exit the M-by-M matrix of<br>
*>                         the left singular vectors, including an ONB<br>
*>                         of the orthogonal complement of the Range(A).<br>
*>          If JOBU = 'W'  .AND. (JOBV.EQ.'V' .AND. JOBT.EQ.'T' .AND. M.EQ.N),<br>
*>                         then U is used as workspace if the procedure<br>
*>                         replaces A with A^t. In that case, [V] is computed<br>
*>                         in U as left singular vectors of A^t and then<br>
*>                         copied back to the V array. This 'W' option is just<br>
*>                         a reminder to the caller that in this case U is<br>
*>                         reserved as workspace of length N*N.<br>
*>          If JOBU = 'N'  U is not referenced, unless JOBT='T'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>          The leading dimension of the array U,  LDU >= 1.<br>
*>          IF  JOBU = 'U' or 'F' or 'W',  then LDU >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension ( LDV, N )<br>
*>          If JOBV = 'V', 'J' then V contains on exit the N-by-N matrix of<br>
*>                         the right singular vectors;<br>
*>          If JOBV = 'W', AND (JOBU.EQ.'U' AND JOBT.EQ.'T' AND M.EQ.N),<br>
*>                         then V is used as workspace if the pprocedure<br>
*>                         replaces A with A^t. In that case, [U] is computed<br>
*>                         in V as right singular vectors of A^t and then<br>
*>                         copied back to the U array. This 'W' option is just<br>
*>                         a reminder to the caller that in this case V is<br>
*>                         reserved as workspace of length N*N.<br>
*>          If JOBV = 'N'  V is not referenced, unless JOBT='T'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V,  LDV >= 1.<br>
*>          If JOBV = 'V' or 'J' or 'W', then LDV >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension at least LWORK.<br>
*>          On exit, if N.GT.0 .AND. M.GT.0 (else not referenced), <br>
*>          WORK(1) = SCALE = WORK(2) / WORK(1) is the scaling factor such<br>
*>                    that SCALE*SVA(1:N) are the computed singular values<br>
*>                    of A. (See the description of SVA().)<br>
*>          WORK(2) = See the description of WORK(1).<br>
*>          WORK(3) = SCONDA is an estimate for the condition number of<br>
*>                    column equilibrated A. (If JOBA .EQ. 'E' or 'G')<br>
*>                    SCONDA is an estimate of DSQRT(||(R^t * R)^(-1)||_1).<br>
*>                    It is computed using DPOCON. It holds<br>
*>                    N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA<br>
*>                    where R is the triangular factor from the QRF of A.<br>
*>                    However, if R is truncated and the numerical rank is<br>
*>                    determined to be strictly smaller than N, SCONDA is<br>
*>                    returned as -1, thus indicating that the smallest<br>
*>                    singular values might be lost.<br>
*><br>
*>          If full SVD is needed, the following two condition numbers are<br>
*>          useful for the analysis of the algorithm. They are provied for<br>
*>          a developer/implementer who is familiar with the details of<br>
*>          the method.<br>
*><br>
*>          WORK(4) = an estimate of the scaled condition number of the<br>
*>                    triangular factor in the first QR factorization.<br>
*>          WORK(5) = an estimate of the scaled condition number of the<br>
*>                    triangular factor in the second QR factorization.<br>
*>          The following two parameters are computed if JOBT .EQ. 'T'.<br>
*>          They are provided for a developer/implementer who is familiar<br>
*>          with the details of the method.<br>
*><br>
*>          WORK(6) = the entropy of A^t*A :: this is the Shannon entropy<br>
*>                    of diag(A^t*A) / Trace(A^t*A) taken as point in the<br>
*>                    probability simplex.<br>
*>          WORK(7) = the entropy of A*A^t.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          Length of WORK to confirm proper allocation of work space.<br>
*>          LWORK depends on the job:<br>
*><br>
*>          If only SIGMA is needed ( JOBU.EQ.'N', JOBV.EQ.'N' ) and<br>
*>            -> .. no scaled condition estimate required (JOBE.EQ.'N'):<br>
*>               LWORK >= max(2*M+N,4*N+1,7). This is the minimal requirement.<br>
*>               ->> For optimal performance (blocked code) the optimal value<br>
*>               is LWORK >= max(2*M+N,3*N+(N+1)*NB,7). Here NB is the optimal<br>
*>               block size for DGEQP3 and DGEQRF.<br>
*>               In general, optimal LWORK is computed as <br>
*>               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DGEQRF), 7).        <br>
*>            -> .. an estimate of the scaled condition number of A is<br>
*>               required (JOBA='E', 'G'). In this case, LWORK is the maximum<br>
*>               of the above and N*N+4*N, i.e. LWORK >= max(2*M+N,N*N+4*N,7).<br>
*>               ->> For optimal performance (blocked code) the optimal value <br>
*>               is LWORK >= max(2*M+N,3*N+(N+1)*NB, N*N+4*N, 7).<br>
*>               In general, the optimal length LWORK is computed as<br>
*>               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DGEQRF), <br>
*>                                                     N+N*N+LWORK(DPOCON),7).<br>
*><br>
*>          If SIGMA and the right singular vectors are needed (JOBV.EQ.'V'),<br>
*>            -> the minimal requirement is LWORK >= max(2*M+N,4*N+1,7).<br>
*>            -> For optimal performance, LWORK >= max(2*M+N,3*N+(N+1)*NB,7),<br>
*>               where NB is the optimal block size for DGEQP3, DGEQRF, DGELQF,<br>
*>               DORMLQ. In general, the optimal length LWORK is computed as<br>
*>               LWORK >= max(2*M+N,N+LWORK(DGEQP3), N+LWORK(DPOCON), <br>
*>                       N+LWORK(DGELQF), 2*N+LWORK(DGEQRF), N+LWORK(DORMLQ)).<br>
*><br>
*>          If SIGMA and the left singular vectors are needed<br>
*>            -> the minimal requirement is LWORK >= max(2*M+N,4*N+1,7).<br>
*>            -> For optimal performance:<br>
*>               if JOBU.EQ.'U' :: LWORK >= max(2*M+N,3*N+(N+1)*NB,7),<br>
*>               if JOBU.EQ.'F' :: LWORK >= max(2*M+N,3*N+(N+1)*NB,N+M*NB,7),<br>
*>               where NB is the optimal block size for DGEQP3, DGEQRF, DORMQR.<br>
*>               In general, the optimal length LWORK is computed as<br>
*>               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DPOCON),<br>
*>                        2*N+LWORK(DGEQRF), N+LWORK(DORMQR)). <br>
*>               Here LWORK(DORMQR) equals N*NB (for JOBU.EQ.'U') or <br>
*>               M*NB (for JOBU.EQ.'F').<br>
*>               <br>
*>          If the full SVD is needed: (JOBU.EQ.'U' or JOBU.EQ.'F') and <br>
*>            -> if JOBV.EQ.'V'  <br>
*>               the minimal requirement is LWORK >= max(2*M+N,6*N+2*N*N). <br>
*>            -> if JOBV.EQ.'J' the minimal requirement is <br>
*>               LWORK >= max(2*M+N, 4*N+N*N,2*N+N*N+6).<br>
*>            -> For optimal performance, LWORK should be additionally<br>
*>               larger than N+M*NB, where NB is the optimal block size<br>
*>               for DORMQR.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension M+3*N.<br>
*>          On exit,<br>
*>          IWORK(1) = the numerical rank determined after the initial<br>
*>                     QR factorization with pivoting. See the descriptions<br>
*>                     of JOBA and JOBR.<br>
*>          IWORK(2) = the number of the computed nonzero singular values<br>
*>          IWORK(3) = if nonzero, a warning message:<br>
*>                     If IWORK(3).EQ.1 then some of the column norms of A<br>
*>                     were denormalized floats. The requested high accuracy<br>
*>                     is not warranted by the data.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           < 0  : if INFO = -i, then the i-th argument had an illegal value.<br>
*>           = 0 :  successfull exit;<br>
*>           > 0 :  DGEJSV  did not converge in the maximal allowed number<br>
*>                  of sweeps. The computed values may be inaccurate.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsing<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  DGEJSV implements a preconditioned Jacobi SVD algorithm. It uses DGEQP3,<br>
*>  DGEQRF, and DGELQF as preprocessors and preconditioners. Optionally, an<br>
*>  additional row pivoting can be used as a preprocessor, which in some<br>
*>  cases results in much higher accuracy. An example is matrix A with the<br>
*>  structure A = D1 * C * D2, where D1, D2 are arbitrarily ill-conditioned<br>
*>  diagonal matrices and C is well-conditioned matrix. In that case, complete<br>
*>  pivoting in the first QR factorizations provides accuracy dependent on the<br>
*>  condition number of C, and independent of D1, D2. Such higher accuracy is<br>
*>  not completely understood theoretically, but it works well in practice.<br>
*>  Further, if A can be written as A = B*D, with well-conditioned B and some<br>
*>  diagonal D, then the high accuracy is guaranteed, both theoretically and<br>
*>  in software, independent of D. For more details see [1], [2].<br>
*>     The computational range for the singular values can be the full range<br>
*>  ( UNDERFLOW,OVERFLOW ), provided that the machine arithmetic and the BLAS<br>
*>  & LAPACK routines called by DGEJSV are implemented to work in that range.<br>
*>  If that is not the case, then the restriction for safe computation with<br>
*>  the singular values in the range of normalized IEEE numbers is that the<br>
*>  spectral condition number kappa(A)=sigma_max(A)/sigma_min(A) does not<br>
*>  overflow. This code (DGEJSV) is best used in this restricted range,<br>
*>  meaning that singular values of magnitude below ||A||_2 / DLAMCH('O') are<br>
*>  returned as zeros. See JOBR for details on this.<br>
*>     Further, this implementation is somewhat slower than the one described<br>
*>  in [1,2] due to replacement of some non-LAPACK components, and because<br>
*>  the choice of some tuning parameters in the iterative part (DGESVJ) is<br>
*>  left to the implementer on a particular machine.<br>
*>     The rank revealing QR factorization (in this code: DGEQP3) should be<br>
*>  implemented as in [3]. We have a new version of DGEQP3 under development<br>
*>  that is more robust than the current one in LAPACK, with a cleaner cut in<br>
*>  rank defficient cases. It will be available in the SIGMA library [4].<br>
*>  If M is much larger than N, it is obvious that the inital QRF with<br>
*>  column pivoting can be preprocessed by the QRF without pivoting. That<br>
*>  well known trick is not used in DGEJSV because in some cases heavy row<br>
*>  weighting can be treated with complete pivoting. The overhead in cases<br>
*>  M much larger than N is then only due to pivoting, but the benefits in<br>
*>  terms of accuracy have prevailed. The implementer/user can incorporate<br>
*>  this extra QRF step easily. The implementer can also improve data movement<br>
*>  (matrix transpose, matrix copy, matrix transposed copy) - this<br>
*>  implementation of DGEJSV uses only the simplest, naive data movement.<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*> \verbatim<br>
*><br>
*> [1] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I.<br>
*>     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342.<br>
*>     LAPACK Working note 169.<br>
*> [2] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II.<br>
*>     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362.<br>
*>     LAPACK Working note 170.<br>
*> [3] Z. Drmac and Z. Bujanovic: On the failure of rank-revealing QR<br>
*>     factorization software - a case study.<br>
*>     ACM Trans. Math. Softw. Vol. 35, No 2 (2008), pp. 1-28.<br>
*>     LAPACK Working note 176.<br>
*> [4] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV,<br>
*>     QSVD, (H,K)-SVD computations.<br>
*>     Department of Mathematics, University of Zagreb, 2008.<br>
*> \endverbatim<br>
*<br>
*>  \par Bugs, examples and comments:<br>
*   =================================<br>
*><br>
*>  Please report all bugs and send interesting examples and/or comments to<br>
*>  drmac@math.hr. Thank you.<br>
*><br>
*  =====================================================================<br>
*/
	public void dgejsv_(char[] JOBA,char[] JOBU,char[] JOBV,char[] JOBR,char[] JOBT,char[] JOBP,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] SVA,double[] U,INTEGER LDU,double[] V,INTEGER LDV,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGELQ2 computes the LQ factorization of a general rectangular matrix using an unblocked algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGELQ2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelq2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelq2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelq2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGELQ2( M, N, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
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
*> DGELQ2 computes an LQ factorization of a real m by n matrix A:<br>
*> A = L * Q.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the m by n matrix A.<br>
*>          On exit, the elements on and below the diagonal of the array<br>
*>          contain the m by min(m,n) lower trapezoidal matrix L (L is<br>
*>          lower triangular if m <= n); the elements above the diagonal,<br>
*>          with the array TAU, represent the orthogonal matrix Q as a<br>
*>          product of elementary reflectors (see Further Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(k) . . . H(2) H(1), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),<br>
*>  and tau in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgelq2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGELQF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGELQF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelqf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelqf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelqf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
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
*> DGELQF computes an LQ factorization of a real M-by-N matrix A:<br>
*> A = L * Q.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the elements on and below the diagonal of the array<br>
*>          contain the m-by-min(m,n) lower trapezoidal matrix L (L is<br>
*>          lower triangular if m <= n); the elements above the diagonal,<br>
*>          with the array TAU, represent the orthogonal matrix Q as a<br>
*>          product of elementary reflectors (see Further Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          The dimension of the array WORK.  LWORK >= max(1,M).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(k) . . . H(2) H(1), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),<br>
*>  and tau in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgelqf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGELS solves overdetermined or underdetermined systems for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGELS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgels.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgels.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgels.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,<br>
*                         INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGELS solves overdetermined or underdetermined real linear systems<br>
*> involving an M-by-N matrix A, or its transpose, using a QR or LQ<br>
*> factorization of A.  It is assumed that A has full rank.<br>
*><br>
*> The following options are provided:<br>
*><br>
*> 1. If TRANS = 'N' and m >= n:  find the least squares solution of<br>
*>    an overdetermined system, i.e., solve the least squares problem<br>
*>                 minimize || B - A*X ||.<br>
*><br>
*> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of<br>
*>    an underdetermined system A * X = B.<br>
*><br>
*> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of<br>
*>    an undetermined system A**T * X = B.<br>
*><br>
*> 4. If TRANS = 'T' and m < n:  find the least squares solution of<br>
*>    an overdetermined system, i.e., solve the least squares problem<br>
*>                 minimize || B - A**T * X ||.<br>
*><br>
*> Several right hand side vectors b and solution vectors x can be<br>
*> handled in a single call; they are stored as the columns of the<br>
*> M-by-NRHS right hand side matrix B and the N-by-NRHS solution<br>
*> matrix X.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': the linear system involves A;<br>
*>          = 'T': the linear system involves A**T.<br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of<br>
*>          columns of the matrices B and X. NRHS >=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit,<br>
*>            if M >= N, A is overwritten by details of its QR<br>
*>                       factorization as returned by DGEQRF;<br>
*>            if M <  N, A is overwritten by details of its LQ<br>
*>                       factorization as returned by DGELQF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the matrix B of right hand side vectors, stored<br>
*>          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS<br>
*>          if TRANS = 'T'.<br>
*>          On exit, if INFO = 0, B is overwritten by the solution<br>
*>          vectors, stored columnwise:<br>
*>          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least<br>
*>          squares solution vectors; the residual sum of squares for the<br>
*>          solution in each column is given by the sum of squares of<br>
*>          elements N+1 to M in that column;<br>
*>          if TRANS = 'N' and m < n, rows 1 to N of B contain the<br>
*>          minimum norm solution vectors;<br>
*>          if TRANS = 'T' and m >= n, rows 1 to M of B contain the<br>
*>          minimum norm solution vectors;<br>
*>          if TRANS = 'T' and m < n, rows 1 to M of B contain the<br>
*>          least squares solution vectors; the residual sum of squares<br>
*>          for the solution in each column is given by the sum of<br>
*>          squares of elements M+1 to N in that column.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= MAX(1,M,N).<br>
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
*>          LWORK >= max( 1, MN + max( MN, NRHS ) ).<br>
*>          For optimal performance,<br>
*>          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).<br>
*>          where MN = min(M,N) and NB is the optimum block size.<br>
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
*>          > 0:  if INFO =  i, the i-th diagonal element of the<br>
*>                triangular factor of A is zero, so that A does not have<br>
*>                full rank; the least squares solution could not be<br>
*>                computed.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgels_(CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGELSD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelsd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelsd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelsd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,<br>
*                          WORK, LWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGELSD computes the minimum-norm solution to a real linear least<br>
*> squares problem:<br>
*>     minimize 2-norm(| b - A*x |)<br>
*> using the singular value decomposition (SVD) of A. A is an M-by-N<br>
*> matrix which may be rank-deficient.<br>
*><br>
*> Several right hand side vectors b and solution vectors x can be<br>
*> handled in a single call; they are stored as the columns of the<br>
*> M-by-NRHS right hand side matrix B and the N-by-NRHS solution<br>
*> matrix X.<br>
*><br>
*> The problem is solved in three steps:<br>
*> (1) Reduce the coefficient matrix A to bidiagonal form with<br>
*>     Householder transformations, reducing the original problem<br>
*>     into a "bidiagonal least squares problem" (BLS)<br>
*> (2) Solve the BLS using a divide and conquer approach.<br>
*> (3) Apply back all the Householder tranformations to solve<br>
*>     the original least squares problem.<br>
*><br>
*> The effective rank of A is determined by treating as zero those<br>
*> singular values which are less than RCOND times the largest singular<br>
*> value.<br>
*><br>
*> The divide and conquer algorithm makes very mild assumptions about<br>
*> floating point arithmetic. It will work on machines with a guard<br>
*> digit in add/subtract, or on those binary machines without guard<br>
*> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or<br>
*> Cray-2. It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of A. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrices B and X. NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, A has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the M-by-NRHS right hand side matrix B.<br>
*>          On exit, B is overwritten by the N-by-NRHS solution<br>
*>          matrix X.  If m >= n and RANK = n, the residual<br>
*>          sum-of-squares for the solution in the i-th column is given<br>
*>          by the sum of squares of elements n+1:m in that column.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,max(M,N)).<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The singular values of A in decreasing order.<br>
*>          The condition number of A in the 2-norm = S(1)/S(min(m,n)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          RCOND is used to determine the effective rank of A.<br>
*>          Singular values S(i) <= RCOND*S(1) are treated as zero.<br>
*>          If RCOND < 0, machine precision is used instead.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RANK<br>
*> \verbatim<br>
*>          RANK is INTEGER<br>
*>          The effective rank of A, i.e., the number of singular values<br>
*>          which are greater than RCOND*S(1).<br>
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
*>          The dimension of the array WORK. LWORK must be at least 1.<br>
*>          The exact minimum amount of workspace needed depends on M,<br>
*>          N and NRHS. As long as LWORK is at least<br>
*>              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,<br>
*>          if M is greater than or equal to N or<br>
*>              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,<br>
*>          if M is less than N, the code will execute correctly.<br>
*>          SMLSIZ is returned by ILAENV and is equal to the maximum<br>
*>          size of the subproblems at the bottom of the computation<br>
*>          tree (usually about 25), and<br>
*>             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )<br>
*>          For good performance, LWORK should generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))<br>
*>          LIWORK >= max(1, 3 * MINMN * NLVL + 11 * MINMN),<br>
*>          where MINMN = MIN( M,N ).<br>
*>          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  the algorithm for computing the SVD failed to converge;<br>
*>                if INFO = i, i off-diagonal elements of an intermediate<br>
*>                bidiagonal form did not converge to zero.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsolve<br>
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
	public void dgelsd_(INTEGER M,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] S,DOUBLE RCOND,INTEGER RANK,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> DGELSS solves overdetermined or underdetermined systems for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGELSS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelss.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelss.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelss.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGELSS computes the minimum norm solution to a real linear least<br>
*> squares problem:<br>
*><br>
*> Minimize 2-norm(| b - A*x |).<br>
*><br>
*> using the singular value decomposition (SVD) of A. A is an M-by-N<br>
*> matrix which may be rank-deficient.<br>
*><br>
*> Several right hand side vectors b and solution vectors x can be<br>
*> handled in a single call; they are stored as the columns of the<br>
*> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix<br>
*> X.<br>
*><br>
*> The effective rank of A is determined by treating as zero those<br>
*> singular values which are less than RCOND times the largest singular<br>
*> value.<br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrices B and X. NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the first min(m,n) rows of A are overwritten with<br>
*>          its right singular vectors, stored rowwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the M-by-NRHS right hand side matrix B.<br>
*>          On exit, B is overwritten by the N-by-NRHS solution<br>
*>          matrix X.  If m >= n and RANK = n, the residual<br>
*>          sum-of-squares for the solution in the i-th column is given<br>
*>          by the sum of squares of elements n+1:m in that column.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,max(M,N)).<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The singular values of A in decreasing order.<br>
*>          The condition number of A in the 2-norm = S(1)/S(min(m,n)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          RCOND is used to determine the effective rank of A.<br>
*>          Singular values S(i) <= RCOND*S(1) are treated as zero.<br>
*>          If RCOND < 0, machine precision is used instead.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RANK<br>
*> \verbatim<br>
*>          RANK is INTEGER<br>
*>          The effective rank of A, i.e., the number of singular values<br>
*>          which are greater than RCOND*S(1).<br>
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
*>          The dimension of the array WORK. LWORK >= 1, and also:<br>
*>          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )<br>
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
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  the algorithm for computing the SVD failed to converge;<br>
*>                if INFO = i, i off-diagonal elements of an intermediate<br>
*>                bidiagonal form did not converge to zero.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgelss_(INTEGER M,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] S,DOUBLE RCOND,INTEGER RANK,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGELSY solves overdetermined or underdetermined systems for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGELSY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelsy.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelsy.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelsy.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            JPVT( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGELSY computes the minimum-norm solution to a real linear least<br>
*> squares problem:<br>
*>     minimize || A * X - B ||<br>
*> using a complete orthogonal factorization of A.  A is an M-by-N<br>
*> matrix which may be rank-deficient.<br>
*><br>
*> Several right hand side vectors b and solution vectors x can be<br>
*> handled in a single call; they are stored as the columns of the<br>
*> M-by-NRHS right hand side matrix B and the N-by-NRHS solution<br>
*> matrix X.<br>
*><br>
*> The routine first computes a QR factorization with column pivoting:<br>
*>     A * P = Q * [ R11 R12 ]<br>
*>                 [  0  R22 ]<br>
*> with R11 defined as the largest leading submatrix whose estimated<br>
*> condition number is less than 1/RCOND.  The order of R11, RANK,<br>
*> is the effective rank of A.<br>
*><br>
*> Then, R22 is considered to be negligible, and R12 is annihilated<br>
*> by orthogonal transformations from the right, arriving at the<br>
*> complete orthogonal factorization:<br>
*>    A * P = Q * [ T11 0 ] * Z<br>
*>                [  0  0 ]<br>
*> The minimum-norm solution is then<br>
*>    X = P * Z**T [ inv(T11)*Q1**T*B ]<br>
*>                 [        0         ]<br>
*> where Q1 consists of the first RANK columns of Q.<br>
*><br>
*> This routine is basically identical to the original xGELSX except<br>
*> three differences:<br>
*>   o The call to the subroutine xGEQPF has been substituted by the<br>
*>     the call to the subroutine xGEQP3. This subroutine is a Blas-3<br>
*>     version of the QR factorization with column pivoting.<br>
*>   o Matrix B (the right hand side) is updated with Blas-3.<br>
*>   o The permutation of matrix B (the right hand side) is faster and<br>
*>     more simple.<br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of<br>
*>          columns of matrices B and X. NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, A has been overwritten by details of its<br>
*>          complete orthogonal factorization.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the M-by-NRHS right hand side matrix B.<br>
*>          On exit, the N-by-NRHS solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,M,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] JPVT<br>
*> \verbatim<br>
*>          JPVT is INTEGER array, dimension (N)<br>
*>          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted<br>
*>          to the front of AP, otherwise column i is a free column.<br>
*>          On exit, if JPVT(i) = k, then the i-th column of AP<br>
*>          was the k-th column of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          RCOND is used to determine the effective rank of A, which<br>
*>          is defined as the order of the largest leading triangular<br>
*>          submatrix R11 in the QR factorization with pivoting of A,<br>
*>          whose estimated condition number < 1/RCOND.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RANK<br>
*> \verbatim<br>
*>          RANK is INTEGER<br>
*>          The effective rank of A, i.e., the order of the submatrix<br>
*>          R11.  This is the same as the order of the submatrix T11<br>
*>          in the complete orthogonal factorization of A.<br>
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
*>          The unblocked strategy requires that:<br>
*>             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),<br>
*>          where MN = min( M, N ).<br>
*>          The block algorithm requires that:<br>
*>             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),<br>
*>          where NB is an upper bound on the blocksize returned<br>
*>          by ILAENV for the routines DGEQP3, DTZRZF, STZRQF, DORMQR,<br>
*>          and DORMRZ.<br>
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
*>          = 0: successful exit<br>
*>          < 0: If INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsolve<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA \n <br>
*>    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n<br>
*>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n<br>
*><br>
*  =====================================================================<br>
*/
	public void dgelsy_(INTEGER M,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,int[] JPVT,DOUBLE RCOND,INTEGER RANK,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGEMQRT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEMQRT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgemqrt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgemqrt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgemqrt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, <br>
*                          C, LDC, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER SIDE, TRANS<br>
*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEMQRT overwrites the general real M-by-N matrix C with<br>
*><br>
*>                 SIDE = 'L'     SIDE = 'R'<br>
*> TRANS = 'N':      Q C            C Q<br>
*> TRANS = 'T':   Q**T C            C Q**T<br>
*><br>
*> where Q is a real orthogonal matrix defined as the product of K<br>
*> elementary reflectors:<br>
*><br>
*>       Q = H(1) H(2) . . . H(K) = I - V T V**T<br>
*><br>
*> generated using the compact WY representation as returned by DGEQRT. <br>
*><br>
*> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.<br>
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
*>          = 'C':  Transpose, apply Q**T.<br>
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
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The block size used for the storage of T.  K >= NB >= 1.<br>
*>          This must be the same value of NB used to generate T<br>
*>          in CGEQRT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension (LDV,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          CGEQRT in the first K columns of its array argument A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V.<br>
*>          If SIDE = 'L', LDA >= max(1,M);<br>
*>          if SIDE = 'R', LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] T<br>
*> \verbatim<br>
*>          T is DOUBLE PRECISION array, dimension (LDT,K)<br>
*>          The upper triangular factors of the block reflectors<br>
*>          as returned by CGEQRT, stored as a NB-by-N matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (LDC,N)<br>
*>          On entry, the M-by-N matrix C.<br>
*>          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q.<br>
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
*>          WORK is DOUBLE PRECISION array. The dimension of<br>
*>          WORK is N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.<br>
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
*> \date November 2013<br>
*<br>
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgemqrt_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER NB,double[] V,INTEGER LDV,double[] T,INTEGER LDT,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGEQL2 computes the QL factorization of a general rectangular matrix using an unblocked algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQL2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeql2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeql2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeql2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQL2( M, N, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
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
*> DGEQL2 computes a QL factorization of a real m by n matrix A:<br>
*> A = Q * L.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the m by n matrix A.<br>
*>          On exit, if m >= n, the lower triangle of the subarray<br>
*>          A(m-n+1:m,1:n) contains the n by n lower triangular matrix L;<br>
*>          if m <= n, the elements on and below the (n-m)-th<br>
*>          superdiagonal contain the m by n lower trapezoidal matrix L;<br>
*>          the remaining elements, with the array TAU, represent the<br>
*>          orthogonal matrix Q as a product of elementary reflectors<br>
*>          (see Further Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(k) . . . H(2) H(1), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in<br>
*>  A(1:m-k+i-1,n-k+i), and tau in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeql2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGEQLF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQLF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqlf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqlf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqlf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
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
*> DGEQLF computes a QL factorization of a real M-by-N matrix A:<br>
*> A = Q * L.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit,<br>
*>          if m >= n, the lower triangle of the subarray<br>
*>          A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L;<br>
*>          if m <= n, the elements on and below the (n-m)-th<br>
*>          superdiagonal contain the M-by-N lower trapezoidal matrix L;<br>
*>          the remaining elements, with the array TAU, represent the<br>
*>          orthogonal matrix Q as a product of elementary reflectors<br>
*>          (see Further Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          The dimension of the array WORK.  LWORK >= max(1,N).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(k) . . . H(2) H(1), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in<br>
*>  A(1:m-k+i-1,n-k+i), and tau in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqlf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGEQP3<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQP3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqp3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqp3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqp3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            JPVT( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEQP3 computes a QR factorization with column pivoting of a<br>
*> matrix A:  A*P = Q*R  using Level 3 BLAS.<br>
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
*>          The number of columns of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the upper triangle of the array contains the<br>
*>          min(M,N)-by-N upper trapezoidal matrix R; the elements below<br>
*>          the diagonal, together with the array TAU, represent the<br>
*>          orthogonal matrix Q as a product of min(M,N) elementary<br>
*>          reflectors.<br>
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
*>          On entry, if JPVT(J).ne.0, the J-th column of A is permuted<br>
*>          to the front of A*P (a leading column); if JPVT(J)=0,<br>
*>          the J-th column of A is a free column.<br>
*>          On exit, if JPVT(J)=K, then the J-th column of A*P was the<br>
*>          the K-th column of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO=0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= 3*N+1.<br>
*>          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB<br>
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
*>          = 0: successful exit.<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real/complex vector<br>
*>  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in<br>
*>  A(i+1:m,i), and tau in TAU(i).<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain<br>
*>    X. Sun, Computer Science Dept., Duke University, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqp3_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,int[] JPVT,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGEQR2 computes the QR factorization of a general rectangular matrix using an unblocked algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQR2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqr2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqr2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqr2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
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
*> DGEQR2 computes a QR factorization of a real m by n matrix A:<br>
*> A = Q * R.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the m by n matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the min(m,n) by n upper trapezoidal matrix R (R is<br>
*>          upper triangular if m >= n); the elements below the diagonal,<br>
*>          with the array TAU, represent the orthogonal matrix Q as a<br>
*>          product of elementary reflectors (see Further Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),<br>
*>  and tau in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqr2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGEQR2P computes the QR factorization of a general rectangular matrix with non-negative diagonal elements using an unblocked algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQR2P + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqr2p.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqr2p.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqr2p.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQR2P( M, N, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
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
*> DGEQR2 computes a QR factorization of a real m by n matrix A:<br>
*> A = Q * R. The diagonal entries of R are nonnegative.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the m by n matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the min(m,n) by n upper trapezoidal matrix R (R is<br>
*>          upper triangular if m >= n). The diagonal entries of R are<br>
*>          nonnegative; the elements below the diagonal,<br>
*>          with the array TAU, represent the orthogonal matrix Q as a<br>
*>          product of elementary reflectors (see Further Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),<br>
*>  and tau in TAU(i).<br>
*><br>
*> See Lapack Working Note 203 for details<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqr2p_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGEQRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
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
*> DGEQRF computes a QR factorization of a real M-by-N matrix A:<br>
*> A = Q * R.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is<br>
*>          upper triangular if m >= n); the elements below the diagonal,<br>
*>          with the array TAU, represent the orthogonal matrix Q as a<br>
*>          product of min(m,n) elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          The dimension of the array WORK.  LWORK >= max(1,N).<br>
*>          For optimum performance LWORK >= N*NB, where NB is<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),<br>
*>  and tau in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqrf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGEQRFP<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQRFP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrfp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrfp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrfp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
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
*> DGEQRFP computes a QR factorization of a real M-by-N matrix A:<br>
*> A = Q * R. The diagonal entries of R are nonnegative.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is<br>
*>          upper triangular if m >= n). The diagonal entries of R<br>
*>          are nonnegative; the elements below the diagonal,<br>
*>          with the array TAU, represent the orthogonal matrix Q as a<br>
*>          product of min(m,n) elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          The dimension of the array WORK.  LWORK >= max(1,N).<br>
*>          For optimum performance LWORK >= N*NB, where NB is<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),<br>
*>  and tau in TAU(i).<br>
*><br>
*> See Lapack Working Note 203 for details<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqrfp_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGEQRT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQRT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INFO, LDA, LDT, M, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A( LDA, * ), T( LDT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEQRT computes a blocked QR factorization of a real M-by-N matrix A<br>
*> using the compact WY representation of Q.  <br>
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
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The block size to be used in the blocked QR.  MIN(M,N) >= NB >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is<br>
*>          upper triangular if M >= N); the elements below the diagonal<br>
*>          are the columns of V.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is DOUBLE PRECISION array, dimension (LDT,MIN(M,N))<br>
*>          The upper triangular block reflectors stored in compact form<br>
*>          as a sequence of upper triangular blocks.  See below<br>
*>          for further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (NB*N)<br>
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
*> \date November 2013<br>
*<br>
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix V stores the elementary reflectors H(i) in the i-th column<br>
*>  below the diagonal. For example, if M=5 and N=3, the matrix V is<br>
*><br>
*>               V = (  1       )<br>
*>                   ( v1  1    )<br>
*>                   ( v1 v2  1 )<br>
*>                   ( v1 v2 v3 )<br>
*>                   ( v1 v2 v3 )<br>
*><br>
*>  where the vi's represent the vectors which define H(i), which are returned<br>
*>  in the matrix A.  The 1's along the diagonal of V are not stored in A.<br>
*><br>
*>  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where each<br>
*>  block is of order NB except for the last block, which is of order <br>
*>  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular block<br>
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB <br>
*>  for the last block) T's are stored in the NB-by-N matrix T as<br>
*><br>
*>               T = (T1 T2 ... TB).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqrt_(INTEGER M,INTEGER N,INTEGER NB,double[] A,INTEGER LDA,double[] T,INTEGER LDT,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGEQRT2 computes a QR factorization of a general real or complex matrix using the compact WY representation of Q.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQRT2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrt2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrt2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrt2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEQRT2( M, N, A, LDA, T, LDT, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER   INFO, LDA, LDT, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), T( LDT, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEQRT2 computes a QR factorization of a real M-by-N matrix A, <br>
*> using the compact WY representation of Q. <br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= N.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the real M-by-N matrix A.  On exit, the elements on and<br>
*>          above the diagonal contain the N-by-N upper triangular matrix R; the<br>
*>          elements below the diagonal are the columns of V.  See below for<br>
*>          further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is DOUBLE PRECISION array, dimension (LDT,N)<br>
*>          The N-by-N upper triangular factor of the block reflector.<br>
*>          The elements on and above the diagonal contain the block<br>
*>          reflector T; the elements below the diagonal are not used.<br>
*>          See below for further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= max(1,N).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix V stores the elementary reflectors H(i) in the i-th column<br>
*>  below the diagonal. For example, if M=5 and N=3, the matrix V is<br>
*><br>
*>               V = (  1       )<br>
*>                   ( v1  1    )<br>
*>                   ( v1 v2  1 )<br>
*>                   ( v1 v2 v3 )<br>
*>                   ( v1 v2 v3 )<br>
*><br>
*>  where the vi's represent the vectors which define H(i), which are returned<br>
*>  in the matrix A.  The 1's along the diagonal of V are not stored in A.  The<br>
*>  block reflector H is then given by<br>
*><br>
*>               H = I - V * T * V**T<br>
*><br>
*>  where V**T is the transpose of V.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqrt2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] T,INTEGER LDT,INTEGER INFO);
/**
*> \brief \b DGEQRT3 recursively computes a QR factorization of a general real or complex matrix using the compact WY representation of Q.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGEQRT3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrt3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrt3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrt3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       RECURSIVE SUBROUTINE DGEQRT3( M, N, A, LDA, T, LDT, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER   INFO, LDA, M, N, LDT<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), T( LDT, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEQRT3 recursively computes a QR factorization of a real M-by-N <br>
*> matrix A, using the compact WY representation of Q. <br>
*><br>
*> Based on the algorithm of Elmroth and Gustavson, <br>
*> IBM J. Res. Develop. Vol 44 No. 4 July 2000.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= N.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the real M-by-N matrix A.  On exit, the elements on and<br>
*>          above the diagonal contain the N-by-N upper triangular matrix R; the<br>
*>          elements below the diagonal are the columns of V.  See below for<br>
*>          further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is DOUBLE PRECISION array, dimension (LDT,N)<br>
*>          The N-by-N upper triangular factor of the block reflector.<br>
*>          The elements on and above the diagonal contain the block<br>
*>          reflector T; the elements below the diagonal are not used.<br>
*>          See below for further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= max(1,N).<br>
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
*> \date June 2016<br>
*<br>
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix V stores the elementary reflectors H(i) in the i-th column<br>
*>  below the diagonal. For example, if M=5 and N=3, the matrix V is<br>
*><br>
*>               V = (  1       )<br>
*>                   ( v1  1    )<br>
*>                   ( v1 v2  1 )<br>
*>                   ( v1 v2 v3 )<br>
*>                   ( v1 v2 v3 )<br>
*><br>
*>  where the vi's represent the vectors which define H(i), which are returned<br>
*>  in the matrix A.  The 1's along the diagonal of V are not stored in A.  The<br>
*>  block reflector H is then given by<br>
*><br>
*>               H = I - V * T * V**T<br>
*><br>
*>  where V**T is the transpose of V.<br>
*><br>
*>  For details of the algorithm, see Elmroth and Gustavson (cited above).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgeqrt3_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] T,INTEGER LDT,INTEGER INFO);
/**
*> \brief \b DGERFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGERFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgerfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgerfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgerfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,<br>
*                          X, LDX, FERR, BERR, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGERFS improves the computed solution to a system of linear<br>
*> equations and provides error bounds and backward error estimates for<br>
*> the solution.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          The original N-by-N matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] AF<br>
*> \verbatim<br>
*>          AF is DOUBLE PRECISION array, dimension (LDAF,N)<br>
*>          The factors L and U from the factorization A = P*L*U<br>
*>          as computed by DGETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>          The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices from DGETRF; for 1<=i<=N, row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          The right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by DGETRS.<br>
*>          On exit, the improved solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The estimated forward error bound for each solution vector<br>
*>          X(j) (the j-th column of the solution matrix X).<br>
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)<br>
*>          is an estimated upper bound for the magnitude of the largest<br>
*>          element in (X(j) - XTRUE) divided by the magnitude of the<br>
*>          largest element in X(j).  The estimate is as reliable as<br>
*>          the estimate for RCOND, and is almost always a slight<br>
*>          overestimate of the true error.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  ITMAX is the maximum number of steps of iterative refinement.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgerfs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGERFSX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGERFSX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgerfsx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgerfsx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgerfsx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGERFSX( TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
*                           R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS,<br>
*                           ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS,<br>
*                           WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS, EQUED<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   X( LDX , * ), WORK( * )<br>
*       DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ),<br>
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
*>    DGERFSX improves the computed solution to a system of linear<br>
*>    equations and provides error bounds and backward error estimates<br>
*>    for the solution.  In addition to normwise error bound, the code<br>
*>    provides maximum componentwise error bound if possible.  See<br>
*>    comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the<br>
*>    error bounds.<br>
*><br>
*>    The original system of linear equations may have been equilibrated<br>
*>    before calling this routine, as described by arguments EQUED, R<br>
*>    and C below. In this case, the solution and error bounds returned<br>
*>    are for the original unequilibrated system.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \verbatim<br>
*>     Some optional parameters are bundled in the PARAMS array.  These<br>
*>     settings determine how refinement is performed, but often the<br>
*>     defaults are acceptable.  If the defaults are acceptable, users<br>
*>     can pass NPARAMS = 0 which prevents the source code from accessing<br>
*>     the PARAMS argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>     Specifies the form of the system of equations:<br>
*>       = 'N':  A * X = B     (No transpose)<br>
*>       = 'T':  A**T * X = B  (Transpose)<br>
*>       = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>     Specifies the form of equilibration that was done to A<br>
*>     before calling this routine. This is needed to compute<br>
*>     the solution and error bounds correctly.<br>
*>       = 'N':  No equilibration<br>
*>       = 'R':  Row equilibration, i.e., A has been premultiplied by<br>
*>               diag(R).<br>
*>       = 'C':  Column equilibration, i.e., A has been postmultiplied<br>
*>               by diag(C).<br>
*>       = 'B':  Both row and column equilibration, i.e., A has been<br>
*>               replaced by diag(R) * A * diag(C).<br>
*>               The right hand side B has been changed accordingly.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>     The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right hand sides, i.e., the number of columns<br>
*>     of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>     The original N-by-N matrix A.<br>
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
*>          AF is DOUBLE PRECISION array, dimension (LDAF,N)<br>
*>     The factors L and U from the factorization A = P*L*U<br>
*>     as computed by DGETRF.<br>
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
*>     The pivot indices from DGETRF; for 1<=i<=N, row i of the<br>
*>     matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (N)<br>
*>     The row scale factors for A.  If EQUED = 'R' or 'B', A is<br>
*>     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R<br>
*>     is not accessed.  <br>
*>     If R is accessed, each element of R should be a power of the radix<br>
*>     to ensure a reliable solution and error estimates. Scaling by<br>
*>     powers of the radix does not cause rounding errors unless the<br>
*>     result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>     The column scale factors for A.  If EQUED = 'C' or 'B', A is<br>
*>     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C<br>
*>     is not accessed. <br>
*>     If C is accessed, each element of C should be a power of the radix<br>
*>     to ensure a reliable solution and error estimates. Scaling by<br>
*>     powers of the radix does not cause rounding errors unless the<br>
*>     result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>     The right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>     On entry, the solution matrix X, as computed by DGETRS.<br>
*>     On exit, the improved solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>     The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>     Componentwise relative backward error.  This is the<br>
*>     componentwise relative backward error of each solution vector X(j)<br>
*>     (i.e., the smallest relative change in any element of A or B that<br>
*>     makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_ERR_BNDS<br>
*> \verbatim<br>
*>          N_ERR_BNDS is INTEGER<br>
*>     Number of error bounds to return for each right hand side<br>
*>     and each type (normwise or componentwise).  See ERR_BNDS_NORM and<br>
*>     ERR_BNDS_COMP below.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_NORM<br>
*> \verbatim<br>
*>          ERR_BNDS_NORM is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>              sqrt(n) * dlamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * dlamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * dlamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_COMP<br>
*> \verbatim<br>
*>          ERR_BNDS_COMP is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>              sqrt(n) * dlamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * dlamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * dlamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NPARAMS<br>
*> \verbatim<br>
*>          NPARAMS is INTEGER<br>
*>     Specifies the number of parameters set in PARAMS.  If .LE. 0, the<br>
*>     PARAMS array is never referenced and default values are used.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] PARAMS<br>
*> \verbatim<br>
*>          PARAMS is DOUBLE PRECISION array, dimension (NPARAMS)<br>
*>     Specifies algorithm parameters.  If an entry is .LT. 0.0, then<br>
*>     that entry will be filled with default value used for that<br>
*>     parameter.  Only positions up to NPARAMS are accessed; defaults<br>
*>     are used for higher-numbered parameters.<br>
*><br>
*>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative<br>
*>            refinement or not.<br>
*>         Default: 1.0D+0<br>
*>            = 0.0 : No refinement is performed, and no error bounds are<br>
*>                    computed.<br>
*>            = 1.0 : Use the double-precision refinement algorithm,<br>
*>                    possibly with doubled-single computations if the<br>
*>                    compilation environment does not support DOUBLE<br>
*>                    PRECISION.<br>
*>              (other values are reserved for future use)<br>
*><br>
*>       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual<br>
*>            computations allowed for refinement.<br>
*>         Default: 10<br>
*>         Aggressive: Set to 100 to permit convergence using approximate<br>
*>                     factorizations or factorizations other than LU. If<br>
*>                     the factorization uses a technique other than<br>
*>                     Gaussian elimination, the guarantees in<br>
*>                     err_bnds_norm and err_bnds_comp may no longer be<br>
*>                     trustworthy.<br>
*><br>
*>       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code<br>
*>            will attempt to find a solution with small componentwise<br>
*>            relative error in the double-precision algorithm.  Positive<br>
*>            is true, 0.0 is false.<br>
*>         Default: 1.0 (attempt componentwise convergence)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit. The solution to every right-hand side is<br>
*>         guaranteed.<br>
*>       < 0:  If INFO = -i, the i-th argument had an illegal value<br>
*>       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization<br>
*>         has been completed, but the factor U is exactly singular, so<br>
*>         the solution and error bounds could not be computed. RCOND = 0<br>
*>         is returned.<br>
*>       = N+J: The solution corresponding to the Jth right-hand side is<br>
*>         not guaranteed. The solutions corresponding to other right-<br>
*>         hand sides K with K > J may not be guaranteed as well, but<br>
*>         only the first such right-hand side is reported. If a small<br>
*>         componentwise error is not requested (PARAMS(3) = 0.0) then<br>
*>         the Jth right-hand side is the first with a normwise error<br>
*>         bound that is not guaranteed (the smallest J such<br>
*>         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)<br>
*>         the Jth right-hand side is the first with either a normwise or<br>
*>         componentwise error bound that is not guaranteed (the smallest<br>
*>         J such that either ERR_BNDS_NORM(J,1) = 0.0 or<br>
*>         ERR_BNDS_COMP(J,1) = 0.0). See the definition of<br>
*>         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information<br>
*>         about all of the right-hand sides check ERR_BNDS_NORM or<br>
*>         ERR_BNDS_COMP.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgerfsx_(CHARACTER TRANS,CHARACTER EQUED,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGERQ2 computes the RQ factorization of a general rectangular matrix using an unblocked algorithm.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGERQ2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgerq2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgerq2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgerq2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGERQ2( M, N, A, LDA, TAU, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
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
*> DGERQ2 computes an RQ factorization of a real m by n matrix A:<br>
*> A = R * Q.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the m by n matrix A.<br>
*>          On exit, if m <= n, the upper triangle of the subarray<br>
*>          A(1:m,n-m+1:n) contains the m by m upper triangular matrix R;<br>
*>          if m >= n, the elements on and above the (m-n)-th subdiagonal<br>
*>          contain the m by n upper trapezoidal matrix R; the remaining<br>
*>          elements, with the array TAU, represent the orthogonal matrix<br>
*>          Q as a product of elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in<br>
*>  A(m-k+i,1:n-k+i-1), and tau in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgerq2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
/**
*> \brief \b DGERQF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGERQF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgerqf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgerqf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgerqf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
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
*> DGERQF computes an RQ factorization of a real M-by-N matrix A:<br>
*> A = R * Q.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit,<br>
*>          if m <= n, the upper triangle of the subarray<br>
*>          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;<br>
*>          if m >= n, the elements on and above the (m-n)-th subdiagonal<br>
*>          contain the M-by-N upper trapezoidal matrix R;<br>
*>          the remaining elements, with the array TAU, represent the<br>
*>          orthogonal matrix Q as a product of min(m,n) elementary<br>
*>          reflectors (see Further Details).<br>
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
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          The dimension of the array WORK.  LWORK >= max(1,M).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in<br>
*>  A(m-k+i,1:n-k+i-1), and tau in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgerqf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGESC2 solves a system of linear equations using the LU factorization with complete pivoting computed by sgetc2.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGESC2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesc2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesc2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesc2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            LDA, N<br>
*       DOUBLE PRECISION   SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), JPIV( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), RHS( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGESC2 solves a system of linear equations<br>
*><br>
*>           A * X = scale* RHS<br>
*><br>
*> with a general N-by-N matrix A using the LU factorization with<br>
*> complete pivoting computed by DGETC2.<br>
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
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the  LU part of the factorization of the n-by-n<br>
*>          matrix A computed by DGETC2:  A = P * L * U * Q<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1, N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RHS<br>
*> \verbatim<br>
*>          RHS is DOUBLE PRECISION array, dimension (N).<br>
*>          On entry, the right hand side vector b.<br>
*>          On exit, the solution vector X.<br>
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
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is DOUBLE PRECISION<br>
*>          On exit, SCALE contains the scale factor. SCALE is chosen<br>
*>          0 <= SCALE <= 1 to prevent owerflow in the solution.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,<br>
*>     Umea University, S-901 87 Umea, Sweden.<br>
*<br>
*  =====================================================================<br>
*/
	public void dgesc2_(INTEGER N,double[] A,INTEGER LDA,double[] RHS,int[] IPIV,int[] JPIV,DOUBLE SCALE);
/**
*> \brief \b DGESDD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGESDD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesdd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesdd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesdd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT,<br>
*                          WORK, LWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ<br>
*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),<br>
*      $                   VT( LDVT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGESDD computes the singular value decomposition (SVD) of a real<br>
*> M-by-N matrix A, optionally computing the left and right singular<br>
*> vectors.  If singular vectors are desired, it uses a<br>
*> divide-and-conquer algorithm.<br>
*><br>
*> The SVD is written<br>
*><br>
*>      A = U * SIGMA * transpose(V)<br>
*><br>
*> where SIGMA is an M-by-N matrix which is zero except for its<br>
*> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and<br>
*> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA<br>
*> are the singular values of A; they are real and non-negative, and<br>
*> are returned in descending order.  The first min(m,n) columns of<br>
*> U and V are the left and right singular vectors of A.<br>
*><br>
*> Note that the routine returns VT = V**T, not V.<br>
*><br>
*> The divide and conquer algorithm makes very mild assumptions about<br>
*> floating point arithmetic. It will work on machines with a guard<br>
*> digit in add/subtract, or on those binary machines without guard<br>
*> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or<br>
*> Cray-2. It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBZ<br>
*> \verbatim<br>
*>          JOBZ is CHARACTER*1<br>
*>          Specifies options for computing all or part of the matrix U:<br>
*>          = 'A':  all M columns of U and all N rows of V**T are<br>
*>                  returned in the arrays U and VT;<br>
*>          = 'S':  the first min(M,N) columns of U and the first<br>
*>                  min(M,N) rows of V**T are returned in the arrays U<br>
*>                  and VT;<br>
*>          = 'O':  If M >= N, the first N columns of U are overwritten<br>
*>                  on the array A and all rows of V**T are returned in<br>
*>                  the array VT;<br>
*>                  otherwise, all columns of U are returned in the<br>
*>                  array U and the first M rows of V**T are overwritten<br>
*>                  in the array A;<br>
*>          = 'N':  no columns of U or rows of V**T are computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the input matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the input matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit,<br>
*>          if JOBZ = 'O',  A is overwritten with the first N columns<br>
*>                          of U (the left singular vectors, stored<br>
*>                          columnwise) if M >= N;<br>
*>                          A is overwritten with the first M rows<br>
*>                          of V**T (the right singular vectors, stored<br>
*>                          rowwise) otherwise.<br>
*>          if JOBZ .ne. 'O', the contents of A are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The singular values of A, sorted so that S(i) >= S(i+1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is DOUBLE PRECISION array, dimension (LDU,UCOL)<br>
*>          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;<br>
*>          UCOL = min(M,N) if JOBZ = 'S'.<br>
*>          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M<br>
*>          orthogonal matrix U;<br>
*>          if JOBZ = 'S', U contains the first min(M,N) columns of U<br>
*>          (the left singular vectors, stored columnwise);<br>
*>          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>          The leading dimension of the array U.  LDU >= 1; if<br>
*>          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VT<br>
*> \verbatim<br>
*>          VT is DOUBLE PRECISION array, dimension (LDVT,N)<br>
*>          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the<br>
*>          N-by-N orthogonal matrix V**T;<br>
*>          if JOBZ = 'S', VT contains the first min(M,N) rows of<br>
*>          V**T (the right singular vectors, stored rowwise);<br>
*>          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>          The leading dimension of the array VT.  LDVT >= 1;<br>
*>          if JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;<br>
*>          if JOBZ = 'S', LDVT >= min(M,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= 1.<br>
*>          If LWORK = -1, a workspace query is assumed.  The optimal<br>
*>          size for the WORK array is calculated and stored in WORK(1),<br>
*>          and no other work except argument checking is performed.<br>
*><br>
*>          Let mx = max(M,N) and mn = min(M,N).<br>
*>          If JOBZ = 'N', LWORK >= 3*mn + max( mx, 7*mn ).<br>
*>          If JOBZ = 'O', LWORK >= 3*mn + max( mx, 5*mn*mn + 4*mn ).<br>
*>          If JOBZ = 'S', LWORK >= 4*mn*mn + 7*mn.<br>
*>          If JOBZ = 'A', LWORK >= 4*mn*mn + 6*mn + mx.<br>
*>          These are not tight minimums in all cases; see comments inside code.<br>
*>          For good performance, LWORK should generally be larger;<br>
*>          a query is recommended.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (8*min(M,N))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  DBDSDC did not converge, updating process failed.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsing<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*> @precisions fortran d -> s<br>
*  =====================================================================<br>
*/
	public void dgesdd_(CHARACTER JOBZ,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] S,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> DGESV computes the solution to system of linear equations A * X = B for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGESV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGESV computes the solution to a real system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.<br>
*><br>
*> The LU decomposition with partial pivoting and row interchanges is<br>
*> used to factor A as<br>
*>    A = P * L * U,<br>
*> where P is a permutation matrix, L is unit lower triangular, and U is<br>
*> upper triangular.  The factored form of A is then used to solve the<br>
*> system of equations A * X = B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of linear equations, i.e., the order of the<br>
*>          matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the N-by-N coefficient matrix A.<br>
*>          On exit, the factors L and U from the factorization<br>
*>          A = P*L*U; the unit diagonal elements of L are not stored.<br>
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
*>          The pivot indices that define the permutation matrix P;<br>
*>          row i of the matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the N-by-NRHS matrix of right hand side matrix B.<br>
*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization<br>
*>                has been completed, but the factor U is exactly<br>
*>                singular, so the solution could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgesv_(INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> DGESVD computes the singular value decomposition (SVD) for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGESVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBU, JOBVT<br>
*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),<br>
*      $                   VT( LDVT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGESVD computes the singular value decomposition (SVD) of a real<br>
*> M-by-N matrix A, optionally computing the left and/or right singular<br>
*> vectors. The SVD is written<br>
*><br>
*>      A = U * SIGMA * transpose(V)<br>
*><br>
*> where SIGMA is an M-by-N matrix which is zero except for its<br>
*> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and<br>
*> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA<br>
*> are the singular values of A; they are real and non-negative, and<br>
*> are returned in descending order.  The first min(m,n) columns of<br>
*> U and V are the left and right singular vectors of A.<br>
*><br>
*> Note that the routine returns V**T, not V.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBU<br>
*> \verbatim<br>
*>          JOBU is CHARACTER*1<br>
*>          Specifies options for computing all or part of the matrix U:<br>
*>          = 'A':  all M columns of U are returned in array U:<br>
*>          = 'S':  the first min(m,n) columns of U (the left singular<br>
*>                  vectors) are returned in the array U;<br>
*>          = 'O':  the first min(m,n) columns of U (the left singular<br>
*>                  vectors) are overwritten on the array A;<br>
*>          = 'N':  no columns of U (no left singular vectors) are<br>
*>                  computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVT<br>
*> \verbatim<br>
*>          JOBVT is CHARACTER*1<br>
*>          Specifies options for computing all or part of the matrix<br>
*>          V**T:<br>
*>          = 'A':  all N rows of V**T are returned in the array VT;<br>
*>          = 'S':  the first min(m,n) rows of V**T (the right singular<br>
*>                  vectors) are returned in the array VT;<br>
*>          = 'O':  the first min(m,n) rows of V**T (the right singular<br>
*>                  vectors) are overwritten on the array A;<br>
*>          = 'N':  no rows of V**T (no right singular vectors) are<br>
*>                  computed.<br>
*><br>
*>          JOBVT and JOBU cannot both be 'O'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the input matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the input matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit,<br>
*>          if JOBU = 'O',  A is overwritten with the first min(m,n)<br>
*>                          columns of U (the left singular vectors,<br>
*>                          stored columnwise);<br>
*>          if JOBVT = 'O', A is overwritten with the first min(m,n)<br>
*>                          rows of V**T (the right singular vectors,<br>
*>                          stored rowwise);<br>
*>          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A<br>
*>                          are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The singular values of A, sorted so that S(i) >= S(i+1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is DOUBLE PRECISION array, dimension (LDU,UCOL)<br>
*>          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.<br>
*>          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;<br>
*>          if JOBU = 'S', U contains the first min(m,n) columns of U<br>
*>          (the left singular vectors, stored columnwise);<br>
*>          if JOBU = 'N' or 'O', U is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>          The leading dimension of the array U.  LDU >= 1; if<br>
*>          JOBU = 'S' or 'A', LDU >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VT<br>
*> \verbatim<br>
*>          VT is DOUBLE PRECISION array, dimension (LDVT,N)<br>
*>          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix<br>
*>          V**T;<br>
*>          if JOBVT = 'S', VT contains the first min(m,n) rows of<br>
*>          V**T (the right singular vectors, stored rowwise);<br>
*>          if JOBVT = 'N' or 'O', VT is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>          The leading dimension of the array VT.  LDVT >= 1; if<br>
*>          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;<br>
*>          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged<br>
*>          superdiagonal elements of an upper bidiagonal matrix B<br>
*>          whose diagonal is in S (not necessarily sorted). B<br>
*>          satisfies A = U * B * VT, so it has the same singular values<br>
*>          as A, and singular vectors related by U and VT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):<br>
*>             - PATH 1  (M much larger than N, JOBU='N') <br>
*>             - PATH 1t (N much larger than M, JOBVT='N')<br>
*>          LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths<br>
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
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if DBDSQR did not converge, INFO specifies how many<br>
*>                superdiagonals of an intermediate bidiagonal form B<br>
*>                did not converge to zero. See the description of WORK<br>
*>                above for details.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsing<br>
*<br>
*  =====================================================================<br>
*/
	public void dgesvd_(CHARACTER JOBU,CHARACTER JOBVT,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] S,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGESVDX computes the singular value decomposition (SVD) for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGESVDX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvdx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvdx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvdx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*     SUBROUTINE DGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, <br>
*    $                    IL, IU, NS, S, U, LDU, VT, LDVT, WORK, <br>
*    $                    LWORK, IWORK, INFO )<br>
*      <br>
*<br>
*     .. Scalar Arguments ..<br>
*      CHARACTER          JOBU, JOBVT, RANGE<br>
*      INTEGER            IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS<br>
*      DOUBLE PRECISION   VL, VU<br>
*     ..<br>
*     .. Array Arguments ..<br>
*     INTEGER            IWORK( * )<br>
*     DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),<br>
*    $                   VT( LDVT, * ), WORK( * )<br>
*     ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>  DGESVDX computes the singular value decomposition (SVD) of a real<br>
*>  M-by-N matrix A, optionally computing the left and/or right singular<br>
*>  vectors. The SVD is written<br>
*> <br>
*>      A = U * SIGMA * transpose(V)<br>
*> <br>
*>  where SIGMA is an M-by-N matrix which is zero except for its<br>
*>  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and<br>
*>  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA<br>
*>  are the singular values of A; they are real and non-negative, and<br>
*>  are returned in descending order.  The first min(m,n) columns of<br>
*>  U and V are the left and right singular vectors of A.<br>
*> <br>
*>  DGESVDX uses an eigenvalue problem for obtaining the SVD, which <br>
*>  allows for the computation of a subset of singular values and <br>
*>  vectors. See DBDSVDX for details.<br>
*> <br>
*>  Note that the routine returns V**T, not V.<br>
*> \endverbatim<br>
*   <br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBU<br>
*> \verbatim<br>
*>          JOBU is CHARACTER*1<br>
*>          Specifies options for computing all or part of the matrix U:<br>
*>          = 'V':  the first min(m,n) columns of U (the left singular<br>
*>                  vectors) or as specified by RANGE are returned in <br>
*>                  the array U;<br>
*>          = 'N':  no columns of U (no left singular vectors) are<br>
*>                  computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVT<br>
*> \verbatim<br>
*>          JOBVT is CHARACTER*1<br>
*>           Specifies options for computing all or part of the matrix<br>
*>           V**T:<br>
*>           = 'V':  the first min(m,n) rows of V**T (the right singular<br>
*>                   vectors) or as specified by RANGE are returned in <br>
*>                   the array VT;<br>
*>           = 'N':  no rows of V**T (no right singular vectors) are<br>
*>                   computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RANGE<br>
*> \verbatim<br>
*>          RANGE is CHARACTER*1<br>
*>          = 'A': all singular values will be found.<br>
*>          = 'V': all singular values in the half-open interval (VL,VU]<br>
*>                 will be found.<br>
*>          = 'I': the IL-th through IU-th singular values will be found. <br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the input matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the input matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the contents of A are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for singular values. VU > VL.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is DOUBLE PRECISION<br>
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
*>          The total number of singular values found,  <br>
*>          0 <= NS <= min(M,N).<br>
*>          If RANGE = 'A', NS = min(M,N); if RANGE = 'I', NS = IU-IL+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The singular values of A, sorted so that S(i) >= S(i+1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is DOUBLE PRECISION array, dimension (LDU,UCOL)<br>
*>          If JOBU = 'V', U contains columns of U (the left singular <br>
*>          vectors, stored columnwise) as specified by RANGE; if <br>
*>          JOBU = 'N', U is not referenced.<br>
*>          Note: The user must ensure that UCOL >= NS; if RANGE = 'V', <br>
*>          the exact value of NS is not known in advance and an upper<br>
*>          bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>          The leading dimension of the array U.  LDU >= 1; if<br>
*>          JOBU = 'V', LDU >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VT<br>
*> \verbatim<br>
*>          VT is DOUBLE PRECISION array, dimension (LDVT,N)<br>
*>          If JOBVT = 'V', VT contains the rows of V**T (the right singular <br>
*>          vectors, stored rowwise) as specified by RANGE; if JOBVT = 'N', <br>
*>          VT is not referenced.<br>
*>          Note: The user must ensure that LDVT >= NS; if RANGE = 'V', <br>
*>          the exact value of NS is not known in advance and an upper <br>
*>          bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVT<br>
*> \verbatim<br>
*>          LDVT is INTEGER<br>
*>          The leading dimension of the array VT.  LDVT >= 1; if<br>
*>          JOBVT = 'V', LDVT >= NS (see above).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          LWORK >= MAX(1,MIN(M,N)*(MIN(M,N)+4)) for the paths (see <br>
*>          comments inside the code):<br>
*>             - PATH 1  (M much larger than N) <br>
*>             - PATH 1t (N much larger than M)<br>
*>          LWORK >= MAX(1,MIN(M,N)*2+MAX(M,N)) for the other paths.<br>
*>          For good performance, LWORK should generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (12*MIN(M,N))<br>
*>          If INFO = 0, the first NS elements of IWORK are zero. If INFO > 0, <br>
*>          then IWORK contains the indices of the eigenvectors that failed <br>
*>          to converge in DBDSVDX/DSTEVX.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>     INFO is INTEGER<br>
*>           = 0:  successful exit<br>
*>           < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>           > 0:  if INFO = i, then i eigenvectors failed to converge<br>
*>                 in DBDSVDX/DSTEVX.<br>
*>                 if INFO = N*2 + 1, an internal error occurred in<br>
*>                 DBDSVDX<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsing<br>
*<br>
*  =====================================================================<br>
*/
	public void dgesvdx_(CHARACTER JOBU,CHARACTER JOBVT,CHARACTER RANGE,INTEGER M,INTEGER N,double[] A,INTEGER LDA,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,INTEGER NS,double[] S,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGESVJ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGESVJ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvj.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvj.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvj.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V,<br>
*                          LDV, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N<br>
*       CHARACTER*1        JOBA, JOBU, JOBV<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), SVA( N ), V( LDV, * ),<br>
*      $                   WORK( LWORK )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGESVJ computes the singular value decomposition (SVD) of a real<br>
*> M-by-N matrix A, where M >= N. The SVD of A is written as<br>
*>                                    [++]   [xx]   [x0]   [xx]<br>
*>              A = U * SIGMA * V^t,  [++] = [xx] * [ox] * [xx]<br>
*>                                    [++]   [xx]<br>
*> where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal<br>
*> matrix, and V is an N-by-N orthogonal matrix. The diagonal elements<br>
*> of SIGMA are the singular values of A. The columns of U and V are the<br>
*> left and the right singular vectors of A, respectively.<br>
*> DGESVJ can sometimes compute tiny singular values and their singular vectors much<br>
*> more accurately than other SVD routines, see below under Further Details.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBA<br>
*> \verbatim<br>
*>          JOBA is CHARACTER* 1<br>
*>          Specifies the structure of A.<br>
*>          = 'L': The input matrix A is lower triangular;<br>
*>          = 'U': The input matrix A is upper triangular;<br>
*>          = 'G': The input matrix A is general M-by-N matrix, M >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBU<br>
*> \verbatim<br>
*>          JOBU is CHARACTER*1<br>
*>          Specifies whether to compute the left singular vectors<br>
*>          (columns of U):<br>
*>          = 'U': The left singular vectors corresponding to the nonzero<br>
*>                 singular values are computed and returned in the leading<br>
*>                 columns of A. See more details in the description of A.<br>
*>                 The default numerical orthogonality threshold is set to<br>
*>                 approximately TOL=CTOL*EPS, CTOL=DSQRT(M), EPS=DLAMCH('E').<br>
*>          = 'C': Analogous to JOBU='U', except that user can control the<br>
*>                 level of numerical orthogonality of the computed left<br>
*>                 singular vectors. TOL can be set to TOL = CTOL*EPS, where<br>
*>                 CTOL is given on input in the array WORK.<br>
*>                 No CTOL smaller than ONE is allowed. CTOL greater<br>
*>                 than 1 / EPS is meaningless. The option 'C'<br>
*>                 can be used if M*EPS is satisfactory orthogonality<br>
*>                 of the computed left singular vectors, so CTOL=M could<br>
*>                 save few sweeps of Jacobi rotations.<br>
*>                 See the descriptions of A and WORK(1).<br>
*>          = 'N': The matrix U is not computed. However, see the<br>
*>                 description of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV<br>
*> \verbatim<br>
*>          JOBV is CHARACTER*1<br>
*>          Specifies whether to compute the right singular vectors, that<br>
*>          is, the matrix V:<br>
*>          = 'V' : the matrix V is computed and returned in the array V<br>
*>          = 'A' : the Jacobi rotations are applied to the MV-by-N<br>
*>                  array V. In other words, the right singular vector<br>
*>                  matrix V is not computed explicitly, instead it is<br>
*>                  applied to an MV-by-N matrix initially stored in the<br>
*>                  first MV rows of V.<br>
*>          = 'N' : the matrix V is not computed and the array V is not<br>
*>                  referenced<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the input matrix A. 1/DLAMCH('E') > M >= 0.  <br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the input matrix A.<br>
*>          M >= N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit :<br>
*>          If JOBU .EQ. 'U' .OR. JOBU .EQ. 'C' :<br>
*>                 If INFO .EQ. 0 :<br>
*>                 RANKA orthonormal columns of U are returned in the<br>
*>                 leading RANKA columns of the array A. Here RANKA <= N<br>
*>                 is the number of computed singular values of A that are<br>
*>                 above the underflow threshold DLAMCH('S'). The singular<br>
*>                 vectors corresponding to underflowed or zero singular<br>
*>                 values are not computed. The value of RANKA is returned<br>
*>                 in the array WORK as RANKA=NINT(WORK(2)). Also see the<br>
*>                 descriptions of SVA and WORK. The computed columns of U<br>
*>                 are mutually numerically orthogonal up to approximately<br>
*>                 TOL=DSQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU.EQ.'C'),<br>
*>                 see the description of JOBU.<br>
*>                 If INFO .GT. 0 :<br>
*>                 the procedure DGESVJ did not converge in the given number<br>
*>                 of iterations (sweeps). In that case, the computed<br>
*>                 columns of U may not be orthogonal up to TOL. The output<br>
*>                 U (stored in A), SIGMA (given by the computed singular<br>
*>                 values in SVA(1:N)) and V is still a decomposition of the<br>
*>                 input matrix A in the sense that the residual<br>
*>                 ||A-SCALE*U*SIGMA*V^T||_2 / ||A||_2 is small.<br>
*><br>
*>          If JOBU .EQ. 'N' :<br>
*>                 If INFO .EQ. 0 :<br>
*>                 Note that the left singular vectors are 'for free' in the<br>
*>                 one-sided Jacobi SVD algorithm. However, if only the<br>
*>                 singular values are needed, the level of numerical<br>
*>                 orthogonality of U is not an issue and iterations are<br>
*>                 stopped when the columns of the iterated matrix are<br>
*>                 numerically orthogonal up to approximately M*EPS. Thus,<br>
*>                 on exit, A contains the columns of U scaled with the<br>
*>                 corresponding singular values.<br>
*>                 If INFO .GT. 0 :<br>
*>                 the procedure DGESVJ did not converge in the given number<br>
*>                 of iterations (sweeps).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SVA<br>
*> \verbatim<br>
*>          SVA is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit :<br>
*>          If INFO .EQ. 0 :<br>
*>          depending on the value SCALE = WORK(1), we have:<br>
*>                 If SCALE .EQ. ONE :<br>
*>                 SVA(1:N) contains the computed singular values of A.<br>
*>                 During the computation SVA contains the Euclidean column<br>
*>                 norms of the iterated matrices in the array A.<br>
*>                 If SCALE .NE. ONE :<br>
*>                 The singular values of A are SCALE*SVA(1:N), and this<br>
*>                 factored representation is due to the fact that some of the<br>
*>                 singular values of A might underflow or overflow.<br>
*>          If INFO .GT. 0 :<br>
*>          the procedure DGESVJ did not converge in the given number of<br>
*>          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MV<br>
*> \verbatim<br>
*>          MV is INTEGER<br>
*>          If JOBV .EQ. 'A', then the product of Jacobi rotations in DGESVJ<br>
*>          is applied to the first MV rows of V. See the description of JOBV.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension (LDV,N)<br>
*>          If JOBV = 'V', then V contains on exit the N-by-N matrix of<br>
*>                         the right singular vectors;<br>
*>          If JOBV = 'A', then V contains the product of the computed right<br>
*>                         singular vector matrix and the initial matrix in<br>
*>                         the array V.<br>
*>          If JOBV = 'N', then V is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V, LDV .GE. 1.<br>
*>          If JOBV .EQ. 'V', then LDV .GE. max(1,N).<br>
*>          If JOBV .EQ. 'A', then LDV .GE. max(1,MV) .<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension max(4,M+N).<br>
*>          On entry :<br>
*>          If JOBU .EQ. 'C' :<br>
*>          WORK(1) = CTOL, where CTOL defines the threshold for convergence.<br>
*>                    The process stops if all columns of A are mutually<br>
*>                    orthogonal up to CTOL*EPS, EPS=DLAMCH('E').<br>
*>                    It is required that CTOL >= ONE, i.e. it is not<br>
*>                    allowed to force the routine to obtain orthogonality<br>
*>                    below EPS.<br>
*>          On exit :<br>
*>          WORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N)<br>
*>                    are the computed singular values of A.<br>
*>                    (See description of SVA().)<br>
*>          WORK(2) = NINT(WORK(2)) is the number of the computed nonzero<br>
*>                    singular values.<br>
*>          WORK(3) = NINT(WORK(3)) is the number of the computed singular<br>
*>                    values that are larger than the underflow threshold.<br>
*>          WORK(4) = NINT(WORK(4)) is the number of sweeps of Jacobi<br>
*>                    rotations needed for numerical convergence.<br>
*>          WORK(5) = max_{i.NE.j} |COS(A(:,i),A(:,j))| in the last sweep.<br>
*>                    This is useful information in cases when DGESVJ did<br>
*>                    not converge, as it can be used to estimate whether<br>
*>                    the output is stil useful and for post festum analysis.<br>
*>          WORK(6) = the largest absolute value over all sines of the<br>
*>                    Jacobi rotation angles in the last sweep. It can be<br>
*>                    useful for a post festum analysis.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          length of WORK, WORK >= MAX(6,M+N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0 : successful exit.<br>
*>          < 0 : if INFO = -i, then the i-th argument had an illegal value<br>
*>          > 0 : DGESVJ did not converge in the maximal allowed number (30)<br>
*>                of sweeps. The output may still be useful. See the<br>
*>                description of WORK.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The orthogonal N-by-N matrix V is obtained as a product of Jacobi plane<br>
*>  rotations. The rotations are implemented as fast scaled rotations of<br>
*>  Anda and Park [1]. In the case of underflow of the Jacobi angle, a<br>
*>  modified Jacobi transformation of Drmac [4] is used. Pivot strategy uses<br>
*>  column interchanges of de Rijk [2]. The relative accuracy of the computed<br>
*>  singular values and the accuracy of the computed singular vectors (in<br>
*>  angle metric) is as guaranteed by the theory of Demmel and Veselic [3].<br>
*>  The condition number that determines the accuracy in the full rank case<br>
*>  is essentially min_{D=diag} kappa(A*D), where kappa(.) is the<br>
*>  spectral condition number. The best performance of this Jacobi SVD<br>
*>  procedure is achieved if used in an  accelerated version of Drmac and<br>
*>  Veselic [5,6], and it is the kernel routine in the SIGMA library [7].<br>
*>  Some tunning parameters (marked with [TP]) are available for the<br>
*>  implementer.<br>
*>  The computational range for the nonzero singular values is the  machine<br>
*>  number interval ( UNDERFLOW , OVERFLOW ). In extreme cases, even<br>
*>  denormalized singular values can be computed with the corresponding<br>
*>  gradual loss of accurate digits.<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>  ============<br>
*><br>
*>  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*> \verbatim<br>
*><br>
*> [1] A. A. Anda and H. Park: Fast plane rotations with dynamic scaling.<br>
*>     SIAM J. matrix Anal. Appl., Vol. 15 (1994), pp. 162-174.<br>
*> [2] P. P. M. De Rijk: A one-sided Jacobi algorithm for computing the<br>
*>     singular value decomposition on a vector computer.<br>
*>     SIAM J. Sci. Stat. Comp., Vol. 10 (1998), pp. 359-371.<br>
*> [3] J. Demmel and K. Veselic: Jacobi method is more accurate than QR.<br>
*> [4] Z. Drmac: Implementation of Jacobi rotations for accurate singular<br>
*>     value computation in floating point arithmetic.<br>
*>     SIAM J. Sci. Comp., Vol. 18 (1997), pp. 1200-1222.<br>
*> [5] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I.<br>
*>     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342.<br>
*>     LAPACK Working note 169.<br>
*> [6] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II.<br>
*>     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362.<br>
*>     LAPACK Working note 170.<br>
*> [7] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV,<br>
*>     QSVD, (H,K)-SVD computations.<br>
*>     Department of Mathematics, University of Zagreb, 2008.<br>
*> \endverbatim<br>
*<br>
*>  \par Bugs, examples and comments:<br>
*   =================================<br>
*><br>
*> \verbatim<br>
*>  ===========================<br>
*>  Please report all bugs and send interesting test examples and comments to<br>
*>  drmac@math.hr. Thank you.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgesvj_(char[] JOBA,char[] JOBU,char[] JOBV,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] SVA,INTEGER MV,double[] V,INTEGER LDV,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGESVX computes the solution to system of linear equations A * X = B for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGESVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
*                          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,<br>
*                          WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, TRANS<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   BERR( * ), C( * ), FERR( * ), R( * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGESVX uses the LU factorization to compute the solution to a real<br>
*> system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.<br>
*><br>
*> Error bounds on the solution and a condition estimate are also<br>
*> provided.<br>
*> \endverbatim<br>
*<br>
*> \par Description:<br>
*  =================<br>
*><br>
*> \verbatim<br>
*><br>
*> The following steps are performed:<br>
*><br>
*> 1. If FACT = 'E', real scaling factors are computed to equilibrate<br>
*>    the system:<br>
*>       TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B<br>
*>       TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B<br>
*>       TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B<br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')<br>
*>    or diag(C)*B (if TRANS = 'T' or 'C').<br>
*><br>
*> 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the<br>
*>    matrix A (after equilibration if FACT = 'E') as<br>
*>       A = P * L * U,<br>
*>    where P is a permutation matrix, L is a unit lower triangular<br>
*>    matrix, and U is upper triangular.<br>
*><br>
*> 3. If some U(i,i)=0, so that U is exactly singular, then the routine<br>
*>    returns with INFO = i. Otherwise, the factored form of A is used<br>
*>    to estimate the condition number of the matrix A.  If the<br>
*>    reciprocal of the condition number is less than machine precision,<br>
*>    INFO = N+1 is returned as a warning, but the routine still goes on<br>
*>    to solve for X and compute error bounds as described below.<br>
*><br>
*> 4. The system of equations is solved for X using the factored form<br>
*>    of A.<br>
*><br>
*> 5. Iterative refinement is applied to improve the computed solution<br>
*>    matrix and calculate error bounds and backward error estimates<br>
*>    for it.<br>
*><br>
*> 6. If equilibration was used, the matrix X is premultiplied by<br>
*>    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so<br>
*>    that it solves the original system before equilibration.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] FACT<br>
*> \verbatim<br>
*>          FACT is CHARACTER*1<br>
*>          Specifies whether or not the factored form of the matrix A is<br>
*>          supplied on entry, and if not, whether the matrix A should be<br>
*>          equilibrated before it is factored.<br>
*>          = 'F':  On entry, AF and IPIV contain the factored form of A.<br>
*>                  If EQUED is not 'N', the matrix A has been<br>
*>                  equilibrated with scaling factors given by R and C.<br>
*>                  A, AF, and IPIV are not modified.<br>
*>          = 'N':  The matrix A will be copied to AF and factored.<br>
*>          = 'E':  The matrix A will be equilibrated if necessary, then<br>
*>                  copied to AF and factored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of linear equations, i.e., the order of the<br>
*>          matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is<br>
*>          not 'N', then A must have been equilibrated by the scaling<br>
*>          factors in R and/or C.  A is not modified if FACT = 'F' or<br>
*>          'N', or if FACT = 'E' and EQUED = 'N' on exit.<br>
*><br>
*>          On exit, if EQUED .ne. 'N', A is scaled as follows:<br>
*>          EQUED = 'R':  A := diag(R) * A<br>
*>          EQUED = 'C':  A := A * diag(C)<br>
*>          EQUED = 'B':  A := diag(R) * A * diag(C).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AF<br>
*> \verbatim<br>
*>          AF is DOUBLE PRECISION array, dimension (LDAF,N)<br>
*>          If FACT = 'F', then AF is an input argument and on entry<br>
*>          contains the factors L and U from the factorization<br>
*>          A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then<br>
*>          AF is the factored form of the equilibrated matrix A.<br>
*><br>
*>          If FACT = 'N', then AF is an output argument and on exit<br>
*>          returns the factors L and U from the factorization A = P*L*U<br>
*>          of the original matrix A.<br>
*><br>
*>          If FACT = 'E', then AF is an output argument and on exit<br>
*>          returns the factors L and U from the factorization A = P*L*U<br>
*>          of the equilibrated matrix A (see the description of A for<br>
*>          the form of the equilibrated matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>          The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          If FACT = 'F', then IPIV is an input argument and on entry<br>
*>          contains the pivot indices from the factorization A = P*L*U<br>
*>          as computed by DGETRF; row i of the matrix was interchanged<br>
*>          with row IPIV(i).<br>
*><br>
*>          If FACT = 'N', then IPIV is an output argument and on exit<br>
*>          contains the pivot indices from the factorization A = P*L*U<br>
*>          of the original matrix A.<br>
*><br>
*>          If FACT = 'E', then IPIV is an output argument and on exit<br>
*>          contains the pivot indices from the factorization A = P*L*U<br>
*>          of the equilibrated matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies the form of equilibration that was done.<br>
*>          = 'N':  No equilibration (always true if FACT = 'N').<br>
*>          = 'R':  Row equilibration, i.e., A has been premultiplied by<br>
*>                  diag(R).<br>
*>          = 'C':  Column equilibration, i.e., A has been postmultiplied<br>
*>                  by diag(C).<br>
*>          = 'B':  Both row and column equilibration, i.e., A has been<br>
*>                  replaced by diag(R) * A * diag(C).<br>
*>          EQUED is an input argument if FACT = 'F'; otherwise, it is an<br>
*>          output argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (N)<br>
*>          The row scale factors for A.  If EQUED = 'R' or 'B', A is<br>
*>          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R<br>
*>          is not accessed.  R is an input argument if FACT = 'F';<br>
*>          otherwise, R is an output argument.  If FACT = 'F' and<br>
*>          EQUED = 'R' or 'B', each element of R must be positive.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>          The column scale factors for A.  If EQUED = 'C' or 'B', A is<br>
*>          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C<br>
*>          is not accessed.  C is an input argument if FACT = 'F';<br>
*>          otherwise, C is an output argument.  If FACT = 'F' and<br>
*>          EQUED = 'C' or 'B', each element of C must be positive.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the N-by-NRHS right hand side matrix B.<br>
*>          On exit,<br>
*>          if EQUED = 'N', B is not modified;<br>
*>          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by<br>
*>          diag(R)*B;<br>
*>          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is<br>
*>          overwritten by diag(C)*B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X<br>
*>          to the original system of equations.  Note that A and B are<br>
*>          modified on exit if EQUED .ne. 'N', and the solution to the<br>
*>          equilibrated system is inv(diag(C))*X if TRANS = 'N' and<br>
*>          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'<br>
*>          and EQUED = 'R' or 'B'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          The estimate of the reciprocal condition number of the matrix<br>
*>          A after equilibration (if done).  If RCOND is less than the<br>
*>          machine precision (in particular, if RCOND = 0), the matrix<br>
*>          is singular to working precision.  This condition is<br>
*>          indicated by a return code of INFO > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The estimated forward error bound for each solution vector<br>
*>          X(j) (the j-th column of the solution matrix X).<br>
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)<br>
*>          is an estimated upper bound for the magnitude of the largest<br>
*>          element in (X(j) - XTRUE) divided by the magnitude of the<br>
*>          largest element in X(j).  The estimate is as reliable as<br>
*>          the estimate for RCOND, and is almost always a slight<br>
*>          overestimate of the true error.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (4*N)<br>
*>          On exit, WORK(1) contains the reciprocal pivot growth<br>
*>          factor norm(A)/norm(U). The "max absolute element" norm is<br>
*>          used. If WORK(1) is much less than 1, then the stability<br>
*>          of the LU factorization of the (equilibrated) matrix A<br>
*>          could be poor. This also means that the solution X, condition<br>
*>          estimator RCOND, and forward error bound FERR could be<br>
*>          unreliable. If factorization fails with 0<INFO<=N, then<br>
*>          WORK(1) contains the reciprocal pivot growth factor for the<br>
*>          leading INFO columns of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is<br>
*>                <= N:  U(i,i) is exactly zero.  The factorization has<br>
*>                       been completed, but the factor U is exactly<br>
*>                       singular, so the solution and error bounds<br>
*>                       could not be computed. RCOND = 0 is returned.<br>
*>                = N+1: U is nonsingular, but RCOND is less than machine<br>
*>                       precision, meaning that the matrix is singular<br>
*>                       to working precision.  Nevertheless, the<br>
*>                       solution and error bounds are computed because<br>
*>                       there are a number of situations where the<br>
*>                       computed solution can be more accurate than the<br>
*>                       value of RCOND would suggest.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgesvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> DGESVXX computes the solution to system of linear equations A * X = B for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGESVXX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvxx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvxx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvxx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGESVXX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
*                           EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW,<br>
*                           BERR, N_ERR_BNDS, ERR_BNDS_NORM,<br>
*                           ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK,<br>
*                           INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, TRANS<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS<br>
*       DOUBLE PRECISION   RCOND, RPVGRW<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   X( LDX , * ),WORK( * )<br>
*       DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ),<br>
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
*>    DGESVXX uses the LU factorization to compute the solution to a<br>
*>    double precision system of linear equations  A * X = B,  where A is an<br>
*>    N-by-N matrix and X and B are N-by-NRHS matrices.<br>
*><br>
*>    If requested, both normwise and maximum componentwise error bounds<br>
*>    are returned. DGESVXX will return a solution with a tiny<br>
*>    guaranteed error (O(eps) where eps is the working machine<br>
*>    precision) unless the matrix is very ill-conditioned, in which<br>
*>    case a warning is returned. Relevant condition numbers also are<br>
*>    calculated and returned.<br>
*><br>
*>    DGESVXX accepts user-provided factorizations and equilibration<br>
*>    factors; see the definitions of the FACT and EQUED options.<br>
*>    Solving with refinement and using a factorization from a previous<br>
*>    DGESVXX call will also produce a solution with either O(eps)<br>
*>    errors or warnings, but we cannot make that claim for general<br>
*>    user-provided factorizations and equilibration factors if they<br>
*>    differ from what DGESVXX would itself produce.<br>
*> \endverbatim<br>
*<br>
*> \par Description:<br>
*  =================<br>
*><br>
*> \verbatim<br>
*><br>
*>    The following steps are performed:<br>
*><br>
*>    1. If FACT = 'E', double precision scaling factors are computed to equilibrate<br>
*>    the system:<br>
*><br>
*>      TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B<br>
*>      TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B<br>
*>      TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B<br>
*><br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')<br>
*>    or diag(C)*B (if TRANS = 'T' or 'C').<br>
*><br>
*>    2. If FACT = 'N' or 'E', the LU decomposition is used to factor<br>
*>    the matrix A (after equilibration if FACT = 'E') as<br>
*><br>
*>      A = P * L * U,<br>
*><br>
*>    where P is a permutation matrix, L is a unit lower triangular<br>
*>    matrix, and U is upper triangular.<br>
*><br>
*>    3. If some U(i,i)=0, so that U is exactly singular, then the<br>
*>    routine returns with INFO = i. Otherwise, the factored form of A<br>
*>    is used to estimate the condition number of the matrix A (see<br>
*>    argument RCOND). If the reciprocal of the condition number is less<br>
*>    than machine precision, the routine still goes on to solve for X<br>
*>    and compute error bounds as described below.<br>
*><br>
*>    4. The system of equations is solved for X using the factored form<br>
*>    of A.<br>
*><br>
*>    5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),<br>
*>    the routine will use iterative refinement to try to get a small<br>
*>    error and error bounds.  Refinement calculates the residual to at<br>
*>    least twice the working precision.<br>
*><br>
*>    6. If equilibration was used, the matrix X is premultiplied by<br>
*>    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so<br>
*>    that it solves the original system before equilibration.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \verbatim<br>
*>     Some optional parameters are bundled in the PARAMS array.  These<br>
*>     settings determine how refinement is performed, but often the<br>
*>     defaults are acceptable.  If the defaults are acceptable, users<br>
*>     can pass NPARAMS = 0 which prevents the source code from accessing<br>
*>     the PARAMS argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in] FACT<br>
*> \verbatim<br>
*>          FACT is CHARACTER*1<br>
*>     Specifies whether or not the factored form of the matrix A is<br>
*>     supplied on entry, and if not, whether the matrix A should be<br>
*>     equilibrated before it is factored.<br>
*>       = 'F':  On entry, AF and IPIV contain the factored form of A.<br>
*>               If EQUED is not 'N', the matrix A has been<br>
*>               equilibrated with scaling factors given by R and C.<br>
*>               A, AF, and IPIV are not modified.<br>
*>       = 'N':  The matrix A will be copied to AF and factored.<br>
*>       = 'E':  The matrix A will be equilibrated if necessary, then<br>
*>               copied to AF and factored.<br>
*> \endverbatim<br>
*><br>
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
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>     The number of right hand sides, i.e., the number of columns<br>
*>     of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>     On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is<br>
*>     not 'N', then A must have been equilibrated by the scaling<br>
*>     factors in R and/or C.  A is not modified if FACT = 'F' or<br>
*>     'N', or if FACT = 'E' and EQUED = 'N' on exit.<br>
*><br>
*>     On exit, if EQUED .ne. 'N', A is scaled as follows:<br>
*>     EQUED = 'R':  A := diag(R) * A<br>
*>     EQUED = 'C':  A := A * diag(C)<br>
*>     EQUED = 'B':  A := diag(R) * A * diag(C).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>     The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AF<br>
*> \verbatim<br>
*>          AF is DOUBLE PRECISION array, dimension (LDAF,N)<br>
*>     If FACT = 'F', then AF is an input argument and on entry<br>
*>     contains the factors L and U from the factorization<br>
*>     A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then<br>
*>     AF is the factored form of the equilibrated matrix A.<br>
*><br>
*>     If FACT = 'N', then AF is an output argument and on exit<br>
*>     returns the factors L and U from the factorization A = P*L*U<br>
*>     of the original matrix A.<br>
*><br>
*>     If FACT = 'E', then AF is an output argument and on exit<br>
*>     returns the factors L and U from the factorization A = P*L*U<br>
*>     of the equilibrated matrix A (see the description of A for<br>
*>     the form of the equilibrated matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>     If FACT = 'F', then IPIV is an input argument and on entry<br>
*>     contains the pivot indices from the factorization A = P*L*U<br>
*>     as computed by DGETRF; row i of the matrix was interchanged<br>
*>     with row IPIV(i).<br>
*><br>
*>     If FACT = 'N', then IPIV is an output argument and on exit<br>
*>     contains the pivot indices from the factorization A = P*L*U<br>
*>     of the original matrix A.<br>
*><br>
*>     If FACT = 'E', then IPIV is an output argument and on exit<br>
*>     contains the pivot indices from the factorization A = P*L*U<br>
*>     of the equilibrated matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>     Specifies the form of equilibration that was done.<br>
*>       = 'N':  No equilibration (always true if FACT = 'N').<br>
*>       = 'R':  Row equilibration, i.e., A has been premultiplied by<br>
*>               diag(R).<br>
*>       = 'C':  Column equilibration, i.e., A has been postmultiplied<br>
*>               by diag(C).<br>
*>       = 'B':  Both row and column equilibration, i.e., A has been<br>
*>               replaced by diag(R) * A * diag(C).<br>
*>     EQUED is an input argument if FACT = 'F'; otherwise, it is an<br>
*>     output argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] R<br>
*> \verbatim<br>
*>          R is DOUBLE PRECISION array, dimension (N)<br>
*>     The row scale factors for A.  If EQUED = 'R' or 'B', A is<br>
*>     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R<br>
*>     is not accessed.  R is an input argument if FACT = 'F';<br>
*>     otherwise, R is an output argument.  If FACT = 'F' and<br>
*>     EQUED = 'R' or 'B', each element of R must be positive.<br>
*>     If R is output, each element of R is a power of the radix.<br>
*>     If R is input, each element of R should be a power of the radix<br>
*>     to ensure a reliable solution and error estimates. Scaling by<br>
*>     powers of the radix does not cause rounding errors unless the<br>
*>     result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (N)<br>
*>     The column scale factors for A.  If EQUED = 'C' or 'B', A is<br>
*>     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C<br>
*>     is not accessed.  C is an input argument if FACT = 'F';<br>
*>     otherwise, C is an output argument.  If FACT = 'F' and<br>
*>     EQUED = 'C' or 'B', each element of C must be positive.<br>
*>     If C is output, each element of C is a power of the radix.<br>
*>     If C is input, each element of C should be a power of the radix<br>
*>     to ensure a reliable solution and error estimates. Scaling by<br>
*>     powers of the radix does not cause rounding errors unless the<br>
*>     result underflows or overflows. Rounding errors during scaling<br>
*>     lead to refining with a matrix that is not equivalent to the<br>
*>     input matrix, producing error estimates that may not be<br>
*>     reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>     On entry, the N-by-NRHS right hand side matrix B.<br>
*>     On exit,<br>
*>     if EQUED = 'N', B is not modified;<br>
*>     if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by<br>
*>        diag(R)*B;<br>
*>     if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is<br>
*>        overwritten by diag(C)*B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>     The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>     If INFO = 0, the N-by-NRHS solution matrix X to the original<br>
*>     system of equations.  Note that A and B are modified on exit<br>
*>     if EQUED .ne. 'N', and the solution to the equilibrated system is<br>
*>     inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or 'B', or<br>
*>     inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>     The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>     Reciprocal scaled condition number.  This is an estimate of the<br>
*>     reciprocal Skeel condition number of the matrix A after<br>
*>     equilibration (if done).  If this is less than the machine<br>
*>     precision (in particular, if it is zero), the matrix is singular<br>
*>     to working precision.  Note that the error may still be small even<br>
*>     if this number is very small and the matrix appears ill-<br>
*>     conditioned.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RPVGRW<br>
*> \verbatim<br>
*>          RPVGRW is DOUBLE PRECISION<br>
*>     Reciprocal pivot growth.  On exit, this contains the reciprocal<br>
*>     pivot growth factor norm(A)/norm(U). The "max absolute element"<br>
*>     norm is used.  If this is much less than 1, then the stability of<br>
*>     the LU factorization of the (equilibrated) matrix A could be poor.<br>
*>     This also means that the solution X, estimated condition numbers,<br>
*>     and error bounds could be unreliable. If factorization fails with<br>
*>     0<INFO<=N, then this contains the reciprocal pivot growth factor<br>
*>     for the leading INFO columns of A.  In DGESVX, this quantity is<br>
*>     returned in WORK(1).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>     Componentwise relative backward error.  This is the<br>
*>     componentwise relative backward error of each solution vector X(j)<br>
*>     (i.e., the smallest relative change in any element of A or B that<br>
*>     makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[in] N_ERR_BNDS<br>
*> \verbatim<br>
*>          N_ERR_BNDS is INTEGER<br>
*>     Number of error bounds to return for each right hand side<br>
*>     and each type (normwise or componentwise).  See ERR_BNDS_NORM and<br>
*>     ERR_BNDS_COMP below.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_NORM<br>
*> \verbatim<br>
*>          ERR_BNDS_NORM is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>              sqrt(n) * dlamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * dlamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated normwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * dlamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*A, where S scales each row by a power of the<br>
*>              radix so all absolute row sums of Z are approximately 1.<br>
*><br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_COMP<br>
*> \verbatim<br>
*>          ERR_BNDS_COMP is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)<br>
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
*>              sqrt(n) * dlamch('Epsilon').<br>
*><br>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,<br>
*>              almost certainly within a factor of 10 of the true error<br>
*>              so long as the next entry is greater than the threshold<br>
*>              sqrt(n) * dlamch('Epsilon'). This error bound should only<br>
*>              be trusted if the previous boolean is true.<br>
*><br>
*>     err = 3  Reciprocal condition number: Estimated componentwise<br>
*>              reciprocal condition number.  Compared with the threshold<br>
*>              sqrt(n) * dlamch('Epsilon') to determine if the error<br>
*>              estimate is "guaranteed". These reciprocal condition<br>
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some<br>
*>              appropriately scaled matrix Z.<br>
*>              Let Z = S*(A*diag(x)), where x is the solution for the<br>
*>              current right-hand side and S scales each row of<br>
*>              A*diag(x) by a power of the radix so all absolute row<br>
*>              sums of Z are approximately 1.<br>
*><br>
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NPARAMS<br>
*> \verbatim<br>
*>          NPARAMS is INTEGER<br>
*>     Specifies the number of parameters set in PARAMS.  If .LE. 0, the<br>
*>     PARAMS array is never referenced and default values are used.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] PARAMS<br>
*> \verbatim<br>
*>          PARAMS is DOUBLE PRECISION array, dimension (NPARAMS)<br>
*>     Specifies algorithm parameters.  If an entry is .LT. 0.0, then<br>
*>     that entry will be filled with default value used for that<br>
*>     parameter.  Only positions up to NPARAMS are accessed; defaults<br>
*>     are used for higher-numbered parameters.<br>
*><br>
*>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative<br>
*>            refinement or not.<br>
*>         Default: 1.0D+0<br>
*>            = 0.0 : No refinement is performed, and no error bounds are<br>
*>                    computed.<br>
*>            = 1.0 : Use the extra-precise refinement algorithm.<br>
*>              (other values are reserved for future use)<br>
*><br>
*>       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual<br>
*>            computations allowed for refinement.<br>
*>         Default: 10<br>
*>         Aggressive: Set to 100 to permit convergence using approximate<br>
*>                     factorizations or factorizations other than LU. If<br>
*>                     the factorization uses a technique other than<br>
*>                     Gaussian elimination, the guarantees in<br>
*>                     err_bnds_norm and err_bnds_comp may no longer be<br>
*>                     trustworthy.<br>
*><br>
*>       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code<br>
*>            will attempt to find a solution with small componentwise<br>
*>            relative error in the double-precision algorithm.  Positive<br>
*>            is true, 0.0 is false.<br>
*>         Default: 1.0 (attempt componentwise convergence)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>       = 0:  Successful exit. The solution to every right-hand side is<br>
*>         guaranteed.<br>
*>       < 0:  If INFO = -i, the i-th argument had an illegal value<br>
*>       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization<br>
*>         has been completed, but the factor U is exactly singular, so<br>
*>         the solution and error bounds could not be computed. RCOND = 0<br>
*>         is returned.<br>
*>       = N+J: The solution corresponding to the Jth right-hand side is<br>
*>         not guaranteed. The solutions corresponding to other right-<br>
*>         hand sides K with K > J may not be guaranteed as well, but<br>
*>         only the first such right-hand side is reported. If a small<br>
*>         componentwise error is not requested (PARAMS(3) = 0.0) then<br>
*>         the Jth right-hand side is the first with a normwise error<br>
*>         bound that is not guaranteed (the smallest J such<br>
*>         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)<br>
*>         the Jth right-hand side is the first with either a normwise or<br>
*>         componentwise error bound that is not guaranteed (the smallest<br>
*>         J such that either ERR_BNDS_NORM(J,1) = 0.0 or<br>
*>         ERR_BNDS_COMP(J,1) = 0.0). See the definition of<br>
*>         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information<br>
*>         about all of the right-hand sides check ERR_BNDS_NORM or<br>
*>         ERR_BNDS_COMP.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgesvxx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE DLA_GERPVGRW);
/**
*> \brief \b DGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGETC2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetc2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetc2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetc2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGETC2( N, A, LDA, IPIV, JPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), JPIV( * )<br>
*       DOUBLE PRECISION   A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGETC2 computes an LU factorization with complete pivoting of the<br>
*> n-by-n matrix A. The factorization has the form A = P * L * U * Q,<br>
*> where P and Q are permutation matrices, L is lower triangular with<br>
*> unit diagonal elements and U is upper triangular.<br>
*><br>
*> This is the Level 2 BLAS algorithm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the n-by-n matrix A to be factored.<br>
*>          On exit, the factors L and U from the factorization<br>
*>          A = P*L*U*Q; the unit diagonal elements of L are not stored.<br>
*>          If U(k, k) appears to be less than SMIN, U(k, k) is given the<br>
*>          value of SMIN, i.e., giving a nonsingular perturbed system.<br>
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
*>          IPIV is INTEGER array, dimension(N).<br>
*>          The pivot indices; for 1 <= i <= N, row i of the<br>
*>          matrix has been interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] JPIV<br>
*> \verbatim<br>
*>          JPIV is INTEGER array, dimension(N).<br>
*>          The pivot indices; for 1 <= j <= N, column j of the<br>
*>          matrix has been interchanged with column JPIV(j).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           = 0: successful exit<br>
*>           > 0: if INFO = k, U(k, k) is likely to produce owerflow if<br>
*>                we try to solve for x in Ax = b. So U is perturbed to<br>
*>                avoid the overflow.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,<br>
*>     Umea University, S-901 87 Umea, Sweden.<br>
*<br>
*  =====================================================================<br>
*/
	public void dgetc2_(INTEGER N,double[] A,INTEGER LDA,int[] IPIV,int[] JPIV,INTEGER INFO);
/**
*> \brief \b DGETF2 computes the LU factorization of a general m-by-n matrix using partial pivoting with row interchanges (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGETF2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetf2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetf2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetf2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGETF2 computes an LU factorization of a general m-by-n matrix A<br>
*> using partial pivoting with row interchanges.<br>
*><br>
*> The factorization has the form<br>
*>    A = P * L * U<br>
*> where P is a permutation matrix, L is lower triangular with unit<br>
*> diagonal elements (lower trapezoidal if m > n), and U is upper<br>
*> triangular (upper trapezoidal if m < n).<br>
*><br>
*> This is the right-looking Level 2 BLAS version of the algorithm.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the m by n matrix to be factored.<br>
*>          On exit, the factors L and U from the factorization<br>
*>          A = P*L*U; the unit diagonal elements of L are not stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (min(M,N))<br>
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -k, the k-th argument had an illegal value<br>
*>          > 0: if INFO = k, U(k,k) is exactly zero. The factorization<br>
*>               has been completed, but the factor U is exactly<br>
*>               singular, and division by zero will occur if it is used<br>
*>               to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgetf2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b DGETRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGETRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGETRF computes an LU factorization of a general M-by-N matrix A<br>
*> using partial pivoting with row interchanges.<br>
*><br>
*> The factorization has the form<br>
*>    A = P * L * U<br>
*> where P is a permutation matrix, L is lower triangular with unit<br>
*> diagonal elements (lower trapezoidal if m > n), and U is upper<br>
*> triangular (upper trapezoidal if m < n).<br>
*><br>
*> This is the right-looking Level 3 BLAS version of the algorithm.<br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix to be factored.<br>
*>          On exit, the factors L and U from the factorization<br>
*>          A = P*L*U; the unit diagonal elements of L are not stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (min(M,N))<br>
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization<br>
*>                has been completed, but the factor U is exactly<br>
*>                singular, and division by zero will occur if it is used<br>
*>                to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgetrf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b DGETRF2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGETRF2 computes an LU factorization of a general M-by-N matrix A<br>
*> using partial pivoting with row interchanges.<br>
*><br>
*> The factorization has the form<br>
*>    A = P * L * U<br>
*> where P is a permutation matrix, L is lower triangular with unit<br>
*> diagonal elements (lower trapezoidal if m > n), and U is upper<br>
*> triangular (upper trapezoidal if m < n).<br>
*><br>
*> This is the recursive version of the algorithm. It divides<br>
*> the matrix into four submatrices:<br>
*>            <br>
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2<br>
*>    A = [ -----|----- ]  with n1 = min(m,n)/2<br>
*>        [  A21 | A22  ]       n2 = n-n1<br>
*>            <br>
*>                                       [ A11 ]<br>
*> The subroutine calls itself to factor [ --- ],<br>
*>                                       [ A12 ]<br>
*>                 [ A12 ]<br>
*> do the swaps on [ --- ], solve A12, update A22,<br>
*>                 [ A22 ]<br>
*><br>
*> then calls itself to factor A22 and do the swaps on A21.<br>
*><br>
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix to be factored.<br>
*>          On exit, the factors L and U from the factorization<br>
*>          A = P*L*U; the unit diagonal elements of L are not stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (min(M,N))<br>
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization<br>
*>                has been completed, but the factor U is exactly<br>
*>                singular, and division by zero will occur if it is used<br>
*>                to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgetrf2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b DGETRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGETRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGETRI computes the inverse of a matrix using the LU factorization<br>
*> computed by DGETRF.<br>
*><br>
*> This method inverts U and then computes inv(A) by solving the system<br>
*> inv(A)*L = inv(U) for inv(A).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the factors L and U from the factorization<br>
*>          A = P*L*U as computed by DGETRF.<br>
*>          On exit, if INFO = 0, the inverse of the original matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices from DGETRF; for 1<=i<=N, row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.  LWORK >= max(1,N).<br>
*>          For optimal performance LWORK >= N*NB, where NB is<br>
*>          the optimal blocksize returned by ILAENV.<br>
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
*>          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is<br>
*>                singular and its inverse could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgetri_(INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGETRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGETRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGETRS solves a system of linear equations<br>
*>    A * X = B  or  A**T * X = B<br>
*> with a general N-by-N matrix A using the LU factorization computed<br>
*> by DGETRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T* X = B  (Transpose)<br>
*>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)<br>
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
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          The factors L and U from the factorization A = P*L*U<br>
*>          as computed by DGETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices from DGETRF; for 1<=i<=N, row i of the<br>
*>          matrix was interchanged with row IPIV(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the right hand side matrix B.<br>
*>          On exit, the solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgetrs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b DGGBAK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGBAK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbak.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbak.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbak.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,<br>
*                          LDV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOB, SIDE<br>
*       INTEGER            IHI, ILO, INFO, LDV, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), V( LDV, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGBAK forms the right or left eigenvectors of a real generalized<br>
*> eigenvalue problem A*x = lambda*B*x, by backward transformation on<br>
*> the computed eigenvectors of the balanced pair of matrices output by<br>
*> DGGBAL.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          Specifies the type of backward transformation required:<br>
*>          = 'N':  do nothing, return immediately;<br>
*>          = 'P':  do backward transformation for permutation only;<br>
*>          = 'S':  do backward transformation for scaling only;<br>
*>          = 'B':  do backward transformations for both permutation and<br>
*>                  scaling.<br>
*>          JOB must be the same as the argument JOB supplied to DGGBAL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'R':  V contains right eigenvectors;<br>
*>          = 'L':  V contains left eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of rows of the matrix V.  N >= 0.<br>
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
*>          The integers ILO and IHI determined by DGGBAL.<br>
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LSCALE<br>
*> \verbatim<br>
*>          LSCALE is DOUBLE PRECISION array, dimension (N)<br>
*>          Details of the permutations and/or scaling factors applied<br>
*>          to the left side of A and B, as returned by DGGBAL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RSCALE<br>
*> \verbatim<br>
*>          RSCALE is DOUBLE PRECISION array, dimension (N)<br>
*>          Details of the permutations and/or scaling factors applied<br>
*>          to the right side of A and B, as returned by DGGBAL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of columns of the matrix V.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension (LDV,M)<br>
*>          On entry, the matrix of right or left eigenvectors to be<br>
*>          transformed, as returned by DTGEVC.<br>
*>          On exit, V is overwritten by the transformed eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the matrix V. LDV >= max(1,N).<br>
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
*> \ingroup doubleGBcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  See R.C. Ward, Balancing the generalized eigenvalue problem,<br>
*>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dggbak_(CHARACTER JOB,CHARACTER SIDE,INTEGER N,INTEGER ILO,INTEGER IHI,double[] LSCALE,double[] RSCALE,INTEGER M,double[] V,INTEGER LDV,INTEGER INFO);
/**
*> \brief \b DGGBAL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGBAL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbal.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbal.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbal.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,<br>
*                          RSCALE, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOB<br>
*       INTEGER            IHI, ILO, INFO, LDA, LDB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ),<br>
*      $                   RSCALE( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGBAL balances a pair of general real matrices (A,B).  This<br>
*> involves, first, permuting A and B by similarity transformations to<br>
*> isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N<br>
*> elements on the diagonal; and second, applying a diagonal similarity<br>
*> transformation to rows and columns ILO to IHI to make the rows<br>
*> and columns as close in norm as possible. Both steps are optional.<br>
*><br>
*> Balancing may reduce the 1-norm of the matrices, and improve the<br>
*> accuracy of the computed eigenvalues and/or eigenvectors in the<br>
*> generalized eigenvalue problem A*x = lambda*B*x.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          Specifies the operations to be performed on A and B:<br>
*>          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0<br>
*>                  and RSCALE(I) = 1.0 for i = 1,...,N.<br>
*>          = 'P':  permute only;<br>
*>          = 'S':  scale only;<br>
*>          = 'B':  both permute and scale.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the input matrix A.<br>
*>          On exit,  A is overwritten by the balanced matrix.<br>
*>          If JOB = 'N', A is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,N)<br>
*>          On entry, the input matrix B.<br>
*>          On exit,  B is overwritten by the balanced matrix.<br>
*>          If JOB = 'N', B is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ILO<br>
*> \verbatim<br>
*>          ILO is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[out] IHI<br>
*> \verbatim<br>
*>          IHI is INTEGER<br>
*>          ILO and IHI are set to integers such that on exit<br>
*>          A(i,j) = 0 and B(i,j) = 0 if i > j and<br>
*>          j = 1,...,ILO-1 or i = IHI+1,...,N.<br>
*>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] LSCALE<br>
*> \verbatim<br>
*>          LSCALE is DOUBLE PRECISION array, dimension (N)<br>
*>          Details of the permutations and scaling factors applied<br>
*>          to the left side of A and B.  If P(j) is the index of the<br>
*>          row interchanged with row j, and D(j)<br>
*>          is the scaling factor applied to row j, then<br>
*>            LSCALE(j) = P(j)    for J = 1,...,ILO-1<br>
*>                      = D(j)    for J = ILO,...,IHI<br>
*>                      = P(j)    for J = IHI+1,...,N.<br>
*>          The order in which the interchanges are made is N to IHI+1,<br>
*>          then 1 to ILO-1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RSCALE<br>
*> \verbatim<br>
*>          RSCALE is DOUBLE PRECISION array, dimension (N)<br>
*>          Details of the permutations and scaling factors applied<br>
*>          to the right side of A and B.  If P(j) is the index of the<br>
*>          column interchanged with column j, and D(j)<br>
*>          is the scaling factor applied to column j, then<br>
*>            LSCALE(j) = P(j)    for J = 1,...,ILO-1<br>
*>                      = D(j)    for J = ILO,...,IHI<br>
*>                      = P(j)    for J = IHI+1,...,N.<br>
*>          The order in which the interchanges are made is N to IHI+1,<br>
*>          then 1 to ILO-1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (lwork)<br>
*>          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and<br>
*>          at least 1 when JOB = 'N' or 'P'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
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
*> \ingroup doubleGBcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  See R.C. WARD, Balancing the generalized eigenvalue problem,<br>
*>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dggbal_(CHARACTER JOB,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER ILO,INTEGER IHI,double[] LSCALE,double[] RSCALE,double[] WORK,INTEGER INFO);
/**
*> \brief <b> DGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGES + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgges.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgges.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgges.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB,<br>
*                         SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR,<br>
*                         LDVSR, WORK, LWORK, BWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBVSL, JOBVSR, SORT<br>
*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            BWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),<br>
*      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ),<br>
*      $                   VSR( LDVSR, * ), WORK( * )<br>
*       ..<br>
*       .. Function Arguments ..<br>
*       LOGICAL            SELCTG<br>
*       EXTERNAL           SELCTG<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGES computes for a pair of N-by-N real nonsymmetric matrices (A,B),<br>
*> the generalized eigenvalues, the generalized real Schur form (S,T),<br>
*> optionally, the left and/or right matrices of Schur vectors (VSL and<br>
*> VSR). This gives the generalized Schur factorization<br>
*><br>
*>          (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )<br>
*><br>
*> Optionally, it also orders the eigenvalues so that a selected cluster<br>
*> of eigenvalues appears in the leading diagonal blocks of the upper<br>
*> quasi-triangular matrix S and the upper triangular matrix T.The<br>
*> leading columns of VSL and VSR then form an orthonormal basis for the<br>
*> corresponding left and right eigenspaces (deflating subspaces).<br>
*><br>
*> (If only the generalized eigenvalues are needed, use the driver<br>
*> DGGEV instead, which is faster.)<br>
*><br>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w<br>
*> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is<br>
*> usually represented as the pair (alpha,beta), as there is a<br>
*> reasonable interpretation for beta=0 or both being zero.<br>
*><br>
*> A pair of matrices (S,T) is in generalized real Schur form if T is<br>
*> upper triangular with non-negative diagonal and S is block upper<br>
*> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond<br>
*> to real generalized eigenvalues, while 2-by-2 blocks of S will be<br>
*> "standardized" by making the corresponding elements of T have the<br>
*> form:<br>
*>         [  a  0  ]<br>
*>         [  0  b  ]<br>
*><br>
*> and the pair of corresponding 2-by-2 blocks in S and T will have a<br>
*> complex conjugate pair of generalized eigenvalues.<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBVSL<br>
*> \verbatim<br>
*>          JOBVSL is CHARACTER*1<br>
*>          = 'N':  do not compute the left Schur vectors;<br>
*>          = 'V':  compute the left Schur vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVSR<br>
*> \verbatim<br>
*>          JOBVSR is CHARACTER*1<br>
*>          = 'N':  do not compute the right Schur vectors;<br>
*>          = 'V':  compute the right Schur vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SORT<br>
*> \verbatim<br>
*>          SORT is CHARACTER*1<br>
*>          Specifies whether or not to order the eigenvalues on the<br>
*>          diagonal of the generalized Schur form.<br>
*>          = 'N':  Eigenvalues are not ordered;<br>
*>          = 'S':  Eigenvalues are ordered (see SELCTG);<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELCTG<br>
*> \verbatim<br>
*>          SELCTG is a LOGICAL FUNCTION of three DOUBLE PRECISION arguments<br>
*>          SELCTG must be declared EXTERNAL in the calling subroutine.<br>
*>          If SORT = 'N', SELCTG is not referenced.<br>
*>          If SORT = 'S', SELCTG is used to select eigenvalues to sort<br>
*>          to the top left of the Schur form.<br>
*>          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if<br>
*>          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either<br>
*>          one of a complex conjugate pair of eigenvalues is selected,<br>
*>          then both complex eigenvalues are selected.<br>
*><br>
*>          Note that in the ill-conditioned case, a selected complex<br>
*>          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j),<br>
*>          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2<br>
*>          in this case.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A, B, VSL, and VSR.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the first of the pair of matrices.<br>
*>          On exit, A has been overwritten by its generalized Schur<br>
*>          form S.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB, N)<br>
*>          On entry, the second of the pair of matrices.<br>
*>          On exit, B has been overwritten by its generalized Schur<br>
*>          form T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SDIM<br>
*> \verbatim<br>
*>          SDIM is INTEGER<br>
*>          If SORT = 'N', SDIM = 0.<br>
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)<br>
*>          for which SELCTG is true.  (Complex conjugate pairs for which<br>
*>          SELCTG is true for either eigenvalue count as 2.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will<br>
*>          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,<br>
*>          and  BETA(j),j=1,...,N are the diagonals of the complex Schur<br>
*>          form (S,T) that would result if the 2-by-2 diagonal blocks of<br>
*>          the real Schur form of (A,B) were further reduced to<br>
*>          triangular form using 2-by-2 complex unitary transformations.<br>
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if<br>
*>          positive, then the j-th and (j+1)-st eigenvalues are a<br>
*>          complex conjugate pair, with ALPHAI(j+1) negative.<br>
*><br>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)<br>
*>          may easily over- or underflow, and BETA(j) may even be zero.<br>
*>          Thus, the user should avoid naively computing the ratio.<br>
*>          However, ALPHAR and ALPHAI will be always less than and<br>
*>          usually comparable with norm(A) in magnitude, and BETA always<br>
*>          less than and usually comparable with norm(B).<br>
*> \endverbatim<br>
*><br>
*> \param[out] VSL<br>
*> \verbatim<br>
*>          VSL is DOUBLE PRECISION array, dimension (LDVSL,N)<br>
*>          If JOBVSL = 'V', VSL will contain the left Schur vectors.<br>
*>          Not referenced if JOBVSL = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVSL<br>
*> \verbatim<br>
*>          LDVSL is INTEGER<br>
*>          The leading dimension of the matrix VSL. LDVSL >=1, and<br>
*>          if JOBVSL = 'V', LDVSL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VSR<br>
*> \verbatim<br>
*>          VSR is DOUBLE PRECISION array, dimension (LDVSR,N)<br>
*>          If JOBVSR = 'V', VSR will contain the right Schur vectors.<br>
*>          Not referenced if JOBVSR = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVSR<br>
*> \verbatim<br>
*>          LDVSR is INTEGER<br>
*>          The leading dimension of the matrix VSR. LDVSR >= 1, and<br>
*>          if JOBVSR = 'V', LDVSR >= N.<br>
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
*>          If N = 0, LWORK >= 1, else LWORK >= 8*N+16.<br>
*>          For good performance , LWORK must generally be larger.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BWORK<br>
*> \verbatim<br>
*>          BWORK is LOGICAL array, dimension (N)<br>
*>          Not referenced if SORT = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          = 1,...,N:<br>
*>                The QZ iteration failed.  (A,B) are not in Schur<br>
*>                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should<br>
*>                be correct for j=INFO+1,...,N.<br>
*>          > N:  =N+1: other than QZ iteration failed in DHGEQZ.<br>
*>                =N+2: after reordering, roundoff changed values of<br>
*>                      some complex eigenvalues so that leading<br>
*>                      eigenvalues in the Generalized Schur form no<br>
*>                      longer satisfy SELCTG=.TRUE.  This could also<br>
*>                      be caused due to scaling.<br>
*>                =N+3: reordering failed in DTGSEN.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dgges_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER SDIM,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VSL,INTEGER LDVSL,double[] VSR,INTEGER LDVSR,double[] WORK,INTEGER LWORK,boolean[] BWORK,INTEGER INFO);
/**
*> \brief <b> DGGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices (blocked algorithm)</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download DGGES3 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgges3.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgges3.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgges3.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB,<br>
*                          SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR,<br>
*                          LDVSR, WORK, LWORK, BWORK, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBVSL, JOBVSR, SORT<br>
*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            BWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),<br>
*      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ),<br>
*      $                   VSR( LDVSR, * ), WORK( * )<br>
*       ..<br>
*       .. Function Arguments ..<br>
*       LOGICAL            SELCTG<br>
*       EXTERNAL           SELCTG<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGES3 computes for a pair of N-by-N real nonsymmetric matrices (A,B),<br>
*> the generalized eigenvalues, the generalized real Schur form (S,T),<br>
*> optionally, the left and/or right matrices of Schur vectors (VSL and<br>
*> VSR). This gives the generalized Schur factorization<br>
*><br>
*>          (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )<br>
*><br>
*> Optionally, it also orders the eigenvalues so that a selected cluster<br>
*> of eigenvalues appears in the leading diagonal blocks of the upper<br>
*> quasi-triangular matrix S and the upper triangular matrix T.The<br>
*> leading columns of VSL and VSR then form an orthonormal basis for the<br>
*> corresponding left and right eigenspaces (deflating subspaces).<br>
*><br>
*> (If only the generalized eigenvalues are needed, use the driver<br>
*> DGGEV instead, which is faster.)<br>
*><br>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w<br>
*> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is<br>
*> usually represented as the pair (alpha,beta), as there is a<br>
*> reasonable interpretation for beta=0 or both being zero.<br>
*><br>
*> A pair of matrices (S,T) is in generalized real Schur form if T is<br>
*> upper triangular with non-negative diagonal and S is block upper<br>
*> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond<br>
*> to real generalized eigenvalues, while 2-by-2 blocks of S will be<br>
*> "standardized" by making the corresponding elements of T have the<br>
*> form:<br>
*>         [  a  0  ]<br>
*>         [  0  b  ]<br>
*><br>
*> and the pair of corresponding 2-by-2 blocks in S and T will have a<br>
*> complex conjugate pair of generalized eigenvalues.<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBVSL<br>
*> \verbatim<br>
*>          JOBVSL is CHARACTER*1<br>
*>          = 'N':  do not compute the left Schur vectors;<br>
*>          = 'V':  compute the left Schur vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVSR<br>
*> \verbatim<br>
*>          JOBVSR is CHARACTER*1<br>
*>          = 'N':  do not compute the right Schur vectors;<br>
*>          = 'V':  compute the right Schur vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SORT<br>
*> \verbatim<br>
*>          SORT is CHARACTER*1<br>
*>          Specifies whether or not to order the eigenvalues on the<br>
*>          diagonal of the generalized Schur form.<br>
*>          = 'N':  Eigenvalues are not ordered;<br>
*>          = 'S':  Eigenvalues are ordered (see SELCTG);<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELCTG<br>
*> \verbatim<br>
*>          SELCTG is a LOGICAL FUNCTION of three DOUBLE PRECISION arguments<br>
*>          SELCTG must be declared EXTERNAL in the calling subroutine.<br>
*>          If SORT = 'N', SELCTG is not referenced.<br>
*>          If SORT = 'S', SELCTG is used to select eigenvalues to sort<br>
*>          to the top left of the Schur form.<br>
*>          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if<br>
*>          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either<br>
*>          one of a complex conjugate pair of eigenvalues is selected,<br>
*>          then both complex eigenvalues are selected.<br>
*><br>
*>          Note that in the ill-conditioned case, a selected complex<br>
*>          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j),<br>
*>          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2<br>
*>          in this case.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A, B, VSL, and VSR.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the first of the pair of matrices.<br>
*>          On exit, A has been overwritten by its generalized Schur<br>
*>          form S.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB, N)<br>
*>          On entry, the second of the pair of matrices.<br>
*>          On exit, B has been overwritten by its generalized Schur<br>
*>          form T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SDIM<br>
*> \verbatim<br>
*>          SDIM is INTEGER<br>
*>          If SORT = 'N', SDIM = 0.<br>
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)<br>
*>          for which SELCTG is true.  (Complex conjugate pairs for which<br>
*>          SELCTG is true for either eigenvalue count as 2.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will<br>
*>          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,<br>
*>          and  BETA(j),j=1,...,N are the diagonals of the complex Schur<br>
*>          form (S,T) that would result if the 2-by-2 diagonal blocks of<br>
*>          the real Schur form of (A,B) were further reduced to<br>
*>          triangular form using 2-by-2 complex unitary transformations.<br>
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if<br>
*>          positive, then the j-th and (j+1)-st eigenvalues are a<br>
*>          complex conjugate pair, with ALPHAI(j+1) negative.<br>
*><br>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)<br>
*>          may easily over- or underflow, and BETA(j) may even be zero.<br>
*>          Thus, the user should avoid naively computing the ratio.<br>
*>          However, ALPHAR and ALPHAI will be always less than and<br>
*>          usually comparable with norm(A) in magnitude, and BETA always<br>
*>          less than and usually comparable with norm(B).<br>
*> \endverbatim<br>
*><br>
*> \param[out] VSL<br>
*> \verbatim<br>
*>          VSL is DOUBLE PRECISION array, dimension (LDVSL,N)<br>
*>          If JOBVSL = 'V', VSL will contain the left Schur vectors.<br>
*>          Not referenced if JOBVSL = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVSL<br>
*> \verbatim<br>
*>          LDVSL is INTEGER<br>
*>          The leading dimension of the matrix VSL. LDVSL >=1, and<br>
*>          if JOBVSL = 'V', LDVSL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VSR<br>
*> \verbatim<br>
*>          VSR is DOUBLE PRECISION array, dimension (LDVSR,N)<br>
*>          If JOBVSR = 'V', VSR will contain the right Schur vectors.<br>
*>          Not referenced if JOBVSR = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVSR<br>
*> \verbatim<br>
*>          LDVSR is INTEGER<br>
*>          The leading dimension of the matrix VSR. LDVSR >= 1, and<br>
*>          if JOBVSR = 'V', LDVSR >= N.<br>
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
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BWORK<br>
*> \verbatim<br>
*>          BWORK is LOGICAL array, dimension (N)<br>
*>          Not referenced if SORT = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          = 1,...,N:<br>
*>                The QZ iteration failed.  (A,B) are not in Schur<br>
*>                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should<br>
*>                be correct for j=INFO+1,...,N.<br>
*>          > N:  =N+1: other than QZ iteration failed in DHGEQZ.<br>
*>                =N+2: after reordering, roundoff changed values of<br>
*>                      some complex eigenvalues so that leading<br>
*>                      eigenvalues in the Generalized Schur form no<br>
*>                      longer satisfy SELCTG=.TRUE.  This could also<br>
*>                      be caused due to scaling.<br>
*>                =N+3: reordering failed in DTGSEN.<br>
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
*> \date January 2015<br>
*<br>
*> \ingroup doubleGEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dgges3_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER SDIM,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VSL,INTEGER LDVSL,double[] VSR,INTEGER LDVSR,double[] WORK,INTEGER LWORK,boolean[] BWORK,INTEGER INFO);
/**
*> \brief <b> DGGESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGESX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggesx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggesx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggesx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA,<br>
*                          B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL,<br>
*                          VSR, LDVSR, RCONDE, RCONDV, WORK, LWORK, IWORK,<br>
*                          LIWORK, BWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBVSL, JOBVSR, SENSE, SORT<br>
*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LIWORK, LWORK, N,<br>
*      $                   SDIM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            BWORK( * )<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),<br>
*      $                   B( LDB, * ), BETA( * ), RCONDE( 2 ),<br>
*      $                   RCONDV( 2 ), VSL( LDVSL, * ), VSR( LDVSR, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*       .. Function Arguments ..<br>
*       LOGICAL            SELCTG<br>
*       EXTERNAL           SELCTG<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGESX computes for a pair of N-by-N real nonsymmetric matrices<br>
*> (A,B), the generalized eigenvalues, the real Schur form (S,T), and,<br>
*> optionally, the left and/or right matrices of Schur vectors (VSL and<br>
*> VSR).  This gives the generalized Schur factorization<br>
*><br>
*>      (A,B) = ( (VSL) S (VSR)**T, (VSL) T (VSR)**T )<br>
*><br>
*> Optionally, it also orders the eigenvalues so that a selected cluster<br>
*> of eigenvalues appears in the leading diagonal blocks of the upper<br>
*> quasi-triangular matrix S and the upper triangular matrix T; computes<br>
*> a reciprocal condition number for the average of the selected<br>
*> eigenvalues (RCONDE); and computes a reciprocal condition number for<br>
*> the right and left deflating subspaces corresponding to the selected<br>
*> eigenvalues (RCONDV). The leading columns of VSL and VSR then form<br>
*> an orthonormal basis for the corresponding left and right eigenspaces<br>
*> (deflating subspaces).<br>
*><br>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w<br>
*> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is<br>
*> usually represented as the pair (alpha,beta), as there is a<br>
*> reasonable interpretation for beta=0 or for both being zero.<br>
*><br>
*> A pair of matrices (S,T) is in generalized real Schur form if T is<br>
*> upper triangular with non-negative diagonal and S is block upper<br>
*> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond<br>
*> to real generalized eigenvalues, while 2-by-2 blocks of S will be<br>
*> "standardized" by making the corresponding elements of T have the<br>
*> form:<br>
*>         [  a  0  ]<br>
*>         [  0  b  ]<br>
*><br>
*> and the pair of corresponding 2-by-2 blocks in S and T will have a<br>
*> complex conjugate pair of generalized eigenvalues.<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBVSL<br>
*> \verbatim<br>
*>          JOBVSL is CHARACTER*1<br>
*>          = 'N':  do not compute the left Schur vectors;<br>
*>          = 'V':  compute the left Schur vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVSR<br>
*> \verbatim<br>
*>          JOBVSR is CHARACTER*1<br>
*>          = 'N':  do not compute the right Schur vectors;<br>
*>          = 'V':  compute the right Schur vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SORT<br>
*> \verbatim<br>
*>          SORT is CHARACTER*1<br>
*>          Specifies whether or not to order the eigenvalues on the<br>
*>          diagonal of the generalized Schur form.<br>
*>          = 'N':  Eigenvalues are not ordered;<br>
*>          = 'S':  Eigenvalues are ordered (see SELCTG).<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELCTG<br>
*> \verbatim<br>
*>          SELCTG is procedure) LOGICAL FUNCTION of three DOUBLE PRECISION arguments<br>
*>          SELCTG must be declared EXTERNAL in the calling subroutine.<br>
*>          If SORT = 'N', SELCTG is not referenced.<br>
*>          If SORT = 'S', SELCTG is used to select eigenvalues to sort<br>
*>          to the top left of the Schur form.<br>
*>          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if<br>
*>          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either<br>
*>          one of a complex conjugate pair of eigenvalues is selected,<br>
*>          then both complex eigenvalues are selected.<br>
*>          Note that a selected complex eigenvalue may no longer satisfy<br>
*>          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) = .TRUE. after ordering,<br>
*>          since ordering may change the value of complex eigenvalues<br>
*>          (especially if the eigenvalue is ill-conditioned), in this<br>
*>          case INFO is set to N+3.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SENSE<br>
*> \verbatim<br>
*>          SENSE is CHARACTER*1<br>
*>          Determines which reciprocal condition numbers are computed.<br>
*>          = 'N' : None are computed;<br>
*>          = 'E' : Computed for average of selected eigenvalues only;<br>
*>          = 'V' : Computed for selected deflating subspaces only;<br>
*>          = 'B' : Computed for both.<br>
*>          If SENSE = 'E', 'V', or 'B', SORT must equal 'S'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A, B, VSL, and VSR.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the first of the pair of matrices.<br>
*>          On exit, A has been overwritten by its generalized Schur<br>
*>          form S.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB, N)<br>
*>          On entry, the second of the pair of matrices.<br>
*>          On exit, B has been overwritten by its generalized Schur<br>
*>          form T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SDIM<br>
*> \verbatim<br>
*>          SDIM is INTEGER<br>
*>          If SORT = 'N', SDIM = 0.<br>
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)<br>
*>          for which SELCTG is true.  (Complex conjugate pairs for which<br>
*>          SELCTG is true for either eigenvalue count as 2.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will<br>
*>          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i<br>
*>          and BETA(j),j=1,...,N  are the diagonals of the complex Schur<br>
*>          form (S,T) that would result if the 2-by-2 diagonal blocks of<br>
*>          the real Schur form of (A,B) were further reduced to<br>
*>          triangular form using 2-by-2 complex unitary transformations.<br>
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if<br>
*>          positive, then the j-th and (j+1)-st eigenvalues are a<br>
*>          complex conjugate pair, with ALPHAI(j+1) negative.<br>
*><br>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)<br>
*>          may easily over- or underflow, and BETA(j) may even be zero.<br>
*>          Thus, the user should avoid naively computing the ratio.<br>
*>          However, ALPHAR and ALPHAI will be always less than and<br>
*>          usually comparable with norm(A) in magnitude, and BETA always<br>
*>          less than and usually comparable with norm(B).<br>
*> \endverbatim<br>
*><br>
*> \param[out] VSL<br>
*> \verbatim<br>
*>          VSL is DOUBLE PRECISION array, dimension (LDVSL,N)<br>
*>          If JOBVSL = 'V', VSL will contain the left Schur vectors.<br>
*>          Not referenced if JOBVSL = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVSL<br>
*> \verbatim<br>
*>          LDVSL is INTEGER<br>
*>          The leading dimension of the matrix VSL. LDVSL >=1, and<br>
*>          if JOBVSL = 'V', LDVSL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VSR<br>
*> \verbatim<br>
*>          VSR is DOUBLE PRECISION array, dimension (LDVSR,N)<br>
*>          If JOBVSR = 'V', VSR will contain the right Schur vectors.<br>
*>          Not referenced if JOBVSR = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVSR<br>
*> \verbatim<br>
*>          LDVSR is INTEGER<br>
*>          The leading dimension of the matrix VSR. LDVSR >= 1, and<br>
*>          if JOBVSR = 'V', LDVSR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCONDE<br>
*> \verbatim<br>
*>          RCONDE is DOUBLE PRECISION array, dimension ( 2 )<br>
*>          If SENSE = 'E' or 'B', RCONDE(1) and RCONDE(2) contain the<br>
*>          reciprocal condition numbers for the average of the selected<br>
*>          eigenvalues.<br>
*>          Not referenced if SENSE = 'N' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCONDV<br>
*> \verbatim<br>
*>          RCONDV is DOUBLE PRECISION array, dimension ( 2 )<br>
*>          If SENSE = 'V' or 'B', RCONDV(1) and RCONDV(2) contain the<br>
*>          reciprocal condition numbers for the selected deflating<br>
*>          subspaces.<br>
*>          Not referenced if SENSE = 'N' or 'E'.<br>
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
*>          If N = 0, LWORK >= 1, else if SENSE = 'E', 'V', or 'B',<br>
*>          LWORK >= max( 8*N, 6*N+16, 2*SDIM*(N-SDIM) ), else<br>
*>          LWORK >= max( 8*N, 6*N+16 ).<br>
*>          Note that 2*SDIM*(N-SDIM) <= N*N/2.<br>
*>          Note also that an error is only returned if<br>
*>          LWORK < max( 8*N, 6*N+16), but if SENSE = 'E' or 'V' or 'B'<br>
*>          this may not be large enough.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the bound on the optimal size of the WORK<br>
*>          array and the minimum size of the IWORK array, returns these<br>
*>          values as the first entries of the WORK and IWORK arrays, and<br>
*>          no error message related to LWORK or LIWORK is issued by<br>
*>          XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))<br>
*>          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK.<br>
*>          If SENSE = 'N' or N = 0, LIWORK >= 1, otherwise<br>
*>          LIWORK >= N+6.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the bound on the optimal size of the<br>
*>          WORK array and the minimum size of the IWORK array, returns<br>
*>          these values as the first entries of the WORK and IWORK<br>
*>          arrays, and no error message related to LWORK or LIWORK is<br>
*>          issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BWORK<br>
*> \verbatim<br>
*>          BWORK is LOGICAL array, dimension (N)<br>
*>          Not referenced if SORT = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          = 1,...,N:<br>
*>                The QZ iteration failed.  (A,B) are not in Schur<br>
*>                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should<br>
*>                be correct for j=INFO+1,...,N.<br>
*>          > N:  =N+1: other than QZ iteration failed in DHGEQZ<br>
*>                =N+2: after reordering, roundoff changed values of<br>
*>                      some complex eigenvalues so that leading<br>
*>                      eigenvalues in the Generalized Schur form no<br>
*>                      longer satisfy SELCTG=.TRUE.  This could also<br>
*>                      be caused due to scaling.<br>
*>                =N+3: reordering failed in DTGSEN.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEeigen<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  An approximate (asymptotic) bound on the average absolute error of<br>
*>  the selected eigenvalues is<br>
*><br>
*>       EPS * norm((A, B)) / RCONDE( 1 ).<br>
*><br>
*>  An approximate (asymptotic) bound on the maximum angular error in<br>
*>  the computed deflating subspaces is<br>
*><br>
*>       EPS * norm((A, B)) / RCONDV( 2 ).<br>
*><br>
*>  See LAPACK User's Guide, section 4.11 for more information.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dggesx_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,CHARACTER SENSE,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER SDIM,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VSL,INTEGER LDVSL,double[] VSR,INTEGER LDVSR,double[] RCONDE,double[] RCONDV,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,boolean[] BWORK,INTEGER INFO);
/**
*> \brief <b> DGGEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,<br>
*                         BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBVL, JOBVR<br>
*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),<br>
*      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),<br>
*      $                   VR( LDVR, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)<br>
*> the generalized eigenvalues, and optionally, the left and/or right<br>
*> generalized eigenvectors.<br>
*><br>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar<br>
*> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is<br>
*> singular. It is usually represented as the pair (alpha,beta), as<br>
*> there is a reasonable interpretation for beta=0, and even for both<br>
*> being zero.<br>
*><br>
*> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)<br>
*> of (A,B) satisfies<br>
*><br>
*>                  A * v(j) = lambda(j) * B * v(j).<br>
*><br>
*> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)<br>
*> of (A,B) satisfies<br>
*><br>
*>                  u(j)**H * A  = lambda(j) * u(j)**H * B .<br>
*><br>
*> where u(j)**H is the conjugate-transpose of u(j).<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBVL<br>
*> \verbatim<br>
*>          JOBVL is CHARACTER*1<br>
*>          = 'N':  do not compute the left generalized eigenvectors;<br>
*>          = 'V':  compute the left generalized eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVR<br>
*> \verbatim<br>
*>          JOBVR is CHARACTER*1<br>
*>          = 'N':  do not compute the right generalized eigenvectors;<br>
*>          = 'V':  compute the right generalized eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A, B, VL, and VR.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the matrix A in the pair (A,B).<br>
*>          On exit, A has been overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB, N)<br>
*>          On entry, the matrix B in the pair (A,B).<br>
*>          On exit, B has been overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will<br>
*>          be the generalized eigenvalues.  If ALPHAI(j) is zero, then<br>
*>          the j-th eigenvalue is real; if positive, then the j-th and<br>
*>          (j+1)-st eigenvalues are a complex conjugate pair, with<br>
*>          ALPHAI(j+1) negative.<br>
*><br>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)<br>
*>          may easily over- or underflow, and BETA(j) may even be zero.<br>
*>          Thus, the user should avoid naively computing the ratio<br>
*>          alpha/beta.  However, ALPHAR and ALPHAI will be always less<br>
*>          than and usually comparable with norm(A) in magnitude, and<br>
*>          BETA always less than and usually comparable with norm(B).<br>
*> \endverbatim<br>
*><br>
*> \param[out] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION array, dimension (LDVL,N)<br>
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one<br>
*>          after another in the columns of VL, in the same order as<br>
*>          their eigenvalues. If the j-th eigenvalue is real, then<br>
*>          u(j) = VL(:,j), the j-th column of VL. If the j-th and<br>
*>          (j+1)-th eigenvalues form a complex conjugate pair, then<br>
*>          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).<br>
*>          Each eigenvector is scaled so the largest component has<br>
*>          abs(real part)+abs(imag. part)=1.<br>
*>          Not referenced if JOBVL = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the matrix VL. LDVL >= 1, and<br>
*>          if JOBVL = 'V', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VR<br>
*> \verbatim<br>
*>          VR is DOUBLE PRECISION array, dimension (LDVR,N)<br>
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one<br>
*>          after another in the columns of VR, in the same order as<br>
*>          their eigenvalues. If the j-th eigenvalue is real, then<br>
*>          v(j) = VR(:,j), the j-th column of VR. If the j-th and<br>
*>          (j+1)-th eigenvalues form a complex conjugate pair, then<br>
*>          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).<br>
*>          Each eigenvector is scaled so the largest component has<br>
*>          abs(real part)+abs(imag. part)=1.<br>
*>          Not referenced if JOBVR = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the matrix VR. LDVR >= 1, and<br>
*>          if JOBVR = 'V', LDVR >= N.<br>
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
*>          The dimension of the array WORK.  LWORK >= max(1,8*N).<br>
*>          For good performance, LWORK must generally be larger.<br>
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
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          = 1,...,N:<br>
*>                The QZ iteration failed.  No eigenvectors have been<br>
*>                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)<br>
*>                should be correct for j=INFO+1,...,N.<br>
*>          > N:  =N+1: other than QZ iteration failed in DHGEQZ.<br>
*>                =N+2: error return from DTGEVC.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dggev_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices (blocked algorithm)</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download DGGEV3 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggev3.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggev3.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggev3.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,<br>
*      $                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK,<br>
*      $                   INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBVL, JOBVR<br>
*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),<br>
*      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),<br>
*      $                   VR( LDVR, * ), WORK( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGEV3 computes for a pair of N-by-N real nonsymmetric matrices (A,B)<br>
*> the generalized eigenvalues, and optionally, the left and/or right<br>
*> generalized eigenvectors.<br>
*><br>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar<br>
*> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is<br>
*> singular. It is usually represented as the pair (alpha,beta), as<br>
*> there is a reasonable interpretation for beta=0, and even for both<br>
*> being zero.<br>
*><br>
*> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)<br>
*> of (A,B) satisfies<br>
*><br>
*>                  A * v(j) = lambda(j) * B * v(j).<br>
*><br>
*> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)<br>
*> of (A,B) satisfies<br>
*><br>
*>                  u(j)**H * A  = lambda(j) * u(j)**H * B .<br>
*><br>
*> where u(j)**H is the conjugate-transpose of u(j).<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBVL<br>
*> \verbatim<br>
*>          JOBVL is CHARACTER*1<br>
*>          = 'N':  do not compute the left generalized eigenvectors;<br>
*>          = 'V':  compute the left generalized eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVR<br>
*> \verbatim<br>
*>          JOBVR is CHARACTER*1<br>
*>          = 'N':  do not compute the right generalized eigenvectors;<br>
*>          = 'V':  compute the right generalized eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A, B, VL, and VR.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the matrix A in the pair (A,B).<br>
*>          On exit, A has been overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB, N)<br>
*>          On entry, the matrix B in the pair (A,B).<br>
*>          On exit, B has been overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will<br>
*>          be the generalized eigenvalues.  If ALPHAI(j) is zero, then<br>
*>          the j-th eigenvalue is real; if positive, then the j-th and<br>
*>          (j+1)-st eigenvalues are a complex conjugate pair, with<br>
*>          ALPHAI(j+1) negative.<br>
*><br>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)<br>
*>          may easily over- or underflow, and BETA(j) may even be zero.<br>
*>          Thus, the user should avoid naively computing the ratio<br>
*>          alpha/beta.  However, ALPHAR and ALPHAI will be always less<br>
*>          than and usually comparable with norm(A) in magnitude, and<br>
*>          BETA always less than and usually comparable with norm(B).<br>
*> \endverbatim<br>
*><br>
*> \param[out] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION array, dimension (LDVL,N)<br>
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one<br>
*>          after another in the columns of VL, in the same order as<br>
*>          their eigenvalues. If the j-th eigenvalue is real, then<br>
*>          u(j) = VL(:,j), the j-th column of VL. If the j-th and<br>
*>          (j+1)-th eigenvalues form a complex conjugate pair, then<br>
*>          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).<br>
*>          Each eigenvector is scaled so the largest component has<br>
*>          abs(real part)+abs(imag. part)=1.<br>
*>          Not referenced if JOBVL = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the matrix VL. LDVL >= 1, and<br>
*>          if JOBVL = 'V', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VR<br>
*> \verbatim<br>
*>          VR is DOUBLE PRECISION array, dimension (LDVR,N)<br>
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one<br>
*>          after another in the columns of VR, in the same order as<br>
*>          their eigenvalues. If the j-th eigenvalue is real, then<br>
*>          v(j) = VR(:,j), the j-th column of VR. If the j-th and<br>
*>          (j+1)-th eigenvalues form a complex conjugate pair, then<br>
*>          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).<br>
*>          Each eigenvector is scaled so the largest component has<br>
*>          abs(real part)+abs(imag. part)=1.<br>
*>          Not referenced if JOBVR = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the matrix VR. LDVR >= 1, and<br>
*>          if JOBVR = 'V', LDVR >= N.<br>
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
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          = 1,...,N:<br>
*>                The QZ iteration failed.  No eigenvectors have been<br>
*>                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)<br>
*>                should be correct for j=INFO+1,...,N.<br>
*>          > N:  =N+1: other than QZ iteration failed in DHGEQZ.<br>
*>                =N+2: error return from DTGEVC.<br>
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
*> \date January 2015<br>
*<br>
*> \ingroup doubleGEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dggev3_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGGEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB,<br>
*                          ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO,<br>
*                          IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE,<br>
*                          RCONDV, WORK, LWORK, IWORK, BWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          BALANC, JOBVL, JOBVR, SENSE<br>
*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N<br>
*       DOUBLE PRECISION   ABNRM, BBNRM<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            BWORK( * )<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),<br>
*      $                   B( LDB, * ), BETA( * ), LSCALE( * ),<br>
*      $                   RCONDE( * ), RCONDV( * ), RSCALE( * ),<br>
*      $                   VL( LDVL, * ), VR( LDVR, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)<br>
*> the generalized eigenvalues, and optionally, the left and/or right<br>
*> generalized eigenvectors.<br>
*><br>
*> Optionally also, it computes a balancing transformation to improve<br>
*> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,<br>
*> LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for<br>
*> the eigenvalues (RCONDE), and reciprocal condition numbers for the<br>
*> right eigenvectors (RCONDV).<br>
*><br>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar<br>
*> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is<br>
*> singular. It is usually represented as the pair (alpha,beta), as<br>
*> there is a reasonable interpretation for beta=0, and even for both<br>
*> being zero.<br>
*><br>
*> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)<br>
*> of (A,B) satisfies<br>
*><br>
*>                  A * v(j) = lambda(j) * B * v(j) .<br>
*><br>
*> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)<br>
*> of (A,B) satisfies<br>
*><br>
*>                  u(j)**H * A  = lambda(j) * u(j)**H * B.<br>
*><br>
*> where u(j)**H is the conjugate-transpose of u(j).<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] BALANC<br>
*> \verbatim<br>
*>          BALANC is CHARACTER*1<br>
*>          Specifies the balance option to be performed.<br>
*>          = 'N':  do not diagonally scale or permute;<br>
*>          = 'P':  permute only;<br>
*>          = 'S':  scale only;<br>
*>          = 'B':  both permute and scale.<br>
*>          Computed reciprocal condition numbers will be for the<br>
*>          matrices after permuting and/or balancing. Permuting does<br>
*>          not change condition numbers (in exact arithmetic), but<br>
*>          balancing does.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVL<br>
*> \verbatim<br>
*>          JOBVL is CHARACTER*1<br>
*>          = 'N':  do not compute the left generalized eigenvectors;<br>
*>          = 'V':  compute the left generalized eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBVR<br>
*> \verbatim<br>
*>          JOBVR is CHARACTER*1<br>
*>          = 'N':  do not compute the right generalized eigenvectors;<br>
*>          = 'V':  compute the right generalized eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SENSE<br>
*> \verbatim<br>
*>          SENSE is CHARACTER*1<br>
*>          Determines which reciprocal condition numbers are computed.<br>
*>          = 'N': none are computed;<br>
*>          = 'E': computed for eigenvalues only;<br>
*>          = 'V': computed for eigenvectors only;<br>
*>          = 'B': computed for eigenvalues and eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A, B, VL, and VR.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the matrix A in the pair (A,B).<br>
*>          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'<br>
*>          or both, then A contains the first part of the real Schur<br>
*>          form of the "balanced" versions of the input A and B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB, N)<br>
*>          On entry, the matrix B in the pair (A,B).<br>
*>          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'<br>
*>          or both, then B contains the second part of the real Schur<br>
*>          form of the "balanced" versions of the input A and B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will<br>
*>          be the generalized eigenvalues.  If ALPHAI(j) is zero, then<br>
*>          the j-th eigenvalue is real; if positive, then the j-th and<br>
*>          (j+1)-st eigenvalues are a complex conjugate pair, with<br>
*>          ALPHAI(j+1) negative.<br>
*><br>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)<br>
*>          may easily over- or underflow, and BETA(j) may even be zero.<br>
*>          Thus, the user should avoid naively computing the ratio<br>
*>          ALPHA/BETA. However, ALPHAR and ALPHAI will be always less<br>
*>          than and usually comparable with norm(A) in magnitude, and<br>
*>          BETA always less than and usually comparable with norm(B).<br>
*> \endverbatim<br>
*><br>
*> \param[out] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION array, dimension (LDVL,N)<br>
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one<br>
*>          after another in the columns of VL, in the same order as<br>
*>          their eigenvalues. If the j-th eigenvalue is real, then<br>
*>          u(j) = VL(:,j), the j-th column of VL. If the j-th and<br>
*>          (j+1)-th eigenvalues form a complex conjugate pair, then<br>
*>          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).<br>
*>          Each eigenvector will be scaled so the largest component have<br>
*>          abs(real part) + abs(imag. part) = 1.<br>
*>          Not referenced if JOBVL = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the matrix VL. LDVL >= 1, and<br>
*>          if JOBVL = 'V', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] VR<br>
*> \verbatim<br>
*>          VR is DOUBLE PRECISION array, dimension (LDVR,N)<br>
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one<br>
*>          after another in the columns of VR, in the same order as<br>
*>          their eigenvalues. If the j-th eigenvalue is real, then<br>
*>          v(j) = VR(:,j), the j-th column of VR. If the j-th and<br>
*>          (j+1)-th eigenvalues form a complex conjugate pair, then<br>
*>          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).<br>
*>          Each eigenvector will be scaled so the largest component have<br>
*>          abs(real part) + abs(imag. part) = 1.<br>
*>          Not referenced if JOBVR = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the matrix VR. LDVR >= 1, and<br>
*>          if JOBVR = 'V', LDVR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ILO<br>
*> \verbatim<br>
*>          ILO is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[out] IHI<br>
*> \verbatim<br>
*>          IHI is INTEGER<br>
*>          ILO and IHI are integer values such that on exit<br>
*>          A(i,j) = 0 and B(i,j) = 0 if i > j and<br>
*>          j = 1,...,ILO-1 or i = IHI+1,...,N.<br>
*>          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] LSCALE<br>
*> \verbatim<br>
*>          LSCALE is DOUBLE PRECISION array, dimension (N)<br>
*>          Details of the permutations and scaling factors applied<br>
*>          to the left side of A and B.  If PL(j) is the index of the<br>
*>          row interchanged with row j, and DL(j) is the scaling<br>
*>          factor applied to row j, then<br>
*>            LSCALE(j) = PL(j)  for j = 1,...,ILO-1<br>
*>                      = DL(j)  for j = ILO,...,IHI<br>
*>                      = PL(j)  for j = IHI+1,...,N.<br>
*>          The order in which the interchanges are made is N to IHI+1,<br>
*>          then 1 to ILO-1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RSCALE<br>
*> \verbatim<br>
*>          RSCALE is DOUBLE PRECISION array, dimension (N)<br>
*>          Details of the permutations and scaling factors applied<br>
*>          to the right side of A and B.  If PR(j) is the index of the<br>
*>          column interchanged with column j, and DR(j) is the scaling<br>
*>          factor applied to column j, then<br>
*>            RSCALE(j) = PR(j)  for j = 1,...,ILO-1<br>
*>                      = DR(j)  for j = ILO,...,IHI<br>
*>                      = PR(j)  for j = IHI+1,...,N<br>
*>          The order in which the interchanges are made is N to IHI+1,<br>
*>          then 1 to ILO-1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ABNRM<br>
*> \verbatim<br>
*>          ABNRM is DOUBLE PRECISION<br>
*>          The one-norm of the balanced matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BBNRM<br>
*> \verbatim<br>
*>          BBNRM is DOUBLE PRECISION<br>
*>          The one-norm of the balanced matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCONDE<br>
*> \verbatim<br>
*>          RCONDE is DOUBLE PRECISION array, dimension (N)<br>
*>          If SENSE = 'E' or 'B', the reciprocal condition numbers of<br>
*>          the eigenvalues, stored in consecutive elements of the array.<br>
*>          For a complex conjugate pair of eigenvalues two consecutive<br>
*>          elements of RCONDE are set to the same value. Thus RCONDE(j),<br>
*>          RCONDV(j), and the j-th columns of VL and VR all correspond<br>
*>          to the j-th eigenpair.<br>
*>          If SENSE = 'N or 'V', RCONDE is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCONDV<br>
*> \verbatim<br>
*>          RCONDV is DOUBLE PRECISION array, dimension (N)<br>
*>          If SENSE = 'V' or 'B', the estimated reciprocal condition<br>
*>          numbers of the eigenvectors, stored in consecutive elements<br>
*>          of the array. For a complex eigenvector two consecutive<br>
*>          elements of RCONDV are set to the same value. If the<br>
*>          eigenvalues cannot be reordered to compute RCONDV(j),<br>
*>          RCONDV(j) is set to 0; this can only occur when the true<br>
*>          value would be very small anyway.<br>
*>          If SENSE = 'N' or 'E', RCONDV is not referenced.<br>
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
*>          The dimension of the array WORK. LWORK >= max(1,2*N).<br>
*>          If BALANC = 'S' or 'B', or JOBVL = 'V', or JOBVR = 'V',<br>
*>          LWORK >= max(1,6*N).<br>
*>          If SENSE = 'E' or 'B', LWORK >= max(1,10*N).<br>
*>          If SENSE = 'V' or 'B', LWORK >= 2*N*N+8*N+16.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N+6)<br>
*>          If SENSE = 'E', IWORK is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BWORK<br>
*> \verbatim<br>
*>          BWORK is LOGICAL array, dimension (N)<br>
*>          If SENSE = 'N', BWORK is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          = 1,...,N:<br>
*>                The QZ iteration failed.  No eigenvectors have been<br>
*>                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)<br>
*>                should be correct for j=INFO+1,...,N.<br>
*>          > N:  =N+1: other than QZ iteration failed in DHGEQZ.<br>
*>                =N+2: error return from DTGEVC.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGEeigen<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Balancing a matrix pair (A,B) includes, first, permuting rows and<br>
*>  columns to isolate eigenvalues, second, applying diagonal similarity<br>
*>  transformation to the rows and columns to make the rows and columns<br>
*>  as close in norm as possible. The computed reciprocal condition<br>
*>  numbers correspond to the balanced matrix. Permuting rows and columns<br>
*>  will not change the condition numbers (in exact arithmetic) but<br>
*>  diagonal scaling will.  For further explanation of balancing, see<br>
*>  section 4.11.1.2 of LAPACK Users' Guide.<br>
*><br>
*>  An approximate error bound on the chordal distance between the i-th<br>
*>  computed generalized eigenvalue w and the corresponding exact<br>
*>  eigenvalue lambda is<br>
*><br>
*>       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)<br>
*><br>
*>  An approximate error bound for the angle between the i-th computed<br>
*>  eigenvector VL(i) or VR(i) is given by<br>
*><br>
*>       EPS * norm(ABNRM, BBNRM) / DIF(i).<br>
*><br>
*>  For further explanation of the reciprocal condition numbers RCONDE<br>
*>  and RCONDV, see section 4.11 of LAPACK User's Guide.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dggevx_(CHARACTER BALANC,CHARACTER JOBVL,CHARACTER JOBVR,CHARACTER SENSE,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER ILO,INTEGER IHI,double[] LSCALE,double[] RSCALE,DOUBLE ABNRM,DOUBLE BBNRM,double[] RCONDE,double[] RCONDV,double[] WORK,INTEGER LWORK,int[] IWORK,boolean[] BWORK,INTEGER INFO);
/**
*> \brief \b DGGGLM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGGLM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggglm.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggglm.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggglm.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), D( * ), WORK( * ),<br>
*      $                   X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGGLM solves a general Gauss-Markov linear model (GLM) problem:<br>
*><br>
*>         minimize || y ||_2   subject to   d = A*x + B*y<br>
*>             x<br>
*><br>
*> where A is an N-by-M matrix, B is an N-by-P matrix, and d is a<br>
*> given N-vector. It is assumed that M <= N <= M+P, and<br>
*><br>
*>            rank(A) = M    and    rank( A B ) = N.<br>
*><br>
*> Under these assumptions, the constrained equation is always<br>
*> consistent, and there is a unique solution x and a minimal 2-norm<br>
*> solution y, which is obtained using a generalized QR factorization<br>
*> of the matrices (A, B) given by<br>
*><br>
*>    A = Q*(R),   B = Q*T*Z.<br>
*>          (0)<br>
*><br>
*> In particular, if matrix B is square nonsingular, then the problem<br>
*> GLM is equivalent to the following weighted linear least squares<br>
*> problem<br>
*><br>
*>              minimize || inv(B)*(d-A*x) ||_2<br>
*>                  x<br>
*><br>
*> where inv(B) denotes the inverse of B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of rows of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of columns of the matrix A.  0 <= M <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of columns of the matrix B.  P >= N-M.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,M)<br>
*>          On entry, the N-by-M matrix A.<br>
*>          On exit, the upper triangular part of the array A contains<br>
*>          the M-by-M upper triangular matrix R.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,P)<br>
*>          On entry, the N-by-P matrix B.<br>
*>          On exit, if N <= P, the upper triangle of the subarray<br>
*>          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;<br>
*>          if N > P, the elements on and above the (N-P)th subdiagonal<br>
*>          contain the N-by-P upper trapezoidal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, D is the left hand side of the GLM equation.<br>
*>          On exit, D is destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (M)<br>
*> \endverbatim<br>
*><br>
*> \param[out] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array, dimension (P)<br>
*><br>
*>          On exit, X and Y are the solutions of the GLM problem.<br>
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
*>          The dimension of the array WORK. LWORK >= max(1,N+M+P).<br>
*>          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB,<br>
*>          where NB is an upper bound for the optimal blocksizes for<br>
*>          DGEQRF, SGERQF, DORMQR and SORMRQ.<br>
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
*>          = 1:  the upper triangular factor R associated with A in the<br>
*>                generalized QR factorization of the pair (A, B) is<br>
*>                singular, so that rank(A) < M; the least squares<br>
*>                solution could not be computed.<br>
*>          = 2:  the bottom (N-M) by (N-M) part of the upper trapezoidal<br>
*>                factor T associated with B in the generalized QR<br>
*>                factorization of the pair (A, B) is singular, so that<br>
*>                rank( A B ) < N; the least squares solution could not<br>
*>                be computed.<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void dggglm_(INTEGER N,INTEGER M,INTEGER P,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] D,double[] X,double[] Y,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGGHD3<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download DGGHD3 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgghd3.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgghd3.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgghd3.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGHD3( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,<br>
*                          LDQ, Z, LDZ, WORK, LWORK, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ, COMPZ<br>
*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),<br>
*      $                   Z( LDZ, * ), WORK( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGHD3 reduces a pair of real matrices (A,B) to generalized upper<br>
*> Hessenberg form using orthogonal transformations, where A is a<br>
*> general matrix and B is upper triangular.  The form of the<br>
*> generalized eigenvalue problem is<br>
*>    A*x = lambda*B*x,<br>
*> and B is typically made upper triangular by computing its QR<br>
*> factorization and moving the orthogonal matrix Q to the left side<br>
*> of the equation.<br>
*><br>
*> This subroutine simultaneously reduces A to a Hessenberg matrix H:<br>
*>    Q**T*A*Z = H<br>
*> and transforms B to another upper triangular matrix T:<br>
*>    Q**T*B*Z = T<br>
*> in order to reduce the problem to its standard form<br>
*>    H*y = lambda*T*y<br>
*> where y = Z**T*x.<br>
*><br>
*> The orthogonal matrices Q and Z are determined as products of Givens<br>
*> rotations.  They may either be formed explicitly, or they may be<br>
*> postmultiplied into input matrices Q1 and Z1, so that<br>
*><br>
*>      Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T<br>
*><br>
*>      Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T<br>
*><br>
*> If Q1 is the orthogonal matrix from the QR factorization of B in the<br>
*> original equation A*x = lambda*B*x, then DGGHD3 reduces the original<br>
*> problem to generalized Hessenberg form.<br>
*><br>
*> This is a blocked variant of DGGHRD, using matrix-matrix<br>
*> multiplications for parts of the computation to enhance performance.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] COMPQ<br>
*> \verbatim<br>
*>          COMPQ is CHARACTER*1<br>
*>          = 'N': do not compute Q;<br>
*>          = 'I': Q is initialized to the unit matrix, and the<br>
*>                 orthogonal matrix Q is returned;<br>
*>          = 'V': Q must contain an orthogonal matrix Q1 on entry,<br>
*>                 and the product Q1*Q is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>          = 'N': do not compute Z;<br>
*>          = 'I': Z is initialized to the unit matrix, and the<br>
*>                 orthogonal matrix Z is returned;<br>
*>          = 'V': Z must contain an orthogonal matrix Z1 on entry,<br>
*>                 and the product Z1*Z is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
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
*>          ILO and IHI mark the rows and columns of A which are to be<br>
*>          reduced.  It is assumed that A is already upper triangular<br>
*>          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are<br>
*>          normally set by a previous call to DGGBAL; otherwise they<br>
*>          should be set to 1 and N respectively.<br>
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the N-by-N general matrix to be reduced.<br>
*>          On exit, the upper triangle and the first subdiagonal of A<br>
*>          are overwritten with the upper Hessenberg matrix H, and the<br>
*>          rest is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB, N)<br>
*>          On entry, the N-by-N upper triangular matrix B.<br>
*>          On exit, the upper triangular matrix T = Q**T B Z.  The<br>
*>          elements below the diagonal are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)<br>
*>          On entry, if COMPQ = 'V', the orthogonal matrix Q1,<br>
*>          typically from the QR factorization of B.<br>
*>          On exit, if COMPQ='I', the orthogonal matrix Q, and if<br>
*>          COMPQ = 'V', the product Q1*Q.<br>
*>          Not referenced if COMPQ='N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.<br>
*>          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)<br>
*>          On entry, if COMPZ = 'V', the orthogonal matrix Z1.<br>
*>          On exit, if COMPZ='I', the orthogonal matrix Z, and if<br>
*>          COMPZ = 'V', the product Z1*Z.<br>
*>          Not referenced if COMPZ='N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.<br>
*>          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in]  LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of the array WORK.  LWORK >= 1.<br>
*>          For optimum performance LWORK >= 6*N*NB, where NB is the<br>
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
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
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
*> \date January 2015<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  This routine reduces A to Hessenberg form and maintains B in<br>
*>  using a blocked variant of Moler and Stewart's original algorithm,<br>
*>  as described by Kagstrom, Kressner, Quintana-Orti, and Quintana-Orti<br>
*>  (BIT 2008).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgghd3_(CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGGHRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGHRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgghrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgghrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgghrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,<br>
*                          LDQ, Z, LDZ, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ, COMPZ<br>
*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGHRD reduces a pair of real matrices (A,B) to generalized upper<br>
*> Hessenberg form using orthogonal transformations, where A is a<br>
*> general matrix and B is upper triangular.  The form of the<br>
*> generalized eigenvalue problem is<br>
*>    A*x = lambda*B*x,<br>
*> and B is typically made upper triangular by computing its QR<br>
*> factorization and moving the orthogonal matrix Q to the left side<br>
*> of the equation.<br>
*><br>
*> This subroutine simultaneously reduces A to a Hessenberg matrix H:<br>
*>    Q**T*A*Z = H<br>
*> and transforms B to another upper triangular matrix T:<br>
*>    Q**T*B*Z = T<br>
*> in order to reduce the problem to its standard form<br>
*>    H*y = lambda*T*y<br>
*> where y = Z**T*x.<br>
*><br>
*> The orthogonal matrices Q and Z are determined as products of Givens<br>
*> rotations.  They may either be formed explicitly, or they may be<br>
*> postmultiplied into input matrices Q1 and Z1, so that<br>
*><br>
*>      Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T<br>
*><br>
*>      Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T<br>
*><br>
*> If Q1 is the orthogonal matrix from the QR factorization of B in the<br>
*> original equation A*x = lambda*B*x, then DGGHRD reduces the original<br>
*> problem to generalized Hessenberg form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] COMPQ<br>
*> \verbatim<br>
*>          COMPQ is CHARACTER*1<br>
*>          = 'N': do not compute Q;<br>
*>          = 'I': Q is initialized to the unit matrix, and the<br>
*>                 orthogonal matrix Q is returned;<br>
*>          = 'V': Q must contain an orthogonal matrix Q1 on entry,<br>
*>                 and the product Q1*Q is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>          = 'N': do not compute Z;<br>
*>          = 'I': Z is initialized to the unit matrix, and the<br>
*>                 orthogonal matrix Z is returned;<br>
*>          = 'V': Z must contain an orthogonal matrix Z1 on entry,<br>
*>                 and the product Z1*Z is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
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
*>          ILO and IHI mark the rows and columns of A which are to be<br>
*>          reduced.  It is assumed that A is already upper triangular<br>
*>          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are<br>
*>          normally set by a previous call to DGGBAL; otherwise they<br>
*>          should be set to 1 and N respectively.<br>
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA, N)<br>
*>          On entry, the N-by-N general matrix to be reduced.<br>
*>          On exit, the upper triangle and the first subdiagonal of A<br>
*>          are overwritten with the upper Hessenberg matrix H, and the<br>
*>          rest is set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB, N)<br>
*>          On entry, the N-by-N upper triangular matrix B.<br>
*>          On exit, the upper triangular matrix T = Q**T B Z.  The<br>
*>          elements below the diagonal are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)<br>
*>          On entry, if COMPQ = 'V', the orthogonal matrix Q1,<br>
*>          typically from the QR factorization of B.<br>
*>          On exit, if COMPQ='I', the orthogonal matrix Q, and if<br>
*>          COMPQ = 'V', the product Q1*Q.<br>
*>          Not referenced if COMPQ='N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.<br>
*>          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)<br>
*>          On entry, if COMPZ = 'V', the orthogonal matrix Z1.<br>
*>          On exit, if COMPZ='I', the orthogonal matrix Z, and if<br>
*>          COMPZ = 'V', the product Z1*Z.<br>
*>          Not referenced if COMPZ='N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.<br>
*>          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  This routine reduces A to Hessenberg and B to triangular form by<br>
*>  an unblocked reduction, as described in _Matrix_Computations_,<br>
*>  by Golub and Van Loan (Johns Hopkins Press.)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgghrd_(CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,INTEGER INFO);
/**
*> \brief <b> DGGLSE solves overdetermined or underdetermined systems for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGLSE + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgglse.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgglse.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgglse.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( * ), D( * ),<br>
*      $                   WORK( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGLSE solves the linear equality-constrained least squares (LSE)<br>
*> problem:<br>
*><br>
*>         minimize || c - A*x ||_2   subject to   B*x = d<br>
*><br>
*> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given<br>
*> M-vector, and d is a given P-vector. It is assumed that<br>
*> P <= N <= M+P, and<br>
*><br>
*>          rank(B) = P and  rank( (A) ) = N.<br>
*>                               ( (B) )<br>
*><br>
*> These conditions ensure that the LSE problem has a unique solution,<br>
*> which is obtained using a generalized RQ factorization of the<br>
*> matrices (B, A) given by<br>
*><br>
*>    B = (0 R)*Q,   A = Z*T*Q.<br>
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
*>          The number of columns of the matrices A and B. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of rows of the matrix B. 0 <= P <= N <= M+P.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the min(M,N)-by-N upper trapezoidal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,N)<br>
*>          On entry, the P-by-N matrix B.<br>
*>          On exit, the upper triangle of the subarray B(1:P,N-P+1:N)<br>
*>          contains the P-by-P upper triangular matrix R.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array, dimension (M)<br>
*>          On entry, C contains the right hand side vector for the<br>
*>          least squares part of the LSE problem.<br>
*>          On exit, the residual sum of squares for the solution<br>
*>          is given by the sum of squares of elements N-P+1 to M of<br>
*>          vector C.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (P)<br>
*>          On entry, D contains the right hand side vector for the<br>
*>          constrained equation.<br>
*>          On exit, D is destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (N)<br>
*>          On exit, X is the solution of the LSE problem.<br>
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
*>          The dimension of the array WORK. LWORK >= max(1,M+N+P).<br>
*>          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,<br>
*>          where NB is an upper bound for the optimal blocksizes for<br>
*>          DGEQRF, SGERQF, DORMQR and SORMRQ.<br>
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
*>          = 1:  the upper triangular factor R associated with B in the<br>
*>                generalized RQ factorization of the pair (B, A) is<br>
*>                singular, so that rank(B) < P; the least squares<br>
*>                solution could not be computed.<br>
*>          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor<br>
*>                T associated with A in the generalized RQ factorization<br>
*>                of the pair (B, A) is singular, so that<br>
*>                rank( (A) ) < N; the least squares solution could not<br>
*>                    ( (B) )<br>
*>                be computed.<br>
*> \endverbatim<br>
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
*> \ingroup doubleOTHERsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgglse_(INTEGER M,INTEGER N,INTEGER P,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] C,double[] D,double[] X,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGGQRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGQRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggqrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggqrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggqrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK,<br>
*                          LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGQRF computes a generalized QR factorization of an N-by-M matrix A<br>
*> and an N-by-P matrix B:<br>
*><br>
*>             A = Q*R,        B = Q*T*Z,<br>
*><br>
*> where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal<br>
*> matrix, and R and T assume one of the forms:<br>
*><br>
*> if N >= M,  R = ( R11 ) M  ,   or if N < M,  R = ( R11  R12 ) N,<br>
*>                 (  0  ) N-M                         N   M-N<br>
*>                    M<br>
*><br>
*> where R11 is upper triangular, and<br>
*><br>
*> if N <= P,  T = ( 0  T12 ) N,   or if N > P,  T = ( T11 ) N-P,<br>
*>                  P-N  N                           ( T21 ) P<br>
*>                                                      P<br>
*><br>
*> where T12 or T21 is upper triangular.<br>
*><br>
*> In particular, if B is square and nonsingular, the GQR factorization<br>
*> of A and B implicitly gives the QR factorization of inv(B)*A:<br>
*><br>
*>              inv(B)*A = Z**T*(inv(T)*R)<br>
*><br>
*> where inv(B) denotes the inverse of the matrix B, and Z**T denotes the<br>
*> transpose of the matrix Z.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of rows of the matrices A and B. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of columns of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of columns of the matrix B.  P >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,M)<br>
*>          On entry, the N-by-M matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the min(N,M)-by-M upper trapezoidal matrix R (R is<br>
*>          upper triangular if N >= M); the elements below the diagonal,<br>
*>          with the array TAUA, represent the orthogonal matrix Q as a<br>
*>          product of min(N,M) elementary reflectors (see Further<br>
*>          Details).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUA<br>
*> \verbatim<br>
*>          TAUA is DOUBLE PRECISION array, dimension (min(N,M))<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix Q (see Further Details).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,P)<br>
*>          On entry, the N-by-P matrix B.<br>
*>          On exit, if N <= P, the upper triangle of the subarray<br>
*>          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;<br>
*>          if N > P, the elements on and above the (N-P)-th subdiagonal<br>
*>          contain the N-by-P upper trapezoidal matrix T; the remaining<br>
*>          elements, with the array TAUB, represent the orthogonal<br>
*>          matrix Z as a product of elementary reflectors (see Further<br>
*>          Details).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUB<br>
*> \verbatim<br>
*>          TAUB is DOUBLE PRECISION array, dimension (min(N,P))<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix Z (see Further Details).<br>
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
*>          The dimension of the array WORK. LWORK >= max(1,N,M,P).<br>
*>          For optimum performance LWORK >= max(N,M,P)*max(NB1,NB2,NB3),<br>
*>          where NB1 is the optimal blocksize for the QR factorization<br>
*>          of an N-by-M matrix, NB2 is the optimal blocksize for the<br>
*>          RQ factorization of an N-by-P matrix, and NB3 is the optimal<br>
*>          blocksize for a call of DORMQR.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(n,m).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - taua * v * v**T<br>
*><br>
*>  where taua is a real scalar, and v is a real vector with<br>
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),<br>
*>  and taua in TAUA(i).<br>
*>  To form Q explicitly, use LAPACK subroutine DORGQR.<br>
*>  To use Q to update another matrix, use LAPACK subroutine DORMQR.<br>
*><br>
*>  The matrix Z is represented as a product of elementary reflectors<br>
*><br>
*>     Z = H(1) H(2) . . . H(k), where k = min(n,p).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - taub * v * v**T<br>
*><br>
*>  where taub is a real scalar, and v is a real vector with<br>
*>  v(p-k+i+1:p) = 0 and v(p-k+i) = 1; v(1:p-k+i-1) is stored on exit in<br>
*>  B(n-k+i,1:p-k+i-1), and taub in TAUB(i).<br>
*>  To form Z explicitly, use LAPACK subroutine DORGRQ.<br>
*>  To use Z to update another matrix, use LAPACK subroutine DORMRQ.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dggqrf_(INTEGER N,INTEGER M,INTEGER P,double[] A,INTEGER LDA,double[] TAUA,double[] B,INTEGER LDB,double[] TAUB,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGGRQF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGRQF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggrqf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggrqf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggrqf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGRQF( M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK,<br>
*                          LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGRQF computes a generalized RQ factorization of an M-by-N matrix A<br>
*> and a P-by-N matrix B:<br>
*><br>
*>             A = R*Q,        B = Z*T*Q,<br>
*><br>
*> where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal<br>
*> matrix, and R and T assume one of the forms:<br>
*><br>
*> if M <= N,  R = ( 0  R12 ) M,   or if M > N,  R = ( R11 ) M-N,<br>
*>                  N-M  M                           ( R21 ) N<br>
*>                                                      N<br>
*><br>
*> where R12 or R21 is upper triangular, and<br>
*><br>
*> if P >= N,  T = ( T11 ) N  ,   or if P < N,  T = ( T11  T12 ) P,<br>
*>                 (  0  ) P-N                         P   N-P<br>
*>                    N<br>
*><br>
*> where T11 is upper triangular.<br>
*><br>
*> In particular, if B is square and nonsingular, the GRQ factorization<br>
*> of A and B implicitly gives the RQ factorization of A*inv(B):<br>
*><br>
*>              A*inv(B) = (R*inv(T))*Z**T<br>
*><br>
*> where inv(B) denotes the inverse of the matrix B, and Z**T denotes the<br>
*> transpose of the matrix Z.<br>
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
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of rows of the matrix B.  P >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrices A and B. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, if M <= N, the upper triangle of the subarray<br>
*>          A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R;<br>
*>          if M > N, the elements on and above the (M-N)-th subdiagonal<br>
*>          contain the M-by-N upper trapezoidal matrix R; the remaining<br>
*>          elements, with the array TAUA, represent the orthogonal<br>
*>          matrix Q as a product of elementary reflectors (see Further<br>
*>          Details).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUA<br>
*> \verbatim<br>
*>          TAUA is DOUBLE PRECISION array, dimension (min(M,N))<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix Q (see Further Details).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,N)<br>
*>          On entry, the P-by-N matrix B.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the min(P,N)-by-N upper trapezoidal matrix T (T is<br>
*>          upper triangular if P >= N); the elements below the diagonal,<br>
*>          with the array TAUB, represent the orthogonal matrix Z as a<br>
*>          product of elementary reflectors (see Further Details).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAUB<br>
*> \verbatim<br>
*>          TAUB is DOUBLE PRECISION array, dimension (min(P,N))<br>
*>          The scalar factors of the elementary reflectors which<br>
*>          represent the orthogonal matrix Z (see Further Details).<br>
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
*>          The dimension of the array WORK. LWORK >= max(1,N,M,P).<br>
*>          For optimum performance LWORK >= max(N,M,P)*max(NB1,NB2,NB3),<br>
*>          where NB1 is the optimal blocksize for the RQ factorization<br>
*>          of an M-by-N matrix, NB2 is the optimal blocksize for the<br>
*>          QR factorization of a P-by-N matrix, and NB3 is the optimal<br>
*>          blocksize for a call of DORMRQ.<br>
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
*>          < 0:  if INF0= -i, the i-th argument had an illegal value.<br>
*> \endverbatim<br>
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
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix Q is represented as a product of elementary reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - taua * v * v**T<br>
*><br>
*>  where taua is a real scalar, and v is a real vector with<br>
*>  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in<br>
*>  A(m-k+i,1:n-k+i-1), and taua in TAUA(i).<br>
*>  To form Q explicitly, use LAPACK subroutine DORGRQ.<br>
*>  To use Q to update another matrix, use LAPACK subroutine DORMRQ.<br>
*><br>
*>  The matrix Z is represented as a product of elementary reflectors<br>
*><br>
*>     Z = H(1) H(2) . . . H(k), where k = min(p,n).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - taub * v * v**T<br>
*><br>
*>  where taub is a real scalar, and v is a real vector with<br>
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:p) is stored on exit in B(i+1:p,i),<br>
*>  and taub in TAUB(i).<br>
*>  To form Z explicitly, use LAPACK subroutine DORGQR.<br>
*>  To use Z to update another matrix, use LAPACK subroutine DORMQR.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dggrqf_(INTEGER M,INTEGER P,INTEGER N,double[] A,INTEGER LDA,double[] TAUA,double[] B,INTEGER LDB,double[] TAUB,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> DGGSVD3 computes the singular value decomposition (SVD) for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGSVD3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggsvd3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggsvd3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggsvd3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,<br>
*                           LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK,<br>
*                           LWORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBQ, JOBU, JOBV<br>
*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), ALPHA( * ), B( LDB, * ),<br>
*      $                   BETA( * ), Q( LDQ, * ), U( LDU, * ),<br>
*      $                   V( LDV, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGSVD3 computes the generalized singular value decomposition (GSVD)<br>
*> of an M-by-N real matrix A and P-by-N real matrix B:<br>
*><br>
*>       U**T*A*Q = D1*( 0 R ),    V**T*B*Q = D2*( 0 R )<br>
*><br>
*> where U, V and Q are orthogonal matrices.<br>
*> Let K+L = the effective numerical rank of the matrix (A**T,B**T)**T,<br>
*> then R is a K+L-by-K+L nonsingular upper triangular matrix, D1 and<br>
*> D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and of the<br>
*> following structures, respectively:<br>
*><br>
*> If M-K-L >= 0,<br>
*><br>
*>                     K  L<br>
*>        D1 =     K ( I  0 )<br>
*>                 L ( 0  C )<br>
*>             M-K-L ( 0  0 )<br>
*><br>
*>                   K  L<br>
*>        D2 =   L ( 0  S )<br>
*>             P-L ( 0  0 )<br>
*><br>
*>                 N-K-L  K    L<br>
*>   ( 0 R ) = K (  0   R11  R12 )<br>
*>             L (  0    0   R22 )<br>
*><br>
*> where<br>
*><br>
*>   C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),<br>
*>   S = diag( BETA(K+1),  ... , BETA(K+L) ),<br>
*>   C**2 + S**2 = I.<br>
*><br>
*>   R is stored in A(1:K+L,N-K-L+1:N) on exit.<br>
*><br>
*> If M-K-L < 0,<br>
*><br>
*>                   K M-K K+L-M<br>
*>        D1 =   K ( I  0    0   )<br>
*>             M-K ( 0  C    0   )<br>
*><br>
*>                     K M-K K+L-M<br>
*>        D2 =   M-K ( 0  S    0  )<br>
*>             K+L-M ( 0  0    I  )<br>
*>               P-L ( 0  0    0  )<br>
*><br>
*>                    N-K-L  K   M-K  K+L-M<br>
*>   ( 0 R ) =     K ( 0    R11  R12  R13  )<br>
*>               M-K ( 0     0   R22  R23  )<br>
*>             K+L-M ( 0     0    0   R33  )<br>
*><br>
*> where<br>
*><br>
*>   C = diag( ALPHA(K+1), ... , ALPHA(M) ),<br>
*>   S = diag( BETA(K+1),  ... , BETA(M) ),<br>
*>   C**2 + S**2 = I.<br>
*><br>
*>   (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored<br>
*>   ( 0  R22 R23 )<br>
*>   in B(M-K+1:L,N+M-K-L+1:N) on exit.<br>
*><br>
*> The routine computes C, S, R, and optionally the orthogonal<br>
*> transformation matrices U, V and Q.<br>
*><br>
*> In particular, if B is an N-by-N nonsingular matrix, then the GSVD of<br>
*> A and B implicitly gives the SVD of A*inv(B):<br>
*>                      A*inv(B) = U*(D1*inv(D2))*V**T.<br>
*> If ( A**T,B**T)**T  has orthonormal columns, then the GSVD of A and B is<br>
*> also equal to the CS decomposition of A and B. Furthermore, the GSVD<br>
*> can be used to derive the solution of the eigenvalue problem:<br>
*>                      A**T*A x = lambda* B**T*B x.<br>
*> In some literature, the GSVD of A and B is presented in the form<br>
*>                  U**T*A*X = ( 0 D1 ),   V**T*B*X = ( 0 D2 )<br>
*> where U and V are orthogonal and X is nonsingular, D1 and D2 are<br>
*> ``diagonal''.  The former GSVD form can be converted to the latter<br>
*> form by taking the nonsingular matrix X as<br>
*><br>
*>                      X = Q*( I   0    )<br>
*>                            ( 0 inv(R) ).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBU<br>
*> \verbatim<br>
*>          JOBU is CHARACTER*1<br>
*>          = 'U':  Orthogonal matrix U is computed;<br>
*>          = 'N':  U is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV<br>
*> \verbatim<br>
*>          JOBV is CHARACTER*1<br>
*>          = 'V':  Orthogonal matrix V is computed;<br>
*>          = 'N':  V is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBQ<br>
*> \verbatim<br>
*>          JOBQ is CHARACTER*1<br>
*>          = 'Q':  Orthogonal matrix Q is computed;<br>
*>          = 'N':  Q is not computed.<br>
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
*>          The number of columns of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of rows of the matrix B.  P >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[out] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*><br>
*>          On exit, K and L specify the dimension of the subblocks<br>
*>          described in Purpose.<br>
*>          K + L = effective numerical rank of (A**T,B**T)**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, A contains the triangular matrix R, or part of R.<br>
*>          See Purpose for details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,N)<br>
*>          On entry, the P-by-N matrix B.<br>
*>          On exit, B contains the triangular matrix R if M-K-L < 0.<br>
*>          See Purpose for details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION array, dimension (N)<br>
*><br>
*>          On exit, ALPHA and BETA contain the generalized singular<br>
*>          value pairs of A and B;<br>
*>            ALPHA(1:K) = 1,<br>
*>            BETA(1:K)  = 0,<br>
*>          and if M-K-L >= 0,<br>
*>            ALPHA(K+1:K+L) = C,<br>
*>            BETA(K+1:K+L)  = S,<br>
*>          or if M-K-L < 0,<br>
*>            ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0<br>
*>            BETA(K+1:M) =S, BETA(M+1:K+L) =1<br>
*>          and<br>
*>            ALPHA(K+L+1:N) = 0<br>
*>            BETA(K+L+1:N)  = 0<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is DOUBLE PRECISION array, dimension (LDU,M)<br>
*>          If JOBU = 'U', U contains the M-by-M orthogonal matrix U.<br>
*>          If JOBU = 'N', U is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>          The leading dimension of the array U. LDU >= max(1,M) if<br>
*>          JOBU = 'U'; LDU >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension (LDV,P)<br>
*>          If JOBV = 'V', V contains the P-by-P orthogonal matrix V.<br>
*>          If JOBV = 'N', V is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V. LDV >= max(1,P) if<br>
*>          JOBV = 'V'; LDV >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension (LDQ,N)<br>
*>          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q.<br>
*>          If JOBQ = 'N', Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q. LDQ >= max(1,N) if<br>
*>          JOBQ = 'Q'; LDQ >= 1 otherwise.<br>
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
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*>          On exit, IWORK stores the sorting information. More<br>
*>          precisely, the following loop will sort ALPHA<br>
*>             for I = K+1, min(M,K+L)<br>
*>                 swap ALPHA(I) and ALPHA(IWORK(I))<br>
*>             endfor<br>
*>          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = 1, the Jacobi-type procedure failed to<br>
*>                converge.  For further details, see subroutine DTGSJA.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  TOLA    DOUBLE PRECISION<br>
*>  TOLB    DOUBLE PRECISION<br>
*>          TOLA and TOLB are the thresholds to determine the effective<br>
*>          rank of (A**T,B**T)**T. Generally, they are set to<br>
*>                   TOLA = MAX(M,N)*norm(A)*MACHEPS,<br>
*>                   TOLB = MAX(P,N)*norm(B)*MACHEPS.<br>
*>          The size of TOLA and TOLB may affect the size of backward<br>
*>          errors of the decomposition.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date August 2015<br>
*<br>
*> \ingroup doubleOTHERsing<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Ming Gu and Huan Ren, Computer Science Division, University of<br>
*>     California at Berkeley, USA<br>
*><br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*>  DGGSVD3 replaces the deprecated subroutine DGGSVD.<br>
*><br>
*  =====================================================================<br>
*/
	public void dggsvd3_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER N,INTEGER P,INTEGER K,INTEGER L,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHA,double[] BETA,double[] U,INTEGER LDU,double[] V,INTEGER LDV,double[] Q,INTEGER LDQ,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGGSVP3<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGGSVP3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggsvp3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggsvp3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggsvp3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,<br>
*                           TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,<br>
*                           IWORK, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBQ, JOBU, JOBV<br>
*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK<br>
*       DOUBLE PRECISION   TOLA, TOLB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),<br>
*      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGGSVP3 computes orthogonal matrices U, V and Q such that<br>
*><br>
*>                    N-K-L  K    L<br>
*>  U**T*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0;<br>
*>                 L ( 0     0   A23 )<br>
*>             M-K-L ( 0     0    0  )<br>
*><br>
*>                  N-K-L  K    L<br>
*>         =     K ( 0    A12  A13 )  if M-K-L < 0;<br>
*>             M-K ( 0     0   A23 )<br>
*><br>
*>                  N-K-L  K    L<br>
*>  V**T*B*Q =   L ( 0     0   B13 )<br>
*>             P-L ( 0     0    0  )<br>
*><br>
*> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular<br>
*> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,<br>
*> otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective<br>
*> numerical rank of the (M+P)-by-N matrix (A**T,B**T)**T. <br>
*><br>
*> This decomposition is the preprocessing step for computing the<br>
*> Generalized Singular Value Decomposition (GSVD), see subroutine<br>
*> DGGSVD3.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBU<br>
*> \verbatim<br>
*>          JOBU is CHARACTER*1<br>
*>          = 'U':  Orthogonal matrix U is computed;<br>
*>          = 'N':  U is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV<br>
*> \verbatim<br>
*>          JOBV is CHARACTER*1<br>
*>          = 'V':  Orthogonal matrix V is computed;<br>
*>          = 'N':  V is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBQ<br>
*> \verbatim<br>
*>          JOBQ is CHARACTER*1<br>
*>          = 'Q':  Orthogonal matrix Q is computed;<br>
*>          = 'N':  Q is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is INTEGER<br>
*>          The number of rows of the matrix B.  P >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, A contains the triangular (or trapezoidal) matrix<br>
*>          described in the Purpose section.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,N)<br>
*>          On entry, the P-by-N matrix B.<br>
*>          On exit, B contains the triangular matrix described in<br>
*>          the Purpose section.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,P).<br>
*> \endverbatim<br>
*><br>
*> \param[in] TOLA<br>
*> \verbatim<br>
*>          TOLA is DOUBLE PRECISION<br>
*> \endverbatim<br>
*><br>
*> \param[in] TOLB<br>
*> \verbatim<br>
*>          TOLB is DOUBLE PRECISION<br>
*><br>
*>          TOLA and TOLB are the thresholds to determine the effective<br>
*>          numerical rank of matrix B and a subblock of A. Generally,<br>
*>          they are set to<br>
*>             TOLA = MAX(M,N)*norm(A)*MACHEPS,<br>
*>             TOLB = MAX(P,N)*norm(B)*MACHEPS.<br>
*>          The size of TOLA and TOLB may affect the size of backward<br>
*>          errors of the decomposition.<br>
*> \endverbatim<br>
*><br>
*> \param[out] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[out] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*><br>
*>          On exit, K and L specify the dimension of the subblocks<br>
*>          described in Purpose section.<br>
*>          K + L = effective numerical rank of (A**T,B**T)**T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] U<br>
*> \verbatim<br>
*>          U is DOUBLE PRECISION array, dimension (LDU,M)<br>
*>          If JOBU = 'U', U contains the orthogonal matrix U.<br>
*>          If JOBU = 'N', U is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDU<br>
*> \verbatim<br>
*>          LDU is INTEGER<br>
*>          The leading dimension of the array U. LDU >= max(1,M) if<br>
*>          JOBU = 'U'; LDU >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension (LDV,P)<br>
*>          If JOBV = 'V', V contains the orthogonal matrix V.<br>
*>          If JOBV = 'N', V is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V. LDV >= max(1,P) if<br>
*>          JOBV = 'V'; LDV >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension (LDQ,N)<br>
*>          If JOBQ = 'Q', Q contains the orthogonal matrix Q.<br>
*>          If JOBQ = 'N', Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q. LDQ >= max(1,N) if<br>
*>          JOBQ = 'Q'; LDQ >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is DOUBLE PRECISION array, dimension (N)<br>
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
*> \date August 2015<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The subroutine uses LAPACK subroutine DGEQP3 for the QR factorization<br>
*>  with column pivoting to detect the effective numerical rank of the<br>
*>  a matrix. It may be replaced by a better rank determination strategy.<br>
*><br>
*>  DGGSVP3 replaces the deprecated subroutine DGGSVP.<br>
*><br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dggsvp3_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER P,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE TOLA,DOUBLE TOLB,INTEGER K,INTEGER L,double[] U,INTEGER LDU,double[] V,INTEGER LDV,double[] Q,INTEGER LDQ,int[] IWORK,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGSVJ0 pre-processor for the routine dgesvj.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGSVJ0 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgsvj0.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgsvj0.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgsvj0.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS,<br>
*                          SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP<br>
*       DOUBLE PRECISION   EPS, SFMIN, TOL<br>
*       CHARACTER*1        JOBV<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), SVA( N ), D( N ), V( LDV, * ),<br>
*      $                   WORK( LWORK )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGSVJ0 is called from DGESVJ as a pre-processor and that is its main<br>
*> purpose. It applies Jacobi rotations in the same way as DGESVJ does, but<br>
*> it does not check convergence (stopping criterion). Few tuning<br>
*> parameters (marked by [TP]) are available for the implementer.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBV<br>
*> \verbatim<br>
*>          JOBV is CHARACTER*1<br>
*>          Specifies whether the output from this procedure is used<br>
*>          to compute the matrix V:<br>
*>          = 'V': the product of the Jacobi rotations is accumulated<br>
*>                 by postmulyiplying the N-by-N array V.<br>
*>                (See the description of V.)<br>
*>          = 'A': the product of the Jacobi rotations is accumulated<br>
*>                 by postmulyiplying the MV-by-N array V.<br>
*>                (See the descriptions of MV and V.)<br>
*>          = 'N': the Jacobi rotations are not accumulated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the input matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the input matrix A.<br>
*>          M >= N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, M-by-N matrix A, such that A*diag(D) represents<br>
*>          the input matrix.<br>
*>          On exit,<br>
*>          A_onexit * D_onexit represents the input matrix A*diag(D)<br>
*>          post-multiplied by a sequence of Jacobi rotations, where the<br>
*>          rotation threshold and the total number of sweeps are given in<br>
*>          TOL and NSWEEP, respectively.<br>
*>          (See the descriptions of D, TOL and NSWEEP.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The array D accumulates the scaling factors from the fast scaled<br>
*>          Jacobi rotations.<br>
*>          On entry, A*diag(D) represents the input matrix.<br>
*>          On exit, A_onexit*diag(D_onexit) represents the input matrix<br>
*>          post-multiplied by a sequence of Jacobi rotations, where the<br>
*>          rotation threshold and the total number of sweeps are given in<br>
*>          TOL and NSWEEP, respectively.<br>
*>          (See the descriptions of A, TOL and NSWEEP.)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SVA<br>
*> \verbatim<br>
*>          SVA is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, SVA contains the Euclidean norms of the columns of<br>
*>          the matrix A*diag(D).<br>
*>          On exit, SVA contains the Euclidean norms of the columns of<br>
*>          the matrix onexit*diag(D_onexit).<br>
*> \endverbatim<br>
*><br>
*> \param[in] MV<br>
*> \verbatim<br>
*>          MV is INTEGER<br>
*>          If JOBV .EQ. 'A', then MV rows of V are post-multipled by a<br>
*>                           sequence of Jacobi rotations.<br>
*>          If JOBV = 'N',   then MV is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension (LDV,N)<br>
*>          If JOBV .EQ. 'V' then N rows of V are post-multipled by a<br>
*>                           sequence of Jacobi rotations.<br>
*>          If JOBV .EQ. 'A' then MV rows of V are post-multipled by a<br>
*>                           sequence of Jacobi rotations.<br>
*>          If JOBV = 'N',   then V is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V,  LDV >= 1.<br>
*>          If JOBV = 'V', LDV .GE. N.<br>
*>          If JOBV = 'A', LDV .GE. MV.<br>
*> \endverbatim<br>
*><br>
*> \param[in] EPS<br>
*> \verbatim<br>
*>          EPS is DOUBLE PRECISION<br>
*>          EPS = DLAMCH('Epsilon')<br>
*> \endverbatim<br>
*><br>
*> \param[in] SFMIN<br>
*> \verbatim<br>
*>          SFMIN is DOUBLE PRECISION<br>
*>          SFMIN = DLAMCH('Safe Minimum')<br>
*> \endverbatim<br>
*><br>
*> \param[in] TOL<br>
*> \verbatim<br>
*>          TOL is DOUBLE PRECISION<br>
*>          TOL is the threshold for Jacobi rotations. For a pair<br>
*>          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is<br>
*>          applied only if DABS(COS(angle(A(:,p),A(:,q)))) .GT. TOL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NSWEEP<br>
*> \verbatim<br>
*>          NSWEEP is INTEGER<br>
*>          NSWEEP is the number of sweeps of Jacobi rotations to be<br>
*>          performed.<br>
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
*>          LWORK is the dimension of WORK. LWORK .GE. M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0 : successful exit.<br>
*>          < 0 : if INFO = -i, then the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> DGSVJ0 is used just to enable DGESVJ to call a simplified version of<br>
*> itself to work on a submatrix of the original matrix.<br>
*><br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)<br>
*><br>
*> \par Bugs, Examples and Comments:<br>
*  =================================<br>
*><br>
*> Please report all bugs and send interesting test examples and comments to<br>
*> drmac@math.hr. Thank you.<br>
*<br>
*  =====================================================================<br>
*/
	public void dgsvj0_(char[] JOBV,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] D,double[] SVA,INTEGER MV,double[] V,INTEGER LDV,DOUBLE EPS,DOUBLE SFMIN,DOUBLE TOL,INTEGER NSWEEP,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGSVJ1 pre-processor for the routine dgesvj, applies Jacobi rotations targeting only particular pivots.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGSVJ1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgsvj1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgsvj1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgsvj1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV,<br>
*                          EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION   EPS, SFMIN, TOL<br>
*       INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP<br>
*       CHARACTER*1        JOBV<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * ), D( N ), SVA( N ), V( LDV, * ),<br>
*      $                   WORK( LWORK )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGSVJ1 is called from DGESVJ as a pre-processor and that is its main<br>
*> purpose. It applies Jacobi rotations in the same way as DGESVJ does, but<br>
*> it targets only particular pivots and it does not check convergence<br>
*> (stopping criterion). Few tunning parameters (marked by [TP]) are<br>
*> available for the implementer.<br>
*><br>
*> Further Details<br>
*> ~~~~~~~~~~~~~~~<br>
*> DGSVJ1 applies few sweeps of Jacobi rotations in the column space of<br>
*> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)<br>
*> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The<br>
*> block-entries (tiles) of the (1,2) off-diagonal block are marked by the<br>
*> [x]'s in the following scheme:<br>
*><br>
*>    | *  *  * [x] [x] [x]|<br>
*>    | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.<br>
*>    | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.<br>
*>    |[x] [x] [x] *  *  * |<br>
*>    |[x] [x] [x] *  *  * |<br>
*>    |[x] [x] [x] *  *  * |<br>
*><br>
*> In terms of the columns of A, the first N1 columns are rotated 'against'<br>
*> the remaining N-N1 columns, trying to increase the angle between the<br>
*> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is<br>
*> tiled using quadratic tiles of side KBL. Here, KBL is a tunning parmeter.<br>
*> The number of sweeps is given in NSWEEP and the orthogonality threshold<br>
*> is given in TOL.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBV<br>
*> \verbatim<br>
*>          JOBV is CHARACTER*1<br>
*>          Specifies whether the output from this procedure is used<br>
*>          to compute the matrix V:<br>
*>          = 'V': the product of the Jacobi rotations is accumulated<br>
*>                 by postmulyiplying the N-by-N array V.<br>
*>                (See the description of V.)<br>
*>          = 'A': the product of the Jacobi rotations is accumulated<br>
*>                 by postmulyiplying the MV-by-N array V.<br>
*>                (See the descriptions of MV and V.)<br>
*>          = 'N': the Jacobi rotations are not accumulated.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the input matrix A.  M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the input matrix A.<br>
*>          M >= N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N1<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*>          N1 specifies the 2 x 2 block partition, the first N1 columns are<br>
*>          rotated 'against' the remaining N-N1 columns of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          On entry, M-by-N matrix A, such that A*diag(D) represents<br>
*>          the input matrix.<br>
*>          On exit,<br>
*>          A_onexit * D_onexit represents the input matrix A*diag(D)<br>
*>          post-multiplied by a sequence of Jacobi rotations, where the<br>
*>          rotation threshold and the total number of sweeps are given in<br>
*>          TOL and NSWEEP, respectively.<br>
*>          (See the descriptions of N1, D, TOL and NSWEEP.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The array D accumulates the scaling factors from the fast scaled<br>
*>          Jacobi rotations.<br>
*>          On entry, A*diag(D) represents the input matrix.<br>
*>          On exit, A_onexit*diag(D_onexit) represents the input matrix<br>
*>          post-multiplied by a sequence of Jacobi rotations, where the<br>
*>          rotation threshold and the total number of sweeps are given in<br>
*>          TOL and NSWEEP, respectively.<br>
*>          (See the descriptions of N1, A, TOL and NSWEEP.)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SVA<br>
*> \verbatim<br>
*>          SVA is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, SVA contains the Euclidean norms of the columns of<br>
*>          the matrix A*diag(D).<br>
*>          On exit, SVA contains the Euclidean norms of the columns of<br>
*>          the matrix onexit*diag(D_onexit).<br>
*> \endverbatim<br>
*><br>
*> \param[in] MV<br>
*> \verbatim<br>
*>          MV is INTEGER<br>
*>          If JOBV .EQ. 'A', then MV rows of V are post-multipled by a<br>
*>                           sequence of Jacobi rotations.<br>
*>          If JOBV = 'N',   then MV is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] V<br>
*> \verbatim<br>
*>          V is DOUBLE PRECISION array, dimension (LDV,N)<br>
*>          If JOBV .EQ. 'V' then N rows of V are post-multipled by a<br>
*>                           sequence of Jacobi rotations.<br>
*>          If JOBV .EQ. 'A' then MV rows of V are post-multipled by a<br>
*>                           sequence of Jacobi rotations.<br>
*>          If JOBV = 'N',   then V is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V,  LDV >= 1.<br>
*>          If JOBV = 'V', LDV .GE. N.<br>
*>          If JOBV = 'A', LDV .GE. MV.<br>
*> \endverbatim<br>
*><br>
*> \param[in] EPS<br>
*> \verbatim<br>
*>          EPS is DOUBLE PRECISION<br>
*>          EPS = DLAMCH('Epsilon')<br>
*> \endverbatim<br>
*><br>
*> \param[in] SFMIN<br>
*> \verbatim<br>
*>          SFMIN is DOUBLE PRECISION<br>
*>          SFMIN = DLAMCH('Safe Minimum')<br>
*> \endverbatim<br>
*><br>
*> \param[in] TOL<br>
*> \verbatim<br>
*>          TOL is DOUBLE PRECISION<br>
*>          TOL is the threshold for Jacobi rotations. For a pair<br>
*>          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is<br>
*>          applied only if DABS(COS(angle(A(:,p),A(:,q)))) .GT. TOL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NSWEEP<br>
*> \verbatim<br>
*>          NSWEEP is INTEGER<br>
*>          NSWEEP is the number of sweeps of Jacobi rotations to be<br>
*>          performed.<br>
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
*>          LWORK is the dimension of WORK. LWORK .GE. M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0 : successful exit.<br>
*>          < 0 : if INFO = -i, then the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)<br>
*<br>
*  =====================================================================<br>
*/
	public void dgsvj1_(char[] JOBV,INTEGER M,INTEGER N,INTEGER N1,double[] A,INTEGER LDA,double[] D,double[] SVA,INTEGER MV,double[] V,INTEGER LDV,DOUBLE EPS,DOUBLE SFMIN,DOUBLE TOL,INTEGER NSWEEP,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DGTCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGTCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND,<br>
*                          WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          NORM<br>
*       INTEGER            INFO, N<br>
*       DOUBLE PRECISION   ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   D( * ), DL( * ), DU( * ), DU2( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGTCON estimates the reciprocal of the condition number of a real<br>
*> tridiagonal matrix A using the LU factorization as computed by<br>
*> DGTTRF.<br>
*><br>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the<br>
*> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] NORM<br>
*> \verbatim<br>
*>          NORM is CHARACTER*1<br>
*>          Specifies whether the 1-norm condition number or the<br>
*>          infinity-norm condition number is required:<br>
*>          = '1' or 'O':  1-norm;<br>
*>          = 'I':         Infinity-norm.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DL<br>
*> \verbatim<br>
*>          DL is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) multipliers that define the matrix L from the<br>
*>          LU factorization of A as computed by DGTTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The n diagonal elements of the upper triangular matrix U from<br>
*>          the LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) elements of the first superdiagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU2<br>
*> \verbatim<br>
*>          DU2 is DOUBLE PRECISION array, dimension (N-2)<br>
*>          The (n-2) elements of the second superdiagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices; for 1 <= i <= n, row i of the matrix was<br>
*>          interchanged with row IPIV(i).  IPIV(i) will always be either<br>
*>          i or i+1; IPIV(i) = i indicates a row interchange was not<br>
*>          required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is DOUBLE PRECISION<br>
*>          If NORM = '1' or 'O', the 1-norm of the original matrix A.<br>
*>          If NORM = 'I', the infinity-norm of the original matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an<br>
*>          estimate of the 1-norm of inv(A) computed in this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup doubleGTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgtcon_(CHARACTER NORM,INTEGER N,double[] DL,double[] D,double[] DU,double[] DU2,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGTRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGTRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2,<br>
*                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   B( LDB, * ), BERR( * ), D( * ), DF( * ),<br>
*      $                   DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ),<br>
*      $                   FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGTRFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is tridiagonal, and provides<br>
*> error bounds and backward error estimates for the solution.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DL<br>
*> \verbatim<br>
*>          DL is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) subdiagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The diagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) superdiagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DLF<br>
*> \verbatim<br>
*>          DLF is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) multipliers that define the matrix L from the<br>
*>          LU factorization of A as computed by DGTTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DF<br>
*> \verbatim<br>
*>          DF is DOUBLE PRECISION array, dimension (N)<br>
*>          The n diagonal elements of the upper triangular matrix U from<br>
*>          the LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DUF<br>
*> \verbatim<br>
*>          DUF is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) elements of the first superdiagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU2<br>
*> \verbatim<br>
*>          DU2 is DOUBLE PRECISION array, dimension (N-2)<br>
*>          The (n-2) elements of the second superdiagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices; for 1 <= i <= n, row i of the matrix was<br>
*>          interchanged with row IPIV(i).  IPIV(i) will always be either<br>
*>          i or i+1; IPIV(i) = i indicates a row interchange was not<br>
*>          required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          The right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by DGTTRS.<br>
*>          On exit, the improved solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The estimated forward error bound for each solution vector<br>
*>          X(j) (the j-th column of the solution matrix X).<br>
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)<br>
*>          is an estimated upper bound for the magnitude of the largest<br>
*>          element in (X(j) - XTRUE) divided by the magnitude of the<br>
*>          largest element in X(j).  The estimate is as reliable as<br>
*>          the estimate for RCOND, and is almost always a slight<br>
*>          overestimate of the true error.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  ITMAX is the maximum number of steps of iterative refinement.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgtrfs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] DLF,double[] DF,double[] DUF,double[] DU2,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> DGTSV computes the solution to system of linear equations A * X = B for GT matrices <b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGTSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGTSV  solves the equation<br>
*><br>
*>    A*X = B,<br>
*><br>
*> where A is an n by n tridiagonal matrix, by Gaussian elimination with<br>
*> partial pivoting.<br>
*><br>
*> Note that the equation  A**T*X = B  may be solved by interchanging the<br>
*> order of the arguments DU and DL.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
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
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DL<br>
*> \verbatim<br>
*>          DL is DOUBLE PRECISION array, dimension (N-1)<br>
*>          On entry, DL must contain the (n-1) sub-diagonal elements of<br>
*>          A.<br>
*><br>
*>          On exit, DL is overwritten by the (n-2) elements of the<br>
*>          second super-diagonal of the upper triangular matrix U from<br>
*>          the LU factorization of A, in DL(1), ..., DL(n-2).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, D must contain the diagonal elements of A.<br>
*><br>
*>          On exit, D is overwritten by the n diagonal elements of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DU<br>
*> \verbatim<br>
*>          DU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          On entry, DU must contain the (n-1) super-diagonal elements<br>
*>          of A.<br>
*><br>
*>          On exit, DU is overwritten by the (n-1) elements of the first<br>
*>          super-diagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the N by NRHS matrix of right hand side matrix B.<br>
*>          On exit, if INFO = 0, the N by NRHS solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, U(i,i) is exactly zero, and the solution<br>
*>               has not been computed.  The factorization has not been<br>
*>               completed unless i = N.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGTsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgtsv_(INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> DGTSVX computes the solution to system of linear equations A * X = B for GT matrices <b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGTSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtsvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtsvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtsvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,<br>
*                          DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,<br>
*                          WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          FACT, TRANS<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       DOUBLE PRECISION   B( LDB, * ), BERR( * ), D( * ), DF( * ),<br>
*      $                   DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ),<br>
*      $                   FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGTSVX uses the LU factorization to compute the solution to a real<br>
*> system of linear equations A * X = B or A**T * X = B,<br>
*> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS<br>
*> matrices.<br>
*><br>
*> Error bounds on the solution and a condition estimate are also<br>
*> provided.<br>
*> \endverbatim<br>
*<br>
*> \par Description:<br>
*  =================<br>
*><br>
*> \verbatim<br>
*><br>
*> The following steps are performed:<br>
*><br>
*> 1. If FACT = 'N', the LU decomposition is used to factor the matrix A<br>
*>    as A = L * U, where L is a product of permutation and unit lower<br>
*>    bidiagonal matrices and U is upper triangular with nonzeros in<br>
*>    only the main diagonal and first two superdiagonals.<br>
*><br>
*> 2. If some U(i,i)=0, so that U is exactly singular, then the routine<br>
*>    returns with INFO = i. Otherwise, the factored form of A is used<br>
*>    to estimate the condition number of the matrix A.  If the<br>
*>    reciprocal of the condition number is less than machine precision,<br>
*>    INFO = N+1 is returned as a warning, but the routine still goes on<br>
*>    to solve for X and compute error bounds as described below.<br>
*><br>
*> 3. The system of equations is solved for X using the factored form<br>
*>    of A.<br>
*><br>
*> 4. Iterative refinement is applied to improve the computed solution<br>
*>    matrix and calculate error bounds and backward error estimates<br>
*>    for it.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] FACT<br>
*> \verbatim<br>
*>          FACT is CHARACTER*1<br>
*>          Specifies whether or not the factored form of A has been<br>
*>          supplied on entry.<br>
*>          = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored<br>
*>                  form of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV<br>
*>                  will not be modified.<br>
*>          = 'N':  The matrix will be copied to DLF, DF, and DUF<br>
*>                  and factored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DL<br>
*> \verbatim<br>
*>          DL is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) subdiagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The n diagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) superdiagonal elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DLF<br>
*> \verbatim<br>
*>          DLF is DOUBLE PRECISION array, dimension (N-1)<br>
*>          If FACT = 'F', then DLF is an input argument and on entry<br>
*>          contains the (n-1) multipliers that define the matrix L from<br>
*>          the LU factorization of A as computed by DGTTRF.<br>
*><br>
*>          If FACT = 'N', then DLF is an output argument and on exit<br>
*>          contains the (n-1) multipliers that define the matrix L from<br>
*>          the LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DF<br>
*> \verbatim<br>
*>          DF is DOUBLE PRECISION array, dimension (N)<br>
*>          If FACT = 'F', then DF is an input argument and on entry<br>
*>          contains the n diagonal elements of the upper triangular<br>
*>          matrix U from the LU factorization of A.<br>
*><br>
*>          If FACT = 'N', then DF is an output argument and on exit<br>
*>          contains the n diagonal elements of the upper triangular<br>
*>          matrix U from the LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DUF<br>
*> \verbatim<br>
*>          DUF is DOUBLE PRECISION array, dimension (N-1)<br>
*>          If FACT = 'F', then DUF is an input argument and on entry<br>
*>          contains the (n-1) elements of the first superdiagonal of U.<br>
*><br>
*>          If FACT = 'N', then DUF is an output argument and on exit<br>
*>          contains the (n-1) elements of the first superdiagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DU2<br>
*> \verbatim<br>
*>          DU2 is DOUBLE PRECISION array, dimension (N-2)<br>
*>          If FACT = 'F', then DU2 is an input argument and on entry<br>
*>          contains the (n-2) elements of the second superdiagonal of<br>
*>          U.<br>
*><br>
*>          If FACT = 'N', then DU2 is an output argument and on exit<br>
*>          contains the (n-2) elements of the second superdiagonal of<br>
*>          U.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          If FACT = 'F', then IPIV is an input argument and on entry<br>
*>          contains the pivot indices from the LU factorization of A as<br>
*>          computed by DGTTRF.<br>
*><br>
*>          If FACT = 'N', then IPIV is an output argument and on exit<br>
*>          contains the pivot indices from the LU factorization of A;<br>
*>          row i of the matrix was interchanged with row IPIV(i).<br>
*>          IPIV(i) will always be either i or i+1; IPIV(i) = i indicates<br>
*>          a row interchange was not required.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          The N-by-NRHS right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)<br>
*>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is DOUBLE PRECISION<br>
*>          The estimate of the reciprocal condition number of the matrix<br>
*>          A.  If RCOND is less than the machine precision (in<br>
*>          particular, if RCOND = 0), the matrix is singular to working<br>
*>          precision.  This condition is indicated by a return code of<br>
*>          INFO > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The estimated forward error bound for each solution vector<br>
*>          X(j) (the j-th column of the solution matrix X).<br>
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)<br>
*>          is an estimated upper bound for the magnitude of the largest<br>
*>          element in (X(j) - XTRUE) divided by the magnitude of the<br>
*>          largest element in X(j).  The estimate is as reliable as<br>
*>          the estimate for RCOND, and is almost always a slight<br>
*>          overestimate of the true error.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is DOUBLE PRECISION array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is<br>
*>                <= N:  U(i,i) is exactly zero.  The factorization<br>
*>                       has not been completed unless i = N, but the<br>
*>                       factor U is exactly singular, so the solution<br>
*>                       and error bounds could not be computed.<br>
*>                       RCOND = 0 is returned.<br>
*>                = N+1: U is nonsingular, but RCOND is less than machine<br>
*>                       precision, meaning that the matrix is singular<br>
*>                       to working precision.  Nevertheless, the<br>
*>                       solution and error bounds are computed because<br>
*>                       there are a number of situations where the<br>
*>                       computed solution can be more accurate than the<br>
*>                       value of RCOND would suggest.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGTsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void dgtsvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] DLF,double[] DF,double[] DUF,double[] DU2,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b DGTTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGTTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgttrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgttrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgttrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGTTRF( N, DL, D, DU, DU2, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   D( * ), DL( * ), DU( * ), DU2( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGTTRF computes an LU factorization of a real tridiagonal matrix A<br>
*> using elimination with partial pivoting and row interchanges.<br>
*><br>
*> The factorization has the form<br>
*>    A = L * U<br>
*> where L is a product of permutation and unit lower bidiagonal<br>
*> matrices and U is upper triangular with nonzeros in only the main<br>
*> diagonal and first two superdiagonals.<br>
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
*> \param[in,out] DL<br>
*> \verbatim<br>
*>          DL is DOUBLE PRECISION array, dimension (N-1)<br>
*>          On entry, DL must contain the (n-1) sub-diagonal elements of<br>
*>          A.<br>
*><br>
*>          On exit, DL is overwritten by the (n-1) multipliers that<br>
*>          define the matrix L from the LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, D must contain the diagonal elements of A.<br>
*><br>
*>          On exit, D is overwritten by the n diagonal elements of the<br>
*>          upper triangular matrix U from the LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DU<br>
*> \verbatim<br>
*>          DU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          On entry, DU must contain the (n-1) super-diagonal elements<br>
*>          of A.<br>
*><br>
*>          On exit, DU is overwritten by the (n-1) elements of the first<br>
*>          super-diagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DU2<br>
*> \verbatim<br>
*>          DU2 is DOUBLE PRECISION array, dimension (N-2)<br>
*>          On exit, DU2 is overwritten by the (n-2) elements of the<br>
*>          second super-diagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices; for 1 <= i <= n, row i of the matrix was<br>
*>          interchanged with row IPIV(i).  IPIV(i) will always be either<br>
*>          i or i+1; IPIV(i) = i indicates a row interchange was not<br>
*>          required.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -k, the k-th argument had an illegal value<br>
*>          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization<br>
*>                has been completed, but the factor U is exactly<br>
*>                singular, and division by zero will occur if it is used<br>
*>                to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \ingroup doubleGTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgttrf_(INTEGER N,double[] DL,double[] D,double[] DU,double[] DU2,int[] IPIV,INTEGER INFO);
/**
*> \brief \b DGTTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGTTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgttrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgttrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgttrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGTTRS solves one of the systems of equations<br>
*>    A*X = B  or  A**T*X = B,<br>
*> with a tridiagonal matrix A using the LU factorization computed<br>
*> by DGTTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations.<br>
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T* X = B  (Transpose)<br>
*>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DL<br>
*> \verbatim<br>
*>          DL is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) multipliers that define the matrix L from the<br>
*>          LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The n diagonal elements of the upper triangular matrix U from<br>
*>          the LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) elements of the first super-diagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU2<br>
*> \verbatim<br>
*>          DU2 is DOUBLE PRECISION array, dimension (N-2)<br>
*>          The (n-2) elements of the second super-diagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices; for 1 <= i <= n, row i of the matrix was<br>
*>          interchanged with row IPIV(i).  IPIV(i) will always be either<br>
*>          i or i+1; IPIV(i) = i indicates a row interchange was not<br>
*>          required.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the matrix of right hand side vectors B.<br>
*>          On exit, B is overwritten by the solution vectors X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup doubleGTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgttrs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] DU2,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b DGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DGTTS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtts2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtts2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtts2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            ITRANS, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGTTS2 solves one of the systems of equations<br>
*>    A*X = B  or  A**T*X = B,<br>
*> with a tridiagonal matrix A using the LU factorization computed<br>
*> by DGTTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITRANS<br>
*> \verbatim<br>
*>          ITRANS is INTEGER<br>
*>          Specifies the form of the system of equations.<br>
*>          = 0:  A * X = B  (No transpose)<br>
*>          = 1:  A**T* X = B  (Transpose)<br>
*>          = 2:  A**T* X = B  (Conjugate transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DL<br>
*> \verbatim<br>
*>          DL is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) multipliers that define the matrix L from the<br>
*>          LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The n diagonal elements of the upper triangular matrix U from<br>
*>          the LU factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU<br>
*> \verbatim<br>
*>          DU is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) elements of the first super-diagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DU2<br>
*> \verbatim<br>
*>          DU2 is DOUBLE PRECISION array, dimension (N-2)<br>
*>          The (n-2) elements of the second super-diagonal of U.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          The pivot indices; for 1 <= i <= n, row i of the matrix was<br>
*>          interchanged with row IPIV(i).  IPIV(i) will always be either<br>
*>          i or i+1; IPIV(i) = i indicates a row interchange was not<br>
*>          required.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)<br>
*>          On entry, the matrix of right hand side vectors B.<br>
*>          On exit, B is overwritten by the solution vectors X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
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
*> \ingroup doubleGTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void dgtts2_(INTEGER ITRANS,INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] DU2,int[] IPIV,double[] B,INTEGER LDB);

}