package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackCP extends Library
{

	public static LapackCP instance = (LapackCP) Native.loadLibrary("liblapack",LapackCP.class);

/**
*> \brief \b CPBCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,<br>
*                          RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
*       REAL               ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               RWORK( * )<br>
*       COMPLEX            AB( LDAB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPBCON estimates the reciprocal of the condition number (in the<br>
*> 1-norm) of a complex Hermitian positive definite band matrix using<br>
*> the Cholesky factorization A = U**H*U or A = L*L**H computed by<br>
*> CPBTRF.<br>
*><br>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the<br>
*> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangular factor stored in AB;<br>
*>          = 'L':  Lower triangular factor stored in AB.<br>
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
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H of the band matrix A, stored in the<br>
*>          first KD+1 rows of the array.  The j-th column of U or L is<br>
*>          stored in the j-th column of the array AB as follows:<br>
*>          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is REAL<br>
*>          The 1-norm (or infinity-norm) of the Hermitian band matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an<br>
*>          estimate of the 1-norm of inv(A) computed in this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpbcon_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,REAL ANORM,REAL RCOND,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPBEQU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBEQU + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbequ.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbequ.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbequ.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
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
*> CPBEQU computes row and column scalings intended to equilibrate a<br>
*> Hermitian positive definite band matrix A and reduce its condition<br>
*> number (with respect to the two-norm).  S contains the scale factors,<br>
*> S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with<br>
*> elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This<br>
*> choice of S puts the condition number of B within a factor N of the<br>
*> smallest possible condition number over all possible diagonal<br>
*> scalings.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangular of A is stored;<br>
*>          = 'L':  Lower triangular of A is stored.<br>
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
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The upper or lower triangle of the Hermitian band matrix A,<br>
*>          stored in the first KD+1 rows of the array.  The j-th column<br>
*>          of A is stored in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array A.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          If INFO = 0, S contains the scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          If INFO = 0, S contains the ratio of the smallest S(i) to<br>
*>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too<br>
*>          large nor too small, it is not worth scaling by S.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
*>          Absolute value of largest matrix element.  If AMAX is very<br>
*>          close to overflow or very close to underflow, the matrix<br>
*>          should be scaled.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpbequ_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
/**
*> \brief \b CPBRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B,<br>
*                          LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPBRFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is Hermitian positive definite<br>
*> and banded, and provides error bounds and backward error estimates<br>
*> for the solution.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
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
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The upper or lower triangle of the Hermitian band matrix A,<br>
*>          stored in the first KD+1 rows of the array.  The j-th column<br>
*>          of A is stored in the j-th column of the array AB as follows:<br>
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
*> \param[in] AFB<br>
*> \verbatim<br>
*>          AFB is COMPLEX array, dimension (LDAFB,N)<br>
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H of the band matrix A as computed by<br>
*>          CPBTRF, in the same storage format as A (see AB).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>          The leading dimension of the array AFB.  LDAFB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by CPBTRS.<br>
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
*>          FERR is REAL array, dimension (NRHS)<br>
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
*>          BERR is REAL array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpbrfs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPBSTF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBSTF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbstf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbstf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbstf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBSTF( UPLO, N, KD, AB, LDAB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPBSTF computes a split Cholesky factorization of a complex<br>
*> Hermitian positive definite band matrix A.<br>
*><br>
*> This routine is designed to be used in conjunction with CHBGST.<br>
*><br>
*> The factorization has the form  A = S**H*S  where S is a band matrix<br>
*> of the same bandwidth as A and the following structure:<br>
*><br>
*>   S = ( U    )<br>
*>       ( M  L )<br>
*><br>
*> where U is upper triangular of order m = (n+kd)/2, and L is lower<br>
*> triangular of order n-m.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first kd+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*><br>
*>          On exit, if INFO = 0, the factor S from the split Cholesky<br>
*>          factorization A = S**H*S. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, the factorization could not be completed,<br>
*>               because the updated element a(i,i) was negative; the<br>
*>               matrix A is not positive definite.<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The band storage scheme is illustrated by the following example, when<br>
*>  N = 7, KD = 2:<br>
*><br>
*>  S = ( s11  s12  s13                     )<br>
*>      (      s22  s23  s24                )<br>
*>      (           s33  s34                )<br>
*>      (                s44                )<br>
*>      (           s53  s54  s55           )<br>
*>      (                s64  s65  s66      )<br>
*>      (                     s75  s76  s77 )<br>
*><br>
*>  If UPLO = 'U', the array AB holds:<br>
*><br>
*>  on entry:                          on exit:<br>
*><br>
*>   *    *   a13  a24  a35  a46  a57   *    *   s13  s24  s53**H s64**H s75**H<br>
*>   *   a12  a23  a34  a45  a56  a67   *   s12  s23  s34  s54**H s65**H s76**H<br>
*>  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55    s66    s77<br>
*><br>
*>  If UPLO = 'L', the array AB holds:<br>
*><br>
*>  on entry:                          on exit:<br>
*><br>
*>  a11  a22  a33  a44  a55  a66  a77  s11    s22    s33    s44  s55  s66  s77<br>
*>  a21  a32  a43  a54  a65  a76   *   s12**H s23**H s34**H s54  s65  s76   *<br>
*>  a31  a42  a53  a64  a64   *    *   s13**H s24**H s53    s64  s75   *    *<br>
*><br>
*>  Array elements marked * are not used by the routine; s12**H denotes<br>
*>  conjg(s12); the diagonal elements of S are real.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cpbstf_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
/**
*> \brief <b> CPBSV computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AB( LDAB, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPBSV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian positive definite band matrix and X<br>
*> and B are N-by-NRHS matrices.<br>
*><br>
*> The Cholesky decomposition is used to factor A as<br>
*>    A = U**H * U,  if UPLO = 'U', or<br>
*>    A = L * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular band matrix, and L is a lower<br>
*> triangular band matrix, with the same number of superdiagonals or<br>
*> subdiagonals as A.  The factored form of A is then used to solve the<br>
*> system of equations A * X = B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of linear equations, i.e., the order of the<br>
*>          matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KD<br>
*> \verbatim<br>
*>          KD is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
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
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).<br>
*>          See below for further details.<br>
*><br>
*>          On exit, if INFO = 0, the triangular factor U or L from the<br>
*>          Cholesky factorization A = U**H*U or A = L*L**H of the band<br>
*>          matrix A, in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          > 0:  if INFO = i, the leading minor of order i of A is not<br>
*>                positive definite, so the factorization could not be<br>
*>                completed, and the solution has not been computed.<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERsolve<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The band storage scheme is illustrated by the following example, when<br>
*>  N = 6, KD = 2, and UPLO = 'U':<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46<br>
*>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56<br>
*>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66<br>
*><br>
*>  Similarly, if UPLO = 'L' the format of A is as follows:<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66<br>
*>     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *<br>
*>     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *<br>
*><br>
*>  Array elements marked * are not used by the routine.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cpbsv_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> CPBSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbsvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbsvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbsvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB,<br>
*                          EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR,<br>
*                          WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               BERR( * ), FERR( * ), RWORK( * ), S( * )<br>
*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPBSVX uses the Cholesky factorization A = U**H*U or A = L*L**H to<br>
*> compute the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian positive definite band matrix and X<br>
*> and B are N-by-NRHS matrices.<br>
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
*>       diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B<br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.<br>
*><br>
*> 2. If FACT = 'N' or 'E', the Cholesky decomposition is used to<br>
*>    factor the matrix A (after equilibration if FACT = 'E') as<br>
*>       A = U**H * U,  if UPLO = 'U', or<br>
*>       A = L * L**H,  if UPLO = 'L',<br>
*>    where U is an upper triangular band matrix, and L is a lower<br>
*>    triangular band matrix.<br>
*><br>
*> 3. If the leading i-by-i principal minor is not positive definite,<br>
*>    then the routine returns with INFO = i. Otherwise, the factored<br>
*>    form of A is used to estimate the condition number of the matrix<br>
*>    A.  If the reciprocal of the condition number is less than machine<br>
*>    precision, INFO = N+1 is returned as a warning, but the routine<br>
*>    still goes on to solve for X and compute error bounds as<br>
*>    described below.<br>
*><br>
*> 4. The system of equations is solved for X using the factored form<br>
*>    of A.<br>
*><br>
*> 5. Iterative refinement is applied to improve the computed solution<br>
*>    matrix and calculate error bounds and backward error estimates<br>
*>    for it.<br>
*><br>
*> 6. If equilibration was used, the matrix X is premultiplied by<br>
*>    diag(S) so that it solves the original system before<br>
*>    equilibration.<br>
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
*>          = 'F':  On entry, AFB contains the factored form of A.<br>
*>                  If EQUED = 'Y', the matrix A has been equilibrated<br>
*>                  with scaling factors given by S.  AB and AFB will not<br>
*>                  be modified.<br>
*>          = 'N':  The matrix A will be copied to AFB and factored.<br>
*>          = 'E':  The matrix A will be equilibrated if necessary, then<br>
*>                  copied to AFB and factored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of linear equations, i.e., the order of the<br>
*>          matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KD<br>
*> \verbatim<br>
*>          KD is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right-hand sides, i.e., the number of columns<br>
*>          of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first KD+1 rows of the array, except<br>
*>          if FACT = 'F' and EQUED = 'Y', then A must contain the<br>
*>          equilibrated matrix diag(S)*A*diag(S).  The j-th column of A<br>
*>          is stored in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).<br>
*>          See below for further details.<br>
*><br>
*>          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by<br>
*>          diag(S)*A*diag(S).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array A.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AFB<br>
*> \verbatim<br>
*>          AFB is COMPLEX array, dimension (LDAFB,N)<br>
*>          If FACT = 'F', then AFB is an input argument and on entry<br>
*>          contains the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H of the band matrix<br>
*>          A, in the same storage format as A (see AB).  If EQUED = 'Y',<br>
*>          then AFB is the factored form of the equilibrated matrix A.<br>
*><br>
*>          If FACT = 'N', then AFB is an output argument and on exit<br>
*>          returns the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H.<br>
*><br>
*>          If FACT = 'E', then AFB is an output argument and on exit<br>
*>          returns the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H of the equilibrated<br>
*>          matrix A (see the description of A for the form of the<br>
*>          equilibrated matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAFB<br>
*> \verbatim<br>
*>          LDAFB is INTEGER<br>
*>          The leading dimension of the array AFB.  LDAFB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies the form of equilibration that was done.<br>
*>          = 'N':  No equilibration (always true if FACT = 'N').<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*>          EQUED is an input argument if FACT = 'F'; otherwise, it is an<br>
*>          output argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A; not accessed if EQUED = 'N'.  S is<br>
*>          an input argument if FACT = 'F'; otherwise, S is an output<br>
*>          argument.  If FACT = 'F' and EQUED = 'Y', each element of S<br>
*>          must be positive.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>          On entry, the N-by-NRHS right hand side matrix B.<br>
*>          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',<br>
*>          B is overwritten by diag(S) * B.<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to<br>
*>          the original system of equations.  Note that if EQUED = 'Y',<br>
*>          A and B are modified on exit, and the solution to the<br>
*>          equilibrated system is inv(diag(S))*X.<br>
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
*>          RCOND is REAL<br>
*>          The estimate of the reciprocal condition number of the matrix<br>
*>          A after equilibration (if done).  If RCOND is less than the<br>
*>          machine precision (in particular, if RCOND = 0), the matrix<br>
*>          is singular to working precision.  This condition is<br>
*>          indicated by a return code of INFO > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is REAL array, dimension (NRHS)<br>
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
*>          BERR is REAL array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, and i is<br>
*>                <= N:  the leading minor of order i of A is<br>
*>                       not positive definite, so the factorization<br>
*>                       could not be completed, and the solution has not<br>
*>                       been computed. RCOND = 0 is returned.<br>
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
*> \ingroup complexOTHERsolve<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The band storage scheme is illustrated by the following example, when<br>
*>  N = 6, KD = 2, and UPLO = 'U':<br>
*><br>
*>  Two-dimensional storage of the Hermitian matrix A:<br>
*><br>
*>     a11  a12  a13<br>
*>          a22  a23  a24<br>
*>               a33  a34  a35<br>
*>                    a44  a45  a46<br>
*>                         a55  a56<br>
*>     (aij=conjg(aji))         a66<br>
*><br>
*>  Band storage of the upper triangle of A:<br>
*><br>
*>      *    *   a13  a24  a35  a46<br>
*>      *   a12  a23  a34  a45  a56<br>
*>     a11  a22  a33  a44  a55  a66<br>
*><br>
*>  Similarly, if UPLO = 'L' the format of A is as follows:<br>
*><br>
*>     a11  a22  a33  a44  a55  a66<br>
*>     a21  a32  a43  a54  a65   *<br>
*>     a31  a42  a53  a64   *    *<br>
*><br>
*>  Array elements marked * are not used by the routine.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cpbsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPBTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite band matrix (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBTF2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbtf2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbtf2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbtf2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBTF2( UPLO, N, KD, AB, LDAB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPBTF2 computes the Cholesky factorization of a complex Hermitian<br>
*> positive definite band matrix A.<br>
*><br>
*> The factorization has the form<br>
*>    A = U**H * U ,  if UPLO = 'U', or<br>
*>    A = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix, U**H is the conjugate transpose<br>
*> of U, and L is lower triangular.<br>
*><br>
*> This is the unblocked version of the algorithm, calling Level 2 BLAS.<br>
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
*>          On entry, the upper or lower triangle of the Hermitian band<br>
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
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -k, the k-th argument had an illegal value<br>
*>          > 0: if INFO = k, the leading minor of order k is not<br>
*>               positive definite, and the factorization could not be<br>
*>               completed.<br>
*> \endverbatim<br>
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
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The band storage scheme is illustrated by the following example, when<br>
*>  N = 6, KD = 2, and UPLO = 'U':<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46<br>
*>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56<br>
*>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66<br>
*><br>
*>  Similarly, if UPLO = 'L' the format of A is as follows:<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66<br>
*>     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *<br>
*>     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *<br>
*><br>
*>  Array elements marked * are not used by the routine.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cpbtf2_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
/**
*> \brief \b CPBTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbtrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbtrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbtrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBTRF( UPLO, N, KD, AB, LDAB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AB( LDAB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPBTRF computes the Cholesky factorization of a complex Hermitian<br>
*> positive definite band matrix A.<br>
*><br>
*> The factorization has the form<br>
*>    A = U**H * U,  if UPLO = 'U', or<br>
*>    A = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is lower triangular.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*><br>
*>          On exit, if INFO = 0, the triangular factor U or L from the<br>
*>          Cholesky factorization A = U**H*U or A = L*L**H of the band<br>
*>          matrix A, in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the leading minor of order i is not<br>
*>                positive definite, and the factorization could not be<br>
*>                completed.<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The band storage scheme is illustrated by the following example, when<br>
*>  N = 6, KD = 2, and UPLO = 'U':<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46<br>
*>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56<br>
*>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66<br>
*><br>
*>  Similarly, if UPLO = 'L' the format of A is as follows:<br>
*><br>
*>  On entry:                       On exit:<br>
*><br>
*>     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66<br>
*>     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *<br>
*>     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *<br>
*><br>
*>  Array elements marked * are not used by the routine.<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989<br>
*<br>
*  =====================================================================<br>
*/
	public void cpbtrf_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
/**
*> \brief \b CPBTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPBTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbtrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbtrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbtrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AB( LDAB, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPBTRS solves a system of linear equations A*X = B with a Hermitian<br>
*> positive definite band matrix A using the Cholesky factorization<br>
*> A = U**H*U or A = L*L**H computed by CPBTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangular factor stored in AB;<br>
*>          = 'L':  Lower triangular factor stored in AB.<br>
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
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
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
*>          AB is COMPLEX array, dimension (LDAB,N)<br>
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H of the band matrix A, stored in the<br>
*>          first KD+1 rows of the array.  The j-th column of U or L is<br>
*>          stored in the j-th column of the array AB as follows:<br>
*>          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpbtrs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b CPFTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPFTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpftrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpftrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpftrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPFTRF( TRANSR, UPLO, N, A, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            N, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( 0: * )<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPFTRF computes the Cholesky factorization of a complex Hermitian<br>
*> positive definite matrix A.<br>
*><br>
*> The factorization has the form<br>
*>    A = U**H * U,  if UPLO = 'U', or<br>
*>    A = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is lower triangular.<br>
*><br>
*> This is the block version of the algorithm, calling Level 3 BLAS.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  The Normal TRANSR of RFP A is stored;<br>
*>          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of RFP A is stored;<br>
*>          = 'L':  Lower triangle of RFP A is stored.<br>
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
*>          A is COMPLEX array, dimension ( N*(N+1)/2 );<br>
*>          On entry, the Hermitian matrix A in RFP format. RFP format is<br>
*>          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'<br>
*>          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is<br>
*>          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is<br>
*>          the Conjugate-transpose of RFP A as defined when<br>
*>          TRANSR = 'N'. The contents of RFP A are defined by UPLO as<br>
*>          follows: If UPLO = 'U' the RFP A contains the nt elements of<br>
*>          upper packed A. If UPLO = 'L' the RFP A contains the elements<br>
*>          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR =<br>
*>          'C'. When TRANSR is 'N' the LDA is N+1 when N is even and N<br>
*>          is odd. See the Note below for more details.<br>
*><br>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky<br>
*>          factorization RFP A = U**H*U or RFP A = L*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the leading minor of order i is not<br>
*>                positive definite, and the factorization could not be<br>
*>                completed.<br>
*><br>
*>  Further Notes on RFP Format:<br>
*>  ============================<br>
*><br>
*>  We first consider Standard Packed Format when N is even.<br>
*>  We give an example where N = 6.<br>
*><br>
*>     AP is Upper             AP is Lower<br>
*><br>
*>   00 01 02 03 04 05       00<br>
*>      11 12 13 14 15       10 11<br>
*>         22 23 24 25       20 21 22<br>
*>            33 34 35       30 31 32 33<br>
*>               44 45       40 41 42 43 44<br>
*>                  55       50 51 52 53 54 55<br>
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
*>           RFP A                   RFP A<br>
*><br>
*>     -- -- -- --                -- -- -- -- -- --<br>
*>     03 13 23 33 00 01 02    33 00 10 20 30 40 50<br>
*>     -- -- -- -- --                -- -- -- -- --<br>
*>     04 14 24 34 44 11 12    43 44 11 21 31 41 51<br>
*>     -- -- -- -- -- --                -- -- -- --<br>
*>     05 15 25 35 45 55 22    53 54 55 22 32 42 52<br>
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
*>           RFP A                   RFP A<br>
*><br>
*>     -- -- --                   -- -- -- -- -- --<br>
*>     02 12 22 00 01             00 10 20 30 40 50<br>
*>     -- -- -- --                   -- -- -- -- --<br>
*>     03 13 23 33 11             33 11 21 31 41 51<br>
*>     -- -- -- -- --                   -- -- -- --<br>
*>     04 14 24 34 44             43 44 22 32 42 52<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpftrf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,INTEGER INFO);
/**
*> \brief \b CPFTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPFTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpftri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpftri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpftri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPFTRI( TRANSR, UPLO, N, A, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPFTRI computes the inverse of a complex Hermitian positive definite<br>
*> matrix A using the Cholesky factorization A = U**H*U or A = L*L**H<br>
*> computed by CPFTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  The Normal TRANSR of RFP A is stored;<br>
*>          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          A is COMPLEX array, dimension ( N*(N+1)/2 );<br>
*>          On entry, the Hermitian matrix A in RFP format. RFP format is<br>
*>          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'<br>
*>          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is<br>
*>          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is<br>
*>          the Conjugate-transpose of RFP A as defined when<br>
*>          TRANSR = 'N'. The contents of RFP A are defined by UPLO as<br>
*>          follows: If UPLO = 'U' the RFP A contains the nt elements of<br>
*>          upper packed A. If UPLO = 'L' the RFP A contains the elements<br>
*>          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR =<br>
*>          'C'. When TRANSR is 'N' the LDA is N+1 when N is even and N<br>
*>          is odd. See the Note below for more details.<br>
*><br>
*>          On exit, the Hermitian inverse of the original matrix, in the<br>
*>          same storage format.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the (i,i) element of the factor U or L is<br>
*>                zero, and the inverse could not be computed.<br>
*> \endverbatim<br>
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
	public void cpftri_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,INTEGER INFO);
/**
*> \brief \b CPFTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPFTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpftrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpftrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpftrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPFTRS( TRANSR, UPLO, N, NRHS, A, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( 0: * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPFTRS solves a system of linear equations A*X = B with a Hermitian<br>
*> positive definite matrix A using the Cholesky factorization<br>
*> A = U**H*U or A = L*L**H computed by CPFTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  The Normal TRANSR of RFP A is stored;<br>
*>          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of RFP A is stored;<br>
*>          = 'L':  Lower triangle of RFP A is stored.<br>
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
*>          A is COMPLEX array, dimension ( N*(N+1)/2 );<br>
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          of RFP A = U**H*U or RFP A = L*L**H, as computed by CPFTRF.<br>
*>          See note below for more details about RFP A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
	public void cpftrs_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b CPOCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpocon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpocon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpocon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, RWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       REAL               ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               RWORK( * )<br>
*       COMPLEX            A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPOCON estimates the reciprocal of the condition number (in the<br>
*> 1-norm) of a complex Hermitian positive definite matrix using the<br>
*> Cholesky factorization A = U**H*U or A = L*L**H computed by CPOTRF.<br>
*><br>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the<br>
*> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H, as computed by CPOTRF.<br>
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
*>          ANORM is REAL<br>
*>          The 1-norm (or infinity-norm) of the Hermitian matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an<br>
*>          estimate of the 1-norm of inv(A) computed in this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
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
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpocon_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,REAL ANORM,REAL RCOND,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPOEQU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOEQU + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpoequ.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpoequ.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpoequ.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOEQU( N, A, LDA, S, SCOND, AMAX, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, N<br>
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
*> CPOEQU computes row and column scalings intended to equilibrate a<br>
*> Hermitian positive definite matrix A and reduce its condition number<br>
*> (with respect to the two-norm).  S contains the scale factors,<br>
*> S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with<br>
*> elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This<br>
*> choice of S puts the condition number of B within a factor N of the<br>
*> smallest possible condition number over all possible diagonal<br>
*> scalings.<br>
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
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The N-by-N Hermitian positive definite matrix whose scaling<br>
*>          factors are to be computed.  Only the diagonal elements of A<br>
*>          are referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          If INFO = 0, S contains the scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          If INFO = 0, S contains the ratio of the smallest S(i) to<br>
*>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too<br>
*>          large nor too small, it is not worth scaling by S.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
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
*>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.<br>
*> \endverbatim<br>
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
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpoequ_(INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
/**
*> \brief \b CPOEQUB<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOEQUB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpoequb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpoequb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpoequb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, N<br>
*       REAL               AMAX, SCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       REAL               S( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPOEQUB computes row and column scalings intended to equilibrate a<br>
*> symmetric positive definite matrix A and reduce its condition number<br>
*> (with respect to the two-norm).  S contains the scale factors,<br>
*> S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with<br>
*> elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This<br>
*> choice of S puts the condition number of B within a factor N of the<br>
*> smallest possible condition number over all possible diagonal<br>
*> scalings.<br>
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
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The N-by-N symmetric positive definite matrix whose scaling<br>
*>          factors are to be computed.  Only the diagonal elements of A<br>
*>          are referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          If INFO = 0, S contains the scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          If INFO = 0, S contains the ratio of the smallest S(i) to<br>
*>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too<br>
*>          large nor too small, it is not worth scaling by S.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
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
*>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.<br>
*> \endverbatim<br>
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
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpoequb_(INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
/**
*> \brief \b CPORFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPORFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cporfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cporfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cporfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X,<br>
*                          LDX, FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPORFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is Hermitian positive definite,<br>
*> and provides error bounds and backward error estimates for the<br>
*> solution.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The Hermitian matrix A.  If UPLO = 'U', the leading N-by-N<br>
*>          upper triangular part of A contains the upper triangular part<br>
*>          of the matrix A, and the strictly lower triangular part of A<br>
*>          is not referenced.  If UPLO = 'L', the leading N-by-N lower<br>
*>          triangular part of A contains the lower triangular part of<br>
*>          the matrix A, and the strictly upper triangular part of A is<br>
*>          not referenced.<br>
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
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H, as computed by CPOTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>          The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by CPOTRS.<br>
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
*>          FERR is REAL array, dimension (NRHS)<br>
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
*>          BERR is REAL array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
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
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cporfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPORFSX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPORFSX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cporfsx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cporfsx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cporfsx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPORFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, S, B,<br>
*                           LDB, X, LDX, RCOND, BERR, N_ERR_BNDS,<br>
*                           ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS,<br>
*                           WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, EQUED<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   X( LDX, * ), WORK( * )<br>
*       REAL               RWORK( * ), S( * ), PARAMS(*), BERR( * ),<br>
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
*>    CPORFSX improves the computed solution to a system of linear<br>
*>    equations when the coefficient matrix is symmetric positive<br>
*>    definite, and provides error bounds and backward error estimates<br>
*>    for the solution.  In addition to normwise error bound, the code<br>
*>    provides maximum componentwise error bound if possible.  See<br>
*>    comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the<br>
*>    error bounds.<br>
*><br>
*>    The original system of linear equations may have been equilibrated<br>
*>    before calling this routine, as described by arguments EQUED and S<br>
*>    below. In this case, the solution and error bounds returned are<br>
*>    for the original unequilibrated system.<br>
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
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>       = 'U':  Upper triangle of A is stored;<br>
*>       = 'L':  Lower triangle of A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>     Specifies the form of equilibration that was done to A<br>
*>     before calling this routine. This is needed to compute<br>
*>     the solution and error bounds correctly.<br>
*>       = 'N':  No equilibration<br>
*>       = 'Y':  Both row and column equilibration, i.e., A has been<br>
*>               replaced by diag(S) * A * diag(S).<br>
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
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     The symmetric matrix A.  If UPLO = 'U', the leading N-by-N<br>
*>     upper triangular part of A contains the upper triangular part<br>
*>     of the matrix A, and the strictly lower triangular part of A<br>
*>     is not referenced.  If UPLO = 'L', the leading N-by-N lower<br>
*>     triangular part of A contains the lower triangular part of<br>
*>     the matrix A, and the strictly upper triangular part of A is<br>
*>     not referenced.<br>
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
*>     A = U**T*U or A = L*L**T, as computed by SPOTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>     The row scale factors for A.  If EQUED = 'Y', A is multiplied on<br>
*>     the left and right by diag(S).  S is an input argument if FACT =<br>
*>     'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED<br>
*>     = 'Y', each element of S must be positive.  If S is output, each<br>
*>     element of S is a power of the radix. If S is input, each element<br>
*>     of S should be a power of the radix to ensure a reliable solution<br>
*>     and error estimates. Scaling by powers of the radix does not cause<br>
*>     rounding errors unless the result underflows or overflows.<br>
*>     Rounding errors during scaling lead to refining with a matrix that<br>
*>     is not equivalent to the input matrix, producing error estimates<br>
*>     that may not be reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>     On entry, the solution matrix X, as computed by SGETRS.<br>
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
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is REAL array, dimension (NRHS)<br>
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
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_COMP<br>
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
*>          PARAMS is REAL array, dimension NPARAMS<br>
*>     Specifies algorithm parameters.  If an entry is .LT. 0.0, then<br>
*>     that entry will be filled with default value used for that<br>
*>     parameter.  Only positions up to NPARAMS are accessed; defaults<br>
*>     are used for higher-numbered parameters.<br>
*><br>
*>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative<br>
*>            refinement or not.<br>
*>         Default: 1.0<br>
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
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (2*N)<br>
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
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cporfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief <b> CPOSV computes the solution to system of linear equations A * X = B for PO matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cposv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cposv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cposv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
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
*> CPOSV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian positive definite matrix and X and B<br>
*> are N-by-NRHS matrices.<br>
*><br>
*> The Cholesky decomposition is used to factor A as<br>
*>    A = U**H* U,  if UPLO = 'U', or<br>
*>    A = L * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and  L is a lower triangular<br>
*> matrix.  The factored form of A is then used to solve the system of<br>
*> equations A * X = B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H.<br>
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
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          > 0:  if INFO = i, the leading minor of order i of A is not<br>
*>                positive definite, so the factorization could not be<br>
*>                completed, and the solution has not been computed.<br>
*> \endverbatim<br>
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
*> \ingroup complexPOsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void cposv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> CPOSVX computes the solution to system of linear equations A * X = B for PO matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cposvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cposvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cposvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED,<br>
*                          S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,<br>
*                          RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               BERR( * ), FERR( * ), RWORK( * ), S( * )<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPOSVX uses the Cholesky factorization A = U**H*U or A = L*L**H to<br>
*> compute the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian positive definite matrix and X and B<br>
*> are N-by-NRHS matrices.<br>
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
*>       diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B<br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.<br>
*><br>
*> 2. If FACT = 'N' or 'E', the Cholesky decomposition is used to<br>
*>    factor the matrix A (after equilibration if FACT = 'E') as<br>
*>       A = U**H* U,  if UPLO = 'U', or<br>
*>       A = L * L**H,  if UPLO = 'L',<br>
*>    where U is an upper triangular matrix and L is a lower triangular<br>
*>    matrix.<br>
*><br>
*> 3. If the leading i-by-i principal minor is not positive definite,<br>
*>    then the routine returns with INFO = i. Otherwise, the factored<br>
*>    form of A is used to estimate the condition number of the matrix<br>
*>    A.  If the reciprocal of the condition number is less than machine<br>
*>    precision, INFO = N+1 is returned as a warning, but the routine<br>
*>    still goes on to solve for X and compute error bounds as<br>
*>    described below.<br>
*><br>
*> 4. The system of equations is solved for X using the factored form<br>
*>    of A.<br>
*><br>
*> 5. Iterative refinement is applied to improve the computed solution<br>
*>    matrix and calculate error bounds and backward error estimates<br>
*>    for it.<br>
*><br>
*> 6. If equilibration was used, the matrix X is premultiplied by<br>
*>    diag(S) so that it solves the original system before<br>
*>    equilibration.<br>
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
*>          = 'F':  On entry, AF contains the factored form of A.<br>
*>                  If EQUED = 'Y', the matrix A has been equilibrated<br>
*>                  with scaling factors given by S.  A and AF will not<br>
*>                  be modified.<br>
*>          = 'N':  The matrix A will be copied to AF and factored.<br>
*>          = 'E':  The matrix A will be equilibrated if necessary, then<br>
*>                  copied to AF and factored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A, except if FACT = 'F' and<br>
*>          EQUED = 'Y', then A must contain the equilibrated matrix<br>
*>          diag(S)*A*diag(S).  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.  A is not modified if<br>
*>          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.<br>
*><br>
*>          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by<br>
*>          diag(S)*A*diag(S).<br>
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
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>          If FACT = 'F', then AF is an input argument and on entry<br>
*>          contains the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H, in the same storage<br>
*>          format as A.  If EQUED .ne. 'N', then AF is the factored form<br>
*>          of the equilibrated matrix diag(S)*A*diag(S).<br>
*><br>
*>          If FACT = 'N', then AF is an output argument and on exit<br>
*>          returns the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H of the original<br>
*>          matrix A.<br>
*><br>
*>          If FACT = 'E', then AF is an output argument and on exit<br>
*>          returns the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H of the equilibrated<br>
*>          matrix A (see the description of A for the form of the<br>
*>          equilibrated matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>          The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies the form of equilibration that was done.<br>
*>          = 'N':  No equilibration (always true if FACT = 'N').<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*>          EQUED is an input argument if FACT = 'F'; otherwise, it is an<br>
*>          output argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A; not accessed if EQUED = 'N'.  S is<br>
*>          an input argument if FACT = 'F'; otherwise, S is an output<br>
*>          argument.  If FACT = 'F' and EQUED = 'Y', each element of S<br>
*>          must be positive.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>          On entry, the N-by-NRHS righthand side matrix B.<br>
*>          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',<br>
*>          B is overwritten by diag(S) * B.<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to<br>
*>          the original system of equations.  Note that if EQUED = 'Y',<br>
*>          A and B are modified on exit, and the solution to the<br>
*>          equilibrated system is inv(diag(S))*X.<br>
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
*>          RCOND is REAL<br>
*>          The estimate of the reciprocal condition number of the matrix<br>
*>          A after equilibration (if done).  If RCOND is less than the<br>
*>          machine precision (in particular, if RCOND = 0), the matrix<br>
*>          is singular to working precision.  This condition is<br>
*>          indicated by a return code of INFO > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is REAL array, dimension (NRHS)<br>
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
*>          BERR is REAL array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, and i is<br>
*>                <= N:  the leading minor of order i of A is<br>
*>                       not positive definite, so the factorization<br>
*>                       could not be completed, and the solution has not<br>
*>                       been computed. RCOND = 0 is returned.<br>
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
*> \ingroup complexPOsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void cposvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief <b> CPOSVXX computes the solution to system of linear equations A * X = B for PO matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOSVXX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cposvxx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cposvxx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cposvxx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED,<br>
*                           S, B, LDB, X, LDX, RCOND, RPVGRW, BERR,<br>
*                           N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP,<br>
*                           NPARAMS, PARAMS, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS<br>
*       REAL               RCOND, RPVGRW<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       REAL               S( * ), PARAMS( * ), BERR( * ), RWORK( * ),<br>
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
*>    CPOSVXX uses the Cholesky factorization A = U**T*U or A = L*L**T<br>
*>    to compute the solution to a complex system of linear equations<br>
*>    A * X = B, where A is an N-by-N symmetric positive definite matrix<br>
*>    and X and B are N-by-NRHS matrices.<br>
*><br>
*>    If requested, both normwise and maximum componentwise error bounds<br>
*>    are returned. CPOSVXX will return a solution with a tiny<br>
*>    guaranteed error (O(eps) where eps is the working machine<br>
*>    precision) unless the matrix is very ill-conditioned, in which<br>
*>    case a warning is returned. Relevant condition numbers also are<br>
*>    calculated and returned.<br>
*><br>
*>    CPOSVXX accepts user-provided factorizations and equilibration<br>
*>    factors; see the definitions of the FACT and EQUED options.<br>
*>    Solving with refinement and using a factorization from a previous<br>
*>    CPOSVXX call will also produce a solution with either O(eps)<br>
*>    errors or warnings, but we cannot make that claim for general<br>
*>    user-provided factorizations and equilibration factors if they<br>
*>    differ from what CPOSVXX would itself produce.<br>
*> \endverbatim<br>
*<br>
*> \par Description:<br>
*  =================<br>
*><br>
*> \verbatim<br>
*><br>
*>    The following steps are performed:<br>
*><br>
*>    1. If FACT = 'E', real scaling factors are computed to equilibrate<br>
*>    the system:<br>
*><br>
*>      diag(S)*A*diag(S)     *inv(diag(S))*X = diag(S)*B<br>
*><br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.<br>
*><br>
*>    2. If FACT = 'N' or 'E', the Cholesky decomposition is used to<br>
*>    factor the matrix A (after equilibration if FACT = 'E') as<br>
*>       A = U**T* U,  if UPLO = 'U', or<br>
*>       A = L * L**T,  if UPLO = 'L',<br>
*>    where U is an upper triangular matrix and L is a lower triangular<br>
*>    matrix.<br>
*><br>
*>    3. If the leading i-by-i principal minor is not positive definite,<br>
*>    then the routine returns with INFO = i. Otherwise, the factored<br>
*>    form of A is used to estimate the condition number of the matrix<br>
*>    A (see argument RCOND).  If the reciprocal of the condition number<br>
*>    is less than machine precision, the routine still goes on to solve<br>
*>    for X and compute error bounds as described below.<br>
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
*>    diag(S) so that it solves the original system before<br>
*>    equilibration.<br>
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
*>       = 'F':  On entry, AF contains the factored form of A.<br>
*>               If EQUED is not 'N', the matrix A has been<br>
*>               equilibrated with scaling factors given by S.<br>
*>               A and AF are not modified.<br>
*>       = 'N':  The matrix A will be copied to AF and factored.<br>
*>       = 'E':  The matrix A will be equilibrated if necessary, then<br>
*>               copied to AF and factored.<br>
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
*>     The number of right hand sides, i.e., the number of columns<br>
*>     of the matrices B and X.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>     On entry, the symmetric matrix A, except if FACT = 'F' and EQUED =<br>
*>     'Y', then A must contain the equilibrated matrix<br>
*>     diag(S)*A*diag(S).  If UPLO = 'U', the leading N-by-N upper<br>
*>     triangular part of A contains the upper triangular part of the<br>
*>     matrix A, and the strictly lower triangular part of A is not<br>
*>     referenced.  If UPLO = 'L', the leading N-by-N lower triangular<br>
*>     part of A contains the lower triangular part of the matrix A, and<br>
*>     the strictly upper triangular part of A is not referenced.  A is<br>
*>     not modified if FACT = 'F' or 'N', or if FACT = 'E' and EQUED =<br>
*>     'N' on exit.<br>
*><br>
*>     On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by<br>
*>     diag(S)*A*diag(S).<br>
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
*>          AF is COMPLEX array, dimension (LDAF,N)<br>
*>     If FACT = 'F', then AF is an input argument and on entry<br>
*>     contains the triangular factor U or L from the Cholesky<br>
*>     factorization A = U**T*U or A = L*L**T, in the same storage<br>
*>     format as A.  If EQUED .ne. 'N', then AF is the factored<br>
*>     form of the equilibrated matrix diag(S)*A*diag(S).<br>
*><br>
*>     If FACT = 'N', then AF is an output argument and on exit<br>
*>     returns the triangular factor U or L from the Cholesky<br>
*>     factorization A = U**T*U or A = L*L**T of the original<br>
*>     matrix A.<br>
*><br>
*>     If FACT = 'E', then AF is an output argument and on exit<br>
*>     returns the triangular factor U or L from the Cholesky<br>
*>     factorization A = U**T*U or A = L*L**T of the equilibrated<br>
*>     matrix A (see the description of A for the form of the<br>
*>     equilibrated matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAF<br>
*> \verbatim<br>
*>          LDAF is INTEGER<br>
*>     The leading dimension of the array AF.  LDAF >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>     Specifies the form of equilibration that was done.<br>
*>       = 'N':  No equilibration (always true if FACT = 'N').<br>
*>       = 'Y':  Both row and column equilibration, i.e., A has been<br>
*>               replaced by diag(S) * A * diag(S).<br>
*>     EQUED is an input argument if FACT = 'F'; otherwise, it is an<br>
*>     output argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>     The row scale factors for A.  If EQUED = 'Y', A is multiplied on<br>
*>     the left and right by diag(S).  S is an input argument if FACT =<br>
*>     'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED<br>
*>     = 'Y', each element of S must be positive.  If S is output, each<br>
*>     element of S is a power of the radix. If S is input, each element<br>
*>     of S should be a power of the radix to ensure a reliable solution<br>
*>     and error estimates. Scaling by powers of the radix does not cause<br>
*>     rounding errors unless the result underflows or overflows.<br>
*>     Rounding errors during scaling lead to refining with a matrix that<br>
*>     is not equivalent to the input matrix, producing error estimates<br>
*>     that may not be reliable.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>     On entry, the N-by-NRHS right hand side matrix B.<br>
*>     On exit,<br>
*>     if EQUED = 'N', B is not modified;<br>
*>     if EQUED = 'Y', B is overwritten by diag(S)*B;<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>     If INFO = 0, the N-by-NRHS solution matrix X to the original<br>
*>     system of equations.  Note that A and B are modified on exit if<br>
*>     EQUED .ne. 'N', and the solution to the equilibrated system is<br>
*>     inv(diag(S))*X.<br>
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
*> \param[out] RPVGRW<br>
*> \verbatim<br>
*>          RPVGRW is REAL<br>
*>     Reciprocal pivot growth.  On exit, this contains the reciprocal<br>
*>     pivot growth factor norm(A)/norm(U). The "max absolute element"<br>
*>     norm is used.  If this is much less than 1, then the stability of<br>
*>     the LU factorization of the (equilibrated) matrix A could be poor.<br>
*>     This also means that the solution X, estimated condition numbers,<br>
*>     and error bounds could be unreliable. If factorization fails with<br>
*>     0<INFO<=N, then this contains the reciprocal pivot growth factor<br>
*>     for the leading INFO columns of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is REAL array, dimension (NRHS)<br>
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
*>     See Lapack Working Note 165 for further details and extra<br>
*>     cautions.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ERR_BNDS_COMP<br>
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
*>          PARAMS is REAL array, dimension NPARAMS<br>
*>     Specifies algorithm parameters.  If an entry is .LT. 0.0, then<br>
*>     that entry will be filled with default value used for that<br>
*>     parameter.  Only positions up to NPARAMS are accessed; defaults<br>
*>     are used for higher-numbered parameters.<br>
*><br>
*>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative<br>
*>            refinement or not.<br>
*>         Default: 1.0<br>
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
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (2*N)<br>
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
*> \ingroup complexPOsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void cposvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL CLA_PORPVGRW);
/**
*> \brief \b CPOTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite matrix (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOTF2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpotf2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpotf2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpotf2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOTF2( UPLO, N, A, LDA, INFO )<br>
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
*> CPOTF2 computes the Cholesky factorization of a complex Hermitian<br>
*> positive definite matrix A.<br>
*><br>
*> The factorization has the form<br>
*>    A = U**H * U ,  if UPLO = 'U', or<br>
*>    A = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is lower triangular.<br>
*><br>
*> This is the unblocked version of the algorithm, calling Level 2 BLAS.<br>
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
*>          On exit, if INFO = 0, the factor U or L from the Cholesky<br>
*>          factorization A = U**H *U  or A = L*L**H.<br>
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
*>          > 0: if INFO = k, the leading minor of order k is not<br>
*>               positive definite, and the factorization could not be<br>
*>               completed.<br>
*> \endverbatim<br>
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
	public void cpotf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b CPOTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpotrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpotrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpotrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOTRF( UPLO, N, A, LDA, INFO )<br>
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
*> CPOTRF computes the Cholesky factorization of a complex Hermitian<br>
*> positive definite matrix A.<br>
*><br>
*> The factorization has the form<br>
*>    A = U**H * U,  if UPLO = 'U', or<br>
*>    A = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is lower triangular.<br>
*><br>
*> This is the block version of the algorithm, calling Level 3 BLAS.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H.<br>
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
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the leading minor of order i is not<br>
*>                positive definite, and the factorization could not be<br>
*>                completed.<br>
*> \endverbatim<br>
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
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpotrf_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b CPOTRF2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       RECURSIVE SUBROUTINE CPOTRF2( UPLO, N, A, LDA, INFO )<br>
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
*> CPOTRF2 computes the Cholesky factorization of a real symmetric<br>
*> positive definite matrix A using the recursive algorithm.<br>
*><br>
*> The factorization has the form<br>
*>    A = U**H * U,  if UPLO = 'U', or<br>
*>    A = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is lower triangular.<br>
*><br>
*> This is the recursive version of the algorithm. It divides<br>
*> the matrix into four submatrices:<br>
*><br>
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2<br>
*>    A = [ -----|----- ]  with n1 = n/2<br>
*>        [  A21 | A22  ]       n2 = n-n1<br>
*><br>
*> The subroutine calls itself to factor A11. Update and scale A21<br>
*> or A12, update A22 then calls itself to factor A22.<br>
*> <br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H.<br>
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
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the leading minor of order i is not<br>
*>                positive definite, and the factorization could not be<br>
*>                completed.<br>
*> \endverbatim<br>
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
	public void cpotrf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b CPOTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpotri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpotri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpotri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOTRI( UPLO, N, A, LDA, INFO )<br>
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
*> CPOTRI computes the inverse of a complex Hermitian positive definite<br>
*> matrix A using the Cholesky factorization A = U**H*U or A = L*L**H<br>
*> computed by CPOTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          On entry, the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H, as computed by<br>
*>          CPOTRF.<br>
*>          On exit, the upper or lower triangle of the (Hermitian)<br>
*>          inverse of A, overwriting the input factor U or L.<br>
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
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the (i,i) element of the factor U or L is<br>
*>                zero, and the inverse could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpotri_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b CPOTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPOTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpotrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpotrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpotrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
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
*> CPOTRS solves a system of linear equations A*X = B with a Hermitian<br>
*> positive definite matrix A using the Cholesky factorization <br>
*> A = U**H*U or A = L*L**H computed by CPOTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H, as computed by CPOTRF.<br>
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
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*> \ingroup complexPOcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpotrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b CPPCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPPCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPPCON( UPLO, N, AP, ANORM, RCOND, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       REAL               ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               RWORK( * )<br>
*       COMPLEX            AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPPCON estimates the reciprocal of the condition number (in the <br>
*> 1-norm) of a complex Hermitian positive definite packed matrix using<br>
*> the Cholesky factorization A = U**H*U or A = L*L**H computed by<br>
*> CPPTRF.<br>
*><br>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the<br>
*> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H, packed columnwise in a linear<br>
*>          array.  The j-th column of U or L is stored in the array AP<br>
*>          as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is REAL<br>
*>          The 1-norm (or infinity-norm) of the Hermitian matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an<br>
*>          estimate of the 1-norm of inv(A) computed in this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cppcon_(CHARACTER UPLO,INTEGER N,float[] AP,REAL ANORM,REAL RCOND,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPPEQU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPPEQU + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppequ.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppequ.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppequ.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
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
*> CPPEQU computes row and column scalings intended to equilibrate a<br>
*> Hermitian positive definite matrix A in packed storage and reduce<br>
*> its condition number (with respect to the two-norm).  S contains the<br>
*> scale factors, S(i)=1/sqrt(A(i,i)), chosen so that the scaled matrix<br>
*> B with elements B(i,j)=S(i)*A(i,j)*S(j) has ones on the diagonal.<br>
*> This choice of S puts the condition number of B within a factor N of<br>
*> the smallest possible condition number over all possible diagonal<br>
*> scalings.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          The upper or lower triangle of the Hermitian matrix A, packed<br>
*>          columnwise in a linear array.  The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          If INFO = 0, S contains the scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCOND<br>
*> \verbatim<br>
*>          SCOND is REAL<br>
*>          If INFO = 0, S contains the ratio of the smallest S(i) to<br>
*>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too<br>
*>          large nor too small, it is not worth scaling by S.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AMAX<br>
*> \verbatim<br>
*>          AMAX is REAL<br>
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
*>          > 0:  if INFO = i, the i-th diagonal element is nonpositive.<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cppequ_(CHARACTER UPLO,INTEGER N,float[] AP,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
/**
*> \brief \b CPPRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPPRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpprfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpprfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpprfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR,<br>
*                          BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX            AFP( * ), AP( * ), B( LDB, * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPPRFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is Hermitian positive definite<br>
*> and packed, and provides error bounds and backward error estimates<br>
*> for the solution.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangle of the Hermitian matrix A, packed<br>
*>          columnwise in a linear array.  The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AFP<br>
*> \verbatim<br>
*>          AFP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H, as computed by SPPTRF/CPPTRF,<br>
*>          packed columnwise in a linear array in the same format as A<br>
*>          (see AP).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by CPPTRS.<br>
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
*>          FERR is REAL array, dimension (NRHS)<br>
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
*>          BERR is REAL array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief <b> CPPSV computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPPSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AP( * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPPSV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian positive definite matrix stored in<br>
*> packed format and X and B are N-by-NRHS matrices.<br>
*><br>
*> The Cholesky decomposition is used to factor A as<br>
*>    A = U**H * U,  if UPLO = 'U', or<br>
*>    A = L * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is a lower triangular<br>
*> matrix.  The factored form of A is then used to solve the system of<br>
*> equations A * X = B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          of the matrix B.  NRHS >= 0.<br>
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
*>          See below for further details.<br>
*><br>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H, in the same storage<br>
*>          format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          > 0:  if INFO = i, the leading minor of order i of A is not<br>
*>                positive definite, so the factorization could not be<br>
*>                completed, and the solution has not been computed.<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERsolve<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The packed storage scheme is illustrated by the following example<br>
*>  when N = 4, UPLO = 'U':<br>
*><br>
*>  Two-dimensional storage of the Hermitian matrix A:<br>
*><br>
*>     a11 a12 a13 a14<br>
*>         a22 a23 a24<br>
*>             a33 a34     (aij = conjg(aji))<br>
*>                 a44<br>
*><br>
*>  Packed storage of the upper triangle of A:<br>
*><br>
*>  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cppsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> CPPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPPSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppsvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppsvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppsvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB,<br>
*                          X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               BERR( * ), FERR( * ), RWORK( * ), S( * )<br>
*       COMPLEX            AFP( * ), AP( * ), B( LDB, * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPPSVX uses the Cholesky factorization A = U**H*U or A = L*L**H to<br>
*> compute the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian positive definite matrix stored in<br>
*> packed format and X and B are N-by-NRHS matrices.<br>
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
*>       diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B<br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.<br>
*><br>
*> 2. If FACT = 'N' or 'E', the Cholesky decomposition is used to<br>
*>    factor the matrix A (after equilibration if FACT = 'E') as<br>
*>       A = U**H * U ,  if UPLO = 'U', or<br>
*>       A = L * L**H,  if UPLO = 'L',<br>
*>    where U is an upper triangular matrix, L is a lower triangular<br>
*>    matrix, and **H indicates conjugate transpose.<br>
*><br>
*> 3. If the leading i-by-i principal minor is not positive definite,<br>
*>    then the routine returns with INFO = i. Otherwise, the factored<br>
*>    form of A is used to estimate the condition number of the matrix<br>
*>    A.  If the reciprocal of the condition number is less than machine<br>
*>    precision, INFO = N+1 is returned as a warning, but the routine<br>
*>    still goes on to solve for X and compute error bounds as<br>
*>    described below.<br>
*><br>
*> 4. The system of equations is solved for X using the factored form<br>
*>    of A.<br>
*><br>
*> 5. Iterative refinement is applied to improve the computed solution<br>
*>    matrix and calculate error bounds and backward error estimates<br>
*>    for it.<br>
*><br>
*> 6. If equilibration was used, the matrix X is premultiplied by<br>
*>    diag(S) so that it solves the original system before<br>
*>    equilibration.<br>
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
*>          = 'F':  On entry, AFP contains the factored form of A.<br>
*>                  If EQUED = 'Y', the matrix A has been equilibrated<br>
*>                  with scaling factors given by S.  AP and AFP will not<br>
*>                  be modified.<br>
*>          = 'N':  The matrix A will be copied to AFP and factored.<br>
*>          = 'E':  The matrix A will be equilibrated if necessary, then<br>
*>                  copied to AFP and factored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array, except if FACT = 'F'<br>
*>          and EQUED = 'Y', then A must contain the equilibrated matrix<br>
*>          diag(S)*A*diag(S).  The j-th column of A is stored in the<br>
*>          array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          See below for further details.  A is not modified if<br>
*>          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.<br>
*><br>
*>          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by<br>
*>          diag(S)*A*diag(S).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AFP<br>
*> \verbatim<br>
*>          AFP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          If FACT = 'F', then AFP is an input argument and on entry<br>
*>          contains the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H, in the same storage<br>
*>          format as A.  If EQUED .ne. 'N', then AFP is the factored<br>
*>          form of the equilibrated matrix A.<br>
*><br>
*>          If FACT = 'N', then AFP is an output argument and on exit<br>
*>          returns the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H * U or A = L * L**H of the original<br>
*>          matrix A.<br>
*><br>
*>          If FACT = 'E', then AFP is an output argument and on exit<br>
*>          returns the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H of the equilibrated<br>
*>          matrix A (see the description of AP for the form of the<br>
*>          equilibrated matrix).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EQUED<br>
*> \verbatim<br>
*>          EQUED is CHARACTER*1<br>
*>          Specifies the form of equilibration that was done.<br>
*>          = 'N':  No equilibration (always true if FACT = 'N').<br>
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by<br>
*>                  diag(S) * A * diag(S).<br>
*>          EQUED is an input argument if FACT = 'F'; otherwise, it is an<br>
*>          output argument.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
*>          The scale factors for A; not accessed if EQUED = 'N'.  S is<br>
*>          an input argument if FACT = 'F'; otherwise, S is an output<br>
*>          argument.  If FACT = 'F' and EQUED = 'Y', each element of S<br>
*>          must be positive.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>          On entry, the N-by-NRHS right hand side matrix B.<br>
*>          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',<br>
*>          B is overwritten by diag(S) * B.<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to<br>
*>          the original system of equations.  Note that if EQUED = 'Y',<br>
*>          A and B are modified on exit, and the solution to the<br>
*>          equilibrated system is inv(diag(S))*X.<br>
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
*>          RCOND is REAL<br>
*>          The estimate of the reciprocal condition number of the matrix<br>
*>          A after equilibration (if done).  If RCOND is less than the<br>
*>          machine precision (in particular, if RCOND = 0), the matrix<br>
*>          is singular to working precision.  This condition is<br>
*>          indicated by a return code of INFO > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is REAL array, dimension (NRHS)<br>
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
*>          BERR is REAL array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is<br>
*>                <= N:  the leading minor of order i of A is<br>
*>                       not positive definite, so the factorization<br>
*>                       could not be completed, and the solution has not<br>
*>                       been computed. RCOND = 0 is returned.<br>
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
*> \ingroup complexOTHERsolve<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The packed storage scheme is illustrated by the following example<br>
*>  when N = 4, UPLO = 'U':<br>
*><br>
*>  Two-dimensional storage of the Hermitian matrix A:<br>
*><br>
*>     a11 a12 a13 a14<br>
*>         a22 a23 a24<br>
*>             a33 a34     (aij = conjg(aji))<br>
*>                 a44<br>
*><br>
*>  Packed storage of the upper triangle of A:<br>
*><br>
*>  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cppsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPPTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPPTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpptrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpptrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpptrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPPTRF( UPLO, N, AP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPPTRF computes the Cholesky factorization of a complex Hermitian<br>
*> positive definite matrix A stored in packed format.<br>
*><br>
*> The factorization has the form<br>
*>    A = U**H * U,  if UPLO = 'U', or<br>
*>    A = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is lower triangular.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*>          See below for further details.<br>
*><br>
*>          On exit, if INFO = 0, the triangular factor U or L from the<br>
*>          Cholesky factorization A = U**H*U or A = L*L**H, in the same<br>
*>          storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the leading minor of order i is not<br>
*>                positive definite, and the factorization could not be<br>
*>                completed.<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The packed storage scheme is illustrated by the following example<br>
*>  when N = 4, UPLO = 'U':<br>
*><br>
*>  Two-dimensional storage of the Hermitian matrix A:<br>
*><br>
*>     a11 a12 a13 a14<br>
*>         a22 a23 a24<br>
*>             a33 a34     (aij = conjg(aji))<br>
*>                 a44<br>
*><br>
*>  Packed storage of the upper triangle of A:<br>
*><br>
*>  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cpptrf_(CHARACTER UPLO,INTEGER N,float[] AP,INTEGER INFO);
/**
*> \brief \b CPPTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPPTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpptri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpptri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpptri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPPTRI( UPLO, N, AP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPPTRI computes the inverse of a complex Hermitian positive definite<br>
*> matrix A using the Cholesky factorization A = U**H*U or A = L*L**H<br>
*> computed by CPPTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangular factor is stored in AP;<br>
*>          = 'L':  Lower triangular factor is stored in AP.<br>
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
*>          On entry, the triangular factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H, packed columnwise as<br>
*>          a linear array.  The j-th column of U or L is stored in the<br>
*>          array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the upper or lower triangle of the (Hermitian)<br>
*>          inverse of A, overwriting the input factor U or L.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the (i,i) element of the factor U or L is<br>
*>                zero, and the inverse could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpptri_(CHARACTER UPLO,INTEGER N,float[] AP,INTEGER INFO);
/**
*> \brief \b CPPTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPPTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpptrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpptrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpptrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            AP( * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPPTRS solves a system of linear equations A*X = B with a Hermitian<br>
*> positive definite matrix A in packed storage using the Cholesky<br>
*> factorization A = U**H*U or A = L*L**H computed by CPPTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored;<br>
*>          = 'L':  Lower triangle of A is stored.<br>
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
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array, dimension (N*(N+1)/2)<br>
*>          The triangular factor U or L from the Cholesky factorization<br>
*>          A = U**H*U or A = L*L**H, packed columnwise in a linear<br>
*>          array.  The j-th column of U or L is stored in the array AP<br>
*>          as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*> \ingroup complexOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b CPSTF2 computes the Cholesky factorization with complete pivoting of complex Hermitian positive semidefinite matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPSTF2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpstf2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpstf2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpstf2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               TOL<br>
*       INTEGER            INFO, LDA, N, RANK<br>
*       CHARACTER          UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       REAL               WORK( 2*N )<br>
*       INTEGER            PIV( N )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPSTF2 computes the Cholesky factorization with complete<br>
*> pivoting of a complex Hermitian positive semidefinite matrix A.<br>
*><br>
*> The factorization has the form<br>
*>    P**T * A * P = U**H * U ,  if UPLO = 'U',<br>
*>    P**T * A * P = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is lower triangular, and<br>
*> P is stored as vector PIV.<br>
*><br>
*> This algorithm does not attempt to check that A is positive<br>
*> semidefinite. This version of the algorithm calls level 2 BLAS.<br>
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
*>          On exit, if INFO = 0, the factor U or L from the Cholesky<br>
*>          factorization as above.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PIV<br>
*> \verbatim<br>
*>          PIV is INTEGER array, dimension (N)<br>
*>          PIV is such that the nonzero entries are P( PIV(K), K ) = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RANK<br>
*> \verbatim<br>
*>          RANK is INTEGER<br>
*>          The rank of A given by the number of steps the algorithm<br>
*>          completed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TOL<br>
*> \verbatim<br>
*>          TOL is REAL<br>
*>          User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) )<br>
*>          will be used. The algorithm terminates at the (K-1)st step<br>
*>          if the pivot <= TOL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
*>          Work space.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          < 0: If INFO = -K, the K-th argument had an illegal value,<br>
*>          = 0: algorithm completed successfully, and<br>
*>          > 0: the matrix A is either rank deficient with computed rank<br>
*>               as returned in RANK, or is not positive semidefinite. See<br>
*>               Section 7 of LAPACK Working Note #161 for further<br>
*>               information.<br>
*> \endverbatim<br>
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
*  =====================================================================<br>
*/
	public void cpstf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] PIV,INTEGER RANK,REAL TOL,float[] WORK,INTEGER INFO);
/**
*> \brief \b CPSTRF computes the Cholesky factorization with complete pivoting of complex Hermitian positive semidefinite matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPSTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpstrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpstrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpstrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               TOL<br>
*       INTEGER            INFO, LDA, N, RANK<br>
*       CHARACTER          UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       REAL               WORK( 2*N )<br>
*       INTEGER            PIV( N )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPSTRF computes the Cholesky factorization with complete<br>
*> pivoting of a complex Hermitian positive semidefinite matrix A.<br>
*><br>
*> The factorization has the form<br>
*>    P**T * A * P = U**H * U ,  if UPLO = 'U',<br>
*>    P**T * A * P = L  * L**H,  if UPLO = 'L',<br>
*> where U is an upper triangular matrix and L is lower triangular, and<br>
*> P is stored as vector PIV.<br>
*><br>
*> This algorithm does not attempt to check that A is positive<br>
*> semidefinite. This version of the algorithm calls level 3 BLAS.<br>
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
*>          On exit, if INFO = 0, the factor U or L from the Cholesky<br>
*>          factorization as above.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] PIV<br>
*> \verbatim<br>
*>          PIV is INTEGER array, dimension (N)<br>
*>          PIV is such that the nonzero entries are P( PIV(K), K ) = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RANK<br>
*> \verbatim<br>
*>          RANK is INTEGER<br>
*>          The rank of A given by the number of steps the algorithm<br>
*>          completed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TOL<br>
*> \verbatim<br>
*>          TOL is REAL<br>
*>          User defined tolerance. If TOL < 0, then N*U*MAX( A(K,K) )<br>
*>          will be used. The algorithm terminates at the (K-1)st step<br>
*>          if the pivot <= TOL.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
*>          Work space.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          < 0: If INFO = -K, the K-th argument had an illegal value,<br>
*>          = 0: algorithm completed successfully, and<br>
*>          > 0: the matrix A is either rank deficient with computed rank<br>
*>               as returned in RANK, or is not positive semidefinite. See<br>
*>               Section 7 of LAPACK Working Note #161 for further<br>
*>               information.<br>
*> \endverbatim<br>
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
*  =====================================================================<br>
*/
	public void cpstrf_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] PIV,INTEGER RANK,REAL TOL,float[] WORK,INTEGER INFO);
/**
*> \brief \b CPTCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPTCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, N<br>
*       REAL               ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), RWORK( * )<br>
*       COMPLEX            E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPTCON computes the reciprocal of the condition number (in the<br>
*> 1-norm) of a complex Hermitian positive definite tridiagonal matrix<br>
*> using the factorization A = L*D*L**H or A = U**H*D*U computed by<br>
*> CPTTRF.<br>
*><br>
*> Norm(inv(A)) is computed by a direct method, and the reciprocal of<br>
*> the condition number is computed as<br>
*>                  RCOND = 1 / (ANORM * norm(inv(A))).<br>
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
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The n diagonal elements of the diagonal matrix D from the<br>
*>          factorization of A, as computed by CPTTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) off-diagonal elements of the unit bidiagonal factor<br>
*>          U or L from the factorization of A, as computed by CPTTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is REAL<br>
*>          The 1-norm of the original matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RCOND<br>
*> \verbatim<br>
*>          RCOND is REAL<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the<br>
*>          1-norm of inv(A) computed in this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
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
*> \ingroup complexPTcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The method used is described in Nicholas J. Higham, "Efficient<br>
*>  Algorithms for Computing the Condition Number of a Tridiagonal<br>
*>  Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cptcon_(INTEGER N,float[] D,float[] E,REAL ANORM,REAL RCOND,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPTEQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPTEQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpteqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpteqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpteqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPZ<br>
*       INTEGER            INFO, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * ), WORK( * )<br>
*       COMPLEX            Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPTEQR computes all eigenvalues and, optionally, eigenvectors of a<br>
*> symmetric positive definite tridiagonal matrix by first factoring the<br>
*> matrix using SPTTRF and then calling CBDSQR to compute the singular<br>
*> values of the bidiagonal factor.<br>
*><br>
*> This routine computes the eigenvalues of the positive definite<br>
*> tridiagonal matrix to high relative accuracy.  This means that if the<br>
*> eigenvalues range over many orders of magnitude in size, then the<br>
*> small eigenvalues and corresponding eigenvectors will be computed<br>
*> more accurately than, for example, with the standard QR method.<br>
*><br>
*> The eigenvectors of a full or band positive definite Hermitian matrix<br>
*> can also be found if CHETRD, CHPTRD, or CHBTRD has been used to<br>
*> reduce this matrix to tridiagonal form.  (The reduction to<br>
*> tridiagonal form, however, may preclude the possibility of obtaining<br>
*> high relative accuracy in the small eigenvalues of the original<br>
*> matrix, if these eigenvalues range over many orders of magnitude.)<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only.<br>
*>          = 'V':  Compute eigenvectors of original Hermitian<br>
*>                  matrix also.  Array Z contains the unitary matrix<br>
*>                  used to reduce the original matrix to tridiagonal<br>
*>                  form.<br>
*>          = 'I':  Compute eigenvectors of tridiagonal matrix also.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the n diagonal elements of the tridiagonal matrix.<br>
*>          On normal exit, D contains the eigenvalues, in descending<br>
*>          order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal<br>
*>          matrix.<br>
*>          On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX array, dimension (LDZ, N)<br>
*>          On entry, if COMPZ = 'V', the unitary matrix used in the<br>
*>          reduction to tridiagonal form.<br>
*>          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the<br>
*>          original Hermitian matrix;<br>
*>          if COMPZ = 'I', the orthonormal eigenvectors of the<br>
*>          tridiagonal matrix.<br>
*>          If INFO > 0 on exit, Z contains the eigenvectors associated<br>
*>          with only the stored eigenvalues.<br>
*>          If  COMPZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          COMPZ = 'V' or 'I', LDZ >= max(1,N).<br>
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
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = i, and i is:<br>
*>                <= N  the Cholesky factorization of the matrix could<br>
*>                      not be performed because the i-th principal minor<br>
*>                      was not positive definite.<br>
*>                > N   the SVD algorithm failed to converge;<br>
*>                      if INFO = N+i, i off-diagonal elements of the<br>
*>                      bidiagonal factor did not converge to zero.<br>
*> \endverbatim<br>
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
*> \ingroup complexPTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpteqr_(CHARACTER COMPZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
/**
*> \brief \b CPTRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPTRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,<br>
*                          FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               BERR( * ), D( * ), DF( * ), FERR( * ),<br>
*      $                   RWORK( * )<br>
*       COMPLEX            B( LDB, * ), E( * ), EF( * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPTRFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is Hermitian positive definite<br>
*> and tridiagonal, and provides error bounds and backward error<br>
*> estimates for the solution.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the superdiagonal or the subdiagonal of the<br>
*>          tridiagonal matrix A is stored and the form of the<br>
*>          factorization:<br>
*>          = 'U':  E is the superdiagonal of A, and A = U**H*D*U;<br>
*>          = 'L':  E is the subdiagonal of A, and A = L*D*L**H.<br>
*>          (The two forms are equivalent if A is real.)<br>
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
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The n real diagonal elements of the tridiagonal matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) off-diagonal elements of the tridiagonal matrix A<br>
*>          (see UPLO).<br>
*> \endverbatim<br>
*><br>
*> \param[in] DF<br>
*> \verbatim<br>
*>          DF is REAL array, dimension (N)<br>
*>          The n diagonal elements of the diagonal matrix D from<br>
*>          the factorization computed by CPTTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] EF<br>
*> \verbatim<br>
*>          EF is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) off-diagonal elements of the unit bidiagonal<br>
*>          factor U or L from the factorization computed by CPTTRF<br>
*>          (see UPLO).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by CPTTRS.<br>
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
*>          FERR is REAL array, dimension (NRHS)<br>
*>          The forward error bound for each solution vector<br>
*>          X(j) (the j-th column of the solution matrix X).<br>
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)<br>
*>          is an estimated upper bound for the magnitude of the largest<br>
*>          element in (X(j) - XTRUE) divided by the magnitude of the<br>
*>          largest element in X(j).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is REAL array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in<br>
*>          any element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
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
*> \ingroup complexPTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cptrfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] DF,float[] EF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief <b> CPTSV computes the solution to system of linear equations A * X = B for PT matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPTSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPTSV( N, NRHS, D, E, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * )<br>
*       COMPLEX            B( LDB, * ), E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPTSV computes the solution to a complex system of linear equations<br>
*> A*X = B, where A is an N-by-N Hermitian positive definite tridiagonal<br>
*> matrix, and X and B are N-by-NRHS matrices.<br>
*><br>
*> A is factored as A = L*D*L**H, and the factored form of A is then<br>
*> used to solve the system of equations.<br>
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
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the n diagonal elements of the tridiagonal matrix<br>
*>          A.  On exit, the n diagonal elements of the diagonal matrix<br>
*>          D from the factorization A = L*D*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is COMPLEX array, dimension (N-1)<br>
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal<br>
*>          matrix A.  On exit, the (n-1) subdiagonal elements of the<br>
*>          unit bidiagonal factor L from the L*D*L**H factorization of<br>
*>          A.  E can also be regarded as the superdiagonal of the unit<br>
*>          bidiagonal factor U from the U**H*D*U factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          > 0:  if INFO = i, the leading minor of order i is not<br>
*>                positive definite, and the solution has not been<br>
*>                computed.  The factorization has not been completed<br>
*>                unless i = N.<br>
*> \endverbatim<br>
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
*> \ingroup complexPTsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void cptsv_(INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> CPTSVX computes the solution to system of linear equations A * X = B for PT matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPTSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptsvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptsvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptsvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,<br>
*                          RCOND, FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          FACT<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               BERR( * ), D( * ), DF( * ), FERR( * ),<br>
*      $                   RWORK( * )<br>
*       COMPLEX            B( LDB, * ), E( * ), EF( * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPTSVX uses the factorization A = L*D*L**H to compute the solution<br>
*> to a complex system of linear equations A*X = B, where A is an<br>
*> N-by-N Hermitian positive definite tridiagonal matrix and X and B<br>
*> are N-by-NRHS matrices.<br>
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
*> 1. If FACT = 'N', the matrix A is factored as A = L*D*L**H, where L<br>
*>    is a unit lower bidiagonal matrix and D is diagonal.  The<br>
*>    factorization can also be regarded as having the form<br>
*>    A = U**H*D*U.<br>
*><br>
*> 2. If the leading i-by-i principal minor is not positive definite,<br>
*>    then the routine returns with INFO = i. Otherwise, the factored<br>
*>    form of A is used to estimate the condition number of the matrix<br>
*>    A.  If the reciprocal of the condition number is less than machine<br>
*>    precision, INFO = N+1 is returned as a warning, but the routine<br>
*>    still goes on to solve for X and compute error bounds as<br>
*>    described below.<br>
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
*>          Specifies whether or not the factored form of the matrix<br>
*>          A is supplied on entry.<br>
*>          = 'F':  On entry, DF and EF contain the factored form of A.<br>
*>                  D, E, DF, and EF will not be modified.<br>
*>          = 'N':  The matrix A will be copied to DF and EF and<br>
*>                  factored.<br>
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
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The n diagonal elements of the tridiagonal matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is COMPLEX array, dimension (N-1)<br>
*>          The (n-1) subdiagonal elements of the tridiagonal matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DF<br>
*> \verbatim<br>
*>          DF is REAL array, dimension (N)<br>
*>          If FACT = 'F', then DF is an input argument and on entry<br>
*>          contains the n diagonal elements of the diagonal matrix D<br>
*>          from the L*D*L**H factorization of A.<br>
*>          If FACT = 'N', then DF is an output argument and on exit<br>
*>          contains the n diagonal elements of the diagonal matrix D<br>
*>          from the L*D*L**H factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] EF<br>
*> \verbatim<br>
*>          EF is COMPLEX array, dimension (N-1)<br>
*>          If FACT = 'F', then EF is an input argument and on entry<br>
*>          contains the (n-1) subdiagonal elements of the unit<br>
*>          bidiagonal factor L from the L*D*L**H factorization of A.<br>
*>          If FACT = 'N', then EF is an output argument and on exit<br>
*>          contains the (n-1) subdiagonal elements of the unit<br>
*>          bidiagonal factor L from the L*D*L**H factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX array, dimension (LDX,NRHS)<br>
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
*>          RCOND is REAL<br>
*>          The reciprocal condition number of the matrix A.  If RCOND<br>
*>          is less than the machine precision (in particular, if<br>
*>          RCOND = 0), the matrix is singular to working precision.<br>
*>          This condition is indicated by a return code of INFO > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] FERR<br>
*> \verbatim<br>
*>          FERR is REAL array, dimension (NRHS)<br>
*>          The forward error bound for each solution vector<br>
*>          X(j) (the j-th column of the solution matrix X).<br>
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)<br>
*>          is an estimated upper bound for the magnitude of the largest<br>
*>          element in (X(j) - XTRUE) divided by the magnitude of the<br>
*>          largest element in X(j).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BERR<br>
*> \verbatim<br>
*>          BERR is REAL array, dimension (NRHS)<br>
*>          The componentwise relative backward error of each solution<br>
*>          vector X(j) (i.e., the smallest relative change in any<br>
*>          element of A or B that makes X(j) an exact solution).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is<br>
*>                <= N:  the leading minor of order i of A is<br>
*>                       not positive definite, so the factorization<br>
*>                       could not be completed, and the solution has not<br>
*>                       been computed. RCOND = 0 is returned.<br>
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
*> \ingroup complexPTsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void cptsvx_(CHARACTER FACT,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] DF,float[] EF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
/**
*> \brief \b CPTTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPTTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpttrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpttrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpttrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPTTRF( N, D, E, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, N<br>
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
*> CPTTRF computes the L*D*L**H factorization of a complex Hermitian<br>
*> positive definite tridiagonal matrix A.  The factorization may also<br>
*> be regarded as having the form A = U**H *D*U.<br>
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
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the n diagonal elements of the tridiagonal matrix<br>
*>          A.  On exit, the n diagonal elements of the diagonal matrix<br>
*>          D from the L*D*L**H factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is COMPLEX array, dimension (N-1)<br>
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal<br>
*>          matrix A.  On exit, the (n-1) subdiagonal elements of the<br>
*>          unit bidiagonal factor L from the L*D*L**H factorization of A.<br>
*>          E can also be regarded as the superdiagonal of the unit<br>
*>          bidiagonal factor U from the U**H *D*U factorization of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -k, the k-th argument had an illegal value<br>
*>          > 0: if INFO = k, the leading minor of order k is not<br>
*>               positive definite; if k < N, the factorization could not<br>
*>               be completed, while if k = N, the factorization was<br>
*>               completed, but D(N) <= 0.<br>
*> \endverbatim<br>
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
*> \ingroup complexPTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpttrf_(INTEGER N,float[] D,float[] E,INTEGER INFO);
/**
*> \brief \b CPTTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPTTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpttrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpttrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpttrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * )<br>
*       COMPLEX            B( LDB, * ), E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPTTRS solves a tridiagonal system of the form<br>
*>    A * X = B<br>
*> using the factorization A = U**H*D*U or A = L*D*L**H computed by CPTTRF.<br>
*> D is a diagonal matrix specified in the vector D, U (or L) is a unit<br>
*> bidiagonal matrix whose superdiagonal (subdiagonal) is specified in<br>
*> the vector E, and X and B are N by NRHS matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies the form of the factorization and whether the<br>
*>          vector E is the superdiagonal of the upper bidiagonal factor<br>
*>          U or the subdiagonal of the lower bidiagonal factor L.<br>
*>          = 'U':  A = U**H*D*U, E is the superdiagonal of U<br>
*>          = 'L':  A = L*D*L**H, E is the subdiagonal of L<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the tridiagonal matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The n diagonal elements of the diagonal matrix D from the<br>
*>          factorization A = U**H*D*U or A = L*D*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is COMPLEX array, dimension (N-1)<br>
*>          If UPLO = 'U', the (n-1) superdiagonal elements of the unit<br>
*>          bidiagonal factor U from the factorization A = U**H*D*U.<br>
*>          If UPLO = 'L', the (n-1) subdiagonal elements of the unit<br>
*>          bidiagonal factor L from the factorization A = L*D*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>          On entry, the right hand side vectors B for the system of<br>
*>          linear equations.<br>
*>          On exit, the solution vectors, X.<br>
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
*> \date June 2016<br>
*<br>
*> \ingroup complexPTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cpttrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b CPTTS2 solves a tridiagonal system of the form AX=B using the L D LH factorization computed by spttrf.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CPTTS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptts2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptts2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptts2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CPTTS2( IUPLO, N, NRHS, D, E, B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IUPLO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * )<br>
*       COMPLEX            B( LDB, * ), E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CPTTS2 solves a tridiagonal system of the form<br>
*>    A * X = B<br>
*> using the factorization A = U**H*D*U or A = L*D*L**H computed by CPTTRF.<br>
*> D is a diagonal matrix specified in the vector D, U (or L) is a unit<br>
*> bidiagonal matrix whose superdiagonal (subdiagonal) is specified in<br>
*> the vector E, and X and B are N by NRHS matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] IUPLO<br>
*> \verbatim<br>
*>          IUPLO is INTEGER<br>
*>          Specifies the form of the factorization and whether the<br>
*>          vector E is the superdiagonal of the upper bidiagonal factor<br>
*>          U or the subdiagonal of the lower bidiagonal factor L.<br>
*>          = 1:  A = U**H *D*U, E is the superdiagonal of U<br>
*>          = 0:  A = L*D*L**H, E is the subdiagonal of L<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the tridiagonal matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NRHS<br>
*> \verbatim<br>
*>          NRHS is INTEGER<br>
*>          The number of right hand sides, i.e., the number of columns<br>
*>          of the matrix B.  NRHS >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The n diagonal elements of the diagonal matrix D from the<br>
*>          factorization A = U**H *D*U or A = L*D*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is COMPLEX array, dimension (N-1)<br>
*>          If IUPLO = 1, the (n-1) superdiagonal elements of the unit<br>
*>          bidiagonal factor U from the factorization A = U**H*D*U.<br>
*>          If IUPLO = 0, the (n-1) subdiagonal elements of the unit<br>
*>          bidiagonal factor L from the factorization A = L*D*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array, dimension (LDB,NRHS)<br>
*>          On entry, the right hand side vectors B for the system of<br>
*>          linear equations.<br>
*>          On exit, the solution vectors, X.<br>
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
*> \date June 2016<br>
*<br>
*> \ingroup complexPTcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void cptts2_(INTEGER IUPLO,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB);

}