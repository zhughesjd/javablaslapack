package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZT extends Library
{

	public static LapackZT instance = (LapackZT) Native.loadLibrary("liblapack",LapackZT.class);

/**
*> \brief \b ZTBCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTBCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztbcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztbcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztbcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK,<br>
*                          RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         AB( LDAB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTBCON estimates the reciprocal of the condition number of a<br>
*> triangular band matrix A, in either the 1-norm or the infinity-norm.<br>
*><br>
*> The norm of A is computed and an estimate is obtained for<br>
*> norm(inv(A)), then the reciprocal of the condition number is<br>
*> computed as<br>
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
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          The number of superdiagonals or subdiagonals of the<br>
*>          triangular band matrix A.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX*16 array, dimension (LDAB,N)<br>
*>          The upper or lower triangular band matrix A, stored in the<br>
*>          first kd+1 rows of the array. The j-th column of A is stored<br>
*>          in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*>          If DIAG = 'U', the diagonal elements of A are not referenced<br>
*>          and are assumed to be 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztbcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTBRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTBRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztbrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztbrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztbrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,<br>
*                          LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX*16         AB( LDAB, * ), B( LDB, * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTBRFS provides error bounds and backward error estimates for the<br>
*> solution to a system of linear equations with a triangular band<br>
*> coefficient matrix.<br>
*><br>
*> The solution matrix X must be computed by ZTBTRS or some other<br>
*> means before entering this routine.  ZTBRFS does not do iterative<br>
*> refinement because doing so cannot improve the backward error.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          The number of superdiagonals or subdiagonals of the<br>
*>          triangular band matrix A.  KD >= 0.<br>
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
*>          AB is COMPLEX*16 array, dimension (LDAB,N)<br>
*>          The upper or lower triangular band matrix A, stored in the<br>
*>          first kd+1 rows of the array. The j-th column of A is stored<br>
*>          in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*>          If DIAG = 'U', the diagonal elements of A are not referenced<br>
*>          and are assumed to be 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
*>          The right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
*>          The solution matrix X.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztbrfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTBTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTBTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztbtrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztbtrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztbtrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,<br>
*                          LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         AB( LDAB, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTBTRS solves a triangular system of the form<br>
*><br>
*>    A * X = B,  A**T * X = B,  or  A**H * X = B,<br>
*><br>
*> where A is a triangular band matrix of order N, and B is an<br>
*> N-by-NRHS matrix.  A check is made to verify that A is nonsingular.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          The number of superdiagonals or subdiagonals of the<br>
*>          triangular band matrix A.  KD >= 0.<br>
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
*>          AB is COMPLEX*16 array, dimension (LDAB,N)<br>
*>          The upper or lower triangular band matrix A, stored in the<br>
*>          first kd+1 rows of AB.  The j-th column of A is stored<br>
*>          in the j-th column of the array AB as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*>          If DIAG = 'U', the diagonal elements of A are not referenced<br>
*>          and are assumed to be 1.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
*>          On entry, the right hand side matrix B.<br>
*>          On exit, if INFO = 0, the solution matrix X.<br>
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
*>          > 0:  if INFO = i, the i-th diagonal element of A is zero,<br>
*>                indicating that the matrix is singular and the<br>
*>                solutions X have not been computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztbtrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZTFSM solves a matrix equation (one operand is a triangular matrix in RFP format).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTFSM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztfsm.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztfsm.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztfsm.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A,<br>
*                         B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO<br>
*       INTEGER            LDB, M, N<br>
*       COMPLEX*16         ALPHA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( 0: * ), B( 0: LDB-1, 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Level 3 BLAS like routine for A in RFP Format.<br>
*><br>
*> ZTFSM  solves the matrix equation<br>
*><br>
*>    op( A )*X = alpha*B  or  X*op( A ) = alpha*B<br>
*><br>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**H.<br>
*><br>
*> A is in Rectangular Full Packed (RFP) Format.<br>
*><br>
*> The matrix X is overwritten on B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  The Normal Form of RFP A is stored;<br>
*>          = 'C':  The Conjugate-transpose Form of RFP A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry, SIDE specifies whether op( A ) appears on the left<br>
*>           or right of X as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the RFP matrix A came from<br>
*>           an upper or lower triangular matrix as follows:<br>
*>           UPLO = 'U' or 'u' RFP A came from an upper triangular matrix<br>
*>           UPLO = 'L' or 'l' RFP A came from a  lower triangular matrix<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS  specifies the form of op( A ) to be used<br>
*>           in the matrix multiplication as follows:<br>
*><br>
*>              TRANS  = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANS  = 'C' or 'c'   op( A ) = conjg( A' ).<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not RFP A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>           NT = N*(N+1)/2. On entry, the matrix A in RFP Format.<br>
*>           RFP Format is described by TRANSR, UPLO and N as follows:<br>
*>           If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;<br>
*>           K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If<br>
*>           TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A as<br>
*>           defined when TRANSR = 'N'. The contents of RFP A are defined<br>
*>           by UPLO as follows: If UPLO = 'U' the RFP A contains the NT<br>
*>           elements of upper packed A either in normal or<br>
*>           conjugate-transpose Format. If UPLO = 'L' the RFP A contains<br>
*>           the NT elements of lower packed A either in normal or<br>
*>           conjugate-transpose Format. The LDA of RFP A is (N+1)/2 when<br>
*>           TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is<br>
*>           even and is N when is odd.<br>
*>           See the Note below for more details. Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>           Before entry,  the leading  m by n part of the array  B must<br>
*>           contain  the  right-hand  side  matrix  B,  and  on exit  is<br>
*>           overwritten by the solution matrix  X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
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
*> \ingroup complex16OTHERcomputational<br>
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
	public void ztfsm_(CHARACTER TRANSR,CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER M,INTEGER N,double[] ALPHA,double[] A,double[] B,INTEGER LDB);
/**
*> \brief \b ZTFTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTFTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztftri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztftri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztftri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTFTRI( TRANSR, UPLO, DIAG, N, A, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO, DIAG<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTFTRI computes the inverse of a triangular matrix A stored in RFP<br>
*> format.<br>
*><br>
*> This is a Level 3 BLAS version of the algorithm.<br>
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
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          A is COMPLEX*16 array, dimension ( N*(N+1)/2 );<br>
*>          On entry, the triangular matrix A in RFP format. RFP format<br>
*>          is described by TRANSR, UPLO, and N as follows: If TRANSR =<br>
*>          'N' then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is<br>
*>          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is<br>
*>          the Conjugate-transpose of RFP A as defined when<br>
*>          TRANSR = 'N'. The contents of RFP A are defined by UPLO as<br>
*>          follows: If UPLO = 'U' the RFP A contains the nt elements of<br>
*>          upper packed A; If UPLO = 'L' the RFP A contains the nt<br>
*>          elements of lower packed A. The LDA of RFP A is (N+1)/2 when<br>
*>          TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is<br>
*>          even and N is odd. See the Note below for more details.<br>
*><br>
*>          On exit, the (triangular) inverse of the original matrix, in<br>
*>          the same storage format.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular<br>
*>               matrix is singular and its inverse can not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
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
	public void ztftri_(CHARACTER TRANSR,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] A,INTEGER INFO);
/**
*> \brief \b ZTFTTP copies a triangular matrix from the rectangular full packed format (TF) to the standard packed format (TP).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTFTTP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztfttp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztfttp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztfttp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTFTTP( TRANSR, UPLO, N, ARF, AP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         AP( 0: * ), ARF( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTFTTP copies a triangular matrix A from rectangular full packed<br>
*> format (TF) to standard packed format (TP).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  ARF is in Normal format;<br>
*>          = 'C':  ARF is in Conjugate-transpose format;<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ARF<br>
*> \verbatim<br>
*>          ARF is COMPLEX*16 array, dimension ( N*(N+1)/2 ),<br>
*>          On entry, the upper or lower triangular matrix A stored in<br>
*>          RFP format. For a further discussion see Notes below.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension ( N*(N+1)/2 ),<br>
*>          On exit, the upper or lower triangular matrix A, packed<br>
*>          columnwise in a linear array. The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
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
*> \ingroup complex16OTHERcomputational<br>
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
*>  We next consider Standard Packed Format when N is odd.<br>
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
	public void ztfttp_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] ARF,double[] AP,INTEGER INFO);
/**
*> \brief \b ZTFTTR copies a triangular matrix from the rectangular full packed format (TF) to the standard full format (TR).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTFTTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztfttr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztfttr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztfttr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( 0: LDA-1, 0: * ), ARF( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTFTTR copies a triangular matrix A from rectangular full packed<br>
*> format (TF) to standard full format (TR).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  ARF is in Normal format;<br>
*>          = 'C':  ARF is in Conjugate-transpose format;<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ARF<br>
*> \verbatim<br>
*>          ARF is COMPLEX*16 array, dimension ( N*(N+1)/2 ),<br>
*>          On entry, the upper or lower triangular matrix A stored in<br>
*>          RFP format. For a further discussion see Notes below.<br>
*> \endverbatim<br>
*><br>
*> \param[out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension ( LDA, N )<br>
*>          On exit, the triangular matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of the array A contains<br>
*>          the upper triangular matrix, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of the array A contains<br>
*>          the lower triangular matrix, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
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
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
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
	public void ztfttr_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] ARF,double[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b ZTGEVC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTGEVC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgevc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgevc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgevc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,<br>
*                          LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, SIDE<br>
*       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         P( LDP, * ), S( LDS, * ), VL( LDVL, * ),<br>
*      $                   VR( LDVR, * ), WORK( * )<br>
*       ..<br>
*  <br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTGEVC computes some or all of the right and/or left eigenvectors of<br>
*> a pair of complex matrices (S,P), where S and P are upper triangular.<br>
*> Matrix pairs of this type are produced by the generalized Schur<br>
*> factorization of a complex matrix pair (A,B):<br>
*> <br>
*>    A = Q*S*Z**H,  B = Q*P*Z**H<br>
*> <br>
*> as computed by ZGGHRD + ZHGEQZ.<br>
*> <br>
*> The right eigenvector x and the left eigenvector y of (S,P)<br>
*> corresponding to an eigenvalue w are defined by:<br>
*> <br>
*>    S*x = w*P*x,  (y**H)*S = w*(y**H)*P,<br>
*> <br>
*> where y**H denotes the conjugate tranpose of y.<br>
*> The eigenvalues are not input to this routine, but are computed<br>
*> directly from the diagonal elements of S and P.<br>
*> <br>
*> This routine returns the matrices X and/or Y of right and left<br>
*> eigenvectors of (S,P), or the products Z*X and/or Q*Y,<br>
*> where Z and Q are input matrices.<br>
*> If Q and Z are the unitary factors from the generalized Schur<br>
*> factorization of a matrix pair (A,B), then Z*X and Q*Y<br>
*> are the matrices of right and left eigenvectors of (A,B).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'R': compute right eigenvectors only;<br>
*>          = 'L': compute left eigenvectors only;<br>
*>          = 'B': compute both right and left eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] HOWMNY<br>
*> \verbatim<br>
*>          HOWMNY is CHARACTER*1<br>
*>          = 'A': compute all right and/or left eigenvectors;<br>
*>          = 'B': compute all right and/or left eigenvectors,<br>
*>                 backtransformed by the matrices in VR and/or VL;<br>
*>          = 'S': compute selected right and/or left eigenvectors,<br>
*>                 specified by the logical array SELECT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          If HOWMNY='S', SELECT specifies the eigenvectors to be<br>
*>          computed.  The eigenvector corresponding to the j-th<br>
*>          eigenvalue is computed if SELECT(j) = .TRUE..<br>
*>          Not referenced if HOWMNY = 'A' or 'B'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices S and P.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is COMPLEX*16 array, dimension (LDS,N)<br>
*>          The upper triangular matrix S from a generalized Schur<br>
*>          factorization, as computed by ZHGEQZ.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDS<br>
*> \verbatim<br>
*>          LDS is INTEGER<br>
*>          The leading dimension of array S.  LDS >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] P<br>
*> \verbatim<br>
*>          P is COMPLEX*16 array, dimension (LDP,N)<br>
*>          The upper triangular matrix P from a generalized Schur<br>
*>          factorization, as computed by ZHGEQZ.  P must have real<br>
*>          diagonal elements.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDP<br>
*> \verbatim<br>
*>          LDP is INTEGER<br>
*>          The leading dimension of array P.  LDP >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is COMPLEX*16 array, dimension (LDVL,MM)<br>
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must<br>
*>          contain an N-by-N matrix Q (usually the unitary matrix Q<br>
*>          of left Schur vectors returned by ZHGEQZ).<br>
*>          On exit, if SIDE = 'L' or 'B', VL contains:<br>
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);<br>
*>          if HOWMNY = 'B', the matrix Q*Y;<br>
*>          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by<br>
*>                      SELECT, stored consecutively in the columns of<br>
*>                      VL, in the same order as their eigenvalues.<br>
*>          Not referenced if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of array VL.  LDVL >= 1, and if<br>
*>          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VR<br>
*> \verbatim<br>
*>          VR is COMPLEX*16 array, dimension (LDVR,MM)<br>
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must<br>
*>          contain an N-by-N matrix Q (usually the unitary matrix Z<br>
*>          of right Schur vectors returned by ZHGEQZ).<br>
*>          On exit, if SIDE = 'R' or 'B', VR contains:<br>
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);<br>
*>          if HOWMNY = 'B', the matrix Z*X;<br>
*>          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by<br>
*>                      SELECT, stored consecutively in the columns of<br>
*>                      VR, in the same order as their eigenvalues.<br>
*>          Not referenced if SIDE = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR.  LDVR >= 1, and if<br>
*>          SIDE = 'R' or 'B', LDVR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MM<br>
*> \verbatim<br>
*>          MM is INTEGER<br>
*>          The number of columns in the arrays VL and/or VR. MM >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of columns in the arrays VL and/or VR actually<br>
*>          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M<br>
*>          is set to N.  Each selected eigenvector occupies one column.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (2*N)<br>
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
*> \ingroup complex16GEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztgevc_(CHARACTER SIDE,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,double[] S,INTEGER LDS,double[] P,INTEGER LDP,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an unitary equivalence transformation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTGEX2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgex2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgex2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgex2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,<br>
*                          LDZ, J1, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            WANTQ, WANTZ<br>
*       INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTGEX2 swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)<br>
*> in an upper triangular matrix pair (A, B) by an unitary equivalence<br>
*> transformation.<br>
*><br>
*> (A, B) must be in generalized Schur canonical form, that is, A and<br>
*> B are both upper triangular.<br>
*><br>
*> Optionally, the matrices Q and Z of generalized Schur vectors are<br>
*> updated.<br>
*><br>
*>        Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H<br>
*>        Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTQ<br>
*> \verbatim<br>
*>          WANTQ is LOGICAL<br>
*>          .TRUE. : update the left transformation matrix Q;<br>
*>          .FALSE.: do not update Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          .TRUE. : update the right transformation matrix Z;<br>
*>          .FALSE.: do not update Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 arrays, dimensions (LDA,N)<br>
*>          On entry, the matrix A in the pair (A, B).<br>
*>          On exit, the updated matrix A.<br>
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
*>          B is COMPLEX*16 arrays, dimensions (LDB,N)<br>
*>          On entry, the matrix B in the pair (A, B).<br>
*>          On exit, the updated matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX*16 array, dimension (LDZ,N)<br>
*>          If WANTQ = .TRUE, on entry, the unitary matrix Q. On exit,<br>
*>          the updated matrix Q.<br>
*>          Not referenced if WANTQ = .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q. LDQ >= 1;<br>
*>          If WANTQ = .TRUE., LDQ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ,N)<br>
*>          If WANTZ = .TRUE, on entry, the unitary matrix Z. On exit,<br>
*>          the updated matrix Z.<br>
*>          Not referenced if WANTZ = .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z. LDZ >= 1;<br>
*>          If WANTZ = .TRUE., LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] J1<br>
*> \verbatim<br>
*>          J1 is INTEGER<br>
*>          The index to the first block (A11, B11).<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           =0:  Successful exit.<br>
*>           =1:  The transformed matrix pair (A, B) would be too far<br>
*>                from generalized Schur form; the problem is ill-<br>
*>                conditioned. <br>
*> \endverbatim<br>
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
*> \ingroup complex16GEauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*>  In the current code both weak and strong stability tests are<br>
*>  performed. The user can omit the strong stability test by changing<br>
*>  the internal logical parameter WANDS to .FALSE.. See ref. [2] for<br>
*>  details.<br>
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
*>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the<br>
*>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in<br>
*>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and<br>
*>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.<br>
*> \n<br>
*>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified<br>
*>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition<br>
*>      Estimation: Theory, Algorithms and Software, Report UMINF-94.04,<br>
*>      Department of Computing Science, Umea University, S-901 87 Umea,<br>
*>      Sweden, 1994. Also as LAPACK Working Note 87. To appear in<br>
*>      Numerical Algorithms, 1996.<br>
*><br>
*  =====================================================================<br>
*/
	public void ztgex2_(LOGICAL WANTQ,LOGICAL WANTZ,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,INTEGER J1,INTEGER INFO);
/**
*> \brief \b ZTGEXC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTGEXC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgexc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgexc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgexc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,<br>
*                          LDZ, IFST, ILST, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            WANTQ, WANTZ<br>
*       INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTGEXC reorders the generalized Schur decomposition of a complex<br>
*> matrix pair (A,B), using an unitary equivalence transformation<br>
*> (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with<br>
*> row index IFST is moved to row ILST.<br>
*><br>
*> (A, B) must be in generalized Schur canonical form, that is, A and<br>
*> B are both upper triangular.<br>
*><br>
*> Optionally, the matrices Q and Z of generalized Schur vectors are<br>
*> updated.<br>
*><br>
*>        Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H<br>
*>        Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] WANTQ<br>
*> \verbatim<br>
*>          WANTQ is LOGICAL<br>
*>          .TRUE. : update the left transformation matrix Q;<br>
*>          .FALSE.: do not update Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          .TRUE. : update the right transformation matrix Z;<br>
*>          .FALSE.: do not update Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the upper triangular matrix A in the pair (A, B).<br>
*>          On exit, the updated matrix A.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>          On entry, the upper triangular matrix B in the pair (A, B).<br>
*>          On exit, the updated matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX*16 array, dimension (LDZ,N)<br>
*>          On entry, if WANTQ = .TRUE., the unitary matrix Q.<br>
*>          On exit, the updated matrix Q.<br>
*>          If WANTQ = .FALSE., Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q. LDQ >= 1;<br>
*>          If WANTQ = .TRUE., LDQ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ,N)<br>
*>          On entry, if WANTZ = .TRUE., the unitary matrix Z.<br>
*>          On exit, the updated matrix Z.<br>
*>          If WANTZ = .FALSE., Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z. LDZ >= 1;<br>
*>          If WANTZ = .TRUE., LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IFST<br>
*> \verbatim<br>
*>          IFST is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ILST<br>
*> \verbatim<br>
*>          ILST is INTEGER<br>
*>          Specify the reordering of the diagonal blocks of (A, B).<br>
*>          The block with row index IFST is moved to row ILST, by a<br>
*>          sequence of swapping between adjacent blocks.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>           =0:  Successful exit.<br>
*>           <0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>           =1:  The transformed matrix pair (A, B) would be too far<br>
*>                from generalized Schur form; the problem is ill-<br>
*>                conditioned. (A, B) may have been partially reordered,<br>
*>                and ILST points to the first row of the current<br>
*>                position of the block being moved.<br>
*> \endverbatim<br>
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
*> \ingroup complex16GEcomputational<br>
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
*>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the<br>
*>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in<br>
*>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and<br>
*>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.<br>
*> \n<br>
*>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified<br>
*>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition<br>
*>      Estimation: Theory, Algorithms and Software, Report<br>
*>      UMINF - 94.04, Department of Computing Science, Umea University,<br>
*>      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.<br>
*>      To appear in Numerical Algorithms, 1996.<br>
*> \n<br>
*>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software<br>
*>      for Solving the Generalized Sylvester Equation and Estimating the<br>
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,<br>
*>      Department of Computing Science, Umea University, S-901 87 Umea,<br>
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK working<br>
*>      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,<br>
*>      1996.<br>
*><br>
*  =====================================================================<br>
*/
	public void ztgexc_(LOGICAL WANTQ,LOGICAL WANTZ,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,INTEGER IFST,INTEGER ILST,INTEGER INFO);
/**
*> \brief \b ZTGSEN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTGSEN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsen.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsen.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsen.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB,<br>
*                          ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF,<br>
*                          WORK, LWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            WANTQ, WANTZ<br>
*       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,<br>
*      $                   M, N<br>
*       DOUBLE PRECISION   PL, PR<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   DIF( * )<br>
*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),<br>
*      $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTGSEN reorders the generalized Schur decomposition of a complex<br>
*> matrix pair (A, B) (in terms of an unitary equivalence trans-<br>
*> formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues<br>
*> appears in the leading diagonal blocks of the pair (A,B). The leading<br>
*> columns of Q and Z form unitary bases of the corresponding left and<br>
*> right eigenspaces (deflating subspaces). (A, B) must be in<br>
*> generalized Schur canonical form, that is, A and B are both upper<br>
*> triangular.<br>
*><br>
*> ZTGSEN also computes the generalized eigenvalues<br>
*><br>
*>          w(j)= ALPHA(j) / BETA(j)<br>
*><br>
*> of the reordered matrix pair (A, B).<br>
*><br>
*> Optionally, the routine computes estimates of reciprocal condition<br>
*> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),<br>
*> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)<br>
*> between the matrix pairs (A11, B11) and (A22,B22) that correspond to<br>
*> the selected cluster and the eigenvalues outside the cluster, resp.,<br>
*> and norms of "projections" onto left and right eigenspaces w.r.t.<br>
*> the selected cluster in the (1,1)-block.<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] IJOB<br>
*> \verbatim<br>
*>          IJOB is integer<br>
*>          Specifies whether condition numbers are required for the<br>
*>          cluster of eigenvalues (PL and PR) or the deflating subspaces<br>
*>          (Difu and Difl):<br>
*>           =0: Only reorder w.r.t. SELECT. No extras.<br>
*>           =1: Reciprocal of norms of "projections" onto left and right<br>
*>               eigenspaces w.r.t. the selected cluster (PL and PR).<br>
*>           =2: Upper bounds on Difu and Difl. F-norm-based estimate<br>
*>               (DIF(1:2)).<br>
*>           =3: Estimate of Difu and Difl. 1-norm-based estimate<br>
*>               (DIF(1:2)).<br>
*>               About 5 times as expensive as IJOB = 2.<br>
*>           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic<br>
*>               version to get it all.<br>
*>           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above)<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTQ<br>
*> \verbatim<br>
*>          WANTQ is LOGICAL<br>
*>          .TRUE. : update the left transformation matrix Q;<br>
*>          .FALSE.: do not update Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WANTZ<br>
*> \verbatim<br>
*>          WANTZ is LOGICAL<br>
*>          .TRUE. : update the right transformation matrix Z;<br>
*>          .FALSE.: do not update Z.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          SELECT specifies the eigenvalues in the selected cluster. To<br>
*>          select an eigenvalue w(j), SELECT(j) must be set to<br>
*>          .TRUE..<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension(LDA,N)<br>
*>          On entry, the upper triangular matrix A, in generalized<br>
*>          Schur canonical form.<br>
*>          On exit, A is overwritten by the reordered matrix A.<br>
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
*>          B is COMPLEX*16 array, dimension(LDB,N)<br>
*>          On entry, the upper triangular matrix B, in generalized<br>
*>          Schur canonical form.<br>
*>          On exit, B is overwritten by the reordered matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16 array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16 array, dimension (N)<br>
*><br>
*>          The diagonal elements of A and B, respectively,<br>
*>          when the pair (A,B) has been reduced to generalized Schur<br>
*>          form.  ALPHA(i)/BETA(i) i=1,...,N are the generalized<br>
*>          eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX*16 array, dimension (LDQ,N)<br>
*>          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix.<br>
*>          On exit, Q has been postmultiplied by the left unitary<br>
*>          transformation matrix which reorder (A, B); The leading M<br>
*>          columns of Q form orthonormal bases for the specified pair of<br>
*>          left eigenspaces (deflating subspaces).<br>
*>          If WANTQ = .FALSE., Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q. LDQ >= 1.<br>
*>          If WANTQ = .TRUE., LDQ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ,N)<br>
*>          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix.<br>
*>          On exit, Z has been postmultiplied by the left unitary<br>
*>          transformation matrix which reorder (A, B); The leading M<br>
*>          columns of Z form orthonormal bases for the specified pair of<br>
*>          left eigenspaces (deflating subspaces).<br>
*>          If WANTZ = .FALSE., Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z. LDZ >= 1.<br>
*>          If WANTZ = .TRUE., LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The dimension of the specified pair of left and right<br>
*>          eigenspaces, (deflating subspaces) 0 <= M <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PL<br>
*> \verbatim<br>
*>          PL is DOUBLE PRECISION<br>
*> \endverbatim<br>
*><br>
*> \param[out] PR<br>
*> \verbatim<br>
*>          PR is DOUBLE PRECISION<br>
*><br>
*>          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the<br>
*>          reciprocal  of the norm of "projections" onto left and right<br>
*>          eigenspace with respect to the selected cluster.<br>
*>          0 < PL, PR <= 1.<br>
*>          If M = 0 or M = N, PL = PR  = 1.<br>
*>          If IJOB = 0, 2 or 3 PL, PR are not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIF<br>
*> \verbatim<br>
*>          DIF is DOUBLE PRECISION array, dimension (2).<br>
*>          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.<br>
*>          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on<br>
*>          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based<br>
*>          estimates of Difu and Difl, computed using reversed<br>
*>          communication with ZLACN2.<br>
*>          If M = 0 or N, DIF(1:2) = F-norm([A, B]).<br>
*>          If IJOB = 0 or 1, DIF is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >=  1<br>
*>          If IJOB = 1, 2 or 4, LWORK >=  2*M*(N-M)<br>
*>          If IJOB = 3 or 5, LWORK >=  4*M*(N-M)<br>
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
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK. LIWORK >= 1.<br>
*>          If IJOB = 1, 2 or 4, LIWORK >=  N+2;<br>
*>          If IJOB = 3 or 5, LIWORK >= MAX(N+2, 2*M*(N-M));<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal size of the IWORK array,<br>
*>          returns this value as the first entry of the IWORK array, and<br>
*>          no error message related to LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>            =0: Successful exit.<br>
*>            <0: If INFO = -i, the i-th argument had an illegal value.<br>
*>            =1: Reordering of (A, B) failed because the transformed<br>
*>                matrix pair (A, B) would be too far from generalized<br>
*>                Schur form; the problem is very ill-conditioned.<br>
*>                (A, B) may have been partially reordered.<br>
*>                If requested, 0 is returned in DIF(*), PL and PR.<br>
*> \endverbatim<br>
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
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  ZTGSEN first collects the selected eigenvalues by computing unitary<br>
*>  U and W that move them to the top left corner of (A, B). In other<br>
*>  words, the selected eigenvalues are the eigenvalues of (A11, B11) in<br>
*><br>
*>              U**H*(A, B)*W = (A11 A12) (B11 B12) n1<br>
*>                              ( 0  A22),( 0  B22) n2<br>
*>                                n1  n2    n1  n2<br>
*><br>
*>  where N = n1+n2 and U**H means the conjugate transpose of U. The first<br>
*>  n1 columns of U and W span the specified pair of left and right<br>
*>  eigenspaces (deflating subspaces) of (A, B).<br>
*><br>
*>  If (A, B) has been obtained from the generalized real Schur<br>
*>  decomposition of a matrix pair (C, D) = Q*(A, B)*Z**H, then the<br>
*>  reordered generalized Schur form of (C, D) is given by<br>
*><br>
*>           (C, D) = (Q*U)*(U**H *(A, B)*W)*(Z*W)**H,<br>
*><br>
*>  and the first n1 columns of Q*U and Z*W span the corresponding<br>
*>  deflating subspaces of (C, D) (Q and Z store Q*U and Z*W, resp.).<br>
*><br>
*>  Note that if the selected eigenvalue is sufficiently ill-conditioned,<br>
*>  then its value may differ significantly from its value before<br>
*>  reordering.<br>
*><br>
*>  The reciprocal condition numbers of the left and right eigenspaces<br>
*>  spanned by the first n1 columns of U and W (or Q*U and Z*W) may<br>
*>  be returned in DIF(1:2), corresponding to Difu and Difl, resp.<br>
*><br>
*>  The Difu and Difl are defined as:<br>
*><br>
*>       Difu[(A11, B11), (A22, B22)] = sigma-min( Zu )<br>
*>  and<br>
*>       Difl[(A11, B11), (A22, B22)] = Difu[(A22, B22), (A11, B11)],<br>
*><br>
*>  where sigma-min(Zu) is the smallest singular value of the<br>
*>  (2*n1*n2)-by-(2*n1*n2) matrix<br>
*><br>
*>       Zu = [ kron(In2, A11)  -kron(A22**H, In1) ]<br>
*>            [ kron(In2, B11)  -kron(B22**H, In1) ].<br>
*><br>
*>  Here, Inx is the identity matrix of size nx and A22**H is the<br>
*>  conjugate transpose of A22. kron(X, Y) is the Kronecker product between<br>
*>  the matrices X and Y.<br>
*><br>
*>  When DIF(2) is small, small changes in (A, B) can cause large changes<br>
*>  in the deflating subspace. An approximate (asymptotic) bound on the<br>
*>  maximum angular error in the computed deflating subspaces is<br>
*><br>
*>       EPS * norm((A, B)) / DIF(2),<br>
*><br>
*>  where EPS is the machine precision.<br>
*><br>
*>  The reciprocal norm of the projectors on the left and right<br>
*>  eigenspaces associated with (A11, B11) may be returned in PL and PR.<br>
*>  They are computed as follows. First we compute L and R so that<br>
*>  P*(A, B)*Q is block diagonal, where<br>
*><br>
*>       P = ( I -L ) n1           Q = ( I R ) n1<br>
*>           ( 0  I ) n2    and        ( 0 I ) n2<br>
*>             n1 n2                    n1 n2<br>
*><br>
*>  and (L, R) is the solution to the generalized Sylvester equation<br>
*><br>
*>       A11*R - L*A22 = -A12<br>
*>       B11*R - L*B22 = -B12<br>
*><br>
*>  Then PL = (F-norm(L)**2+1)**(-1/2) and PR = (F-norm(R)**2+1)**(-1/2).<br>
*>  An approximate (asymptotic) bound on the average absolute error of<br>
*>  the selected eigenvalues is<br>
*><br>
*>       EPS * norm((A, B)) / PL.<br>
*><br>
*>  There are also global error bounds which valid for perturbations up<br>
*>  to a certain restriction:  A lower bound (x) on the smallest<br>
*>  F-norm(E,F) for which an eigenvalue of (A11, B11) may move and<br>
*>  coalesce with an eigenvalue of (A22, B22) under perturbation (E,F),<br>
*>  (i.e. (A + E, B + F), is<br>
*><br>
*>   x = min(Difu,Difl)/((1/(PL*PL)+1/(PR*PR))**(1/2)+2*max(1/PL,1/PR)).<br>
*><br>
*>  An approximate bound on x can be computed from DIF(1:2), PL and PR.<br>
*><br>
*>  If y = ( F-norm(E,F) / x) <= 1, the angles between the perturbed<br>
*>  (L', R') and unperturbed (L, R) left and right deflating subspaces<br>
*>  associated with the selected cluster in the (1,1)-blocks can be<br>
*>  bounded as<br>
*><br>
*>   max-angle(L, L') <= arctan( y * PL / (1 - y * (1 - PL * PL)**(1/2))<br>
*>   max-angle(R, R') <= arctan( y * PR / (1 - y * (1 - PR * PR)**(1/2))<br>
*><br>
*>  See LAPACK User's Guide section 4.11 or the following references<br>
*>  for more information.<br>
*><br>
*>  Note that if the default method for computing the Frobenius-norm-<br>
*>  based estimate DIF is not wanted (see ZLATDF), then the parameter<br>
*>  IDIFJB (see below) should be changed from 3 to 4 (routine ZLATDF<br>
*>  (IJOB = 2 will be used)). See ZTGSYL for more details.<br>
*> \endverbatim<br>
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
*>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the<br>
*>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in<br>
*>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and<br>
*>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.<br>
*> \n<br>
*>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified<br>
*>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition<br>
*>      Estimation: Theory, Algorithms and Software, Report<br>
*>      UMINF - 94.04, Department of Computing Science, Umea University,<br>
*>      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.<br>
*>      To appear in Numerical Algorithms, 1996.<br>
*> \n<br>
*>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software<br>
*>      for Solving the Generalized Sylvester Equation and Estimating the<br>
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,<br>
*>      Department of Computing Science, Umea University, S-901 87 Umea,<br>
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK working<br>
*>      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,<br>
*>      1996.<br>
*><br>
*  =====================================================================<br>
*/
	public void ztgsen_(INTEGER IJOB,LOGICAL WANTQ,LOGICAL WANTZ,boolean[] SELECT,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHA,double[] BETA,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,INTEGER M,DOUBLE PL,DOUBLE PR,double[] DIF,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b ZTGSJA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTGSJA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsja.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsja.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsja.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,<br>
*                          LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,<br>
*                          Q, LDQ, WORK, NCYCLE, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBQ, JOBU, JOBV<br>
*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N,<br>
*      $                   NCYCLE, P<br>
*       DOUBLE PRECISION   TOLA, TOLB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   ALPHA( * ), BETA( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),<br>
*      $                   U( LDU, * ), V( LDV, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTGSJA computes the generalized singular value decomposition (GSVD)<br>
*> of two complex upper triangular (or trapezoidal) matrices A and B.<br>
*><br>
*> On entry, it is assumed that matrices A and B have the following<br>
*> forms, which may be obtained by the preprocessing subroutine ZGGSVP<br>
*> from a general M-by-N matrix A and P-by-N matrix B:<br>
*><br>
*>              N-K-L  K    L<br>
*>    A =    K ( 0    A12  A13 ) if M-K-L >= 0;<br>
*>           L ( 0     0   A23 )<br>
*>       M-K-L ( 0     0    0  )<br>
*><br>
*>            N-K-L  K    L<br>
*>    A =  K ( 0    A12  A13 ) if M-K-L < 0;<br>
*>       M-K ( 0     0   A23 )<br>
*><br>
*>            N-K-L  K    L<br>
*>    B =  L ( 0     0   B13 )<br>
*>       P-L ( 0     0    0  )<br>
*><br>
*> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular<br>
*> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,<br>
*> otherwise A23 is (M-K)-by-L upper trapezoidal.<br>
*><br>
*> On exit,<br>
*><br>
*>        U**H *A*Q = D1*( 0 R ),    V**H *B*Q = D2*( 0 R ),<br>
*><br>
*> where U, V and Q are unitary matrices.<br>
*> R is a nonsingular upper triangular matrix, and D1<br>
*> and D2 are ``diagonal'' matrices, which are of the following<br>
*> structures:<br>
*><br>
*> If M-K-L >= 0,<br>
*><br>
*>                     K  L<br>
*>        D1 =     K ( I  0 )<br>
*>                 L ( 0  C )<br>
*>             M-K-L ( 0  0 )<br>
*><br>
*>                    K  L<br>
*>        D2 = L   ( 0  S )<br>
*>             P-L ( 0  0 )<br>
*><br>
*>                N-K-L  K    L<br>
*>   ( 0 R ) = K (  0   R11  R12 ) K<br>
*>             L (  0    0   R22 ) L<br>
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
*>                K M-K K+L-M<br>
*>     D1 =   K ( I  0    0   )<br>
*>          M-K ( 0  C    0   )<br>
*><br>
*>                  K M-K K+L-M<br>
*>     D2 =   M-K ( 0  S    0   )<br>
*>          K+L-M ( 0  0    I   )<br>
*>            P-L ( 0  0    0   )<br>
*><br>
*>                N-K-L  K   M-K  K+L-M<br>
*> ( 0 R ) =    K ( 0    R11  R12  R13  )<br>
*>           M-K ( 0     0   R22  R23  )<br>
*>         K+L-M ( 0     0    0   R33  )<br>
*><br>
*> where<br>
*> C = diag( ALPHA(K+1), ... , ALPHA(M) ),<br>
*> S = diag( BETA(K+1),  ... , BETA(M) ),<br>
*> C**2 + S**2 = I.<br>
*><br>
*> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored<br>
*>     (  0  R22 R23 )<br>
*> in B(M-K+1:L,N+M-K-L+1:N) on exit.<br>
*><br>
*> The computation of the unitary transformation matrices U, V or Q<br>
*> is optional.  These matrices may either be formed explicitly, or they<br>
*> may be postmultiplied into input matrices U1, V1, or Q1.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBU<br>
*> \verbatim<br>
*>          JOBU is CHARACTER*1<br>
*>          = 'U':  U must contain a unitary matrix U1 on entry, and<br>
*>                  the product U1*U is returned;<br>
*>          = 'I':  U is initialized to the unit matrix, and the<br>
*>                  unitary matrix U is returned;<br>
*>          = 'N':  U is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV<br>
*> \verbatim<br>
*>          JOBV is CHARACTER*1<br>
*>          = 'V':  V must contain a unitary matrix V1 on entry, and<br>
*>                  the product V1*V is returned;<br>
*>          = 'I':  V is initialized to the unit matrix, and the<br>
*>                  unitary matrix V is returned;<br>
*>          = 'N':  V is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBQ<br>
*> \verbatim<br>
*>          JOBQ is CHARACTER*1<br>
*>          = 'Q':  Q must contain a unitary matrix Q1 on entry, and<br>
*>                  the product Q1*Q is returned;<br>
*>          = 'I':  Q is initialized to the unit matrix, and the<br>
*>                  unitary matrix Q is returned;<br>
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
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*><br>
*>          K and L specify the subblocks in the input matrices A and B:<br>
*>          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,,N-L+1:N)<br>
*>          of A and B, whose GSVD is going to be computed by ZTGSJA.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the M-by-N matrix A.<br>
*>          On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular<br>
*>          matrix R or part of R.  See Purpose for details.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>          On entry, the P-by-N matrix B.<br>
*>          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains<br>
*>          a part of R.  See Purpose for details.<br>
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
*>          TOLA and TOLB are the convergence criteria for the Jacobi-<br>
*>          Kogbetliantz iteration procedure. Generally, they are the<br>
*>          same as used in the preprocessing step, say<br>
*>              TOLA = MAX(M,N)*norm(A)*MAZHEPS,<br>
*>              TOLB = MAX(P,N)*norm(B)*MAZHEPS.<br>
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
*>            ALPHA(K+1:K+L) = diag(C),<br>
*>            BETA(K+1:K+L)  = diag(S),<br>
*>          or if M-K-L < 0,<br>
*>            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0<br>
*>            BETA(K+1:M) = S, BETA(M+1:K+L) = 1.<br>
*>          Furthermore, if K+L < N,<br>
*>            ALPHA(K+L+1:N) = 0 and<br>
*>            BETA(K+L+1:N)  = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] U<br>
*> \verbatim<br>
*>          U is COMPLEX*16 array, dimension (LDU,M)<br>
*>          On entry, if JOBU = 'U', U must contain a matrix U1 (usually<br>
*>          the unitary matrix returned by ZGGSVP).<br>
*>          On exit,<br>
*>          if JOBU = 'I', U contains the unitary matrix U;<br>
*>          if JOBU = 'U', U contains the product U1*U.<br>
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
*> \param[in,out] V<br>
*> \verbatim<br>
*>          V is COMPLEX*16 array, dimension (LDV,P)<br>
*>          On entry, if JOBV = 'V', V must contain a matrix V1 (usually<br>
*>          the unitary matrix returned by ZGGSVP).<br>
*>          On exit,<br>
*>          if JOBV = 'I', V contains the unitary matrix V;<br>
*>          if JOBV = 'V', V contains the product V1*V.<br>
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
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX*16 array, dimension (LDQ,N)<br>
*>          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually<br>
*>          the unitary matrix returned by ZGGSVP).<br>
*>          On exit,<br>
*>          if JOBQ = 'I', Q contains the unitary matrix Q;<br>
*>          if JOBQ = 'Q', Q contains the product Q1*Q.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] NCYCLE<br>
*> \verbatim<br>
*>          NCYCLE is INTEGER<br>
*>          The number of cycles required for convergence.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          = 1:  the procedure does not converge after MAXIT cycles.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  MAXIT   INTEGER<br>
*>          MAXIT specifies the total loops that the iterative procedure<br>
*>          may take. If after MAXIT cycles, the routine fails to<br>
*>          converge, we return INFO = 1.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  ZTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce<br>
*>  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L<br>
*>  matrix B13 to the form:<br>
*><br>
*>           U1**H *A13*Q1 = C1*R1; V1**H *B13*Q1 = S1*R1,<br>
*><br>
*>  where U1, V1 and Q1 are unitary matrix.<br>
*>  C1 and S1 are diagonal matrices satisfying<br>
*><br>
*>                C1**2 + S1**2 = I,<br>
*><br>
*>  and R1 is an L-by-L nonsingular upper triangular matrix.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztgsja_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER P,INTEGER N,INTEGER K,INTEGER L,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE TOLA,DOUBLE TOLB,double[] ALPHA,double[] BETA,double[] U,INTEGER LDU,double[] V,INTEGER LDV,double[] Q,INTEGER LDQ,double[] WORK,INTEGER NCYCLE,INTEGER INFO);
/**
*> \brief \b ZTGSNA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTGSNA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsna.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsna.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsna.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,<br>
*                          LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, JOB<br>
*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   DIF( * ), S( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), VL( LDVL, * ),<br>
*      $                   VR( LDVR, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTGSNA estimates reciprocal condition numbers for specified<br>
*> eigenvalues and/or eigenvectors of a matrix pair (A, B).<br>
*><br>
*> (A, B) must be in generalized Schur canonical form, that is, A and<br>
*> B are both upper triangular.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          Specifies whether condition numbers are required for<br>
*>          eigenvalues (S) or eigenvectors (DIF):<br>
*>          = 'E': for eigenvalues only (S);<br>
*>          = 'V': for eigenvectors only (DIF);<br>
*>          = 'B': for both eigenvalues and eigenvectors (S and DIF).<br>
*> \endverbatim<br>
*><br>
*> \param[in] HOWMNY<br>
*> \verbatim<br>
*>          HOWMNY is CHARACTER*1<br>
*>          = 'A': compute condition numbers for all eigenpairs;<br>
*>          = 'S': compute condition numbers for selected eigenpairs<br>
*>                 specified by the array SELECT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          If HOWMNY = 'S', SELECT specifies the eigenpairs for which<br>
*>          condition numbers are required. To select condition numbers<br>
*>          for the corresponding j-th eigenvalue and/or eigenvector,<br>
*>          SELECT(j) must be set to .TRUE..<br>
*>          If HOWMNY = 'A', SELECT is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the square matrix pair (A, B). N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The upper triangular matrix A in the pair (A,B).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>          The upper triangular matrix B in the pair (A, B).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is COMPLEX*16 array, dimension (LDVL,M)<br>
*>          IF JOB = 'E' or 'B', VL must contain left eigenvectors of<br>
*>          (A, B), corresponding to the eigenpairs specified by HOWMNY<br>
*>          and SELECT.  The eigenvectors must be stored in consecutive<br>
*>          columns of VL, as returned by ZTGEVC.<br>
*>          If JOB = 'V', VL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL. LDVL >= 1; and<br>
*>          If JOB = 'E' or 'B', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VR<br>
*> \verbatim<br>
*>          VR is COMPLEX*16 array, dimension (LDVR,M)<br>
*>          IF JOB = 'E' or 'B', VR must contain right eigenvectors of<br>
*>          (A, B), corresponding to the eigenpairs specified by HOWMNY<br>
*>          and SELECT.  The eigenvectors must be stored in consecutive<br>
*>          columns of VR, as returned by ZTGEVC.<br>
*>          If JOB = 'V', VR is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR. LDVR >= 1;<br>
*>          If JOB = 'E' or 'B', LDVR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (MM)<br>
*>          If JOB = 'E' or 'B', the reciprocal condition numbers of the<br>
*>          selected eigenvalues, stored in consecutive elements of the<br>
*>          array.<br>
*>          If JOB = 'V', S is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIF<br>
*> \verbatim<br>
*>          DIF is DOUBLE PRECISION array, dimension (MM)<br>
*>          If JOB = 'V' or 'B', the estimated reciprocal condition<br>
*>          numbers of the selected eigenvectors, stored in consecutive<br>
*>          elements of the array.<br>
*>          If the eigenvalues cannot be reordered to compute DIF(j),<br>
*>          DIF(j) is set to 0; this can only occur when the true value<br>
*>          would be very small anyway.<br>
*>          For each eigenvalue/vector specified by SELECT, DIF stores<br>
*>          a Frobenius norm-based estimate of Difl.<br>
*>          If JOB = 'E', DIF is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MM<br>
*> \verbatim<br>
*>          MM is INTEGER<br>
*>          The number of elements in the arrays S and DIF. MM >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of elements of the arrays S and DIF used to store<br>
*>          the specified condition numbers; for each selected eigenvalue<br>
*>          one element is used. If HOWMNY = 'A', M is set to N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,N).<br>
*>          If JOB = 'V' or 'B', LWORK >= max(1,2*N*N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N+2)<br>
*>          If JOB = 'E', IWORK is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: Successful exit<br>
*>          < 0: If INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The reciprocal of the condition number of the i-th generalized<br>
*>  eigenvalue w = (a, b) is defined as<br>
*><br>
*>          S(I) = (|v**HAu|**2 + |v**HBu|**2)**(1/2) / (norm(u)*norm(v))<br>
*><br>
*>  where u and v are the right and left eigenvectors of (A, B)<br>
*>  corresponding to w; |z| denotes the absolute value of the complex<br>
*>  number, and norm(u) denotes the 2-norm of the vector u. The pair<br>
*>  (a, b) corresponds to an eigenvalue w = a/b (= v**HAu/v**HBu) of the<br>
*>  matrix pair (A, B). If both a and b equal zero, then (A,B) is<br>
*>  singular and S(I) = -1 is returned.<br>
*><br>
*>  An approximate error bound on the chordal distance between the i-th<br>
*>  computed generalized eigenvalue w and the corresponding exact<br>
*>  eigenvalue lambda is<br>
*><br>
*>          chord(w, lambda) <=   EPS * norm(A, B) / S(I),<br>
*><br>
*>  where EPS is the machine precision.<br>
*><br>
*>  The reciprocal of the condition number of the right eigenvector u<br>
*>  and left eigenvector v corresponding to the generalized eigenvalue w<br>
*>  is defined as follows. Suppose<br>
*><br>
*>                   (A, B) = ( a   *  ) ( b  *  )  1<br>
*>                            ( 0  A22 ),( 0 B22 )  n-1<br>
*>                              1  n-1     1 n-1<br>
*><br>
*>  Then the reciprocal condition number DIF(I) is<br>
*><br>
*>          Difl[(a, b), (A22, B22)]  = sigma-min( Zl )<br>
*><br>
*>  where sigma-min(Zl) denotes the smallest singular value of<br>
*><br>
*>         Zl = [ kron(a, In-1) -kron(1, A22) ]<br>
*>              [ kron(b, In-1) -kron(1, B22) ].<br>
*><br>
*>  Here In-1 is the identity matrix of size n-1 and X**H is the conjugate<br>
*>  transpose of X. kron(X, Y) is the Kronecker product between the<br>
*>  matrices X and Y.<br>
*><br>
*>  We approximate the smallest singular value of Zl with an upper<br>
*>  bound. This is done by ZLATDF.<br>
*><br>
*>  An approximate error bound for a computed eigenvector VL(i) or<br>
*>  VR(i) is given by<br>
*><br>
*>                      EPS * norm(A, B) / DIF(i).<br>
*><br>
*>  See ref. [2-3] for more details and further references.<br>
*> \endverbatim<br>
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
*>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the<br>
*>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in<br>
*>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and<br>
*>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.<br>
*><br>
*>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified<br>
*>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition<br>
*>      Estimation: Theory, Algorithms and Software, Report<br>
*>      UMINF - 94.04, Department of Computing Science, Umea University,<br>
*>      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.<br>
*>      To appear in Numerical Algorithms, 1996.<br>
*><br>
*>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software<br>
*>      for Solving the Generalized Sylvester Equation and Estimating the<br>
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,<br>
*>      Department of Computing Science, Umea University, S-901 87 Umea,<br>
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working<br>
*>      Note 75.<br>
*>      To appear in ACM Trans. on Math. Software, Vol 22, No 1, 1996.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztgsna_(CHARACTER JOB,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] S,double[] DIF,INTEGER MM,INTEGER M,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b ZTGSY2 solves the generalized Sylvester equation (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTGSY2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsy2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsy2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsy2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,<br>
*                          LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N<br>
*       DOUBLE PRECISION   RDSCAL, RDSUM, SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),<br>
*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTGSY2 solves the generalized Sylvester equation<br>
*><br>
*>             A * R - L * B = scale * C               (1)<br>
*>             D * R - L * E = scale * F<br>
*><br>
*> using Level 1 and 2 BLAS, where R and L are unknown M-by-N matrices,<br>
*> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,<br>
*> N-by-N and M-by-N, respectively. A, B, D and E are upper triangular<br>
*> (i.e., (A,D) and (B,E) in generalized Schur form).<br>
*><br>
*> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output<br>
*> scaling factor chosen to avoid overflow.<br>
*><br>
*> In matrix notation solving equation (1) corresponds to solve<br>
*> Zx = scale * b, where Z is defined as<br>
*><br>
*>        Z = [ kron(In, A)  -kron(B**H, Im) ]             (2)<br>
*>            [ kron(In, D)  -kron(E**H, Im) ],<br>
*><br>
*> Ik is the identity matrix of size k and X**H is the conjuguate transpose of X.<br>
*> kron(X, Y) is the Kronecker product between the matrices X and Y.<br>
*><br>
*> If TRANS = 'C', y in the conjugate transposed system Z**H*y = scale*b<br>
*> is solved for, which is equivalent to solve for R and L in<br>
*><br>
*>             A**H * R  + D**H * L   = scale * C           (3)<br>
*>             R  * B**H + L  * E**H  = scale * -F<br>
*><br>
*> This case is used to compute an estimate of Dif[(A, D), (B, E)] =<br>
*> = sigma_min(Z) using reverse communicaton with ZLACON.<br>
*><br>
*> ZTGSY2 also (IJOB >= 1) contributes to the computation in ZTGSYL<br>
*> of an upper bound on the separation between to matrix pairs. Then<br>
*> the input (A, D), (B, E) are sub-pencils of two matrix pairs in<br>
*> ZTGSYL.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N', solve the generalized Sylvester equation (1).<br>
*>          = 'T': solve the 'transposed' system (3).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IJOB<br>
*> \verbatim<br>
*>          IJOB is INTEGER<br>
*>          Specifies what kind of functionality to be performed.<br>
*>          =0: solve (1) only.<br>
*>          =1: A contribution from this subsystem to a Frobenius<br>
*>              norm-based estimate of the separation between two matrix<br>
*>              pairs is computed. (look ahead strategy is used).<br>
*>          =2: A contribution from this subsystem to a Frobenius<br>
*>              norm-based estimate of the separation between two matrix<br>
*>              pairs is computed. (DGECON on sub-systems is used.)<br>
*>          Not referenced if TRANS = 'T'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          On entry, M specifies the order of A and D, and the row<br>
*>          dimension of C, F, R and L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          On entry, N specifies the order of B and E, and the column<br>
*>          dimension of C, F, R and L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA, M)<br>
*>          On entry, A contains an upper triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the matrix A. LDA >= max(1, M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB, N)<br>
*>          On entry, B contains an upper triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the matrix B. LDB >= max(1, N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array, dimension (LDC, N)<br>
*>          On entry, C contains the right-hand-side of the first matrix<br>
*>          equation in (1).<br>
*>          On exit, if IJOB = 0, C has been overwritten by the solution<br>
*>          R.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the matrix C. LDC >= max(1, M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is COMPLEX*16 array, dimension (LDD, M)<br>
*>          On entry, D contains an upper triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDD<br>
*> \verbatim<br>
*>          LDD is INTEGER<br>
*>          The leading dimension of the matrix D. LDD >= max(1, M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is COMPLEX*16 array, dimension (LDE, N)<br>
*>          On entry, E contains an upper triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDE<br>
*> \verbatim<br>
*>          LDE is INTEGER<br>
*>          The leading dimension of the matrix E. LDE >= max(1, N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] F<br>
*> \verbatim<br>
*>          F is COMPLEX*16 array, dimension (LDF, N)<br>
*>          On entry, F contains the right-hand-side of the second matrix<br>
*>          equation in (1).<br>
*>          On exit, if IJOB = 0, F has been overwritten by the solution<br>
*>          L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDF<br>
*> \verbatim<br>
*>          LDF is INTEGER<br>
*>          The leading dimension of the matrix F. LDF >= max(1, M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is DOUBLE PRECISION<br>
*>          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions<br>
*>          R and L (C and F on entry) will hold the solutions to a<br>
*>          slightly perturbed system but the input matrices A, B, D and<br>
*>          E have not been changed. If SCALE = 0, R and L will hold the<br>
*>          solutions to the homogeneous system with C = F = 0.<br>
*>          Normally, SCALE = 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RDSUM<br>
*> \verbatim<br>
*>          RDSUM is DOUBLE PRECISION<br>
*>          On entry, the sum of squares of computed contributions to<br>
*>          the Dif-estimate under computation by ZTGSYL, where the<br>
*>          scaling factor RDSCAL (see below) has been factored out.<br>
*>          On exit, the corresponding sum of squares updated with the<br>
*>          contributions from the current sub-system.<br>
*>          If TRANS = 'T' RDSUM is not touched.<br>
*>          NOTE: RDSUM only makes sense when ZTGSY2 is called by<br>
*>          ZTGSYL.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] RDSCAL<br>
*> \verbatim<br>
*>          RDSCAL is DOUBLE PRECISION<br>
*>          On entry, scaling factor used to prevent overflow in RDSUM.<br>
*>          On exit, RDSCAL is updated w.r.t. the current contributions<br>
*>          in RDSUM.<br>
*>          If TRANS = 'T', RDSCAL is not touched.<br>
*>          NOTE: RDSCAL only makes sense when ZTGSY2 is called by<br>
*>          ZTGSYL.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          On exit, if INFO is set to<br>
*>            =0: Successful exit<br>
*>            <0: If INFO = -i, input argument number i is illegal.<br>
*>            >0: The matrix pairs (A, D) and (B, E) have common or very<br>
*>                close eigenvalues.<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,<br>
*>     Umea University, S-901 87 Umea, Sweden.<br>
*<br>
*  =====================================================================<br>
*/
	public void ztgsy2_(CHARACTER TRANS,INTEGER IJOB,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] C,INTEGER LDC,double[] D,INTEGER LDD,double[] E,INTEGER LDE,double[] F,INTEGER LDF,DOUBLE SCALE,DOUBLE RDSUM,DOUBLE RDSCAL,INTEGER INFO);
/**
*> \brief \b ZTGSYL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTGSYL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsyl.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsyl.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsyl.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,<br>
*                          LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF,<br>
*      $                   LWORK, M, N<br>
*       DOUBLE PRECISION   DIF, SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),<br>
*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTGSYL solves the generalized Sylvester equation:<br>
*><br>
*>             A * R - L * B = scale * C            (1)<br>
*>             D * R - L * E = scale * F<br>
*><br>
*> where R and L are unknown m-by-n matrices, (A, D), (B, E) and<br>
*> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,<br>
*> respectively, with complex entries. A, B, D and E are upper<br>
*> triangular (i.e., (A,D) and (B,E) in generalized Schur form).<br>
*><br>
*> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1<br>
*> is an output scaling factor chosen to avoid overflow.<br>
*><br>
*> In matrix notation (1) is equivalent to solve Zx = scale*b, where Z<br>
*> is defined as<br>
*><br>
*>        Z = [ kron(In, A)  -kron(B**H, Im) ]        (2)<br>
*>            [ kron(In, D)  -kron(E**H, Im) ],<br>
*><br>
*> Here Ix is the identity matrix of size x and X**H is the conjugate<br>
*> transpose of X. Kron(X, Y) is the Kronecker product between the<br>
*> matrices X and Y.<br>
*><br>
*> If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b<br>
*> is solved for, which is equivalent to solve for R and L in<br>
*><br>
*>             A**H * R + D**H * L = scale * C           (3)<br>
*>             R * B**H + L * E**H = scale * -F<br>
*><br>
*> This case (TRANS = 'C') is used to compute an one-norm-based estimate<br>
*> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)<br>
*> and (B,E), using ZLACON.<br>
*><br>
*> If IJOB >= 1, ZTGSYL computes a Frobenius norm-based estimate of<br>
*> Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the<br>
*> reciprocal of the smallest singular value of Z.<br>
*><br>
*> This is a level-3 BLAS algorithm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': solve the generalized sylvester equation (1).<br>
*>          = 'C': solve the "conjugate transposed" system (3).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IJOB<br>
*> \verbatim<br>
*>          IJOB is INTEGER<br>
*>          Specifies what kind of functionality to be performed.<br>
*>          =0: solve (1) only.<br>
*>          =1: The functionality of 0 and 3.<br>
*>          =2: The functionality of 0 and 4.<br>
*>          =3: Only an estimate of Dif[(A,D), (B,E)] is computed.<br>
*>              (look ahead strategy is used).<br>
*>          =4: Only an estimate of Dif[(A,D), (B,E)] is computed.<br>
*>              (ZGECON on sub-systems is used).<br>
*>          Not referenced if TRANS = 'C'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The order of the matrices A and D, and the row dimension of<br>
*>          the matrices C, F, R and L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices B and E, and the column dimension<br>
*>          of the matrices C, F, R and L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA, M)<br>
*>          The upper triangular matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1, M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB, N)<br>
*>          The upper triangular matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1, N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array, dimension (LDC, N)<br>
*>          On entry, C contains the right-hand-side of the first matrix<br>
*>          equation in (1) or (3).<br>
*>          On exit, if IJOB = 0, 1 or 2, C has been overwritten by<br>
*>          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,<br>
*>          the solution achieved during the computation of the<br>
*>          Dif-estimate.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDC >= max(1, M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is COMPLEX*16 array, dimension (LDD, M)<br>
*>          The upper triangular matrix D.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDD<br>
*> \verbatim<br>
*>          LDD is INTEGER<br>
*>          The leading dimension of the array D. LDD >= max(1, M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is COMPLEX*16 array, dimension (LDE, N)<br>
*>          The upper triangular matrix E.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDE<br>
*> \verbatim<br>
*>          LDE is INTEGER<br>
*>          The leading dimension of the array E. LDE >= max(1, N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] F<br>
*> \verbatim<br>
*>          F is COMPLEX*16 array, dimension (LDF, N)<br>
*>          On entry, F contains the right-hand-side of the second matrix<br>
*>          equation in (1) or (3).<br>
*>          On exit, if IJOB = 0, 1 or 2, F has been overwritten by<br>
*>          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,<br>
*>          the solution achieved during the computation of the<br>
*>          Dif-estimate.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDF<br>
*> \verbatim<br>
*>          LDF is INTEGER<br>
*>          The leading dimension of the array F. LDF >= max(1, M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIF<br>
*> \verbatim<br>
*>          DIF is DOUBLE PRECISION<br>
*>          On exit DIF is the reciprocal of a lower bound of the<br>
*>          reciprocal of the Dif-function, i.e. DIF is an upper bound of<br>
*>          Dif[(A,D), (B,E)] = sigma-min(Z), where Z as in (2).<br>
*>          IF IJOB = 0 or TRANS = 'C', DIF is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is DOUBLE PRECISION<br>
*>          On exit SCALE is the scaling factor in (1) or (3).<br>
*>          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,<br>
*>          to a slightly perturbed system but the input matrices A, B,<br>
*>          D and E have not been changed. If SCALE = 0, R and L will<br>
*>          hold the solutions to the homogenious system with C = F = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK > = 1.<br>
*>          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= max(1,2*M*N).<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (M+N+2)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>            =0: successful exit<br>
*>            <0: If INFO = -i, the i-th argument had an illegal value.<br>
*>            >0: (A, D) and (B, E) have common or very close<br>
*>                eigenvalues.<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYcomputational<br>
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
*>  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software<br>
*>      for Solving the Generalized Sylvester Equation and Estimating the<br>
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,<br>
*>      Department of Computing Science, Umea University, S-901 87 Umea,<br>
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working<br>
*>      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,<br>
*>      No 1, 1996.<br>
*> \n<br>
*>  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester<br>
*>      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.<br>
*>      Appl., 15(4):1045-1060, 1994.<br>
*> \n<br>
*>  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with<br>
*>      Condition Estimators for Solving the Generalized Sylvester<br>
*>      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,<br>
*>      July 1989, pp 745-751.<br>
*><br>
*  =====================================================================<br>
*/
	public void ztgsyl_(CHARACTER TRANS,INTEGER IJOB,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] C,INTEGER LDC,double[] D,INTEGER LDD,double[] E,INTEGER LDE,double[] F,INTEGER LDF,DOUBLE SCALE,DOUBLE DIF,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b ZTPCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, RWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            INFO, N<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPCON estimates the reciprocal of the condition number of a packed<br>
*> triangular matrix A, in either the 1-norm or the infinity-norm.<br>
*><br>
*> The norm of A is computed and an estimate is obtained for<br>
*> norm(inv(A)), then the reciprocal of the condition number is<br>
*> computed as<br>
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
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangular matrix A, packed columnwise in<br>
*>          a linear array.  The j-th column of A is stored in the array<br>
*>          AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          If DIAG = 'U', the diagonal elements of A are not referenced<br>
*>          and are assumed to be 1.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztpcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] AP,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTPMQRT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPMQRT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpmqrt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpmqrt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpmqrt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,<br>
*                           A, LDA, B, LDB, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER SIDE, TRANS<br>
*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), <br>
*      $          WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPMQRT applies a complex orthogonal matrix Q obtained from a <br>
*> "triangular-pentagonal" complex block reflector H to a general<br>
*> complex matrix C, which consists of two blocks A and B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q**H from the Left;<br>
*>          = 'R': apply Q or Q**H from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'C':  Transpose, apply Q**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix B. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix B. N >= 0.<br>
*> \endverbatim<br>
*> <br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The number of elementary reflectors whose product defines<br>
*>          the matrix Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The order of the trapezoidal part of V.  <br>
*>          K >= L >= 0.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The block size used for the storage of T.  K >= NB >= 1.<br>
*>          This must be the same value of NB used to generate T<br>
*>          in CTPQRT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is COMPLEX*16 array, dimension (LDA,K)<br>
*>          The i-th column must contain the vector which defines the<br>
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by<br>
*>          CTPQRT in B.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDV<br>
*> \verbatim<br>
*>          LDV is INTEGER<br>
*>          The leading dimension of the array V.<br>
*>          If SIDE = 'L', LDV >= max(1,M);<br>
*>          if SIDE = 'R', LDV >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] T<br>
*> \verbatim<br>
*>          T is COMPLEX*16 array, dimension (LDT,K)<br>
*>          The upper triangular factors of the block reflectors<br>
*>          as returned by CTPQRT, stored as a NB-by-K matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension<br>
*>          (LDA,N) if SIDE = 'L' or <br>
*>          (LDA,K) if SIDE = 'R'<br>
*>          On entry, the K-by-N or M-by-K matrix A.<br>
*>          On exit, A is overwritten by the corresponding block of <br>
*>          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. <br>
*>          If SIDE = 'L', LDC >= max(1,K);<br>
*>          If SIDE = 'R', LDC >= max(1,M). <br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>          On entry, the M-by-N matrix B.<br>
*>          On exit, B is overwritten by the corresponding block of<br>
*>          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. <br>
*>          LDB >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array. The dimension of WORK is<br>
*>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The columns of the pentagonal matrix V contain the elementary reflectors<br>
*>  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a <br>
*>  trapezoidal block V2:<br>
*><br>
*>        V = [V1]<br>
*>            [V2].<br>
*><br>
*>  The size of the trapezoidal block V2 is determined by the parameter L, <br>
*>  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L<br>
*>  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular;<br>
*>  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.<br>
*><br>
*>  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K. <br>
*>                      [B]   <br>
*>  <br>
*>  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K.<br>
*><br>
*>  The complex orthogonal matrix Q is formed from V and T.<br>
*><br>
*>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.<br>
*><br>
*>  If TRANS='C' and SIDE='L', C is on exit replaced with Q**H * C.<br>
*><br>
*>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.<br>
*><br>
*>  If TRANS='C' and SIDE='R', C is on exit replaced with C * Q**H.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztpmqrt_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER L,INTEGER NB,double[] V,INTEGER LDV,double[] T,INTEGER LDT,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZTPQRT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPQRT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpqrt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpqrt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpqrt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPQRT computes a blocked QR factorization of a complex <br>
*> "triangular-pentagonal" matrix C, which is composed of a <br>
*> triangular block A and pentagonal block B, using the compact <br>
*> WY representation for Q.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix B.  <br>
*>          M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix B, and the order of the<br>
*>          triangular matrix A.<br>
*>          N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The number of rows of the upper trapezoidal part of B.<br>
*>          MIN(M,N) >= L >= 0.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          The block size to be used in the blocked QR.  N >= NB >= 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the upper triangular N-by-N matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the upper triangular matrix R.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>          On entry, the pentagonal M-by-N matrix B.  The first M-L rows <br>
*>          are rectangular, and the last L rows are upper trapezoidal.<br>
*>          On exit, B contains the pentagonal matrix V.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is COMPLEX*16 array, dimension (LDT,N)<br>
*>          The upper triangular block reflectors stored in compact form<br>
*>          as a sequence of upper triangular blocks.  See Further Details.<br>
*> \endverbatim<br>
*>          <br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= NB.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (NB*N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The input matrix C is a (N+M)-by-N matrix  <br>
*><br>
*>               C = [ A ]<br>
*>                   [ B ]        <br>
*><br>
*>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal<br>
*>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N<br>
*>  upper trapezoidal matrix B2:<br>
*><br>
*>               B = [ B1 ]  <- (M-L)-by-N rectangular<br>
*>                   [ B2 ]  <-     L-by-N upper trapezoidal.<br>
*><br>
*>  The upper trapezoidal matrix B2 consists of the first L rows of a<br>
*>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, <br>
*>  B is rectangular M-by-N; if M=L=N, B is upper triangular.  <br>
*><br>
*>  The matrix W stores the elementary reflectors H(i) in the i-th column<br>
*>  below the diagonal (of A) in the (N+M)-by-N input matrix C<br>
*><br>
*>               C = [ A ]  <- upper triangular N-by-N<br>
*>                   [ B ]  <- M-by-N pentagonal<br>
*><br>
*>  so that W can be represented as<br>
*><br>
*>               W = [ I ]  <- identity, N-by-N<br>
*>                   [ V ]  <- M-by-N, same form as B.<br>
*><br>
*>  Thus, all of information needed for W is contained on exit in B, which<br>
*>  we call V above.  Note that V has the same form as B; that is, <br>
*><br>
*>               V = [ V1 ] <- (M-L)-by-N rectangular<br>
*>                   [ V2 ] <-     L-by-N upper trapezoidal.<br>
*><br>
*>  The columns of V represent the vectors which define the H(i)'s.  <br>
*><br>
*>  The number of blocks is B = ceiling(N/NB), where each<br>
*>  block is of order NB except for the last block, which is of order <br>
*>  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular block<br>
*>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB-by-IB <br>
*>  for the last block) T's are stored in the NB-by-N matrix T as<br>
*><br>
*>               T = [T1 T2 ... TB].<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztpqrt_(INTEGER M,INTEGER N,INTEGER L,INTEGER NB,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] T,INTEGER LDT,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZTPQRT2 computes a QR factorization of a real or complex "triangular-pentagonal" matrix, which is composed of a triangular block and a pentagonal block, using the compact WY representation for Q.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPQRT2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpqrt2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpqrt2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpqrt2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER   INFO, LDA, LDB, LDT, N, M, L<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16   A( LDA, * ), B( LDB, * ), T( LDT, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPQRT2 computes a QR factorization of a complex "triangular-pentagonal"<br>
*> matrix C, which is composed of a triangular block A and pentagonal block B, <br>
*> using the compact WY representation for Q.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The total number of rows of the matrix B.  <br>
*>          M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix B, and the order of<br>
*>          the triangular matrix A.<br>
*>          N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The number of rows of the upper trapezoidal part of B.  <br>
*>          MIN(M,N) >= L >= 0.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the upper triangular N-by-N matrix A.<br>
*>          On exit, the elements on and above the diagonal of the array<br>
*>          contain the upper triangular matrix R.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>          On entry, the pentagonal M-by-N matrix B.  The first M-L rows <br>
*>          are rectangular, and the last L rows are upper trapezoidal.<br>
*>          On exit, B contains the pentagonal matrix V.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] T<br>
*> \verbatim<br>
*>          T is COMPLEX*16 array, dimension (LDT,N)<br>
*>          The N-by-N upper triangular factor T of the block reflector.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= max(1,N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The input matrix C is a (N+M)-by-N matrix  <br>
*><br>
*>               C = [ A ]<br>
*>                   [ B ]        <br>
*><br>
*>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal<br>
*>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N<br>
*>  upper trapezoidal matrix B2:<br>
*><br>
*>               B = [ B1 ]  <- (M-L)-by-N rectangular<br>
*>                   [ B2 ]  <-     L-by-N upper trapezoidal.<br>
*><br>
*>  The upper trapezoidal matrix B2 consists of the first L rows of a<br>
*>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0, <br>
*>  B is rectangular M-by-N; if M=L=N, B is upper triangular.  <br>
*><br>
*>  The matrix W stores the elementary reflectors H(i) in the i-th column<br>
*>  below the diagonal (of A) in the (N+M)-by-N input matrix C<br>
*><br>
*>               C = [ A ]  <- upper triangular N-by-N<br>
*>                   [ B ]  <- M-by-N pentagonal<br>
*><br>
*>  so that W can be represented as<br>
*><br>
*>               W = [ I ]  <- identity, N-by-N<br>
*>                   [ V ]  <- M-by-N, same form as B.<br>
*><br>
*>  Thus, all of information needed for W is contained on exit in B, which<br>
*>  we call V above.  Note that V has the same form as B; that is, <br>
*><br>
*>               V = [ V1 ] <- (M-L)-by-N rectangular<br>
*>                   [ V2 ] <-     L-by-N upper trapezoidal.<br>
*><br>
*>  The columns of V represent the vectors which define the H(i)'s.  <br>
*>  The (M+N)-by-(M+N) block reflector H is then given by<br>
*><br>
*>               H = I - W * T * W**H<br>
*><br>
*>  where W**H is the conjugate transpose of W and T is the upper triangular<br>
*>  factor of the block reflector.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztpqrt2_(INTEGER M,INTEGER N,INTEGER L,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] T,INTEGER LDT,INTEGER INFO);
/**
*> \brief \b ZTPRFB applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex matrix, which is composed of two blocks.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPRFB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztprfb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztprfb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztprfb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, <br>
*                          V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER DIRECT, SIDE, STOREV, TRANS<br>
*       INTEGER   K, L, LDA, LDB, LDT, LDV, LDWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16   A( LDA, * ), B( LDB, * ), T( LDT, * ), <br>
*      $          V( LDV, * ), WORK( LDWORK, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPRFB applies a complex "triangular-pentagonal" block reflector H or its <br>
*> conjugate transpose H**H to a complex matrix C, which is composed of two <br>
*> blocks A and B, either from the left or right.<br>
*> <br>
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
*>          = 'C': Columns<br>
*>          = 'R': Rows<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix B.  <br>
*>          M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix B.  <br>
*>          N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>          The order of the matrix T, i.e. the number of elementary<br>
*>          reflectors whose product defines the block reflector.  <br>
*>          K >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] L<br>
*> \verbatim<br>
*>          L is INTEGER<br>
*>          The order of the trapezoidal part of V.  <br>
*>          K >= L >= 0.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] V<br>
*> \verbatim<br>
*>          V is COMPLEX*16 array, dimension<br>
*>                                (LDV,K) if STOREV = 'C'<br>
*>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'<br>
*>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'<br>
*>          The pentagonal matrix V, which contains the elementary reflectors<br>
*>          H(1), H(2), ..., H(K).  See Further Details.<br>
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
*>          T is COMPLEX*16 array, dimension (LDT,K)<br>
*>          The triangular K-by-K matrix T in the representation of the<br>
*>          block reflector.  <br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. <br>
*>          LDT >= K.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension<br>
*>          (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R'<br>
*>          On entry, the K-by-N or M-by-K matrix A.<br>
*>          On exit, A is overwritten by the corresponding block of <br>
*>          H*C or H**H*C or C*H or C*H**H.  See Futher Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. <br>
*>          If SIDE = 'L', LDC >= max(1,K);<br>
*>          If SIDE = 'R', LDC >= max(1,M). <br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>          On entry, the M-by-N matrix B.<br>
*>          On exit, B is overwritten by the corresponding block of<br>
*>          H*C or H**H*C or C*H or C*H**H.  See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. <br>
*>          LDB >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension<br>
*>          (LDWORK,N) if SIDE = 'L',<br>
*>          (LDWORK,K) if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDWORK<br>
*> \verbatim<br>
*>          LDWORK is INTEGER<br>
*>          The leading dimension of the array WORK.<br>
*>          If SIDE = 'L', LDWORK >= K; <br>
*>          if SIDE = 'R', LDWORK >= M.<br>
*> \endverbatim<br>
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
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The matrix C is a composite matrix formed from blocks A and B.<br>
*>  The block B is of size M-by-N; if SIDE = 'R', A is of size M-by-K, <br>
*>  and if SIDE = 'L', A is of size K-by-N.<br>
*><br>
*>  If SIDE = 'R' and DIRECT = 'F', C = [A B].<br>
*><br>
*>  If SIDE = 'L' and DIRECT = 'F', C = [A] <br>
*>                                      [B].<br>
*><br>
*>  If SIDE = 'R' and DIRECT = 'B', C = [B A].<br>
*><br>
*>  If SIDE = 'L' and DIRECT = 'B', C = [B]<br>
*>                                      [A]. <br>
*><br>
*>  The pentagonal matrix V is composed of a rectangular block V1 and a <br>
*>  trapezoidal block V2.  The size of the trapezoidal block is determined by <br>
*>  the parameter L, where 0<=L<=K.  If L=K, the V2 block of V is triangular;<br>
*>  if L=0, there is no trapezoidal block, thus V = V1 is rectangular.<br>
*><br>
*>  If DIRECT = 'F' and STOREV = 'C':  V = [V1]<br>
*>                                         [V2]<br>
*>     - V2 is upper trapezoidal (first L rows of K-by-K upper triangular)<br>
*><br>
*>  If DIRECT = 'F' and STOREV = 'R':  V = [V1 V2]<br>
*><br>
*>     - V2 is lower trapezoidal (first L columns of K-by-K lower triangular)<br>
*><br>
*>  If DIRECT = 'B' and STOREV = 'C':  V = [V2]<br>
*>                                         [V1]<br>
*>     - V2 is lower trapezoidal (last L rows of K-by-K lower triangular)<br>
*><br>
*>  If DIRECT = 'B' and STOREV = 'R':  V = [V2 V1]<br>
*>    <br>
*>     - V2 is upper trapezoidal (last L columns of K-by-K upper triangular)<br>
*><br>
*>  If STOREV = 'C' and SIDE = 'L', V is M-by-K with V2 L-by-K.<br>
*><br>
*>  If STOREV = 'C' and SIDE = 'R', V is N-by-K with V2 L-by-K.<br>
*><br>
*>  If STOREV = 'R' and SIDE = 'L', V is K-by-M with V2 K-by-L.<br>
*><br>
*>  If STOREV = 'R' and SIDE = 'R', V is K-by-N with V2 K-by-L.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztprfb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,INTEGER L,double[] V,INTEGER LDV,double[] T,INTEGER LDT,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] WORK,INTEGER LDWORK);
/**
*> \brief \b ZTPRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztprfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztprfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztprfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX,<br>
*                          FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX*16         AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPRFS provides error bounds and backward error estimates for the<br>
*> solution to a system of linear equations with a triangular packed<br>
*> coefficient matrix.<br>
*><br>
*> The solution matrix X must be computed by ZTPTRS or some other<br>
*> means before entering this routine.  ZTPRFS does not do iterative<br>
*> refinement because doing so cannot improve the backward error.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangular matrix A, packed columnwise in<br>
*>          a linear array.  The j-th column of A is stored in the array<br>
*>          AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          If DIAG = 'U', the diagonal elements of A are not referenced<br>
*>          and are assumed to be 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
*>          The right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
*>          The solution matrix X.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztprfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,double[] AP,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTPTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztptri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztptri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztptri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPTRI( UPLO, DIAG, N, AP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPTRI computes the inverse of a complex upper or lower triangular<br>
*> matrix A stored in packed format.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangular matrix A, stored<br>
*>          columnwise in a linear array.  The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          See below for further details.<br>
*>          On exit, the (triangular) inverse of the original matrix, in<br>
*>          the same packed storage format.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular<br>
*>                matrix is singular and its inverse can not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  A triangular matrix A can be transferred to packed storage using one<br>
*>  of the following program segments:<br>
*><br>
*>  UPLO = 'U':                      UPLO = 'L':<br>
*><br>
*>        JC = 1                           JC = 1<br>
*>        DO 2 J = 1, N                    DO 2 J = 1, N<br>
*>           DO 1 I = 1, J                    DO 1 I = J, N<br>
*>              AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J)<br>
*>      1    CONTINUE                    1    CONTINUE<br>
*>           JC = JC + J                      JC = JC + N - J + 1<br>
*>      2 CONTINUE                       2 CONTINUE<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztptri_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] AP,INTEGER INFO);
/**
*> \brief \b ZTPTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztptrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztptrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztptrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         AP( * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPTRS solves a triangular system of the form<br>
*><br>
*>    A * X = B,  A**T * X = B,  or  A**H * X = B,<br>
*><br>
*> where A is a triangular matrix of order N stored in packed format,<br>
*> and B is an N-by-NRHS matrix.  A check is made to verify that A is<br>
*> nonsingular.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangular matrix A, packed columnwise in<br>
*>          a linear array.  The j-th column of A is stored in the array<br>
*>          AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
*>          On entry, the right hand side matrix B.<br>
*>          On exit, if INFO = 0, the solution matrix X.<br>
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
*>          > 0:  if INFO = i, the i-th diagonal element of A is zero,<br>
*>                indicating that the matrix is singular and the<br>
*>                solutions X have not been computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztptrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,double[] AP,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZTPTTF copies a triangular matrix from the standard packed format (TP) to the rectangular full packed format (TF).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPTTF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpttf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpttf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpttf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPTTF( TRANSR, UPLO, N, AP, ARF, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         AP( 0: * ), ARF( 0: * )<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPTTF copies a triangular matrix A from standard packed format (TP)<br>
*> to rectangular full packed format (TF).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  ARF in Normal format is wanted;<br>
*>          = 'C':  ARF in Conjugate-transpose format is wanted.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
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
*>          AP is COMPLEX*16 array, dimension ( N*(N+1)/2 ),<br>
*>          On entry, the upper or lower triangular matrix A, packed<br>
*>          columnwise in a linear array. The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ARF<br>
*> \verbatim<br>
*>          ARF is COMPLEX*16 array, dimension ( N*(N+1)/2 ),<br>
*>          On exit, the upper or lower triangular matrix A stored in<br>
*>          RFP format. For a further discussion see Notes below.<br>
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
*> \ingroup complex16OTHERcomputational<br>
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
	public void ztpttf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] AP,double[] ARF,INTEGER INFO);
/**
*> \brief \b ZTPTTR copies a triangular matrix from the standard packed format (TP) to the standard full format (TR).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTPTTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpttr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpttr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpttr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPTTR( UPLO, N, AP, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPTTR copies a triangular matrix A from standard packed format (TP)<br>
*> to standard full format (TR).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular.<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension ( N*(N+1)/2 ),<br>
*>          On entry, the upper or lower triangular matrix A, packed<br>
*>          columnwise in a linear array. The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension ( LDA, N )<br>
*>          On exit, the triangular matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
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
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztpttr_(CHARACTER UPLO,INTEGER N,double[] AP,double[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b ZTRCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK,<br>
*                          RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRCON estimates the reciprocal of the condition number of a<br>
*> triangular matrix A, in either the 1-norm or the infinity-norm.<br>
*><br>
*> The norm of A is computed and an estimate is obtained for<br>
*> norm(inv(A)), then the reciprocal of the condition number is<br>
*> computed as<br>
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
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The triangular matrix A.  If UPLO = 'U', the leading N-by-N<br>
*>          upper triangular part of the array A contains the upper<br>
*>          triangular matrix, and the strictly lower triangular part of<br>
*>          A is not referenced.  If UPLO = 'L', the leading N-by-N lower<br>
*>          triangular part of the array A contains the lower triangular<br>
*>          matrix, and the strictly upper triangular part of A is not<br>
*>          referenced.  If DIAG = 'U', the diagonal elements of A are<br>
*>          also not referenced and are assumed to be 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztrcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTREVC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTREVC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrevc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrevc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrevc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,<br>
*                          LDVR, MM, M, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, SIDE<br>
*       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTREVC computes some or all of the right and/or left eigenvectors of<br>
*> a complex upper triangular matrix T.<br>
*> Matrices of this type are produced by the Schur factorization of<br>
*> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.<br>
*> <br>
*> The right eigenvector x and the left eigenvector y of T corresponding<br>
*> to an eigenvalue w are defined by:<br>
*> <br>
*>              T*x = w*x,     (y**H)*T = w*(y**H)<br>
*> <br>
*> where y**H denotes the conjugate transpose of the vector y.<br>
*> The eigenvalues are not input to this routine, but are read directly<br>
*> from the diagonal of T.<br>
*> <br>
*> This routine returns the matrices X and/or Y of right and left<br>
*> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an<br>
*> input matrix.  If Q is the unitary factor that reduces a matrix A to<br>
*> Schur form T, then Q*X and Q*Y are the matrices of right and left<br>
*> eigenvectors of A.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'R':  compute right eigenvectors only;<br>
*>          = 'L':  compute left eigenvectors only;<br>
*>          = 'B':  compute both right and left eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] HOWMNY<br>
*> \verbatim<br>
*>          HOWMNY is CHARACTER*1<br>
*>          = 'A':  compute all right and/or left eigenvectors;<br>
*>          = 'B':  compute all right and/or left eigenvectors,<br>
*>                  backtransformed using the matrices supplied in<br>
*>                  VR and/or VL;<br>
*>          = 'S':  compute selected right and/or left eigenvectors,<br>
*>                  as indicated by the logical array SELECT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be<br>
*>          computed.<br>
*>          The eigenvector corresponding to the j-th eigenvalue is<br>
*>          computed if SELECT(j) = .TRUE..<br>
*>          Not referenced if HOWMNY = 'A' or 'B'.<br>
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
*>          T is COMPLEX*16 array, dimension (LDT,N)<br>
*>          The upper triangular matrix T.  T is modified, but restored<br>
*>          on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. LDT >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is COMPLEX*16 array, dimension (LDVL,MM)<br>
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must<br>
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of<br>
*>          Schur vectors returned by ZHSEQR).<br>
*>          On exit, if SIDE = 'L' or 'B', VL contains:<br>
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;<br>
*>          if HOWMNY = 'B', the matrix Q*Y;<br>
*>          if HOWMNY = 'S', the left eigenvectors of T specified by<br>
*>                           SELECT, stored consecutively in the columns<br>
*>                           of VL, in the same order as their<br>
*>                           eigenvalues.<br>
*>          Not referenced if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL.  LDVL >= 1, and if<br>
*>          SIDE = 'L' or 'B', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VR<br>
*> \verbatim<br>
*>          VR is COMPLEX*16 array, dimension (LDVR,MM)<br>
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must<br>
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of<br>
*>          Schur vectors returned by ZHSEQR).<br>
*>          On exit, if SIDE = 'R' or 'B', VR contains:<br>
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;<br>
*>          if HOWMNY = 'B', the matrix Q*X;<br>
*>          if HOWMNY = 'S', the right eigenvectors of T specified by<br>
*>                           SELECT, stored consecutively in the columns<br>
*>                           of VR, in the same order as their<br>
*>                           eigenvalues.<br>
*>          Not referenced if SIDE = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR.  LDVR >= 1, and if<br>
*>          SIDE = 'R' or 'B'; LDVR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MM<br>
*> \verbatim<br>
*>          MM is INTEGER<br>
*>          The number of columns in the arrays VL and/or VR. MM >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of columns in the arrays VL and/or VR actually<br>
*>          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M<br>
*>          is set to N.  Each selected eigenvector occupies one<br>
*>          column.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The algorithm used in this program is basically backward (forward)<br>
*>  substitution, with scaling to make the the code robust against<br>
*>  possible overflow.<br>
*><br>
*>  Each eigenvector is normalized so that the element of largest<br>
*>  magnitude has magnitude 1; here the magnitude of a complex number<br>
*>  (x,y) is taken to be |x| + |y|.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztrevc_(CHARACTER SIDE,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,double[] T,INTEGER LDT,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTREVC3<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTREVC3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrevc3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrevc3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrevc3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL,<br>
*                           VR, LDVR, MM, M, WORK, LWORK, RWORK, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, SIDE<br>
*       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTREVC3 computes some or all of the right and/or left eigenvectors of<br>
*> a complex upper triangular matrix T.<br>
*> Matrices of this type are produced by the Schur factorization of<br>
*> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.<br>
*><br>
*> The right eigenvector x and the left eigenvector y of T corresponding<br>
*> to an eigenvalue w are defined by:<br>
*><br>
*>              T*x = w*x,     (y**H)*T = w*(y**H)<br>
*><br>
*> where y**H denotes the conjugate transpose of the vector y.<br>
*> The eigenvalues are not input to this routine, but are read directly<br>
*> from the diagonal of T.<br>
*><br>
*> This routine returns the matrices X and/or Y of right and left<br>
*> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an<br>
*> input matrix. If Q is the unitary factor that reduces a matrix A to<br>
*> Schur form T, then Q*X and Q*Y are the matrices of right and left<br>
*> eigenvectors of A.<br>
*><br>
*> This uses a Level 3 BLAS version of the back transformation.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'R':  compute right eigenvectors only;<br>
*>          = 'L':  compute left eigenvectors only;<br>
*>          = 'B':  compute both right and left eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] HOWMNY<br>
*> \verbatim<br>
*>          HOWMNY is CHARACTER*1<br>
*>          = 'A':  compute all right and/or left eigenvectors;<br>
*>          = 'B':  compute all right and/or left eigenvectors,<br>
*>                  backtransformed using the matrices supplied in<br>
*>                  VR and/or VL;<br>
*>          = 'S':  compute selected right and/or left eigenvectors,<br>
*>                  as indicated by the logical array SELECT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be<br>
*>          computed.<br>
*>          The eigenvector corresponding to the j-th eigenvalue is<br>
*>          computed if SELECT(j) = .TRUE..<br>
*>          Not referenced if HOWMNY = 'A' or 'B'.<br>
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
*>          T is COMPLEX*16 array, dimension (LDT,N)<br>
*>          The upper triangular matrix T.  T is modified, but restored<br>
*>          on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. LDT >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is COMPLEX*16 array, dimension (LDVL,MM)<br>
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must<br>
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of<br>
*>          Schur vectors returned by ZHSEQR).<br>
*>          On exit, if SIDE = 'L' or 'B', VL contains:<br>
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;<br>
*>          if HOWMNY = 'B', the matrix Q*Y;<br>
*>          if HOWMNY = 'S', the left eigenvectors of T specified by<br>
*>                           SELECT, stored consecutively in the columns<br>
*>                           of VL, in the same order as their<br>
*>                           eigenvalues.<br>
*>          Not referenced if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL.<br>
*>          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VR<br>
*> \verbatim<br>
*>          VR is COMPLEX*16 array, dimension (LDVR,MM)<br>
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must<br>
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of<br>
*>          Schur vectors returned by ZHSEQR).<br>
*>          On exit, if SIDE = 'R' or 'B', VR contains:<br>
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;<br>
*>          if HOWMNY = 'B', the matrix Q*X;<br>
*>          if HOWMNY = 'S', the right eigenvectors of T specified by<br>
*>                           SELECT, stored consecutively in the columns<br>
*>                           of VR, in the same order as their<br>
*>                           eigenvalues.<br>
*>          Not referenced if SIDE = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR.<br>
*>          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MM<br>
*> \verbatim<br>
*>          MM is INTEGER<br>
*>          The number of columns in the arrays VL and/or VR. MM >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of columns in the arrays VL and/or VR actually<br>
*>          used to store the eigenvectors.<br>
*>          If HOWMNY = 'A' or 'B', M is set to N.<br>
*>          Each selected eigenvector occupies one column.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of array WORK. LWORK >= max(1,2*N).<br>
*>          For optimum performance, LWORK >= N + 2*N*NB, where NB is<br>
*>          the optimal blocksize.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (LRWORK)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LRWORK<br>
*> \verbatim<br>
*>          LRWORK is INTEGER<br>
*>          The dimension of array RWORK. LRWORK >= max(1,N).<br>
*><br>
*>          If LRWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the RWORK array, returns<br>
*>          this value as the first entry of the RWORK array, and no error<br>
*>          message related to LRWORK is issued by XERBLA.<br>
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
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2011<br>
*<br>
*  @precisions fortran z -> c<br>
*<br>
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The algorithm used in this program is basically backward (forward)<br>
*>  substitution, with scaling to make the the code robust against<br>
*>  possible overflow.<br>
*><br>
*>  Each eigenvector is normalized so that the element of largest<br>
*>  magnitude has magnitude 1; here the magnitude of a complex number<br>
*>  (x,y) is taken to be |x| + |y|.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztrevc3_(CHARACTER SIDE,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,double[] T,INTEGER LDT,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,INTEGER INFO);
/**
*> \brief \b ZTREXC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTREXC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrexc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrexc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrexc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ<br>
*       INTEGER            IFST, ILST, INFO, LDQ, LDT, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         Q( LDQ, * ), T( LDT, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTREXC reorders the Schur factorization of a complex matrix<br>
*> A = Q*T*Q**H, so that the diagonal element of T with row index IFST<br>
*> is moved to row ILST.<br>
*><br>
*> The Schur form T is reordered by a unitary similarity transformation<br>
*> Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by<br>
*> postmultplying it with Z.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] COMPQ<br>
*> \verbatim<br>
*>          COMPQ is CHARACTER*1<br>
*>          = 'V':  update the matrix Q of Schur vectors;<br>
*>          = 'N':  do not update Q.<br>
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
*>          T is COMPLEX*16 array, dimension (LDT,N)<br>
*>          On entry, the upper triangular matrix T.<br>
*>          On exit, the reordered upper triangular matrix.<br>
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
*>          Q is COMPLEX*16 array, dimension (LDQ,N)<br>
*>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.<br>
*>          On exit, if COMPQ = 'V', Q has been postmultiplied by the<br>
*>          unitary transformation matrix Z which reorders T.<br>
*>          If COMPQ = 'N', Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IFST<br>
*> \verbatim<br>
*>          IFST is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILST<br>
*> \verbatim<br>
*>          ILST is INTEGER<br>
*><br>
*>          Specify the reordering of the diagonal elements of T:<br>
*>          The element with row index IFST is moved to row ILST by a<br>
*>          sequence of transpositions between adjacent elements.<br>
*>          1 <= IFST <= N; 1 <= ILST <= N.<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztrexc_(CHARACTER COMPQ,INTEGER N,double[] T,INTEGER LDT,double[] Q,INTEGER LDQ,INTEGER IFST,INTEGER ILST,INTEGER INFO);
/**
*> \brief \b ZTRRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,<br>
*                          LDX, FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, LDA, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRRFS provides error bounds and backward error estimates for the<br>
*> solution to a system of linear equations with a triangular<br>
*> coefficient matrix.<br>
*><br>
*> The solution matrix X must be computed by ZTRTRS or some other<br>
*> means before entering this routine.  ZTRRFS does not do iterative<br>
*> refinement because doing so cannot improve the backward error.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The triangular matrix A.  If UPLO = 'U', the leading N-by-N<br>
*>          upper triangular part of the array A contains the upper<br>
*>          triangular matrix, and the strictly lower triangular part of<br>
*>          A is not referenced.  If UPLO = 'L', the leading N-by-N lower<br>
*>          triangular part of the array A contains the lower triangular<br>
*>          matrix, and the strictly upper triangular part of A is not<br>
*>          referenced.  If DIAG = 'U', the diagonal elements of A are<br>
*>          also not referenced and are assumed to be 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
*>          The right hand side matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
*>          The solution matrix X.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztrrfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTRSEN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRSEN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsen.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsen.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsen.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S,<br>
*                          SEP, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ, JOB<br>
*       INTEGER            INFO, LDQ, LDT, LWORK, M, N<br>
*       DOUBLE PRECISION   S, SEP<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       COMPLEX*16         Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRSEN reorders the Schur factorization of a complex matrix<br>
*> A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in<br>
*> the leading positions on the diagonal of the upper triangular matrix<br>
*> T, and the leading columns of Q form an orthonormal basis of the<br>
*> corresponding right invariant subspace.<br>
*><br>
*> Optionally the routine computes the reciprocal condition numbers of<br>
*> the cluster of eigenvalues and/or the invariant subspace.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          Specifies whether condition numbers are required for the<br>
*>          cluster of eigenvalues (S) or the invariant subspace (SEP):<br>
*>          = 'N': none;<br>
*>          = 'E': for eigenvalues only (S);<br>
*>          = 'V': for invariant subspace only (SEP);<br>
*>          = 'B': for both eigenvalues and invariant subspace (S and<br>
*>                 SEP).<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPQ<br>
*> \verbatim<br>
*>          COMPQ is CHARACTER*1<br>
*>          = 'V': update the matrix Q of Schur vectors;<br>
*>          = 'N': do not update Q.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          SELECT specifies the eigenvalues in the selected cluster. To<br>
*>          select the j-th eigenvalue, SELECT(j) must be set to .TRUE..<br>
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
*>          T is COMPLEX*16 array, dimension (LDT,N)<br>
*>          On entry, the upper triangular matrix T.<br>
*>          On exit, T is overwritten by the reordered matrix T, with the<br>
*>          selected eigenvalues as the leading diagonal elements.<br>
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
*>          Q is COMPLEX*16 array, dimension (LDQ,N)<br>
*>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.<br>
*>          On exit, if COMPQ = 'V', Q has been postmultiplied by the<br>
*>          unitary transformation matrix which reorders T; the leading M<br>
*>          columns of Q form an orthonormal basis for the specified<br>
*>          invariant subspace.<br>
*>          If COMPQ = 'N', Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.<br>
*>          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX*16 array, dimension (N)<br>
*>          The reordered eigenvalues of T, in the same order as they<br>
*>          appear on the diagonal of T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The dimension of the specified invariant subspace.<br>
*>          0 <= M <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION<br>
*>          If JOB = 'E' or 'B', S is a lower bound on the reciprocal<br>
*>          condition number for the selected cluster of eigenvalues.<br>
*>          S cannot underestimate the true reciprocal condition number<br>
*>          by more than a factor of sqrt(N). If M = 0 or N, S = 1.<br>
*>          If JOB = 'N' or 'V', S is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SEP<br>
*> \verbatim<br>
*>          SEP is DOUBLE PRECISION<br>
*>          If JOB = 'V' or 'B', SEP is the estimated reciprocal<br>
*>          condition number of the specified invariant subspace. If<br>
*>          M = 0 or N, SEP = norm(T).<br>
*>          If JOB = 'N' or 'E', SEP is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If JOB = 'N', LWORK >= 1;<br>
*>          if JOB = 'E', LWORK = max(1,M*(N-M));<br>
*>          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  ZTRSEN first collects the selected eigenvalues by computing a unitary<br>
*>  transformation Z to move them to the top left corner of T. In other<br>
*>  words, the selected eigenvalues are the eigenvalues of T11 in:<br>
*><br>
*>          Z**H * T * Z = ( T11 T12 ) n1<br>
*>                         (  0  T22 ) n2<br>
*>                            n1  n2<br>
*><br>
*>  where N = n1+n2. The first<br>
*>  n1 columns of Z span the specified invariant subspace of T.<br>
*><br>
*>  If T has been obtained from the Schur factorization of a matrix<br>
*>  A = Q*T*Q**H, then the reordered Schur factorization of A is given by<br>
*>  A = (Q*Z)*(Z**H*T*Z)*(Q*Z)**H, and the first n1 columns of Q*Z span the<br>
*>  corresponding invariant subspace of A.<br>
*><br>
*>  The reciprocal condition number of the average of the eigenvalues of<br>
*>  T11 may be returned in S. S lies between 0 (very badly conditioned)<br>
*>  and 1 (very well conditioned). It is computed as follows. First we<br>
*>  compute R so that<br>
*><br>
*>                         P = ( I  R ) n1<br>
*>                             ( 0  0 ) n2<br>
*>                               n1 n2<br>
*><br>
*>  is the projector on the invariant subspace associated with T11.<br>
*>  R is the solution of the Sylvester equation:<br>
*><br>
*>                        T11*R - R*T22 = T12.<br>
*><br>
*>  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote<br>
*>  the two-norm of M. Then S is computed as the lower bound<br>
*><br>
*>                      (1 + F-norm(R)**2)**(-1/2)<br>
*><br>
*>  on the reciprocal of 2-norm(P), the true reciprocal condition number.<br>
*>  S cannot underestimate 1 / 2-norm(P) by more than a factor of<br>
*>  sqrt(N).<br>
*><br>
*>  An approximate error bound for the computed average of the<br>
*>  eigenvalues of T11 is<br>
*><br>
*>                         EPS * norm(T) / S<br>
*><br>
*>  where EPS is the machine precision.<br>
*><br>
*>  The reciprocal condition number of the right invariant subspace<br>
*>  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.<br>
*>  SEP is defined as the separation of T11 and T22:<br>
*><br>
*>                     sep( T11, T22 ) = sigma-min( C )<br>
*><br>
*>  where sigma-min(C) is the smallest singular value of the<br>
*>  n1*n2-by-n1*n2 matrix<br>
*><br>
*>     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )<br>
*><br>
*>  I(m) is an m by m identity matrix, and kprod denotes the Kronecker<br>
*>  product. We estimate sigma-min(C) by the reciprocal of an estimate of<br>
*>  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)<br>
*>  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).<br>
*><br>
*>  When SEP is small, small changes in T can cause large changes in<br>
*>  the invariant subspace. An approximate bound on the maximum angular<br>
*>  error in the computed right invariant subspace is<br>
*><br>
*>                      EPS * norm(T) / SEP<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztrsen_(CHARACTER JOB,CHARACTER COMPQ,boolean[] SELECT,INTEGER N,double[] T,INTEGER LDT,double[] Q,INTEGER LDQ,double[] W,INTEGER M,DOUBLE S,DOUBLE SEP,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZTRSNA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRSNA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsna.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsna.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsna.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,<br>
*                          LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, JOB<br>
*       INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       DOUBLE PRECISION   RWORK( * ), S( * ), SEP( * )<br>
*       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   WORK( LDWORK, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRSNA estimates reciprocal condition numbers for specified<br>
*> eigenvalues and/or right eigenvectors of a complex upper triangular<br>
*> matrix T (or of any matrix Q*T*Q**H with Q unitary).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          Specifies whether condition numbers are required for<br>
*>          eigenvalues (S) or eigenvectors (SEP):<br>
*>          = 'E': for eigenvalues only (S);<br>
*>          = 'V': for eigenvectors only (SEP);<br>
*>          = 'B': for both eigenvalues and eigenvectors (S and SEP).<br>
*> \endverbatim<br>
*><br>
*> \param[in] HOWMNY<br>
*> \verbatim<br>
*>          HOWMNY is CHARACTER*1<br>
*>          = 'A': compute condition numbers for all eigenpairs;<br>
*>          = 'S': compute condition numbers for selected eigenpairs<br>
*>                 specified by the array SELECT.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          If HOWMNY = 'S', SELECT specifies the eigenpairs for which<br>
*>          condition numbers are required. To select condition numbers<br>
*>          for the j-th eigenpair, SELECT(j) must be set to .TRUE..<br>
*>          If HOWMNY = 'A', SELECT is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix T. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] T<br>
*> \verbatim<br>
*>          T is COMPLEX*16 array, dimension (LDT,N)<br>
*>          The upper triangular matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T. LDT >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is COMPLEX*16 array, dimension (LDVL,M)<br>
*>          If JOB = 'E' or 'B', VL must contain left eigenvectors of T<br>
*>          (or of any Q*T*Q**H with Q unitary), corresponding to the<br>
*>          eigenpairs specified by HOWMNY and SELECT. The eigenvectors<br>
*>          must be stored in consecutive columns of VL, as returned by<br>
*>          ZHSEIN or ZTREVC.<br>
*>          If JOB = 'V', VL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL.<br>
*>          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VR<br>
*> \verbatim<br>
*>          VR is COMPLEX*16 array, dimension (LDVR,M)<br>
*>          If JOB = 'E' or 'B', VR must contain right eigenvectors of T<br>
*>          (or of any Q*T*Q**H with Q unitary), corresponding to the<br>
*>          eigenpairs specified by HOWMNY and SELECT. The eigenvectors<br>
*>          must be stored in consecutive columns of VR, as returned by<br>
*>          ZHSEIN or ZTREVC.<br>
*>          If JOB = 'V', VR is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR.<br>
*>          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (MM)<br>
*>          If JOB = 'E' or 'B', the reciprocal condition numbers of the<br>
*>          selected eigenvalues, stored in consecutive elements of the<br>
*>          array. Thus S(j), SEP(j), and the j-th columns of VL and VR<br>
*>          all correspond to the same eigenpair (but not in general the<br>
*>          j-th eigenpair, unless all eigenpairs are selected).<br>
*>          If JOB = 'V', S is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SEP<br>
*> \verbatim<br>
*>          SEP is DOUBLE PRECISION array, dimension (MM)<br>
*>          If JOB = 'V' or 'B', the estimated reciprocal condition<br>
*>          numbers of the selected eigenvectors, stored in consecutive<br>
*>          elements of the array.<br>
*>          If JOB = 'E', SEP is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MM<br>
*> \verbatim<br>
*>          MM is INTEGER<br>
*>          The number of elements in the arrays S (if JOB = 'E' or 'B')<br>
*>           and/or SEP (if JOB = 'V' or 'B'). MM >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of elements of the arrays S and/or SEP actually<br>
*>          used to store the estimated condition numbers.<br>
*>          If HOWMNY = 'A', M is set to N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (LDWORK,N+6)<br>
*>          If JOB = 'E', WORK is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDWORK<br>
*> \verbatim<br>
*>          LDWORK is INTEGER<br>
*>          The leading dimension of the array WORK.<br>
*>          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
*>          If JOB = 'E', RWORK is not referenced.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The reciprocal of the condition number of an eigenvalue lambda is<br>
*>  defined as<br>
*><br>
*>          S(lambda) = |v**H*u| / (norm(u)*norm(v))<br>
*><br>
*>  where u and v are the right and left eigenvectors of T corresponding<br>
*>  to lambda; v**H denotes the conjugate transpose of v, and norm(u)<br>
*>  denotes the Euclidean norm. These reciprocal condition numbers always<br>
*>  lie between zero (very badly conditioned) and one (very well<br>
*>  conditioned). If n = 1, S(lambda) is defined to be 1.<br>
*><br>
*>  An approximate error bound for a computed eigenvalue W(i) is given by<br>
*><br>
*>                      EPS * norm(T) / S(i)<br>
*><br>
*>  where EPS is the machine precision.<br>
*><br>
*>  The reciprocal of the condition number of the right eigenvector u<br>
*>  corresponding to lambda is defined as follows. Suppose<br>
*><br>
*>              T = ( lambda  c  )<br>
*>                  (   0    T22 )<br>
*><br>
*>  Then the reciprocal condition number is<br>
*><br>
*>          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )<br>
*><br>
*>  where sigma-min denotes the smallest singular value. We approximate<br>
*>  the smallest singular value by the reciprocal of an estimate of the<br>
*>  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is<br>
*>  defined to be abs(T(1,1)).<br>
*><br>
*>  An approximate error bound for a computed right eigenvector VR(i)<br>
*>  is given by<br>
*><br>
*>                      EPS * norm(T) / SEP(i)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztrsna_(CHARACTER JOB,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,double[] T,INTEGER LDT,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] S,double[] SEP,INTEGER MM,INTEGER M,double[] WORK,INTEGER LDWORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZTRSYL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRSYL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsyl.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsyl.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsyl.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,<br>
*                          LDC, SCALE, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANA, TRANB<br>
*       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N<br>
*       DOUBLE PRECISION   SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRSYL solves the complex Sylvester matrix equation:<br>
*><br>
*>    op(A)*X + X*op(B) = scale*C or<br>
*>    op(A)*X - X*op(B) = scale*C,<br>
*><br>
*> where op(A) = A or A**H, and A and B are both upper triangular. A is<br>
*> M-by-M and B is N-by-N; the right hand side C and the solution X are<br>
*> M-by-N; and scale is an output scale factor, set <= 1 to avoid<br>
*> overflow in X.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANA<br>
*> \verbatim<br>
*>          TRANA is CHARACTER*1<br>
*>          Specifies the option op(A):<br>
*>          = 'N': op(A) = A    (No transpose)<br>
*>          = 'C': op(A) = A**H (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANB<br>
*> \verbatim<br>
*>          TRANB is CHARACTER*1<br>
*>          Specifies the option op(B):<br>
*>          = 'N': op(B) = B    (No transpose)<br>
*>          = 'C': op(B) = B**H (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] ISGN<br>
*> \verbatim<br>
*>          ISGN is INTEGER<br>
*>          Specifies the sign in the equation:<br>
*>          = +1: solve op(A)*X + X*op(B) = scale*C<br>
*>          = -1: solve op(A)*X - X*op(B) = scale*C<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The order of the matrix A, and the number of rows in the<br>
*>          matrices X and C. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix B, and the number of columns in the<br>
*>          matrices X and C. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,M)<br>
*>          The upper triangular matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,N)<br>
*>          The upper triangular matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array, dimension (LDC,N)<br>
*>          On entry, the M-by-N right hand side matrix C.<br>
*>          On exit, C is overwritten by the solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>          The leading dimension of the array C. LDC >= max(1,M)<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is DOUBLE PRECISION<br>
*>          The scale factor, scale, set <= 1 to avoid overflow in X.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          = 1: A and B have common or very close eigenvalues; perturbed<br>
*>               values were used to solve the equation (but the matrices<br>
*>               A and B are unchanged).<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztrsyl_(CHARACTER TRANA,CHARACTER TRANB,INTEGER ISGN,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] C,INTEGER LDC,DOUBLE SCALE,INTEGER INFO);
/**
*> \brief \b ZTRTI2 computes the inverse of a triangular matrix (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRTI2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrti2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrti2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrti2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRTI2 computes the inverse of a complex upper or lower triangular<br>
*> matrix.<br>
*><br>
*> This is the Level 2 BLAS version of the algorithm.<br>
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
*>          The order of the matrix A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the triangular matrix A.  If UPLO = 'U', the<br>
*>          leading n by n upper triangular part of the array A contains<br>
*>          the upper triangular matrix, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n by n lower triangular part of the array A contains<br>
*>          the lower triangular matrix, and the strictly upper<br>
*>          triangular part of A is not referenced.  If DIAG = 'U', the<br>
*>          diagonal elements of A are also not referenced and are<br>
*>          assumed to be 1.<br>
*><br>
*>          On exit, the (triangular) inverse of the original matrix, in<br>
*>          the same storage format.<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztrti2_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b ZTRTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrtri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrtri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrtri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRTRI computes the inverse of a complex upper or lower triangular<br>
*> matrix A.<br>
*><br>
*> This is the Level 3 BLAS version of the algorithm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the triangular matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of the array A contains<br>
*>          the upper triangular matrix, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of the array A contains<br>
*>          the lower triangular matrix, and the strictly upper<br>
*>          triangular part of A is not referenced.  If DIAG = 'U', the<br>
*>          diagonal elements of A are also not referenced and are<br>
*>          assumed to be 1.<br>
*>          On exit, the (triangular) inverse of the original matrix, in<br>
*>          the same storage format.<br>
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
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular<br>
*>               matrix is singular and its inverse can not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztrtri_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b ZTRTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrtrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrtrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrtrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRTRS solves a triangular system of the form<br>
*><br>
*>    A * X = B,  A**T * X = B,  or  A**H * X = B,<br>
*><br>
*> where A is a triangular matrix of order N, and B is an N-by-NRHS<br>
*> matrix.  A check is made to verify that A is nonsingular.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          Specifies the form of the system of equations:<br>
*>          = 'N':  A * X = B     (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>          = 'N':  A is non-unit triangular;<br>
*>          = 'U':  A is unit triangular.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The triangular matrix A.  If UPLO = 'U', the leading N-by-N<br>
*>          upper triangular part of the array A contains the upper<br>
*>          triangular matrix, and the strictly lower triangular part of<br>
*>          A is not referenced.  If UPLO = 'L', the leading N-by-N lower<br>
*>          triangular part of the array A contains the lower triangular<br>
*>          matrix, and the strictly upper triangular part of A is not<br>
*>          referenced.  If DIAG = 'U', the diagonal elements of A are<br>
*>          also not referenced and are assumed to be 1.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
*>          On entry, the right hand side matrix B.<br>
*>          On exit, if INFO = 0, the solution matrix X.<br>
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
*>          > 0: if INFO = i, the i-th diagonal element of A is zero,<br>
*>               indicating that the matrix is singular and the solutions<br>
*>               X have not been computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztrtrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZTRTTF copies a triangular matrix from the standard full format (TR) to the rectangular full packed format (TF).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRTTF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrttf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrttf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrttf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRTTF( TRANSR, UPLO, N, A, LDA, ARF, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( 0: LDA-1, 0: * ), ARF( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRTTF copies a triangular matrix A from standard full format (TR)<br>
*> to rectangular full packed format (TF) .<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  ARF in Normal mode is wanted;<br>
*>          = 'C':  ARF in Conjugate Transpose mode is wanted;<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
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
*>          A is COMPLEX*16 array, dimension ( LDA, N )<br>
*>          On entry, the triangular matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of the array A contains<br>
*>          the upper triangular matrix, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of the array A contains<br>
*>          the lower triangular matrix, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the matrix A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ARF<br>
*> \verbatim<br>
*>          ARF is COMPLEX*16 array, dimension ( N*(N+1)/2 ),<br>
*>          On exit, the upper or lower triangular matrix A stored in<br>
*>          RFP format. For a further discussion see Notes below.<br>
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
*> \ingroup complex16OTHERcomputational<br>
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
	public void ztrttf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] ARF,INTEGER INFO);
/**
*> \brief \b ZTRTTP copies a triangular matrix from the standard full format (TR) to the standard packed format (TP).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTRTTP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrttp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrttp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrttp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRTTP( UPLO, N, A, LDA, AP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRTTP copies a triangular matrix A from full format (TR) to standard<br>
*> packed format (TP).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  A is upper triangular;<br>
*>          = 'L':  A is lower triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices AP and A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the triangular matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension ( N*(N+1)/2 ),<br>
*>          On exit, the upper or lower triangular matrix A, packed<br>
*>          columnwise in a linear array. The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ztrttp_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] AP,INTEGER INFO);
/**
*> \brief \b ZTZRZF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZTZRZF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztzrzf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztzrzf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztzrzf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTZRZF reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A<br>
*> to upper triangular form by means of unitary transformations.<br>
*><br>
*> The upper trapezoidal matrix A is factored as<br>
*><br>
*>    A = ( R  0 ) * Z,<br>
*><br>
*> where Z is an N-by-N unitary matrix and R is an M-by-M upper<br>
*> triangular matrix.<br>
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
*>          The number of columns of the matrix A.  N >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the leading M-by-N upper trapezoidal part of the<br>
*>          array A must contain the matrix to be factorized.<br>
*>          On exit, the leading M-by-M upper triangular part of A<br>
*>          contains the upper triangular matrix R, and elements M+1 to<br>
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
*>          TAU is COMPLEX*16 array, dimension (M)<br>
*>          The scalar factors of the elementary reflectors.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
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
*> \date April 2012<br>
*<br>
*> \ingroup complex16OTHERcomputational<br>
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
*>  The N-by-N matrix Z can be computed by<br>
*><br>
*>     Z =  Z(1)*Z(2)* ... *Z(M)<br>
*><br>
*>  where each N-by-N Z(k) is given by<br>
*><br>
*>     Z(k) = I - tau(k)*v(k)*v(k)**H<br>
*><br>
*>  with v(k) is the kth row vector of the M-by-N matrix<br>
*><br>
*>     V = ( I   A(:,M+1:N) )<br>
*><br>
*>  I is the M-by-M identity matrix, A(:,M+1:N) <br>
*>  is the output stored in A on exit from DTZRZF,<br>
*>  and tau(k) is the kth element of the array TAU.<br>
*><br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztzrzf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);

}